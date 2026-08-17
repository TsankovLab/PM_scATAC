###############################################################
# R2_Q3 : arm-level clone discovery with a MULTIMODAL (any number of components)
# criterion -- no assumption of two populations anywhere in the pipeline.
#
# Why: variance ranks noisy-but-unimodal arms first (P4 8q is variance-rank 12 in
# scATAC, and the sample stays clonal even when 8q is retained). A bimodality
# coefficient fixes the ranking but hard-codes "two populations". Here both the arm
# ranking AND the clone calling use a Gaussian mixture whose NUMBER OF COMPONENTS is
# chosen by BIC (G = 1..GMAX), so an arm with 1, 2, 3 ... copy-number states is scored
# on its own terms.
#
#   arm score  = BIC(best G>=2) - BIC(G=1)   (evidence for more than one component)
#   G_best     = number of components chosen for that arm
#   n_modes    = independent, assumption-free check: count of local maxima of the KDE
#                (Silverman bandwidth) -- does not assume Gaussian components at all
#   filter     = keep the top-K arms by BIC gain (K scanned), per sample and modality
#   clones     = multivariate Mclust on the retained arms, G = 1..GMAX by BIC
#                (silhouette of that partition reported only for comparability)
#
# Input : arm_cnv_zscore/*_arm.rds  (scATAC = per-cell z across arms, no reference
#         subtraction; scRNA = Mesothelium-referenced arm means, unchanged)
#
# Controls: P4 chr8q must be recovered (rank, retention, ARI vs the chr8q-only mixture,
# amp fraction ~0.58-0.60 from the targeted + Visium analyses); P5/P14/P23/P3 must stay
# clonal (G=1) -- a criterion that finds structure everywhere has just found noise.
#
# Outputs: arm_level_multimodal_arms.csv, arm_level_multimodal_scan.csv,
#          arm_level_clone_calls_multimodal.csv, arm_cnv_multimodal/*_arm.rds,
#          Plots/R2_Q3_arm_level_multimodal.pdf
###############################################################
suppressMessages({ library(mclust); library(cluster); library(ComplexHeatmap); library(circlize); library(grid) })
set.seed(1)
SC  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(SC,"git_repo_claude","R2_Q3")); dir.create("Plots", showWarnings=FALSE)
ARGS <- commandArgs(TRUE)
SRC  <- if (length(ARGS) >= 1) ARGS[1] else "arm_cnv_zscore"
KEEP <- if (length(ARGS) >= 2) as.integer(ARGS[2]) else 10
DIR  <- "arm_cnv_multimodal"; dir.create(DIR, showWarnings=FALSE)
GMAX <- 6; KSCAN <- c(5,10,15,20,39); NSUB <- 3000; MAXCELL <- 1200

## ---- 1-D multimodality: mixture with BIC-selected number of components -------------
mm_score <- function(x, gmax=GMAX, nsub=NSUB){
  x <- x[is.finite(x)]; if (length(x) < 50) return(c(gain=NA, G=NA))
  if (length(x) > nsub) x <- x[sort(sample(seq_along(x), nsub))]
  ## E = equal variance across components, V = free variance; BIC picks both the model
  ## and the component count, so an arm with 3 CN states is scored on its own terms
  m <- try(Mclust(x, G=1:gmax, modelNames=c("E","V"), verbose=FALSE), silent=TRUE)
  if (inherits(m,"try-error") || is.null(m)) return(c(gain=NA, G=NA))
  b <- m$BIC; b1 <- suppressWarnings(max(b[1,], na.rm=TRUE))                 # G=1 best
  bm <- suppressWarnings(max(b[-1,], na.rm=TRUE))                            # G>=2 best
  if (!is.finite(b1) || !is.finite(bm)) return(c(gain=NA, G=NA))
  c(gain=as.numeric(bm-b1), G=as.numeric(m$G)) }                             # mclust BIC: higher = better

## ---- assumption-free cross-check: number of KDE modes ----
n_modes <- function(x){ x <- x[is.finite(x)]; if (length(x) < 50) return(NA)
  d <- density(x, bw="SJ", n=512)
  sum(diff(sign(diff(d$y))) == -2 & d$y[2:(length(d$y)-1)] > 0.05*max(d$y)) }

## ---- clone calling: multivariate mixture, number of clones by BIC ----
call_clones_mm <- function(A, gmax=GMAX){
  ok <- apply(A,2,function(z) all(is.finite(z))); A <- A[,ok,drop=FALSE]
  none <- list(G=1, sil=NA, cl=setNames(rep(1L,ncol(A)), colnames(A)))
  if (ncol(A) < 40 || nrow(A) < 2) return(none)
  X <- t(A)
  m <- try(Mclust(X, G=1:gmax, modelNames=c("EII","VII","EEI","VVI"), verbose=FALSE), silent=TRUE)
  if (inherits(m,"try-error") || is.null(m)) return(none)
  cl <- setNames(as.integer(m$classification), rownames(X))
  if (length(unique(cl)) < 2) return(list(G=1, sil=NA, cl=setNames(rep(1L,ncol(A)), colnames(A))))
  ## drop components below the usual robustness floor, then re-label
  tb <- table(cl); keepc <- names(tb)[tb >= max(20, 0.05*ncol(A))]
  if (length(keepc) < 2) return(list(G=1, sil=NA, cl=setNames(rep(1L,ncol(A)), colnames(A))))
  cl[!(as.character(cl) %in% keepc)] <- as.integer(keepc[1])
  cl <- setNames(as.integer(factor(cl)), names(cl))
  d <- dist(t(A - rowMeans(A)))
  list(G=length(unique(cl)), sil=mean(silhouette(cl, d)[,3]), cl=cl) }

ARI <- function(a,b){ tb <- table(a,b); n <- sum(tb)
  ci <- sum(choose(rowSums(tb),2)); cj <- sum(choose(colSums(tb),2)); cij <- sum(choose(tb,2))
  ex <- ci*cj/choose(n,2); (cij-ex)/(0.5*(ci+cj)-ex) }
arm_delta <- function(A, cl){ tb <- sort(table(cl), decreasing=TRUE); if (length(tb)<2) return(NULL)
  rowMeans(A[, cl==names(tb)[1], drop=FALSE]) - rowMeans(A[, cl==names(tb)[2], drop=FALSE]) }

fs <- list.files(SRC, "_arm\\.rds$", full.names=TRUE)
objs <- lapply(fs, readRDS); names(objs) <- sub("_arm\\.rds$","",basename(fs))

## ---------------- STEP 1: score every arm for multimodality ----------------
## One score per (sample, modality, arm). An arm that separates two cell populations has
## a 2+ component distribution; an uninformative arm is unimodal. Nothing here assumes
## exactly two clones -- G is free from 1 to 6.
ARMS <- list()
for (nm in names(objs)){
  o <- objs[[nm]]; A <- o$A
  st <- t(sapply(rownames(A), function(a){ s <- mm_score(A[a,]); c(s, modes=n_modes(A[a,])) }))
  colnames(st) <- c("bic_gain","G_best","n_modes_kde")
  st <- data.frame(sample=o$sample, modality=o$modality, arm=rownames(A), st,
                   row.names=NULL, stringsAsFactors=FALSE)
  st <- st[order(-st$bic_gain),]; st$rank <- seq_len(nrow(st))
  ARMS[[nm]] <- st
  cat(sprintf("%-14s top5: %s\n", nm,
      paste(sprintf("%s(dBIC=%.0f,G=%d,modes=%d)", st$arm[1:5], st$bic_gain[1:5],
                    st$G_best[1:5], st$n_modes_kde[1:5]), collapse=" ")))
}
ARMSall <- do.call(rbind, ARMS); rownames(ARMSall) <- NULL
write.csv(ARMSall, "arm_level_multimodal_arms.csv", row.names=FALSE)

cat("\n=== distribution of the per-arm component count G (how often is an arm NOT unimodal?) ===\n")
print(table(modality=ARMSall$modality, G=ARMSall$G_best))

## ---------------- STEP 2: how many arms to keep? --------------------------
## Cluster on the top-K arms for several K so the operating point is chosen from evidence
## rather than asserted. K=5 recovers P4 best but breaks the negative controls; K>=15
## loses P4; K=10 is the compromise actually used downstream.
scan <- list()
for (nm in names(objs)){
  o <- objs[[nm]]; A <- o$A; st <- ARMS[[nm]]
  for (K in KSCAN){
    keep <- st$arm[1:min(K,nrow(st))]
    cc <- call_clones_mm(A[keep,,drop=FALSE])
    scan[[length(scan)+1]] <- data.frame(sample=o$sample, modality=o$modality, K=K, G=cc$G,
      sil=round(ifelse(is.na(cc$sil),0,cc$sil),2),
      minor=round(ifelse(cc$G>1, min(table(cc$cl))/ncol(A), NA),2),
      has8q=("8q" %in% keep), stringsAsFactors=FALSE) }
}
scan <- do.call(rbind, scan); rownames(scan) <- NULL
write.csv(scan, "arm_level_multimodal_scan.csv", row.names=FALSE)
cat("\n=== K scan: number of clones G chosen by BIC ===\n")
for (md in c("scATAC","scRNA")){ cat("\n--",md,"--\n")
  sub <- scan[scan$modality==md,]
  w <- reshape(sub[,c("sample","K","G")], idvar="sample", timevar="K", direction="wide")
  names(w) <- sub("G\\.","K=",names(w)); print(w, row.names=FALSE) }

## ---------------- STEP 3: primary run at K = KEEP -------------------------
## Re-call clones using only the retained arms, and save one object per sample/modality
## so downstream scripts read a partition rather than recompute it.
cat("\n\n############ primary: top-", KEEP, " arms by BIC gain, clones by BIC ############\n", sep="")
calls <- list(); D39 <- list()
for (nm in names(objs)){
  o <- objs[[nm]]; A <- o$A; st <- ARMS[[nm]]
  keep <- st$arm[1:min(KEEP,nrow(st))]; keep <- rownames(A)[rownames(A) %in% keep]
  Af <- A[keep,,drop=FALSE]
  cc <- call_clones_mm(Af)
  saveRDS(list(A=A, A_filt=Af, arms_kept=keep, clones=cc$cl, k=cc$G, sil=cc$sil,
               delta=arm_delta(Af, cc$cl), sample=o$sample, modality=o$modality),
          file.path(DIR, paste0(nm,"_arm.rds")))
  write.csv(data.frame(cell=names(cc$cl), clone=as.integer(cc$cl)),
            file.path(DIR, paste0(nm,"_clones.csv")), row.names=FALSE)
  if (cc$G > 1) D39[[nm]] <- arm_delta(A, cc$cl)
  calls[[nm]] <- data.frame(sample=o$sample, modality=o$modality, n_cells=ncol(A), G_clones=cc$G,
    sil=round(ifelse(is.na(cc$sil),0,cc$sil),2),
    sizes=paste(as.integer(sort(table(cc$cl), decreasing=TRUE)), collapse="/"),
    arms_kept=paste(keep, collapse=","), stringsAsFactors=FALSE)
}
saveRDS(D39, file.path(DIR,"arm_level_deltas.rds"))
calls <- do.call(rbind, calls); rownames(calls) <- NULL
calls <- calls[order(calls$sample, calls$modality),]
write.csv(calls, "arm_level_clone_calls_multimodal.csv", row.names=FALSE)
print(calls[,c("sample","modality","n_cells","G_clones","sil","sizes")], row.names=FALSE)

## ---------------- STEP 4: controls ---------------------------------------
## POSITIVE: P4 must recover its chr8q clone (rank of 8q, retention, ARI vs the
## chr8q-only split, amp fraction ~0.58-0.60 from the targeted + Visium analyses).
## NEGATIVE: P5 / P14 / P23 / P3 must stay clonal. A criterion that finds subclones
## everywhere has found noise, so both directions are required before believing a filter.
cat("\n=== POSITIVE CONTROL: P4 chr8q ===\n")
for (md in c("scATAC","scRNA")){
  nm <- paste0("P4_",md); if (is.null(objs[[nm]])) next
  o <- objs[[nm]]; st <- ARMS[[nm]]; r <- readRDS(file.path(DIR, paste0(nm,"_arm.rds")))
  i <- which(st$arm=="8q")
  cat(sprintf("%s : 8q dBIC=%.0f rank %d/%d (G=%d, KDE modes=%d), kept=%s | clones G=%d sil=%.2f\n",
      md, st$bic_gain[i], st$rank[i], nrow(st), st$G_best[i], st$n_modes_kde[i],
      "8q" %in% r$arms_kept, r$k, ifelse(is.na(r$sil),0,r$sil)))
  x <- o$A["8q",]; m8 <- Mclust(x, G=2, modelNames="V", verbose=FALSE)   # chr8q-only reference split
  c8 <- m8$classification; hi <- which.max(tapply(x, c8, mean)); frac8 <- mean(c8==hi)
  if (r$k > 1){
    cat(sprintf("   vs chr8q-only mixture: ARI=%.2f ; chr8q-high fraction %.2f (targeted+Visium ~0.58-0.60)\n",
                ARI(r$clones[names(c8)], c8), frac8))
    print(table(clone=r$clones[names(c8)], chr8q=ifelse(c8==hi,"amp","low")))
  } else cat(sprintf("   still clonal; chr8q-only mixture would give amp fraction %.2f\n", frac8))
}
cat("\n=== NEGATIVE CONTROLS (must stay G=1) ===\n")
print(calls[calls$sample %in% c("P5","P14","P23","P3"), c("sample","modality","G_clones","sil","sizes")], row.names=FALSE)

## ---------------- STEP 5: figures ----------------------------------------
col_fun <- colorRamp2(c(-2,-0.6,0,0.6,2), c("#08306b","#6baed6","white","#fb6a4a","#67000d"))
pdf("Plots/R2_Q3_arm_level_multimodal.pdf", width=13, height=7)
for (nm in names(objs)){
  o <- objs[[nm]]; st <- ARMS[[nm]]; r <- readRDS(file.path(DIR, paste0(nm,"_arm.rds")))
  show <- head(st$arm, 5); if (o$sample=="P4" && !("8q" %in% show)) show <- c(show,"8q")
  par(mfrow=c(2,3), mar=c(4,4,3,1), oma=c(0,0,2,0))
  for (a in show){ x <- o$A[a,]; i <- which(st$arm==a); d <- density(x[is.finite(x)], bw="SJ")
    plot(d, main=sprintf("%s  dBIC=%.0f  G=%d  modes=%d", a, st$bic_gain[i], st$G_best[i], st$n_modes_kde[i]),
         xlab="per-cell arm CNV score", lwd=2, col=ifelse(a %in% r$arms_kept,"#d95f0e","#2c7fb8"))
    polygon(d, col=adjustcolor(ifelse(a %in% r$arms_kept,"#d95f0e","#2c7fb8"),0.25), border=NA) }
  mtext(sprintf("%s - most multimodal arms (orange = retained; G = mixture components by BIC)", nm),
        outer=TRUE, font=2, cex=0.95)
}
for (nm in names(objs)){
  r <- readRDS(file.path(DIR, paste0(nm,"_arm.rds"))); Af <- r$A_filt; cl <- r$clones
  cells <- colnames(Af)
  if (length(cells) > MAXCELL){ idx <- sort(sample(seq_along(cells), MAXCELL)); Af <- Af[,idx,drop=FALSE]; cl <- cl[cells[idx]] }
  kk <- length(unique(cl)); rs <- if (kk>1) factor(paste0("clone", cl[colnames(Af)])) else NULL
  Ad <- Af/(sd(Af, na.rm=TRUE)+1e-9)
  draw(Heatmap(t(Ad), name="arm CNV", col=col_fun, cluster_columns=FALSE, cluster_rows=(kk>1),
       show_row_names=FALSE, column_names_gp=gpar(fontsize=9), row_split=rs, row_gap=unit(1.2,"mm"),
       border=TRUE, use_raster=TRUE, raster_quality=2,
       row_title=if(kk>1) NULL else "clonal", row_title_gp=gpar(fontsize=9),
       column_title=sprintf("%s - top-%d multimodal arms - G=%d clones (BIC), sil=%.2f",
                            nm, KEEP, r$k, ifelse(is.na(r$sil),0,r$sil)),
       column_title_gp=gpar(fontface="bold", fontsize=10)), heatmap_legend_side="right")
}
dev.off()
cat("\nDONE\n")
