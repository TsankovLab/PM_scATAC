###############################################################
# R2_Q3 : ARM-LEVEL scATAC CNV *without reference subtraction* -- per-cell Z-SCORE
# ACROSS ARMS instead (same logic as R2_Q3_GL_cnv_heatmaps.R, but per arm).
#
# Two readings of "no reference subtraction", both built here:
#   A "zscore" : GL log2FC (GC-matched-background kept, which is what makes the windows
#                comparable) -> window QC -> arm mean -> PER-CELL Z ACROSS THE 39 ARMS.
#                Replaces the previous per-cell median centring; also removes per-cell
#                amplitude (depth) differences, so clustering sees profile SHAPE.
#   B "noref"  : GL raw fragment counts -> log1p -> per-cell z across windows -> window
#                QC -> arm mean -> per-cell z across arms. NO GC-matched background at
#                all. (Recorded caveat: at bin level this reintroduced the full
#                gene-density/accessibility bias -- 17p/q,19p/q,20q,22q lit up in every
#                sample. Included so the comparison is explicit rather than assumed.)
#
# scRNA arm matrices are NOT changed (only scATAC was asked about): the existing
# Mesothelium-referenced arm matrices are copied into each variant directory so the
# downstream tree / validation / bimodality scripts can run unchanged on each variant.
#
# Outputs per variant V in {zscore, noref}:
#   arm_cnv_<V>/<sample>_scATAC_arm.rds  (+ scRNA copies), arm_cnv_<V>/arm_level_deltas.rds
#   arm_level_clone_calls_<V>.csv, arm_level_top_arms_<V>.csv, arm_level_bad_windows_<V>.csv
#   Plots/R2_Q3_arm_level_heatmaps_<V>.pdf
# Then run R2_Q3_arm_level_tree.R / R2_Q3_arm_level_selected.R with args <dir> <suffix>.
###############################################################
suppressMessages({ library(SummarizedExperiment); library(GenomicRanges); library(cluster)
                   library(ComplexHeatmap); library(circlize); library(grid); library(gridExtra) })
set.seed(1)
SC  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
CND <- file.path(SC,"per_sample_QC_signac","CNV_analysis")
OUT <- file.path(SC,"git_repo_claude","R2_Q3"); setwd(OUT); dir.create("Plots", showWarnings=FALSE)
AUT <- paste0("chr",1:22); MAXCELL <- 1200; MINWIN <- 2
CEN <- c(1,123.6,2,93.1,3,91.1,4,50.4,5,48.3,6,59.2,7,59.5,8,44.9,9,44.4,10,40.6,11,52.7,12,35.9,
         13,17.0,14,17.1,15,19.1,16,36.3,17,23.8,18,17.7,19,25.8,20,28.0,21,11.9,22,14.5)
CENmb <- setNames(CEN[seq(2,length(CEN),2)], paste0("chr",CEN[seq(1,length(CEN),2)]))
ACRO  <- c("13p","14p","15p","21p","22p")
ARMS  <- setdiff(unlist(lapply(1:22, function(i) paste0(i, c("p","q")))), ACRO)
armof <- function(chr,mb) paste0(sub("chr","",chr), ifelse(mb < CENmb[chr],"p","q"))
ATAC_SAMP <- c("P1","P4","P5","P8","P10","P11","P12","P14","P23")

clean_windows <- function(M, wid, arm, tag, clip=2){
  sdw <- apply(M, 1, function(z) mad(z[is.finite(z)])); sdw[!is.finite(sdw)] <- 0
  medw <- apply(M, 1, function(z) median(z[is.finite(z)]))
  thr <- max(1.5, 4*median(sdw[sdw>0]))
  bad <- which(sdw > thr | !is.finite(medw) | abs(medw) > 1.5)
  why <- ifelse(!is.finite(medw[bad]), "nonfinite",
         ifelse(abs(medw[bad]) > 1.5, "extreme_median", "extreme_mad"))
  bd <- if (!length(bad)) data.frame(source=character(0), window=character(0), arm=character(0),
                                     median=numeric(0), mad=numeric(0), threshold=numeric(0),
                                     rule=character(0), stringsAsFactors=FALSE) else
        data.frame(source=tag, window=wid[bad], arm=arm[bad], median=round(medw[bad],2),
                   mad=round(sdw[bad],2), threshold=round(thr,2), rule=why, stringsAsFactors=FALSE)
  if (length(bad)) { M <- M[-bad,,drop=FALSE]; arm <- arm[-bad] }
  M[!is.finite(M)] <- NA; M <- pmax(pmin(M, clip), -clip)
  list(M=M, arm=arm, bad=bd) }

call_clones_arm <- function(A, mode=c("euclid","pearson"), kmax=6){
  mode <- match.arg(mode)
  ok <- apply(A,2,function(z) all(is.finite(z))); A <- A[,ok,drop=FALSE]
  none <- list(k=1, sil=NA, cl=setNames(rep(1L,ncol(A)), colnames(A)))
  if (ncol(A) < 40) return(none)
  if (mode=="euclid"){ Ac <- A - rowMeans(A); d <- dist(t(Ac)) } else { d <- as.dist(1-cor(A)) }
  hc <- hclust(d, "ward.D2"); best <- none
  for (k in 2:min(kmax, ncol(A)-1)){
    cl <- cutree(hc,k); if (length(unique(cl))<k) next
    if (any(table(cl) < max(20, 0.05*ncol(A)))) next
    s <- mean(silhouette(cl,d)[,3])
    if (is.na(best$sil) || s > best$sil) best <- list(k=k, sil=s, cl=cl) }
  if (is.na(best$sil) || best$sil < 0.10) best <- list(k=1, sil=best$sil, cl=setNames(rep(1L,ncol(A)), colnames(A)))
  best }

arm_delta <- function(A, cl){
  tb <- sort(table(cl), decreasing=TRUE); if (length(tb)<2) return(NULL)
  rowMeans(A[, cl==names(tb)[1], drop=FALSE]) - rowMeans(A[, cl==names(tb)[2], drop=FALSE]) }

## ---- the whole point of this script -----------------------------------------------
## Per-cell z across the 39 arms. Unlike median centring this also divides out per-cell
## AMPLITUDE, so the clustering compares profile SHAPE rather than how strong a cell's
## overall CNV signal happens to be. Adopted variant; see header for the metrics.
zrows <- function(A){                      # per-cell z ACROSS ARMS (columns = cells)
  mu <- colMeans(A, na.rm=TRUE); sdv <- apply(A, 2, sd, na.rm=TRUE); sdv[!is.finite(sdv)|sdv==0] <- 1
  sweep(sweep(A, 2, mu, "-"), 2, sdv, "/") }

ann   <- read.csv(file.path(SC,"tumor_compartment","scatac_ArchR","barcode_annotation.csv"), stringsAsFactors=FALSE)
malig <- ann$barcode[ann$celltype_lv1 %in% c("Tumor","Malignant")]

## Two readings of "no reference subtraction", built side by side so the comparison is
## explicit: zscore KEEPS the GC-matched background (log2FC) and only changes the
## per-cell normalisation; noref discards the background entirely and uses raw counts.
for (V in c("zscore","noref")){
  DIR <- paste0("arm_cnv_", V); dir.create(DIR, showWarnings=FALSE)
  cat("\n##################### variant:", V, "#####################\n")
  ARM <- list(); BAD <- list()
  for (s in ATAC_SAMP){
    o  <- readRDS(file.path(CND, sprintf("CNV_LFC_GC_%s_ws_1e+07_ss_5e+06.rds", s)))
    M  <- if (V=="zscore") as.matrix(assays(o)$log2FC) else as.matrix(assays(o)$counts)
    cn <- colnames(o); colnames(M) <- ifelse(grepl("#",cn), cn, paste0(s,"#",cn))
    rr <- rowRanges(o); wchr <- as.character(seqnames(rr)); wmb <- (start(rr)+end(rr))/2/1e6
    keep <- wchr %in% AUT; M <- M[keep,,drop=FALSE]; wchr <- wchr[keep]; wmb <- wmb[keep]
    cells <- intersect(malig, colnames(M)); M <- M[,cells,drop=FALSE]
    ## noref: log1p the raw counts and z-score per cell across windows. Recorded result --
    ## this reintroduces the full gene-density/accessibility bias (17p/q, 19p/q, 20q, 22q
    ## light up in every sample), which is why zscore is the adopted variant.
    if (V=="noref"){                                   # raw counts, no background at all
      M[!is.finite(M)] <- 0; M <- log1p(M); M <- scale(M); M[!is.finite(M)] <- 0 }
    aw <- armof(wchr,wmb); wid <- sprintf("%s:%.0fMb", wchr, wmb)
    cw <- clean_windows(M, wid, aw, paste0("scATAC_",s)); M <- cw$M; aw <- cw$arm
    BAD[[s]] <- cw$bad
    use <- names(which(table(aw) >= MINWIN)); use <- intersect(ARMS, use)
    A <- t(sapply(use, function(a) colMeans(M[aw==a,,drop=FALSE], na.rm=TRUE)))
    A <- zrows(A)                                      # <-- per-cell z across arms, no reference subtraction
    ARM[[s]] <- A
    cat(sprintf("scATAC %-4s : %d cells, %d windows kept (%d dropped) -> %d arms\n",
                s, ncol(A), nrow(M), nrow(cw$bad), nrow(A)))
  }
  armset <- Reduce(intersect, lapply(ARM, rownames))
  armlv  <- unlist(lapply(1:22, function(i) intersect(paste0(i,c("p","q")), armset)))
  for (s in names(ARM)) ARM[[s]] <- ARM[[s]][armlv,,drop=FALSE]
  cat("common arms:", length(armlv), "\n")

  ## clone discovery, identical caller to R2_Q3_arm_level_cnv.R, now on the z-scored arms
  calls <- list(); tops <- list(); DELTA <- list()
  for (s in names(ARM)){
    A  <- ARM[[s]]
    cc <- call_clones_arm(A, "euclid"); lg <- call_clones_arm(A, "pearson")
    d  <- arm_delta(A, cc$cl); key <- paste0(s,"_scATAC")
    saveRDS(list(A=A, clones=cc$cl, k=cc$k, sil=cc$sil, delta=d, sample=s, modality="scATAC"),
            file.path(DIR, paste0(key,"_arm.rds")))
    write.csv(data.frame(cell=names(cc$cl), clone=as.integer(cc$cl)),
              file.path(DIR, paste0(key,"_clones.csv")), row.names=FALSE)
    calls[[s]] <- data.frame(sample=s, modality="scATAC", n_cells=ncol(A),
                             k_arm=cc$k, sil_arm=round(ifelse(is.na(cc$sil),0,cc$sil),2),
                             k_armPearson=lg$k, sil_armPearson=round(ifelse(is.na(lg$sil),0,lg$sil),2),
                             minor_frac=round(ifelse(cc$k>1, min(table(cc$cl))/ncol(A), NA),2))
    if (!is.null(d)){ DELTA[[key]] <- d; o5 <- order(-abs(d))[1:5]
      tops[[s]] <- data.frame(sample=s, modality="scATAC", rank=1:5, arm=names(d)[o5],
                              delta=round(as.numeric(d)[o5],3)) }
  }
  ## only the scATAC normalisation is being varied, so the scRNA arm matrices from
  ## R2_Q3_arm_level_cnv.R are copied in unchanged -- this keeps each variant directory
  ## self-contained and lets the downstream scripts run on it without edits
  for (f in list.files("arm_cnv", "_scRNA_arm\\.rds$", full.names=TRUE)) file.copy(f, DIR, overwrite=TRUE)
  D0 <- readRDS("arm_cnv/arm_level_deltas.rds")
  for (k in grep("_scRNA$", names(D0), value=TRUE)) DELTA[[k]] <- D0[[k]][armlv]
  saveRDS(DELTA, file.path(DIR,"arm_level_deltas.rds"))
  callsR <- read.csv("arm_level_clone_calls.csv", stringsAsFactors=FALSE)
  calls  <- rbind(do.call(rbind, calls), callsR[callsR$modality=="scRNA",])
  rownames(calls) <- NULL; calls <- calls[order(calls$sample, calls$modality),]
  write.csv(calls, sprintf("arm_level_clone_calls_%s.csv", V), row.names=FALSE)
  tops <- do.call(rbind, tops); rownames(tops) <- NULL
  write.csv(tops, sprintf("arm_level_top_arms_%s.csv", V), row.names=FALSE)
  BAD <- do.call(rbind, BAD); rownames(BAD) <- NULL
  write.csv(BAD, sprintf("arm_level_bad_windows_%s.csv", V), row.names=FALSE)

  cat("\n=== calls (", V, ") ===\n", sep=""); print(calls, row.names=FALSE)
  cat("\n=== top discriminating arms (", V, ") ===\n", sep=""); print(tops[tops$rank<=3,], row.names=FALSE)
  cat("\n=== top-arm recurrence across unrelated tumours (artefact guard) ===\n")
  t1 <- tops[tops$rank==1,]; print(sort(table(t1$arm), decreasing=TRUE))
  cat("\n=== P4 chr8q check ===\n")
  if (is.null(DELTA[["P4_scATAC"]])) cat("P4 called clonal\n") else {
    d <- DELTA[["P4_scATAC"]]; r <- rank(-abs(d))[["8q"]]
    cat(sprintf("8q delta=%+.3f, rank %d/%d\n", d[["8q"]], r, length(d))) }

  ## heatmaps
  col_fun <- colorRamp2(c(-2,-0.6,0,0.6,2), c("#08306b","#6baed6","white","#fb6a4a","#67000d"))
  pdf(sprintf("Plots/R2_Q3_arm_level_heatmaps_%s.pdf", V), width=9, height=5.6)
  for (s in names(ARM)){
    A <- ARM[[s]]; cc <- readRDS(file.path(DIR, paste0(s,"_scATAC_arm.rds")))
    cl <- cc$clones; cells <- colnames(A)
    if (length(cells) > MAXCELL){ idx <- sort(sample(seq_along(cells), MAXCELL)); A <- A[,idx,drop=FALSE]; cl <- cl[cells[idx]] }
    kk <- length(unique(cl)); rs <- if (kk>1) factor(paste0("clone", cl[colnames(A)])) else NULL
    draw(Heatmap(t(A), name="arm z", col=col_fun, cluster_columns=FALSE, cluster_rows=(kk>1),
      show_row_names=FALSE, column_names_gp=gpar(fontsize=7), row_split=rs, row_gap=unit(1.2,"mm"),
      border=TRUE, use_raster=TRUE, raster_quality=2,
      row_title=if(kk>1) NULL else "clonal", row_title_gp=gpar(fontsize=9),
      column_title=sprintf("%s scATAC arm-level CNV (%s: per-cell z across arms, no reference) - k=%d, sil=%.2f",
                           s, V, cc$k, ifelse(is.na(cc$sil),0,cc$sil)),
      column_title_gp=gpar(fontface="bold", fontsize=10)), heatmap_legend_side="right")
  }
  dev.off()
}
cat("\nDONE\n")
