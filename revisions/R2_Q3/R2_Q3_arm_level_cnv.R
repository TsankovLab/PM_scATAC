###############################################################
# R2_Q3 : ARM-LEVEL clone discovery, both modalities.
#
# Rationale (from the P4 chr8q lesson): genome-wide 1-Mb / 10-Mb clustering dilutes
# a real arm-level subclonal event across ~2000 shared bins, and single-bin noise
# dominates the cell-cell distance. Averaging every bin of a chromosome arm into ONE
# score cancels noise ~sqrt(N) (that is exactly how the P4 chr8q clone was recovered).
# Here we do that systematically: per cell, aggregate the CNV signal per chromosome
# arm in BOTH modalities, then re-run clone discovery in the 39-arm feature space.
#
#   scATAC : GL scCNA log2FC (GC-corrected 10-Mb windows, 5-Mb step)  -> arm means
#   scRNA  : binned-expression CNV (1-Mb bins vs Mesothelium ref, no inferCNV) -> arm means
#   acrocentric p-arms (13p/14p/15p/21p/22p) excluded (satellite/rDNA artefacts)
#
# WINDOW/BIN QC (essential -- found the hard way, see arm_level_bad_windows.csv):
#   the GL CNV_LFC_GC objects contain pericentromeric / acrocentric heterochromatin
#   windows whose log2FC collapses to ~-8..-24 in a variable fraction of cells
#   (chr16:40Mb at the 16q centromere edge, chr1:125-140Mb, chr9:60-65Mb, chr14:15Mb,
#   chr15:15-20Mb). Their cross-cell SD is 2-9 vs ~0.35 for a normal window, so a naive
#   arm MEAN is entirely driven by them: without this QC, 16q is the top "clone"-
#   discriminating arm in 7/9 scATAC samples with an identical delta ~1.8 (i.e. the
#   split is dropout, not CNV). We therefore drop windows whose robust cross-cell SD
#   exceeds max(1.5, 4x the median window SD), then winsorise to +-2 before averaging.
#
# Two clusterings per sample/modality:
#   (a) PRIMARY  = arm matrix CENTRED PER ARM across the sample's cells (removes the
#       truncal/shared CNV that makes every malignant cell look alike) + Euclidean
#       ward.D2 -> subclonal variation only.
#   (b) legacy   = 1-Pearson on the uncentred arm profiles (comparable to the previous
#       bin-level calls) -- reported side by side so the change is attributable.
# k chosen by silhouette (k=2..6) with the same robustness rule used throughout R2_Q3:
#   every cluster >= max(20, 5%) cells and best silhouette >= 0.10, else clonal.
#
# Outputs: arm_cnv/<sample>_<modality>_arm.rds (arm matrix + labels),
#          arm_level_clone_calls.csv, arm_level_top_arms.csv,
#          Plots/R2_Q3_arm_level_heatmaps.pdf
###############################################################
suppressMessages({ library(Seurat); library(Matrix); library(SummarizedExperiment); library(GenomicRanges)
                   library(cluster); library(ComplexHeatmap); library(circlize); library(grid); library(gridExtra) })
set.seed(1)
SC  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
CND <- file.path(SC,"per_sample_QC_signac","CNV_analysis")
OUT <- file.path(SC,"git_repo_claude","R2_Q3"); setwd(OUT)
dir.create("Plots", showWarnings=FALSE); dir.create("arm_cnv", showWarnings=FALSE)
BIN <- 1e6; AUT <- paste0("chr",1:22); MAXCELL <- 1200; MINWIN <- 2

## centromere positions (Mb, hg38) -> p/q assignment
CEN <- c(1,123.6,2,93.1,3,91.1,4,50.4,5,48.3,6,59.2,7,59.5,8,44.9,9,44.4,10,40.6,11,52.7,12,35.9,
         13,17.0,14,17.1,15,19.1,16,36.3,17,23.8,18,17.7,19,25.8,20,28.0,21,11.9,22,14.5)
CENmb <- setNames(CEN[seq(2,length(CEN),2)], paste0("chr",CEN[seq(1,length(CEN),2)]))
ACRO  <- c("13p","14p","15p","21p","22p")
ARMS  <- setdiff(unlist(lapply(1:22, function(i) paste0(i, c("p","q")))), ACRO)
armof <- function(chr,mb) paste0(sub("chr","",chr), ifelse(mb < CENmb[chr],"p","q"))

## ---- STEP 1 (helper): window QC -- kill heterochromatin-dropout windows, winsorise --
## Two independent drop rules, BOTH required:
##   (a) mad(window) > max(1.5, 4*median mad)  -> grossly variable across cells
##   (b) |median(window)| > 1.5 or non-finite  -> systematically collapsed/blown up
## Rule (b) is what catches chr16:40Mb, which is CONSISTENTLY broken rather than
## variable and slips past the MAD rule. Survivors are winsorised to +-2 so a residual
## outlier window cannot dominate the arm mean it feeds into.
## returns list(M=cleaned matrix, bad=data.frame of dropped windows)
clean_windows <- function(M, wid, arm, tag, clip=2){
  sdw <- apply(M, 1, function(z) mad(z[is.finite(z)]))
  sdw[!is.finite(sdw)] <- 0
  medw <- apply(M, 1, function(z) median(z[is.finite(z)]))
  thr <- max(1.5, 4*median(sdw[sdw>0]))
  ## (a) grossly variable across cells, (b) systematically collapsed/blown-up window
  bad <- which(sdw > thr | !is.finite(medw) | abs(medw) > 1.5)
  why <- ifelse(!is.finite(medw[bad]), "nonfinite",
         ifelse(abs(medw[bad]) > 1.5, "extreme_median", "extreme_mad"))
  bd <- if (!length(bad)) data.frame(source=character(0), window=character(0), arm=character(0),
                                     median=numeric(0), mad=numeric(0), threshold=numeric(0),
                                     rule=character(0), stringsAsFactors=FALSE) else
        data.frame(source=tag, window=wid[bad], arm=arm[bad],
                   median=round(medw[bad],2), mad=round(sdw[bad],2), threshold=round(thr,2),
                   rule=why, stringsAsFactors=FALSE)
  if (length(bad)) { M <- M[-bad,,drop=FALSE]; arm <- arm[-bad] }
  M[!is.finite(M)] <- NA
  M <- pmax(pmin(M, clip), -clip)                       # winsorise residual outliers
  list(M=M, arm=arm, bad=bd) }

smooth_chr <- function(v, chr, k=5){ o<-v; for(ch in unique(chr)){ i<-which(chr==ch); x<-v[i]; n<-length(x); s<-x
  if(n>=k) for(j in seq_len(n)){ lo<-max(1,j-k%/%2); hi<-min(n,j+k%/%2); s[j]<-mean(x[lo:hi]) }; o[i]<-s }; o }

## ---- STEP 2 (helper): clone discovery in arm space --------------------------------
## "euclid" (primary): centre PER ARM across the sample's cells first. That removes the
## truncal CNV shared by every malignant cell, leaving only between-cell variation --
## without it the distance is dominated by what the clones have in common.
## k is chosen by silhouette over 2..6, a k is admissible only if every cluster holds
## >= max(20, 5%) cells, and a best silhouette < 0.10 returns CLONAL (k=1). The caller
## can therefore answer "no subclones", which is what makes the negative controls work.
## d = "euclid" on per-arm-centred matrix (primary) or "pearson" on raw arm profiles (legacy)
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

## per-arm difference between the two largest clones (signed, clone1 = larger)
arm_delta <- function(A, cl){
  tb <- sort(table(cl), decreasing=TRUE); if (length(tb)<2) return(NULL)
  g1 <- names(tb)[1]; g2 <- names(tb)[2]
  d <- rowMeans(A[, cl==g1, drop=FALSE]) - rowMeans(A[, cl==g2, drop=FALSE])
  attr(d,"n") <- c(as.integer(tb[1]), as.integer(tb[2])); d }

######################## STEP 3. scATAC : GL log2FC -> arm scores ########################
## Precomputed Granja-Lareau scCNA objects: 10 Mb windows / 5 Mb step, Tn5 insertions
## normalised against the 100 nearest-GC windows. We take log2FC (NOT the z assay: z
## divides by the background SD and is dominated by low-variance-background windows).
## Per sample: autosomes -> malignant cells -> window QC -> mean per arm -> per-cell
## median centring.
ann   <- read.csv(file.path(SC,"tumor_compartment","scatac_ArchR","barcode_annotation.csv"), stringsAsFactors=FALSE)
malig <- ann$barcode[ann$celltype_lv1 %in% c("Tumor","Malignant")]
ATAC_SAMP <- c("P1","P4","P5","P8","P10","P11","P12","P14","P23")
ARM <- list()   # ARM[[modality]][[sample]] = arms x cells
BAD <- list()   # dropped windows/bins
ARM$scATAC <- list()
for (s in ATAC_SAMP){
  o  <- readRDS(file.path(CND, sprintf("CNV_LFC_GC_%s_ws_1e+07_ss_5e+06.rds", s)))
  M  <- as.matrix(assays(o)$log2FC)
  cn <- colnames(o); colnames(M) <- ifelse(grepl("#",cn), cn, paste0(s,"#",cn))
  rr <- rowRanges(o); wchr <- as.character(seqnames(rr)); wmb <- (start(rr)+end(rr))/2/1e6
  keep <- wchr %in% AUT; M <- M[keep,,drop=FALSE]; wchr <- wchr[keep]; wmb <- wmb[keep]
  cells <- intersect(malig, colnames(M)); M <- M[,cells,drop=FALSE]
  aw <- armof(wchr,wmb); wid <- sprintf("%s:%.0fMb", wchr, wmb)
  cw <- clean_windows(M, wid, aw, paste0("scATAC_",s)); M <- cw$M; aw <- cw$arm
  BAD[[paste0("scATAC_",s)]] <- cw$bad
  use <- names(which(table(aw) >= MINWIN)); use <- intersect(ARMS, use)
  ## arm score = mean of that arm's surviving windows. Averaging N windows cancels noise
  ## ~sqrt(N); this is exactly why a single-arm event (P4 chr8q) becomes visible at all.
  A <- t(sapply(use, function(a) colMeans(M[aw==a,,drop=FALSE], na.rm=TRUE)))   # ARM AGGREGATION
  A <- sweep(A, 2, apply(A,2,median,na.rm=TRUE), "-")                           # per-cell centring
  ARM$scATAC[[s]] <- A
  cat(sprintf("scATAC %-4s : %d cells, %d windows kept (%d dropped) -> %d arms\n",
              s, ncol(A), nrow(M), nrow(cw$bad), nrow(A)))
}
rm(M); gc()

######################## STEP 4. scRNA : binned expression -> arm scores ########################
## Deliberately NOT inferCNV values (its tumour_subclusters are over-segmented). Genes ->
## 1 Mb bins -> bin mean per cell -> subtract the same bin's mean over benign Mesothelium
## cells -> per-cell median centring -> 5-bin running mean within each chromosome -> the
## same window QC as scATAC -> mean per arm.
ic <- readRDS(file.path(SC,"tumor_compartment","scrna","infercnv","run.final.infercnv_obj")); go <- ic@gene_order; rm(ic); gc()
gpos <- data.frame(gene=rownames(go), chr=as.character(go$chr), start=go$start, stringsAsFactors=FALSE)
gpos <- gpos[gpos$chr %in% AUT,]; gpos$bin <- paste0(gpos$chr,":",floor(gpos$start/BIN))
srt <- readRDS(file.path(SC,"tumor_compartment","scrna","scRNA_meso.rds")); srt <- NormalizeData(srt, verbose=FALSE)
Mx  <- GetAssayData(srt, slot="data")
mal <- colnames(srt)[srt$celltype=="Malignant"]; refc <- colnames(srt)[srt$celltype=="Mesothelium"]
g   <- intersect(gpos$gene, rownames(Mx)); gp <- gpos[match(g,gpos$gene),]; Mg <- Mx[g,]
refm <- Matrix::rowMeans(Mg[,refc,drop=FALSE]); Aagg <- fac2sparse(factor(gp$bin)); gcb <- as.numeric(Matrix::rowSums(Aagg))
RNA <- sweep(as.matrix(Aagg %*% Mg[,mal,drop=FALSE])/gcb, 1, as.numeric(Aagg %*% refm)/gcb, "-")
rownames(RNA) <- levels(factor(gp$bin))
rchr <- sub(":.*","",rownames(RNA)); rmb <- as.integer(sub(".*:","",rownames(RNA)))
oo <- order(match(rchr,AUT), rmb); RNA <- RNA[oo,]; rchr <- rchr[oo]; rmb <- rmb[oo]
RNA <- sweep(RNA, 2, apply(RNA,2,median), "-")
RNA <- apply(RNA, 2, function(z) smooth_chr(z, rchr))
rna_arm <- armof(rchr, rmb)
cwR <- clean_windows(RNA, rownames(RNA), rna_arm, "scRNA_allbins")   # same QC as scATAC
RNA <- cwR$M; rna_arm <- cwR$arm; BAD[["scRNA"]] <- cwR$bad
cat(sprintf("scRNA bin QC : %d bins kept, %d dropped\n", nrow(RNA), nrow(cwR$bad)))
rm(Mx, Mg, srt); gc()

RNA_SAMP <- c("P1","P3","P4","P5","P8","P11","P12","P13","P14")
ARM$scRNA <- list()
for (s in RNA_SAMP){
  cells <- intersect(read.csv(file.path("scRNA_cnv",paste0(s,"_subclones.csv")), stringsAsFactors=FALSE)$cell, colnames(RNA))
  if (length(cells) < 40) { cat(sprintf("scRNA  %-4s : only %d cells -> skipped\n", s, length(cells))); next }
  use <- names(which(table(rna_arm) >= MINWIN)); use <- intersect(ARMS, use)
  A <- t(sapply(use, function(a) colMeans(RNA[rna_arm==a, cells, drop=FALSE], na.rm=TRUE)))
  A <- sweep(A, 2, apply(A,2,median,na.rm=TRUE), "-")
  ARM$scRNA[[s]] <- A
  cat(sprintf("scRNA  %-4s : %d cells, %d bins -> %d arms\n", s, ncol(A), nrow(RNA), nrow(A)))
}
rm(RNA); gc()

## intersect the arm sets so every sample and both modalities share one feature space --
## required before any cross-modal comparison of arm profiles or deltas
armset <- Reduce(intersect, lapply(unlist(ARM, recursive=FALSE), rownames))
armlv  <- unlist(lapply(1:22, function(i) intersect(paste0(i,c("p","q")), armset)))
cat("\ncommon arm features:", length(armlv), "->", paste(armlv, collapse=","), "\n\n")
for (md in names(ARM)) for (s in names(ARM[[md]])) ARM[[md]][[s]] <- ARM[[md]][[s]][armlv,,drop=FALSE]

######################## STEP 5. clone discovery, per sample and modality ########################
## Primary = euclid on the per-arm-centred matrix; a legacy 1-Pearson call on the
## uncentred profiles is computed alongside so the change from the bin-level era is
## attributable rather than assumed.
calls <- list(); tops <- list(); DELTA <- list()
for (md in names(ARM)) for (s in names(ARM[[md]])){
  A  <- ARM[[md]][[s]]
  cc <- call_clones_arm(A, "euclid")        # primary
  lg <- call_clones_arm(A, "pearson")       # legacy comparison
  d  <- arm_delta(A, cc$cl)
  key <- paste0(s,"_",md)
  saveRDS(list(A=A, clones=cc$cl, k=cc$k, sil=cc$sil, delta=d, sample=s, modality=md),
          file.path("arm_cnv", paste0(key,"_arm.rds")))
  write.csv(data.frame(cell=names(cc$cl), clone=as.integer(cc$cl)),
            file.path("arm_cnv", paste0(key,"_clones.csv")), row.names=FALSE)
  calls[[key]] <- data.frame(sample=s, modality=md, n_cells=ncol(A),
                             k_arm=cc$k, sil_arm=round(ifelse(is.na(cc$sil),0,cc$sil),2),
                             k_armPearson=lg$k, sil_armPearson=round(ifelse(is.na(lg$sil),0,lg$sil),2),
                             minor_frac=round(ifelse(cc$k>1, min(table(cc$cl))/ncol(A), NA),2))
  if (!is.null(d)){
    DELTA[[key]] <- d
    o3 <- order(-abs(d))[1:5]
    tops[[key]] <- data.frame(sample=s, modality=md, rank=1:5, arm=names(d)[o3], delta=round(as.numeric(d)[o3],3))
  }
}
calls <- do.call(rbind, calls); rownames(calls) <- NULL
calls <- calls[order(calls$sample, calls$modality),]
write.csv(calls, "arm_level_clone_calls.csv", row.names=FALSE)
tops  <- do.call(rbind, tops); rownames(tops) <- NULL
write.csv(tops, "arm_level_top_arms.csv", row.names=FALSE)
saveRDS(DELTA, "arm_cnv/arm_level_deltas.rds")

BAD <- do.call(rbind, BAD); rownames(BAD) <- NULL
write.csv(BAD, "arm_level_bad_windows.csv", row.names=FALSE)
cat("=== dropped windows/bins (heterochromatin dropout QC) ===\n")
cat(nrow(BAD), "dropped in total; recurrent offenders across samples:\n")
print(head(sort(table(BAD$window), decreasing=TRUE), 15))

cat("\n=== ARM-LEVEL clone calls (primary = per-arm-centred Euclidean) ===\n"); print(calls, row.names=FALSE)
cat("\n=== top discriminating arms per split (clone1 - clone2, clone1 = larger) ===\n")
if (!is.null(tops)) print(tops[tops$rank<=3,], row.names=FALSE, digits=3)

## artefact guard: an arm that is the TOP discriminator in many unrelated tumours is
## technical, not CNV (this is exactly how the 16q dropout artefact was caught).
cat("\n=== how often is each arm the top discriminator across (unrelated) tumours? ===\n")
for (md in c("scATAC","scRNA")){
  t1 <- tops[tops$modality==md & tops$rank==1,]
  if (!nrow(t1)) next
  cat(md, ": ", paste(sprintf("%s x%d", names(sort(table(t1$arm),decreasing=TRUE)),
                              sort(table(t1$arm),decreasing=TRUE)), collapse=", "), "\n", sep="")
}

## sanity check: does the arm-level scan now recover the spatially-validated P4 chr8q clone?
cat("\n=== P4 chr8q recovery check ===\n")
for (md in c("scATAC","scRNA")){
  key <- paste0("P4_",md); if (is.null(DELTA[[key]])) { cat(md,": P4 called clonal\n"); next }
  d <- DELTA[[key]]; r <- rank(-abs(d))[["8q"]]
  cat(sprintf("%s : 8q delta=%+.3f, rank %d/%d by |delta| (pct %.2f)\n",
              md, d[["8q"]], r, length(d), 1-(r-1)/length(d)))
}

######################## STEP 6. per-sample arm heatmaps ########################
col_fun <- colorRamp2(c(-2,-0.6,0,0.6,2), c("#08306b","#6baed6","white","#fb6a4a","#67000d"))
mk_ht <- function(A, cl, title){
  Ad <- A/(sd(A)+1e-9)
  cells <- colnames(Ad)
  if (length(cells) > MAXCELL){ idx <- sort(sample(seq_along(cells), MAXCELL)); Ad <- Ad[,idx,drop=FALSE]; cl <- cl[cells[idx]] }
  kk <- length(unique(cl)); rs <- if (kk>1) factor(paste0("clone", cl[colnames(Ad)])) else NULL
  Heatmap(t(Ad), name="arm CNV", col=col_fun, cluster_columns=FALSE, cluster_rows=(kk>1),
    show_row_names=FALSE, show_column_names=TRUE, column_names_gp=gpar(fontsize=7),
    row_split=rs, row_gap=unit(1.2,"mm"), border=TRUE, use_raster=TRUE, raster_quality=2,
    row_title=if(kk>1) NULL else "clonal", row_title_gp=gpar(fontsize=9),
    column_title=title, column_title_gp=gpar(fontface="bold", fontsize=10),
    heatmap_legend_param=list(direction="horizontal")) }

pdf("Plots/R2_Q3_arm_level_heatmaps.pdf", width=15, height=5.6)
for (s in sort(unique(calls$sample))){
  panels <- list()
  for (md in c("scATAC","scRNA")){
    if (is.null(ARM[[md]][[s]])) next
    A <- ARM[[md]][[s]]; r <- calls[calls$sample==s & calls$modality==md,]
    cl <- readRDS(file.path("arm_cnv", paste0(s,"_",md,"_arm.rds")))$clones
    ttl <- sprintf("%s  %s  (arm-level CNV)  -  k=%d, sil=%.2f, n=%d", s, md, r$k_arm, r$sil_arm, r$n_cells)
    panels[[md]] <- grid.grabExpr(draw(mk_ht(A, cl, ttl), heatmap_legend_side="bottom"))
  }
  grid.newpage(); grid.draw(arrangeGrob(grobs=panels, ncol=length(panels)))
}
dev.off()

cat("\nDONE\n")
