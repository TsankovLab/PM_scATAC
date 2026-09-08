###############################################################################
# STEP 4 -- do the two annotations describe the same cells, without any integration?
#
# Step 3 depends on a CCA anchor transfer; if the transfer is wrong, the agreement
# number is wrong.  This step asks the same question with no integration at all.
# The only thing shared between the two modalities here is the gene identity.
#
#   A. profile correlation
#      average the ArchR gene score of each ATAC cell type, and the log-normalised
#      expression of each RNA cell type, over the cell types present in BOTH; z-score
#      each gene across those same cell types in each modality (removes the modality
#      baseline -- gene scores carry gene length and promoter accessibility, and the
#      z must be taken over an identical column set or the profiles are not
#      comparable); correlate every ATAC profile against every RNA profile.
#
#   B. marker enrichment  (the primary test -- more robust than A)
#      take the top NMARK genes that define each cell type IN THE RNA DATA, then ask
#      how specifically those genes are accessible across the ATAC cell types.  If
#      the two annotations name the same biology, the RNA markers of cell type k must
#      score highest in the ATAC cells the chromatin annotation called k.
#      A permutation over random gene sets of the same size gives the null.
#
# Cell types are restricted to those with >= MINCELLS cells in both modalities.
# Glia is therefore dropped: the chromatin annotation has no Glia cluster at all
# (the RNA object has 55 Glia cells).  Plasma and Alveolar sit just above the floor
# in both assays and their profiles are correspondingly noisy -- flagged, not hidden.
#
# Input : main/scatac_ArchR, main/scrna/srt.rds
# Output: pseudobulk_atac.csv, pseudobulk_rna.csv, pseudobulk_correlation.csv,
#         pseudobulk_bestmatch.csv, marker_enrichment.csv, marker_bestmatch.csv,
#         rna_markers.csv, marker_genescore_by_atac_label.csv
###############################################################################
source("00_common.R")
suppressMessages({ library(ArchR); library(Seurat); library(SummarizedExperiment) })
addArchRThreads(4)

MINCELLS <- 50     # a cell type must have this many cells in BOTH modalities
NGENES   <- 2000   # genes used for the profile correlation
NMARK    <- 100    # RNA marker genes per cell type
NPERM    <- 1000   # random gene sets for the marker-enrichment null

## ---- ATAC pseudobulk ---------------------------------------------------------
proj <- loadArchRProject(ARCHR, showLogo = FALSE)
proj@projectMetadata$outputDirectory <- file.path(ROOT, "archr_scratch")
dir.create(file.path(ROOT, "archr_scratch"), showWarnings = FALSE)
gse <- getGroupSE(proj, useMatrix = "GeneScoreMatrix", groupBy = ATAC_LABEL,
                  divideN = TRUE, scaleTo = 10000)
GA <- log2(as.matrix(assay(gse)) + 1)
rownames(GA) <- rowData(gse)$name
GA <- GA[!duplicated(rownames(GA)), , drop = FALSE]
nA <- table(as.character(getCellColData(proj)[[ATAC_LABEL]]))
rm(proj); invisible(gc(FALSE))

## ---- RNA pseudobulk ----------------------------------------------------------
s   <- readRDS(SRNA)
lab <- as.character(s@meta.data[[RNA_LABEL]])
D   <- SeuratObject::LayerData(s[["RNA"]], layer = "data")   # log-normalised
nR  <- table(lab)
GR  <- vapply(names(nR), function(k) Matrix::rowMeans(D[, lab == k, drop = FALSE]),
              numeric(nrow(D)))
rownames(GR) <- rownames(D)
rm(s, D); invisible(gc(FALSE))

## ---- the cell types the two annotations share --------------------------------
CT <- intersect(names(nA)[nA >= MINCELLS], names(nR)[nR >= MINCELLS])
CT <- CELLTYPES[CELLTYPES %in% CT]
GA <- GA[, CT, drop = FALSE]; GR <- GR[, CT, drop = FALSE]
cat("cell types profiled in both modalities:", length(CT), "->",
    paste(CT, collapse = ", "), "\n")
cat("dropped:", paste(setdiff(union(names(nA), names(nR)), CT), collapse = ", "), "\n")
cat("ATAC pseudobulk:", nrow(GA), "genes | RNA pseudobulk:", nrow(GR), "genes\n")
write.csv(round(GA, 4), "pseudobulk_atac.csv")
write.csv(round(GR, 4), "pseudobulk_rna.csv")

## ---- z within modality, over the SAME columns --------------------------------
g  <- intersect(rownames(GA), rownames(GR))
zA <- t(scale(t(GA[g, , drop = FALSE])))
zR <- t(scale(t(GR[g, , drop = FALSE])))
ok <- complete.cases(zA) & complete.cases(zR)
zA <- zA[ok, ]; zR <- zR[ok, ]
cat("genes measured in both:", nrow(zA), "\n")

## ---- A. profile correlation ---------------------------------------------------
## Pearson, not Spearman: the z profiles are near zero for most genes, and ranking
## 2000 mostly-flat genes lets the noise dominate (Spearman puts the diagonal first
## for only 6 of 12 cell types, Pearson for 9; the Spearman matrix is written too).
sel <- head(order(apply(zR, 1, var) + apply(zA, 1, var), decreasing = TRUE), NGENES)
CM  <- cor(zA[sel, ], zR[sel, ], method = "pearson")
write.csv(round(CM, 4), "pseudobulk_correlation.csv")
write.csv(round(cor(zA[sel, ], zR[sel, ], method = "spearman"), 4),
          "pseudobulk_correlation_spearman.csv")
cat("\n=== A. Pearson correlation of z-scored profiles",
    sprintf("(%d genes) ===\n", length(sel)))
cat("    rows = chromatin-defined cell type, columns = expression-defined cell type\n")
print(round(CM, 2))

bm <- function(M, tag){
  d <- data.frame(atac_label = rownames(M),
                  best_rna = colnames(M)[apply(M, 1, which.max)],
                  best = apply(M, 1, max),
                  same = M[cbind(seq_len(nrow(M)), match(rownames(M), colnames(M)))],
                  runner_up = apply(M, 1, function(x) colnames(M)[order(-x)[2]]),
                  runner_up_v = apply(M, 1, function(x) sort(x, decreasing = TRUE)[2]),
                  stringsAsFactors = FALSE)
  d$margin  <- d$best - d$runner_up_v
  d$correct <- d$best_rna == d$atac_label
  d$rank_of_same <- vapply(seq_len(nrow(M)), function(i)
    which(order(-M[i, ]) == match(rownames(M)[i], colnames(M))), integer(1))
  cat(sprintf("\n=== %s: best RNA match for each chromatin-defined cell type ===\n", tag))
  print(transform(d, best = round(best, 3), same = round(same, 3),
                  runner_up_v = round(runner_up_v, 3), margin = round(margin, 3)),
        row.names = FALSE)
  cat(sprintf("diagonal is the best match for %d of %d cell types; median rank of the matching type = %g\n",
              sum(d$correct), nrow(d), median(d$rank_of_same)))
  d
}
BM <- bm(CM, "A. profile correlation")
write.csv(BM, "pseudobulk_bestmatch.csv", row.names = FALSE)

## ---- B. RNA markers scored in ATAC chromatin ----------------------------------
## marker = highest z in that RNA cell type, i.e. defined without touching the ATAC.
MK <- lapply(CT, function(k) rownames(zR)[head(order(zR[, k], decreasing = TRUE), NMARK)])
names(MK) <- CT
write.csv(data.frame(celltype = rep(CT, each = NMARK), gene = unlist(MK)),
          "rna_markers.csv", row.names = FALSE)

ME <- t(vapply(CT, function(k) colMeans(zA[MK[[k]], , drop = FALSE]), numeric(length(CT))))
ME <- t(ME)                                  # rows = ATAC cell type, cols = marker set
write.csv(round(ME, 4), "marker_enrichment.csv")
cat(sprintf("\n=== B. mean ATAC gene-score z of each RNA marker set (%d genes/set) ===\n", NMARK))
cat("    rows = chromatin-defined cell type, columns = the RNA cell type the markers came from\n")
print(round(ME, 2))
BM2 <- bm(ME, "B. marker enrichment")
write.csv(BM2, "marker_bestmatch.csv", row.names = FALSE)

## permutation null: random gene sets of the same size
set.seed(1)
nullmax <- replicate(NPERM, {
  r <- sample(rownames(zA), NMARK); max(abs(colMeans(zA[r, , drop = FALSE]))) })
cat(sprintf("\nnull (%d random %d-gene sets): max |mean z| = %.3f (95th pct %.3f)\n",
            NPERM, NMARK, max(nullmax), quantile(nullmax, .95)))
diagv <- diag(ME)
cat(sprintf("matched-diagonal z: %s\n",
            paste(sprintf("%s %.2f", names(diagv), diagv), collapse = " | ")))
cat(sprintf("cell types whose matched z exceeds the 95th-percentile null: %d of %d\n",
            sum(diagv > quantile(nullmax, .95)), length(diagv)))

## ---- provenance: the markers the chromatin annotation was made from -----------
## (git_repo/main_analysis/scatac_main_ArchR.R, "Plot gene score of cell type markers")
MARKERS <- c("KRT19","AMOTL2","EGFR","HP","BDKRB1","ZBTB7C","SFTA3","SFTPB","LINC00261",
             "COL1A1","FBN1","COL6A2","MYH11","COL4A2","COL4A1","CLDN5","VWF","CDH5",
             "LYZ","IL1B","CD83","CD3E","BCL11B","RUNX3","GNLY","KLRD1","ITGAE",
             "CD79A","PAX5","BLK","IGLL5","ELL2","VASH2","MAD1L1","ZFAT")
mk <- intersect(MARKERS, rownames(GA))
write.csv(round(GA[mk, , drop = FALSE], 4), "marker_genescore_by_atac_label.csv")
cat("\ncanonical annotation markers found in the gene score matrix:",
    length(mk), "of", length(MARKERS), "\n")

cat("\nDONE\n")
