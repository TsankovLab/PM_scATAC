###############################################################################
# STEP 1 -- export both annotations and document where each one came from.
#
# This is the "document the data" step.  It writes flat tables of every labelled
# cell in each modality, plus a provenance summary, so the comparison in steps 3-6
# never has to reach back into the ArchR project or the Seurat object.
#
# Input : main/scatac_ArchR (ArchR project), main/scrna/srt.rds (Seurat)
# Output: atac_cells.csv, rna_cells.csv, annotation_inventory.csv
###############################################################################
source("00_common.R")
suppressMessages({ library(ArchR); library(Seurat) }); addArchRThreads(1)

## ---- scATAC ------------------------------------------------------------------
proj <- loadArchRProject(ARCHR, showLogo = FALSE)
cd   <- getCellColData(proj)
A <- data.frame(cell    = rownames(cd),
                sample  = as.character(cd$Sample),
                cluster = as.character(cd$Clusters),
                atac_label = as.character(cd[[ATAC_LABEL]]),
                nFrags  = cd$nFrags,
                TSSEnrichment = cd$TSSEnrichment,
                stringsAsFactors = FALSE)
## the UMAP the labels were drawn on, so the figures can reuse it
um <- getEmbedding(proj, embedding = "UMAP", returnDF = TRUE)
A$UMAP1 <- um[A$cell, 1]; A$UMAP2 <- um[A$cell, 2]
write.csv(A, "atac_cells.csv", row.names = FALSE)
cat("scATAC:", nrow(A), "cells |", length(unique(A$sample)), "patients |",
    length(unique(A$cluster)), "LSI clusters |",
    length(unique(A$atac_label)), "cell types\n")
print(table(A$atac_label))

## ---- scRNA -------------------------------------------------------------------
s <- readRDS(SRNA); m <- s@meta.data
R <- data.frame(cell = rownames(m),
                sample = as.character(m$sampleID),
                rna_label = as.character(m[[RNA_LABEL]]),
                rna_label_fine = as.character(m$celltype),
                nCount_RNA = m$nCount_RNA, nFeature_RNA = m$nFeature_RNA,
                stringsAsFactors = FALSE)
write.csv(R, "rna_cells.csv", row.names = FALSE)
cat("\nscRNA :", nrow(R), "cells |", length(unique(R$sample)), "patients |",
    length(unique(R$rna_label)), "cell types\n")
print(table(R$rna_label))

## ---- what overlaps -----------------------------------------------------------
cat("\npatients:  ATAC only:", paste(setdiff(A$sample, R$sample), collapse = ", "),
    "| RNA only:", paste(setdiff(R$sample, A$sample), collapse = ", "),
    "| shared:", length(intersect(unique(A$sample), unique(R$sample))), "\n")
cat("labels  :  ATAC only:", paste(setdiff(A$atac_label, R$rna_label), collapse = ", "),
    "| RNA only:", paste(setdiff(R$rna_label, A$atac_label), collapse = ", "), "\n")
## the two modalities are separate assays on the same tumours, NOT a multiome:
cat("cell barcodes shared between the modalities:",
    length(intersect(sub(".*#", "", A$cell), R$cell)), "\n")

inv <- merge(
  aggregate(cell ~ sample, A, length), aggregate(cell ~ sample, R, length),
  by = "sample", all = TRUE, suffixes = c("_atac", "_rna"))
names(inv) <- c("sample", "n_atac", "n_rna")
inv$modalities <- ifelse(is.na(inv$n_rna), "ATAC only",
                  ifelse(is.na(inv$n_atac), "RNA only", "both"))
write.csv(inv, "annotation_inventory.csv", row.names = FALSE)
cat("\n=== cells per patient ===\n"); print(inv, row.names = FALSE)

cat("\nDONE\n")
