#!/usr/bin/env Rscript
# diag_archr.R — inspect main ArchR project: cellColData cols, Sample values,
# celltype annotation presence, and whether WT1 is in GeneScoreMatrix.
suppressPackageStartupMessages(library(ArchR))
addArchRGenome("hg38"); addArchRThreads(1)

SRC <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR"
proj <- loadArchRProject(SRC, showLogo = FALSE)

message("nCells: ", nCells(proj))
message("cellColData cols: ", paste(colnames(proj@cellColData), collapse=", "))
message("=== Sample table ===");        print(table(proj$Sample))
ct <- grep("celltype|cell_type|annot", colnames(proj@cellColData), value=TRUE, ignore.case=TRUE)
message("celltype-like cols: ", paste(ct, collapse=", "))
for (c in ct) { message("--- ", c, " ---"); print(table(proj@cellColData[[c]])) }
message("=== available matrices ===")
print(getAvailableMatrices(proj))
gsm <- tryCatch(getFeatures(proj, useMatrix="GeneScoreMatrix"), error=function(e) NULL)
if (!is.null(gsm)) message("WT1 in GeneScoreMatrix: ", "WT1" %in% gsm)
