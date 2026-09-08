###############################################################################
# STEP 2 -- transfer the transcriptomic labels onto the scATAC cells.
#
# The two modalities are separate assays on the same tumours (step 1 confirms the
# barcodes do not match, so this is not a multiome).  The only way to ask whether
# the chromatin annotation and the transcriptomic annotation agree CELL BY CELL is
# to bring the RNA labels across:
#
#   ArchR GeneScoreMatrix  --(Seurat CCA anchors)-->  every ATAC cell gets
#       predictedGroup  the RNA celltype_lv1 it looks most like
#       predictedScore  prediction.score.max, the anchor-weighted confidence (0-1)
#
# The transfer is UNCONSTRAINED: it is never told the chromatin label of the ATAC
# cell, nor the patient of either cell.  That is the point -- a constrained transfer
# would be handed the answer and could not test agreement.
#
# This is the same procedure as ArchR's addGeneIntegrationMatrix (gene score as the
# query feature space, CCA anchors, LSI as the transfer weight reduction), written
# out directly because ArchR 1.0.2 predates Seurat 5 and its internal Seurat object
# handling fails against SeuratObject 5.x ("more elements supplied than there are to
# replace", from an Assay5 created inside the block loop).  Writing it here also
# makes every parameter visible instead of buried in package defaults.
#
# NOTHING IS WRITTEN BACK INTO THE SHARED ArchR PROJECT -- the project is only read.
#
# Runtime ~1-3 h.  Submit with submit_02_label_transfer.sh.
#
# Input : main/scatac_ArchR, main/scrna/srt.rds
# Output: atac_label_transfer.csv (one row per ATAC cell)
###############################################################################
source("00_common.R")
suppressMessages({ library(ArchR); library(Seurat); library(Matrix)
                   library(SummarizedExperiment) })
addArchRThreads(1)
t0 <- Sys.time()

NREF_PER_CT <- 3000    # reference cells kept per RNA cell type (all of them if fewer)
NFEAT       <- 2000    # variable genes defining the shared feature space
NDIM        <- 30      # CCA / LSI dimensions
BLOCK       <- 10000   # ATAC cells anchored at a time, to cap peak memory

## ---- ATAC gene score + LSI ----------------------------------------------------
## NB on the gene score matrix: ArchR returns a dgCMatrix whose row names are not
## the gene symbols (they are in rowData).  Do NOT rebuild it with sparseMatrix()
## to attach them -- that produces an object whose `[` no longer dispatches
## ("object of type 'S4' is not subsettable").  Index the original BY POSITION,
## carry the symbols in a plain character vector, and name the small per-block
## matrix instead.  Matrix must be attached for the S4 subsetting methods.
proj <- loadArchRProject(ARCHR, showLogo = FALSE)
LSI  <- getReducedDims(proj, "IterativeLSI", corCutOff = 0.75,
                       dimsToUse = seq_len(NDIM))
cat("ATAC:", nCells(proj), "cells | LSI", ncol(LSI), "dims kept\n")
se <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", verbose = FALSE)
rm(proj); invisible(gc(FALSE))
GS     <- assay(se)
gname  <- as.character(rowData(se)$name)
cells  <- colnames(se)
gkeep  <- which(!duplicated(gname))          # integer positions, first row per symbol
gname  <- gname[gkeep]
cat("gene score:", length(gkeep), "genes x", length(cells), "cells\n")

## ---- RNA reference, subsampled per cell type ----------------------------------
## Subsampling is stratified so the rare types (Alveolar 147, Glia 55, Plasma 179)
## keep every cell while the large ones are capped -- an unstratified subsample
## would drop the rare types below the anchor threshold entirely.
s   <- readRDS(SRNA)
lab <- as.character(s@meta.data[[RNA_LABEL]])
set.seed(1)
keep <- unlist(lapply(split(seq_along(lab), lab), function(i)
  if (length(i) <= NREF_PER_CT) i else sample(i, NREF_PER_CT)))
ref <- CreateSeuratObject(counts = SeuratObject::LayerData(s[["RNA"]], layer = "counts")[, keep],
                          meta.data = data.frame(celltype = lab[keep],
                                                 row.names = colnames(s)[keep]))
rm(s); invisible(gc(FALSE))
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, nfeatures = NFEAT, verbose = FALSE)
cat("RNA reference:", ncol(ref), "cells\n"); print(table(ref$celltype))

## ---- shared feature space -------------------------------------------------------
feats <- intersect(VariableFeatures(ref), gname)
cat("shared variable features:", length(feats), "of", NFEAT, "\n")
ref <- ScaleData(ref, features = feats, verbose = FALSE)

## ---- transfer, in blocks ---------------------------------------------------------
blocks <- split(seq_along(cells), ceiling(seq_along(cells) / BLOCK))
cat("transferring in", length(blocks), "blocks of up to", BLOCK, "cells\n")
res <- lapply(seq_along(blocks), function(b){
  i <- blocks[[b]]
  m <- GS[gkeep, i, drop = FALSE]            # positional only
  m <- as(as.matrix(m), "CsparseMatrix")     # a clean matrix Seurat can name
  dimnames(m) <- list(gname, cells[i])
  q <- CreateSeuratObject(counts = m, assay = "GeneScore"); rm(m)
  q <- NormalizeData(q, verbose = FALSE)
  q <- ScaleData(q, features = feats, verbose = FALSE)
  an <- FindTransferAnchors(reference = ref, query = q, reduction = "cca",
                            features = feats, dims = seq_len(NDIM), verbose = FALSE)
  wr <- CreateDimReducObject(embeddings = LSI[cells[i], , drop = FALSE],
                             key = "LSI_", assay = "GeneScore")
  pr <- TransferData(anchorset = an, refdata = ref$celltype,
                     weight.reduction = wr, dims = seq_len(ncol(LSI)), verbose = FALSE)
  cat(sprintf("  block %d/%d: %d cells, %d anchors, %.1f min elapsed\n", b, length(blocks),
              length(i), nrow(an@anchors), difftime(Sys.time(), t0, units = "mins")))
  data.frame(cell = cells[i], predictedGroup = as.character(pr$predicted.id),
             predictedScore = pr$prediction.score.max, stringsAsFactors = FALSE)
})
TR <- do.call(rbind, res)
write.csv(TR, "atac_label_transfer.csv", row.names = FALSE)

cat("\ntransferred labels:\n"); print(table(TR$predictedGroup))
cat("\npredictedScore quantiles:\n")
print(round(quantile(TR$predictedScore, c(0, .1, .25, .5, .75, .9, 1)), 3))
cat("\nelapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
cat("DONE\n")
