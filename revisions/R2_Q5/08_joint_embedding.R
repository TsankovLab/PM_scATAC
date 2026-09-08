###############################################################################
# STEP 8 -- a joint CCA embedding of the scATAC and scRNA cells.
#
# Steps 3-5 score the agreement numerically.  This puts both modalities in ONE
# space so the agreement can be seen, and so it can be measured in a way the label
# transfer cannot fake: whether cells of the two assays actually mix.
#
# Procedure (the Seurat/Signac cross-modality co-embedding):
#   1. CCA anchors between the RNA reference (expression) and the ATAC cells
#      (gene score) -- the same anchors step 2 uses for the label transfer.
#   2. TransferData with the RNA EXPRESSION MATRIX as refdata, weighted by the ATAC
#      LSI, giving every ATAC cell an imputed expression profile over the shared
#      variable genes.  The ATAC cell is now described in RNA units.
#   3. Merge the imputed ATAC profiles with the real RNA profiles, centre, PCA, UMAP.
#      One embedding, both assays.
#
# What the embedding can and cannot show.  Imputation pulls ATAC cells toward the
# reference by construction, so "the two modalities overlap" is NOT evidence on its
# own -- a co-embedding will look mixed even for a bad transfer.  The measurements
# below are therefore about STRUCTURE, not overlap:
#   modality mixing   fraction of each cell's k nearest neighbours drawn from the
#                     other assay, per cell type.  Compared against the fraction
#                     expected from the cohort composition alone.
#   cross-modality    for each ATAC cell, the majority CHROMATIN-independent label
#   label agreement   among its nearest RNA neighbours, scored against its own
#                     chromatin label.  This is the joint-space analogue of step 3
#                     and should reproduce it if the embedding is honest.
#   centroid geometry per cell type, the distance between the two assays' centroids
#                     relative to the spread between cell types.  Small within,
#                     large between = the joint space is organised by cell type and
#                     not by assay.
#   dispersion        the spread of each assay's cells around its own centroid.
#                     Imputation is a weighted average over anchors, so imputed
#                     profiles are SHRUNK toward the local reference mean and the
#                     ATAC cloud is denser than the RNA one.  This is measured
#                     because it is the reason the neighbour-mixing statistic and
#                     the centroid statistic disagree, and reporting one without
#                     the other would misrepresent the embedding.
#
# Runtime ~1-2 h.  Submit with submit_08_joint_embedding.sh.
#
# Input : main/scatac_ArchR, main/scrna/srt.rds, atac_cells.csv, rna_cells.csv
# Output: joint_embedding.csv, joint_mixing_by_celltype.csv,
#         joint_label_agreement.csv, joint_centroid_distance.csv,
#         joint_dispersion.csv
###############################################################################
source("00_common.R")
suppressMessages({ library(ArchR); library(Seurat); library(Matrix)
                   library(SummarizedExperiment) })
addArchRThreads(1)
t0 <- Sys.time()

NREF_PER_CT <- 3000   # reference cells per RNA cell type for the anchors (as step 2)
NFEAT       <- 2000   # shared variable genes
NDIM        <- 30     # CCA / LSI / PCA dimensions
BLOCK       <- 10000  # ATAC cells imputed at a time
KNN         <- 30     # neighbours for the mixing and label statistics

## ---- ATAC gene score + LSI ----------------------------------------------------
## See step 2 on why the gene score matrix is only ever indexed by position.
proj <- loadArchRProject(ARCHR, showLogo = FALSE)
LSI  <- getReducedDims(proj, "IterativeLSI", corCutOff = 0.75, dimsToUse = seq_len(NDIM))
se   <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", verbose = FALSE)
rm(proj); invisible(gc(FALSE))
GS    <- assay(se)
gname <- as.character(rowData(se)$name)
acell <- colnames(se)
gkeep <- which(!duplicated(gname)); gname <- gname[gkeep]
cat("ATAC:", length(acell), "cells |", length(gkeep), "genes | LSI", ncol(LSI), "dims\n")

## ---- RNA: full object for the embedding, subsample for the anchors -------------
s     <- readRDS(SRNA)
rlab  <- as.character(s@meta.data[[RNA_LABEL]])
rsamp <- as.character(s@meta.data$sampleID)
RNAc  <- SeuratObject::LayerData(s[["RNA"]], layer = "counts")
rcell <- colnames(s); rm(s); invisible(gc(FALSE))

set.seed(1)
sub <- unlist(lapply(split(seq_along(rlab), rlab), function(i)
  if (length(i) <= NREF_PER_CT) i else sample(i, NREF_PER_CT)))
ref <- CreateSeuratObject(counts = RNAc[, sub],
                          meta.data = data.frame(celltype = rlab[sub],
                                                 row.names = rcell[sub]))
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, nfeatures = NFEAT, verbose = FALSE)
feats <- intersect(VariableFeatures(ref), gname)
ref <- ScaleData(ref, features = feats, verbose = FALSE)
cat("RNA reference for anchors:", ncol(ref), "cells |", length(feats),
    "shared variable genes\n")

## the RNA half of the co-embedding: every RNA cell, log-normalised, shared genes
rna_all <- CreateSeuratObject(counts = RNAc)
rna_all <- NormalizeData(rna_all, verbose = FALSE)
RD <- SeuratObject::LayerData(rna_all[["RNA"]], layer = "data")[feats, , drop = FALSE]
rm(rna_all, RNAc); invisible(gc(FALSE))
cat("RNA for the embedding:", ncol(RD), "cells\n")

## ---- CCA anchors + imputation, in blocks --------------------------------------
refdata <- SeuratObject::LayerData(ref[["RNA"]], layer = "data")[feats, , drop = FALSE]
blocks  <- split(seq_along(acell), ceiling(seq_along(acell) / BLOCK))
AD <- matrix(0, nrow = length(feats), ncol = length(acell),
             dimnames = list(feats, acell))
for (b in seq_along(blocks)){
  i <- blocks[[b]]
  m <- GS[gkeep, i, drop = FALSE]
  m <- as(as.matrix(m), "CsparseMatrix"); dimnames(m) <- list(gname, acell[i])
  q <- CreateSeuratObject(counts = m, assay = "GeneScore"); rm(m)
  q <- NormalizeData(q, verbose = FALSE)
  q <- ScaleData(q, features = feats, verbose = FALSE)
  an <- FindTransferAnchors(reference = ref, query = q, reduction = "cca",
                            features = feats, dims = seq_len(NDIM), verbose = FALSE)
  wr <- CreateDimReducObject(embeddings = LSI[acell[i], , drop = FALSE],
                             key = "LSI_", assay = "GeneScore")
  im <- TransferData(anchorset = an, refdata = refdata, weight.reduction = wr,
                     dims = seq_len(ncol(LSI)), verbose = FALSE)
  AD[, i] <- as.matrix(SeuratObject::LayerData(im, layer = "data"))[feats, ]
  cat(sprintf("  block %d/%d imputed, %d anchors, %.1f min\n", b, length(blocks),
              nrow(an@anchors), difftime(Sys.time(), t0, units = "mins")))
  rm(q, an, wr, im); invisible(gc(FALSE))
}
rm(GS, se); invisible(gc(FALSE))

## ---- co-embed -------------------------------------------------------------------
## Centre only (do.scale = FALSE): imputed values are already on the reference's
## log-normalised scale, and unit-variance scaling would inflate genes the
## imputation left nearly constant.
J <- cbind(as.matrix(RD), AD)
modality <- c(rep("scRNA", ncol(RD)), rep("scATAC", ncol(AD)))
rm(RD, AD); invisible(gc(FALSE))
cat("joint matrix:", nrow(J), "genes x", ncol(J), "cells\n")
J <- J - rowMeans(J)
pca <- irlba::prcomp_irlba(t(J), n = NDIM, center = FALSE, scale. = FALSE)
PC  <- pca$x; rownames(PC) <- colnames(J); colnames(PC) <- paste0("PC", seq_len(NDIM))
rm(J); invisible(gc(FALSE))
cat("PCA done,", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

um <- uwot::umap(PC, n_neighbors = 30, min_dist = 0.3, metric = "cosine",
                 n_threads = 4, verbose = FALSE)
E <- data.frame(cell = rownames(PC), modality = modality,
                UMAP1 = um[, 1], UMAP2 = um[, 2], stringsAsFactors = FALSE)

## labels: chromatin label for ATAC cells, expression label for RNA cells
A <- read.csv("atac_cells.csv", stringsAsFactors = FALSE)
R <- read.csv("rna_cells.csv",  stringsAsFactors = FALSE)
E$celltype <- ifelse(E$modality == "scATAC",
                     A$atac_label[match(E$cell, A$cell)],
                     R$rna_label[match(E$cell, R$cell)])
E$sample   <- ifelse(E$modality == "scATAC",
                     A$sample[match(E$cell, A$cell)],
                     R$sample[match(E$cell, R$cell)])
write.csv(E, "joint_embedding.csv", row.names = FALSE)
cat("joint embedding written,", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

## ---- 1. modality mixing ---------------------------------------------------------
nn  <- FNN::get.knn(PC, k = KNN)$nn.index
oth <- matrix(modality[nn], nrow = nrow(nn)) != modality
E$frac_other_modality <- rowMeans(oth)
exp_atac <- mean(modality == "scATAC")   # what a perfectly mixed space would give
MX <- do.call(rbind, lapply(split(E, list(E$celltype, E$modality), drop = TRUE), function(d)
  data.frame(celltype = d$celltype[1], modality = d$modality[1], n = nrow(d),
             frac_other = mean(d$frac_other_modality),
             expected = if (d$modality[1] == "scATAC") 1 - exp_atac else exp_atac,
             stringsAsFactors = FALSE)))
MX$mixing_ratio <- MX$frac_other / MX$expected
MX <- MX[order(MX$celltype, MX$modality), ]
write.csv(MX, "joint_mixing_by_celltype.csv", row.names = FALSE)
cat(sprintf("\n=== modality mixing (%d nearest neighbours in the joint PCA) ===\n", KNN))
cat(sprintf("cohort is %.1f%% scATAC / %.1f%% scRNA; ratio 1.0 = as mixed as composition allows\n",
            100 * exp_atac, 100 * (1 - exp_atac)))
print(transform(MX, frac_other = round(frac_other, 3), expected = round(expected, 3),
                mixing_ratio = round(mixing_ratio, 2)), row.names = FALSE)
cat(sprintf("\noverall: %.3f of neighbours from the other assay (expected %.3f)\n",
            mean(E$frac_other_modality),
            2 * exp_atac * (1 - exp_atac) / (exp_atac^2 + (1 - exp_atac)^2 +
                                             2 * exp_atac * (1 - exp_atac))))

## ---- 2. cross-modality label agreement in the joint space -------------------------
## For each ATAC cell, the commonest RNA cell type among its nearest RNA neighbours.
ri <- which(modality == "scRNA"); ai <- which(modality == "scATAC")
rnn <- FNN::get.knnx(PC[ri, ], PC[ai, ], k = KNN)$nn.index
rl  <- E$celltype[ri]
maj <- apply(matrix(rl[rnn], nrow = nrow(rnn)), 1,
             function(x) names(sort(table(x), decreasing = TRUE))[1])
JA <- data.frame(cell = E$cell[ai], atac_label = E$celltype[ai],
                 joint_rna_label = maj, stringsAsFactors = FALSE)
JA$agree <- JA$atac_label == JA$joint_rna_label
write.csv(JA, "joint_label_agreement.csv", row.names = FALSE)
cat(sprintf("\n=== joint-space label agreement: %.3f of ATAC cells (kappa %.3f, ARI %.3f) ===\n",
            mean(JA$agree), kappa2(JA$atac_label, JA$joint_rna_label),
            ARI(JA$atac_label, JA$joint_rna_label)))
print(round(100 * prop.table(table(JA$atac_label, JA$joint_rna_label), 1), 1))

## ---- 3. centroid geometry ----------------------------------------------------------
cen <- do.call(rbind, lapply(split(seq_len(nrow(PC)), paste(E$celltype, E$modality)),
                             function(i) colMeans(PC[i, , drop = FALSE])))
key <- do.call(rbind, strsplit(rownames(cen), " (?=[^ ]+$)", perl = TRUE))
DM  <- as.matrix(dist(cen))
cts <- intersect(CELLTYPES, unique(key[, 1]))
CD  <- do.call(rbind, lapply(cts, function(k){
  a <- which(key[, 1] == k & key[, 2] == "scATAC")
  r <- which(key[, 1] == k & key[, 2] == "scRNA")
  if (!length(a) || !length(r)) return(NULL)
  oth <- which(key[, 1] != k)
  data.frame(celltype = k, within_celltype_across_assay = DM[a, r],
             median_to_other_celltypes = median(DM[a, oth]),
             ratio = DM[a, r] / median(DM[a, oth]), stringsAsFactors = FALSE) }))
write.csv(CD, "joint_centroid_distance.csv", row.names = FALSE)
cat("\n=== centroid geometry: same cell type across assays vs across cell types ===\n")
cat("ratio < 1 means the two assays' cells of a type sit closer to each other than to other types\n")
print(transform(CD, within_celltype_across_assay = round(within_celltype_across_assay, 2),
                median_to_other_celltypes = round(median_to_other_celltypes, 2),
                ratio = round(ratio, 3)), row.names = FALSE)
cat(sprintf("\ncell types where the assays sit closer to each other than to other types: %d of %d\n",
            sum(CD$ratio < 1), nrow(CD)))

## ---- 4. dispersion: why 2 and 3 disagree ------------------------------------------
## RMS distance of each assay's cells to their own (cell type, assay) centroid.
disp <- function(i) sqrt(mean(rowSums((PC[i, , drop = FALSE] -
                       rep(colMeans(PC[i, , drop = FALSE]), each = length(i)))^2)))
SP <- do.call(rbind, lapply(cts, function(k){
  a <- which(E$celltype == k & modality == "scATAC")
  r <- which(E$celltype == k & modality == "scRNA")
  if (length(a) < 5 || length(r) < 5) return(NULL)
  data.frame(celltype = k, n_atac = length(a), n_rna = length(r),
             spread_atac = disp(a), spread_rna = disp(r),
             shrinkage = disp(a) / disp(r), stringsAsFactors = FALSE) }))
write.csv(SP, "joint_dispersion.csv", row.names = FALSE)
cat("\n=== dispersion of each assay around its own centroid ===\n")
cat("shrinkage < 1 means the imputed ATAC cells are more tightly packed than the real RNA cells\n")
print(transform(SP, spread_atac = round(spread_atac, 2), spread_rna = round(spread_rna, 2),
                shrinkage = round(shrinkage, 3)), row.names = FALSE)
cat(sprintf("median shrinkage across cell types: %.3f\n", median(SP$shrinkage)))

## how far is an ATAC cell from the RNA manifold, in units of RNA-to-RNA spacing?
d_ar <- FNN::get.knnx(PC[ri, ], PC[ai, ], k = 1)$nn.dist[, 1]
d_rr <- FNN::get.knnx(PC[ri, ], PC[ri, ], k = 2)$nn.dist[, 2]
cat(sprintf("\nnearest RNA cell: median distance from an ATAC cell %.2f, from another RNA cell %.2f (ratio %.2f)\n",
            median(d_ar), median(d_rr), median(d_ar) / median(d_rr)))

## mixing on a size-balanced subsample, so cohort composition is not a confound
set.seed(1)
bal <- unlist(lapply(cts, function(k){
  a <- which(E$celltype == k & modality == "scATAC")
  r <- which(E$celltype == k & modality == "scRNA")
  n <- min(length(a), length(r)); if (n < 20) return(integer(0))
  c(sample(a, n), sample(r, n)) }))
if (length(bal) > 2 * KNN){
  nb <- FNN::get.knn(PC[bal, ], k = KNN)$nn.index
  mb <- modality[bal]
  cat(sprintf("balanced subsample (%d cells, 50/50 per cell type): %.3f of neighbours from the other assay (0.5 = perfect)\n",
              length(bal), mean(matrix(mb[nb], nrow = nrow(nb)) != mb)))
}

cat("\nelapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\nDONE\n")
