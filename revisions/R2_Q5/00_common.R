###############################################################################
# R2_Q5 -- shared settings for the ATAC-vs-RNA cell type annotation comparison.
#
# Reviewer 2, Question 5: "Does cell type annotation clustered with open chromatin
# accessibility overlap with cell annotation using transcriptomic data?"
#
# Nothing here runs an analysis; it fixes the paths, the label vocabulary and the
# agreement statistics so the same definitions are used by every step.
###############################################################################

## ---- paths ------------------------------------------------------------------
SC    <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
ROOT  <- file.path(SC, "git_repo_claude", "R2_Q5")
ARCHR <- file.path(SC, "main", "scatac_ArchR")      # whole-cohort scATAC project
SRNA  <- file.path(SC, "main", "scrna", "srt.rds")  # whole-cohort scRNA Seurat object
setwd(ROOT)
dir.create("Plots", showWarnings = FALSE)

## ArchR 1.0.2 predates Seurat 5 and mixes CreateSeuratObject() with the v3-only
## CreateAssayObject().  Forcing v3 assays keeps the two consistent; without this
## addGeneIntegrationMatrix() fails inside FindTransferAnchors.
options(Seurat.object.assay.version = "v3")

## ---- what is being compared -------------------------------------------------
## ATAC labels: annotated from CHROMATIN ONLY.  IterativeLSI on the TileMatrix
##   (25 000 variable features, LSIMethod 2) -> addClusters(resolution 3) -> 59
##   clusters -> each cluster named from the gene score of canonical markers
##   (git_repo/main_analysis/scatac_main_ArchR.R).  The scRNA object was never used.
## RNA labels : annotated from EXPRESSION ONLY, in the matched scRNA cohort.
## The two are therefore independent annotations of the same patients, and the
## question is how far they agree.
ATAC_LABEL <- "celltype_lv1"
RNA_LABEL  <- "celltype_lv1"

## the 12 ATAC / 13 RNA level-1 types; Glia exists only in the RNA object
CELLTYPES <- c("Malignant", "Mesothelium", "Alveolar", "Fibroblasts", "SmoothMuscle",
               "Endothelial", "Myeloid", "pDCs", "T_cells", "NK", "B_cells",
               "Plasma", "Glia")

## broad compartments -- used for the "is the disagreement between neighbours or
## across the tissue?" question, never for the headline agreement number
COMPARTMENT <- c(Malignant = "malignant", Mesothelium = "malignant",
                 Alveolar = "epithelial",
                 Fibroblasts = "stromal", SmoothMuscle = "stromal",
                 Endothelial = "stromal", Glia = "stromal",
                 Myeloid = "myeloid", pDCs = "myeloid",
                 T_cells = "lymphoid", NK = "lymphoid",
                 B_cells = "lymphoid", Plasma = "lymphoid")

## ---- samples ----------------------------------------------------------------
## ATAC has 11 patients, RNA has 10; P23 was profiled by ATAC only.  Every
## patient-level statistic is restricted to the 10 shared patients.
SAMPLES_ATAC   <- c("P1","P3","P4","P5","P8","P10","P11","P12","P13","P14","P23")
SAMPLES_SHARED <- c("P1","P3","P4","P5","P8","P10","P11","P12","P13","P14")

## ---- palette (git_repo/utils/palettes.R) ------------------------------------
CTCOL <- c(Malignant = "plum4", Mesothelium = "olivedrab1", Alveolar = "black",
           Glia = "lawngreen", Fibroblasts = "azure4", SmoothMuscle = "blueviolet",
           Endothelial = "brown", Myeloid = "cornflowerblue", pDCs = "tomato",
           T_cells = "firebrick1", NK = "gold1", B_cells = "magenta2",
           Plasma = "lightsalmon1")

## ---- agreement statistics ---------------------------------------------------
## Adjusted Rand index -- partition agreement, ignores the label names entirely.
ARI <- function(a, b){
  tb <- table(a, b); n <- sum(tb)
  ci <- sum(choose(rowSums(tb), 2)); cj <- sum(choose(colSums(tb), 2))
  ex <- ci * cj / choose(n, 2)
  (sum(choose(tb, 2)) - ex) / (0.5 * (ci + cj) - ex)
}

## Cohen's kappa -- agreement of the NAMED labels corrected for chance.  Reported
## alongside raw agreement because the cohort is 43% malignant, so a labeller that
## said "Malignant" every time would already look 43% right.
kappa2 <- function(a, b){
  lv <- union(levels(factor(a)), levels(factor(b)))
  tb <- table(factor(a, lv), factor(b, lv)); n <- sum(tb)
  po <- sum(diag(tb)) / n
  pe <- sum(rowSums(tb) * colSums(tb)) / n^2
  (po - pe) / (1 - pe)
}

## Per-class recall / precision / F1 of the RNA-transferred label against the
## chromatin label, treating the chromatin label as the reference.
class_metrics <- function(ref, pred){
  lv <- sort(union(unique(ref), unique(pred)))
  do.call(rbind, lapply(lv, function(k){
    tp <- sum(ref == k & pred == k); fn <- sum(ref == k & pred != k)
    fp <- sum(ref != k & pred == k)
    rec <- if (tp + fn > 0) tp / (tp + fn) else NA
    pre <- if (tp + fp > 0) tp / (tp + fp) else NA
    data.frame(celltype = k, n_atac = tp + fn, n_rna = tp + fp, tp = tp,
               recall = rec, precision = pre,
               f1 = if (is.na(rec) || is.na(pre) || rec + pre == 0) NA
                    else 2 * rec * pre / (rec + pre), stringsAsFactors = FALSE)
  }))
}

set.seed(1)












