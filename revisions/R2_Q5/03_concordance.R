###############################################################################
# STEP 3 -- how far do the two annotations agree, cell by cell?
#
# Joins the chromatin label (step 1) to the transferred transcriptomic label
# (step 2) and scores the agreement several ways, because "do they overlap?" has
# more than one honest answer:
#
#   raw agreement   fraction of ATAC cells whose chromatin label equals the RNA
#                   label transferred onto them.  Easy to read, but flattered by
#                   the fact that 43% of the cohort is malignant.
#   Cohen's kappa   the same thing corrected for what agreement by chance alone
#                   would give at this composition.  The number to quote.
#   ARI             partition agreement ignoring the names -- would still be high
#                   if the two annotations split the cells identically but called
#                   the groups different things.
#   per class       recall and precision per cell type: which types agree and which
#                   do not, which is what a reader needs.
#   compartment     agreement after collapsing to malignant / epithelial / stromal /
#                   myeloid / lymphoid, i.e. how much disagreement is between
#                   neighbouring subtypes rather than across the tissue.
#   by score        agreement as a function of the transfer's own confidence.  If
#                   the disagreements are concentrated at low predictedScore they
#                   are transfer failures, not annotation failures.
#   by cluster      the chromatin annotation was assigned to 59 LSI CLUSTERS, not
#                   to cells.  Cluster-level purity of the transferred labels shows
#                   whether disagreement is scattered or whole clusters.
#
# P23 is reported separately throughout: it has no matched scRNA, so its cells can
# only borrow a label from another patient's tumour.
#
# Input : atac_cells.csv, atac_label_transfer.csv
# Output: concordance_summary.csv, confusion_matrix.csv, concordance_by_celltype.csv,
#         concordance_by_score.csv, concordance_by_cluster.csv, atac_joined.csv
###############################################################################
source("00_common.R")

A  <- read.csv("atac_cells.csv", stringsAsFactors = FALSE)
TR <- read.csv("atac_label_transfer.csv", stringsAsFactors = FALSE)
X  <- merge(A, TR, by = "cell")
stopifnot(nrow(X) == nrow(A))
X$agree <- X$atac_label == X$predictedGroup
X$comp_atac <- COMPARTMENT[X$atac_label]
X$comp_rna  <- COMPARTMENT[X$predictedGroup]
X$comp_agree <- X$comp_atac == X$comp_rna
X$matched_rna <- X$sample %in% SAMPLES_SHARED
write.csv(X, "atac_joined.csv", row.names = FALSE)

## ---- headline ----------------------------------------------------------------
sets <- list(`all ATAC cells` = rep(TRUE, nrow(X)),
             `patients with matched scRNA` = X$matched_rna,
             `P23 (ATAC only)` = !X$matched_rna,
             `predictedScore >= 0.5` = X$predictedScore >= 0.5)
S <- do.call(rbind, lapply(names(sets), function(k){
  i <- sets[[k]]
  data.frame(subset = k, n_cells = sum(i),
             agreement = mean(X$agree[i]),
             kappa = kappa2(X$atac_label[i], X$predictedGroup[i]),
             ARI = ARI(X$atac_label[i], X$predictedGroup[i]),
             compartment_agreement = mean(X$comp_agree[i], na.rm = TRUE),
             median_score = median(X$predictedScore[i]),
             stringsAsFactors = FALSE) }))
write.csv(S, "concordance_summary.csv", row.names = FALSE)
cat("=== agreement between the chromatin annotation and the transferred RNA annotation ===\n")
print(transform(S, agreement = round(agreement, 4), kappa = round(kappa, 4),
                ARI = round(ARI, 4), compartment_agreement = round(compartment_agreement, 4),
                median_score = round(median_score, 3)), row.names = FALSE)

## ---- confusion matrix ----------------------------------------------------------
lv <- intersect(CELLTYPES, union(X$atac_label, X$predictedGroup))
CF <- table(factor(X$atac_label, lv), factor(X$predictedGroup, lv))
write.csv(as.data.frame.matrix(CF), "confusion_matrix.csv")
cat("\n=== confusion matrix: chromatin label (rows) vs transferred RNA label (cols), row % ===\n")
print(round(100 * prop.table(CF, 1), 1))

## ---- per cell type --------------------------------------------------------------
CT <- class_metrics(X$atac_label, X$predictedGroup)
CT$median_score <- vapply(CT$celltype, function(k)
  median(X$predictedScore[X$atac_label == k]), numeric(1))
CT$top_disagreement <- vapply(CT$celltype, function(k){
  d <- X$predictedGroup[X$atac_label == k & !X$agree]
  if (!length(d)) NA_character_ else names(sort(table(d), decreasing = TRUE))[1] },
  character(1))
CT <- CT[order(-CT$n_atac), ]
write.csv(CT, "concordance_by_celltype.csv", row.names = FALSE)
cat("\n=== per chromatin-defined cell type ===\n")
print(transform(CT, recall = round(recall, 3), precision = round(precision, 3),
                f1 = round(f1, 3), median_score = round(median_score, 3)), row.names = FALSE)

## ---- agreement vs the transfer's own confidence ----------------------------------
X$score_bin <- cut(X$predictedScore, c(-Inf, .2, .4, .5, .6, .8, Inf),
                   labels = c("<0.2","0.2-0.4","0.4-0.5","0.5-0.6","0.6-0.8",">=0.8"))
SB <- do.call(rbind, lapply(split(X, X$score_bin), function(d)
  if (!nrow(d)) NULL else data.frame(score_bin = as.character(d$score_bin[1]),
    n_cells = nrow(d), pct_of_cells = 100 * nrow(d) / nrow(X),
    agreement = mean(d$agree), compartment_agreement = mean(d$comp_agree, na.rm = TRUE),
    stringsAsFactors = FALSE)))
write.csv(SB, "concordance_by_score.csv", row.names = FALSE)
cat("\n=== agreement by transfer confidence ===\n")
print(transform(SB, pct_of_cells = round(pct_of_cells, 1), agreement = round(agreement, 3),
                compartment_agreement = round(compartment_agreement, 3)), row.names = FALSE)

## ---- cluster level ----------------------------------------------------------------
## The chromatin annotation was made per LSI cluster, so this is the level at which
## the two annotations can actually be said to disagree.
top <- function(v) names(sort(table(v), decreasing = TRUE))[1]
pur <- function(v) max(table(v)) / length(v)
CL <- do.call(rbind, lapply(split(X, X$cluster), function(d)
  data.frame(cluster = d$cluster[1], n_cells = nrow(d),
             atac_label = top(d$atac_label), atac_purity = pur(d$atac_label),
             rna_label = top(d$predictedGroup), rna_purity = pur(d$predictedGroup),
             agree = top(d$atac_label) == top(d$predictedGroup),
             cell_agreement = mean(d$agree), median_score = median(d$predictedScore),
             n_samples = length(unique(d$sample)), stringsAsFactors = FALSE)))
CL <- CL[order(-CL$n_cells), ]
write.csv(CL, "concordance_by_cluster.csv", row.names = FALSE)
cat(sprintf("\n=== LSI clusters: majority label agrees in %d of %d clusters (%.1f%% of cells) ===\n",
            sum(CL$agree), nrow(CL), 100 * sum(CL$n_cells[CL$agree]) / sum(CL$n_cells)))
cat("\nclusters where the majority labels disagree:\n")
print(transform(CL[!CL$agree, ], atac_purity = round(atac_purity, 2),
                rna_purity = round(rna_purity, 2), cell_agreement = round(cell_agreement, 2),
                median_score = round(median_score, 3)), row.names = FALSE)

cat("\nDONE\n")
