###############################################################################
# R1_Q8 -- peak-overlap heatmap for the 20 sarcomatoid-associated TFs
#
# The boxplot figure ranks 20 TFs by how well their motif activity / gene score tracks the
# cNMF20 sarcomatoid score, but says nothing about whether those TFs are independent. Many
# of them bind similar sequences, so a "hit" can simply be the same peaks counted twice.
# This quantifies that directly from the ArchR motif matches (peaks x motifs, binary):
#
#   Jaccard   |A n B| / |A u B|   -- symmetric, the headline panel
#   fraction  |A n B| / |A|       -- asymmetric, row TF's peaks that also carry the column TF
#
# Peaks are the project peak set; a peak "carries" a TF if ArchR's motif scan matched it.
# TF order is the published top_sarc_TF order on the main panel, and the clustered order on
# the second, so co-binding blocks are visible without losing the figure's ordering.
#
# Output: sarc_TF_peak_overlap_jaccard.csv, sarc_TF_peak_overlap_counts.csv,
#         Plots/sarc_TF_peak_overlap_heatmap.pdf (+ .png)
###############################################################################
suppressPackageStartupMessages({ library(ArchR); library(ComplexHeatmap); library(circlize)
                                 library(RColorBrewer); library(paletteer) })
addArchRThreads(8); addArchRGenome("hg38")
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
TUM  <- file.path(ROOT, "tumor_compartment", "scatac_ArchR")
OUT  <- file.path(ROOT, "git_repo_claude", "R1_Q8")
setwd(OUT); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "palettes.R"))

top_sarc_TF <- readRDS(file.path(TUM, "top_sarc_TF.rds"))
proj <- loadArchRProject(TUM, showLogo = FALSE)
M    <- getMatches(proj, "Motif")
cn   <- colnames(M)
## motif columns carry a numeric suffix the TF names do not
mcol <- sapply(top_sarc_TF, function(tf) { i <- grep(paste0("^", tf, "_[0-9]+$"), cn); cn[i[1]] })
stopifnot(!any(is.na(mcol)))

B <- assay(M)[, mcol, drop = FALSE]          # peaks x 20, logical
colnames(B) <- top_sarc_TF
B <- as(B, "dgCMatrix")
npk <- Matrix::colSums(B)
cat("peaks in the project:", nrow(B), "\n")
cat("peaks carrying each motif:\n"); print(sort(npk, decreasing = TRUE))

inter <- as.matrix(Matrix::crossprod(B))      # |A n B|
union <- outer(npk, npk, "+") - inter
J     <- inter / union
F     <- inter / npk                          # row-normalised: fraction of ROW's peaks
diag(J) <- 1
write.csv(J,     "sarc_TF_peak_overlap_jaccard.csv")
write.csv(inter, "sarc_TF_peak_overlap_counts.csv")

offdiag <- J[upper.tri(J)]
cat("\nJaccard off-diagonal: median", round(median(offdiag), 3),
    "| max", round(max(offdiag), 3), "\n")
pr <- which(J == max(offdiag), arr.ind = TRUE)[1, ]
cat("most overlapping pair:", rownames(J)[pr[1]], "/", colnames(J)[pr[2]],
    sprintf("(J = %.3f, %d shared peaks)\n", max(offdiag), inter[pr[1], pr[2]]))
cat("\npairs with Jaccard > 0.3:\n")
ij <- which(J > 0.3 & upper.tri(J), arr.ind = TRUE)
if (nrow(ij)) {
  d <- data.frame(TF1 = rownames(J)[ij[, 1]], TF2 = colnames(J)[ij[, 2]],
                  jaccard = J[ij], shared = inter[ij])
  print(d[order(-d$jaccard), ], row.names = FALSE, digits = 3)
} else cat("  none\n")

col_j <- colorRamp2(seq(0, max(offdiag), length.out = 9), brewer.pal(9, "PuBu"))
col_f <- colorRamp2(seq(0, 1, length.out = 9), brewer.pal(9, "OrRd"))
cell_fun <- function(m, digits) function(j, i, x, y, w, h, fill) {
  v <- m[i, j]
  if (i != j) grid.text(sprintf(paste0("%.", digits, "f"), v), x, y,
                        gp = gpar(fontsize = 5, col = if (v > max(m[upper.tri(m)]) * .6) "white" else "grey20"))
}
ha <- rowAnnotation(`peaks (k)` = anno_barplot(npk[top_sarc_TF] / 1000, gp = gpar(fill = "grey60", col = NA),
                                               width = unit(1.6, "cm")),
                    annotation_name_gp = gpar(fontsize = 7))
h1 <- Heatmap(J[top_sarc_TF, top_sarc_TF], name = "Jaccard", col = col_j,
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 7.5), column_names_gp = gpar(fontsize = 7.5),
              column_names_rot = 45, right_annotation = ha,
              cell_fun = cell_fun(J[top_sarc_TF, top_sarc_TF], 2),
              column_title = "Motif peak overlap (Jaccard), published TF order",
              column_title_gp = gpar(fontsize = 9))
h2 <- Heatmap(J, name = "Jaccard ", col = col_j,
              clustering_distance_rows = function(m) as.dist(1 - J),
              clustering_distance_columns = function(m) as.dist(1 - J),
              row_names_gp = gpar(fontsize = 7.5), column_names_gp = gpar(fontsize = 7.5),
              column_names_rot = 45, cell_fun = cell_fun(J, 2),
              column_title = "same, clustered", column_title_gp = gpar(fontsize = 9))
h3 <- Heatmap(F[top_sarc_TF, top_sarc_TF], name = "fraction of\nrow TF peaks", col = col_f,
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 7.5), column_names_gp = gpar(fontsize = 7.5),
              column_names_rot = 45, cell_fun = cell_fun(F[top_sarc_TF, top_sarc_TF], 2),
              column_title = "Asymmetric: fraction of the row TF's peaks carrying the column TF",
              column_title_gp = gpar(fontsize = 9))

pdf("Plots/sarc_TF_peak_overlap_heatmap.pdf", width = 8.5, height = 7.5)
draw(h1); draw(h2); draw(h3)
dev.off()
png("Plots/sarc_TF_peak_overlap_heatmap.png", width = 8.5, height = 7.5, units = "in", res = 200)
draw(h1); dev.off()
cat("\nDONE\n")
