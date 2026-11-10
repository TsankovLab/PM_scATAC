###############################################################################
# R2_Q10 -- keratin expression in epithelioid tumour cells compared with normal mesothelium (scRNA-seq)
#
# Data: main/scrna/srt.rds (58,950 cells, 10 patients). Malignant cells and non-malignant mesothelium are the
# celltype_lv1 labels "Malignant" and "Mesothelium".
#
# Epithelioid vs sarcomatoid is decided with the scS-score, as elsewhere in the manuscript: the sarcomatoid
# programme cNMF20 of the 25-programme malignant cNMF (git_repo/files/malignant_cnmf_genelist_25_nfeat_5000.rds),
# top 20 genes (as in bulk_analysis.R), scored per cell with Seurat AddModuleScore. No keratin is among its genes,
# so the grouping is not built from the genes being tested.
#   primary     tumour-level call: a tumour is epithelioid when the median scS of its malignant cells is below
#               SCS_CUT. Five tumours form a tight low-scS cluster (medians <= 0.07); the next ones start at 0.29.
#               SCS_CUT is the middle of that first break. Epithelioid cells = malignant cells of those tumours.
#               The tumours above the cut are a mix of intermediate and high scS and are labelled "higher scS".
#   sensitivity cell-level call: any malignant cell with scS < SCS_CUT, whatever its tumour.
# The per-cell scS distribution is unimodal, so there is no natural per-cell threshold; that is why the tumour-level
# call is the primary one.
#
# Keratins: every KRT gene detected in >= 5% of the cells of at least one group.
# Comparison: epithelioid malignant cells vs mesothelium -- per gene the % of cells expressing, the mean
# log-normalised expression, log2 fold change of the means, and a Wilcoxon test over cells (BH over the keratins).
# IMPORTANT LIMIT: the mesothelium cells come from very few patients (see the printed table), and those patients are
# not the epithelioid ones, so the comparison is between patients and cells are not independent replicates. A
# pseudobulk per patient is therefore shown next to the cell-level result, and the sample-level test is reported
# for what it is.
#
# Output: R2_Q10_scS_per_sample.csv, R2_Q10_keratin_stats.csv, R2_Q10_keratin_pseudobulk_per_sample.csv
#         (the first two are inputs of R2_Q10_keratins_scATAC_genescore.R, run it after this script)
###############################################################################
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
options(width = 220); set.seed(1)
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"; OUT <- file.path(ROOT, "git_repo_claude", "R2_Q10")
setwd(OUT)
N_SCS <- 20; MIN_PCT <- 0.05; MIN_CELLS <- 20

srt <- readRDS(file.path(ROOT, "main", "scrna", "srt.rds"))
srt <- srt[, srt$celltype_lv1 %in% c("Malignant", "Mesothelium")]
mods <- readRDS(file.path(ROOT, "git_repo", "files", "malignant_cnmf_genelist_25_nfeat_5000.rds")); sg <- intersect(head(mods[["cNMF20"]], N_SCS), rownames(srt))
cat("scS genes used (", length(sg), "of", N_SCS, "):", paste(sg, collapse = ", "), "\n")
srt <- AddModuleScore(srt, features = list(sg), name = "scS", seed = 1); md <- srt@meta.data; md$scS <- md$scS1; md$sample <- as.character(md$sampleID)

## ---- tumour-level epithelioid call -----------------------------------------------------------
mal <- md[md$celltype_lv1 == "Malignant", ]; ss <- data.frame(sample = names(table(mal$sample)), n_malignant = as.integer(table(mal$sample)))
ss$median_scS <- tapply(mal$scS, mal$sample, median)[ss$sample]; ss <- ss[order(ss$median_scS), ]
## cut at the FIRST clear break above the low-scS cluster: going up the ordered medians, the first gap more than twice as wide as any
## gap below it. (The single largest gap is not used: it only isolates the most sarcomatoid tumour.)
gaps <- diff(ss$median_scS); brk <- which(sapply(seq_along(gaps), function(i) i > 1 && gaps[i] > 2 * max(gaps[seq_len(i - 1)])))[1]
SCS_CUT <- mean(ss$median_scS[brk + 0:1]); ss$call <- ifelse(ss$median_scS < SCS_CUT, "epithelioid", "higher scS")
ss$n_mesothelium <- as.integer(table(factor(md$sample[md$celltype_lv1 == "Mesothelium"], levels = ss$sample)))
ss$pct_cells_below_cut <- round(100 * tapply(mal$scS < SCS_CUT, mal$sample, mean)[ss$sample], 1)
write.csv(ss, "R2_Q10_scS_per_sample.csv", row.names = FALSE)
cat(sprintf("\nscS cut-off (middle of the first clear break in the ordered tumour medians): %.3f\n", SCS_CUT)); print(ss, row.names = FALSE, digits = 3)
EPI <- ss$sample[ss$call == "epithelioid"]
md$group <- ifelse(md$celltype_lv1 == "Mesothelium", "Normal mesothelium", ifelse(md$sample %in% EPI, "Epithelioid tumour cells", "Higher-scS tumour cells"))
md$group_cell <- ifelse(md$celltype_lv1 == "Mesothelium", "Normal mesothelium", ifelse(md$scS < SCS_CUT, "Epithelioid tumour cells", "Higher-scS tumour cells"))
GL <- c("Normal mesothelium", "Epithelioid tumour cells", "Higher-scS tumour cells"); cat("\ncells per group (tumour-level call):\n"); print(table(md$group)[GL]); cat("cells per group (cell-level call):\n"); print(table(md$group_cell)[GL])
cat("scS of the mesothelium cells: median", round(median(md$scS[md$group == "Normal mesothelium"]), 3), "(they express the mesenchymal genes of the score, so the score is only used to split MALIGNANT cells)\n")

## ---- keratins ---------------------------------------------------------------------------------
X <- GetAssayData(srt, layer = "data"); kr <- sort(grep("^KRT[0-9]", rownames(X), value = TRUE))
pct <- sapply(GL, function(g) Matrix::rowMeans(X[kr, md$group == g, drop = FALSE] > 0)); kr <- kr[apply(pct, 1, max) >= MIN_PCT]
cat("\nkeratin genes detected in >=", 100 * MIN_PCT, "% of the cells of a group:", length(kr), "\n")
stat <- function(gcol) do.call(rbind, lapply(kr, function(g) { a <- X[g, md[[gcol]] == "Epithelioid tumour cells"]; b <- X[g, md[[gcol]] == "Normal mesothelium"]
  data.frame(gene = g, pct_epithelioid = 100 * mean(a > 0), pct_mesothelium = 100 * mean(b > 0), mean_epithelioid = mean(a), mean_mesothelium = mean(b),
             log2FC = log2((mean(expm1(a)) + 0.01) / (mean(expm1(b)) + 0.01)), p_wilcox = suppressWarnings(wilcox.test(a, b)$p.value)) }))
ST <- stat("group"); ST$fdr <- p.adjust(ST$p_wilcox, "BH"); S2 <- stat("group_cell"); ST$log2FC_cell_level_call <- S2$log2FC; ST$fdr_cell_level_call <- p.adjust(S2$p_wilcox, "BH")
## pseudobulk per patient x compartment: mean of normalised counts, patients with >= MIN_CELLS cells
md$unit <- paste(md$sample, ifelse(md$celltype_lv1 == "Mesothelium", "Mesothelium", "Malignant"), sep = "|"); tu <- table(md$unit); units <- names(tu)[tu >= MIN_CELLS]
PB <- sapply(units, function(u) Matrix::rowMeans(expm1(X[kr, md$unit == u, drop = FALSE]))); PBl <- log2(PB + 0.01)
ui <- data.frame(unit = units, sample = sub("\\|.*", "", units), compartment = sub(".*\\|", "", units), n_cells = as.integer(tu[units])); ui$call <- ifelse(ui$compartment == "Mesothelium", "Normal mesothelium", ss$call[match(ui$sample, ss$sample)]); ui$scS <- ss$median_scS[match(ui$sample, ss$sample)]
write.csv(cbind(ui, t(round(PBl, 3))), "R2_Q10_keratin_pseudobulk_per_sample.csv", row.names = FALSE)
ue <- ui$unit[ui$call == "epithelioid"]; um <- ui$unit[ui$call == "Normal mesothelium"]
ST$pseudobulk_log2FC <- rowMeans(PBl[, ue, drop = FALSE]) - rowMeans(PBl[, um, drop = FALSE]); ST$n_epi_tumours_above_all_meso <- rowSums(PBl[, ue, drop = FALSE] > apply(PBl[, um, drop = FALSE], 1, max)); ST$n_epi_tumours_below_all_meso <- rowSums(PBl[, ue, drop = FALSE] < apply(PBl[, um, drop = FALSE], 1, min))
ST$p_pseudobulk <- sapply(kr, function(g) suppressWarnings(wilcox.test(PBl[g, ue], PBl[g, um])$p.value)); ST <- ST[order(-ST$log2FC), ]
write.csv(ST, "R2_Q10_keratin_stats.csv", row.names = FALSE)
cat(sprintf("\n=== keratins: epithelioid tumour cells (%d cells, %d tumours) vs normal mesothelium (%d cells, %d patients) ===\n", sum(md$group == GL[2]), length(ue), sum(md$group == GL[1]), length(um))); print(format(ST, digits = 3), row.names = FALSE)
cat(sprintf("\ncell level: %d of %d keratins higher in epithelioid cells (FDR < 0.05), %d lower; sample level: smallest possible two-sided Wilcoxon p with %d vs %d units = %.3f\n",
            sum(ST$fdr < 0.05 & ST$log2FC > 0), nrow(ST), sum(ST$fdr < 0.05 & ST$log2FC < 0), length(ue), length(um), 2 / choose(length(ue) + length(um), length(um))))

cat("\nDONE\n")
