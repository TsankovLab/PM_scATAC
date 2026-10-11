# Spec 40: repeat of spec39 (screen-to-bulk agreement across all 30 cNMF modules, k=30), but with cell cycle regressed out of the
# BULK expression matrix first (same per-line reg() convention as spec36/37/38's bulk half), before scoring modules with
# Seurat::AddModuleScore. The screen side is unchanged (cNMF usages are a precomputed decomposition, not something to cc-regress
# post hoc). 6/76 bulk samples lack a precomputed cc score (bulk_cellcycle_per_sample_RUN1.csv) and are dropped for this variant.
# needs: same as spec39, + bulk_cellcycle_per_sample_RUN1.csv
# outputs: spec40_all_modules_agreement_ccregressed.csv, spec40_top_KO_confusion_ccregressed.csv,
#          fig2X_ccregressed_all_modules_correlation.{png,pdf}, fig2X_ccregressed_top_bottom_KO_confusion.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat))
CNMF_D <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CRISPR_cropseq_analysis/_cellranger_raw_Filter_400_800_25/no_harmony/cNMF/cNMF_crispr_20_30_vf5000/cnmf"
K <- "30"; MODULES <- as.character(1:30)
GROUPS <- c("NTC", ORDER)   # the 9 groups shared by screen and bulk

## ---------------------------------------------------------------- screen: per-module top50 genes + per-KO usage (unchanged from spec39, not cc-regressed)
gss <- read.delim(file.path(CNMF_D, sprintf("cnmf.gene_spectra_score.k_%s.dt_0_3.txt", K)), row.names = 1, check.names = FALSE)
usages <- read.delim(file.path(CNMF_D, sprintf("cnmf.usages.k_%s.dt_0_3.consensus.txt", K)), row.names = 1, check.names = FALSE)
srt <- load_screen(); grp_scr <- setNames(as.character(srt$merged_call), colnames(srt))
common <- intersect(rownames(usages), names(grp_scr)); grp_scr <- grp_scr[common]; usages <- usages[common, ]
cat(sprintf("screen: %d cells (overlap cNMF usages / filtered screen)\n", length(common)))

top_genes_list <- lapply(MODULES, function(m) { x <- sort(setNames(as.numeric(gss[m, ]), colnames(gss)), decreasing = TRUE); names(x)[1:50] })
names(top_genes_list) <- paste0("m", MODULES)

screen_usage <- sapply(MODULES, function(m) sapply(GROUPS, function(grp) mean(usages[grp_scr == grp, m])))
rownames(screen_usage) <- GROUPS

## ---------------------------------------------------------------- bulk: regress cell cycle out per line first, then score all 30 modules
fab <- load_full_audited_bulk(); meta_all <- fab$meta; lg_all <- fab$lg
ccb <- read.csv(file.path(OUT, "bulk_cellcycle_per_sample_RUN1.csv"))
has_cc <- rownames(meta_all)[rownames(meta_all) %in% ccb$sample]
cat(sprintf("dropping %d/%d bulk samples with no precomputed cc score: %s\n", nrow(meta_all) - length(has_cc), nrow(meta_all),
            paste(meta_all[setdiff(rownames(meta_all), has_cc), "Sample_ID"], collapse = ", ")))
meta <- meta_all[has_cc, ]
lg_cc <- lg_all[, rownames(meta), drop = FALSE]
for (cl in LINES) {
  s <- rownames(meta)[meta$Cell_Line == cl]; cc <- ccb$cc[match(s, ccb$sample)]; ccc <- cc - mean(cc)
  M <- lg_cc[, s, drop = FALSE]; b <- as.numeric((M - rowMeans(M)) %*% ccc) / sum(ccc^2)
  lg_cc[, s] <- M - outer(b, ccc)
}

options(Seurat.object.assay.version = "v3")   # AddModuleScore errors ("subscript out of bounds") on the default Assay5 class here; classic Assay works
srt_bulk <- CreateSeuratObject(counts = as.matrix(fab$cs[, rownames(meta), drop = FALSE]))
srt_bulk <- SetAssayData(srt_bulk, slot = "data", new.data = as.matrix(lg_cc))
srt_bulk <- AddModuleScore(srt_bulk, features = top_genes_list, name = "m", seed = 1)
bulk_score_cells <- as.matrix(srt_bulk@meta.data[, paste0("m", seq_along(MODULES))]); colnames(bulk_score_cells) <- MODULES
bulk_score_cells <- bulk_score_cells[rownames(meta), ]

bulk_score <- sapply(MODULES, function(m) sapply(GROUPS, function(grp) mean(bulk_score_cells[meta$Gene_Control == grp, m])))
rownames(bulk_score) <- GROUPS

## ---------------------------------------------------------------- per-module agreement: correlation + top/bottom match
agree <- do.call(rbind, lapply(MODULES, function(m) {
  s <- screen_usage[, m]; b <- bulk_score[, m]
  pr <- cor(s, b, method = "pearson"); sr <- suppressWarnings(cor(s, b, method = "spearman"))
  pt <- suppressWarnings(cor.test(s, b, method = "pearson")$p.value)
  data.frame(module = m, screen_top = names(which.max(s)), screen_bottom = names(which.min(s)),
             bulk_top = names(which.max(b)), bulk_bottom = names(which.min(b)),
             top_match = names(which.max(s)) == names(which.max(b)), bottom_match = names(which.min(s)) == names(which.min(b)),
             pearson_r = pr, spearman_rho = sr, pearson_p = pt)
}))
agree <- agree[order(-agree$pearson_r), ]
write.csv(agree, file.path(OUT, "spec40_all_modules_agreement_ccregressed.csv"), row.names = FALSE)

options(width = 200)
cat("\n=== per-module screen-vs-bulk(cc-regressed) agreement (all 30 cNMF modules, k=30), sorted by Pearson r ===\n")
print(format(agree, digits = 3), row.names = FALSE)

cat(sprintf("\n=== summary across all %d modules (bulk cc-regressed) ===\n", nrow(agree)))
cat(sprintf("top-KO match (screen argmax == bulk argmax):    %d/%d (%.0f%%) [chance level with %d groups = %.0f%%]\n",
            sum(agree$top_match), nrow(agree), 100 * mean(agree$top_match), length(GROUPS), 100 / length(GROUPS)))
cat(sprintf("bottom-KO match (screen argmin == bulk argmin): %d/%d (%.0f%%) [chance level = %.0f%%]\n",
            sum(agree$bottom_match), nrow(agree), 100 * mean(agree$bottom_match), 100 / length(GROUPS)))
cat(sprintf("mean Pearson r = %.3f | median = %.3f | modules with r > 0: %d/%d | r > 0.3: %d/%d | significant (p<0.05): %d/%d\n",
            mean(agree$pearson_r), median(agree$pearson_r), sum(agree$pearson_r > 0), nrow(agree),
            sum(agree$pearson_r > 0.3), nrow(agree), sum(agree$pearson_p < 0.05), nrow(agree)))
cat(sprintf("one-sample t-test (is mean r > 0?): t=%.2f, p=%.4f\n", t.test(agree$pearson_r)$statistic, t.test(agree$pearson_r)$p.value))

## ---------------------------------------------------------------- confusion matrix
conf_top <- table(factor(agree$screen_top, levels = GROUPS), factor(agree$bulk_top, levels = GROUPS))
conf_bottom <- table(factor(agree$screen_bottom, levels = GROUPS), factor(agree$bulk_bottom, levels = GROUPS))
write.csv(as.data.frame.matrix(conf_top), file.path(OUT, "spec40_top_KO_confusion_ccregressed.csv"))
cat("\n=== confusion: screen's top-KO (rows) vs bulk's top-KO (cols), tallied over 30 modules ===\n"); print(conf_top)
cat("\n=== how often is each KO the SCREEN top module owner? ===\n"); print(table(agree$screen_top))
cat("\n=== how often is each KO the BULK top module owner? ===\n"); print(table(agree$bulk_top))

## ---------------------------------------------------------------- fig: per-module Pearson r, sorted, colour = significant
save_fig("fig2X_ccregressed_all_modules_correlation", 11, 6, function() {
  ord <- agree[order(-agree$pearson_r), ]
  cols <- ifelse(ord$pearson_p < 0.05 & ord$pearson_r > 0, "#1baf7a", ifelse(ord$pearson_p < 0.05 & ord$pearson_r < 0, "#eb6834", adjustcolor(INK2, 0.4)))
  par(mar = c(5, 4.5, 3.6, 1), family = "sans")
  bp <- barplot(ord$pearson_r, names.arg = paste0("m", ord$module), col = cols, border = INK2, las = 2, cex.names = 0.65,
                ylab = "Pearson r (screen vs bulk, cc-regressed)", ylim = c(-1, 1))
  abline(h = 0, col = INK)
  legend("topright", legend = c("p<0.05, r>0", "p<0.05, r<0", "n.s."), fill = c("#1baf7a", "#eb6834", adjustcolor(INK2, 0.4)), bty = "n", cex = 0.7)
  mtext("Screen-to-bulk agreement across all 30 cNMF modules (k=30), bulk cc-regressed", side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext(sprintf("per module: Pearson r between screen per-KO usage and bulk per-KO AddModuleScore (NTC + 8 shared KOs); top-KO match %d/30, bottom-KO match %d/30",
                sum(agree$top_match), sum(agree$bottom_match)), side = 3, adj = 0, line = 0.9, cex = 0.58, col = INK2)
})

## ---------------------------------------------------------------- fig: confusion heatmap
save_fig("fig2X_ccregressed_top_bottom_KO_confusion", 12, 6, function() {
  layout(matrix(1:2, nrow = 1)); par(oma = c(0, 0, 3, 0), family = "sans")
  draw_conf <- function(cm, title) {
    n <- length(GROUPS); par(mar = c(6.5, 6.5, 2.5, 1))
    image(seq_len(n), seq_len(n), t(cm)[, n:1], col = colorRampPalette(c("white", "#4a3aa7"))(101), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(n)) for (j in seq_len(n)) if (cm[i, j] > 0) text(j, n - i + 1, cm[i, j], cex = 0.9, col = ifelse(cm[i, j] > max(cm) / 2, "white", INK))
    axis(1, at = seq_len(n), labels = GROUPS, las = 2, tick = FALSE, cex.axis = 0.75, col.axis = INK2)
    axis(2, at = seq_len(n), labels = rev(GROUPS), las = 2, tick = FALSE, cex.axis = 0.75, col.axis = INK2)
    box(col = MUTED, lwd = 0.6); mtext(title, side = 3, line = 0.6, font = 2, cex = 0.85)
    mtext("bulk", side = 1, line = 5, cex = 0.7); mtext("screen", side = 2, line = 5, cex = 0.7)
  }
  draw_conf(conf_top, "top-KO owner: screen (rows) vs bulk cc-regressed (cols)")
  draw_conf(conf_bottom, "bottom-KO owner: screen (rows) vs bulk cc-regressed (cols)")
  mtext("Which KO group 'owns' each module's screen/bulk extreme, tallied across 30 modules (bulk cc-regressed)", outer = TRUE, side = 3, adj = 0, line = 1.2, font = 2, cex = 1.0)
})
