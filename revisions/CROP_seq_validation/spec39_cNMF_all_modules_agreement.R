# Spec 39: for ALL 30 cNMF modules (k=30 solution, CROP-seq screen's own cNMF sweep), check general screen-to-bulk agreement --
# not just "does the single highest/lowest KO match" (spec38, done module-by-module) but a systematic sweep: for each module,
# take its top-50 defining genes, get (a) mean per-cell cNMF usage per KO group in the screen and (b) mean Seurat::AddModuleScore
# per KO group in bulk RNA-seq (matched to the 9 groups both datasets share: NTC + the 8 CROP-seq/bulk KOs), then compute the
# Pearson/Spearman correlation between the screen and bulk per-KO profiles, and whether the top/bottom KO matches between the two.
# needs: CRISPR_cropseq_analysis/.../cNMF_crispr_20_30_vf5000/cnmf/cnmf.usages.k_30.dt_0_3.consensus.txt,
#        cnmf.gene_spectra_score.k_30.dt_0_3.txt, spec_common.R (load_screen()), repro_common.R (load_full_audited_bulk())
# outputs: spec39_all_modules_agreement.csv, spec39_top_KO_confusion.csv,
#          fig2X_all_modules_correlation.{png,pdf}, fig2X_top_bottom_KO_confusion.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat))
CNMF_D <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CRISPR_cropseq_analysis/_cellranger_raw_Filter_400_800_25/no_harmony/cNMF/cNMF_crispr_20_30_vf5000/cnmf"
K <- "30"; MODULES <- as.character(1:30)
GROUPS <- c("NTC", ORDER)   # the 9 groups shared by screen and bulk

## ---------------------------------------------------------------- screen: per-module top50 genes + per-KO usage (restricted to shared groups)
gss <- read.delim(file.path(CNMF_D, sprintf("cnmf.gene_spectra_score.k_%s.dt_0_3.txt", K)), row.names = 1, check.names = FALSE)
usages <- read.delim(file.path(CNMF_D, sprintf("cnmf.usages.k_%s.dt_0_3.consensus.txt", K)), row.names = 1, check.names = FALSE)
srt <- load_screen(); grp_scr <- setNames(as.character(srt$merged_call), colnames(srt))
common <- intersect(rownames(usages), names(grp_scr)); grp_scr <- grp_scr[common]; usages <- usages[common, ]
cat(sprintf("screen: %d cells (overlap cNMF usages / filtered screen)\n", length(common)))

top_genes_list <- lapply(MODULES, function(m) { x <- sort(setNames(as.numeric(gss[m, ]), colnames(gss)), decreasing = TRUE); names(x)[1:50] })
names(top_genes_list) <- paste0("m", MODULES)

screen_usage <- sapply(MODULES, function(m) sapply(GROUPS, function(grp) mean(usages[grp_scr == grp, m])))
rownames(screen_usage) <- GROUPS   # rows = groups, cols = modules

## ---------------------------------------------------------------- bulk: score all 30 modules at once via AddModuleScore
fab <- load_full_audited_bulk(); meta <- fab$meta; lg <- fab$lg
options(Seurat.object.assay.version = "v3")   # AddModuleScore errors ("subscript out of bounds") on the default Assay5 class here; classic Assay works
srt_bulk <- CreateSeuratObject(counts = as.matrix(fab$cs)); srt_bulk <- SetAssayData(srt_bulk, slot = "data", new.data = as.matrix(lg))
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
write.csv(agree, file.path(OUT, "spec39_all_modules_agreement.csv"), row.names = FALSE)

options(width = 200)
cat("\n=== per-module screen-vs-bulk agreement (all 30 cNMF modules, k=30), sorted by Pearson r ===\n")
print(format(agree, digits = 3), row.names = FALSE)

cat(sprintf("\n=== summary across all %d modules ===\n", nrow(agree)))
cat(sprintf("top-KO match (screen argmax == bulk argmax):    %d/%d (%.0f%%) [chance level with %d groups = %.0f%%]\n",
            sum(agree$top_match), nrow(agree), 100 * mean(agree$top_match), length(GROUPS), 100 / length(GROUPS)))
cat(sprintf("bottom-KO match (screen argmin == bulk argmin): %d/%d (%.0f%%) [chance level = %.0f%%]\n",
            sum(agree$bottom_match), nrow(agree), 100 * mean(agree$bottom_match), 100 / length(GROUPS)))
cat(sprintf("mean Pearson r = %.3f | median = %.3f | modules with r > 0: %d/%d | r > 0.3: %d/%d | significant (p<0.05): %d/%d\n",
            mean(agree$pearson_r), median(agree$pearson_r), sum(agree$pearson_r > 0), nrow(agree),
            sum(agree$pearson_r > 0.3), nrow(agree), sum(agree$pearson_p < 0.05), nrow(agree)))
cat(sprintf("one-sample t-test (is mean r > 0?): t=%.2f, p=%.4f\n", t.test(agree$pearson_r)$statistic, t.test(agree$pearson_r)$p.value))

## ---------------------------------------------------------------- confusion matrix: which KO is "top" in screen vs "top" in bulk, tallied across the 30 modules
conf_top <- table(factor(agree$screen_top, levels = GROUPS), factor(agree$bulk_top, levels = GROUPS))
conf_bottom <- table(factor(agree$screen_bottom, levels = GROUPS), factor(agree$bulk_bottom, levels = GROUPS))
write.csv(as.data.frame.matrix(conf_top), file.path(OUT, "spec39_top_KO_confusion.csv"))
cat("\n=== confusion: screen's top-KO (rows) vs bulk's top-KO (cols), tallied over 30 modules ===\n"); print(conf_top)
cat("\n=== how often is each KO the SCREEN top module owner? ===\n"); print(table(agree$screen_top))
cat("\n=== how often is each KO the BULK top module owner? ===\n"); print(table(agree$bulk_top))

## ---------------------------------------------------------------- fig: per-module Pearson r, sorted, colour = significant
save_fig("fig2X_all_modules_correlation", 11, 6, function() {
  ord <- agree[order(-agree$pearson_r), ]
  cols <- ifelse(ord$pearson_p < 0.05 & ord$pearson_r > 0, "#1baf7a", ifelse(ord$pearson_p < 0.05 & ord$pearson_r < 0, "#eb6834", adjustcolor(INK2, 0.4)))
  par(mar = c(5, 4.5, 3.6, 1), family = "sans")
  bp <- barplot(ord$pearson_r, names.arg = paste0("m", ord$module), col = cols, border = INK2, las = 2, cex.names = 0.65,
                ylab = "Pearson r (screen usage vs bulk AddModuleScore, across 9 groups)", ylim = c(-1, 1))
  abline(h = 0, col = INK)
  legend("topright", legend = c("p<0.05, r>0", "p<0.05, r<0", "n.s."), fill = c("#1baf7a", "#eb6834", adjustcolor(INK2, 0.4)), bty = "n", cex = 0.7)
  mtext("Screen-to-bulk agreement across all 30 cNMF modules (k=30)", side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext(sprintf("per module: Pearson r between screen per-KO usage and bulk per-KO AddModuleScore (NTC + 8 shared KOs); top-KO match %d/30, bottom-KO match %d/30",
                sum(agree$top_match), sum(agree$bottom_match)), side = 3, adj = 0, line = 0.9, cex = 0.58, col = INK2)
})

## ---------------------------------------------------------------- fig: confusion heatmap (screen top-KO vs bulk top-KO)
save_fig("fig2X_top_bottom_KO_confusion", 12, 6, function() {
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
  draw_conf(conf_top, "top-KO owner: screen (rows) vs bulk (cols)")
  draw_conf(conf_bottom, "bottom-KO owner: screen (rows) vs bulk (cols)")
  mtext("Which KO group 'owns' each module's screen/bulk extreme, tallied across 30 modules", outer = TRUE, side = 3, adj = 0, line = 1.2, font = 2, cex = 1.0)
})

## ---------------------------------------------------------------- alternate view: per-KO correlation ACROSS the 30 modules, using RAW group
## values (no NTC differencing at all -- literally screen_usage[k,] vs bulk_score[k,], the KO group's own raw profile in each dataset)
raw_ko <- do.call(rbind, lapply(c("NTC", ORDER), function(k) {
  s <- screen_usage[k, ]; b <- bulk_score[k, ]
  pr <- cor.test(s, b, method = "pearson"); sr <- suppressWarnings(cor.test(s, b, method = "spearman", exact = FALSE))
  data.frame(KO = k, pearson_r = pr$estimate, pearson_p = pr$p.value, spearman_rho = sr$estimate, spearman_p = sr$p.value)
}))
raw_ko <- raw_ko[order(-raw_ko$pearson_r), ]
write.csv(raw_ko, file.path(OUT, "spec39_raw_KO_correlation_no_NTC_diff.csv"), row.names = FALSE)
cat("\n=== per-KO correlation, RAW group profiles across 30 modules (no NTC differencing) -- NTC itself included as a baseline ===\n")
print(format(raw_ko, digits = 3), row.names = FALSE)
ko_only <- raw_ko[raw_ko$KO != "NTC", ]
cat(sprintf("\nKOs only: mean Pearson r = %.3f | median = %.3f | KOs with r > 0: %d/%d | significant (p<0.05): %d/%d\n",
            mean(ko_only$pearson_r), median(ko_only$pearson_r), sum(ko_only$pearson_r > 0), nrow(ko_only), sum(ko_only$pearson_p < 0.05), nrow(ko_only)))
cat(sprintf("CAUTION: NTC alone (no perturbation at all) correlates at r=%.3f (p=%.2e) -- essentially as high as most KOs. This raw,\n", raw_ko$pearson_r[raw_ko$KO == "NTC"], raw_ko$pearson_p[raw_ko$KO == "NTC"]))
cat("non-differenced correlation is dominated by which modules are generically high/low-expressed in this cell-type context in BOTH\n")
cat("datasets, not by KO-specific perturbation signal -- it is NOT a valid measure of screen-to-bulk KO validation by itself.\n")

save_fig("fig2X_raw_KO_correlation_no_NTC_diff", 8.5, 5.5, function() {
  ord <- raw_ko[order(-raw_ko$pearson_r), ]
  cols <- ifelse(ord$KO == "NTC", "grey40", ifelse(ord$pearson_p < 0.05, KO_COL[ord$KO], adjustcolor(KO_COL[ord$KO], 0.4)))
  par(mar = c(6, 4.8, 4.2, 1), family = "sans")
  bp <- barplot(ord$pearson_r, names.arg = ord$KO, col = cols, border = INK2, las = 2, cex.names = 0.85,
                ylab = "Pearson r (raw screen vs raw bulk, 30 modules)", ylim = c(-1, 1))
  abline(h = 0, col = INK)
  abline(h = raw_ko$pearson_r[raw_ko$KO == "NTC"], col = "grey40", lty = 2)
  mtext("Raw screen KO profile vs raw bulk KO profile (no NTC differencing)", side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext("dashed line = NTC's own baseline correlation (grey bar); most KOs don't clearly exceed it -- this raw comparison mostly reflects", side = 3, adj = 0, line = 0.9, cex = 0.58, col = INK2)
  mtext("generic module-expression-level similarity between the two datasets, NOT KO-specific perturbation signal", side = 3, adj = 0, line = 0.0, cex = 0.58, col = INK2)
})
