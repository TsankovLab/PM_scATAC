# Spec 46: use NCI-H2052 -- the ONE cell line present in BOTH the CROP-seq screen and the bulk RNA-seq -- to model and remove
# platform/technology-specific noise before comparing screen to the other 3 bulk lines. Logic: since NCI-H2052 is the exact same
# genetic background in both experiments, any systematic screen-vs-H2052-bulk relationship that deviates from a 1:1 line must be a
# technology effect (different platform dynamic range/compression/sensitivity per Cm module), not a biological cross-line effect --
# because there IS no cross-line biology in this one comparison. For each of the 20 Cm modules separately, fit a linear "technology
# transfer function" (H2052 bulk KO-NTC diff ~ screen KO-NTC diff, across the 8 shared KOs) and use it to rescale the screen's own
# diff profile onto the bulk scale, per module, BEFORE correlating with the other 3 (true cross-line) bulk lines. A plain single
# global rescaling would not change Pearson r at all (correlation is invariant to one common affine transform) -- fitting PER MODULE
# and applying different module-specific corrections is what can actually change the result, since it reweights each module's
# contribution to the cross-module correlation based on how reliably it transfers between platforms in the one line where the
# biology is guaranteed identical.
# needs: git_repo/files/scRNA_PM_modules.rds, spec_common.R (load_screen()), repro_common.R (load_bulk())
# outputs: spec46_H2052_technology_transfer_fits.csv, spec46_before_after_correction.csv,
#          fig3E_H2052_technology_transfer_fits.{png,pdf}, fig3E_before_after_correction.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat))

cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds")
top50 <- lapply(cm_genes, head, 50); CM <- names(top50)
options(Seurat.object.assay.version = "v3")

## ---------------------------------------------------------------- screen: KO-NTC diff (median), same as spec45
srt <- load_screen(); srt <- AddModuleScore(srt, features = top50, name = "Cm_", seed = 1)
screen_cell_scores <- as.matrix(srt@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(screen_cell_scores) <- CM
grp_scr <- as.character(srt$merged_call)
screen_median <- t(sapply(unique(grp_scr), function(g) apply(screen_cell_scores[grp_scr == g, , drop = FALSE], 2, median)))
rownames(screen_median) <- unique(grp_scr)
common_ko <- intersect(rownames(screen_median), ORDER)
screen_diff <- sweep(screen_median[common_ko, CM, drop = FALSE], 2, screen_median["NTC", CM], "-")   # KO - NTC

## ---------------------------------------------------------------- bulk: same top-50 gene sets, all 4 lines
b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1)
srt_bulk <- CreateSeuratObject(counts = as.matrix(b$counts)); srt_bulk <- SetAssayData(srt_bulk, slot = "data", new.data = as.matrix(lg))
srt_bulk <- AddModuleScore(srt_bulk, features = top50, name = "Cm_", seed = 1)
bulk_cell_scores <- as.matrix(srt_bulk@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(bulk_cell_scores) <- CM
bulk_cell_scores <- bulk_cell_scores[rownames(meta), ]

bulk_diff_by_line <- function(cl) {
  mm <- meta[meta$Cell_Line == cl, ]; ntc_s <- rownames(mm)[mm$Gene_Control == "NTC"]
  t(sapply(common_ko, function(k) { ko_s <- rownames(mm)[mm$Gene_Control == k]
    if (!length(ko_s) || !length(ntc_s)) return(setNames(rep(NA, length(CM)), CM))
    apply(bulk_cell_scores[ko_s, CM, drop = FALSE], 2, median) - apply(bulk_cell_scores[ntc_s, CM, drop = FALSE], 2, median) }))
}
h2052_diff <- bulk_diff_by_line("NCI-H2052"); rownames(h2052_diff) <- common_ko   # calibration set: n=1 NTC sample, used only for fitting, not evaluation
OTHER_LINES <- setdiff(LINES, "NCI-H2052")
other_diff <- lapply(OTHER_LINES, function(cl) { d <- bulk_diff_by_line(cl); rownames(d) <- common_ko; d })
names(other_diff) <- OTHER_LINES

## ---------------------------------------------------------------- fit ONE technology-transfer regression per Cm module: H2052_bulk ~ screen (n=8 KOs)
fits <- lapply(CM, function(m) lm(y ~ x, data = data.frame(y = h2052_diff[, m], x = screen_diff[, m])))
names(fits) <- CM
fit_summary <- do.call(rbind, lapply(CM, function(m) {
  f <- fits[[m]]; s <- summary(f)
  data.frame(Cm = m, intercept = coef(f)[1], slope = coef(f)[2], r2 = s$r.squared, p_slope = coef(s)[2, 4])
}))
write.csv(fit_summary, file.path(OUT, "spec46_H2052_technology_transfer_fits.csv"), row.names = FALSE)
options(width = 200)
cat("=== per-module technology-transfer fit (H2052 bulk KO-NTC diff ~ screen KO-NTC diff, n=8 KOs) ===\n")
print(format(fit_summary[order(-fit_summary$r2), ], digits = 3), row.names = FALSE)
cat(sprintf("\nmedian R2 = %.3f (how well screen predicts THIS SAME line's bulk, per module -- the technology ceiling)\n", median(fit_summary$r2)))

## ---------------------------------------------------------------- apply the per-module fit to rescale the screen profile onto the bulk scale, then re-test vs the OTHER 3 lines
screen_corrected <- screen_diff
for (m in CM) screen_corrected[, m] <- predict(fits[[m]], newdata = data.frame(x = screen_diff[, m]))

before_after <- do.call(rbind, lapply(OTHER_LINES, function(cl) do.call(rbind, lapply(common_ko, function(k) {
  bl <- other_diff[[cl]][k, ]; if (all(is.na(bl))) return(NULL)
  r_before <- cor(screen_diff[k, ], bl, use = "complete.obs")
  r_after <- cor(screen_corrected[k, ], bl, use = "complete.obs")
  data.frame(line = cl, KO = k, r_before = r_before, r_after = r_after, delta = r_after - r_before)
}))))
write.csv(before_after, file.path(OUT, "spec46_before_after_correction.csv"), row.names = FALSE)
cat("\n=== screen-vs-OTHER-line correlation, before vs after H2052-calibrated per-module technology correction ===\n")
print(format(before_after, digits = 3), row.names = FALSE)
cat(sprintf("\nmedian r_before = %.3f | median r_after = %.3f | median delta = %+.3f | improved in %d/%d line-KO pairs\n",
            median(before_after$r_before), median(before_after$r_after), median(before_after$delta), sum(before_after$delta > 0), nrow(before_after)))

## ---------------------------------------------------------------- fig: per-module technology-transfer fit quality (R2), sorted
save_fig("fig3E_H2052_technology_transfer_fits", 9, 5.5, function() {
  ord <- fit_summary[order(-fit_summary$r2), ]
  cols <- ifelse(ord$p_slope < 0.05 & ord$slope > 0, "#1baf7a", adjustcolor(INK2, 0.4))
  par(mar = c(5, 4.8, 3.6, 1), family = "sans")
  bp <- barplot(ord$r2, names.arg = ord$Cm, col = cols, border = INK2, las = 2, cex.names = 0.7, ylab = "R2 (H2052 bulk ~ screen, n=8 KOs)", ylim = c(0, 1))
  mtext("Per-module technology-transfer fit quality (screen vs H2052's OWN bulk)", side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext("green = significant positive slope (p<0.05); this R2 is the best-case ceiling for screen-to-bulk agreement (same line, cross-platform only)", side = 3, adj = 0, line = 0.9, cex = 0.6, col = INK2)
})

## ---------------------------------------------------------------- fig: before vs after correction, per TF boxplot pair
save_fig("fig3E_before_after_correction", 8.5, 5.5, function() {
  ko_ord <- names(sort(sapply(split(before_after$r_before, before_after$KO), median), decreasing = TRUE))
  n_ko <- length(ko_ord)
  par(mar = c(5, 4.8, 3.6, 1), family = "sans")
  plot(NA, xlim = c(0.5, n_ko + 0.5), ylim = c(-1, 1), xaxt = "n", xlab = "", ylab = "Pearson r (screen vs other-line bulk, KO-NTC diff)")
  axis(1, at = seq_len(n_ko), labels = FALSE)
  text(seq_len(n_ko), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), ko_ord, srt = 45, adj = 1, xpd = NA)
  for (i in seq_len(n_ko)) {
    d0 <- before_after$r_before[before_after$KO == ko_ord[i]]
    boxplot(d0, at = i - 0.18, add = TRUE, boxwex = 0.3, col = "grey85", border = "grey40", outline = FALSE, axes = FALSE)
    d1 <- before_after$r_after[before_after$KO == ko_ord[i]]
    boxplot(d1, at = i + 0.18, add = TRUE, boxwex = 0.3, col = adjustcolor(KO_COL[ko_ord[i]], 0.5), border = INK2, outline = FALSE, axes = FALSE)
  }
  abline(h = 0, col = INK, lty = 2)
  legend("bottomleft", legend = c("before correction", "after H2052 correction"), fill = c("grey85", adjustcolor("grey40", 0.5)), border = c("grey40", INK2), bty = "n", cex = 0.65)
  mtext("Effect of H2052-calibrated technology correction on screen-to-bulk agreement", side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext("grey = raw screen-vs-other-line correlation; coloured = after per-module rescaling using NCI-H2052 as the same-line technology calibrator", side = 3, adj = 0, line = 0.9, cex = 0.6, col = INK2)
})
