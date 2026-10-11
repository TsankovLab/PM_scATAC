# Spec 63: boxplot of intra-KO (same-KO replicate pairs) vs cross-KO (different-KO pairs) sample-sample Spearman correlation, bulk NCI-H2052 only, using the
# spec62 sample correlation matrices (20 Cm modules, top 1000/2000/3000 variable genes). One panel per feature space; NTC excluded from "intra" (only 1 NTC sample
# in the curated sheet, so no same-KO pair exists for it) but shown as a labelled reference point (NTC vs each KO's replicate mean).
# needs: spec62_H2052_bulk_KO_correlation_heatmap.R (spec62_H2052_sample_correlation.csv)
# outputs: spec63_intra_vs_cross_summary.csv, fig3Y_H2052_intra_vs_cross_replicate_boxplot.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
long <- read.csv(file.path(OUT, "spec62_H2052_sample_correlation_raw.csv"), stringsAsFactors = FALSE)
ko_of <- function(s) sub("^H2052_", "", sub("_R[0-9]+$", "", s))
long$ko_a <- ko_of(long$sample_a); long$ko_b <- ko_of(long$sample_b)
long$pair_type <- ifelse(long$ko_a == long$ko_b, "intra-KO (replicates)", "cross-KO"); long$pair_type[long$ko_a == "NTC" | long$ko_b == "NTC"] <- "vs NTC"
spaces <- unique(long$space)

summ <- do.call(rbind, lapply(spaces, function(sp) { d <- long[long$space == sp, ]
  do.call(rbind, lapply(split(d, d$pair_type), function(x) data.frame(space = sp, pair_type = x$pair_type[1], n = nrow(x), mean_rho = mean(x$rho), median_rho = median(x$rho), sd_rho = sd(x$rho)))) }))
write.csv(summ, file.path(OUT, "spec63_intra_vs_cross_summary.csv"), row.names = FALSE)
options(width = 200); cat("\n=== intra-KO (replicate pairs) vs cross-KO vs vs-NTC, per feature space (NCI-H2052 bulk, Spearman) ===\n"); print(format(summ, digits = 3), row.names = FALSE)
for (sp in spaces) { d <- long[long$space == sp, ]; iw <- d$rho[d$pair_type == "intra-KO (replicates)"]; cr <- d$rho[d$pair_type == "cross-KO"]
  wt <- wilcox.test(iw, cr); cat(sprintf("%-24s intra-KO median %.3f vs cross-KO median %.3f (Wilcoxon p=%.3g)\n", sp, median(iw), median(cr), wt$p.value)) }

save_fig("fig3Y_H2052_intra_vs_cross_replicate_boxplot", 10.5, 5.6, function() {
  layout(matrix(1:4, nrow = 1)); par(oma = c(0, 0, 3.4, 0), family = "sans")
  ord <- c("intra-KO (replicates)", "cross-KO", "vs NTC"); cols <- c("#1baf7a", "grey75", "#eda100")
  for (sp in spaces) { d <- long[long$space == sp, ]; d$pair_type <- factor(d$pair_type, levels = ord); par(mar = c(7, 4.4, 2.6, 0.8))
    boxplot(rho ~ pair_type, d, col = cols, border = INK2, outline = FALSE, las = 2, xlab = "", ylab = "Spearman rho (sample pair)", names = c("intra-KO\n(replicates)", "cross-KO", "vs NTC"))
    set.seed(1); stripchart(rho ~ pair_type, d, vertical = TRUE, add = TRUE, method = "jitter", jitter = 0.12, pch = 16, cex = 0.55, col = adjustcolor(INK, 0.45))
    iw <- d$rho[d$pair_type == "intra-KO (replicates)"]; cr <- d$rho[d$pair_type == "cross-KO"]; wt <- wilcox.test(iw, cr)
    mtext(sp, side = 3, adj = 0, line = 0.3, font = 2, cex = 0.82)
    mtext(sprintf("intra %.3f vs cross %.3f (%s)", median(iw), median(cr), if (wt$p.value < 0.001) "p<0.001" else sprintf("p=%.2g", wt$p.value)), side = 3, adj = 0, line = -0.7, cex = 0.55, col = INK2) }
  mtext("NCI-H2052 bulk: intra-KO (replicate) vs cross-KO vs vs-NTC sample-sample Spearman correlation", outer = TRUE, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext("intra-KO = same-KO replicate pairs (n=8); cross-KO = different-KO pairs excluding NTC (n=120); vs NTC = any KO sample vs the single NTC sample (n=16)", outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.55, col = INK2)
})
