# Spec 57: per-TF view of spec56: Spearman correlation (NTC - KO) between the CROP-seq screen and bulk NCI-H2052 for every TF (columns) under every one of the 15
# NTC combinations (rows; H2052_NTC_R3 + the QC-dropped H28_NTC_R1/2/3), one panel per feature space (20 Cm modules, top 1000/2000/3000 variable genes).
# Figure A = own-KO rho; figure B = separation (own-KO rho minus median rho of the screen KO vs the other H2052 KOs).
# needs: spec56_H2052_NTC_subset_search.R (spec56_NTC_subset_scores.csv)
# outputs: fig3O_H2052_NTC_combos_rho_per_TF.{png,pdf}, fig3P_H2052_NTC_combos_separation_per_TF.{png,pdf}, spec57_rho_wide.csv
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
R <- read.csv(file.path(OUT, "spec56_NTC_subset_scores.csv"), stringsAsFactors = FALSE)
spaces <- unique(R$space); subs <- unique(R[, c("subset", "n_ntc")]); subs <- subs[order(subs$n_ntc, subs$subset), ]
kos <- names(sort(tapply(R$own, R$KO, mean), decreasing = TRUE))   # TFs ordered by mean own-KO rho
write.csv(reshape(R[, c("space", "subset", "KO", "own")], idvar = c("space", "subset"), timevar = "KO", direction = "wide"), file.path(OUT, "spec57_rho_wide.csv"), row.names = FALSE)
mk <- function(name, col, lab, ttl) save_fig(name, 14, 12.5, function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.4, 0), family = "sans"); vmax <- max(abs(R[[col]]), na.rm = TRUE)
  for (sp in spaces) { d <- R[R$space == sp, ]; M <- sapply(kos, function(k) sapply(subs$subset, function(s) d[[col]][d$KO == k & d$subset == s])); rownames(M) <- subs$subset
    par(mar = c(4.2, 12, 2.4, 0.6)); n <- nrow(M)
    image(seq_along(kos), seq_len(n), t(M)[, n:1], col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(n)) for (j in seq_along(kos)) text(j, n - i + 1, sprintf("%.2f", M[i, j]), cex = 0.62, col = INK)
    axis(1, at = seq_along(kos), labels = FALSE); text(seq_along(kos), 0.2, kos, srt = 45, adj = 1, xpd = NA, cex = 0.8)
    text(0.4, n:1, rownames(M), adj = 1, xpd = NA, cex = 0.62, col = ifelse(rownames(M) == "H2052R3", "#eda100", ifelse(subs$n_ntc == 4, "#4a3aa7", INK2)), font = ifelse(rownames(M) %in% c("H2052R3") | subs$n_ntc == 4, 2, 1))
    gb <- which(diff(subs$n_ntc) != 0); for (g in gb) abline(h = n - g + 0.5, col = "white", lwd = 2.5); box(col = MUTED, lwd = 0.6)
    mtext(sprintf("%s (mean %s = %.3f)", sp, lab, mean(d[[col]])), side = 3, adj = 0, line = 0.3, font = 2, cex = 0.85) }
  mtext(ttl, outer = TRUE, side = 3, adj = 0, line = 1.9, font = 2, cex = 1.0)
  mtext("rows = H2052 NTC combination (orange = H2052_NTC_R3 only, purple = all 4; white lines separate 1/2/3/4-sample sets); columns = TF; Spearman, NTC-KO", outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.6, col = INK2) })
mk("fig3O_H2052_NTC_combos_rho_per_TF", "own", "rho", "Screen vs bulk H2052, own-KO Spearman rho per TF for all 15 H2052 NTC combinations")
mk("fig3P_H2052_NTC_combos_separation_per_TF", "sep", "separation", "Screen vs bulk H2052, separation (own-KO rho minus other-KO rho) per TF for all 15 H2052 NTC combinations")
options(width = 200); cat("mean own-KO rho per TF over all 15 combinations x 4 spaces:\n"); print(round(sort(tapply(R$own, R$KO, mean), decreasing = TRUE), 3))
cat("\nmean separation per TF:\n"); print(round(sort(tapply(R$sep, R$KO, mean), decreasing = TRUE), 3))
cat("\nfraction of (combination x space) cells with own rho > 0, per TF:\n"); print(round(tapply(R$own > 0, R$KO, mean), 2)[kos])
