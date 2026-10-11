# Spec 44: screen-to-bulk agreement for the ORIGINAL/established 20 "Cm" modules (the ones used throughout this whole pipeline --
# cm_sample_scores_AddModuleScore_bulk_RUN1.csv, screen_cm_diff_matrix.csv -- predating and distinct from this session's own cNMF
# k=20/30 sweep exploration). Mirrors spec39/43's methodology (per-KO Pearson r across modules; "down in screen -> down in bulk?"
# agreement) but using these pre-existing, already-computed Cm scores rather than recomputing from raw cNMF output.
# needs: screen_cm_diff_matrix.csv (pre-existing, screen KO-vs-NTC diff, NCI-H2052 only, all 20 Cm), cm_sample_scores_AddModuleScore_bulk_RUN1.csv
#        (pre-existing, per-bulk-sample Cm scores), repro_common.R (load_bulk(), the curated sample sheet these Cm scores were built against)
# outputs: spec44_20Cm_agreement.csv, fig3C_20Cm_screen_bulk_correlation.{png,pdf}, fig3C_20Cm_down_agreement.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
PP <- file.path(PROC, "de", "perturbation_panel")

screen_diff <- read.csv(file.path(PP, "screen_cm_diff_matrix.csv"), row.names = 1, check.names = FALSE)
# sign convention check: this pre-existing file's own companion figure (screen_vs_bulk_cm_heatmaps.png) is captioned "NTC - KO
# module score" -- i.e. stored value = NTC minus KO, the OPPOSITE of the KO-minus-NTC convention used everywhere else in this
# pipeline (spec34-43, cm_diff below). Flip it here so both sides of every comparison mean the same thing (positive = KO > NTC).
screen_diff <- -screen_diff
CM <- colnames(screen_diff)
b <- load_bulk(); meta <- b$meta
sc_df <- read.csv(file.path(PP, "cm_sample_scores_AddModuleScore_bulk_RUN1.csv"), check.names = FALSE)
rownames(sc_df) <- sc_df$sample
common_samples <- intersect(rownames(meta), rownames(sc_df))
cat(sprintf("bulk Cm scores available for %d/%d curated samples\n", length(common_samples), nrow(meta)))
meta <- meta[common_samples, ]; sc_df <- sc_df[common_samples, ]

ntc_s <- rownames(meta)[meta$Gene_Control == "NTC"]
bulk_diff <- t(sapply(ORDER, function(k) {
  ko_s <- rownames(meta)[meta$Gene_Control == k]
  colMeans(sc_df[ko_s, CM, drop = FALSE]) - colMeans(sc_df[ntc_s, CM, drop = FALSE])
}))
rownames(bulk_diff) <- ORDER

common_ko <- intersect(rownames(screen_diff), ORDER)
res <- do.call(rbind, lapply(common_ko, function(k) {
  s <- as.numeric(screen_diff[k, CM]); b_ <- bulk_diff[k, CM]
  pr <- cor.test(s, b_, method = "pearson")
  sr <- suppressWarnings(cor.test(s, b_, method = "spearman", exact = FALSE))
  down_cm <- CM[s < 0]
  down_agree <- mean(b_[down_cm] < 0)
  data.frame(KO = k, pearson_r = pr$estimate, pearson_p = pr$p.value, spearman_rho = sr$estimate, spearman_p = sr$p.value,
             screen_top_Cm = CM[which.max(s)], bulk_top_Cm = CM[which.max(b_)],
             n_screen_down_Cm = length(down_cm), n_bulk_down_agree = sum(b_[down_cm] < 0), frac_down_agree = down_agree)
}))
res <- res[order(-res$spearman_rho), ]
write.csv(res, file.path(OUT, "spec44_20Cm_agreement.csv"), row.names = FALSE)

options(width = 200)
cat("\n=== per-KO screen-vs-bulk agreement, the 20 established Cm modules ===\n")
print(format(res, digits = 3), row.names = FALSE)
cat(sprintf("\nPearson:  mean r = %.3f | median = %.3f | KOs with r > 0: %d/%d | significant (p<0.05): %d/%d\n",
            mean(res$pearson_r), median(res$pearson_r), sum(res$pearson_r > 0), nrow(res), sum(res$pearson_p < 0.05), nrow(res)))
cat(sprintf("Spearman: mean rho = %.3f | median = %.3f | KOs with rho > 0: %d/%d | significant (p<0.05): %d/%d\n",
            mean(res$spearman_rho), median(res$spearman_rho), sum(res$spearman_rho > 0), nrow(res), sum(res$spearman_p < 0.05), nrow(res)))
all_pairs <- do.call(rbind, lapply(common_ko, function(k) data.frame(KO = k, Cm = CM, screen = as.numeric(screen_diff[k, CM]), bulk = bulk_diff[k, CM])))
down_pairs <- all_pairs[all_pairs$screen < 0, ]
cat(sprintf("overall: of %d KO x Cm pairs down in screen, %d (%.0f%%) also down in bulk [chance = 50%%]\n",
            nrow(down_pairs), sum(down_pairs$bulk < 0), 100 * mean(down_pairs$bulk < 0)))
bt <- binom.test(sum(down_pairs$bulk < 0), nrow(down_pairs), p = 0.5)
cat(sprintf("binomial test vs 50%% chance: p = %.4g\n", bt$p.value))

## ---------------------------------------------------------------- per-line correlation: same Pearson r, computed separately within each cell line's own bulk samples
bulk_diff_line <- lapply(LINES, function(cl) {
  mm <- meta[meta$Cell_Line == cl, ]; ntc_l <- rownames(mm)[mm$Gene_Control == "NTC"]
  t(sapply(ORDER, function(k) { ko_l <- rownames(mm)[mm$Gene_Control == k]
    if (!length(ko_l) || length(ntc_l) < 2) return(setNames(rep(NA, length(CM)), CM))   # NCI-H2052 excluded: only 1 NTC sample in the curated sheet, no replication
    colMeans(sc_df[ko_l, CM, drop = FALSE]) - colMeans(sc_df[ntc_l, CM, drop = FALSE]) }))
})
names(bulk_diff_line) <- LINES; for (cl in LINES) rownames(bulk_diff_line[[cl]]) <- ORDER

r_per_line <- do.call(rbind, lapply(LINES, function(cl) do.call(rbind, lapply(common_ko, function(k) {
  b_ <- bulk_diff_line[[cl]][k, CM]; if (all(is.na(b_))) return(NULL)
  s <- as.numeric(screen_diff[k, CM])
  data.frame(line = cl, KO = k, pearson_r = cor(s, b_, use = "complete.obs"),
             spearman_rho = suppressWarnings(cor(s, b_, method = "spearman", use = "complete.obs")))
}))))
write.csv(r_per_line, file.path(OUT, "spec44_20Cm_agreement_per_line.csv"), row.names = FALSE)
cat("\n=== per-line, per-KO Pearson r / Spearman rho (screen vs that line's own bulk, across 20 Cm modules) ===\n")
print(format(r_per_line, digits = 3), row.names = FALSE)
cat("\nmedian Pearson r per line:\n"); print(round(sapply(split(r_per_line$pearson_r, r_per_line$line), median), 3))
cat("median Spearman rho per line:\n"); print(round(sapply(split(r_per_line$spearman_rho, r_per_line$line), median), 3))

save_fig("fig3C_20Cm_correlation_per_line_boxplot", 7.5, 5.5, function() {
  ko_ord <- names(sort(sapply(split(r_per_line$spearman_rho, r_per_line$KO), median), decreasing = TRUE))
  r_per_line$KO <- factor(r_per_line$KO, levels = ko_ord)
  par(mar = c(5, 4.8, 3.6, 1), family = "sans")
  bx <- boxplot(spearman_rho ~ KO, data = r_per_line, col = adjustcolor(KO_COL[ko_ord], 0.35), border = INK2,
                outline = FALSE, las = 2, ylab = "Spearman rho (screen vs bulk, across 20 Cm modules)", xlab = "", ylim = c(-1, 1))
  set.seed(1); for (i in seq_along(ko_ord)) { d <- r_per_line[r_per_line$KO == ko_ord[i], ]
    points(jitter(rep(i, nrow(d)), amount = 0.15), d$spearman_rho, pch = 16, cex = 1, col = LINE_COL[d$line]) }
  abline(h = 0, col = INK, lty = 2)
  present_lines <- sort(unique(r_per_line$line))
  legend("bottomleft", legend = present_lines, col = LINE_COL[present_lines], pch = 16, bty = "n", cex = 0.65)
  mtext("Screen-to-bulk correlation (20 Cm modules), per TF (Spearman)", side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext("sorted by median (high to low); one point per cell line (NCI-H2052 excluded, only 1 NTC sample); boxplot = distribution across the 3 remaining lines", side = 3, adj = 0, line = 0.9, cex = 0.6, col = INK2)
})

## ---------------------------------------------------------------- fig: per-KO correlation bar
save_fig("fig3C_20Cm_screen_bulk_correlation", 8, 5.5, function() {
  ord <- res[order(-res$spearman_rho), ]
  cols <- ifelse(ord$spearman_p < 0.05, KO_COL[ord$KO], adjustcolor(KO_COL[ord$KO], 0.4))
  par(mar = c(6, 4.8, 3.6, 1), family = "sans")
  bp <- barplot(ord$spearman_rho, names.arg = ord$KO, col = cols, border = INK2, las = 2, cex.names = 0.85,
                ylab = "Spearman rho (screen vs bulk, across 20 Cm modules)", ylim = c(-1, 1))
  abline(h = 0, col = INK)
  mtext("Screen-to-bulk agreement, the 20 established Cm modules (Spearman)", side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext("solid = p<0.05; faded = n.s.; screen = NCI-H2052 only (Perturb-seq line); bulk = pooled across all 4 lines", side = 3, adj = 0, line = 0.9, cex = 0.6, col = INK2)
})

## ---------------------------------------------------------------- fig: down-in-screen -> down-in-bulk agreement, per KO
save_fig("fig3C_20Cm_down_agreement", 8, 5.5, function() {
  ord <- res[order(-res$frac_down_agree), ]
  par(mar = c(6, 4.8, 3.6, 1), family = "sans")
  bp <- barplot(ord$frac_down_agree, names.arg = ord$KO, col = adjustcolor(KO_COL[ord$KO], 0.85), border = INK2, las = 2, cex.names = 0.85,
                ylab = "fraction of screen-down Cm modules also down in bulk", ylim = c(0, 1))
  abline(h = 0.5, col = INK, lty = 2)
  text(bp, ord$frac_down_agree, sprintf("%d/%d", ord$n_bulk_down_agree, ord$n_screen_down_Cm), pos = 3, cex = 0.7, xpd = NA)
  mtext("The 20 Cm modules: of those down vs NTC in screen, how many also down in bulk?", side = 3, adj = 0, line = 2.2, font = 2, cex = 0.95)
  mtext("dashed line = 50% chance; numbers = Cm modules agreeing / Cm modules down in screen for that KO", side = 3, adj = 0, line = 0.9, cex = 0.6, col = INK2)
})
