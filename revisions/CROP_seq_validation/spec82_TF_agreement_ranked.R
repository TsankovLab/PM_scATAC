# Spec 82: TFs ranked by how well their bulk KO profile agrees across the 4 cell lines (S70 primary configuration, top-500 KO-signal genes, from spec76 outputs).
# Bar = mean same-KO Spearman over the 6 line pairs, dots = the 6 pairs; grey dashed tick = median correlation with the OTHER KOs (null level); label = permutation p (BH FDR stars).
# Shown without covariate and with cell cycle regressed; TFs ordered by their best (smallest) p over the two.
# needs: spec76_final_figures.R (SET=S70)
# outputs: spec82_TF_agreement_ranked.csv, fig9D_TF_agreement_ranked.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
COVL <- c(none = "no covariate", cc = "cell cycle regressed"); CC2 <- c(none = "#2a78d6", cc = "#eb6834")
pk <- read.csv(file.path(OUT, "spec76_primary_per_KO_S70.csv")); pk <- pk[pk$space == "signal500" & pk$cov %in% names(COVL), ]
intra <- lapply(setNames(names(COVL), names(COVL)), function(cv) { f <- file.path(OUT, sprintf("spec76_KOline_correlation_signal500_%s_S70.csv", cv)); x <- read.csv(f, check.names = FALSE); M <- as.matrix(x[, -1]); dimnames(M) <- list(x[, 1], x[, 1])
  sl <- sub("NCI-|-211H", "", LINES); setNames(lapply(ORDER, function(k) sapply(PAIRS, function(p) M[paste(k, sub("NCI-|-211H", "", p[1])), paste(k, sub("NCI-|-211H", "", p[2]))])), ORDER) })
best <- tapply(pk$perm_p, pk$KO, min); o <- names(sort(best)); o <- o[order(best[o], -tapply(pk$mean_intra, pk$KO, max)[o])]
tab <- do.call(rbind, lapply(o, function(k) { d <- pk[pk$KO == k, ]; data.frame(rank = match(k, o), KO = k, p_none = d$perm_p[d$cov == "none"], fdr_none = d$perm_fdr[d$cov == "none"], rho_none = d$mean_intra[d$cov == "none"],
  p_cc = d$perm_p[d$cov == "cc"], fdr_cc = d$perm_fdr[d$cov == "cc"], rho_cc = d$mean_intra[d$cov == "cc"], pairs_pos_none = d$n_pos_pairs[d$cov == "none"], pairs_pos_cc = d$n_pos_pairs[d$cov == "cc"]) }))
write.csv(tab, file.path(OUT, "spec82_TF_agreement_ranked.csv"), row.names = FALSE); print(tab, digits = 3, row.names = FALSE)

save_fig("fig9D_TF_agreement_ranked", 9.6, 5.6, function() {
  par(mar = c(5.6, 4.6, 3.6, 1), family = "sans"); n <- length(o); at <- outer(c(-0.19, 0.19), seq_len(n), "+")
  yl <- range(c(unlist(intra), 0)); yl <- yl + c(-0.02, 0.16) * diff(yl)
  plot(NA, xlim = c(0.5, n + 0.5), ylim = yl, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "same-KO Spearman rho between two cell lines", cex.axis = 0.8); abline(h = 0, col = INK, lwd = 0.8)
  rect(seq_len(n) - 0.48, yl[1], seq_len(n) + 0.48, yl[2], col = ifelse(seq_len(n) <= 3, adjustcolor("#1baf7a", 0.08), NA), border = NA)
  for (j in seq_along(COVL)) { cv <- names(COVL)[j]
    for (i in seq_len(n)) { k <- o[i]; r <- pk[pk$KO == k & pk$cov == cv, ]; x <- at[j, i]
      rect(x - 0.17, 0, x + 0.17, r$mean_intra, col = adjustcolor(CC2[[cv]], 0.55), border = CC2[[cv]], lwd = 0.8)
      segments(x - 0.17, r$median_cross, x + 0.17, r$median_cross, col = "grey40", lty = 2, lwd = 1)
      points(x + seq(-0.1, 0.1, length.out = 6), intra[[cv]][[k]], pch = 16, cex = 0.65, col = INK2)
      text(x, yl[2] - (j - 1) * 0.06 * diff(yl), sprintf("%s%s", formatC(r$perm_p, format = "g", digits = 2), stars(r$perm_fdr) |> (\(s) if (s == "ns") "" else paste0(" ", s))()),
           cex = 0.62, col = CC2[[cv]], font = ifelse(r$perm_p < 0.05, 2, 1), srt = 0) } }
  text(seq_len(n), par("usr")[3] - 0.04 * diff(yl), o, srt = 0, adj = c(0.5, 1), xpd = NA, cex = 0.85, font = 2, col = KO_COL[o])
  legend("bottom", inset = c(0, -0.2), xpd = NA, legend = c(COVL, "median vs other KOs (null)"), fill = c(adjustcolor(CC2, 0.55), NA), border = c(CC2, NA), lty = c(NA, NA, 2), col = c(NA, NA, "grey40"), bty = "n", cex = 0.65, ncol = 3)
  mtext("TFs ranked by cross-line agreement of their bulk KO profile (S70)", side = 3, adj = 0, line = 2.6, font = 2, cex = 0.95)
  mtext("bar = mean same-KO rho over the 6 line pairs (dots); numbers = permutation p (blue: no covariate, orange: cell cycle regressed), * FDR<0.05 ** <0.01; shaded = reproducible TFs", side = 3, adj = 0, line = 1.5, cex = 0.58, col = INK2)
})
