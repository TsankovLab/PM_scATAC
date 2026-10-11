# Spec 84: per-TF agreement on two axes, side by side, same TF order (ranked by cross-line agreement, spec82):
#   left  bulk vs CROP-seq: Spearman over the TF's 20 Cm module effects, screen vs bulk (4 lines combined, spec75 S70, cell cycle adjusted on both sides);
#         KO vs NTC (screen median difference vs bulk Stouffer z) and KO vs the other KOs (KO-specific, no NTC); p = cor.test (n = 20 modules)
#   right bulk vs bulk: mean same-KO Spearman over the 6 line pairs (KO-centred profiles, top-500 KO-signal genes; spec76 S70), dots = line pairs; p = KO-label permutation
# needs: spec75_Cm_preserved.R (SET=S70), spec76_final_figures.R (SET=S70), spec82_TF_agreement_ranked.R
# outputs: spec84_TF_agreement_two_axes.csv, fig9F_TF_agreement_two_axes.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
rd <- function(f) read.csv(file.path(OUT, f), stringsAsFactors = FALSE, check.names = FALSE)
tab <- rd("spec75_Cm_KO_table_S70.csv"); pk <- rd("spec76_primary_per_KO_S70.csv"); pk <- pk[pk$space == "signal500", ]; o <- rd("spec82_TF_agreement_ranked.csv")$KO
cp <- do.call(rbind, lapply(o, function(k) { d <- tab[tab$KO == k, ]; a <- cor.test(d$screen_effect, d$bulk_z, method = "spearman", exact = FALSE); b <- cor.test(d$screen_rest_effect, d$bulk_rest_z, method = "spearman", exact = FALSE)
  data.frame(KO = k, crop_rho_vsNTC = unname(a$estimate), crop_p_vsNTC = a$p.value, crop_rho_specific = unname(b$estimate), crop_p_specific = b$p.value,
             bulk_rho_none = pk$mean_intra[pk$KO == k & pk$cov == "none"], bulk_p_none = pk$perm_p[pk$KO == k & pk$cov == "none"], bulk_fdr_none = pk$perm_fdr[pk$KO == k & pk$cov == "none"],
             bulk_rho_cc = pk$mean_intra[pk$KO == k & pk$cov == "cc"], bulk_p_cc = pk$perm_p[pk$KO == k & pk$cov == "cc"], bulk_fdr_cc = pk$perm_fdr[pk$KO == k & pk$cov == "cc"]) }))
cp$crop_fdr_vsNTC <- p.adjust(cp$crop_p_vsNTC, "BH"); cp$crop_fdr_specific <- p.adjust(cp$crop_p_specific, "BH")
write.csv(cp, file.path(OUT, "spec84_TF_agreement_two_axes.csv"), row.names = FALSE); options(width = 200); print(cp, digits = 3, row.names = FALSE)
intra <- lapply(c(none = "none", cc = "cc"), function(cv) { x <- rd(sprintf("spec76_KOline_correlation_signal500_%s_S70.csv", cv)); M <- as.matrix(x[, -1]); dimnames(M) <- list(x[, 1], x[, 1])
  setNames(lapply(o, function(k) sapply(PAIRS, function(p) M[paste(k, sub("NCI-|-211H", "", p[1])), paste(k, sub("NCI-|-211H", "", p[2]))])), o) })

fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
bars <- function(V, P, Fd, cols, ylab, ylim, dots = NULL) { n <- length(o); at <- outer(c(-0.19, 0.19), seq_len(n), "+")
  plot(NA, xlim = c(0.5, n + 0.5), ylim = ylim, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = ylab, cex.axis = 0.8); abline(h = 0, col = INK, lwd = 0.8)
  for (j in 1:2) for (i in seq_len(n)) { x <- at[j, i]; rect(x - 0.17, 0, x + 0.17, V[i, j], col = adjustcolor(cols[j], 0.6), border = cols[j], lwd = 0.8)
    if (!is.null(dots)) points(x + seq(-0.1, 0.1, length.out = 6), dots[[j]][[o[i]]], pch = 16, cex = 0.55, col = INK2)
    text(x, ylim[2] - (j - 1) * 0.07 * diff(ylim), paste0(fp(P[i, j]), ifelse(Fd[i, j] < 0.01, " **", ifelse(Fd[i, j] < 0.05, " *", ""))), cex = 0.6, col = cols[j], font = ifelse(P[i, j] < 0.05, 2, 1)) }
  text(seq_len(n), ylim[1] - 0.05 * diff(ylim), o, xpd = NA, cex = 0.82, font = 2, col = KO_COL[o], adj = c(0.5, 1)) }
save_fig("fig9F_TF_agreement_two_axes", 15, 5.6, function() {
  layout(matrix(1:2, nrow = 1)); par(oma = c(2.2, 0, 3.4, 0.4), mar = c(3.4, 4.6, 2.6, 1), family = "sans")
  bars(as.matrix(cp[, c("crop_rho_vsNTC", "crop_rho_specific")]), as.matrix(cp[, c("crop_p_vsNTC", "crop_p_specific")]), as.matrix(cp[, c("crop_fdr_vsNTC", "crop_fdr_specific")]),
       c("grey50", "#eb6834"), "Spearman, screen vs bulk (over the 20 Cm modules)", c(-1, 1.15))
  mtext("A  bulk vs CROP-seq", side = 3, adj = 0, line = 1, font = 2, cex = 0.9)
  legend("bottom", inset = c(0, -0.22), xpd = NA, c("KO vs NTC", "KO vs other KOs (KO-specific)"), fill = adjustcolor(c("grey50", "#eb6834"), 0.6), border = c("grey50", "#eb6834"), bty = "n", horiz = TRUE, cex = 0.72)
  yl <- range(c(unlist(intra), 0)); yl <- yl + c(-0.02, 0.2) * diff(yl)
  bars(as.matrix(cp[, c("bulk_rho_none", "bulk_rho_cc")]), as.matrix(cp[, c("bulk_p_none", "bulk_p_cc")]), as.matrix(cp[, c("bulk_fdr_none", "bulk_fdr_cc")]),
       c("#2a78d6", "#1baf7a"), "same-KO Spearman between two bulk lines", yl, dots = intra)
  mtext("B  bulk vs bulk (across the 4 cell lines)", side = 3, adj = 0, line = 1, font = 2, cex = 0.9)
  legend("bottom", inset = c(0, -0.22), xpd = NA, c("no covariate", "cell cycle regressed"), fill = adjustcolor(c("#2a78d6", "#1baf7a"), 0.6), border = c("#2a78d6", "#1baf7a"), bty = "n", horiz = TRUE, cex = 0.72)
  mtext("Agreement of each TF's KO effect: with the CROP-seq screen (A) and between bulk cell lines (B); TFs ordered by B", outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 1)
  mtext("numbers = p (A: Spearman test over 20 modules; B: KO-label permutation), bold = p < 0.05, * / ** = BH FDR < 0.05 / 0.01 over the 8 TFs; B dots = the 6 line pairs; S70, NTC - KO", outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.6, col = INK2)
})
