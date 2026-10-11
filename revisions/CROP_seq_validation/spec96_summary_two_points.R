# Spec 96: the two headline results in one figure, from saved outputs (no recomputation).
# Row 1 -- NCI-H2052, CROP-seq vs bulk (same cell line):
#   A shared response: module average over the 8 KOs, screen vs bulk (20 Cm; spec81)      B KO-specific part: each KO x Cm cell minus its module average (160 cells)
#   C per TF, gene level: own-KO Spearman (t vs other KOs, genes in the top 2000 KO-signal genes of both data sets; spec95), dash = median vs the other bulk KOs,
#     # = rank of the own KO among the 8 bulk KOs (per-TF empirical p = rank / 8); overall p = exact permutation of the bulk KO labels
# Row 2 -- bulk, across the 4 cell lines (S70, KO-centred NTC - KO t, top-500 KO-signal genes; spec76):
#   D / E per TF: same KO in two lines (6 pairs) vs other KOs, permutation p (BH stars), no covariate / cell cycle regressed      F all KOs pooled, global permutation p
# needs: spec81, spec95, spec76 (SET=S70), spec82 (TF order)
# outputs: fig10_summary_crop_vs_bulk_and_across_lines.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
rd <- function(f) read.csv(file.path(OUT, f), stringsAsFactors = FALSE, check.names = FALSE); fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
x <- rd("spec81_Cm_KO_NTCminusKO_screen_and_bulk.csv"); W <- function(s) { d <- x[x$source == s, ]; M <- tapply(d$NTC_minus_KO, list(d$KO, d$Cm), mean); M[ORDER, paste0("Cm", 1:20)] }
S <- W("screen"); B <- W("NCI-H2052"); sh <- cor.test(colMeans(S), colMeans(B), method = "spearman", exact = FALSE); Ss <- sweep(S, 2, colMeans(S)); Bs <- sweep(B, 2, colMeans(B)); spc <- cor.test(as.vector(Ss), as.vector(Bs), method = "spearman", exact = FALSE)
SPACE <- "signal_top2000"; k <- rd("spec95_shared_space_per_KO.csv"); k <- k[k$space == SPACE & k$view == "t", ]; g <- rd("spec95_shared_space.csv"); g <- g[g$space == SPACE & g$view == "t", ]
pk <- rd("spec76_primary_per_KO_S70.csv"); pk <- pk[pk$space == "signal500", ]; gl <- rd("spec76_primary_global_S70.csv"); gl <- gl[gl$space == "signal500", ]
pairs_of <- function(cv) { z <- rd(sprintf("spec76_KOline_correlation_signal500_%s_S70.csv", cv)); M <- as.matrix(z[, -1]); rn <- z[, 1]; dimnames(M) <- list(rn, rn); ko <- sub(" .*", "", rn); ln <- sub(".* ", "", rn)
  ut <- which(upper.tri(M) & outer(ln, ln, "!="), arr.ind = TRUE); data.frame(a = ko[ut[, 1]], b = ko[ut[, 2]], rho = M[ut]) }
PR <- list(none = pairs_of("none"), cc = pairs_of("cc")); COVL <- c(none = "no covariate", cc = "cell cycle regressed")

save_fig("fig10_summary_crop_vs_bulk_and_across_lines", 15, 10, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE), widths = c(1, 1, 1.25)); par(oma = c(0, 0, 4.6, 0.4), mar = c(4.6, 4.6, 3.4, 1), family = "sans")
  ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.82); mtext(b, side = 3, adj = 0, line = 0.4, cex = 0.6, col = INK2) }
  ## A
  a <- colMeans(S); b <- colMeans(B); plot(a, b, pch = 16, col = INK2, cex = 1.4, bty = "n", las = 1, xlab = "CROP-seq: module average over the 8 KOs (NTC - KO)", ylab = "bulk H2052: module average over the 8 KOs", cex.axis = 0.8); abline(h = 0, v = 0, col = GRID); abline(lm(b ~ a), col = "#eb6834", lwd = 1.5)
  text(a, b, names(a), pos = 3, cex = 0.55, col = INK2, xpd = NA); ttl("A  shared response to any KO: reproduced", sprintf("20 Cm modules: Spearman %.2f, p = %s", sh$estimate, fp(sh$p.value)))
  ## B
  plot(as.vector(Ss), as.vector(Bs), pch = 16, col = adjustcolor(KO_COL[rep(ORDER, 20)], 0.75), cex = 1, bty = "n", las = 1, xlab = "CROP-seq: KO-specific module effect", ylab = "bulk H2052: KO-specific module effect", cex.axis = 0.8); abline(h = 0, v = 0, col = GRID)
  legend("topleft", ORDER, pch = 16, col = KO_COL[ORDER], bty = "n", cex = 0.55, ncol = 2); ttl("B  KO-specific module effects: not reproduced", sprintf("160 KO x Cm cells minus module average: Spearman %.2f, p = %s", spc$estimate, fp(spc$p.value)))
  ## C
  o <- k$KO[order(-k$own)]; d <- k[match(o, k$KO), ]; par(mar = c(4.6, 4.6, 3.4, 1)); bp <- barplot(d$own, col = adjustcolor(KO_COL[o], 0.8), border = NA, las = 1, ylim = c(min(d$own) - 0.08, max(d$own) + 0.22), ylab = "Spearman, CROP-seq KO vs bulk H2052 KO", names.arg = rep("", 8), cex.axis = 0.8)
  abline(h = 0, col = INK); segments(bp - 0.4, d$median_other, bp + 0.4, d$median_other, lty = 2, col = "grey25", lwd = 1.3)
  text(bp, pmax(d$own, 0) + 0.02, sprintf("#%d\np=%s", d$rank, fp(d$rank / 8)), cex = 0.62, adj = c(0.5, 0), font = ifelse(d$rank == 1, 2, 1)); text(bp, par("usr")[3] - 0.02 * diff(par("usr")[3:4]), o, srt = 45, adj = 1, xpd = NA, cex = 0.82, font = 2, col = KO_COL[o])
  legend("topright", "median vs the other 7 bulk KOs", lty = 2, col = "grey25", bty = "n", cex = 0.65)
  ttl("C  per TF, gene level: only TEAD4 and PITX1 match", sprintf("%d shared KO-signal genes; all 8 KOs: permutation p = %s; # = rank of own KO among 8", g$n_genes, fp(g$perm_p)))
  ## D, E
  o2 <- rd("spec82_TF_agreement_ranked.csv")$KO
  for (cv in names(COVL)) { d <- PR[[cv]]; kk <- pk[pk$cov == cv, ]; yr <- range(d$rho) + c(-0.02, 0.16) * diff(range(d$rho)); par(mar = c(4.6, 4.6, 3.4, 1))
    plot(NA, xlim = c(0.5, 8.5), ylim = yr, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "Spearman between two bulk cell lines", cex.axis = 0.8); abline(h = 0, lty = 2, col = INK)
    for (i in seq_along(o2)) { w <- d$rho[d$a == o2[i] & d$b == o2[i]]; z <- d$rho[(d$a == o2[i]) != (d$b == o2[i])]
      boxplot(z, at = i - 0.2, add = TRUE, boxwex = 0.32, col = "#d5d4cf", border = INK2, outline = FALSE, axes = FALSE); boxplot(w, at = i + 0.2, add = TRUE, boxwex = 0.32, col = adjustcolor(KO_COL[[o2[i]]], 0.65), border = INK2, outline = FALSE, axes = FALSE)
      points(i + 0.2 + seq(-0.1, 0.1, length.out = length(w)), w, pch = 16, cex = 0.75, col = INK2); r <- kk[kk$KO == o2[i], ]
      text(i, yr[2], sprintf("p=%s\n%s", fp(r$perm_p), ifelse(r$perm_fdr < 0.01, "**", ifelse(r$perm_fdr < 0.05, "*", ""))), cex = 0.68, adj = c(0.5, 1), font = ifelse(r$perm_p < 0.05, 2, 1)) }
    text(1:8, par("usr")[3] - 0.03 * diff(par("usr")[3:4]), o2, srt = 45, adj = 1, xpd = NA, cex = 0.82, font = 2, col = KO_COL[o2]); gg <- gl[gl$cov == cv, ]
    ttl(sprintf("%s  per TF across the 4 bulk lines (%s)", if (cv == "none") "D" else "E", COVL[[cv]]), sprintf("colour = same KO (6 line pairs), grey = other KOs; global p = %s; * / ** FDR < 0.05 / 0.01", fp(gg$perm_p))) }
  ## F
  v <- list(PR$none$rho[PR$none$a == PR$none$b], PR$none$rho[PR$none$a != PR$none$b], PR$cc$rho[PR$cc$a == PR$cc$b], PR$cc$rho[PR$cc$a != PR$cc$b]); at <- c(1, 1.8, 3.2, 4)
  boxplot(v, at = at, col = c("#eb6834", "#d5d4cf", "#eb6834", "#d5d4cf"), border = INK2, outline = FALSE, xaxt = "n", las = 1, ylab = "Spearman between two bulk cell lines", cex.axis = 0.8, boxwex = 0.6, ylim = range(unlist(v)) + c(0, 0.08)); abline(h = 0, lty = 2, col = INK)
  axis(1, at = at, labels = rep(c("same KO", "other KO"), 2), tick = FALSE, cex.axis = 0.72); mtext(COVL, side = 1, at = c(1.4, 3.6), line = 2.4, cex = 0.62, font = 2)
  text(c(1.4, 3.6), par("usr")[4] - 0.02, sprintf("p = %s", fp(gl$perm_p[match(c("none", "cc"), gl$cov)])), cex = 0.8, font = 2)
  ttl("F  all KOs pooled", sprintf("same-KO %.3f vs other-KO %.3f (no covariate)", gl$diag_mean[gl$cov == "none"], gl$off_mean[gl$cov == "none"]))
  mtext("Top: NCI-H2052, CROP-seq vs bulk RNA-seq (same cell line).   Bottom: bulk RNA-seq, reproducibility of each TF KO across the 4 cell lines", outer = TRUE, side = 3, adj = 0, line = 3, font = 2, cex = 1.05)
  mtext("NTC - KO effects; A-B: Cm module scores; C: t of each KO vs the other KOs; D-F: KO-centred limma t, top-500 KO-signal genes chosen from the 2 lines not compared, 20,000 label permutations (S70)", outer = TRUE, side = 3, adj = 0, line = 1.7, cex = 0.65, col = INK2)
})
