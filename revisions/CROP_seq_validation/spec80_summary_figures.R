# Spec 80: two summary figures of the S70 primary analysis, from saved outputs (no recomputation).
#   fig9A  bulk vs CROP-seq on the 20 Cm modules: (A) screen effect vs bulk combined z for all 160 KO x Cm pairs, KO-specific preserved pairs labelled;
#          (B) per KO, Spearman over its 20 Cm between screen and bulk, vs NTC and vs the other KOs (KO-specific, no NTC);
#          (C) the preserved pairs: screen and every bulk line (signed z), with robustness over sample-set / covariate variants.
#   fig9B  bulk cell lines: same KO in two lines vs different KOs (KO-centred NTC - KO profiles, top-500 KO-signal genes), no covariate and cell cycle regressed;
#          per KO boxplots (points = the 6 line pairs), permutation p (BH FDR) per KO, and the global comparison.
# needs: spec75_Cm_preserved.R (S70, S70N4, S75, S70 noCC), spec76_final_figures.R (S70)
# outputs: fig9A_bulk_vs_cropseq.{png,pdf}, fig9B_bulk_within_vs_outside_KO.{png,pdf}, spec80_preserved_Cm.csv
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
rd <- function(f) read.csv(file.path(OUT, f), stringsAsFactors = FALSE, check.names = FALSE)
tab <- rd("spec75_Cm_KO_table_S70.csv"); pl <- rd("spec75_Cm_KO_per_line_S70.csv"); cmg <- readRDS(CM_RDS)
VAR <- c(S70 = "S70", S70N4 = "+3 H2052 NTC", S75 = "75 samples", S70_noCC = "no cell cycle")
vt <- lapply(names(VAR), function(v) { y <- rd(sprintf("spec75_Cm_KO_table_%s.csv", v)); setNames(y$preserved | y$preserved_specific, paste(y$KO, y$Cm)) }); names(vt) <- names(VAR)
key <- paste(tab$KO, tab$Cm); tab$n_variants <- rowSums(sapply(vt, function(v) v[key]))
sel <- tab[tab$preserved | (tab$KO == "TEAD4" & tab$Cm == "Cm20"), ]; sel <- sel[order(-sel$n_variants, sel$bulk_p), ]
sel$top_genes <- sapply(sel$Cm, function(c) paste(head(cmg[[c]], 5), collapse = ", ")); sel$after_KO <- ifelse(sel$screen_effect > 0, "down", "up")
write.csv(sel[, c("KO", "Cm", "top_genes", "after_KO", "screen_effect", "screen_fdr", "bulk_z", "bulk_fdr", "n_lines_same_sign", "bulk_rest_z", "bulk_rest_p", "preserved_specific", "n_variants")], file.path(OUT, "spec80_preserved_Cm.csv"), row.names = FALSE)

## ================================================================ fig9A bulk vs CROP-seq
save_fig("fig9A_bulk_vs_cropseq", 15.5, 5.6, function() {
  layout(matrix(1:3, nrow = 1), widths = c(1.05, 0.9, 1.35)); par(oma = c(0, 0, 3.2, 0.3), family = "sans"); ttl <- function(x) mtext(x, side = 3, adj = 0, line = 0.6, font = 2, cex = 0.8)
  ## A scatter
  par(mar = c(4.4, 4.4, 2.6, 0.8)); sig <- tab$screen_fdr < 0.05; ct <- cor.test(tab$screen_effect, tab$bulk_z, method = "spearman", exact = FALSE)
  plot(tab$screen_effect, tab$bulk_z, pch = 16, cex = ifelse(sig, 1.15, 0.7), col = adjustcolor(ifelse(sig, KO_COL[tab$KO], "grey75"), 0.8), bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8,
       xlab = "CROP-seq: module score NTC - KO (median)", ylab = "bulk: combined z over 4 lines (NTC - KO)"); abline(h = 0, v = 0, col = GRID)
  ps <- tab$preserved_specific; points(tab$screen_effect[ps], tab$bulk_z[ps], pch = 1, cex = 2.1, lwd = 1.5); text(tab$screen_effect[ps], tab$bulk_z[ps], paste(tab$KO[ps], tab$Cm[ps]), pos = 4, cex = 0.55, xpd = NA)
  legend("topleft", c(ORDER, "not significant in screen"), pch = 16, col = c(KO_COL[ORDER], "grey75"), bty = "n", cex = 0.58, ncol = 2)
  ttl("A  all 160 KO x Cm effects"); mtext(sprintf("Spearman %.2f (p %.2g); same sign in %.0f%% of %d screen-significant pairs; circled = preserved & KO-specific", ct$estimate, ct$p.value, 100 * mean(sign(tab$screen_effect[sig]) == sign(tab$bulk_z[sig])), sum(sig)), side = 3, adj = 0, line = -0.4, cex = 0.5, col = INK2)
  ## B per KO
  par(mar = c(5.2, 4.2, 2.6, 0.8)); pk <- do.call(rbind, lapply(ORDER, function(k) { d <- tab[tab$KO == k, ]; a <- cor.test(d$screen_effect, d$bulk_z, method = "spearman", exact = FALSE); b <- cor.test(d$screen_rest_effect, d$bulk_rest_z, method = "spearman", exact = FALSE)
    data.frame(KO = k, vs_NTC = a$estimate, p_NTC = a$p.value, vs_rest = b$estimate, p_rest = b$p.value) })); o <- order(-pk$vs_rest)
  bp <- barplot(t(as.matrix(pk[o, c("vs_NTC", "vs_rest")])), beside = TRUE, col = c("grey70", "#eb6834"), border = NA, las = 1, ylim = c(-1, 1), ylab = "Spearman over the 20 Cm (screen vs bulk)", names.arg = rep("", 8), cex.axis = 0.75); abline(h = 0, col = INK)
  pv <- t(as.matrix(pk[o, c("p_NTC", "p_rest")])); rv <- t(as.matrix(pk[o, c("vs_NTC", "vs_rest")])); text(bp, rv + ifelse(rv >= 0, 0.03, -0.03), ifelse(pv < 0.001, "***", ifelse(pv < 0.01, "**", ifelse(pv < 0.05, "*", ""))), adj = c(0.5, ifelse(rv >= 0, 0, 1)), cex = 0.8)
  text(colMeans(bp), -1.08, pk$KO[o], srt = 45, adj = 1, xpd = NA, cex = 0.78, font = 2, col = KO_COL[pk$KO[o]])
  legend("topright", c("KO vs NTC", "KO vs other KOs (KO-specific)"), fill = c("grey70", "#eb6834"), border = NA, bty = "n", cex = 0.65); ttl("B  per KO, over its 20 modules")
  ## C preserved pairs: screen + 4 lines
  par(mar = c(4.2, 13.5, 2.6, 5.5)); n <- nrow(sel); cols <- c("screen", LINES)
  Z <- cbind(screen = sign(sel$screen_effect) * pmin(-log10(sel$screen_fdr), 10) / 10 * 3, sapply(LINES, function(cl) { d <- pl[pl$line == cl, ]; d$z[match(paste(sel$KO, sel$Cm), paste(d$KO, d$Cm))] }))
  Zc <- pmax(pmin(Z, 3), -3); image(seq_len(ncol(Z)), seq_len(n), t(Zc)[, n:1, drop = FALSE], col = colorRampPalette(c("#2a78d6", "#f0efec", "#eb6834"))(101), zlim = c(-3, 3), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(n)) for (j in seq_len(ncol(Z))) if (is.finite(Z[i, j])) text(j, n - i + 1, if (j == 1) formatC(sel$screen_fdr[i], format = "e", digits = 0) else sprintf("%+.1f", Z[i, j]), cex = 0.55, col = INK)
  abline(v = 1.5, col = "white", lwd = 3); axis(1, at = seq_len(ncol(Z)), labels = c("CROP-seq\n(FDR)", sub("NCI-|-211H", "", LINES)), tick = FALSE, cex.axis = 0.7, padj = 0.5)
  labs <- sprintf("%s  %s  (%s after KO)", sel$KO, sel$Cm, sel$after_KO); text(0.4, n:1, labs, adj = 1, xpd = NA, cex = 0.72, font = 2, col = KO_COL[sel$KO]); text(0.4, n:1 - 0.3, sel$top_genes, adj = 1, xpd = NA, cex = 0.5, col = INK2)
  text(ncol(Z) + 0.6, n:1, sprintf("%d / 4 variants", sel$n_variants), adj = 0, xpd = NA, cex = 0.6, font = ifelse(sel$n_variants == 4, 2, 1)); box(col = MUTED, lwd = 0.6)
  ttl("C  preserved pairs: screen and each bulk line"); mtext("bulk per-line z (orange = module lower after KO); right: variants where preserved", side = 3, adj = 0, line = -0.4, cex = 0.5, col = INK2)
  mtext("How well bulk RNA-seq reproduces the CROP-seq KO effects on the 20 Cm modules (S70, NTC - KO, cell cycle adjusted)", outer = TRUE, side = 3, adj = 0, line = 1.6, font = 2, cex = 1)
  mtext("screen: AddModuleScore on cells, Wilcoxon; bulk: module score ~ group + cell cycle + % assigned (limma) per line, Stouffer over lines; KO-specific = same test against the other KOs (no NTC)", outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2)
})

## ================================================================ fig9B bulk cell lines: within KO vs outside KO
pk76 <- rd("spec76_primary_per_KO_S70.csv"); gl76 <- rd("spec76_primary_global_S70.csv"); COVL <- c(none = "no covariate", cc = "cell cycle regressed")
save_fig("fig9B_bulk_within_vs_outside_KO", 14, 5.4, function() {
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 0.42)); par(oma = c(0, 0, 3.2, 0.3), family = "sans"); ttl <- function(x) mtext(x, side = 3, adj = 0, line = 0.6, font = 2, cex = 0.8); G <- list()
  for (cv in names(COVL)) { M <- as.matrix(rd(sprintf("spec76_KOline_correlation_signal500_%s_S70.csv", cv))[, -1]); rn <- rd(sprintf("spec76_KOline_correlation_signal500_%s_S70.csv", cv))[, 1]; dimnames(M) <- list(rn, rn)
    ko <- sub(" .*", "", rn); ln <- sub(".* ", "", rn); ut <- which(upper.tri(M) & outer(ln, ln, "!="), arr.ind = TRUE)
    d <- data.frame(a = ko[ut[, 1]], b = ko[ut[, 2]], rho = M[ut]); G[[cv]] <- d
    k <- pk76[pk76$cov == cv & pk76$space == "signal500", ]; o <- k$KO[order(k$perm_p, -k$mean_intra)]
    par(mar = c(4.6, 4.4, 2.6, 0.8)); yr <- range(d$rho) + c(-0.02, 0.1) * diff(range(d$rho))
    plot(NA, xlim = c(0.5, 8.5), ylim = yr, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "Spearman rho between two cell lines", cex.axis = 0.75); abline(h = 0, lty = 2, col = INK)
    for (i in seq_along(o)) { w <- d$rho[d$a == o[i] & d$b == o[i]]; x <- d$rho[(d$a == o[i]) != (d$b == o[i])]
      boxplot(x, at = i - 0.2, add = TRUE, boxwex = 0.32, col = "#d5d4cf", border = INK2, outline = FALSE, axes = FALSE)
      boxplot(w, at = i + 0.2, add = TRUE, boxwex = 0.32, col = adjustcolor(KO_COL[[o[i]]], 0.55), border = INK2, outline = FALSE, axes = FALSE)
      points(i + 0.2 + seq(-0.1, 0.1, length.out = length(w)), w, pch = 16, cex = 0.9, col = INK2)
      r <- k[k$KO == o[i], ]; text(i, yr[2], sprintf("p=%s\n%s", formatC(r$perm_p, format = "g", digits = 2), stars(r$perm_fdr)), cex = 0.6, adj = c(0.5, 1), font = ifelse(r$perm_fdr < 0.05, 2, 1)) }
    text(1:8, par("usr")[3] - 0.04 * diff(par("usr")[3:4]), o, srt = 45, adj = 1, xpd = NA, cex = 0.78, font = 2, col = KO_COL[o])
    g <- gl76[gl76$cov == cv & gl76$space == "signal500", ]; ttl(sprintf("%s  %s", if (cv == "none") "A" else "B", COVL[[cv]]))
    mtext(sprintf("global: same-KO %.3f vs other-KO %.3f, permutation p %s (family-wise %s)", g$diag_mean, g$off_mean, formatC(g$perm_p, format = "g", digits = 2), formatC(g$familywise_p, format = "g", digits = 2)), side = 3, adj = 0, line = -0.4, cex = 0.52, col = INK2) }
  par(mar = c(4.6, 4.2, 2.6, 0.6)); v <- list(G$none$rho[G$none$a == G$none$b], G$none$rho[G$none$a != G$none$b], G$cc$rho[G$cc$a == G$cc$b], G$cc$rho[G$cc$a != G$cc$b])
  boxplot(v, at = c(1, 1.7, 3, 3.7), col = c("#eb6834", "#d5d4cf", "#eb6834", "#d5d4cf"), border = INK2, outline = FALSE, xaxt = "n", las = 1, ylab = "Spearman rho between two cell lines", cex.axis = 0.75, boxwex = 0.55); abline(h = 0, lty = 2, col = INK)
  axis(1, at = c(1, 1.7, 3, 3.7), labels = rep(c("same\nKO", "other\nKO"), 2), tick = FALSE, cex.axis = 0.65, padj = 0.5); mtext(c("no cov.", "cc regr."), side = 1, at = c(1.35, 3.35), line = 2.6, cex = 0.6, font = 2)
  ttl("C  all KOs")
  mtext("Do bulk KO profiles reproduce across cell lines? Same KO in two lines (colour, n = 6 line pairs) vs different KOs (grey) -- S70, top-500 KO-signal genes", outer = TRUE, side = 3, adj = 0, line = 1.6, font = 2, cex = 1)
  mtext("TMM logCPM, limma moderated t of NTC - KO, profiles centred over the 8 KOs of each line, Spearman; genes chosen from the 2 lines not being compared; p = KO-label permutation (20,000), stars = BH FDR over the 8 KOs", outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2)
})
