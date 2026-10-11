# Spec 78: side-by-side comparison of the RUN3 sample set (S70 minus RUN3_EXCLUDE, 60 samples) against S70 (primary, 70 samples), from the saved outputs of spec73-76.
#   A global cross-line agreement (permutation z) per feature space, no covariate and cell cycle regressed (spec76 primary configuration)
#   B per-TF cross-line agreement (permutation z, top-500 KO-signal genes, no covariate), paired S70 -> RUN3
#   C whole benchmark: global permutation z of every spec73 configuration (centred / z-scored profiles, t statistics)
#   D same KO in another line vs replicate ceiling (spec74 per-sample profiles, top-500 signal genes)
#   E-F CROP-seq vs bulk (4 lines combined): screen effect vs bulk z over the 160 KO x Cm pairs, preserved pairs outlined (spec75, cell cycle adjusted)
# needs: spec73 (S70 and RUN3), spec74 (S70 and SET=RUN3), spec75 (SET=S70/RUN3), spec76 (SET=S70/RUN3)
# outputs: fig8E_RUN3_vs_S70.{png,pdf}, spec78_RUN3_vs_S70_summary.csv
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
SETS <- c("S70", "RUN3"); SCOL <- c(S70 = "#2a78d6", RUN3 = "#eb6834"); rd <- function(f) read.csv(file.path(OUT, f), stringsAsFactors = FALSE)
suf <- function(s, base) if (s == "S70") base else sub("\\.csv$", "_RUN3.csv", base)
glob <- lapply(setNames(SETS, SETS), function(s) rd(sprintf("spec76_primary_global_%s.csv", s)))
perko <- lapply(setNames(SETS, SETS), function(s) rd(sprintf("spec76_primary_per_KO_%s.csv", s)))
g73 <- list(S70 = subset(rd("spec73_grid.csv"), set == "S70"), RUN3 = rd("spec73_grid_RUN3.csv")); g73 <- lapply(g73, function(d) d[d$prep != "none" & d$stat %in% c("t", "t_rest"), ])
rel <- lapply(setNames(SETS, SETS), function(s) { r <- rd(suf(s, "spec74_replicate_reliability.csv")); r[r$space == "signal500", ] })
cm <- lapply(setNames(SETS, SETS), function(s) rd(sprintf("spec75_Cm_KO_table_%s.csv", s)))
SPL <- c(signal500 = "KO-signal genes", screen = "CROP-seq genes", Cm20_genes = "Cm modules", Hallmark = "Hallmark")

summ <- do.call(rbind, lapply(SETS, function(s) { x <- cm[[s]]; sig <- x$screen_fdr < 0.05; ct <- cor.test(x$screen_effect, x$bulk_z, method = "spearman", exact = FALSE); r <- rel[[s]]; g <- glob[[s]]
  data.frame(set = s, n_samples = c(S70 = 70, RUN3 = 60)[[s]], global_p_signal_none = g$perm_p[g$cov == "none" & g$space == "signal500"], global_p_signal_cc = g$perm_p[g$cov == "cc" & g$space == "signal500"],
             frac_configs_p05 = mean(g73[[s]]$perm_p < 0.05), n_KO_perm05_signal_none = sum(perko[[s]]$perm_p[perko[[s]]$cov == "none" & perko[[s]]$space == "signal500"] < 0.05),
             same_KO_other_line_rho = median(r$rho[r$class == "same KO, other line"]), replicate_rho = median(r$rho[r$class == "replicates (same KO, same line)"]),
             screen_vs_bulk_rho = unname(ct$estimate), screen_vs_bulk_p = ct$p.value, sign_agree_screen_sig = mean(sign(x$screen_effect[sig]) == sign(x$bulk_z[sig])), n_preserved = sum(x$preserved), n_preserved_specific = sum(x$preserved_specific)) }))
write.csv(summ, file.path(OUT, "spec78_RUN3_vs_S70_summary.csv"), row.names = FALSE); options(width = 220); print(format(summ, digits = 3), row.names = FALSE)

save_fig("fig8E_RUN3_vs_S70", 15, 8.8, function() {
  layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE), widths = c(1.25, 1.1, 0.8)); par(oma = c(0, 0, 3.4, 0.4), family = "sans")
  ttl <- function(x) mtext(x, side = 3, adj = 0, line = 0.6, font = 2, cex = 0.8); z05 <- qnorm(0.95)
  ## A global agreement
  par(mar = c(6.6, 4.2, 2.6, 0.6)); ids <- expand.grid(space = names(SPL), cov = c("none", "cc"), stringsAsFactors = FALSE)
  Z <- sapply(SETS, function(s) mapply(function(sp, cv) glob[[s]]$perm_z[glob[[s]]$space == sp & glob[[s]]$cov == cv], ids$space, ids$cov))
  P <- sapply(SETS, function(s) mapply(function(sp, cv) glob[[s]]$perm_p[glob[[s]]$space == sp & glob[[s]]$cov == cv], ids$space, ids$cov))
  bp <- barplot(t(Z), beside = TRUE, col = SCOL, border = NA, las = 1, ylab = "global agreement (permutation z)", names.arg = rep("", nrow(ids)), ylim = range(0, Z) * c(1.1, 1.18), cex.axis = 0.75)
  abline(h = z05, lty = 3, col = MUTED); text(bp, t(Z), ifelse(t(P) < 0.001, "***", ifelse(t(P) < 0.01, "**", ifelse(t(P) < 0.05, "*", ""))), pos = 3, cex = 0.8, xpd = NA)
  text(colMeans(bp), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), paste0(SPL[ids$space], "\n", ifelse(ids$cov == "none", "no covariate", "cell cycle regr.")), srt = 45, adj = 1, xpd = NA, cex = 0.62)
  legend("topright", c("S70 (70 samples)", "RUN3 (60 samples)"), fill = SCOL, border = NA, bty = "n", cex = 0.72); ttl("A  cross-line agreement, primary configuration")
  ## B per TF
  par(mar = c(4.6, 4.2, 2.6, 0.6)); pk <- lapply(SETS, function(s) { d <- perko[[s]]; d <- d[d$cov == "none" & d$space == "signal500", ]; setNames(d$perm_z, d$KO)[ORDER] }); names(pk) <- SETS
  o <- order(-pk$S70); plot(NA, xlim = c(0.5, 8.5), ylim = range(c(unlist(pk), z05)) + c(-0.3, 0.3), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "same-KO agreement (permutation z)", cex.axis = 0.75)
  abline(h = c(0, z05), lty = c(2, 3), col = c(INK, MUTED)); segments(seq_len(8), pk$S70[o], seq_len(8), pk$RUN3[o], col = "grey60", lwd = 1.5)
  points(seq_len(8), pk$S70[o], pch = 16, col = SCOL[["S70"]], cex = 1.6); points(seq_len(8), pk$RUN3[o], pch = 16, col = SCOL[["RUN3"]], cex = 1.6)
  np <- perko$RUN3; np <- setNames(np$n_pairs[np$cov == "none" & np$space == "signal500"], np$KO[np$cov == "none" & np$space == "signal500"])[ORDER][o]
  text(seq_len(8), par("usr")[3] + 0.02 * diff(par("usr")[3:4]), ifelse(np < 6, sprintf("%d/6 pairs", np), ""), cex = 0.55, col = SCOL[["RUN3"]], adj = c(0.5, 0))
  text(seq_len(8), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), ORDER[o], srt = 45, adj = 1, xpd = NA, cex = 0.75, font = 2, col = KO_COL[ORDER[o]]); ttl("B  per TF (KO-signal genes, no covariate)")
  ## C whole benchmark
  par(mar = c(4.6, 4.2, 2.6, 0.6)); boxplot(lapply(SETS, function(s) g73[[s]]$perm_z), names = sprintf("%s\n%.0f%% p<0.05", SETS, 100 * sapply(SETS, function(s) mean(g73[[s]]$perm_p < 0.05))), col = adjustcolor(SCOL, 0.5), border = INK2,
    outline = FALSE, las = 1, ylab = "global permutation z (every configuration)", cex.axis = 0.75); abline(h = c(0, z05), lty = c(2, 3), col = c(INK, MUTED)); ttl(sprintf("C  whole benchmark (%d configs each)", nrow(g73$S70)))
  ## D reliability
  par(mar = c(5.6, 4.2, 2.6, 0.6)); cls <- c("replicates (same KO, same line)", "same KO, other line", "different KO, other line"); lab <- c("replicates\n(same line)", "same KO\n(other line)", "other KO\n(other line)")
  v <- unlist(lapply(cls, function(cl) lapply(SETS, function(s) rel[[s]]$rho[rel[[s]]$class == cl])), recursive = FALSE); at <- c(1, 1.7, 3, 3.7, 5, 5.7)
  boxplot(v, at = at, col = adjustcolor(rep(SCOL, 3), 0.5), border = INK2, outline = FALSE, xaxt = "n", las = 1, ylab = "Spearman rho (per-sample profiles)", cex.axis = 0.75, boxwex = 0.55); abline(h = 0, lty = 2, col = INK)
  text(c(1.35, 3.35, 5.35), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), lab, xpd = NA, cex = 0.68, adj = c(0.5, 1))
  md <- sapply(v, median); text(at, md, sprintf("%.2f", md), pos = 3, cex = 0.55, col = INK); ttl("D  replicate ceiling vs cross-line")
  ## E-F CROP-seq vs bulk
  for (s in SETS) { x <- cm[[s]]; par(mar = c(4.6, 4.2, 2.6, 0.6)); sig <- x$screen_fdr < 0.05; ct <- cor.test(x$screen_effect, x$bulk_z, method = "spearman", exact = FALSE)
    plot(x$screen_effect, x$bulk_z, pch = 16, cex = ifelse(sig, 1.1, 0.7), col = adjustcolor(ifelse(sig, KO_COL[x$KO], "grey70"), 0.75), bty = "n", las = 1, xlab = sprintf("CROP-seq effect (%s, module score)", "NTC - KO"), ylab = "bulk combined z (4 lines)", cex.axis = 0.75, cex.lab = 0.8)
    abline(h = 0, v = 0, col = GRID); pr <- x$preserved; points(x$screen_effect[pr], x$bulk_z[pr], pch = 1, cex = 2, lwd = 1.6, col = INK); if (any(pr)) text(x$screen_effect[pr], x$bulk_z[pr], paste(x$KO[pr], x$Cm[pr]), pos = 4, cex = 0.5, col = INK)
    ttl(sprintf("%s  CROP-seq vs bulk, %s", if (s == "S70") "E" else "F", s)); mtext(sprintf("rho %.2f (p %.2g); same sign %.0f%% of %d screen-sig.; %d preserved (circled)", ct$estimate, ct$p.value, 100 * mean(sign(x$screen_effect[sig]) == sign(x$bulk_z[sig])), sum(sig), sum(pr)), side = 3, adj = 0, line = -0.4, cex = 0.55, col = SCOL[[s]], font = 2) }
  mtext("RUN3 (60 samples) vs S70 (70 samples): removing the flagged samples lowers agreement across cell lines and with CROP-seq", outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 1.0)
  mtext("TMM logCPM, limma moderated t of NTC - KO, KO-centred profiles, Spearman; p / z from KO-label permutation within lines; * p<0.05 ** p<0.01 *** p<0.001; dotted = one-sided 5 %; RUN3 lacks BPTF and HMGA1 in MSTO-211H and TEAD4 in NCI-H2452", outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.58, col = INK2)
})
