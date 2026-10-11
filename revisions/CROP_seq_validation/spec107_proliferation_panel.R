# Spec 107: proliferation panel. Proliferation score of a bulk sample = Seurat CellCycleScoring S.Score + G2M.Score (cc.genes; the "cc" column of the spec71 audit).
# For every KO sample: mean score of the NTC samples of its cell line minus the score of the sample (NTC - KO; > 0 = less proliferation after the KO). Boxplots per KO, one point per KO sample (S70).
#   A  the 4 lines pooled, KOs ordered by the median; p = lm(NTC - KO ~ 1) per KO with the line as fixed effect dropped (one-sample t over the samples) and, below, the number of lines with a positive mean
#   B  one panel per cell line
# CROP-seq (NCI-H2052 screen): the same score per cell (spec75_screen_scores_by_cell.rds, cc), mean NTC - mean KO, shown as a diamond on its own axis in C with the Wilcoxon p
# needs: spec71 audit (cc), spec72 (S70 sample list), spec75 (screen scores)
# outputs: spec107_proliferation_NTC_minus_KO.csv, spec107_proliferation_summary.csv, fig22_proliferation_panel.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
fab <- load_full_audited_bulk(); meta <- fab$meta; q <- load_audit(); pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); ids <- pr$samples; SH <- function(l) sub("NCI-|-211H", "", l); fp <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2)))
d <- data.frame(sample = ids, line = meta[ids, "Cell_Line"], group = meta[ids, "Gene_Control"], cc = q[ids, "cc"]); ntc <- tapply(d$cc[d$group == "NTC"], d$line[d$group == "NTC"], mean); nn <- table(d$line[d$group == "NTC"])
D <- d[d$group %in% KOS, ]; D$NTC_mean <- as.numeric(ntc[D$line]); D$NTC_minus_KO <- D$NTC_mean - D$cc; write.csv(D, file.path(OUT, "spec107_proliferation_NTC_minus_KO.csv"), row.names = FALSE)
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); gs <- as.character(sc$group)
SM <- do.call(rbind, lapply(KOS, function(k) { v <- D$NTC_minus_KO[D$group == k]; pl <- tapply(D$NTC_minus_KO[D$group == k], D$line[D$group == k], mean)[LINES]
  data.frame(KO = k, n = length(v), median = median(v), mean = mean(v), p = t.test(v)$p.value, n_lines_pos = sum(pl > 0), t(setNames(round(pl, 3), SH(LINES))), screen_NTC_minus_KO = mean(sc$cc[gs == "NTC"]) - mean(sc$cc[gs == k]), screen_p = wilcox.test(sc$cc[gs == "NTC"], sc$cc[gs == k])$p.value) }))
SM$fdr <- p.adjust(SM$p, "BH"); SM <- SM[order(-SM$median), ]; write.csv(SM, file.path(OUT, "spec107_proliferation_summary.csv"), row.names = FALSE); options(width = 200); print(SM, digits = 2, row.names = FALSE)
cat("NTC samples per line:", paste(SH(names(nn)), nn, collapse = ", "), "\n"); o <- SM$KO; PCH <- setNames(c(21, 22, 24, 23), LINES)
save_fig("fig22_proliferation_panel", 15, 9.2, function() {
  layout(rbind(c(1, 1, 2, 2), 3:6), heights = c(1.15, 1)); par(oma = c(0, 0, 3, 0.4), family = "sans"); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.85); mtext(b, side = 3, adj = 0, line = 0.4, cex = 0.62, col = INK2) }
  bx <- function(dd, ord, main, sub, stats = NULL, cexp = 1.2) { yr <- range(dd$NTC_minus_KO) + c(-0.05, 0.22) * diff(range(dd$NTC_minus_KO)); plot(NA, xlim = c(0.5, 8.5), ylim = yr, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "proliferation score, NTC - KO", cex.axis = 0.8); abline(h = 0, lty = 2, col = INK)
    for (i in seq_along(ord)) { e <- dd[dd$group == ord[i], ]; if (nrow(e) >= 4) boxplot(e$NTC_minus_KO, at = i, add = TRUE, boxwex = 0.6, col = adjustcolor(KO_COL[[ord[i]]], 0.45), border = INK2, outline = FALSE, axes = FALSE) else if (nrow(e)) segments(i - 0.25, mean(e$NTC_minus_KO), i + 0.25, mean(e$NTC_minus_KO), lwd = 2, col = KO_COL[[ord[i]]])
      if (nrow(e)) points(i + seq(-0.2, 0.2, length.out = nrow(e))[order(order(e$line))], e$NTC_minus_KO, pch = PCH[e$line], bg = LINE_COL[e$line], col = INK2, cex = cexp)
      if (!is.null(stats)) { s <- stats[stats$KO == ord[i], ]; text(i, yr[2], sprintf("p=%s%s\n%d/4 lines", fp(s$p), ifelse(s$fdr < 0.05, " *", ""), s$n_lines_pos), cex = 0.62, adj = c(0.5, 1), font = ifelse(s$p < 0.05, 2, 1)) } }
    text(1:8, par("usr")[3] - 0.03 * diff(par("usr")[3:4]), ord, srt = 45, adj = 1, xpd = NA, cex = 0.85, font = 2, col = KO_COL[ord]); ttl(main, sub) }
  par(mar = c(4.6, 4.8, 3.4, 1)); bx(D, o, "A  bulk, 4 cell lines pooled", "one point per KO sample; > 0 = less proliferation after the KO; one-sample t-test, * FDR < 0.05; n/4 = lines with a positive mean", SM); legend("bottomleft", SH(LINES), pch = PCH, pt.bg = LINE_COL[LINES], col = INK2, bty = "n", cex = 0.7, horiz = TRUE)
  v <- SM$screen_NTC_minus_KO; yr <- range(c(0, v)) + c(-0.1, 0.25) * diff(range(c(0, v))); bp <- barplot(v, col = adjustcolor(KO_COL[o], 0.8), border = NA, las = 1, ylim = yr, names.arg = rep("", 8), ylab = "proliferation score, mean NTC - mean KO (cells)", cex.axis = 0.8); abline(h = 0, col = INK)
  text(bp, pmax(v, 0) + 0.02 * diff(yr), sprintf("p=%s", fp(SM$screen_p)), cex = 0.62, adj = c(0.5, 0), font = ifelse(SM$screen_p < 0.05, 2, 1)); text(bp, par("usr")[3] - 0.03 * diff(par("usr")[3:4]), o, srt = 45, adj = 1, xpd = NA, cex = 0.85, font = 2, col = KO_COL[o])
  ttl("C  CROP-seq (NCI-H2052 screen), same KO order", "S.Score + G2M.Score per cell; Wilcoxon NTC vs KO cells")
  for (cl in LINES) { par(mar = c(4.4, 4.4, 3.2, 0.6)); bx(D[D$line == cl, ], o, sprintf("B  %s", SH(cl)), sprintf("%d NTC sample%s; bar = mean of the replicates", nn[[cl]], ifelse(nn[[cl]] > 1, "s", "")), cexp = 1.3) }
  mtext("Proliferation after each KO: NTC - KO of the cell-cycle score (S.Score + G2M.Score), relative to the NTC of the same cell line", outer = TRUE, side = 3, adj = 0, line = 1.5, font = 2, cex = 1)
  mtext("bulk RNA-seq, S70; each KO sample subtracted from the mean of the NTC samples of its line; KOs ordered by the pooled median", outer = TRUE, side = 3, adj = 0, line = 0.3, cex = 0.65, col = INK2)
})
