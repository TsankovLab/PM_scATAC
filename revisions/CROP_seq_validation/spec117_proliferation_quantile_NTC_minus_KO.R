# Spec 117: proliferation score, NTC - KO in each bulk cell line, with QUANTILE NORMALISATION applied before scoring.
# Proliferation score = Seurat CellCycleScoring S.Score + G2M.Score (cc.genes), exactly as the "cc" of the spec71 audit used in spec107/108. The only change here:
#   current   Seurat NormalizeData (log1p CP10K) -> CellCycleScoring                                                   [= spec71 / spec107]
#   quantile  the same log-normalised matrix, limma::normalizeQuantiles across the 70 S70 samples, -> CellCycleScoring
#   (tmmq, TMM logCPM per line + quantile as in spec103, is scored too and reported as a cross-check of the input scale)
# NTC - KO per cell line = mean score of the line's NTC samples minus the mean score of the KO's replicates; > 0 = less proliferation after the KO.
# Replicates are averaged first, so each line contributes one number per KO. NCI-H2052 has a single NTC sample: its whole column rests on that one sample.
# Per KO over the 4 lines: number of lines with a drop and a one-sample t-test on the 4 line values (n = 4, so this is weak by construction).
# CROP-seq (NCI-H2052 screen): the same score per cell, mean NTC - mean KO (spec75), for comparison; Spearman over the 8 KOs against each line.
# needs: spec71 audit (cc), spec72 (S70 sample list), spec75 (screen scores)
# fig32: the scores themselves with NTC as a group, z-scored within each cell line (bulk) / over the cells (screen); KO medians against the NTC median.
# outputs: spec117_proliferation_scaled_medians.csv, fig32_proliferation_scaled_with_NTC.{png,pdf}, spec117_proliferation_scores_per_sample.csv, spec117_proliferation_NTC_minus_KO_per_line.csv, spec117_proliferation_summary.csv, fig31_proliferation_quantile_NTC_minus_KO.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3", width = 220); SH <- function(l) sub("NCI-|-211H", "", l)
fp <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2)))
fab <- load_full_audited_bulk(); meta <- fab$meta; cs <- fab$cs; q <- load_audit(); pr <- readRDS(file.path(OUT, "spec72_profiles_S70.rds")); ids <- pr$samples; line <- meta[ids, "Cell_Line"]; grp <- meta[ids, "Gene_Control"]
det <- rownames(cs)[rowSums(cs[, ids]) > 0]; qn <- function(y) { z <- limma::normalizeQuantiles(y); dimnames(z) <- dimnames(y); z }
## Seurat replaces "_" by "-" in feature names, so counts and data are given the same sanitised names (rows are in the same order: all matrices derive from cs[det, ids])
cnt <- cs[det, ids]; rownames(cnt) <- gsub("_", "-", det)
ccs <- function(y) { y <- y - min(y); rownames(y) <- rownames(cnt); o <- CreateSeuratObject(counts = cnt); o <- SetAssayData(o, slot = "data", new.data = as.matrix(y))
  o <- suppressWarnings(CellCycleScoring(o, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes, set.ident = FALSE)); setNames(as.numeric(o$S.Score + o$G2M.Score), ids) }
o0 <- NormalizeData(CreateSeuratObject(counts = cnt), verbose = FALSE); yL <- as.matrix(GetAssayData(o0, slot = "data"))[rownames(cnt), ids]                       # log1p CP10K, as spec71
yT <- do.call(cbind, lapply(LINES, function(cl) normalise(cs[det, ids[line == cl]], "tmm")$y))[, ids]
SC <- data.frame(sample = ids, line = line, group = grp, audit_cc = q[ids, "cc"], current = ccs(yL), quantile = ccs(qn(yL)), tmmq = ccs(qn(yT)))
write.csv(SC, file.path(OUT, "spec117_proliferation_scores_per_sample.csv"), row.names = FALSE)
cat("check: 'current' recomputed on the 70 samples vs the audit cc (scored on all 76): r =", round(cor(SC$current, SC$audit_cc), 4), "\n")
cat("score before vs after quantile normalisation, over the 70 samples: r =", round(cor(SC$current, SC$quantile), 4), "; within line:", paste(sprintf("%s %.3f", SH(LINES), sapply(LINES, function(cl) cor(SC$current[line == cl], SC$quantile[line == cl]))), collapse = "  "), "\n")
cat("quantile (log1p CP10K input) vs tmmq (TMM logCPM input): r =", round(cor(SC$quantile, SC$tmmq), 4), "\n")

dlt <- function(col) sapply(LINES, function(cl) sapply(KOS, function(k) mean(SC[[col]][line == cl & grp == "NTC"]) - mean(SC[[col]][line == cl & grp == k])))
DQ <- dlt("quantile"); DC <- dlt("current"); DT <- dlt("tmmq")
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); gs <- as.character(sc$group); scr <- sapply(KOS, function(k) mean(sc$cc[gs == "NTC"]) - mean(sc$cc[gs == k])); scr_p <- sapply(KOS, function(k) wilcox.test(sc$cc[gs == "NTC"], sc$cc[gs == k])$p.value)
L <- do.call(rbind, lapply(LINES, function(cl) data.frame(line = cl, KO = KOS, n_NTC = sum(line == cl & grp == "NTC"), n_KO = sapply(KOS, function(k) sum(line == cl & grp == k)),
  NTC_minus_KO_quantile = DQ[, cl], NTC_minus_KO_current = DC[, cl], NTC_minus_KO_tmmq = DT[, cl]))); write.csv(L, file.path(OUT, "spec117_proliferation_NTC_minus_KO_per_line.csv"), row.names = FALSE)
SM <- data.frame(KO = KOS, round(DQ, 3), check.names = FALSE); colnames(SM)[-1] <- SH(LINES); SM$mean_4_lines <- rowMeans(DQ); SM$n_lines_drop <- rowSums(DQ > 0); SM$p_t <- apply(DQ, 1, function(v) t.test(v)$p.value); SM$fdr <- p.adjust(SM$p_t, "BH")
SM$mean_current <- rowMeans(DC); SM$n_lines_drop_current <- rowSums(DC > 0); SM$screen_NTC_minus_KO <- scr; SM$screen_p <- scr_p; SM <- SM[order(-SM$mean_4_lines), ]
write.csv(SM, file.path(OUT, "spec117_proliferation_summary.csv"), row.names = FALSE)
cat("\n=== proliferation score, NTC - KO per cell line, QUANTILE-normalised (> 0 = less proliferation after the KO) ===\n"); print(format(SM, digits = 2), row.names = FALSE)
cat("NTC samples per line:", paste(SH(LINES), sapply(LINES, function(cl) sum(line == cl & grp == "NTC")), collapse = ", "), "\n")
cat("\nsame, without quantile normalisation (current):\n"); print(round(DC[SM$KO, ], 3))
cat("\nchange caused by the quantile normalisation (quantile - current):\n"); print(round((DQ - DC)[SM$KO, ], 3)); cat("largest absolute change:", round(max(abs(DQ - DC)), 3), "; sign changes:", sum(sign(DQ) != sign(DC)), "of 32\n")
cat("\nall KOs vs NTC, per line (mean NTC - KO over the 8 KOs):", paste(sprintf("%s %+.3f", SH(LINES), colMeans(DQ)), collapse = "  "), "; screen", sprintf("%+.3f", mean(scr)), "\n")
ag <- sapply(LINES, function(cl) cor(DQ[, cl], scr[KOS], method = "spearman")); xl <- combn(LINES, 2, function(p) cor(DQ[, p[1]], DQ[, p[2]], method = "spearman"))
cat("KO ranking, Spearman over the 8 KOs -- screen vs each line:", paste(sprintf("%s %+.2f", SH(LINES), ag), collapse = "  "), "; screen vs 4-line mean", sprintf("%+.2f", cor(rowMeans(DQ), scr[KOS], method = "spearman")), "\n")
cat("bulk line vs bulk line:", paste(sprintf("%s %+.2f", combn(SH(LINES), 2, paste, collapse = "-"), xl), collapse = "  "), "\n")

## KOs ordered by the CROP-seq effect (largest drop in proliferation first), so every panel reads against the screen
SM <- SM[order(-SM$screen_NTC_minus_KO), ]; write.csv(SM, file.path(OUT, "spec117_proliferation_summary.csv"), row.names = FALSE)
o <- SM$KO; PCH <- setNames(c(21, 22, 24, 23), LINES); cat("\nKO order (CROP-seq NTC - KO, decreasing):", paste(o, collapse = " > "), "\n")
## CROP-seq drawn the same way as the bulk panel: NTC - KO per unit of observation (here the KO cell: mean score of the NTC cells minus the cell's score), one box per KO
ntcS <- mean(sc$cc[gs == "NTC"]); CELL <- lapply(setNames(KOS, KOS), function(k) ntcS - sc$cc[gs == k])
save_fig("fig31_proliferation_quantile_NTC_minus_KO", 16, 9.2, function() {
  layout(matrix(c(1, 1, 1, 2, 2, 2, 3, 4, 5, 6, 7, 7), nrow = 2, byrow = TRUE), heights = c(1.15, 1)); par(family = "sans", oma = c(0, 0.5, 3.4, 0.4)); set.seed(1)
  par(mar = c(4.2, 4.8, 3.6, 1)); st <- sapply(CELL[o], function(v) boxplot.stats(v)$stats); yr <- range(st) + c(-0.08, 0.3) * diff(range(st))
  plot(NA, xlim = c(0.5, 8.5), ylim = yr, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "proliferation score, NTC - KO (per cell)", cex.axis = 0.85); abline(h = 0, lty = 2, col = INK)
  for (i in seq_along(o)) { v <- CELL[[o[i]]]; pv <- v[v >= yr[1] & v <= yr[2]]; pv <- if (length(pv) > 150) sample(pv, 150) else pv
    points(jitter(rep(i, length(pv)), amount = 0.22), pv, pch = 16, cex = 0.45, col = adjustcolor(KO_COL[[o[i]]], 0.35))
    boxplot(v, at = i, add = TRUE, boxwex = 0.6, col = adjustcolor(KO_COL[[o[i]]], 0.6), border = INK2, outline = FALSE, axes = FALSE); points(i, mean(v), pch = 23, bg = "white", col = INK, cex = 1.4)
    text(i, yr[2], sprintf("%d cells\np %s", length(v), fp(scr_p[o[i]])), cex = 0.7, adj = c(0.5, 1), font = ifelse(scr_p[o[i]] < 0.05, 2, 1)) }
  text(seq_along(o), yr[1] - 0.06 * diff(yr), o, xpd = NA, font = 2, cex = 0.95, col = KO_COL[o])
  mtext("A  CROP-seq screen (NCI-H2052), one point per KO cell", side = 3, adj = 0, line = 1.6, font = 2, cex = 0.85); mtext("mean of the NTC cells minus the cell; diamond = mean; Wilcoxon NTC vs KO cells; up to 150 cells drawn per KO; KOs ordered by this mean", side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2)
  yr <- range(DQ) + c(-0.08, 0.3) * diff(range(DQ)); plot(NA, xlim = c(0.5, 8.5), ylim = yr, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "proliferation score, NTC - KO (per cell line)", cex.axis = 0.85); abline(h = 0, lty = 2, col = INK)
  for (i in seq_along(o)) { v <- DQ[o[i], ]; boxplot(v, at = i, add = TRUE, boxwex = 0.6, col = adjustcolor(KO_COL[[o[i]]], 0.6), border = INK2, outline = FALSE, axes = FALSE)
    points(i + seq(-0.2, 0.2, length.out = 4), v, pch = PCH[LINES], bg = LINE_COL[LINES], col = INK2, cex = 1.5); points(i, mean(v), pch = 23, bg = "white", col = INK, cex = 1.4)
    text(i, yr[2], sprintf("%d/4 lines\np %s", SM$n_lines_drop[i], fp(SM$p_t[i])), cex = 0.7, adj = c(0.5, 1), font = ifelse(SM$p_t[i] < 0.05, 2, 1)) }
  text(seq_along(o), yr[1] - 0.06 * diff(yr), o, xpd = NA, font = 2, cex = 0.95, col = KO_COL[o]); legend("bottomleft", SH(LINES), pch = PCH, pt.bg = LINE_COL[LINES], col = INK2, bty = "n", cex = 0.75, horiz = TRUE)
  mtext("B  bulk RNA-seq, 4 cell lines, one point per line (replicates averaged), same KO order", side = 3, adj = 0, line = 1.6, font = 2, cex = 0.85); mtext("mean of the line's NTC samples minus the KO; diamond = mean; lines with a drop and one-sample t-test over the 4 line values", side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2)
  for (cl in LINES) { par(mar = c(5, 4.4, 3.2, 0.6)); v <- DQ[o, cl]; bp <- barplot(v, col = adjustcolor(KO_COL[o], 0.85), border = NA, las = 2, names.arg = o, cex.names = 0.7, cex.axis = 0.75, ylim = range(c(0, DQ)) * 1.15, ylab = "NTC - KO")
    abline(h = 0, col = INK); points(bp, DC[o, cl], pch = 4, cex = 0.8, col = INK2); mtext(sprintf("%s  (%d NTC)", SH(cl), sum(line == cl & grp == "NTC")), side = 3, adj = 0, line = 1.3, font = 2, cex = 0.78, col = LINE_COL[[cl]])
    mtext("x = without quantile norm.", side = 3, adj = 0, line = 0.3, cex = 0.5, col = INK2) }
  par(mar = c(5, 4.8, 3.2, 1)); plot(scr[KOS], rowMeans(DQ), pch = 21, bg = KO_COL[KOS], col = INK2, cex = 2, bty = "n", las = 1, xlab = "CROP-seq: NTC - KO", ylab = "bulk: NTC - KO, mean of 4 lines", cex.axis = 0.8); abline(h = 0, v = 0, lty = 2, col = GRID)
  text(scr[KOS], rowMeans(DQ), KOS, pos = 3, cex = 0.7, font = 2, col = KO_COL[KOS], xpd = NA); mtext("C  screen vs bulk", side = 3, adj = 0, line = 1.3, font = 2, cex = 0.82)
  mtext(sprintf("Spearman %.2f over the 8 KOs; per line: %s", cor(rowMeans(DQ), scr[KOS], method = "spearman"), paste(sprintf("%s %.2f", SH(LINES), ag), collapse = ", ")), side = 3, adj = 0, line = 0.3, cex = 0.55, col = INK2)
  mtext("Proliferation after each KO, CROP-seq and bulk side by side: NTC - KO of the cell-cycle score", outer = TRUE, side = 3, adj = 0.01, line = 1.7, font = 2, cex = 1)
  mtext("proliferation = S.Score + G2M.Score (Seurat CellCycleScoring); > 0 = less proliferation after the KO; bulk: log-normalised counts, quantile normalisation across the 70 samples, then scored", outer = TRUE, side = 3, adj = 0.01, line = 0.5, cex = 0.6, col = INK2)
})
cat("\nCROP-seq, NTC - KO per cell: median / mean / fraction of cells below the NTC mean\n"); print(round(t(sapply(CELL[o], function(v) c(n = length(v), median = median(v), mean = mean(v), frac_lower = mean(v > 0)))), 3))

## ---------------------------------------------------------------- fig32: the scores themselves, NTC shown as a group, scaled within each data set
## Bulk: the quantile-normalised score is z-scored WITHIN each cell line across that line's NTC + KO samples (removes the line's baseline and spread), then the 4 lines are
## pooled, one point per sample. CROP-seq: z-scored across the cells of the same 9 groups. Same group order in both (NTC, then the KOs by the screen effect).
## Reference line = the NTC median; every KO is annotated with its median minus the NTC median. The 8 comparisons share one NTC group, so "all 8 below" is a description
## of a shared shift (any KO vs NTC), not 8 independent tests.
GRO <- c("NTC", o); SC$z <- ave(SC$quantile, SC$line, FUN = function(v) (v - mean(v)) / sd(v)); keepS <- gs %in% GRO; zS <- as.numeric(scale(sc$cc[keepS])); gK <- gs[keepS]
BZ <- lapply(setNames(GRO, GRO), function(g) SC$z[SC$group == g]); BL <- lapply(setNames(GRO, GRO), function(g) SC$line[SC$group == g]); CZ <- lapply(setNames(GRO, GRO), function(g) zS[gK == g])
tab <- function(Z) data.frame(group = GRO, n = lengths(Z), median = sapply(Z, median), median_minus_NTC = sapply(Z, median) - median(Z$NTC), p_vs_NTC = sapply(GRO, function(g) if (g == "NTC") NA else wilcox.test(Z[[g]], Z$NTC)$p.value))
TB <- tab(BZ); TC <- tab(CZ); write.csv(cbind(data = rep(c("bulk", "CROP-seq"), each = length(GRO)), rbind(TB, TC)), file.path(OUT, "spec117_proliferation_scaled_medians.csv"), row.names = FALSE)
cat("\n=== scaled proliferation score, medians per group ===\nbulk (z within line, 4 lines pooled):\n"); print(format(TB, digits = 3), row.names = FALSE); cat("CROP-seq (z over the cells of the 9 groups):\n"); print(format(TC, digits = 3), row.names = FALSE)
## the screen also carries TF KOs that were not taken to bulk (MEF2A, MEF2D): panel A shows every screen KO, z-scored over the cells of all its groups; the 8 shared KOs keep their relative order
allS <- setdiff(sort(unique(gs)), "NTC"); oS <- allS[order(-sapply(allS, function(k) mean(sc$cc[gs == "NTC"]) - mean(sc$cc[gs == k])))]; GROS <- c("NTC", oS)
zA <- as.numeric(scale(sc$cc)); CZ <- lapply(setNames(GROS, GROS), function(g) zA[gs == g])
TC <- data.frame(group = GROS, n = lengths(CZ), median = sapply(CZ, median), median_minus_NTC = sapply(CZ, median) - median(CZ$NTC), p_vs_NTC = sapply(GROS, function(g) if (g == "NTC") NA else wilcox.test(CZ[[g]], CZ$NTC)$p.value))
TC$in_bulk <- TC$group %in% c("NTC", KOS); cat("\nCROP-seq, ALL TF KOs of the screen (z over all", length(zA), "cells):\n"); print(format(TC, digits = 3), row.names = FALSE)
write.csv(TC, file.path(OUT, "spec117_proliferation_scaled_medians_screen_all_KOs.csv"), row.names = FALSE)
nb <- sum(TB$median_minus_NTC[-1] < 0); nc <- sum(TC$median_minus_NTC[-1] < 0)
cat(sprintf("KOs with a median below NTC: bulk %d of 8, CROP-seq %d of %d\n", nb, nc, length(oS)))
pl <- sapply(LINES, function(cl) sapply(o, function(k) median(SC$z[SC$line == cl & SC$group == k]) - median(SC$z[SC$line == cl & SC$group == "NTC"])))
cat("per line, KO median minus NTC median (z):\n"); print(round(pl, 2)); cat("KOs below NTC per line:", paste(sprintf("%s %d/8", SH(LINES), colSums(pl < 0)), collapse = "  "), "\n")
allB <- unlist(BZ[o]); allC <- unlist(CZ[oS]); cat(sprintf("all KO samples vs NTC, bulk: median %.2f vs %.2f, Wilcoxon p = %s; CROP-seq: all KO cells %.2f vs %.2f, p = %s\n", median(allB), median(BZ$NTC), fp(wilcox.test(allB, BZ$NTC)$p.value), median(allC), median(CZ$NTC), fp(wilcox.test(allC, CZ$NTC)$p.value)))
## compact layout: A screen | B bulk (same axis style, NTC median as the reference, zone below it shaded) | C the shift of every KO in both technologies on one axis
save_fig("fig32_proliferation_scaled_with_NTC", 7.4, 2.9, res = 300, fun = function() {
  layout(matrix(1:3, nrow = 1), widths = c(1.18, 1, 0.85)); par(family = "sans", oma = c(0, 0, 1.2, 0.1), mgp = c(1.45, 0.35, 0), tcl = -0.2); set.seed(1); gc <- c(NTC = "#7a7a7a", KO_COL[oS]); BLUE <- "#2a78d6"; ORA <- "#eb6834"
  pan <- function(Z, T, ttl, pts) { GRO <- names(Z); nG <- length(GRO); par(mar = c(3.0, 2.6, 2.1, 0.2)); bs <- sapply(Z, function(v) boxplot.stats(v)$stats); dr <- if (is.null(pts)) range(unlist(Z)) else range(bs); h <- diff(dr); m0 <- T$median[1]
    sg <- which(T$p_vs_NTC < 0.05); step <- 0.085 * h; top <- dr[2] + 0.09 * h; yr <- c(dr[1] - 0.03 * h, top + max(length(sg), 1) * step + 0.03 * h)   # room for one bracket per significant KO
    plot(NA, xlim = c(0.45, nG + 0.55), ylim = yr, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "proliferation score (z)", cex.lab = 0.72); axis(2, at = pretty(dr), las = 1, cex.axis = 0.62)
    rect(1.5, yr[1], nG + 0.6, m0, col = adjustcolor(BLUE, 0.09), border = NA); segments(0.45, m0, nG + 0.6, m0, col = "#7a7a7a", lwd = 1, lty = 2)
    for (i in seq_along(GRO)) { g <- GRO[i]; v <- Z[[g]]; low <- T$median_minus_NTC[i] < 0
      if (is.null(pts)) points(jitter(rep(i, length(v)), amount = 0.2), v, pch = 16, cex = 0.45, col = adjustcolor(gc[[g]], 0.8)) else { pv <- v[v >= dr[1] & v <= dr[2]]; if (length(pv) > pts) pv <- sample(pv, pts); points(jitter(rep(i, length(pv)), amount = 0.28), pv, pch = 16, cex = 0.18, col = adjustcolor(gc[[g]], 0.28)) }
      boxplot(v, at = i, add = TRUE, boxwex = 0.86, col = adjustcolor(gc[[g]], 0.55), border = INK2, outline = FALSE, axes = FALSE, medlwd = 2.2, lwd = 0.6)
      if (i > 1) text(i, dr[2] + 0.045 * h, ifelse(low, "\u25BC", "\u25B2"), col = ifelse(low, BLUE, ORA), cex = 0.6) }
    for (k in seq_along(sg)) { y <- top + (k - 1) * step; x1 <- sg[k]; segments(c(1, 1, x1), c(y, y - 0.3 * step, y - 0.3 * step), c(x1, 1, x1), c(y, y, y), lwd = 0.6, col = INK)   # bracket NTC -> KO
      text((1 + x1) / 2, y + 0.08 * step, sprintf("p = %s", fp(T$p_vs_NTC[x1])), cex = 0.5, adj = c(0.5, 0), col = INK) }
    text(seq_along(GRO), yr[1] - 0.04 * diff(yr), GRO, xpd = NA, font = 2, cex = 0.55, col = gc[GRO], srt = 45, adj = c(1, 1))
    mtext(ttl, side = 3, adj = 0, line = 1.1, font = 2, cex = 0.55); mtext(sprintf("%d of %d KOs below the NTC median", sum(T$median_minus_NTC[-1] < 0), nG - 1), side = 3, adj = 0, line = 0.3, font = 2, cex = 0.52, col = BLUE) }
  pan(CZ, TC, "A  CROP-seq (NCI-H2052), all TF KOs, cells", 120); pan(BZ, TB, "B  bulk RNA-seq, 4 lines, samples", NULL)
  par(mar = c(3.0, 2.6, 2.1, 0.3)); nS <- length(oS); m <- rbind(TC$median_minus_NTC[-1], TB$median_minus_NTC[match(oS, TB$group)]); yr <- range(c(m, 0), na.rm = TRUE) + c(-0.06, 0.22) * diff(range(c(m, 0), na.rm = TRUE))
  plot(NA, xlim = c(0.4, nS + 0.6), ylim = yr, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "KO median - NTC median (z)", cex.axis = 0.62, cex.lab = 0.72); rect(0.4, yr[1], nS + 0.7, 0, col = adjustcolor(BLUE, 0.09), border = NA); segments(0.4, 0, nS + 0.7, 0, col = "#7a7a7a", lwd = 1, lty = 2)
  for (i in seq_len(nS)) { hb <- is.finite(m[2, i]); dx <- if (hb) 0.18 else 0; segments(i - dx, 0, i - dx, m[1, i], col = KO_COL[[oS[i]]], lwd = 1.1); points(i - dx, m[1, i], pch = 21, bg = KO_COL[[oS[i]]], col = INK, cex = 0.95)
    if (hb) { segments(i + dx, 0, i + dx, m[2, i], col = KO_COL[[oS[i]]], lwd = 1.1); points(i + dx, m[2, i], pch = 24, bg = "white", col = KO_COL[[oS[i]]], cex = 0.85, lwd = 1.1) } }
  text(seq_len(nS), yr[1] - 0.04 * diff(yr), oS, xpd = NA, font = 2, cex = 0.55, col = KO_COL[oS], srt = 45, adj = c(1, 1)); legend("topright", c("CROP-seq", "bulk"), pch = c(21, 24), pt.bg = c("grey40", "white"), col = c(INK, "grey40"), bty = "n", cex = 0.5, pt.cex = 0.75, horiz = TRUE, inset = c(0, -0.02))
  mtext("C  shift of each KO", side = 3, adj = 0, line = 1.1, font = 2, cex = 0.55); mtext(sprintf("%d of %d below NTC", sum(m < 0, na.rm = TRUE), sum(is.finite(m))), side = 3, adj = 0, line = 0.3, font = 2, cex = 0.52, col = BLUE)
  mtext("Proliferation score (S + G2M), z-scored within each data set (bulk: within cell line, quantile-normalised); dashed = NTC median; \u25BC below / \u25B2 above NTC; brackets = Wilcoxon vs NTC, p < 0.05 only", outer = TRUE, side = 3, adj = 0.01, line = 0.35, cex = 0.42, col = INK2)
})
sh <- TC$median_minus_NTC[match(o, TC$group)] < 0 & TB$median_minus_NTC[match(o, TB$group)] < 0; cat(sprintf("\nbelow NTC in BOTH technologies for %d of the 8 shared KOs (%s)\n", sum(sh), paste(o[sh], collapse = ", ")))
