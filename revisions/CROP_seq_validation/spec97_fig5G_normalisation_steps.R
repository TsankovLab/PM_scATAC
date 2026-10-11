# Spec 97: step-by-step illustration of how fig5G_KO_x_target_per_line_rings (spec69) is normalised, on real data of one line (env LINE, default NCI-H2052).
# The steps are computed with the SAME expressions as spec69 (lines 8, 46-49, 99) and the result is checked against the values spec69 saved
# (spec69_KO_x_target_per_line_rowz.csv); the script stops if they differ.
#   1 per-sample log2(CP10K + 1) of the 8 target genes          2 median per group (NTC + 8 KOs)                    3 minus the gene's mean over the 9 group medians
#   4 divided by the gene's SD over the 9 group medians (= z)   5 transposed, each group row z-scored over 8 targets   6 the rings: one ring per line, coloured by step 5
# needs: repro_common.R (load_bulk()), spec69 outputs
# outputs: spec97_fig5G_steps_<line>.csv, fig5H_fig5G_normalisation_steps.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
CLX <- Sys.getenv("LINE", unset = "NCI-H2052"); CROP_PAL <- rev(as.character(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 100)))
## ---- identical to spec69 ----
KOS <- ORDER; b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1); GR <- c("NTC", KOS)
Mfull <- lapply(setNames(LINES, LINES), function(cl) { mm <- meta[meta$Cell_Line == cl, ]
  t(sapply(KOS, function(gene) sapply(GR, function(g) { s <- rownames(mm)[mm$Gene_Control == g]; if (!length(s) || !gene %in% rownames(lg)) NA else median(lg[gene, s]) }))) })
Zfull <- lapply(Mfull, function(M) { Z <- (M - rowMeans(M, na.rm = TRUE)) / matrixStats::rowSds(M, na.rm = TRUE); Z[!is.finite(Z)] <- NA; Z })
Kline <- lapply(Zfull, function(Z) { K <- t(Z[KOS, GR]); Kz <- (K - rowMeans(K, na.rm = TRUE)) / matrixStats::rowSds(K, na.rm = TRUE); Kz[!is.finite(Kz)] <- NA; Kz })
## ---- check against what spec69 saved ----
sv <- read.csv(file.path(OUT, "spec69_KO_x_target_per_line_rowz.csv"), check.names = FALSE)
dmax <- max(sapply(LINES, function(cl) { x <- sv[sv$line == cl, ]; max(abs(as.matrix(x[match(GR, x$KO), KOS]) - round(Kline[[cl]][GR, KOS], 3))) }))
cat(sprintf("max |difference| between these steps and the values plotted in fig5G (all 4 lines): %.4f\n", dmax)); stopifnot(dmax < 1e-9)
## ---- the intermediate tables of the example line ----
mm <- meta[meta$Cell_Line == CLX, ]; smp <- unlist(lapply(GR, function(g) sort(rownames(mm)[mm$Gene_Control == g]))); sg <- mm[smp, "Gene_Control"]
S1 <- lg[KOS, smp]; S2 <- Mfull[[CLX]]; mu <- rowMeans(S2); S3 <- S2 - mu; sdv <- matrixStats::rowSds(S2); S4 <- S3 / sdv; S5a <- t(S4[KOS, GR]); rmu <- rowMeans(S5a); rsd <- matrixStats::rowSds(S5a); S5 <- (S5a - rmu) / rsd
stopifnot(max(abs(S4 - Zfull[[CLX]])) < 1e-12, max(abs(S5 - Kline[[CLX]])) < 1e-12)
long <- rbind(data.frame(step = "2 group median log2(CP10K+1)", row = rownames(S2)[row(S2)], col = colnames(S2)[col(S2)], value = as.vector(S2)), data.frame(step = "3 minus gene mean over 9 groups", row = rownames(S3)[row(S3)], col = colnames(S3)[col(S3)], value = as.vector(S3)),
              data.frame(step = "4 divided by gene SD over 9 groups", row = rownames(S4)[row(S4)], col = colnames(S4)[col(S4)], value = as.vector(S4)), data.frame(step = "5 row z over 8 targets (plotted)", row = rownames(S5)[row(S5)], col = colnames(S5)[col(S5)], value = as.vector(S5)))
write.csv(long, file.path(OUT, sprintf("spec97_fig5G_steps_%s.csv", gsub("[^A-Za-z0-9]", "", CLX))), row.names = FALSE)
ex <- "TCF3"; cat(sprintf("\nworked example, %s, target %s in the %s KO:\n  samples: %s\n  step 2 median = %.3f | gene mean over the 9 group medians = %.3f, SD = %.3f\n  step 3 = %.3f - %.3f = %.3f | step 4 = %.3f / %.3f = %.3f\n  %s KO row over the 8 targets: mean %.3f, SD %.3f | step 5 = (%.3f - %.3f) / %.3f = %.3f\n",
  CLX, ex, ex, paste(sprintf("%s %.3f", smp[sg == ex], S1[ex, sg == ex]), collapse = ", "), S2[ex, ex], mu[ex], sdv[ex], S2[ex, ex], mu[ex], S3[ex, ex], S3[ex, ex], sdv[ex], S4[ex, ex], ex, rmu[ex], rsd[ex], S4[ex, ex], rmu[ex], rsd[ex], S5[ex, ex]))

SEQ <- colorRampPalette(c("#f7f7f5", "#eda100", "#b3400f"))(100)
hm <- function(M, ttl, sub, pal, sym = TRUE, digits = 2, rcol = NULL, ccol = NULL, box = TRUE, side = NULL, side_lab = NULL, cex_num = 0.6, vsep = NULL) { n <- nrow(M); m <- ncol(M)
  zl <- if (sym) c(-1, 1) * max(abs(M), na.rm = TRUE) else range(M, na.rm = TRUE); par(mar = c(5.2, 5.6, 4, if (is.null(side)) 1 else 6.5))
  image(seq_len(m), seq_len(n), t(M)[, n:1, drop = FALSE], col = pal, zlim = zl, axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(n)) for (j in seq_len(m)) text(j, n - i + 1, formatC(M[i, j], format = "f", digits = digits), cex = cex_num)
  if (box) for (k in KOS) { i <- match(k, rownames(M)); j <- match(k, colnames(M)); if (!is.na(i) && !is.na(j)) rect(j - 0.5, n - i + 0.5, j + 0.5, n - i + 1.5, border = "black", lwd = 2) }
  if (!is.null(vsep)) abline(v = vsep, col = "white", lwd = 2)
  text(seq_len(m), 0.35, colnames(M), srt = 45, adj = 1, xpd = NA, cex = 0.7, font = 2, col = if (is.null(ccol)) INK2 else ccol); text(0.4, n:1, rownames(M), adj = 1, xpd = NA, cex = 0.75, font = 2, col = if (is.null(rcol)) INK2 else rcol); box(col = MUTED)
  if (!is.null(side)) { text(m + 0.7, n + 0.75, side_lab, adj = 0, xpd = NA, cex = 0.58, font = 2, col = INK2); text(m + 0.7, n:1, side, adj = 0, xpd = NA, cex = 0.6, col = INK2) }
  mtext(ttl, side = 3, adj = 0, line = 2.2, font = 2, cex = 0.8); mtext(sub, side = 3, adj = 0, line = 0.3, cex = 0.55, col = INK2) }
save_fig("fig5H_fig5G_normalisation_steps", 17, 11.5, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.6, 0.4), family = "sans"); gc <- c("grey40", KO_COL[KOS])
  S1p <- S1; colnames(S1p) <- sub("^[^_]+_", "", smp); hm(S1p, sprintf("1  per-sample log2(CP10K + 1), %s", CLX), "rows = the 8 target genes; columns = every sample of the line, grouped NTC then the 8 KOs", SEQ, sym = FALSE, rcol = KO_COL[KOS], ccol = gc[match(sg, GR)], box = FALSE, cex_num = 0.42, vsep = which(sg[-1] != sg[-length(sg)]) + 0.5)
  hm(S2, "2  median per group", "median of the samples of each group (NTC + 8 KOs); boxed = the gene's own KO", SEQ, sym = FALSE, rcol = KO_COL[KOS], ccol = gc, side = sprintf("%.2f   %.2f", mu, sdv), side_lab = "mean    SD\nover the 9 groups")
  hm(S3, "3  minus the gene's mean over the 9 groups", "step 2 minus its row mean (unweighted mean of the 9 group medians, NTC counted once)", CROP_PAL, rcol = KO_COL[KOS], ccol = gc)
  hm(S4, "4  divided by the gene's SD over the 9 groups (= fig5E)", "step 3 / row SD: every gene on the same scale (z-score across the 9 groups)", CROP_PAL, rcol = KO_COL[KOS], ccol = gc)
  hm(S5, "5  transposed, each group row z-scored over the 8 targets", "rows = groups, columns = targets; (step 4 - row mean) / row SD = value drawn in fig5G", CROP_PAL, rcol = gc, ccol = KO_COL[KOS], side = sprintf("%+.2f   %.2f", rmu, rsd), side_lab = "mean    SD\nof the step-4 row")
  n <- length(GR); m <- length(KOS); vm <- max(abs(unlist(Kline)), na.rm = TRUE); RAD <- c(0.47, 0.36, 0.25, 0.14)       # identical to spec69
  colfun <- function(v) if (is.na(v)) "grey85" else CROP_PAL[max(1, min(100, round((v + vm) / (2 * vm) * 99) + 1))]
  circ <- function(x, y, r, col, border = NA, lwd = 1) { a <- seq(0, 2 * pi, length.out = 60); polygon(x + r * cos(a), y + r * sin(a), col = col, border = border, lwd = lwd, xpd = NA) }
  par(mar = c(5.2, 5.6, 4, 6.5)); plot(NA, xlim = c(0.5, m + 0.5), ylim = c(0.5, n + 0.5), asp = 1, axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i"); li0 <- match(CLX, LINES)
  for (i in seq_len(n)) for (j in seq_len(m)) { y <- n - i + 1; for (li in seq_along(LINES)) circ(j, y, RAD[li], colfun(Kline[[LINES[li]]][GR[i], KOS[j]])); for (li in seq_along(LINES)[-1]) circ(j, y, RAD[li], NA, border = "white", lwd = 0.4) }
  for (k in KOS) { i <- match(k, GR); j <- match(k, KOS); rect(j - 0.5, n - i + 0.5, j + 0.5, n - i + 1.5, border = "black", lwd = 1.6) }
  text(seq_len(m), 0.35, KOS, srt = 45, adj = 1, xpd = NA, cex = 0.7, font = 2, col = KO_COL[KOS]); text(0.35, n:1, GR, adj = 1, xpd = NA, cex = 0.75, font = 2, col = gc)
  text(m + 0.8, n:(n - 3), sprintf("ring %d: %s", 1:4, sub("NCI-|-211H", "", LINES)), adj = 0, xpd = NA, cex = 0.62, font = ifelse(seq_along(LINES) == li0, 2, 1), col = INK2); text(m + 0.8, n - 4.2, "(1 = outer)", adj = 0, xpd = NA, cex = 0.55, col = INK2)
  text(m + 0.8, n - 5.6, sprintf("colour scale\n-%.2f .. +%.2f\n(max |value|,\nall 4 lines)", vm, vm), adj = c(0, 1), xpd = NA, cex = 0.55, col = INK2)
  mtext("6  fig5G: steps 1-5 done in every line, one ring per line", side = 3, adj = 0, line = 2.2, font = 2, cex = 0.8); mtext(sprintf("ring %d (%s) is coloured with the values of panel 5; same palette, one symmetric scale for all rings", li0, CLX), side = 3, adj = 0, line = 0.3, cex = 0.55, col = INK2)
  mtext(sprintf("How fig5G_KO_x_target_per_line_rings is normalised, step by step (real values, %s; 70-sample set)", CLX), outer = TRUE, side = 3, adj = 0, line = 2, font = 2, cex = 1.05)
  mtext(sprintf("computed with the same code as spec69 and checked against the values it plotted (max difference %.1e over all 4 lines); steps 3-5 are done separately within each cell line", dmax), outer = TRUE, side = 3, adj = 0, line = 0.7, cex = 0.65, col = INK2)
})
