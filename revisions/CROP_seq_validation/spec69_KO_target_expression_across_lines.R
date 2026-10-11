# Spec 69: knockdown check -- expression of each KO's OWN target gene (bulk only), NTC vs its own KO, across all 4 cell lines, as ONE heatmap. Group-median
# log2 CP10K per (line, NTC/KO) column, z-scored PER ROW (per TF) across all 8 columns jointly -- so, unlike spec27 (which z-scores separately within each
# line's own samples), knockdown depth is directly comparable across lines on the same colour scale.
# needs: repro_common.R (load_bulk())
# outputs: spec69_KO_target_expression.csv, fig5D_KO_target_expression_heatmap.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
KOS <- ORDER; b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1)
get_vals <- function(cl, gene) { if (!gene %in% rownames(lg)) return(NULL); mm <- meta[meta$Cell_Line == cl, ]; ntc <- rownames(mm)[mm$Gene_Control == "NTC"]; ko <- rownames(mm)[mm$Gene_Control == gene]
  if (!length(ntc) || !length(ko)) return(NULL); list(NTC = lg[gene, ntc], KO = lg[gene, ko]) }
lin <- function(x) pmax(2^x - 1, 0)   # bulk lg = log2(CP10K+1); invert to linear CP10K for the %% knockdown ratio

rows <- list()
for (k in KOS) for (cl in LINES) { v <- get_vals(cl, k); if (is.null(v)) next
  p <- if (length(v$NTC) > 1 && length(v$KO) > 1) tryCatch(t.test(v$NTC, v$KO)$p.value, error = function(e) NA) else NA
  rows[[length(rows) + 1]] <- data.frame(KO = k, line = cl, n_NTC = length(v$NTC), n_KO = length(v$KO), median_NTC = median(v$NTC), median_KO = median(v$KO),
                                          mean_NTC = mean(v$NTC), mean_KO = mean(v$KO), pct_knockdown = 1 - (mean(lin(v$KO)) / mean(lin(v$NTC))), p = p) }
res <- do.call(rbind, rows); write.csv(res, file.path(OUT, "spec69_KO_target_expression.csv"), row.names = FALSE)
options(width = 200); cat("\n=== bulk KO target gene expression: NTC vs own KO, per line (%% knockdown = 1 - mean linear KO / mean linear NTC) ===\n")
for (k in KOS) { cat(sprintf("\n-- %s --\n", k)); print(format(res[res$KO == k, c("line", "n_NTC", "n_KO", "mean_NTC", "mean_KO", "pct_knockdown", "p")], digits = 3), row.names = FALSE) }
cat("\nmean %% knockdown per line (across the 8 KOs):\n"); print(round(tapply(res$pct_knockdown, res$line, mean, na.rm = TRUE), 3))
cat("\nmean %% knockdown per KO (across the 4 lines):\n"); print(round(sort(tapply(res$pct_knockdown, res$KO, mean, na.rm = TRUE), decreasing = TRUE), 3))

## ---------------------------------------------------------------- heatmap: rows = 8 target genes, columns = NTC/KO x 4 lines (8 cols), z-scored per row across all 8
short_l <- sub("NCI-|-211H", "", LINES)
M <- do.call(rbind, lapply(KOS, function(k) { v <- unlist(lapply(LINES, function(cl) { d <- get_vals(cl, k); if (is.null(d)) c(NA, NA) else c(median(d$NTC), median(d$KO)) }))
  setNames(v, paste(rep(short_l, each = 2), c("NTC", "KO"))) }))
rownames(M) <- KOS
Z <- (M - rowMeans(M, na.rm = TRUE)) / matrixStats::rowSds(M, na.rm = TRUE); Z[!is.finite(Z)] <- NA

save_fig("fig5D_KO_target_expression_heatmap", 8.5, 6.4, function() {
  n <- nrow(Z); m <- ncol(Z); par(mar = c(4.2, 5.5, 3.6, 1), family = "sans")
  image(seq_len(m), seq_len(n), t(Z)[, n:1], col = DIV, zlim = c(-max(abs(Z), na.rm = TRUE), max(abs(Z), na.rm = TRUE)), axes = FALSE, xlab = "", ylab = "")
  na_i <- which(is.na(Z), arr.ind = TRUE); if (nrow(na_i)) rect(na_i[, 2] - 0.5, n - na_i[, 1] + 0.5, na_i[, 2] + 0.5, n - na_i[, 1] + 1.5, col = "#efefef", border = NA)
  for (i in seq_len(n)) for (j in seq_len(m)) if (!is.na(M[i, j])) text(j, n - i + 1, sprintf("%.2f", M[i, j]), cex = 0.62, col = INK)
  vb <- seq(2.5, m - 0.5, by = 2); abline(v = vb, col = "white", lwd = 2.5)
  axis(1, at = seq_len(m), labels = FALSE); text(seq_len(m), 0.35, sub(".* ", "", colnames(Z)), srt = 0, adj = 0.5, xpd = NA, cex = 0.65, col = INK2)
  text((seq_len(m / 2) * 2 - 0.5), -0.55, short_l, xpd = NA, cex = 0.82, font = 2)
  for (ri in seq_len(n)) text(0.35, n - ri + 1, rownames(Z)[ri], adj = 1, xpd = NA, cex = 0.78, font = 2, col = KO_COL[rownames(Z)[ri]]); box(col = MUTED, lwd = 0.6)
  mtext("Bulk knockdown check: NTC vs each KO's own target gene, all 4 lines", side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext("colour = z-score per row (TF), scaled jointly across all 8 columns; numbers = median log2 CP10K; white lines separate lines", side = 3, adj = 0, line = 0.6, cex = 0.58, col = INK2)
})

## ---------------------------------------------------------------- fig5E: each target gene across ALL KO groups (NTC + 8 KOs) in each line, z-scored per gene WITHIN each line
## so the reference is the rest of the KOs (not a single NTC) and line baseline differences are removed. Own-KO cell boxed; rank 1 = lowest of the 9 groups.
GR <- c("NTC", KOS)
Mfull <- lapply(setNames(LINES, LINES), function(cl) { mm <- meta[meta$Cell_Line == cl, ]
  t(sapply(KOS, function(gene) sapply(GR, function(g) { s <- rownames(mm)[mm$Gene_Control == g]; if (!length(s) || !gene %in% rownames(lg)) NA else median(lg[gene, s]) }))) })
Zfull <- lapply(Mfull, function(M) { Z <- (M - rowMeans(M, na.rm = TRUE)) / matrixStats::rowSds(M, na.rm = TRUE); Z[!is.finite(Z)] <- NA; Z })
own <- do.call(rbind, lapply(LINES, function(cl) { Z <- Zfull[[cl]]; data.frame(line = cl, KO = KOS, own_z = sapply(KOS, function(k) Z[k, k]),
  rank_of_own = sapply(KOS, function(k) rank(Z[k, ], na.last = "keep")[k]), n_groups = sapply(KOS, function(k) sum(!is.na(Z[k, ])))) }))
write.csv(own, file.path(OUT, "spec69_KO_target_vs_other_KOs.csv"), row.names = FALSE)
cat("\n=== own-KO z-score of each target gene vs all 9 groups in its line (negative = lower than the other KOs); rank 1 = lowest ===\n")
print(round(tapply(own$own_z, list(own$KO, own$line), mean)[KOS, LINES], 2)); print(tapply(own$rank_of_own, list(own$KO, own$line), mean)[KOS, LINES])
cat("\nmean own-KO z per KO (across lines):\n"); print(round(sort(tapply(own$own_z, own$KO, mean)), 2))

save_fig("fig5E_KO_target_vs_other_KOs_heatmap", 15, 5.2, function() {
  layout(matrix(1:4, nrow = 1)); par(oma = c(0, 0, 3.4, 0), family = "sans"); vm <- 2.5
  for (cl in LINES) { Z <- Zfull[[cl]]; M <- Mfull[[cl]]; n <- nrow(Z); m <- ncol(Z); par(mar = c(4.6, 5, 2.4, 0.6))
    image(seq_len(m), seq_len(n), t(pmax(pmin(Z, vm), -vm))[, n:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = "")
    na_i <- which(is.na(Z), arr.ind = TRUE); if (nrow(na_i)) rect(na_i[, 2] - 0.5, n - na_i[, 1] + 0.5, na_i[, 2] + 0.5, n - na_i[, 1] + 1.5, col = "#efefef", border = NA)
    for (i in seq_len(n)) for (j in seq_len(m)) if (!is.na(Z[i, j])) text(j, n - i + 1, sprintf("%.1f", Z[i, j]), cex = 0.5, col = INK)
    for (i in seq_len(n)) { j <- match(rownames(Z)[i], colnames(Z)); rect(j - 0.5, n - i + 0.5, j + 0.5, n - i + 1.5, border = "black", lwd = 2) }
    abline(v = 1.5, col = "white", lwd = 3)
    axis(1, at = seq_len(m), labels = FALSE); text(seq_len(m), 0.3, colnames(Z), srt = 60, adj = 1, xpd = NA, cex = 0.7, col = c("grey40", KO_COL[KOS]))
    for (i in seq_len(n)) text(0.4, n - i + 1, rownames(Z)[i], adj = 1, xpd = NA, cex = 0.75, font = 2, col = KO_COL[rownames(Z)[i]]); box(col = MUTED, lwd = 0.6)
    mtext(sprintf("%s  (mean own-KO z %.2f)", cl, mean(diag(Z[KOS, KOS]), na.rm = TRUE)), side = 3, adj = 0, line = 0.4, font = 2, cex = 0.82, col = LINE_COL[[cl]]) }
  mtext("KO target gene expression across all KO groups, per cell line (rows = target gene, columns = KO/NTC group)", outer = TRUE, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext("median log2 CP10K, z-scored per gene within each line across NTC + 8 KOs; boxed = the gene's own KO (blue = knocked down relative to the other KOs)", outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.58, col = INK2)
})

## ---------------------------------------------------------------- fig5F: median across the 4 lines, rows = KO groups (NTC + 8 KOs), columns = TF targets, scaled per ROW (KO)
## Input = the per-line gene z-scores of fig5E (each target centred/scaled across the 9 groups within its line, so gene abundance -- e.g. HMGA1's much higher
## level -- does not dominate the row scaling). Median across lines, transpose, then z-score each KO row across the 8 targets: a specific knockdown shows as
## the boxed (own-target) cell being the most negative in its row.
Zmed <- apply(simplify2array(lapply(Zfull, function(Z) Z[KOS, GR])), c(1, 2), median, na.rm = TRUE)   # targets x groups
Km <- t(Zmed)                                                                                         # groups (rows) x targets (cols)
Kr <- (Km - rowMeans(Km)) / matrixStats::rowSds(Km); Kr[!is.finite(Kr)] <- NA
write.csv(data.frame(KO = rownames(Km), round(Km, 3), check.names = FALSE), file.path(OUT, "spec69_KO_x_target_median_across_lines.csv"), row.names = FALSE)
spec <- data.frame(KO = KOS, own_row_z = sapply(KOS, function(k) Kr[k, k]), own_rank_in_row = sapply(KOS, function(k) rank(Kr[k, ])[k]))
cat("\n=== KO x target (median across lines), row-scaled: own-target z and rank within its KO row (1 = most reduced of the 8 targets) ===\n"); print(format(spec, digits = 2), row.names = FALSE)

CROP_PAL <- rev(as.character(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 100)))   # same palette as the CROP-seq knockdown heatmap (CRISPR_cropseq.R)
save_fig("fig5F_KO_x_target_median_heatmap", 2.6, 2.45, res = 300, fun = function() {
  n <- nrow(Kr); m <- ncol(Kr); vm <- max(abs(Kr), na.rm = TRUE); par(mar = c(2.3, 2.9, 1.0, 1.5), family = "sans", mgp = c(1, 0.2, 0))
  image(seq_len(m), seq_len(n), t(Kr)[, n:1], col = CROP_PAL, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = "")
  for (k in KOS) { i <- match(k, rownames(Kr)); j <- match(k, colnames(Kr)); rect(j - 0.5, n - i + 0.5, j + 0.5, n - i + 1.5, border = "black", lwd = 0.9) }
  abline(h = n - 0.5, col = "white", lwd = 1.4)
  text(seq_len(m), 0.4, colnames(Kr), srt = 45, adj = 1, xpd = NA, cex = 0.42, font = 2, col = KO_COL[colnames(Kr)])
  text(0.42, n:1, rownames(Kr), adj = 1, xpd = NA, cex = 0.42, font = 2, col = c("grey40", KO_COL[KOS])); box(col = INK2, lwd = 0.4)
  mtext("KO target expression (median 4 lines, row z)", side = 3, adj = 0, line = 0.25, font = 2, cex = 0.45)
  x0 <- m + 0.7; x1 <- m + 0.95; ys <- seq(n * 0.3, n * 0.7, length.out = 101)
  rect(x0, ys[-101], x1, ys[-1], col = CROP_PAL, border = NA, xpd = NA); rect(x0, ys[1], x1, ys[101], border = INK2, lwd = 0.3, xpd = NA)
  text(x1 + 0.1, c(ys[1], ys[51], ys[101]), sprintf("%.1f", c(-vm, 0, vm)), adj = 0, cex = 0.34, xpd = NA, col = INK2); text((x0 + x1) / 2, ys[101] + 0.35, "z", cex = 0.36, xpd = NA, col = INK2)
})

## ---------------------------------------------------------------- fig5G: same KO x target layout as fig5F but no median -- each cell = 4 concentric rings, one per line
## (outer -> inner: H2052, MSTO, H2452, H28), each coloured by that line's own value: per-line gene z, transposed, z-scored per KO row within the line.
##
## This panel extends the 8 screen KOs with RUNX1, RUNX2 and SNAI2. Those three were knocked down in the bulk RNA-seq but are NOT in the CROP-seq
## screen, which is the only reason ORDER (and so fig5D-F) leaves them out -- see the note in repro_common.R. Coverage is uneven and is left as grey
## rings rather than imputed: RUNX1 exists only in MSTO-211H and NCI-H2052, SNAI2 is missing from NCI-H28. The matrices are therefore rebuilt here
## over the extended group set instead of reusing Zfull, so each target gene is centred across all 12 KO groups + NTC rather than the original 9;
## fig5D/E/F above are untouched and still use the 8-KO set.
EXTRA  <- c("RUNX1", "RUNX2", "SNAI2")   # bulk-only KOs, absent from the CROP-seq screen
KOS_X  <- c(KOS, EXTRA); GR_X <- c("NTC", KOS_X)
Mfull_x <- lapply(setNames(LINES, LINES), function(cl) { mm <- meta[meta$Cell_Line == cl, ]
  t(sapply(KOS_X, function(gene) sapply(GR_X, function(g) { s <- rownames(mm)[mm$Gene_Control == g]; if (!length(s) || !gene %in% rownames(lg)) NA else median(lg[gene, s]) }))) })
Zfull_x <- lapply(Mfull_x, function(M) { Z <- (M - rowMeans(M, na.rm = TRUE)) / matrixStats::rowSds(M, na.rm = TRUE); Z[!is.finite(Z)] <- NA; Z })
Kline <- lapply(Zfull_x, function(Z) { K <- t(Z[KOS_X, GR_X]); Kz <- (K - rowMeans(K, na.rm = TRUE)) / matrixStats::rowSds(K, na.rm = TRUE); Kz[!is.finite(Kz)] <- NA; Kz })
write.csv(do.call(rbind, lapply(LINES, function(cl) data.frame(line = cl, KO = rownames(Kline[[cl]]), round(Kline[[cl]], 3), check.names = FALSE))), file.path(OUT, "spec69_KO_x_target_per_line_rowz.csv"), row.names = FALSE)
cat("\n=== groups available per line (n samples; 0 = ring left grey) ===\n")
print(sapply(LINES, function(cl) sapply(GR_X, function(g) sum(meta$Cell_Line == cl & meta$Gene_Control == g))))
cat("\n=== per line: own-target rank within KO row (1 = most reduced of the", length(KOS_X), "targets) ===\n")
print(sapply(LINES, function(cl) sapply(KOS_X, function(k) rank(Kline[[cl]][k, ], na.last = "keep")[k])))
cat("\nmean own-target row z per KO (across lines with data):\n")
print(round(sort(sapply(KOS_X, function(k) mean(sapply(LINES, function(cl) Kline[[cl]][k, k]), na.rm = TRUE))), 2))
cat("\nexpression level of the added targets (mean log2 CP10K over all samples):\n")
print(round(sapply(EXTRA, function(g) mean(lg[g, ])), 2))

save_fig("fig5G_KO_x_target_per_line_rings", 4.6, 4.2, res = 300, fun = function() {
  n <- length(GR_X); m <- length(KOS_X); vm <- max(abs(unlist(Kline)), na.rm = TRUE); RAD <- c(0.47, 0.36, 0.25, 0.14)
  colfun <- function(v) if (is.na(v)) "grey85" else CROP_PAL[max(1, min(100, round((v + vm) / (2 * vm) * 99) + 1))]
  circ <- function(x, y, r, col) { a <- seq(0, 2 * pi, length.out = 60); polygon(x + r * cos(a), y + r * sin(a), col = col, border = NA, xpd = NA) }
  par(mar = c(2.6, 3.2, 2.0, 3.0), family = "sans"); plot(NA, xlim = c(0.5, m + 0.5), ylim = c(0.5, n + 0.5), asp = 1, axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  for (i in seq_len(n)) for (j in seq_len(m)) { y <- n - i + 1
    for (li in seq_along(LINES)) circ(j, y, RAD[li], colfun(Kline[[LINES[li]]][GR_X[i], KOS_X[j]]))
    for (li in seq_along(LINES)[-1]) { a <- seq(0, 2 * pi, length.out = 60); lines(j + RAD[li] * cos(a), y + RAD[li] * sin(a), col = "white", lwd = 0.25) } }
  for (k in KOS_X) { i <- match(k, GR_X); j <- match(k, KOS_X); rect(j - 0.5, n - i + 0.5, j + 0.5, n - i + 1.5, border = "black", lwd = 0.9) }
  abline(h = n - 0.5, col = "grey60", lwd = 0.5)
  ## dashed separator: screen KOs (left/top) vs the bulk-only additions
  abline(v = length(KOS) + 0.5, col = "grey35", lwd = 0.6, lty = 2); abline(h = n - length(KOS) - 0.5, col = "grey35", lwd = 0.6, lty = 2)
  text(seq_len(m), 0.35, KOS_X, srt = 45, adj = 1, xpd = NA, cex = 0.42, font = 2, col = KO_COL[KOS_X])
  text(0.35, n:1, GR_X, adj = 1, xpd = NA, cex = 0.42, font = 2, col = c("grey40", KO_COL[KOS_X]))
  mtext("KO target expression per cell line (row z)", side = 3, adj = 0, line = 1.05, font = 2, cex = 0.45)
  mtext("dashed = CROP-seq screen KOs | RUNX1/RUNX2/SNAI2 bulk-only; grey ring = line not assayed", side = 3, adj = 0, line = 0.35, cex = 0.3, col = INK2)
  x0 <- m + 0.75; x1 <- m + 1.0; ys <- seq(n * 0.45, n * 0.85, length.out = 101)
  rect(x0, ys[-101], x1, ys[-1], col = CROP_PAL, border = NA, xpd = NA); rect(x0, ys[1], x1, ys[101], border = INK2, lwd = 0.3, xpd = NA)
  text(x1 + 0.1, c(ys[1], ys[51], ys[101]), sprintf("%.1f", c(-vm, 0, vm)), adj = 0, cex = 0.34, xpd = NA, col = INK2); text((x0 + x1) / 2, ys[101] + 0.4, "z", cex = 0.36, xpd = NA, col = INK2)
  cx <- m + 0.95; cy <- n * 0.10; for (li in seq_along(LINES)) circ(cx, cy, RAD[li], c("grey25", "grey50", "grey72", "grey90")[li])
  for (li in seq_along(LINES)) { yy <- cy + RAD[li] - 0.055; segments(cx, yy, cx + 0.62, cy + 0.55 - (li - 1) * 0.32, col = INK2, lwd = 0.3, xpd = NA)
    text(cx + 0.66, cy + 0.55 - (li - 1) * 0.32, sub("NCI-|-211H", "", LINES[li]), adj = 0, cex = 0.3, xpd = NA, col = INK2) }
  text(cx, cy + 0.72, "rings (out -> in)", cex = 0.3, xpd = NA, col = INK2)
})
