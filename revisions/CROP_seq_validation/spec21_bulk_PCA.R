# Spec 21: PCA of the bulk RNA-seq RUN1 samples, one panel per cell line (a combined PCA across lines is dominated by cell-line differences, not KO),
# coloured by Gene_Control with the SAME palette as the CROP-seq guide PCA (KO_COL, shared via repro_common.R -> also used by spec17/fig21 and the cross-line heatmaps).
# Per line: top 3000 variable genes (log2 CP10K, mean CP10K >= 0.1), raw and after regressing the cell-cycle index (cc = S.Score + G2M.Score) out per gene, pooled over the line's samples.
# needs: bulk_cellcycle_RUN1.R
# outputs: fig25_bulk_PCA_by_KO.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R"))
b <- load_bulk(); cp10k <- b$cp10k; meta <- b$meta; lg <- log2(cp10k + 1)
ccb <- read.csv(file.path(OUT, "bulk_cellcycle_per_sample_RUN1.csv"))
meta$rep <- sub(".*_(R[0-9]+)$", "\\1", meta$Sample_ID)
NTOP <- 3000

pca_line <- function(cl, cc_adj) {
  s <- rownames(meta)[meta$Cell_Line == cl]; ok <- rowMeans(cp10k[, s]) >= 0.1; y <- lg[ok, s, drop = FALSE]
  if (cc_adj) { cc <- ccb$cc[match(s, ccb$sample)]; ccc <- cc - mean(cc); bb <- as.numeric((y - rowMeans(y)) %*% ccc) / sum(ccc^2); y <- y - outer(bb, ccc) }
  top <- names(sort(matrixStats::rowVars(y), decreasing = TRUE))[seq_len(min(NTOP, nrow(y)))]
  yt <- y[top, , drop = FALSE]; yt <- yt[apply(yt, 1, sd) > 0, , drop = FALSE]; p <- prcomp(t(yt), center = TRUE, scale. = TRUE)
  list(p = p, ve = p$sdev^2 / sum(p$sdev^2), samples = s)
}

save_fig("fig25_bulk_PCA_by_KO", 15, 8, function() {
  layout(matrix(1:8, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4.5, 0.5), mar = c(4, 4.2, 3.6, 1), family = "sans")
  for (cc_adj in c(FALSE, TRUE)) for (cl in LINES) {
    r <- pca_line(cl, cc_adj); X <- r$p$x[, 1:2]; s <- r$samples; ko <- meta[s, "Gene_Control"]; rg <- apply(X, 2, range); pad <- 0.18 * apply(rg, 2, diff)
    plot(X[, 1], X[, 2], xlim = rg[, 1] + c(-1, 1) * pad[1], ylim = rg[, 2] + c(-1, 1) * pad[2], pch = 21, bg = adjustcolor(KO_COL[ko], 0.9), col = "white", cex = 1.9, bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8,
         xlab = sprintf("PC1 (%.0f%%)", 100 * r$ve[1]), ylab = sprintf("PC2 (%.0f%%)", 100 * r$ve[2])); abline(h = 0, v = 0, col = GRID)
    text(X[, 1], X[, 2], paste(ko, meta[s, "rep"]), pos = 3, cex = 0.5, offset = 0.4)
    mtext(sprintf("%s (RUN1)", cl), side = 3, adj = 0, line = 1.9, font = 2, cex = 0.85, col = LINE_COL[[cl]])
    mtext(sprintf("%s, %d genes, %d samples", ifelse(cc_adj, "cell-cycle regressed", "raw"), length(r$p$center), length(s)), side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2)
    if (cl == LINES[length(LINES)]) legend("topright", names(KO_COL)[names(KO_COL) %in% ko], pt.bg = KO_COL[names(KO_COL) %in% ko], pch = 21, col = "white", bty = "n", cex = 0.55, pt.cex = 1.1, ncol = 2)
  }
  mtext("PCA of bulk RNA-seq KO samples (RUN1), coloured by target gene (same palette as the CROP-seq guide PCA)", outer = TRUE, side = 3, adj = 0, line = 2.8, font = 2, cex = 1.05)
  mtext("one panel per cell line (a combined PCA is dominated by cell-line differences); top 3000 variable genes per line, mean CP10K >= 0.1; top row raw, bottom row cell-cycle regressed", outer = TRUE, side = 3, adj = 0, line = 1.4, cex = 0.6, col = INK2)
})
