# Spec 24: QC of the CROP-seq screen dataset (same cell filtering as CRISPR_cropseq.R, via load_screen()).
#   Analogous to CRISPR_cropseq.R's "QC per guide" panel (cells per KO + VlnPlot of nFeature/nCount/percent.mt), but coloured
#   with the same KO palette (KO_COL) used throughout this pipeline (bulk PCA, guide PCA, cross-line heatmaps) instead of ggplot defaults.
# A: cells per KO; B: log10(genes detected); C: log10(UMI count); D: % mitochondrial reads -- one violin per KO (tick = median)
# needs: nothing beyond the screen Seurat object (load_screen(), in spec_common.R)
# outputs: spec24_screen_QC.csv, fig28_screen_QC.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))

srt <- load_screen(); md <- srt@meta.data
KO_LEVELS <- c("NTC", ORDER, EXTRA); KO_LEVELS <- KO_LEVELS[KO_LEVELS %in% md$merged_call]
grp <- factor(md$merged_call, KO_LEVELS)
qc <- data.frame(KO = grp, log10_nFeature = log10(md$nFeature_RNA), log10_nCount = log10(md$nCount_RNA), pct_mito = md$percent.mt)
write.csv(qc, file.path(OUT, "spec24_screen_QC.csv"), row.names = FALSE)
n_cells <- table(grp)
options(width = 200); cat(sprintf("\nCROP-seq screen QC (%d cells total, %d KO groups):\n", nrow(qc), length(KO_LEVELS)))
print(do.call(rbind, lapply(KO_LEVELS, function(k) { d <- qc[qc$KO == k, ]
  data.frame(KO = k, n = nrow(d), median_genes = round(median(10^d$log10_nFeature)), median_UMI = round(median(10^d$log10_nCount)), median_pct_mito = round(median(d$pct_mito), 2)) })), row.names = FALSE)

## ---------------------------------------------------------------- figure (violins, coloured by KO_COL to match the rest of the pipeline)
k <- length(KO_LEVELS)
ko_axis <- function(at) { for (i in seq_along(at)) text(at[i], par("usr")[3] - diff(par("usr")[3:4]) * 0.03, KO_LEVELS[i], srt = 45, adj = 1, xpd = NA, font = 2, cex = 0.6, col = KO_COL[KO_LEVELS[i]]) }

panel_violin <- function(vals, ylab, main) {
  yl <- range(vals); yl <- yl + c(-1, 1) * 0.06 * diff(yl)
  plot(NA, xlim = c(0.3, k + 0.7), ylim = yl, xaxt = "n", xlab = "", ylab = ylab, bty = "n", las = 1, cex.axis = 0.65, cex.lab = 0.7)
  abline(h = pretty(yl), col = GRID, lwd = 0.6)
  for (i in seq_len(k)) {
    v <- vals[qc$KO == KO_LEVELS[i]]
    if (length(v) > 1 && sd(v) > 0) violin_at(v, i, KO_COL[[KO_LEVELS[i]]]) else points(i, v[1], pch = 16, cex = 0.4, col = KO_COL[[KO_LEVELS[i]]])
    q <- quantile(v, c(0.25, 0.5, 0.75)); segments(i - 0.16, q[2], i + 0.16, q[2], lwd = 1.6, col = INK)
  }
  ko_axis(seq_len(k)); mtext(main, side = 3, adj = 0, line = 0.5, font = 2, cex = 0.72)
}

save_fig("fig28_screen_QC", 9.5, 3.3, function() {
  layout(matrix(1:4, nrow = 1), widths = c(0.8, 1, 1, 1)); par(oma = c(3.4, 0, 2.3, 0.2), family = "sans")

  par(mar = c(0.2, 3.6, 1.2, 0.4)); bp <- barplot(as.numeric(n_cells), col = KO_COL[KO_LEVELS], border = INK, lwd = 0.5, space = 0.3, ylab = "cells", las = 1, cex.axis = 0.65, cex.lab = 0.7, xaxt = "n")
  ko_axis(bp); mtext("A  Cells per KO", side = 3, adj = 0, line = 0.5, font = 2, cex = 0.72)

  par(mar = c(0.2, 3.6, 1.2, 0.4)); panel_violin(qc$log10_nFeature, "log10(genes detected)", "B  Gene complexity")
  par(mar = c(0.2, 3.6, 1.2, 0.4)); panel_violin(qc$log10_nCount, "log10(UMI count)", "C  UMI count")
  par(mar = c(0.2, 3.6, 1.2, 0.4)); panel_violin(qc$pct_mito, "% mitochondrial reads", "D  Mitochondrial %")

  mtext(sprintf("CROP-seq screen QC (%d cells, %d KO groups, %s)", nrow(qc), k, SCREEN_LINE), outer = TRUE, side = 3, adj = 0, line = 1.1, font = 2, cex = 0.68)
  mtext("colour = KO, matches the bulk RNA-seq KO palette used throughout this pipeline; tick = median", outer = TRUE, side = 3, adj = 0, line = 0.2, cex = 0.42, col = INK2)
})
