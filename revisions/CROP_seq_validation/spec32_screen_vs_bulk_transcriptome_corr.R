# Spec 32: whole-transcriptome correspondence between the CROP-seq screen and bulk RNA-seq, for the SAME KO, per cell line --
# gene-level Spearman correlation of screen's log-norm mean diff-vs-NTC (spec_screen_stats.rds$diff0) against bulk's DESeq2
# log2FC-vs-NTC (spec_bulk_stats.rds$lfc0), swept over 5 feature-space sizes (top500/1000/3000/5000 variable genes, whole transcriptome).
# needs: spec1_prepare_stats.R (spec_screen_stats.rds, spec_bulk_stats.rds)
# outputs: spec_screen_vs_bulk_transcriptome_corr.csv, fig2K_screen_vs_bulk_transcriptome_corr_summary.{png,pdf},
#          fig2M_screen_vs_bulk_transcriptome_corr_top1000.{png,pdf} (same heatmap, fixed at top1000 for every line instead of each line's own best space),
#          fig2L_screen_vs_bulk_scatter_PITX1_HMGA1.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))

ss <- readRDS(file.path(OUT, "spec_screen_stats.rds")); bs <- readRDS(file.path(OUT, "spec_bulk_stats.rds"))
b <- load_bulk(); lg <- log2(b$cp10k + 1); meta <- b$meta
SPACES <- c("top500", "top1000", "top3000", "top5000", "all_genes"); KOS <- ORDER

feats_cache <- list(); genes_cache <- list()
res <- list()
for (cl in LINES) {
  common_genes <- intersect(ss$genes, bs[[cl]]$genes); genes_cache[[cl]] <- common_genes
  s <- rownames(meta)[meta$Cell_Line == cl]
  vord <- names(sort(setNames(matrixStats::rowVars(lg[common_genes, s]), common_genes), decreasing = TRUE))
  feats_cache[[cl]] <- list(top500 = vord[1:500], top1000 = vord[1:1000], top3000 = vord[1:3000], top5000 = vord[1:min(5000, length(vord))], all_genes = common_genes)
  for (sp in SPACES) {
    feats <- feats_cache[[cl]][[sp]]
    for (k in KOS) {
      if (!(k %in% colnames(ss$diff0)) || !(k %in% colnames(bs[[cl]]$lfc0))) next
      x <- ss$diff0[feats, k]; y <- bs[[cl]]$lfc0[feats, k]; ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 10) next
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      res[[length(res) + 1]] <- data.frame(line = cl, space = sp, KO = k, n_genes = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
    }
  }
}
res <- do.call(rbind, res); write.csv(res, file.path(OUT, "spec_screen_vs_bulk_transcriptome_corr.csv"), row.names = FALSE)
summ <- aggregate(rho ~ line + space, res, mean)
best_space <- setNames(sapply(LINES, function(cl) { d <- summ[summ$line == cl, ]; d$space[which.max(d$rho)] }), LINES)
options(width = 200); cat("best space per line:\n"); print(best_space)

## ---------------------------------------------------------------- fig2K: heatmap (KO x line, at each line's best space) + mean-rho-by-space summary
best <- do.call(rbind, lapply(LINES, function(cl) res[res$line == cl & res$space == best_space[[cl]], ]))
ko_ord <- names(sort(tapply(best$rho, best$KO, mean), decreasing = TRUE))
M <- matrix(NA, length(ko_ord), length(LINES), dimnames = list(ko_ord, LINES))
P <- M
for (i in seq_len(nrow(best))) M[best$KO[i], best$line[i]] <- best$rho[i]
for (i in seq_len(nrow(best))) P[best$KO[i], best$line[i]] <- best$p[i]

save_fig("fig2K_screen_vs_bulk_transcriptome_corr_summary", 11, 5.6, function() {
  layout(matrix(1:2, nrow = 1), widths = c(1.3, 1)); par(oma = c(0, 0, 3.6, 0), family = "sans")
  par(mar = c(7, 7, 2.5, 4.5))
  vmax <- max(abs(M), na.rm = TRUE)
  image(seq_len(ncol(M)), seq_len(nrow(M)), t(M)[, nrow(M):1, drop = FALSE], col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
  axis(1, at = seq_len(ncol(M)), labels = colnames(M), las = 2, tick = FALSE, cex.axis = 0.75)
  for (ri in seq_len(nrow(M))) text(0.15, nrow(M) - ri + 1, rownames(M)[ri], adj = 1, xpd = NA, cex = 0.85, font = 2, col = KO_COL[rownames(M)[ri]])
  for (ri in seq_len(nrow(M))) for (ci in seq_len(ncol(M))) {
    v <- M[ri, ci]; pv <- P[ri, ci]
    if (is.na(v)) next
    lab <- sprintf("%.2f%s", v, ifelse(is.na(pv), "", ifelse(pv < 0.001, "***", ifelse(pv < 0.01, "**", ifelse(pv < 0.05, "*", "")))))
    text(ci, nrow(M) - ri + 1, lab, cex = 0.62, col = if (abs(v) > vmax * 0.55) "white" else INK)
  }
  box(col = MUTED, lwd = 0.6)
  mtext("A  Screen vs bulk correlation (KO signature, gene level)", side = 3, adj = 0, line = 1.3, font = 2, cex = 0.85)
  mtext("rho at each line's best-performing feature space; *** p<0.001 ** p<0.01 * p<0.05", side = 3, adj = 0, line = 0.2, cex = 0.55, col = INK2)

  par(mar = c(5, 4.5, 2.5, 1))
  cols <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"), LINES)
  plot(NA, xlim = c(0.7, length(SPACES) + 0.3), ylim = range(summ$rho) + c(-0.02, 0.02), xaxt = "n", xlab = "", ylab = "mean rho across the 8 KOs", bty = "n", las = 1, cex.axis = 0.8, cex.lab = 0.8)
  abline(h = 0, col = GRID); axis(1, at = seq_along(SPACES), labels = SPACES, las = 2, tick = FALSE, cex.axis = 0.75)
  for (cl in LINES) { d <- summ[summ$line == cl, ][match(SPACES, summ$space[summ$line == cl]), ]; lines(seq_along(SPACES), d$rho, col = cols[cl], lwd = 2); points(seq_along(SPACES), d$rho, pch = 21, bg = cols[cl], col = "white", cex = 1.3) }
  legend("topright", LINES, col = cols[LINES], lwd = 2, bty = "n", cex = 0.65)
  mtext("B  Mean rho by feature-space size", side = 3, adj = 0, line = 1.3, font = 2, cex = 0.85)
  mtext("smaller/more-variable gene sets modestly outperform whole transcriptome", side = 3, adj = 0, line = 0.2, cex = 0.55, col = INK2)
  mtext("Whole-transcriptome correspondence between CROP-seq screen and bulk RNA-seq, per cell line, same KO", outer = TRUE, side = 3, adj = 0, line = 1.6, font = 2, cex = 1.0)
  draw_colorbar(c(0.03, 0.30), c(-vmax, vmax), "Spearman rho (screen diff vs bulk log2FC, panel A)", cex = 0.55)
})

## ---------------------------------------------------------------- fig2M: same heatmap as fig2K panel A, but at a FIXED feature space (top1000) for every line, instead of each line's own best space
fixed_space <- "top1000"
dF <- res[res$space == fixed_space, ]
ko_ordF <- names(sort(tapply(dF$rho, dF$KO, mean), decreasing = TRUE))
MF <- matrix(NA, length(ko_ordF), length(LINES), dimnames = list(ko_ordF, LINES)); PF <- MF
for (i in seq_len(nrow(dF))) { MF[dF$KO[i], dF$line[i]] <- dF$rho[i]; PF[dF$KO[i], dF$line[i]] <- dF$p[i] }
save_fig("fig2M_screen_vs_bulk_transcriptome_corr_top1000", 7, 5.6, function() {
  par(mar = c(7, 7, 4.5, 4.5), family = "sans")
  vmax <- max(abs(MF), na.rm = TRUE)
  image(seq_len(ncol(MF)), seq_len(nrow(MF)), t(MF)[, nrow(MF):1, drop = FALSE], col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
  axis(1, at = seq_len(ncol(MF)), labels = colnames(MF), las = 2, tick = FALSE, cex.axis = 0.8)
  for (ri in seq_len(nrow(MF))) text(0.15, nrow(MF) - ri + 1, rownames(MF)[ri], adj = 1, xpd = NA, cex = 0.9, font = 2, col = KO_COL[rownames(MF)[ri]])
  for (ri in seq_len(nrow(MF))) for (ci in seq_len(ncol(MF))) {
    v <- MF[ri, ci]; pv <- PF[ri, ci]; if (is.na(v)) next
    lab <- sprintf("%.2f%s", v, ifelse(is.na(pv), "", ifelse(pv < 0.001, "***", ifelse(pv < 0.01, "**", ifelse(pv < 0.05, "*", "")))))
    text(ci, nrow(MF) - ri + 1, lab, cex = 0.7, col = if (abs(v) > vmax * 0.55) "white" else INK)
  }
  box(col = MUTED, lwd = 0.6)
  mtext("Screen vs bulk correlation, fixed at top1000 variable genes (every line)", side = 3, adj = 0, line = 2.6, font = 2, cex = 0.88)
  mtext("Spearman rho, screen Δ vs bulk log2FC, per KO per line; *** p<0.001 ** p<0.01 * p<0.05", side = 3, adj = 0, line = 1.4, cex = 0.55, col = INK2)
  draw_colorbar(c(0.60, 0.98), c(-vmax, vmax), "Spearman rho (top1000 space)", cex = 0.55)
})

## ---------------------------------------------------------------- fig2L: gene-level scatter, PITX1 (best) vs HMGA1 (worst/inconsistent), all 4 lines
save_fig("fig2L_screen_vs_bulk_scatter_PITX1_HMGA1", 15, 7, function() {
  layout(matrix(1:8, nrow = 2, byrow = TRUE)); par(oma = c(0, 4, 4, 0.5), family = "sans")
  for (k in c("PITX1", "HMGA1")) {
    for (cl in LINES) {
      sp <- best_space[[cl]]; feats <- feats_cache[[cl]][[sp]]
      x <- ss$diff0[feats, k]; y <- bs[[cl]]$lfc0[feats, k]; ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
      r <- res[res$line == cl & res$space == sp & res$KO == k, ]
      par(mar = c(4, 4.2, 3.6, 1))
      plot(x, y, pch = 16, cex = 0.35, col = adjustcolor(LINE_COL[[cl]], 0.35), bty = "n", las = 1, cex.axis = 0.7, cex.lab = 0.75,
           xlab = "screen Δ (log-norm, vs NTC)", ylab = "bulk log2FC (vs NTC)")
      abline(h = 0, v = 0, col = GRID); abline(lm(y ~ x), col = INK, lwd = 1.3)
      mtext(sprintf("%s, %s (%s space)", k, cl, sp), side = 3, adj = 0, line = 1.7, font = 2, cex = 0.75, col = LINE_COL[[cl]])
      mtext(sprintf("Spearman ρ=%+.3f (p=%.2g), n=%d genes", r$rho, r$p, r$n_genes), side = 3, adj = 0, line = 0.5, cex = 0.6, col = INK2)
    }
  }
  mtext("Gene-level screen vs bulk correspondence: PITX1 (top row, best-reproducing KO) vs HMGA1 (bottom row, inconsistent/wrong-sign)", outer = TRUE, side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext("each dot = one gene; line = OLS fit; feature space = that line's best-performing space from fig2K panel B", outer = TRUE, side = 3, adj = 0, line = 0.9, cex = 0.6, col = INK2)
})
