# Spec 47: repeat the screen-vs-bulk agreement check (own-KO vs cross-KO null, per TF) but using raw gene-level KO-NTC diffs across
# the top-N most variable genes (N = 1000/3000/5000, set via TOPN env var) INSTEAD OF the 20 Cm modules -- a much larger, independent
# feature space, following this pipeline's established "topN variable genes" convention (spec22/29/30: genes ranked by bulk rowVar,
# averaged across lines). No "technology correction" here (dropped per feedback -- it inflated matched/mismatched pairs equally and
# didn't survive its own specificity check). Instead, NCI-H2052 (the one line shared between screen and bulk) is shown as individual
# reference points overlaid on the plot -- a direct "same line, cross-platform only" noise/ceiling reference -- NOT folded into either
# the own-KO or cross-KO boxplot distributions.
# needs: spec_common.R (load_screen()), repro_common.R (load_bulk())
# outputs: spec47_topN_genes_agreement.csv (TOPN in filename), fig3F_topN_genes_own_vs_null_H2052ref.{png,pdf} (TOPN in filename)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
TOPN <- as.integer(Sys.getenv("TOPN", "3000"))
# DIFF_DIR=NTC_KO -> NTC - KO on both sides (Pearson r is invariant to flipping both); default KO_NTC
DIFF_DIR <- Sys.getenv("DIFF_DIR", "KO_NTC"); SGN <- if (DIFF_DIR == "NTC_KO") -1 else 1; DLAB <- if (DIFF_DIR == "NTC_KO") "NTC-KO" else "KO-NTC"; DSUF <- if (DIFF_DIR == "NTC_KO") "_NTCminusKO" else ""

## ---------------------------------------------------------------- screen: per-gene, per-KO-group median expression (log-norm.)
srt <- load_screen(); Y <- GetAssayData(srt, layer = "data"); grp_scr <- as.character(srt$merged_call)
common_ko <- intersect(unique(grp_scr), ORDER)

## ---------------------------------------------------------------- bulk: variable-gene ranking (bulk rowVar, averaged across lines -- spec22/29/30 convention)
b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1)
genes_shared <- intersect(rownames(Y), rownames(lg))
lv <- sapply(LINES, function(cl) matrixStats::rowVars(as.matrix(lg[genes_shared, rownames(meta)[meta$Cell_Line == cl], drop = FALSE])))
vord <- genes_shared[order(rowMeans(lv), decreasing = TRUE)]
topgenes <- vord[1:min(TOPN, length(vord))]
cat(sprintf("TOPN=%d requested, %d genes used (shared between screen and bulk, ranked by mean bulk rowVar across %d lines)\n", TOPN, length(topgenes), length(LINES)))

## ---------------------------------------------------------------- screen KO-NTC diff (median), top genes only
screen_median <- t(sapply(c("NTC", common_ko), function(g) apply(as.matrix(Y[topgenes, grp_scr == g, drop = FALSE]), 1, median)))
rownames(screen_median) <- c("NTC", common_ko)
screen_diff <- SGN * sweep(screen_median[common_ko, , drop = FALSE], 2, screen_median["NTC", ], "-")

## ---------------------------------------------------------------- bulk KO-NTC diff (median), per line, top genes only
bulk_diff_by_line <- function(cl) {
  mm <- meta[meta$Cell_Line == cl, ]; ntc_s <- rownames(mm)[mm$Gene_Control == "NTC"]
  t(sapply(common_ko, function(k) { ko_s <- rownames(mm)[mm$Gene_Control == k]
    if (!length(ko_s) || !length(ntc_s)) return(setNames(rep(NA, length(topgenes)), topgenes))
    SGN * (apply(lg[topgenes, ko_s, drop = FALSE], 1, median) - apply(lg[topgenes, ntc_s, drop = FALSE], 1, median)) }))
}
h2052_diff <- bulk_diff_by_line("NCI-H2052"); rownames(h2052_diff) <- common_ko   # reference only, n=1 NTC -- not used in boxplot distributions
OTHER_LINES <- setdiff(LINES, "NCI-H2052")
other_diff <- lapply(OTHER_LINES, function(cl) { d <- bulk_diff_by_line(cl); rownames(d) <- common_ko; d })
names(other_diff) <- OTHER_LINES

## ---------------------------------------------------------------- own-KO, cross-KO (null), and H2052-reference correlations
own_ko <- do.call(rbind, lapply(OTHER_LINES, function(cl) do.call(rbind, lapply(common_ko, function(k) {
  bl <- other_diff[[cl]][k, ]; if (all(is.na(bl))) return(NULL)
  data.frame(line = cl, KO = k, pearson_r = cor(screen_diff[k, ], bl, use = "complete.obs"))
}))))
cross_ko <- do.call(rbind, lapply(OTHER_LINES, function(cl) do.call(rbind, lapply(common_ko, function(screen_k) {
  do.call(rbind, lapply(setdiff(common_ko, screen_k), function(bulk_k) {
    bl <- other_diff[[cl]][bulk_k, ]; if (all(is.na(bl))) return(NULL)
    data.frame(line = cl, screen_KO = screen_k, bulk_KO = bulk_k, pearson_r = cor(screen_diff[screen_k, ], bl, use = "complete.obs"))
  }))
}))))
h2052_ref <- data.frame(KO = common_ko, pearson_r = sapply(common_ko, function(k) cor(screen_diff[k, ], h2052_diff[k, ], use = "complete.obs")))

write.csv(own_ko, file.path(OUT, sprintf("spec47_top%d_genes_own_KO%s.csv", TOPN, DSUF)), row.names = FALSE)
write.csv(cross_ko, file.path(OUT, sprintf("spec47_top%d_genes_cross_KO%s.csv", TOPN, DSUF)), row.names = FALSE)
write.csv(h2052_ref, file.path(OUT, sprintf("spec47_top%d_genes_H2052_reference%s.csv", TOPN, DSUF)), row.names = FALSE)

options(width = 200)
cat(sprintf("\n=== top%d genes: own-KO vs cross-KO (null) vs H2052 reference ===\n", TOPN))
cat(sprintf("own-KO median r    = %.3f\n", median(own_ko$pearson_r)))
cat(sprintf("cross-KO median r  = %.3f (separation = %+.3f)\n", median(cross_ko$pearson_r), median(own_ko$pearson_r) - median(cross_ko$pearson_r)))
cat(sprintf("H2052 (same-line, cross-platform only) median r = %.3f\n", median(h2052_ref$pearson_r)))
print(format(h2052_ref, digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- fig: per-TF, own-KO (coloured) vs cross-KO (grey, null) vs H2052 (reference points)
save_fig(sprintf("fig3F_top%d_genes_own_vs_null_H2052ref%s", TOPN, DSUF), 9.5, 5.5, function() {
  ko_ord <- names(sort(sapply(split(own_ko$pearson_r, own_ko$KO), median), decreasing = TRUE))
  n_ko <- length(ko_ord)
  par(mar = c(5, 4.8, 3.6, 1), family = "sans")
  plot(NA, xlim = c(0.5, n_ko + 0.5), ylim = c(-1, 1), xaxt = "n", xlab = "",
       ylab = sprintf("Pearson r (screen vs bulk, top%d genes, %s diff)", TOPN, DLAB))
  axis(1, at = seq_len(n_ko), labels = FALSE)
  text(seq_len(n_ko), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), ko_ord, srt = 45, adj = 1, xpd = NA)
  for (i in seq_len(n_ko)) {
    d_null <- cross_ko$pearson_r[cross_ko$screen_KO == ko_ord[i]]
    boxplot(d_null, at = i - 0.18, add = TRUE, boxwex = 0.3, col = "grey85", border = "grey40", outline = FALSE, axes = FALSE)
  }
  for (i in seq_len(n_ko)) {
    d_own <- own_ko$pearson_r[own_ko$KO == ko_ord[i]]
    boxplot(d_own, at = i + 0.18, add = TRUE, boxwex = 0.3, col = adjustcolor(KO_COL[ko_ord[i]], 0.5), border = INK2, outline = FALSE, axes = FALSE)
  }
  set.seed(1); for (i in seq_len(n_ko)) { d <- own_ko[own_ko$KO == ko_ord[i], ]
    points(jitter(rep(i + 0.18, nrow(d)), amount = 0.06), d$pearson_r, pch = 16, cex = 1, col = LINE_COL[d$line]) }
  # H2052 reference: individual points, distinct marker, deliberately OUTSIDE/ABOVE the two boxes at each TF's position
  for (i in seq_len(n_ko)) points(i, h2052_ref$pearson_r[h2052_ref$KO == ko_ord[i]], pch = 18, cex = 1.8, col = LINE_COL[["NCI-H2052"]])
  abline(h = 0, col = INK, lty = 2)
  present_lines <- sort(unique(own_ko$line))
  legend("bottomleft", legend = c(present_lines, "NCI-H2052 (same-line ref.)", "own KO", "other KOs (null)"),
         col = c(LINE_COL[present_lines], LINE_COL[["NCI-H2052"]], "grey30", "grey60"),
         pch = c(rep(16, length(present_lines)), 18, 15, 15), bty = "n", cex = 0.6)
  mtext(sprintf("Screen-to-bulk correlation, top%d variable genes (%s diff), per TF", TOPN, DLAB), side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext("coloured = own bulk KO; grey = other bulk KOs (null); diamond = NCI-H2052 (same line as screen, cross-platform-only reference/noise ceiling)", side = 3, adj = 0, line = 0.9, cex = 0.58, col = INK2)
})
