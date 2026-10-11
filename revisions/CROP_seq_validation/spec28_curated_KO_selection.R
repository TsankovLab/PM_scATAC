# Spec 28: pick a curated sample set (NTC + KO replicates) from the 76-sample audited universe that minimizes the number of KO replicates
# whose own target-gene expression is HIGHER than the NTC reference for that line, while dropping as few samples as possible.
# NTC reference = median log2(CP10K+1) of ALL NTC replicates in that line (all NTC samples are always kept: dropping an NTC never fixes a
# KO>NTC violation, since the comparison is per KO-replicate vs the pooled NTC level, so removing NTCs only destabilises the reference).
# Per (line, KO) group: if not every replicate violates, drop only the violating ones (0 violations left, minimal sample loss). If EVERY
# replicate of a KO violates, dropping all of them would remove that KO entirely from that line -- instead keep the single least-bad
# replicate (smallest excess over the NTC reference) so the KO group survives with 1 residual, unavoidable violation.
# needs: spec27_knockdown_expression_heatmap.R (spec27_knockdown_expression_full76_<line>.csv)
# outputs: spec28_ko_vs_ntc_violations.csv (every KO replicate x gene check), spec28_curated_KO_selection.csv (keep/drop decision + reason),
#          fig2J_knockdown_expression_heatmap_curated.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))

fab <- load_full_audited_bulk(); meta <- fab$meta

res <- list()
for (cl in LINES) {
  f <- file.path(OUT, sprintf("spec27_knockdown_expression_full76_%s.csv", gsub("[^A-Za-z0-9]", "", cl)))
  M <- read.csv(f, row.names = 1, check.names = FALSE)
  grp <- meta[colnames(M), "Gene_Control"]
  for (g in rownames(M)) {
    ntc_med <- median(as.numeric(M[g, grp == "NTC"]))
    ko_cols <- which(grp == g); if (!length(ko_cols)) next
    for (cidx in ko_cols) res[[length(res) + 1]] <- data.frame(line = cl, KO = g, sample = colnames(M)[cidx], rep = meta[colnames(M)[cidx], "rep"],
      expr = M[g, cidx], ntc_median = ntc_med, violation = M[g, cidx] > ntc_med)
  }
}
res <- do.call(rbind, res); res$excess <- res$expr - res$ntc_median
write.csv(res, file.path(OUT, "spec28_ko_vs_ntc_violations.csv"), row.names = FALSE)

decide <- do.call(rbind, lapply(split(res, list(res$line, res$KO), drop = TRUE), function(d) {
  if (!any(d$violation)) { d$keep <- TRUE; d$reason <- "no violation"; return(d) }
  if (all(d$violation)) { best <- which.min(d$excess); d$keep <- FALSE; d$keep[best] <- TRUE
    d$reason <- ifelse(d$keep, "all reps violate - kept least-bad rep", "all reps violate - dropped (worse than kept rep)"); return(d) }
  d$keep <- !d$violation; d$reason <- ifelse(d$keep, "kept (non-violating)", "dropped (violates NTC)"); d
}))
rownames(decide) <- NULL; decide <- decide[order(decide$line, decide$KO, decide$rep), ]
write.csv(decide, file.path(OUT, "spec28_curated_KO_selection.csv"), row.names = FALSE)

options(width = 200)
cat(sprintf("KO replicates: %d total, %d kept, %d dropped\n", nrow(decide), sum(decide$keep), sum(!decide$keep)))
cat(sprintf("violations: %d before curation, %d remaining among kept samples (unavoidable: every replicate of that KO/line violated)\n", sum(decide$violation), sum(decide$violation & decide$keep)))
cat("\ndropped samples:\n"); print(decide[!decide$keep, c("line", "KO", "sample", "expr", "ntc_median", "reason")], row.names = FALSE)

## ---------------------------------------------------------------- fig2J: same layout as fig2H/fig2I, restricted to the curated sample set
kept_samples <- decide$sample[decide$keep]
screen_order <- c(ORDER, "MEF2A", "MEF2D")
col_order_for <- function(meta, samples, cl) {
  m <- meta[samples[samples %in% rownames(meta) & meta[samples, "Cell_Line"] == cl], , drop = FALSE]
  kos <- setdiff(sort(unique(m$Gene_Control), method = "radix"), "NTC")
  ko_ord <- c(intersect(screen_order, kos), setdiff(kos, screen_order))
  s <- unlist(lapply(c("NTC", ko_ord), function(g) sort(rownames(m)[m$Gene_Control == g])))
  list(samples = s, kos = ko_ord, meta = m)
}
all_kept_samples <- c(kept_samples, rownames(meta)[meta$Cell_Line %in% LINES & meta$Gene_Control == "NTC"])   # NTC always kept in full

save_fig("fig2J_knockdown_expression_heatmap_curated", 17, 5.4, function() {
  layout(matrix(1:4, nrow = 1)); par(oma = c(0, 0, 4.6, 0.5), mar = c(6.5, 5.5, 2.2, 1), family = "sans")
  for (cl in LINES) {
    co <- col_order_for(meta, all_kept_samples, cl); s <- co$samples; kos <- co$kos; m <- co$meta
    genes <- kos[kos %in% rownames(fab$lg)]
    M <- fab$lg[genes, s, drop = FALSE]; Z <- (M - rowMeans(M)) / matrixStats::rowSds(M); Z[!is.finite(Z)] <- 0
    grp <- m[s, "Gene_Control"]; n <- length(s); k <- length(genes)
    image(seq_len(n), seq_len(k), t(Z)[, k:1, drop = FALSE], col = DIV, zlim = c(-2.5, 2.5), axes = FALSE, xlab = "", ylab = "")
    edges <- which(grp[-1] != grp[-n]); for (e in edges) abline(v = e + 0.5, col = "white", lwd = 2)
    for (gi in seq_len(k)) text(par("usr")[1] - 0.015 * diff(par("usr")[1:2]), k - gi + 1, genes[gi], adj = 1, xpd = NA, cex = 0.7, font = 2, col = KO_COL[genes[gi]])
    axis(1, at = seq_len(n), labels = paste(grp, m[s, "rep"]), las = 2, tick = FALSE, cex.axis = 0.45, col.axis = INK2, line = -0.4)
    box(col = MUTED, lwd = 0.6)
    mtext(sprintf("%s (%d samples)", cl, n), side = 3, adj = 0, line = 0.7, font = 2, cex = 0.9, col = LINE_COL[[cl]])
  }
  mtext("Expression of each KO's own target gene, CURATED sample set (violating replicates dropped)", outer = TRUE, side = 3, adj = 0, line = 2.9, font = 2, cex = 1.05)
  mtext(sprintf("%d of the original 76 samples kept (%d KO replicates dropped, all NTC kept); a KO replicate was dropped when its own target-gene expression exceeded the NTC median for that line",
                length(all_kept_samples), sum(!decide$keep)), outer = TRUE, side = 3, adj = 0, line = 1.5, cex = 0.6, col = INK2)
  mtext("values = log2(CP10K+1), z-scored per gene across the CURATED samples of that line; white lines = KO group boundaries", outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.6, col = INK2)
  draw_colorbar(c(0.40, 0.60), c(-2.5, 2.5), "gene expression z-score (within cell line, curated set)")
})
