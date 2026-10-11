# Spec 27: expression of each KO's own target gene, every individual replicate, plus all NTC samples -- one heatmap per cell line.
# Rows = KO target genes (all bulk KOs present in that line, incl. SNAI2/RUNX1/RUNX2 where present); columns = every sample (KO reps + NTC),
# grouped by KO (CROP-seq screen order first, non-screen KOs after, NTC first), individual replicates within a group ordered by rep number.
# Values = log2(CP10K+1), z-scored per gene (row) across that line's own samples -- so a knockdown shows as a dark column-block in its own row.
# A raw (non-z-scored) log2 CP10K version of the same matrix is also written out per line as CSV.
# needs: nothing beyond load_bulk() / load_full_audited_bulk() (repro_common.R)
# outputs: spec27_knockdown_expression_<line>.csv (raw log2 CP10K, 70-sample set), fig2H_knockdown_expression_heatmap_by_line.{png,pdf} (70 samples),
#          spec27_knockdown_expression_full76_<line>.csv, fig2I_knockdown_expression_heatmap_by_line_full76.{png,pdf} (all 76 audited samples)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))

screen_order <- c(ORDER, "MEF2A", "MEF2D")   # CROP-seq screen order reference; MEF2A/MEF2D aren't in bulk, harmless if absent

col_order_for <- function(meta, cl) {
  m <- meta[meta$Cell_Line == cl, ]
  kos <- setdiff(sort(unique(m$Gene_Control), method = "radix"), "NTC")
  ko_ord <- c(intersect(screen_order, kos), setdiff(kos, screen_order))   # screen-shared KOs first (screen order), then any bulk-only KOs, alphabetical
  grp_order <- c("NTC", ko_ord)
  s <- unlist(lapply(grp_order, function(g) sort(rownames(m)[m$Gene_Control == g])))
  list(samples = s, kos = ko_ord, meta = m)
}

draw_knockdown_heatmap <- function(meta, lg, csv_prefix, fig_title) {
  layout(matrix(1:4, nrow = 1)); par(oma = c(0, 0, 4.6, 0.5), mar = c(6.5, 5.5, 2.2, 1), family = "sans")
  for (cl in LINES) {
    co <- col_order_for(meta, cl); s <- co$samples; kos <- co$kos; m <- co$meta
    genes <- kos[kos %in% rownames(lg)]
    M <- lg[genes, s, drop = FALSE]
    Z <- (M - rowMeans(M)) / matrixStats::rowSds(M); Z[!is.finite(Z)] <- 0
    write.csv(round(M, 3), file.path(OUT, sprintf("%s_%s.csv", csv_prefix, gsub("[^A-Za-z0-9]", "", cl))))

    grp <- m[s, "Gene_Control"]; n <- length(s); k <- length(genes)
    image(seq_len(n), seq_len(k), t(Z)[, k:1, drop = FALSE], col = DIV, zlim = c(-2.5, 2.5), axes = FALSE, xlab = "", ylab = "")
    edges <- which(grp[-1] != grp[-n]); for (e in edges) abline(v = e + 0.5, col = "white", lwd = 2)
    for (gi in seq_len(k)) text(par("usr")[1] - 0.015 * diff(par("usr")[1:2]), k - gi + 1, genes[gi], adj = 1, xpd = NA, cex = 0.7, font = 2, col = KO_COL[genes[gi]])
    axis(1, at = seq_len(n), labels = paste(grp, m[s, "rep"]), las = 2, tick = FALSE, cex.axis = 0.45, col.axis = INK2, line = -0.4)
    box(col = MUTED, lwd = 0.6)
    mtext(sprintf("%s (%d samples)", cl, n), side = 3, adj = 0, line = 0.7, font = 2, cex = 0.9, col = LINE_COL[[cl]])
  }
  mtext(fig_title, outer = TRUE, side = 3, adj = 0, line = 2.9, font = 2, cex = 1.05)
  mtext("rows = KO target genes (row label coloured by gene); columns = individual samples, grouped by KO (NTC first, then CROP-seq screen order, then bulk-only KOs); white lines = group boundaries",
        outer = TRUE, side = 3, adj = 0, line = 1.5, cex = 0.6, col = INK2)
  mtext("values = log2(CP10K+1), z-scored per gene across that line's own samples; a successful knockdown shows as a dark block sitting in its own gene's row, under its own group of columns",
        outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.6, col = INK2)
  draw_colorbar(c(0.40, 0.60), c(-2.5, 2.5), "gene expression z-score (within cell line)")
}

## ---------------------------------------------------------------- fig2H: the 70-sample set (load_bulk(), same scope as spec22's KO-reproducibility analysis)
b <- load_bulk()
save_fig("fig2H_knockdown_expression_heatmap_by_line", 17, 5.4, function() {
  draw_knockdown_heatmap(b$meta, log2(b$cp10k + 1), "spec27_knockdown_expression", "Expression of each KO's own target gene, every replicate + all NTC, by cell line")
})

## ---------------------------------------------------------------- fig2I: the full 76-sample audited universe (load_full_audited_bulk()) -- adds the 6 orphaned/undocumented
## samples (with their corrected Cell_Line); note this universe excludes SNAI2/RUNX1/RUNX2 by construction (see load_full_audited_bulk(), repro_common.R), unlike fig2H
fab <- load_full_audited_bulk()
save_fig("fig2I_knockdown_expression_heatmap_by_line_full76", 17, 5.4, function() {
  draw_knockdown_heatmap(fab$meta, fab$lg, "spec27_knockdown_expression_full76", "Expression of each KO's own target gene, all 76 audited samples (incl. orphaned), by cell line")
})
