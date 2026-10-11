# Spec 35: for each KO's top50-upregulated-vs-NTC gene set (from spec34, screen-derived), show the MEDIAN expression of that gene
# set, KO vs NTC, side by side across the CROP-seq screen and each of the 4 bulk cell lines. Raw (non-z-scored) log expression --
# this is the direct, literal "does the screen's median shift also show up in each bulk line" comparison.
# needs: spec34_screen_onevsrest_module_bulk.R (spec34_screen_top50_upregulated_genes.csv), spec_common.R (load_screen())
# outputs: spec35_median_expr_vs_ntc.csv, fig2R_top50_median_expr_vs_ntc.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))

top50 <- read.csv(file.path(OUT, "spec34_screen_top50_upregulated_genes.csv"))
top50 <- top50[top50$comparison == "vsNTC", ]

fab <- load_full_audited_bulk(); meta <- fab$meta; lg <- fab$lg
srt <- load_screen(); Yscr <- GetAssayData(srt, layer = "data"); grp_scr <- as.character(srt$merged_call)

SOURCES <- c("CROP-seq", LINES)
res <- list()
for (k in ORDER) {
  g <- top50$gene[top50$KO == k]; g_scr <- g[g %in% rownames(Yscr)]; g_bulk <- g[g %in% rownames(lg)]

  # CROP-seq screen: per-cell median across the gene set, then median across cells in each group
  cell_med <- Matrix::colMeans(Yscr[g_scr, , drop = FALSE])   # mean across genes per cell (median across ~50 genes per cell is mostly 0s at single-cell depth; mean is the standard choice here)
  m_ko <- median(cell_med[grp_scr == k]); m_ntc <- median(cell_med[grp_scr == "NTC"])
  res[[length(res) + 1]] <- data.frame(KO = k, source = "CROP-seq", median_ko = m_ko, median_ntc = m_ntc, delta = m_ko - m_ntc, n_ko = sum(grp_scr == k), n_ntc = sum(grp_scr == "NTC"))

  # bulk, per line: per-sample mean across the gene set, then median across that line's KO replicates / NTC replicates
  for (cl in LINES) {
    m <- meta[meta$Cell_Line == cl, ]; s_ko <- rownames(m)[m$Gene_Control == k]; s_ntc <- rownames(m)[m$Gene_Control == "NTC"]
    samp_med <- colMeans(lg[g_bulk, c(s_ko, s_ntc), drop = FALSE])
    m_ko <- median(samp_med[s_ko]); m_ntc <- median(samp_med[s_ntc])
    res[[length(res) + 1]] <- data.frame(KO = k, source = cl, median_ko = m_ko, median_ntc = m_ntc, delta = m_ko - m_ntc, n_ko = length(s_ko), n_ntc = length(s_ntc))
  }
}
res <- do.call(rbind, res)
write.csv(res, file.path(OUT, "spec35_median_expr_vs_ntc.csv"), row.names = FALSE)
options(width = 200); cat("=== median expression of each KO's top50 (screen, vs-NTC) genes: KO vs NTC, per source ===\n")
print(format(res, digits = 3), row.names = FALSE)
cat("\nsign agreement with CROP-seq (delta > 0 in both): \n")
for (k in ORDER) { scr_sign <- sign(res$delta[res$KO == k & res$source == "CROP-seq"])
  bulk_signs <- sign(res$delta[res$KO == k & res$source %in% LINES]); cat(sprintf("  %-8s screen delta=%.3f | agree in %d/4 lines\n", k, res$delta[res$KO==k & res$source=="CROP-seq"], sum(bulk_signs == scr_sign))) }

## ---------------------------------------------------------------- figure: one panel per KO, bars = delta (median KO - median NTC) per source
save_fig("fig2R_top50_median_expr_vs_ntc", 14, 6.2, function() {
  layout(matrix(1:8, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4, 0), family = "sans")
  cols <- c("CROP-seq" = LINE_COL[["CROP-seq"]], setNames(LINE_COL[LINES], LINES))
  for (k in ORDER) {
    d <- res[res$KO == k, ]; d <- d[match(SOURCES, d$source), ]
    par(mar = c (4.2, 4.2, 2.6, 0.8))
    ylim <- range(c(0, d$delta)) + c(-1, 1) * 0.15 * diff(range(c(0, d$delta)))
    bp <- barplot(d$delta, names.arg = SOURCES, col = adjustcolor(cols[SOURCES], 0.85), border = INK2, las = 2, cex.names = 0.68, cex.axis = 0.75,
                  ylim = ylim, ylab = "median expr, KO - NTC")
    abline(h = 0, col = INK)
    mtext(k, side = 3, adj = 0, line = 0.6, font = 2, cex = 0.9, col = KO_COL[k])
  }
  mtext("Median expression of each KO's top50 screen-upregulated (vs NTC) genes: KO minus NTC, per source", outer = TRUE, side = 3, adj = 0, line = 2.4, font = 2, cex = 1.0)
  mtext("bars = CROP-seq screen (single line) then each of the 4 bulk cell lines; positive = KO higher than NTC (matching direction found in the screen)", outer = TRUE, side = 3, adj = 0, line = 1.1, cex = 0.6, col = INK2)
})

## ---------------------------------------------------------------- same data as a heatmap: rows = KO, columns = source (screen + 4 lines)
save_fig("fig2S_top50_median_expr_vs_ntc_heatmap", 7, 5.2, function() {
  M <- matrix(NA_real_, nrow = length(ORDER), ncol = length(SOURCES), dimnames = list(ORDER, SOURCES))
  for (i in seq_len(nrow(res))) M[res$KO[i], res$source[i]] <- res$delta[i]
  n_r <- nrow(M); n_c <- ncol(M)
  par(mar = c(6.5, 7, 3.4, 1), family = "sans")
  vmax <- max(abs(M)); image(seq_len(n_c), seq_len(n_r), heat_z(M), col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
  abline(v = 1.5, col = "white", lwd = 2.2)
  for (i in seq_len(n_r)) for (j in seq_len(n_c)) text(j, n_r - i + 1, sprintf("%.2f", M[i, j]), cex = 0.62, col = ifelse(abs(M[i, j]) > 0.55 * vmax, "white", INK))
  axis(1, at = seq_len(n_c), labels = SOURCES, las = 2, tick = FALSE, cex.axis = 0.78, col.axis = INK2)
  for (ri in seq_len(n_r)) text(0.35, n_r - ri + 1, rownames(M)[ri], adj = 1, xpd = NA, cex = 0.85, font = 2, col = KO_COL[rownames(M)[ri]])
  box(col = MUTED, lwd = 0.6)
  mtext("Median expr of KO's top50 screen-upreg. genes: KO minus NTC", side = 3, adj = 0, line = 2.0, font = 2, cex = 0.95)
  mtext("rows = KO; columns = CROP-seq screen then each bulk cell line; positive = KO higher than NTC (screen direction)", side = 3, adj = 0, line = 0.8, cex = 0.58, col = INK2)
  draw_colorbar(c(0.06, 0.30), c(-vmax, vmax), "median expr, KO - NTC", cex = 0.55)
})

## ---------------------------------------------------------------- full KO-signature x bulk-KO-condition matrix: for each screen-derived
## top50 (vsNTC) gene set, score its ACTIVITY in bulk (all 76 samples) via Seurat::AddModuleScore (the same control-gene-corrected module
## score convention used throughout this pipeline for the 20 Cm modules, cm_sample_scores_AddModuleScore_bulk_RUN1.csv) -- NOT a plain
## z-score mean -- then take the MEDIAN score across every bulk KO condition's samples (pooled across all 4 cell lines) -- shows whether
## a KO's own signature is specific to its own bulk KO, or is elevated broadly across other bulk KOs too (NTC included as baseline column)
options(Seurat.object.assay.version = "v3")   # AddModuleScore errors ("subscript out of bounds") on the default Assay5 class here; classic Assay works
srt_bulk <- CreateSeuratObject(counts = as.matrix(fab$cs)); srt_bulk <- SetAssayData(srt_bulk, slot = "data", new.data = as.matrix(lg))
feats <- lapply(ORDER, function(k) top50$gene[top50$KO == k]); names(feats) <- NULL
srt_bulk <- AddModuleScore(srt_bulk, features = feats, name = "sig", seed = 1)
sig_score <- as.matrix(srt_bulk@meta.data[, paste0("sig", seq_along(ORDER))]); colnames(sig_score) <- ORDER; rownames(sig_score) <- colnames(srt_bulk)
sig_score <- sig_score[rownames(meta), ]
COLS_KO <- c("NTC", ORDER)
Mko <- matrix(NA_real_, nrow = length(ORDER), ncol = length(COLS_KO), dimnames = list(ORDER, COLS_KO))
for (i in ORDER) for (j in COLS_KO) Mko[i, j] <- median(sig_score[rownames(meta)[meta$Gene_Control == j], i])
write.csv(cbind(signature_KO = rownames(Mko), Mko), file.path(OUT, "spec35_KO_signature_x_bulk_KO_matrix.csv"), row.names = FALSE)
cat("\n=== median AddModuleScore activity of each screen KO's top50 (vsNTC) signature, across every bulk KO condition's pooled samples ===\n")
print(round(Mko, 3))

# row-scale (z-score each signature's own row across the 9 bulk-KO columns) -- each raw row is nearly flat (dominated by the signature's
# own baseline expression level), so scaling row-wise is what actually reveals the relative column-wise (bulk-KO) pattern within each row
row_z <- function(M) { mu <- rowMeans(M, na.rm = TRUE); sdv <- apply(M, 1, sd, na.rm = TRUE); sdv[sdv == 0] <- 1; sweep(sweep(M, 1, mu, "-"), 1, sdv, "/") }

save_fig("fig2T_KO_signature_x_bulk_KO_matrix", 7.4, 5.6, function() {
  Mz <- row_z(Mko); n_r <- nrow(Mz); n_c <- ncol(Mz)
  par(mar = c(6.5, 7, 3.4, 1), family = "sans")
  vmax <- max(abs(Mz), na.rm = TRUE)
  image(seq_len(n_c), seq_len(n_r), heat_z(Mz), col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
  abline(v = 1.5, col = "white", lwd = 2.2)
  for (i in seq_len(n_r)) { j <- match(rownames(Mz)[i], COLS_KO); rect(j - 0.5, n_r - i + 0.5, j + 0.5, n_r - i + 1.5, border = INK, lwd = 2) }
  for (i in seq_len(n_r)) for (j in seq_len(n_c)) text(j, n_r - i + 1, sprintf("%.2f", Mz[i, j]), cex = 0.6, col = ifelse(abs(Mz[i, j]) > 0.55 * vmax, "white", INK))
  axis(1, at = seq_len(n_c), labels = COLS_KO, las = 2, tick = FALSE, cex.axis = 0.78, col.axis = INK2)
  for (ri in seq_len(n_r)) text(0.35, n_r - ri + 1, rownames(Mz)[ri], adj = 1, xpd = NA, cex = 0.85, font = 2, col = KO_COL[rownames(Mz)[ri]])
  box(col = MUTED, lwd = 0.6)
  mtext("Screen KO signature (top50 vs NTC) scored in bulk, row-scaled across bulk KO", side = 3, adj = 0, line = 2.0, font = 2, cex = 0.9)
  mtext("rows = screen-derived gene-set KO (z-scored across its own row); columns = bulk KO condition (all 4 lines pooled); black outline = matching KO", side = 3, adj = 0, line = 0.8, cex = 0.58, col = INK2)
  draw_colorbar(c(0.06, 0.30), c(-vmax, vmax), "row z-score (AddModuleScore)", cex = 0.55)
})

## ---------------------------------------------------------------- same KO x KO matrix, split per cell line -- is one line responsible
## for the low overall level and/or the poor diagonal specificity seen in the pooled fig2T?
Mko_line <- list(); line_summary <- list()
for (cl in LINES) {
  m <- meta[meta$Cell_Line == cl, ]
  M <- matrix(NA_real_, nrow = length(ORDER), ncol = length(COLS_KO), dimnames = list(ORDER, COLS_KO))
  for (i in ORDER) for (j in COLS_KO) { s <- rownames(m)[m$Gene_Control == j]; if (length(s)) M[i, j] <- median(sig_score[s, i]) }
  Mko_line[[cl]] <- M
  diag_v <- sapply(ORDER, function(k) M[k, k]); offdiag_v <- M[, ORDER][row(M[, ORDER]) != col(M[, ORDER])]
  line_summary[[cl]] <- data.frame(line = cl, overall_mean = mean(M, na.rm = TRUE), diag_mean = mean(diag_v, na.rm = TRUE),
                                    offdiag_mean = mean(offdiag_v, na.rm = TRUE), specificity = mean(diag_v, na.rm = TRUE) - mean(offdiag_v, na.rm = TRUE),
                                    NTC_col_mean = mean(M[, "NTC"], na.rm = TRUE))
}
line_summary <- do.call(rbind, line_summary)
write.csv(line_summary, file.path(OUT, "spec35_per_line_specificity_summary.csv"), row.names = FALSE)
cat("\n=== per-line breakdown: overall level and diagonal (own-KO) specificity of the KO signature x bulk-KO matrix ===\n")
cat("(overall_mean: is this line's bulk generally pulled positive/negative regardless of KO identity?\n")
cat(" specificity = diag_mean - offdiag_mean: is a KO's own signature higher in its own bulk KO than in other bulk KOs?)\n")
print(format(line_summary, digits = 3), row.names = FALSE)

save_fig("fig2U_KO_signature_x_bulk_KO_matrix_per_line", 13, 11, function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.5, 0), family = "sans")
  Mko_line_z <- lapply(Mko_line, row_z); vmax <- max(abs(unlist(Mko_line_z)), na.rm = TRUE)
  for (cl in LINES) {
    M <- Mko_line_z[[cl]]; n_r <- nrow(M); n_c <- ncol(M)
    par(mar = c(6, 6.5, 3, 1))
    image(seq_len(n_c), seq_len(n_r), heat_z(M), col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
    abline(v = 1.5, col = "white", lwd = 2)
    for (i in seq_len(n_r)) { j <- match(rownames(M)[i], COLS_KO); rect(j - 0.5, n_r - i + 0.5, j + 0.5, n_r - i + 1.5, border = INK, lwd = 2) }
    for (i in seq_len(n_r)) for (j in seq_len(n_c)) text(j, n_r - i + 1, sprintf("%.2f", M[i, j]), cex = 0.55, col = ifelse(abs(M[i, j]) > 0.55 * vmax, "white", INK))
    axis(1, at = seq_len(n_c), labels = COLS_KO, las = 2, tick = FALSE, cex.axis = 0.68, col.axis = INK2)
    for (ri in seq_len(n_r)) text(0.35, n_r - ri + 1, rownames(M)[ri], adj = 1, xpd = NA, cex = 0.75, font = 2, col = KO_COL[rownames(M)[ri]])
    box(col = MUTED, lwd = 0.6)
    ls <- line_summary[line_summary$line == cl, ]
    mtext(sprintf("%s   (overall=%.2f, specificity=%.2f)", cl, ls$overall_mean, ls$specificity), side = 3, adj = 0, line = 0.8, font = 2, cex = 0.85, col = LINE_COL[[cl]])
  }
  mtext("Screen KO signature x bulk KO condition, per cell line, row-scaled -- which line drags down the pooled matrix?", outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 1.05)
  mtext("black outline = matching KO (row=col); each row z-scored within its own line; overall/specificity computed on the RAW (unscaled) values", outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.62, col = INK2)
})

## ---------------------------------------------------------------- same matrix, one full-size standalone heatmap per line (same layout/size as fig2T)
draw_line_matrix <- function(cl) {
  M <- row_z(Mko_line[[cl]]); n_r <- nrow(M); n_c <- ncol(M)
  par(mar = c(6.5, 7, 3.4, 1), family = "sans")
  vmax <- max(abs(M), na.rm = TRUE)
  image(seq_len(n_c), seq_len(n_r), heat_z(M), col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
  abline(v = 1.5, col = "white", lwd = 2.2)
  for (i in seq_len(n_r)) { j <- match(rownames(M)[i], COLS_KO); rect(j - 0.5, n_r - i + 0.5, j + 0.5, n_r - i + 1.5, border = INK, lwd = 2) }
  for (i in seq_len(n_r)) for (j in seq_len(n_c)) text(j, n_r - i + 1, sprintf("%.2f", M[i, j]), cex = 0.6, col = ifelse(abs(M[i, j]) > 0.55 * vmax, "white", INK))
  axis(1, at = seq_len(n_c), labels = COLS_KO, las = 2, tick = FALSE, cex.axis = 0.78, col.axis = INK2)
  for (ri in seq_len(n_r)) text(0.35, n_r - ri + 1, rownames(M)[ri], adj = 1, xpd = NA, cex = 0.85, font = 2, col = KO_COL[rownames(M)[ri]])
  box(col = MUTED, lwd = 0.6)
  mtext(sprintf("%s: screen KO signature (top50 vs NTC) scored in bulk, row-scaled", cl), side = 3, adj = 0, line = 2.0, font = 2, cex = 0.9, col = LINE_COL[[cl]])
  mtext("rows = screen-derived gene-set KO (z-scored across its own row); columns = bulk KO condition (this line only); black outline = matching KO", side = 3, adj = 0, line = 0.8, cex = 0.58, col = INK2)
  draw_colorbar(c(0.06, 0.30), c(-vmax, vmax), "row z-score (AddModuleScore)", cex = 0.55)
}
## ---------------------------------------------------------------- same 4 heatmaps bundled into ONE multi-page PDF (one line per page) -- no separate per-line files
pdf_path <- file.path(OUT, "fig2T_per_line_KO_signature_x_bulk_KO_matrix_allLines.pdf")
cairo_pdf(pdf_path, width = 7.4, height = 5.6, onefile = TRUE)
for (cl in LINES) draw_line_matrix(cl)
dev.off()
message("saved fig2T_per_line_KO_signature_x_bulk_KO_matrix_allLines (multi-page pdf, ", length(LINES), " pages)")
