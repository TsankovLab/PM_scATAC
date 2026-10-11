# Spec 36: repeat of spec35's KO-signature x bulk-KO AddModuleScore analysis (fig2T/fig2U), but with each cell line's own
# cell-cycle score regressed out of the log2(CP10K+1) bulk expression matrix FIRST (same per-line linear-regression-residual
# convention used throughout this pipeline for the cc-adjusted t statistic -- spec20/spec22/spec31's `reg()`), before scoring the
# screen-derived top50 (vs NTC) gene sets with Seurat::AddModuleScore. 6/76 bulk samples have no precomputed cc score
# (bulk_cellcycle_per_sample_RUN1.csv) and are dropped for this variant only (H28_TCF3_R1, H28_HMGA1_R2, H28_NTC_R1/2/3,
# H2452_TWIST1_R2) -- 70 samples remain.
# needs: spec34_screen_onevsrest_module_bulk.R (spec34_screen_top50_upregulated_genes.csv), bulk_cellcycle_per_sample_RUN1.csv
# outputs: spec36_KO_signature_x_bulk_KO_matrix_ccregressed.csv, spec36_per_line_specificity_summary_ccregressed.csv,
#          fig2T_ccregressed_KO_signature_x_bulk_KO_matrix.{png,pdf}, fig2U_ccregressed_KO_signature_x_bulk_KO_matrix_per_line.{png,pdf},
#          fig2T_ccregressed_<line>_KO_signature_x_bulk_KO_matrix.{png,pdf}, fig2T_ccregressed_per_line_allLines.pdf (multi-page)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat))

top50 <- read.csv(file.path(OUT, "spec34_screen_top50_upregulated_genes.csv")); top50 <- top50[top50$comparison == "vsNTC", ]
fab <- load_full_audited_bulk(); meta <- fab$meta; lg <- fab$lg
ccb <- read.csv(file.path(OUT, "bulk_cellcycle_per_sample_RUN1.csv"))
has_cc <- rownames(meta)[rownames(meta) %in% ccb$sample]
cat(sprintf("dropping %d/%d bulk samples with no precomputed cc score: %s\n", nrow(meta) - length(has_cc), nrow(meta),
            paste(meta[setdiff(rownames(meta), has_cc), "Sample_ID"], collapse = ", ")))
meta <- meta[has_cc, ]

## ---------------------------------------------------------------- regress cc out of lg, per line (same convention as spec20/22/31's reg())
lg_cc <- lg[, rownames(meta), drop = FALSE]
for (cl in LINES) {
  s <- rownames(meta)[meta$Cell_Line == cl]; cc <- ccb$cc[match(s, ccb$sample)]; ccc <- cc - mean(cc)
  M <- lg_cc[, s, drop = FALSE]; b <- as.numeric((M - rowMeans(M)) %*% ccc) / sum(ccc^2)
  lg_cc[, s] <- M - outer(b, ccc)
}

## ---------------------------------------------------------------- score the screen top50 (vsNTC) signatures in the cc-regressed bulk matrix via AddModuleScore
options(Seurat.object.assay.version = "v3")   # AddModuleScore errors ("subscript out of bounds") on the default Assay5 class here; classic Assay works
srt_bulk <- CreateSeuratObject(counts = as.matrix(fab$cs[, rownames(meta), drop = FALSE]))
srt_bulk <- SetAssayData(srt_bulk, slot = "data", new.data = as.matrix(lg_cc))
feats <- lapply(ORDER, function(k) top50$gene[top50$KO == k]); names(feats) <- NULL
srt_bulk <- AddModuleScore(srt_bulk, features = feats, name = "sig", seed = 1)
sig_score <- as.matrix(srt_bulk@meta.data[, paste0("sig", seq_along(ORDER))]); colnames(sig_score) <- ORDER; rownames(sig_score) <- colnames(srt_bulk)
sig_score <- sig_score[rownames(meta), ]

COLS_KO <- c("NTC", ORDER)
Mko <- matrix(NA_real_, nrow = length(ORDER), ncol = length(COLS_KO), dimnames = list(ORDER, COLS_KO))
for (i in ORDER) for (j in COLS_KO) { s <- rownames(meta)[meta$Gene_Control == j]; if (length(s)) Mko[i, j] <- median(sig_score[s, i]) }
write.csv(cbind(signature_KO = rownames(Mko), Mko), file.path(OUT, "spec36_KO_signature_x_bulk_KO_matrix_ccregressed.csv"), row.names = FALSE)
options(width = 200); cat("\n=== median AddModuleScore activity (cc-regressed bulk) of each screen KO's top50 (vsNTC) signature, per bulk KO condition ===\n")
print(round(Mko, 3))

row_z <- function(M) { mu <- rowMeans(M, na.rm = TRUE); sdv <- apply(M, 1, sd, na.rm = TRUE); sdv[sdv == 0] <- 1; sweep(sweep(M, 1, mu, "-"), 1, sdv, "/") }

save_fig("fig2T_ccregressed_KO_signature_x_bulk_KO_matrix", 7.4, 5.6, function() {
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
  mtext("Screen KO signature scored in bulk (cc regressed), row-scaled across bulk KO", side = 3, adj = 0, line = 2.0, font = 2, cex = 0.9)
  mtext("cell-cycle score regressed out of bulk expr per line first; rows z-scored across own row; columns = bulk KO (4 lines pooled, 70 samples w/ cc score)", side = 3, adj = 0, line = 0.8, cex = 0.58, col = INK2)
  draw_colorbar(c(0.06, 0.30), c(-vmax, vmax), "row z-score (AddModuleScore, cc-regressed)", cex = 0.55)
})

## ---------------------------------------------------------------- per-line matrices + summary
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
write.csv(line_summary, file.path(OUT, "spec36_per_line_specificity_summary_ccregressed.csv"), row.names = FALSE)
cat("\n=== per-line breakdown (cc-regressed bulk): overall level and diagonal (own-KO) specificity ===\n")
print(format(line_summary, digits = 3), row.names = FALSE)

save_fig("fig2U_ccregressed_KO_signature_x_bulk_KO_matrix_per_line", 13, 11, function() {
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
  mtext("Screen KO signature x bulk KO condition, per cell line, cc-regressed + row-scaled", outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 1.05)
  mtext("black outline = matching KO (row=col); cc regressed out of bulk expr per line first; each row z-scored within its own line; overall/specificity on RAW (unscaled) values", outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.62, col = INK2)
})

## ---------------------------------------------------------------- standalone per-line heatmaps + combined multi-page pdf
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
  mtext(sprintf("%s: screen KO signature scored in bulk (cc regressed), row-scaled", cl), side = 3, adj = 0, line = 2.0, font = 2, cex = 0.9, col = LINE_COL[[cl]])
  mtext("cc regressed out per line first; rows z-scored across own row; columns = bulk KO condition (this line only)", side = 3, adj = 0, line = 0.8, cex = 0.58, col = INK2)
  draw_colorbar(c(0.06, 0.30), c(-vmax, vmax), "row z-score (AddModuleScore, cc-regressed)", cex = 0.55)
}
pdf_path <- file.path(OUT, "fig2T_ccregressed_per_line_KO_signature_x_bulk_KO_matrix_allLines.pdf")
cairo_pdf(pdf_path, width = 7.4, height = 5.6, onefile = TRUE)
for (cl in LINES) draw_line_matrix(cl)
dev.off()
message("saved fig2T_ccregressed_per_line_KO_signature_x_bulk_KO_matrix_allLines (multi-page pdf, ", length(LINES), " pages)")
