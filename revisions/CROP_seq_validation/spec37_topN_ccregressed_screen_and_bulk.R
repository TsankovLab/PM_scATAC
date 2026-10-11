# Spec 37: like spec36 (bulk KO-signature x bulk-KO AddModuleScore matrix, cc regressed out of bulk expr first), but ALSO regress
# cell cycle out of the CROP-seq SCREEN before identifying the top upregulated DE genes per KO (uses ss$tcc from spec1_prepare_stats.R
# -- Welch t vs NTC with the screen's cc index regressed out pooled over all cells, instead of ss$t0's plain t), and tries TWO gene-set
# sizes (top50 and top25) rather than just top50. No individual per-line PDF/PNG files -- only the combined 2x2 grid and one
# multi-page PDF per topN.
# needs: spec1_prepare_stats.R (spec_screen_stats.rds, for ss$tcc), spec34...(bulk_genes via load_full_audited_bulk), bulk_cellcycle_per_sample_RUN1.csv
# outputs (per N in {50,25}): spec37_screen_topN_ccregressed_upregulated_genes.csv, spec37_KO_signature_x_bulk_KO_matrix_topN_ccregressed.csv,
#          spec37_per_line_specificity_summary_topN_ccregressed.csv, fig2V_topN_ccregressed_KO_signature_x_bulk_KO_matrix.{png,pdf},
#          fig2V_topN_ccregressed_KO_signature_x_bulk_KO_matrix_per_line.{png,pdf}, fig2V_topN_ccregressed_per_line_allLines.pdf (multi-page)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat))

fab <- load_full_audited_bulk(); meta_all <- fab$meta; lg_all <- fab$lg; bulk_genes <- rownames(lg_all)
ccb <- read.csv(file.path(OUT, "bulk_cellcycle_per_sample_RUN1.csv"))
has_cc <- rownames(meta_all)[rownames(meta_all) %in% ccb$sample]
cat(sprintf("dropping %d/%d bulk samples with no precomputed cc score: %s\n", nrow(meta_all) - length(has_cc), nrow(meta_all),
            paste(meta_all[setdiff(rownames(meta_all), has_cc), "Sample_ID"], collapse = ", ")))
meta <- meta_all[has_cc, ]

## ---------------------------------------------------------------- regress cc out of bulk lg, per line (same convention as spec20/22/31/36's reg())
lg_cc <- lg_all[, rownames(meta), drop = FALSE]
for (cl in LINES) {
  s <- rownames(meta)[meta$Cell_Line == cl]; cc <- ccb$cc[match(s, ccb$sample)]; ccc <- cc - mean(cc)
  M <- lg_cc[, s, drop = FALSE]; b <- as.numeric((M - rowMeans(M)) %*% ccc) / sum(ccc^2)
  lg_cc[, s] <- M - outer(b, ccc)
}

## ---------------------------------------------------------------- screen DE genes, WITH cell cycle regressed out (ss$tcc instead of spec34's ss$t0)
ss <- readRDS(file.path(OUT, "spec_screen_stats.rds"))
row_z <- function(M) { mu <- rowMeans(M, na.rm = TRUE); sdv <- apply(M, 1, sd, na.rm = TRUE); sdv[sdv == 0] <- 1; sweep(sweep(M, 1, mu, "-"), 1, sdv, "/") }
COLS_KO <- c("NTC", ORDER)
options(Seurat.object.assay.version = "v3")   # AddModuleScore errors ("subscript out of bounds") on the default Assay5 class here; classic Assay works
options(width = 200)

for (topN in c(50, 25)) {
  cat(sprintf("\n\n########## topN = %d, screen DE cc-regressed, bulk cc-regressed ##########\n", topN))
  top <- list()
  for (k in ORDER) {
    if (!k %in% colnames(ss$tcc)) next
    x <- ss$tcc[, k]; x <- x[names(x) %in% bulk_genes & is.finite(x)]
    top[[k]] <- data.frame(KO = k, gene = names(sort(x, decreasing = TRUE))[1:topN], t_cc = sort(x, decreasing = TRUE)[1:topN])
  }
  top <- do.call(rbind, top); rownames(top) <- NULL
  write.csv(top, file.path(OUT, sprintf("spec37_screen_top%d_ccregressed_upregulated_genes.csv", topN)), row.names = FALSE)

  srt_bulk <- CreateSeuratObject(counts = as.matrix(fab$cs[, rownames(meta), drop = FALSE]))
  srt_bulk <- SetAssayData(srt_bulk, slot = "data", new.data = as.matrix(lg_cc))
  feats <- lapply(ORDER, function(k) top$gene[top$KO == k]); names(feats) <- NULL
  srt_bulk <- AddModuleScore(srt_bulk, features = feats, name = "sig", seed = 1)
  sig_score <- as.matrix(srt_bulk@meta.data[, paste0("sig", seq_along(ORDER))]); colnames(sig_score) <- ORDER; rownames(sig_score) <- colnames(srt_bulk)
  sig_score <- sig_score[rownames(meta), ]

  Mko <- matrix(NA_real_, nrow = length(ORDER), ncol = length(COLS_KO), dimnames = list(ORDER, COLS_KO))
  for (i in ORDER) for (j in COLS_KO) { s <- rownames(meta)[meta$Gene_Control == j]; if (length(s)) Mko[i, j] <- median(sig_score[s, i]) }
  write.csv(cbind(signature_KO = rownames(Mko), Mko), file.path(OUT, sprintf("spec37_KO_signature_x_bulk_KO_matrix_top%d_ccregressed.csv", topN)), row.names = FALSE)
  cat(sprintf("=== median AddModuleScore activity (screen DE + bulk both cc-regressed, top%d) ===\n", topN)); print(round(Mko, 3))

  save_fig(sprintf("fig2V_top%d_ccregressed_KO_signature_x_bulk_KO_matrix", topN), 7.4, 5.6, function() {
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
    mtext(sprintf("Screen top%d DE genes (cc-regressed) scored in bulk (cc-regressed), row-scaled", topN), side = 3, adj = 0, line = 2.0, font = 2, cex = 0.9)
    mtext("cell cycle regressed out in BOTH the screen (gene selection) and bulk (scoring); rows z-scored across own row; columns = bulk KO (4 lines pooled)", side = 3, adj = 0, line = 0.8, cex = 0.58, col = INK2)
    draw_colorbar(c(0.06, 0.30), c(-vmax, vmax), "row z-score (AddModuleScore)", cex = 0.55)
  })

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
  write.csv(line_summary, file.path(OUT, sprintf("spec37_per_line_specificity_summary_top%d_ccregressed.csv", topN)), row.names = FALSE)
  cat(sprintf("\n=== per-line breakdown (top%d, screen+bulk cc-regressed) ===\n", topN)); print(format(line_summary, digits = 3), row.names = FALSE)

  save_fig(sprintf("fig2V_top%d_ccregressed_KO_signature_x_bulk_KO_matrix_per_line", topN), 13, 11, function() {
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
    mtext(sprintf("Screen top%d DE genes (cc-regressed) scored in bulk (cc-regressed), per line, row-scaled", topN), outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 1.05)
    mtext("black outline = matching KO (row=col); cc regressed out in both screen gene selection and bulk scoring; each row z-scored within its own line", outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.62, col = INK2)
  })

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
    mtext(sprintf("%s: screen top%d DE (cc-regressed) scored in bulk (cc-regressed), row-scaled", cl, topN), side = 3, adj = 0, line = 2.0, font = 2, cex = 0.9, col = LINE_COL[[cl]])
    mtext("cc regressed out in both screen gene selection and bulk scoring; rows z-scored across own row; columns = bulk KO (this line only)", side = 3, adj = 0, line = 0.8, cex = 0.58, col = INK2)
    draw_colorbar(c(0.06, 0.30), c(-vmax, vmax), "row z-score (AddModuleScore)", cex = 0.55)
  }
  pdf_path <- file.path(OUT, sprintf("fig2V_top%d_ccregressed_per_line_KO_signature_x_bulk_KO_matrix_allLines.pdf", topN))
  cairo_pdf(pdf_path, width = 7.4, height = 5.6, onefile = TRUE)
  for (cl in LINES) draw_line_matrix(cl)
  dev.off()
  message(sprintf("saved fig2V_top%d_ccregressed_per_line_KO_signature_x_bulk_KO_matrix_allLines (multi-page pdf, %d pages)", topN, length(LINES)))
}
