# Spec 25: PCA of the bulk RNA-seq samples, using the FULL, AUDITED sample universe (load_full_audited_bulk(), repro_common.R) -- sample list built
#   directly from counts/ (every *_ReadsPerGene.out.tab file), not just the curated sheet, catching samples entirely missing from the sheet
#   (H2452_TWIST1_R2, H28_HMGA1_R2, H28_NTC_R1/2/3). Every KO except SNAI2/RUNX1/RUNX2 (not in the CROP-seq screen), the 4 real lines, ALL statuses
#   kept (meta_all$removed_run1 flags status != "OK"). meta_all$Cell_Line is the corrected/genotyped line -- including the registered NTC
#   reassignment (H28_NTC_R1/2/3 -> NCI-H2052, H2052_NTC_R1/2 -> NCI-H28; spec25_NTC_reassignment_table.csv) -- so the per-line panels below group
#   samples by their ACTUAL identity, not the original sequencing-batch label.
# fig35 -- ONE combined PCA (top 5000 variable genes, whole transcriptome, raw), coloured by AS-NAMED/original cell line, no other classification --
#          plain dots, every one labelled with its Sample_ID (de-overlapped via maptools::pointLabel, thin leader line back to each point).
# fig36 -- per-line PCA with the REGISTERED NTC correction applied, each line's own top 5000 variable genes, raw, highlighting the reassigned samples.
# fig37 -- same as fig36, cell-cycle regressed, no colour/shape classification -- coloured by cell line only, every point labelled.
# fig38 -- same sample set/PCA as fig37, highlighting samples NOT present in RUN1 (status != OK).
# fig39 -- same sample set/PCA as fig37/38, highlighting VISUAL_OUTLIERS (repro_common.R): manually-flagged outliers from this PCA (NCI-H2052:
#          PC2 > 20; plus specific named samples in the other 3 lines). Also used to highlight fig27_QC_RUN1_screenKOs_byKO.
# fig40 -- same sample set/PCA as fig37/38/39, highlighting which of the 76 belong to the 70-sample set used by spec22's KO-reproducibility
#          analysis (load_bulk(), screen KOSET: NTC + the 8 CROP-seq-matching KOs) -- the 6 NOT in that 70 are the orphaned/undocumented samples
#          (H28_NTC_R1/2/3, H2452_TWIST1_R2, H28_HMGA1_R2, H28_TCF3_R1), all of which are also VISUAL_OUTLIERS (fig39).
# fig41 -- same idea as fig40, inverted: triangle = EXCLUDED (not in the 70) instead of in the 70, and uses the RAW (non-cc-regressed) per-line
#          PCA (pca_line4, same as fig36) instead of the cc-regressed one.
# needs: spec1_prepare_stats.R (not strictly required here, but keeps the CROP_OUT dir consistent with the rest of the pipeline)
# outputs: fig35_bulk_PCA_full_audit_plain_labelled.{png,pdf}, fig36_bulk_PCA_by_line_full_NTC_correction.{png,pdf},
#          fig37_bulk_PCA_by_line_full_NTC_correction_ccregressed.{png,pdf}, fig38_bulk_PCA_by_line_ccregressed_highlight_not_in_RUN1.{png,pdf},
#          fig39_bulk_PCA_by_line_ccregressed_highlight_visual_outliers.{png,pdf}, fig40_bulk_PCA_by_line_ccregressed_highlight_in_70_KOrepro.{png,pdf},
#          fig41_bulk_PCA_by_line_raw_highlight_not_in_70_KOrepro.{png,pdf}, spec25_full_audit_samples.csv, spec25_NTC_reassignment_table.csv
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))

## ==================================================================== the FULL, AUDITED sample universe -- the ONE canonical sample set used by every figure below and by spec23_QC_RUN1.R
fab <- load_full_audited_bulk(); meta_all <- fab$meta; cp10k_all <- fab$cp10k; lg_all <- fab$lg; cs <- fab$cs
write.csv(meta_all, file.path(OUT, "spec25_full_audit_samples.csv"), row.names = FALSE)
cat(sprintf("\nfull audited sample universe (every KO except SNAI2/RUNX1/RUNX2, 4 lines): %d samples, %d swapped cell line, %d removed in RUN1\n", nrow(meta_all), sum(meta_all$swapped), sum(meta_all$removed_run1)))
print(table(meta_all$Cell_Line_original, meta_all$flag))
cat("\nNTC per AS-NAMED/original cell line (ground truth for 'how many NTC were sequenced under this label'):\n"); print(table(meta_all$Cell_Line_original[meta_all$Gene_Control == "NTC"]))

## ---------------------------------------------------------------- fig35 PCA: top 5000 variable genes (whole transcriptome), raw
ok34 <- rowMeans(cp10k_all) >= 0.1; y34 <- lg_all[ok34, , drop = FALSE]
top34 <- names(sort(matrixStats::rowVars(y34), decreasing = TRUE))[seq_len(min(5000, nrow(y34)))]
yt34 <- y34[top34, , drop = FALSE]; yt34 <- yt34[apply(yt34, 1, sd) > 0, , drop = FALSE]
p34 <- prcomp(t(yt34), center = TRUE, scale. = TRUE); ve34 <- p34$sdev^2 / sum(p34$sdev^2)

## ---------------------------------------------------------------- fig35: PCA of the full audited sample set, coloured by AS-NAMED cell line, large canvas; labels de-overlapped (maptools::pointLabel) with leader lines back to each point
suppressPackageStartupMessages(library(maptools))
save_fig("fig35_bulk_PCA_full_audit_plain_labelled", 24, 20, function() {
  par(mar = c(4.5, 4.8, 3.8, 1), family = "sans")
  X <- p34$x[, 1:2]; cl <- meta_all$Cell_Line_original
  rg <- apply(X, 2, range); pad <- 0.12 * apply(rg, 2, diff)
  plot(X[, 1], X[, 2], xlim = rg[, 1] + c(-1, 1) * pad[1], ylim = rg[, 2] + c(-1, 1) * pad[2], type = "n", bty = "n", las = 1, cex.axis = 1.1, cex.lab = 1.2,
       xlab = sprintf("PC1 (%.0f%%)", 100 * ve34[1]), ylab = sprintf("PC2 (%.0f%%)", 100 * ve34[2])); abline(h = 0, v = 0, col = GRID)
  points(X[, 1], X[, 2], pch = 21, bg = adjustcolor(LINE_COL[cl], 0.85), col = "white", cex = 2.2)
  lbl <- pointLabel(X[, 1], X[, 2], meta_all$Sample_ID, cex = 0.55, doPlot = FALSE)
  segments(X[, 1], X[, 2], lbl$x, lbl$y, col = adjustcolor(INK2, 0.45), lwd = 0.4)
  text(lbl$x, lbl$y, meta_all$Sample_ID, cex = 0.55, col = INK2)
  legend("bottomright", LINES, pt.bg = LINE_COL[LINES], pch = 21, col = "white", bty = "n", pt.cex = 2.2, cex = 1.1)
  mtext(sprintf("PCA of all %d bulk RNA-seq samples (meso lines only, SNAI2/RUNX1/RUNX2 excluded), coloured by AS-NAMED/original cell line", nrow(meta_all)), side = 3, adj = 0, line = 2.3, font = 2, cex = 1.0)
  mtext("top 5000 variable genes (mean CP10K >= 0.1), raw; every point labelled with its own Sample_ID (de-overlapped; thin line = actual point)", side = 3, adj = 0, line = 1.1, cex = 0.75, col = INK2)
})

## ==================================================================== REGISTERED NTC reassignment (based on the fig35 clustering evidence)
#   H28_NTC_R1, H28_NTC_R2, H28_NTC_R3      ->  reassigned to NCI-H2052   (registered in H28_NTC_NATIVE, repro_common.R)
#   H2052_NTC_R1, H2052_NTC_R2              ->  reassigned to NCI-H28    (already registered in sample_metadata_RUN1.tsv itself: label_corrected = True)
# This table is the single source of truth for both corrections; meta_all (built from H28_NTC_NATIVE + the sheet) already reflects it automatically.
ntc_reassign <- data.frame(
  Sample_ID = c("H28_NTC_R1", "H28_NTC_R2", "H28_NTC_R3", "H2052_NTC_R1", "H2052_NTC_R2"),
  Cell_Line_original_asnamed = c("NCI-H28", "NCI-H28", "NCI-H28", "NCI-H2052", "NCI-H2052"),
  Cell_Line_corrected = c("NCI-H2052", "NCI-H2052", "NCI-H2052", "NCI-H28", "NCI-H28"),
  registered_in = c(rep("H28_NTC_NATIVE (repro_common.R)", 3), rep("sample_metadata_RUN1.tsv (label_corrected)", 2)),
  evidence = "PCA clustering, fig35 (spec25_bulk_PCA_prefilter.R)", stringsAsFactors = FALSE)
write.csv(ntc_reassign, file.path(OUT, "spec25_NTC_reassignment_table.csv"), row.names = FALSE)
cat("\n=== registered NTC reassignment table ===\n"); print(ntc_reassign, row.names = FALSE)
cat("\nNTC per CORRECTED Cell_Line after both reassignments:\n"); print(table(meta_all$Cell_Line[meta_all$Gene_Control == "NTC"]))

## ---------------------------------------------------------------- fig36: per-line PCA with the FULL correction applied (both NTC reassignments), each line's own top 5000 variable genes, raw
NTOP2 <- 5000
pca_line4 <- function(cl) {
  s <- rownames(meta_all)[meta_all$Cell_Line == cl]     # fully corrected Cell_Line (both reassignments applied via H28_NTC_NATIVE + the sheet)
  ok <- rowMeans(cp10k_all[, s]) >= 0.1; y <- lg_all[ok, s, drop = FALSE]
  top <- names(sort(matrixStats::rowVars(y), decreasing = TRUE))[seq_len(min(NTOP2, nrow(y)))]
  yt <- y[top, , drop = FALSE]; yt <- yt[apply(yt, 1, sd) > 0, , drop = FALSE]; p <- prcomp(t(yt), center = TRUE, scale. = TRUE)
  list(p = p, ve = p$sdev^2 / sum(p$sdev^2), samples = s, n_expr = sum(ok))
}

save_fig("fig36_bulk_PCA_by_line_full_NTC_correction", 15, 4.6, function() {
  layout(matrix(1:4, nrow = 1)); par(oma = c(0, 0, 4.4, 0.5), mar = c(4, 4.2, 3.6, 1), family = "sans")
  for (cl in LINES) {
    r <- pca_line4(cl); X <- r$p$x[, 1:2]; s <- r$samples; reassigned <- s %in% ntc_reassign$Sample_ID
    rg <- apply(X, 2, range); pad <- 0.18 * apply(rg, 2, diff)
    plot(X[, 1], X[, 2], xlim = rg[, 1] + c(-1, 1) * pad[1], ylim = rg[, 2] + c(-1, 1) * pad[2], type = "n", bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8,
         xlab = sprintf("PC1 (%.0f%%)", 100 * r$ve[1]), ylab = sprintf("PC2 (%.0f%%)", 100 * r$ve[2])); abline(h = 0, v = 0, col = GRID)
    points(X[!reassigned, 1], X[!reassigned, 2], pch = 21, bg = adjustcolor(LINE_COL[[cl]], 0.85), col = "white", cex = 1.9)
    if (any(reassigned)) points(X[reassigned, 1], X[reassigned, 2], pch = 24, bg = adjustcolor(LINE_COL[[cl]], 0.9), col = INK, cex = 2.3, lwd = 1.6)
    text(X[, 1], X[, 2], paste0(meta_all[s, "Gene_Control"], " ", meta_all[s, "rep"], ifelse(reassigned, " (reassigned)", "")), pos = 3, cex = 0.5, offset = 0.4, font = ifelse(reassigned, 2, 1))
    mtext(sprintf("%s (corrected)", cl), side = 3, adj = 0, line = 1.9, font = 2, cex = 0.85, col = LINE_COL[[cl]])
    mtext(sprintf("top %d variable genes (of %d expressed), %d samples (%d reassigned)", length(r$p$center), r$n_expr, length(s), sum(reassigned)), side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2)
  }
  mtext("PCA of bulk RNA-seq samples per cell line, with BOTH NTC reassignments registered (see spec25_NTC_reassignment_table.csv)", outer = TRUE, side = 3, adj = 0, line = 2.8, font = 2, cex = 1.0)
  mtext("every KO except SNAI2/RUNX1/RUNX2; triangle outline = reassigned sample (H28_NTC_R1/2/3 -> NCI-H2052; H2052_NTC_R1/2 -> NCI-H28); each line's own top 5000 variable genes, raw", outer = TRUE, side = 3, adj = 0, line = 1.4, cex = 0.58, col = INK2)
})

## ---------------------------------------------------------------- fig37/38/39 shared setup: single Seurat object over ALL samples -> cell-cycle index (no per-line name-prefixing to worry about)
suppressPackageStartupMessages(library(Seurat))
obj <- CreateSeuratObject(counts = cs[, rownames(meta_all)], min.cells = 0, min.features = 0)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- CellCycleScoring(obj, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes, set.ident = FALSE)
cc_full <- setNames(obj$S.Score + obj$G2M.Score, colnames(obj))

pca_line5 <- function(cl) {
  s <- rownames(meta_all)[meta_all$Cell_Line == cl]
  ok <- rowMeans(cp10k_all[, s]) >= 0.1; y <- lg_all[ok, s, drop = FALSE]
  cc <- cc_full[s]; ccc <- cc - mean(cc); bb <- as.numeric((y - rowMeans(y)) %*% ccc) / sum(ccc^2); y <- y - outer(bb, ccc)
  top <- names(sort(matrixStats::rowVars(y), decreasing = TRUE))[seq_len(min(NTOP2, nrow(y)))]
  yt <- y[top, , drop = FALSE]; yt <- yt[apply(yt, 1, sd) > 0, , drop = FALSE]; p <- prcomp(t(yt), center = TRUE, scale. = TRUE)
  list(p = p, ve = p$sdev^2 / sum(p$sdev^2), samples = s, n_expr = sum(ok))
}

save_fig("fig37_bulk_PCA_by_line_full_NTC_correction_ccregressed", 15, 4.6, function() {
  layout(matrix(1:4, nrow = 1)); par(oma = c(0, 0, 4.4, 0.5), mar = c(4, 4.2, 3.6, 1), family = "sans")
  n_tot <- 0
  for (cl in LINES) {
    r <- pca_line5(cl); X <- r$p$x[, 1:2]; s <- r$samples; n_tot <- n_tot + length(s)
    rg <- apply(X, 2, range); pad <- 0.18 * apply(rg, 2, diff)
    plot(X[, 1], X[, 2], xlim = rg[, 1] + c(-1, 1) * pad[1], ylim = rg[, 2] + c(-1, 1) * pad[2], type = "n", bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8,
         xlab = sprintf("PC1 (%.0f%%)", 100 * r$ve[1]), ylab = sprintf("PC2 (%.0f%%)", 100 * r$ve[2])); abline(h = 0, v = 0, col = GRID)
    points(X[, 1], X[, 2], pch = 21, bg = adjustcolor(LINE_COL[[cl]], 0.85), col = "white", cex = 1.9)
    text(X[, 1], X[, 2], paste(meta_all[s, "Gene_Control"], meta_all[s, "rep"]), pos = 3, cex = 0.5, offset = 0.4, col = INK2)
    mtext(cl, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.85, col = LINE_COL[[cl]])
    mtext(sprintf("top %d variable genes (of %d expressed), %d samples", length(r$p$center), r$n_expr, length(s)), side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2)
  }
  mtext("PCA of bulk RNA-seq samples per cell line, cell-cycle regressed -- coloured by cell line, no other classification", outer = TRUE, side = 3, adj = 0, line = 2.8, font = 2, cex = 1.0)
  mtext(sprintf("%d samples total (same audited sample set as fig31, the across-cell-line PCA); every KO except SNAI2/RUNX1/RUNX2; each line's own top 5000 variable genes, cell-cycle regressed; every point labelled with its own Gene_Control + replicate", n_tot), outer = TRUE, side = 3, adj = 0, line = 1.4, cex = 0.58, col = INK2)
})

## ---------------------------------------------------------------- fig38: same sample set/PCA as fig37, highlighting samples NOT present in RUN1 (status != OK)
save_fig("fig38_bulk_PCA_by_line_ccregressed_highlight_not_in_RUN1", 15, 4.6, function() {
  layout(matrix(1:4, nrow = 1)); par(oma = c(0, 0, 4.4, 0.5), mar = c(4, 4.2, 3.6, 1), family = "sans")
  n_tot <- 0; n_rm <- 0
  for (cl in LINES) {
    r <- pca_line5(cl); X <- r$p$x[, 1:2]; s <- r$samples; rm_i <- meta_all[s, "removed_run1"]; n_tot <- n_tot + length(s); n_rm <- n_rm + sum(rm_i)
    rg <- apply(X, 2, range); pad <- 0.18 * apply(rg, 2, diff)
    plot(X[, 1], X[, 2], xlim = rg[, 1] + c(-1, 1) * pad[1], ylim = rg[, 2] + c(-1, 1) * pad[2], type = "n", bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8,
         xlab = sprintf("PC1 (%.0f%%)", 100 * r$ve[1]), ylab = sprintf("PC2 (%.0f%%)", 100 * r$ve[2])); abline(h = 0, v = 0, col = GRID)
    points(X[!rm_i, 1], X[!rm_i, 2], pch = 21, bg = adjustcolor(LINE_COL[[cl]], 0.85), col = "white", cex = 1.9)
    if (any(rm_i)) points(X[rm_i, 1], X[rm_i, 2], pch = 24, bg = adjustcolor(LINE_COL[[cl]], 0.9), col = INK, cex = 2.3, lwd = 1.6)
    text(X[, 1], X[, 2], paste0(meta_all[s, "Gene_Control"], " ", meta_all[s, "rep"], ifelse(rm_i, " (not in RUN1)", "")), pos = 3, cex = 0.5, offset = 0.4, font = ifelse(rm_i, 2, 1))
    mtext(cl, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.85, col = LINE_COL[[cl]])
    mtext(sprintf("top %d variable genes (of %d expressed), %d samples (%d not in RUN1)", length(r$p$center), r$n_expr, length(s), sum(rm_i)), side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2)
  }
  mtext("PCA of bulk RNA-seq samples per cell line, cell-cycle regressed -- highlighting samples NOT present in RUN1", outer = TRUE, side = 3, adj = 0, line = 2.8, font = 2, cex = 1.0)
  mtext(sprintf("%d samples total, %d not in RUN1 (status != OK; same audited sample set as fig37); every KO except SNAI2/RUNX1/RUNX2; each line's own top 5000 variable genes, cell-cycle regressed; triangle outline = excluded from RUN1", n_tot, n_rm), outer = TRUE, side = 3, adj = 0, line = 1.4, cex = 0.58, col = INK2)
})

## ---------------------------------------------------------------- fig39: same sample set/PCA as fig37/38, highlighting the VISUAL_OUTLIERS (repro_common.R) -- manually flagged from this same PCA: NCI-H2052 samples with PC2 > 20, plus specific named samples in the other 3 lines
save_fig("fig39_bulk_PCA_by_line_ccregressed_highlight_visual_outliers", 15, 4.6, function() {
  layout(matrix(1:4, nrow = 1)); par(oma = c(0, 0, 4.4, 0.5), mar = c(4, 4.2, 3.6, 1), family = "sans")
  n_tot <- 0; n_out <- 0
  for (cl in LINES) {
    r <- pca_line5(cl); X <- r$p$x[, 1:2]; s <- r$samples; out_i <- meta_all[s, "Sample_ID"] %in% VISUAL_OUTLIERS; n_tot <- n_tot + length(s); n_out <- n_out + sum(out_i)
    rg <- apply(X, 2, range); pad <- 0.18 * apply(rg, 2, diff)
    plot(X[, 1], X[, 2], xlim = rg[, 1] + c(-1, 1) * pad[1], ylim = rg[, 2] + c(-1, 1) * pad[2], type = "n", bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8,
         xlab = sprintf("PC1 (%.0f%%)", 100 * r$ve[1]), ylab = sprintf("PC2 (%.0f%%)", 100 * r$ve[2])); abline(h = 0, v = 0, col = GRID)
    points(X[!out_i, 1], X[!out_i, 2], pch = 21, bg = adjustcolor(LINE_COL[[cl]], 0.85), col = "white", cex = 1.9)
    if (any(out_i)) points(X[out_i, 1], X[out_i, 2], pch = 24, bg = adjustcolor(LINE_COL[[cl]], 0.9), col = INK, cex = 2.3, lwd = 1.6)
    text(X[, 1], X[, 2], paste0(meta_all[s, "Sample_ID"], ifelse(out_i, " (excluded)", "")), pos = 3, cex = 0.5, offset = 0.4, font = ifelse(out_i, 2, 1))
    mtext(cl, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.85, col = LINE_COL[[cl]])
    mtext(sprintf("top %d variable genes (of %d expressed), %d samples (%d excluded)", length(r$p$center), r$n_expr, length(s), sum(out_i)), side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2)
  }
  mtext("PCA of bulk RNA-seq samples per cell line, cell-cycle regressed -- highlighting the visually-identified outlier samples", outer = TRUE, side = 3, adj = 0, line = 2.8, font = 2, cex = 1.0)
  mtext(sprintf("%d samples total, %d flagged (NCI-H2052: PC2 > 20; MSTO_TWIST1_R2, MSTO_SOX9_R3; H2452_BPTF_R1, H2452_TWIST1_R2; H28_TCF3_R1, H28_HMGA1_R2); each line's own top 5000 variable genes, cell-cycle regressed; triangle outline = flagged for exclusion", n_tot, n_out), outer = TRUE, side = 3, adj = 0, line = 1.4, cex = 0.58, col = INK2)
})

## ---------------------------------------------------------------- fig40: same sample set/PCA as fig37/38/39, highlighting which of the 76 are part of the 70-sample set used by spec22's KO-reproducibility analysis (load_bulk(), screen KOSET)
b70 <- load_bulk(); ids70 <- b70$meta$Sample_ID[b70$meta$Gene_Control %in% c("NTC", ORDER)]
cat(sprintf("\nfig40: %d/%d full-audit samples are part of spec22's 70-sample KO-reproducibility set; %d are not:\n", sum(meta_all$Sample_ID %in% ids70), nrow(meta_all), sum(!meta_all$Sample_ID %in% ids70)))
print(meta_all$Sample_ID[!meta_all$Sample_ID %in% ids70])

save_fig("fig40_bulk_PCA_by_line_ccregressed_highlight_in_70_KOrepro", 15, 4.6, function() {
  layout(matrix(1:4, nrow = 1)); par(oma = c(0, 0, 4.4, 0.5), mar = c(4, 4.2, 3.6, 1), family = "sans")
  n_tot <- 0; n_in <- 0
  for (cl in LINES) {
    r <- pca_line5(cl); X <- r$p$x[, 1:2]; s <- r$samples; in70_i <- meta_all[s, "Sample_ID"] %in% ids70; n_tot <- n_tot + length(s); n_in <- n_in + sum(in70_i)
    rg <- apply(X, 2, range); pad <- 0.18 * apply(rg, 2, diff)
    plot(X[, 1], X[, 2], xlim = rg[, 1] + c(-1, 1) * pad[1], ylim = rg[, 2] + c(-1, 1) * pad[2], type = "n", bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8,
         xlab = sprintf("PC1 (%.0f%%)", 100 * r$ve[1]), ylab = sprintf("PC2 (%.0f%%)", 100 * r$ve[2])); abline(h = 0, v = 0, col = GRID)
    points(X[!in70_i, 1], X[!in70_i, 2], pch = 21, bg = adjustcolor(LINE_COL[[cl]], 0.85), col = "white", cex = 1.9)
    if (any(in70_i)) points(X[in70_i, 1], X[in70_i, 2], pch = 24, bg = adjustcolor(LINE_COL[[cl]], 0.9), col = INK, cex = 2.3, lwd = 1.6)
    text(X[, 1], X[, 2], paste0(meta_all[s, "Gene_Control"], " ", meta_all[s, "rep"], ifelse(in70_i, "", " (not in 70)")), pos = 3, cex = 0.5, offset = 0.4, font = ifelse(in70_i, 2, 1))
    mtext(cl, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.85, col = LINE_COL[[cl]])
    mtext(sprintf("top %d variable genes (of %d expressed), %d samples (%d in the 70)", length(r$p$center), r$n_expr, length(s), sum(in70_i)), side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2)
  }
  mtext("PCA of bulk RNA-seq samples per cell line, cell-cycle regressed -- highlighting the 70 samples used by spec22's KO-reproducibility analysis (fig26)", outer = TRUE, side = 3, adj = 0, line = 2.8, font = 2, cex = 1.0)
  mtext(sprintf("%d samples total (full audited universe), %d in the fig26/spec22 70-sample set (load_bulk(), NTC + the 8 CROP-seq KOs); triangle outline = in the 70; plain circle = the 6 orphaned/undocumented samples only present here (all also VISUAL_OUTLIERS, fig39)", n_tot, n_in), outer = TRUE, side = 3, adj = 0, line = 1.4, cex = 0.58, col = INK2)
})

## ---------------------------------------------------------------- fig41: same as fig40 but (a) triangle = EXCLUDED (not in the 70), plain circle = in the 70, and (b) the raw (non-cc-regressed) per-line PCA (pca_line4, same as fig36) instead of the cc-regressed one
save_fig("fig41_bulk_PCA_by_line_raw_highlight_not_in_70_KOrepro", 15, 4.6, function() {
  layout(matrix(1:4, nrow = 1)); par(oma = c(0, 0, 4.4, 0.5), mar = c(4, 4.2, 3.6, 1), family = "sans")
  n_tot <- 0; n_out <- 0
  for (cl in LINES) {
    r <- pca_line4(cl); X <- r$p$x[, 1:2]; s <- r$samples; notin70_i <- !meta_all[s, "Sample_ID"] %in% ids70; n_tot <- n_tot + length(s); n_out <- n_out + sum(notin70_i)
    rg <- apply(X, 2, range); pad <- 0.18 * apply(rg, 2, diff)
    plot(X[, 1], X[, 2], xlim = rg[, 1] + c(-1, 1) * pad[1], ylim = rg[, 2] + c(-1, 1) * pad[2], type = "n", bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8,
         xlab = sprintf("PC1 (%.0f%%)", 100 * r$ve[1]), ylab = sprintf("PC2 (%.0f%%)", 100 * r$ve[2])); abline(h = 0, v = 0, col = GRID)
    points(X[!notin70_i, 1], X[!notin70_i, 2], pch = 21, bg = adjustcolor(LINE_COL[[cl]], 0.85), col = "white", cex = 1.9)
    if (any(notin70_i)) points(X[notin70_i, 1], X[notin70_i, 2], pch = 24, bg = adjustcolor(LINE_COL[[cl]], 0.9), col = INK, cex = 2.3, lwd = 1.6)
    text(X[, 1], X[, 2], paste0(meta_all[s, "Gene_Control"], " ", meta_all[s, "rep"], ifelse(notin70_i, " (not in 70)", "")), pos = 3, cex = 0.5, offset = 0.4, font = ifelse(notin70_i, 2, 1))
    mtext(cl, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.85, col = LINE_COL[[cl]])
    mtext(sprintf("top %d variable genes (of %d expressed), %d samples (%d not in the 70)", length(r$p$center), r$n_expr, length(s), sum(notin70_i)), side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2)
  }
  mtext("PCA of bulk RNA-seq samples per cell line, raw (not cc-regressed) -- highlighting samples NOT in spec22's 70-sample KO-reproducibility set", outer = TRUE, side = 3, adj = 0, line = 2.8, font = 2, cex = 1.0)
  mtext(sprintf("%d samples total (full audited universe), %d NOT in the fig26/spec22 70-sample set; triangle outline = excluded (the 6 orphaned/undocumented samples, all also VISUAL_OUTLIERS, fig39); plain circle = in the 70", n_tot, n_out), outer = TRUE, side = 3, adj = 0, line = 1.4, cex = 0.58, col = INK2)
})
