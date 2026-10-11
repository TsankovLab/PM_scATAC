# Spec 3: does any processing make bulk KO profiles match their CROP-seq KO better than the other KOs?
# Variants: distance-based (KO-NTC log2 diff, RMS-scaled or not), the original CROP-seq heatmap processing (z-score across KO groups),
#           test statistics (screen Welch t vs bulk DESeq2 Wald), all expressed genes vs the screen's variable genes, Spearman,
#           cell-cycle regressed (cc = S.Score + G2M.Score, both datasets), and screen-defined KO signatures scored in bulk.
# Each variant gives a rows (bulk line x KO) x 8 screen-KO correlation/score matrix; specificity = matched entry after row z-scoring, exact permutation p over all 8! KO assignments.
# Row sets: ALL3 (3 lines), H2052 (the screen's own line), EFF (bulk on-target knockdown lfc < -0.5, padj < 0.05), EFF2 (EFF and the KO also knocked down in the screen, t < -2).
# needs: spec1_prepare_stats.R, screen_export_figs.R, bulk_cellcycle_RUN1.R
# outputs: spec_method_comparison.csv, fig6_specificity_methods.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_variants.R"))

ctx <- load_spec_context(THREE); ss <- ctx$ss
cat("features: variable genes", length(ctx$varg), "(in all stat sets:", length(ctx$F_var), ") | all expressed genes", length(ctx$F_all), "\n")
row_sets <- make_row_sets(ctx); V <- make_core_variants(ctx)
cat("EFF rows:", paste(paste(attr(row_sets, "eff")$cell_line, attr(row_sets, "eff")$KO), collapse = ", "), "\nEFF2 rows:", paste(paste(attr(row_sets, "eff2")$cell_line, attr(row_sets, "eff2")$KO), collapse = ", "), "\n")

# column sets: all 8 shared screen KOs, or only screen KOs with >= 50 cells (drops HMGA1 n=13 and TEAD2 n=20, whose profiles are noisy hubs)
POWERED <- intersect(ORDER, names(ss$n)[ss$n >= 50])
col_sets <- list(all8 = ORDER, powered = POWERED); cat("screen KOs with >= 50 cells:", paste(POWERED, collapse = ", "), "\n")
res <- evaluate_variants(V, row_sets, col_sets); write.csv(res, file.path(OUT, "spec_method_comparison.csv"), row.names = FALSE)
options(width = 220)
for (cs in names(col_sets)) for (rs in names(row_sets)) { cat("\n===== columns:", cs, "| rows:", rs, "=====\n"); r <- res[res$row_set == rs & res$col_set == cs, c("id", "n_rows", "matched_r", "other_r", "z_matched", "perm_p", "top1", "top3", "median_rank")]; print(format(r, digits = 3), row.names = FALSE) }
cat("\nSmallest permutation p per column/row set:\n"); print(aggregate(perm_p ~ col_set + row_set, res, min), row.names = FALSE)

## ---------------------------------------------------------------- figure: row-standardised matched score with exact permutation p
for (cs in names(col_sets)) save_fig(paste0("fig6_specificity_methods_", cs), 15, 8.5, function() {
  ids <- names(V); k <- length(ids); sets <- c("ALL3", "H2052", "EFF", "EFF2")
  layout(matrix(1:4, nrow = 1), widths = c(1.9, 1, 1, 1)); par(oma = c(3, 0, 4.5, 1), family = "sans")
  for (si in seq_along(sets)) {
    d <- res[res$row_set == sets[si] & res$col_set == cs, ]; par(mar = c(2.5, if (si == 1) 22 else 0.5, 2.6, 0.5))
    z <- d$z_matched[match(ids, d$id)]; p <- d$perm_p[match(ids, d$id)]; n <- d$n_rows[1]
    plot(NA, xlim = c(-1.1, 1.9), ylim = c(k + 0.5, 0.5), axes = FALSE, xlab = "", ylab = "")
    abline(v = 0, col = INK2); abline(h = seq_len(k), col = GRID, lwd = 0.6)
    rect(0, seq_len(k) - 0.35, ifelse(is.na(z), 0, z), seq_len(k) + 0.35, col = ifelse(!is.na(p) & p < 0.05, "#eb6834", "#9db7d6"), border = NA)
    text(ifelse(is.na(z), 0, pmax(z, 0)) + 0.03, seq_len(k), ifelse(is.na(p), "", ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p))), adj = 0, cex = 0.55, col = INK2)
    axis(1, cex.axis = 0.7, lwd = 0.5, mgp = c(1, 0.3, 0)); mtext(sprintf("%s (%d rows)", c("all 3 lines", "NCI-H2052 (screen line)", "effective bulk knockdown", "knockdown in bulk and screen")[si], n), side = 3, adj = 0, line = 0.6, font = 2, cex = 0.8)
    if (si == 1) axis(2, at = seq_len(k), labels = paste0(ids, "  ", sapply(V, `[[`, "desc")), las = 1, tick = FALSE, cex.axis = 0.62, line = -0.5)
  }
  mtext("Specificity of bulk KO vs CROP-seq KO correspondence across processing variants", outer = TRUE, side = 3, adj = 0, line = 2.6, font = 2, cex = 1.1)
  mtext(sprintf("columns = %s; bar = matched screen KO's score after standardising each bulk row across those screen KOs (0 = no better than the other KOs; ~1.4 = perfect); orange = exact permutation p < 0.05", if (cs == "all8") "all 8 shared screen KOs" else paste0("screen KOs with >= 50 cells (", paste(POWERED, collapse = ", "), ")")),
        outer = TRUE, side = 3, adj = 0, line = 1.3, cex = 0.6, col = INK2)
  mtext("row-standardised matched score (z)", outer = TRUE, side = 1, line = 1.2, cex = 0.7, col = INK2)
})
