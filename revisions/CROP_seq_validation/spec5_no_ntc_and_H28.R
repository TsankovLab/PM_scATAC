# Spec 5: (a) compare KO groups directly WITHOUT subtracting NTC (absolute expression; raw / centred across KOs / z-scored, with and without cell-cycle regression,
#             on the screen's variable genes or the 3000 most expressed genes) and (b) add NCI-H28 as a fourth bulk line (rows ALL4) and on its own (H28).
# The core KO-NTC variants of spec3 are re-run on the same 4-line gene set for reference.
# needs: spec1_prepare_stats.R
# outputs: spec_no_ntc_H28_comparison.csv, fig7_specificity_no_ntc_H28_{all8,powered}.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_variants.R"))

ctx <- load_spec_context(LINES); ss <- ctx$ss
cat("genes: variable", length(ctx$F_var), "| all expressed", length(ctx$F_all), "(common to the 4 bulk lines and the screen)\n")
row_sets <- make_row_sets(ctx)
V <- c(make_core_variants(ctx), make_abs_variants(ctx))
POWERED <- intersect(ORDER, names(ss$n)[ss$n >= 50]); col_sets <- list(all8 = ORDER, powered = POWERED)
res <- evaluate_variants(V, row_sets, col_sets); write.csv(res, file.path(OUT, "spec_no_ntc_H28_comparison.csv"), row.names = FALSE)
res$family <- ifelse(startsWith(res$id, "N"), "absolute (no NTC)", "KO - NTC")
options(width = 230)
cat("\n=== absolute (no NTC subtraction) variants ===\n")
for (cs in names(col_sets)) for (rs in c("ALL3", "ALL4", "H2052", "H28")) {
  cat("\n--- columns:", cs, "| rows:", rs, "---\n"); r <- res[res$family == "absolute (no NTC)" & res$col_set == cs & res$row_set == rs, c("id", "description", "n_rows", "matched_r", "other_r", "z_matched", "perm_p", "top1")]
  print(format(r, digits = 3), row.names = FALSE)
}
cat("\n=== mean z_matched by family and row set ===\n"); print(round(tapply(res$z_matched, list(paste(res$family, res$col_set), res$row_set), mean), 3))
cat("\ntests:", nrow(res), "| p < 0.05:", sum(res$perm_p < 0.05), "(expected by chance ~", round(0.05 * nrow(res), 1), ")\n")
cat("\ntests with p < 0.05:\n"); print(format(res[res$perm_p < 0.05, c("id", "description", "col_set", "row_set", "n_rows", "matched_r", "other_r", "z_matched", "perm_p", "top1")], digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- figure: one panel per row set, all variants
for (cs in names(col_sets)) save_fig(paste0("fig7_specificity_no_ntc_H28_", cs), 15, 12.5, function() {
  ids <- names(V); k <- length(ids); sets <- c("ALL3", "ALL4", "H2052", "H28")
  layout(matrix(1:4, nrow = 1), widths = c(2.2, 1, 1, 1)); par(oma = c(3, 0, 4.5, 1), family = "sans")
  for (si in seq_along(sets)) {
    d <- res[res$row_set == sets[si] & res$col_set == cs, ]; par(mar = c(2.5, if (si == 1) 27 else 0.5, 2.6, 0.5))
    z <- d$z_matched[match(ids, d$id)]; p <- d$perm_p[match(ids, d$id)]; n <- d$n_rows[1]
    plot(NA, xlim = c(-1.1, 1.9), ylim = c(k + 0.5, 0.5), axes = FALSE, xlab = "", ylab = "")
    abline(v = 0, col = INK2); abline(h = seq_len(k), col = GRID, lwd = 0.6); abline(h = sum(!startsWith(ids, "N")) + 0.5, col = INK2, lwd = 0.8, lty = 2)
    rect(0, seq_len(k) - 0.35, ifelse(is.na(z), 0, z), seq_len(k) + 0.35, col = ifelse(!is.na(p) & p < 0.05, "#eb6834", ifelse(startsWith(ids, "N"), "#8fa88a", "#9db7d6")), border = NA)
    text(ifelse(is.na(z), 0, pmax(z, 0)) + 0.03, seq_len(k), ifelse(is.na(p), "", ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p))), adj = 0, cex = 0.5, col = INK2)
    axis(1, cex.axis = 0.7, lwd = 0.5, mgp = c(1, 0.3, 0)); mtext(sprintf("%s (%d rows)", c("3 lines (no NCI-H28)", "all 4 lines", "NCI-H2052 (screen line)", "NCI-H28 only")[si], n), side = 3, adj = 0, line = 0.6, font = 2, cex = 0.8)
    if (si == 1) axis(2, at = seq_len(k), labels = paste0(ids, "  ", sapply(V, `[[`, "desc")), las = 1, tick = FALSE, cex.axis = 0.55, line = -0.5)
  }
  mtext("Bulk KO vs CROP-seq KO specificity: KO - NTC variants (blue) vs absolute KO profiles without NTC subtraction (green), incl. NCI-H28", outer = TRUE, side = 3, adj = 0, line = 2.6, font = 2, cex = 1)
  mtext(sprintf("columns = %s; bar = matched screen KO's row-standardised score (0 = no better than other KOs); orange = exact permutation p < 0.05", if (cs == "all8") "all 8 shared screen KOs" else paste0("screen KOs with >= 50 cells (", paste(POWERED, collapse = ", "), ")")),
        outer = TRUE, side = 3, adj = 0, line = 1.3, cex = 0.6, col = INK2)
  mtext("row-standardised matched score (z)", outer = TRUE, side = 1, line = 1.2, cex = 0.7, col = INK2)
})
