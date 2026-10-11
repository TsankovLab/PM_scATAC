# Spec 4: are KO effects conserved between the bulk cell lines themselves? Same specificity test, but one bulk line plays the role of the "screen".
# If the same KO does not resemble itself across lines, KO effects are context (cell line) specific and a mismatch with the screen is expected.
# needs: spec1_prepare_stats.R
# outputs: spec_cross_line.csv
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
bs <- readRDS(file.path(OUT, "spec_bulk_stats.rds")); ot <- read.csv(file.path(OUT, "spec_ontarget_knockdown.csv"))
common <- Reduce(intersect, lapply(bs[LINES], `[[`, "genes"))
cat("genes expressed in all four lines:", length(common), "\n")
effective <- function(cl, ko) any(ot$cell_line == cl & ot$KO == ko & ot$lfc < -0.5 & ot$padj < 0.05)
res <- list()
for (what in c("stat0", "statcc")) for (centre in c(FALSE, TRUE)) for (a in LINES) for (bl in setdiff(LINES, a)) {
  S <- bs[[bl]][[what]][common, ]; A <- bs[[a]][[what]][common, ]
  cm <- cross_matrix(S, setNames(list(A), a), centre, FALSE)
  cols <- intersect(ORDER, colnames(S)); cm$C <- cm$C[, cols, drop = FALSE]
  both <- sapply(cm$lab$ko, function(k) effective(a, k) && effective(bl, k))
  for (rs in c("all", "both_effective")) {
    sel <- if (rs == "all") rep(TRUE, nrow(cm$C)) else both; sel <- sel & cm$lab$ko %in% cols; if (sum(sel) < 3) next
    res[[length(res) + 1]] <- data.frame(stat = what, centred = centre, rows_line = a, reference_line = bl, row_set = rs, evaluate_specificity(cm$C[sel, , drop = FALSE], cm$lab$ko[sel], cols = cols))
  }
}
res <- do.call(rbind, res); write.csv(res, file.path(OUT, "spec_cross_line.csv"), row.names = FALSE)
options(width = 220)
s <- aggregate(cbind(z_matched, matched_r, other_r, top1 / n_rows, perm_p < 0.05) ~ stat + centred + row_set, res, mean)
names(s) <- c("stat", "centred", "row_set", "mean_z_matched", "mean_matched_r", "mean_other_r", "mean_top1_frac", "frac_pairs_p<0.05"); print(format(s, digits = 3), row.names = FALSE)
cat("\nper line pair (t-stat, not cc-regressed, centred, all KOs):\n"); print(format(res[res$stat == "stat0" & res$centred & res$row_set == "all", c("rows_line", "reference_line", "n_rows", "matched_r", "other_r", "z_matched", "perm_p", "top1", "median_rank")], digits = 3), row.names = FALSE)
