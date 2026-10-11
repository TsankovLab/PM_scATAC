# Spec 7: which (cell line, TF) pairs reproduce the CROP-seq screen consistently across methods, with selection-aware inference.
# For a family of methods, every bulk row (line x TF) gets a consistency score = fraction of methods in which its matched screen KO ranks in the top 2 of the 6
# well-powered screen KOs (>= 50 cells; chance = 1/3). The same statistics are recomputed for all 720 relabelings of the screen KOs, so p-values include the
# selection of the best rows / TFs / lines (family-wise). Two families: cell-cycle regressed variants and variants without regression.
# needs: spec1_prepare_stats.R
# outputs: spec_best_set.csv, spec_best_set_summary.csv, fig9_best_set.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_variants.R"))

ctx <- load_spec_context(LINES); ss <- ctx$ss
V <- c(make_core_variants(ctx), make_abs_variants(ctx))
POWERED <- intersect(ORDER, names(ss$n)[ss$n >= 50]); P <- get_perms(length(POWERED))
fam <- list(cc_regressed = c("F1", "F2", "F3", "F4", "G3", "G4", "H1", "H2", "N05", "N06", "N07", "N08", "N13", "N14", "N15", "N16"),
            no_cc = c("A1", "A2", "B1", "B2", "C1", "D1", "D2", "D3", "D4", "E1", "E2", "G1", "G2", "H0", "N01", "N02", "N03", "N04", "N09", "N10", "N11", "N12"))
mats <- lapply(unlist(fam), function(id) { cm <- V[[id]]$fun(); keep <- cm$lab$ko %in% POWERED; list(C = cm$C[keep, POWERED, drop = FALSE], lab = cm$lab[keep, ]) }); names(mats) <- unlist(fam)
lab <- mats[[1]]$lab; stopifnot(all(sapply(mats, function(m) identical(m$lab, lab))))
n <- nrow(lab); j <- match(lab$ko, POWERED); cat("rows:", n, "(", paste(unique(lab$line), collapse = ", "), ") | screen KOs:", paste(POWERED, collapse = ", "), "\n")

# top-2 indicator per method x row x permutation (perm 0 = observed)
top2 <- function(ids) {
  arr <- array(NA, c(length(ids), n, nrow(P))); for (m in seq_along(ids)) { C <- mats[[ids[m]]]$C
    for (i in seq_len(n)) { matched <- C[i, P[, j[i]]]; arr[m, i, ] <- rowSums(sapply(seq_len(ncol(C)), function(k) C[i, k] > matched)) + 1 <= 2 } }
  arr }
ident <- which(apply(P, 1, function(p) all(p == seq_along(p))))
summ <- list(); rows_out <- list()
for (fn in names(fam)) {
  arr <- top2(fam[[fn]]); cons <- apply(arr, c(2, 3), mean)               # rows x perms
  obs <- cons[, ident]
  tf_of <- lab$ko; ln_of <- lab$line; tfs <- POWERED; lines <- unique(ln_of)
  stat_row <- function(cv) c(mean_all = mean(cv), n_ge75 = sum(cv >= 0.75), max_row = max(cv),
                              max_TF = max(sapply(tfs, function(t) mean(cv[tf_of == t]))), max_line = max(sapply(lines, function(l) mean(cv[ln_of == l]))))
  null <- t(apply(cons, 2, stat_row)); ob <- stat_row(obs)
  pv <- sapply(names(ob), function(s) mean(null[, s] >= ob[s]))
  summ[[fn]] <- data.frame(family = fn, n_methods = length(fam[[fn]]), statistic = names(ob), observed = as.numeric(ob), null_mean = colMeans(null), perm_p = as.numeric(pv))
  rp <- sapply(seq_len(n), function(i) mean(cons[i, ] >= obs[i])); rp_fw <- sapply(seq_len(n), function(i) mean(null[, "max_row"] >= obs[i]))
  rows_out[[fn]] <- data.frame(family = fn, line = ln_of, TF = tf_of, consistency = obs, chance = 1 / 3, row_perm_p = rp, familywise_p = rp_fw,
                               mean_matched_rank = sapply(seq_len(n), function(i) mean(sapply(fam[[fn]], function(id) { C <- mats[[id]]$C; sum(C[i, ] > C[i, j[i]]) + 1 }))))
}
summ <- do.call(rbind, summ); rows_out <- do.call(rbind, rows_out); rownames(summ) <- NULL; rownames(rows_out) <- NULL
write.csv(summ, file.path(OUT, "spec_best_set_summary.csv"), row.names = FALSE); write.csv(rows_out, file.path(OUT, "spec_best_set.csv"), row.names = FALSE)
options(width = 200)
cat("\n=== selection-aware permutation tests (720 relabelings) ===\n"); print(format(summ, digits = 3), row.names = FALSE)
for (fn in names(fam)) { cat("\n=== (line, TF) consistency, family:", fn, "===\n"); d <- rows_out[rows_out$family == fn, ]
  print(round(tapply(d$consistency, list(d$TF, d$line), mean), 2)); cat("mean over lines:\n"); print(round(sort(tapply(d$consistency, d$TF, mean), decreasing = TRUE), 2)); cat("mean over TFs:\n"); print(round(sort(tapply(d$consistency, d$line, mean), decreasing = TRUE), 2)) }
cat("\nrows with consistency >= 0.75 (cc_regressed):\n"); print(format(rows_out[rows_out$family == "cc_regressed" & rows_out$consistency >= 0.75, ], digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- figure
save_fig("fig9_best_set", 12.5, 5.6, function() {
  layout(matrix(1:2, nrow = 1)); par(oma = c(0, 0, 4, 0.5), family = "sans")
  for (fn in names(fam)) {
    d <- rows_out[rows_out$family == fn, ]; tfs <- POWERED; lines <- LINES
    M <- sapply(lines, function(l) sapply(tfs, function(t) { v <- d$consistency[d$line == l & d$TF == t]; if (length(v)) v else NA }))
    par(mar = c(2, 5, 6, 1)); image(seq_along(lines), seq_along(tfs), t(M)[, length(tfs):1], col = colorRampPalette(c("#f0efec", "#eb6834"))(50), zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_along(tfs)) for (k in seq_along(lines)) { v <- M[i, k]; text(k, length(tfs) - i + 1, sprintf("%.2f", v), cex = 0.85, font = ifelse(!is.na(v) && v >= 0.75, 2, 1)) }
    axis(2, at = length(tfs):1, labels = tfs, las = 1, tick = FALSE, cex.axis = 0.9, line = -0.5); mtext(sub("NCI-", "", lines), side = 3, at = seq_along(lines), cex = 0.85, line = 0.4)
    s <- summ[summ$family == fn, ]
    mtext(sprintf("%s family (%d methods)", ifelse(fn == "cc_regressed", "cell-cycle regressed", "no regression"), length(fam[[fn]])), side = 3, adj = 0, line = 3.5, font = 2, cex = 0.85)
    mtext(sprintf("mean consistency %.2f (chance 0.33; permutation p = %.3f) | rows >= 0.75: %d (p = %.3f)", s$observed[s$statistic == "mean_all"], s$perm_p[s$statistic == "mean_all"], s$observed[s$statistic == "n_ge75"], s$perm_p[s$statistic == "n_ge75"]), side = 3, adj = 0, line = 2.2, cex = 0.6, col = INK2)
    mtext(sprintf("best TF (mean over lines) p = %.3f | best line p = %.3f", s$perm_p[s$statistic == "max_TF"], s$perm_p[s$statistic == "max_line"]), side = 3, adj = 0, line = 1.2, cex = 0.6, col = INK2)
  }
  mtext("Which bulk (line, TF) pairs reproduce the CROP-seq screen consistently across methods?", outer = TRUE, side = 3, adj = 0, line = 2.5, font = 2, cex = 1)
  mtext("cell = fraction of methods in which the matched screen KO ranks in the top 2 of the 6 screen KOs with >= 50 cells; p-values from all 720 relabelings of the screen KOs (selection included)", outer = TRUE, side = 3, adj = 0, line = 1.2, cex = 0.6, col = INK2)
})
