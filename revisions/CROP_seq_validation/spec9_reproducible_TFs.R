# Spec 9: drop TFs whose own KO signature is not reproducible, then redraw the best-method bulk-vs-CROP-seq matrix.
# Reproducibility is judged WITHOUT looking at the bulk-vs-screen match (no selection on the outcome), on the same cc-adjusted, centred profiles as variant F4:
#   bulk:   for each line, replicate 1 vs replicate 2 of the same KO (per-sample log2 CP10K minus NTC, cell-cycle regressed, centred across all KO samples);
#           the KO's own replicate must rank 1st or 2nd among the KOs of that line ("hit"), in at least half of the lines where it has 2 replicates
#   screen: random halves of the cells (100 splits, cc regressed, centred across KOs); the KO's own half must rank 1st in >= 50% of splits
# A TF is kept only if it passes both. Rule fixed before running.
# needs: spec1_prepare_stats.R
# outputs: spec_TF_reproducibility.csv, spec_filtered_matrix.csv, fig11a_TF_reproducibility, fig11b_reproducible_TF_matrix
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_variants.R"))
set.seed(11)
ctx <- load_spec_context(LINES); ss <- ctx$ss; bs <- ctx$bs
genes <- Reduce(intersect, c(list(ss$genes), lapply(bs[LINES], `[[`, "genes")))

## ---------------------------------------------------------------- bulk replicate reproducibility
bulk_rows <- list()
for (cl in LINES) {
  m <- ctx$meta[ctx$meta$Cell_Line == cl, ]; s <- rownames(m); y <- ctx$lg[genes, s, drop = FALSE]; cc <- ctx$ccb$cc[match(s, ctx$ccb$sample)]
  bco <- as.numeric((y - rowMeans(y)) %*% (cc - mean(cc))) / sum((cc - mean(cc))^2); y <- y - outer(bco, cc - mean(cc))
  ntc <- rowMeans(y[, m$Gene_Control == "NTC", drop = FALSE]); ko_s <- s[m$Gene_Control != "NTC"]; D <- y[, ko_s, drop = FALSE] - ntc; D <- D - rowMeans(D)
  gc <- setNames(m[ko_s, "Gene_Control"], ko_s); kos <- sort(names(which(table(gc) == 2)), method = "radix")
  A <- sapply(kos, function(k) D[, sort(ko_s[gc == k])[1]]); B <- sapply(kos, function(k) D[, sort(ko_s[gc == k])[2]]); C <- cor(A, B)
  bulk_rows[[cl]] <- data.frame(line = cl, TF = kos, r = diag(C), rank = sapply(seq_along(kos), function(i) sum(C[i, ] > C[i, i]) + 1), n_kos = length(kos),
                                z = sapply(seq_along(kos), function(i) (C[i, i] - mean(C[i, ])) / sd(C[i, ])))
}
br <- do.call(rbind, bulk_rows); br$hit <- br$rank <= 2

## ---------------------------------------------------------------- screen split-half reproducibility
srt <- load_screen(); Y <- GetAssayData(srt, layer = "data")[genes, ]; grp <- srt$merged_call; ccs <- srt$cc
bco <- (as.numeric(Y %*% (ccs - mean(ccs))) / (ncol(Y) - 1)) / var(ccs); kos_all <- sort(setdiff(unique(grp), "NTC"), method = "radix")
acc <- lapply(1:100, function(i) {
  h <- ave(seq_along(grp), grp, FUN = function(x) sample(rep(1:2, length.out = length(x)))); key <- paste(grp, h); lev <- unique(key)
  G <- sparseMatrix(i = seq_along(key), j = match(key, lev), x = 1, dims = c(length(key), length(lev)), dimnames = list(NULL, lev)); cnt <- colSums(G)
  M <- sweep(as.matrix(Y %*% G), 2, cnt, "/"); mc <- as.numeric(crossprod(G, ccs)) / cnt; Mr <- M - outer(bco, mc)
  A <- Mr[, paste(kos_all, 1)] - Mr[, "NTC 1"]; B <- Mr[, paste(kos_all, 2)] - Mr[, "NTC 2"]; A <- A - rowMeans(A); B <- B - rowMeans(B); C <- cor(A, B)
  data.frame(TF = kos_all, r = diag(C), top1 = sapply(seq_along(kos_all), function(i) sum(C[i, ] > C[i, i]) == 0), z = sapply(seq_along(kos_all), function(i) (C[i, i] - mean(C[i, ])) / sd(C[i, ])))
})
sr <- aggregate(cbind(r, top1, z) ~ TF, do.call(rbind, acc), mean)

## ---------------------------------------------------------------- combine, apply rule
tab <- do.call(rbind, lapply(ORDER, function(k) { b <- br[br$TF == k, ]; s <- sr[sr$TF == k, ]
  data.frame(TF = k, screen_cells = ss$n[[k]], bulk_lines_with_2reps = nrow(b), bulk_hit_lines = sum(b$hit), bulk_hit_frac = if (nrow(b)) mean(b$hit) else NA, bulk_mean_r = mean(b$r), bulk_mean_z = mean(b$z),
             screen_top1_frac = s$top1, screen_mean_r = s$r, screen_mean_z = s$z) }))
tab$bulk_pass <- tab$bulk_hit_frac >= 0.5; tab$screen_pass <- tab$screen_top1_frac >= 0.5; tab$keep <- tab$bulk_pass & tab$screen_pass
write.csv(tab, file.path(OUT, "spec_TF_reproducibility.csv"), row.names = FALSE); options(width = 220)
cat("\n=== per-line bulk replicate reproducibility (rank of own replicate among the line's KOs) ===\n"); print(format(br, digits = 3), row.names = FALSE)
cat("\n=== TF reproducibility and rule ===\n"); print(format(tab, digits = 3), row.names = FALSE)
KEEP <- tab$TF[tab$keep]; cat("\nkept TFs:", paste(KEEP, collapse = ", "), "| removed:", paste(setdiff(ORDER, KEEP), collapse = ", "), "\n")

## ---------------------------------------------------------------- F4 matrix restricted to reproducible TFs
COLS <- KEEP; rows_lines <- c("NCI-H2052", "MSTO-211H", "NCI-H2452", "NCI-H28")
S <- prep_dataset(ss$tcc[genes, COLS, drop = FALSE], TRUE, FALSE); S_all <- prep_dataset(ss$tcc[genes, c(ORDER, EXTRA)], TRUE, FALSE)
rows <- list(); lab <- list()
for (cl in rows_lines) { mm <- bs[[cl]]$statcc[genes, ]; mm[!is.finite(mm)] <- 0; mm <- prep_dataset(mm, TRUE, FALSE)      # centring across ALL KOs of the line and the screen, as in F4
  for (ko in COLS) if (ko %in% colnames(mm)) { rows[[length(rows) + 1]] <- sapply(COLS, function(cc) cor(mm[, ko], S_all[, cc])); lab[[length(lab) + 1]] <- c(cl, ko) } }
C <- do.call(rbind, rows); lab <- setNames(as.data.frame(do.call(rbind, lab), stringsAsFactors = FALSE), c("line", "ko")); rownames(C) <- paste(lab$line, lab$ko); write.csv(C, file.path(OUT, "spec_filtered_matrix.csv"))
sets <- list(`3 lines` = lab$line %in% THREE, `4 lines` = rep(TRUE, nrow(lab)), `NCI-H2052 only` = lab$line == "NCI-H2052")
st <- do.call(rbind, lapply(sets, function(sel) evaluate_specificity(C[sel, , drop = FALSE], lab$ko[sel], cols = COLS))); cat("\n=== best method restricted to reproducible TFs ===\n"); print(format(st, digits = 3))
rank_k <- sapply(seq_len(nrow(lab)), function(i) sum(C[i, ] > C[i, lab$ko[i]]) + 1)

## ---------------------------------------------------------------- figures
save_fig("fig11a_TF_reproducibility", 12, 4.6, function() {
  layout(matrix(1:2, nrow = 1), widths = c(1.35, 1)); par(oma = c(0, 0, 3.5, 0.5), family = "sans"); ord <- rev(ORDER); n <- length(ord)
  par(mar = c(3.5, 6, 3, 0.5)); plot(NA, xlim = c(-1.2, 3.3), ylim = c(0.5, n + 0.5), axes = FALSE, xlab = "", ylab = ""); abline(v = 0, col = INK2); abline(h = seq_len(n), col = GRID, lwd = 0.6)
  for (i in seq_len(n)) { b <- br[br$TF == ord[i], ]; if (nrow(b)) points(b$z, rep(i, nrow(b)) + seq(-0.15, 0.15, length.out = nrow(b)), pch = 21, bg = LINE_COL[b$line], col = ifelse(b$hit, INK, "white"), cex = 1.5, lwd = 1.2) }
  axis(1, cex.axis = 0.8, lwd = 0.5); kp <- ord %in% KEEP; axis(2, at = which(kp), labels = ord[kp], las = 1, tick = FALSE, cex.axis = 0.95, font.axis = 2, col.axis = INK); axis(2, at = which(!kp), labels = ord[!kp], las = 1, tick = FALSE, cex.axis = 0.95, font.axis = 1, col.axis = MUTED)
  mtext("replicate 1 vs 2 of the same KO: z vs other KOs of the line", side = 1, line = 2.3, cex = 0.7); mtext("Bulk reproducibility (black ring = own replicate ranks top 2)", side = 3, adj = 0, line = 0.7, font = 2, cex = 0.85)
  legend("topright", sub("NCI-", "", LINES), pt.bg = LINE_COL[LINES], pch = 21, col = "white", bty = "n", cex = 0.7, pt.cex = 1.3)
  par(mar = c(3.5, 0.5, 3, 0.5)); s <- tab$screen_top1_frac[match(ord, tab$TF)]; plot(NA, xlim = c(0, 1.6), ylim = c(0.5, n + 0.5), axes = FALSE, xlab = "", ylab = ""); abline(v = 0.5, lty = 2, col = INK2)
  rect(0, seq_len(n) - 0.32, s, seq_len(n) + 0.32, col = ifelse(s >= 0.5, "#8fa88a", "#c9c9c9"), border = NA); text(pmax(s, 0) + 0.02, seq_len(n), sprintf("%.0f%%  (n=%d cells)", 100 * s, ss$n[ord]), adj = 0, cex = 0.7, col = INK2)
  axis(1, at = c(0, 0.5, 1), labels = c("0", "50%", "100%"), cex.axis = 0.8, lwd = 0.5); mtext("share of random cell splits in which the KO matches its own other half", side = 1, line = 2.3, cex = 0.7); mtext("Screen reproducibility", side = 3, adj = 0, line = 0.7, font = 2, cex = 0.85)
  mtext("Which TF knockouts have a reproducible signature?", outer = TRUE, side = 3, adj = 0, line = 1.7, font = 2, cex = 1)
  mtext(sprintf("kept (bold): %s | removed: %s. Rule: own replicate in the top 2 in >= 50%% of lines AND own screen half first in >= 50%% of splits", paste(KEEP, collapse = ", "), paste(setdiff(ORDER, KEEP), collapse = ", ")), outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.6, col = INK2)
})
n <- nrow(C); k <- ncol(C); vmax <- ceiling(quantile(abs(C), 0.99, names = FALSE) * 20) / 20
save_fig("fig11b_reproducible_TF_matrix", 8.2, 9.5, function() {
  par(mar = c(5, 9, 12, 6.5), family = "sans"); image(seq_len(k), seq_len(n), heat_z(pmax(pmin(C, vmax), -vmax)), col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
  axis(3, at = seq_len(k), labels = sprintf("%s (%d)", COLS, ss$n[COLS]), las = 2, tick = FALSE, cex.axis = 0.95, line = -0.4, font.axis = 2); axis(2, at = n:1, labels = lab$ko, las = 1, tick = FALSE, cex.axis = 0.85, line = -0.4)
  for (i in seq_len(n)) for (j in seq_len(k)) text(j, n - i + 1, sub("^([+-])0\\.", "\\1.", sprintf("%+.2f", C[i, j])), cex = 0.75, col = if (abs(C[i, j]) < 0.6 * vmax) INK2 else "white")
  for (i in seq_len(n)) { j <- match(lab$ko[i], COLS); rect(j - 0.5, n - i + 0.5, j + 0.5, n - i + 1.5, border = INK, lwd = 2 ) }
  edges <- which(lab$line[-1] != lab$line[-n]); abline(h = n - edges + 0.5, col = "white", lwd = 5); bounds <- c(0, edges, n); par(xpd = NA)
  for (b in seq_along(rows_lines)) { y1 <- n - bounds[b + 1] + 0.5; y2 <- n - bounds[b] + 0.5; rect(-0.85, y1, -0.6, y2, col = LINE_COL[[rows_lines[b]]], border = NA); text(-1.25, (y1 + y2) / 2, rows_lines[b], srt = 90, font = 2, cex = 0.9, col = LINE_COL[[rows_lines[b]]]) }
  text(k + 0.7, n:1, sprintf("rank %d/%d", rank_k, k), adj = 0, cex = 0.75, col = ifelse(rank_k <= 1, "#b3400f", MUTED), font = ifelse(rank_k <= 1, 2, 1))
  at <- -1.35
  mtext(sprintf("Best method, reproducible TFs only (%s)", paste(KEEP, collapse = ", ")), side = 3, adj = 0, line = 10.4, font = 2, cex = 0.95, at = at)
  mtext("cell-cycle adjusted test statistics, all expressed genes, centred across KOs; rows = bulk KO per line, columns = CROP-seq KO (cells)", side = 3, adj = 0, line = 9.2, cex = 0.58, col = INK2, at = at)
  mtext(sprintf("exact permutation test: 3 lines z = %.2f (p = %.3f) | 4 lines z = %.2f (p = %.3f) | NCI-H2052 only z = %.2f (p = %.3f)", st$z_matched[1], st$perm_p[1], st$z_matched[2], st$perm_p[2], st$z_matched[3], st$perm_p[3]), side = 3, adj = 0, line = 8.2, cex = 0.58, col = INK2, at = at)
  mtext(sprintf("matched KO ranks 1st in %d of %d rows (chance %.1f); mean r matched %+.3f vs other KOs %+.3f (3 lines)", st$top1[2], st$n_rows[2], st$n_rows[2] / k, st$matched_r[1], st$other_r[1]), side = 3, adj = 0, line = 7.2, cex = 0.58, col = INK2, at = at)
  par(xpd = FALSE); draw_colorbar(c(0.3, 0.7), c(-vmax, vmax), "Pearson r")
})
