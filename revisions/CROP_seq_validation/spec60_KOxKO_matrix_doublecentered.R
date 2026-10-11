# Spec 60: two-way (double) centering of the spec58 KO x KO matrices. Each 8x8 Spearman matrix (rows = screen KO, columns = bulk H2052 KO, NTC - KO) has the row means and the
# column means removed (resid = C - rowMean - colMean + grandMean = the KO-by-KO INTERACTION part), which drops both the screen-KO offset (row effect) and the bulk-KO offset
# (column effect) that make all NTC combinations look alike. Diagonal metrics on the residual matrix: mean diagonal residual, top1 (# screen KOs whose max residual is on the
# diagonal), mean rank of the diagonal in its row, exact permutation p (all 8! relabelings of bulk KOs) for the mean diagonal residual.
# needs: spec58_H2052_NTC_combos_KOxKO_matrix.R (spec58_KOxKO_matrices.rds)
# outputs: spec60_doublecentered_diagonal_metrics.csv, fig3U_KOxKO_matrix_dc_{Cm,top1000,top2000,top3000}.{png,pdf}, fig3V_dc_diagonal_summary.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
mats <- readRDS(file.path(OUT, "spec58_KOxKO_matrices.rds")); spaces <- names(mats); sp_tag <- c("Cm", "top1000", "top2000", "top3000"); names(sp_tag) <- spaces
lab <- names(mats[[1]]); K <- nrow(mats[[1]][[1]]); KOS <- rownames(mats[[1]][[1]]); nn <- sapply(strsplit(lab, "+", fixed = TRUE), length)
perms <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(seq_len(K))) })
dc <- function(C) { R <- C - rowMeans(C)[row(C)] - colMeans(C)[col(C)] + mean(C); dimnames(R) <- dimnames(C); R }
met <- list(); Z <- list()
for (sp in spaces) for (s in lab) { C <- dc(mats[[sp]][[s]]); Z[[sp]][[s]] <- C; dm <- mean(diag(C)); pd <- apply(perms, 1, function(p) mean(C[cbind(seq_len(K), p)]))
  met[[length(met) + 1]] <- data.frame(space = sp, subset = s, n_ntc = length(strsplit(s, "+", fixed = TRUE)[[1]]), diag_mean_resid = dm, top1_rows = sum(apply(C, 1, which.max) == seq_len(K)),
       top1_cols = sum(apply(C, 2, which.max) == seq_len(K)), mean_row_rank = mean(sapply(seq_len(K), function(r) rank(-C[r, ])[r])), perm_p = mean(pd >= dm - 1e-12)) }
met <- do.call(rbind, met); write.csv(met, file.path(OUT, "spec60_doublecentered_diagonal_metrics.csv"), row.names = FALSE)
summ <- do.call(rbind, lapply(split(met, met$subset), function(d) data.frame(subset = d$subset[1], n_ntc = d$n_ntc[1], mean_diag_resid = mean(d$diag_mean_resid), mean_top1 = mean(d$top1_rows), mean_row_rank = mean(d$mean_row_rank), min_perm_p = min(d$perm_p),
     Cm = d$diag_mean_resid[d$space == spaces[1]], top1000 = d$diag_mean_resid[d$space == spaces[2]], top2000 = d$diag_mean_resid[d$space == spaces[3]], top3000 = d$diag_mean_resid[d$space == spaces[4]])))
summ <- summ[order(-summ$mean_diag_resid), ]; options(width = 230)
cat("\n=== double-centered (row+column effects removed): diagonal quality per H2052 NTC combination; chance: resid 0, top1 1/8, mean row rank 4.5 ===\n"); print(format(summ, digits = 2), row.names = FALSE)
cat("\n=== per space: best combination, its top1, row rank, permutation p; R3-only and all-4 for reference ===\n")
for (sp in spaces) { d <- met[met$space == sp, ]; b <- d[which.max(d$diag_mean_resid), ]; r3 <- d[d$subset == "H2052R3", ]; a4 <- d[d$n_ntc == 4, ]
  cat(sprintf("%-24s best %-24s resid %+.3f top1 %d/8 rank %.1f p %.3f | R3 only %+.3f top1 %d rank %.1f p %.3f | all4 %+.3f top1 %d rank %.1f p %.3f\n", sp, b$subset, b$diag_mean_resid, b$top1_rows, b$mean_row_rank, b$perm_p, r3$diag_mean_resid, r3$top1_rows, r3$mean_row_rank, r3$perm_p, a4$diag_mean_resid, a4$top1_rows, a4$mean_row_rank, a4$perm_p)) }
cat("\nnumber of (combination x space) matrices with permutation p < 0.05:", sum(met$perm_p < 0.05), "of", nrow(met), "| < 0.10:", sum(met$perm_p < 0.10), "\n")
cat("\nmean per-KO diagonal residual over all combinations (rows = KO; columns = spaces):\n"); print(round(sapply(spaces, function(sp) rowMeans(sapply(lab, function(s) diag(Z[[sp]][[s]])))), 3))
cat("\ncorrelation between the 64 residual entries of different combinations (median | min) per space:\n")
for (sp in spaces) { V <- sapply(lab, function(s) as.vector(Z[[sp]][[s]])); ut <- cor(V)[upper.tri(diag(length(lab)))]; cat(sprintf("  %-24s %.2f | %.2f\n", sp, median(ut), min(ut))) }

for (sp in spaces) save_fig(paste0("fig3U_KOxKO_matrix_dc_", sp_tag[[sp]]), 15, 15.5, function() {
  layout(matrix(1:16, nrow = 4, byrow = TRUE)); par(oma = c(0, 0, 3.4, 0), family = "sans"); d <- met[met$space == sp, ]; best <- d$subset[which.max(d$diag_mean_resid)]
  vm <- max(abs(unlist(lapply(Z[[sp]], as.vector))))
  for (s in lab) { C <- Z[[sp]][[s]]; m <- d[d$subset == s, ]; par(mar = c(3.6, 3.6, 3.4, 0.5))
    image(seq_len(K), seq_len(K), t(C)[, K:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(K)) for (j in seq_len(K)) text(j, K - i + 1, sprintf("%.2f", C[i, j]), cex = 0.42, col = INK)
    rect(seq_len(K) - 0.5, K - seq_len(K) + 0.5, seq_len(K) + 0.5, K - seq_len(K) + 1.5, border = "black", lwd = 2)
    axis(1, at = seq_len(K), labels = KOS, las = 2, tick = FALSE, cex.axis = 0.62, line = -0.6); axis(2, at = K:1, labels = KOS, las = 1, tick = FALSE, cex.axis = 0.62, line = -0.6)
    box(col = if (s == best) "#d62728" else MUTED, lwd = if (s == best) 3 else 0.6)
    mtext(s, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.62, col = if (s == best) "#d62728" else INK)
    mtext(sprintf("diag resid %+.3f | top1 %d/8 | p=%.3f", m$diag_mean_resid, m$top1_rows, m$perm_p), side = 3, adj = 0, line = 0.4, cex = 0.55, col = INK2) }
  par(mar = c(3.6, 1, 3.4, 1)); plot.new(); text(0.02, 0.9, "rows = CROP-seq screen KO\ncolumns = bulk H2052 KO\nSpearman rho with row and column\nmeans removed (interaction part)\nblack squares = diagonal\nred frame = best diagonal", adj = c(0, 1), cex = 0.85)
  xs <- seq(0.05, 0.95, length.out = 42); rect(xs[-42], 0.12, xs[-1], 0.2, col = DIV[round(seq(1, 101, length.out = 41))], border = NA); text(c(0.05, 0.5, 0.95), 0.05, sprintf("%+.2f", c(-vm, 0, vm)), cex = 0.7); text(0.5, 0.26, "residual", cex = 0.7)
  mtext(sprintf("Double-centered KO x KO correlation (row + column effects removed), screen vs bulk H2052, every H2052 NTC combination: %s", sp), outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 0.95) })

save_fig("fig3V_dc_diagonal_summary", 11, 6.5, function() {
  W <- sapply(spaces, function(sp) sapply(summ$subset, function(s) met$diag_mean_resid[met$space == sp & met$subset == s])); rownames(W) <- summ$subset; n <- nrow(W); vm <- max(abs(W))
  T1 <- sapply(spaces, function(sp) sapply(summ$subset, function(s) met$top1_rows[met$space == sp & met$subset == s]))
  par(mar = c(4.5, 15, 3.6, 1), family = "sans"); image(seq_len(ncol(W)), seq_len(n), t(W)[, n:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(n)) for (j in seq_len(ncol(W))) text(j, n - i + 1, sprintf("%+.3f  (%d/8)", W[i, j], T1[i, j]), cex = 0.68, col = INK)
  axis(1, at = seq_len(ncol(W)), labels = FALSE); text(seq_len(ncol(W)), 0.2, colnames(W), srt = 20, adj = 1, xpd = NA, cex = 0.8)
  text(0.35, n:1, rownames(W), adj = 1, xpd = NA, cex = 0.7, col = ifelse(rownames(W) == "H2052R3", "#eda100", ifelse(summ$n_ntc == 4, "#4a3aa7", INK)), font = ifelse(rownames(W) == "H2052R3" | summ$n_ntc == 4, 2, 1)); box(col = MUTED, lwd = 0.6)
  mtext("Double-centered diagonal quality per H2052 NTC combination (rows ranked by mean over spaces)", side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext("cell = mean diagonal residual of the KO x KO matrix after removing row and column means (chance = 0) and (number of screen KOs whose best bulk match is the same KO)", side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2) })
