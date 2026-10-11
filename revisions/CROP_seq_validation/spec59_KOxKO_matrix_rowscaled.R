# Spec 59: row-wise scaling of the spec58 KO x KO matrices. Each row (one screen KO) is z-scored across the 8 bulk-H2052 KOs, which removes the column/offset
# effect (a bulk KO profile that correlates with every screen KO) and asks only: does this screen KO rank its own bulk KO high? Diagonal metrics on the scaled matrices:
# mean diagonal z, top1 (# screen KOs whose max z is on the diagonal), mean rank of the diagonal within its row, exact permutation p (all 8! relabelings) for the mean diagonal z.
# needs: spec58_H2052_NTC_combos_KOxKO_matrix.R (spec58_KOxKO_matrices.rds)
# outputs: spec59_rowscaled_diagonal_metrics.csv, fig3S_KOxKO_matrix_rowz_{Cm,top1000,top2000,top3000}.{png,pdf}, fig3T_rowz_diagonal_summary.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
mats <- readRDS(file.path(OUT, "spec58_KOxKO_matrices.rds")); spaces <- names(mats); sp_tag <- c("Cm", "top1000", "top2000", "top3000"); names(sp_tag) <- spaces
lab <- names(mats[[1]]); K <- nrow(mats[[1]][[1]]); KOS <- rownames(mats[[1]][[1]]); nn <- sapply(strsplit(lab, "+", fixed = TRUE), length)
perms <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(seq_len(K))) })
rz <- function(C) { Z <- t(scale(t(C))); dimnames(Z) <- dimnames(C); Z }   # row-wise z-score
met <- list(); Z <- list()
for (sp in spaces) for (s in lab) { C <- rz(mats[[sp]][[s]]); Z[[sp]][[s]] <- C; dm <- mean(diag(C)); pd <- apply(perms, 1, function(p) mean(C[cbind(seq_len(K), p)]))
  met[[length(met) + 1]] <- data.frame(space = sp, subset = s, n_ntc = length(strsplit(s, "+", fixed = TRUE)[[1]]), diag_mean_z = dm, top1_rows = sum(apply(C, 1, which.max) == seq_len(K)),
       mean_row_rank = mean(sapply(seq_len(K), function(r) rank(-C[r, ])[r])), perm_p = mean(pd >= dm - 1e-12)) }
met <- do.call(rbind, met); write.csv(met, file.path(OUT, "spec59_rowscaled_diagonal_metrics.csv"), row.names = FALSE)
summ <- do.call(rbind, lapply(split(met, met$subset), function(d) data.frame(subset = d$subset[1], n_ntc = d$n_ntc[1], mean_diag_z = mean(d$diag_mean_z), mean_top1 = mean(d$top1_rows), mean_row_rank = mean(d$mean_row_rank), min_perm_p = min(d$perm_p),
     Cm = d$diag_mean_z[d$space == spaces[1]], top1000 = d$diag_mean_z[d$space == spaces[2]], top2000 = d$diag_mean_z[d$space == spaces[3]], top3000 = d$diag_mean_z[d$space == spaces[4]])))
summ <- summ[order(-summ$mean_diag_z), ]; options(width = 230)
cat("\n=== row-scaled (z per screen KO across bulk KOs): diagonal quality per H2052 NTC combination; chance mean_diag_z = 0, top1 = 1/8, mean row rank = 4.5 ===\n"); print(format(summ, digits = 2), row.names = FALSE)
cat("\n=== per space: best combination, its top1, mean diagonal row rank, permutation p; R3-only and all-4 for reference ===\n")
for (sp in spaces) { d <- met[met$space == sp, ]; b <- d[which.max(d$diag_mean_z), ]; r3 <- d[d$subset == "H2052R3", ]; a4 <- d[d$n_ntc == 4, ]
  cat(sprintf("%-24s best %-24s z %+.2f top1 %d/8 rank %.1f p %.3f | R3 only z %+.2f top1 %d rank %.1f p %.3f | all4 z %+.2f top1 %d rank %.1f p %.3f\n", sp, b$subset, b$diag_mean_z, b$top1_rows, b$mean_row_rank, b$perm_p, r3$diag_mean_z, r3$top1_rows, r3$mean_row_rank, r3$perm_p, a4$diag_mean_z, a4$top1_rows, a4$mean_row_rank, a4$perm_p)) }
cat("\nper-KO diagonal z in the best combination overall (", summ$subset[1], ") and in R3 only, per space:\n")
for (s in c(summ$subset[1], "H2052R3")) { cat(" ", s, "\n"); print(round(sapply(spaces, function(sp) diag(Z[[sp]][[s]])), 2)) }

for (sp in spaces) save_fig(paste0("fig3S_KOxKO_matrix_rowz_", sp_tag[[sp]]), 15, 15.5, function() {
  layout(matrix(1:16, nrow = 4, byrow = TRUE)); par(oma = c(0, 0, 3.4, 0), family = "sans"); d <- met[met$space == sp, ]; best <- d$subset[which.max(d$diag_mean_z)]
  for (s in lab) { C <- Z[[sp]][[s]]; m <- d[d$subset == s, ]; par(mar = c(3.6, 3.6, 3.4, 0.5))
    image(seq_len(K), seq_len(K), t(C)[, K:1], col = DIV, zlim = c(-2.5, 2.5), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(K)) for (j in seq_len(K)) text(j, K - i + 1, sprintf("%.1f", C[i, j]), cex = 0.5, col = INK)
    rect(seq_len(K) - 0.5, K - seq_len(K) + 0.5, seq_len(K) + 0.5, K - seq_len(K) + 1.5, border = "black", lwd = 2)
    axis(1, at = seq_len(K), labels = KOS, las = 2, tick = FALSE, cex.axis = 0.62, line = -0.6); axis(2, at = K:1, labels = KOS, las = 1, tick = FALSE, cex.axis = 0.62, line = -0.6)
    box(col = if (s == best) "#d62728" else MUTED, lwd = if (s == best) 3 else 0.6)
    mtext(s, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.62, col = if (s == best) "#d62728" else INK)
    mtext(sprintf("diag z %+.2f | top1 %d/8 | p=%.3f", m$diag_mean_z, m$top1_rows, m$perm_p), side = 3, adj = 0, line = 0.4, cex = 0.55, col = INK2) }
  par(mar = c(3.6, 1, 3.4, 1)); plot.new(); text(0.02, 0.9, "rows = CROP-seq screen KO\ncolumns = bulk H2052 KO\nSpearman rho, z-scored per row\n(across the 8 bulk KOs)\nblack squares = diagonal\nred frame = best diagonal z", adj = c(0, 1), cex = 0.85)
  xs <- seq(0.05, 0.95, length.out = 42); rect(xs[-42], 0.12, xs[-1], 0.2, col = DIV[round(seq(1, 101, length.out = 41))], border = NA); text(c(0.05, 0.5, 0.95), 0.05, c("-2.5", "0", "+2.5"), cex = 0.7); text(0.5, 0.26, "row z-score", cex = 0.7)
  mtext(sprintf("Row-scaled KO x KO correlation, screen vs bulk H2052, every H2052 NTC combination: %s", sp), outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 1.0) })

save_fig("fig3T_rowz_diagonal_summary", 11, 6.5, function() {
  W <- sapply(spaces, function(sp) sapply(summ$subset, function(s) met$diag_mean_z[met$space == sp & met$subset == s])); rownames(W) <- summ$subset; n <- nrow(W); vm <- max(abs(W))
  T1 <- sapply(spaces, function(sp) sapply(summ$subset, function(s) met$top1_rows[met$space == sp & met$subset == s]))
  par(mar = c(4.5, 15, 3.6, 1), family = "sans"); image(seq_len(ncol(W)), seq_len(n), t(W)[, n:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(n)) for (j in seq_len(ncol(W))) text(j, n - i + 1, sprintf("%+.2f  (%d/8)", W[i, j], T1[i, j]), cex = 0.68, col = INK)
  axis(1, at = seq_len(ncol(W)), labels = FALSE); text(seq_len(ncol(W)), 0.2, colnames(W), srt = 20, adj = 1, xpd = NA, cex = 0.8)
  text(0.35, n:1, rownames(W), adj = 1, xpd = NA, cex = 0.7, col = ifelse(rownames(W) == "H2052R3", "#eda100", ifelse(summ$n_ntc == 4, "#4a3aa7", INK)), font = ifelse(rownames(W) == "H2052R3" | summ$n_ntc == 4, 2, 1)); box(col = MUTED, lwd = 0.6)
  mtext("Row-scaled diagonal quality per H2052 NTC combination (rows ranked by mean over spaces)", side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext("cell = mean diagonal z-score (each screen KO z-scored across bulk KOs; chance = 0) and (number of screen KOs whose best bulk match is the same KO)", side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2) })
