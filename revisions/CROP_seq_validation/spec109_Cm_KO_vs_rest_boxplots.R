# Spec 109: specific Cm module scores in one KO vs the rest (NTC + the other 7 KOs), with REPLICATES AVERAGED and the 4 CELL LINES AGGREGATED.
# Pairs: Cm20/TEAD4, Cm8+Cm9/TWIST1, Cm10/TCF3, Cm10+Cm1/TEAD2.
#
# The unit of observation is the (cell line, group) replicate-mean, not the individual sample. Within each line the module is z-scored
# across that line's samples, then the replicates of each of the 9 groups are averaged. The target KO therefore contributes 4 values
# (one per line) and the rest contributes 8 non-target groups x 4 lines = 32. Averaging first matters: the groups have 1-3 replicates,
# so pooling raw samples would weight a 3-replicate NTC more than a 1-replicate KO, and would also treat technical replicates as if they
# were independent observations of the KO effect.
# z within line is a linear transform inside a line, so it changes nothing about the comparison; it only puts the 4 lines on one axis.
#
# Significance: group-label permutation WITHIN each line (10,000) -- in every line one of the 9 groups is relabelled "the KO" and the
# statistic mean(KO) - mean(rest) recomputed, so the null preserves both the per-line structure and the 4-vs-32 design. Wilcoxon and
# Welch on the same 4 vs 32 values are reported alongside. BH across the 6 pairs.
# fig21b keeps the per-line view (KO vs the rest of its own line) for reference.
# needs: repro_common.R, profiles.R, agree_common.R (KOS), spec72_profiles_S70.rds
# outputs: spec109_Cm_KO_vs_rest_stats.csv, spec109_Cm_KO_vs_rest_group_means.csv,
#          fig21_Cm_KO_vs_rest_boxplots.{png,pdf}, fig21b_Cm_KO_vs_rest_per_line.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "10000")); GR <- c("NTC", KOS)
PAIRS <- list(c("TEAD4", "Cm20"), c("TWIST1", "Cm8"), c("TWIST1", "Cm9"), c("TCF3", "Cm10"), c("TEAD2", "Cm10"), c("TEAD2", "Cm1"))

pr <- readRDS(file.path(OUT, "spec72_profiles_S70.rds")); fab <- load_full_audited_bulk(); meta <- fab$meta
ids <- pr$samples; line <- meta[ids, "Cell_Line"]; grp <- meta[ids, "Gene_Control"]
S <- pr$cm_scores$tmm[, ids, drop = FALSE]
Z <- S; for (cl in LINES) { j <- which(line == cl); Z[, j] <- t(scale(t(S[, j, drop = FALSE]))) }; Z[!is.finite(Z)] <- NA

## G[module, group, line] = replicate mean of the within-line z
G <- array(NA_real_, c(nrow(Z), length(GR), length(LINES)), dimnames = list(rownames(Z), GR, LINES))
nrep <- matrix(0L, length(GR), length(LINES), dimnames = list(GR, LINES))
for (cl in LINES) for (g in GR) { j <- which(line == cl & grp == g); nrep[g, cl] <- length(j)
  if (length(j)) G[, g, cl] <- rowMeans(Z[, j, drop = FALSE], na.rm = TRUE) }
cat("replicates averaged per group x line:\n"); print(nrep)
write.csv(do.call(rbind, lapply(LINES, function(cl) data.frame(line = cl, group = GR, t(round(G[, , cl], 4))))),
          file.path(OUT, "spec109_Cm_KO_vs_rest_group_means.csv"), row.names = FALSE)

set.seed(1)
rows <- list()
for (pp in PAIRS) { ko <- pp[1]; cm <- pp[2]
  a <- G[cm, ko, ]; b <- as.vector(G[cm, setdiff(GR, ko), ])          # 4 vs 32 replicate-means
  a <- a[is.finite(a)]; b <- b[is.finite(b)]; obs <- mean(a) - mean(b)
  null <- replicate(N_PERM, { d <- sapply(LINES, function(cl) { v <- G[cm, , cl]; v <- v[is.finite(v)]
      k <- sample(length(v), 1); c(v[k], mean(v[-k])) })
    mean(d[1, ]) - mean(d[2, ]) })
  rows[[length(rows) + 1]] <- data.frame(KO = ko, Cm = cm, n_KO = length(a), n_rest = length(b),
    mean_KO = mean(a), mean_rest = mean(b), delta = obs,
    perm_p = (1 + sum(abs(null) >= abs(obs))) / (N_PERM + 1),
    wilcox_p = wilcox.test(a, b)$p.value, welch_p = t.test(a, b)$p.value) }
res <- do.call(rbind, rows); res$perm_fdr <- p.adjust(res$perm_p, "BH")
write.csv(res, file.path(OUT, "spec109_Cm_KO_vs_rest_stats.csv"), row.names = FALSE)
options(width = 200)
cat("\n=== KO vs rest, replicates averaged and the 4 lines aggregated (unit = line x group mean) ===\n")
print(format(res, digits = 3), row.names = FALSE)
cat("\nthe KO's replicate-mean in each line (the 4 points in the coloured box):\n")
for (pp in PAIRS) cat(sprintf("%-7s %-5s : %s\n", pp[1], pp[2],
  paste(sprintf("%s %+.2f", sub("NCI-|-211H", "", LINES), G[pp[2], pp[1], ]), collapse = "  ")))

fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
LCOL <- LINE_COL[LINES]; short <- sub("NCI-|-211H", "", LINES)
save_fig("fig21_Cm_KO_vs_rest_boxplots", 11, 6.8, function() {
  par(mfrow = c(2, 3), oma = c(3.0, 1.2, 3.8, 0.6), mar = c(2.8, 4.3, 3.6, 0.8), family = "sans")
  for (i in seq_len(nrow(res))) { ko <- res$KO[i]; cm <- res$Cm[i]
    a <- G[cm, ko, ]; b <- as.vector(G[cm, setdiff(GR, ko), ]); bl <- rep(LINES, each = length(GR) - 1)
    yr <- range(c(a, b), na.rm = TRUE) + c(-0.12, 0.42) * diff(range(c(a, b), na.rm = TRUE))
    boxplot(list(rest = b, KO = a), ylim = yr, col = c("#d5d4cf", adjustcolor(KO_COL[[ko]], 0.75)), border = INK2,
            outline = FALSE, boxwex = 0.55, las = 1, bty = "n", cex.axis = 0.8,
            names = c(sprintf("rest\n(8 groups x 4 lines, n=%d)", res$n_rest[i]), sprintf("%s\n(1 per line, n=%d)", ko, res$n_KO[i])),
            ylab = sprintf("%s score (z within line, replicate mean)", cm))
    abline(h = 0, lty = 3, col = GRID)
    points(jitter(rep(1, length(b)), amount = 0.16), b, pch = 21, bg = adjustcolor(LCOL[bl], 0.6), col = INK2, cex = 0.7)
    points(seq(1.88, 2.12, length.out = length(a)), a, pch = 21, bg = LCOL[LINES], col = INK2, cex = 1.5)
    st <- ifelse(res$perm_fdr[i] < 0.01, " **", ifelse(res$perm_fdr[i] < 0.05, " *", ""))
    mtext(sprintf("%s in %s KO", cm, ko), side = 3, adj = 0, line = 2.0, font = 2, cex = 0.82, col = KO_COL[[ko]])
    mtext(sprintf("delta = %+.2f  |  perm p = %s%s", res$delta[i], fp(res$perm_p[i]), st), side = 3, adj = 0, line = 1.0,
          cex = 0.6, col = INK2, font = ifelse(res$perm_p[i] < 0.05, 2, 1))
    mtext(sprintf("wilcox %s | welch %s", fp(res$wilcox_p[i]), fp(res$welch_p[i])), side = 3, adj = 0, line = 0.2, cex = 0.5, col = MUTED) }
  mtext("Cm module score: the KO vs the rest (NTC + the other 7 KOs), replicates averaged, 4 cell lines aggregated", outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 0.95)
  mtext(sprintf("unit = one (cell line, group) replicate-mean of the within-line z; the KO contributes one point per line; group-label permutation within line (%s perms), two-sided; * / ** BH FDR < 0.05 / 0.01", format(N_PERM, big.mark = ",")),
        outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.57, col = INK2)
  legend("bottom", short, pt.bg = LCOL, pch = 21, col = INK2, horiz = TRUE, bty = "n", xpd = NA, inset = c(0, -0.30), cex = 0.78, pt.cex = 1.2)
})

## ---------------------------------------------------------------- fig21b: the per-line view, kept for reference
pl <- do.call(rbind, lapply(PAIRS, function(pp) { ko <- pp[1]; cm <- pp[2]
  do.call(rbind, lapply(LINES, function(cl) data.frame(KO = ko, Cm = cm, line = cl,
    KO_mean = G[cm, ko, cl], rest_mean = mean(G[cm, setdiff(GR, ko), cl], na.rm = TRUE),
    delta = G[cm, ko, cl] - mean(G[cm, setdiff(GR, ko), cl], na.rm = TRUE))) ) }))
write.csv(pl, file.path(OUT, "spec109_Cm_KO_vs_rest_per_line.csv"), row.names = FALSE)
cat("\n=== per line: KO replicate-mean minus the mean of the other 8 group means ===\n")
print(format(pl, digits = 3), row.names = FALSE)
save_fig("fig21b_Cm_KO_vs_rest_per_line", 11, 6.8, function() {
  par(mfrow = c(2, 3), oma = c(2.6, 1.2, 3.6, 0.6), mar = c(3.0, 4.3, 3.4, 0.8), family = "sans")
  for (i in seq_len(nrow(res))) { ko <- res$KO[i]; cm <- res$Cm[i]; d <- pl[pl$KO == ko & pl$Cm == cm, ]
    bp <- barplot(d$delta, names.arg = short, col = adjustcolor(KO_COL[[ko]], 0.8), border = INK2, las = 1,
                  ylim = range(c(0, d$delta)) * 1.25, ylab = sprintf("%s: KO - rest (z)", cm), cex.axis = 0.8, cex.names = 0.72)
    abline(h = 0, col = INK, lwd = 0.6)
    mtext(sprintf("%s in %s KO", cm, ko), side = 3, adj = 0, line = 1.4, font = 2, cex = 0.8, col = KO_COL[[ko]])
    mtext(sprintf("%d of 4 lines agree in sign", max(sum(d$delta > 0), sum(d$delta < 0))), side = 3, adj = 0, line = 0.4, cex = 0.58, col = INK2) }
  mtext("Per cell line: the KO's replicate-mean minus the mean of the other 8 groups of that line", outer = TRUE, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.92)
  mtext("same quantity as the left-hand figure before the lines are aggregated; a pair only survives aggregation if the 4 bars share a sign", outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2)
})
