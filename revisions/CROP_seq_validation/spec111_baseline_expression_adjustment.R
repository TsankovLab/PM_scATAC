# Spec 111: does the KO effect on a Cm module scale with how much of the target TF the line expressed to begin with?
# Rationale: a line that barely expresses TEAD4 has little to lose when TEAD4 is knocked out, so its delta is mostly noise and should not
# be averaged in with equal weight. If that is right, (a) |delta| should grow with baseline target expression, and (b) down-weighting or
# dropping the low-baseline lines should sharpen the aggregated KO-vs-rest separation.
#
# baseline(line, KO) = mean log2 CP10K of the target gene across that line's NON-target samples (NTC + the other 7 KOs). The rest-mean is
#   used rather than NTC alone because H2052 has a single NTC; NTC-only is computed too and reported for comparison.
# delta(line, Cm, KO)  = the KO's replicate-mean minus the mean of the other 8 group means, on the within-line z (as in spec109).
# Three aggregations of the 4 lines are compared per pair, each with the same within-line group-label permutation (10,000):
#   equal      the spec109 result, every line weight 1
#   weighted   line weight = baseline, so a line that does not express the target contributes little
#   top2       only the 2 lines with the highest baseline
# "Separation" is judged on the permutation p and on how many lines agree in sign -- not on delta alone, which weighting can inflate.
# needs: repro_common.R, profiles.R, agree_common.R, spec72_profiles_S70.rds
# outputs: spec111_baseline_vs_effect.csv, spec111_adjustment_comparison.csv, fig23_baseline_expression_adjustment.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "10000")); GR <- c("NTC", KOS); CM <- paste0("Cm", 1:20)
PAIRS <- list(c("TEAD4", "Cm20"), c("TWIST1", "Cm8"), c("TWIST1", "Cm9"), c("TCF3", "Cm10"), c("TEAD2", "Cm10"), c("TEAD2", "Cm1"))

pr <- readRDS(file.path(OUT, "spec72_profiles_S70.rds")); fab <- load_full_audited_bulk(); meta <- fab$meta
ids <- pr$samples; line <- meta[ids, "Cell_Line"]; grp <- meta[ids, "Gene_Control"]
lg <- log2(sweep(fab$cs[, ids, drop = FALSE], 2, colSums(fab$cs[, ids, drop = FALSE]), "/") * 1e4 + 1)
S <- pr$cm_scores$tmm[, ids, drop = FALSE]
Z <- S; for (cl in LINES) { j <- which(line == cl); Z[, j] <- t(scale(t(S[, j, drop = FALSE]))) }; Z[!is.finite(Z)] <- NA
G <- array(NA_real_, c(nrow(Z), length(GR), length(LINES)), dimnames = list(rownames(Z), GR, LINES))
for (cl in LINES) for (g in GR) { j <- which(line == cl & grp == g); if (length(j)) G[, g, cl] <- rowMeans(Z[, j, drop = FALSE], na.rm = TRUE) }

## baseline expression of each target TF in each line
base <- sapply(LINES, function(cl) sapply(KOS, function(k) mean(lg[k, line == cl & grp != k])))
base_ntc <- sapply(LINES, function(cl) sapply(KOS, function(k) mean(lg[k, line == cl & grp == "NTC"])))
cat("=== baseline target expression (mean log2 CP10K in the non-target samples of each line) ===\n"); print(round(base, 2))
cat("\nsame from NTC samples only:\n"); print(round(base_ntc, 2))
cat("\ncorrelation between the two definitions:", round(cor(as.vector(base), as.vector(base_ntc)), 3), "\n")

## (a) does baseline predict |delta|? over all 20 Cm x 8 KOs x 4 lines
dl <- array(NA_real_, c(length(CM), length(KOS), length(LINES)), dimnames = list(CM, KOS, LINES))
for (cl in LINES) for (k in KOS) dl[, k, cl] <- G[CM, k, cl] - rowMeans(G[CM, setdiff(GR, k), cl], na.rm = TRUE)
A <- data.frame(Cm = rep(CM, length(KOS) * length(LINES)), KO = rep(rep(KOS, each = length(CM)), length(LINES)),
                line = rep(LINES, each = length(CM) * length(KOS)), delta = as.vector(dl),
                baseline = as.vector(sapply(LINES, function(cl) rep(base[, cl], each = length(CM)))))
A$abs_delta <- abs(A$delta); write.csv(A, file.path(OUT, "spec111_baseline_vs_effect.csv"), row.names = FALSE)
r_all <- suppressWarnings(cor.test(A$baseline, A$abs_delta, method = "spearman", exact = FALSE))
cat(sprintf("\n=== |delta| vs baseline target expression, all %d (Cm x KO x line) observations ===\nSpearman rho = %.3f, p = %.3g\n", nrow(A), r_all$estimate, r_all$p.value))
cat("\nper KO (over 20 modules x 4 lines):\n")
print(round(t(sapply(KOS, function(k) { s <- A[A$KO == k, ]; r <- suppressWarnings(cor.test(s$baseline, s$abs_delta, method = "spearman", exact = FALSE))
  c(rho = unname(r$estimate), p = r$p.value) })), 3))
cat("\nmean |delta| by baseline quartile:\n")
A$q <- cut(A$baseline, quantile(A$baseline, 0:4 / 4), include.lowest = TRUE, labels = paste0("Q", 1:4))
print(round(tapply(A$abs_delta, A$q, mean), 3))

## (b) does weighting / dropping low-baseline lines improve the 6 pairs?
set.seed(1)
agg <- function(cm, k, w) { d <- dl[cm, k, ]; ok <- is.finite(d) & w > 0; sum(d[ok] * w[ok]) / sum(w[ok]) }
perm_agg <- function(cm, k, w) {
  obs <- agg(cm, k, w)
  null <- replicate(N_PERM, { d <- sapply(LINES, function(cl) { v <- G[cm, , cl]; v <- v[is.finite(v)]
      j <- sample(length(v), 1); v[j] - mean(v[-j]) })
    ok <- is.finite(d) & w > 0; sum(d[ok] * w[ok]) / sum(w[ok]) })
  c(delta = obs, p = (1 + sum(abs(null) >= abs(obs))) / (N_PERM + 1))
}
cmp <- do.call(rbind, lapply(PAIRS, function(pp) { k <- pp[1]; cm <- pp[2]; b <- base[k, ]
  w_eq <- setNames(rep(1, length(LINES)), LINES); w_wt <- b
  w_t2 <- setNames(as.numeric(rank(-b) <= 2), LINES)
  e <- perm_agg(cm, k, w_eq); w <- perm_agg(cm, k, w_wt); t2 <- perm_agg(cm, k, w_t2)
  data.frame(KO = k, Cm = cm, n_lines_same_sign = max(sum(dl[cm, k, ] > 0), sum(dl[cm, k, ] < 0)),
             baseline_range = sprintf("%.2f-%.2f", min(b), max(b)),
             delta_equal = e[1], p_equal = e[2], delta_weighted = w[1], p_weighted = w[2],
             delta_top2 = t2[1], p_top2 = t2[2], top2_lines = paste(sub("NCI-|-211H", "", LINES[rank(-b) <= 2]), collapse = "+")) }))
write.csv(cmp, file.path(OUT, "spec111_adjustment_comparison.csv"), row.names = FALSE)
cat("\n=== aggregation of the 4 lines: equal vs baseline-weighted vs the 2 highest-baseline lines ===\n")
print(format(cmp, digits = 3), row.names = FALSE)
cat("\nmedian permutation p:  equal", round(median(cmp$p_equal), 3), " weighted", round(median(cmp$p_weighted), 3),
    " top2", round(median(cmp$p_top2), 3), "\n")
cat("pairs with p < 0.05:   equal", sum(cmp$p_equal < 0.05), " weighted", sum(cmp$p_weighted < 0.05), " top2", sum(cmp$p_top2 < 0.05), "\n")

fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
save_fig("fig23_baseline_expression_adjustment", 12.5, 5.6, function() {
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1.15, 1, 1.1)); par(family = "sans", oma = c(0, 0, 3.4, 0.4), mar = c(4.4, 4.6, 3.0, 1))
  plot(A$baseline, A$abs_delta, pch = 21, bg = adjustcolor(KO_COL[A$KO], 0.6), col = NA, cex = 0.8, bty = "n", las = 1,
       xlab = "baseline target expression in that line (log2 CP10K)", ylab = "|KO effect| on the module (z)", cex.axis = 0.85)
  abline(lm(abs_delta ~ baseline, A), col = INK, lty = 2)
  mtext("all 20 Cm x 8 KO x 4 lines", side = 3, adj = 0, line = 1.5, font = 2, cex = 0.85)
  mtext(sprintf("Spearman %.2f, p = %s (%d points)", r_all$estimate, fp(r_all$p.value), nrow(A)), side = 3, adj = 0, line = 0.4, cex = 0.62, col = INK2)
  bx <- boxplot(abs_delta ~ q, A, col = "#d5d4cf", border = INK2, outline = FALSE, las = 1, bty = "n", cex.axis = 0.85,
                xlab = "baseline expression quartile", ylab = "|KO effect| (z)")
  mtext("binned", side = 3, adj = 0, line = 1.5, font = 2, cex = 0.85)
  mtext(sprintf("mean |delta|: %s", paste(sprintf("%s %.2f", levels(A$q), tapply(A$abs_delta, A$q, mean)), collapse = "  ")), side = 3, adj = 0, line = 0.4, cex = 0.58, col = INK2)
  m <- t(as.matrix(cmp[, c("p_equal", "p_weighted", "p_top2")])); colnames(m) <- paste(cmp$Cm, cmp$KO)
  bp <- barplot(m, beside = TRUE, col = c("#d5d4cf", "#8fa88a", "#4a3aa7"), border = INK2, las = 2, cex.names = 0.6,
                ylab = "permutation p", ylim = c(0, max(m) * 1.25))
  abline(h = 0.05, lty = 2, col = "red"); legend("topright", c("equal", "baseline-weighted", "top-2 baseline lines"),
         fill = c("#d5d4cf", "#8fa88a", "#4a3aa7"), border = NA, bty = "n", cex = 0.7)
  mtext("the 6 pairs: does adjusting help?", side = 3, adj = 0, line = 1.5, font = 2, cex = 0.85)
  mtext("lower = better separation; red = 0.05", side = 3, adj = 0, line = 0.4, cex = 0.6, col = INK2)
  mtext("Does baseline expression of the target TF explain the KO effect, and does conditioning on it improve separation?", outer = TRUE, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.95)
  mtext("baseline = mean log2 CP10K of the target in that line's non-target samples; effect = KO replicate-mean minus the mean of the other 8 group means (within-line z)", outer = TRUE, side = 3, adj = 0, line = 0.4, cex = 0.6, col = INK2)
})
