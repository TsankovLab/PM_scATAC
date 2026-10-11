# Spec 70: systematic benchmark of how to match CROP-seq screen KO signal to bulk RNA-seq KO signal. Instead of testing one normalization at a time, every
# analysis choice varied across spec44-69 is crossed in one grid and scored with ONE metric, per bulk target:
#   statistic   : diff (screen mean log-norm diff vs NTC | bulk DESeq2 log2FC)  or  t (screen Welch t | bulk DESeq2 Wald stat); with/without cell-cycle regression
#   features    : top-N genes ranked by (a) bulk variance across lines [current default], (b) screen signal = max |t| over screen KOs,
#                 (c) bulk signal = mean over lines of max |Wald| over KOs, (d) both = rank-sum of (b)+(c); plus all shared genes
#   N           : 200, 500, 1000, 3000
#   centring    : none  or  KO-centred (subtract each gene's mean across the 8 KO profiles, per platform -> removes the effect shared by all KOs vs the common NTC)
#   correlation : Pearson or Spearman
#   bulk target : each of the 4 lines, and POOLED = mean of the per-line statistics over MSTO/H2452/H28 (the lines with >=2 NTC)
# Metric: 8x8 screen-KO x bulk-KO correlation matrix -> diag_minus_off (mean diagonal - mean off-diagonal), top1 (# screen KOs whose best bulk match is itself),
# exact permutation p (all 8! bulk-label relabelings). NTC-subtracted by construction (both platforms' statistics are KO vs NTC).
# Overfitting guard: leave-one-line-out -- pick the config maximising mean diag_minus_off on 3 lines, score it on the held-out line; compared with the
# fixed default (diff, bulk-variance top1000, no centring, Spearman).
# needs: spec1_prepare_stats.R (spec_screen_stats.rds, spec_bulk_stats.rds)
# outputs: spec70_benchmark_grid.csv, spec70_benchmark_LOLO.csv, spec70_best_config_per_KO.csv, fig6A_matching_benchmark.{png,pdf}, fig6B_best_config_matrices.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
ss <- readRDS(file.path(OUT, "spec_screen_stats.rds")); bs <- readRDS(file.path(OUT, "spec_bulk_stats.rds"))
KOS <- ORDER; K <- length(KOS); POOL_LINES <- c("MSTO-211H", "NCI-H2452", "NCI-H28")
genes <- Reduce(intersect, c(list(ss$genes), lapply(LINES, function(cl) bs[[cl]]$genes)))
cat(sprintf("shared genes (screen-detected & tested in all 4 bulk lines): %d\n", length(genes)))

STAT <- list(diff = c(s = "diff0", b = "lfc0"), t = c(s = "t0", b = "stat0"), diff_cc = c(s = "diffcc", b = "lfccc"), t_cc = c(s = "tcc", b = "statcc"))
clean <- function(M) { M[!is.finite(M)] <- 0; M }
S_all <- lapply(STAT, function(x) clean(ss[[x[["s"]]]][genes, KOS]))
B_all <- lapply(STAT, function(x) { L <- lapply(setNames(LINES, LINES), function(cl) clean(bs[[cl]][[x[["b"]]]][genes, KOS]))
  L[["pooled"]] <- Reduce(`+`, L[POOL_LINES]) / length(POOL_LINES); L })
TARGETS <- c(LINES, "pooled")

## feature rankings (computed once, from the uncorrected t / Wald statistics; bulk variance from log2 CP10K as in spec22/47)
b <- load_bulk(); lg <- log2(b$cp10k + 1); meta <- b$meta
bvar <- rowMeans(sapply(LINES, function(cl) matrixStats::rowVars(as.matrix(lg[genes, rownames(meta)[meta$Cell_Line == cl]]))))
rank_by <- list(bulk_var = order(-bvar),
                screen_signal = order(-apply(abs(S_all$t), 1, max)),
                bulk_signal = order(-rowMeans(sapply(POOL_LINES, function(cl) apply(abs(B_all$t[[cl]]), 1, max)))))
rank_by$both_signal <- order(rank(-apply(abs(S_all$t), 1, max)) + rank(-rowMeans(sapply(POOL_LINES, function(cl) apply(abs(B_all$t[[cl]]), 1, max)))))
NS <- c(200, 500, 1000, 3000)

perms <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(seq_len(K))) })
pidx <- cbind(rep(seq_len(K), times = nrow(perms)), as.vector(t(perms)))
score <- function(C) { dm <- mean(diag(C)); om <- mean(C[row(C) != col(C)]); pd <- colMeans(matrix(C[pidx], nrow = K))
  c(diag_minus_off = dm - om, top1 = sum(apply(C, 1, which.max) == seq_len(K)), perm_p = mean(pd >= dm - 1e-12)) }
kcentre <- function(M) M - rowMeans(M)

grid <- expand.grid(stat = names(STAT), feat = c(names(rank_by), "all"), N = NS, centre = c(FALSE, TRUE), method = c("pearson", "spearman"), stringsAsFactors = FALSE)
grid <- unique(transform(grid, N = ifelse(feat == "all", length(genes), N)))
cat(sprintf("configurations: %d x %d targets\n", nrow(grid), length(TARGETS)))
res <- vector("list", nrow(grid) * length(TARGETS)); ii <- 0; diagv <- list()
for (gi in seq_len(nrow(grid))) { g <- grid[gi, ]
  f <- if (g$feat == "all") seq_along(genes) else rank_by[[g$feat]][seq_len(g$N)]
  S <- S_all[[g$stat]][f, , drop = FALSE]; if (g$centre) S <- kcentre(S)
  for (tg in TARGETS) { B <- B_all[[g$stat]][[tg]][f, , drop = FALSE]; if (g$centre) B <- kcentre(B)
    C <- cor(S, B, method = g$method); ii <- ii + 1
    res[[ii]] <- data.frame(g, target = tg, t(score(C))); diagv[[ii]] <- data.frame(g, target = tg, KO = KOS, own = diag(C), null = sapply(seq_len(K), function(i) median(C[i, -i]))) } }
res <- do.call(rbind, res); diagv <- do.call(rbind, diagv)
res$config <- with(res, sprintf("%s | %s%s | %s | %s", stat, feat, ifelse(feat == "all", "", paste0(" top", N)), ifelse(centre, "KO-centred", "uncentred"), method))
diagv$config <- with(diagv, sprintf("%s | %s%s | %s | %s", stat, feat, ifelse(feat == "all", "", paste0(" top", N)), ifelse(centre, "KO-centred", "uncentred"), method))
write.csv(res, file.path(OUT, "spec70_benchmark_grid.csv"), row.names = FALSE)
diagv$sep <- diagv$own - diagv$null; write.csv(diagv, file.path(OUT, "spec70_benchmark_per_KO.csv"), row.names = FALSE)
cat("\n=== per-KO robustness across ALL configurations: fraction of configs with own > median(other), per target ===\n")
print(round(tapply(diagv$sep > 0, list(diagv$KO, diagv$target), mean)[KOS, TARGETS], 2))
cat("\nmedian separation across all configurations, per KO x target:\n"); print(round(tapply(diagv$sep, list(diagv$KO, diagv$target), median)[KOS, TARGETS], 3))

DEFAULT <- "diff | bulk_var top1000 | uncentred | spearman"
agg <- aggregate(cbind(diag_minus_off, top1) ~ config, res[res$target %in% LINES, ], mean)
agg$lines_positive <- sapply(agg$config, function(cf) sum(res$diag_minus_off[res$config == cf & res$target %in% LINES] > 0))
agg$pooled_dmo <- res$diag_minus_off[res$target == "pooled"][match(agg$config, res$config[res$target == "pooled"])]
agg$pooled_p <- res$perm_p[res$target == "pooled"][match(agg$config, res$config[res$target == "pooled"])]
agg <- agg[order(-agg$diag_minus_off), ]
options(width = 230)
cat("\n=== top 20 configurations by mean diag_minus_off over the 4 lines (in-sample; see LOLO below for an unbiased estimate) ===\n")
print(format(head(agg, 20), digits = 3), row.names = FALSE)
cat(sprintf("\ndefault config (%s): rank %d/%d, mean diag-off %.3f, pooled %.3f\n", DEFAULT, match(DEFAULT, agg$config), nrow(agg), agg$diag_minus_off[agg$config == DEFAULT], agg$pooled_dmo[agg$config == DEFAULT]))

## main-effect summary: average diag_minus_off (over lines) for each level of each factor, marginalising the others
cat("\n=== main effects (mean diag_minus_off over lines and all other factor levels) ===\n")
rl <- res[res$target %in% LINES, ]
for (f in c("stat", "feat", "N", "centre", "method")) { m <- tapply(rl$diag_minus_off, rl[[f]], mean); cat(sprintf("%-7s %s\n", f, paste(sprintf("%s=%+.3f", names(m), m), collapse = "  "))) }
cat("\nper-target main effect of centring and of feature ranking:\n"); print(round(tapply(res$diag_minus_off, list(res$target, res$centre), mean), 3)); print(round(tapply(res$diag_minus_off, list(res$target, res$feat), mean), 3))

## leave-one-line-out selection
lolo <- do.call(rbind, lapply(LINES, function(ho) { tr <- aggregate(diag_minus_off ~ config, res[res$target %in% setdiff(LINES, ho), ], mean); best <- tr$config[which.max(tr$diag_minus_off)]
  data.frame(held_out = ho, picked = best, heldout_dmo_picked = res$diag_minus_off[res$target == ho & res$config == best], heldout_p_picked = res$perm_p[res$target == ho & res$config == best],
             heldout_dmo_default = res$diag_minus_off[res$target == ho & res$config == DEFAULT], heldout_dmo_median_config = median(res$diag_minus_off[res$target == ho])) }))
write.csv(lolo, file.path(OUT, "spec70_benchmark_LOLO.csv"), row.names = FALSE)
cat("\n=== leave-one-line-out: config picked on 3 lines, scored on the held-out line ===\n"); print(format(lolo, digits = 3), row.names = FALSE)
cat(sprintf("mean held-out diag-off: picked %.3f | default %.3f | median config %.3f\n", mean(lolo$heldout_dmo_picked), mean(lolo$heldout_dmo_default), mean(lolo$heldout_dmo_median_config)))

## best overall config: per-KO own vs null, per line and pooled
BEST <- agg$config[1]
bk <- diagv[diagv$config == BEST, ]; bk$sep <- bk$own - bk$null
write.csv(bk, file.path(OUT, "spec70_best_config_per_KO.csv"), row.names = FALSE)
cat(sprintf("\n=== best config: %s -- per-KO separation (own - median other) ===\n", BEST)); print(round(tapply(bk$sep, list(bk$KO, bk$target), mean)[KOS, TARGETS], 3))
cat("\nsame for default config:\n"); dk <- diagv[diagv$config == DEFAULT, ]; print(round(tapply(dk$own - dk$null, list(dk$KO, dk$target), mean)[KOS, TARGETS], 3))
cat("\nbest config per target (diag-off, top1, perm p):\n"); print(format(res[res$config == BEST, c("target", "diag_minus_off", "top1", "perm_p")], digits = 3), row.names = FALSE)

## figures
save_fig("fig6A_matching_benchmark", 14, 8.8, function() {
  layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.4, 0), family = "sans")
  lab <- list(stat = "statistic", feat = "feature ranking", N = "number of genes", centre = "KO-centred", method = "correlation")
  for (f in names(lab)) { par(mar = c(6.5, 4.4, 2.4, 0.8)); lv <- sort(unique(rl[[f]])); lv <- if (f == "feat") c("bulk_var", "screen_signal", "bulk_signal", "both_signal", "all") else lv
    boxplot(lapply(lv, function(v) rl$diag_minus_off[rl[[f]] == v]), names = if (f == "centre") c("uncentred", "KO-centred") else as.character(lv), las = 2, col = "grey85", border = INK2, outline = FALSE, ylab = "diag - off-diag (per line)", cex.axis = 0.75)
    abline(h = 0, lty = 2, col = INK); mtext(lab[[f]], side = 3, adj = 0, line = 0.3, font = 2, cex = 0.85) }
  par(mar = c(6.5, 4.4, 2.4, 0.8)); v <- rbind(picked = lolo$heldout_dmo_picked, default = lolo$heldout_dmo_default); colnames(v) <- sub("NCI-|-211H", "", lolo$held_out)
  barplot(v, beside = TRUE, col = c("#1baf7a", "grey70"), border = NA, ylab = "held-out diag - off-diag", las = 2); abline(h = 0, col = INK)
  legend("topleft", c("config picked on other 3 lines", "default"), fill = c("#1baf7a", "grey70"), border = NA, bty = "n", cex = 0.7)
  mtext("leave-one-line-out", side = 3, adj = 0, line = 0.3, font = 2, cex = 0.85)
  mtext("Which way of matching screen KO signal to bulk KO signal works best? Grid of every analysis choice, scored per bulk line", outer = TRUE, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext(sprintf("%d configurations x 4 lines; metric = mean diagonal - mean off-diagonal of the 8x8 screen-KO x bulk-KO correlation matrix (0 = no KO-specific match)", nrow(grid)), outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.58, col = INK2)
})

mat_for <- function(cf, tg) { g <- res[res$config == cf & res$target == tg, ][1, ]; f <- if (g$feat == "all") seq_along(genes) else rank_by[[g$feat]][seq_len(g$N)]
  S <- S_all[[g$stat]][f, , drop = FALSE]; B <- B_all[[g$stat]][[tg]][f, , drop = FALSE]; if (g$centre) { S <- kcentre(S); B <- kcentre(B) }; cor(S, B, method = g$method) }
save_fig("fig6B_best_config_matrices", 16, 7.6, function() {
  layout(matrix(1:10, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.4, 0), family = "sans")
  for (cf in c(DEFAULT, BEST)) for (tg in TARGETS) { C <- mat_for(cf, tg); r <- res[res$config == cf & res$target == tg, ]; par(mar = c(3.8, 3.8, 3.2, 0.4)); vm <- max(abs(C))
    image(seq_len(K), seq_len(K), t(C)[, K:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(K)) for (j in seq_len(K)) text(j, K - i + 1, sprintf("%.2f", C[i, j]), cex = 0.45, col = INK)
    rect(seq_len(K) - 0.5, K - seq_len(K) + 0.5, seq_len(K) + 0.5, K - seq_len(K) + 1.5, border = "black", lwd = 1.6)
    axis(1, at = seq_len(K), labels = KOS, las = 2, tick = FALSE, cex.axis = 0.55, line = -0.7); axis(2, at = K:1, labels = KOS, las = 1, tick = FALSE, cex.axis = 0.55, line = -0.7); box(col = MUTED, lwd = 0.6)
    mtext(sprintf("%s: %s", if (cf == DEFAULT) "default" else "best", sub("NCI-|-211H", "", tg)), side = 3, adj = 0, line = 1.3, font = 2, cex = 0.72)
    mtext(sprintf("diag-off %+.3f | top1 %d/8 | p=%.3f", r$diag_minus_off, r$top1, r$perm_p), side = 3, adj = 0, line = 0.35, cex = 0.55, col = INK2) }
  mtext(sprintf("Screen KO (rows) x bulk KO (columns) correlation: default [%s] (top) vs best [%s] (bottom)", DEFAULT, BEST), outer = TRUE, side = 3, adj = 0, line = 1.7, font = 2, cex = 0.85)
})
