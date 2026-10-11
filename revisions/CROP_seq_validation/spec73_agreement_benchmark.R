# Spec 73: benchmark of the cross-cell-line agreement of the bulk KO profiles (NTC - KO) over every analysis choice, scored with ONE metric.
#   from spec72 : sample set (S70 | S70N4 | S75) x normalisation (cp10k | tmm | vst | voomw) x covariates (none | cc | cc_pa | cc_pa_depth) x statistic (med | coef | t | t_rest)
#   feature space: Cm20 (AddModuleScore) | Cm20_genes (mean gene statistic over the top-50 genes of each Cm) | Hallmark (mean gene statistic per Hallmark set) |
#                  top 500/1000/2000/3000 variable genes (mean within-line variance) | all expressed genes |
#                  signal500/1000 = genes with the largest KO signal (max |t| over KOs) in the 2 lines NOT in the pair being compared (no circularity) |
#                  screen = union of the top-100 genes (|Welch t|, cc regressed) of each KO in the CROP-seq screen (independent of the bulk data)
#   profile prep : none | centred (each feature minus its mean over the 8 KOs of the line: removes what all KOs share, i.e. the shared-NTC term) | z (centred + scaled)
#   correlation  : pearson | spearman
# Metric (agree_common.R agreement()): per KO mean same-KO correlation over the 6 line pairs vs other-KO correlations; Wilcoxon p and label-permutation p (N_PERM); global =
# mean diagonal - mean off-diagonal (dmo) with permutation p. Overfitting guard: leave-one-line-out -- config chosen on the 3 pairs without line h, scored on the 3 pairs with h.
# needs: spec72_build_profiles.R, spec1_prepare_stats.R (spec_screen_stats.rds)
# outputs: spec73_grid.csv, spec73_grid_per_KO.csv, spec73_LOLO.csv, spec73_axis_effects.csv, spec73_reference_configs_per_KO.csv, fig7D_benchmark_axes.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages(library(parallel)); NC <- as.integer(Sys.getenv("NCORES", unset = "8")); N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "2000"))
SETN <- strsplit(Sys.getenv("SETS", unset = "S70,S70N4,S75"), ",")[[1]]; perms <- make_perms(N_PERM)
OSUF <- if (identical(SETN, c("S70", "S70N4", "S75"))) "" else paste0("_", paste(SETN, collapse = "_")); REFSET <- if ("S70" %in% SETN) "S70" else SETN[1]   # e.g. SETS=RUN3 -> spec73_*_RUN3.csv
top50 <- cm_top50(); ss <- readRDS(file.path(OUT, "spec_screen_stats.rds"))
hall <- { suppressPackageStartupMessages(library(msigdbr)); h <- msigdbr(species = "Homo sapiens", category = "H"); split(h$gene_symbol, sub("HALLMARK_", "", h$gs_name)) }
PREPS <- list(none = c(FALSE, FALSE), centred = c(TRUE, FALSE), z = c(TRUE, TRUE)); METHODS <- c("pearson", "spearman")

eval_set <- function(sn) { pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_%s%s.rds", sn, DSUF))); genes <- pr$genes
  vord <- genes[order(rowMeans(pr$gene_var), decreasing = TRUE)]
  scr <- intersect(unique(unlist(lapply(KOS, function(k) ss$genes[order(abs(ss$tcc[, k]), decreasing = TRUE)[1:100]]))), genes)
  sig_pair <- lapply(PAIRS, function(p) genes[order(rowMeans(pr$gene_signal[, setdiff(LINES, p), drop = FALSE]), decreasing = TRUE)])
  setmean <- function(A, sets) { sets <- lapply(sets, intersect, genes); sets <- sets[sapply(sets, length) >= 10]; out <- array(NA_real_, c(length(sets), dim(A)[2:3]), dimnames = c(list(names(sets)), dimnames(A)[2:3]))
    for (i in seq_along(sets)) out[i, , ] <- colMeans(A[sets[[i]], , , drop = FALSE]); out }
  cfgs <- expand.grid(norm = names(pr$P), cov = names(pr$P[[1]]), stat = names(pr$P[[1]][[1]]), stringsAsFactors = FALSE)
  res <- mclapply(seq_len(nrow(cfgs)), function(i) { cf <- cfgs[i, ]; A <- pr$P[[cf$norm]][[cf$cov]][[cf$stat]]; Mo <- pr$M[[cf$norm]][[cf$cov]][[cf$stat]]
    spaces <- list(Cm20 = function(pi) Mo, Cm20_genes = local({ x <- setmean(A, top50); function(pi) x }), Hallmark = local({ x <- setmean(A, hall); function(pi) x }),
                   top500 = function(pi) A[vord[1:500], , ], top1000 = function(pi) A[vord[1:1000], , ], top2000 = function(pi) A[vord[1:2000], , ], top3000 = function(pi) A[vord[1:3000], , ],
                   all_genes = function(pi) A, signal500 = function(pi) A[sig_pair[[pi]][1:500], , ], signal1000 = function(pi) A[sig_pair[[pi]][1:1000], , ], screen = function(pi) A[scr, , ])
    g <- list(); k <- list()
    for (sp in names(spaces)) for (pp in names(PREPS)) for (me in METHODS) {
      CL <- lapply(seq_along(PAIRS), function(pi) { X <- spaces[[sp]](pi); p <- PAIRS[[pi]]; C <- pair_cor(prep(X[, , p[1]], PREPS[[pp]][1], PREPS[[pp]][2]), prep(X[, , p[2]], PREPS[[pp]][1], PREPS[[pp]][2]), me); dimnames(C) <- list(KOS, KOS); C })
      a <- agreement(CL, perms); pd <- sapply(CL, function(C) mean(diag(C), na.rm = TRUE) - mean(C[row(C) != col(C)], na.rm = TRUE)); names(pd) <- paste0("dmo_", PLAB)
      id <- data.frame(set = sn, cf, space = sp, prep = pp, method = me, stringsAsFactors = FALSE)
      g[[length(g) + 1]] <- cbind(id, t(a$global), n_KO_perm05 = sum(a$per_KO$perm_p < 0.05), n_KO_fdr10 = sum(a$per_KO$perm_fdr < 0.1), n_KO_wilcox05 = sum(a$per_KO$wilcox_p < 0.05), t(pd))
      k[[length(k) + 1]] <- cbind(id, a$per_KO, row.names = NULL) }
    list(g = do.call(rbind, g), k = do.call(rbind, k)) }, mc.cores = NC)
  bad <- sapply(res, inherits, "try-error"); if (any(bad)) stop(res[[which(bad)[1]]])
  list(g = do.call(rbind, lapply(res, `[[`, "g")), k = do.call(rbind, lapply(res, `[[`, "k"))) }
R <- lapply(SETN, eval_set); grid <- do.call(rbind, lapply(R, `[[`, "g")); perko <- do.call(rbind, lapply(R, `[[`, "k")); rownames(grid) <- NULL
grid$config <- with(grid, paste(set, norm, cov, stat, space, prep, method, sep = " | ")); perko$config <- with(perko, paste(set, norm, cov, stat, space, prep, method, sep = " | "))
write.csv(grid, file.path(OUT, paste0("spec73_grid", OSUF, ".csv")), row.names = FALSE); write.csv(perko, file.path(OUT, paste0("spec73_grid_per_KO", OSUF, ".csv")), row.names = FALSE)

## ---------------------------------------------------------------- marginal effect of every axis (median dmo / perm z over all configs sharing the level)
AX <- c("set", "norm", "cov", "stat", "space", "prep", "method")
axis_eff <- do.call(rbind, lapply(AX, function(a) do.call(rbind, lapply(split(grid, grid[[a]]), function(d) data.frame(axis = a, level = d[[a]][1], n_configs = nrow(d), median_dmo = median(d$dmo), median_perm_z = median(d$perm_z),
  frac_p05 = mean(d$perm_p < 0.05), median_n_KO_perm05 = median(d$n_KO_perm05), median_off_mean = median(d$off_mean))))))
write.csv(axis_eff, file.path(OUT, paste0("spec73_axis_effects", OSUF, ".csv")), row.names = FALSE)

## ---------------------------------------------------------------- leave-one-line-out selection (within S70)
pc <- grep("^dmo_", names(grid), value = TRUE); lolo <- do.call(rbind, lapply(SETN, function(sn) do.call(rbind, lapply(LINES, function(h) { d <- grid[grid$set == sn, ]; short <- sub("NCI-|-211H", "", h)
  with_h <- pc[sapply(strsplit(sub("dmo_", "", pc), "-"), function(x) short %in% x)]; tr <- rowMeans(d[, setdiff(pc, with_h)]); te <- rowMeans(d[, with_h]); b <- which.max(tr)
  leg <- which(d$norm == "cp10k" & d$cov == "none" & d$stat == "med" & d$space == "top1000" & d$prep == "none" & d$method == "spearman")
  pri <- which(d$norm == "voomw" & d$cov == "cc_pa" & d$stat == "t" & d$space == "top2000" & d$prep == "centred" & d$method == "spearman")
  data.frame(set = sn, held_out = h, picked = d$config[b], train_dmo_picked = tr[b], heldout_dmo_picked = te[b], heldout_dmo_legacy = te[leg], heldout_dmo_principled = te[pri], heldout_dmo_median_config = median(te),
             heldout_rank_pct_of_picked = 100 * mean(te <= te[b])) }))))
write.csv(lolo, file.path(OUT, paste0("spec73_LOLO", OSUF, ".csv")), row.names = FALSE)

## ---------------------------------------------------------------- reference configurations, per KO
refs <- rbind(data.frame(name = "legacy (spec68)", norm = "cp10k", cov = "none", stat = "med", prep = "none", method = "spearman"),
              data.frame(name = "principled", norm = "voomw", cov = "cc_pa", stat = "t", prep = "centred", method = "spearman"))
ref_k <- do.call(rbind, lapply(seq_len(nrow(refs)), function(i) { r <- refs[i, ]; d <- perko[perko$norm == r$norm & perko$cov == r$cov & perko$stat == r$stat & perko$prep == r$prep & perko$method == r$method & perko$space %in% c("Cm20", "top1000", "top2000", "all_genes", "screen"), ]
  cbind(ref = r$name, d) })); write.csv(ref_k, file.path(OUT, paste0("spec73_reference_configs_per_KO", OSUF, ".csv")), row.names = FALSE)

options(width = 250); cat(sprintf("\n%d configurations scored (%d permutations each)\n", nrow(grid), N_PERM))
cat("\n=== marginal effect of each axis ===\n"); print(axis_eff, digits = 3, row.names = FALSE)
cat(sprintf("\n=== top 25 configurations by global permutation z (%s only) ===\n", REFSET)); d <- grid[grid$set == REFSET, ]; print(head(d[order(-d$perm_z), c("config", "diag_mean", "off_mean", "dmo", "top1", "perm_p", "perm_z", "n_KO_perm05", "n_KO_fdr10")], 25), digits = 3, row.names = FALSE)
cat("\n=== reference configurations: global ===\n"); for (i in seq_len(nrow(refs))) { r <- refs[i, ]; d <- grid[grid$norm == r$norm & grid$cov == r$cov & grid$stat == r$stat & grid$prep == r$prep & grid$method == r$method, ]
  cat("\n", r$name, "\n"); print(d[order(d$set, d$space), c("set", "space", "diag_mean", "off_mean", "dmo", "top1", "perm_p", "perm_z", "n_KO_perm05", "n_KO_fdr10")], digits = 3, row.names = FALSE) }
cat("\n=== leave-one-line-out ===\n"); print(lolo, digits = 3, row.names = FALSE)
if ("S70" %in% SETN) { old <- read.csv(file.path(OUT, "spec68_bulk_pairwise_NTC_KO_stats.csv")); new <- ref_k[ref_k$ref == "legacy (spec68)" & ref_k$set == "S70" & ref_k$space %in% c("Cm20", "top1000", "top2000"), ]
old$space2 <- c("20 Cm modules" = "Cm20", "top 1000 variable genes" = "top1000", "top 2000 variable genes" = "top2000", "top 3000 variable genes" = "top3000")[old$space]
mm <- merge(old, new, by.x = c("space2", "KO"), by.y = c("space", "KO")); cat(sprintf("\nregression check vs spec68 (legacy config; gene filter and module-score background differ slightly): cor(median intra) = %.3f, cor(sep) = %.3f, max |diff median intra| = %.3f\n",
  cor(mm$median_intra.x, mm$median_intra.y), cor(mm$sep.x, mm$sep.y), max(abs(mm$median_intra.x - mm$median_intra.y)))) }

## ---------------------------------------------------------------- figure: distribution of the global permutation z per level of every axis
save_fig(paste0("fig7D_benchmark_axes", OSUF), 15, 7.6, function() {
  layout(matrix(1:8, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.2, 0.5), mar = c(7, 4.2, 2.2, 0.6), family = "sans")
  for (a in AX) { lv <- names(sort(tapply(grid$perm_z, grid[[a]], median), decreasing = TRUE)); plot(NA, xlim = c(0.4, length(lv) + 0.6), ylim = range(grid$perm_z), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "global agreement (permutation z)", cex.axis = 0.75, cex.lab = 0.8)
    abline(h = 0, col = INK, lty = 2); abline(h = qnorm(0.95), col = MUTED, lty = 3)
    for (i in seq_along(lv)) boxplot(grid$perm_z[grid[[a]] == lv[i]], at = i, add = TRUE, boxwex = 0.6, col = adjustcolor("#2a78d6", 0.35), border = INK2, outline = FALSE, axes = FALSE)
    text(seq_along(lv), par("usr")[3] - 0.05 * diff(par("usr")[3:4]), lv, srt = 45, adj = 1, xpd = NA, cex = 0.75); mtext(a, side = 3, adj = 0, line = 0.4, font = 2, cex = 0.85) }
  d <- grid[grid$set == REFSET, ]; plot(d$off_mean, d$diag_mean, pch = 16, cex = 0.4, col = adjustcolor(ifelse(d$prep == "none", "#eb6834", "#2a78d6"), 0.35), bty = "n", las = 1, xlab = "mean other-KO correlation (off-diagonal)", ylab = "mean same-KO correlation (diagonal)", cex.axis = 0.75, cex.lab = 0.8)
  abline(0, 1, col = INK, lty = 2); legend("topleft", c("NTC - KO as is", "centred / z"), col = c("#eb6834", "#2a78d6"), pch = 16, bty = "n", cex = 0.7); mtext(sprintf("shared component (%s)", REFSET), side = 3, adj = 0, line = 0.4, font = 2, cex = 0.85)
  mtext(sprintf("Benchmark of the cross-line agreement of bulk KO profiles (%s): %d configurations", DLAB, nrow(grid)), outer = TRUE, side = 3, adj = 0, line = 1.7, font = 2, cex = 1)
  mtext("each box = all configurations sharing that level; z = (mean same-KO correlation - permutation mean) / permutation sd, KO labels permuted within every line; dotted = one-sided 5%", outer = TRUE, side = 3, adj = 0, line = 0.4, cex = 0.62, col = INK2)
})
