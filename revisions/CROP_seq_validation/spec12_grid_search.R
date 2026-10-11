# Spec 12: exhaustive exploratory search for the strongest bulk-KO vs CROP-seq-KO correspondence (KNOWN TO BE OPTIMISTIC; a search-wide null is provided).
# Dimensions searched:
#   feature space : 20 Cm modules | screen's variable genes | bulk-derived top-3000 variable genes | top-3000 expressed genes | whole transcriptome (all expressed genes)
#   profile / normalisation : KO-NTC log2 CP10K difference | + cell-cycle regression | DESeq2 VST difference | + cc | pooled-variance t statistic | + cc | absolute expression (no NTC) | + cc
#   scaling (differences) : none | equal-RMS ; centring across KOs : no | yes ; correlation : Pearson | Spearman
#   replicate handling : average of replicates | ONLY the replicate best matching the screen KO (selection-biased; see control)
#   TF removal : every subset of the 8 shared TFs with >= 4 TFs ; cell-line removal : every non-empty subset of the 4 lines (>= 12 bulk rows overall)
# Score = mean over bulk rows of the matched screen KO's score after standardising each row across the selected screen KOs (z; 0 = no better than the other KOs).
# Search-wide null: the SAME search on bulk rows with randomly permuted KO labels; the max z over the whole grid is recorded for each permutation.
# needs: spec1_prepare_stats.R
# outputs: spec12_full_set.csv, spec12_top_results.csv, spec12_null_max.csv, fig15_grid_search.{png,pdf}, fig16_best_combination_matrix.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_variants.R"))
suppressPackageStartupMessages({ library(DESeq2); library(parallel) })
set.seed(12); NPERM <- as.integer(Sys.getenv("SPEC12_NPERM", "100")); MINROWS <- 12; MINTF <- 4; NCORE <- 6
ctx <- load_spec_context(LINES); ss <- ctx$ss; lg <- ctx$lg; meta <- ctx$meta; ccb <- ctx$ccb; counts <- load_bulk()$counts
genes_all <- ctx$F_all; ALLC <- c(ORDER, EXTRA)

## ---------------------------------------------------------------- feature spaces
lv <- sapply(LINES, function(cl) matrixStats::rowVars(lg[genes_all, rownames(meta)[meta$Cell_Line == cl]])); rownames(lv) <- genes_all
F_sets <- list(all_genes = genes_all, var_screen = ctx$F_var, var_bulk = names(sort(rowMeans(lv), decreasing = TRUE))[1:3000],
               top3000_expr = genes_all[order(ss$abs0[genes_all, "NTC"], decreasing = TRUE)][1:3000])
cat("feature sets:", paste(names(F_sets), sapply(F_sets, length), collapse = " | "), "\n")

## ---------------------------------------------------------------- screen profiles (genes from spec1; modules recomputed from the screen object)
sc_df <- read.csv(file.path(PP, "cm_sample_scores_AddModuleScore_bulk_RUN1.csv"), check.names = FALSE); cm_names <- setdiff(colnames(sc_df), c("sample", "cell_line"))
mods <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scrna/scRNA_PM_modules.rds"); gl <- lapply(mods[cm_names], head, 50)
srt <- load_screen(); srt <- AddModuleScore(srt, features = gl); cl_cols <- grep("^Cluster", colnames(srt@meta.data), value = TRUE)
Mscr <- t(as.matrix(srt@meta.data[, cl_cols])); rownames(Mscr) <- cm_names
Smod <- list(reg0 = screen_ko_stats(Mscr, srt$merged_call, srt$cc, FALSE), regcc = screen_ko_stats(Mscr, srt$merged_call, srt$cc, TRUE)); rm(srt)
TYPE_SRC <- c(diff = "diff0", diff_cc = "diffcc", vst = "diff0", vst_cc = "diffcc", t = "t0", t_cc = "tcc", abs = "abs0", abs_cc = "abscc")
screen_mat <- function(space, type) {
  if (space == "modules") { s <- if (type %in% c("diff_cc", "vst_cc", "t_cc", "abs_cc")) Smod$regcc else Smod$reg0
    m <- switch(sub("_cc$", "", type), diff = s$diff, t = s$t, abs = s$abs, vst = s$diff) } else m <- ss[[TYPE_SRC[[type]]]][F_sets[[space]], ]
  m[, colnames(m) != "NTC", drop = FALSE] }

## ---------------------------------------------------------------- bulk replicate-level profiles per line
make_types <- function(y, cc, grp, vstm = NULL) {
  ko <- grp != "NTC"; ntc <- grp == "NTC"; nn <- sum(ntc); ccc <- cc - mean(cc)
  reg <- function(M) { b <- as.numeric((M - rowMeans(M)) %*% ccc) / sum(ccc^2); M - outer(b, ccc) }
  D <- function(M) M[, ko, drop = FALSE] - rowMeans(M[, ntc, drop = FALSE])
  tstat <- function(M, with_cc) { g <- factor(grp, levels = c("NTC", sort(setdiff(unique(grp), "NTC")))); X <- if (with_cc) model.matrix(~ cc + g) else model.matrix(~ g)
    q <- qr(X); res <- qr.resid(q, t(M)); s2 <- pmax(colSums(res^2) / (ncol(M) - q$rank), 1e-4); Ma <- if (with_cc) M - outer(qr.coef(q, t(M))["cc", ], ccc) else M
    (Ma[, ko, drop = FALSE] - rowMeans(Ma[, ntc, drop = FALSE])) / sqrt(s2 * (1 + 1 / nn)) }
  out <- list(diff = D(y), diff_cc = D(reg(y)), t = tstat(y, FALSE), t_cc = tstat(y, TRUE), abs = y[, ko, drop = FALSE], abs_cc = reg(y)[, ko, drop = FALSE])
  if (!is.null(vstm)) { out$vst <- D(vstm); out$vst_cc <- D(reg(vstm)) }
  lapply(out, function(m) { m[!is.finite(m)] <- 0; m }) }
B <- list()
for (cl in LINES) { m <- meta[meta$Cell_Line == cl, ]; s <- rownames(m); cc <- ccb$cc[match(s, ccb$sample)]; grp <- as.character(m$Gene_Control)
  x <- round(counts[, s]); storage.mode(x) <- "integer"; x <- x[rowSums(x >= 10) >= 2, ]; vs <- vst(x, blind = TRUE)[genes_all, ]
  B[[cl]] <- list(gene = make_types(lg[genes_all, s], cc, grp, vs), modules = make_types(t(as.matrix(sc_df[match(s, sc_df$sample), cm_names])) |> `colnames<-`(s), cc, grp), ko_samples = s[grp != "NTC"], ko = grp[grp != "NTC"]) }
cat("bulk profiles built for", length(LINES), "lines\n")

## ---------------------------------------------------------------- configurations
cfgs <- list()
for (space in c("modules", names(F_sets))) { types <- if (space == "modules") c("diff", "diff_cc", "t", "t_cc", "abs", "abs_cc") else names(TYPE_SRC)
  for (ty in types) for (scl in if (ty %in% c("diff", "diff_cc", "vst", "vst_cc")) c("none", "rms") else "none") for (cen in c(FALSE, TRUE)) for (met in c("pearson", "spearman"))
    cfgs[[length(cfgs) + 1]] <- data.frame(space = space, type = ty, scale = scl, centred = cen, cor = met) }
cfgs <- do.call(rbind, cfgs); cfgs$id <- seq_len(nrow(cfgs)); cat("configurations:", nrow(cfgs), "\n")
prep_mat <- function(M, centre, scale) { if (centre) M <- M - rowMeans(M); if (scale == "rms") M <- M / rms_scale(M); M[!is.finite(M)] <- 0; M }
samp <- do.call(rbind, lapply(LINES, function(cl) { k <- B[[cl]]$ko; data.frame(line = cl, ko = k, sample = B[[cl]]$ko_samples)[k %in% ORDER, ] })); rownames(samp) <- NULL
grows <- unique(samp[, c("line", "ko")]); rownames(grows) <- NULL; grows$key <- paste(grows$line, grows$ko); samp$key <- paste(samp$line, samp$ko)
grp_idx <- lapply(grows$key, function(k) which(samp$key == k)); ko_idx <- match(grows$ko, ORDER)
run_config <- function(i) { cf <- cfgs[i, ]; feats <- if (cf$space == "modules") cm_names else F_sets[[cf$space]]; src <- if (cf$space == "modules") "modules" else "gene"
  S <- prep_mat(screen_mat(cf$space, cf$type), cf$centred, cf$scale)[, ALLC]; if (cf$cor == "spearman") S <- apply(S, 2, rank)
  Rrep <- matrix(NA, nrow(samp), length(ALLC)); Cavg <- matrix(NA, nrow(grows), length(ALLC))
  for (cl in LINES) { M <- B[[cl]][[src]][[cf$type]]; if (src == "gene") M <- M[feats, , drop = FALSE]; M <- prep_mat(M, cf$centred, cf$scale)
    idx <- which(samp$line == cl); Mr <- M[, samp$sample[idx], drop = FALSE]; Rrep[idx, ] <- cor(if (cf$cor == "spearman") apply(Mr, 2, rank) else Mr, S)
    gi <- which(grows$line == cl); G <- sapply(gi, function(g) rowMeans(M[, samp$sample[grp_idx[[g]]], drop = FALSE])); Cavg[gi, ] <- cor(if (cf$cor == "spearman") apply(G, 2, rank) else G, S) }
  list(R = Rrep, Cavg = Cavg) }
t0 <- Sys.time(); RES <- mclapply(seq_len(nrow(cfgs)), run_config, mc.cores = NCORE); cat("profiles/correlations done in", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
stopifnot(!any(sapply(RES, function(r) inherits(r, "try-error"))))

## ---------------------------------------------------------------- grid evaluation
TS <- unlist(lapply(MINTF:8, function(k) combn(8, k, simplify = FALSE)), recursive = FALSE)              # TF subsets (column indexes into ORDER)
subsets_lines <- unlist(lapply(1:4, function(k) combn(LINES, k, simplify = FALSE)), recursive = FALSE)
LS <- sapply(subsets_lines, function(ls) as.numeric(grows$line %in% ls))                                     # rows x line-subsets
best_rows <- function(R, match_col) t(sapply(seq_along(grp_idx), function(g) { ix <- grp_idx[[g]]; R[ix[which.max(R[ix, match_col[g]])], ] }))
ctl_rows <- function(R) t(sapply(seq_along(grp_idx), function(g) apply(R[grp_idx[[g]], , drop = FALSE], 2, max)))
grid_eval <- function(C, lab) {          # C rows x 8 ; lab = label index (1..8) each row is treated as
  Z <- T1 <- N <- matrix(NA, ncol(LS), length(TS))
  for (t in seq_along(TS)) { cols <- TS[[t]]; sel <- lab %in% cols; if (sum(sel) < MINROWS) next
    sub <- C[sel, cols, drop = FALSE]; mu <- rowMeans(sub); sdv <- matrixStats::rowSds(sub); j <- match(lab[sel], cols); mv <- sub[cbind(seq_len(nrow(sub)), j)]
    mz <- (mv - mu) / sdv; mz[!is.finite(mz)] <- 0; top <- as.numeric(mv >= matrixStats::rowMaxs(sub) - 1e-12)
    zf <- numeric(nrow(C)); zf[sel] <- mz; tf <- numeric(nrow(C)); tf[sel] <- top; cnt <- as.numeric(sel); n <- as.numeric(crossprod(LS, cnt)); ok <- n >= MINROWS
    Z[ok, t] <- (as.numeric(crossprod(LS, zf)) / n)[ok]; T1[ok, t] <- (as.numeric(crossprod(LS, tf)) / n)[ok]; N[ok, t] <- n[ok] }
  list(Z = Z, T1 = T1, N = N) }
mode_C <- function(r, mode, match_col) switch(mode, avg = r$Cavg, best = best_rows(r$R, match_col), ctl = ctl_rows(r$R))
full_ts <- length(TS); full_ls <- ncol(LS)                                                                      # index of the all-TF subset and the all-line subset
t0 <- Sys.time(); obs <- list(); k <- 0
for (i in seq_len(nrow(cfgs))) for (mode in c("avg", "best")) { g <- grid_eval(mode_C(RES[[i]], mode, ko_idx), ko_idx); k <- k + 1; obs[[k]] <- c(list(cfg = i, mode = mode), g) }
cat("observed grid:", length(obs), "config-modes x", ncol(LS), "line subsets x", length(TS), "TF subsets in", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")

## observed results (long table of all valid cells is large; keep the full-set table and the top results)
full <- do.call(rbind, lapply(obs, function(o) cbind(cfgs[o$cfg, ], mode = o$mode, z_full = o$Z[full_ls, full_ts], top1_full = o$T1[full_ls, full_ts], z_noH28 = o$Z[match("MSTO-211H|NCI-H2052|NCI-H2452", sapply(subsets_lines, function(x) paste(sort(x, method = "radix"), collapse = "|"))), full_ts])))
write.csv(full, file.path(OUT, "spec12_full_set.csv"), row.names = FALSE)
long <- do.call(rbind, lapply(obs, function(o) { ii <- which(!is.na(o$Z), arr.ind = TRUE); data.frame(cfg = o$cfg, mode = o$mode, ls = ii[, 1], ts = ii[, 2], z = o$Z[ii], top1 = o$T1[ii], n_rows = o$N[ii]) }))
long <- long[order(-long$z), ]; cat("valid grid cells:", nrow(long), "\n")
top <- head(long, 200); top$description <- with(cfgs[top$cfg, ], sprintf("%s | %s | scale %s | %s | %s", space, type, scale, ifelse(centred, "centred", "uncentred"), cor))
top$lines <- sapply(top$ls, function(i) paste(sub("NCI-", "", subsets_lines[[i]]), collapse = "+")); top$TFs <- sapply(top$ts, function(i) paste(ORDER[TS[[i]]], collapse = ","))
write.csv(top, file.path(OUT, "spec12_top_results.csv"), row.names = FALSE)

## ---------------------------------------------------------------- search-wide null: same search with permuted KO labels
null_one <- function(p) { set.seed(1000 + p); perm <- sample(8); lab <- perm[ko_idx]; mx <- c(all = -Inf, avg = -Inf, best = -Inf, full_sets = -Inf)
  for (i in seq_len(nrow(cfgs))) for (mode in c("avg", "best")) { g <- grid_eval(mode_C(RES[[i]], mode, lab), lab); z <- g$Z; m <- max(z, na.rm = TRUE)
    mx["all"] <- max(mx["all"], m); mx[mode] <- max(mx[mode], m); if (!is.na(z[full_ls, full_ts])) mx["full_sets"] <- max(mx["full_sets"], z[full_ls, full_ts]) }
  mx }
t0 <- Sys.time(); nl <- mclapply(seq_len(NPERM), null_one, mc.cores = NCORE); nullmax <- do.call(rbind, nl); cat("search-wide null:", NPERM, "permutations in", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
write.csv(nullmax, file.path(OUT, "spec12_null_max.csv"), row.names = FALSE)
obs_max <- c(all = max(long$z), avg = max(long$z[long$mode == "avg"]), best = max(long$z[long$mode == "best"]), full_sets = max(full$z_full, na.rm = TRUE))
fw <- sapply(names(obs_max), function(n) (1 + sum(nullmax[, n] >= obs_max[n])) / (NPERM + 1))
cat("\n=== search-wide comparison: best z observed vs the max z the same search finds on permuted labels ===\n"); print(round(rbind(observed_max = obs_max, null_median = apply(nullmax, 2, median), null_95th = apply(nullmax, 2, quantile, 0.95), family_wise_p = fw), 3))

## ---------------------------------------------------------------- summaries
options(width = 250); cat("\n=== 25 best combinations (>=", MINTF, "TFs, >=", MINROWS, "rows) ===\n")
print(format(top[1:25, c("mode", "description", "lines", "TFs", "n_rows", "z", "top1")], digits = 3), row.names = FALSE)
cat("\n=== best combination with ALL TFs and ALL lines (no subsetting) ===\n"); f <- full[order(-full$z_full), ][1:10, ]; f$description <- with(f, sprintf("%s | %s | scale %s | %s | %s", space, type, scale, ifelse(centred, "centred", "uncentred"), cor)); print(format(f[, c("mode", "description", "z_full", "top1_full", "z_noH28")], digits = 3), row.names = FALSE)
cat("\n=== main effects (mean z over all configs; ALL TFs, ALL lines) ===\n")
for (d in c("space", "type", "scale", "centred", "cor", "mode")) { a <- tapply(full$z_full, full[[d]], mean); cat(sprintf("%-8s", d), paste(sprintf("%s %.3f", names(a), a), collapse = " | "), "\n") }
saveRDS(list(cfgs = cfgs, RES = RES, grows = grows, samp = samp, obs = obs, top = top, subsets_lines = subsets_lines, TS = TS, nullmax = nullmax, obs_max = obs_max, fw = fw, full = full), file.path(OUT, "spec12_grid_search.rds"))

## ---------------------------------------------------------------- exact permutation p and selection-bias control for the best hits
top_u <- top[!duplicated(paste(top$mode, round(top$z, 6), top$ls, top$ts)), ]
cfg_desc <- function(i) with(cfgs[i, ], sprintf("%s | %s | scale %s | %s | %s", space, type, scale, ifelse(centred, "centred", "uncentred"), cor))
detail <- function(cfg, mode, ls_i, ts_i) { cols <- ORDER[TS[[ts_i]]]; sel <- grows$line %in% subsets_lines[[ls_i]] & grows$ko %in% cols; r <- RES[[cfg]]
  C <- mode_C(r, mode, ko_idx); ct <- mode_C(r, "ctl", ko_idx); ev <- function(M) { colnames(M) <- ALLC; evaluate_specificity(M[sel, , drop = FALSE], grows$ko[sel], cols = cols) }
  list(sel = sel, cols = cols, C = C, obs = ev(C), avg = ev(r$Cavg), ctl = ev(ct)) }
b1 <- top_u[1, ]; d1 <- detail(b1$cfg, b1$mode, b1$ls, b1$ts)
fb <- full[order(-full$z_full), ][1, ]; fcfg <- which(cfgs$space == fb$space & cfgs$type == fb$type & cfgs$scale == fb$scale & cfgs$centred == fb$centred & cfgs$cor == fb$cor); d2 <- detail(fcfg, fb$mode, full_ls, full_ts)
cat(sprintf("\nBEST OVERALL: %s [%s] lines %s TFs %s: z=%.2f, exact permutation p (this subset only) = %.4f; replicates averaged z=%.2f; selection-bias control z=%.2f\n", cfg_desc(b1$cfg), b1$mode, b1$lines, b1$TFs, d1$obs$z_matched, d1$obs$perm_p, d1$avg$z_matched, d1$ctl$z_matched))
cat(sprintf("BEST WITH ALL TFs AND LINES: %s [%s]: z=%.2f, exact permutation p = %.4f; averaged z=%.2f; control z=%.2f\n", cfg_desc(fcfg), fb$mode, d2$obs$z_matched, d2$obs$perm_p, d2$avg$z_matched, d2$ctl$z_matched))

## ---------------------------------------------------------------- figure 15: search summary
save_fig("fig15_grid_search", 16, 9.6, function() {
  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE), widths = c(1.7, 1), heights = c(1, 0.9)); par(oma = c(0, 0, 4, 0.5), family = "sans")
  dims <- c("space", "type", "scale", "centred", "cor", "mode"); lev <- lapply(dims, function(d) sort(unique(as.character(full[[d]])))); n <- sum(lengths(lev)) + length(dims) - 1
  par(mar = c(3.6, 9.5, 3.4, 0.8)); plot(NA, xlim = c(-0.6, 1.0), ylim = c(n + 0.5, 0.5), axes = FALSE, xlab = "", ylab = ""); abline(v = 0, col = INK2); y <- 0; cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7")
  for (di in seq_along(dims)) { for (l in lev[[di]]) { y <- y + 1; v <- full$z_full[as.character(full[[dims[di]]]) == l]; points(v, y + runif(length(v), -0.28, 0.28), pch = 16, cex = 0.35, col = adjustcolor(INK2, 0.35))
      points(mean(v), y, pch = 23, bg = cols[di], col = INK, cex = 1.5); text(-0.62, y, sprintf("%s: %s", dims[di], l), adj = 1, xpd = NA, cex = 0.65) }; y <- y + 1; abline(h = y - 0.5, col = GRID, lwd = 0.6) }
  axis(1, cex.axis = 0.7, lwd = 0.5); mtext("z (all TFs, all 4 lines)", side = 1, line = 2.2, cex = 0.7)
  mtext("Main effects across all configurations", side = 3, adj = 0, line = 1, font = 2, cex = 0.85); mtext("diamond = mean; dots = individual configurations", side = 3, adj = 0, line = 0.1, cex = 0.55, col = INK2)
  tp <- top_u[1:20, ]; par(mar = c(3.4, 30, 3.4, 1.5)); plot(NA, xlim = c(0, max(tp$z, quantile(nullmax[, "all"], 0.95)) * 1.05), ylim = c(20.5, 0.5), axes = FALSE, xlab = "", ylab = "")
  rect(0, seq_len(20) - 0.36, tp$z, seq_len(20) + 0.36, col = ifelse(tp$mode == "best", "#eb6834", "#2a78d6"), border = NA); abline(v = c(median(nullmax[, "all"]), quantile(nullmax[, "all"], 0.95)), lty = c(2, 3), col = INK)
  axis(1, cex.axis = 0.7, lwd = 0.5); axis(2, at = 1:20, labels = sprintf("%s | %s | %s", ifelse(tp$mode == "best", "best rep", "avg"), sapply(tp$cfg, function(i) sub(" \\| scale none", "", cfg_desc(i))), paste(tp$lines, tp$TFs, sep = ": ")), las = 1, tick = FALSE, cex.axis = 0.5, line = -0.5)
  mtext("Top 20 combinations (any TF / line subset)", side = 3, adj = 0, line = 1, font = 2, cex = 0.85); mtext(sprintf("dashed / dotted = median / 95th percentile of the best z the same search finds on permuted labels (%d permutations)", NPERM), side = 3, adj = 0, line = 0.1, cex = 0.55, col = INK2)
  par(mar = c(3.6, 4, 3.4, 0.8)); rg <- range(c(nullmax, obs_max)); brk <- seq(floor(rg[1] * 10) / 10, ceiling(rg[2] * 10) / 10, by = 0.05)
  h1 <- hist(nullmax[, "all"], breaks = brk, plot = FALSE); h2 <- hist(nullmax[, "full_sets"], breaks = brk, plot = FALSE); plot(NA, xlim = range(brk), ylim = c(0, max(h1$counts, h2$counts) * 1.15), xlab = "", ylab = "permutations", bty = "n", cex.axis = 0.7, cex.lab = 0.75)
  rect(head(brk, -1), 0, brk[-1], h1$counts, col = adjustcolor(GREY, 0.9), border = "white"); rect(head(brk, -1), 0, brk[-1], h2$counts, col = adjustcolor("#2a78d6", 0.6), border = "white")
  abline(v = obs_max["all"], lwd = 2, col = "#b3400f"); abline(v = obs_max["full_sets"], lwd = 2, col = "#2a78d6"); mtext("best z found by the search", side = 1, line = 2.2, cex = 0.7)
  legend("topleft", c(sprintf("any subset: observed %.2f (family-wise p = %.2f)", obs_max["all"], fw["all"]), sprintf("all TFs & lines: observed %.2f (p = %.2f)", obs_max["full_sets"], fw["full_sets"]), "grey = null, any subset", "blue = null, all TFs & lines"), lty = c(1, 1, NA, NA), lwd = 2, col = c("#b3400f", "#2a78d6", NA, NA), fill = c(NA, NA, GREY, adjustcolor("#2a78d6", 0.6)), border = NA, bty = "n", cex = 0.65, x.intersp = 0.5)
  mtext("Is the best hit better than what the search finds by chance?", side = 3, adj = 0, line = 0.9, font = 2, cex = 0.85)
  mtext("Exploratory search for the strongest bulk-vs-CROP-seq KO correspondence", outer = TRUE, side = 3, adj = 0, line = 2.2, font = 2, cex = 1.05)
  mtext(sprintf("%d configurations x 2 replicate modes x %d TF subsets x %d line subsets; z = matched screen KO's row-standardised score (0 = no better than other KOs)", nrow(cfgs), length(TS), ncol(LS)), outer = TRUE, side = 3, adj = 0, line = 0.9, cex = 0.62, col = INK2)
})

## ---------------------------------------------------------------- figure 16: matrices of the best combinations
draw_best <- function(name, d, cfg, mode, lines_txt, sub_txt, subset_tfs) {
  lines_in <- unique(grows$line[d$sel]); sel <- grows$line %in% lines_in                                   # all 8 TFs of the selected lines
  Call <- d$C[sel, , drop = FALSE]; colnames(Call) <- ALLC; lab <- grows[sel, ]; jm <- match(lab$ko, ORDER)
  Z <- t(scale(t(Call[, ORDER]))); Z[!is.finite(Z)] <- 0; tf_perf <- sapply(ORDER, function(k) mean(Z[cbind(seq_len(nrow(Z)), jm)][lab$ko == k]))      # matched row-z per TF
  cols <- names(sort(tf_perf, decreasing = TRUE)); lab$ko <- factor(lab$ko, cols); o <- order(match(lab$line, lines_in), as.integer(lab$ko)); lab <- lab[o, ]; lab$ko <- as.character(lab$ko)
  C <- Call[o, cols, drop = FALSE]; Cx <- Call[o, EXTRA, drop = FALSE]; ev <- evaluate_specificity(C, lab$ko, cols = cols); part <- length(subset_tfs) < length(ORDER)
  rl <- if (mode == "best") sapply(seq_len(nrow(C)), function(i) { g <- which(grows$key == lab$key[i]); ix <- grp_idx[[g]]; sprintf("%s (%s)", lab$ko[i], sub(".*_", "", samp$sample[ix[which.max(RES[[cfg]]$R[ix, ko_idx[g]])]])) }) else lab$ko
  draw_ordered_matrix(name, C, Cx, lab, ss$n, rl, cfg_desc(cfg),
    sprintf("%s; lines: %s; all 8 shared TFs ordered best to worst (mean matched row-z), rows ordered the same way in each line; screen-only KOs MEF2A / MEF2D at the right", ifelse(mode == "best", "best replicate per KO", "replicates averaged"), lines_txt),
    sprintf("all-TF score z = %.2f (exact permutation p = %.3f) | %s", ev$z_matched, ev$perm_p, sub_txt), if (part) "black dot after a column name = TF in the subset that gave the best score in the search" else NULL, if (part) subset_tfs else NULL)
}
draw_best("fig16a_best_combination_any_subset", d1, b1$cfg, b1$mode, b1$lines, sprintf("best 4-TF subset z = %.2f, search-wide p = %.2f", d1$obs$z_matched, fw["all"]), d1$cols)
draw_best("fig16b_best_combination_all_TFs_and_lines", d2, fcfg, fb$mode, "all 4", sprintf("search-wide p (all-TF/all-line configs) = %.2f", fw["full_sets"]), d2$cols)
