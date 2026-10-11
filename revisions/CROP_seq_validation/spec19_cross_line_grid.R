# Spec 19: cross-cell-line reproducibility of each bulk KO over a grid of choices.
#   feature space: 20 Cm modules | top 500 / 1000 / 3000 variable genes (mean within-line variance) | all genes expressed in the 4 lines
#   profile: pooled-variance t | cc-adjusted t | KO - NTC difference (RMS-scaled) | cc-regressed difference          (centred on the line's mean KO profile in every case)
#   replicates: averaged | ONLY the replicate best matching the same KO in the other line | control (best replicate picked separately for every cell)
# For an ordered line pair a -> b and KO k: row = profile of KO k in a (average / chosen replicate), correlated (Pearson) with all KO profiles of b; score z = matched (same-KO) r standardised across
# the KOs of b. KO score = mean z over the directions/line pairs where k exists in both lines. Significance: permutation of KO labels within every line (same shuffles for every configuration; the replicate
# choice is redone under each shuffle), statistic = mean over KOs of the KO score; family-wise p = P(max over all configurations/modes >= observed max).
# needs: spec1_prepare_stats.R
# outputs: spec19_config_summary.csv, spec19_KO_scores.csv, fig23_cross_line_grid.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_variants.R"))
suppressPackageStartupMessages(library(parallel)); set.seed(19); NPERM <- as.integer(Sys.getenv("SPEC19_NPERM", "2000")); NCORE <- as.integer(Sys.getenv("SPEC19_NCORE", "6")); options(warn = as.integer(Sys.getenv("SPEC19_WARN", "0")))
ctx <- load_spec_context(LINES); bs <- ctx$bs; lg <- ctx$lg; meta <- ctx$meta; ccb <- ctx$ccb
genes <- Reduce(intersect, lapply(bs[LINES], `[[`, "genes"))
GENESET <- Sys.getenv("CROP_GENESET", "screen"); KOSET <- Sys.getenv("CROP_KOSET", "screen"); SUF <- paste0(if (GENESET == "screen") "_screengenes" else "", if (KOSET == "screen") "_screenKOs" else "")    # KOSET "screen": only the 8 KOs that are also in the CROP-seq screen    # "screen": only genes measured in the CROP-seq screen (detected in >= 1% of cells)
if (GENESET == "screen") genes <- intersect(genes, ctx$ss$genes); cat("genes used:", length(genes), "(gene set:", GENESET, "| KO set:", KOSET, ")\n")
sc_df <- read.csv(file.path(PP, "cm_sample_scores_AddModuleScore_bulk_RUN1.csv"), check.names = FALSE); cm_names <- setdiff(colnames(sc_df), c("sample", "cell_line"))
lv <- sapply(LINES, function(cl) matrixStats::rowVars(lg[genes, rownames(meta)[meta$Cell_Line == cl]])); rownames(lv) <- genes; vord <- names(sort(rowMeans(lv), decreasing = TRUE))
SP <- list(modules = NULL, top500 = vord[1:500], top1000 = vord[1:1000], top3000 = vord[1:3000], all_genes = genes)
make_types <- function(y, cc, grp) {
  ko <- grp != "NTC"; ntc <- grp == "NTC"; nn <- sum(ntc); ccc <- cc - mean(cc)
  reg <- function(M) { b <- as.numeric((M - rowMeans(M)) %*% ccc) / sum(ccc^2); M - outer(b, ccc) }
  D <- function(M) M[, ko, drop = FALSE] - rowMeans(M[, ntc, drop = FALSE])
  tstat <- function(M, with_cc) { g <- factor(grp, levels = c("NTC", sort(setdiff(unique(grp), "NTC")))); X <- if (with_cc) model.matrix(~ cc + g) else model.matrix(~ g)
    q <- qr(X); res <- qr.resid(q, t(M)); s2 <- pmax(colSums(res^2) / (ncol(M) - q$rank), 1e-4); Ma <- if (with_cc) M - outer(qr.coef(q, t(M))["cc", ], ccc) else M
    (Ma[, ko, drop = FALSE] - rowMeans(Ma[, ntc, drop = FALSE])) / sqrt(s2 * (1 + 1 / nn)) }
  out <- list(t = tstat(y, FALSE), t_cc = tstat(y, TRUE), diff = D(y), diff_cc = D(reg(y))); lapply(out, function(m) { m[!is.finite(m)] <- 0; m }) }
TYPES <- c("t", "t_cc", "diff", "diff_cc"); TLAB <- c(t = "t", t_cc = "t + cc", diff = "diff", diff_cc = "diff + cc")

## ---------------------------------------------------------------- per line: replicate profiles and group profiles for every (space, type)
PR <- list()
for (cl in LINES) { m <- meta[meta$Cell_Line == cl, ]; s <- rownames(m); cc <- ccb$cc[match(s, ccb$sample)]; grp <- as.character(m$Gene_Control)
  for (sp in names(SP)) { y <- if (sp == "modules") t(as.matrix(sc_df[match(s, sc_df$sample), cm_names])) |> `colnames<-`(s) else lg[SP[[sp]], s, drop = FALSE]; tt <- make_types(y, cc, grp); ko <- grp[grp != "NTC"]
    for (ty in TYPES) { M <- tt[[ty]]; ko_use <- ko; if (KOSET == "screen") { keepr <- ko %in% ORDER; M <- M[, keepr, drop = FALSE]; ko_use <- ko[keepr] }
      kos <- sort(unique(ko_use), method = "radix"); G <- sapply(setNames(kos, kos), function(k) rowMeans(M[, ko_use == k, drop = FALSE])); mu <- rowMeans(G); Rr <- M - mu; G <- G - mu
      if (grepl("diff", ty)) { sc <- sqrt(rowMeans(G^2)); sc[sc < 1e-12] <- NA; Rr <- Rr / sc; G <- G / sc; Rr[!is.finite(Rr)] <- 0; G[!is.finite(G)] <- 0 }
      PR[[paste(cl, sp, ty)]] <- list(R = Rr, G = G, ko = ko_use, kos = kos) } } }
cfgs <- expand.grid(ty = TYPES, sp = names(SP), stringsAsFactors = FALSE)[, 2:1]; cfgs$id <- paste(cfgs$sp, cfgs$ty); cat("configurations:", nrow(cfgs), "\n")

## ---------------------------------------------------------------- Z matrices per configuration and direction (independent of the permutation)
dirs <- do.call(rbind, lapply(LINES, function(a) data.frame(a = a, b = setdiff(LINES, a)))); zrow <- function(r, j) (r[j] - mean(r)) / sd(r)
build_cfg <- function(sp, ty) lapply(seq_len(nrow(dirs)), function(d) { a <- dirs$a[d]; b <- dirs$b[d]; A <- PR[[paste(a, sp, ty)]]; B <- PR[[paste(b, sp, ty)]]
  Cavg <- cor(A$G, B$G); Crep <- cor(A$R, B$G); nA <- length(A$kos); nB <- length(B$kos); Za <- Zb <- Zc <- matrix(NA, nA, nB)
  for (ga in seq_len(nA)) { rep_i <- which(A$ko == A$kos[ga]); for (cb in seq_len(nB)) {
      Za[ga, cb] <- zrow(Cavg[ga, ], cb); rr <- Crep[rep_i, , drop = FALSE]; Zb[ga, cb] <- zrow(rr[which.max(rr[, cb]), ], cb); Zc[ga, cb] <- zrow(apply(rr, 2, max), cb) } }
  list(avg = Za, best = Zb, ctl = Zc, ia = match(intersect(A$kos, B$kos), A$kos), ib = match(intersect(A$kos, B$kos), B$kos), shared = intersect(A$kos, B$kos), nA = nA, nB = nB, a = a, b = b) })
CF <- mclapply(seq_len(nrow(cfgs)), function(i) build_cfg(cfgs$sp[i], cfgs$ty[i]), mc.cores = NCORE)
kos_all <- sort(unique(unlist(lapply(LINES, function(l) PR[[paste(l, "all_genes t")]]$kos))), method = "radix"); modes <- c("avg", "best", "ctl")
ko_scores <- function(cf, mode, pa = NULL) {                                    # named vector: mean z per KO; pa = list of permutations per line (or NULL = identity)
  tot <- setNames(numeric(length(kos_all)), kos_all); cnt <- tot
  for (d in seq_along(cf)) { x <- cf[[d]]; pia <- if (is.null(pa)) seq_len(x$nA) else pa[[x$a]]; pib <- if (is.null(pa)) seq_len(x$nB) else pa[[x$b]]
    v <- x[[mode]][cbind(pia[x$ia], pib[x$ib])]; tot[x$shared] <- tot[x$shared] + v; cnt[x$shared] <- cnt[x$shared] + 1 }
  ifelse(cnt > 0, tot / pmax(cnt, 1), NA) }
nk <- sapply(LINES, function(l) length(PR[[paste(l, "all_genes t")]]$kos))
obs <- lapply(seq_len(nrow(cfgs)), function(i) sapply(modes, function(m) ko_scores(CF[[i]], m))); glob_obs <- t(sapply(obs, function(M) colMeans(M, na.rm = TRUE)))
perms <- lapply(seq_len(NPERM), function(p) setNames(lapply(LINES, function(l) sample(nk[[l]])), LINES))
t0 <- Sys.time(); nullres <- mclapply(perms, function(pa) { g <- t(sapply(seq_len(nrow(cfgs)), function(i) sapply(modes, function(m) mean(ko_scores(CF[[i]], m, pa), na.rm = TRUE)))); g }, mc.cores = NCORE)
NG <- array(unlist(nullres), dim = c(nrow(cfgs), length(modes), NPERM)); cat("permutations done:", NPERM, "in", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
pv <- sapply(seq_along(modes), function(m) sapply(seq_len(nrow(cfgs)), function(i) (1 + sum(NG[i, m, ] >= glob_obs[i, m])) / (NPERM + 1)))
maxnull <- apply(NG, 3, function(g) c(all = max(g), avg = max(g[, 1]), best = max(g[, 2]))); maxobs <- c(all = max(glob_obs), avg = max(glob_obs[, 1]), best = max(glob_obs[, 2]))
fw <- sapply(names(maxobs), function(n) (1 + sum(maxnull[n, ] >= maxobs[n])) / (NPERM + 1))
summ <- do.call(rbind, lapply(seq_along(modes), function(m) data.frame(cfgs, mode = modes[m], mean_z = glob_obs[, m], null_mean = apply(NG[, m, , drop = FALSE], 1, mean), null_sd = apply(NG[, m, , drop = FALSE], 1, sd), perm_p = pv[, m]))); summ <- summ[order(-summ$mean_z), ]; write.csv(summ, file.path(OUT, paste0("spec19_config_summary", SUF, ".csv")), row.names = FALSE)
kosc <- do.call(rbind, lapply(seq_len(nrow(cfgs)), function(i) do.call(rbind, lapply(seq_along(modes), function(m) data.frame(cfgs[i, ], mode = modes[m], KO = kos_all, z = unname(obs[[i]][, m])))))); write.csv(kosc, file.path(OUT, paste0("spec19_KO_scores", SUF, ".csv")), row.names = FALSE)
options(width = 200); cat("\n=== 15 best configurations (mean KO z; permutation p for that configuration only) ===\n"); print(format(head(summ, 15), digits = 3), row.names = FALSE)
cat("\n=== search-wide: best observed vs the maximum the same grid finds on shuffled labels ===\n"); print(round(rbind(observed_max = maxobs, null_median = apply(maxnull, 1, median), null_95th = apply(maxnull, 1, quantile, 0.95), family_wise_p = fw), 3))
cat("\n=== main effects: mean of the KO-mean z over configurations ===\n"); for (d in c("sp", "ty", "mode")) { a <- tapply(summ$mean_z, summ[[d]], mean); cat(sprintf("%-5s", d), paste(sprintf("%s %.3f", names(a), a), collapse = " | "), "\n") }
cat("\nmean over configurations: observed z vs mean z under shuffled labels (selection alone):\n"); a <- aggregate(cbind(mean_z, null_mean) ~ mode, summ, mean); a[, -1] <- round(a[, -1], 3); print(a, row.names = FALSE)
cat("\nmean z by cc adjustment (t vs t+cc, diff vs diff+cc):\n"); print(round(tapply(summ$mean_z, list(summ$ty, summ$mode), mean), 3))
bestcfg <- summ[summ$mode == "avg", ][1, ]; cat(sprintf("\nbest replicate-averaged configuration: %s | %s (mean z %.2f, p = %.3f)\n", bestcfg$sp, bestcfg$ty, bestcfg$mean_z, bestcfg$perm_p))
avgm <- do.call(cbind, lapply(cfgs$id, function(id) kosc$z[kosc$id == id & kosc$mode == "avg"])); rownames(avgm) <- kos_all; colnames(avgm) <- cfgs$id
bestm <- do.call(cbind, lapply(cfgs$id, function(id) kosc$z[kosc$id == id & kosc$mode == "best"])); rownames(bestm) <- kos_all; colnames(bestm) <- cfgs$id; ctlm <- do.call(cbind, lapply(cfgs$id, function(id) kosc$z[kosc$id == id & kosc$mode == "ctl"])); rownames(ctlm) <- kos_all
ord <- names(sort(rowMeans(avgm, na.rm = TRUE), decreasing = TRUE)); cat("\nper-KO mean z over all configurations (averaged replicates | best replicate | control):\n"); print(round(cbind(avg = rowMeans(avgm, na.rm = TRUE), best = rowMeans(bestm, na.rm = TRUE), control = rowMeans(ctlm, na.rm = TRUE))[ord, ], 2))

## ---------------------------------------------------------------- figure
save_fig(paste0("fig23_cross_line_grid", SUF), 16, 10.5, function() {
  layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE), heights = c(1, 1)); par(oma = c(0, 0, 4, 0.5), family = "sans")
  hm <- function(M, ttl) { par(mar = c(6.5, 5, 3.6, 0.5)); nn <- ncol(M); image(seq_len(nn), seq_along(ord), t(pmax(pmin(M[ord, ], 2), -2))[, length(ord):1], col = DIV, zlim = c(-2, 2), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_along(ord)) for (j in seq_len(nn)) text(j, length(ord) - i + 1, ifelse(is.na(M[ord[i], j]), "", sprintf("%.1f", M[ord[i], j])), cex = 0.45); axis(2, at = length(ord):1, labels = ord, las = 1, tick = FALSE, cex.axis = 0.7)
    axis(1, at = seq_len(nn), labels = sub(" ", "\n", cfgs$id), las = 2, tick = FALSE, cex.axis = 0.45, line = -0.4); sp_edges <- seq(4, nn - 1, 4) + 0.5; abline(v = sp_edges, col = "white", lwd = 3); mtext(ttl, side = 3, adj = 0, line = 0.6, font = 2, cex = 0.85) }
  hm(avgm, "A  KO score (mean z) per configuration: replicates averaged"); hm(bestm, "B  best replicate (selection-biased)")
  par(mar = c(3.6, 8, 3.6, 1)); dims <- c("sp", "ty", "mode"); lev <- lapply(dims, function(d) sort(unique(as.character(summ[[d]])))); nn <- sum(lengths(lev)) + length(dims) - 1; plot(NA, xlim = c(-0.3, 0.9), ylim = c(nn + 0.5, 0.5), axes = FALSE, xlab = "", ylab = ""); abline(v = 0, col = INK2); y <- 0
  for (di in seq_along(dims)) { for (l in lev[[di]]) { y <- y + 1; v <- summ$mean_z[as.character(summ[[dims[di]]]) == l]; points(v, y + runif(length(v), -0.25, 0.25), pch = 16, cex = 0.4, col = adjustcolor(INK2, 0.4)); points(mean(v), y, pch = 23, bg = c("#2a78d6", "#eb6834", "#4a3aa7")[di], col = INK, cex = 1.4)
      text(-0.32, y, sprintf("%s: %s", c("features", "profile", "replicates")[di], ifelse(dims[di] == "ty", TLAB[l], l)), adj = 1, xpd = NA, cex = 0.65) }; y <- y + 1; abline(h = y - 0.5, col = GRID) }
  axis(1, cex.axis = 0.7, lwd = 0.5); mtext("mean KO score z", side = 1, line = 2.2, cex = 0.7); mtext("C  Main effects over configurations", side = 3, adj = 0, line = 0.6, font = 2, cex = 0.85)
  par(mar = c(3.6, 4.4, 3.6, 1)); rg <- range(c(maxnull["all", ], maxobs["all"])); brk <- seq(floor(rg[1] * 20) / 20, ceiling(rg[2] * 20) / 20 + 0.05, by = 0.025); h <- hist(maxnull["all", ], breaks = brk, plot = FALSE)
  plot(NA, xlim = range(brk), ylim = c(0, max(h$counts) * 1.2), xlab = "", ylab = "shuffles", bty = "n", cex.axis = 0.7, cex.lab = 0.75); rect(head(brk, -1), 0, brk[-1], h$counts, col = adjustcolor(GREY, 0.9), border = "white")
  abline(v = maxobs["all"], lwd = 2, col = "#b3400f"); abline(v = maxobs["avg"], lwd = 2, col = "#2a78d6", lty = 2); mtext("best mean z found by the grid", side = 1, line = 2.2, cex = 0.7)
  legend("topright", c(sprintf("observed best, any mode %.2f (family-wise p = %.2f)", maxobs["all"], fw["all"]), sprintf("observed best, replicates averaged %.2f (p = %.2f)", maxobs["avg"], fw["avg"]), "grey: same grid on shuffled KO labels"), lty = c(1, 2, NA), lwd = 2, col = c("#b3400f", "#2a78d6", NA), fill = c(NA, NA, GREY), border = NA, bty = "n", cex = 0.6, x.intersp = 0.5)
  mtext("D  Is the best configuration better than chance?", side = 3, adj = 0, line = 0.6, font = 2, cex = 0.85)
  mtext("Cross-cell-line reproducibility of each bulk KO across feature spaces, cc adjustment and replicate handling", outer = TRUE, side = 3, adj = 0, line = 2.4, font = 2, cex = 1.05)
  mtext(paste0(if (GENESET == "screen") sprintf("only the %d genes measured in the CROP-seq screen; ", length(genes)) else "", if (KOSET == "screen") "only the 8 KOs shared with the CROP-seq screen; " else "", sprintf("%d configurations x 3 replicate modes; %d label shuffles within lines", nrow(cfgs), NPERM)), outer = TRUE, side = 3, adj = 0, line = 1.3, cex = 0.58, col = INK2)
  mtext("KO score z = same-KO correlation between two lines standardised against the other KOs (0 = no better than other KOs); columns of A/B grouped by feature space (modules, top500, top1000, top3000, all genes), then profile type", outer = TRUE, side = 3, adj = 0, line = 0.4, cex = 0.55, col = INK2)
})
