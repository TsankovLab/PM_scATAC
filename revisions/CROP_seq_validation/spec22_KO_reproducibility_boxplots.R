# Spec 22: standalone figure with ONLY the boxplot panels from fig24 (cross-line correlations of each KO with itself vs with the other KOs), no heatmap.
# Two panels: A replicates averaged, B best replicate per KO/line (selection-biased). Same computation as spec20_best_config_heatmap.R; default configuration =
# the one used for fig24_best_config_heatmap_modules_screenKOs_t_cc (20 Cm modules, cc-adjusted t statistic, CROP-seq genes/KOs only). Significance shown as
# asterisks from the one-sided Wilcoxon p (pairs share samples, so this is optimistic, especially for panel B); KOs ordered most- to least-significant by that p, per panel.
# Set SPEC22_SPACE / SPEC22_TYPE to show another configuration (space: modules, top500, top1000, top3000, all_genes; type: t, t_cc, diff, diff_cc).
# Set SPEC22_COR_METHOD=spearman for Spearman rho instead of Pearson r (comparison run only, does not affect the default fig26).
# Set SPEC22_ZSCORE=1 to z-score each feature (module/gene) across samples within each line before the stat step (comparison run only, does not affect the default fig26).
# needs: spec1_prepare_stats.R
# outputs: spec22_KO_reproducibility.csv, fig26_KO_reproducibility_boxplots.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_variants.R"))
SPACE <- Sys.getenv("SPEC22_SPACE", "modules"); TYPE <- Sys.getenv("SPEC22_TYPE", "t_cc")
EXCL_OUTLIERS <- Sys.getenv("SPEC22_EXCLUDE_VISUAL_OUTLIERS", "0") == "1"     # drop VISUAL_OUTLIERS (repro_common.R) -- comparison run only, does not affect the default fig26
METHOD <- Sys.getenv("SPEC22_COR_METHOD", "pearson")                        # pearson (default) or spearman -- comparison run only, does not affect the default fig26
ZSCORE <- Sys.getenv("SPEC22_ZSCORE", "0") == "1"                           # z-score features per line before the stat step -- comparison run only, does not affect the default fig26
ctx <- load_spec_context(LINES); bs <- ctx$bs; lg <- ctx$lg; meta <- ctx$meta; ccb <- ctx$ccb; genes <- Reduce(intersect, lapply(bs[LINES], `[[`, "genes"))
if (EXCL_OUTLIERS) { dropped <- rownames(meta)[meta$Sample_ID %in% VISUAL_OUTLIERS]; cat("excluding VISUAL_OUTLIERS present in this sample set:", paste(meta[dropped, "Sample_ID"], collapse = ", "), "\n"); meta <- meta[!meta$Sample_ID %in% VISUAL_OUTLIERS, ] }
GENESET <- Sys.getenv("CROP_GENESET", "screen"); KOSET <- Sys.getenv("CROP_KOSET", "screen")
SUF <- paste0(if (SPACE == "modules") "_modules" else if (GENESET == "screen") "_screengenes" else paste0("_", SPACE, "_allgenes"), if (KOSET == "screen") "_screenKOs" else "_allKOs", if (TYPE != "t") paste0("_", TYPE) else "", if (EXCL_OUTLIERS) "_noOutliers" else "", if (METHOD == "spearman") "_spearman" else "", if (ZSCORE) "_zscore" else "")
if (GENESET == "screen") genes <- intersect(genes, ctx$ss$genes); cat("genes used:", length(genes), "(gene set:", GENESET, ")\n")
sc_df <- read.csv(file.path(PP, "cm_sample_scores_AddModuleScore_bulk_RUN1.csv"), check.names = FALSE); cm_names <- setdiff(colnames(sc_df), c("sample", "cell_line"))
lv <- sapply(LINES, function(cl) matrixStats::rowVars(lg[genes, rownames(meta)[meta$Cell_Line == cl]])); rownames(lv) <- genes; vord <- names(sort(rowMeans(lv), decreasing = TRUE))
feats <- switch(SPACE, modules = NULL, top500 = vord[1:500], top1000 = vord[1:1000], top3000 = vord[1:3000], top5000 = vord[1:min(5000, length(vord))], all_genes = genes)
make_types <- function(y, cc, grp) {
  ko <- grp != "NTC"; ntc <- grp == "NTC"; nn <- sum(ntc); ccc <- cc - mean(cc); reg <- function(M) { b <- as.numeric((M - rowMeans(M)) %*% ccc) / sum(ccc^2); M - outer(b, ccc) }
  D <- function(M) M[, ko, drop = FALSE] - rowMeans(M[, ntc, drop = FALSE])
  tstat <- function(M, with_cc) { g <- factor(grp, levels = c("NTC", sort(setdiff(unique(grp), "NTC")))); X <- if (with_cc) model.matrix(~ cc + g) else model.matrix(~ g); q <- qr(X); res <- qr.resid(q, t(M))
    s2 <- pmax(colSums(res^2) / (ncol(M) - q$rank), 1e-4); Ma <- if (with_cc) M - outer(qr.coef(q, t(M))["cc", ], ccc) else M; (Ma[, ko, drop = FALSE] - rowMeans(Ma[, ntc, drop = FALSE])) / sqrt(s2 * (1 + 1 / nn)) }
  out <- list(t = tstat(y, FALSE), t_cc = tstat(y, TRUE), diff = D(y), diff_cc = D(reg(y))); lapply(out, function(m) { m[!is.finite(m)] <- 0; m }) }
P <- list(); for (cl in LINES) { m <- meta[meta$Cell_Line == cl, ]; s <- rownames(m); cc <- ccb$cc[match(s, ccb$sample)]; grp <- as.character(m$Gene_Control)
  y <- if (SPACE == "modules") t(as.matrix(sc_df[match(s, sc_df$sample), cm_names])) |> `colnames<-`(s) else lg[feats, s, drop = FALSE]
  if (ZSCORE) { y <- (y - rowMeans(y)) / matrixStats::rowSds(y); y[!is.finite(y)] <- 0 }
  M <- make_types(y, cc, grp)[[TYPE]]; ko <- grp[grp != "NTC"]; if (KOSET == "screen") { keepr <- ko %in% ORDER; M <- M[, keepr, drop = FALSE]; ko <- ko[keepr] }
  kos <- sort(unique(ko), method = "radix"); G <- sapply(setNames(kos, kos), function(k) rowMeans(M[, ko == k, drop = FALSE])); mu <- rowMeans(G); R <- M - mu; G <- G - mu
  if (grepl("diff", TYPE)) { sc <- sqrt(rowMeans(G^2)); sc[sc < 1e-12] <- NA; R <- R / sc; G <- G / sc; R[!is.finite(R)] <- 0; G[!is.finite(G)] <- 0 }; P[[cl]] <- list(R = R, G = G, ko = ko, kos = kos) }
kos_all <- sort(unique(unlist(lapply(P, `[[`, "kos"))), method = "radix")
units <- do.call(rbind, lapply(kos_all, function(k) do.call(rbind, lapply(LINES, function(l) if (k %in% P[[l]]$kos) data.frame(KO = k, line = l, id = paste(k, l), stringsAsFactors = FALSE))))); n <- nrow(units)
avgP <- sapply(seq_len(n), function(i) P[[units$line[i]]]$G[, units$KO[i]]); colnames(avgP) <- units$id
pickrep <- vapply(seq_len(n), function(i) { l <- units$line[i]; k <- units$KO[i]; reps <- colnames(P[[l]]$R)[P[[l]]$ko == k]; others <- which(units$KO == k & units$line != l)
  if (!length(others) || length(reps) == 1) return(reps[1]); sc <- sapply(reps, function(r) mean(cor(P[[l]]$R[, r], avgP[, others, drop = FALSE], method = METHOD))); unname(reps[which.max(sc)]) }, character(1))
NPERM <- as.integer(Sys.getenv("SPEC22_NPERM", "5000")); set.seed(22)
Rall <- do.call(cbind, lapply(LINES, function(l) P[[l]]$R)); rep_group <- unlist(lapply(LINES, function(l) paste(P[[l]]$ko, l)))
CorRR <- cor(Rall, method = METHOD); CorRG <- cor(Rall, avgP, method = METHOD); CG <- cor(avgP, method = METHOD); same_line <- outer(units$line, units$line, "==")
lines_of <- setNames(lapply(LINES, function(l) which(units$line == l)), LINES)
mats_for <- function(lab) {
  CAl <- CG; CAl[same_line] <- NA; chosen <- integer(n)
  for (u in seq_len(n)) { reps <- which(rep_group == units$id[u]); k <- lab[u]; others <- which(lab == k & units$line != units$line[u])
    chosen[u] <- if (!length(others) || length(reps) == 1) reps[1] else reps[which.max(rowMeans(CorRG[reps, others, drop = FALSE]))] }
  CBl <- CorRR[chosen, chosen]; CBl[same_line] <- NA; list(A = CAl, B = CBl, chosen = chosen) }
ko_stat <- function(C, lab, k) { i <- which(lab == k); o <- which(lab != k); sm <- C[i, i][upper.tri(C[i, i])]; sm <- sm[!is.na(sm)]; ot <- as.vector(C[i, o]); ot <- ot[!is.na(ot)]
  list(same = sm, other = ot, diff = if (length(sm)) median(sm) - median(ot) else NA_real_) }
lab0 <- units$KO; M0 <- mats_for(lab0); kos_stat <- kos_all[kos_all %in% units$KO]
tab <- do.call(rbind, lapply(c("A", "B"), function(m) do.call(rbind, lapply(kos_stat, function(k) { st <- ko_stat(M0[[m]], lab0, k)
  data.frame(version = ifelse(m == "A", "averaged", "best replicate"), KO = k, n_same_pairs = length(st$same), median_same = if (length(st$same)) median(st$same) else NA, median_other = median(st$other),
             p_wilcoxon = if (length(st$same) >= 2) wilcox.test(st$same, st$other, alternative = "greater", exact = FALSE)$p.value else NA, diff_obs = st$diff) }))))
DIST <- lapply(c(A = "A", B = "B"), function(m) setNames(lapply(kos_stat, function(k) ko_stat(M0[[m]], lab0, k)), kos_stat))
perm_diff <- replicate(NPERM, { lab <- lab0; for (l in LINES) { ix <- lines_of[[l]]; lab[ix] <- lab0[ix][sample(length(ix))] }; M <- mats_for(lab)
  c(sapply(kos_stat, function(k) ko_stat(M$A, lab, k)$diff), sapply(kos_stat, function(k) ko_stat(M$B, lab, k)$diff)) })
nA <- length(kos_stat); tab$p_perm <- NA_real_
for (i in seq_len(nrow(tab))) { j <- if (tab$version[i] == "averaged") match(tab$KO[i], kos_stat) else nA + match(tab$KO[i], kos_stat); tab$p_perm[i] <- if (is.na(tab$diff_obs[i])) NA else (1 + sum(perm_diff[j, ] >= tab$diff_obs[i], na.rm = TRUE)) / (NPERM + 1) }
tab$p_perm_BH <- ave(tab$p_perm, tab$version, FUN = function(p) p.adjust(p, "BH")); write.csv(tab, file.path(OUT, paste0("spec22_KO_reproducibility", SUF, ".csv")), row.names = FALSE)
options(width = 200); cat(sprintf("configuration: %s | %s | %d genes/features | %d permutations\n", SPACE, TYPE, nrow(avgP), NPERM)); print(format(tab[, c("version", "KO", "n_same_pairs", "median_same", "median_other", "p_wilcoxon", "p_perm", "p_perm_BH")], digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- figure: boxplots only, ordered most- to least-significant (permutation p) per panel
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
MOD <- SPACE == "modules"; RLO <- if (MOD) -1 else -0.5; RHI <- if (MOD) 1 else 0.6
save_fig(paste0("fig26_KO_reproducibility_boxplots", SUF), 9.5, 5.4, function() {
  layout(matrix(1:2, nrow = 1)); par(oma = c(3.4, 0, 5.4, 0.3), family = "sans")
  panel <- function(ver, m, ttl) {
    tb <- tab[tab$version == ver, ]; o <- tb[order(ifelse(is.na(tb$p_wilcoxon), Inf, tb$p_wilcoxon)), ]; k <- nrow(o)
    par(mar = c(2.4, 6.6, 3.4, 2.6)); plot(NA, xlim = c(RLO, RHI), ylim = c(k + 0.6, 0.4), axes = FALSE, xlab = "", ylab = "")
    axis(1, at = pretty(c(RLO, RHI)), cex.axis = 0.75, lwd = 0.6, mgp = c(1.4, 0.35, 0)); abline(v = 0, col = GRID, lwd = 0.8)
    bx <- function(v, yc, hh, col) { if (!length(v)) return(invisible()); q <- quantile(v, c(0.25, 0.5, 0.75)); iq <- q[3] - q[1]; lo <- min(v[v >= q[1] - 1.5 * iq]); hi <- max(v[v <= q[3] + 1.5 * iq])
      segments(lo, yc, hi, yc, col = INK2, lwd = 0.9); rect(q[1], yc - hh, q[3], yc + hh, col = adjustcolor(col, 0.9), border = INK2, lwd = 0.9); segments(q[2], yc - hh, q[2], yc + hh, col = INK, lwd = 2)
      points(v, yc + runif(length(v), -hh * 0.7, hh * 0.7), pch = 16, cex = if (length(v) > 20) 0.35 else 0.6, col = adjustcolor(INK, 0.45)) }
    for (i in seq_len(k)) { yc <- i; st <- DIST[[m]][[o$KO[i]]]
      bx(st$other, yc + 0.24, 0.17, GREY); bx(st$same, yc - 0.24, 0.17, KO_COL[o$KO[i]])
      text(RHI + (RHI - RLO) * 0.04, yc, stars(o$p_wilcoxon[i]), adj = 0, cex = 0.85, xpd = NA, font = 2) }
    for (i in seq_len(k)) text(RLO - (RHI - RLO) * 0.03, i, o$KO[i], adj = 1, xpd = NA, cex = 0.78, font = 2, col = KO_COL[o$KO[i]])
    mtext(ttl, side = 3, adj = 0, line = 2.3, font = 2, cex = 0.85)
    mtext(sprintf("ordered by Wilcoxon p (most to least significant); n=%d line pairs", o$n_same_pairs[1]), side = 3, adj = 0, line = 1.1, cex = 0.55, col = INK2)
  }
  panel("averaged", "A", "A  Replicates averaged"); panel("best replicate", "B", "B  Best replicate per KO and line (selection-biased)")
  legend("bottomright", inset = c(0.02, -0.34), xpd = NA, horiz = TRUE, legend = c("same KO, other line", "other KOs"), fill = c(adjustcolor(INK, 0.9), GREY), border = INK2, bty = "n", cex = 0.65)
  mtext("How well does each bulk KO reproduce across cell lines?", outer = TRUE, side = 3, adj = 0, line = 2.9, font = 2, cex = 1)
  mtext(sprintf("%s features, %s statistic, centred on each line%s%s", if (MOD) "20 Cm module" else ifelse(SPACE == "all_genes", ifelse(GENESET == "screen", "CROP-seq-panel", "all"), SPACE), TYPE, if (KOSET == "screen") "; the 8 CROP-seq KOs" else "", if (ZSCORE) "; z-scored per feature" else ""), outer = TRUE, side = 3, adj = 0, line = 1.9, cex = 0.6, col = INK2)
  corlab <- if (METHOD == "spearman") "Spearman rho" else "Pearson r"
  mtext(sprintf("boxplots: cross-line %s of each KO with itself (colour) vs with the other KOs (grey), median/IQR/range, dots = individual line pairs", corlab), outer = TRUE, side = 3, adj = 0, line = 1.0, cex = 0.55, col = INK2)
  mtext("*** p<0.001, ** p<0.01, * p<0.05, ns = not significant (one-sided Wilcoxon test on the pooled cross-line correlations; pairs share samples so this is optimistic, esp. for B)", outer = TRUE, side = 3, adj = 0, line = 0.1, cex = 0.55, col = INK2)
  mtext(sprintf("cross-line %s", corlab), side = 1, outer = TRUE, line = 2.1, cex = 0.7, col = INK2)
})
