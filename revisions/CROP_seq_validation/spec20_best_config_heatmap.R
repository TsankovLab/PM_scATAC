# Spec 20: cross-cell-line correlation heatmap (blue-white-red diverging palette; side table: per KO median same-KO vs other-KO cross-line r with one-sided Wilcoxon and label-permutation p)
# (original description:) cross-cell-line correlation heatmap of the bulk KO profiles for the best configuration of spec19 (default: all expressed genes, pooled-variance t statistic, centred on each line's mean KO profile).
# A: replicate-averaged KO profiles; B: one replicate per KO/line (the replicate best correlated with the same KO in the other lines; selection-biased). Rows/columns = (KO, line), ordered by KO then line;
# same-line cells are greyed (shared NTC / background); boxed blocks = the same KO in different lines. C: per KO, mean same-KO r (across lines) vs mean r with other KOs.
# Set SPEC20_SPACE / SPEC20_TYPE to show another configuration (space: modules, top500, top1000, top3000, all_genes; type: t, t_cc, diff, diff_cc).
# Set SPEC20_COR_METHOD=spearman / SPEC20_ZSCORE=1 / SPEC20_EXCLUDE_VISUAL_OUTLIERS=1 for the same comparison knobs as spec22 (does not affect the default fig24).
# needs: spec1_prepare_stats.R
# outputs: spec20_block_summary.csv (medians and p-values), fig24_best_config_heatmap.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_variants.R"))
SPACE <- Sys.getenv("SPEC20_SPACE", "all_genes"); TYPE <- Sys.getenv("SPEC20_TYPE", "t")
METHOD <- Sys.getenv("SPEC20_COR_METHOD", "pearson"); ZSCORE <- Sys.getenv("SPEC20_ZSCORE", "0") == "1"; EXCL_OUTLIERS <- Sys.getenv("SPEC20_EXCLUDE_VISUAL_OUTLIERS", "0") == "1"
ctx <- load_spec_context(LINES); bs <- ctx$bs; lg <- ctx$lg; meta <- ctx$meta; ccb <- ctx$ccb; genes <- Reduce(intersect, lapply(bs[LINES], `[[`, "genes"))
if (EXCL_OUTLIERS) { dropped <- rownames(meta)[meta$Sample_ID %in% VISUAL_OUTLIERS]; cat("excluding VISUAL_OUTLIERS present in this sample set:", paste(meta[dropped, "Sample_ID"], collapse = ", "), "\n"); meta <- meta[!meta$Sample_ID %in% VISUAL_OUTLIERS, ] }
GENESET <- Sys.getenv("CROP_GENESET", "screen"); KOSET <- Sys.getenv("CROP_KOSET", "screen"); SUF <- paste0(if (SPACE == "modules") "_modules" else if (GENESET == "screen") "_screengenes" else paste0("_", SPACE, "_allgenes"), if (KOSET == "screen") "_screenKOs" else "_allKOs", if (TYPE != "t") paste0("_", TYPE) else "", if (EXCL_OUTLIERS) "_noOutliers" else "", if (METHOD == "spearman") "_spearman" else "", if (ZSCORE) "_zscore" else "")    # "screen": only genes measured in the CROP-seq screen (detected in >= 1% of cells)
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
# units and their average / chosen-replicate profiles
units <- do.call(rbind, lapply(kos_all, function(k) do.call(rbind, lapply(LINES, function(l) if (k %in% P[[l]]$kos) data.frame(KO = k, line = l, id = paste(k, l), stringsAsFactors = FALSE))))); n <- nrow(units)
avgP <- sapply(seq_len(n), function(i) P[[units$line[i]]]$G[, units$KO[i]]); colnames(avgP) <- units$id
pickrep <- vapply(seq_len(n), function(i) { l <- units$line[i]; k <- units$KO[i]; reps <- colnames(P[[l]]$R)[P[[l]]$ko == k]; others <- which(units$KO == k & units$line != l)
  if (!length(others) || length(reps) == 1) return(reps[1]); sc <- sapply(reps, function(r) mean(cor(P[[l]]$R[, r], avgP[, others, drop = FALSE], method = METHOD))); unname(reps[which.max(sc)]) }, character(1))
bestP <- sapply(seq_len(n), function(i) P[[units$line[i]]]$R[, pickrep[i]]); colnames(bestP) <- units$id
PAL <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100); DIV <- PAL       # blue (negative) - white - red (positive) diverging
NPERM <- as.integer(Sys.getenv("SPEC20_NPERM", "5000")); set.seed(20)
# ---- all replicate profiles and their correlations (label independent) ----
Rall <- do.call(cbind, lapply(LINES, function(l) P[[l]]$R)); rep_line <- unlist(lapply(LINES, function(l) rep(l, ncol(P[[l]]$R)))); rep_group <- unlist(lapply(LINES, function(l) paste(P[[l]]$ko, l)))
CorRR <- cor(Rall, method = METHOD); CorRG <- cor(Rall, avgP, method = METHOD); CG <- cor(avgP, method = METHOD); same_line <- outer(units$line, units$line, "==")
lines_of <- setNames(lapply(LINES, function(l) which(units$line == l)), LINES)
# labelling: lab = KO label of every unit. Builds the average and best-replicate correlation matrices for that labelling.
mats_for <- function(lab) {
  CAl <- CG; CAl[same_line] <- NA
  chosen <- integer(n)
  for (u in seq_len(n)) { reps <- which(rep_group == units$id[u]); k <- lab[u]; others <- which(lab == k & units$line != units$line[u])
    chosen[u] <- if (!length(others) || length(reps) == 1) reps[1] else reps[which.max(rowMeans(CorRG[reps, others, drop = FALSE]))] }
  CBl <- CorRR[chosen, chosen]; CBl[same_line] <- NA; list(A = CAl, B = CBl, chosen = chosen) }
ko_stat <- function(C, lab, k) { i <- which(lab == k); o <- which(lab != k); sm <- C[i, i][upper.tri(C[i, i])]; sm <- sm[!is.na(sm)]; ot <- as.vector(C[i, o]); ot <- ot[!is.na(ot)]
  list(same = sm, other = ot, diff = if (length(sm)) median(sm) - median(ot) else NA_real_) }
lab0 <- units$KO; M0 <- mats_for(lab0); CA <- M0$A; CB <- M0$B; pick_lab <- sub(".*_", "", colnames(Rall)[M0$chosen])
kos_stat <- kos_all[kos_all %in% units$KO]
tab <- do.call(rbind, lapply(c("A", "B"), function(m) do.call(rbind, lapply(kos_stat, function(k) { st <- ko_stat(M0[[m]], lab0, k)
  data.frame(version = ifelse(m == "A", "averaged", "best replicate"), KO = k, n_same_pairs = length(st$same), median_same = if (length(st$same)) median(st$same) else NA, median_other = median(st$other), mean_same = if (length(st$same)) mean(st$same) else NA,
             p_wilcoxon = if (length(st$same) >= 2) wilcox.test(st$same, st$other, alternative = "greater", exact = FALSE)$p.value else NA, diff_obs = st$diff) }))))
DIST <- lapply(c(A = "A", B = "B"), function(m) setNames(lapply(kos_stat, function(k) ko_stat(M0[[m]], lab0, k)), kos_stat))     # per KO: same-KO and other-KO correlations (for the boxplots)
# label permutation: KO labels shuffled among the KOs present in each line; replicate choice redone under every shuffle
perm_diff <- replicate(NPERM, { lab <- lab0; for (l in LINES) { ix <- lines_of[[l]]; lab[ix] <- lab0[ix][sample(length(ix))] }; M <- mats_for(lab)
  c(sapply(kos_stat, function(k) ko_stat(M$A, lab, k)$diff), sapply(kos_stat, function(k) ko_stat(M$B, lab, k)$diff)) })
nA <- length(kos_stat); tab$p_perm <- NA_real_
for (i in seq_len(nrow(tab))) { j <- if (tab$version[i] == "averaged") match(tab$KO[i], kos_stat) else nA + match(tab$KO[i], kos_stat); tab$p_perm[i] <- if (is.na(tab$diff_obs[i])) NA else (1 + sum(perm_diff[j, ] >= tab$diff_obs[i], na.rm = TRUE)) / (NPERM + 1) }
tab$p_perm_BH <- ave(tab$p_perm, tab$version, FUN = function(p) p.adjust(p, "BH")); write.csv(tab, file.path(OUT, paste0("spec20_block_summary", SUF, ".csv")), row.names = FALSE)
options(width = 200); cat(sprintf("configuration: %s | %s | %d genes/features | %d permutations\n", SPACE, TYPE, nrow(avgP), NPERM)); print(format(tab[, c("version", "KO", "n_same_pairs", "median_same", "median_other", "p_wilcoxon", "p_perm", "p_perm_BH")], digits = 3), row.names = FALSE)

# ---- order KOs from top to bottom by the one-sided Wilcoxon p (smallest first), separately in each panel ----
n <- nrow(units); kcol <- KO_COL; short <- function(l) sub("NCI-", "", sub("-211H", "", l)); fp <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
ordv <- function(ver) { tb <- tab[tab$version == ver, ]; o <- tb$KO[order(ifelse(is.na(tb$p_wilcoxon), Inf, tb$p_wilcoxon), tb$p_perm)]; ou <- unlist(lapply(o, function(k) which(units$KO == k)))
  list(ordk = o, ordu = ou, U = units[ou, ], st = cumsum(c(1, table(factor(units$KO[ou], o))))[seq_along(o)], en = cumsum(table(factor(units$KO[ou], o)))) }
OA <- ordv("averaged"); OB <- ordv("best replicate"); CAo <- CA[OA$ordu, OA$ordu]; CBo <- CB[OB$ordu, OB$ordu]
lab_A <- paste(OA$U$KO, short(OA$U$line)); lab_B <- paste0(OB$U$KO, " ", short(OB$U$line), " (", pick_lab[OB$ordu], ")")
save_fig(paste0("fig24_best_config_heatmap", SUF), 17, 5.4, function() {
  layout(matrix(1:2, nrow = 1)); par(oma = c(3.2, 0, 3.4, 0.3), family = "sans"); MOD <- SPACE == "modules"; VM <- if (MOD) 0.8 else 0.5; RLO <- if (MOD) -1 else -0.5; RHI <- if (MOD) 1 else 0.6
  panel <- function(C, ver, m, ttl, labs, O) {
    ordk <- O$ordk; starts <- O$st; ends <- O$en
    par(mar = c(2.4, 6.3, 3.4, 20.5)); vmax <- VM; image(seq_len(n), seq_len(n), heat_z(pmax(pmin(C, vmax), -vmax)), col = PAL, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
    na <- which(is.na(C), arr.ind = TRUE); rect(na[, 2] - 0.5, n - na[, 1] + 0.5, na[, 2] + 0.5, n - na[, 1] + 1.5, col = "#efefef", border = NA)
    axis(2, at = n:1, labels = labs, las = 1, tick = FALSE, cex.axis = 0.4, line = -0.6)
    unit_in <- par("pin")[1] / n; xin <- function(inch) n + 0.5 + inch / unit_in; xr <- function(r) xin(0.95 + (pmax(pmin(r, RHI), RLO) - RLO) / (RHI - RLO) * 1.75); tb <- tab[tab$version == ver, ]
    segments(xr(0), 0.5, xr(0), n + 0.5, col = GRID, lwd = 0.8, xpd = NA)
    bx <- function(v, yc, hh, col) { if (!length(v)) return(invisible()); q <- quantile(v, c(0.25, 0.5, 0.75)); iq <- q[3] - q[1]; lo <- min(v[v >= q[1] - 1.5 * iq]); hi <- max(v[v <= q[3] + 1.5 * iq])
      segments(xr(lo), yc, xr(hi), yc, col = INK2, lwd = 0.7, xpd = NA); rect(xr(q[1]), yc - hh, xr(q[3]), yc + hh, col = adjustcolor(col, 0.9), border = INK2, lwd = 0.7, xpd = NA); segments(xr(q[2]), yc - hh, xr(q[2]), yc + hh, col = INK, lwd = 1.6, xpd = NA)
      points(xr(v), yc + runif(length(v), -hh * 0.75, hh * 0.75), pch = 16, cex = if (length(v) > 20) 0.16 else 0.32, col = adjustcolor(INK, 0.45), xpd = NA) }
    for (i in seq_along(ordk)) { ytop <- n - starts[i] + 1.5; ybot <- n - ends[i] + 0.5; h <- ytop - ybot; k <- ordk[i]; r <- tb[tb$KO == k, ]; st <- DIST[[m]][[k]]
      rect(starts[i] - 0.5, ybot, ends[i] + 0.5, ytop, border = INK, lwd = 1.0); rect(xin(0.02), ybot + 0.05, xin(0.09), ytop - 0.05, col = kcol[k], border = NA, xpd = NA); text(xin(0.13), (ytop + ybot) / 2, k, adj = 0, cex = 0.62, font = 2, xpd = NA)
      bx(st$same, ytop - h * 0.27, h * 0.15, kcol[k]); bx(st$other, ybot + h * 0.27, h * 0.15, GREY)
      text(xin(3.4), (ytop + ybot) / 2, fp(r$p_wilcoxon), adj = 1, cex = 0.6, font = ifelse(!is.na(r$p_wilcoxon) & r$p_wilcoxon < 0.05, 2, 1), xpd = NA); text(xin(3.98), (ytop + ybot) / 2, fp(r$p_perm), adj = 1, cex = 0.6, font = ifelse(!is.na(r$p_perm) & r$p_perm < 0.05, 2, 1), xpd = NA) }
    tk <- if (MOD) c(-0.8, 0, 0.8) else c(-0.4, 0, 0.4); segments(xr(tk), 0.5, xr(tk), 0.1, col = INK2, lwd = 0.6, xpd = NA); text(xr(tk), -0.9, tk, cex = 0.5, col = INK2, xpd = NA); text(xr(0.1), -2.3, "cross-line r (colour = same KO, grey = other KOs)", cex = 0.5, col = INK2, xpd = NA)
    text(xin(c(3.4, 3.98)), n + 1.6, c("p\nWilcoxon", "p\npermutation"), adj = 1, cex = 0.5, col = INK2, xpd = NA); text(xr(0.05), n + 1.6, "median and IQR of the\ncorrelations between lines", cex = 0.5, col = INK2, xpd = NA)
    mtext(ttl, side = 3, adj = 0, line = 1.7, font = 2, cex = 0.85) }
  panel(CAo, "averaged", "A", "A  Replicates averaged", lab_A, OA); panel(CBo, "best replicate", "B", "B  Best replicate per KO and line (selection-biased)", lab_B, OB)
  corlab <- if (METHOD == "spearman") "Spearman rho" else "Pearson r"
  draw_colorbar(c(0.38, 0.62), c(-VM, VM), sprintf("%s between the KO profiles of two cell lines", corlab), cex = 0.55)
  mtext(sprintf("Cross-cell-line reproducibility of bulk KO profiles (%s features, %s statistic, centred on each line%s%s)", if (MOD) "20 Cm module" else ifelse(SPACE == "all_genes", ifelse(GENESET == "screen", "CROP-seq-panel", "all"), SPACE), TYPE, if (KOSET == "screen") "; the 8 CROP-seq KOs" else "", if (ZSCORE) "; z-scored per feature" else ""), outer = TRUE, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext(sprintf("%d features; heatmap boxes = same KO in different lines, grey cells = same line; boxplots: cross-line correlations of the KO with itself (colour) vs with the other KOs (grey); KOs ordered by the Wilcoxon p (smallest at the top) within each panel; one-sided Wilcoxon p (pairs share samples) and label-permutation p (%d shuffles; for B the replicate choice is redone each time)", nrow(avgP), NPERM), outer = TRUE, side = 3, adj = 0, line = 0.7, cex = 0.5, col = INK2) })
