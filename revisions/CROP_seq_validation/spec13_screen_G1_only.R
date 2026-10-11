# Spec 13: remove ALL proliferating cells from the CROP-seq screen (Seurat Phase S or G2M; keep G1 only), recompute the screen KO profiles from the remaining cells,
# and correlate them with the cell-cycle-adjusted bulk (DESeq2 ~ cc + condition statistics / cc-regressed log differences), as in the best variants (F1-F4, H1-H2).
# Screen versions compared: (A) all cells, cc regressed (previous best)  (B) all cells, no regression  (C) G1 cells only  (D) G1 cells only + residual cc regressed.
# Bulk is cc-adjusted in every version. Columns: the 8 shared screen KOs, or those with >= 50 cells IN THE VERSION USED (so the powered set can differ).
# needs: spec1_prepare_stats.R
# outputs: spec13_phase_counts.csv, spec13_G1_only_results.csv, fig17_G1_only_screen.{png,pdf}, fig17b_G1_only_matrix.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_variants.R"))
ctx <- load_spec_context(LINES); ss <- ctx$ss; bs <- ctx$bs
srt <- load_screen(); md <- srt@meta.data
ph <- table(md$merged_call, md$Phase); phc <- data.frame(KO = rownames(ph), G1 = ph[, "G1"], S = ph[, "S"], G2M = ph[, "G2M"], total = rowSums(ph)); phc$frac_proliferating <- 1 - phc$G1 / phc$total
write.csv(phc, file.path(OUT, "spec13_phase_counts.csv"), row.names = FALSE); options(width = 200); print(format(phc, digits = 3), row.names = FALSE)
keep <- md$Phase == "G1"; cat(sprintf("\nkeeping %d of %d cells (G1); removing %d proliferating cells (%.0f%%)\n", sum(keep), length(keep), sum(!keep), 100 * mean(!keep)))
Y <- GetAssayData(srt, layer = "data")[ss$genes, ]; grp <- md$merged_call; ccv <- md$cc
G1 <- screen_ko_stats(Y[, keep], grp[keep], ccv[keep], FALSE); G1r <- screen_ko_stats(Y[, keep], grp[keep], ccv[keep], TRUE)
n_g1 <- G1$n; cat("cells per KO after removal:", paste(names(n_g1), n_g1, collapse = " | "), "\n")
SV <- list(`A all cells, cc regressed` = list(t = ss$tcc, diff = ss$diffcc, n = ss$n), `B all cells, no regression` = list(t = ss$t0, diff = ss$diff0, n = ss$n),
           `C G1 cells only` = list(t = G1$t, diff = G1$diff, n = n_g1), `D G1 only + cc regressed` = list(t = G1r$t, diff = G1r$diff, n = n_g1))
genes_all <- ctx$F_all; genes_var <- ctx$F_var

## ---------------------------------------------------------------- bulk cc-adjusted inputs
bulk_stat <- function(feats) setNames(lapply(LINES, function(cl) { x <- bs[[cl]]$statcc[feats, , drop = FALSE]; x[!is.finite(x)] <- 0; x }), LINES)
bulk_diff_cc <- function(feats) setNames(lapply(LINES, function(cl) { m <- ctx$meta[ctx$meta$Cell_Line == cl, ]; y <- ctx$lg[feats, rownames(m), drop = FALSE]; cc <- ctx$ccb$cc[match(rownames(m), ctx$ccb$sample)]
  b <- as.numeric((y - rowMeans(y)) %*% (cc - mean(cc))) / sum((cc - mean(cc))^2); y <- y - outer(b, cc); ntc <- rowMeans(y[, m$Gene_Control == "NTC", drop = FALSE]); kos <- sort(setdiff(unique(m$Gene_Control), "NTC"), method = "radix")
  sapply(setNames(kos, kos), function(k) rowMeans(y[, m$Gene_Control == k, drop = FALSE]) - ntc) }), LINES)
row_sets <- list(ALL3 = function(lab) lab$line %in% THREE, ALL4 = function(lab) rep(TRUE, nrow(lab)), H2052 = function(lab) lab$line == "NCI-H2052")
res <- list(); mats <- list()
for (sv in names(SV)) for (fs in c("all_genes", "var_genes")) for (ty in c("t", "diff")) for (cen in c(FALSE, TRUE)) {
  feats <- if (fs == "all_genes") genes_all else genes_var; S <- SV[[sv]][[ty]][feats, ]; Bl <- if (ty == "t") bulk_stat(feats) else bulk_diff_cc(feats)
  if (ty == "diff") S <- S[, colnames(S) != "NTC", drop = FALSE]
  cm <- cross_matrix(S, Bl, cen, ty == "diff"); mats[[paste(sv, fs, ty, cen)]] <- cm
  for (cs in c("all8", "powered")) { cols <- if (cs == "all8") ORDER else intersect(ORDER, names(SV[[sv]]$n)[SV[[sv]]$n >= 50])
    for (rs in names(row_sets)) { sel <- row_sets[[rs]](cm$lab) & cm$lab$ko %in% cols; if (sum(sel) < 3) next; Cm <- cm$C; colnames(Cm) <- ORDER
      res[[length(res) + 1]] <- data.frame(screen = sv, features = fs, stat = ty, centred = cen, cols = cs, n_cols = length(cols), rows = rs, evaluate_specificity(Cm[sel, , drop = FALSE], cm$lab$ko[sel], cols = cols)) } } }
res <- do.call(rbind, res); rownames(res) <- NULL; write.csv(res, file.path(OUT, "spec13_G1_only_results.csv"), row.names = FALSE)
cat("\n=== best per screen version (max z over features/stat/centring), by column set and row set ===\n")
best <- do.call(rbind, lapply(split(res, list(res$screen, res$cols, res$rows), drop = TRUE), function(d) d[which.max(d$z_matched), ])); best <- best[order(best$cols, best$rows, best$screen), ]
print(format(best[, c("screen", "cols", "rows", "features", "stat", "centred", "n_rows", "matched_r", "other_r", "z_matched", "perm_p", "top1", "median_rank")], digits = 3), row.names = FALSE)
cat("\n=== mean z by screen version (over all feature/stat/centring variants) ===\n"); print(round(tapply(res$z_matched, list(res$screen, paste(res$cols, res$rows)), mean), 3))
cat("\ntests:", nrow(res), "| p<0.05:", sum(res$perm_p < 0.05), "(expected ~", round(0.05 * nrow(res), 1), ")\n")

## ---------------------------------------------------------------- figure: z by screen version
save_fig("fig17_G1_only_screen", 13, 8.4, function() {
  layout(matrix(1:2, nrow = 1), widths = c(1, 1)); par(oma = c(0, 0, 4, 0.5), family = "sans"); vcol <- c("#2a78d6", GREY, "#eb6834", "#b3400f")
  for (cs in c("all8", "powered")) { d <- res[res$cols == cs & res$rows %in% c("ALL3", "ALL4"), ]; svs <- names(SV); rsn <- c("ALL3", "ALL4")
    par(mar = c(7.5, 4.6, 3.6, 0.8)); plot(NA, xlim = c(0.4, length(svs) * 2 + 0.6), ylim = range(c(d$z_matched, 0)) + c(-0.05, 0.12), xaxt = "n", xlab = "", ylab = "row-standardised matched score z", bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8); abline(h = 0, col = INK2)
    for (i in seq_along(svs)) for (j in seq_along(rsn)) { v <- d$z_matched[d$screen == svs[i] & d$rows == rsn[j]]; at <- (i - 1) * 2 + j; boxplot(v, at = at, add = TRUE, axes = FALSE, boxwex = 0.7, outline = FALSE, col = adjustcolor(vcol[i], ifelse(j == 1, 0.85, 0.45)), border = INK2, medcol = INK, medlwd = 1.6, whisklty = 1, staplelty = 1)
      points(at + runif(length(v), -0.2, 0.2), v, pch = 16, cex = 0.4, col = adjustcolor(INK, 0.5)); pv <- d$perm_p[d$screen == svs[i] & d$rows == rsn[j]]; points(at, max(v), pch = ifelse(min(pv) < 0.05, 8, NA), cex = 0.9) }
    axis(1, at = seq(1.5, length(svs) * 2, 2), labels = sub(" ", "\n", svs), tick = FALSE, cex.axis = 0.62, line = 0.6)
    mtext(sprintf("%s (n cols: %s)", ifelse(cs == "all8", "all 8 shared screen KOs", "screen KOs with >= 50 cells in that version"), paste(sapply(svs, function(s) length(if (cs == "all8") ORDER else intersect(ORDER, names(SV[[s]]$n)[SV[[s]]$n >= 50]))), collapse = "/")), side = 3, adj = 0, line = 0.7, font = 2, cex = 0.8) }
  mtext("Removing proliferating cells from the CROP-seq screen before correlating with cc-adjusted bulk", outer = TRUE, side = 3, adj = 0, line = 2.4, font = 2, cex = 1)
  mtext(sprintf("each box = %d feature/statistic/centring variants (dark = 3 lines, light = 4 lines); G1-only keeps %d of %d cells; * = the best variant has exact permutation p < 0.05", length(unique(paste(res$features, res$stat, res$centred))), sum(keep), length(keep)), outer = TRUE, side = 3, adj = 0, line = 1.1, cex = 0.6, col = INK2)
})
## matrix of the best G1-only variant (all 8 shared TFs ordered best -> worst + screen-only KOs MEF2A / MEF2D, diagonal ordered)
ALLC <- c(ORDER, EXTRA)
g <- res[res$screen == "C G1 cells only" & res$cols == "all8" & res$rows == "ALL4", ]; g <- g[which.max(g$z_matched), ]
feats <- if (g$features == "all_genes") genes_all else genes_var; S <- SV[["C G1 cells only"]][[g$stat]][feats, ]; Bl <- if (g$stat == "t") bulk_stat(feats) else bulk_diff_cc(feats); if (g$stat == "diff") S <- S[, colnames(S) != "NTC", drop = FALSE]
cm <- cross_matrix(S, Bl, g$centred, g$stat == "diff", cols = ALLC); Call <- cm$C; colnames(Call) <- ALLC; lab <- cm$lab
Z <- t(scale(t(Call[, ORDER]))); Z[!is.finite(Z)] <- 0; jm <- match(lab$ko, ORDER); perf <- sapply(ORDER, function(k) mean(Z[cbind(seq_len(nrow(Z)), jm)][lab$ko == k])); cols <- names(sort(perf, decreasing = TRUE))
lab$ko <- factor(lab$ko, cols); o <- order(match(lab$line, LINES), as.integer(lab$ko)); lab <- lab[o, ]; lab$ko <- as.character(lab$ko); C <- Call[o, cols]; Cx <- Call[o, EXTRA, drop = FALSE]
draw_ordered_matrix("fig17b_G1_only_matrix", C, Cx, lab, n_g1, lab$ko, "CROP-seq G1 cells only vs cc-adjusted bulk (replicates averaged)",
  sprintf("%s genes, %s statistic, %s; all 8 shared TFs ordered best to worst, then the screen-only KOs MEF2A / MEF2D; cells per KO in brackets (G1 only)", g$features, ifelse(g$stat == "t", "test", "log-difference"), ifelse(g$centred, "centred across KOs", "uncentred")),
  sprintf("all 4 lines, 8 TFs: z = %.2f (exact permutation p = %.3f); best variant of this screen version, so optimistic", g$z_matched, g$perm_p))
