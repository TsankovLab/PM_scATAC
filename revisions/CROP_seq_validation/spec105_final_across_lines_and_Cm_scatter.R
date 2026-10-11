# Spec 105: (1) the retained across-line plot: each KO vs the mean of the 9 groups (NTC + 8 KOs) of its line (limma moderated t, no NTC subtraction; spec104 "g9_t"), top-500 KO-signal genes, Spearman,
#               TFs as rows sorted from the most to the least significant (KO-label permutation p)
#           (2) NCI-H2052, CROP-seq vs bulk: mean Cm module score of each KO (AddModuleScore, top-50 genes; z within data set), relative to the mean of the 9 groups, 20 Cm x 8 KOs; all points and one panel per KO
# needs: spec72 (S70), spec75 (screen scores)
# outputs: spec105_across_lines_9group_t_per_KO.csv, spec105_H2052_Cm_KO_means.csv, fig19_across_lines_9group_sorted.{png,pdf}, fig20_H2052_Cm_KO_scatter.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000")); permsL <- make_perms(N_PERM); pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); fab <- load_full_audited_bulk(); meta <- fab$meta; ids <- pr$samples; GR <- c("NTC", KOS); CL <- "NCI-H2052"
fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
T9 <- array(NA_real_, c(length(pr$genes), 8, 4), dimnames = list(pr$genes, KOS, LINES))
for (cl in LINES) { s <- ids[meta[ids, "Cell_Line"] == cl]; y <- normalise(fab$cs[pr$genes, s], "tmm")$y; g <- factor(meta[s, "Gene_Control"], levels = GR); des <- model.matrix(~ 0 + g); colnames(des) <- GR
  L <- sapply(KOS, function(k) { v <- setNames(rep(1 / 9, 9), GR); v[k] <- v[k] - 1; v }); T9[, , cl] <- limma::eBayes(limma::contrasts.fit(limma::lmFit(y, des), L))$t }
a <- agreement(cl_from_array(T9, "signal", pr$gene_signal, centre = FALSE), permsL); k <- a$per_KO[order(a$per_KO$perm_p), ]; write.csv(k, file.path(OUT, "spec105_across_lines_9group_t_per_KO.csv"), row.names = FALSE); options(width = 200); print(k[, c("KO", "mean_intra", "median_cross", "perm_p", "perm_fdr")], digits = 3, row.names = FALSE)
save_fig("fig19_across_lines_9group_sorted", 7.4, 8.6, function() { o <- k$KO; n <- length(o); par(mar = c(4.4, 6.2, 4.6, 7.4), family = "sans"); xr <- range(c(unlist(a$cross), a$intra))
  plot(NA, xlim = xr, ylim = c(n + 0.6, 0.4), yaxt = "n", bty = "n", las = 1, ylab = "", xlab = "Spearman between two bulk cell lines", cex.axis = 0.85); abline(v = 0, lty = 2, col = INK)
  ## boxwex 0.5 at +/-0.25 is the geometric maximum: the two boxes meet in the middle and their
  ## outer edges reach the next row. Extra thickness comes from the taller canvas, not a bigger boxwex.
  for (i in seq_len(n)) { j <- match(o[i], a$per_KO$KO); v <- a$intra[o[i], ]; boxplot(a$cross[[j]], at = i + 0.25, add = TRUE, horizontal = TRUE, boxwex = 0.5, col = "#d5d4cf", border = INK2, outline = FALSE, axes = FALSE)
    boxplot(v, at = i - 0.25, add = TRUE, horizontal = TRUE, boxwex = 0.5, col = adjustcolor(KO_COL[[o[i]]], 0.7), border = INK2, outline = FALSE, axes = FALSE); points(v, i - 0.25 + seq(-0.15, 0.15, length.out = length(v)), pch = 16, cex = 0.75, col = INK2)
    text(xr[2] + 0.03 * diff(xr), i, sprintf("p = %s%s", fp(k$perm_p[i]), ifelse(k$perm_fdr[i] < 0.01, " **", ifelse(k$perm_fdr[i] < 0.05, " *", ""))), adj = 0, xpd = NA, cex = 0.8, font = ifelse(k$perm_p[i] < 0.05, 2, 1)) }
  text(xr[1] - 0.04 * diff(xr), seq_len(n), o, adj = 1, xpd = NA, font = 2, cex = 0.95, col = KO_COL[o]); legend(xr[1], 0.55, c("same KO (6 line pairs; TF colour)", "other KOs"), fill = c("grey45", "#d5d4cf"), border = NA, bty = "n", cex = 0.68, xpd = NA, horiz = TRUE, yjust = 0.5)
  mtext("Reproducibility of each TF KO across the 4 bulk cell lines", side = 3, adj = 0, line = 3, font = 2, cex = 1, at = xr[1] - 0.2 * diff(xr))
  mtext(sprintf("each KO vs the mean of the 9 groups (NTC + 8 KOs) of its line, limma t; top-500 KO-signal genes; all KOs: p = %s", fp(a$global[["perm_p"]])), side = 3, adj = 0, line = 1.8, cex = 0.62, col = INK2, at = xr[1] - 0.2 * diff(xr))
  mtext("KO-label permutation p (20,000 permutations within lines); * / ** FDR < 0.05 / 0.01; S70", side = 3, adj = 0, line = 0.8, cex = 0.62, col = INK2, at = xr[1] - 0.2 * diff(xr)) })

## ---------------------------------------------------------------- (2) H2052: Cm x KO means, CROP-seq vs bulk
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); CM <- paste0("Cm", 1:20); sH <- ids[meta[ids, "Cell_Line"] == CL]; gH <- meta[sH, "Gene_Control"]; gs <- as.character(sc$group); cs_ <- gs %in% GR
zS <- scale(sc$score[cs_, CM]); zB <- t(scale(t(pr$cm_scores$tmm[CM, sH])))                                    # z per module within the data set (screen: cells of the 9 groups; bulk: H2052 samples)
mS <- sapply(GR, function(g) colMeans(zS[gs[cs_] == g, , drop = FALSE])); mB <- sapply(GR, function(g) rowMeans(zB[, gH == g, drop = FALSE])); S <- (mS - rowMeans(mS))[, KOS]; B <- (mB - rowMeans(mB))[, KOS]   # KO mean minus mean of the 9 groups
D <- data.frame(Cm = rep(CM, 8), KO = rep(KOS, each = 20), cropseq = as.vector(S), bulk = as.vector(B)); write.csv(D, file.path(OUT, "spec105_H2052_Cm_KO_means.csv"), row.names = FALSE)
ct <- function(x, y) { r <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)); c(r$estimate, r$p.value) }; all_r <- ct(D$cropseq, D$bulk); perKO <- t(sapply(KOS, function(kk) ct(S[, kk], B[, kk]))); perCm <- t(sapply(CM, function(m) ct(S[m, ], B[m, ])))
cat(sprintf("\nH2052 Cm x KO (160 points): Spearman %.2f, p = %.3g\nper KO (over 20 Cm):\n", all_r[1], all_r[2])); print(round(perKO, 3)); cat("per Cm (over 8 KOs): median rho", round(median(perCm[, 1]), 2), "; modules with rho > 0.5:", paste(CM[perCm[, 1] > 0.5], collapse = " "), "\n")
save_fig("fig20_H2052_Cm_KO_scatter", 14, 7.6, function() {
  layout(cbind(matrix(1, 2, 2), matrix(2:9, nrow = 2, byrow = TRUE)), widths = c(1, 1, rep(0.82, 4))); par(oma = c(0, 0, 3.2, 0.4), family = "sans", mar = c(4.6, 4.8, 3.2, 1))
  plot(D$cropseq, D$bulk, pch = 21, bg = adjustcolor(KO_COL[D$KO], 0.85), col = INK2, cex = 1.4, bty = "n", las = 1, xlab = "CROP-seq: Cm score of the KO, relative to the mean of the 9 groups (z)", ylab = "bulk H2052: Cm score of the KO, relative to the mean of the 9 groups (z)", cex.axis = 0.85)
  abline(h = 0, v = 0, lty = 2, col = GRID); abline(lm(bulk ~ cropseq, D), col = INK2, lty = 3); far <- order(-abs(scale(D$cropseq)) - abs(scale(D$bulk)))[1:14]; text(D$cropseq[far], D$bulk[far], D$Cm[far], pos = 3, cex = 0.6, col = KO_COL[D$KO[far]], xpd = NA)
  legend("topleft", KOS, pch = 21, pt.bg = KO_COL[KOS], col = INK2, bty = "n", cex = 0.75, pt.cex = 1.2); mtext("all 20 Cm x 8 KOs", side = 3, adj = 0, line = 1.5, font = 2, cex = 0.85); mtext(sprintf("Spearman %.2f, p = %s (160 points)", all_r[1], fp(all_r[2])), side = 3, adj = 0, line = 0.4, cex = 0.65, col = INK2)
  o <- KOS[order(-perKO[, 1])]; for (kk in o) { par(mar = c(3.8, 3.6, 3.2, 0.6)); plot(S[, kk], B[, kk], pch = 21, bg = adjustcolor(KO_COL[[kk]], 0.85), col = INK2, cex = 1.2, bty = "n", las = 1, xlab = "", ylab = "", cex.axis = 0.7); abline(h = 0, v = 0, lty = 2, col = GRID); abline(lm(B[, kk] ~ S[, kk]), col = INK2, lty = 3)
    text(S[, kk], B[, kk], sub("Cm", "", CM), pos = 3, cex = 0.5, col = INK2, xpd = NA); mtext(kk, side = 3, adj = 0, line = 1.4, font = 2, cex = 0.8, col = KO_COL[[kk]]); mtext(sprintf("rho %.2f, p %s", perKO[kk, 1], fp(perKO[kk, 2])), side = 3, adj = 0, line = 0.4, cex = 0.6, col = INK2, font = ifelse(perKO[kk, 2] < 0.05, 2, 1))
    mtext("CROP-seq", side = 1, line = 2.4, cex = 0.55); mtext("bulk H2052", side = 2, line = 2.4, cex = 0.55) }
  mtext("NCI-H2052: mean Cm module score of each KO, CROP-seq vs bulk RNA-seq", outer = TRUE, side = 3, adj = 0, line = 1.7, font = 2, cex = 1)
  mtext("AddModuleScore (top-50 module genes), z-scored per module within each data set, KO mean minus the mean of the 9 groups (NTC + 8 KOs); > 0 = module higher in that KO; right: one panel per KO, numbers = Cm", outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.65, col = INK2)
})
