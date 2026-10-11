# Spec 100: figures for the settings retained after specs 98-99, and Cm20 in the TEAD4 KO.
# fig13  A  NCI-H2052, CROP-seq (rows) vs bulk (columns): KO x KO Spearman on the genes in the top 2000 variable genes of BOTH data sets, profiles centred on the mean of the 9 groups (NTC + 8 KOs)
#        B  own-KO Spearman per TF for the top-1000 and top-2000 intersections, # = rank of the own KO among the 8 bulk KOs, exact permutation p of the bulk KO labels (8!) for all KOs
#        C / D  bulk across the 4 lines, primary configuration (S70, TMM + limma t of NTC - KO, KO-centred, top-500 KO-signal genes from the 2 lines not compared, Spearman; spec76): same KO vs other KOs
#               per TF with KO-label permutation p, no covariate / cell cycle regressed
# fig14  Cm20 (AddModuleScore, top-50 genes; z-scored within each data set) in NTC, TEAD4 KO and the other KOs: CROP-seq cells and the bulk samples of each line; NTC - KO per KO and data set
# needs: spec72 (S70), spec75 (screen scores), spec76 (SET=S70), spec82, spec_common.R
# outputs: fig13_final_settings.{png,pdf}, fig14_Cm20_TEAD4.{png,pdf}, spec100_H2052_KOxKO_int2000.csv, spec100_Cm20_TEAD4.csv, spec100_Cm20_NTCminusKO.csv
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages(library(Seurat)); CL <- "NCI-H2052"; GR <- c("NTC", KOS); fab <- load_full_audited_bulk(); meta <- fab$meta; pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); ids <- pr$samples
rd <- function(f) read.csv(file.path(OUT, f), stringsAsFactors = FALSE, check.names = FALSE); fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2)); SH <- function(l) sub("NCI-|-211H", "", l)
srt <- load_screen(); grp <- as.character(srt$merged_call); X <- GetAssayData(srt, layer = "counts")
PB <- sapply(GR, function(k) Matrix::rowSums(X[, grp == k, drop = FALSE])); lcS <- edgeR::cpm(edgeR::DGEList(PB), log = TRUE, prior.count = 2); expS <- rownames(lcS)[rowSums(edgeR::cpm(edgeR::DGEList(PB)) >= 10) >= 1]
G <- intersect(pr$genes, expS); hvS <- intersect(VariableFeatures(FindVariableFeatures(srt, selection.method = "vst", nfeatures = nrow(srt), verbose = FALSE)), G)
sH <- ids[meta[ids, "Cell_Line"] == CL]; yH <- normalise(fab$cs[pr$genes, sH], "tmm")$y; bH <- sapply(GR, function(k) rowMeans(yH[, meta[sH, "Gene_Control"] == k, drop = FALSE]))
c9 <- function(M) (M[, GR] - rowMeans(M[, GR]))[, KOS]; vH <- G[order(pr$gene_var[G, CL], decreasing = TRUE)]
perms8 <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(1:8)) })
RES <- lapply(c(int1000 = 1000, int2000 = 2000), function(n) { f <- intersect(hvS[seq_len(n)], vH[seq_len(n)]); C <- cor(c9(lcS[f, ]), c9(bH[f, ]), method = "spearman"); dm <- mean(diag(C))
  list(n = length(f), C = C, p = mean(apply(perms8, 1, function(p) mean(C[cbind(1:8, p)])) >= dm - 1e-12), rank = sapply(1:8, function(i) rank(-C[i, ])[i]), dmo = dm - mean(C[row(C) != col(C)])) })
write.csv(RES$int2000$C, file.path(OUT, "spec100_H2052_KOxKO_int2000.csv"))
pk <- rd("spec76_primary_per_KO_S70.csv"); pk <- pk[pk$space == "signal500", ]; gl <- rd("spec76_primary_global_S70.csv"); gl <- gl[gl$space == "signal500", ]
pairs_of <- function(cv) { z <- rd(sprintf("spec76_KOline_correlation_signal500_%s_S70.csv", cv)); M <- as.matrix(z[, -1]); rn <- z[, 1]; ko <- sub(" .*", "", rn); ln <- sub(".* ", "", rn); ut <- which(upper.tri(M) & outer(ln, ln, "!="), arr.ind = TRUE); data.frame(a = ko[ut[, 1]], b = ko[ut[, 2]], rho = M[ut]) }
PR <- list(none = pairs_of("none"), cc = pairs_of("cc")); COVL <- c(none = "no covariate", cc = "cell cycle regressed"); o2 <- rd("spec82_TF_agreement_ranked.csv")$KO

save_fig("fig13_final_settings", 13, 10.6, function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4.2, 0.4), family = "sans"); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.82); mtext(b, side = 3, adj = 0, line = 0.4, cex = 0.6, col = INK2) }
  ## A heatmap
  C <- RES$int2000$C; zl <- max(abs(C)); par(mar = c(5.4, 6.4, 3.4, 5)); pal <- colorRampPalette(c("#2a78d6", "#f7f6f2", "#eb6834"))(101)
  image(1:8, 1:8, t(C[8:1, ]), col = pal, zlim = c(-zl, zl), axes = FALSE, xlab = "", ylab = ""); text(rep(1:8, each = 8), rep(8:1, 8), sprintf("%.2f", as.vector(C)), cex = 0.62, font = ifelse(as.vector(row(C) == col(C)), 2, 1)); rect(1:8 - 0.5, 8:1 - 0.5, 1:8 + 0.5, 8:1 + 0.5, border = INK, lwd = 1.4)
  text(1:8, 0.35, KOS, srt = 45, adj = 1, xpd = NA, cex = 0.8, font = 2, col = KO_COL[KOS]); text(0.4, 8:1, KOS, adj = 1, xpd = NA, cex = 0.8, font = 2, col = KO_COL[KOS]); mtext("bulk H2052 KO", side = 1, line = 4, cex = 0.7); mtext("CROP-seq KO", side = 2, line = 5, cex = 0.7)
  yy <- seq(2, 7, length.out = 102); rect(9, yy[-102], 9.35, yy[-1], col = pal, border = NA, xpd = NA); text(9.45, c(2, 4.5, 7), sprintf("%.1f", c(-zl, 0, zl)), adj = 0, xpd = NA, cex = 0.62); text(9, 7.4, "Spearman", adj = 0, xpd = NA, cex = 0.62)
  ttl("A  H2052: CROP-seq vs bulk, KO x KO", sprintf("%d genes in the top 2000 variable of both; 9-group centring; own KO on the diagonal", RES$int2000$n))
  ## B bars
  o <- KOS[order(-diag(RES$int2000$C))]; M <- rbind(diag(RES$int1000$C)[o], diag(RES$int2000$C)[o]); R <- rbind(RES$int1000$rank[match(o, KOS)], RES$int2000$rank[match(o, KOS)]); par(mar = c(4.6, 4.6, 3.4, 1))
  bp <- barplot(M, beside = TRUE, col = rbind(adjustcolor(KO_COL[o], 0.45), adjustcolor(KO_COL[o], 0.9)), border = NA, las = 1, names.arg = rep("", 8), ylab = "own-KO Spearman, CROP-seq vs bulk H2052", cex.axis = 0.8, ylim = range(M) + c(-0.06, 0.16)); abline(h = 0, col = INK)
  text(bp, pmax(M, 0) + 0.015, sprintf("#%s", format(R)), cex = 0.6, font = ifelse(R <= 1, 2, 1)); text(colMeans(bp), par("usr")[3] - 0.02 * diff(par("usr")[3:4]), o, srt = 45, adj = 1, xpd = NA, cex = 0.82, font = 2, col = KO_COL[o])
  legend("topright", sprintf("top %s in both (%d genes): all KOs p = %s", c("1000", "2000"), c(RES$int1000$n, RES$int2000$n), fp(c(RES$int1000$p, RES$int2000$p))), fill = c("grey75", "grey35"), border = NA, bty = "n", cex = 0.62)
  ttl("B  H2052 per TF: TEAD4, PITX1 and TEAD2 match", "light / dark = top-1000 / top-2000 intersection; # = rank of the own KO among the 8 bulk KOs")
  ## C, D
  for (cv in names(COVL)) { d <- PR[[cv]]; kk <- pk[pk$cov == cv, ]; yr <- range(d$rho) + c(-0.02, 0.16) * diff(range(d$rho)); par(mar = c(4.6, 4.6, 3.4, 1))
    plot(NA, xlim = c(0.5, 8.5), ylim = yr, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "Spearman between two bulk cell lines", cex.axis = 0.8); abline(h = 0, lty = 2, col = INK)
    for (i in seq_along(o2)) { w <- d$rho[d$a == o2[i] & d$b == o2[i]]; z <- d$rho[(d$a == o2[i]) != (d$b == o2[i])]
      boxplot(z, at = i - 0.2, add = TRUE, boxwex = 0.32, col = "#d5d4cf", border = INK2, outline = FALSE, axes = FALSE); boxplot(w, at = i + 0.2, add = TRUE, boxwex = 0.32, col = adjustcolor(KO_COL[[o2[i]]], 0.65), border = INK2, outline = FALSE, axes = FALSE)
      points(i + 0.2 + seq(-0.1, 0.1, length.out = length(w)), w, pch = 16, cex = 0.75, col = INK2); r <- kk[kk$KO == o2[i], ]
      text(i, yr[2], sprintf("p=%s\n%s", fp(r$perm_p), ifelse(r$perm_fdr < 0.01, "**", ifelse(r$perm_fdr < 0.05, "*", ""))), cex = 0.68, adj = c(0.5, 1), font = ifelse(r$perm_p < 0.05, 2, 1)) }
    text(1:8, par("usr")[3] - 0.03 * diff(par("usr")[3:4]), o2, srt = 45, adj = 1, xpd = NA, cex = 0.82, font = 2, col = KO_COL[o2]); gg <- gl[gl$cov == cv, ]
    ttl(sprintf("%s  across the 4 bulk lines (%s)", if (cv == "none") "C" else "D", COVL[[cv]]), sprintf("colour = same KO (6 line pairs), grey = other KOs; global p = %s; * / ** FDR < 0.05 / 0.01", fp(gg$perm_p))) }
  mtext("Retained settings.  Top: NCI-H2052 CROP-seq vs bulk (genes variable in both, 9-group centring).  Bottom: bulk across lines (top-500 KO-signal genes, KO-centred limma t)", outer = TRUE, side = 3, adj = 0, line = 2.6, font = 2, cex = 0.9)
  mtext("Spearman; S70 sample set; permutation of KO labels (exact 8! for A-B, 20,000 within-line permutations for C-D)", outer = TRUE, side = 3, adj = 0, line = 1.4, cex = 0.65, col = INK2)
})

## ---------------------------------------------------------------- Cm20 in TEAD4
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); DS <- c("CROP-seq", LINES)
dat <- rbind(data.frame(ds = "CROP-seq", group = as.character(sc$group), v = sc$score[, "Cm20"]), data.frame(ds = meta[ids, "Cell_Line"], group = meta[ids, "Gene_Control"], v = pr$cm_scores$tmm["Cm20", ids])); dat <- dat[dat$group %in% GR, ]
dat$z <- ave(dat$v, dat$ds, FUN = function(x) (x - mean(x)) / sd(x)); dat$cls <- ifelse(dat$group == "NTC", "NTC", ifelse(dat$group == "TEAD4", "TEAD4", "other"))
st <- do.call(rbind, lapply(DS, function(d) { a <- dat$z[dat$ds == d & dat$cls == "NTC"]; b <- dat$z[dat$ds == d & dat$cls == "TEAD4"]; o <- dat$z[dat$ds == d & dat$cls == "other"]
  data.frame(ds = d, n_NTC = length(a), n_TEAD4 = length(b), NTC_minus_TEAD4 = mean(a) - mean(b), p_vs_NTC = if (d == "CROP-seq") wilcox.test(a, b)$p.value else if (min(length(a), length(b)) >= 2) t.test(a, b)$p.value else NA,
             otherKO_minus_TEAD4 = mean(o) - mean(b), p_vs_otherKO = if (d == "CROP-seq") wilcox.test(o, b)$p.value else t.test(o, b)$p.value) }))
bk <- dat[dat$ds != "CROP-seq", ]; f1 <- summary(lm(z ~ ds + I(group == "TEAD4"), bk[bk$cls != "other", ]))$coef; f2 <- summary(lm(z ~ ds + I(group == "TEAD4"), bk[bk$cls != "NTC", ]))$coef
comb <- c(vs_NTC = f1[nrow(f1), 4], vs_other = f2[nrow(f2), 4]); est <- c(vs_NTC = -f1[nrow(f1), 1], vs_other = -f2[nrow(f2), 1])
E <- do.call(rbind, lapply(DS, function(d) data.frame(ds = d, KO = KOS, NTC_minus_KO = sapply(KOS, function(k) mean(dat$z[dat$ds == d & dat$group == "NTC"]) - mean(dat$z[dat$ds == d & dat$group == k]))))); E$rank_most_down <- ave(-E$NTC_minus_KO, E$ds, FUN = rank)
write.csv(st, file.path(OUT, "spec100_Cm20_TEAD4.csv"), row.names = FALSE); write.csv(E, file.path(OUT, "spec100_Cm20_NTCminusKO.csv"), row.names = FALSE)
options(width = 200); cat("\n=== Cm20, TEAD4 KO (z within data set) ===\n"); print(st, digits = 3, row.names = FALSE); cat(sprintf("bulk, 4 lines combined (lm z ~ line + TEAD4): NTC - TEAD4 = %.2f, p = %.3g; other KOs - TEAD4 = %.2f, p = %.3g\n", est[1], comb[1], est[2], comb[2]))
cat("\n-- NTC - KO of Cm20 per KO (rank 1 = most reduced by the KO)\n"); print(reshape(transform(E[, c("KO", "ds")], v = sprintf("%.2f (%g)", E$NTC_minus_KO, E$rank_most_down)), idvar = "KO", timevar = "ds", direction = "wide"), row.names = FALSE)

save_fig("fig14_Cm20_TEAD4", 13, 5.6, function() {
  layout(matrix(1:2, nrow = 1), widths = c(1.25, 1)); par(oma = c(0, 0, 2.6, 0.4), family = "sans", mar = c(4.4, 4.6, 3.6, 1)); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.85); mtext(b, side = 3, adj = 0, line = 0.5, cex = 0.62, col = INK2) }
  CC <- c(NTC = "grey45", TEAD4 = KO_COL[["TEAD4"]], other = "#d5d4cf"); yr <- c(-3, 4.3); plot(NA, xlim = c(0.4, 5.6), ylim = yr, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "Cm20 score (z within data set)", cex.axis = 0.8); abline(h = 0, lty = 2, col = GRID)
  for (i in seq_along(DS)) for (j in 1:3) { cl <- names(CC)[j]; v <- dat$z[dat$ds == DS[i] & dat$cls == cl]; at <- i + (j - 2) * 0.27
    if (length(v) >= 4) boxplot(v, at = at, add = TRUE, boxwex = 0.24, col = adjustcolor(CC[[cl]], 0.7), border = INK2, outline = FALSE, axes = FALSE)
    if (DS[i] != "CROP-seq") points(at + seq(-0.07, 0.07, length.out = length(v)), v, pch = 21, bg = CC[[cl]], col = INK2, cex = 0.9) else segments(at - 0.1, mean(v), at + 0.1, mean(v), lwd = 2, col = INK) }
  axis(1, at = 1:5, labels = c("CROP-seq\n(cells)", SH(LINES)), tick = FALSE, cex.axis = 0.78, padj = 0.5)
  text(1:5, 4.3, sprintf("NTC - TEAD4 %+.2f\np %s\nother - TEAD4 %+.2f\np %s", st$NTC_minus_TEAD4, ifelse(is.na(st$p_vs_NTC), "n/a (1 NTC)", fp(st$p_vs_NTC)), st$otherKO_minus_TEAD4, fp(st$p_vs_otherKO)), cex = 0.56, adj = c(0.5, 1))
  legend("bottomright", c("NTC", "TEAD4 KO", "other 7 KOs"), pt.bg = CC, pch = 22, pt.cex = 1.3, bty = "n", cex = 0.68)
  ttl("A  Cm20 in NTC, TEAD4 KO and the other KOs", sprintf("bulk, 4 lines combined: NTC - TEAD4 %+.2f SD, p = %s; other KOs - TEAD4 %+.2f SD, p = %s", est[1], fp(comb[1]), est[2], fp(comb[2])))
  plot(NA, xlim = c(0.5, 5.5), ylim = range(E$NTC_minus_KO) + c(-0.15, 0.3), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "Cm20, NTC - KO (z units; > 0 = lower in the KO)", cex.axis = 0.8); abline(h = 0, lty = 2, col = INK)
  for (i in seq_along(DS)) { e <- E[E$ds == DS[i], ]; points(rep(i, 8) + ifelse(e$KO == "TEAD4", 0, seq(-0.18, 0.18, length.out = 8)), e$NTC_minus_KO, pch = 21, bg = KO_COL[e$KO], col = INK2, cex = ifelse(e$KO == "TEAD4", 2, 1.1)); text(i, max(E$NTC_minus_KO) + 0.25, sprintf("TEAD4 #%g of 8", e$rank_most_down[e$KO == "TEAD4"]), cex = 0.62, font = 2, col = KO_COL[["TEAD4"]]) }
  axis(1, at = 1:5, labels = c("CROP-seq", SH(LINES)), tick = FALSE, cex.axis = 0.78); legend("bottomleft", KOS, pt.bg = KO_COL[KOS], pch = 21, bty = "n", cex = 0.55, ncol = 4)
  ttl("B  Cm20 effect of every KO", "large dot = TEAD4; # = rank of TEAD4 among the 8 KOs (1 = largest drop)")
  mtext("Cm20 in the TEAD4 KO: CROP-seq and bulk RNA-seq across the 4 cell lines (AddModuleScore, top-50 module genes; S70)", outer = TRUE, side = 3, adj = 0, line = 1.2, font = 2, cex = 0.95)
})
