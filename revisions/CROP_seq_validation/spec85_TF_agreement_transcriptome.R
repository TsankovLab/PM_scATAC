# Spec 85: per-TF agreement at the TRANSCRIPTOME level (genes, not Cm modules), with rescaling, on the whole transcriptome and on top features.
#   profiles : NTC - KO t statistic per gene; screen = Welch t of cells (spec1; no / cell cycle regressed), bulk = limma moderated t per line (spec72 S70, TMM; none / cc)
#   rescaling: each gene z-scored across the 8 KOs within each data set (screen, each bulk line) -> only KO-specific differences are compared
#   features : all expressed genes shared by screen and bulk | top 2000 variable genes (mean within-line bulk variance) |
#              top 500 KO-signal genes (screen vs bulk: rank sum of screen max|t| and bulk mean max|t|; bulk vs bulk: signal in the 2 lines NOT compared, as spec76)
#   A bulk vs CROP-seq: per TF, mean same-KO Spearman over the 4 (screen, line) pairs; p = KO-label permutation within the screen and within every line
#   B bulk vs bulk     : per TF, mean same-KO Spearman over the 6 line pairs; same permutation
# needs: spec1_prepare_stats.R (spec_screen_stats.rds), spec72_build_profiles.R (S70), spec82 (TF order)
# outputs: spec85_TF_agreement_transcriptome.csv, spec85_global.csv, fig9G_TF_agreement_transcriptome.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000"))
ss <- readRDS(file.path(OUT, "spec_screen_stats.rds")); pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF)))
G <- intersect(pr$genes, ss$genes); cat("genes shared by screen and bulk:", length(G), "\n")
o <- read.csv(file.path(OUT, "spec82_TF_agreement_ranked.csv"))$KO
COVL <- c(none = "no covariate", cc = "cell cycle regressed")
SCR <- list(none = -ss$t0[G, KOS], cc = -ss$tcc[G, KOS])                 # screen t is KO - NTC -> NTC - KO
BLK <- list(none = pr$P$tmm$none$t[G, , ], cc = pr$P$tmm$cc$t[G, , ])
for (cv in names(SCR)) SCR[[cv]][!is.finite(SCR[[cv]])] <- 0
vord <- G[order(rowMeans(pr$gene_var[G, ]), decreasing = TRUE)]
sig_scr <- apply(abs(SCR$none), 1, max); sig_blk <- rowMeans(pr$gene_signal[G, ])
top_sig <- G[order(rank(-sig_scr) + rank(-sig_blk))][1:500]
SP <- c(all = "whole transcriptome", top2000 = "top 2000 variable genes", signal500 = "top 500 KO-signal genes")
zs <- function(P) prep(P, TRUE, TRUE)
DS <- c("screen", LINES); agree5 <- agreement; environment(agree5) <- list2env(list(LINES = DS), parent = globalenv())   # permutations over 5 data sets
p5 <- make_perms(N_PERM, L = 5, seed = 2); p4 <- make_perms(N_PERM)

res <- list(); glob <- list()
for (cv in names(COVL)) for (sp in names(SP)) {
  f <- switch(sp, all = G, top2000 = vord[1:2000], signal500 = top_sig)
  ## A screen vs each bulk line (4 pairs)
  CLa <- lapply(LINES, function(cl) { C <- pair_cor(zs(SCR[[cv]][f, ]), zs(BLK[[cv]][f, , cl]), "spearman"); dimnames(C) <- list(KOS, KOS); C })
  a <- agree5(CLa, p5, pairs = lapply(LINES, function(cl) c("screen", cl)))
  ## B bulk vs bulk (6 pairs)
  CLb <- if (sp == "signal500") cl_from_array(BLK[[cv]], "signal", pr$gene_signal[G, ], zscore = TRUE) else cl_from_array(BLK[[cv]], f, zscore = TRUE)
  b <- agreement(CLb, p4)
  for (ax in c("A", "B")) { r <- if (ax == "A") a else b
    res[[length(res) + 1]] <- data.frame(axis = ax, cov = cv, space = sp, r$per_KO[, c("KO", "mean_intra", "median_cross", "n_pos_pairs", "perm_p", "perm_fdr", "perm_z")])
    glob[[length(glob) + 1]] <- data.frame(axis = ax, cov = cv, space = sp, t(r$global)) } }
res <- do.call(rbind, res); glob <- do.call(rbind, glob)
write.csv(res, file.path(OUT, "spec85_TF_agreement_transcriptome.csv"), row.names = FALSE); write.csv(glob, file.path(OUT, "spec85_global.csv"), row.names = FALSE)
options(width = 230); cat("\n=== global (perm p) ===\n"); print(glob[, c("axis", "cov", "space", "diag_mean", "off_mean", "perm_p", "perm_z")], digits = 3, row.names = FALSE)
for (ax in c("A", "B")) { cat(sprintf("\n=== %s per TF: perm p ===\n", c(A = "bulk vs CROP-seq", B = "bulk vs bulk")[ax])); x <- res[res$axis == ax, ]
  print(reshape(transform(x[, c("KO", "cov", "space")], v = sprintf("%.3f (r %.2f)", x$perm_p, x$mean_intra), id = paste(x$cov, x$space))[, c("KO", "id", "v")], idvar = "KO", timevar = "id", direction = "wide")[match(o, unique(x$KO)), ], row.names = FALSE) }

## ---------------------------------------------------------------- figure: rows = feature spaces, columns = A / B; 2 bars per TF (no covariate, cell cycle)
fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2)); COLS <- list(A = c("grey50", "#eb6834"), B = c("#2a78d6", "#1baf7a"))
save_fig("fig9G_TF_agreement_transcriptome", 15, 11.5, function() {
  layout(matrix(1:6, nrow = 3, byrow = TRUE)); par(oma = c(0.5, 0, 3.4, 0.4), mar = c(3.2, 4.6, 2.6, 1), family = "sans"); n <- length(o)
  yl <- range(c(res$mean_intra, 0)); yl <- yl + c(-0.03, 0.25) * diff(yl)
  for (sp in names(SP)) for (ax in c("A", "B")) { x <- res[res$axis == ax & res$space == sp, ]; cols <- COLS[[ax]]; at <- outer(c(-0.19, 0.19), seq_len(n), "+")
    plot(NA, xlim = c(0.5, n + 0.5), ylim = yl, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = if (ax == "A") "same-KO Spearman, screen vs bulk line" else "same-KO Spearman, bulk line vs bulk line", cex.axis = 0.75, cex.lab = 0.8)
    abline(h = 0, col = INK, lwd = 0.8)
    for (j in 1:2) { cv <- names(COVL)[j]; for (i in seq_len(n)) { r <- x[x$cov == cv & x$KO == o[i], ]; xx <- at[j, i]
      rect(xx - 0.17, 0, xx + 0.17, r$mean_intra, col = adjustcolor(cols[j], 0.6), border = cols[j], lwd = 0.8); segments(xx - 0.17, r$median_cross, xx + 0.17, r$median_cross, lty = 2, col = "grey35")
      text(xx, yl[2] - (j - 1) * 0.08 * diff(yl), paste0(fp(r$perm_p), ifelse(r$perm_fdr < 0.01, " **", ifelse(r$perm_fdr < 0.05, " *", ""))), cex = 0.58, col = cols[j], font = ifelse(r$perm_p < 0.05, 2, 1)) } }
    text(seq_len(n), yl[1] - 0.05 * diff(yl), o, xpd = NA, cex = 0.8, font = 2, col = KO_COL[o], adj = c(0.5, 1))
    g <- glob[glob$axis == ax & glob$space == sp, ]
    mtext(sprintf("%s: %s", c(A = "bulk vs CROP-seq", B = "bulk vs bulk")[ax], SP[[sp]]), side = 3, adj = 0, line = 1.1, font = 2, cex = 0.82)
    mtext(sprintf("global p: %s", paste(sprintf("%s %s", COVL[g$cov], fp(g$perm_p)), collapse = " | ")), side = 3, adj = 0, line = 0.25, cex = 0.55, col = INK2) }
  mtext("Agreement of each TF's KO effect on the transcriptome: with CROP-seq (left) and between bulk cell lines (right); TFs ordered by cross-line agreement", outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 1)
  mtext("NTC - KO t per gene, z-scored across the 8 KOs within every data set (rescaled); bar = mean same-KO Spearman (left: 4 screen-line pairs, right: 6 line pairs); left bars grey/orange and right bars blue/green = no covariate / cell cycle regressed; dashed = median vs other KOs; numbers = KO-label permutation p, bold < 0.05, * / ** BH FDR < 0.05 / 0.01", outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.6, col = INK2)
})
