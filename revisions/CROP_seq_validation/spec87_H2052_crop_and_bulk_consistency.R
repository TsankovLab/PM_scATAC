# Spec 87: per feature space (whole transcriptome | top 1000 / 2000 / 3000 variable genes | 20 Cm modules), two panels:
#   left  bar per TF = Spearman between the CROP-seq KO profile and the bulk NCI-H2052 KO profile (same KO, same cell line);
#         p = Spearman test over the features (genes are not independent -> nominal), label = rank of the own KO among the 8 bulk KOs (1 = best match)
#   right per TF: same KO in two bulk lines (6 line pairs, coloured box + dots) vs different KOs (grey), p = KO-label permutation within lines (BH stars)
# Profiles: NTC - KO t (genes) / module effect (Cm), rescaled = each feature z-scored across the 8 KOs within each data set (KO-specific part); no covariate.
# Variable genes ranked by mean within-line bulk variance (spec72 gene_var).
# needs: spec1_prepare_stats.R, spec72_build_profiles.R (S70), spec75 (spec75_screen_scores_by_cell.rds)
# outputs: spec87_crop_vs_H2052.csv, spec87_bulk_consistency.csv, fig9I_H2052_crop_and_bulk_consistency.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000")); perms <- make_perms(N_PERM)
ss <- readRDS(file.path(OUT, "spec_screen_stats.rds")); pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); CM <- names(cm_top50())
G <- intersect(pr$genes, ss$genes); vord <- G[order(rowMeans(pr$gene_var[G, ]), decreasing = TRUE)]
S_gene <- -ss$t0[G, KOS]; S_gene[!is.finite(S_gene)] <- 0                                   # screen Welch t is KO - NTC -> NTC - KO
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); grp <- as.character(sc$group); n0 <- apply(sc$score[grp == "NTC", CM], 2, median)
S_mod <- sapply(KOS, function(k) SGN * (n0 - apply(sc$score[grp == k, CM, drop = FALSE], 2, median)))
SPACES <- c(all = "whole transcriptome", top1000 = "top 1000 variable genes", top2000 = "top 2000 variable genes", top3000 = "top 3000 variable genes", Cm = "20 Cm modules")
z <- function(P) prep(P, TRUE, TRUE)
crop <- list(); cons <- list(); GL <- list(); CLs <- list()
for (sp in names(SPACES)) {
  if (sp == "Cm") { S <- S_mod[CM, ]; A <- pr$M$tmm$none$t[CM, , ] } else { f <- if (sp == "all") G else vord[1:as.integer(sub("top", "", sp))]; S <- S_gene[f, ]; A <- pr$P$tmm$none$t[f, , ] }
  Sz <- z(S); Bz <- z(A[, , "NCI-H2052"]); C <- cor(Sz, Bz, method = "spearman")
  crop[[sp]] <- data.frame(space = sp, KO = KOS, rho = diag(C), p = sapply(KOS, function(k) cor.test(Sz[, k], Bz[, k], method = "spearman", exact = FALSE)$p.value),
                           rank_own = sapply(seq_along(KOS), function(i) rank(-C[i, ])[i]), median_other = sapply(seq_along(KOS), function(i) median(C[i, -i])))
  CL <- cl_from_array(A, NULL, zscore = TRUE); a <- agreement(CL, perms); CLs[[sp]] <- a
  cons[[sp]] <- cbind(space = sp, a$per_KO); GL[[sp]] <- data.frame(space = sp, t(a$global)) }
crop <- do.call(rbind, crop); cons <- do.call(rbind, cons); GL <- do.call(rbind, GL)
write.csv(crop, file.path(OUT, "spec87_crop_vs_H2052.csv"), row.names = FALSE); write.csv(cons, file.path(OUT, "spec87_bulk_consistency.csv"), row.names = FALSE)
options(width = 220); cat("CROP-seq vs bulk H2052, own-KO rho (rank among 8):\n"); print(reshape(transform(crop[, c("KO", "space")], v = sprintf("%.2f (%d)", crop$rho, crop$rank_own)), idvar = "KO", timevar = "space", direction = "wide"), row.names = FALSE)
cat("\nbulk across lines, permutation p:\n"); print(reshape(transform(cons[, c("KO", "space")], v = sprintf("%.3f (r %.2f)", cons$perm_p, cons$mean_intra)), idvar = "KO", timevar = "space", direction = "wide"), row.names = FALSE); print(GL[, c("space", "diag_mean", "off_mean", "perm_p")], digits = 3, row.names = FALSE)

fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
save_fig("fig9I_H2052_crop_and_bulk_consistency", 14, 17, function() {
  layout(matrix(1:10, nrow = 5, byrow = TRUE), widths = c(0.85, 1.15)); par(oma = c(0, 0, 3.6, 0.4), mar = c(3.4, 4.6, 3, 1), family = "sans")
  for (sp in names(SPACES)) {
    d <- crop[crop$space == sp, ]; o <- d$KO[order(-d$rho)]; d <- d[match(o, d$KO), ]
    bp <- barplot(d$rho, col = adjustcolor(KO_COL[o], 0.7), border = NA, las = 1, ylim = c(min(-0.3, d$rho) - 0.05, max(0.7, d$rho) + 0.15), ylab = "Spearman, CROP-seq vs bulk H2052", names.arg = rep("", 8), cex.axis = 0.75)
    abline(h = 0, col = INK); points(bp, d$median_other, pch = "-", cex = 2.2, col = "grey30")
    text(bp, pmax(d$rho, 0) + 0.03, sprintf("p=%s\n#%d", fp(d$p), d$rank_own), cex = 0.55, adj = c(0.5, 0), font = ifelse(d$rank_own == 1, 2, 1))
    text(bp, par("usr")[3] - 0.03 * diff(par("usr")[3:4]), o, srt = 45, adj = 1, xpd = NA, cex = 0.75, font = 2, col = KO_COL[o])
    mtext(sprintf("CROP-seq vs bulk NCI-H2052 -- %s", SPACES[[sp]]), side = 3, adj = 0, line = 1.2, font = 2, cex = 0.78)
    mtext("bar = own KO; dash = median vs the other bulk KOs; # = rank of own KO among 8", side = 3, adj = 0, line = 0.3, cex = 0.5, col = INK2)
    a <- CLs[[sp]]; k <- a$per_KO; ob <- k$KO[order(k$perm_p, -k$mean_intra)]; yr <- range(c(unlist(a$cross), a$intra), na.rm = TRUE); yr <- yr + c(-0.03, 0.18) * diff(yr)
    plot(NA, xlim = c(0.5, 8.5), ylim = yr, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "Spearman between two bulk lines", cex.axis = 0.75); abline(h = 0, lty = 2, col = INK)
    for (i in seq_along(ob)) { j <- match(ob[i], k$KO); boxplot(a$cross[[j]], at = i - 0.2, add = TRUE, boxwex = 0.32, col = "#d5d4cf", border = INK2, outline = FALSE, axes = FALSE)
      boxplot(a$intra[ob[i], ], at = i + 0.2, add = TRUE, boxwex = 0.32, col = adjustcolor(KO_COL[[ob[i]]], 0.6), border = INK2, outline = FALSE, axes = FALSE)
      points(i + 0.2 + seq(-0.1, 0.1, length.out = 6), a$intra[ob[i], ], pch = 16, cex = 0.7, col = INK2)
      text(i, yr[2], sprintf("p=%s%s", fp(k$perm_p[j]), ifelse(k$perm_fdr[j] < 0.01, "\n**", ifelse(k$perm_fdr[j] < 0.05, "\n*", ""))), cex = 0.6, adj = c(0.5, 1), font = ifelse(k$perm_p[j] < 0.05, 2, 1)) }
    text(1:8, par("usr")[3] - 0.03 * diff(par("usr")[3:4]), ob, srt = 45, adj = 1, xpd = NA, cex = 0.75, font = 2, col = KO_COL[ob])
    g <- GL[GL$space == sp, ]; mtext(sprintf("bulk across 4 lines -- %s", SPACES[[sp]]), side = 3, adj = 0, line = 1.2, font = 2, cex = 0.78)
    mtext(sprintf("colour = same KO in two lines (dots = 6 pairs), grey = other KOs; global p = %s", fp(g$perm_p)), side = 3, adj = 0, line = 0.3, cex = 0.5, col = INK2) }
  mtext("Left: CROP-seq vs bulk in the same cell line (NCI-H2052). Right: consistency of each TF KO across the 4 bulk cell lines", outer = TRUE, side = 3, adj = 0, line = 2, font = 2, cex = 1)
  mtext("NTC - KO profiles rescaled (z across the 8 KOs per data set), no covariate, S70; left p = Spearman test over features (nominal: genes are correlated); right p = KO-label permutation (20,000), * / ** BH FDR < 0.05 / 0.01", outer = TRUE, side = 3, adj = 0, line = 0.8, cex = 0.6, col = INK2)
})
