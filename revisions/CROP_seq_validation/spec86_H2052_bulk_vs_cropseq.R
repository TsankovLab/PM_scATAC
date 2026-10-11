# Spec 86: NCI-H2052 only -- bulk (S70) vs CROP-seq (same cell line), KO x KO Spearman matrices on: whole transcriptome | top 1000 / 2000 / 3000 variable genes | 20 Cm modules.
#   genes  : NTC - KO t per gene (screen Welch t on cells; bulk limma moderated t, TMM) ; modules: screen median(NTC) - median(KO) of the cell module score vs bulk module t
#   prep   : rescaled = each feature z-scored across the 8 KOs within each data set (KO-specific part only) | uncentred = NTC - KO as is (includes the response shared by all KOs)
#   variable genes ranked by bulk H2052 variance (TMM logCPM, spec72 gene_var)
#   global p = exact permutation of the bulk KO labels (all 8! = 40320), statistic = mean diagonal; per TF = rank of the own KO in its row (1 = best of 8) -- one pair only, so per-TF
#   significance cannot be finer than 1/8
# Shown without covariate and with cell cycle regressed (screen: cell cc score; bulk: cc covariate).
# needs: spec1_prepare_stats.R, spec72_build_profiles.R (S70), spec75 (spec75_screen_scores_by_cell.rds)
# outputs: spec86_H2052_global.csv, spec86_H2052_per_KO.csv, fig9H_H2052_bulk_vs_cropseq.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(grid) })
CL <- "NCI-H2052"; ss <- readRDS(file.path(OUT, "spec_screen_stats.rds")); pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); CM <- names(cm_top50())
G <- intersect(pr$genes, ss$genes); vord <- G[order(pr$gene_var[G, CL], decreasing = TRUE)]
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); grp <- as.character(sc$group)
scr_mod <- function(cc) { S <- sc$score[, CM]; if (cc) S <- apply(S, 2, function(y) resid(lm(y ~ sc$cc)) + mean(y)); n <- apply(S[grp == "NTC", ], 2, median); sapply(KOS, function(k) SGN * (n - apply(S[grp == k, , drop = FALSE], 2, median))) }
SCR <- list(none = -ss$t0[G, KOS], cc = -ss$tcc[G, KOS]); for (cv in names(SCR)) SCR[[cv]][!is.finite(SCR[[cv]])] <- 0
SPACES <- c(all = "whole transcriptome", top1000 = "top 1000 variable genes", top2000 = "top 2000 variable genes", top3000 = "top 3000 variable genes", Cm = "20 Cm modules")
PREP <- c(rescaled = "rescaled (z across KOs)", uncentred = "uncentred (NTC - KO)"); COVL <- c(none = "no covariate", cc = "cell cycle regressed")
perms <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(1:8)) })
mats <- list(); glob <- list(); pk <- list()
for (cv in names(COVL)) for (sp in names(SPACES)) for (pp in names(PREP)) {
  if (sp == "Cm") { S <- scr_mod(cv == "cc")[CM, ]; B <- pr$M$tmm[[cv]]$t[CM, , CL] } else { f <- if (sp == "all") G else vord[1:as.integer(sub("top", "", sp))]; S <- SCR[[cv]][f, ]; B <- pr$P$tmm[[cv]]$t[f, , CL] }
  if (pp == "rescaled") { S <- prep(S, TRUE, TRUE); B <- prep(B, TRUE, TRUE) }
  C <- pair_cor(S, B, "spearman"); dimnames(C) <- list(KOS, KOS); mats[[paste(cv, sp, pp)]] <- C
  dm <- mean(diag(C)); nul <- apply(perms, 1, function(p) mean(C[cbind(1:8, p)]))
  glob[[length(glob) + 1]] <- data.frame(cov = cv, space = sp, prep = pp, n_features = nrow(S), diag_mean = dm, off_mean = mean(C[row(C) != col(C)]), perm_p = mean(nul >= dm - 1e-12), top1 = sum(apply(C, 1, which.max) == 1:8))
  pk[[length(pk) + 1]] <- data.frame(cov = cv, space = sp, prep = pp, KO = KOS, own = diag(C), median_other = sapply(1:8, function(i) median(C[i, -i])), rank_in_row = sapply(1:8, function(i) rank(-C[i, ])[i])) }
glob <- do.call(rbind, glob); pk <- do.call(rbind, pk)
write.csv(glob, file.path(OUT, "spec86_H2052_global.csv"), row.names = FALSE); write.csv(pk, file.path(OUT, "spec86_H2052_per_KO.csv"), row.names = FALSE)
options(width = 220); print(glob, digits = 3, row.names = FALSE)
cat("\nper TF own-KO rho (rank in row), rescaled:\n"); for (cv in names(COVL)) { x <- pk[pk$prep == "rescaled" & pk$cov == cv, ]; cat(" ", cv, "\n"); print(reshape(transform(x[, c("KO", "space")], v = sprintf("%.2f (%d)", x$own, x$rank_in_row)), idvar = "KO", timevar = "space", direction = "wide"), row.names = FALSE) }

## ---------------------------------------------------------------- figure: row 1 = rescaled heatmaps (no covariate) per space; row 2 = per-TF own rho, 4 variants, per space
PCOL <- c("none rescaled" = "#2a78d6", "cc rescaled" = "#1baf7a", "none uncentred" = "grey65", "cc uncentred" = "grey35")
save_fig("fig9H_H2052_bulk_vs_cropseq", 17, 8.6, function() {
  grid.newpage(); pushViewport(viewport(layout = grid.layout(3, 5, heights = unit(c(1.2, 1, 1.05), c("cm", "null", "null")))))
  grid.text("NCI-H2052: bulk RNA-seq vs CROP-seq (same cell line), KO x KO Spearman of NTC - KO profiles", x = unit(4, "mm"), y = unit(0.65, "npc"), just = "left", gp = gpar(fontsize = 11, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:5))
  grid.text("top: rescaled profiles (each feature z-scored across the 8 KOs), no covariate; rows = CROP-seq KO, columns = bulk KO, diagonal boxed. bottom: own-KO rho per TF, 4 variants; global p = exact permutation of bulk KO labels (8!)", x = unit(4, "mm"), y = unit(0.2, "npc"), just = "left", gp = gpar(fontsize = 7, col = INK2), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:5))
  for (j in seq_along(SPACES)) { sp <- names(SPACES)[j]; C <- mats[[paste("none", sp, "rescaled")]]; g <- glob[glob$cov == "none" & glob$space == sp & glob$prep == "rescaled", ]; lim <- max(0.25, max(abs(C)))
    h <- Heatmap(C, name = sp, col = colorRamp2(c(-lim, 0, lim), c("#2a78d6", "#f0efec", "#eb6834")), cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = j == 5,
      heatmap_legend_param = list(title = "rho", title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6)), row_names_gp = gpar(fontsize = 7, col = KO_COL[KOS], fontface = "bold"), row_names_side = "left",
      column_names_gp = gpar(fontsize = 7, col = KO_COL[KOS], fontface = "bold"), column_names_rot = 45, rect_gp = gpar(col = "white", lwd = 0.5),
      column_title = sprintf("%s (%d)\ndiag %.3f vs off %.3f, p = %s", SPACES[[sp]], g$n_features, g$diag_mean, g$off_mean, formatC(g$perm_p, format = "g", digits = 2)), column_title_gp = gpar(fontsize = 7.5, fontface = "bold"),
      cell_fun = function(jj, ii, x, y, w, hh, fill) { grid.text(sprintf("%.2f", C[ii, jj]), x, y, gp = gpar(fontsize = 5)); if (ii == jj) grid.rect(x, y, w, hh, gp = gpar(fill = NA, col = INK, lwd = 1.4)) })
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = j)); draw(h, newpage = FALSE, padding = unit(c(1, 1, 1, 1), "mm")); popViewport()
    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = j)); pushViewport(plotViewport(c(3.2, 3.4, 2.2, 0.6))); x <- pk[pk$space == sp, ]
    yr <- c(-0.6, 0.85); pushViewport(dataViewport(xscale = c(0.5, 8.5), yscale = yr)); grid.yaxis(gp = gpar(fontsize = 6)); grid.lines(unit(c(0.5, 8.5), "native"), unit(c(0, 0), "native"), gp = gpar(col = INK2, lty = 2))
    off <- c(-0.3, -0.1, 0.1, 0.3); for (v in seq_along(PCOL)) { pv <- strsplit(names(PCOL)[v], " ")[[1]]; d <- x[x$cov == pv[1] & x$prep == pv[2], ]
      grid.rect(unit(1:8 + off[v], "native"), unit(0, "native"), unit(0.18, "native"), unit(d$own, "native"), just = c("centre", "bottom"), gp = gpar(fill = PCOL[v], col = NA))
      if (pv[2] == "rescaled") grid.text(ifelse(d$rank_in_row == 1, "1st", ""), unit(1:8 + off[v], "native"), unit(pmax(d$own, 0) + 0.04, "native"), gp = gpar(fontsize = 5, col = PCOL[v], fontface = "bold")) }
    grid.text(KOS, unit(1:8, "native"), unit(-1.1, "lines"), rot = 45, just = "right", gp = gpar(fontsize = 6.5, col = KO_COL[KOS], fontface = "bold"))
    gg <- glob[glob$space == sp, ]; grid.text(paste(sprintf("%s: p %s", c("none rescaled" = "rescaled", "cc rescaled" = "rescaled+cc", "none uncentred" = "uncentred", "cc uncentred" = "uncentred+cc")[paste(gg$cov, gg$prep)], formatC(gg$perm_p, format = "g", digits = 2)), collapse = "; "),
      x = unit(0, "npc"), y = unit(1, "npc") + unit(0.6, "lines"), just = "left", gp = gpar(fontsize = 5.6, col = INK2))
    if (j == 1) { grid.text("own-KO Spearman (screen vs H2052 bulk)", x = unit(-2.6, "lines"), rot = 90, gp = gpar(fontsize = 6.5))
      grid.legend(c("rescaled", "rescaled + cc", "uncentred", "uncentred + cc"), pch = 15, gp = gpar(col = PCOL, fontsize = 6), vp = viewport(x = unit(0.78, "npc"), y = unit(0.82, "npc"), width = unit(0.3, "npc"))) }
    popViewport(3) }
})
