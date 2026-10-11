# Spec 81: the CROP-seq "knock_downs_heatmap2" (git_repo/tumor_analysis/CRISPR_cropseq.R) reproduced for bulk RNA-seq, side by side with the screen.
# Value = median(NTC) - median(KO) of the Cm module score (AddModuleScore, top-50 genes) -- exactly the screen recipe; rows = KO, columns = the 20 Cm; palette
# rev(ggthemes::Orange-Blue-White Diverging), colour range capped symmetrically at min(|min|, |max|) per panel as in the original; stars = p < 0.05 / 0.01 / 0.001.
#   CROP-seq : cells, no cell-cycle regression (as the original figure), Wilcoxon NTC vs KO cells; the 8 KOs that are also in bulk (MEF2A / MEF2D not profiled in bulk)
#   bulk     : S70, module score on TMM logCPM (spec72), each line: median over NTC samples - median over KO samples, p = limma (score ~ group + % assigned; spec75, no cell cycle);
#              "4 lines" = mean of the per-line values, p = Stouffer over the lines (spec75)
# Row and column order: hierarchical clustering of the CROP-seq panel (ComplexHeatmap default), kept fixed in every panel so cells can be compared one to one.
# needs: spec72 (spec72_profiles_S70.rds), spec75 (spec75_screen_scores_by_cell.rds, spec75_Cm_KO_per_line_S70_noCC.csv, spec75_Cm_KO_table_S70_noCC.csv)
# outputs: spec81_Cm_KO_NTCminusKO_screen_and_bulk.csv, fig9C_Cm_KO_heatmap_cropseq_vs_bulk.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages({ library(ComplexHeatmap); library(grid) })
PAL <- rev(as.character(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 100)))
CM <- names(cm_top50()); meta <- load_full_audited_bulk()$meta
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); B <- pr$cm_scores$tmm[CM, pr$samples]
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); S <- sc$score[, CM]; grp <- as.character(sc$group)
pl <- read.csv(file.path(OUT, "spec75_Cm_KO_per_line_S70_noCC.csv")); tb <- read.csv(file.path(OUT, "spec75_Cm_KO_table_S70_noCC.csv"))

## screen: identical recipe to CRISPR_cropseq.R
Dm <- list(); Pm <- list()
Dm$screen <- t(sapply(ORDER, function(k) SGN * (apply(S[grp == "NTC", ], 2, median) - apply(S[grp == k, , drop = FALSE], 2, median))))
Pm$screen <- t(sapply(ORDER, function(k) sapply(CM, function(m) wilcox.test(S[grp == "NTC", m], S[grp == k, m], exact = FALSE)$p.value)))
## bulk per line
for (cl in LINES) { s <- pr$samples[meta[pr$samples, "Cell_Line"] == cl]; g <- meta[s, "Gene_Control"]
  Dm[[cl]] <- t(sapply(ORDER, function(k) SGN * (apply(B[, s[g == "NTC"], drop = FALSE], 1, median) - apply(B[, s[g == k], drop = FALSE], 1, median))))
  d <- pl[pl$line == cl, ]; Pm[[cl]] <- t(sapply(ORDER, function(k) d$p[d$KO == k][match(CM, d$Cm[d$KO == k])])) }
Dm$bulk4 <- Reduce(`+`, Dm[LINES]) / length(LINES); Pm$bulk4 <- t(sapply(ORDER, function(k) tb$bulk_p[tb$KO == k][match(CM, tb$Cm[tb$KO == k])]))
Dm$bulk4med <- apply(simplify2array(Dm[LINES]), c(1, 2), median); Pm$bulk4med <- Pm$bulk4                       # median of the 4 lines (same combined p)
for (n in names(Dm)) dimnames(Dm[[n]]) <- dimnames(Pm[[n]]) <- list(ORDER, CM)
long <- do.call(rbind, lapply(names(Dm), function(n) data.frame(source = n, KO = rep(ORDER, length(CM)), Cm = rep(CM, each = length(ORDER)), NTC_minus_KO = as.vector(Dm[[n]]), p = as.vector(Pm[[n]]))))
write.csv(long, file.path(OUT, "spec81_Cm_KO_NTCminusKO_screen_and_bulk.csv"), row.names = FALSE)

## agreement of each bulk panel with the screen over the 160 cells
agr <- sapply(setdiff(names(Dm), "screen"), function(n) { ct <- cor.test(as.vector(Dm$screen), as.vector(Dm[[n]]), method = "spearman", exact = FALSE); sg <- Pm$screen < 0.05
  c(rho = unname(ct$estimate), p = ct$p.value, sign = mean(sign(Dm$screen[sg]) == sign(Dm[[n]][sg]))) })
options(width = 200); cat("agreement with the CROP-seq panel (160 KO x Cm cells):\n"); print(round(agr, 3))

## fixed order from the CROP-seq clustering
ro <- hclust(dist(Dm$screen))$order; co <- hclust(dist(t(Dm$screen)))$order
TITLE <- c(screen = "CROP-seq (NCI-H2052 cells)", bulk4 = "bulk, mean of 4 lines", bulk4med = "bulk, median of 4 lines", setNames(paste("bulk", LINES), LINES))
hm <- function(n) { D <- Dm[[n]][ro, co]; P <- Pm[[n]][ro, co]; lim <- min(abs(range(D))); D <- pmax(pmin(D, lim), -lim)
  sub <- if (n == "screen") sprintf("%d-%d cells per KO", min(table(grp)[ORDER]), max(table(grp)[ORDER])) else sprintf("vs screen: rho %.2f, same sign %.0f%%", agr["rho", n], 100 * agr["sign", n])
  Heatmap(D, name = n, col = PAL, cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize = 7, col = KO_COL[rownames(D)], fontface = "bold"),
          column_names_rot = 45, column_names_gp = gpar(fontsize = 6.5), column_title = paste0(TITLE[[n]], "\n", sub), column_title_gp = gpar(fontsize = 7.5, fontface = "bold"),
          heatmap_legend_param = list(title = "NTC - KO", title_gp = gpar(fontsize = 6.5), labels_gp = gpar(fontsize = 6), grid_width = unit(2.5, "mm")),
          cell_fun = function(j, i, x, y, w, h, fill) { p <- P[i, j]; if (is.finite(p) && p < 0.05) grid.text(if (p < 0.001) "***" else if (p < 0.01) "**" else "*", x, y, vjust = 0.8, gp = gpar(fontsize = 5)) }) }
save_fig("fig9C_Cm_KO_heatmap_cropseq_vs_bulk", 15, 6.2, function() {
  grid.newpage(); pushViewport(viewport(layout = grid.layout(3, 3, heights = unit(c(1.1, 1, 1), c("cm", "null", "null")))))
  grid.text("Effect of each KO on the 20 Cm modules, NTC - KO (orange = module lower after KO): CROP-seq recipe applied to bulk RNA-seq", x = unit(4, "mm"), y = unit(0.62, "npc"), just = "left", gp = gpar(fontsize = 11, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  grid.text("median(NTC) - median(KO) of the module score; same row / column order in every panel (clustering of the CROP-seq panel); colour capped per panel at min(|min|, |max|) as in the original; stars: CROP-seq Wilcoxon, bulk limma per line / Stouffer over lines", x = unit(4, "mm"), y = unit(0.2, "npc"), just = "left", gp = gpar(fontsize = 7, col = INK2), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  pos <- list(screen = c(2, 1), bulk4 = c(3, 1), "NCI-H2052" = c(2, 2), "MSTO-211H" = c(2, 3), "NCI-H2452" = c(3, 2), "NCI-H28" = c(3, 3))
  for (n in names(pos)) { pushViewport(viewport(layout.pos.row = pos[[n]][1], layout.pos.col = pos[[n]][2])); draw(hm(n), newpage = FALSE, padding = unit(c(1, 2, 1, 2), "mm")); popViewport() }
})

## fig9C2: CROP-seq vs bulk mean vs bulk median of the 4 lines, same order and recipe
save_fig("fig9C2_Cm_KO_heatmap_cropseq_vs_bulk_mean_median", 15, 3.6, function() {
  grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 3, heights = unit(c(0.9, 1), c("cm", "null")))))
  grid.text("Effect of each KO on the 20 Cm modules (NTC - KO, orange = module lower after KO): CROP-seq vs bulk summarised over the 4 lines as mean or median", x = unit(4, "mm"), y = unit(0.6, "npc"), just = "left", gp = gpar(fontsize = 10, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  grid.text("same row / column order (clustering of the CROP-seq panel); colour capped per panel at min(|min|, |max|); stars: CROP-seq Wilcoxon, bulk Stouffer over lines (same for mean and median)", x = unit(4, "mm"), y = unit(0.15, "npc"), just = "left", gp = gpar(fontsize = 6.5, col = INK2), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  for (j in 1:3) { n <- c("screen", "bulk4", "bulk4med")[j]; pushViewport(viewport(layout.pos.row = 2, layout.pos.col = j)); draw(hm(n), newpage = FALSE, padding = unit(c(1, 2, 1, 2), "mm")); popViewport() }
})
