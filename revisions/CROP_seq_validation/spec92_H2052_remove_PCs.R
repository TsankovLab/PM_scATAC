# Spec 92: CROP-seq vs bulk restricted to NCI-H2052 (the screen's own cell line) -- effect of removing the top k KO axes (PCs of the KO-centred NTC - KO profiles, k = 0..3)
# from bulk only, or from bulk AND screen (spec91 outputs; no recomputation).
#   A per feature space: mean own-KO minus mean other-KO Spearman (exact permutation p of the bulk KO labels, 8!)
#   B per TF: own-KO Spearman, coloured; number = rank of the own KO among the 8 bulk KOs (1 = best match)
# needs: spec91_remove_bulk_PC_axes.R run twice (STRIP_SCREEN=0 and STRIP_SCREEN=1)
# outputs: spec92_H2052_remove_PCs.csv, fig9N_H2052_remove_PCs.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(grid) })
MODES <- c(bulk = "axes removed from bulk only", both = "axes removed from bulk and screen"); SPL <- c(hvg = "screen top-5000 HVGs", top2000 = "top 2000 bulk variable genes", all = "all shared genes")
rd <- function(f) read.csv(file.path(OUT, f))
VS <- rbind(cbind(mode = "bulk", rd("spec91_vs_screen.csv")), cbind(mode = "both", rd("spec91_vs_screen_both.csv"))); VS <- VS[VS$line == "NCI-H2052", ]; VS$dmo <- VS$diag - VS$off
VK <- rbind(cbind(mode = "bulk", rd("spec91_vs_screen_per_KO.csv")), cbind(mode = "both", rd("spec91_vs_screen_per_KO_both.csv"))); VK <- VK[VK$line == "NCI-H2052", ]
write.csv(merge(VK, VS[, c("mode", "space", "k", "dmo", "perm_p")], by = c("mode", "space", "k")), file.path(OUT, "spec92_H2052_remove_PCs.csv"), row.names = FALSE)
options(width = 200); print(VS[, c("mode", "space", "k", "diag", "off", "dmo", "perm_p", "top1")], digits = 3, row.names = FALSE)

save_fig("fig9N_H2052_remove_PCs", 15, 8.6, function() {
  grid.newpage(); pushViewport(viewport(layout = grid.layout(3, 3, heights = unit(c(1.2, 1, 1.15), c("cm", "null", "null")))))
  grid.text("CROP-seq vs bulk, NCI-H2052 only (same cell line): does removing the dominant KO axes improve the match?", x = unit(4, "mm"), y = unit(0.65, "npc"), just = "left", gp = gpar(fontsize = 11, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  grid.text("NTC - KO profiles centred over the 8 KOs, top k PCs projected out (PC1 = proliferation vs inflammation axis); top: own-KO minus other-KO Spearman, exact permutation p; bottom: own-KO rho per TF, number = rank of own KO among 8", x = unit(4, "mm"), y = unit(0.2, "npc"), just = "left", gp = gpar(fontsize = 7, col = INK2), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  for (j in seq_along(SPL)) { sp <- names(SPL)[j]; pushViewport(viewport(layout.pos.row = 2, layout.pos.col = j)); pushViewport(plotViewport(c(3.2, 4.2, 2.2, 1))); d <- VS[VS$space == sp, ]
    yr <- range(c(d$dmo, 0)) + c(-0.02, 0.03); pushViewport(dataViewport(xscale = c(-0.3, 3.3), yscale = yr)); grid.xaxis(at = 0:3, gp = gpar(fontsize = 7)); grid.yaxis(gp = gpar(fontsize = 7))
    grid.lines(unit(c(-0.3, 3.3), "native"), unit(c(0, 0), "native"), gp = gpar(lty = 2, col = INK2)); mc <- c(bulk = "#2a78d6", both = "#eb6834")
    for (m in names(MODES)) { e <- d[d$mode == m, ]; e <- e[order(e$k), ]; grid.lines(unit(e$k, "native"), unit(e$dmo, "native"), gp = gpar(col = mc[m], lwd = 2))
      grid.points(unit(e$k, "native"), unit(e$dmo, "native"), pch = ifelse(e$perm_p < 0.05, 16, 1), size = unit(2.6, "mm"), gp = gpar(col = mc[m]))
      grid.text(sprintf("p %.2f", e$perm_p), unit(e$k, "native"), unit(e$dmo, "native") + unit(ifelse(m == "bulk", 3, -3), "mm"), gp = gpar(fontsize = 5.5, col = mc[m])) }
    grid.text("top k KO axes removed", y = unit(-2.4, "lines"), gp = gpar(fontsize = 7)); grid.text("own-KO - other-KO rho", x = unit(-3, "lines"), rot = 90, gp = gpar(fontsize = 7))
    grid.text(SPL[[sp]], y = unit(1, "npc") + unit(0.8, "lines"), gp = gpar(fontsize = 8, fontface = "bold"))
    if (j == 1) grid.legend(MODES, pch = 16, gp = gpar(col = mc, fontsize = 6.5), vp = viewport(x = unit(0.7, "npc"), y = unit(0.88, "npc")))
    popViewport(3) }
  for (j in seq_along(MODES)) { m <- names(MODES)[j]; x <- VK[VK$mode == m & VK$space == "hvg", ]
    M <- sapply(0:3, function(k) x$own[x$k == k][match(KOS, x$KO[x$k == k])]); R <- sapply(0:3, function(k) x$rank[x$k == k][match(KOS, x$KO[x$k == k])]); dimnames(M) <- dimnames(R) <- list(KOS, paste(0:3, "removed"))
    h <- Heatmap(M, name = m, col = colorRamp2(c(-0.35, 0, 0.35), c("#2a78d6", "#f0efec", "#eb6834")), cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = j == 2,
      heatmap_legend_param = list(title = "own-KO rho", title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6)), row_names_side = "left", row_names_gp = gpar(fontsize = 8, col = KO_COL[KOS], fontface = "bold"),
      column_names_gp = gpar(fontsize = 7.5), column_names_rot = 0, column_names_centered = TRUE, rect_gp = gpar(col = "white", lwd = 0.8),
      column_title = sprintf("per TF, screen HVGs -- %s", MODES[[m]]), column_title_gp = gpar(fontsize = 8, fontface = "bold"),
      cell_fun = function(jj, ii, xx, yy, w, hh, fill) grid.text(sprintf("%.2f\n#%d", M[ii, jj], R[ii, jj]), xx, yy, gp = gpar(fontsize = 6, fontface = ifelse(R[ii, jj] == 1, "bold", "plain"))))
    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = j)); draw(h, newpage = FALSE, padding = unit(c(2, 4, 2, 4), "mm")); popViewport() }
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 3)); pushViewport(plotViewport(c(3, 1, 2.5, 1)))
  grid.text(paste("How to read", "", "#1 = the screen KO matches its own bulk KO", "better than any of the 7 other bulk KOs.", "", "Without removal: PITX1 #2, TEAD4 #3.", "After removing the axes: PITX1 / TEAD4", "drop, BPTF becomes #1.", sep = "\n"), x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 7.5, col = INK2))
  popViewport(2)
})
