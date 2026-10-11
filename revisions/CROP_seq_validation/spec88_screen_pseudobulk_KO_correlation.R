# Spec 88: CROP-seq screen pseudobulk per KO group (NTC + every KO incl. MEF2A / MEF2D): counts summed over the cells of the group -> logCPM (edgeR, prior 2);
# top 5000 highly variable genes of the cells (Seurat vst) that have CPM >= 10 in at least one group. KO x KO correlation heatmaps:
#   A raw logCPM, Spearman (dominated by the shared baseline, all ~1)
#   B gene-centred (each gene's mean over the groups subtracted, fig13c / spec62 convention), Spearman
#   C NTC - KO log fold change (KOs only), centred over KOs, Spearman -- same view as the bulk KO profiles
# Hierarchical clustering (complete, 1 - rho) on the centred matrix; numbers in cells; cell counts as side bars.
# needs: spec_common.R (load_screen())
# outputs: spec88_screen_pseudobulk_logCPM.csv, spec88_screen_pseudobulk_cor.csv, fig9J_screen_pseudobulk_KO_correlation.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages({ library(Seurat); library(edgeR); library(ComplexHeatmap); library(circlize); library(grid) })
TOPN <- as.integer(Sys.getenv("TOPN", unset = "5000"))
srt <- load_screen(); X <- GetAssayData(srt, layer = "counts"); grp <- as.character(srt$merged_call)
GR <- c("NTC", ORDER, setdiff(sort(unique(grp)), c("NTC", ORDER))); ncell <- table(grp)[GR]
PB <- sapply(GR, function(g) Matrix::rowSums(X[, grp == g, drop = FALSE]))
lc <- cpm(DGEList(PB), log = TRUE, prior.count = 2); keep <- rowSums(cpm(DGEList(PB)) >= 10) >= 1; lc <- lc[keep, ]
# top variable genes from the CELLS (Seurat vst), not from the pseudobulks: variance across pseudobulks is dominated by drop-outs of the small groups (HMGA1 13, TEAD2 20 cells)
hv <- VariableFeatures(FindVariableFeatures(srt, selection.method = "vst", nfeatures = TOPN, verbose = FALSE)); top <- intersect(hv, rownames(lc)); L <- lc[top, ]
write.csv(round(L, 4), file.path(OUT, "spec88_screen_pseudobulk_logCPM.csv"))
cat(sprintf("groups: %d, cells per group %d-%d, expressed genes %d, top %d used\n", length(GR), min(ncell), max(ncell), nrow(lc), length(top)))
C_raw <- cor(L, method = "spearman"); Lc <- L - rowMeans(L); C_cen <- cor(Lc, method = "spearman")
KOg <- setdiff(GR, "NTC"); FC <- -(L[, KOg] - L[, "NTC"]); FC <- FC - rowMeans(FC); C_fc <- cor(FC, method = "spearman")   # NTC - KO, centred over KOs
write.csv(rbind(cbind(type = "raw", as.data.frame(round(C_raw, 4))), cbind(type = "gene_centred", as.data.frame(round(C_cen, 4)))), file.path(OUT, "spec88_screen_pseudobulk_cor.csv"))
options(width = 200); cat("\ngene-centred Spearman:\n"); print(round(C_cen, 2)); cat("\nNTC - KO (centred over KOs) Spearman:\n"); print(round(C_fc, 2))

GCOL <- c(NTC = "grey45", KO_COL[ORDER], MEF2A = "#8c6d31", MEF2D = "#bd9e39"); GCOL <- GCOL[GR]
hm <- function(C, ttl, lim, clus, ann = TRUE, legend_title = "Spearman") { n <- rownames(C)
  Heatmap(C, name = ttl, col = colorRamp2(c(-lim, 0, lim), c("#2a78d6", "#f0efec", "#eb6834")), cluster_rows = clus, cluster_columns = clus, show_heatmap_legend = TRUE,
    heatmap_legend_param = list(title = legend_title, title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6)), row_names_gp = gpar(fontsize = 8, col = GCOL[n], fontface = "bold"),
    column_names_gp = gpar(fontsize = 8, col = GCOL[n], fontface = "bold"), column_names_rot = 45, rect_gp = gpar(col = "white", lwd = 0.6), row_dend_width = unit(6, "mm"), column_dend_height = unit(6, "mm"),
    right_annotation = if (ann) rowAnnotation(cells = anno_barplot(as.numeric(ncell[n]), gp = gpar(fill = GCOL[n], col = NA), width = unit(1.4, "cm"), axis_param = list(gp = gpar(fontsize = 5.5))), annotation_name_gp = gpar(fontsize = 6.5)) else NULL,
    column_title = ttl, column_title_gp = gpar(fontsize = 8.5, fontface = "bold"), cell_fun = function(j, i, x, y, w, h, fill) grid.text(sprintf("%.2f", C[i, j]), x, y, gp = gpar(fontsize = 5.2))) }
clus <- hclust(as.dist(1 - C_cen), method = "complete")
save_fig("fig9J_screen_pseudobulk_KO_correlation", 17, 6.4, function() {
  grid.newpage(); pushViewport(viewport(layout = grid.layout(2, 3, heights = unit(c(1.2, 1), c("cm", "null")), widths = unit(c(1, 1.1, 1), "null"))))
  grid.text(sprintf("CROP-seq pseudobulk per KO: correlation over the top %d variable genes (NCI-H2052 cells)", length(top)), x = unit(4, "mm"), y = unit(0.65, "npc"), just = "left", gp = gpar(fontsize = 11, fontface = "bold"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  grid.text("pseudobulk = counts summed over the cells of each guide group, logCPM; A raw; B each gene centred over the groups; C NTC - KO log fold change centred over KOs (the view used for bulk KO profiles); right bars = cells per group", x = unit(4, "mm"), y = unit(0.2, "npc"), just = "left", gp = gpar(fontsize = 7, col = INK2), vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  dr <- function(h, j) { pushViewport(viewport(layout.pos.row = 2, layout.pos.col = j)); draw(h, newpage = FALSE, padding = unit(c(2, 2, 2, 2), "mm")); popViewport() }
  lr <- min(C_raw[row(C_raw) != col(C_raw)]); dr(Heatmap(C_raw, name = "raw", col = colorRamp2(c(lr, 1), c("#f0efec", "#eb6834")), cluster_rows = clus, cluster_columns = clus, row_names_gp = gpar(fontsize = 8, col = GCOL, fontface = "bold"),
    column_names_gp = gpar(fontsize = 8, col = GCOL, fontface = "bold"), column_names_rot = 45, rect_gp = gpar(col = "white", lwd = 0.6), row_dend_width = unit(6, "mm"), column_dend_height = unit(6, "mm"),
    heatmap_legend_param = list(title = "Spearman", title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6)), column_title = "A  raw logCPM", column_title_gp = gpar(fontsize = 8.5, fontface = "bold"),
    cell_fun = function(j, i, x, y, w, h, fill) grid.text(sprintf("%.3f", C_raw[i, j]), x, y, gp = gpar(fontsize = 4.6))), 1)
  dr(hm(C_cen, "B  gene-centred", max(abs(C_cen[row(C_cen) != col(C_cen)])), clus), 2)
  dr(hm(C_fc, "C  NTC - KO, centred over KOs", max(abs(C_fc[row(C_fc) != col(C_fc)])), hclust(as.dist(1 - C_fc), "complete"), ann = FALSE), 3)
})
