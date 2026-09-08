###############################################################################
# STEP 4 -- overview figure: the subclones epiAneufinder finds in every tumour.
#
#  A  CNV landscape of every clone the ARM_MIN criterion calls, rows grouped by tumour,
#     columns = the shared 5 Mb bins in genome order, split by chromosome.
#     Colour = (fraction of cells called GAIN) - (fraction called LOSS).
#     A tumour whose primary split fell below ARM_MIN contributes a single row.
#  B  What separates the two clones of the tumours that DID split: the per-bin
#     difference (minor clone c2 - major clone c1), with the size of the largest
#     difference as a bar. A real focal event shows as a solid block on one arm
#     (P4, chr8q); the tumours left as one clone had nothing of that size to show.
#
# Input : epi_clone_profiles.rds (step 3)
# Output: Plots/epiAneufinder_subclones_overview.pdf
###############################################################################
suppressMessages({ library(ComplexHeatmap); library(circlize); library(grid); library(gridExtra) })
source("00_common.R")
pdf(NULL)          # keep grid.grabExpr() from leaving a stray Rplots.pdf behind
P <- readRDS("epi_clone_profiles.rds")
Z <- P$Z; DL <- P$DL; meta <- P$meta; co <- P$coord; SPLIT <- P$split
chrf <- factor(co$chr, levels = paste0("chr", 1:22))

col_fun <- colorRamp2(c(-0.6, -0.15, 0, 0.15, 0.6),
                      c("#08519c", "#9ecae1", "white", "#fb9a99", "#a50f15"))
star <- ifelse(meta$leaf == P$chr8q_leaf, "<- chr8q clone", "")

## ---- panel A ----------------------------------------------------------------
## NOTE: do not pass a scalar `column_title` to Heatmap() when column_split is set --
## it replaces the per-slice chromosome labels. Panel titles go in via draw() below.
htA <- Heatmap(Z, name = "gain - loss", col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = factor(meta$sample, levels = SAMPLES),
  row_title_rot = 0, row_title_gp = gpar(fontsize = 8, fontface = "bold"),
  row_gap = unit(0.8, "mm"),
  column_split = chrf, column_gap = unit(0.4, "mm"), column_title_gp = gpar(fontsize = 6.5),
  show_column_names = FALSE,
  row_labels = sprintf("%s  (n=%s)%s", meta$leaf, format(meta$n_cells, big.mark = ","),
                       ifelse(meta$driver == "unsplit", "  one clone", "")),
  row_names_gp = gpar(fontsize = 6.5),
  left_annotation = rowAnnotation(
    tumour = factor(meta$sample, levels = SAMPLES),   # factor -> legend follows row order
    cells  = anno_barplot(meta$n_cells, width = unit(11, "mm"), gp = gpar(fill = "grey65", col = NA)),
    col = list(tumour = SAMPCOL),
    annotation_name_gp = gpar(fontsize = 6.5), simple_anno_size = unit(3.5, "mm")),
  right_annotation = rowAnnotation(
    hl = anno_text(star, gp = gpar(fontsize = 6.5, col = COL_HI, fontface = "bold"))),
  border = TRUE,
  heatmap_legend_param = list(legend_height = unit(20, "mm"),
                              title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6)))

## ---- panel B ----------------------------------------------------------------
dlab <- sprintf("%s   %s %s", SPLIT, P$amax[, "arm"],
                ifelse(as.numeric(P$amax[, "delta"]) < 0, "loss", "gain"))
htB <- Heatmap(DL, name = "clone2 - clone1",
  col = colorRamp2(c(-0.6, 0, 0.6), c("#08519c", "white", "#a50f15")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = chrf, column_gap = unit(0.4, "mm"), column_title_gp = gpar(fontsize = 6.5),
  show_column_names = FALSE, row_labels = dlab, row_names_gp = gpar(fontsize = 7),
  left_annotation = rowAnnotation(tumour = factor(SPLIT, levels = SAMPLES),
    col = list(tumour = SAMPCOL),
    show_legend = FALSE, annotation_name_gp = gpar(fontsize = 6.5),
    simple_anno_size = unit(3.5, "mm")),
  right_annotation = rowAnnotation(
    `max|d|` = anno_barplot(P$strength, width = unit(13, "mm"), ylim = c(0, 1),
      gp = gpar(fill = ifelse(SPLIT == "P4", COL_HI, "grey65"), col = NA)),
    annotation_name_gp = gpar(fontsize = 6.5)),
  border = TRUE,
  heatmap_legend_param = list(legend_height = unit(16, "mm"),
                              title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6)))

gA <- grid.grabExpr(draw(htA, heatmap_legend_side = "right", annotation_legend_side = "right",
        column_title = sprintf("A   CNV landscape of every clone (%d clones in %d tumours)",
                               nrow(Z), length(unique(meta$sample))),
        column_title_gp = gpar(fontface = "bold", fontsize = 10)))
gB <- grid.grabExpr(draw(htB, heatmap_legend_side = "right", annotation_legend_side = "right",
        column_title = sprintf(paste("B   What separates the two clones of the %d tumours",
                                     "that split (per-bin difference, c2 - c1)"), length(SPLIT)),
        column_title_gp = gpar(fontface = "bold", fontsize = 10)))
pdf("Plots/epiAneufinder_subclones_overview.pdf", width = 13.5, height = 10)
grid.arrange(gA, gB, ncol = 1, heights = c(1.9, 1),
  top = textGrob(sprintf(paste("epiAneufinder clones across all scATAC tumours (5 Mb windows,",
                               "malignant cells; split accepted at arm-level separation >= %.2f)"), ARM_MIN),
                 gp = gpar(fontface = "bold", fontsize = 13)))
invisible(dev.off())
cat("DONE -> Plots/epiAneufinder_subclones_overview.pdf\n")
