###############################################################################
# R2_Q14 -- heatmap of the significant BAP1-associated TFs across compartments.
#
# Rows    TF motifs that are a hit in at least one compartment (motif P < 0.01 in the
#         histology-adjusted model AND the TF gene score detectably accessible and moving in
#         the same direction; R2_Q14_BAP1_TF_by_compartment.R), split by where the hit sits
#         (tumour-intrinsic / shared / TME with tumour trend / TME-only).
# Columns Malignant, Myeloid, Stroma, T/NK, B/plasma (T/NK and B/plasma flagged: their hit
#         counts are not above chance by label permutation).
# Cells   moderated t, BAP1-lost vs retained (red = higher in BAP1-lost); "*" = hit there.
# Right   TCGA MESO inferred TF activity (Hmeljak 2018 tab 2B): colour = direction when
#         FDR < 0.05, pale = same/opposite but not significant, grey = TF not in the table.
#         Row names in bold = TCGA-concordant (tab 2B FDR < 0.05, same direction as the hit).
# Output: Plots/BAP1_TF_compartment_heatmap.pdf (horizontal: TFs as columns, compartments as rows)
#         Plots/BAP1_TF_compartment_heatmap_vertical.pdf (TFs as rows)
###############################################################################
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(RColorBrewer); library(paletteer); library(grid) })
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R2_Q14")); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "palettes.R"))
COMPS <- c("Malignant", "Myeloid", "Stroma", "TNK", "B_Plasma")
NAME  <- c(Malignant = "Malignant", Myeloid = "Myeloid", Stroma = "Stroma", TNK = "T/NK", B_Plasma = "B/plasma")
COL   <- c(Malignant = palette_celltype_lv1[["Malignant"]], Myeloid = palette_celltype_lv1[["Myeloid"]],
           Stroma = palette_celltype_lv1[["Fibroblasts"]], TNK = palette_celltype_lv1[["T_cells"]],
           B_Plasma = palette_celltype_lv1[["B_cells"]])
TLIM <- 5

A   <- read.csv("BAP1_TF_by_compartment_all.csv", stringsAsFactors = FALSE)
HT  <- read.csv("BAP1_TF_by_compartment_hits.csv", stringsAsFactors = FALSE)
FDR <- read.csv("BAP1_TF_by_compartment_permFDR.csv", stringsAsFactors = FALSE)
P2B <- read.csv("affinity_regression/prepared/paper_2B.csv", stringsAsFactors = FALSE)
uncal <- FDR$compartment[FDR$empirical_FDR > 0.25]

tf <- HT$motif
Tm <- sapply(COMPS, function(cp) { d <- A[A$compartment == cp, ]; d$motif_t_adj[match(tf, d$motif)] }); rownames(Tm) <- tf
Hm <- sapply(COMPS, function(cp) { d <- A[A$compartment == cp, ]; d$hit[match(tf, d$motif)] %in% TRUE }); rownames(Hm) <- tf
colnames(Tm) <- colnames(Hm) <- COMPS

loc_levels <- c("tumour-intrinsic", "shared", "TME, tumour trend", "TME-only", "discordant")
loc <- factor(HT$localisation, levels = intersect(loc_levels, unique(HT$localisation)))
ord <- order(loc, -apply(abs(Tm * Hm), 1, max, na.rm = TRUE))
Tm <- Tm[ord, ]; Hm <- Hm[ord, ]; HT <- HT[ord, ]; loc <- droplevels(loc[ord])

## TCGA tab 2B annotation relative to the hit direction
sym <- sub("\\.[0-9]+$", "", rownames(Tm))
est <- P2B$estimate_tf[match(sym, P2B$TF)]; pfdr <- P2B$p_adj_tf[match(sym, P2B$TF)]
hit_dir <- ifelse(HT$direction == "up in BAP1-lost", 1, ifelse(HT$direction == "down in BAP1-lost", -1, NA))
tcga <- ifelse(is.na(est), "not in table",
         ifelse(pfdr < 0.05, ifelse(est > 0, "higher, FDR<0.05", "lower, FDR<0.05"),
                ifelse(est > 0, "higher, ns", "lower, ns")))
concordant <- !is.na(est) & pfdr < 0.05 & sign(est) == hit_dir
cat("rows:", nrow(Tm), "| TCGA-concordant:", paste(rownames(Tm)[concordant], collapse = ", "), "\n")
print(table(localisation = loc)); print(table(tcga))

## project deviation palette (palettes.R: palette_deviation, ggthemes Red-Black-White Diverging),
## stops laid out as in palette_deviation_centered, rescaled from +-4 to +-TLIM
col_fun <- colorRamp2(c(TLIM, TLIM / 2, 0, -TLIM / 2, -TLIM),
                      c(palette_deviation[1], palette_deviation[27], "white", palette_deviation[65], palette_deviation[100]))
TCGA_COL <- c(`higher, FDR<0.05` = col_fun(TLIM), `lower, FDR<0.05` = col_fun(-TLIM),
              `higher, ns` = col_fun(TLIM * 0.3), `lower, ns` = col_fun(-TLIM * 0.3), `not in table` = "white")
## asterisk colour: white on dark cells, black on light ones
star_col <- function(fill) { rgb <- col2rgb(fill) / 255; if (sum(c(.299, .587, .114) * rgb) < .5) "white" else "black" }
col_lab <- ifelse(COMPS %in% uncal, paste0(NAME[COMPS], "*"), NAME[COMPS])
top <- HeatmapAnnotation(compartment = COMPS, col = list(compartment = COL), show_legend = FALSE,
                         show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"))
right <- rowAnnotation(`TCGA inferred\nTF activity` = tcga,
  col = list(`TCGA inferred\nTF activity` = TCGA_COL),
  simple_anno_size = unit(4, "mm"), annotation_name_gp = gpar(fontsize = 7), annotation_name_rot = 90, gp = gpar(col = "grey80", lwd = .4),
  annotation_legend_param = list(`TCGA inferred\nTF activity` = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                                     labels_gp = gpar(fontsize = 7), border = "grey70")))
ht <- Heatmap(Tm, name = "t", col = col_fun, na_col = "grey85",
  cluster_rows = FALSE, cluster_columns = FALSE, row_split = loc, row_gap = unit(1.5, "mm"),
  row_title_gp = gpar(fontsize = 8, fontface = "bold"), row_title_rot = 0,
  column_labels = col_lab, column_names_rot = 40, column_names_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 7, fontface = ifelse(concordant, "bold", "plain")),
  top_annotation = top, right_annotation = right,
  rect_gp = gpar(col = "white", lwd = .6),
  cell_fun = function(j, i, x, y, w, h, fill) if (Hm[i, j]) grid.text("*", x, y - unit(.6, "mm"), gp = gpar(fontsize = 10, fontface = "bold", col = star_col(fill))),
  heatmap_legend_param = list(title = "t (BAP1-lost\nvs retained)", at = c(-TLIM, 0, TLIM),
                              title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)),
  column_title = "BAP1-associated TF activity by compartment (histology-adjusted)",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  width = unit(5 * 7, "mm"), height = unit(nrow(Tm) * 3.4, "mm"))
pdf("Plots/BAP1_TF_compartment_heatmap_vertical.pdf", width = 5.6, height = nrow(Tm) * 3.4 / 25.4 + 3.2)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(24, 4, 4, 4), "mm"))          # bottom, left, top, right: room for the footnote
grid.text(paste("* in cell: hit (motif P < 0.01, gene score concordant)",
                "* after compartment: hit count not above chance",
                "bold TF: same direction, FDR < 0.05 in TCGA inferred TF activity (tab 2B)", sep = "\n"),
          x = unit(6, "mm"), y = unit(7, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 6.5, col = "grey30"))
dev.off()
cat("DONE\n")

## ---- horizontal version: TFs as columns, compartments as rows -------------------------
Tt <- t(Tm); Ht <- t(Hm)
loc_title <- c(`tumour-intrinsic` = "tumour-intrinsic", shared = "shared",
               `TME, tumour trend` = "TME, tumour trend", `TME-only` = "TME-only", discordant = "discordant")
top_h <- HeatmapAnnotation(`TCGA inferred TF activity` = tcga,
  col = list(`TCGA inferred TF activity` = TCGA_COL),
  simple_anno_size = unit(4, "mm"), annotation_name_gp = gpar(fontsize = 7), annotation_name_side = "left", gp = gpar(col = "grey80", lwd = .4),
  annotation_legend_param = list(`TCGA inferred TF activity` = list(title = "TCGA inferred\nTF activity", title_gp = gpar(fontsize = 8, fontface = "bold"),
                                                                   labels_gp = gpar(fontsize = 7), border = "grey70")))
left_h <- rowAnnotation(compartment = COMPS, col = list(compartment = COL), show_legend = FALSE,
                        show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"))
hh <- Heatmap(Tt, name = "t_h", col = col_fun, na_col = "grey85",
  cluster_rows = FALSE, cluster_columns = FALSE, column_split = loc, column_gap = unit(1.5, "mm"),
  column_title = loc_title[levels(loc)], column_title_gp = gpar(fontsize = 8, fontface = "bold"),
  row_labels = col_lab, row_names_side = "left", row_names_gp = gpar(fontsize = 9),
  column_names_rot = 45, column_names_gp = gpar(fontsize = 6.5, fontface = ifelse(concordant, "bold", "plain")),
  top_annotation = top_h, left_annotation = left_h,
  rect_gp = gpar(col = "white", lwd = .6),
  cell_fun = function(j, i, x, y, w, h, fill) if (Ht[i, j]) grid.text("*", x, y - unit(.6, "mm"), gp = gpar(fontsize = 10, fontface = "bold", col = star_col(fill))),
  heatmap_legend_param = list(title = "t (BAP1-lost\nvs retained)", at = c(-TLIM, 0, TLIM),
                              title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)),
  width = unit(ncol(Tt) * 4.4, "mm"), height = unit(5 * 6, "mm"))
pdf("Plots/BAP1_TF_compartment_heatmap.pdf", width = ncol(Tt) * 4.4 / 25.4 + 4.6, height = 3.6)
draw(hh, heatmap_legend_side = "right", annotation_legend_side = "right",
     column_title = "BAP1-associated TF activity by compartment (histology-adjusted)",
     column_title_gp = gpar(fontsize = 10, fontface = "bold"), padding = unit(c(8, 4, 2, 10), "mm"))
grid.text(paste("* in cell: hit (motif P < 0.01, gene score concordant)   |   * after compartment: hit count not above chance",
                "   |   bold TF: same direction, FDR < 0.05 in TCGA inferred TF activity (tab 2B)"),
          x = unit(6, "mm"), y = unit(4, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 6.5, col = "grey30"))
dev.off()
cat("horizontal DONE\n")
