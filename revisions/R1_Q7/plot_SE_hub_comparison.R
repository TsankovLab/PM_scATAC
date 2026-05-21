# plot_SE_hub_comparison.R
#
# Six side-by-side heatmaps (one per hub condition: CoAcc r0.2-0.6 + Stitching).
# Rows = SE databases, columns = cell types.
# Row + column order fixed from r=0.3 reference (hclust, diagonal arrangement).
# Same colour scale across all panels.
#
# Output: plot_SE_hub_comparison.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

script_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude"

results <- read.csv(file.path(script_dir, "myeloid_SE_hub_enrichment.csv"),
                    stringsAsFactors = FALSE)

# ---- short SE labels ---------------------------------------------------------
se_short <- c(
  "AT1"                                                                    = "AT1",
  "AT2"                                                                    = "AT2",
  "B.cells"                                                                = "B cells",
  "ENCODE_Blood_primary_cell_Bcell"                                        = "Blood B cell",
  "ENCODE_Blood_primary_cell_CD14-positive-monocyte"                       = "Blood Monocyte",
  "ENCODE_Bone_primary_cell_osteoblast"                                    = "Bone Osteoblast",
  "ENCODE_Brain_and_spinal_cord_primary_cell_astrocyte"                    = "Brain Astrocyte",
  "ENCODE_Cervix_cell_line_HeLa-S3"                                        = "HeLa-S3",
  "ENCODE_Embryo_in_vitro_differentiated_cells_cardiac-muscle-cell"        = "Cardiac Muscle",
  "ENCODE_Embryo_primary_cell_neural-progenitor-cell"                      = "Neural Progenitor",
  "ENCODE_Embryo_stem_cell_H9"                                             = "H9 ESC",
  "ENCODE_Lung_primary_cell_fibroblast-of-lung"                            = "Lung Fibroblast",
  "ENCODE_Mammary_gland_primary_cell_mammary-epithelial-cell"              = "Mammary Epithelial",
  "ENCODE_Prostate_cell_line_C4-2B"                                        = "Prostate C4-2B",
  "ENCODE_Skin_Induced_pluripotent_stem_cell_line_GM23338"                 = "iPSC GM23338",
  "ENCODE_Skin_in_vitro_differentiated_cells_bipolar-neuron"               = "Bipolar Neuron",
  "ENCODE_Skin_primary_cell_keratinocyte"                                  = "Keratinocyte",
  "ENCODE_Smooth_muscle_in_vitro_differentiated_cells_smooth-muscle-cell"  = "Smooth Muscle",
  "ENCODE_Submandibular_salivary_gland_cell_line_ACC112"                   = "Salivary Gland",
  "ENCODE_Umbilical_vein_primary_cell_endothelial-cell-of-umbilical-vein"  = "HUVEC",
  "Endothelial"                                                            = "Endothelial",
  "NK"                                                                     = "NK",
  "Smoothmuscle"                                                           = "Smooth Muscle*",
  "Tcell"                                                                  = "T cell"
)

# ---- ordering from r=0.3 reference -------------------------------------------
# fixed column order (ct1) and row order (ct2) as specified
ct_order <- c("Malignant","Mesothelium","Alveolar",
               "Fibroblasts","SmoothMuscle","Endothelial",
               "Myeloid","T_cells",
               "B_cells","Plasma","pDCs")

se_order <- c(
  "ENCODE_Mammary_gland_primary_cell_mammary-epithelial-cell",
  "ENCODE_Skin_primary_cell_keratinocyte",
  "AT1", "AT2",
  "ENCODE_Lung_primary_cell_fibroblast-of-lung",
  "Smoothmuscle",
  "Endothelial",
  "ENCODE_Blood_primary_cell_CD14-positive-monocyte",
  "Tcell",
  "ENCODE_Blood_primary_cell_Bcell"
)

ct_order <- intersect(ct_order, unique(results$celltype))
se_order <- intersect(se_order, unique(results$se_db))

# ---- shared colour scale (global max across all conditions) ------------------
global_max <- quantile(results$neg_log10_fdr, 0.99, na.rm = TRUE)
col_fun    <- colorRamp2(c(0, global_max / 2, global_max),
                         c("#440154", "#31688e", "#fde725"))

# ---- helper: build matrix in fixed order ------------------------------------
.make_mat <- function(cond) {
  df  <- results[results$condition == cond &
                   results$se_db %in% se_order &
                   results$celltype %in% ct_order, ]
  mat <- matrix(0, nrow = length(se_order), ncol = length(ct_order),
                dimnames = list(se_order, ct_order))
  for (i in seq_len(nrow(df)))
    mat[df$se_db[i], df$celltype[i]] <- df$neg_log10_fdr[i]
  mat
}

# ---- build 7 heatmaps --------------------------------------------------------
all_conds   <- c("CoAccess_r0.2","CoAccess_r0.3","CoAccess_r0.4",
                 "CoAccess_r0.5","CoAccess_r0.6","Stitching","Stitching_CT")
all_labels  <- c("CoAcc r=0.2","CoAcc r=0.3","CoAcc r=0.4",
                 "CoAcc r=0.5","CoAcc r=0.6","Stitch (union)","Stitch (CT)")
cond_order  <- intersect(all_conds, unique(results$condition))
cond_labels <- all_labels[match(cond_order, all_conds)]

row_labels_short <- ifelse(se_order %in% names(se_short),
                           se_short[se_order], se_order)

sig_thresh <- -log10(0.05)   # dotted significance line (for colour, not drawn)

make_ht <- function(cond, label, show_row = FALSE, show_legend = FALSE, idx = 1) {
  mat <- .make_mat(cond)
  Heatmap(
    mat,
    name                 = paste0("h", idx),   # unique internal name per panel
    col                  = col_fun,
    cluster_rows         = FALSE,
    cluster_columns      = FALSE,
    row_labels           = if (show_row) row_labels_short else rep("", nrow(mat)),
    row_names_gp         = gpar(fontsize = 6.5),
    row_names_side       = "left",
    column_labels        = ct_order,
    column_names_gp      = gpar(fontsize = 7),
    column_names_rot     = 45,
    column_title         = label,
    column_title_gp      = gpar(fontsize = 8.5, fontface = "bold"),
    column_title_side    = "top",
    width                = unit(length(ct_order)  * 3.5, "mm"),
    height               = unit(length(se_order)  * 3.2, "mm"),
    border               = TRUE,
    rect_gp              = gpar(col = NA),
    show_heatmap_legend  = show_legend,
    heatmap_legend_param = list(
      title         = "-log10(FDR)",
      title_gp      = gpar(fontsize = 8, fontface = "bold"),
      labels_gp     = gpar(fontsize = 7),
      legend_height = unit(3, "cm"),
      at            = c(0, round(global_max / 2, 1), round(global_max, 1))
    )
  )
}

ht_list <- NULL
for (k in seq_along(cond_order)) {
  ht <- make_ht(cond_order[k], cond_labels[k],
                show_row    = (k == 1),
                show_legend = (k == length(cond_order)),
                idx         = k)
  ht_list <- if (is.null(ht_list)) ht else ht_list + ht
}

# ---- save --------------------------------------------------------------------
pdf(file.path(script_dir, "plot_SE_hub_comparison.pdf"), width = 16, height = 5)
draw(ht_list,
     column_title     = "SE enrichment in DA cHub regions  |  FDR-based score  |  rows/cols ordered by r=0.3 reference",
     column_title_gp  = gpar(fontsize = 10, fontface = "bold"),
     ht_gap           = unit(3, "mm"),
     padding          = unit(c(3, 3, 6, 3), "mm"))
dev.off()

message("Saved plot_SE_hub_comparison.pdf")
