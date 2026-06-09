# plot_SE_hub_comparison.R
#
# Heatmaps of SE enrichment across the parameter sweep:
#   CoAcc r=0.2-0.6, Stitching (consensus), Stitching_CT (per-CT peaks)
#   min_peaks = 3/4/5, TSS_dist = 0 / 2000 bp
#
# One PDF page per (tss_dist, min_peaks) combination, 7 heatmaps per row.
# Condition labels use the form: CoAccess_r0.3_mp3_tss0, Stitching_mp4_tss2000, etc.
# Also handles legacy condition labels without _tss suffix.
#
# Output: plot_SE_hub_comparison.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

script_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q7"

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
  "Smoothmuscle"                                                           = "Smooth Muscle",
  "Tcell"                                                                  = "T cell"
)

# ---- fixed row/column ordering (from r=0.3 reference) -----------------------
ct_order <- c("Malignant", "Mesothelium", "Alveolar",
               "Fibroblasts", "SmoothMuscle", "Endothelial",
               "Myeloid", "T_cells",
               "B_cells", "Plasma", "pDCs")

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

# ---- shared colour scale (99th-percentile global max) -----------------------
global_max <- quantile(results$neg_log10_fdr, 0.99, na.rm = TRUE)
col_fun    <- colorRamp2(c(0, global_max / 2, global_max),
                         c("#440154", "#31688e", "#fde725"))

# ---- helper: build score matrix in fixed order ------------------------------
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

row_labels_short <- ifelse(se_order %in% names(se_short), se_short[se_order], se_order)

# ---- detect parameter sweep structure ----------------------------------------
avail_conds <- unique(results$condition)

has_tss   <- any(grepl("_tss\\d+$", avail_conds))
tss_vals  <- if (has_tss)
               sort(unique(as.integer(sub(".*_tss(\\d+)$", "\\1",
                 avail_conds[grepl("_tss\\d+$", avail_conds)])))) else 0L

has_mp    <- any(grepl("_mp\\d+", avail_conds))
mp_values <- if (has_mp)
               sort(unique(as.integer(sub(".*_mp(\\d+).*", "\\1",
                 avail_conds[grepl("_mp\\d+", avail_conds)])))) else integer(0)

cor_base   <- c(paste0("CoAccess_r", c("0.2", "0.3", "0.4", "0.5", "0.6")),
                "Stitching", "Stitching_CT")
cor_labels <- c("CoAcc r=0.2", "CoAcc r=0.3", "CoAcc r=0.4",
                "CoAcc r=0.5", "CoAcc r=0.6", "Stitch", "Stitch CT")

# ---- heatmap factory (global idx ensures unique ComplexHeatmap names) --------
ht_idx <- 0L

make_ht <- function(cond, label, show_row = FALSE, show_legend = FALSE) {
  ht_idx <<- ht_idx + 1L
  mat <- .make_mat(cond)
  Heatmap(
    mat,
    name                 = paste0("h", ht_idx),
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

# ---- main: one page per (tss_dist, min_peaks) --------------------------------
if (has_mp) {

  pdf(file.path(script_dir, "plot_SE_hub_comparison.pdf"), width = 19, height = 5)

  page_i <- 0L

  for (tss_dist in tss_vals) {
    tss_label <- if (tss_dist == 0) "no TSS filter" else
                   sprintf("TSS_dist >= %d bp", tss_dist)

    for (mp in mp_values) {
      page_i <- page_i + 1L

      if (has_tss) {
        conds_all <- paste0(cor_base, sprintf("_mp%d_tss%d", mp, tss_dist))
      } else {
        conds_all <- paste0(cor_base, sprintf("_mp%d", mp))
      }
      if (length(intersect(conds_all, avail_conds)) == 0) next

      ht_row <- NULL
      for (k in seq_along(conds_all)) {
        ht <- make_ht(conds_all[k], cor_labels[k],
                      show_row    = (k == 1),
                      show_legend = (k == length(conds_all)))
        ht_row <- if (is.null(ht_row)) ht else ht_row + ht
      }

      draw(ht_row,
           column_title     = sprintf("SE enrichment  |  %s  |  min_peaks = %d",
                                      tss_label, mp),
           column_title_gp  = gpar(fontsize = 10, fontface = "bold"),
           ht_gap           = unit(3, "mm"),
           padding          = unit(c(3, 3, 6, 3), "mm"),
           newpage          = (page_i > 1))
    }
  }

  dev.off()

} else {
  # Legacy: no min_peaks sweep in CSV
  cond_order  <- intersect(c(cor_base, "Stitching_CT"), avail_conds)
  cond_labels <- c(cor_labels, "Stitch (CT)")[match(cond_order,
                   c(cor_base, "Stitching_CT"))]
  ht_list <- NULL
  for (k in seq_along(cond_order)) {
    ht <- make_ht(cond_order[k], cond_labels[k],
                  show_row    = (k == 1),
                  show_legend = (k == length(cond_order)))
    ht_list <- if (is.null(ht_list)) ht else ht_list + ht
  }
  pdf(file.path(script_dir, "plot_SE_hub_comparison.pdf"), width = 16, height = 5)
  draw(ht_list,
       column_title     = "SE enrichment in DA cHub regions  |  FDR-based score",
       column_title_gp  = gpar(fontsize = 10, fontface = "bold"),
       ht_gap           = unit(3, "mm"),
       padding          = unit(c(3, 3, 6, 3), "mm"))
  dev.off()
}

message("Saved plot_SE_hub_comparison.pdf")
