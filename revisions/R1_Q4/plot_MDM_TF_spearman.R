# plot_MDM_TF_spearman.R  (R1_Q4)
#
# Pairwise Spearman correlations between all PBMC and meso cell types
# using the 14 MDM/inflammatory-axis TF chromVAR z-scores.
# Input: MDM_TF_chromvar_combined.csv (from run_chromvar_combined.R)
# Run on login node: Rscript plot_MDM_TF_spearman.R

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(reshape2)
})

out_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4"
csv_in  <- file.path(out_dir, "MDM_TF_chromvar_combined.csv")

# ---- Load data ---------------------------------------------------------------
df <- read.csv(csv_in, row.names = 1, stringsAsFactors = FALSE)
# df: rows = 14 TFs, cols = cell types

# ---- Pairwise Spearman correlation -------------------------------------------
# Each column is the TF activity profile of one cell type (14-dimensional vector)
# Correlate all pairs: cor(df) operates on columns
cor_mat <- cor(df, method = "spearman")

# ---- Annotations -------------------------------------------------------------
ct_names <- colnames(cor_mat)

source_vec <- ifelse(grepl("^PBMC", ct_names), "Blood (NSCLC)", "Mesothelioma")
source_pal <- c("Blood (NSCLC)" = "#1f78b4", "Mesothelioma" = "#e31a1c")

ct_pal <- c(
  PBMC_CD14_Mono          = "#d6604d",
  PBMC_CD16_Mono          = "#f4a582",
  PBMC_DC                 = "#8e7cb8",
  Mesothelioma_CD14_Mono  = "#b2182b",
  Mesothelioma_CD16_Mono  = "#fc8d59",
  Mesothelioma_DC         = "#542788",
  TAM_CXCLs               = "#4393c3",
  TAM_interstitial        = "#99d8c9",
  TAM_MARCO               = "#2ca25f",
  TAM_TREM2               = "#74c476"
)

# Pretty labels
pretty_labels <- gsub("PBMC_",         "Blood ",      ct_names)
pretty_labels <- gsub("Mesothelioma_", "Meso ",       pretty_labels)
pretty_labels <- gsub("_", " ",                       pretty_labels)

ha_col <- HeatmapAnnotation(
  Source    = source_vec,
  col       = list(Source = source_pal),
  annotation_name_gp = gpar(fontsize = 7.5),
  show_legend = TRUE
)
ha_row <- rowAnnotation(
  Source    = source_vec,
  col       = list(Source = source_pal),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

col_fun <- colorRamp2(c(-1, 0, 1), c("#2166ac", "white", "#d6604d"))

# ---- Heatmap -----------------------------------------------------------------
ht <- Heatmap(
  cor_mat,
  name              = "Spearman r",
  col               = col_fun,
  top_annotation    = ha_col,
  left_annotation   = ha_row,
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  clustering_method_rows    = "complete",
  clustering_method_columns = "complete",
  clustering_distance_rows    = "euclidean",
  clustering_distance_columns = "euclidean",
  row_labels        = pretty_labels,
  column_labels     = pretty_labels,
  row_names_gp      = gpar(fontsize = 9),
  column_names_gp   = gpar(fontsize = 9),
  column_names_rot  = 45,
  border            = TRUE,
  rect_gp           = gpar(col = "white", lwd = 0.5),
  cell_fun          = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", cor_mat[i, j]), x, y,
              gp = gpar(fontsize = 7,
                        col = ifelse(abs(cor_mat[i, j]) > 0.5, "white", "grey20")))
  },
  heatmap_legend_param = list(
    title_gp      = gpar(fontsize = 8, fontface = "bold"),
    labels_gp     = gpar(fontsize = 7),
    legend_height = unit(3, "cm"),
    at            = c(-1, -0.5, 0, 0.5, 1)
  )
)

pdf(file.path(out_dir, "MDM_TF_spearman_heatmap.pdf"), width = 7, height = 6)
draw(ht,
     column_title    = "Pairwise Spearman correlation of 14 MDM TF activity profiles",
     column_title_gp = gpar(fontsize = 10, fontface = "bold"),
     padding         = unit(c(4, 4, 8, 4), "mm"))
dev.off()
message("Saved MDM_TF_spearman_heatmap.pdf")

# ---- Save correlation table --------------------------------------------------
cor_df        <- as.data.frame(cor_mat)
cor_df$cell_type <- rownames(cor_df)
write.csv(cor_df[, c("cell_type", colnames(cor_mat))],
          file.path(out_dir, "MDM_TF_spearman_matrix.csv"),
          row.names = FALSE, quote = FALSE)

# ---- Print summary -----------------------------------------------------------
message("\nSpearman correlations with PBMC CD14 Mono (sorted):")
ref <- cor_mat["PBMC_CD14_Mono", ]
ref <- sort(ref[names(ref) != "PBMC_CD14_Mono"], decreasing = TRUE)
print(round(ref, 3))

message("\nDone.")
