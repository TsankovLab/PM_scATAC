# plot_MDM_TF_activity.R  (R1_Q4)
#
# Visualises the 14 MDM/inflammatory-axis TFs across:
#   - NSCLC blood monocytes (Mye_CD14.Mono, Mye_CD16.Mono)
#   - Meso myeloid cell types (Mono_CD14, Mono_CD16, TAM_CXCLs,
#                              TAM_interstitial, TAM_MARCO, TAM_TREM2, cDCs)
#
# TF activity = chromVAR z-scores averaged per cell type.
# Blood source: average_TF_activity_PBMC.csv (pbmc_myeloid project)
# Meso source : MotifMatrix from meso myeloid ArchR project
#
# Outputs:
#   plot_MDM_TF_heatmap.pdf        -- heatmap rows=TFs, cols=all cell types
#   plot_MDM_TF_mean_barplot.pdf   -- mean activity across 14 TFs per cell type
#
# Runs on cluster (~5 min); submit via submit_meso_blood_comparison.sh
# Or run interactively: OPENBLAS_NUM_THREADS=1 Rscript plot_MDM_TF_activity.R

suppressPackageStartupMessages({
  library(ArchR)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
addArchRGenome("hg38")
addArchRThreads(threads = 4)

script_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4"
meso_dir     <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR"
blood_tf_csv <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/pbmc_myeloid/average_TF_activity_PBMC.csv"

MDM_TFs <- c("FOSL1","FOSL2","BACH1","PPARG","NFKB2","KLF12",
             "HIVEP3","SMAD1","NFKB1","REL","RUNX1","SNAI1","RUNX2","NFAT5")

meso_ct_order  <- c("Mono_CD14","Mono_CD16","TAM_CXCLs",
                    "TAM_interstitial","TAM_MARCO","TAM_TREM2","cDCs")
blood_ct_order <- c("Mye_CD14.Mono","Mye_CD16.Mono")

palette_myeloid <- c(
  Mono_CD14        = "#d6604d",
  Mono_CD16        = "#f4a582",
  TAM_CXCLs        = "#4393c3",
  TAM_interstitial = "#99d8c9",
  TAM_MARCO        = "#2ca25f",
  TAM_TREM2        = "#74c476",
  cDCs             = "#8e7cb8"
)

# ---- 1. Load meso TF activity ------------------------------------------------
message("Loading meso ArchR project...")
archp <- loadArchRProject(meso_dir, showLogo = FALSE)

# Restore celltype_lv3 if absent
ct_map <- c(cnmf_cluster_1="Mono_CD16", cnmf_cluster_2="TAM_MARCO",
            cnmf_cluster_3="Mono_CD14", cnmf_cluster_4="TAM_CXCLs",
            cnmf_cluster_5="TAM_TREM2", cnmf_cluster_6="TAM_interstitial",
            cnmf_cluster_7="cDCs")
if (!"celltype_lv3" %in% colnames(archp@cellColData))
  archp$celltype_lv3 <- ct_map[as.character(archp$cnmf_cluster)]
archp$celltype_lv3 <- factor(as.character(archp$celltype_lv3), levels = meso_ct_order)

message("Averaging MotifMatrix per celltype_lv3...")
mSE  <- getGroupSE(archp, useMatrix = "MotifMatrix", groupBy = "celltype_lv3",
                   divideN = TRUE, scaleTo = NULL, threads = getArchRThreads())
mMat <- assays(mSE)[[1]]
rownames(mMat) <- rowData(mSE)$name
rownames(mMat) <- gsub("_.*", "", rownames(mMat))   # strip motif ID suffix

# ---- 2. Load blood TF activity -----------------------------------------------
message("Loading blood TF activity...")
blood_tf <- read.csv(blood_tf_csv, stringsAsFactors = FALSE)
blood_tf <- blood_tf[seq_len(min(870, nrow(blood_tf))), ]
rownames(blood_tf) <- blood_tf[, 1]
blood_tf <- blood_tf[, -1, drop = FALSE]

# ---- 3. Intersect on MDM_TFs -------------------------------------------------
meso_tfs  <- intersect(MDM_TFs, rownames(mMat))
blood_tfs <- intersect(MDM_TFs, rownames(blood_tf))
common_tfs <- intersect(meso_tfs, blood_tfs)

missing <- setdiff(MDM_TFs, common_tfs)
if (length(missing))
  message("  TFs not found in both datasets: ", paste(missing, collapse = ", "))
message(sprintf("  Using %d / %d MDM TFs", length(common_tfs), length(MDM_TFs)))

# Sub-matrices for the common TFs
meso_sub  <- mMat[common_tfs,
                  meso_ct_order[meso_ct_order %in% colnames(mMat)], drop = FALSE]
blood_sub <- as.matrix(blood_tf[common_tfs,
                  blood_ct_order[blood_ct_order %in% colnames(blood_tf)], drop = FALSE])

# Pretty column labels
colnames(blood_sub) <- c("Mye_CD14.Mono" = "Blood\nCD14 Mono",
                          "Mye_CD16.Mono" = "Blood\nCD16 Mono")[colnames(blood_sub)]

# Combined matrix: blood | meso
combined <- cbind(blood_sub, meso_sub)

# ---- 4. Heatmap --------------------------------------------------------------
message("Building heatmap...")

# Row-scale (z-score per TF across all cell types so blood and meso are comparable)
combined_z <- t(scale(t(combined)))

# Colour scale capped at ±2 sd
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#d6604d"))

# Column split: Blood | Meso
col_split <- factor(c(rep("Blood (NSCLC)", ncol(blood_sub)),
                      rep("Mesothelioma",  ncol(meso_sub))),
                    levels = c("Blood (NSCLC)", "Mesothelioma"))

# Column annotation: cell type colour bar
ct_col_vec <- c(
  setNames(rep("#1f78b4", ncol(blood_sub)), colnames(blood_sub)),
  palette_myeloid[meso_ct_order[meso_ct_order %in% colnames(meso_sub)]]
)

ha_col <- HeatmapAnnotation(
  `Cell type` = colnames(combined_z),
  col = list(`Cell type` = ct_col_vec),
  show_legend  = FALSE,
  annotation_name_gp = gpar(fontsize = 7)
)

ht <- Heatmap(
  combined_z,
  name              = "z-score",
  col               = col_fun,
  top_annotation    = ha_col,
  column_split      = col_split,
  cluster_rows      = TRUE,
  cluster_columns   = FALSE,
  clustering_distance_rows = "euclidean",
  row_names_gp      = gpar(fontsize = 9, fontface = "italic"),
  row_names_side    = "left",
  column_names_gp   = gpar(fontsize = 8),
  column_names_rot  = 45,
  column_title_gp   = gpar(fontsize = 9, fontface = "bold"),
  border            = TRUE,
  rect_gp           = gpar(col = "white", lwd = 0.5),
  width             = unit(ncol(combined_z) * 6, "mm"),
  height            = unit(length(common_tfs) * 6, "mm"),
  heatmap_legend_param = list(
    title      = "z-score\n(row-scaled)",
    title_gp   = gpar(fontsize = 8, fontface = "bold"),
    labels_gp  = gpar(fontsize = 7),
    legend_height = unit(3, "cm"),
    at = c(-2, -1, 0, 1, 2)
  )
)

pdf(file.path(script_dir, "plot_MDM_TF_heatmap.pdf"), width = 8, height = 5.5)
draw(ht,
     column_title    = "MDM / inflammatory-axis TF activity (chromVAR z-score)",
     column_title_gp = gpar(fontsize = 10, fontface = "bold"),
     padding         = unit(c(4, 4, 8, 4), "mm"))
dev.off()
message("Saved plot_MDM_TF_heatmap.pdf")

# ---- 5. Mean MDM TF activity barplot ----------------------------------------
# Average z-score across the 14 TFs per cell type — single summary score
mean_act <- colMeans(combined_z, na.rm = TRUE)
bar_df <- data.frame(
  celltype = names(mean_act),
  mean_z   = as.numeric(mean_act),
  source   = ifelse(grepl("Blood", names(mean_act)), "Blood (NSCLC)", "Mesothelioma"),
  stringsAsFactors = FALSE
)
bar_df$celltype <- factor(bar_df$celltype,
  levels = bar_df$celltype[order(bar_df$mean_z, decreasing = TRUE)])

bar_fill <- c(
  "Blood\nCD14 Mono"   = "#1f78b4",
  "Blood\nCD16 Mono"   = "#a6cee3",
  palette_myeloid
)

p_bar <- ggplot(bar_df, aes(x = celltype, y = mean_z, fill = celltype)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey50") +
  scale_fill_manual(values = bar_fill) +
  labs(
    x     = NULL,
    y     = "Mean z-score across 14 MDM TFs",
    title = "MDM / inflammatory TF activity"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

pdf(file.path(script_dir, "plot_MDM_TF_mean_barplot.pdf"), width = 5, height = 4)
print(p_bar)
dev.off()
message("Saved plot_MDM_TF_mean_barplot.pdf")

# ---- 6. Save numeric table ---------------------------------------------------
out_df <- as.data.frame(combined_z)
out_df$TF <- rownames(out_df)
write.csv(out_df[, c("TF", colnames(combined_z))],
          file.path(script_dir, "MDM_TF_activity_table.csv"),
          row.names = FALSE, quote = FALSE)

message("\nDone. Outputs in ", script_dir)
