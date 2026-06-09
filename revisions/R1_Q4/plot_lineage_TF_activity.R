# plot_lineage_TF_activity.R  (R1_Q4)
#
# chromVAR TF activity (MotifMatrix z-scores) for known myeloid lineage TFs
# across meso myeloid subtypes (celltype_lv3).
#
# Plots:
#   1. Heatmap  — average z-score per celltype_lv3 × TF (row-scaled)
#   2. Violins  — per-cell z-score distributions per TF, split by celltype_lv3
#
# Outputs:
#   R1_Q4/output/lineage_TF_activity_heatmap.pdf
#   R1_Q4/output/lineage_TF_activity_violins.pdf

suppressPackageStartupMessages({
  library(ArchR)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(paletteer)
  library(SummarizedExperiment)
})

addArchRGenome("hg19")
addArchRThreads(4)

utils_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils"
source(file.path(utils_dir, "ggplot_aestetics.R"))

palette_myeloid <- rev(paletteer::paletteer_d("MoMAColors::VanGogh"))
palette_myeloid <- setNames(as.character(palette_myeloid),
  c("Mono_CD14","TAM_CXCLs","TAM_MARCO","TAM_TREM2","cDCs","Mono_CD16","TAM_interstitial"))

# ── Known myeloid lineage TFs ──────────────────────────────────────────────────
lineage_TFs <- c(
  # Master regulators
  "SPI1", "CEBPA", "CEBPB", "CEBPD",
  # Macrophage differentiation
  "MAFB", "MAF", "KLF4", "EGR1", "EGR2", "NR4A1",
  # DC / pDC
  "IRF4", "IRF8", "IRF1",
  # AP-1 / inflammatory
  "FOS", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND",
  # NF-kB
  "NFKB1", "NFKB2", "RELA", "REL",
  # STAT
  "STAT1", "STAT3",
  # Other myeloid
  "RUNX1", "MYC", "PPARG", "MAFB"
)
lineage_TFs <- unique(lineage_TFs)

ct_order <- c("Mono_CD14","Mono_CD16","TAM_CXCLs",
              "TAM_interstitial","TAM_MARCO","TAM_TREM2","cDCs")

ct_map <- c(cnmf_cluster_1="Mono_CD16", cnmf_cluster_2="TAM_MARCO",
            cnmf_cluster_3="Mono_CD14", cnmf_cluster_4="TAM_CXCLs",
            cnmf_cluster_5="TAM_TREM2", cnmf_cluster_6="TAM_interstitial",
            cnmf_cluster_7="cDCs")

out_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output"
proj_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR"

# ── Load project ───────────────────────────────────────────────────────────────
message("Loading ArchR project...")
archp <- loadArchRProject(proj_dir, showLogo = FALSE)

if (!"celltype_lv3" %in% colnames(archp@cellColData))
  archp$celltype_lv3 <- ct_map[as.character(archp$cnmf_cluster)]
archp$celltype_lv3 <- factor(as.character(archp$celltype_lv3), levels = ct_order)

# ── Average MotifMatrix per celltype_lv3 ───────────────────────────────────────
message("Averaging MotifMatrix per celltype_lv3...")
mSE  <- getGroupSE(archp, useMatrix = "MotifMatrix", groupBy = "celltype_lv3",
                   divideN = TRUE, scaleTo = NULL, threads = getArchRThreads())
mMat <- assays(mSE)[[1]]
rownames(mMat) <- rowData(mSE)$name
rownames(mMat) <- gsub("_.*", "", rownames(mMat))

# ── Subset to lineage TFs ──────────────────────────────────────────────────────
found <- intersect(lineage_TFs, rownames(mMat))
missing <- setdiff(lineage_TFs, rownames(mMat))
if (length(missing) > 0) message("Not found in MotifMatrix: ", paste(missing, collapse = ", "))
message("Plotting ", length(found), " TFs")

avg_mat <- mMat[found, ct_order, drop = FALSE]

# Row-scale (z-score across cell types per TF)
avg_scaled <- t(scale(t(avg_mat)))

# ── 1. Heatmap ─────────────────────────────────────────────────────────────────
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#D73027"))
ct_cols <- palette_myeloid[ct_order]

ht <- Heatmap(
  avg_scaled,
  name            = "z-score\n(row-scaled)",
  col             = col_fun,
  cluster_columns = FALSE,
  cluster_rows    = TRUE,
  show_row_names  = TRUE,
  show_column_names = TRUE,
  row_names_gp    = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  top_annotation  = HeatmapAnnotation(
    celltype = ct_order,
    col      = list(celltype = ct_cols),
    show_legend = FALSE,
    show_annotation_name = FALSE
  ),
  column_title    = "Myeloid lineage TF activity (chromVAR z-score)",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  width  = unit(5, "cm"),
  height = unit(length(found) * 0.4, "cm")
)

pdf(file.path(out_dir, "lineage_TF_activity_heatmap.pdf"),
    width = 6, height = max(4, length(found) * 0.35 + 2))
draw(ht)
dev.off()
message("Heatmap saved.")

# ── 2. Per-cell violin plots ───────────────────────────────────────────────────
message("Extracting per-cell MotifMatrix...")
motif_se   <- getMatrixFromProject(archp, useMatrix = "MotifMatrix")
motif_mat  <- assays(motif_se)[[1]]
rownames(motif_mat) <- gsub("_.*", "", rowData(motif_se)$name)

found_rows <- intersect(found, rownames(motif_mat))
cell_df <- as.data.frame(t(as.matrix(motif_mat[found_rows, , drop = FALSE])))
cell_df$celltype_lv3 <- archp$celltype_lv3[match(rownames(cell_df), getCellNames(archp))]
cell_df$celltype_lv3 <- factor(as.character(cell_df$celltype_lv3), levels = ct_order)

long_df <- cell_df |>
  pivot_longer(cols = all_of(found_rows), names_to = "TF", values_to = "z_score") |>
  filter(!is.na(celltype_lv3))

# TF order: by variance across cell types (most differential first)
tf_var <- tapply(long_df$z_score, long_df$TF, var, na.rm = TRUE)
tf_order <- names(sort(tf_var, decreasing = TRUE))
long_df$TF <- factor(long_df$TF, levels = tf_order)

p_vln <- ggplot(long_df, aes(x = celltype_lv3, y = z_score, fill = celltype_lv3)) +
  vlp + bxpv +
  facet_wrap(~ TF, ncol = 6, scales = "free_y") +
  scale_fill_manual(values = palette_myeloid, guide = "none") +
  labs(x = NULL, y = "chromVAR z-score",
       title = "Myeloid lineage TF activity by subtype") +
  gtheme +
  theme(
    plot.title       = element_text(size = 9),
    axis.text.x      = element_text(size = 5),
    axis.text.y      = element_text(size = 5),
    axis.title.y     = element_text(size = 7),
    strip.text       = element_text(size = 7, face = "bold.italic"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    panel.spacing    = unit(1.5, "mm")
  )

nrows_vln <- ceiling(length(found_rows) / 6)
ggsave(file.path(out_dir, "lineage_TF_activity_violins.pdf"),
       p_vln, width = 13, height = nrows_vln * 2.5)
message("Violins saved.")
