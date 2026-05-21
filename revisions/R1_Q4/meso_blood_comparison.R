# meso_blood_comparison.R  (R1_Q4)
#
# Epigenomic comparison of meso myeloid cell types with NSCLC blood monocytes.
#
# Analysis 1 — chromVAR deviations using blood peak sets as reference:
#   For each meso cell, compute how much it opens blood CD14/CD16 Mono peaks
#   above background. Compares inflamed (high mod_2 = MDM) vs non_inflamed
#   (low mod_2 = resident macrophage) cells.
#
# Analysis 2 — TF activity Spearman correlation:
#   Average chromVAR TF z-scores per meso cell type correlated with
#   blood monocyte TF activity (from average_TF_activity_PBMC.csv).
#
# Outputs:
#   plot_chromVAR_blood_peaks_heatmap.pdf
#   plot_chromVAR_blood_peaks_violin.pdf
#   plot_TF_correlation_blood.pdf
#   TF_correlation_blood_results.csv
#
# Requires cluster: 4 CPUs, 64 GB, ~1.5 h
# Run via submit_meso_blood_comparison.sh

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})
addArchRGenome("hg38")
addArchRThreads(threads = 4)

script_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4"
meso_dir     <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR"
pbmc_dir     <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/pbmc_myeloid"
blood_tf_csv <- file.path(pbmc_dir, "average_TF_activity_PBMC.csv")
pbmc_peaks   <- file.path(pbmc_dir, "PeakCalls")

# ---- Cell-type palette -------------------------------------------------------
cell_order <- c("Mono_CD14","Mono_CD16","TAM_CXCLs","TAM_TREM2",
                "TAM_MARCO","TAM_interstitial","cDCs")

palette_myeloid <- c(
  Mono_CD14        = "#d6604d",
  Mono_CD16        = "#f4a582",
  TAM_CXCLs        = "#4393c3",
  TAM_TREM2        = "#92c5de",
  TAM_MARCO        = "#2ca25f",
  TAM_interstitial = "#99d8c9",
  cDCs             = "#8e7cb8"
)

# cnmf_cluster -> celltype_lv3 mapping (scatac_myeloid_ArchR.R l228-235)
ct_map <- c(
  cnmf_cluster_1 = "Mono_CD16",
  cnmf_cluster_2 = "TAM_MARCO",
  cnmf_cluster_3 = "Mono_CD14",
  cnmf_cluster_4 = "TAM_CXCLs",
  cnmf_cluster_5 = "TAM_TREM2",
  cnmf_cluster_6 = "TAM_interstitial",
  cnmf_cluster_7 = "cDCs"
)

# ---- Load meso ArchR project -------------------------------------------------
message("Loading meso myeloid ArchR project...")
archp <- loadArchRProject(meso_dir, showLogo = FALSE)

# Restore celltype_lv3 if absent
if (!"celltype_lv3" %in% colnames(archp@cellColData)) {
  message("  Restoring celltype_lv3 from cnmf_cluster...")
  archp$celltype_lv3 <- ct_map[as.character(archp$cnmf_cluster)]
}
archp$celltype_lv3 <- factor(as.character(archp$celltype_lv3), levels = cell_order)

# Restore inflamed label (mod_2 must be in cellColData)
if (!"inflamed" %in% colnames(archp@cellColData)) {
  if ("mod_2" %in% colnames(archp@cellColData)) {
    message("  Restoring inflamed label from mod_2...")
    archp$inflamed <- ifelse(archp$mod_2 > 0, "inflamed", "non_inflamed")
  } else {
    message("  WARNING: mod_2 not found; inflamed label unavailable.")
  }
}

message(sprintf("  %d cells; cell type distribution:", ncol(archp)))
print(table(archp$celltype_lv3))
if ("inflamed" %in% colnames(archp@cellColData))
  print(table(archp$celltype_lv3, archp$inflamed))

# ---- Analysis 1: chromVAR deviations with blood peaks -----------------------
PBMC_names <- c("Mye_CD14.Mono", "Mye_CD16.Mono")

if (!"NSCLC_PBMCsMatrix" %in% getAvailableMatrices(archp)) {
  message("\nComputing chromVAR deviations for PBMC monocyte peaks...")
  PBMC_peaks <- lapply(PBMC_names, function(f)
    readRDS(file.path(pbmc_peaks, paste0(f, "-reproduciblePeaks.gr.rds"))))
  names(PBMC_peaks) <- PBMC_names

  archp <- addBgdPeaks(archp, force = FALSE)
  archp <- addPeakAnnotations(
    ArchRProj = archp,
    regions   = GRangesList(PBMC_peaks),
    name      = "NSCLC_PBMCs",
    force     = FALSE
  )
  archp <- addDeviationsMatrix(
    ArchRProj       = archp,
    peakAnnotation  = "NSCLC_PBMCs",
    force           = FALSE
  )
  message("  Done.")
} else {
  message("\nNSCLC_PBMCsMatrix already present.")
}

# --- 1a. Per-cell z-scores ---
message("Extracting per-cell deviation z-scores...")
devRaw  <- getMatrixFromProject(archp, useMatrix = "NSCLC_PBMCsMatrix")
dev_z   <- assays(devRaw)[[1]]           # rows = annotations, cols = cells
feat_nm <- rowData(devRaw)$name
rownames(dev_z) <- feat_nm

# Keep only z-score rows for CD14 and CD16 (ignore Mye_DC if present)
keep_rows <- intersect(PBMC_names, rownames(dev_z))
dev_z <- dev_z[keep_rows, , drop = FALSE]

df_dev <- do.call(rbind, lapply(keep_rows, function(pbmc) {
  data.frame(
    z_score  = as.numeric(dev_z[pbmc, ]),
    celltype = as.character(archp$celltype_lv3),
    inflamed = if ("inflamed" %in% colnames(archp@cellColData))
                 as.character(archp$inflamed) else "unknown",
    pbmc_ref = pbmc,
    stringsAsFactors = FALSE
  )
}))
df_dev$celltype <- factor(df_dev$celltype, levels = cell_order)
df_dev$pbmc_ref <- recode(df_dev$pbmc_ref,
  "Mye_CD14.Mono" = "Blood CD14 Mono",
  "Mye_CD16.Mono" = "Blood CD16 Mono")

# --- 1b. Grouped mean z-scores heatmap ---
message("Building grouped deviation heatmap...")
pbmcSE <- getGroupSE(
  ArchRProj = archp,
  useMatrix  = "NSCLC_PBMCsMatrix",
  groupBy    = "celltype_lv3",
  divideN    = TRUE,
  scaleTo    = NULL,
  threads    = getArchRThreads()
)
pbmc_grp <- assays(pbmcSE)[[1]]
rownames(pbmc_grp) <- rowData(pbmcSE)$name
pbmc_grp <- pbmc_grp[keep_rows, cell_order[cell_order %in% colnames(pbmc_grp)],
                     drop = FALSE]

col_dev <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#d6604d"))

ht_dev <- Heatmap(
  pbmc_grp,
  name             = "z-score",
  col              = col_dev,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  row_labels       = c("Mye_CD14.Mono" = "Blood CD14 Mono",
                       "Mye_CD16.Mono" = "Blood CD16 Mono")[rownames(pbmc_grp)],
  row_names_gp     = gpar(fontsize = 8),
  column_names_gp  = gpar(fontsize = 9),
  column_names_rot = 45,
  border           = TRUE,
  rect_gp          = gpar(col = NA),
  height           = unit(1.6, "cm"),
  cell_fun         = function(j, i, x, y, width, height, fill)
    grid.text(sprintf("%.2f", pbmc_grp[i, j]), x, y,
              gp = gpar(fontsize = 7, col = "white")),
  heatmap_legend_param = list(
    title    = "mean\nz-score",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7)
  )
)

pdf(file.path(script_dir, "plot_chromVAR_blood_peaks_heatmap.pdf"), width = 6, height = 3)
draw(ht_dev,
     column_title    = "chromVAR: accessibility at NSCLC blood monocyte peaks",
     column_title_gp = gpar(fontsize = 9, fontface = "bold"),
     padding         = unit(c(4, 4, 8, 4), "mm"))
dev.off()
message("Saved plot_chromVAR_blood_peaks_heatmap.pdf")

# --- 1c. Per-cell violin split by celltype and inflamed ----------------------
p_violin <- ggplot(df_dev, aes(x = celltype, y = z_score, fill = celltype)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey60", linetype = "dashed") +
  geom_violin(scale = "width", alpha = 0.85, colour = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white",
               colour = "grey30", linewidth = 0.4) +
  facet_wrap(~pbmc_ref, ncol = 1) +
  scale_fill_manual(values = palette_myeloid) +
  labs(x = NULL,
       y = "chromVAR z-score",
       title = "Accessibility deviation at NSCLC blood monocyte peaks\nin meso myeloid cell types") +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "none",
        panel.grid.minor   = element_blank(),
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold", size = 9))

# If inflamed label available, add violin split by inflamed vs non_inflamed
if (all(df_dev$inflamed != "unknown")) {
  df_mac <- df_dev[df_dev$celltype %in% c("Mono_CD14","TAM_MARCO","TAM_TREM2"), ]
  df_mac$inflamed <- factor(df_mac$inflamed, levels = c("inflamed","non_inflamed"))

  # Compute Wilcoxon p-values manually per celltype x pbmc_ref
  sig_df <- do.call(rbind, lapply(unique(df_mac$pbmc_ref), function(pr) {
    do.call(rbind, lapply(levels(df_mac$celltype), function(ct) {
      sub <- df_mac[df_mac$pbmc_ref == pr & df_mac$celltype == ct, ]
      if (length(unique(sub$inflamed)) < 2) return(NULL)
      pv <- wilcox.test(z_score ~ inflamed, data = sub)$p.value
      stars <- dplyr::case_when(
        pv < 0.001 ~ "***", pv < 0.01 ~ "**", pv < 0.05 ~ "*", TRUE ~ "ns")
      y_pos <- quantile(sub$z_score, 0.99, na.rm = TRUE) + 0.15
      data.frame(pbmc_ref = pr, celltype = ct, label = stars,
                 y = y_pos, stringsAsFactors = FALSE)
    }))
  }))
  sig_df$celltype <- factor(sig_df$celltype, levels = cell_order)

  p_inflamed <- ggplot(df_mac, aes(x = celltype, y = z_score, fill = inflamed)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey60", linetype = "dashed") +
    geom_violin(scale = "width", alpha = 0.8, position = position_dodge(0.9), colour = NA) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white",
                 colour = "grey30", linewidth = 0.4,
                 position = position_dodge(0.9)) +
    geom_text(data = sig_df, aes(x = celltype, y = y, label = label),
              inherit.aes = FALSE, size = 3.5, colour = "grey20") +
    facet_wrap(~pbmc_ref, ncol = 1) +
    scale_fill_manual(values = c(inflamed = "#d6604d", non_inflamed = "#4393c3"),
                      labels  = c(inflamed = "Inflamed (high mod_2)",
                                  non_inflamed = "Non-inflamed (low mod_2)")) +
    labs(x = NULL, y = "chromVAR z-score", fill = NULL,
         title = "Inflamed vs non-inflamed at blood monocyte peaks") +
    theme_bw(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
          legend.position  = "bottom",
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold", size = 9))

  pdf(file.path(script_dir, "plot_chromVAR_blood_peaks_violin.pdf"), width = 7, height = 9)
  print(p_violin / p_inflamed + plot_layout(heights = c(1.2, 1)))
  dev.off()
} else {
  pdf(file.path(script_dir, "plot_chromVAR_blood_peaks_violin.pdf"), width = 6, height = 7)
  print(p_violin)
  dev.off()
}
message("Saved plot_chromVAR_blood_peaks_violin.pdf")

# ---- Analysis 2: TF activity Spearman correlation ---------------------------
message("\nLoading blood TF activity...")
blood_tf <- read.csv(blood_tf_csv, stringsAsFactors = FALSE)
blood_tf <- blood_tf[seq_len(min(870, nrow(blood_tf))), ]  # TF rows only
rownames(blood_tf) <- blood_tf[, 1]
blood_tf <- blood_tf[, -1, drop = FALSE]

# Top 300 TFs most active in each blood monocyte type
blood_cd14 <- head(blood_tf[order(-blood_tf$Mye_CD14.Mono), ], 300)
blood_cd16 <- head(blood_tf[order(-blood_tf$Mye_CD16.Mono), ], 300)

message("Loading meso TF activity (MotifMatrix)...")
mSE <- getGroupSE(
  ArchRProj = archp,
  useMatrix  = "MotifMatrix",
  groupBy    = "celltype_lv3",
  divideN    = TRUE,
  scaleTo    = NULL,
  threads    = getArchRThreads()
)
mMat <- assays(mSE)[[1]]
rownames(mMat) <- rowData(mSE)$name
rownames(mMat) <- gsub("_.*", "", rownames(mMat))   # strip motif ID suffix

# Spearman r: mean TF profile per meso cell type vs blood monocyte ranking
cor_df <- do.call(rbind, lapply(cell_order[cell_order %in% colnames(mMat)], function(ct) {
  tfs14 <- intersect(rownames(blood_cd14), rownames(mMat))
  tfs16 <- intersect(rownames(blood_cd16), rownames(mMat))
  r_cd14 <- cor(mMat[tfs14, ct], blood_cd14[tfs14, "Mye_CD14.Mono"],
                method = "spearman", use = "complete.obs")
  r_cd16 <- cor(mMat[tfs16, ct], blood_cd16[tfs16, "Mye_CD16.Mono"],
                method = "spearman", use = "complete.obs")
  data.frame(celltype = ct, r_CD14 = r_cd14, r_CD16 = r_cd16,
             stringsAsFactors = FALSE)
}))
cor_df$celltype <- factor(cor_df$celltype, levels = cell_order)
write.csv(cor_df, file.path(script_dir, "TF_correlation_blood_results.csv"),
          row.names = FALSE, quote = FALSE)
message("TF correlation results:")
print(cor_df, row.names = FALSE)

cor_long <- pivot_longer(cor_df, cols = c(r_CD14, r_CD16),
                         names_to = "blood_ref", values_to = "r") %>%
  mutate(blood_ref = recode(blood_ref,
    r_CD14 = "vs Blood CD14 Mono",
    r_CD16 = "vs Blood CD16 Mono"))

p_cor <- ggplot(cor_long, aes(x = celltype, y = r, fill = celltype)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey50") +
  facet_wrap(~blood_ref, ncol = 1) +
  scale_fill_manual(values = palette_myeloid) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = NULL,
       y = "Spearman r (TF activity profile)",
       title = "TF activity correlation between meso myeloid\ncell types and NSCLC blood monocytes") +
  theme_bw(base_size = 11) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "none",
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold", size = 9))

pdf(file.path(script_dir, "plot_TF_correlation_blood.pdf"), width = 5, height = 6)
print(p_cor)
dev.off()
message("Saved plot_TF_correlation_blood.pdf")

# ---- Optional: inflamed vs non-inflamed TF correlation ----------------------
if ("inflamed" %in% colnames(archp@cellColData)) {
  message("\nComputing TF correlation split by inflamed/non_inflamed...")
  mSE_inf <- getGroupSE(
    ArchRProj = archp[archp$celltype_lv3 %in% c("Mono_CD14","TAM_MARCO","TAM_TREM2")],
    useMatrix  = "MotifMatrix",
    groupBy    = "inflamed",
    divideN    = TRUE, scaleTo = NULL, threads = getArchRThreads()
  )
  mMat_inf <- assays(mSE_inf)[[1]]
  rownames(mMat_inf) <- rowData(mSE_inf)$name
  rownames(mMat_inf) <- gsub("_.*", "", rownames(mMat_inf))

  cor_inf <- do.call(rbind, lapply(colnames(mMat_inf), function(gr) {
    tfs14 <- intersect(rownames(blood_cd14), rownames(mMat_inf))
    r_cd14 <- cor(mMat_inf[tfs14, gr], blood_cd14[tfs14, "Mye_CD14.Mono"],
                  method = "spearman", use = "complete.obs")
    data.frame(group = gr, r_CD14 = r_cd14, stringsAsFactors = FALSE)
  }))

  p_inf_cor <- ggplot(cor_inf, aes(x = group, y = r_CD14,
                                    fill = group)) +
    geom_col(width = 0.6, colour = "grey30", linewidth = 0.3) +
    scale_fill_manual(values = c(inflamed = "#d6604d", non_inflamed = "#4393c3")) +
    labs(x = NULL, y = "Spearman r (vs Blood CD14 Mono)",
         title = "TF activity correlation\nby inflammation score") +
    theme_bw(base_size = 11) +
    theme(legend.position  = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())

  pdf(file.path(script_dir, "plot_TF_correlation_inflamed.pdf"), width = 3.5, height = 4)
  print(p_inf_cor)
  dev.off()
  message("Saved plot_TF_correlation_inflamed.pdf")
}

message("\nAll outputs saved to ", script_dir)
