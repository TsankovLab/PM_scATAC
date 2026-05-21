# plot_peak_similarity.R  (R1_Q4)
#
# Peak-level overlap analysis to establish epigenomic similarity between
# NSCLC blood CD14 monocytes and meso myeloid cell types.
#
# Comparisons:
#   Blood  : Mye_CD14.Mono, Mye_CD16.Mono  (NSCLC PBMC project)
#   Meso   : Mono_CD14 (=cnmf_cluster_3, high inflammation / MDM),
#             Mono_CD16 (=cnmf_cluster_1), TAM_CXCLs (=cnmf_cluster_4),
#             TAM_TREM2 (=cnmf_cluster_5), TAM_MARCO (=cnmf_cluster_2),
#             TAM_interstitial (=cnmf_cluster_6)
#   cnmf assignments from scatac_myeloid_ArchR.R lines 228-235
#
# Outputs:
#   plot_peak_jaccard_heatmap.pdf  -- pairwise Jaccard similarity heatmap
#   plot_peak_overlap_barplot.pdf  -- % blood CD14 Mono peaks per meso type
#
# Runs on login node (<2 min, no fragment / ArchR project loading)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(dplyr)
})

script_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4"
pbmc_peaks  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/pbmc_myeloid/PeakCalls"
meso_peaks  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/PeakCalls"

# ---- peak set definitions ----------------------------------------------------
# Display name -> peak file basename (without -reproduciblePeaks.gr.rds)
blood_map <- c(
  "Blood CD14 Mono"  = "Mye_CD14.Mono",
  "Blood CD16 Mono"  = "Mye_CD16.Mono"
)

# Exact celltype_lv3 named peaksets
meso_map <- c(
  "Meso Mono_CD14"         = "Mono_CD14",
  "Meso Mono_CD16"         = "Mono_CD16",
  "Meso TAM_CXCLs"         = "TAM_CXCLs",
  "Meso TAM_interstitial"  = "TAM_interstitial",
  "Meso TAM_MARCO"         = "TAM_MARCO",
  "Meso TAM_TREM2"         = "TAM_TREM2",
  "Meso cDCs"              = "cDCs"
)

# ---- load peaks --------------------------------------------------------------
message("Loading peak sets...")
load_peaks <- function(dir, basenames) {
  lapply(basenames, function(f) {
    path <- file.path(dir, paste0(f, "-reproduciblePeaks.gr.rds"))
    if (!file.exists(path)) stop("Peak file not found: ", path)
    readRDS(path)
  })
}

blood_list <- load_peaks(pbmc_peaks, blood_map)
names(blood_list) <- names(blood_map)

meso_list  <- load_peaks(meso_peaks, meso_map)
names(meso_list) <- names(meso_map)

peak_sets <- c(blood_list, meso_list)

for (nm in names(peak_sets))
  message(sprintf("  %-28s  %d peaks", nm, length(peak_sets[[nm]])))

# ---- Jaccard similarity ------------------------------------------------------
jaccard <- function(a, b) {
  n_ab <- sum(overlapsAny(a, b, ignore.strand = TRUE))
  n_ab / (length(a) + length(b) - n_ab)
}

message("\nComputing pairwise Jaccard similarity (", length(peak_sets), "x", length(peak_sets), ")...")
n <- length(peak_sets)
sim_mat <- matrix(0, n, n, dimnames = list(names(peak_sets), names(peak_sets)))
for (i in seq_len(n))
  for (j in seq_len(n))
    sim_mat[i, j] <- jaccard(peak_sets[[i]], peak_sets[[j]])

# Print diagonal-free summary
message("Jaccard similarity with Blood CD14 Mono (sorted):")
jac_cd14 <- sort(sim_mat["Blood CD14 Mono", names(peak_sets) != "Blood CD14 Mono"],
                 decreasing = TRUE)
print(round(jac_cd14, 4))

# ---- Annotations for heatmap -------------------------------------------------
source_ann <- ifelse(grepl("^Blood", rownames(sim_mat)), "Blood (NSCLC)", "Mesothelioma")

celltype_ann <- dplyr::case_when(
  grepl("Mono_CD14",                        rownames(sim_mat)) ~ "CD14 Mono / MDM",
  grepl("Mono_CD16|Blood CD16",             rownames(sim_mat)) ~ "CD16 Mono",
  grepl("TAM_MARCO|TAM_TREM2|TAM_inter",   rownames(sim_mat)) ~ "Resident Macrophage",
  grepl("TAM_CXCLs",                        rownames(sim_mat)) ~ "Inflammatory TAM",
  TRUE                                                          ~ "Other"
)

source_col <- c("Blood (NSCLC)" = "#1f78b4", "Mesothelioma" = "#e31a1c")
ct_col <- c(
  "CD14 Mono / MDM"   = "#ff7f00",
  "CD16 Mono"         = "#fdbf6f",
  "Resident Macrophage" = "#33a02c",
  "Inflammatory TAM"  = "#b2df8a",
  "Other"             = "#cab2d6"
)

ann_top  <- HeatmapAnnotation(
  Source    = source_ann,
  Cell_type = celltype_ann,
  col       = list(Source = source_col, Cell_type = ct_col),
  annotation_name_gp = gpar(fontsize = 7.5),
  show_legend = TRUE
)
ann_left <- rowAnnotation(
  Source    = source_ann,
  Cell_type = celltype_ann,
  col       = list(Source = source_col, Cell_type = ct_col),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

max_off <- max(sim_mat[upper.tri(sim_mat)])
col_fun <- colorRamp2(c(0, max_off / 2, max_off),
                      c("#f7fbff", "#6baed6", "#08306b"))

ht <- Heatmap(
  sim_mat,
  name             = "Jaccard",
  col              = col_fun,
  top_annotation   = ann_top,
  left_annotation  = ann_left,
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  clustering_distance_rows    = "euclidean",
  clustering_distance_columns = "euclidean",
  row_names_gp     = gpar(fontsize = 8),
  column_names_gp  = gpar(fontsize = 8),
  column_names_rot = 45,
  border           = TRUE,
  rect_gp          = gpar(col = NA),
  cell_fun         = function(j, i, x, y, width, height, fill) {
    if (abs(sim_mat[i, j]) > 0.03)
      grid.text(sprintf("%.3f", sim_mat[i, j]), x, y,
                gp = gpar(fontsize = 5.5, col = "white"))
  },
  heatmap_legend_param = list(
    title      = "Jaccard\nsimilarity",
    title_gp   = gpar(fontsize = 8, fontface = "bold"),
    labels_gp  = gpar(fontsize = 7),
    legend_height = unit(3, "cm")
  )
)

pdf(file.path(script_dir, "plot_peak_jaccard_heatmap.pdf"), width = 7.5, height = 6.5)
draw(ht,
     column_title     = "Pairwise peak overlap (Jaccard) — blood vs meso myeloid",
     column_title_gp  = gpar(fontsize = 10, fontface = "bold"),
     padding          = unit(c(4, 4, 8, 4), "mm"))
dev.off()
message("Saved plot_peak_jaccard_heatmap.pdf")

# ---- Overlap fraction barplot ------------------------------------------------
blood_cd14 <- peak_sets[["Blood CD14 Mono"]]
frac_df <- data.frame(
  celltype = names(peak_sets),
  n_peaks  = sapply(peak_sets, length),
  n_overlap = sapply(peak_sets, function(x)
    sum(overlapsAny(blood_cd14, x, ignore.strand = TRUE))),
  source   = source_ann,
  stringsAsFactors = FALSE
)
frac_df$frac_pct <- frac_df$n_overlap / length(blood_cd14) * 100

frac_df$celltype <- factor(frac_df$celltype,
  levels = frac_df$celltype[order(frac_df$frac_pct, decreasing = TRUE)])

ct_fill <- c(
  "Blood CD14 Mono"        = "#1f78b4",
  "Blood CD16 Mono"        = "#a6cee3",
  "Meso Mono_CD14"         = "#d6604d",
  "Meso Mono_CD16"         = "#f4a582",
  "Meso TAM_CXCLs"         = "#4393c3",
  "Meso TAM_interstitial"  = "#99d8c9",
  "Meso TAM_MARCO"         = "#2ca25f",
  "Meso TAM_TREM2"         = "#74c476",
  "Meso cDCs"              = "#8e7cb8"
)

p_bar <- ggplot(frac_df, aes(x = celltype, y = frac_pct, fill = celltype)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = ct_fill) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(frac_df$frac_pct) * 1.08)) +
  labs(
    x     = NULL,
    y     = "Blood CD14 Mono peaks accessible (%)",
    title = "Fraction of blood CD14 monocyte peaks\noverlapping each myeloid cell type's peak set"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

pdf(file.path(script_dir, "plot_peak_overlap_barplot.pdf"), width = 5.5, height = 4.5)
print(p_bar)
dev.off()
message("Saved plot_peak_overlap_barplot.pdf")

# ---- Summary table -----------------------------------------------------------
message("\nOverlap summary table:")
print(frac_df[order(-frac_df$frac_pct),
              c("celltype","n_peaks","n_overlap","frac_pct")], row.names = FALSE)

write.csv(frac_df, file.path(script_dir, "peak_overlap_summary.csv"),
          row.names = FALSE, quote = FALSE)
message("Done.")
