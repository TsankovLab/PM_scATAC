# plot_pbmc_projection_umap.R  (R1_Q4)
#
# Maps Greenleaf PBMC/BM BioClassification labels and prediction confidence
# onto the myeloid query UMAP_H embedding.
#
# Plots (one PDF, 4 panels):
#   1. UMAP coloured by celltype_lv3  (existing myeloid annotation, reference)
#   2. UMAP coloured by pred_label    (KNN BioClassification from reference)
#   3. UMAP coloured by pred_conf     (KNN confidence 0–1, continuous)
#   4. UMAP coloured by mean_knn_dist (KNN distance, QC)
#
# Input  : R1_Q4/output/pbmc_projection/myeloid_all_knn_labels.csv
# Output : R1_Q4/output/pbmc_projection/umap_pbmc_projection.pdf
# Run via: bsub < submit_plot_pbmc_projection.sh

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(patchwork)
})

addArchRThreads(1)
addArchRGenome("hg19")

# ── Paths ──────────────────────────────────────────────────────────────────────
proj_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR"
knn_file  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection/myeloid_all_knn_labels.csv"
out_file  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection/umap_pbmc_projection.pdf"

# ── cnmf_cluster → celltype_lv3 map ───────────────────────────────────────────
ct_map <- c(
  cnmf_cluster_1 = "Mono_CD16",
  cnmf_cluster_2 = "TAM_MARCO",
  cnmf_cluster_3 = "Mono_CD14",
  cnmf_cluster_4 = "TAM_CXCLs",
  cnmf_cluster_5 = "TAM_TREM2",
  cnmf_cluster_6 = "TAM_interstitial",
  cnmf_cluster_7 = "cDCs"
)

# ── Load project ───────────────────────────────────────────────────────────────
message("Loading ArchR project...")
archp <- loadArchRProject(proj_dir, showLogo = FALSE)

if (!"celltype_lv3" %in% colnames(archp@cellColData))
  archp$celltype_lv3 <- ct_map[as.character(archp$cnmf_cluster)]

# ── UMAP coordinates ───────────────────────────────────────────────────────────
emb <- getEmbedding(archp, embedding = "UMAP_H", returnDF = TRUE)
colnames(emb) <- c("UMAP1", "UMAP2")
emb$cell <- rownames(emb)

# ── KNN labels ─────────────────────────────────────────────────────────────────
knn <- read.csv(knn_file)

# ── Merge ──────────────────────────────────────────────────────────────────────
meta <- data.frame(
  cell          = getCellNames(archp),
  celltype_lv3  = as.character(archp$celltype_lv3),
  stringsAsFactors = FALSE
)
df <- merge(emb, meta, by = "cell")
df <- merge(df,  knn,  by = "cell", all.x = TRUE)

set.seed(42)
df <- df[sample(nrow(df)), ]

message("Cells in plot: ", nrow(df))

# ── Shared theme ───────────────────────────────────────────────────────────────
umap_theme <- function() {
  theme_classic(base_size = 10) +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    legend.key.height = unit(5, "mm"),
    legend.text  = element_text(size = 7),
    legend.title = element_text(size = 8)
  )
}

pt <- 0.15

# ── Panel 1: celltype_lv3 ──────────────────────────────────────────────────────
ct_cols <- c(
  Mono_CD16        = "#E41A1C",
  Mono_CD14        = "#FF7F00",
  TAM_MARCO        = "#4DAF4A",
  TAM_CXCLs        = "#984EA3",
  TAM_TREM2        = "#377EB8",
  TAM_interstitial = "#A65628",
  cDCs             = "#F781BF"
)

p1 <- ggplot(df, aes(UMAP1, UMAP2, colour = celltype_lv3)) +
  geom_point(size = pt, stroke = 0) +
  scale_colour_manual(values = ct_cols, name = "celltype_lv3") +
  guides(colour = guide_legend(override.aes = list(size = 2.5))) +
  labs(title = "Myeloid annotation (celltype_lv3)") +
  umap_theme()

# ── Panel 2: pred_label ────────────────────────────────────────────────────────
pred_labels <- sort(unique(df$pred_label[!is.na(df$pred_label)]))
n_labels    <- length(pred_labels)
label_cols  <- setNames(
  colorRampPalette(c(
    "#1F78B4","#33A02C","#E31A1C","#FF7F00","#6A3D9A",
    "#B15928","#A6CEE3","#B2DF8A","#FB9A99","#FDBF6F",
    "#CAB2D6","#808080"
  ))(n_labels),
  pred_labels
)

p2 <- ggplot(df, aes(UMAP1, UMAP2, colour = pred_label)) +
  geom_point(size = pt, stroke = 0) +
  scale_colour_manual(values = label_cols, name = "BioClassification\n(pred_label)", na.value = "grey80") +
  guides(colour = guide_legend(override.aes = list(size = 2.5))) +
  labs(title = "Predicted BioClassification (KNN)") +
  umap_theme()

# ── Panel 3: pred_conf ─────────────────────────────────────────────────────────
p3 <- ggplot(df, aes(UMAP1, UMAP2, colour = pred_conf)) +
  geom_point(size = pt, stroke = 0) +
  scale_colour_gradientn(
    colours  = c("#f7f7f7", "#fdae61", "#d73027"),
    limits   = c(0, 1),
    name     = "Confidence\n(pred_conf)"
  ) +
  labs(title = "Prediction confidence") +
  umap_theme()

# ── Panel 4: mean_knn_dist ─────────────────────────────────────────────────────
dist_cap <- quantile(df$mean_knn_dist, 0.99, na.rm = TRUE)
df$dist_capped <- pmin(df$mean_knn_dist, dist_cap)

p4 <- ggplot(df, aes(UMAP1, UMAP2, colour = dist_capped)) +
  geom_point(size = pt, stroke = 0) +
  scale_colour_gradientn(
    colours = c("#2166ac", "#f7f7f7", "#d73027"),
    name    = "Mean KNN\ndist (capped\n99th pct)"
  ) +
  labs(title = "Mean KNN distance (QC)") +
  umap_theme()

# ── Save ───────────────────────────────────────────────────────────────────────
combined <- (p1 | p2) / (p3 | p4)

ggsave(out_file, combined, width = 14, height = 11)
message("Saved: ", out_file)
