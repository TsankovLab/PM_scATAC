# plot_per_label_confidence_umap.R  (R1_Q4)
#
# Computes per-reference-label KNN confidence scores for each query cell
# and plots one UMAP panel per BioClassification label on the myeloid UMAP_H.
#
# For each query cell, finds k=25 nearest reference cells in LSI space and
# computes the fraction of those neighbors belonging to each BioClassification
# type — giving a soft per-label score in [0,1] for every cell.
#
# Inputs:
#   myeloid_all_lsi_projection.rds   query cells in reference LSI space
#   greenleaf_PBMC_BM ArchR project  reference LSI + BioClassification labels
#   myeloid scatac_ArchR project     UMAP_H coordinates
#
# Output:
#   R1_Q4/output/pbmc_projection/umap_per_label_confidence.pdf
#
# Run via: bsub < submit_plot_per_label_confidence.sh

suppressPackageStartupMessages({
  library(ArchR)
  library(FNN)
  library(ggplot2)
  library(patchwork)
})

addArchRThreads(1)
addArchRGenome("hg19")

K <- 25

# ── Paths ──────────────────────────────────────────────────────────────────────
query_proj_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR"
ref_proj_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Xmen/greenleaf_PBMC_BM"
lsi_rds        <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection/myeloid_all_lsi_projection.rds"
out_file       <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection/umap_per_label_confidence.pdf"

# ── cnmf_cluster → celltype_lv3 map ───────────────────────────────────────────
ct_map <- c(
  cnmf_cluster_1 = "Mono_CD16",   cnmf_cluster_2 = "TAM_MARCO",
  cnmf_cluster_3 = "Mono_CD14",   cnmf_cluster_4 = "TAM_CXCLs",
  cnmf_cluster_5 = "TAM_TREM2",   cnmf_cluster_6 = "TAM_interstitial",
  cnmf_cluster_7 = "cDCs"
)

# ── Load query LSI projection ──────────────────────────────────────────────────
message("Loading query LSI projection...")
lsi_q <- as.matrix(readRDS(lsi_rds))   # cells × dims
message("  Query cells: ", nrow(lsi_q), "  dims: ", ncol(lsi_q))

# ── Load reference LSI + BioClassification ─────────────────────────────────────
message("Loading reference ArchR project...")
ref_archp  <- loadArchRProject(ref_proj_dir, showLogo = FALSE)
lsi_ref_obj <- getReducedDims(ref_archp, reducedDims = "IterativeLSI", returnMatrix = FALSE)
svd_ref_all <- as.matrix(lsi_ref_obj$matSVD)

ref_meta   <- as.data.frame(getCellColData(ref_archp, select = "BioClassification"))
ref_labels <- ref_meta[["BioClassification"]]
names(ref_labels) <- rownames(ref_meta)

common_ref        <- intersect(rownames(svd_ref_all), names(ref_labels))
svd_ref           <- svd_ref_all[common_ref, , drop = FALSE]
ref_labels_aligned <- ref_labels[common_ref]
message("  Reference cells used: ", length(common_ref))

# ── KNN: query → reference ────────────────────────────────────────────────────
message("Computing KNN (k=", K, ")...")
knn_res <- FNN::get.knnx(data = svd_ref, query = lsi_q, k = K)
# knn_res$nn.index: nQuery × K  (indices into svd_ref rows)

# ── Per-label fraction matrix ──────────────────────────────────────────────────
all_labels <- sort(unique(ref_labels_aligned))
message("BioClassification labels: ", paste(all_labels, collapse = ", "))

score_mat <- matrix(0.0,
  nrow = nrow(lsi_q), ncol = length(all_labels),
  dimnames = list(rownames(lsi_q), all_labels)
)

for (i in seq_len(nrow(lsi_q))) {
  nb_labels <- ref_labels_aligned[knn_res$nn.index[i, ]]
  tab <- table(nb_labels)
  score_mat[i, names(tab)] <- as.numeric(tab) / K
}

message("Per-label score matrix: ", nrow(score_mat), " cells x ", ncol(score_mat), " labels")

# ── Load query UMAP_H ──────────────────────────────────────────────────────────
message("Loading query UMAP_H...")
query_archp <- loadArchRProject(query_proj_dir, showLogo = FALSE)
if (!"celltype_lv3" %in% colnames(query_archp@cellColData))
  query_archp$celltype_lv3 <- ct_map[as.character(query_archp$cnmf_cluster)]

emb <- getEmbedding(query_archp, embedding = "UMAP_H", returnDF = TRUE)
colnames(emb) <- c("UMAP1", "UMAP2")
emb$cell <- rownames(emb)

# ── Build plot data frame ──────────────────────────────────────────────────────
score_df      <- as.data.frame(score_mat, check.names = FALSE)
score_df$cell <- rownames(score_mat)

df <- merge(emb, score_df, by = "cell")
df$celltype_lv3 <- as.character(query_archp$celltype_lv3)[match(df$cell, getCellNames(query_archp))]

set.seed(42)
df <- df[sample(nrow(df)), ]
message("Final plot cells: ", nrow(df))

# ── Shared theme ───────────────────────────────────────────────────────────────
umap_theme <- function() {
  theme_classic(base_size = 9) +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    legend.key.height = unit(5, "mm"),
    legend.text  = element_text(size = 7),
    legend.title = element_text(size = 7),
    plot.title   = element_text(size = 9, face = "bold")
  )
}

pt <- 0.15

# ── Reference annotation panel ─────────────────────────────────────────────────
ct_cols <- c(
  Mono_CD16 = "#E41A1C", Mono_CD14 = "#FF7F00", TAM_MARCO = "#4DAF4A",
  TAM_CXCLs = "#984EA3", TAM_TREM2 = "#377EB8", TAM_interstitial = "#A65628",
  cDCs      = "#F781BF"
)
p_ref <- ggplot(df[order(df$celltype_lv3), ], aes(UMAP1, UMAP2, colour = celltype_lv3)) +
  geom_point(size = pt, stroke = 0) +
  scale_colour_manual(values = ct_cols, name = "celltype_lv3") +
  guides(colour = guide_legend(override.aes = list(size = 2.5))) +
  labs(title = "celltype_lv3 (reference)") +
  umap_theme()

# ── Build per-label panels ──────────────────────────────────────────────────────
label_panels <- lapply(seq_along(all_labels), function(i) {
  lbl <- all_labels[i]
  df$score <- as.numeric(df[[lbl]])
  df <- df[order(df$score), ]
  ggplot(df, aes(UMAP1, UMAP2, colour = score)) +
    geom_point(size = pt, stroke = 0) +
    scale_colour_gradientn(
      colours = c("#f0f0f0", "#fdae61", "#d73027"),
      limits  = c(0, 1),
      name    = "KNN\nfraction"
    ) +
    labs(title = lbl) +
    umap_theme()
})

# ── Layout: reference + one panel per label ─────────────────────────────────────
all_panels <- c(list(p_ref), label_panels)
ncols      <- 4
nrows      <- ceiling(length(all_panels) / ncols)

combined <- wrap_plots(all_panels, ncol = ncols)

ggsave(out_file, combined,
       width  = ncols * 3.5,
       height = nrows * 3.2)
message("Saved: ", out_file)

# ── Save score matrix + celltype_lv3 as CSV for downstream use ────────────────
score_csv <- sub("umap_per_label_confidence.pdf", "per_label_scores.csv", out_file)
write.csv(df[, c("cell", "celltype_lv3", all_labels)], score_csv,
          row.names = FALSE, quote = FALSE)
message("Score CSV: ", score_csv)
