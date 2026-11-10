# compute_per_label_scores.R  (R1_Q4)
#
# Per-reference-label KNN scores for every myeloid scATAC cell projected onto the Greenleaf PBMC/BM reference.
# For each query cell, finds k=25 nearest reference cells in LSI space and computes the fraction of those neighbours
# belonging to each BioClassification type (soft per-label score in [0,1]).
#
# Inputs:
#   myeloid_all_lsi_projection.rds   query cells in reference LSI space (project_myeloid_on_pbmc.sh)
#   greenleaf_PBMC_BM ArchR project  reference LSI + BioClassification labels
#   myeloid scatac_ArchR project     celltype_lv3 and UMAP_H (cells kept = cells with an embedding)
#
# Output:
#   R1_Q4/output/pbmc_projection/per_label_scores.csv   (input of the three plot_* scripts)
#
# Run via: bsub < submit_compute_per_label_scores.sh

suppressPackageStartupMessages({
  library(ArchR)
  library(FNN)
})

addArchRThreads(1)
addArchRGenome("hg19")

K <- 25

# ── Paths ──────────────────────────────────────────────────────────────────────
query_proj_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR"
ref_proj_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Xmen/greenleaf_PBMC_BM"
lsi_rds        <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection/myeloid_all_lsi_projection.rds"
score_csv      <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection/per_label_scores.csv"

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

# ── Build score table ──────────────────────────────────────────────────────
score_df      <- as.data.frame(score_mat, check.names = FALSE)
score_df$cell <- rownames(score_mat)

df <- merge(emb, score_df, by = "cell")
df$celltype_lv3 <- as.character(query_archp$celltype_lv3)[match(df$cell, getCellNames(query_archp))]

set.seed(42)
df <- df[sample(nrow(df)), ]
message("Cells: ", nrow(df))

# ── Save score matrix + celltype_lv3 ───────────────────────────────────────────
write.csv(df[, c("cell", "celltype_lv3", all_labels)], score_csv,
          row.names = FALSE, quote = FALSE)
message("Score CSV: ", score_csv)
