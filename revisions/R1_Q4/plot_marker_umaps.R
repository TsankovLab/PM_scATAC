# plot_marker_umaps.R  (R1_Q4)
#
# 1. Volcano plots for CD14_Mono, CD16_Mono, DC DA peaks (fixed rendering)
# 2. Harmony UMAP coloured by gene-activity scores for key myeloid markers
#
# Run on login node: Rscript plot_marker_umaps.R

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(patchwork)
  library(SummarizedExperiment)
})
addArchRGenome("hg38")
addArchRThreads(4)

proj_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/meso_vs_pbmc"
out_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4"

message("Loading project...")
proj <- loadArchRProject(proj_dir, showLogo = FALSE)

# ---- 1. Volcano plots --------------------------------------------------------
message("Generating volcano plots...")

marker_targets <- c("CD14_Mono", "CD16_Mono", "DC")

volcano_plot <- function(mkr, tgt, fdr_cut = 0.05, lfc_cut = 1) {
  col <- colnames(mkr)[1]
  df <- data.frame(
    Log2FC = as.numeric(assays(mkr)[["Log2FC"]][[col]]),
    FDR    = as.numeric(assays(mkr)[["FDR"]][[col]]),
    stringsAsFactors = FALSE
  )
  df$neglog10FDR <- -log10(pmax(df$FDR, 1e-300))
  df$sig <- df$FDR < fdr_cut & df$Log2FC > lfc_cut
  df$col <- ifelse(df$sig, "#d6604d", "#aaaaaa")
  n_sig  <- sum(df$sig, na.rm = TRUE)

  ggplot(df, aes(x = Log2FC, y = neglog10FDR)) +
    geom_point(colour = df$col, size = 0.4, alpha = 0.5) +
    geom_vline(xintercept = lfc_cut,  linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("n=", n_sig), size = 3.5, colour = "#d6604d") +
    labs(
      x     = expression(log[2]~FC),
      y     = expression(-log[10]~FDR),
      title = paste0(tgt, " vs rest")
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 10, face = "bold"))
}

vol_list <- lapply(marker_targets, function(tgt) {
  rds <- file.path(proj_dir, paste0("DAP_", tgt, "_vs_rest.rds"))
  if (!file.exists(rds)) { message("  Missing: ", rds); return(NULL) }
  mkr <- readRDS(rds)
  volcano_plot(mkr, tgt)
})
names(vol_list) <- marker_targets

# Individual PDFs
for (tgt in marker_targets) {
  if (!is.null(vol_list[[tgt]]))
    ggsave(file.path(out_dir, paste0("volcano_", tgt, ".pdf")),
           vol_list[[tgt]], width = 4.5, height = 4)
}

# Combined panel
vol_panel <- wrap_plots(vol_list[!sapply(vol_list, is.null)], nrow = 1)
ggsave(file.path(out_dir, "volcano_all_targets.pdf"), vol_panel,
       width = 13, height = 4)
message("Saved volcano plots.")

# ---- 2. Gene-score UMAPs -----------------------------------------------------
message("Extracting gene scores...")

markers <- c(
  "CD14", "FCGR3A",           # monocyte
  "TREM2", "MARCO", "MRC1",   # resident / TREM2 macrophage
  "SPP1", "CXCL10",           # inflammatory TAM
  "CD68", "CSF1R",            # pan-macrophage
  "ITGAX", "CLEC9A"           # DC
)

# Pull GeneScoreMatrix for these genes
gsm <- getMatrixFromProject(
  ArchRProj  = proj,
  useMatrix  = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose    = FALSE
)
gene_names <- rowData(gsm)$name
keep_idx   <- which(gene_names %in% markers)
found      <- gene_names[keep_idx]
message("  Found: ", paste(found, collapse = ", "))

score_mat <- as.matrix(assay(gsm)[keep_idx, , drop = FALSE])
rownames(score_mat) <- gene_names[keep_idx]

# UMAP coordinates (Harmony)
emb <- getEmbedding(proj, embedding = "UMAP_H", returnDF = TRUE)
colnames(emb) <- c("UMAP1", "UMAP2")
emb$cell <- rownames(emb)

# Merge
score_df <- as.data.frame(t(score_mat))
score_df$cell <- colnames(score_mat)
plot_df <- merge(emb, score_df, by = "cell")

# Shuffle rows for unbiased overplotting
set.seed(42)
plot_df <- plot_df[sample(nrow(plot_df)), ]

# Plot function
gene_umap <- function(df, gene, pt_size = 0.1) {
  if (!gene %in% colnames(df)) return(NULL)
  vals <- df[[gene]]
  # Winsorise top 1% for display
  cap <- quantile(vals, 0.99, na.rm = TRUE)
  vals_capped <- pmin(vals, cap)
  df$score <- vals_capped

  ggplot(df, aes(x = UMAP1, y = UMAP2, colour = score)) +
    geom_point(size = pt_size, stroke = 0) +
    scale_colour_gradientn(
      colours = c("#f7f7f7", "#fee090", "#fc8d59", "#d73027"),
      name    = "Gene\nscore"
    ) +
    labs(title = gene) +
    theme_classic(base_size = 10) +
    theme(
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.line  = element_blank(),
      plot.title = element_text(size = 10, face = "bold.italic"),
      legend.key.height = unit(8, "mm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7)
    )
}

gene_plots <- lapply(found, function(g) gene_umap(plot_df, g))
names(gene_plots) <- found
gene_plots <- gene_plots[!sapply(gene_plots, is.null)]

# Grid layout
n_genes <- length(gene_plots)
ncols   <- 4
nrows   <- ceiling(n_genes / ncols)
p_grid  <- wrap_plots(gene_plots, ncol = ncols)

ggsave(file.path(out_dir, "umap_gene_scores.pdf"), p_grid,
       width = ncols * 3.5, height = nrows * 3.2)
message("Saved gene-score UMAPs: ", file.path(out_dir, "umap_gene_scores.pdf"))

message("\nDone.")
