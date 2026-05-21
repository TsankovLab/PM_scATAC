# plot_meso_pbmc_umaps.R  (R1_Q4)
# Loads the combined meso+PBMC ArchR project and generates UMAP plots
# using ggplot2 directly (bypasses ArchR plotEmbedding rendering issues).
# Run on login node: Rscript plot_meso_pbmc_umaps.R

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(patchwork)
})
addArchRGenome("hg38")
addArchRThreads(1)

proj_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/meso_vs_pbmc"
plot_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4"

message("Loading project...")
proj <- loadArchRProject(proj_dir, showLogo = FALSE)

# ---- Extract UMAP coordinates ------------------------------------------------
get_umap_df <- function(proj, embedding) {
  emb <- getEmbedding(proj, embedding = embedding, returnDF = TRUE)
  colnames(emb) <- c("UMAP1", "UMAP2")
  emb$cell <- rownames(emb)
  emb$tissue          <- proj$tissue[match(emb$cell, rownames(proj@cellColData))]
  emb$celltype_unified <- proj$celltype_unified[match(emb$cell, rownames(proj@cellColData))]
  emb$Clusters_H      <- proj$Clusters_H[match(emb$cell, rownames(proj@cellColData))]
  emb
}

df_h   <- get_umap_df(proj, "UMAP_H")
df_lsi <- get_umap_df(proj, "UMAP")

# ---- Palettes ----------------------------------------------------------------
pal_tissue <- c(PBMC = "#1f78b4", Mesothelioma = "#e31a1c")

pal_ct <- c(
  CD14_Mono        = "#d6604d",
  CD16_Mono        = "#f4a582",
  DC               = "#8e7cb8",
  TAM_CXCLs        = "#4393c3",
  TAM_interstitial = "#99d8c9",
  TAM_MARCO        = "#2ca25f",
  TAM_TREM2        = "#74c476",
  Unknown          = "#cccccc"
)

# ---- Plot helper -------------------------------------------------------------
umap_plot <- function(df, colour_col, pal, title, pt_size = 0.15, alpha = 0.6) {
  # Shuffle so neither group occludes the other
  df <- df[sample(nrow(df)), ]
  aes_col <- if (is.null(pal)) aes(colour = .data[[colour_col]])
             else               aes(colour = .data[[colour_col]])
  p <- ggplot(df, aes(x = UMAP1, y = UMAP2, colour = .data[[colour_col]])) +
    geom_point(size = pt_size, alpha = alpha, stroke = 0) +
    labs(title = title, colour = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.line  = element_blank(),
      plot.title = element_text(size = 10, face = "bold"),
      legend.key.size = unit(4, "mm"),
      legend.text = element_text(size = 7)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
  if (!is.null(pal))
    p <- p + scale_colour_manual(values = pal, na.value = "#cccccc")
  p
}

# ---- Harmony UMAP ------------------------------------------------------------
p1 <- umap_plot(df_h,   "tissue",           pal_tissue, "Tissue (Harmony UMAP)")
p2 <- umap_plot(df_h,   "celltype_unified",  pal_ct,    "Cell type (Harmony UMAP)")
p3 <- umap_plot(df_h,   "Clusters_H",        NULL,      "Clusters (Harmony UMAP)")

# ---- LSI UMAP ----------------------------------------------------------------
p4 <- umap_plot(df_lsi, "tissue",            pal_tissue, "Tissue (LSI UMAP)")
p5 <- umap_plot(df_lsi, "celltype_unified",  pal_ct,     "Cell type (LSI UMAP)")

# ---- Save individual PDFs ----------------------------------------------------
ggsave(file.path(plot_dir, "umap_harmony_tissue.pdf"),    p1, width = 5.5, height = 4.5)
ggsave(file.path(plot_dir, "umap_harmony_celltype.pdf"),  p2, width = 6,   height = 4.5)
ggsave(file.path(plot_dir, "umap_harmony_clusters.pdf"),  p3, width = 5.5, height = 4.5)
ggsave(file.path(plot_dir, "umap_lsi_tissue.pdf"),        p4, width = 5.5, height = 4.5)
ggsave(file.path(plot_dir, "umap_lsi_celltype.pdf"),      p5, width = 6,   height = 4.5)

# ---- Combined panel ----------------------------------------------------------
p_panel <- (p1 | p2 | p3) / (p4 | p5 | plot_spacer()) +
  plot_annotation(
    title   = "Meso + PBMC myeloid — combined ArchR project",
    theme   = theme(plot.title = element_text(size = 12, face = "bold"))
  )
ggsave(file.path(plot_dir, "umap_combined_panel.pdf"), p_panel,
       width = 15, height = 9)

message("Saved UMAPs to ", plot_dir)
