# plot_gene_scores_imputed.R  (R1_Q4)
# Plots imputed gene-score UMAPs for monocyte/macrophage markers on the
# non-harmonized UMAP of the combined meso+PBMC ArchR project.
# Run on login node: Rscript plot_gene_scores_imputed.R

suppressPackageStartupMessages({
  library(ArchR)
  library(patchwork)
  library(viridis)
})
addArchRGenome("hg38")
addArchRThreads(4)

proj_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/meso_vs_pbmc"
out_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4"

message("Loading project...")
proj <- loadArchRProject(proj_dir, showLogo = FALSE)

# ---- Immune marker list ------------------------------------------------------
immune_markers <- list(
  monocyte   = c("CD14", "FCGR3A", "CCR2", "CX3CR1", "S100A8", "S100A9",
                 "LYZ", "VCAN", "FCN1", "CDKN1C"),
  macrophage = c("CD68", "CSF1R", "MRC1", "MARCO", "TREM2", "C1QA", "C1QB",
                 "APOE", "FOLR2", "TIMD4", "LYVE1", "SELENOP"),
  TAM        = c("SPP1", "CXCL10", "CXCL9", "IL1B", "CCL2", "CCL3",
                 "FN1", "VEGFA"),
  DC         = c("CLEC9A", "ITGAX", "XCR1", "CD1C", "FCER1A", "LILRA4")
)

# Filter to genes present in GeneScoreMatrix
gene_features  <- getFeatures(proj, "GeneScoreMatrix")
immune_markers <- unlist(immune_markers)[unlist(immune_markers) %in% gene_features]
message(sprintf("Plotting %d marker genes: %s",
                length(immune_markers), paste(immune_markers, collapse = ", ")))

# ---- Impute weights ----------------------------------------------------------
if (is.null(getImputeWeights(proj))) {
  message("Adding impute weights...")
  proj <- addImputeWeights(proj)
}

# ---- plotEmbedding (non-harmonized UMAP) ------------------------------------
message("Generating feature plots...")

pdf(file.path(tempdir(), "tmp_archr.pdf"))   # active device for rasterization
p_genes <- plotEmbedding(
  ArchRProj     = proj,
  colorBy       = "GeneScoreMatrix",
  name          = immune_markers,
  embedding     = "UMAP",
  pal           = viridis::plasma(100),
  imputeWeights = getImputeWeights(proj)
)
dev.off()

# Tidy up each plot
p_genes <- lapply(seq_along(p_genes), function(i) {
  p_genes[[i]] +
    labs(title = names(p_genes)[i], x = NULL, y = NULL) +
    theme(
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      legend.position  = "none",
      plot.subtitle    = element_blank(),
      plot.caption     = element_blank()
    )
})

# ---- Save --------------------------------------------------------------------
pdf(file.path(out_dir, "umap_gene_scores_imputed.pdf"), width = 25, height = 25)
print(wrap_plots(p_genes))
dev.off()
message("Saved: umap_gene_scores_imputed.pdf")
