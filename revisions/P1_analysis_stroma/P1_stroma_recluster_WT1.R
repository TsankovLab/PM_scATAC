#!/usr/bin/env Rscript
# P1_stroma_recluster_WT1.R
#
# Revision analysis: thorough re-analysis of the normal stromal compartment
# (Fibroblasts + Mesothelium) in patient P1.
#
# Steps:
#   1. Load the main ArchR project.
#   2. Subset to P1 cells annotated as Fibroblasts or Mesothelium (celltype_lv1),
#      writing an isolated sub-project into this folder.
#   3. Regenerate dimensionality reduction (IterativeLSI), clustering and UMAP
#      de novo on the subset.
#   4. Plot WT1 (mesothelial marker) gene-activity score on the new UMAP, plus
#      UMAPs coloured by the original annotation and by the new clusters.
#
# Outputs -> P1_analysis_stroma/
#   scatac_ArchR_P1_stroma/            (subset + re-embedded ArchR project)
#   P1_stroma_umap_panels.pdf          (annotation / clusters / WT1 UMAPs)
#   P1_stroma_WT1_violin.pdf           (WT1 gene score by annotation & cluster)
#   P1_stroma_cluster_composition.csv  (cluster x celltype_lv1 cross-tab)

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
})

addArchRGenome("hg38")
addArchRThreads(4)
set.seed(1)

SRC     <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR"
OUTDIR  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/P1_analysis_stroma"
SUBPROJ <- file.path(OUTDIR, "scatac_ArchR_P1_stroma")

SAMPLE      <- "P1"
CELLTYPES   <- c("Fibroblasts", "Mesothelium")
CT_COL      <- "celltype_lv1"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load full project ─────────────────────────────────────────────────────
message("Loading main ArchR project: ", SRC)
proj <- loadArchRProject(SRC, showLogo = FALSE)
message(sprintf("  Loaded %d cells", nCells(proj)))

# ── 2. Subset to P1 stroma (Fibroblasts + Mesothelium) ───────────────────────
keep <- as.logical(as.vector(proj$Sample) == SAMPLE &
                   as.vector(proj@cellColData[[CT_COL]]) %in% CELLTYPES)
message(sprintf("  P1 %s cells: %d", paste(CELLTYPES, collapse="+"), sum(keep)))
print(table(proj@cellColData[[CT_COL]][keep]))
stopifnot(sum(keep) > 0)

cells_keep <- rownames(proj@cellColData)[keep]

message("Writing isolated sub-project: ", SUBPROJ)
sub <- subsetArchRProject(
  ArchRProj       = proj,
  cells           = cells_keep,
  outputDirectory = SUBPROJ,
  dropCells       = TRUE,
  force           = TRUE,
  logFile         = NULL
)
# carry a clean annotation column
sub$celltype <- sub@cellColData[[CT_COL]]
message(sprintf("  Sub-project: %d cells", nCells(sub)))

# ── 3. Regenerate LSI / clusters / UMAP de novo ──────────────────────────────
message("Re-running IterativeLSI on the subset ...")
sub <- addIterativeLSI(
  ArchRProj      = sub,
  useMatrix      = "TileMatrix",
  name           = "IterativeLSI",
  iterations     = 2,
  clusterParams  = list(resolution = c(0.2), sampleCells = NULL, n.start = 10),
  varFeatures    = 25000,
  dimsToUse      = 1:30,
  force          = TRUE
)

message("Clustering ...")
sub <- addClusters(
  input       = sub,
  reducedDims = "IterativeLSI",
  method      = "Seurat",
  name        = "Clusters",
  resolution  = 0.8,
  force       = TRUE
)
message("Cluster sizes:"); print(table(sub$Clusters))

message("Computing UMAP ...")
sub <- addUMAP(
  ArchRProj   = sub,
  reducedDims = "IterativeLSI",
  name        = "UMAP",
  nNeighbors  = 30,
  minDist     = 0.4,
  metric      = "cosine",
  force       = TRUE
)

message("Adding gene-score impute weights ...")
sub <- addImputeWeights(sub, reducedDims = "IterativeLSI")

saveArchRProject(ArchRProj = sub, outputDirectory = SUBPROJ, overwrite = TRUE, load = FALSE)

# Cluster x annotation composition
comp <- as.data.frame.matrix(table(sub$Clusters, sub$celltype))
write.csv(comp, file.path(OUTDIR, "P1_stroma_cluster_composition.csv"))
message("Cluster x celltype composition:"); print(comp)

# ── 4. Plots ─────────────────────────────────────────────────────────────────
message("Plotting UMAP panels ...")
p_ann <- plotEmbedding(sub, colorBy = "cellColData", name = "celltype",
                       embedding = "UMAP", size = 1.2, labelAsFactors = FALSE,
                       labelMeans = TRUE) +
  ggtitle("P1 stroma — original annotation")

p_cl  <- plotEmbedding(sub, colorBy = "cellColData", name = "Clusters",
                       embedding = "UMAP", size = 1.2, labelAsFactors = FALSE,
                       labelMeans = TRUE) +
  ggtitle("P1 stroma — de novo clusters")

p_wt1 <- plotEmbedding(sub, colorBy = "GeneScoreMatrix", name = "WT1",
                       embedding = "UMAP", size = 1.2,
                       imputeWeights = getImputeWeights(sub)) +
  ggtitle("P1 stroma — WT1 gene score")

pdf(file.path(OUTDIR, "P1_stroma_umap_panels.pdf"), width = 16, height = 5)
ggAlignPlots(p_ann, p_cl, p_wt1, type = "h")
dev.off()
message("Saved: P1_stroma_umap_panels.pdf")

# WT1 by annotation and by cluster (violin)
vln_ann <- plotGroups(sub, groupBy = "celltype", colorBy = "GeneScoreMatrix",
                      name = "WT1", plotAs = "violin", alpha = 0.5,
                      imputeWeights = getImputeWeights(sub)) +
  ggtitle("WT1 gene score by annotation")
vln_cl  <- plotGroups(sub, groupBy = "Clusters", colorBy = "GeneScoreMatrix",
                      name = "WT1", plotAs = "violin", alpha = 0.5,
                      imputeWeights = getImputeWeights(sub)) +
  ggtitle("WT1 gene score by cluster")

pdf(file.path(OUTDIR, "P1_stroma_WT1_violin.pdf"), width = 11, height = 5)
ggAlignPlots(vln_ann, vln_cl, type = "h")
dev.off()
message("Saved: P1_stroma_WT1_violin.pdf")

message("Done.")
