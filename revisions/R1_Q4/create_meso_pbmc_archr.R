# create_meso_pbmc_archr.R  (R1_Q4)
#
# Combines NSCLC blood (PBMC) and mesothelioma myeloid ArchR Arrow files into
# a single project, runs dimensionality reduction + Harmony-corrected UMAP,
# and identifies DA peaks for CD14 Monocytes, CD16 Monocytes, and DCs.
#
# PBMC cells: onlyATAC (multiome_ann2 == "onlyATAC") from pbmc_myeloid project
# Meso cells: all cells from myeloid_cells/scatac_ArchR project
#
# Harmony correction: Sample + tissue (PBMC vs Mesothelioma)
# Marker analysis  : CD14_Mono, CD16_Mono, DC  (one-vs-rest Wilcoxon)
#
# Output: /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/
#          myeloid_cells/meso_vs_pbmc/
#
# Requires: 8 CPUs, 128 GB RAM, 6 h  -->  submit via submit_create_meso_pbmc.sh

suppressPackageStartupMessages({
  library(ArchR)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
addArchRGenome("hg38")
addArchRThreads(threads = 8)

set.seed(1234)

# ---- Paths -------------------------------------------------------------------
pbmc_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/pbmc_myeloid"
meso_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR"
out_dir    <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/meso_vs_pbmc"
plot_dir   <- file.path(out_dir, "Plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load source projects and collect metadata ----------------------------
message("Loading PBMC project...")
pbmc <- loadArchRProject(pbmc_dir, showLogo = FALSE)

message("Loading meso myeloid project...")
meso <- loadArchRProject(meso_dir, showLogo = FALSE)

# PBMC: filter to onlyATAC cells
onlyATAC <- rownames(pbmc@cellColData)[pbmc$multiome_ann2 == "onlyATAC"]
message(sprintf("  PBMC onlyATAC cells: %d", length(onlyATAC)))

# Meso: restore celltype_lv3 if absent
ct_map <- c(cnmf_cluster_1 = "Mono_CD16", cnmf_cluster_2 = "TAM_MARCO",
            cnmf_cluster_3 = "Mono_CD14", cnmf_cluster_4 = "TAM_CXCLs",
            cnmf_cluster_5 = "TAM_TREM2", cnmf_cluster_6 = "TAM_interstitial",
            cnmf_cluster_7 = "cDCs")
if (!"celltype_lv3" %in% colnames(meso@cellColData))
  meso$celltype_lv3 <- ct_map[as.character(meso$cnmf_cluster)]
meso_cells <- rownames(meso@cellColData)
message(sprintf("  Meso cells: %d", length(meso_cells)))

# Save metadata tables for later transfer
pbmc_meta <- as.data.frame(pbmc@cellColData[onlyATAC,
               c("Sample","celltype_simplified","tissue"), drop = FALSE])
meso_meta <- as.data.frame(meso@cellColData[meso_cells,
               c("Sample","celltype_lv3"), drop = FALSE])
meso_meta$tissue <- "Mesothelioma"

# ---- 2. Symlink Arrow files to output dir ------------------------------------
arrow_out <- file.path(out_dir, "ArrowFiles")
dir.create(arrow_out, showWarnings = FALSE)

# Only the 11 PBMC samples that contain onlyATAC cells
onlyATAC_samples <- c("Lu895P","Lu898P","Lu929P","Lu941P","Lu954P",
                      "Lu934P","Lu936P","Lu912P","Lu931P","Lu979P","Lu927P")
pbmc_arrows <- file.path(pbmc_dir, "ArrowFiles",
                         paste0(onlyATAC_samples, ".arrow"))
meso_arrows <- list.files(file.path(meso_dir, "ArrowFiles"),
                          pattern = "^P.*\\.arrow$", full.names = TRUE)
all_arrows  <- c(pbmc_arrows, meso_arrows)

for (f in all_arrows) {
  lnk <- file.path(arrow_out, basename(f))
  if (!file.exists(lnk)) file.symlink(f, lnk)
}
arrow_links <- list.files(arrow_out, pattern = "\\.arrow$", full.names = TRUE)
message(sprintf("  Arrow files linked: %d", length(arrow_links)))

# ---- 3. Create combined ArchR project ----------------------------------------
message("Creating combined ArchR project...")
if (!file.exists(file.path(out_dir, "Save-ArchR-Project.rds"))) {
  combined <- ArchRProject(
    ArrowFiles       = arrow_links,
    outputDirectory  = out_dir,
    copyArrows       = FALSE,
    showLogo         = FALSE
  )
  combined <- saveArchRProject(combined, overwrite = TRUE)
} else {
  message("  Loading existing project...")
  combined <- loadArchRProject(out_dir, showLogo = FALSE)
}

message(sprintf("  Combined project: %d cells", nrow(combined@cellColData)))

# ---- 4. Subset to desired cells ----------------------------------------------
keep <- intersect(c(onlyATAC, meso_cells), rownames(combined@cellColData))
message(sprintf("  Keeping %d cells (PBMC onlyATAC + meso)", length(keep)))
combined <- combined[keep]

# ---- 5. Add harmonised metadata ----------------------------------------------
# tissue
combined$tissue <- "Mesothelioma"
combined$tissue[rownames(combined@cellColData) %in% onlyATAC] <- "PBMC"

# Unified cell type: map both namespaces to a common label
# PBMC: "Mye_CD14 Mono" -> "CD14_Mono", "Mye_CD16 Mono" -> "CD16_Mono", "Mye_DC" -> "DC"
# Meso: "Mono_CD14" -> "CD14_Mono", "Mono_CD16" -> "CD16_Mono", "cDCs" -> "DC",
#        TAMs -> kept with their names
pbmc_ctmap <- c(
  "Mye_CD14 Mono" = "CD14_Mono",
  "Mye_CD16 Mono" = "CD16_Mono",
  "Mye_DC"        = "DC"
)
meso_ctmap <- c(
  "Mono_CD14"        = "CD14_Mono",
  "Mono_CD16"        = "CD16_Mono",
  "cDCs"             = "DC",
  "TAM_CXCLs"        = "TAM_CXCLs",
  "TAM_interstitial" = "TAM_interstitial",
  "TAM_MARCO"        = "TAM_MARCO",
  "TAM_TREM2"        = "TAM_TREM2"
)

cells_all <- rownames(combined@cellColData)

# PBMC cell type
pbmc_ct <- pbmc_meta$celltype_simplified[match(cells_all, rownames(pbmc_meta))]
pbmc_ct_unified <- pbmc_ctmap[pbmc_ct]

# Meso cell type
meso_ct <- meso_meta$celltype_lv3[match(cells_all, rownames(meso_meta))]
meso_ct_unified <- meso_ctmap[meso_ct]

combined$celltype_unified <- ifelse(combined$tissue == "PBMC",
                                    pbmc_ct_unified, meso_ct_unified)
combined$celltype_unified[is.na(combined$celltype_unified)] <- "Unknown"

# Original cell type (with source namespace)
combined$celltype_orig <- ifelse(combined$tissue == "PBMC", pbmc_ct, as.character(meso_ct))

message("Cell type × tissue distribution:")
print(table(combined$celltype_unified, combined$tissue))

# ---- 6. Dimensionality reduction ---------------------------------------------
message("\nRunning IterativeLSI...")
combined <- addIterativeLSI(
  ArchRProj    = combined,
  useMatrix    = "TileMatrix",
  name         = "IterativeLSI",
  iterations   = 2,
  clusterParams = list(resolution = 0.2, sampleCells = 10000, n.start = 10),
  varFeatures  = 25000,
  LSIMethod    = 2,
  force        = TRUE
)

# ---- 7. Harmony correction: Sample + tissue ----------------------------------
message("Running Harmony (Sample + tissue)...")
combined <- addHarmony(
  ArchRProj   = combined,
  reducedDims = "IterativeLSI",
  name        = "Harmony",
  groupBy     = c("Sample", "tissue"),
  force       = TRUE
)

# ---- 8. UMAP -----------------------------------------------------------------
message("Running UMAP (LSI)...")
combined <- addUMAP(
  ArchRProj   = combined,
  reducedDims = "IterativeLSI",
  name        = "UMAP",
  force       = TRUE
)

message("Running UMAP (Harmony)...")
combined <- addUMAP(
  ArchRProj   = combined,
  reducedDims = "Harmony",
  name        = "UMAP_H",
  force       = TRUE
)

# ---- 9. Clustering -----------------------------------------------------------
message("Running clustering (Harmony)...")
combined <- addClusters(
  input       = combined,
  reducedDims = "Harmony",
  name        = "Clusters_H",
  resolution  = 0.8,
  force       = TRUE
)

# ---- 10. Imputation weights --------------------------------------------------
combined <- addImputeWeights(combined)

# ---- 11. Save project --------------------------------------------------------
combined <- saveArchRProject(combined, overwrite = TRUE)
message("Project saved.")

# ---- 12. UMAP plots ----------------------------------------------------------
message("Plotting UMAPs...")

palette_tissue   <- c(PBMC = "#1f78b4", Mesothelioma = "#e31a1c")
palette_celltype <- c(
  CD14_Mono        = "#d6604d",
  CD16_Mono        = "#f4a582",
  DC               = "#8e7cb8",
  TAM_CXCLs        = "#4393c3",
  TAM_interstitial = "#99d8c9",
  TAM_MARCO        = "#2ca25f",
  TAM_TREM2        = "#74c476",
  Unknown          = "#cccccc"
)

plot_embed <- function(proj, embedding, col_by, pal, fname, ...) {
  p <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData",
                     name = col_by, embedding = embedding,
                     pal = pal, labelMeans = FALSE, ...)
  plotPDF(p, name = fname, ArchRProj = proj, addDOC = FALSE,
          width = 5, height = 5)
}

# Harmony UMAP coloured by tissue, cell type, cluster, sample
plot_embed(combined, "UMAP_H", "tissue",          palette_tissue,   "umap_harmony_tissue.pdf")
plot_embed(combined, "UMAP_H", "celltype_unified", palette_celltype, "umap_harmony_celltype.pdf")
plot_embed(combined, "UMAP_H", "Clusters_H",      NULL,              "umap_harmony_clusters.pdf")

# LSI UMAP for comparison
plot_embed(combined, "UMAP",   "tissue",           palette_tissue,   "umap_lsi_tissue.pdf")
plot_embed(combined, "UMAP",   "celltype_unified",  palette_celltype, "umap_lsi_celltype.pdf")

message("UMAP plots saved to Plots/")

# ---- 13. Marker peaks: CD14 Mono, CD16 Mono, DC vs rest ---------------------
message("\nRunning marker peak analysis...")

# Only use cells with defined unified cell type
cells_for_markers <- rownames(combined@cellColData)[
  combined$celltype_unified %in% c("CD14_Mono","CD16_Mono","DC","TAM_CXCLs",
                                   "TAM_interstitial","TAM_MARCO","TAM_TREM2")]
proj_markers <- combined[cells_for_markers]

# Call peaks if not present
if (!"PeakMatrix" %in% getAvailableMatrices(proj_markers)) {
  message("  Calling peaks per celltype_unified...")
  proj_markers <- addGroupCoverages(
    ArchRProj  = proj_markers,
    groupBy    = "celltype_unified",
    force      = FALSE
  )
  proj_markers <- addReproduciblePeakSet(
    ArchRProj    = proj_markers,
    groupBy      = "celltype_unified",
    pathToMacs2  = findMacs2(),
    force        = FALSE
  )
  proj_markers <- addPeakMatrix(proj_markers, force = FALSE)
  proj_markers <- saveArchRProject(proj_markers, overwrite = TRUE)
}

# Marker features: each cell type vs all others (Wilcoxon)
marker_targets <- c("CD14_Mono", "CD16_Mono", "DC")

for (tgt in marker_targets) {
  rds_out <- file.path(out_dir, paste0("DAP_", tgt, "_vs_rest.rds"))
  if (!file.exists(rds_out)) {
    message(sprintf("  Markers: %s vs rest...", tgt))
    bgd <- setdiff(unique(proj_markers$celltype_unified), tgt)
    # ArchR getMarkerFeatures: useGroups = target, bgdGroups = rest
    mkr <- getMarkerFeatures(
      ArchRProj  = proj_markers,
      useMatrix  = "PeakMatrix",
      groupBy    = "celltype_unified",
      useGroups  = tgt,
      bgdGroups  = bgd,
      testMethod = "wilcoxon",
      bias       = c("TSSEnrichment", "log10(nFrags)"),
      k          = 100
    )
    saveRDS(mkr, rds_out)
    message(sprintf("    Saved %s", basename(rds_out)))
  } else {
    message(sprintf("  %s already computed, skipping.", tgt))
  }
}

# ---- 14. Marker volcano / MA plots -------------------------------------------
message("Plotting marker results...")
for (tgt in marker_targets) {
  rds_out <- file.path(out_dir, paste0("DAP_", tgt, "_vs_rest.rds"))
  if (file.exists(rds_out)) {
    mkr <- readRDS(rds_out)
    p_ma <- markerPlot(seMarker = mkr, name = tgt,
                       cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1",
                       plotAs = "MA")
    p_vol <- markerPlot(seMarker = mkr, name = tgt,
                        cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1",
                        plotAs = "Volcano")
    plotPDF(p_ma, p_vol,
            name    = paste0("markers_", tgt, "_MA_Volcano.pdf"),
            ArchRProj = proj_markers,
            addDOC  = FALSE, width = 5, height = 5)
  }
}

# Summary table of top markers per target
message("\nTop DA peaks per target (FDR < 0.05, Log2FC > 1):")
for (tgt in marker_targets) {
  rds_out <- file.path(out_dir, paste0("DAP_", tgt, "_vs_rest.rds"))
  if (file.exists(rds_out)) {
    mkr <- readRDS(rds_out)
    col <- colnames(mkr)[1]
    rd  <- as.data.frame(rowData(mkr))
    res <- data.frame(
      Log2FC = as.numeric(assays(mkr)[["Log2FC"]][[col]]),
      FDR    = as.numeric(assays(mkr)[["FDR"]][[col]]),
      region = paste0(rd$seqnames, ":", rd$start, "-", rd$end),
      stringsAsFactors = FALSE
    )
    sig <- res[!is.na(res$FDR) & res$FDR < 0.05 & res$Log2FC > 1, ]
    sig <- sig[order(-sig$Log2FC), ]
    message(sprintf("  %s: %d significant peaks", tgt, nrow(sig)))
    if (nrow(sig) > 0)
      write.csv(sig, file.path(out_dir, paste0("DAP_", tgt, "_sig.csv")),
                row.names = FALSE, quote = FALSE)
  }
}

message("\nDone. Project at: ", out_dir)
