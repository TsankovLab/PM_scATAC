# run_chromvar_combined.R  (R1_Q4)
#
# 1. Adds cisbp motif annotations to the combined meso+PBMC ArchR project
# 2. Runs chromVAR (addDeviationsMatrix) on the combined peak set
# 3. Compares the 14 MDM/inflammatory-axis TFs across cell types × tissue
#    producing a heatmap and barplot saved to git_repo_claude/R1_Q4/
#
# Submit via submit_chromvar_combined.sh (36 CPUs, 576 GB, 8 h)

suppressPackageStartupMessages({
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(patchwork)
})
addArchRGenome("hg38")
addArchRThreads(threads = 36)
set.seed(1234)

proj_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/meso_vs_pbmc"
out_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4"

MDM_TFs <- c("FOSL1","FOSL2","BACH1","PPARG","NFKB2","KLF12",
             "HIVEP3","SMAD1","NFKB1","REL","RUNX1","SNAI1","RUNX2","NFAT5")

message("Loading project...")
archp <- loadArchRProject(proj_dir, showLogo = FALSE)

# ---- 0. Unified peak set (callPeaks.R logic) ----------------------------------
# Arrow files from PBMC and meso individual projects have different peak sets
# causing rowSums mismatch. Call a unified peak set on the combined project.
macs2_path       <- "/sc/arion/work/giottb01/conda/envs/meso_scatac/bin/macs2"
metaGroupName    <- "celltype_unified"
peak_reproducibility <- 2

message("Calling unified peaks on combined project (metaGroupName = ", metaGroupName, ")...")
archp <- addGroupCoverages(
  ArchRProj     = archp,
  groupBy       = metaGroupName,
  force         = TRUE,
  minCells      = 20,
  maxCells      = 500,
  minReplicates = 2,
  sampleRatio   = 0.8,
  useLabels     = TRUE
)

pdf(file.path(tempdir(), "tmp_peaks.pdf"))
archp <- addReproduciblePeakSet(
  ArchRProj       = archp,
  groupBy         = metaGroupName,
  peakMethod      = "Macs2",
  pathToMacs2     = macs2_path,
  reproducibility = as.character(peak_reproducibility),
  maxPeaks        = 500000,
  minCells        = 20,
  force           = TRUE
)
dev.off()

archp <- addPeakMatrix(archp)
archp <- saveArchRProject(ArchRProj = archp, load = TRUE)
message(sprintf("  Unified peak set: %d peaks", length(archp@peakSet)))

# ---- 1. Motif annotations ----------------------------------------------------
message("Adding motif annotations (cisbp)...")
archp <- addMotifAnnotations(
  ArchRProj = archp,
  motifSet  = "cisbp",
  name      = "Motif",
  force     = FALSE
)

# ---- 2. chromVAR deviations --------------------------------------------------
message("Running chromVAR (addDeviationsMatrix)...")
archp <- addDeviationsMatrix(
  ArchRProj       = archp,
  peakAnnotation  = "Motif",
  force           = TRUE
)

# ---- 3. Save -----------------------------------------------------------------
archp <- saveArchRProject(ArchRProj = archp, load = TRUE)
message("Project saved.")

# ---- 4. Build tissue × celltype grouping ------------------------------------
# Shared types (CD14_Mono, CD16_Mono, DC) exist in both tissues;
# TAM subtypes are meso-only. Create a combined label for getGroupSE.
ct   <- as.character(archp$celltype_unified)
tiss <- as.character(archp$tissue)

archp$tissue_ct <- ifelse(
  ct %in% c("CD14_Mono", "CD16_Mono", "DC"),
  paste0(tiss, "_", ct),
  ct
)

# Define display order: Blood first, then Meso mono/DC, then TAMs
ct_order <- c(
  "PBMC_CD14_Mono", "PBMC_CD16_Mono", "PBMC_DC",
  "Mesothelioma_CD14_Mono", "Mesothelioma_CD16_Mono", "Mesothelioma_DC",
  "TAM_CXCLs", "TAM_interstitial", "TAM_MARCO", "TAM_TREM2"
)
ct_order <- ct_order[ct_order %in% unique(archp$tissue_ct)]
archp$tissue_ct <- factor(archp$tissue_ct, levels = ct_order)

# ---- 5. Average deviations per group ----------------------------------------
message("Averaging MotifMatrix per tissue_ct group...")
devSE <- getGroupSE(
  ArchRProj = archp,
  useMatrix  = "MotifMatrix",
  groupBy    = "tissue_ct",
  divideN    = TRUE,
  scaleTo    = NULL
)

devMat <- assays(devSE)[[1]]
rownames(devMat) <- rowData(devSE)$name
rownames(devMat) <- gsub("_.*", "", rownames(devMat))   # strip motif ID suffix

# Keep only columns in defined order
devMat <- devMat[, ct_order[ct_order %in% colnames(devMat)], drop = FALSE]

# ---- 6. Extract MDM TFs ------------------------------------------------------
found_tfs <- intersect(MDM_TFs, rownames(devMat))
missing   <- setdiff(MDM_TFs, found_tfs)
if (length(missing))
  message("  TFs not found: ", paste(missing, collapse = ", "))
message(sprintf("  Using %d / %d MDM TFs", length(found_tfs), length(MDM_TFs)))

sub_mat <- devMat[found_tfs, , drop = FALSE]

# Save numeric table
out_tab <- as.data.frame(sub_mat)
out_tab$TF <- rownames(out_tab)
write.csv(out_tab[, c("TF", colnames(sub_mat))],
          file.path(out_dir, "MDM_TF_chromvar_combined.csv"),
          row.names = FALSE, quote = FALSE)

# ---- 7. Heatmap --------------------------------------------------------------
message("Building heatmap...")

# Row z-score across all groups
sub_z <- t(scale(t(sub_mat)))
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#d6604d"))

# Column split: Blood vs Meso
col_split <- factor(
  ifelse(grepl("^PBMC", colnames(sub_z)), "Blood (NSCLC)", "Mesothelioma"),
  levels = c("Blood (NSCLC)", "Mesothelioma")
)

# Column colour annotation
source_pal <- c("Blood (NSCLC)" = "#1f78b4", "Mesothelioma" = "#e31a1c")
ha <- HeatmapAnnotation(
  Source = col_split,
  col    = list(Source = source_pal),
  annotation_name_gp = gpar(fontsize = 7.5),
  show_legend = TRUE
)

# Prettier column labels
col_labels <- colnames(sub_z)
col_labels <- gsub("PBMC_",          "Blood ",         col_labels)
col_labels <- gsub("Mesothelioma_",  "Meso ",          col_labels)
col_labels <- gsub("_", " ",         col_labels)

ht <- Heatmap(
  sub_z,
  name              = "z-score",
  col               = col_fun,
  top_annotation    = ha,
  column_split      = col_split,
  cluster_rows      = TRUE,
  cluster_columns   = FALSE,
  clustering_distance_rows = "euclidean",
  row_names_gp      = gpar(fontsize = 9, fontface = "italic"),
  row_names_side    = "left",
  column_labels     = col_labels,
  column_names_gp   = gpar(fontsize = 8),
  column_names_rot  = 45,
  column_title_gp   = gpar(fontsize = 9, fontface = "bold"),
  border            = TRUE,
  rect_gp           = gpar(col = "white", lwd = 0.5),
  width             = unit(ncol(sub_z) * 8, "mm"),
  height            = unit(length(found_tfs) * 7, "mm"),
  heatmap_legend_param = list(
    title          = "z-score\n(row-scaled)",
    title_gp       = gpar(fontsize = 8, fontface = "bold"),
    labels_gp      = gpar(fontsize = 7),
    legend_height  = unit(3, "cm"),
    at             = c(-2, -1, 0, 1, 2)
  )
)

pdf(file.path(out_dir, "MDM_TF_chromvar_heatmap.pdf"), width = 9, height = 5.5)
draw(ht,
     column_title    = "MDM / inflammatory-axis TF activity (chromVAR z-score)",
     column_title_gp = gpar(fontsize = 10, fontface = "bold"),
     padding         = unit(c(4, 4, 8, 4), "mm"))
dev.off()
message("Saved MDM_TF_chromvar_heatmap.pdf")

# ---- 8. Mean MDM TF activity barplot ----------------------------------------
mean_act <- colMeans(sub_z, na.rm = TRUE)
bar_df <- data.frame(
  celltype = names(mean_act),
  mean_z   = as.numeric(mean_act),
  source   = ifelse(grepl("^PBMC", names(mean_act)), "Blood (NSCLC)", "Mesothelioma"),
  stringsAsFactors = FALSE
)
bar_df$celltype <- factor(bar_df$celltype,
  levels = bar_df$celltype[order(bar_df$mean_z, decreasing = TRUE)])
bar_df$label <- col_labels[match(as.character(bar_df$celltype), colnames(sub_z))]
bar_df$label <- factor(bar_df$label, levels = bar_df$label[order(bar_df$mean_z, decreasing=TRUE)])

bar_fill <- setNames(
  c(rep("#1f78b4", 3), rep("#e31a1c", 3),
    "#4393c3", "#99d8c9", "#2ca25f", "#74c476")[seq_along(ct_order)],
  ct_order
)

p_bar <- ggplot(bar_df, aes(x = label, y = mean_z,
                             fill = celltype)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey50") +
  scale_fill_manual(values = bar_fill) +
  labs(
    x     = NULL,
    y     = "Mean z-score across MDM TFs",
    title = "MDM / inflammatory TF activity — combined project"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

pdf(file.path(out_dir, "MDM_TF_chromvar_barplot.pdf"), width = 6, height = 4)
print(p_bar)
dev.off()
message("Saved MDM_TF_chromvar_barplot.pdf")

message("\nDone. Outputs in ", out_dir)
