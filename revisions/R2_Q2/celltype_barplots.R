# celltype_barplots.R  (R2_Q2)
#
# Two-panel barplot comparing cell-type composition across samples:
#   Panel 1: proportions per sample (replicates main analysis celltype_barplots.pdf)
#   Panel 2: raw cell counts per sample, annotated with per-sample total
#
# Output: celltype_barplots.pdf  (height=6, width=11)
#
# Requires cluster: 2 CPUs, 32 GB, ~20 min
# Run via submit_celltype_barplots.sh

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(patchwork)
  library(paletteer)
  library(circlize)
  library(RColorBrewer)
})

addArchRGenome("hg38")
addArchRThreads(threads = 2)

script_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q2"
utils_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils"
main_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR"

source(file.path(utils_dir, "ggplot_aestetics.R"))
source(file.path(utils_dir, "useful_functions.R"))
source(file.path(utils_dir, "palettes.R"))

# ---- Load project ------------------------------------------------------------
message("Loading main ArchR project...")
archp <- loadArchRProject(main_dir, showLogo = FALSE)

# ---- Prepare metadata --------------------------------------------------------
archp_meta <- as.data.frame(archp@cellColData)

# Align SmoothMuscle label with palette key (palette.R renames it to 'Smooth Muscle')
archp_meta$celltype_lv1 <- as.character(archp_meta$celltype_lv1)
archp_meta$celltype_lv1[archp_meta$celltype_lv1 == "SmoothMuscle"] <- "Smooth Muscle"
archp_meta$celltype_lv1 <- factor(archp_meta$celltype_lv1, levels = names(palette_celltype_lv1))

sample_order <- rev(c("P1","P3","P4","P5","P8","P10","P11","P12","P13","P14","P23"))
archp_meta$Sample2 <- factor(archp_meta$Sample, levels = sample_order)

# ---- Panel 1: proportions ----------------------------------------------------
bp_prop <- cellComp(
  seurat_obj  = archp_meta,
  metaGroups  = c("Sample2", "celltype_lv1"),
  plot_as     = "bar",
  pal         = palette_celltype_lv1
) + gtheme + coord_flip() +
  labs(x = NULL, y = "Proportion", fill = "Cell type") +
  theme(legend.position = "none")

# ---- Panel 2: raw counts with total annotation -------------------------------
counts_df <- cellComp(
  seurat_obj  = archp_meta,
  metaGroups  = c("Sample2", "celltype_lv1"),
  plot_as     = "bar",
  pal         = palette_celltype_lv1,
  prop        = FALSE,
  returnDF    = TRUE
)

totals_df <- aggregate(Freq ~ Sample2, data = counts_df, FUN = sum)
colnames(totals_df)[2] <- "total"

bp_counts <- ggplot(counts_df, aes(x = Sample2, y = Freq, fill = celltype_lv1)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_text(
    data        = totals_df,
    aes(x = Sample2, y = total, label = paste0("n=", total)),
    hjust       = -0.15,
    size        = 2.8,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = palette_celltype_lv1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  gtheme + coord_flip() +
  labs(x = NULL, y = "Cell count", fill = "Cell type")

# ---- Combine and save --------------------------------------------------------
combined <- bp_prop + bp_counts + plot_layout(guides = "collect") &
  theme(legend.position = "right")

pdf(file.path(script_dir, "celltype_barplots.pdf"), height = 6, width = 11)
print(combined)
dev.off()

message("Done. Output: ", file.path(script_dir, "celltype_barplots.pdf"))
