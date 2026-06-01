# plot_mod2_vs_cd14mono2.R  (R1_Q4)
#
# mod_2 vs 12_CD14.Mono.2 scatter with top + right density marginals.
# Layout (left→right): scatter | right density | legend
#                      top density above scatter, spacer above right/legend
#
# Inputs:
#   myeloid scatac_ArchR project   mod_2, celltype_lv3
#   per_label_scores.csv           12_CD14.Mono.2 prediction score
#
# Output: R1_Q4/output/pbmc_projection/scatter_mod2_cd14mono2.pdf

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(cowplot)
  library(paletteer)
})

addArchRThreads(1)
addArchRGenome("hg19")

utils_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils"
source(file.path(utils_dir, "ggplot_aestetics.R"))

palette_myeloid <- rev(paletteer::paletteer_d("MoMAColors::VanGogh"))
palette_myeloid <- setNames(as.character(palette_myeloid),
  c("Mono_CD14","TAM_CXCLs","TAM_MARCO","TAM_TREM2","cDCs","Mono_CD16","TAM_interstitial"))

# ── Paths ──────────────────────────────────────────────────────────────────────
proj_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR"
score_csv <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection/per_label_scores.csv"
out_file  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection/scatter_mod2_cd14mono2.pdf"

ct_map <- c(
  cnmf_cluster_1 = "Mono_CD16",   cnmf_cluster_2 = "TAM_MARCO",
  cnmf_cluster_3 = "Mono_CD14",   cnmf_cluster_4 = "TAM_CXCLs",
  cnmf_cluster_5 = "TAM_TREM2",   cnmf_cluster_6 = "TAM_interstitial",
  cnmf_cluster_7 = "cDCs"
)

# ── Load ──────────────────────────────────────────────────────────────────────
message("Loading ArchR project...")
archp <- loadArchRProject(proj_dir, showLogo = FALSE)
if (!"celltype_lv3" %in% colnames(archp@cellColData))
  archp$celltype_lv3 <- ct_map[as.character(archp$cnmf_cluster)]

meta <- data.frame(
  cell         = getCellNames(archp),
  mod_2        = as.numeric(archp$mod_2),
  celltype_lv3 = as.character(archp$celltype_lv3),
  stringsAsFactors = FALSE
)
scores <- read.csv(score_csv, check.names = FALSE)[, c("cell","12_CD14.Mono.2")]
df     <- merge(meta, scores, by = "cell")
df     <- df[!is.na(df$mod_2) & !is.na(df[["12_CD14.Mono.2"]]), ]
df$celltype_lv3 <- factor(df$celltype_lv3, levels = names(palette_myeloid))
set.seed(42); df <- df[sample(nrow(df)), ]
message("Cells: ", nrow(df))

# ── Spearman ───────────────────────────────────────────────────────────────────
ct        <- cor.test(df$mod_2, df[["12_CD14.Mono.2"]], method = "spearman")
pv_label  <- if (ct$p.value < 2.2e-16) "p < 2.2e-16" else sprintf("p = %.2e", ct$p.value)
cor_label <- sprintf("rho = %.2f\n%s", ct$estimate, pv_label)

# ── Shared density theme ───────────────────────────────────────────────────────
dens_theme <- gtheme_no_text + theme(
  axis.text    = element_text(size = 6),
  axis.title   = element_text(size = 7),
  legend.position = "none"
)

# ── Top density: mod_2 ────────────────────────────────────────────────────────
p_top <- ggplot(df, aes(x = mod_2, fill = celltype_lv3, colour = celltype_lv3)) +
  geom_density(alpha = 0.25, linewidth = 0.4) +
  scale_fill_manual(values = palette_myeloid, guide = "none") +
  scale_colour_manual(values = palette_myeloid, guide = "none") +
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(x = NULL, y = "Density") +
  dens_theme

# ── Right density: 12_CD14.Mono.2 (rotated) ───────────────────────────────────
p_right <- ggplot(df, aes(x = `12_CD14.Mono.2`, fill = celltype_lv3, colour = celltype_lv3)) +
  geom_density(alpha = 0.25, linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = palette_myeloid, guide = "none") +
  scale_colour_manual(values = palette_myeloid, guide = "none") +
  scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1)) +
  labs(x = NULL, y = "Density") +
  dens_theme

# ── Scatter (no legend) ────────────────────────────────────────────────────────
p_scatter <- ggplot(df, aes(x = mod_2, y = `12_CD14.Mono.2`, colour = celltype_lv3)) +
  geom_point(size = 0.3, stroke = 0, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30",
              linewidth = 0.5, fill = "grey80") +
  annotate("text", x = Inf, y = Inf, label = cor_label,
           hjust = 1.1, vjust = 1.5, size = 2.8, colour = "grey20") +
  scale_colour_manual(
    values = palette_myeloid, name = "celltype_lv3",
    guide  = guide_legend(override.aes = list(size = 2.5, alpha = 1))
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Inflammatory epigenomic signature (mod_2)",
       y = "CD14 Mono 2 prediction score",
       title = "mod_2 vs 12_CD14.Mono.2 KNN score") +
  gtheme +
  theme(plot.title   = element_text(size = 8),
        axis.title   = element_text(size = 8),
        axis.text    = element_text(size = 7),
        legend.position = "right")

# ── Extract legend ─────────────────────────────────────────────────────────────
leg <- cowplot::get_legend(
  p_scatter +
  theme(legend.text  = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.size = unit(3, "mm"))
)
p_scatter_nl <- p_scatter + theme(legend.position = "none")

# ── Assemble: cowplot::plot_grid ───────────────────────────────────────────────
# Widths:  scatter=4, right_dens=0.8, legend=1.2
# Heights: top_dens=1, scatter=3
bottom_row <- plot_grid(p_scatter_nl, p_right, leg,  nrow = 1, rel_widths = c(4, 1.2, 1.2),
                        align = "h", axis = "tb")
top_row    <- plot_grid(p_top,        NULL,    NULL, nrow = 1, rel_widths = c(4, 1.2, 1.2),
                        align = "h", axis = "tb")
final      <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 3),
                        align = "v", axis = "lr")

ggsave(out_file, final, width = 5.5, height = 3.5)
message("Saved: ", out_file)
