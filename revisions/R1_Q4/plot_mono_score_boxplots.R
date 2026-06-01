# plot_mono_score_boxplots.R  (R1_Q4)
#
# Boxplots of KNN prediction score for 12_CD14.Mono.2 grouped by celltype_lv3,
# ordered by median score, coloured by myeloid_palette.
#
# Input  : R1_Q4/output/pbmc_projection/per_label_scores.csv
# Output : R1_Q4/output/pbmc_projection/boxplot_mono_scores.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(paletteer)
})

utils_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils"
source(file.path(utils_dir, "ggplot_aestetics.R"))

palette_myeloid <- rev(paletteer::paletteer_d("MoMAColors::VanGogh"))
palette_myeloid <- setNames(as.character(palette_myeloid),
  c("Mono_CD14", "TAM_CXCLs", "TAM_MARCO", "TAM_TREM2", "cDCs", "Mono_CD16", "TAM_interstitial"))

# ── Paths ──────────────────────────────────────────────────────────────────────
score_csv <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection/per_label_scores.csv"
out_file  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection/boxplot_mono_scores.pdf"

# ── Load ───────────────────────────────────────────────────────────────────────
df <- read.csv(score_csv, check.names = FALSE)
message("Cells: ", nrow(df))

# ── Order celltype_lv3 by median 12_CD14.Mono.2 score ─────────────────────────
med_order <- df |>
  (\(d) tapply(d[["12_CD14.Mono.2"]], d$celltype_lv3, median, na.rm = TRUE))() |>
  sort(decreasing = TRUE) |>
  names()

df$celltype_lv3 <- factor(df$celltype_lv3, levels = med_order)

# ── Plot ───────────────────────────────────────────────────────────────────────
p <- ggplot(df, aes(x = celltype_lv3, y = `12_CD14.Mono.2`, fill = celltype_lv3)) +
  vlp +
  bxpv +
  scale_fill_manual(values = palette_myeloid, guide = "none") +
  scale_y_continuous(
    limits = c(0, 1),
    name   = "KNN prediction score\n(fraction of 25 neighbors)"
  ) +
  labs(
    x     = NULL,
    title = "CD14 Monocyte 2 (12_CD14.Mono.2) prediction score"
  ) +
  gtheme +
  theme(plot.title = element_text(size = 8))

ggsave(out_file, p, width = 4, height = 3)
message("Saved: ", out_file)
