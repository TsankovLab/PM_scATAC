# plot_all_label_scores_boxplots.R  (R1_Q4)
#
# Violin + boxplot of KNN prediction score for every reference BioClassification
# label, split by meso myeloid celltype_lv3. One facet per reference label.
#
# Input  : R1_Q4/output/pbmc_projection/per_label_scores.csv
# Output : R1_Q4/output/pbmc_projection/all_label_scores_boxplots.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(paletteer)
})

utils_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils"
source(file.path(utils_dir, "ggplot_aestetics.R"))

palette_myeloid <- rev(paletteer::paletteer_d("MoMAColors::VanGogh"))
palette_myeloid <- setNames(as.character(palette_myeloid),
  c("Mono_CD14","TAM_CXCLs","TAM_MARCO","TAM_TREM2","cDCs","Mono_CD16","TAM_interstitial"))

# ── Paths ──────────────────────────────────────────────────────────────────────
score_csv <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection/per_label_scores.csv"
out_file  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection/all_label_scores_boxplots.pdf"

# ── Load ───────────────────────────────────────────────────────────────────────
df <- read.csv(score_csv, check.names = FALSE)

# all BioClassification columns (everything except cell and celltype_lv3)
label_cols <- setdiff(colnames(df), c("cell", "celltype_lv3"))
message("Reference labels: ", length(label_cols))

# ── celltype_lv3 order: by median 12_CD14.Mono.2 score (consistent) ───────────
ct_order <- names(sort(tapply(df[["12_CD14.Mono.2"]], df$celltype_lv3,
                              median, na.rm = TRUE), decreasing = TRUE))
df$celltype_lv3 <- factor(df$celltype_lv3, levels = ct_order)

# ── Pivot to long ──────────────────────────────────────────────────────────────
long_df <- df |>
  pivot_longer(cols      = all_of(label_cols),
               names_to  = "ref_label",
               values_to = "score")

# order facets: most enriched labels first (by overall mean score)
label_order <- long_df |>
  group_by(ref_label) |>
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(mean_score)) |>
  pull(ref_label)

long_df$ref_label <- factor(long_df$ref_label, levels = label_order)

# ── Plot ───────────────────────────────────────────────────────────────────────
p <- ggplot(long_df, aes(x = celltype_lv3, y = score, fill = celltype_lv3)) +
  vlp +
  bxpv +
  facet_wrap(~ ref_label, ncol = 8, scales = "fixed") +
  scale_fill_manual(values = palette_myeloid, guide = "none") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0, 0.5, 1),
                     name   = "KNN prediction score") +
  labs(x = NULL, title = "Reference BioClassification prediction scores by myeloid subtype") +
  gtheme +
  theme(
    plot.title      = element_text(size = 9),
    axis.title.y    = element_text(size = 7),
    axis.text.x     = element_text(size = 5, angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 5),
    strip.text      = element_text(size = 6, face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    panel.spacing   = unit(1.5, "mm")
  )

nrows  <- ceiling(length(label_cols) / 8)
ggsave(out_file, p, width = 14, height = nrows * 2.2)
message("Saved: ", out_file)
