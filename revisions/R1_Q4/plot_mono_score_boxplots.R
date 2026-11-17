# plot_mono_score_boxplots.R  (R1_Q4)
#
# Boxplots of KNN prediction score for 12_CD14.Mono.2 grouped by celltype_lv3,
# ordered by median score, coloured by myeloid_palette.
#
# Input  : R1_Q4/output/pbmc_projection/per_label_scores.csv
# Output : R1_Q4/output/pbmc_projection/boxplot_mono_scores.pdf

source(local({a <- grep("^--file=", commandArgs(FALSE), value = TRUE); for (d in c(if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))), getwd())) { while (!file.exists(file.path(d, "config.R")) && d != dirname(d)) d <- dirname(d); if (file.exists(file.path(d, "config.R"))) break }; Sys.getenv("MESO_CONFIG", file.path(d, "config.R"))}))   # data / results / reference roots: see config.R at the repo root
suppressPackageStartupMessages({
  library(ggplot2)
  library(paletteer)
})

utils_dir <- paste0(REPO_DIR, "/utils")
source(file.path(utils_dir, "ggplot_aestetics.R"))

palette_myeloid <- rev(paletteer::paletteer_d("MoMAColors::VanGogh"))
palette_myeloid <- setNames(as.character(palette_myeloid),
  c("Mono_CD14", "TAM_CXCLs", "TAM_MARCO", "TAM_TREM2", "cDCs", "Mono_CD16", "TAM_interstitial"))

# ── Paths ──────────────────────────────────────────────────────────────────────
score_csv <- paste0(RESULTS_DIR, "/R1_Q4/output/pbmc_projection/per_label_scores.csv")
out_file  <- paste0(RESULTS_DIR, "/R1_Q4/output/pbmc_projection/boxplot_mono_scores.pdf")

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
