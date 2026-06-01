# plot_sample_cd14mono2_enrichment.R  (R1_Q4)
#
# Checks which samples are enriched for CD14 Mono 2 prediction score.
# Sample is extracted from the cell name (Sample#barcode format).
#
# Outputs:
#   sample_cd14mono2_boxplot.pdf   score distribution per sample (ordered by median)
#   sample_cd14mono2_summary.csv   per-sample stats (n, mean, median, pct_high)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(paletteer)
})

utils_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils"
source(file.path(utils_dir, "ggplot_aestetics.R"))

# ── Paths ──────────────────────────────────────────────────────────────────────
score_csv <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection/per_label_scores.csv"
out_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/pbmc_projection"

# ── Load & extract sample ──────────────────────────────────────────────────────
df <- read.csv(score_csv, check.names = FALSE)
df$sample <- sub("#.*", "", df$cell)

message("Samples: ", paste(sort(unique(df$sample)), collapse = ", "))
message("Cells per sample:"); print(table(df$sample))

# ── Order by median score ──────────────────────────────────────────────────────
med_order <- df |>
  group_by(sample) |>
  summarise(median_score = median(`12_CD14.Mono.2`, na.rm = TRUE),
            mean_score   = mean(`12_CD14.Mono.2`,   na.rm = TRUE),
            pct_high     = mean(`12_CD14.Mono.2` > 0.4, na.rm = TRUE),
            n_cells      = n(),
            .groups = "drop") |>
  arrange(desc(median_score))

message("\nPer-sample summary:")
print(as.data.frame(med_order))

write.csv(med_order, file.path(out_dir, "sample_cd14mono2_summary.csv"),
          row.names = FALSE)

df$sample <- factor(df$sample, levels = med_order$sample)

# ── Sample colour palette ──────────────────────────────────────────────────────
n_samples   <- length(levels(df$sample))
sample_cols <- setNames(
  as.character(paletteer::paletteer_d("ggsci::default_igv")[seq_len(n_samples)]),
  levels(df$sample)
)

# ── Boxplot ────────────────────────────────────────────────────────────────────
p <- ggplot(df, aes(x = sample, y = `12_CD14.Mono.2`, fill = sample)) +
  bxp +
  scale_fill_manual(values = sample_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 1),
                     name   = "CD14 Mono 2 prediction score") +
  labs(x = NULL, title = "CD14 Mono 2 score by sample (ordered by median)") +
  gtheme +
  theme(plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text  = element_text(size = 8))

ggsave(file.path(out_dir, "sample_cd14mono2_boxplot.pdf"), p,
       width = 5, height = 3.5)
message("Saved boxplot and summary CSV.")
