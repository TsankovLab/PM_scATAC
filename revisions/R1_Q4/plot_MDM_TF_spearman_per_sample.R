# plot_MDM_TF_spearman_per_sample.R  (R1_Q4)
#
# Per-sample Spearman correlations of 14 MDM TF chromVAR z-scores between:
#   - Each meso cell type (per meso sample, >= 10 cells required)
#   - Each blood reference cell type (PBMC_CD14_Mono, PBMC_CD16_Mono, PBMC_DC)
#     averaged across all PBMC cells (stable reference)
#
# Output: one Spearman r per (meso sample × meso cell type × blood cell type)
# Boxplot: x = meso cell type, y = r, colour = blood reference cell type
#
# Run on login node: Rscript plot_MDM_TF_spearman_per_sample.R

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(SummarizedExperiment)
  library(ggpubr)
})
addArchRGenome("hg38")
addArchRThreads(4)

proj_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/meso_vs_pbmc"
out_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4"

MDM_TFs <- c("FOSL1","FOSL2","BACH1","PPARG","NFKB2","KLF12",
             "HIVEP3","SMAD1","NFKB1","REL","RUNX1","SNAI1","RUNX2","NFAT5")

MIN_CELLS <- 50   # minimum cells per sample × cell type group

message("Loading project...")
proj <- loadArchRProject(proj_dir, showLogo = FALSE)

# Recreate tissue_ct column (same logic as run_chromvar_combined.R)
ct   <- as.character(proj$celltype_unified)
tiss <- as.character(proj$tissue)
proj$tissue_ct <- ifelse(
  ct %in% c("CD14_Mono", "CD16_Mono", "DC"),
  paste0(tiss, "_", ct),
  ct
)

# Sample × tissue_ct grouping label
proj$sample_tct <- paste0(proj$Sample, "__", proj$tissue_ct)

# ---- Average MotifMatrix per sample × cell type ------------------------------
message("Averaging MotifMatrix per sample x tissue_ct...")
devSE <- getGroupSE(proj,
                    useMatrix = "MotifMatrix",
                    groupBy   = "sample_tct",
                    divideN   = TRUE,
                    scaleTo   = NULL)

# Tidy row names: strip motif-ID suffix (e.g. "FOSL1_1" -> "FOSL1")
rownames(devSE) <- gsub("_.*", "", rowData(devSE)$name)

# Keep only MDM TFs present in the matrix
found_tfs <- intersect(MDM_TFs, rownames(devSE))
message(sprintf("  Found %d / %d MDM TFs", length(found_tfs), length(MDM_TFs)))
devSE <- devSE[found_tfs, ]

# Convert to data.frame: rows = TFs, cols = sample__tissue_ct groups
mat <- assays(devSE)[[1]]

# Parse column names back to sample + tissue_ct
col_meta <- data.frame(
  col_name  = colnames(mat),
  sample    = sub("__.*", "", colnames(mat)),
  tissue_ct = sub(".*__", "", colnames(mat)),
  stringsAsFactors = FALSE
)
col_meta$tissue <- ifelse(grepl("^PBMC", col_meta$tissue_ct), "PBMC", "Mesothelioma")

# ---- Cell count filter -------------------------------------------------------
cell_counts <- table(proj$sample_tct)
col_meta$n_cells <- as.integer(cell_counts[col_meta$col_name])
col_meta$n_cells[is.na(col_meta$n_cells)] <- 0

message(sprintf("  Cell count range across groups: %d – %d",
                min(col_meta$n_cells), max(col_meta$n_cells)))
message(sprintf("  Groups with < %d cells (excluded): %d",
                MIN_CELLS, sum(col_meta$n_cells < MIN_CELLS)))
if (any(col_meta$n_cells < MIN_CELLS))
  message("    ", paste(col_meta$col_name[col_meta$n_cells < MIN_CELLS], collapse=", "))

col_meta_pass <- col_meta[col_meta$n_cells >= MIN_CELLS, ]
message(sprintf("  Groups passing filter (>= %d cells): %d / %d",
                MIN_CELLS, nrow(col_meta_pass), nrow(col_meta)))

# ---- Blood reference: average TF profile per blood cell type -----------------
blood_types <- c("PBMC_CD14_Mono")
blood_ref <- sapply(blood_types, function(bct) {
  cols <- col_meta_pass$col_name[col_meta_pass$tissue_ct == bct]
  if (length(cols) == 0) return(rep(NA, length(found_tfs)))
  # Weight by cell count for proper average
  wts <- col_meta_pass$n_cells[col_meta_pass$tissue_ct == bct]
  rowSums(mat[, cols, drop = FALSE] * rep(wts, each = nrow(mat))) / sum(wts)
})
rownames(blood_ref) <- found_tfs
message("Blood reference cell types: ", paste(colnames(blood_ref), collapse = ", "))

# ---- Per-meso-sample Spearman correlations -----------------------------------
meso_types <- c("Mesothelioma_CD14_Mono", "Mesothelioma_CD16_Mono",
                "Mesothelioma_DC", "TAM_CXCLs", "TAM_interstitial",
                "TAM_MARCO", "TAM_TREM2")

meso_meta <- col_meta_pass[col_meta_pass$tissue_ct %in% meso_types, ]
message(sprintf("  Meso sample x cell type groups: %d", nrow(meso_meta)))

cor_rows <- lapply(seq_len(nrow(meso_meta)), function(i) {
  meso_col  <- meso_meta$col_name[i]
  meso_prof <- mat[, meso_col]
  lapply(blood_types, function(bct) {
    blood_prof <- blood_ref[, bct]
    valid <- !is.na(meso_prof) & !is.na(blood_prof)
    if (sum(valid) < 5) return(NULL)
    r <- cor(meso_prof[valid], blood_prof[valid], method = "spearman")
    data.frame(
      meso_sample   = meso_meta$sample[i],
      meso_ct       = meso_meta$tissue_ct[i],
      blood_ct      = bct,
      spearman_r    = r,
      n_cells       = meso_meta$n_cells[i],
      stringsAsFactors = FALSE
    )
  })
})
cor_df <- do.call(rbind, Filter(Negate(is.null), unlist(cor_rows, recursive = FALSE)))

# Clean up labels for plotting
cor_df$meso_ct_label <- gsub("Mesothelioma_", "Meso ", cor_df$meso_ct)
cor_df$meso_ct_label <- gsub("_", " ", cor_df$meso_ct_label)
cor_df$blood_ct_label <- c(
  PBMC_CD14_Mono = "Blood CD14 Mono",
  PBMC_CD16_Mono = "Blood CD16 Mono",
  PBMC_DC        = "Blood DC"
)[cor_df$blood_ct]

# Order meso cell types by descending median Spearman r vs Blood CD14 Mono
ct_median_order <- cor_df %>%
  group_by(meso_ct_label) %>%
  summarise(med = median(spearman_r, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med))
cor_df$meso_ct_label <- factor(cor_df$meso_ct_label,
                                levels = ct_median_order$meso_ct_label)
cor_df$blood_ct_label <- factor(cor_df$blood_ct_label,
                                 levels = c("Blood CD14 Mono","Blood CD16 Mono","Blood DC"))

blood_pal <- c(
  "Blood CD14 Mono" = "#d6604d",
  "Blood CD16 Mono" = "#f4a582",
  "Blood DC"        = "#8e7cb8"
)

# Save table
write.csv(cor_df, file.path(out_dir, "MDM_TF_spearman_per_sample.csv"),
          row.names = FALSE, quote = FALSE)

# ---- Pairwise t-tests --------------------------------------------------------
ct_levels <- levels(cor_df$meso_ct_label)
pairs_list <- combn(ct_levels, 2, simplify = FALSE)

ttest_rows <- lapply(pairs_list, function(pair) {
  g1 <- cor_df$spearman_r[cor_df$meso_ct_label == pair[1]]
  g2 <- cor_df$spearman_r[cor_df$meso_ct_label == pair[2]]
  if (length(g1) < 2 || length(g2) < 2) return(NULL)
  tt <- t.test(g1, g2)
  data.frame(
    group1   = pair[1],
    group2   = pair[2],
    n1       = length(g1),
    n2       = length(g2),
    mean1    = round(mean(g1), 3),
    mean2    = round(mean(g2), 3),
    t_stat   = round(tt$statistic, 3),
    df       = round(tt$parameter, 1),
    p_value  = tt$p.value,
    stringsAsFactors = FALSE
  )
})
ttest_df <- do.call(rbind, Filter(Negate(is.null), ttest_rows))
ttest_df$p_adj <- p.adjust(ttest_df$p_value, method = "BH")
ttest_df$sig   <- ifelse(ttest_df$p_adj < 0.001, "***",
                  ifelse(ttest_df$p_adj < 0.01,  "**",
                  ifelse(ttest_df$p_adj < 0.05,  "*", "ns")))

write.csv(ttest_df, file.path(out_dir, "MDM_TF_spearman_ttest_pairs.csv"),
          row.names = FALSE, quote = FALSE)

# ---- Boxplot with significance brackets (raw p < 0.05 only) ------------------
sig_raw <- ttest_df[ttest_df$p_value < 0.05, ]
sig_comparisons <- lapply(seq_len(nrow(sig_raw)), function(i)
  c(as.character(sig_raw$group1[i]), as.character(sig_raw$group2[i])))

message(sprintf("  Significant pairs (raw p < 0.05): %d", length(sig_comparisons)))

jitter_pos <- position_jitter(width = 0.15, seed = 42)

p <- ggplot(cor_df, aes(x = meso_ct_label, y = spearman_r)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_boxplot(fill = "#d6604d", colour = "grey20", outlier.shape = NA,
               width = 0.55, linewidth = 0.4, alpha = 0.35) +
  geom_point(colour = "#d6604d", size = 1.8, alpha = 0.8, position = jitter_pos) +
  geom_text(aes(label = n_cells), position = jitter_pos,
            size = 2.5, vjust = -0.7, colour = "grey30") +
  labs(
    x       = NULL,
    y       = "Spearman r (14 MDM TFs)",
    title   = "Per-sample TF profile correlation: meso vs Blood CD14 Mono",
    caption = "Brackets: t-test, raw p; * <0.05  ** <0.01  *** <0.001"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption     = element_text(size = 7, colour = "grey50")
  )

if (length(sig_comparisons) > 0) {
  p <- p + stat_compare_means(
    comparisons  = sig_comparisons,
    method       = "t.test",
    label        = "p.signif",
    tip.length   = 0.01,
    size         = 3,
    bracket.size = 0.3
  )
}

ggsave(file.path(out_dir, "MDM_TF_spearman_per_sample_boxplot.pdf"),
       p, width = 7, height = 6)
message("Saved MDM_TF_spearman_per_sample_boxplot.pdf")

# ---- Print summary -----------------------------------------------------------
message("\nMedian Spearman r per meso cell type (sorted):")
summ <- cor_df %>%
  group_by(meso_ct_label) %>%
  summarise(n = n(), median_r = round(median(spearman_r, na.rm=TRUE), 3),
            mean_r = round(mean(spearman_r, na.rm=TRUE), 3),
            .groups = "drop") %>%
  arrange(desc(median_r))
print(as.data.frame(summ), row.names = FALSE)

message("\nSignificant pairwise t-tests (BH-adjusted p < 0.05):")
sig_out <- ttest_df[ttest_df$p_adj < 0.05, c("group1","group2","n1","n2","t_stat","p_value","p_adj","sig")]
sig_out$p_value <- signif(sig_out$p_value, 3)
sig_out$p_adj   <- signif(sig_out$p_adj,   3)
print(as.data.frame(sig_out), row.names = FALSE)

message("\nDone.")
