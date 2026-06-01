# plot_allTF_spearman_per_sample.R  (R1_Q4)
#
# Per-sample Spearman correlations of ALL chromVAR TF z-scores between:
#   - Each meso cell type (per meso sample, >= 50 cells required)
#   - Blood CD14 Mono averaged across all PBMC cells (stable reference)
#
# Identical pipeline to plot_MDM_TF_spearman_per_sample.R but uses the full
# MotifMatrix (all TF motifs) instead of the 14 MDM-specific TFs.
#
# Output: all_TF_spearman_per_sample_boxplot.pdf

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(SummarizedExperiment)
  library(ggpubr)
})
addArchRGenome("hg38")
addArchRThreads(4)

proj_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/meso_vs_pbmc"
out_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4"

MIN_CELLS <- 50

message("Loading project...")
proj <- loadArchRProject(proj_dir, showLogo = FALSE)

ct   <- as.character(proj$celltype_unified)
tiss <- as.character(proj$tissue)
proj$tissue_ct <- ifelse(
  ct %in% c("CD14_Mono", "CD16_Mono", "DC"),
  paste0(tiss, "_", ct),
  ct
)
proj$sample_tct <- paste0(proj$Sample, "__", proj$tissue_ct)

# ---- Average MotifMatrix per sample × cell type ------------------------------
message("Averaging MotifMatrix per sample x tissue_ct...")
devSE <- getGroupSE(proj,
                    useMatrix = "MotifMatrix",
                    groupBy   = "sample_tct",
                    divideN   = TRUE,
                    scaleTo   = NULL)

# Use full motif names to avoid duplicate rownames across TF families
rownames(devSE) <- rowData(devSE)$name
n_motifs <- nrow(devSE)
message(sprintf("  Total motifs in MotifMatrix: %d", n_motifs))

mat <- assays(devSE)[[1]]

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

message(sprintf("  Groups with < %d cells (excluded): %d",
                MIN_CELLS, sum(col_meta$n_cells < MIN_CELLS)))

col_meta_pass <- col_meta[col_meta$n_cells >= MIN_CELLS, ]
message(sprintf("  Groups passing filter: %d / %d",
                nrow(col_meta_pass), nrow(col_meta)))

# ---- Blood reference: weighted average across PBMC_CD14_Mono samples ---------
blood_types <- "PBMC_CD14_Mono"
blood_cols  <- col_meta_pass$col_name[col_meta_pass$tissue_ct == blood_types]
blood_wts   <- col_meta_pass$n_cells[col_meta_pass$tissue_ct == blood_types]
blood_ref   <- rowSums(mat[, blood_cols, drop=FALSE] * rep(blood_wts, each=nrow(mat))) / sum(blood_wts)
message(sprintf("  Blood reference built from %d PBMC_CD14_Mono groups", length(blood_cols)))

# ---- Per-meso-sample Spearman correlations -----------------------------------
meso_types <- c("Mesothelioma_CD14_Mono", "Mesothelioma_CD16_Mono",
                "Mesothelioma_DC", "TAM_CXCLs", "TAM_interstitial",
                "TAM_MARCO", "TAM_TREM2")

meso_meta <- col_meta_pass[col_meta_pass$tissue_ct %in% meso_types, ]
message(sprintf("  Meso sample x cell type groups: %d", nrow(meso_meta)))

cor_rows <- lapply(seq_len(nrow(meso_meta)), function(i) {
  meso_prof <- mat[, meso_meta$col_name[i]]
  valid     <- !is.na(meso_prof) & !is.na(blood_ref)
  if (sum(valid) < 10) return(NULL)
  r <- cor(meso_prof[valid], blood_ref[valid], method = "spearman")
  data.frame(
    meso_sample = meso_meta$sample[i],
    meso_ct     = meso_meta$tissue_ct[i],
    spearman_r  = r,
    n_cells     = meso_meta$n_cells[i],
    stringsAsFactors = FALSE
  )
})
cor_df <- do.call(rbind, Filter(Negate(is.null), cor_rows))

# Clean labels
cor_df$meso_ct_label <- gsub("Mesothelioma_", "Meso ", cor_df$meso_ct)
cor_df$meso_ct_label <- gsub("_", " ", cor_df$meso_ct_label)

ct_median_order <- cor_df %>%
  group_by(meso_ct_label) %>%
  summarise(med = median(spearman_r, na.rm=TRUE), .groups="drop") %>%
  arrange(desc(med))
cor_df$meso_ct_label <- factor(cor_df$meso_ct_label, levels=ct_median_order$meso_ct_label)

write.csv(cor_df, file.path(out_dir, "all_TF_spearman_per_sample.csv"),
          row.names=FALSE, quote=FALSE)

# ---- Pairwise t-tests --------------------------------------------------------
ct_levels <- levels(cor_df$meso_ct_label)
pairs_list <- combn(ct_levels, 2, simplify=FALSE)

ttest_rows <- lapply(pairs_list, function(pair) {
  g1 <- cor_df$spearman_r[cor_df$meso_ct_label == pair[1]]
  g2 <- cor_df$spearman_r[cor_df$meso_ct_label == pair[2]]
  if (length(g1) < 2 || length(g2) < 2) return(NULL)
  tt <- t.test(g1, g2)
  data.frame(group1=pair[1], group2=pair[2],
             n1=length(g1), n2=length(g2),
             mean1=round(mean(g1),3), mean2=round(mean(g2),3),
             t_stat=round(tt$statistic,3), df=round(tt$parameter,1),
             p_value=tt$p.value, stringsAsFactors=FALSE)
})
ttest_df <- do.call(rbind, Filter(Negate(is.null), ttest_rows))
ttest_df$p_adj <- p.adjust(ttest_df$p_value, method="BH")
ttest_df$sig   <- ifelse(ttest_df$p_adj < 0.001, "***",
                  ifelse(ttest_df$p_adj < 0.01,  "**",
                  ifelse(ttest_df$p_adj < 0.05,  "*", "ns")))

write.csv(ttest_df, file.path(out_dir, "all_TF_spearman_ttest_pairs.csv"),
          row.names=FALSE, quote=FALSE)

# ---- Boxplot -----------------------------------------------------------------
sig_raw <- ttest_df[ttest_df$p_value < 0.05, ]
sig_comparisons <- lapply(seq_len(nrow(sig_raw)), function(i)
  c(as.character(sig_raw$group1[i]), as.character(sig_raw$group2[i])))

message(sprintf("  Significant pairs (raw p < 0.05): %d", length(sig_comparisons)))

jitter_pos <- position_jitter(width=0.15, seed=42)

p <- ggplot(cor_df, aes(x=meso_ct_label, y=spearman_r)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey50", linewidth=0.4) +
  geom_boxplot(fill="#4DBBD5", colour="grey20", outlier.shape=NA,
               width=0.55, linewidth=0.4, alpha=0.35) +
  geom_point(colour="#4DBBD5", size=1.8, alpha=0.8, position=jitter_pos) +
  geom_text(aes(label=n_cells), position=jitter_pos,
            size=2.5, vjust=-0.7, colour="grey30") +
  labs(
    x       = NULL,
    y       = sprintf("Spearman r (%d TF motifs)", n_motifs),
    title   = "Per-sample TF profile correlation: meso vs Blood CD14 Mono",
    subtitle = "All chromVAR TF motifs",
    caption = "Brackets: t-test, raw p; * <0.05  ** <0.01  *** <0.001"
  ) +
  theme_bw(base_size=11) +
  theme(
    axis.text.x        = element_text(angle=45, hjust=1, size=9),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption       = element_text(size=7, colour="grey50")
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

ggsave(file.path(out_dir, "all_TF_spearman_per_sample_boxplot.pdf"),
       p, width=7, height=6)
message("Saved all_TF_spearman_per_sample_boxplot.pdf")

# ---- Summary -----------------------------------------------------------------
message("\nMedian Spearman r per meso cell type (sorted):")
summ <- cor_df %>%
  group_by(meso_ct_label) %>%
  summarise(n=n(), median_r=round(median(spearman_r,na.rm=TRUE),3),
            mean_r=round(mean(spearman_r,na.rm=TRUE),3), .groups="drop") %>%
  arrange(desc(median_r))
print(as.data.frame(summ), row.names=FALSE)

message("\nDone.")
