# R1_Q3 – F14: SOX9 SCENIC regulon AUC vs scS-score per cell
#
# The SCENIC AUC matrix contains per-cell SOX9(+) regulon activity.
# We correlate it with the scS-score (cNMF20 module), computed from the same
# expression matrix used for SCENIC (expr_vg_5000.txt).
# Tumor samples: p11, p12, p13, p811, p826, p848
#
# Also plots all 46 regulons vs scS-score to contextualize SOX9
#
# Outputs:
#   F14a_SOX9_regulon_vs_scS_scatter.pdf   – scatter + per-sample
#   F14b_all_regulons_vs_scS_heatmap.pdf   – heatmap: regulon x sample correlation
#   F14c_SOX9_regulon_umap_proxy.pdf       – violin per sample (SOX9 AUC + scS)
#   F14_regulon_scS_correlation.csv        – full correlation table

set.seed(1234)

projdir   <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
scenic_dir <- file.path(projdir, 'tumor_compartment', 'scrna', 'SCENIC',
                         'vg_5000_mw_tss500bp', 'tumor_programs')
scrna_dir  <- file.path(projdir, 'tumor_compartment', 'scrna')
outdir     <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))

library(patchwork)
library(ggpubr)
library(tidyr)
library(dplyr)

# Tumor samples in SCENIC
tumor_prefixes <- c('p11', 'p12', 'p13', 'p811', 'p826', 'p848')

# Sample color palette (use a sarcomatoid gradient)
pal_samples <- setNames(
  paletteer::paletteer_c("grDevices::RdYlBu", length(tumor_prefixes)),
  tumor_prefixes
)

##############################################################################
# 1. Load SCENIC AUC matrix
##############################################################################
message('Loading AUC matrix...')
auc <- read.csv(file.path(scenic_dir, 'auc_mtx.csv'), row.names = 1,
                check.names = FALSE)
colnames(auc) <- gsub('\\(\\+\\)', '', colnames(auc))
colnames(auc) <- gsub('\\(\\-\\)', '', colnames(auc))

# Filter to tumor cells only
auc$sample <- sapply(rownames(auc), function(x) strsplit(x, '_')[[1]][1])
auc <- auc[auc$sample %in% tumor_prefixes, ]
message('Tumor cells in AUC: ', nrow(auc))
message('Regulons: ', paste(colnames(auc)[colnames(auc) != 'sample'], collapse=', '))

##############################################################################
# 2. Compute scS-score per cell from expression matrix
##############################################################################
message('Loading cNMF20 gene list...')
cnmf_rds <- file.path(scrna_dir, 'cnmf_genelist_25_nfeat_5000.rds')
cnmf_spectra <- readRDS(cnmf_rds)
sarc_genes   <- head(cnmf_spectra[['cNMF20']], 20)

message('Loading expression matrix (may take a minute)...')
expr <- read.table(file.path(scenic_dir, 'expr_vg_5000.txt'),
                   header = TRUE, row.names = 1, sep = '\t', check.names = FALSE)

# Filter to tumor cells
expr <- expr[rownames(expr) %in% rownames(auc), , drop = FALSE]
message('Tumor cells in expr: ', nrow(expr))

# Compute scS-score per cell = mean expression of cNMF20 top genes present
sarc_genes_present <- sarc_genes[sarc_genes %in% colnames(expr)]
message('cNMF20 genes present in expr: ', length(sarc_genes_present))
scS_score <- rowMeans(expr[, sarc_genes_present, drop = FALSE], na.rm = TRUE)

##############################################################################
# 3. Build combined data frame
##############################################################################
common_cells <- intersect(rownames(auc), names(scS_score))
df <- data.frame(
  cell     = common_cells,
  scS      = scS_score[common_cells],
  SOX9_auc = auc[common_cells, 'SOX9'],
  sample   = auc[common_cells, 'sample'],
  stringsAsFactors = FALSE
)
df$sample <- factor(df$sample, levels = tumor_prefixes)
message('Cells in combined df: ', nrow(df))

##############################################################################
# F14a – SOX9 regulon AUC vs scS-score: global + per sample
##############################################################################
global_cor <- cor(df$SOX9_auc, df$scS, method = 'spearman', use = 'complete.obs')
message('Global Spearman r (SOX9 AUC vs scS-score): ', round(global_cor, 3))

p_global <- ggplot(df, aes(x = scS, y = SOX9_auc, color = sample)) +
  geom_point(size = 0.5, alpha = 0.4) +
  geom_smooth(aes(group = 1), method = 'lm', se = TRUE,
    color = 'black', linewidth = 0.7) +
  stat_cor(aes(group = 1), method = 'spearman', size = 3.5,
    label.x.npc = 'left', label.y.npc = 'top', color = 'black') +
  scale_color_manual(values = pal_samples) +
  labs(x = 'scS-score (cNMF20)', y = 'SOX9 regulon AUC',
       title = 'SOX9 SCENIC regulon vs scS-score',
       color = 'Sample') +
  gtheme

# Per-sample scatter
per_sam_plots <- lapply(tumor_prefixes, function(sam) {
  d <- df[df$sample == sam, ]
  if (nrow(d) < 10) return(NULL)
  ggplot(d, aes(x = scS, y = SOX9_auc)) +
    geom_point(size = 0.6, alpha = 0.5, color = pal_samples[sam]) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.5) +
    stat_cor(method = 'spearman', size = 3, label.x.npc = 'left', label.y.npc = 'top') +
    labs(x = 'scS-score', y = 'SOX9 AUC', title = sam) +
    gtheme
})
per_sam_plots <- Filter(Negate(is.null), per_sam_plots)

pdf(file.path(outdir, 'Plots', 'F14a_SOX9_regulon_vs_scS_scatter.pdf'),
    width = 5 * (1 + ceiling(length(per_sam_plots) / 2)), height = 8)
print(wrap_plots(c(list(p_global), per_sam_plots),
                 ncol = 1 + ceiling(length(per_sam_plots) / 2)))
dev.off()

##############################################################################
# F14b – All regulons vs scS-score: per-sample correlation heatmap
##############################################################################
regulon_cols <- colnames(auc)[colnames(auc) != 'sample']

cor_all <- do.call(rbind, lapply(tumor_prefixes, function(sam) {
  idx  <- df$sample == sam
  if (sum(idx) < 10) return(NULL)
  scS_s <- df$scS[idx]
  sapply_df <- sapply(regulon_cols, function(reg) {
    auc_s <- auc[common_cells[idx], reg]
    if (var(auc_s, na.rm = TRUE) == 0) return(NA)
    cor(scS_s, auc_s, method = 'spearman', use = 'complete.obs')
  })
  data.frame(sample = sam, t(sapply_df), check.names = FALSE)
}))

# Also compute global correlations
global_cors <- sapply(regulon_cols, function(reg) {
  auc_v <- auc[common_cells, reg]
  if (var(auc_v, na.rm = TRUE) == 0) return(NA)
  cor(df$scS, auc_v, method = 'spearman', use = 'complete.obs')
})
global_row <- data.frame(sample = 'All', t(global_cors), check.names = FALSE)
cor_all <- rbind(global_row, cor_all)

# Order regulons by global correlation
reg_order <- names(sort(global_cors, decreasing = TRUE, na.last = TRUE))
cor_long <- pivot_longer(cor_all, -sample, names_to = 'regulon', values_to = 'rho')
cor_long$regulon <- factor(cor_long$regulon, levels = reg_order)
cor_long$sample  <- factor(cor_long$sample,
  levels = c('All', tumor_prefixes))

f14b <- ggplot(cor_long, aes(x = sample, y = regulon, fill = rho)) +
  geom_tile(color = 'white', linewidth = 0.3) +
  scale_fill_gradientn(
    colours = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100)),
    limits = c(-1, 1), na.value = 'grey80', name = 'Spearman r') +
  geom_text(data = cor_long[cor_long$sample == 'All' & !is.na(cor_long$rho) & abs(cor_long$rho) > 0.1, ],
    aes(label = sprintf('%.2f', rho)), size = 2.5) +
  labs(x = NULL, y = NULL,
       title = 'SCENIC regulon AUC vs scS-score\n(Spearman r per sample)') +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(angle = 30, hjust = 1))

pdf(file.path(outdir, 'Plots', 'F14b_all_regulons_vs_scS_heatmap.pdf'),
    width = 5, height = 8)
print(f14b)
dev.off()

##############################################################################
# F14c – Violin: SOX9 AUC per sample (ordered by mean scS-score)
##############################################################################
# Order samples by mean scS-score
sample_means <- df %>% group_by(sample) %>%
  summarise(mean_scS = mean(scS), mean_auc = mean(SOX9_auc), .groups = 'drop') %>%
  arrange(mean_scS)
df$sample_ord <- factor(df$sample, levels = sample_means$sample)

v1 <- ggplot(df, aes(x = sample_ord, y = SOX9_auc, fill = sample_ord)) +
  geom_violin(trim = TRUE, scale = 'width', linewidth = 0.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, linewidth = 0.3, fill = 'white', alpha = 0.8) +
  scale_fill_manual(values = pal_samples) +
  labs(x = 'Sample (ordered by mean scS-score)', y = 'SOX9 regulon AUC',
       title = 'SOX9 regulon activity per sample') +
  gtheme + NoLegend()

v2 <- ggplot(df, aes(x = sample_ord, y = scS, fill = sample_ord)) +
  geom_violin(trim = TRUE, scale = 'width', linewidth = 0.2, alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA, linewidth = 0.3, fill = 'white', alpha = 0.8) +
  scale_fill_manual(values = pal_samples) +
  labs(x = 'Sample (ordered by mean scS-score)', y = 'scS-score (cNMF20)',
       title = 'scS-score per sample') +
  gtheme + NoLegend()

pdf(file.path(outdir, 'Plots', 'F14c_SOX9_regulon_scS_violins.pdf'), width = 8, height = 7)
print(wrap_plots(v1, v2, ncol = 1))
dev.off()

##############################################################################
# Save full correlation table
##############################################################################
write.csv(cor_all, file.path(outdir, 'F14_regulon_scS_correlation.csv'),
          row.names = FALSE)

message('F14 done. Global SOX9 AUC vs scS-score Spearman r = ', round(global_cor, 3))
