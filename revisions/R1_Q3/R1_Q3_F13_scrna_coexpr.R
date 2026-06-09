# R1_Q3 – F13: SOX9 / RUNX2 / SNAI2 co-expression in mesothelioma scRNA-seq
#
# Analyses:
#   1. Feature UMAPs – expression of each gene on the global UMAP
#   2. Co-expression UMAP – joint SOX9+RUNX2+SNAI2 score per cell
#   3. Per-sample feature UMAPs (split by patient)
#   4. Metacell scatter plots – SOX9 vs RUNX2, SOX9 vs SNAI2, RUNX2 vs SNAI2
#      per sample (metacells smooth technical noise for reliable correlations)
#   5. Per-sample correlation heatmap summary
#
# Outputs:
#   F13a_feature_umaps.pdf
#   F13b_coexpr_umap.pdf
#   F13c_per_sample_feature_umaps.pdf
#   F13d_metacell_coexpr_scatter.pdf
#   F13e_per_sample_cor_heatmap.pdf

set.seed(1234)

projdir  <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
scrna_dir <- file.path(projdir, 'tumor_compartment', 'scrna')
outdir    <- file.path(projdir, 'git_repo_claude', 'R1_Q3')
dir.create(file.path(outdir, 'Plots'), showWarnings = FALSE)

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))

library(patchwork)
library(ggpubr)

genes <- c('SOX9', 'RUNX2', 'SNAI2')
pairs <- combn(genes, 2, simplify = FALSE)   # 3 pairwise combinations

# Samples ordered by sarcomatoid score (low → high)
sarc_order <- read.csv(file.path(scrna_dir, 'cnmf20_sarcomatoid_sample_order.csv'),
                        row.names = 1)
tumor_sams <- sarc_order$sampleID[!sarc_order$sampleID %in% c('HU37','HU62','normal_pleura')]
sample_levels <- tumor_sams[order(sarc_order$x[sarc_order$sampleID %in% tumor_sams])]

pal_sarc <- colorRampPalette(c('#2166AC','#FFFFBF','#D6604D'))(length(sample_levels))
palette_sample <- setNames(pal_sarc, sample_levels)

##############################################################################
# Load data
##############################################################################
message('Loading srt...')
srt <- readRDS(file.path(scrna_dir, 'srt.rds'))

# Filter to tumor samples only
srt <- srt[, srt$sampleID %in% tumor_sams]
srt$sampleID <- factor(srt$sampleID, levels = sample_levels)

# Identify UMAP reduction name
umap_name <- grep('umap|UMAP', names(srt@reductions), ignore.case = TRUE, value = TRUE)[1]
message('Using reduction: ', umap_name)

# Check genes present
genes_present <- genes[genes %in% rownames(srt)]
message('Genes found: ', paste(genes_present, collapse = ', '))

##############################################################################
# F13a – Feature UMAPs: one panel per gene
##############################################################################
feat_plots <- FeaturePlot(srt,
  features   = genes_present,
  reduction  = umap_name,
  order      = TRUE,
  pt.size    = 0.3,
  cols       = c('lightgrey', '#D6604D'),
  combine    = FALSE
)
feat_plots <- lapply(feat_plots, function(p)
  p + gtheme + theme(legend.position = 'right'))

pdf(file.path(outdir, 'Plots', 'F13a_feature_umaps.pdf'), width = 5 * length(genes_present), height = 5)
print(wrap_plots(feat_plots, ncol = length(genes_present)))
dev.off()

##############################################################################
# F13b – Co-expression UMAP: mean z-score of all three genes per cell
##############################################################################
expr_mat <- as.matrix(srt@assays$RNA@data[genes_present, , drop = FALSE])
expr_z   <- t(apply(expr_mat, 1, function(x) {
  s <- sd(x); if (s == 0) return(x); (x - mean(x)) / s
}))
srt$coexpr_score <- colMeans(expr_z, na.rm = TRUE)

coexpr_umap <- FeaturePlot(srt,
  features  = 'coexpr_score',
  reduction = umap_name,
  order     = TRUE,
  pt.size   = 0.3,
  cols      = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100))
) + labs(title = 'SOX9 / RUNX2 / SNAI2 co-expression score') + gtheme

pdf(file.path(outdir, 'Plots', 'F13b_coexpr_umap.pdf'), width = 5, height = 5)
print(coexpr_umap)
dev.off()

##############################################################################
# F13c – Per-sample feature UMAPs (split by patient, ordered by scS-score)
##############################################################################
pdf(file.path(outdir, 'Plots', 'F13c_per_sample_feature_umaps.pdf'),
    width = 4 * length(genes_present), height = 4 * length(sample_levels))

for (sam in sample_levels) {
  srt_sub <- srt[, srt$sampleID == sam]
  if (ncol(srt_sub) < 20) next
  p_list <- FeaturePlot(srt_sub,
    features  = genes_present,
    reduction = umap_name,
    order     = TRUE,
    pt.size   = 0.5,
    cols      = c('lightgrey', '#D6604D'),
    combine   = FALSE
  )
  p_list <- lapply(p_list, function(p)
    p + gtheme + labs(subtitle = sam) + theme(legend.position = 'right'))
  print(wrap_plots(p_list, ncol = length(genes_present)) +
    plot_annotation(title = sam))
}
dev.off()

##############################################################################
# Load metacells
##############################################################################
message('Loading metacells...')
mc <- readRDS(file.path(scrna_dir, 'metacells.rds'))

# Filter to tumor samples
sam_col <- grep('sample|Sample|sampleID', colnames(mc@meta.data), value = TRUE)[1]
mc <- mc[, mc@meta.data[[sam_col]] %in% tumor_sams]
mc@meta.data$sampleID <- factor(mc@meta.data[[sam_col]], levels = sample_levels)

mc_genes <- genes[genes %in% rownames(mc)]
message('Metacell genes found: ', paste(mc_genes, collapse = ', '))

# Extract log-normalized expression
mc_expr <- as.data.frame(t(as.matrix(mc@assays$RNA@data[mc_genes, , drop = FALSE])))
mc_expr$sampleID <- mc@meta.data$sampleID

##############################################################################
# F13d – Metacell pairwise scatter per sample
##############################################################################
make_mc_scatter <- function(gx, gy, df) {
  ggplot(df, aes_string(x = gx, y = gy, color = 'sampleID')) +
    geom_point(size = 1.2, alpha = 0.7) +
    geom_smooth(aes(group = sampleID, color = sampleID),
      method = 'lm', se = FALSE, linewidth = 0.5) +
    scale_color_manual(values = palette_sample) +
    stat_cor(aes(group = 1), method = 'pearson',
      label.x.npc = 'left', label.y.npc = 'top', size = 3,
      color = 'black') +
    labs(x = paste0(gx, ' (metacell)'), y = paste0(gy, ' (metacell)'),
         title = paste0(gx, ' vs ', gy, ' – all samples'),
         color = 'Sample (low→high scS)') +
    gtheme + theme(legend.position = 'right')
}

# Global scatter (all samples together, colored by sample)
global_scatter <- lapply(pairs, function(p) make_mc_scatter(p[1], p[2], mc_expr))

pdf(file.path(outdir, 'Plots', 'F13d_metacell_coexpr_scatter.pdf'),
    width = 6 * length(pairs), height = 6)
print(wrap_plots(global_scatter, ncol = length(pairs)))
dev.off()

# Per-sample scatter (one row per sample, one column per pair)
pdf(file.path(outdir, 'Plots', 'F13d2_metacell_coexpr_per_sample.pdf'),
    width = 4 * length(pairs), height = 4 * length(sample_levels))

for (sam in sample_levels) {
  df_sam <- mc_expr[mc_expr$sampleID == sam, ]
  if (nrow(df_sam) < 5) next
  p_list <- lapply(pairs, function(p) {
    gx <- p[1]; gy <- p[2]
    if (!all(c(gx, gy) %in% colnames(df_sam))) return(NULL)
    ggplot(df_sam, aes_string(x = gx, y = gy)) +
      geom_point(size = 1.5, alpha = 0.8,
        color = palette_sample[as.character(sam)]) +
      geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.5) +
      stat_cor(method = 'pearson', size = 3,
        label.x.npc = 'left', label.y.npc = 'top') +
      labs(x = gx, y = gy, title = paste0(sam, ': ', gx, ' vs ', gy)) +
      gtheme
  })
  p_list <- Filter(Negate(is.null), p_list)
  if (length(p_list) > 0)
    print(wrap_plots(p_list, ncol = length(pairs)) +
      plot_annotation(title = sam))
}
dev.off()

##############################################################################
# F13e – Per-sample Pearson correlation summary heatmap
##############################################################################
cor_per_sam <- do.call(rbind, lapply(sample_levels, function(sam) {
  df_sam <- mc_expr[mc_expr$sampleID == sam, mc_genes, drop = FALSE]
  if (nrow(df_sam) < 5) return(NULL)
  cor_mat <- cor(df_sam, method = 'pearson', use = 'complete.obs')
  pair_names <- sapply(pairs, function(p) paste0(p[1], '_', p[2]))
  pair_cors  <- sapply(pairs, function(p) {
    if (!all(p %in% rownames(cor_mat))) return(NA)
    cor_mat[p[1], p[2]]
  })
  data.frame(sample = sam, t(setNames(pair_cors, pair_names)))
}))

cor_long <- pivot_longer(cor_per_sam, -sample,
  names_to = 'pair', values_to = 'r')
cor_long$sample <- factor(cor_long$sample, levels = sample_levels)
cor_long$pair   <- gsub('_', ' vs ', cor_long$pair)

f13e <- ggplot(cor_long, aes(x = pair, y = sample, fill = r)) +
  geom_tile(color = 'white', linewidth = 0.5) +
  geom_text(aes(label = sprintf('%.2f', r)), size = 3) +
  scale_fill_gradientn(
    colours = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100)),
    limits = c(-1, 1), name = 'Pearson r') +
  labs(x = NULL, y = 'Sample (ordered by scS-score)',
       title = 'Per-sample co-expression (metacells)') +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(face = 'italic'))

pdf(file.path(outdir, 'Plots', 'F13e_per_sample_cor_heatmap.pdf'), width = 6, height = 5)
print(f13e)
dev.off()

message('F13 done.')
