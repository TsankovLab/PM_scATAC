# F13 metacell-only: run after F13a-c already succeeded
set.seed(1234)

projdir  <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
scrna_dir <- file.path(projdir, 'tumor_compartment', 'scrna')
outdir    <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))

library(patchwork)
library(ggpubr)
library(tidyr)

genes <- c('SOX9', 'RUNX2', 'SNAI2')
pairs <- combn(genes, 2, simplify = FALSE)

sarc_order <- read.csv(file.path(scrna_dir, 'cnmf20_sarcomatoid_sample_order.csv'),
                        row.names = 1)
tumor_sams    <- sarc_order$sampleID[!sarc_order$sampleID %in% c('HU37','HU62','normal_pleura')]
sample_levels <- tumor_sams[order(sarc_order$x[sarc_order$sampleID %in% tumor_sams])]
palette_sample <- setNames(
  colorRampPalette(c('#2166AC','#FFFFBF','#D6604D'))(length(sample_levels)),
  sample_levels)

message('Loading metacells...')
mc <- readRDS(file.path(scrna_dir, 'metacells.rds'))

# Print class and slot info for debugging
cat('Metacell class:', class(mc), '\n')
cat('Meta cols:', paste(head(colnames(mc@meta.data), 15), collapse = ', '), '\n')
cat('Assay class:', class(mc@assays$RNA), '\n')

sam_col <- grep('sample|Sample|sampleID', colnames(mc@meta.data), value = TRUE)[1]
cat('Sample column:', sam_col, '\n')
cat('Sample values:', paste(head(unique(mc@meta.data[[sam_col]]), 5), collapse = ', '), '\n')

mc <- mc[, mc@meta.data[[sam_col]] %in% tumor_sams]
mc@meta.data$sampleID <- factor(mc@meta.data[[sam_col]], levels = sample_levels)

mc_genes <- genes[genes %in% rownames(mc)]
cat('Genes in metacells:', paste(mc_genes, collapse = ', '), '\n')

# Seurat v5 (Assay5) vs v4 (@data) compatibility
mc_mat <- tryCatch(
  GetAssayData(mc, assay = 'RNA', layer = 'data')[mc_genes, , drop = FALSE],
  error = function(e) {
    message('GetAssayData layer failed, trying slot...')
    tryCatch(
      GetAssayData(mc, assay = 'RNA', slot = 'data')[mc_genes, , drop = FALSE],
      error = function(e2) mc@assays$RNA@data[mc_genes, , drop = FALSE]
    )
  }
)
mc_expr <- as.data.frame(t(as.matrix(mc_mat)))
mc_expr$sampleID <- mc@meta.data$sampleID
cat('Metacell expression dim:', dim(mc_expr), '\n')

##############################################################################
# F13d – Global metacell scatter (all samples, colored by sample)
##############################################################################
make_mc_scatter <- function(gx, gy, df) {
  ggplot(df, aes_string(x = gx, y = gy, color = 'sampleID')) +
    geom_point(size = 1.2, alpha = 0.7) +
    geom_smooth(aes(group = sampleID, color = sampleID),
      method = 'lm', se = FALSE, linewidth = 0.5) +
    scale_color_manual(values = palette_sample) +
    stat_cor(aes(group = 1), method = 'pearson',
      label.x.npc = 'left', label.y.npc = 'top', size = 3, color = 'black') +
    labs(x = paste0(gx, ' (metacell)'), y = paste0(gy, ' (metacell)'),
         title = paste0(gx, ' vs ', gy),
         color = 'Sample') +
    gtheme + theme(legend.position = 'right')
}

global_scatter <- lapply(pairs, function(p) make_mc_scatter(p[1], p[2], mc_expr))

pdf(file.path(outdir, 'Plots', 'F13d_metacell_coexpr_scatter.pdf'),
    width = 6 * length(pairs), height = 6)
print(wrap_plots(global_scatter, ncol = length(pairs)))
dev.off()

##############################################################################
# F13d2 – Per-sample scatter
##############################################################################
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
# F13e – Per-sample correlation heatmap
##############################################################################
cor_per_sam <- do.call(rbind, lapply(sample_levels, function(sam) {
  df_sam <- mc_expr[mc_expr$sampleID == sam, mc_genes, drop = FALSE]
  if (nrow(df_sam) < 5) return(NULL)
  cor_mat    <- cor(df_sam, method = 'pearson', use = 'complete.obs')
  pair_names <- sapply(pairs, function(p) paste0(p[1], '_', p[2]))
  pair_cors  <- sapply(pairs, function(p) {
    if (!all(p %in% rownames(cor_mat))) return(NA)
    cor_mat[p[1], p[2]]
  })
  data.frame(sample = sam, t(setNames(pair_cors, pair_names)))
}))

cor_long <- pivot_longer(cor_per_sam, -sample, names_to = 'pair', values_to = 'r')
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
  theme(panel.grid = element_blank(), axis.text.x = element_text(face = 'italic'))

pdf(file.path(outdir, 'Plots', 'F13e_per_sample_cor_heatmap.pdf'), width = 6, height = 5)
print(f13e)
dev.off()

message('F13 metacell done.')
