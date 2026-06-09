# R1_Q3 addon – Bueno cohort: axis score by pathological histology + % sarcomatoid
#
# Uses HistologyReduced (Epithelioid / Biphasic / Sarcomatoid) from bueno_meta_cli.rds,
# which is the actual pathologic diagnosis — less granular and more clinically direct
# than the 4-class molecular ConsensusCluster subtype already shown in F1-F3.
# Also adds axis score vs % sarcomatoid component scatter (mirrors F4 for Mesomics).
#
# Outputs:
#   F9_bueno_axis_by_path_histology.pdf  – boxplots by pathological subtype
#   F9b_bueno_axis_vs_sarc_pct.pdf       – scatter vs % sarcomatoid component

set.seed(1234)

projdir  <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
bulk_dir <- file.path(projdir, 'bulkRNA_meso')
outdir   <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))

library(ggpubr)
library(rstatix)
library(patchwork)

axis_tfs <- c('SOX9', 'RUNX1', 'RUNX2', 'SNAI2')

hist_palette <- setNames(
  as.character(paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1, 3, 7)]),
  c('Epithelioid', 'Biphasic', 'Sarcomatoid')
)

# ---- Load data -----------------------------------------------------------
bueno_expr     <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))[['bueno']]
bueno_meta_mol <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))[['bueno']]
bueno_meta_cli <- readRDS(file.path(bulk_dir, 'bueno_meta_cli.rds'))

# Compute composite regulatory axis score (same z-score approach as F3)
genes_present <- axis_tfs[axis_tfs %in% rownames(bueno_expr)]
axis_mat      <- scale(t(bueno_expr[genes_present, , drop = FALSE]))  # samples x genes
axis_score    <- rowMeans(axis_mat, na.rm = TRUE)

# Build data frame with pathological histology from bueno_meta_cli
# Note: X..Sarcomatoid contains ">99" as a string — coerce to numeric (>99 → 99.5)
clean_pct <- function(x) as.numeric(sub('^>', '', as.character(x)))

bueno_df <- data.frame(
  axis_score     = axis_score[rownames(bueno_meta_cli)],
  mol_subtype    = bueno_meta_cli$ConsensusCluster,
  path_histology = bueno_meta_cli$HistologyReduced,
  sarc_pct       = clean_pct(bueno_meta_cli$X..Sarcomatoid),
  epit_pct       = clean_pct(bueno_meta_cli$X..Epithelioid),
  row.names      = rownames(bueno_meta_cli)
)
# Drop Desmoplastic (n=1) and missing
bueno_df <- bueno_df[!is.na(bueno_df$path_histology) &
                       bueno_df$path_histology != 'Desmoplastic', ]
bueno_df$path_histology <- factor(bueno_df$path_histology,
  levels = c('Epithelioid', 'Biphasic', 'Sarcomatoid'))

##############################################################################
# F9 – Composite axis score by pathological histology (Bueno)
##############################################################################

# Individual gene boxplots
gene_plots <- lapply(axis_tfs, function(gene) {
  if (!gene %in% rownames(bueno_expr)) return(NULL)
  df <- data.frame(
    expression    = as.numeric(bueno_expr[gene, rownames(bueno_df)]),
    path_histology = bueno_df$path_histology
  )
  stat.test <- tryCatch(df %>%
    t_test(expression ~ path_histology) %>% adjust_pvalue(method = 'none') %>%
    add_significance() %>% add_xy_position(x = 'path_histology', step.increase = 0.45),
    error = function(e) NULL)
  if (!is.null(stat.test))
    stat.test <- stat.test[
      stat.test$group1 %in% c('Epithelioid', 'Sarcomatoid') &
      stat.test$group2 %in% c('Epithelioid', 'Sarcomatoid'), ]

  p <- ggplot(df, aes(x = path_histology, y = expression)) +
    geom_boxplot(aes(fill = path_histology),
      outlier.shape = 16, outlier.size = 1, outlier.alpha = 0.3,
      lwd = 0.3, alpha = 0.8, width = 0.55) +
    scale_fill_manual(values = hist_palette) +
    labs(title = gene, y = 'log2 expression', x = NULL) +
    gtheme + NoLegend()
  if (!is.null(stat.test) && nrow(stat.test) > 0)
    p <- p + stat_pvalue_manual(stat.test, label = 'p.adj.signif',
      hide.ns = TRUE, bracket.nudge.y = 0)
  p
})
gene_plots <- Filter(Negate(is.null), gene_plots)

# Composite axis score boxplot
stat_axis <- tryCatch(bueno_df %>%
  t_test(axis_score ~ path_histology) %>% adjust_pvalue(method = 'none') %>%
  add_significance() %>% add_xy_position(x = 'path_histology', step.increase = 0.45),
  error = function(e) NULL)
if (!is.null(stat_axis))
  stat_axis <- stat_axis[
    stat_axis$group1 %in% c('Epithelioid', 'Sarcomatoid') &
    stat_axis$group2 %in% c('Epithelioid', 'Sarcomatoid'), ]

p_composite <- ggplot(bueno_df, aes(x = path_histology, y = axis_score)) +
  geom_boxplot(aes(fill = path_histology),
    outlier.shape = 16, outlier.size = 1, outlier.alpha = 0.3,
    lwd = 0.3, alpha = 0.8, width = 0.55) +
  scale_fill_manual(values = hist_palette) +
  labs(title = 'Composite axis score',
       y = 'Regulatory axis score\n(mean z-score SOX9/RUNX1/RUNX2/SNAI2)',
       x = NULL) +
  gtheme + NoLegend()
if (!is.null(stat_axis) && nrow(stat_axis) > 0)
  p_composite <- p_composite + stat_pvalue_manual(stat_axis, label = 'p.adj.signif',
    hide.ns = TRUE, bracket.nudge.y = 0)

pdf(file.path(outdir, 'Plots', 'F9_bueno_axis_by_path_histology.pdf'), width = 10, height = 4)
print(wrap_plots(c(gene_plots, list(p_composite)), ncol = length(axis_tfs) + 1))
dev.off()

##############################################################################
# F9b – Composite axis score vs % sarcomatoid component (Bueno) – standalone
##############################################################################
bueno_sarc <- bueno_df[!is.na(bueno_df$sarc_pct), ]
bueno_epit <- bueno_df[!is.na(bueno_df$epit_pct), ]

f9b_sarc <- ggplot(bueno_sarc, aes(x = sarc_pct, y = axis_score)) +
  geom_point(aes(fill = path_histology), shape = 21, size = 2.5, alpha = 0.8, stroke = 0.5) +
  geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.6) +
  stat_cor(method = 'pearson', label.x.npc = 'left', label.y.npc = 'top', size = 3.5) +
  scale_fill_manual(values = hist_palette) +
  labs(x = '% Sarcomatoid component (Bueno et al.)',
       y = 'Regulatory axis score\n(SOX9/RUNX1/RUNX2/SNAI2)',
       fill = 'Histology',
       title = 'Axis score vs % sarcomatoid') +
  gtheme

f9b_epit <- ggplot(bueno_epit, aes(x = epit_pct, y = axis_score)) +
  geom_point(aes(fill = path_histology), shape = 21, size = 2.5, alpha = 0.8, stroke = 0.5) +
  geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.6) +
  stat_cor(method = 'pearson', label.x.npc = 'left', label.y.npc = 'top', size = 3.5) +
  scale_fill_manual(values = hist_palette) +
  labs(x = '% Epithelioid component (Bueno et al.)',
       y = 'Regulatory axis score\n(SOX9/RUNX1/RUNX2/SNAI2)',
       fill = 'Histology',
       title = 'Axis score vs % epithelioid') +
  gtheme

pdf(file.path(outdir, 'Plots', 'F9b_bueno_axis_vs_pct.pdf'), width = 9, height = 4)
print(wrap_plots(f9b_sarc, f9b_epit, ncol = 2))
dev.off()

##############################################################################
# F9b2 – Individual gene expression vs % sarcomatoid OR % epithelioid (Bueno)
##############################################################################
make_pct_scatter <- function(gene, df, expr_mat, pct_col, xlab) {
  if (!gene %in% rownames(expr_mat)) return(NULL)
  df$expr    <- as.numeric(expr_mat[gene, rownames(df)])
  df$pct_var <- df[[pct_col]]
  ggplot(df, aes(x = pct_var, y = expr)) +
    geom_point(aes(fill = path_histology), shape = 21, size = 2, alpha = 0.8, stroke = 0.4) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.6) +
    stat_cor(method = 'pearson', label.x.npc = 'left', label.y.npc = 'top', size = 3) +
    scale_fill_manual(values = hist_palette) +
    labs(x = xlab, y = paste0(gene, ' expression (log2)'),
         fill = 'Histology', title = gene) +
    gtheme + NoLegend()
}

make_composite_scatter <- function(df, pct_col, xlab, title) {
  df$pct_var <- df[[pct_col]]
  ggplot(df, aes(x = pct_var, y = axis_score)) +
    geom_point(aes(fill = path_histology), shape = 21, size = 2, alpha = 0.8, stroke = 0.4) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.6) +
    stat_cor(method = 'pearson', label.x.npc = 'left', label.y.npc = 'top', size = 3) +
    scale_fill_manual(values = hist_palette) +
    labs(x = xlab,
         y = 'Regulatory axis score\n(SOX9/RUNX1/RUNX2/SNAI2)',
         fill = 'Histology', title = title) +
    gtheme + NoLegend()
}

bueno_sarc <- bueno_df[!is.na(bueno_df$sarc_pct), ]
bueno_epit <- bueno_df[!is.na(bueno_df$epit_pct), ]

# % Sarcomatoid panels
gene_sarc_plots <- lapply(axis_tfs, make_pct_scatter,
  df = bueno_sarc, expr_mat = bueno_expr,
  pct_col = 'sarc_pct', xlab = '% Sarcomatoid component')
gene_sarc_plots <- Filter(Negate(is.null), gene_sarc_plots)
gene_sarc_plots[['Composite']] <- make_composite_scatter(
  bueno_sarc, 'sarc_pct', '% Sarcomatoid component', 'Composite axis score')

# % Epithelioid panels
gene_epit_plots <- lapply(axis_tfs, make_pct_scatter,
  df = bueno_epit, expr_mat = bueno_expr,
  pct_col = 'epit_pct', xlab = '% Epithelioid component')
gene_epit_plots <- Filter(Negate(is.null), gene_epit_plots)
gene_epit_plots[['Composite']] <- make_composite_scatter(
  bueno_epit, 'epit_pct', '% Epithelioid component', 'Composite axis score')

ncols <- length(gene_sarc_plots)

pdf(file.path(outdir, 'Plots', 'F9b2_bueno_genes_vs_sarc_pct.pdf'),
    width = 4 * ncols, height = 4)
print(wrap_plots(gene_sarc_plots, ncol = ncols))
dev.off()

pdf(file.path(outdir, 'Plots', 'F9b3_bueno_genes_vs_epit_pct.pdf'),
    width = 4 * ncols, height = 4)
print(wrap_plots(gene_epit_plots, ncol = ncols))
dev.off()

##############################################################################
# F9c – Molecular subtype vs pathological histology concordance (Bueno)
# Shows how the 4-class molecular phenotype maps onto 3-class pathology
##############################################################################
cross_df <- as.data.frame(table(
  Molecular  = bueno_df$mol_subtype,
  Pathological = bueno_df$path_histology
))
cross_df$Molecular <- factor(cross_df$Molecular,
  levels = c('Epithelioid', 'Biphasic-E', 'Biphasic-S', 'Sarcomatoid'))

f9c <- ggplot(cross_df, aes(x = Molecular, y = Freq, fill = Pathological)) +
  geom_bar(stat = 'identity', position = 'fill', alpha = 0.85) +
  scale_fill_manual(values = hist_palette) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = 'Molecular subtype (ConsensusCluster)',
       y = 'Proportion',
       fill = 'Pathological histology',
       title = 'Molecular vs pathological subtype (Bueno et al.)') +
  gtheme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf(file.path(outdir, 'Plots', 'F9c_bueno_mol_vs_path_histology.pdf'), width = 5, height = 4)
print(f9c)
dev.off()

##############################################################################
# F9d – Pairwise co-expression of SOX9 / RUNX2 / SNAI2 across 3 bulk cohorts
##############################################################################
coexp_genes  <- c('SOX9', 'RUNX2', 'SNAI2')
coexp_pairs  <- combn(coexp_genes, 2, simplify = FALSE)   # 3 pairs

meso_bulk_l    <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))
meso_bulk_meta <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))
studies        <- c('bueno', 'tcga', 'mesomics')
study_labels   <- c(bueno = 'Bueno et al.', tcga = 'TCGA', mesomics = 'Mesomics')

bulk_palette <- setNames(
  as.character(paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1, 2, 5, 7, 3)]),
  c('Epithelioid', 'Biphasic-E', 'Biphasic-S', 'Biphasic', 'Sarcomatoid')
)

make_coexp_scatter <- function(gx, gy, study) {
  bulk <- meso_bulk_l[[study]]
  meta <- meso_bulk_meta[[study]]
  if (!gx %in% rownames(bulk) || !gy %in% rownames(bulk)) return(NULL)
  df <- data.frame(
    x       = as.numeric(bulk[gx, ]),
    y       = as.numeric(bulk[gy, ]),
    subtype = meta$subtype
  )
  df <- df[!is.na(df$subtype), ]
  ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(fill = subtype), shape = 21, size = 1.5, alpha = 0.7, stroke = 0.3) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.5) +
    stat_cor(method = 'pearson', label.x.npc = 'left', label.y.npc = 'top', size = 3) +
    scale_fill_manual(values = bulk_palette) +
    labs(x = paste0(gx, ' (log2)'), y = paste0(gy, ' (log2)'),
         title = study_labels[study]) +
    gtheme + NoLegend()
}

# One row per pair, one column per cohort (3 x 3 grid)
coexp_plots <- lapply(coexp_pairs, function(pair) {
  lapply(studies, function(study) make_coexp_scatter(pair[1], pair[2], study))
})

for (i in seq_along(coexp_pairs)) {
  pair    <- coexp_pairs[[i]]
  plots_i <- Filter(Negate(is.null), coexp_plots[[i]])
  pdf(file.path(outdir, 'Plots',
        paste0('F9d_coexpr_', pair[1], '_', pair[2], '_bulk.pdf')),
      width = 4 * length(plots_i), height = 4)
  print(wrap_plots(plots_i, ncol = length(plots_i)))
  dev.off()
}

# Also one combined figure: all pairs x all cohorts (3 rows x 3 cols)
all_coexp <- do.call(c, coexp_plots)
all_coexp <- Filter(Negate(is.null), all_coexp)
pdf(file.path(outdir, 'Plots', 'F9d_coexpr_SOX9_RUNX2_SNAI2_all_cohorts.pdf'),
    width = 12, height = 4 * length(coexp_pairs))
print(wrap_plots(all_coexp, ncol = length(studies)))
dev.off()

message('F9 done.')
