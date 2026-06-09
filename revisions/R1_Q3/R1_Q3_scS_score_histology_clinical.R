##############################################################################
# R1_Q3 – Clinical significance of the sarcomatoid-like (scS) regulatory axis
#
# Reviewer request: show how the SOX9/RUNX/SNAI2 regulatory program relates
# to histologic subtype at the patient level, varies within tumors, and
# associates with clinical outcomes.
#
# Metadata column names (from bulk_RNA_studies_metadata.rds):
#   subtype  – harmonized histologic category (Epithelioid / Biphasic / Sarcomatoid)
#   census   – survival time (years for Bueno; months for TCGA, Mesomics)
#   status   – death event (1=dead, 0=alive)
#   sarc_score – bulk RNA-derived scS-score (mean expression of top cNMF20 genes)
#   Sarcomatoid.Component – % sarcomatoid component (Mesomics only)
#
# Figures produced:
#   F1  – scS-score (bulk RNA) by histologic subtype (3 cohorts)
#   F2  – SOX9/RUNX1/RUNX2/SNAI2 expression by histologic subtype (3 cohorts)
#   F3  – Composite regulatory axis score by histologic subtype (3 cohorts)
#   F4  – Axis score vs % sarcomatoid component (Mesomics continuous measure)
#   F5  – KM survival by high/low regulatory axis score (Mesomics + Bueno + TCGA)
#   F6  – Forest plot: Cox HRs of SOX9/RUNX1/RUNX2/SNAI2 across 3 cohorts
#   F7  – TCGA bulk-ATAC scS-score by histology (chromatin-level validation)
#   F8  – Within-tumor heterogeneity: per-cell scS-score violin per sample (scRNA)
##############################################################################

set.seed(1234)

# ---- Paths ---------------------------------------------------------------
projdir <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
bulk_dir <- file.path(projdir, 'bulkRNA_meso')
outdir   <- file.path(projdir, 'git_repo_claude', 'R1_Q3')
dir.create(file.path(outdir, 'Plots'), recursive = TRUE, showWarnings = FALSE)

# ---- Packages ------------------------------------------------------------
source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))

library(survival)
library(survminer)
library(ggpubr)
library(tidyr)
library(dplyr)
library(patchwork)
library(rstatix)

# ---- Color palettes ------------------------------------------------------
bulk_palette <- setNames(
  as.character(paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1, 2, 5, 7, 3)]),
  c('Epithelioid', 'Biphasic-E', 'Biphasic-S', 'Sarcomatoid', 'Biphasic')
)
km_palette <- c('High axis score' = '#8B1A1A', 'Low axis score' = '#1A3A5C')

# TFs of the SOX9/RUNX/SNAI2 regulatory axis
axis_tfs <- c('SOX9', 'RUNX1', 'RUNX2', 'SNAI2')

##############################################################################
# 1. Load bulk RNA data (already harmonized: subtype, census, status, sarc_score)
##############################################################################
meso_bulk_l    <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))
meso_bulk_meta <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))
studies        <- c('bueno', 'tcga', 'mesomics')
study_labels   <- c(bueno = 'Bueno et al.', tcga = 'TCGA', mesomics = 'Mesomics')

# ---- Helper: harmonize subtype factor across cohorts --------------------
harmonize_subtype <- function(meta, study) {
  df <- data.frame(subtype = meta$subtype)
  # Collapse Biphasic subtypes when missing the -E/-S distinction
  df$subtype[df$subtype %in% c('Biphasic-E', 'Biphasic-S')] <-
    df$subtype[df$subtype %in% c('Biphasic-E', 'Biphasic-S')]
  df$subtype[is.na(df$subtype)] <- NA
  df$subtype <- factor(df$subtype,
    levels = c('Epithelioid', 'Biphasic-E', 'Biphasic-S', 'Biphasic', 'Sarcomatoid'))
  df
}

##############################################################################
# F1 – scS-score (bulk RNA, previously derived) by histologic subtype
##############################################################################
f1_plots <- lapply(studies, function(study) {
  meta <- meso_bulk_meta[[study]]
  df   <- data.frame(sarc_score = meta$sarc_score, subtype = meta$subtype)
  df   <- df[!is.na(df$subtype) & !is.na(df$sarc_score), ]
  df$subtype <- factor(df$subtype,
    levels = intersect(c('Epithelioid','Biphasic-E','Biphasic-S','Biphasic','Sarcomatoid'),
                       unique(df$subtype)))

  stat.test <- tryCatch(df %>%
    t_test(sarc_score ~ subtype) %>% adjust_pvalue(method = 'none') %>%
    add_significance() %>% add_xy_position(x = 'subtype', step.increase = 0.45),
    error = function(e) NULL)
  if (!is.null(stat.test))
    stat.test <- stat.test[
      stat.test$group1 %in% c('Epithelioid','Sarcomatoid') &
      stat.test$group2 %in% c('Epithelioid','Sarcomatoid'), ]

  p <- ggplot(df, aes(x = subtype, y = sarc_score)) +
    geom_boxplot(aes(fill = subtype),
      outlier.shape = 16, outlier.size = 1, outlier.alpha = 0.3,
      lwd = 0.3, alpha = 0.8, width = 0.6) +
    scale_fill_manual(values = bulk_palette) +
    labs(title = study_labels[study],
         y = 'scS-score (bulk RNA)', x = NULL) +
    gtheme + NoLegend()
  if (!is.null(stat.test) && nrow(stat.test) > 0)
    p <- p + stat_pvalue_manual(stat.test, label = 'p.adj.signif',
      hide.ns = TRUE, bracket.nudge.y = 0)
  p
})

pdf(file.path(outdir, 'Plots', 'F1_scS_score_by_histology.pdf'),
    width = 10, height = 4)
print(wrap_plots(f1_plots, ncol = 3))
dev.off()

##############################################################################
# F2 – SOX9/RUNX1/RUNX2/SNAI2 individual gene expression by histologic subtype
##############################################################################
make_gene_bp <- function(gene, study) {
  bulk <- meso_bulk_l[[study]]
  meta <- meso_bulk_meta[[study]]
  if (!gene %in% rownames(bulk)) return(NULL)
  df <- data.frame(
    expression = as.numeric(bulk[gene, ]),
    subtype    = meta$subtype
  )
  df <- df[!is.na(df$subtype), ]
  df$subtype <- factor(df$subtype,
    levels = intersect(c('Epithelioid','Biphasic-E','Biphasic-S','Biphasic','Sarcomatoid'),
                       unique(df$subtype)))

  stat.test <- tryCatch(df %>%
    t_test(expression ~ subtype) %>% adjust_pvalue(method = 'none') %>%
    add_significance() %>% add_xy_position(x = 'subtype', step.increase = 0.45),
    error = function(e) NULL)
  if (!is.null(stat.test))
    stat.test <- stat.test[
      stat.test$group1 %in% c('Epithelioid','Sarcomatoid') &
      stat.test$group2 %in% c('Epithelioid','Sarcomatoid'), ]

  p <- ggplot(df, aes(x = subtype, y = expression)) +
    geom_boxplot(aes(fill = subtype),
      outlier.shape = 16, outlier.size = 1, outlier.alpha = 0.3,
      lwd = 0.3, alpha = 0.8, width = 0.6) +
    scale_fill_manual(values = bulk_palette) +
    labs(title = paste0(gene, ' – ', study_labels[study]),
         y = 'log2 expression', x = NULL) +
    gtheme + NoLegend()
  if (!is.null(stat.test) && nrow(stat.test) > 0)
    p <- p + stat_pvalue_manual(stat.test, label = 'p.adj.signif',
      hide.ns = TRUE, bracket.nudge.y = 0)
  p
}

f2_plots <- lapply(axis_tfs, function(gene)
  lapply(studies, function(study) make_gene_bp(gene, study)))
f2_flat  <- Filter(Negate(is.null), do.call(c, f2_plots))

pdf(file.path(outdir, 'Plots', 'F2_axis_TFs_by_histology.pdf'),
    width = 12, height = 3.5 * length(axis_tfs))
print(wrap_plots(f2_flat, ncol = length(studies)))
dev.off()

##############################################################################
# F3 – Composite regulatory axis score by histologic subtype (3 cohorts)
##############################################################################
compute_axis_score <- function(study) {
  bulk  <- meso_bulk_l[[study]]
  meta  <- meso_bulk_meta[[study]]
  genes <- axis_tfs[axis_tfs %in% rownames(bulk)]
  if (length(genes) == 0) return(NULL)
  mat   <- scale(t(bulk[genes, , drop = FALSE]))  # samples x genes
  score <- rowMeans(mat, na.rm = TRUE)
  df    <- data.frame(axis_score = score, subtype = meta$subtype,
                      census = meta$census, status = meta$status)
  df$subtype <- factor(df$subtype,
    levels = intersect(c('Epithelioid','Biphasic-E','Biphasic-S','Biphasic','Sarcomatoid'),
                       unique(df$subtype)))
  df$study <- study_labels[study]
  df
}

axis_df_l <- setNames(lapply(studies, compute_axis_score), studies)

f3_plots <- lapply(studies, function(study) {
  df <- axis_df_l[[study]]
  df <- df[!is.na(df$subtype), ]
  stat.test <- tryCatch(df %>%
    t_test(axis_score ~ subtype) %>% adjust_pvalue(method = 'none') %>%
    add_significance() %>% add_xy_position(x = 'subtype', step.increase = 0.45),
    error = function(e) NULL)
  if (!is.null(stat.test))
    stat.test <- stat.test[
      stat.test$group1 %in% c('Epithelioid','Sarcomatoid') &
      stat.test$group2 %in% c('Epithelioid','Sarcomatoid'), ]

  p <- ggplot(df, aes(x = subtype, y = axis_score)) +
    geom_boxplot(aes(fill = subtype),
      outlier.shape = 16, outlier.size = 1, outlier.alpha = 0.3,
      lwd = 0.3, alpha = 0.8, width = 0.6) +
    scale_fill_manual(values = bulk_palette) +
    labs(title = study_labels[study],
         y = 'Regulatory axis score\n(mean z-score SOX9/RUNX1/RUNX2/SNAI2)',
         x = NULL) +
    gtheme + NoLegend()
  if (!is.null(stat.test) && nrow(stat.test) > 0)
    p <- p + stat_pvalue_manual(stat.test, label = 'p.adj.signif',
      hide.ns = TRUE, bracket.nudge.y = 0)
  p
})

pdf(file.path(outdir, 'Plots', 'F3_composite_axis_by_histology.pdf'),
    width = 10, height = 4)
print(wrap_plots(f3_plots, ncol = 3))
dev.off()

##############################################################################
# F4 – Axis score vs % sarcomatoid component (Mesomics continuous measure)
##############################################################################
# Mesomics metadata retains Sarcomatoid.Component in meso_bulk_meta_l
msm_df <- axis_df_l[['mesomics']]
# Attach % sarcomatoid component from stored metadata
msm_meta <- meso_bulk_meta[['mesomics']]
msm_df$sarc_pct  <- msm_meta$Sarcomatoid.Component
msm_df$epit_pct  <- msm_meta$Epithelioid.Component
msm_df <- msm_df[!is.na(msm_df$sarc_pct) & !is.na(msm_df$axis_score), ]

f4 <- ggplot(msm_df, aes(x = sarc_pct, y = axis_score)) +
  geom_point(aes(fill = subtype), shape = 21, size = 2.5, alpha = 0.85, stroke = 0.5) +
  geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.6) +
  stat_cor(method = 'pearson', label.x.npc = 'left', label.y.npc = 'top', size = 3.5) +
  scale_fill_manual(values = bulk_palette, na.value = 'grey70') +
  labs(x = '% Sarcomatoid component (Mesomics)',
       y = 'Regulatory axis score\n(SOX9/RUNX1/RUNX2/SNAI2)',
       fill = 'Histology',
       title = 'Axis score correlates with sarcomatoid component (Mesomics)') +
  gtheme

pdf(file.path(outdir, 'Plots', 'F4_axis_score_vs_sarc_pct_mesomics.pdf'),
    width = 5, height = 4)
print(f4)
dev.off()

##############################################################################
# F5 – Kaplan-Meier survival by high/low regulatory axis score (3 cohorts)
##############################################################################
km_one_cohort <- function(study, time_unit = '') {
  df <- axis_df_l[[study]]
  df <- df[!is.na(df$census) & !is.na(df$status) & !is.na(df$axis_score), ]
  df$status  <- as.integer(df$status)
  df$census  <- as.numeric(as.character(df$census))
  df <- df[!is.na(df$census) & !is.na(df$status), ]
  df$score_group <- ifelse(df$axis_score >= median(df$axis_score),
    'High axis score', 'Low axis score')
  df$score_group <- factor(df$score_group, levels = c('High axis score', 'Low axis score'))

  fit <- survfit(Surv(census, status) ~ score_group, data = df)
  km  <- ggsurvplot(fit, data = df,
    palette      = km_palette,
    pval         = TRUE,
    pval.size    = 3.5,
    risk.table   = TRUE,
    risk.table.fontsize = 3,
    legend.title = '',
    legend.labs  = c('High axis score', 'Low axis score'),
    xlab         = paste0('Time (', time_unit, ')'),
    ylab         = 'Overall survival probability',
    title        = paste0(study_labels[study], ': SOX9/RUNX1/RUNX2/SNAI2 axis'),
    ggtheme      = theme_classic(base_size = 10)
  )
  km
}

km_mesomics <- km_one_cohort('mesomics', 'months')
km_bueno    <- km_one_cohort('bueno', 'years')
km_tcga     <- km_one_cohort('tcga', 'months')

pdf(file.path(outdir, 'Plots', 'F5_KM_axis_score_mesomics.pdf'), width = 6, height = 6)
print(km_mesomics)
dev.off()
pdf(file.path(outdir, 'Plots', 'F5_KM_axis_score_bueno.pdf'), width = 6, height = 6)
print(km_bueno)
dev.off()
pdf(file.path(outdir, 'Plots', 'F5_KM_axis_score_tcga.pdf'), width = 6, height = 6)
print(km_tcga)
dev.off()

##############################################################################
# F6 – Forest plot: Cox HRs for SOX9/RUNX1/RUNX2/SNAI2 across 3 cohorts
##############################################################################
cox_files <- list(
  bueno    = file.path(bulk_dir, 'cox_regression_results_bueno.csv'),
  tcga     = file.path(bulk_dir, 'cox_regression_results_tcga.csv'),
  mesomics = file.path(bulk_dir, 'cox_regression_results_mesomics.csv')
)

cox_df <- do.call(rbind, lapply(names(cox_files), function(s) {
  df        <- read.csv(cox_files[[s]], row.names = 1)
  df$gene   <- rownames(df)
  df$cohort <- study_labels[s]
  df
}))
cox_df <- cox_df[cox_df$gene %in% axis_tfs, ]
cox_df$gene   <- factor(cox_df$gene, levels = rev(axis_tfs))
cox_df$cohort <- factor(cox_df$cohort, levels = unname(study_labels))
cox_df$sig    <- cox_df$P_value_C < 0.05

f6 <- ggplot(cox_df, aes(x = HR, y = gene, color = cohort, shape = sig)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey50', linewidth = 0.4) +
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.2, linewidth = 0.5,
    position = position_dodge(width = 0.6)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
    labels = c('p≥0.05', 'p<0.05'), name = '') +
  scale_color_manual(values = c('Bueno et al.' = '#1B4F72',
                                 'TCGA'         = '#7D3C98',
                                 'Mesomics'     = '#1E8449'), name = 'Cohort') +
  labs(x = 'Hazard Ratio (95% CI)', y = NULL,
       title = 'Cox PH: SOX9/RUNX/SNAI2 regulatory axis') +
  gtheme

pdf(file.path(outdir, 'Plots', 'F6_forest_plot_axis_TFs.pdf'), width = 6, height = 3.5)
print(f6)
dev.off()

##############################################################################
# F7 – TCGA bulk-ATAC scS-score by histology (chromatin-level validation)
##############################################################################
tcga_atac <- read.csv(file.path(bulk_dir, 'TCGA_ATAC_ID_sarc_score_histology.csv'),
                       row.names = 1, stringsAsFactors = FALSE)
tcga_atac <- tcga_atac[!is.na(tcga_atac$histology), ]
tcga_atac_uniq <- tcga_atac[!duplicated(tcga_atac$TCGA_ID), ]
tcga_atac_uniq$histology <- factor(tcga_atac_uniq$histology,
  levels = c('Epithelioid', 'Biphasic', 'Sarcomatoid'))

stat_f7 <- tryCatch(tcga_atac_uniq %>%
  t_test(sarc_score ~ histology) %>% adjust_pvalue(method = 'none') %>%
  add_significance() %>% add_xy_position(x = 'histology', step.increase = 0.45),
  error = function(e) NULL)
if (!is.null(stat_f7))
  stat_f7 <- stat_f7[
    stat_f7$group1 %in% c('Epithelioid','Sarcomatoid') &
    stat_f7$group2 %in% c('Epithelioid','Sarcomatoid'), ]

f7 <- ggplot(tcga_atac_uniq, aes(x = histology, y = sarc_score)) +
  geom_boxplot(aes(fill = histology),
    outlier.shape = 16, outlier.size = 1.5, outlier.alpha = 0.4,
    lwd = 0.4, alpha = 0.8, width = 0.55) +
  scale_fill_manual(values = c(
    Epithelioid = unname(bulk_palette['Epithelioid']),
    Biphasic    = unname(bulk_palette['Biphasic']),
    Sarcomatoid = unname(bulk_palette['Sarcomatoid'])
  )) +
  labs(x = NULL, y = 'Bulk ATAC scS-score (TCGA)',
       title = 'Chromatin-level scS-score by histologic subtype\n(TCGA bulk ATAC)') +
  gtheme + NoLegend()
if (!is.null(stat_f7) && nrow(stat_f7) > 0)
  f7 <- f7 + stat_pvalue_manual(stat_f7, label = 'p.adj.signif',
    hide.ns = TRUE, bracket.nudge.y = 0)

pdf(file.path(outdir, 'Plots', 'F7_TCGA_ATAC_scS_by_histology.pdf'), width = 4, height = 4)
print(f7)
dev.off()

##############################################################################
# F8 – Within-tumor heterogeneity: per-cell scS-score violin per patient (scRNA)
##############################################################################
srt_path <- file.path(projdir, 'tumor_compartment', 'scrna', 'scRNA_meso.rds')
if (file.exists(srt_path)) {
  srt <- readRDS(srt_path)

  sarc_order_df <- read.csv(
    file.path(projdir, 'tumor_compartment', 'scrna', 'cnmf20_sarcomatoid_sample_order.csv'),
    row.names = 1
  )
  tumor_sams    <- sarc_order_df$sampleID[
    !sarc_order_df$sampleID %in% c('HU37', 'HU62', 'normal_pleura')]
  sample_levels <- tumor_sams[order(sarc_order_df$x[sarc_order_df$sampleID %in% tumor_sams])]

  meta <- srt@meta.data[srt@meta.data$sampleID %in% tumor_sams, ]
  meta$sampleID <- factor(meta$sampleID, levels = sample_levels)

  pal_samples <- setNames(
    colorRampPalette(c('#2166AC', '#FFFFBF', '#D6604D'))(length(sample_levels)),
    sample_levels
  )

  f8 <- ggplot(meta, aes(x = sampleID, y = cNMF20, fill = sampleID)) +
    geom_violin(trim = TRUE, scale = 'width', linewidth = 0.2, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.3, alpha = 0.9, fill = 'white') +
    scale_fill_manual(values = pal_samples) +
    labs(x = 'Patient sample (ordered by mean scS-score → sarcomatoid)',
         y = 'scS-score per cell (cNMF20)',
         title = 'Within-tumor heterogeneity of the sarcomatoid-like program') +
    gtheme + NoLegend() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  pdf(file.path(outdir, 'Plots', 'F8_within_tumor_scS_score_violins.pdf'),
      width = 7, height = 4)
  print(f8)
  dev.off()

  # F8b – Correlation between sample-level mean scS-score and SOX9/SNAI2 expression
  if (all(c('SOX9', 'SNAI2') %in% rownames(srt))) {
    gene_mean <- sapply(tumor_sams, function(s) {
      cells <- rownames(meta)[meta$sampleID == s]
      c(SOX9  = mean(srt@assays$RNA@data['SOX9',  cells], na.rm = TRUE),
        SNAI2 = mean(srt@assays$RNA@data['SNAI2', cells], na.rm = TRUE),
        scS   = mean(meta$cNMF20[meta$sampleID == s], na.rm = TRUE))
    })
    gene_mean_df <- as.data.frame(t(gene_mean))
    gene_mean_df$sampleID <- rownames(gene_mean_df)

    f8b_sox9 <- ggplot(gene_mean_df, aes(x = scS, y = SOX9)) +
      geom_point(shape = 21, fill = 'steelblue', size = 3, stroke = 0.5) +
      geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.6) +
      stat_cor(method = 'pearson', size = 3.5) +
      labs(x = 'Mean scS-score (cNMF20)', y = 'Mean SOX9 expression',
           title = 'SOX9 vs scS-score per patient') +
      gtheme

    f8b_snai2 <- ggplot(gene_mean_df, aes(x = scS, y = SNAI2)) +
      geom_point(shape = 21, fill = 'firebrick', size = 3, stroke = 0.5) +
      geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.6) +
      stat_cor(method = 'pearson', size = 3.5) +
      labs(x = 'Mean scS-score (cNMF20)', y = 'Mean SNAI2 expression',
           title = 'SNAI2 vs scS-score per patient') +
      gtheme

    pdf(file.path(outdir, 'Plots', 'F8b_SOX9_SNAI2_vs_scS_per_patient.pdf'),
        width = 8, height = 4)
    print(wrap_plots(f8b_sox9, f8b_snai2, ncol = 2))
    dev.off()
  }
}

##############################################################################
# Summary table – Cox HR for axis TFs (write CSV for supplementary table)
##############################################################################
cox_summary <- cox_df[, c('gene', 'cohort', 'HR', 'CI', 'P_value_C')]
colnames(cox_summary) <- c('Gene', 'Cohort', 'HR', '95% CI', 'p-value (univariate Cox)')
write.csv(cox_summary,
  file.path(outdir, 'axis_TF_cox_summary.csv'),
  row.names = FALSE
)

message('Done. Figures written to: ', file.path(outdir, 'Plots'))
