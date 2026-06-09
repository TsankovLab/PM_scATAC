# R1_Q3 – F20: SOX9 expression vs scS-score in bulk RNA, corrected for inflammation
#
# Approach: partial correlation — regress inflammation score out of BOTH
# SOX9 expression and the bulk scS-score, then correlate residuals.
# Compared side-by-side with the uncorrected correlation.
#
# Inflammation score: tumor-intrinsic NF-kB/inflammatory signature
#   IL1A, IL1B, CDKN1A, KDM6B, NFKBIA — expressed by tumor cells (not immune cells)
#   consistent with HALLMARK_INFLAMMATORY_RESPONSE enrichment in SOX9 co-expression
#
# Also shown:
#   - SOX9 residuals by histological subtype (does subtype signal survive correction?)
#   - Scatter: corrected SOX9 vs corrected scS-score (colored by histology)
#
# Outputs:
#   F20a_SOX9_scS_raw_vs_corrected_scatter.pdf
#   F20b_SOX9_corrected_by_histology.pdf
#   F20c_inflammation_score_diagnostics.pdf
#   F20_partial_correlation_table.csv

set.seed(1234)

projdir  <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
bulk_dir <- file.path(projdir, 'bulkRNA_meso')
outdir   <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))
library(patchwork); library(ggpubr); library(tidyr); library(dplyr)

studies      <- c('bueno', 'tcga', 'mesomics')
study_labels <- c(bueno = 'Bueno et al.', tcga = 'TCGA', mesomics = 'Mesomics')

bulk_palette <- setNames(
  as.character(paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1, 2, 5, 7, 3)]),
  c('Epithelioid', 'Biphasic-E', 'Biphasic-S', 'Sarcomatoid', 'Biphasic'))

# Tumor-intrinsic inflammation gene signature (NF-kB / inflammatory response)
# These genes are expressed BY TUMOR CELLS, not immune infiltrates
inflam_genes <- c('IL1A', 'IL1B', 'CDKN1A', 'KDM6B', 'NFKBIA')

##############################################################################
# Load data
##############################################################################
meso_bulk_l    <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))
meso_bulk_meta <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))

##############################################################################
# Compute inflammation score per cohort
##############################################################################
get_inflam_score <- function(study) {
  bulk  <- meso_bulk_l[[study]]
  genes <- inflam_genes[inflam_genes %in% rownames(bulk)]
  message(study, ': tumor-intrinsic inflammation score from genes: ',
          paste(genes, collapse = ', '))
  if (length(genes) == 0) return(setNames(rep(0, ncol(bulk)), colnames(bulk)))
  # Mean log2 expression of the signature genes (already log2 in bulk matrix)
  score <- colMeans(bulk[genes, , drop = FALSE], na.rm = TRUE)
  setNames(as.numeric(score), colnames(bulk))
}

##############################################################################
# Build per-cohort data frames and regress out inflammation
##############################################################################
results <- list()

for (study in studies) {
  bulk <- meso_bulk_l[[study]]
  meta <- meso_bulk_meta[[study]]

  if (!'SOX9' %in% rownames(bulk)) next
  if (!'sarc_score' %in% colnames(meta)) next

  inflam <- get_inflam_score(study)
  sams   <- colnames(bulk)
  n      <- length(sams)

  df <- data.frame(
    sample    = sams,
    SOX9      = as.numeric(bulk['SOX9', sams]),
    scS_score = as.numeric(meta$sarc_score),
    inflam    = as.numeric(inflam[sams]),
    subtype   = meta$subtype,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$scS_score) & !is.na(df$inflam), ]
  df$subtype <- factor(df$subtype,
    levels = intersect(c('Epithelioid','Biphasic-E','Biphasic-S','Biphasic','Sarcomatoid'),
                       unique(df$subtype)))

  # Regress inflammation out of SOX9 and scS-score
  lm_sox9  <- lm(SOX9      ~ inflam, data = df)
  lm_scS   <- lm(scS_score ~ inflam, data = df)
  df$SOX9_corr  <- residuals(lm_sox9)
  df$scS_corr   <- residuals(lm_scS)

  # Correlations: raw and inflammation-corrected
  r_raw    <- as.numeric(cor(df$SOX9,      df$scS_score, method = 'pearson',  use = 'complete.obs'))
  r_corr   <- as.numeric(cor(df$SOX9_corr, df$scS_corr,  method = 'pearson',  use = 'complete.obs'))
  rho_raw  <- as.numeric(cor(df$SOX9,      df$scS_score, method = 'spearman', use = 'complete.obs'))
  rho_corr <- as.numeric(cor(df$SOX9_corr, df$scS_corr,  method = 'spearman', use = 'complete.obs'))

  r2_inflam_sox9 <- summary(lm_sox9)$r.squared
  r2_inflam_scS  <- summary(lm_scS)$r.squared

  message(sprintf('%s  raw r=%.3f rho=%.3f | corrected r=%.3f rho=%.3f | R2_inflam_SOX9=%.3f R2_inflam_scS=%.3f',
    study_labels[study], r_raw, rho_raw, r_corr, rho_corr,
    r2_inflam_sox9, r2_inflam_scS))

  df$study   <- study_labels[study]
  df$r_raw   <- r_raw;   df$r_corr  <- r_corr
  df$rho_raw <- rho_raw; df$rho_corr <- rho_corr
  df$r2_inflam_sox9 <- r2_inflam_sox9
  df$r2_inflam_scS  <- r2_inflam_scS
  results[[study]] <- df
}

all_df <- do.call(rbind, results)
all_df$study <- factor(all_df$study, levels = unname(study_labels))

# Summary table
cor_summary <- all_df %>%
  group_by(study) %>%
  slice(1) %>%
  select(study, r_raw, r_corr, rho_raw, rho_corr, r2_inflam_sox9, r2_inflam_scS) %>%
  ungroup()
write.csv(cor_summary, file.path(outdir, 'F20_partial_correlation_table.csv'),
          row.names = FALSE)

cat('\nCorrelation of SOX9 vs scS-score (raw vs inflammation-corrected):\n')
cat(sprintf('%-15s %10s %10s %10s %10s %14s %14s\n',
    'Cohort', 'Pearson r', 'r (corr)', 'Spearman r', 'rho (corr)',
    'R2 inflam-SOX9', 'R2 inflam-scS'))
cat(strrep('-', 85), '\n')
for (i in seq_len(nrow(cor_summary)))
  cat(sprintf('%-15s %10.3f %10.3f %10.3f %10.3f %14.3f %14.3f\n',
    cor_summary$study[i], cor_summary$r_raw[i], cor_summary$r_corr[i],
    cor_summary$rho_raw[i], cor_summary$rho_corr[i],
    cor_summary$r2_inflam_sox9[i], cor_summary$r2_inflam_scS[i]))

##############################################################################
# F20a – Scatter: raw vs corrected (2 columns per cohort)
##############################################################################
make_scatter <- function(df, x, y, xlabel, ylabel, title, corr_val) {
  ggplot(df, aes_string(x = x, y = y)) +
    geom_point(aes(fill = subtype), shape = 21, size = 1.5, alpha = 0.7, stroke = 0.3) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.6) +
    scale_fill_manual(values = bulk_palette, na.value = 'grey70') +
    annotate('text', x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3,
      label = sprintf('r = %.3f', corr_val)) +
    labs(x = xlabel, y = ylabel, title = title) +
    gtheme + NoLegend()
}

scatter_plots <- list()
for (study in studies) {
  df   <- results[[study]]
  lab  <- study_labels[study]
  scatter_plots[[paste(study,'raw')]] <- make_scatter(
    df, 'scS_score', 'SOX9',
    'scS-score (bulk RNA)', 'SOX9 expression (log2)',
    paste0(lab, '\nRaw'), df$r_raw[1])
  scatter_plots[[paste(study,'corr')]] <- make_scatter(
    df, 'scS_corr', 'SOX9_corr',
    'scS-score (inflam-corrected)', 'SOX9 (inflam-corrected)',
    paste0(lab, '\nInflammation-corrected'), df$r_corr[1])
}

pdf(file.path(outdir, 'Plots', 'F20a_SOX9_scS_raw_vs_corrected_scatter.pdf'),
    width = 4 * length(studies) * 2, height = 4)
print(wrap_plots(scatter_plots, ncol = length(studies) * 2))
dev.off()

##############################################################################
# F20b – SOX9 (raw and corrected) by histological subtype
##############################################################################
subtype_plots <- list()
for (study in studies) {
  df  <- results[[study]]
  lab <- study_labels[study]
  df2 <- df[!is.na(df$subtype), ]
  if (nrow(df2) < 10) next

  make_box <- function(yvar, ytitle, corr_label) {
    stat.test <- tryCatch(df2 %>%
      t_test(reformulate('subtype', yvar)) %>%
      adjust_pvalue(method = 'none') %>% add_significance() %>%
      add_xy_position(x = 'subtype', step.increase = 0.4),
      error = function(e) NULL)
    if (!is.null(stat.test))
      stat.test <- stat.test[
        stat.test$group1 %in% c('Epithelioid','Sarcomatoid') &
        stat.test$group2 %in% c('Epithelioid','Sarcomatoid'), ]

    p <- ggplot(df2, aes_string(x = 'subtype', y = yvar)) +
      geom_boxplot(aes(fill = subtype),
        outlier.shape = 16, outlier.size = 1, outlier.alpha = 0.3,
        lwd = 0.3, alpha = 0.8, width = 0.6) +
      scale_fill_manual(values = bulk_palette) +
      labs(x = NULL, y = ytitle, title = paste0(lab, '\n', corr_label)) +
      gtheme + NoLegend()
    if (!is.null(stat.test) && nrow(stat.test) > 0)
      p <- p + stat_pvalue_manual(stat.test, label = 'p.adj.signif',
                                   hide.ns = TRUE, bracket.nudge.y = 0)
    p
  }
  subtype_plots[[paste(study,'raw')]]  <- make_box('SOX9',      'SOX9 expression', 'Raw')
  subtype_plots[[paste(study,'corr')]] <- make_box('SOX9_corr', 'SOX9 (inflam-corrected)', 'Corrected')
}

pdf(file.path(outdir, 'Plots', 'F20b_SOX9_corrected_by_histology.pdf'),
    width = 4 * length(studies) * 2, height = 4)
print(wrap_plots(subtype_plots, ncol = length(studies) * 2))
dev.off()

##############################################################################
# F20c – Inflammation score diagnostics (distribution + vs SOX9 + vs scS)
##############################################################################
diag_plots <- list()
for (study in studies) {
  df  <- results[[study]]
  lab <- study_labels[study]
  diag_plots[[paste(study,'dist')]] <- ggplot(df, aes(x = inflam)) +
    geom_histogram(fill = 'steelblue', bins = 30, alpha = 0.8) +
    labs(x = 'Inflammation score', y = 'N samples',
         title = paste0(lab, '\nInflammation score dist.')) +
    gtheme
  diag_plots[[paste(study,'vs_sox9')]] <- ggplot(df, aes(x = inflam, y = SOX9)) +
    geom_point(aes(fill = subtype), shape = 21, size = 1.5, alpha = 0.7, stroke = 0.3) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.5) +
    stat_cor(method = 'pearson', size = 3, label.x.npc = 'left', label.y.npc = 'top') +
    scale_fill_manual(values = bulk_palette, na.value = 'grey70') +
    labs(x = 'Inflammation score', y = 'SOX9 expression',
         title = paste0(lab, '\nSOX9 vs inflammation')) +
    gtheme + NoLegend()
  diag_plots[[paste(study,'vs_scS')]] <- ggplot(df, aes(x = inflam, y = scS_score)) +
    geom_point(aes(fill = subtype), shape = 21, size = 1.5, alpha = 0.7, stroke = 0.3) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.5) +
    stat_cor(method = 'pearson', size = 3, label.x.npc = 'left', label.y.npc = 'top') +
    scale_fill_manual(values = bulk_palette, na.value = 'grey70') +
    labs(x = 'Inflammation score', y = 'scS-score (bulk)',
         title = paste0(lab, '\nscS-score vs inflammation')) +
    gtheme + NoLegend()
}

pdf(file.path(outdir, 'Plots', 'F20c_inflammation_score_diagnostics.pdf'),
    width = 12, height = 4 * length(studies))
print(wrap_plots(diag_plots, ncol = 3))
dev.off()

message('F20 done.')
