# R1_Q3 – F12: Immune-corrected expression of SOX9 / RUNX2 / SNAI2 by histologic subtype
#
# Approach: regress out total immune content from each gene's expression
# (lm residuals), then plot corrected expression by histologic subtype.
# Shows side-by-side: raw vs immune-corrected, to isolate the tumor-intrinsic signal.
#
# Immune score used:
#   Mesomics  – Total.Immune.cells.quanTIseq (formally deconvolved)
#   Bueno/TCGA – mean expression of pan-immune marker genes (PTPRC, CD45)
#                supplemented by CD3D, CD19, CD68 for robustness
#
# Outputs:
#   F12a_raw_vs_corrected_boxplots.pdf      – raw and corrected side-by-side per cohort
#   F12b_immune_corrected_combined.pdf      – corrected only, all 3 genes x 3 cohorts
#   F12_immune_correction_summary.csv       – R² of immune regression per gene/cohort

set.seed(1234)

projdir  <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
bulk_dir <- file.path(projdir, 'bulkRNA_meso')
outdir   <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))

library(tidyr)
library(dplyr)
library(patchwork)
library(rstatix)

anchor_genes <- c('SOX9', 'RUNX2', 'SNAI2')
studies      <- c('bueno', 'tcga', 'mesomics')
study_labels <- c(bueno = 'Bueno et al.', tcga = 'TCGA', mesomics = 'Mesomics')

bulk_palette <- setNames(
  as.character(paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1, 2, 5, 7, 3)]),
  c('Epithelioid', 'Biphasic-E', 'Biphasic-S', 'Sarcomatoid', 'Biphasic')
)

meso_bulk_l    <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))
meso_bulk_meta <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))

# Broad immune gene set for PC1-based immune score (covers T, B, NK, myeloid lineages)
# Derived from published pan-cancer immune gene signatures (Bindea 2013, TIMER)
broad_immune_genes <- c(
  # T cells
  'CD3D','CD3E','CD3G','CD2','CD7','CD28','ITK',
  # CD8 cytotoxic
  'CD8A','CD8B','GZMB','GZMA','GZMK','PRF1','NKG7',
  # CD4 helper
  'CD4','IL7R','TCF7','CCR7',
  # Tregs
  'FOXP3','IL2RA','CTLA4','IKZF2',
  # B cells
  'CD19','MS4A1','CD79A','CD79B','IGHM','IGKC',
  # NK cells
  'KLRD1','KLRB1','NCR1','GNLY','NCAM1',
  # Myeloid / macrophage
  'CD68','CD163','MRC1','CSF1R','CD14','FCGR3A','ITGAM',
  # Dendritic cells
  'ITGAX','CLEC9A','FLT3','THBD',
  # General leukocyte
  'PTPRC','INFLAMMATORY'='CXCR4','LCP2','LCK','ZAP70'
)
broad_immune_genes <- unique(broad_immune_genes[!is.na(broad_immune_genes)])

##############################################################################
# Helper: compute robust PC1-based immune score per cohort
##############################################################################
get_immune_score <- function(study) {
  bulk <- meso_bulk_l[[study]]
  meta <- meso_bulk_meta[[study]]

  if (study == 'mesomics') {
    # For Mesomics: mean z-score of individual quanTIseq immune cell fractions
    # (excludes Total immune to avoid double-counting)
    qt_cols <- grep('quanTIseq', colnames(meta), value = TRUE,  ignore.case = TRUE)
    qt_cols <- qt_cols[!grepl('Total', qt_cols, ignore.case = TRUE)]
    if (length(qt_cols) >= 3) {
      qt_mat <- as.data.frame(lapply(meta[, qt_cols, drop = FALSE], as.numeric))
      qt_mat_z <- as.data.frame(lapply(qt_mat, function(x) {
        s <- sd(x, na.rm = TRUE)
        if (is.na(s) || s == 0) return(rep(0, length(x)))
        z <- (x - mean(x, na.rm = TRUE)) / s
        z[!is.finite(z)] <- 0
        z
      }))
      score <- rowMeans(qt_mat_z, na.rm = TRUE)
      return(setNames(score, colnames(bulk)))
    }
  }

  # For all cohorts (and Mesomics fallback): mean of z-scored broad immune genes
  # (robust alternative to PCA — equivalent to PC1 when genes are co-regulated)
  genes <- broad_immune_genes[broad_immune_genes %in% rownames(bulk)]
  message(study, ': using ', length(genes), ' immune genes')
  mat   <- bulk[genes, , drop = FALSE]   # genes x samples
  # Z-score each gene across samples, replace any NA/Inf with 0
  mat_z <- t(apply(mat, 1, function(x) {
    x <- as.numeric(x)
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    z <- (x - mean(x, na.rm = TRUE)) / s
    z[!is.finite(z)] <- 0
    z
  }))
  # Mean z-score across genes = composite immune score per sample
  score <- colMeans(mat_z, na.rm = TRUE)
  # Orient: higher = more immune (check against PTPRC)
  if ('PTPRC' %in% rownames(bulk)) {
    if (cor(score, as.numeric(bulk['PTPRC', ]), method = 'spearman', use = 'complete.obs') < 0)
      score <- -score
  }
  setNames(score, colnames(bulk))
}

##############################################################################
# Helper: immune-correct a gene, return data.frame with raw + residuals
##############################################################################
correct_immune <- function(gene, study) {
  bulk   <- meso_bulk_l[[study]]
  meta   <- meso_bulk_meta[[study]]
  if (!gene %in% rownames(bulk)) return(NULL)

  immune_score <- get_immune_score(study)

  # Align: both immune_score and bulk are indexed by colnames(bulk)
  sams    <- colnames(bulk)
  expr    <- as.numeric(bulk[gene, sams])
  imm     <- as.numeric(immune_score[sams])
  subtype <- meta$subtype   # already aligned positionally to bulk columns

  df <- data.frame(
    sample   = sams,
    raw      = expr,
    immune   = imm,
    subtype  = subtype,
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$subtype) & !is.na(df$immune), ]

  if (nrow(df) < 10) {
    message('Skipping ', gene, ' in ', study, ': too few samples after filtering')
    return(NULL)
  }

  # Regress out immune score
  fit <- lm(raw ~ immune, data = df)
  df$corrected <- residuals(fit)
  df$r2        <- summary(fit)$r.squared
  df$gene      <- gene
  df$study     <- study
  df$cohort    <- study_labels[study]

  # Harmonize subtype factor
  df$subtype <- factor(df$subtype,
    levels = intersect(c('Epithelioid','Biphasic-E','Biphasic-S','Biphasic','Sarcomatoid'),
                       unique(df$subtype)))
  df
}

##############################################################################
# Compute corrected data for all genes x studies
##############################################################################
corrected_l <- list()
r2_records  <- list()

for (gene in anchor_genes) {
  for (study in studies) {
    df <- correct_immune(gene, study)
    if (is.null(df)) next
    key <- paste(gene, study, sep = '|')
    corrected_l[[key]] <- df
    r2_records[[key]] <- data.frame(gene = gene, cohort = study_labels[study],
                                     R2 = unique(df$r2))
  }
}

r2_df <- do.call(rbind, r2_records)
write.csv(r2_df, file.path(outdir, 'F12_immune_correction_summary.csv'), row.names = FALSE)
message('Immune R² summary:')
print(r2_df)

##############################################################################
# Helper: make one boxplot panel (raw or corrected)
##############################################################################
make_panel <- function(df, yvar, ytitle, show_stats = TRUE) {
  p <- ggplot(df, aes_string(x = 'subtype', y = yvar)) +
    geom_boxplot(aes(fill = subtype),
      outlier.shape = 16, outlier.size = 1, outlier.alpha = 0.3,
      lwd = 0.3, alpha = 0.8, width = 0.6) +
    scale_fill_manual(values = bulk_palette) +
    labs(x = NULL, y = ytitle) +
    gtheme + NoLegend()

  if (show_stats) {
    stat.test <- tryCatch(df %>%
      t_test(reformulate('subtype', yvar)) %>%
      adjust_pvalue(method = 'none') %>% add_significance() %>%
      add_xy_position(x = 'subtype', step.increase = 0.45),
      error = function(e) NULL)
    if (!is.null(stat.test)) {
      stat.test <- stat.test[
        stat.test$group1 %in% c('Epithelioid','Sarcomatoid') &
        stat.test$group2 %in% c('Epithelioid','Sarcomatoid'), ]
      if (nrow(stat.test) > 0)
        p <- p + stat_pvalue_manual(stat.test, label = 'p.adj.signif',
          hide.ns = TRUE, bracket.nudge.y = 0)
    }
  }
  p
}

##############################################################################
# F12a – Raw vs immune-corrected side-by-side, per gene per cohort
##############################################################################
pdf(file.path(outdir, 'Plots', 'F12a_raw_vs_corrected_boxplots.pdf'),
    width = 14, height = 4 * length(anchor_genes))

for (gene in anchor_genes) {
  panels <- list()
  for (study in studies) {
    key <- paste(gene, study, sep = '|')
    if (!key %in% names(corrected_l)) next
    df <- corrected_l[[key]]
    r2 <- round(unique(df$r2), 3)

    p_raw  <- make_panel(df, 'raw', paste0(gene, ' expression (log2)')) +
      labs(title = paste0(study_labels[study], ' — raw'))
    p_corr <- make_panel(df, 'corrected',
      paste0(gene, ' (immune-corrected, R²=', r2, ')')) +
      labs(title = paste0(study_labels[study], ' — immune-corrected'))
    panels <- c(panels, list(p_raw, p_corr))
  }
  print(wrap_plots(panels, ncol = length(studies) * 2) +
    plot_annotation(title = gene, theme = theme(plot.title = element_text(face = 'bold'))))
}
dev.off()

##############################################################################
# F12b – Immune-corrected only, all genes x cohorts (compact summary)
##############################################################################
corr_plots <- list()
for (gene in anchor_genes) {
  for (study in studies) {
    key <- paste(gene, study, sep = '|')
    if (!key %in% names(corrected_l)) next
    df <- corrected_l[[key]]
    r2 <- round(unique(df$r2), 3)
    corr_plots[[key]] <- make_panel(df, 'corrected',
      'Immune-corrected expression') +
      labs(title = paste0(gene, ' – ', study_labels[study],
                          '\n(R²=', r2, ')'))
  }
}

pdf(file.path(outdir, 'Plots', 'F12b_immune_corrected_combined.pdf'),
    width = 4 * length(studies), height = 4 * length(anchor_genes))
print(wrap_plots(corr_plots, ncol = length(studies)))
dev.off()

message('F12 done.')
