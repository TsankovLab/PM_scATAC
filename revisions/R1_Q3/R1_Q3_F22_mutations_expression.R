# R1_Q3 – F22: SOX9 / RUNX2 / SNAI2 expression by genetic mutation status
#
# Genetic data available:
#   Mesomics  – IHC.BAP1 (YES=loss, NO=retained, formal IHC)
#   Bueno     – NF2..FISH. (FISH-based NF2 status), NF2.status.2.5M.array
#   All cohorts – BAP1 / NF2 / CDKN2A RNA expression as proxy for loss
#
# Analyses:
#   F22a – Boxplots: anchor expression by BAP1-IHC status (Mesomics)
#   F22b – Boxplots: anchor expression by NF2 FISH status (Bueno)
#   F22c – Scatter: anchor vs BAP1 expression (all 3 cohorts)
#   F22d – Scatter: anchor vs NF2 expression (all 3 cohorts)
#   F22e – Heatmap: anchor correlation with driver gene expression (all cohorts)

set.seed(1234)

projdir  <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
bulk_dir <- file.path(projdir, 'bulkRNA_meso')
outdir   <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))
library(patchwork); library(ggpubr); library(rstatix); library(tidyr); library(dplyr)

anchor_genes <- c('SOX9', 'RUNX2', 'SNAI2')
# Key driver genes in mesothelioma
driver_genes <- c('BAP1', 'NF2', 'CDKN2A', 'CDKN2B', 'TP53', 'SETD2', 'LATS2')

studies      <- c('bueno', 'tcga', 'mesomics')
study_labels <- c(bueno = 'Bueno et al.', tcga = 'TCGA', mesomics = 'Mesomics')

bulk_palette <- setNames(
  as.character(paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1, 2, 5, 7, 3)]),
  c('Epithelioid', 'Biphasic-E', 'Biphasic-S', 'Sarcomatoid', 'Biphasic'))

mut_palette <- c('BAP1 loss' = '#D6604D', 'BAP1 retained' = '#4393C3',
                 'NF2 loss' = '#E08214', 'NF2 intact' = '#542788')

##############################################################################
# Load bulk RNA data
##############################################################################
meso_bulk_l    <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))
meso_bulk_meta <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))
bueno_meta_cli <- readRDS(file.path(bulk_dir, 'bueno_meta_cli.rds'))
mesomics_meta  <- read.csv(file.path(bulk_dir, 'mesomics_SamplesOverview.csv'),
                            stringsAsFactors = FALSE)

##############################################################################
# F22a – Mesomics: anchor expression by BAP1-IHC status
##############################################################################
bulk_msm <- meso_bulk_l[['mesomics']]

# Build Mesomics df with IHC.BAP1 — match on Sample column
msm_meta_std <- meso_bulk_meta[['mesomics']]
# IHC.BAP1: YES = BAP1 loss, NO = BAP1 retained
bap1_ihc <- mesomics_meta$IHC.BAP1[match(colnames(bulk_msm), mesomics_meta$Sample)]
bap1_ihc[bap1_ihc %in% c('NA', '', 'YES in MMS/NO in MME')] <- NA
bap1_ihc[bap1_ihc == 'YES'] <- 'BAP1 loss'
bap1_ihc[bap1_ihc == 'NO']  <- 'BAP1 retained'

msm_df <- data.frame(
  BAP1_status = bap1_ihc,
  subtype     = msm_meta_std$subtype,
  stringsAsFactors = FALSE
)
for (g in anchor_genes[anchor_genes %in% rownames(bulk_msm)])
  msm_df[[g]] <- as.numeric(bulk_msm[g, ])
msm_df <- msm_df[!is.na(msm_df$BAP1_status), ]

f22a_plots <- lapply(anchor_genes, function(anc) {
  if (!anc %in% colnames(msm_df)) return(NULL)
  stat.test <- tryCatch(
    msm_df %>% t_test(reformulate('BAP1_status', anc)) %>%
      add_significance() %>% add_xy_position(x = 'BAP1_status'),
    error = function(e) NULL)
  p <- ggplot(msm_df, aes_string(x = 'BAP1_status', y = anc)) +
    geom_boxplot(aes(fill = BAP1_status),
      outlier.shape = 16, outlier.size = 1, outlier.alpha = 0.4,
      lwd = 0.3, alpha = 0.85, width = 0.55) +
    scale_fill_manual(values = mut_palette) +
    labs(x = NULL, y = paste0(anc, ' (log2)'),
         title = paste0(anc, ' – Mesomics (BAP1 IHC)'),
         subtitle = paste0('n=', sum(msm_df$BAP1_status=='BAP1 loss'),
                           ' loss | ', sum(msm_df$BAP1_status=='BAP1 retained'), ' retained')) +
    gtheme + NoLegend()
  if (!is.null(stat.test))
    p <- p + tryCatch(stat_pvalue_manual(stat.test, label = 'p', hide.ns = FALSE), error = function(e) NULL)
  p
})

pdf(file.path(outdir, 'Plots', 'F22a_anchors_by_BAP1_IHC_mesomics.pdf'),
    width = 4 * length(anchor_genes), height = 4)
print(wrap_plots(Filter(Negate(is.null), f22a_plots), ncol = length(anchor_genes)))
dev.off()
message('F22a done')

##############################################################################
# F22b – Bueno: anchor expression by NF2 FISH status
##############################################################################
bulk_bueno <- meso_bulk_l[['bueno']]

# NF2..FISH. column from bueno_meta_cli
nf2_fish <- bueno_meta_cli$NF2..FISH.[match(colnames(bulk_bueno), rownames(bueno_meta_cli))]
cat('NF2 FISH values (Bueno):', paste(sort(unique(nf2_fish)), collapse=', '), '\n')

bueno_df <- data.frame(
  NF2_fish = as.character(nf2_fish),
  subtype  = meso_bulk_meta[['bueno']]$subtype,
  stringsAsFactors = FALSE
)
for (g in anchor_genes[anchor_genes %in% rownames(bulk_bueno)])
  bueno_df[[g]] <- as.numeric(bulk_bueno[g, ])
# Simplify NF2 status to loss vs intact
bueno_df$NF2_status <- ifelse(
  grepl('loss|delete|absent|mono|null|0', bueno_df$NF2_fish, ignore.case = TRUE),
  'NF2 loss', ifelse(
  grepl('normal|intact|present|di|2|wt', bueno_df$NF2_fish, ignore.case = TRUE),
  'NF2 intact', NA))
bueno_df2 <- bueno_df[!is.na(bueno_df$NF2_status), ]
cat('NF2 status (Bueno): loss=', sum(bueno_df2$NF2_status=='NF2 loss'),
    ' | intact=', sum(bueno_df2$NF2_status=='NF2 intact'), '\n')

if (nrow(bueno_df2) > 10) {
  f22b_plots <- lapply(anchor_genes, function(anc) {
    if (!anc %in% colnames(bueno_df2)) return(NULL)
    stat.test <- tryCatch(
      bueno_df2 %>% t_test(reformulate('NF2_status', anc)) %>%
        add_significance() %>% add_xy_position(x = 'NF2_status'),
      error = function(e) NULL)
    p <- ggplot(bueno_df2, aes_string(x = 'NF2_status', y = anc)) +
      geom_boxplot(aes(fill = NF2_status),
        outlier.shape = 16, outlier.size = 1, outlier.alpha = 0.4,
        lwd = 0.3, alpha = 0.85, width = 0.55) +
      scale_fill_manual(values = mut_palette) +
      labs(x = NULL, y = paste0(anc, ' (log2)'),
           title = paste0(anc, ' – Bueno (NF2 FISH)')) +
      gtheme + NoLegend()
    if (!is.null(stat.test))
      p <- p + tryCatch(stat_pvalue_manual(stat.test, label = 'p', hide.ns = FALSE), error = function(e) NULL)
    p
  })
  pdf(file.path(outdir, 'Plots', 'F22b_anchors_by_NF2_FISH_bueno.pdf'),
      width = 4 * length(anchor_genes), height = 4)
  print(wrap_plots(Filter(Negate(is.null), f22b_plots), ncol = length(anchor_genes)))
  dev.off()
  message('F22b done')
} else {
  message('F22b: insufficient NF2 FISH data after filtering')
}

##############################################################################
# F22c/d/e – All cohorts: anchor expression vs driver gene RNA expression
##############################################################################
make_cor_scatter <- function(anchor, driver, study) {
  bulk <- meso_bulk_l[[study]]
  meta <- meso_bulk_meta[[study]]
  if (!anchor %in% rownames(bulk) || !driver %in% rownames(bulk)) return(NULL)
  df <- data.frame(
    anchor_expr = as.numeric(bulk[anchor, ]),
    driver_expr = as.numeric(bulk[driver, ]),
    subtype     = meta$subtype,
    stringsAsFactors = FALSE
  )
  df$subtype <- factor(df$subtype,
    levels = intersect(c('Epithelioid','Biphasic-E','Biphasic-S','Biphasic','Sarcomatoid'),
                       unique(df$subtype)))
  r   <- as.numeric(cor(df$anchor_expr, df$driver_expr, method='pearson',  use='complete.obs'))
  rho <- as.numeric(cor(df$anchor_expr, df$driver_expr, method='spearman', use='complete.obs'))
  ggplot(df, aes(x = driver_expr, y = anchor_expr)) +
    geom_point(aes(fill = subtype), shape = 21, size = 1.5, alpha = 0.7, stroke = 0.3) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.5) +
    annotate('text', x=-Inf, y=Inf, hjust=-0.1, vjust=1.5, size=2.8,
      label = sprintf('r=%.3f  rho=%.3f', r, rho)) +
    scale_fill_manual(values = bulk_palette, na.value = 'grey70') +
    labs(x = paste0(driver, ' (log2)'), y = paste0(anchor, ' (log2)'),
         title = paste0(anchor, ' vs ', driver, '\n', study_labels[study])) +
    gtheme + NoLegend()
}

# F22c – vs BAP1
for (driver in c('BAP1', 'NF2')) {
  plots <- lapply(anchor_genes, function(anc)
    lapply(studies, function(study) make_cor_scatter(anc, driver, study)))
  plots_flat <- Filter(Negate(is.null), do.call(c, plots))
  n_valid <- length(plots_flat)
  if (n_valid > 0) {
    pdf(file.path(outdir, 'Plots', paste0('F22_anchors_vs_', driver, '_expr.pdf')),
        width = 4 * length(studies), height = 4 * length(anchor_genes))
    print(wrap_plots(plots_flat, ncol = length(studies)))
    dev.off()
    message('F22 ', driver, ' scatter done')
  }
}

# F22e – Correlation heatmap: anchor vs all driver genes
cor_hm_df <- do.call(rbind, lapply(studies, function(study) {
  bulk <- meso_bulk_l[[study]]
  do.call(rbind, lapply(anchor_genes, function(anc) {
    if (!anc %in% rownames(bulk)) return(NULL)
    do.call(rbind, lapply(driver_genes, function(drv) {
      if (!drv %in% rownames(bulk)) return(NULL)
      data.frame(
        cohort  = study_labels[study],
        anchor  = anc,
        driver  = drv,
        pearson = as.numeric(cor(bulk[anc,], bulk[drv,], method='pearson',  use='complete.obs')),
        spearman= as.numeric(cor(bulk[anc,], bulk[drv,], method='spearman', use='complete.obs'))
      )
    }))
  }))
}))

write.csv(cor_hm_df, file.path(outdir, 'F22_anchor_driver_correlation.csv'), row.names=FALSE)

cat('\nAnchor vs driver gene correlations (Pearson):\n')
cat(sprintf('%-12s %-8s %-10s %10s %12s\n', 'Cohort','Anchor','Driver','Pearson r','Spearman r'))
cat(strrep('-', 55),'\n')
for (i in seq_len(nrow(cor_hm_df)))
  cat(sprintf('%-12s %-8s %-10s %10.3f %12.3f\n',
    cor_hm_df$cohort[i], cor_hm_df$anchor[i], cor_hm_df$driver[i],
    cor_hm_df$pearson[i], cor_hm_df$spearman[i]))

# Plot heatmap
cor_hm_df$cohort <- factor(cor_hm_df$cohort, levels = unname(study_labels))
cor_hm_df$anchor <- factor(cor_hm_df$anchor, levels = anchor_genes)

f22e <- ggplot(cor_hm_df, aes(x = anchor, y = driver, fill = pearson)) +
  geom_tile(color = 'white', linewidth = 0.4) +
  geom_text(aes(label = sprintf('%.2f', pearson)), size = 2.8) +
  scale_fill_gradientn(
    colours = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100)),
    limits = c(-1, 1), name = 'Pearson r') +
  facet_wrap(~ cohort, ncol = length(levels(cor_hm_df$cohort))) +
  labs(x = NULL, y = NULL,
       title = 'Anchor genes vs driver gene expression') +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = 'bold'), panel.grid = element_blank(),
        axis.text.y = element_text(face = 'italic'),
        axis.text.x = element_text(face = 'italic'))

pdf(file.path(outdir, 'Plots', 'F22e_anchor_vs_driver_heatmap.pdf'),
    width = 4 * length(levels(cor_hm_df$cohort)), height = 5)
print(f22e)
dev.off()

message('F22 done.')
