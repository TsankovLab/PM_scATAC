# R1_Q3 – F21: Correlation of SOX9 / RUNX2 / SNAI2 with TGF-β1 signaling in bulk RNA
#
# Analysis:
#   1. Anchor genes vs TGFB1 expression (direct)
#   2. Anchor genes vs HALLMARK_TGF_BETA_SIGNALING score (mean of gene set)
#   3. Anchor genes vs individual TGF-β pathway genes (heatmap per cohort)
#
# All three bulk cohorts; points colored by histological subtype.

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

##############################################################################
# Load bulk RNA data
##############################################################################
meso_bulk_l    <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))
meso_bulk_meta <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))

##############################################################################
# Load Hallmark TGF-beta signaling gene set
##############################################################################
lines    <- readLines(file.path(projdir, 'git_repo', 'files', 'h.all.v7.4.symbols.gmt'))
tgfb_line <- lines[grep('HALLMARK_TGF_BETA_SIGNALING', lines)]
tgfb_genes <- strsplit(tgfb_line, '\t')[[1]][-(1:2)]
message('HALLMARK_TGF_BETA_SIGNALING: ', length(tgfb_genes), ' genes')
message('Genes: ', paste(tgfb_genes, collapse = ', '))

anchor_genes <- c('SOX9', 'RUNX2', 'SNAI2')

##############################################################################
# Build per-cohort data frames
##############################################################################
make_df <- function(study) {
  bulk <- meso_bulk_l[[study]]
  meta <- meso_bulk_meta[[study]]

  # TGF-β gene set score
  tgfb_present <- tgfb_genes[tgfb_genes %in% rownames(bulk)]
  tgfb_score   <- colMeans(bulk[tgfb_present, , drop = FALSE], na.rm = TRUE)

  # Individual key TGF-β genes
  key_tgfb <- c('TGFB1','TGFB2','TGFB3','TGFBR1','TGFBR2','SMAD2','SMAD3',
                 'SMAD4','SMAD7','LTBP1','FN1','COL1A1','COL3A1','ACTA2',
                 'IL1B','IL1A')
  key_present <- key_tgfb[key_tgfb %in% rownames(bulk)]

  df <- data.frame(
    TGFB_score = as.numeric(tgfb_score),
    subtype    = meta$subtype,
    stringsAsFactors = FALSE
  )
  # Add anchor genes
  for (g in anchor_genes[anchor_genes %in% rownames(bulk)])
    df[[g]] <- as.numeric(bulk[g, ])
  # Add individual TGF-β genes
  for (g in key_present) df[[g]] <- as.numeric(bulk[g, ])

  df$subtype <- factor(df$subtype,
    levels = intersect(c('Epithelioid','Biphasic-E','Biphasic-S','Biphasic','Sarcomatoid'),
                       unique(df$subtype)))
  df$study <- study_labels[study]
  df
}

dfs <- setNames(lapply(studies, make_df), studies)
dfs <- Filter(Negate(is.null), dfs)

##############################################################################
# Correlation summary: each anchor vs TGF-β features
##############################################################################
tgfb_features <- c('TGFB_score','TGFB1','TGFB2','TGFB3','TGFBR1','TGFBR2',
                    'SMAD2','SMAD3','SMAD4','SMAD7','FN1','COL1A1','COL3A1','ACTA2')

cor_summary <- do.call(rbind, lapply(names(dfs), function(study) {
  df <- dfs[[study]]
  do.call(rbind, lapply(anchor_genes, function(anc) {
    if (!anc %in% colnames(df)) return(NULL)
    feats <- intersect(tgfb_features, colnames(df))
    do.call(rbind, lapply(feats, function(g) {
      if (sum(!is.na(df[[g]])) < 10) return(NULL)
      data.frame(
        cohort   = study_labels[study],
        anchor   = anc,
        feature  = g,
        pearson  = as.numeric(cor(df[[anc]], df[[g]], method = 'pearson',  use = 'complete.obs')),
        spearman = as.numeric(cor(df[[anc]], df[[g]], method = 'spearman', use = 'complete.obs')),
        n        = sum(!is.na(df[[g]]) & !is.na(df[[anc]]))
      )
    }))
  }))
}))

write.csv(cor_summary, file.path(outdir, 'F21_SOX9_TGFB_correlation.csv'), row.names = FALSE)

cat('\nAnchor genes vs TGF-β signaling:\n')
cat(sprintf('%-15s %-8s %-18s %10s %12s\n', 'Cohort', 'Anchor', 'Feature', 'Pearson r', 'Spearman r'))
cat(strrep('-', 65), '\n')
for (i in seq_len(nrow(cor_summary)))
  cat(sprintf('%-15s %-8s %-18s %10.3f %12.3f\n',
    cor_summary$cohort[i], cor_summary$anchor[i], cor_summary$feature[i],
    cor_summary$pearson[i], cor_summary$spearman[i]))

##############################################################################
# F21a – Anchor genes vs TGF-β score scatter (anchor x cohort grid)
##############################################################################
scatter_score <- lapply(anchor_genes, function(anc) {
  lapply(names(dfs), function(study) {
    df  <- dfs[[study]]
    lab <- study_labels[study]
    if (!anc %in% colnames(df)) return(NULL)
    sub <- cor_summary[cor_summary$cohort == lab & cor_summary$anchor == anc &
                       cor_summary$feature == 'TGFB_score', ]
    r <- if (nrow(sub) > 0) sub$pearson[1] else NA
    rho <- if (nrow(sub) > 0) sub$spearman[1] else NA
    ggplot(df, aes_string(x = 'TGFB_score', y = anc)) +
      geom_point(aes(fill = subtype), shape = 21, size = 1.8, alpha = 0.75, stroke = 0.3) +
      geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.6) +
      annotate('text', x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3,
        label = sprintf('r=%.3f  ρ=%.3f', r, rho)) +
      scale_fill_manual(values = bulk_palette, na.value = 'grey70') +
      labs(x = 'TGF-β signaling score', y = paste0(anc, ' (log2)'),
           title = paste0(anc, ' – ', lab)) +
      gtheme + NoLegend()
  })
})
scatter_score_flat <- Filter(Negate(is.null), do.call(c, scatter_score))

pdf(file.path(outdir, 'Plots', 'F21a_anchors_vs_TGFbeta_score.pdf'),
    width = 4 * length(names(dfs)), height = 4 * length(anchor_genes))
print(wrap_plots(scatter_score_flat, ncol = length(names(dfs))))
dev.off()

##############################################################################
# F21b – Anchor genes vs TGFB1 gene expression scatter
##############################################################################
scatter_tgfb1 <- lapply(anchor_genes, function(anc) {
  lapply(names(dfs), function(study) {
    df  <- dfs[[study]]
    lab <- study_labels[study]
    if (!anc %in% colnames(df) || !'TGFB1' %in% colnames(df)) return(NULL)
    sub <- cor_summary[cor_summary$cohort == lab & cor_summary$anchor == anc &
                       cor_summary$feature == 'TGFB1', ]
    r <- if (nrow(sub) > 0) sub$pearson[1] else NA
    rho <- if (nrow(sub) > 0) sub$spearman[1] else NA
    ggplot(df, aes_string(x = 'TGFB1', y = anc)) +
      geom_point(aes(fill = subtype), shape = 21, size = 1.8, alpha = 0.75, stroke = 0.3) +
      geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.6) +
      annotate('text', x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3,
        label = sprintf('r=%.3f  ρ=%.3f', r, rho)) +
      scale_fill_manual(values = bulk_palette, na.value = 'grey70') +
      labs(x = 'TGFB1 (log2)', y = paste0(anc, ' (log2)'),
           title = paste0(anc, ' vs TGFB1 – ', lab)) +
      gtheme + NoLegend()
  })
})
scatter_tgfb1_flat <- Filter(Negate(is.null), do.call(c, scatter_tgfb1))

pdf(file.path(outdir, 'Plots', 'F21b_anchors_vs_TGFB1.pdf'),
    width = 4 * length(names(dfs)), height = 4 * length(anchor_genes))
print(wrap_plots(scatter_tgfb1_flat, ncol = length(names(dfs))))
dev.off()

##############################################################################
# F21c – Heatmap: anchor x TGF-β feature, per cohort (faceted)
##############################################################################
cor_hm <- cor_summary[cor_summary$feature != 'TGFB_score', ]
cor_hm$cohort <- factor(cor_hm$cohort, levels = unname(study_labels))
cor_hm$anchor <- factor(cor_hm$anchor, levels = anchor_genes)

feat_order <- cor_hm %>%
  group_by(feature) %>%
  summarise(mean_r = mean(pearson, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(mean_r)) %>% pull(feature)
cor_hm$feature <- factor(cor_hm$feature, levels = feat_order)

f21c <- ggplot(cor_hm, aes(x = anchor, y = feature, fill = pearson)) +
  geom_tile(color = 'white', linewidth = 0.4) +
  geom_text(aes(label = sprintf('%.2f', pearson)), size = 2.5) +
  scale_fill_gradientn(
    colours = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100)),
    limits = c(-1, 1), name = 'Pearson r') +
  facet_wrap(~ cohort, ncol = length(levels(cor_hm$cohort))) +
  labs(x = NULL, y = NULL,
       title = 'SOX9 / RUNX2 / SNAI2 correlation with TGF-beta pathway genes') +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = 'bold'), panel.grid = element_blank(),
        axis.text.y = element_text(face = 'italic', size = 7),
        axis.text.x = element_text(face = 'italic'))

pdf(file.path(outdir, 'Plots', 'F21c_anchors_TGFB_genes_heatmap.pdf'),
    width = 4 * length(levels(cor_hm$cohort)), height = 7)
print(f21c)
dev.off()

##############################################################################
# F21d – SOX9 / RUNX2 / SNAI2 vs IL1B scatter (3 anchors x 3 cohorts)
##############################################################################
scatter_il1b <- lapply(anchor_genes, function(anc) {
  lapply(names(dfs), function(study) {
    df  <- dfs[[study]]
    lab <- study_labels[study]
    if (!anc %in% colnames(df) || !'IL1B' %in% colnames(df)) return(NULL)
    sub <- cor_summary[cor_summary$cohort == lab & cor_summary$anchor == anc &
                       cor_summary$feature == 'IL1B', ]
    r   <- if (nrow(sub) > 0) sub$pearson[1]  else NA
    rho <- if (nrow(sub) > 0) sub$spearman[1] else NA
    ggplot(df, aes_string(x = 'IL1B', y = anc)) +
      geom_point(aes(fill = subtype), shape = 21, size = 1.8, alpha = 0.75, stroke = 0.3) +
      geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.6) +
      annotate('text', x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3,
        label = sprintf('r=%.3f  rho=%.3f', r, rho)) +
      scale_fill_manual(values = bulk_palette, na.value = 'grey70') +
      labs(x = 'IL1B (log2)', y = paste0(anc, ' (log2)'),
           title = paste0(anc, ' vs IL1B - ', lab)) +
      gtheme + NoLegend()
  })
})
scatter_il1b_flat <- Filter(Negate(is.null), do.call(c, scatter_il1b))

pdf(file.path(outdir, 'Plots', 'F21d_anchors_vs_IL1B.pdf'),
    width = 4 * length(names(dfs)), height = 4 * length(anchor_genes))
print(wrap_plots(scatter_il1b_flat, ncol = length(names(dfs))))
dev.off()
message('F21 done.')
