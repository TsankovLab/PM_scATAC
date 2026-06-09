# R1_Q3 – F10: Top co-expressed TFs with SOX9 / RUNX2 / SNAI2 in bulk RNA cohorts
#
# For each anchor gene (SOX9, RUNX2, SNAI2), compute Spearman correlation against
# all active TFs identified in the scATAC analysis across 3 bulk cohorts.
# Outputs:
#   F10_coexpr_top_TFs_heatmap.pdf     – heatmap: TF x cohort, faceted by anchor gene
#   F10_coexpr_top_TFs_dotplot.pdf     – dot/lollipop: top correlated TFs per anchor
#   coexpr_top_TFs_table.csv           – full ranked table

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

anchor_genes <- c('SOX9', 'RUNX2', 'SNAI2')
studies      <- c('bueno', 'tcga', 'mesomics')
study_labels <- c(bueno = 'Bueno et al.', tcga = 'TCGA', mesomics = 'Mesomics')

# TF universe: active TFs from scATAC analysis (110 TFs)
tf_universe <- rownames(read.csv(
  file.path(bulk_dir, 'activeTF_sarcomatoid_correlation.csv'), row.names = 1))

meso_bulk_l    <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))
meso_bulk_meta <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))

##############################################################################
# Compute Spearman correlations: anchor gene vs all TFs, per cohort
##############################################################################
cor_results <- list()

for (study in studies) {
  bulk <- meso_bulk_l[[study]]
  tfs_present <- tf_universe[tf_universe %in% rownames(bulk)]

  for (anchor in anchor_genes) {
    if (!anchor %in% rownames(bulk)) next
    anchor_expr <- as.numeric(bulk[anchor, ])
    # exclude anchor itself from TF list
    tfs_query <- tfs_present[tfs_present != anchor]

    cor_vec <- apply(bulk[tfs_query, , drop = FALSE], 1, function(tf_expr)
      cor(anchor_expr, as.numeric(tf_expr), method = 'spearman', use = 'complete.obs')
    )
    cor_results[[paste(anchor, study, sep = '|')]] <- data.frame(
      anchor  = anchor,
      TF      = names(cor_vec),
      rho     = cor_vec,
      cohort  = study_labels[study],
      stringsAsFactors = FALSE
    )
  }
}

cor_df <- do.call(rbind, cor_results)
rownames(cor_df) <- NULL

# Save full table
write.csv(cor_df, file.path(outdir, 'coexpr_top_TFs_table.csv'), row.names = FALSE)

##############################################################################
# Select top N TFs per anchor: consistent across cohorts (median rho)
##############################################################################
top_n <- 15

top_tfs_per_anchor <- cor_df %>%
  group_by(anchor, TF) %>%
  summarise(median_rho = median(rho), .groups = 'drop') %>%
  group_by(anchor) %>%
  slice_max(order_by = abs(median_rho), n = top_n) %>%
  arrange(anchor, desc(median_rho))

# Wide matrix for heatmap: rows = TF, cols = cohort, facet = anchor
hm_df <- cor_df %>%
  filter(paste(anchor, TF) %in%
    paste(top_tfs_per_anchor$anchor, top_tfs_per_anchor$TF))

##############################################################################
# F10a – Heatmap: top TFs per anchor gene (median rho across cohorts)
##############################################################################
hm_plots <- lapply(anchor_genes, function(anc) {
  sub <- top_tfs_per_anchor[top_tfs_per_anchor$anchor == anc, ]
  sub$TF <- factor(sub$TF, levels = rev(sub$TF[order(sub$median_rho)]))

  ggplot(sub, aes(x = 1, y = TF, fill = median_rho)) +
    geom_tile(color = 'white', linewidth = 0.4) +
    scale_fill_gradientn(
      colours = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100)),
      limits  = c(-1, 1), name = 'Median\nSpearman r') +
    labs(x = NULL, y = NULL,
         title = paste0('Co-expressed with ', anc)) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank())
})

pdf(file.path(outdir, 'Plots', 'F10a_coexpr_top_TFs_median_heatmap.pdf'),
    width = 3 * length(anchor_genes), height = 6)
print(wrap_plots(hm_plots, ncol = length(anchor_genes)))
dev.off()

##############################################################################
# F10b – Per-cohort heatmap: TF x cohort, one panel per anchor gene
##############################################################################
cohort_hm_plots <- lapply(anchor_genes, function(anc) {
  top_tfs <- top_tfs_per_anchor$TF[top_tfs_per_anchor$anchor == anc]
  tf_order <- top_tfs_per_anchor$TF[top_tfs_per_anchor$anchor == anc][
    order(top_tfs_per_anchor$median_rho[top_tfs_per_anchor$anchor == anc])]

  sub <- cor_df[cor_df$anchor == anc & cor_df$TF %in% top_tfs, ]
  sub$TF     <- factor(sub$TF, levels = tf_order)
  sub$cohort <- factor(sub$cohort, levels = unname(study_labels))

  ggplot(sub, aes(x = cohort, y = TF, fill = rho)) +
    geom_tile(color = 'white', linewidth = 0.4) +
    scale_fill_gradientn(
      colours = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100)),
      limits  = c(-1, 1), name = 'Spearman r') +
    labs(x = NULL, y = NULL, title = paste0('Co-expressed with ', anc)) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank())
})

pdf(file.path(outdir, 'Plots', 'F10b_coexpr_top_TFs_per_cohort_heatmap.pdf'),
    width = 3.5 * length(anchor_genes), height = 6)
print(wrap_plots(cohort_hm_plots, ncol = length(anchor_genes), guides = 'collect'))
dev.off()

##############################################################################
# F10c – Lollipop plot: top positive and negative TFs per anchor (median rho)
##############################################################################
top_pos_neg <- cor_df %>%
  group_by(anchor, TF) %>%
  summarise(median_rho = median(rho), .groups = 'drop') %>%
  group_by(anchor) %>%
  arrange(desc(median_rho)) %>%
  slice(c(1:10, (n()-9):n())) %>%   # top 10 positive + top 10 negative
  ungroup()

lollipop_plots <- lapply(anchor_genes, function(anc) {
  sub <- top_pos_neg[top_pos_neg$anchor == anc, ]
  sub$TF <- factor(sub$TF, levels = sub$TF[order(sub$median_rho)])
  sub$direction <- ifelse(sub$median_rho > 0, 'Positive', 'Negative')

  ggplot(sub, aes(x = median_rho, y = TF, color = direction)) +
    geom_segment(aes(x = 0, xend = median_rho, y = TF, yend = TF),
      linewidth = 0.6) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey50', linewidth = 0.4) +
    scale_color_manual(values = c(Positive = '#B2182B', Negative = '#2166AC')) +
    labs(x = 'Median Spearman r (3 cohorts)', y = NULL,
         title = paste0('Top co-expressed TFs with ', anc),
         color = NULL) +
    xlim(-1, 1) +
    gtheme + theme(legend.position = 'none')
})

pdf(file.path(outdir, 'Plots', 'F10c_coexpr_top_TFs_lollipop.pdf'),
    width = 4 * length(anchor_genes), height = 5)
print(wrap_plots(lollipop_plots, ncol = length(anchor_genes)))
dev.off()

message('F10 done.')
