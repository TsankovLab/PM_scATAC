# R1_Q3 – F17: Top genes co-expressed with SOX9 / RUNX2 / SNAI2 per sample
#            Using metacells for robust per-sample correlation
#
# For each anchor gene x sample: Spearman correlation vs all genes in metacell matrix
# Outputs:
#   F17a_coexpr_top_genes_lollipop.pdf   – top 20 pos/neg per anchor (median across samples)
#   F17b_coexpr_heatmap_per_sample.pdf   – top genes x sample heatmap
#   F17c_coexpr_consistent_genes.pdf     – genes consistently co-expressed across ≥5 samples
#   F17_coexpr_top_genes_table.csv       – full ranked table per anchor x sample

set.seed(1234)

projdir   <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
scrna_dir <- file.path(projdir, 'tumor_compartment', 'scrna')
outdir    <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))
library(patchwork); library(tidyr); library(dplyr)

anchor_genes <- c('SOX9', 'RUNX2', 'SNAI2')
top_n        <- 20   # top N positive + top N negative per anchor per sample
min_mc       <- 5    # minimum metacells to attempt correlation

# Sample order
sarc_order    <- read.csv(file.path(scrna_dir, 'cnmf20_sarcomatoid_sample_order.csv'),
                           row.names = 1)
tumor_sams    <- sarc_order$sampleID[
  !sarc_order$sampleID %in% c('HU37','HU62','normal_pleura')]
sample_levels <- tumor_sams[order(sarc_order$x[sarc_order$sampleID %in% tumor_sams])]

pal_samples <- setNames(
  colorRampPalette(c('#2166AC','#FFFFBF','#D6604D'))(length(sample_levels)),
  sample_levels)

##############################################################################
# Load metacells
##############################################################################
message('Loading metacells...')
mc <- readRDS(file.path(scrna_dir, 'metacells.rds'))

# Detect sample column
sam_col <- grep('sample|sampleID', colnames(mc@meta.data),
                ignore.case = TRUE, value = TRUE)[1]
message('Sample column: ', sam_col)
mc <- mc[, mc@meta.data[[sam_col]] %in% tumor_sams]
mc@meta.data$sampleID <- factor(mc@meta.data[[sam_col]], levels = sample_levels)
message('Metacells: ', ncol(mc), ' | Samples: ', nlevels(mc@meta.data$sampleID))

# Extract full expression matrix (Seurat v5 / v4 compatible)
message('Extracting expression matrix...')
expr_mat <- tryCatch(
  as.matrix(GetAssayData(mc, assay = 'RNA', layer = 'data')),
  error = function(e)
    tryCatch(as.matrix(GetAssayData(mc, assay = 'RNA', slot = 'data')),
             error = function(e2) as.matrix(mc@assays$RNA@data))
)
message('Matrix dim: ', nrow(expr_mat), ' genes x ', ncol(expr_mat), ' metacells')

# Check anchor genes
anchors_ok <- anchor_genes[anchor_genes %in% rownames(expr_mat)]
message('Anchor genes found: ', paste(anchors_ok, collapse = ', '))

# Gene filter: keep genes expressed in ≥20% of metacells
expr_rate <- rowMeans(expr_mat > 0)
genes_keep <- rownames(expr_mat)[expr_rate >= 0.2 & !rownames(expr_mat) %in% anchors_ok]
message('Genes passing expression filter: ', length(genes_keep))

##############################################################################
# Compute per-sample Spearman correlation: anchor vs all genes
##############################################################################
message('Computing per-sample correlations...')

cor_results <- list()

for (sam in sample_levels) {
  idx  <- mc@meta.data$sampleID == sam
  n_mc <- sum(idx)
  if (n_mc < min_mc) {
    message('  ', sam, ': skipping (', n_mc, ' metacells)')
    next
  }
  message('  ', sam, ': n_mc=', n_mc)
  mat_s <- expr_mat[, idx, drop = FALSE]   # genes x metacells (this sample)

  for (anchor in anchors_ok) {
    anchor_vec <- as.numeric(mat_s[anchor, ])
    if (var(anchor_vec) == 0) next

    # Correlate anchor against all filtered genes
    cor_vec <- apply(mat_s[genes_keep, , drop = FALSE], 1, function(g) {
      if (var(g) == 0) return(NA)
      cor(anchor_vec, g, method = 'spearman', use = 'complete.obs')
    })

    cor_results[[paste(anchor, sam, sep = '|')]] <- data.frame(
      anchor  = anchor,
      sample  = sam,
      gene    = names(cor_vec),
      rho     = cor_vec,
      n_mc    = n_mc,
      stringsAsFactors = FALSE
    )
    message('  ', sam, ' / ', anchor, ' : n_mc=', n_mc,
            ' top_gene=', names(which.max(cor_vec)))
  }
}

cor_df <- do.call(rbind, cor_results)
cor_df <- cor_df[!is.na(cor_df$rho), ]
cor_df$sample <- factor(cor_df$sample, levels = sample_levels)

# Save ALL genes per anchor x sample (needed for F18 full-gene agreement analysis)
write.csv(cor_df[order(cor_df$anchor, cor_df$sample, -cor_df$rho), ],
          file.path(outdir, 'F17_coexpr_all_genes_table.csv'),
          row.names = FALSE)
# Also save top 50 for supplementary
top_table <- cor_df %>%
  group_by(anchor, sample) %>%
  slice_max(order_by = abs(rho), n = 50) %>%
  arrange(anchor, sample, desc(rho))
write.csv(top_table,
          file.path(outdir, 'F17_coexpr_top_genes_table.csv'),
          row.names = FALSE)

##############################################################################
# Median rho across samples per anchor → select top genes for plotting
##############################################################################
median_rho <- cor_df %>%
  group_by(anchor, gene) %>%
  summarise(
    median_rho  = median(rho, na.rm = TRUE),
    n_samples   = sum(!is.na(rho)),
    .groups = 'drop'
  )

# Top 20 pos + 20 neg per anchor by median rho
top_genes_per_anchor <- median_rho %>%
  group_by(anchor) %>%
  arrange(desc(median_rho)) %>%
  slice(c(1:top_n, (n()-top_n+1):n())) %>%
  ungroup()

##############################################################################
# F17a – Lollipop: top co-expressed genes per anchor (median rho)
##############################################################################
lollipop_plots <- lapply(anchors_ok, function(anc) {
  sub <- top_genes_per_anchor[top_genes_per_anchor$anchor == anc, ]
  sub$gene      <- factor(sub$gene, levels = sub$gene[order(sub$median_rho)])
  sub$direction <- ifelse(sub$median_rho > 0, 'Positive', 'Negative')

  ggplot(sub, aes(x = median_rho, y = gene, color = direction)) +
    geom_segment(aes(x = 0, xend = median_rho, y = gene, yend = gene),
                 linewidth = 0.6) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey60', linewidth = 0.4) +
    scale_color_manual(values = c(Positive = '#B2182B', Negative = '#2166AC')) +
    labs(x = 'Median Spearman r (across samples)', y = NULL,
         title = paste0('Top genes co-expressed\nwith ', anc),
         color = NULL) +
    gtheme + NoLegend()
})

pdf(file.path(outdir, 'Plots', 'F17a_coexpr_top_genes_lollipop.pdf'),
    width = 4.5 * length(anchors_ok), height = 8)
print(wrap_plots(lollipop_plots, ncol = length(anchors_ok)))
dev.off()

##############################################################################
# F17b – Per-sample heatmap: top genes x sample (one panel per anchor)
##############################################################################
hm_plots <- lapply(anchors_ok, function(anc) {
  # Use top 30 by abs median rho
  top_genes <- median_rho %>%
    filter(anchor == anc) %>%
    slice_max(order_by = abs(median_rho), n = 30) %>%
    arrange(median_rho) %>%
    pull(gene)

  hm_df <- cor_df %>%
    filter(anchor == anc, gene %in% top_genes) %>%
    select(gene, sample, rho)
  hm_df$gene <- factor(hm_df$gene, levels = top_genes)

  ggplot(hm_df, aes(x = sample, y = gene, fill = rho)) +
    geom_tile(color = 'white', linewidth = 0.3) +
    scale_fill_gradientn(
      colours = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100)),
      limits = c(-1, 1), na.value = 'grey85', name = 'Spearman r') +
    labs(x = NULL, y = NULL, title = paste0('Co-expressed with ', anc)) +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 1),
          axis.text.y = element_text(size = 7))
})

pdf(file.path(outdir, 'Plots', 'F17b_coexpr_heatmap_per_sample.pdf'),
    width = 5 * length(anchors_ok), height = 9)
print(wrap_plots(hm_plots, ncol = length(anchors_ok), guides = 'collect'))
dev.off()

##############################################################################
# F17c – Genes consistently co-expressed across ≥ half the samples
##############################################################################
n_min_samples <- max(3, floor(length(sample_levels) / 2))

consistent <- median_rho %>%
  filter(n_samples >= n_min_samples) %>%
  group_by(anchor) %>%
  slice_max(order_by = abs(median_rho), n = 25) %>%
  arrange(anchor, desc(median_rho))

consistent_plots <- lapply(anchors_ok, function(anc) {
  sub <- consistent[consistent$anchor == anc, ]
  sub$gene      <- factor(sub$gene, levels = sub$gene[order(sub$median_rho)])
  sub$direction <- ifelse(sub$median_rho > 0, 'Positive', 'Negative')
  ggplot(sub, aes(x = median_rho, y = gene, color = direction,
                  size = n_samples)) +
    geom_segment(aes(x = 0, xend = median_rho, y = gene, yend = gene),
                 linewidth = 0.5) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey60', linewidth = 0.4) +
    scale_color_manual(values = c(Positive = '#B2182B', Negative = '#2166AC')) +
    scale_size_continuous(range = c(1.5, 4), name = 'N samples') +
    labs(x = 'Median Spearman r', y = NULL,
         title = paste0(anc, ': consistent genes\n(≥', n_min_samples, ' samples)'),
         color = NULL) +
    gtheme
})

pdf(file.path(outdir, 'Plots', 'F17c_coexpr_consistent_genes.pdf'),
    width = 5 * length(anchors_ok), height = 7)
print(wrap_plots(consistent_plots, ncol = length(anchors_ok)))
dev.off()

message('F17 done.')
