# R1_Q3 – F18: Agreement between metacell co-expression and bulk RNA co-expression
#
# For each anchor (SOX9/RUNX2/SNAI2):
#   1. Load median Spearman r per gene from metacell analysis (F17 output)
#   2. Compute Spearman r per gene in each bulk cohort (Bueno, TCGA, Mesomics)
#   3. Correlate the metacell ranking vs the bulk ranking — measures reproducibility
#   4. Highlight top consistently co-expressed genes in both modalities
#
# Outputs:
#   F18a_metacell_bulk_scatter.pdf       – scatter: metacell r vs bulk r per anchor x cohort
#   F18b_top_genes_agreement_heatmap.pdf – top genes in both modalities
#   F18_metacell_bulk_cor_table.csv      – full gene-level table

set.seed(1234)

projdir  <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
bulk_dir <- file.path(projdir, 'bulkRNA_meso')
outdir   <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))
library(patchwork); library(ggpubr); library(tidyr); library(dplyr)

anchor_genes  <- c('SOX9', 'RUNX2', 'SNAI2')
studies       <- c('bueno', 'tcga', 'mesomics')
study_labels  <- c(bueno = 'Bueno et al.', tcga = 'TCGA', mesomics = 'Mesomics')

##############################################################################
# 1. Load metacell median correlations (F17 output)
##############################################################################
# Use full gene table if available, otherwise fall back to top-50
f17_path <- file.path(outdir, 'F17_coexpr_all_genes_table.csv')
if (!file.exists(f17_path))
  f17_path <- file.path(outdir, 'F17_coexpr_top_genes_table.csv')
if (!file.exists(f17_path)) stop('F17 table not found — run F17 first')
message('Loading F17 from: ', basename(f17_path))

f17 <- read.csv(f17_path, stringsAsFactors = FALSE)

# Compute median rho across samples per anchor x gene (all genes in F17, not just top 50)
metacell_median <- f17 %>%
  group_by(anchor, gene) %>%
  summarise(
    mc_rho       = median(rho, na.rm = TRUE),
    mc_n_samples = sum(!is.na(rho)),
    .groups = 'drop'
  )
message('Metacell genes per anchor: ',
        paste(tapply(metacell_median$gene, metacell_median$anchor, length), collapse=', '))

##############################################################################
# 2. Compute bulk RNA Spearman correlations per anchor per cohort
##############################################################################
meso_bulk_l <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))

bulk_cors <- do.call(rbind, lapply(anchor_genes, function(anchor) {
  do.call(rbind, lapply(studies, function(study) {
    bulk <- meso_bulk_l[[study]]
    if (!anchor %in% rownames(bulk)) return(NULL)
    anchor_expr <- as.numeric(bulk[anchor, ])
    # All genes except the anchor itself
    other_genes <- rownames(bulk)[rownames(bulk) != anchor]
    cor_vec <- apply(bulk[other_genes, , drop = FALSE], 1, function(g)
      tryCatch(cor(anchor_expr, as.numeric(g), method = 'spearman', use = 'complete.obs'),
               error = function(e) NA_real_))
    data.frame(
      anchor = anchor,
      gene   = names(cor_vec),
      bulk_rho  = cor_vec,
      cohort = study_labels[study],
      stringsAsFactors = FALSE
    )
  }))
}))

# Average bulk rho across cohorts per anchor x gene
bulk_cors   <- bulk_cors[!is.na(bulk_cors$bulk_rho), ]
bulk_median <- bulk_cors %>%
  group_by(anchor, gene) %>%
  summarise(bulk_rho_median = median(bulk_rho, na.rm = TRUE),
            bulk_n_cohorts  = sum(!is.na(bulk_rho)),
            .groups = 'drop')

##############################################################################
# 3. Join metacell and bulk correlations
##############################################################################
combined <- inner_join(metacell_median, bulk_median, by = c('anchor','gene'))
combined <- combined[!is.na(combined$mc_rho) & !is.na(combined$bulk_rho_median), ]

# Agreement correlation: how well do metacell ranks predict bulk ranks?
safe_cor <- function(x, y, method)
  tryCatch(cor(x, y, method = method, use = 'complete.obs'), error = function(e) NA_real_)

agreement <- combined %>%
  group_by(anchor) %>%
  summarise(
    r_pearson  = safe_cor(mc_rho, bulk_rho_median, 'pearson'),
    r_spearman = safe_cor(mc_rho, bulk_rho_median, 'spearman'),
    n_genes    = n(),
    .groups = 'drop'
  )
message('\nAgreement between metacell and bulk co-expression:')
print(as.data.frame(agreement))

# Save full table
write.csv(combined, file.path(outdir, 'F18_metacell_bulk_cor_table.csv'),
          row.names = FALSE)

##############################################################################
# F18a – Scatter: metacell median r vs bulk median r (one panel per anchor)
##############################################################################
# Label top 15 genes in each quadrant (high in both, or discordant)
label_top <- function(df, n = 15) {
  df %>% mutate(
    quadrant = case_when(
      mc_rho > 0 & bulk_rho_median > 0 ~ 'both_pos',
      mc_rho < 0 & bulk_rho_median < 0 ~ 'both_neg',
      TRUE ~ 'discordant'
    )
  ) %>%
  group_by(quadrant) %>%
  slice_max(order_by = abs(mc_rho) + abs(bulk_rho_median), n = n) %>%
  ungroup()
}

scatter_plots <- lapply(anchor_genes, function(anc) {
  sub      <- combined[combined$anchor == anc, ]
  agr      <- agreement[agreement$anchor == anc, ]
  top_labs <- label_top(sub, n = 12)

  ggplot(sub, aes(x = mc_rho, y = bulk_rho_median)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey70', linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey70', linewidth = 0.3) +
    geom_point(alpha = 0.3, size = 0.8, color = 'grey40') +
    geom_point(data = top_labs, aes(color = quadrant), size = 1.5, alpha = 0.8) +
    ggrepel::geom_text_repel(data = top_labs, aes(label = gene, color = quadrant),
      size = 2.2, max.overlaps = 25, segment.size = 0.2,
      min.segment.length = 0) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.6) +
    annotate('text', x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3,
      label = sprintf('r=%.3f  ρ=%.3f  n=%d',
        agr$r_pearson, agr$r_spearman, agr$n_genes)) +
    scale_color_manual(values = c(both_pos = '#B2182B', both_neg = '#2166AC',
                                   discordant = '#F4A582')) +
    labs(x = 'Metacell median Spearman r (scRNA)',
         y = 'Bulk RNA median Spearman r (3 cohorts)',
         title = paste0(anc, ' co-expression: scRNA vs bulk')) +
    gtheme + NoLegend()
})

pdf(file.path(outdir, 'Plots', 'F18a_metacell_bulk_scatter.pdf'),
    width = 5 * length(anchor_genes), height = 5)
print(wrap_plots(scatter_plots, ncol = length(anchor_genes)))
dev.off()

##############################################################################
# F18b – Heatmap: top consistently co-expressed genes (high in BOTH modalities)
##############################################################################
# Select genes with |mc_rho| > 0.2 AND |bulk_rho| > 0.1 (both modalities agree)
concordant <- combined %>%
  filter(abs(mc_rho) > 0.15 & abs(bulk_rho_median) > 0.1 &
           sign(mc_rho) == sign(bulk_rho_median)) %>%
  group_by(anchor) %>%
  slice_max(order_by = (abs(mc_rho) + abs(bulk_rho_median)) / 2, n = 25) %>%
  arrange(anchor, desc(mc_rho))

# Long format for heatmap
concordant_long <- concordant %>%
  select(anchor, gene, mc_rho, bulk_rho_median) %>%
  pivot_longer(c(mc_rho, bulk_rho_median), names_to = 'modality', values_to = 'rho') %>%
  mutate(modality = recode(modality,
    mc_rho = 'scRNA (metacell)', bulk_rho_median = 'Bulk RNA (3 cohorts)'))

hm_plots <- lapply(anchor_genes, function(anc) {
  sub <- concordant_long[concordant_long$anchor == anc, ]
  if (nrow(sub) == 0) return(NULL)
  gene_order <- concordant$gene[concordant$anchor == anc]
  sub$gene     <- factor(sub$gene, levels = gene_order)
  sub$modality <- factor(sub$modality,
    levels = c('scRNA (metacell)', 'Bulk RNA (3 cohorts)'))

  ggplot(sub, aes(x = modality, y = gene, fill = rho)) +
    geom_tile(color = 'white', linewidth = 0.4) +
    geom_text(aes(label = sprintf('%.2f', rho)), size = 2.5) +
    scale_fill_gradientn(
      colours = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100)),
      limits = c(-1, 1), name = 'Spearman r') +
    labs(x = NULL, y = NULL,
         title = paste0('Concordant co-expression\nwith ', anc)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          axis.text.y = element_text(size = 7, face = 'italic'))
})
concordant_plots <- Filter(Negate(is.null), hm_plots)

pdf(file.path(outdir, 'Plots', 'F18b_top_genes_agreement_heatmap.pdf'),
    width = 4 * length(concordant_plots), height = 9)
print(wrap_plots(concordant_plots, ncol = length(concordant_plots), guides = 'collect'))
dev.off()

# Also per-cohort scatter (bulk cohort breakdown)
scatter_cohort <- lapply(anchor_genes, function(anc) {
  mc_sub   <- metacell_median[metacell_median$anchor == anc, ]
  bulk_sub <- bulk_cors[bulk_cors$anchor == anc, ]
  df       <- inner_join(mc_sub, bulk_sub, by = c('anchor','gene'))
  df       <- df[!is.na(df$mc_rho) & !is.na(df$bulk_rho), ]
  df$cohort <- factor(df$cohort, levels = unname(study_labels))

  agr_cohort <- df %>%
    group_by(cohort) %>%
    summarise(r = safe_cor(mc_rho, bulk_rho, 'spearman'),
              n = n(), .groups = 'drop') %>%
    mutate(label = ifelse(is.na(r), 'ρ=NA', sprintf('ρ=%.2f', r)))

  ggplot(df, aes(x = mc_rho, y = bulk_rho)) +
    geom_point(alpha = 0.2, size = 0.6, aes(color = cohort)) +
    geom_smooth(method = 'lm', se = FALSE, aes(color = cohort), linewidth = 0.6) +
    geom_text(data = agr_cohort, aes(label = label, color = cohort),
      x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 2.8) +
    scale_color_manual(values = c('Bueno et al.'='#1B4F72','TCGA'='#7D3C98',
                                   'Mesomics'='#1E8449')) +
    labs(x = 'Metacell median Spearman r',
         y = 'Bulk Spearman r (per cohort)',
         title = paste0(anc, ': scRNA vs bulk\n(per cohort)'), color = NULL) +
    facet_wrap(~ cohort, ncol = 1) +
    gtheme + NoLegend()
})

pdf(file.path(outdir, 'Plots', 'F18c_metacell_vs_bulk_per_cohort.pdf'),
    width = 4 * length(anchor_genes), height = 10)
print(wrap_plots(scatter_cohort, ncol = length(anchor_genes)))
dev.off()

message('F18 done.')
