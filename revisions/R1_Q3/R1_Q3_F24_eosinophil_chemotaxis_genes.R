# R1_Q3 – F24: Contribution of genes to SOX9 GO EOSINOPHIL CHEMOTAXIS enrichment
#
# Shows:
#   1. Heatmap: per-sample Spearman r of each eosinophil chemotaxis gene with SOX9
#   2. GSEA enrichment plot (running score) on the global median-ranked list
#   3. Scatter plots of top contributing genes vs SOX9 (metacells, all samples)
#   4. Validation in bulk RNA: eosinophil chemotaxis genes vs SOX9

set.seed(1234)

projdir  <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
bulk_dir <- file.path(projdir, 'bulkRNA_meso')
outdir   <- file.path(projdir, 'git_repo_claude', 'R1_Q3')
scrna_dir <- file.path(projdir, 'tumor_compartment', 'scrna')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))
library(fgsea); library(patchwork); library(ggpubr); library(tidyr); library(dplyr)

# Gene set
eos_genes <- c('DAPK2','HRH1','LGALS3','CCL3L3','CCL1','CCL2','CCL3','CCL3L1',
               'CCL5','CCL7','CCL8','CCL11','CCL13','CCL15','CCL18','CCL23',
               'CCL24','XCL1','CX3CL1','XCL2','SCG2')

# Sample order
sarc_order <- read.csv(file.path(scrna_dir, 'cnmf20_sarcomatoid_sample_order.csv'),
                        row.names = 1)
tumor_sams    <- sarc_order$sampleID[!sarc_order$sampleID %in% c('HU37','HU62','normal_pleura')]
sample_levels <- tumor_sams[order(sarc_order$x[sarc_order$sampleID %in% tumor_sams])]

pal_div <- rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100))
pal_samples <- setNames(
  colorRampPalette(c('#2166AC','#FFFFBF','#D6604D'))(length(sample_levels)),
  sample_levels)

bulk_palette <- setNames(
  as.character(paletteer::paletteer_d("rcartocolor::ArmyRose")[c(1,2,5,7,3)]),
  c('Epithelioid','Biphasic-E','Biphasic-S','Sarcomatoid','Biphasic'))

##############################################################################
# Load F17 per-sample SOX9 correlations
##############################################################################
f17 <- read.csv(file.path(outdir, 'F17_coexpr_all_genes_table.csv'),
                stringsAsFactors = FALSE)
sox9_f17 <- f17[f17$anchor == 'SOX9', ]

# Median rho across samples (global ranking)
sox9_global <- sox9_f17 %>%
  group_by(gene) %>%
  summarise(median_rho = median(rho, na.rm = TRUE), .groups = 'drop')
ranked_global <- setNames(sox9_global$median_rho, sox9_global$gene)
ranked_global <- sort(ranked_global[!duplicated(names(ranked_global))],
                      decreasing = TRUE)

##############################################################################
# F24a – Heatmap: per-sample Spearman r of eos genes with SOX9
##############################################################################
# Filter to detected genes only
eos_detected <- eos_genes[eos_genes %in% sox9_f17$gene]
message('Detected eosinophil chemotaxis genes: ', paste(eos_detected, collapse=', '))

hm_df <- sox9_f17[sox9_f17$gene %in% eos_detected, ] %>%
  filter(sample %in% sample_levels) %>%
  mutate(sample = factor(sample, levels = sample_levels))

# Order genes by median rho
gene_order <- hm_df %>%
  group_by(gene) %>%
  summarise(med = median(rho, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(med)) %>% pull(gene)
hm_df$gene <- factor(hm_df$gene, levels = rev(gene_order))

f24a <- ggplot(hm_df, aes(x = sample, y = gene, fill = rho)) +
  geom_tile(color = 'white', linewidth = 0.5) +
  geom_text(aes(label = sprintf('%.2f', rho)), size = 2.8) +
  scale_fill_gradientn(colours = pal_div, limits = c(-1,1),
    oob = scales::squish, name = 'Spearman r') +
  labs(x = 'Sample (low→high scS-score)', y = NULL,
       title = 'SOX9 co-expression: GO Eosinophil Chemotaxis genes',
       subtitle = paste0('Per-sample metacell Spearman r (', length(eos_detected),
                         '/', length(eos_genes), ' genes detected)')) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = 'italic'),
        axis.text.x = element_text(angle = 30, hjust = 1))

pdf(file.path(outdir, 'Plots', 'F24a_eos_genes_SOX9_heatmap.pdf'), width = 7, height = 5)
print(f24a)
dev.off()

##############################################################################
# F24b – GSEA enrichment plot on global ranked list
##############################################################################
eos_in_ranked <- eos_genes[eos_genes %in% names(ranked_global)]
message('Eos genes in global ranked list: ', length(eos_in_ranked))

gsea_res <- fgseaMultilevel(
  pathways   = list(GO_EOSINOPHIL_CHEMOTAXIS = eos_genes),
  stats      = ranked_global,
  minSize    = 5, maxSize    = 500,
  nPermSimple = 10000, eps = 0
)
message(sprintf('GSEA: NES=%.3f  padj=%.4f  leadingEdge: %s',
  gsea_res$NES, gsea_res$padj,
  paste(unlist(gsea_res$leadingEdge), collapse=', ')))

leading_edge <- unlist(gsea_res$leadingEdge)

# Enrichment plot
f24b <- plotEnrichment(list(GO_EOSINOPHIL_CHEMOTAXIS = eos_genes)[[1]], ranked_global) +
  labs(title = sprintf('GSEA: GO EOSINOPHIL CHEMOTAXIS ~ SOX9\nNES=%.2f  padj=%.4f',
                       gsea_res$NES, gsea_res$padj),
       x = 'Gene rank (by Spearman r with SOX9)',
       y = 'Enrichment score') +
  theme_classic(base_size = 10)

pdf(file.path(outdir, 'Plots', 'F24b_eos_GSEA_enrichment_plot.pdf'), width = 6, height = 4)
print(f24b)
dev.off()

##############################################################################
# F24c – Scatter: top leading-edge genes vs SOX9 in metacells (all samples)
##############################################################################
message('Leading edge genes: ', paste(leading_edge, collapse=', '))

mc <- readRDS(file.path(scrna_dir, 'metacells.rds'))
sam_col <- grep('sample|sampleID', colnames(mc@meta.data), ignore.case=TRUE, value=TRUE)[1]
mc <- mc[, mc@meta.data[[sam_col]] %in% tumor_sams]
mc@meta.data$sampleID <- factor(mc@meta.data[[sam_col]], levels = sample_levels)

mc_mat <- tryCatch(
  GetAssayData(mc, assay='RNA', layer='data'),
  error = function(e) tryCatch(
    GetAssayData(mc, assay='RNA', slot='data'),
    error = function(e2) mc@assays$RNA@data))

genes_plot <- intersect(leading_edge, rownames(mc_mat))
if ('SOX9' %in% rownames(mc_mat) && length(genes_plot) > 0) {
  mc_df <- as.data.frame(t(as.matrix(mc_mat[c('SOX9', genes_plot), , drop=FALSE])))
  mc_df$sampleID <- mc@meta.data$sampleID

  scatter_plots <- lapply(genes_plot, function(gene) {
    ggplot(mc_df, aes_string(x='SOX9', y=gene, color='sampleID')) +
      geom_point(size=0.6, alpha=0.5) +
      geom_smooth(aes(group=1), method='lm', se=TRUE,
        color='black', linewidth=0.6) +
      stat_cor(aes(group=1), method='spearman', size=2.8,
        label.x.npc='left', label.y.npc='top', color='black') +
      scale_color_manual(values=pal_samples) +
      labs(x='SOX9', y=gene,
           title=paste0(gene, '\n(leading edge)')) +
      gtheme + NoLegend()
  })

  pdf(file.path(outdir, 'Plots', 'F24c_leading_edge_vs_SOX9_metacells.pdf'),
      width=4*min(length(scatter_plots),4), height=4*ceiling(length(scatter_plots)/4))
  print(wrap_plots(scatter_plots, ncol=min(length(scatter_plots),4)))
  dev.off()
  message('F24c done (', length(scatter_plots), ' genes)')
}

##############################################################################
# F24d – Validate in bulk RNA: eos leading-edge genes vs SOX9
##############################################################################
studies      <- c('bueno','tcga','mesomics')
study_labels <- c(bueno='Bueno et al.', tcga='TCGA', mesomics='Mesomics')
meso_bulk_l  <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))
meso_bulk_meta <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))

bulk_scatter <- lapply(genes_plot, function(gene) {
  lapply(studies, function(study) {
    bulk <- meso_bulk_l[[study]]
    meta <- meso_bulk_meta[[study]]
    if (!gene %in% rownames(bulk) || !'SOX9' %in% rownames(bulk)) return(NULL)
    df <- data.frame(
      SOX9    = as.numeric(bulk['SOX9',]),
      gene_y  = as.numeric(bulk[gene,]),
      subtype = meta$subtype, stringsAsFactors=FALSE)
    df$subtype <- factor(df$subtype,
      levels=intersect(c('Epithelioid','Biphasic-E','Biphasic-S','Biphasic','Sarcomatoid'),
                       unique(df$subtype)))
    r <- as.numeric(cor(df$SOX9, df$gene_y, method='pearson', use='complete.obs'))
    ggplot(df, aes(x=SOX9, y=gene_y)) +
      geom_point(aes(fill=subtype), shape=21, size=1.2, alpha=0.7, stroke=0.3) +
      geom_smooth(method='lm', se=TRUE, color='black', linewidth=0.5) +
      annotate('text', x=-Inf, y=Inf, hjust=-0.1, vjust=1.5, size=2.8,
        label=sprintf('r=%.3f', r)) +
      scale_fill_manual(values=bulk_palette, na.value='grey70') +
      labs(x='SOX9', y=gene,
           title=paste0(gene,' – ', study_labels[study])) +
      gtheme + NoLegend()
  })
})
bulk_scatter_flat <- Filter(Negate(is.null), do.call(c, bulk_scatter))

if (length(bulk_scatter_flat) > 0) {
  pdf(file.path(outdir, 'Plots', 'F24d_leading_edge_vs_SOX9_bulk.pdf'),
      width=4*length(studies), height=4*length(genes_plot))
  print(wrap_plots(bulk_scatter_flat, ncol=length(studies)))
  dev.off()
  message('F24d done')
}

message('F24 done. Leading edge: ', paste(leading_edge, collapse=', '))
