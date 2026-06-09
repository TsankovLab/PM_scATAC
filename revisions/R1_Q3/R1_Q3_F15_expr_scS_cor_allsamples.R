# R1_Q3 – F15: SOX9 / RUNX2 / SNAI2 expression vs scS-score per cell
# All 9 tumor scRNA samples, no metacells
#
# Outputs:
#   F15a_expr_scS_scatter_global.pdf       – global scatter (all cells, colored by sample)
#   F15b_expr_scS_scatter_per_sample.pdf   – per-sample scatter (3×3 grid per gene)
#   F15c_expr_scS_cor_heatmap.pdf          – Pearson + Spearman r heatmap
#   F15_expr_scS_correlation.csv           – full table

set.seed(1234)

projdir   <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
scrna_dir <- file.path(projdir, 'tumor_compartment', 'scrna')
outdir    <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))
library(patchwork); library(ggpubr); library(tidyr); library(dplyr)

anchor_genes <- c('SOX9', 'RUNX2', 'SNAI2')

# Sample order by scS-score (low → high = epithelioid → sarcomatoid)
sarc_order    <- read.csv(file.path(scrna_dir, 'cnmf20_sarcomatoid_sample_order.csv'),
                           row.names = 1)
tumor_sams    <- sarc_order$sampleID[!sarc_order$sampleID %in% c('HU37','HU62','normal_pleura')]
sample_levels <- tumor_sams[order(sarc_order$x[sarc_order$sampleID %in% tumor_sams])]

pal_samples <- setNames(
  colorRampPalette(c('#2166AC','#FFFFBF','#D6604D'))(length(sample_levels)),
  sample_levels)

##############################################################################
# Load srt and compute scS-score per cell
##############################################################################
message('Loading srt...')
srt <- readRDS(file.path(scrna_dir, 'srt.rds'))
srt <- srt[, srt$sampleID %in% tumor_sams]
srt$sampleID <- factor(srt$sampleID, levels = sample_levels)
message('Cells: ', ncol(srt), ' | Samples: ', paste(levels(srt$sampleID), collapse=', '))

# Compute scS-score (cNMF20 top 20 genes)
cnmf_spectra  <- readRDS(file.path(scrna_dir, 'cnmf_genelist_25_nfeat_5000.rds'))
sarc_genes    <- head(cnmf_spectra[['cNMF20']], 20)
sarc_genes_ok <- sarc_genes[sarc_genes %in% rownames(srt)]
message('cNMF20 genes present: ', length(sarc_genes_ok))

srt <- AddModuleScore(srt, features = list(sarc_genes_ok), name = 'scS_score')
srt$scS_score <- srt$scS_score1
srt$scS_score1 <- NULL

# Check anchor genes
anchor_ok <- anchor_genes[anchor_genes %in% rownames(srt)]
message('Anchor genes present: ', paste(anchor_ok, collapse=', '))

##############################################################################
# Build per-cell data frame
##############################################################################
expr_mat <- as.matrix(GetAssayData(srt, assay='RNA', layer='data')[anchor_ok,,drop=FALSE])
df <- data.frame(
  cell     = colnames(srt),
  sampleID = srt$sampleID,
  scS      = srt$scS_score,
  t(expr_mat),
  check.names = FALSE
)

samples_all <- c('All', levels(srt$sampleID))

##############################################################################
# Compute Pearson + Spearman per gene per sample
##############################################################################
cor_results <- do.call(rbind, lapply(anchor_ok, function(gene) {
  do.call(rbind, lapply(samples_all, function(sam) {
    sub <- if (sam == 'All') df else df[df$sampleID == sam, ]
    n   <- nrow(sub)
    if (n < 10 || var(sub[[gene]], na.rm=TRUE) == 0) {
      return(data.frame(gene=gene, sample=sam, n_cells=n,
                        pearson=NA, spearman=NA))
    }
    data.frame(
      gene     = gene,
      sample   = sam,
      n_cells  = n,
      pearson  = cor(sub[[gene]], sub$scS, method='pearson',  use='complete.obs'),
      spearman = cor(sub[[gene]], sub$scS, method='spearman', use='complete.obs')
    )
  }))
}))

write.csv(cor_results, file.path(outdir, 'F15_expr_scS_correlation.csv'), row.names=FALSE)

# Print table
for (gene in anchor_ok) {
  sub <- cor_results[cor_results$gene == gene, ]
  cat(sprintf('\n%s expression vs scS-score:\n', gene))
  cat(sprintf('%-10s %8s %12s %12s\n', 'Sample','N cells','Pearson r','Spearman r'))
  cat(strrep('-', 46), '\n')
  for (i in seq_len(nrow(sub)))
    cat(sprintf('%-10s %8d %12.3f %12.3f\n',
      sub$sample[i], sub$n_cells[i],
      ifelse(is.na(sub$pearson[i]),  NA, sub$pearson[i]),
      ifelse(is.na(sub$spearman[i]), NA, sub$spearman[i])))
}

##############################################################################
# F15a – Global scatter per gene (all cells, colored by sample)
##############################################################################
f15a_plots <- lapply(anchor_ok, function(gene) {
  ggplot(df, aes_string(x='scS', y=gene, color='sampleID')) +
    geom_point(size=0.3, alpha=0.3) +
    geom_smooth(aes(group=1), method='lm', se=TRUE,
      color='black', linewidth=0.7) +
    stat_cor(aes(group=1), method='pearson', size=3,
      label.x.npc='left', label.y.npc='top', color='black') +
    scale_color_manual(values=pal_samples) +
    labs(x='scS-score (cNMF20)', y=paste0(gene,' expression'),
         title=paste0(gene,' vs scS (all samples)'), color='Sample') +
    gtheme
})

pdf(file.path(outdir,'Plots','F15a_expr_scS_scatter_global.pdf'),
    width=5*length(anchor_ok), height=5)
print(wrap_plots(f15a_plots, ncol=length(anchor_ok)))
dev.off()

##############################################################################
# F15b – Per-sample scatter (one page per gene, grid of samples)
##############################################################################
pdf(file.path(outdir,'Plots','F15b_expr_scS_scatter_per_sample.pdf'),
    width=4*3, height=4*ceiling(length(sample_levels)/3))

for (gene in anchor_ok) {
  p_list <- lapply(sample_levels, function(sam) {
    sub <- df[df$sampleID == sam, ]
    if (nrow(sub) < 10) return(NULL)
    r_p <- cor(sub[[gene]], sub$scS, method='pearson',  use='complete.obs')
    r_s <- cor(sub[[gene]], sub$scS, method='spearman', use='complete.obs')
    ggplot(sub, aes_string(x='scS', y=gene)) +
      geom_point(size=0.4, alpha=0.4, color=pal_samples[sam]) +
      geom_smooth(method='lm', se=TRUE, color='black', linewidth=0.5) +
      labs(x='scS-score', y=paste0(gene,' expr'),
           title=sam,
           subtitle=sprintf('r=%.2f  ρ=%.2f', r_p, r_s)) +
      gtheme
  })
  p_list <- Filter(Negate(is.null), p_list)
  print(wrap_plots(p_list, ncol=3) +
    plot_annotation(title=paste0(gene,' expression vs scS-score per sample')))
}
dev.off()

##############################################################################
# F15c – Heatmap: Pearson and Spearman r (gene × sample)
##############################################################################
cor_long <- pivot_longer(cor_results, c(pearson, spearman),
  names_to='method', values_to='r')
cor_long$sample <- factor(cor_long$sample, levels=samples_all)
cor_long$gene   <- factor(cor_long$gene,   levels=anchor_ok)
cor_long$method <- factor(cor_long$method, levels=c('pearson','spearman'))

f15c <- ggplot(cor_long, aes(x=sample, y=gene, fill=r)) +
  geom_tile(color='white', linewidth=0.4) +
  geom_text(aes(label=ifelse(!is.na(r), sprintf('%.2f',r), '')), size=2.8) +
  scale_fill_gradientn(
    colours=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging",100)),
    limits=c(-1,1), na.value='grey85', name='r') +
  facet_wrap(~method, ncol=1) +
  labs(x=NULL, y=NULL,
       title='SOX9 / RUNX2 / SNAI2 expression vs scS-score') +
  theme_minimal(base_size=10) +
  theme(strip.text=element_text(face='bold'), panel.grid=element_blank(),
        axis.text.x=element_text(angle=35, hjust=1),
        axis.text.y=element_text(face='italic'))

pdf(file.path(outdir,'Plots','F15c_expr_scS_cor_heatmap.pdf'), width=8, height=5)
print(f15c)
dev.off()

message('F15 done.')
