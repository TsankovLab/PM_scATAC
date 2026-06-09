# R1_Q3 – F16: SOX9 regulon activity vs scS-score (top 50 cNMF20 genes)
#            All 9 PM tumor scRNA samples, per-cell, no metacells
#
# SOX9 regulon activity = mean expression of SCENIC target genes:
#   ACTBL2, COL6A3, CXCR4, FN1, GDF15, GPC6, ITGA5, UPP1
# (SCENIC AUC only covers 3 of our 9 samples, so we use target gene expression
#  as a robust proxy that works across all samples)
#
# Outputs:
#   F16a_SOX9regulon_scS50_scatter_global.pdf
#   F16b_SOX9regulon_scS50_scatter_per_sample.pdf
#   F16c_SOX9regulon_scS50_cor_heatmap.pdf
#   F16_SOX9regulon_scS50_correlation.csv

set.seed(1234)

projdir   <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
scrna_dir <- file.path(projdir, 'tumor_compartment', 'scrna')
outdir    <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))
library(patchwork); library(ggpubr); library(tidyr); library(dplyr)

# SOX9 SCENIC regulon target genes
sox9_targets <- c('ACTBL2', 'COL6A3', 'CXCR4', 'FN1', 'GDF15', 'GPC6', 'ITGA5', 'UPP1')

# Top 50 cNMF20 genes for scS-score
cnmf_spectra <- readRDS(file.path(scrna_dir, 'cnmf_genelist_25_nfeat_5000.rds'))
sarc_genes50 <- head(cnmf_spectra[['cNMF20']], 50)

# Sample order by scS-score
sarc_order    <- read.csv(file.path(scrna_dir, 'cnmf20_sarcomatoid_sample_order.csv'),
                           row.names = 1)
tumor_sams    <- sarc_order$sampleID[!sarc_order$sampleID %in% c('HU37','HU62','normal_pleura')]
sample_levels <- tumor_sams[order(sarc_order$x[sarc_order$sampleID %in% tumor_sams])]

pal_samples <- setNames(
  colorRampPalette(c('#2166AC','#FFFFBF','#D6604D'))(length(sample_levels)),
  sample_levels)

##############################################################################
# Load srt
##############################################################################
message('Loading srt...')
srt <- readRDS(file.path(scrna_dir, 'srt.rds'))
srt <- srt[, srt$sampleID %in% tumor_sams]
srt$sampleID <- factor(srt$sampleID, levels = sample_levels)
message('Cells: ', ncol(srt), ' | Samples: ', nlevels(srt$sampleID))

# Verify genes
sarc_genes_ok  <- sarc_genes50[sarc_genes50 %in% rownames(srt)]
sox9_targs_ok  <- sox9_targets[sox9_targets %in% rownames(srt)]
message('Top-50 cNMF20 genes present: ', length(sarc_genes_ok))
message('SOX9 target genes present:   ', length(sox9_targs_ok), ' -> ',
        paste(sox9_targs_ok, collapse=', '))

##############################################################################
# Compute scores with AddModuleScore
##############################################################################
message('Computing scS-score (top 50 cNMF20)...')
srt <- AddModuleScore(srt, features = list(sarc_genes_ok), name = 'scS50_')
srt$scS50 <- srt$scS50_1; srt$scS50_1 <- NULL

message('Computing SOX9 regulon score (target genes)...')
srt <- AddModuleScore(srt, features = list(sox9_targs_ok), name = 'SOX9reg_')
srt$SOX9_regulon <- srt$SOX9reg_1; srt$SOX9reg_1 <- NULL

##############################################################################
# Build per-cell data frame
##############################################################################
df <- data.frame(
  sampleID     = srt$sampleID,
  scS50        = srt$scS50,
  SOX9_regulon = srt$SOX9_regulon,
  stringsAsFactors = FALSE
)

samples_all <- c('All', levels(srt$sampleID))

##############################################################################
# Compute Pearson + Spearman per sample
##############################################################################
cor_results <- do.call(rbind, lapply(samples_all, function(sam) {
  sub <- if (sam == 'All') df else df[df$sampleID == sam, ]
  n   <- nrow(sub)
  if (n < 10 || var(sub$SOX9_regulon, na.rm=TRUE) == 0) {
    return(data.frame(sample=sam, n_cells=n, pearson=NA, spearman=NA))
  }
  data.frame(
    sample   = sam,
    n_cells  = n,
    pearson  = cor(sub$SOX9_regulon, sub$scS50, method='pearson',  use='complete.obs'),
    spearman = cor(sub$SOX9_regulon, sub$scS50, method='spearman', use='complete.obs')
  )
}))

write.csv(cor_results, file.path(outdir, 'F16_SOX9regulon_scS50_correlation.csv'),
          row.names = FALSE)

# Print
cat('\nSOX9 regulon score (target genes) vs scS-score (top 50 cNMF20):\n')
cat(sprintf('%-12s %8s %12s %12s\n', 'Sample', 'N cells', 'Pearson r', 'Spearman r'))
cat(strrep('-', 48), '\n')
for (i in seq_len(nrow(cor_results)))
  cat(sprintf('%-12s %8d %12.3f %12.3f\n',
    cor_results$sample[i], cor_results$n_cells[i],
    cor_results$pearson[i], cor_results$spearman[i]))

##############################################################################
# F16a – Global scatter (all cells, colored by sample)
##############################################################################
r_g  <- cor_results$pearson[cor_results$sample  == 'All']
rho_g <- cor_results$spearman[cor_results$sample == 'All']

f16a <- ggplot(df, aes(x = scS50, y = SOX9_regulon, color = sampleID)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_smooth(aes(group = 1), method = 'lm', se = TRUE,
    color = 'black', linewidth = 0.7) +
  annotate('text', x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3.5,
    label = sprintf('r = %.3f  ρ = %.3f', r_g, rho_g)) +
  scale_color_manual(values = pal_samples) +
  labs(x = 'scS-score (top 50 cNMF20 genes)',
       y = 'SOX9 regulon score\n(mean expr of SCENIC targets)',
       title = 'SOX9 regulon activity vs scS-score — all 9 samples',
       color = 'Sample') +
  gtheme

pdf(file.path(outdir, 'Plots', 'F16a_SOX9regulon_scS50_scatter_global.pdf'),
    width = 6, height = 5)
print(f16a)
dev.off()

##############################################################################
# F16b – Per-sample scatter
##############################################################################
per_sam_plots <- lapply(sample_levels, function(sam) {
  sub <- df[df$sampleID == sam, ]
  if (nrow(sub) < 10) return(NULL)
  r  <- cor_results$pearson[cor_results$sample  == sam]
  rho <- cor_results$spearman[cor_results$sample == sam]
  ggplot(sub, aes(x = scS50, y = SOX9_regulon)) +
    geom_point(size = 0.4, alpha = 0.4, color = pal_samples[sam]) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.5) +
    labs(x = 'scS-score (top 50)', y = 'SOX9 regulon',
         title = sam,
         subtitle = sprintf('r=%.3f  ρ=%.3f', r, rho)) +
    gtheme
})
per_sam_plots <- Filter(Negate(is.null), per_sam_plots)
n_cols <- 3

pdf(file.path(outdir, 'Plots', 'F16b_SOX9regulon_scS50_scatter_per_sample.pdf'),
    width = 4 * n_cols, height = 4 * ceiling(length(per_sam_plots) / n_cols))
print(wrap_plots(per_sam_plots, ncol = n_cols) +
  plot_annotation(title = 'SOX9 regulon vs scS-score per sample'))
dev.off()

##############################################################################
# F16c – Correlation heatmap: sample × method
##############################################################################
cor_long <- cor_results[cor_results$sample != 'All', ] %>%
  pivot_longer(c(pearson, spearman), names_to = 'method', values_to = 'r')
cor_long$sample <- factor(cor_long$sample, levels = sample_levels)
cor_long$method <- factor(cor_long$method, levels = c('pearson', 'spearman'))

f16c <- ggplot(cor_long, aes(x = sample, y = method, fill = r)) +
  geom_tile(color = 'white', linewidth = 0.5) +
  geom_text(aes(label = ifelse(!is.na(r), sprintf('%.3f', r), 'NA')), size = 3) +
  scale_fill_gradientn(
    colours = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100)),
    limits = c(-1, 1), na.value = 'grey85', name = 'r') +
  labs(x = 'Sample (ordered low→high scS-score)', y = NULL,
       title = 'SOX9 regulon vs scS-score (top 50 cNMF20) — per sample') +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1))

pdf(file.path(outdir, 'Plots', 'F16c_SOX9regulon_scS50_cor_heatmap.pdf'),
    width = 9, height = 3)
print(f16c)
dev.off()

message('F16 done.')
