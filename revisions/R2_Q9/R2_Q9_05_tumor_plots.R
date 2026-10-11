# R2_Q9 -- tumour-compartment plots coloured by sample (scatac_tumor_SOX9.R), improved palette:
#   Sarcomatoid_signatures_boxplot   violin + boxplot of the scS-score (cNMF20) of malignant cells per tumour
#   SOX6_9_scSscore_scatterplot      SOX6 / SOX9 expression vs mean scS-score per sample (scRNA-seq)
source ('R2_Q9_common.R')
suppressPackageStartupMessages ({ library(Seurat); library(ggpubr) })

# scS violin ####
archp = load_project ('tumor_compartment/scatac_ArchR')
archp_meta = as.data.frame (archp@cellColData)
sams = as.character(unique(archp_meta$Sample3))
tumor_sams = sams[!sams %in% c('normal1','P11_HOX')] # remove normal,low cell numbers and outlier samples
archp_meta = archp_meta[archp_meta$Sample3 %in% tumor_sams,]
average_by_group <- as.data.frame(archp_meta) %>%
  group_by (Sample3) %>%
  summarise(Average = median(cNMF20))
average_by_group = average_by_group[order(-average_by_group$Average),]
sample_sarc_order = factor (archp_meta$Sample3, levels = average_by_group$Sample3)

ccomp_df = as.data.frame (archp@cellColData)
ccomp_df = ccomp_df[ccomp_df$Sample3 %in% tumor_sams,]
ccomp_df$Sample3 = factor (ccomp_df$Sample3, levels = levels (sample_sarc_order))
mk_box = function (pal) ggplot (ccomp_df, aes_string (x= 'Sample3', y= 'cNMF20')) +
  geom_violin (trim=TRUE, aes_string (fill = 'Sample3'),size=2,
    width=1,
    scale='width',
    linewidth = .2, alpha=0.7) +
  geom_boxplot (aes_string(fill = 'Sample3'),
    linewidth = .2,
    width=0.2,
    outlier.alpha = 0.2,
    outlier.size = 1,
    outlier.shape = NA,
     size=0.3, alpha=0.7
     ) + ylim (c(-0.4,.8))+
  gtheme +
  scale_fill_manual (values= pal) +
  NoLegend()
pdf (file.path (outdir, 'tumor_Sarcomatoid_signatures_boxplot.pdf'),width=3,2)
print (mk_box (palette_sample))
print (mk_box (palette_sample_old))
dev.off()

# SOX6/9 sc-Score scatterplot ####
srt = readRDS (file.path (root, 'tumor_compartment', 'scrna', 'srt.rds'))
avg_exp = log2(AverageExpression (srt, group.by = 'sampleID', feature = c('SOX6','SOX9'))[[1]]+1)
avg_exp = t(avg_exp)
# x = mean scS-score (cNMF20 module score) of the malignant cells of each sample. The original recomputes the
# module scores on the fly (ModScoreCor); the saved scRNA object does not carry them, so the per-sample means
# written by the scRNA analysis (cnmf20_sarcomatoid_sample_order.csv) are used.
ccomp_df = read.csv (file.path (root, 'tumor_compartment', 'scrna', 'cnmf20_sarcomatoid_sample_order.csv'), row.names = 1)
ccomp_df = ccomp_df[order (ccomp_df$sampleID), ]
ccomp_df = cbind (ccomp_df, avg_exp[ccomp_df$sampleID,])
# as in the original, the row names are not sample names, so this line does not drop HU37 / HU62 (normal samples)
ccomp_df = ccomp_df[!rownames (ccomp_df) %in% c('HU37','HU62'),]
print (ccomp_df)
mk_sp = function (gene, pal) ggplot(ccomp_df, aes(x = x, y = .data[[gene]])) +
  geom_point(alpha = .8, shape = 21, stroke = 1, aes (fill = sampleID)) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson",
           label.x.npc = "left",  # place in left side of plot
           label.y.npc = "top") + # place near top
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  gtheme
# the original passes rev(palette_sample): with a named palette the order has no effect
pdf (file.path (outdir, 'tumor_SOX6_9_scSscore_scatterplot.pdf'), height=3,width=3.5)
print (mk_sp ('SOX6', palette_sample)); print (mk_sp ('SOX9', palette_sample))
print (mk_sp ('SOX6', palette_sample_old)); print (mk_sp ('SOX9', palette_sample_old))
dev.off()
