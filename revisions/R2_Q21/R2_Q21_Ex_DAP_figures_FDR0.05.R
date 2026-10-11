# R2_Q21 -- figures of git_repo/tnk_analysis/scatac_NKTcells_ArchR.R that depend on the CD8 exhausted DAP cut-off,
# regenerated with FDR < 0.05 & Log2FC > 0 (884 peaks) instead of FDR < 0.01 (418 peaks), so that they use the
# same peak set as the labelled heatmap (tnk_exhaustion_annotation) and the table of R2_Q21.
#
#   Exhausted_peaks_MA_plot                 MA plot of CD8_exhausted vs CD8
#   TFmodule_enrichments                    hypergeometric enrichment of the 3 TF modules in the DAP
#   TF_in_module_enrichments                same per TF of module 1
#   density_coverage_celltype_lv2_CD8_ext   summed coverage profile around the DAP per cell type
#   CD8_ext_peakset_enrichments             coverage heatmap around the DAP per cell type
#
# Code is kept as in the original script; only the cut-off changes. Saved results of the T/NK project are read
# (test result, motif matches, TF modules, coverage cache); the project is not modified.
# The GSEA plot (DAP_exhaustion_enrichment_plot) ranks all peaks and does not depend on the cut-off.

suppressPackageStartupMessages({
  library(ArchR); library(SummarizedExperiment); library(ComplexHeatmap)
  library(ggplot2); library(patchwork); library(paletteer)
})
root    = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
tnk_dir = file.path (root, 'NKT_cells', 'scatac_ArchR')
source (file.path (root, 'git_repo', 'utils', 'ggplot_aestetics.R'))
palette_fragments = paletteer_c("ggthemes::Classic Orange-White-Blue",40)   # as in git_repo/utils/palettes.R
dir.create ('Plots', showWarnings = FALSE)

EXH_FDR = 0.05
metaGroupName = 'celltype_lv2'

DAP_list = readRDS (file.path (tnk_dir, 'DAP_CD8_CD8_ext_pairwise.rds'))
DAP_res = do.call (cbind, (assays(DAP_list)))
colnames (DAP_res) = names(assays(DAP_list))
DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames(DAP_res) = as.character(DAP_res_regions)
DAP_res_sig = DAP_res[DAP_res$FDR < EXH_FDR & DAP_res$Log2FC > 0, ]
message (nrow(DAP_res_sig), ' peaks at FDR < ', EXH_FDR)

# MA plot ####
pdf(file.path('Plots','Exhausted_peaks_MA_plot_FDR0.05.pdf'), width=5,height=5)
pma <- markerPlot(seMarker = DAP_list, name = 'CD8_exhausted', cutOff = "FDR <= 0.05", plotAs = "MA")
print (pma)
dev.off()

# Calculate hypergeometric test of peaks with motifs ####
mm = readRDS (file.path (tnk_dir, 'Annotations', 'Motif-Matches-In-Peaks.rds'))
mmMat = assay(mm)
colnames(mmMat) = gsub ('_.*','',colnames(mmMat))
colnames(mmMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames(mmMat))
n_peaks = nrow (mm)

# TF modules (kmeans of the TF activity correlation, saved by the original script)
mods = read.csv (file.path (tnk_dir, 'regTNK_modules.csv'), row.names = 1, check.names = FALSE)
mods = lapply (mods, function(x) x[!is.na(x) & x != ''])
stopifnot (all (unlist (mods) %in% colnames(mmMat)))

hyper_at = function (cols, sig)
  {
  module_peaks = rowRanges(mm)[which (rowSums (mmMat[, cols, drop=F]) > 0)]
  q = length (queryHits (findOverlaps (module_peaks, GRanges (rownames(sig)))))
  n = n_peaks - nrow (sig)
  m = nrow (sig)
  k = length (module_peaks)
  phyper (q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  }
sig01 = DAP_res[DAP_res$FDR < .01 & DAP_res$Log2FC > 0, ]
hyper_res = sapply (mods, hyper_at, sig = DAP_res_sig)

# Barplot of module enrichment in CD8 ext peaks ####
mod_df = data.frame (module = paste0('module',names(mods)), pval = -log10(hyper_res+1e-9))
bp = ggplot (mod_df, aes (x = module, y = pval, fill = module)) +
geom_bar (stat = 'identity',color='grey22') +
scale_fill_manual (values = c(module1 = 'darkred',module2='grey',module3='grey')) +
gtheme
pdf (file.path ('Plots','TFmodule_enrichments_FDR0.05.pdf'), width=4,height=3)
print (bp)
dev.off()

# Module 1 is enriched in CD8 exhausted peaks. Find which TF in it is mostly enriched ####
ext_module = '1'
hyper_res_TF_df = data.frame (
  TF = mods[[ext_module]],
  hyper = p.adjust (sapply (mods[[ext_module]], hyper_at, sig = DAP_res_sig)),
  hyper_FDR0.01_set = p.adjust (sapply (mods[[ext_module]], hyper_at, sig = sig01)))
hyper_res_TF_df = hyper_res_TF_df[order(hyper_res_TF_df$hyper), ]
write.csv (hyper_res_TF_df, 'TF_in_module_enrichments_FDR0.05.csv', row.names = FALSE)
write.csv (data.frame (module = names(mods), n_TF = lengths(mods), p_FDR0.05_set = hyper_res,
  p_FDR0.01_set = sapply (mods, hyper_at, sig = sig01)), 'TFmodule_enrichments_FDR0.05.csv', row.names = FALSE)

hyper_res_TF_df$pval= -log10(hyper_res_TF_df$hyper)
hyper_res_TF_df$TF = factor (hyper_res_TF_df$TF, levels = rev(hyper_res_TF_df$TF))
top_TF = 13
hyper_res_TF_df$label = ifelse (hyper_res_TF_df$TF %in% head(hyper_res_TF_df$TF, top_TF), as.character(hyper_res_TF_df$TF), '')

bp = ggplot (hyper_res_TF_df[hyper_res_TF_df$pval > -log10(0.05),], aes (x = TF, y = pval, fill = pval)) +
geom_bar (stat = 'identity') +
scale_fill_gradient(low = "white", high = "darkred") +
geom_text(aes(label = label, y = pval + .1), hjust = -0.05,na.rm = TRUE, angle=0,size=3) +
gtheme_no_rot +
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
coord_flip()
pdf (file.path ('Plots','TF_in_module_enrichments_FDR0.05.pdf'), width=10,height=3)
print (bp)
dev.off()

# Coverage around the DAP per cell type ####
# window coverage at FDR < 0.05 was computed from the fragments by tnk_exhaustion_annotation/annotate_CD8ext_heatmap.R
# with the code of the original script (peak +- 3 kb, 50 bp windows, step 25, normalised by reads in TSS)
ext_l = readRDS (file.path (tnk_dir, paste0('ext_peaks_',metaGroupName,'_FDR',EXH_FDR,'.rds')))
stopifnot (ncol (ext_l[[1]]) == nrow (DAP_res_sig))

ext_den = lapply (ext_l, function(x) rowSums(x))
ext_den = as.data.frame (do.call (cbind, ext_den))
ext_den$bin = seq_len (nrow (ext_den))
den = lapply (c('FGFBP2_NK','CD4','CD8','KLRC1_NK','CD8_exhausted','Tregs'), function(x) ggplot (ext_den[,c(x,'bin')], aes_string (x= 'bin', y = x)) +
geom_density(color="navyblue", fill="navyblue", alpha=0.5, stat='identity') +
ylim (c(0, max(ext_den[, names(ext_l)]))) +
gtheme_no_rot)
pdf (file.path ('Plots',paste0('density_coverage_',metaGroupName,'_CD8_ext_FDR0.05.pdf')), height=4,width=12)
print (wrap_plots (den, ncol= length(unique(names (ext_l)))))
dev.off()

median_celltype = lapply (ext_l, function(x) colMeans (x))
median_celltype_order = sapply (median_celltype, mean)
print (round (sort (median_celltype_order), 3))
median_celltype_order = names(median_celltype_order)[order(median_celltype_order)]
ext_l2 = lapply (seq_along(ext_l), function(x) ext_l[[x]][,order(-median_celltype[[x]])])
ext_df = do.call (rbind, ext_l2)
ext_df = as.data.frame (t(log10(ext_df+1)))
colnames (ext_df) = rep (names(ext_l), each=dim(ext_l[[1]])[1])
ha = HeatmapAnnotation (bar1 = anno_barplot(colSums(ext_df),gp = gpar(fill = "azure4",border =NA,lty='blank'),border =FALSE, baseline=200,lty='blank'))

pdf (file.path ('Plots','CD8_ext_peakset_enrichments_FDR0.05.pdf'), height=4,width=6)
draw (Heatmap (ext_df,
  top_annotation = ha,
  column_split=factor(colnames(ext_df), levels=median_celltype_order),
  cluster_rows=F,
  column_title_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 0),
  row_names_gp = gpar(fontsize = 0),
  col = rev(palette_fragments),
  cluster_columns=F,
  border=T))
dev.off()

print (read.csv ('TFmodule_enrichments_FDR0.05.csv'))
print (head (read.csv ('TF_in_module_enrichments_FDR0.05.csv'), 20))
