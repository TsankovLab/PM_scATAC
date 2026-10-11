# R2_Q21 -- add the nearest gene to the table of peaks more accessible in CD8 exhausted than in CD8 T cells
#
# The table is the one exported as DAP_CD8_Ex.csv by git_repo/tnk_analysis/scatac_NKTcells_ArchR.R
# ("Differential Peaks in CD8 exhausted vs CD8"): getMarkerFeatures on the PeakMatrix, CD8_exhausted vs CD8,
# Wilcoxon, k = 100, not binarised; significant there = FDR < 0.01 & Log2FC > 0 (418 peaks).
# Here the cut-off is FDR < 0.05 & Log2FC > 0 (884 peaks), the peak set of the labelled heatmap.
# Here the saved test result is read back (nothing is recomputed, the T/NK project is not modified) and each
# peak gets a gene by the rule of the labelled heatmap (Peak2Gene link, else nearest gene) plus the ArchR
# peak-set annotation (nearest gene, distance to its TSS, peak type).

suppressPackageStartupMessages({
  library(ArchR)
  library(SummarizedExperiment)
})

tnk_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR'

# Same steps as the original script ####
DAP_list = readRDS (file.path (tnk_dir, 'DAP_CD8_CD8_ext_pairwise.rds'))
DAP_res = do.call (cbind, (assays(DAP_list)))
colnames (DAP_res) = names(assays(DAP_list))
DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames(DAP_res) = as.character(DAP_res_regions)

# Cut-off harmonised with the labelled heatmap: FDR < 0.05 (884 peaks). The original script used FDR < 0.01 (418 peaks).
EXH_FDR = 0.05
DAP_res_sig = DAP_res[DAP_res$FDR < EXH_FDR & DAP_res$Log2FC > 0, ]

# Add nearest gene from the peak set ####
archp_rds = readRDS (file.path (tnk_dir, 'Save-ArchR-Project.rds'))   # read only, the project is not modified
ps = archp_rds@peakSet
names (ps) = NULL
ps_ann = as.data.frame (mcols(ps)[, c('nearestGene','distToTSS','peakType','distToGeneStart')])
rownames (ps_ann) = as.character (ps)
stopifnot (all (rownames(DAP_res_sig) %in% rownames(ps_ann)))

# Gene annotation consistent with the labelled heatmap ####
# tnk_exhaustion_annotation/Plots/CD8_ext_peakset_focused_FGFBP2_CD8__KLRC1_CD8ex_FDR0.05.pdf labels peaks with
# the gene they are linked to by Peak2Gene (correlation >= 0.35), not with the nearest gene, and falls back to
# proximity only for a few genes without a link. The same rule is applied to every peak of the table:
#   gene = heatmap label if the peak is labelled in the heatmap,
#          else the Peak2Gene-linked gene with the highest correlation (>= 0.35),
#          else the nearest gene.
hm_dir  = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/tnk_exhaustion_annotation'
COR_CUT = 0.35   # as in annotate_CD8ext_heatmap.R

# Peak2Gene links stored in the T/NK project by scatac_NKTcells_ArchR.R (addPeak2GeneLinks, GeneScoreMatrix,
# IterativeLSI, maxDist 250 kb), filtered as getPeak2GeneLinks (corCutOff = 0.2) does: correlation >= 0.2,
# FDR <= 1e-4, peak and gene in the top 75% by variance. This is the content of the p2g_links.rds used for the heatmap.
p2g_all_links = metadata (archp_rds@peakSet)$Peak2GeneLinks
keep = which (p2g_all_links$Correlation >= 0.2 & p2g_all_links$FDR <= 1e-4 & p2g_all_links$VarQATAC >= 0.25 & p2g_all_links$VarQRNA >= 0.25)
p2g = data.frame (
  peak = as.character (metadata (p2g_all_links)$peakSet)[p2g_all_links$idxATAC[keep]],
  gene = metadata (p2g_all_links)$geneSet$name[p2g_all_links$idxRNA[keep]],
  cor  = p2g_all_links$Correlation[keep],
  fdr  = p2g_all_links$FDR[keep])
p2g = p2g[p2g$cor >= COR_CUT & p2g$peak %in% rownames(DAP_res_sig), ]
p2g = p2g[order (p2g$peak, -p2g$cor), ]
p2g_top = p2g[!duplicated (p2g$peak), ]
rownames (p2g_top) = p2g_top$peak
p2g_all = sapply (split (paste0 (p2g$gene, ' (', round (p2g$cor, 2), ')'), p2g$peak), paste, collapse = '; ')

# heatmap labels: region_idx indexes the FDR < 0.05 gained peaks in test order
sig05 = rownames(DAP_res)[DAP_res$FDR < .05 & DAP_res$Log2FC > 0]
lab = read.csv (file.path (hm_dir, 'exhaustion_gene_labels_FDR0.05.csv'))
lab$peak = sig05[lab$region_idx]
stopifnot (all (abs (DAP_res[lab$peak, 'Log2FC'] - lab$Log2FC) < 1e-9))
stopifnot (all (lab$peak %in% rownames(DAP_res_sig)))   # every heatmap label is a row of the table
rownames (lab) = lab$peak
# every Peak2Gene label of the heatmap must be recovered from the links used here
stopifnot (all (mapply (function(p, g) g %in% p2g$gene[p2g$peak == p], lab$peak[lab$method == 'P2G'], lab$gene[lab$method == 'P2G'])))

pk = rownames (DAP_res_sig)
ann = ps_ann[pk, ]
ann$heatmap_label = lab[pk, 'gene']
ann$p2g_gene      = p2g_top[pk, 'gene']
ann$p2g_cor       = p2g_top[pk, 'cor']
ann$p2g_genes_all = unname (p2g_all[pk])
ann$gene = ifelse (!is.na (ann$heatmap_label), ann$heatmap_label, ifelse (!is.na (ann$p2g_gene), ann$p2g_gene, ann$nearestGene))
ann$gene_source = ifelse (!is.na (ann$heatmap_label), paste0 ('heatmap_', lab[pk, 'method']), ifelse (!is.na (ann$p2g_gene), 'P2G', 'nearestGene'))
ann = ann[, c('gene','gene_source','heatmap_label','p2g_gene','p2g_cor','p2g_genes_all','nearestGene','distToTSS','peakType','distToGeneStart')]

DAP_res_sig = cbind (ann, DAP_res_sig)
DAP_res_sig = DAP_res_sig[order (-DAP_res_sig$Log2FC), ]
# exported columns, in the order of the supplementary table
out_cols = c('Log2FC','Mean','FDR','Pval','MeanDiff','AUC','MeanBGD','gene','gene_source','p2g_gene','p2g_cor',
  'p2g_genes_all','nearestGene','distToTSS','peakType','distToGeneStart')
write.csv (cbind (region = rownames(DAP_res_sig), DAP_res_sig[, out_cols]), 'DAP_CD8_Ex_FDR0.05_gene_annotation.csv', row.names = FALSE)

message (nrow(DAP_res_sig), ' peaks; ', length(unique(DAP_res_sig$gene)), ' distinct genes')
print (table (DAP_res_sig$gene_source))
message ('gene differs from nearest gene in ', sum (DAP_res_sig$gene != DAP_res_sig$nearestGene | is.na(DAP_res_sig$nearestGene)), ' peaks')
hl = DAP_res_sig[!is.na(DAP_res_sig$heatmap_label), c('gene','gene_source','p2g_gene','p2g_cor','nearestGene')]
print (hl)
message ('heatmap labels where the top Peak2Gene gene is a different gene: ', sum (hl$gene_source == 'heatmap_P2G' & hl$gene != hl$p2g_gene))
print (head (sort (table (DAP_res_sig$gene), decreasing = TRUE), 15))
