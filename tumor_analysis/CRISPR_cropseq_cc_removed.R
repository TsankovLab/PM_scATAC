conda activate scrnatools
R
source ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scrna_pipeline/load_libraries.R')
source ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scrna_pipeline/ggplot_aestetics.R')
source ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scrna_pipeline/useful_functions.R')
setwd ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CRISPR_cropseq_analysis/_cellranger_raw_Filter_400_800_25/no_harmony')

setwd("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CRISPR_cropseq_analysis/_cellranger_raw_Filter_400_800_25/no_harmony")
srt = readRDS ('srt_filtered.rds')

scrna_pipeline_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scrna_pipeline'

### Remove proliferating cells ####
cc_threshold = -0
srt$cycling = ifelse (srt$cc > cc_threshold, 'proliferating','rest')
pdf (file.path ('Plots','proliferation_distribution.pdf'))
hist (srt$cc, breaks=100)
abline (v = cc_threshold, col='red')
dev.off()

bp = cellComp(
  seurat_obj = srt,
  metaGroups = c('crispr_calls','cycling'),
  plot_as = 'bar',
  prop = F
  #pal = palette_celltype_lv1
  ) + gtheme
bp1 = cellComp(
  seurat_obj = srt,
  metaGroups = c('crispr_calls','cycling'),
  plot_as = 'bar',
  prop = T
  #pal = palette_celltype_lv1
  ) + gtheme
bp$data$crispr_calls = factor(bp$data$crispr_calls, levels = unique(bp1$data$crispr_calls[order(-bp1$data$Freq[bp1$data$cycling == 'proliferating'])]))
pdf (file.path('Plots',paste0('cellcycle_barplots.pdf')), height=6, width=12)
bp
dev.off()


srt = srt[,srt$cc < cc_threshold]
#srt = srt[, srt$crispr_calls %in% names(table (srt$crispr_calls)[table (srt$crispr_calls) > 9])]


# Import sarcomatoid module ####
library(readxl)
cnmf_spectra_unique_comb_full = as.list (read_excel( "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/MPM_naive_study/reproduction2/PM_scRNA_atlas/data/cnmf_per_compartment.xlsx", sheet = "Cms_full"))
cnmf_spectra_unique_comb_full = lapply (cnmf_spectra_unique_comb_full, function(x) na.omit (x[x != 'NA']))
# compute module scores ####
if (!all (names(cnmf_spectra_unique_comb_full) %in% colnames(srt@meta.data))) {
srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = lapply(cnmf_spectra_unique_comb_full, function(x) head(x,50)),
        cor_threshold = NULL, 
        pos_threshold = NULL, # 
        listName = 'Cm', outdir = paste0(projdir,'Plots/'))
}

library(tidyverse)
library(pheatmap)

df = as.data.frame (srt@meta.data) 
#df$merged_call[is.na (df$merged_call)] = 'NT'
#df = df[!is.na(df$merged_call),]
#top_kds = names(head(table (df$merged_call)[order(-table(df$merged_call))],50))
#df = df[df$merged_call %in% top_kds,]

# Extract Cm columns
cm_cols <- grep("^Cm\\d+", colnames(df), value = TRUE)

# Get mean expression for NTC_1
ntc_df1 = df[df$merged_call == "NTC", ]
# ntc_df3 = df[df$merged_call == "NTC_3", ]
ntc_means1 = colMeans(ntc_df1[, cm_cols], na.rm = TRUE)
# ntc_means3 = colMeans(ntc_df3[, cm_cols], na.rm = TRUE)
#ntc_means = rowMeans (cbind(ntc_means1, ntc_means3))
ntc_means = ntc_means1

# Get all unique crispr_calls (including NTC_1)
all_calls <- unique(df$merged_call)
# Initialize matrix
diff_matrix <- matrix(nrow = length(all_calls), ncol = length(cm_cols))
rownames(diff_matrix) <- all_calls
colnames(diff_matrix) <- cm_cols

# Compute difference: NTC_1 mean - crispr_call mean
for (call in all_calls) {
  group_df <- df[df$merged_call == call, ]
  group_means <- colMeans(group_df[, cm_cols], na.rm = TRUE)
  
  diff_matrix[call, ] <- ntc_means - group_means
}

# Initialize p-value matrix
pval_matrix <- matrix(NA, nrow = length(all_calls), ncol = length(cm_cols))
rownames(pval_matrix) <- all_calls
colnames(pval_matrix) <- cm_cols

# Loop over calls and Cm columns
for (call in all_calls) {
  
  group_df <- df[df$merged_call == call, ]
  
  for (cm in cm_cols) {
    control_vals <- df[df$merged_call %in% c("NTC"), cm]
    group_vals   <- group_df[[cm]]
    
    # Check if both groups have >1 unique values
    if (length(unique(control_vals)) > 1 && length(unique(group_vals)) > 1) {
      pval_matrix[call, cm] <- wilcox.test(control_vals, group_vals, exact = FALSE)$p.value
    } else {
      pval_matrix[call, cm] <- NA
    }
  }
}
#diff_matrix[diff_matrix > 0.2] = 0.2
#diff_matrix[diff_matrix < -0.2] = -0.2
ha = HeatmapAnnotation (cells = anno_barplot(as.vector(unname(table (srt$merged_call)[rownames(diff_matrix)]))))
pdf (file.path ('Plots','prolif_removed_knock_downs_heatmap3.pdf'), height=5, width=4)
Heatmap (t(diff_matrix), row_names_gp=gpar(fontsize = 7), top_annotation=ha,
  column_split = ifelse (rownames (diff_matrix) == 'NTC','control','KOs'),
  column_names_rot=45,
  border=T,
  column_names_gp = gpar (fontsize=7),
  cell_fun = function (j, i, x, y, width, height, fill) 
            {
           if (t(pval_matrix)[i, j] < 0.001)
              {
               grid.text("***", x, y, just='center', vjust=.8,
                gp = gpar(fontsize = 5, col='black'))
              } else {
              if(t(pval_matrix)[i, j] < 0.01)
                  {
                  grid.text("**", x, y, just='center', vjust=.8,
                  gp = gpar(fontsize = 5, col='black'))   
                  } else {
                  if(t(pval_matrix)[i, j] < 0.05)
                    {
                    grid.text("*", x, y, just='center', vjust=.8,
                    gp = gpar(fontsize = 5, col='black'))         
                    }}}
      })
dev.off()





#kds = kds[30:length (kds)]
srt$celltype3 = 'celltype3'
top_kds = unique(srt$merged_call)
top_kds2 = top_kds[!top_kds %in% c('NTC')]
for (kd in top_kds2)
  {
  force = T
  do.fgsea = TRUE
  rankby = 'LFC' # Ranking to input in fgsea can be 'LFC' or 'pval_signedLFC'
  logfcThreshold = 0.25 # min logFC threshold to consider for testing. Should be set on 0
  pvalAdjTrheshold = 0.05
  topGenes = 20
  addGene=NULL # Add gene(s) to include in the DEG heatmap
  FeatureSets = list (all = NULL)  # Specify subset of features to run the test on
  metaGroupName1 = 'celltype3' # Specify metaGroup for celltype / clustering
  metaGroupName2 = 'merged_call' # Specify metaGroup of condition to compare
  deg2Ident = c('NTC', kd) #vector of 2 elements specifying the comparison between two groups from metaGroupName2
  top_pathways = Inf
  source (file.path(scrna_pipeline_dir,'DEG2.R'))
  }

rankby = 'LFC'

gmt_annotations = c(
#'c2.cp.kegg.v7.1.symbol.gmt',
'c2.cp.reactome.v7.1.symbol.gmt',
'c5.bp.v7.1.symbol.gmt',
'h.all.v7.1.symbol.gmt'
)

fgseaL = list()
for (kd in top_kds)
  {
  projdir_deg2 = paste0 ('deg2_',metaGroupName2,'_', metaGroupName1,'_',deg2Ident[[1]],'_vs_',kd,'_lfc_',logfcThreshold)
  if (file.exists(file.path(projdir_deg2,paste0('fGSEA_annotation_',paste(gmt_annotations,collapse='_'),'_rankby_',rankby,'.rds'))))
  fgseaL[[kd]] = readRDS (file.path(projdir_deg2, paste0('fGSEA_annotation_',paste(gmt_annotations, collapse='_'),'_rankby_',rankby,'.rds')))  
  }

#ann
#ann_terms = unique(unlist(lapply (fgseaL, function(x) x[[ann]][[1]]$pathway)))

for (ann in gmt_annotations)
  {
  #ann = 'c5.bp.v7.1.symbol.gmt'
  fgseaL_ann = lapply (fgseaL, function(x) x[[ann]])
  fgseaL_ann = unlist(fgseaL_ann, recursive=F)
  pvalAdjTrheshold = 0.05
  top_pathways = Inf
  fgseaResAll_dp = dotGSEA (fgseaL_ann, padj_threshold = pvalAdjTrheshold, 
      type = 'fgsea',top_pathways = top_pathways,
      cluster_rows=T,
      cluster_cols=T)
    
  pdf (file.path('Plots',paste0('fGSEA_annotation_',ann,'_filtered_proliferation_dotplot.pdf')), width=6.5, height=5)
  print(fgseaResAll_dp)
  dev.off()
  }  




guides = unique(srt$merged_call)
guides = guides[!guides %in% 'NTC']
srt$celltype = 'crispr'
gdot_p2 = geneDot (
  seurat_obj = srt, 
  gene = unname(unlist(guides)),
  y = 'celltype',
  x = 'merged_call',
  scale.data = T,
  include_NA = FALSE,
  min_expression = 0,
  returnDF = F,
  plotcol = as.character(rev(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 100))))
gdot_p2$data = gdot_p2$data %>%
  filter(x_axis == "NTC" | gene == x_axis)
gdot_p2$data$x_axis = ifelse(gdot_p2$data$x_axis == 'NTC','control', 'guide')

dp = DotPlot (srt, group.by = 'merged_call', features = c('SOX9'), scale =F)
dp_data1 = dp$data[dp$data$id %in% c('NTC', 'SOX9'),]
dp_data2 = gdot_p2$data[gdot_p2$data$y_axis %in% c('NTC', 'SOX9'),]
pdf (file.path('Plots','guides_expression_cc_removed_dotplot.pdf'), height=3.5, width=3)
#gdot_p2
gdot_p2
dev.off()




# import TF targets MsigDB database ####
tf_db =  read.gmt('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/c3.tft.v7.1.symbol.gmt')
tf_db = split (tf_db, tf_db$term)
guides_targets = lapply (unique(srt$merged_call), function(x) which (grepl (x, names (tf_db))))
guides_targets = lapply (guides_targets, function(x) x[1])
guides_targets = unlist (guides_targets[!is.na(guides_targets)])
tf_db = tf_db[guides_targets]
tf_db = lapply (tf_db, function(x) x[,2])
tf_db = tf_db[names (tf_db) != 'ATCMNTCCGY_UNKNOWN']
# Import Harmonizome databases for the missing TCF3 TEAD4 TWIST1 and MEF2A
d1 = read.table ('../../../git_repo/Harmonizome_TF_databases/CHEA Transcription Factor Targets/gene_attribute_edges.txt', header=T)
# d11 = read.table ('../../../git_repo/Harmonizome_TF_databases/CHEA Transcription Factor Targets/attribute_list_entries.txt', header=T)
d1 = read.table ('../../../git_repo/Harmonizome_TF_databases/ENCODE Transcription Factor Targets/gene_attribute_edges.txt', header=T)
d1 = read.table ('../../../git_repo/Harmonizome_TF_databases/JASPAR Predicted Transcription Factor Targets/gene_attribute_edges.txt', header=T)
#tf_db = c(tf_db, split(d1$source[d1$target %in% c('TCF3','TWIST1','MEF2A','TEAD4')], d1$target[d1$target %in% c('TCF3','TWIST1','MEF2A','TEAD4')]))
#vf = VariableFeatures(FindVariableFeatures(srt, nfeat= 10000))
tf_db = lapply (tf_db, function(x) x[x %in% vf])
srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = tf_db,
        cor_threshold = NULL, 
        pos_threshold = NULL, # 
        listName = 'TF_targets', outdir = paste0(projdir,'Plots/'))

ccomp = srt@meta.data
head (ccomp)
ccomp = ccomp[, c(names (tf_db), 'merged_call')]
ccomp = gather (ccomp, TF, module, 1:(ncol(ccomp)-1))
ccomp$guide_control = ifelse (ccomp$merged_call == 'NTC', 'control','guide')
ccomp$TF = sapply (ccomp$TF, function(x) unlist(strsplit(x, '_'))[1])
ccomp = ccomp[ccomp$merged_call == ccomp$TF | ccomp$merged_call == 'NTC',]
ccomp$merged_call[ccomp$merged_call == 'NTC'] = ccomp$TF[ccomp$merged_call == 'NTC']

bp = ggplot (ccomp, aes (x = merged_call, y = module)) + geom_boxplot (aes (fill = guide_control)) + gtheme

pdf (file.path ('Plots','TF_targets_cc_removed_boxplots.pdf'),5,5)
bp
dev.off()






#### Show expression of guides in relative genes ####
guides_mod = unique (srt$merged_call)
guides_mod = guides_mod[guides_mod != 'NTC']
srt$control_guide = ifelse (srt$merged_call == 'NTC','control','guide')
ccomp = srt@meta.data
ccomp = ccomp[, !colnames(ccomp) %in% guides_mod]
sub_data = t(srt@assays$RNA@layers$data[rownames(srt) %in% unlist(guides_mod),])
colnames (sub_data) = rownames (srt) [rownames(srt) %in% unlist(guides_mod)]
ccomp = cbind (ccomp[,c('merged_call','control_guide')], sub_data)
ccomp2 = gather (ccomp, guide, expression, 3:(ncol(ccomp)))
ccomp2 = ccomp2[ccomp2$merged_call == ccomp2$guide | ccomp2$merged_call == 'NTC', ]
ccomp2$merged_call[ccomp2$merged_call == 'NTC'] = ccomp2$guide[ccomp2$merged_call == 'NTC']
gp1 = ggplot (ccomp2, aes (x = merged_call, y = expression)) + geom_violin (
  aes (fill = control_guide),
  trim=TRUE,size=2,
    width=1,
    scale='width',
    linewidth = .3, alpha=.6) + 
geom_boxplot (aes (fill = control_guide),
    linewidth = .3,
    width=1,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.8, alpha=.6
     )  + gtheme +
scale_fill_manual (values = c(control = 'grey', guide='brown'))

stat.test <- ccomp2 %>%
  group_by (merged_call) %>%
  wilcox_test(reformulate('control_guide', 'expression')) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
#stat.test = stat.test[stat.test$group1 == 'NTC',]  
stat.test <- stat.test %>% add_xy_position (x = "merged_call", step.increase=0.3)
gp1 = gp1 + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
  bracket.nudge.y = -1, hide.ns = TRUE,
   label = "p.adj.signif") 
# gp2 = ggplot (ccomp, aes (x = merged_call, y = G2M.Score)) + geom_boxplot()
# gp3 = ggplot (ccomp, aes (x = merged_call, y = S.Score)) + geom_boxplot()

gene = c('SOX9','SNAI2','RUNX2','RUNX1','SOX6','TWIST1')

pdf (file.path ('Plots','target_genes_cc_removed_violinplot.pdf'), height=4, width=5)
gp1
#DotPlot (srt, features = gene, group.by = 'merged_call')
#DotPlot (srt, feature = unname (unlist(guides_mod)), group.by = 'control_guide')
dev.off()







# ### Check expression of TF target genes inferred from scATAC-seq ####
# target_genes = readRDS ('../../../tumor_compartment/scatac_ArchR/genes_with_tf_peaks.rds')
# var_feat = VariableFeatures (srt)
# target_genes = lapply (target_genes, function(x) x[x %in% var_feat])
# guides = names (target_genes)
# names (target_genes) = paste0(names (target_genes), '_target_genes')

# srt = ModScoreCor (
#         seurat_obj = srt, 
#         geneset_list = target_genes,
#         cor_threshold = NULL, 
#         pos_threshold = NULL, # 
#         listName = 'tg', outdir = paste0(projdir,'Plots/'))

# ccomp = srt@meta.data
# bp = lapply (guides, function(x)
# {
# bp = ggplot (ccomp[ccomp$merged_call %in% c('NTC',x),], aes_string (
#   x = 'merged_call', y = paste0(x,'_target_genes'))) + 
# geom_boxplot () + gtheme
# bp$data$merged_call = factor(bp$data$merged_call, levels = c('NTC',x))
# bp
# })


# pdf (file.path ('Plots','target_gene_modules_dotplot.pdf'))
# wrap_plots (bp)
# dev.off()

srt$comb = combgenes (srt, genes = c('RUNX2','SNAI2'))

bp = cellComp(
  seurat_obj = srt,
  metaGroups = c('merged_call','comb'),
  plot_as = 'bar',
  prop = T
  #pal = palette_celltype_lv1
  ) + gtheme


pdf (file.path ('Plots','donstream_cc_removed_TF_SOX9_prop.pdf'))
bp
dev.off()






gene = c('SOX9','SNAI2','RUNX2','SNAI3','NFYC','NFIC')

pdf (file.path ('Plots','SOX9_target_genes_dotplot.pdf'), height=4, width=6)
DotPlot (srt[,srt$merged_call %in% c('SOX9','NTC')], features = gene, group.by = 'merged_call') + gtheme
dev.off()

### Assess enrichment of Sox9 targets from Fuchs paper ####

srt$celltype_cc_filtered = 'celltype_cc_filtered'
top_kds = 'SOX9'
top_kds2 = top_kds[!top_kds %in% c('NTC')]
#srt$crispr_calls2[srt$crispr_calls2 %in% c('NTC')] = 'NTC'
for (kd in top_kds2)
  {
  force = F
  do.fgsea = TRUE
  rankby = 'LFC' # Ranking to input in fgsea can be 'LFC' or 'pval_signedLFC'
  rankby = 'pval_signedLFC'
  logfcThreshold = 0.25 # min logFC threshold to consider for testing. Should be set on 0
  pvalAdjTrheshold = 0.05
  topGenes = 20
  addGene=NULL # Add gene(s) to include in the DEG heatmap
  FeatureSets = list (all = NULL)  # Specify subset of features to run the test on
  metaGroupName1 = 'celltype_cc_filtered' # Specify metaGroup for celltype / clustering
  metaGroupName2 = 'merged_call' # Specify metaGroup of condition to compare
  deg2Ident = c('NTC', kd) #vector of 2 elements specifying the comparison between two groups from metaGroupName2
  top_pathways = Inf
  source (file.path(scrna_pipeline_dir,'DEG2.R'))
  }

rankby = 'pval_signedLFC'

sox9_gs = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/gene_sets/Sox9_targets_mouse.csv')
sox9_gs_W1 = sox9_gs$symbol[sox9_gs$W1 > 0]
sox9_gs_W2 = sox9_gs$symbol[sox9_gs$W2 > 0]
sox9_gs_W6 = sox9_gs$symbol[sox9_gs$W6 > 0]
sox9_gs_W12 = sox9_gs$symbol[sox9_gs$W12 > 0]
sox9_gs_class = split (sox9_gs$symbol, sox9_gs$class)
pathways = sox9_gs_class
pathways = lapply (pathways, function(x) mouseMan (x))
pathways = lapply (pathways, function(x) x$HGNC.symbol)

deg2Cluster = deg2
if (rankby == 'pval_signedLFC') fgsea_ranks = -log10 (deg2Cluster$p_val + 1e-300) * sign (deg2Cluster$avg_log2FC)
if (rankby == 'LFC') fgsea_ranks = deg2Cluster$avg_log2FC
#fgsea_ranks = -log10 (deg2Cluster$p_val + 1e-300) * sign (deg2Cluster$avg_log2FC)

fgsea_ranks = setNames (fgsea_ranks, deg2Cluster$gene)    
fgsea_ranks = fgsea_ranks[fgsea_ranks != 0]
#degCluster = degCluster[degCluster$p_val_adj < 01,]
#fgsea_ranks = -log10(degCluster$p_val_adj) * sign (degCluster$avg_log2FC)
#fgsea_ranks = (-log10(degCluster$p_val_adj) + 1e-10) + abs(degCluster$avg_log2FC) * sign(degCluster$avg_log2FC)
#fgsea_ranks = degCluster$avg_log2FC
fgseaRes = fgseaMultilevel (pathways, 
    fgsea_ranks, 
#    minSize=15, 
#    maxSize=500, # changed this from 1500 to 1000 cause it generated an error
    BPPARAM = NULL)
#fgseaResCol = collapsePathways (fgseaRes, stats = fgsea_ranks, pathway = pathways)
#mainPathways = fgseaRes[fgseaRes$pathway %in% fgseaResCol$mainPathways]

fgsea_ranks

message ('Print enrichment plots for each signficant pathway and cell type')
  
        sig_pathways = fgseaRes$pathway[fgseaRes$padj < pvalAdjTrheshold]
        #gmt.file = paste0 ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/',org,'/',names(fgseaResAll)[x])
        #pathways = gmtPathways (gmt.file)
  pdf()
  ep = lapply (sig_pathways, function(z) {
           print (plotEnrichment(pathways[[z]],
           fgsea_ranks) + 
           labs(title=z))
        })
  dev.off()
  pdf (file.path('Plots','SOX9_Fuchs_enrichment_plots.pdf'),5,3)
  print(ep)
  dev.off() 



pdf (file.path ('Plots','SOX9_target_genes_cc_removed_W2to6_dotplot.pdf'), height=4, width=506)
DotPlot (srt[,srt$merged_call %in% c('SOX9','NTC')], features = unique(pathways$W12), group.by = 'merged_call') + gtheme
dev.off()






#### Run cNMF ####
nfeat = 3000
nfeat = 5000
force=F
k_list = c(5:10)
k_selections = c(20:30)
k_list = c(20:30)
cores= 100
cnmf_name = 'crispr'

### RUN consensus NMF ####
source (file.path (scrna_pipeline_dir, 'cnmf_prepare_inputs.R'))
  
### Import and format spectra files ####
k_selection = 30
cnmf_name = 'crispr'
source (file.path (scrna_pipeline_dir,'cnmf_format_spectra_files.R')) 

#lapply (cnmf_spectra_unique, function(x) head(x, 20))

srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = lapply(cnmf_spectra_unique, function(x) head(x,50)),
        cor_threshold = NULL, 
        pos_threshold = NULL, # 
        listName = 'Cm', outdir = paste0(projdir,'Plots/'))

#ccomp = srt@meta.data[,c(names (cnmf_spectra_unique),names (cnmf_spectra_unique_comb_full))]
ccomp = srt@meta.data[,names (cnmf_spectra_unique)]
ccomp_cor = cor (ccomp)
pdf (file.path ('Plots','cnmf_cor_cc_removed_heatmap.pdf'))
Heatmap (ccomp_cor)
dev.off()

ccomp = srt@meta.data
bp = lapply (names (cnmf_spectra_unique), function(x) 
  ggplot (ccomp, aes_string (x = 'merged_call', y = x)) + geom_boxplot () + gtheme)

pdf (file.path ('Plots','cNMF_cc_removed_cc_removed_boxplots.pdf'),20,20)
wrap_plots (bp)
dev.off()

sapply (unique (srt$merged_call), function(x) sapply (cnmf_spectra_unique, function(y) x %in% y))



# Compare with cnmf from patient samples
ov = ovmat (c(cnmf_spectra_unique_comb_full, cnmf_spectra_unique), compare_lists=list(names(
  cnmf_spectra_unique_comb_full), names (cnmf_spectra_unique)))

pdf (file.path ('Plots','cellcycle_overlap_cc_removed_cnmf_heatmap.pdf'))
ov
dev.off()



### Check distributions of Cms after cc correction
ccomp = srt@meta.data
ccomp2 = srt2@meta.data

bp = lapply (names (cnmf_spectra_unique_comb_full), function(x) 
  ggplot (ccomp, aes_string (x = 'crispr_calls', y = x)) + geom_boxplot () + gtheme)

bp2 = lapply (names (cnmf_spectra_unique_comb_full), function(x) 
  ggplot (ccomp2, aes_string (x = 'crispr_calls', y = x)) + geom_boxplot () + gtheme)

pdf (file.path ('Plots','Cm_cc_removed_boxplots.pdf'),20,20)
wrap_plots (bp)
wrap_plots (bp2)
dev.off()








# import TF targets MsigDB database ####
tf_db = read.gmt('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/c3.tft.v7.1.symbol.gmt')
tf_db = split (tf_db, tf_db$term)
guides_targets = lapply (unique(srt$merged_call), function(x) which (grepl (x, names (tf_db))))
guides_targets = lapply (guides_targets, function(x) x[1])
guides_targets = unlist (guides_targets[!is.na(guides_targets)])
tf_db = tf_db[guides_targets]
tf_db = lapply (tf_db, function(x) x[,2])
tf_db = tf_db[names (tf_db) != 'ATCMNTCCGY_UNKNOWN']
# Import Harmonizome databases for the missing TCF3 TEAD4 TWIST1 and MEF2A
d1 = read.table ('../../../git_repo/Harmonizome_TF_databases/CHEA Transcription Factor Targets/gene_attribute_edges.txt', header=T)
# d11 = read.table ('../../../git_repo/Harmonizome_TF_databases/CHEA Transcription Factor Targets/attribute_list_entries.txt', header=T)
d1 = read.table ('../../../git_repo/Harmonizome_TF_databases/JASPAR Predicted Transcription Factor Targets/gene_attribute_edges.txt', header=T)
d1 = read.table ('../../../git_repo/Harmonizome_TF_databases/ENCODE Transcription Factor Targets/gene_attribute_edges.txt', header=T)
tf_db = c(tf_db, split(d1$source[d1$target %in% c('TCF3','TWIST1','MEF2A','TEAD4')], d1$target[d1$target %in% c('TCF3','TWIST1','MEF2A','TEAD4')]))

# Add signature from Fuchs paper 
sox9_gs = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/gene_sets/Sox9_targets_mouse.csv')
sox9_gs_W1 = sox9_gs$symbol[sox9_gs$W1 > 1]
sox9_gs_W2 = sox9_gs$symbol[sox9_gs$W2 > 1]
sox9_gs_W6 = sox9_gs$symbol[sox9_gs$W6 > 1]
sox9_gs_W12 = sox9_gs$symbol[sox9_gs$W12 > 1]
sox9_gs_class = split (sox9_gs$symbol, sox9_gs$class)
pathways = sox9_gs_class
pathways = lapply (pathways, function(x) mouseMan (x))
pathways = lapply (pathways, function(x) x$HGNC.symbol)
names (pathways)[2] = 'D0toW2'
tf_db = c (tf_db, pathways)

### Get TWIST1 targets from paper
twist1 = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/TWIST1_targets_PMID27546376.csv')
twist1 = mouseMan (twist1[,1])
tf_db = c (tf_db, list(TWIST1 = twist1$HGNC.symbol))

### Get TEAD4 targets from paper
tead4 = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/TEAD4_targets_PMID21701496.csv')
tead4 = mouseMan (tead4$Primary.Annotation)
tf_db[['TEAD4']] = tead4$HGNC.symbol

vf = VariableFeatures(FindVariableFeatures(srt, nfeat= 10000))
tf_db = lapply (tf_db, function(x) x[x %in% vf])

### Remove genes associated to cell cycle in each regulon ####
norm_mat = srt@assays$RNA@layers$data
cc_cor = cor (srt$Cm16, t(as.matrix(norm_mat)))
colnames (cc_cor) = rownames(srt)

tf_db = lapply (tf_db, function(x) x[x%in%rownames (srt)])
cor_res = sapply (tf_db, function(x) 
  {
  cor_res = sapply (x, function(y) cor.test (norm_mat[rownames(srt) == y,], unname(srt$Cm16))$p.value)
  names(cor_res)[cor_res > 0.01]
  })
cor_res = lapply (cor_res, function(x) x[!is.na(x)])
tf_db = cor_res

srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = tf_db,
        cor_threshold = NULL, 
        pos_threshold = NULL, # 
        listName = 'TF_targets', outdir = paste0(projdir,'Plots/'))

ccomp = srt@meta.data
head (ccomp)
ccomp = ccomp[, c(names (tf_db), 'merged_call')]
ccomp = gather (ccomp, TF, module, 1:(ncol(ccomp)-1))
ccomp$guide_control = ifelse (ccomp$merged_call == 'NTC', 'control','guide')
ccomp$TF = sapply (ccomp$TF, function(x) unlist(strsplit(x, '_'))[1])
ccomp$merged_call[ccomp$TF == 'D0' & ccomp$merged_call %in% c('SOX9','NTC')] = 'SOX9.1' 
ccomp$merged_call[ccomp$TF == 'D0toW2' & ccomp$merged_call %in% c('SOX9','NTC')] = 'SOX9.0' 
ccomp$merged_call[ccomp$TF == 'W1to6' & ccomp$merged_call %in% c('SOX9','NTC')] = 'SOX9.2' 
ccomp$merged_call[ccomp$TF == 'W2to6' & ccomp$merged_call %in% c('SOX9','NTC')] = 'SOX9.3' 
ccomp$merged_call[ccomp$TF == 'W12' & ccomp$merged_call %in% c('SOX9','NTC')] = 'SOX9.4'
ccomp$TF[ccomp$TF == 'D0' & ccomp$merged_call == 'SOX9.1'] = 'SOX9.1' 
ccomp$TF[ccomp$TF == 'D0toW2' & ccomp$merged_call == 'SOX9.0'] = 'SOX9.0' 
ccomp$TF[ccomp$TF == 'W1to6' & ccomp$merged_call == 'SOX9.2'] = 'SOX9.2' 
ccomp$TF[ccomp$TF == 'W2to6' & ccomp$merged_call == 'SOX9.3'] = 'SOX9.3' 
ccomp$TF[ccomp$TF == 'W12' & ccomp$merged_call == 'SOX9.4'] = 'SOX9.4'
ccomp = ccomp[ccomp$merged_call == ccomp$TF | ccomp$merged_call == 'NTC',]
ccomp$merged_call[ccomp$merged_call == 'NTC'] = ccomp$TF[ccomp$merged_call == 'NTC']
#ccomp = ccomp[!ccomp$merged_call %in% c('SOX9','SOX9')]
bp = ggplot (ccomp, aes (x = merged_call, y = module)) + geom_boxplot (aes (fill = guide_control)) + gtheme

pdf (file.path ('Plots','TF_targets_cc_genes_removed_boxplots.pdf'),5,5)
bp
dev.off()






