conda activate scrnatools 
R

#data.dir = '/broad/hptmp/bgiotti/JiaMoPool/cellranger_output/RAPA05_JiaPool1_1_v1/filtered_feature_bc_matrix/'
#samples = c(paste0('GFP_',1:3), paste0('CB2_PTENL_',1:3), paste0('CB2_PTENL_C124S_',1:3))
#samples_path1 = paste0('/broad/hptmp/bgiotti/JiaMoPool/cellranger_output/RAPA05_JiaPool1_1_v1/',samples,'/gex/')
#samples_path2 = paste0('/broad/hptmp/bgiotti/JiaMoPool/cellranger_output/RAPA05_JiaPool1_2_v1/',samples,'/gex/')
samples_path =
  c('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/cropseq/raw/Sample1Cropseq/cellranger_output/ALTS06_sample1_0_v1/multi/count/raw_feature_bc_matrix.h5'
    )

meta = data.frame (
  sampleID = 'cropseq'
  )
meta = cbind (meta, sample_path = samples_path)

names(samples_path) = meta$sampleID
#meta = meta[!duplicated(meta$assignment),]
#samples_path = '/broad/hptmp/bgiotti/JiaMoPool/cellranger_output/RAPA05_JiaPool1_1_v1/raw_feature_bc_matrix/'
# Set project directory
proj_name = 'CRISPR_cropseq'
projdir_init = paste0("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/",proj_name,"_analysis")

# Set cellranger_to_seurat parameters
cr_to_seurat = list(
  run_cellbender = FALSE, # if this is set to TRUE then input raw (and not filtered) cellranger count matrices
  cellbender_samples = NULL,
  #cellbender_parameters = list(p811 = c(expected_cells = 12000,total_droplets_included = 100000, low_count_threshold = 5)),
  cellbender_parameters = NULL,
  #cellbender_dir = '/broad/hptmp/bgiotti/meso_cellbender/',
  org = 'human',
  datatype = 'RNA',
  cr_output = 'raw',
  samples_path = samples_path,
  meta = meta, 
  is.hashed = FALSE
  )

# Set QC parameters
qc_params = list(
  filtering = 'hard', # 'emptyDrops' or 'hard' filtering
  nFeat = 400,#400, # Number of features per cells. default 400
  nCounts = 800, #1000, # Number of UMI per cell. Default 800
  pchM = 25, # Percent mitochondrial genes. Default 25 
  remove.samples = NULL, # Remove bad samples. Takes vector of sampleIDs 
  processInd = FALSE # run preprocessing per sample before running it on the merged data
  )

### Data processing and clustering variables ###
harmony_params = list(
  batch = c('no')
  #batch = 'no'
  )

data_processing_param = list(
  variablefeatures = 'scran', # options are 'scran', 'seurat','manual'
  nfeat = 3000, # number of variable genes to consider for dimentionality reduction
  sigPCs = 15,
  vars_to_regress=NULL,
  metaGroupNames = c('sampleID'),
  res = c(0.2, 0.8, 2, 5, 10) # denovo cluster resolutions 
  )
    
# Initiate pipeline ####
force=FALSE # re run pipeline from the beginning no matter if objects are found
subclustername = NULL
scrna_pipeline_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scrna_pipeline'
source (file.path(scrna_pipeline_dir,'master_scrna.R'))


if (!file.exists ('../../_cellranger_raw_Filter_400_800_25/no_harmony/infercnv/all_samples_subsampled_Inf_ref_broad_meso/infercnv.results.obj.Rds'))
  {
  ### INFERCNV ####
  # Get mesothelium ref from broad normal lung data
  projdir_ref = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Normal_lung/normal_lung_scrna/'
  ref = readRDS (paste0(projdir_ref, 'mesothelium_srt.rds'))
  #saveRDS (ref, file.path (projdir_ref, 'mesothelium_srt.rds'))
  #subsample_ct = 1000
  #ref_bc = lapply (unique(ref$celltype), function(x) sample(colnames(ref)[ref$celltype == x], subsample_ct, replace=T))
  #names (ref_bc) = unique(ref$celltype)
  #cts = c('Mesothelium','Fibroblast')
  #ref = ref[,colnames (ref) %in% unlist(ref_bc[cts])]
  ref = ref[,ref$celltype == 'Mesothelium']
  ref$sampleID2 = 'normal_lung'
  ref_name = 'broad_meso'
  # srt$sampleID = paste0(srt$sampleID2, '_broad_ref')
  srt$celltype = srt$sampleID
  #srt$celltype2[srt$celltype2 == 'Fibroblasts'] = 'Malignant'
  metaGroupName1 = 'celltype'
  metaGroupName1Element = unique (srt$celltype)
  metaGroupName2 = 'sampleID'
  per_sample = FALSE
  # cancer_clusters = unique (srt$celltype_sample[grep ('malignant', srt$celltype_sample)])
  #subsample=min (table(srt$sampleID2))
  subsample=Inf
  force=F
  source (file.path(scrna_pipeline_dir, 'inferCNV_run.R'))
  } else {
  library (infercnv) 
  icnv = readRDS ('../../_cellranger_raw_Filter_400_800_25/no_harmony/infercnv/all_samples_subsampled_Inf_ref_broad_meso/infercnv.results.obj.Rds') # output folder
  }

# Filter cells based on correlation to average CNV profile ####
lexpr.data = log2(icnv@expr.data) 
lexpr.data = lexpr.data[,-icnv@reference_grouped_cell_indices[[1]]]
avg_cnv = rowMeans (lexpr.data)
cnv_cor = cor (lexpr.data, avg_cnv)
cnv_cor_filt = rownames(cnv_cor)[cnv_cor[,1] > 0.5]
srt = srt[, colnames (srt) %in% cnv_cor_filt]

pdf (file.path ('Plots','cnv_cor_distribution.pdf'))
hist (cnv_cor)
dev.off()


#### Map guides detected to cells ####
crispr_calls = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/cropseq/raw/Sample1Cropseq/cellranger_output/ALTS06_sample1_0_v1/per_sample_outs/ALTS06_sample1_0/count/crispr_analysis/protospacer_calls_per_cell.csv')
#crispr_ref = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/cropseq/raw/Sample1Cropseq/cellranger_output/ALTS06_sample1_0_v1/per_sample_outs/ALTS06_sample1_0/count/crispr_analysis/feature_reference.csv')
srt$crispr_calls = crispr_calls$feature_call[match (colnames(srt), crispr_calls$cell_barcode)]

# Identify cells with one gRNA called or in combination with NTC ####
x <- srt$crispr_calls

keep <- sapply(strsplit(x, "\\|"), function(parts) {
  # Rule 1: single gene
  if (length(parts) == 1) return(TRUE)
  
  # Check if NTC is present
  ntc_present <- any(grepl("^NTC", parts))
  
  if (ntc_present) {
    # Remove NTC from list
    other_genes <- parts[!grepl("^NTC", parts)]
    # After removing NTC, get unique base gene names
    gene_names <- sub("_\\d+$", "", other_genes)
    # Keep only if there is exactly ONE unique non-NTC gene
    return(length(unique(gene_names)) == 1)
  }
  
  # Rule 3: all same gene name ignoring _number
  gene_names <- sub("_\\d+$", "", parts)
  if (length(unique(gene_names)) == 1) return(TRUE)
  
  return(FALSE)
})
srt$unique_guide = keep
ccomp = as.data.frame (table (srt$unique_guide))
bp = ggplot (ccomp, aes (x = Var1, y = Freq, fill = Var1)) + geom_bar (stat = 'identity') + 
scale_fill_manual (values = c(`TRUE` = 'brown',`FALSE` = 'grey')) + gtheme
bp$data$Var1 = factor (bp$data$Var1, levels = c('TRUE','FALSE'))
pdf (file.path ('Plots','fraction_unique_guides_barplot.pdf'), width=3, height=3)
bp
dev.off()


### Filter non-unique guides ####
srt <- srt[,keep]


# # Remove cells with no guides ####
# srt = srt[,!is.na(srt$crispr_calls)]

res2 = as.data.frame (table (srt$crispr_calls))
res2 = res2[order(-res2$Freq),]
res2$Var1 = factor (res2$Var1, levels = unique (res2$Var1))
bp = ggplot (res2, aes (x = Var1, y = Freq)) + geom_bar (stat = 'identity') + gtheme

pdf (file.path ('Plots','gRNA_abundances_not_filtered_barplot.pdf'), width=93, height=20)
bp
dev.off()


res2 = as.data.frame (table (srt$crispr_calls))
res2 = res2[order(-res2$Freq),]
res2$Var1 = factor (res2$Var1, levels = unique (res2$Var1))
bp = ggplot (res2, aes (x = Var1, y = Freq)) + geom_bar (stat = 'identity') + gtheme

pdf (file.path ('Plots','gRNA_abundances_barplot.pdf'), width=13, height=5)
bp
dev.off()

### Further remove guides with less than n cells ####
# Show genes left in after also filtering by lowest amount of cells 
cells_threshold = 9
gene_left = as.character(res2$Var1[res2$Freq > cells_threshold])
srt = srt[,srt$crispr_calls %in% gene_left]

pdf (file.path ('Plots','gRNA_abundances_barplot.pdf'), width=10, height=5)
bp + geom_hline(yintercept = cells_threshold)
dev.off()

### Check consistency of guides ####
var_feat_avg = log2(AverageExpression (srt, group.by = 'crispr_calls', features = VariableFeatures (srt))[[1]]+1)
#var_feat_avg = t(scale(t(var_feat_avg)))
var_feat_avg_cor = cor(as.matrix(var_feat_avg), method='pearson')
colnames (var_feat_avg_cor) = gsub ('-','_', colnames(var_feat_avg_cor))
#ha = HeatmapAnnotation (cells = anno_barplot(as.numeric(table (srt$crispr_calls)[colnames(var_feat_avg_cor)])))
hm = Heatmap (var_feat_avg_cor[rownames(var_feat_avg_cor) != 'BPTF-2|BPTF-3',colnames(var_feat_avg_cor) != 'BPTF-2|BPTF-3'], 
 # top_annotation = ha, 
  clustering_distance_columns = 'pearson',
    clustering_distance_rows = 'pearson',
    row_names_gp = gpar (fontsize = 5), 
    column_names_gp = gpar (fontsize = 5),
    border=T
    )

pdf (file.path ('Plots','guides_cor_heatmap.pdf'),height = 4,width = 4.5)
hm
dev.off()

# Merge guides ####
merged_guides = list (
  MEF2D = c('MEF2D-1','MEF2D-3','MEF2D−3|NTC−3','MEF2D−1|MEF2D−3','MEF2D−1|NTC−3','MEF2D−1|MEF2D−3|NTC−3','MEF2D−1|NTC−1'),
  MEF2A = c('MEF2A−3','MEF2A−1','MEF2A−3|NTC−3','MEF2A−1|NTC−3','MEF2A−3|NTC−1'),
  PITX1 = c('PITX1−2','PITX1−3','PITX1−1','PITX1−2|PITX1−3','PITX1−1|PITX1−2','PITX1−1|PITX1−3','NTC-1|PITX1-2','NTC-3|PITX1-2','NTC-3|PITX1-3','NTC-1|PITX1-3'),
  BPTF = c('BPTF−3','BPTF−2'),
  HMGA1 = 'HMGA1−3',
  TEAD2 = 'TEAD2-1',
  NTC = c('NTC-1','NTC-3'),
  TCF3 = c('TCF3-1','TCF3-3'),
  SOX9 = 'SOX9-3',
  TWIST1 = c('TWIST1-1','TWIST1-3'),
  TEAD4 = c('TEAD4-2','TEAD4-3'))
merged_guides = lapply (merged_guides, function(x) gsub("[−–-]", "_", x))
srt$merged_call = srt$crispr_calls

for (i in names(merged_guides)) srt$merged_call[srt$merged_call %in% merged_guides[[i]]] = i
table (srt$merged_call)  

srt = srt[, !grepl ('\\|', srt$merged_call)]


### QC per guide ####
ccomp_df = as.data.frame (table(srt$merged_call))
cc_p2 = ggplot (ccomp_df, aes (x= Var1, y= Freq)) +
        geom_bar (position="stack", stat="identity") +
        theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1)) + 
        ggtitle (paste('Tot cells',ncol(srt))) + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
srt$nFeature_RNAL = log10 (srt$nFeature_RNA)
srt$nCount_RNAL = log10 (srt$nCount_RNA)
vln_p = VlnPlot (srt, features = c("nFeature_RNAL", "nCount_RNAL", "percent.mt"), combine=T,group.by = 'merged_call',pt.size = 0, ncol = 3)
png (file.path("Plots","QC_nFeat_nCount_m.percent_vlnPlot.png"), 3000, 1000, res=300)
print (cc_p2 | vln_p[[1]] | vln_p[[2]] | vln_p[[3]]) + plot_layout (widths=c(1,2,2,2))
dev.off()




### Show changes in variable genes ####
var_feat_avg = log2(AverageExpression (srt, group.by = 'merged_call', features = VariableFeatures(FindVariableFeatures(srt, nfeat = 5000)))[[1]]+1)
#var_feat_avg = t(scale(t(var_feat_avg)))
#var_feat_avg_cor = cor(as.matrix(var_feat_avg), method='pearson')
colnames (var_feat_avg) = gsub ('-','_', colnames(var_feat_avg))
var_scaled = scale(t(as.matrix(var_feat_avg)))
var_scaled[is.na(var_scaled)] = 0
var_scaled = var_scaled[,apply(var_scaled,2,var) > 0]
#var_scaled = var_scaled[!rownames(var_scaled) %in% c('HMGA1','TEAD2'),]
ncells = table (srt$merged_call)
#ncells = ncells[!names(ncells) %in% c('HMGA1','TEAD2')]

ha = HeatmapAnnotation (cells = anno_barplot(as.numeric(ncells[rownames(var_scaled)])), which='row')
hm = Heatmap (var_scaled, right_annotation = ha, 
  row_split = ifelse (rownames(var_scaled) == 'NTC','control','KO'),
  clustering_distance_columns = 'pearson',
  , clustering_distance_rows = 'pearson', 
  column_names_gp=gpar (fontsize=0),
  row_names_gp=gpar (fontsize=8),
  col = rev(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 100)),
  border=T)
pdf (file.path ('Plots','most_variable_genes_heatmap.pdf'),width=4,height=2.5)
hm
dev.off()


#### Show expression of guides in relative genes ####
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


### Try with dotplot
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
pdf (file.path('Plots','guides_expression_dotplot.pdf'), height=3.5, width=3)
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

# Import Harmonizome databases for the missing TCF3 TEAD4 TWIST1 and MEF2A
d1 = read.table ('../../../git_repo/Harmonizome_TF_databases/CHEA Transcription Factor Targets/gene_list_terms.txt', header=T)
d1 = read.table ('../../../git_repo/Harmonizome_TF_databases/CHEA Transcription Factor Targets/gene_list_terms.txt', header=T)
d1 = read.table ('../../../git_repo/Harmonizome_TF_databases/CHEA Transcription Factor Targets/gene_list_terms.txt', header=T)

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
ccomp = ccomp[ccomp$merged_call != 'ATCMNTCCGY',]
bp = ggplot (ccomp, aes (x = merged_call, y = module)) + geom_boxplot (aes (fill = guide_control)) + gtheme

pdf (file.path ('Plots','TF_targets_boxplots.pdf'),5,5)
bp
dev.off()







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




# Show lower Cm17 and higher Cm2 in SOX9 KD compared to control ####
module_genes = head(cnmf_spectra_unique_comb_full[['Cm17']],50)
module_genes = c('RUNX2','RUNX1','SNAI2')
#module_genes = rownames(srt)[grep('TGFB', rownames(srt))]
#module_genes = 'SNAI2'
module_genes_exp = t(srt@assays$RNA@layers$data[rownames(srt) %in% module_genes,])
rownames(module_genes_exp) = colnames(srt)
colnames(module_genes_exp) = rownames(srt)[rownames(srt) %in% module_genes]
ccomp = srt@meta.data[, c('merged_call'),drop=F]
ccomp = cbind (ccomp, module_genes_exp)
head (ccomp)
ccomp = ccomp[ccomp$merged_call %in% c('NTC','SOX9'),]
ccomp = gather (ccomp, gene, expression, 2:ncol(ccomp))
#ccomp[ccomp$crispr_calls == 'NTC_3'] = 'NTC_1'
gp1 = ggplot (ccomp, aes (x = gene, y = expression)) + 
 #geom_point (shape = 21, alpha=.5) +
  # geom_violin (trim=TRUE,size=2,
  #   width=1,
  #   scale='width',
  #   linewidth = .2, alpha=0.7) + 
  geom_boxplot (aes (fill = merged_call),
    linewidth = .2,
    width=.3,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.7
     ) + 
  scale_fill_manual (values = c(SOX9 = 'brown', NTC = 'grey60')) + 
  gtheme
  
# gp2 = ggplot (ccomp, aes (x = merged_call, y = Cm17, fill = merged_call)) + 
#   #geom_point (shape = 21, alpha=.5) +
#   # geom_violin (trim=TRUE,size=2,
#   #   width=1,
#   #   scale='width',
#   #   linewidth = .2, alpha=0.7) + 
#   geom_boxplot (
#     linewidth = .2,
#     width=.3,
#     outlier.alpha = 0.2,
#     outlier.size = 1,
#      size=0.6, alpha=0.7
#      ) + 
#   scale_fill_manual (values = c(SOX9 = 'brown', NTC = 'grey60')) + 
#   gtheme
  
  library (rstatix)
stat.test <- ccomp %>%
  group_by(gene) %>%
  filter(n_distinct(merged_call) == 2) %>%      # need two groups
  filter(n_distinct(expression) > 1) %>%        # need variation in expression
  wilcox_test(expression ~ merged_call) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")


stat.test <- stat.test %>%
  add_xy_position(x = "gene", dodge = 0.8)
gp1 <- gp1 +
  stat_pvalue_manual(
    stat.test,
    remove.bracket = FALSE,
    bracket.nudge.y = 0.1,   # nudge a bit above your boxes
    hide.ns = FALSE,
    label = "p.adj.signif"
  ) +
  NoLegend()
        
  #stat_ellipse(aes(color = crispr_calls), type = "norm", size = 1) +
# srt$SNAI2_prop = ifelse (srt$SNAI2 > 0, 'pos','neg')
# bp = cellComp (srt[,srt$merged_call %in% c('NTC','SOX9')],
#    metaGroups = c('merged_call','SNAI2_prop'),
#    #subset = 'CD8_exhausted',
#    prop=T,
#    plot_as = 'bar',
#    #pal = palette_tnk_cells
#    ) + gtheme
pdf (file.path ('Plots','Cm17_genes_in_SOX9.pdf'),width = 14,3)
#DotPlot (srt, group.by = 'merged_call', features= module_genes, scale=F) + gtheme
wrap_plots (gp1)
#bp
dev.off()

srt$comb = combgenes (srt, genes = c('RUNX2','SNAI2'))

bp = cellComp(
  seurat_obj = srt,
  metaGroups = c('merged_call','comb'),
  plot_as = 'bar',
  prop = T
  #pal = palette_celltype_lv1
  ) + gtheme


pdf (file.path ('Plots','donstream_TF_SOX9_prop.pdf'))
bp
dev.off()



### Check cell cycle genes ####
var_feat_avg = log2(AverageExpression (srt, group.by = 'merged_call', features = VariableFeatures(FindVariableFeatures(srt, nfeat = 5000)))[[1]]+1)
cc.genes <- readLines('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/gene_sets/regev_lab_cell_cycle_genes.txt')
  s.genes <- cc.genes[1:43]
  g2m.genes <- cc.genes[44:97]

var_scaled = scale(t(as.matrix(var_feat_avg)))
var_scaled[is.na(var_scaled)] = 0
var_scaled = var_scaled[,apply(var_scaled,2,var) > 0]
var_scaled = var_scaled[, colnames(var_scaled) %in% cc.genes]
hm = Heatmap (var_scaled, right_annotation = ha, 
  row_split = ifelse (rownames(var_scaled) == 'NTC','control','KO'),
  clustering_distance_columns = 'pearson',
  , clustering_distance_rows = 'pearson', 
  column_names_gp=gpar (fontsize=0),
  row_names_gp=gpar (fontsize=8),
  col = rev(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 100)),
  border=T)
pdf (file.path ('Plots','cc_genes_heatmap.pdf'),width=4,height=2.5)
hm
dev.off()

### Boxplot cc genes ####
ccomp = srt@meta.data[,c('cc','merged_call','S.Score','G2M.Score'), drop=F]
bxpv = geom_boxplot (
    linewidth = .3,
    width=.1,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.8, alpha=.6
     ) 

vlp = geom_violin (trim=TRUE,size=2,
    width=1,
    scale='width',
    linewidth = .3, alpha=.6)

ccomp2 = aggregate (ccomp[,'cc', drop=F], by = list(merged_call = srt$merged_call), mean)
ccomp$merged_call = factor (ccomp$merged_call, levels = ccomp2$merged_call[order(-ccomp2$cc)])

  

gp1 = ggplot (ccomp, aes (x = merged_call, y = cc)) + geom_violin (
  aes (fill = merged_call),
  trim=TRUE,size=2,
    width=1,
    scale='width',
    linewidth = .3, alpha=.6) + 
geom_boxplot (
    linewidth = .3,
    width=.1,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.8, alpha=.6
     )  + gtheme +
scale_fill_manual (values = c(NTC = 'grey', TWIST1='brown',
  TCF3 = 'brown', PITX1 = 'brown', TEAD4 = 'brown', TEAD2 = 'brown',
  MEF2A = 'brown', MEF2D = 'brown', SOX9 = 'brown', BPTF = 'brown', HMGA1 = 'brown'))

stat.test <- ccomp %>%
  wilcox_test(reformulate('merged_call', 'cc')) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test = stat.test[stat.test$group1 == 'NTC',]  
stat.test <- stat.test %>% add_xy_position (x = "merged_call", step.increase=0.3)
gp1 = gp1 + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
  bracket.nudge.y = -1, hide.ns = TRUE,
   label = "p.adj.signif") 
# gp2 = ggplot (ccomp, aes (x = merged_call, y = G2M.Score)) + geom_boxplot()
# gp3 = ggplot (ccomp, aes (x = merged_call, y = S.Score)) + geom_boxplot()
pdf (file.path ('Plots','cellcycle_boxplots.pdf'), height=3, width = 4)
gp1
dev.off()

### Make heatmap of all KDs against reference for all cNMFs
# Load libraries

library(tidyverse)

df = as.data.frame (srt@meta.data) 
df$merged_call[is.na (df$merged_call)] = 'NT'
df = df[!is.na(df$merged_call),]
#top_kds = names(head(table (df$crispr_calls)[order(-table(df$crispr_calls))],50))
#df = df[df$crispr_calls %in% top_kds,]

# Extract Cm columns
cm_cols <- grep("^Cm\\d+", colnames(df), value = TRUE)

# Get mean expression for NTC_1
ntc_df1 <- df[df$merged_call == "NTC", ]
#ntc_df3 <- df[df$merged_call == "NTC_3", ]
ntc_means1 <- colMeans(ntc_df1[, cm_cols], na.rm = TRUE)
ntc_means1 = apply (ntc_df1[, cm_cols], 2, median)
#ntc_means3 <- colMeans(ntc_df3[, cm_cols], na.rm = TRUE)
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
  #group_means <- colMeans(group_df[, cm_cols], na.rm = TRUE)
  group_means = apply (group_df[, cm_cols], 2, median)
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
pdf (file.path ('Plots','knock_downs_heatmap2.pdf'), height=5, width=5)
Heatmap (t(diff_matrix), row_names_gp=gpar(fontsize = 7), top_annotation=ha,
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


# # De novo marker discovery
# srt2 = srt
# top_kds = names(head(table (df$merged_call)[order(-table(df$merged_call))],30))
# srt = srt[,srt$merged_call %in% top_kds]
# enricher_universe = 'all'
# logfcThreshold = .25
# pvalAdjTrheshold = 0.01
# metaGroupName = 'merged_call'
# top_pathways = 10
# top_genes = 5
# force = F
# source (file.path(scrna_pipeline_dir,'DEG_standard.R'))


#DEG2

# kds = names(top_kds)
# kd_to_remove = c('NTC_1|NTC_3','NTC_1')
# kds = kds[!kds %in% kd_to_remove]

#kds = kds[30:length (kds)]
srt$celltype2 = 'celltype2'
top_kds = unique(srt$merged_call)
top_kds2 = top_kds[!top_kds %in% c('NTC')]
#srt$crispr_calls2[srt$crispr_calls2 %in% c('NTC')] = 'NTC'
for (kd in top_kds2)
  {
  force = T
  do.fgsea = TRUE
  rankby = 'LFC' # Ranking to input in fgsea can be 'LFC' or 'pval_signedLFC'
  rankby = 'pval_signedLFC'
  logfcThreshold = 0.25 # min logFC threshold to consider for testing. Should be set on 0
  pvalAdjTrheshold = 0.05
  topGenes = 20
  addGene=NULL # Add gene(s) to include in the DEG heatmap
  FeatureSets = list (all = NULL)  # Specify subset of features to run the test on
  metaGroupName1 = 'celltype2' # Specify metaGroup for celltype / clustering
  metaGroupName2 = 'merged_call' # Specify metaGroup of condition to compare
  deg2Ident = c('NTC', kd) #vector of 2 elements specifying the comparison between two groups from metaGroupName2
  top_pathways = Inf
  source (file.path(scrna_pipeline_dir,'DEG2.R'))
  }

rankby = 'pval_signedLFC'

gmt_annotations = c(
#'c2.cp.kegg.v7.1.symbol.gmt',
# 'c2.cp.reactome.v7.1.symbol.gmt',
# 'c5.bp.v7.1.symbol.gmt',
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
top_pathways = 10
fgseaResAll_dp = dotGSEA (fgseaL_ann, padj_threshold = pvalAdjTrheshold, 
    type = 'fgsea',top_pathways = top_pathways,
    cluster_rows=T,
    cluster_cols=T)
  
pdf (file.path('Plots',paste0('fGSEA_annotation_',ann,'_proliferation_cells_included_dotplot.pdf')), width=8, height=6)
print(fgseaResAll_dp)
dev.off()
}  


#### Run SCENIC ####
force = T
org = 'human'
motif_window = 'tss500bp'#'10kbp'
scenic_name = 'crispr_tfs'
genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=7000))
source (file.path(scrna_pipeline_dir, 'SCENIC.R'))

# Run SCENIC plots ####
srt$scenic_mods = 'crispr_tfs'
motif_window = 'tss500bp'#'10kbp'
genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=5000))
metaGroupNames = c('merged_call','scenic_mods','merged_call')
reductionName = 'umap'
source (file.path(scrna_pipeline_dir, 'SCENIC_plots.R'))

#auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
# rownames (auc_mtx) = auc_mtx[,1]
# auc_mtx = auc_mtx[,-1]
# colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

ccomp = srt@meta.data
#ccomp = cbind (ccomp, auc_mtx[rownames(ccomp),])

tpn = c('NFKB2','OSR1','PBX1','PBX3','SOX18','GLIS2','GATA6')
tpn = paste0('SCENIC_',tpn)
gpl = lapply (tpn, function(x) ggplot (ccomp, aes_string (x = 'merged_call', y = x)) + geom_boxplot())

pdf (file.path ('Plots','scenic_modules_boxplots.pdf'))
gpl
dev.off()






km = readRDS (file.path ('..','scatac_ArchR','TF_activity_modules.rds'))
regulon_TFs_in_modules = list(
  km1 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 1])],
  km2 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 2])]
  )
auc_mtx = auc_mtx[, colnames(auc_mtx) %in% regulon_TFs_in_modules$km2]

srt$mod_2 = unname(rowMeans (auc_mtx))






#saveRDS (srt, 'srt_filtered.rds')














# Load cnmfs ####
library(readxl)
# Bueno ####
top_bueno_genes = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/MPM_naive_study/reproduction2/PM_scRNA_atlas/data/bueno_molecular_subtype_deg.rds')
if (!all (names(top_bueno_genes)  %in% colnames(srt@meta.data))) {

srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = top_bueno_genes, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # 
        listName = 'Bueno', outdir = paste0(projdir,'Plots/'))
}

### FIGURE 2B - Plot boxplots of sarcomatoid cnmf ####
sarc_nmf = 'Cm17'
ccomp_df = srt@meta.data
#ccomp_df = ccomp_df[ccomp_df$celltype == 'Malignant',]
ccomp_df = aggregate (ccomp_df[,sarc_nmf,drop=F], by=as.list(ccomp_df[,'crispr_calls',drop=F]), 'mean')
rownames(ccomp_df) = ccomp_df[,1]
ccomp_df = ccomp_df[,-1, drop=F]
sample_order = rownames(ccomp_df)[order(-ccomp_df$Cm17)]

ccomp_df = srt@meta.data
#ccomp_df = ccomp_df[ccomp_df$celltype == 'Malignant',]
ccomp_df$crispr_calls = factor (ccomp_df$crispr_calls, levels = sample_order)

box = ggplot (ccomp_df, aes_string (x= 'crispr_calls', y= sarc_nmf)) +
  geom_violin (trim=TRUE, aes_string (fill = 'crispr_calls'),size=2,
    width=1,
    scale='width',
    linewidth = .2, alpha=0.7) +
  geom_boxplot (aes_string(fill = 'crispr_calls'),
    linewidth = .2,
    width=0.2,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.3, alpha=0.7
     ) +
  gtheme +
#  scale_fill_manual (values= palette_paired) +
  NoLegend()
  
pdf(file.path('Plots',paste0('FIGURE_2B_',sarc_nmf,'_signatures_boxplot.pdf')),width=6,2) #width = 10, height = 11,
print (box)
dev.off()


#### FIGURE 2A - Generate Neftel diagram using the four subtypes from bueno ####
max_sarc = pmax(srt$Sarcomatoid, srt$`Biphasic-S`)
max_epit = pmax(srt$Epithelioid, srt$`Biphasic-E`)
srt$y_axis <- log2(abs(max_sarc - max_epit) + 1)
srt$y_axis[max_epit > max_sarc] <- -1 * srt$y_axis[max_epit > max_sarc]

srt$x_axis = 0
srt$x_axis[srt$y_axis > 0] = log2(abs(srt$Sarcomatoid - srt$`Biphasic-S`) + 1)[srt$y_axis > 0]
srt$x_axis[srt$y_axis > 0 & srt$Sarcomatoid > srt$`Biphasic-S`] <- -1 * srt$x_axis[srt$y_axis > 0 & srt$Sarcomatoid > srt$`Biphasic-S`]

srt$x_axis[srt$y_axis < 0] = log2(abs(srt$Epithelioid - srt$`Biphasic-E`) + 1)[srt$y_axis < 0]
srt$x_axis[srt$y_axis < 0 & srt$`Biphasic-E` > srt$Epithelioid] <- -1 * srt$x_axis[srt$y_axis < 0 & srt$`Biphasic-E` > srt$Epithelioid]
srt$bueno_color = 0
srt$bueno_color[srt$x_axis < 0 & srt$y_axis < 0] = 'Biphasic-E'
srt$bueno_color[srt$x_axis > 0 & srt$y_axis > 0] = 'Biphasic-S'
srt$bueno_color[srt$x_axis < 0 & srt$y_axis > 0] = 'Sarcomatoid'
srt$bueno_color[srt$x_axis > 0 & srt$y_axis < 0] = 'Epithelioid'

p2l = lapply (unique(srt$crispr_calls), function(x) 
  {    
  df = srt@meta.data[srt$crispr_calls == x,]
  tot_cells = nrow(df)    
  ggplot(df, aes(x=x_axis, y=y_axis), color='white') +
  geom_point(alpha=1, shape=21, stroke=.25, size=1) +
  xlim (c(-2,2)) + ylim (c(-2,2)) +
  geom_vline(xintercept = 0,linetype = 'dashed', size=.1) +
  geom_hline(yintercept = 0,linetype = 'dashed', size=.1) +
  #scale_fill_manual(values = palette_bulk) +
  annotate("text", x = -1.1, y = 1.8, label = paste0("Sarco ",round(sum(df$x_axis < 0 & df$y_axis > 0) / tot_cells *100,1),'%') , size=3.5) +
  annotate("text", x = 1.2, y = 1.8, label = paste0("Bi-S ",round(sum(df$x_axis > 0 & df$y_axis > 0)/ tot_cells *100,1),'%'), size=3.5) +
  annotate("text", x = -1.2, y = -1.8, label = paste0("Bi-E ",round(sum(df$x_axis < 0 & df$y_axis < 0)/ tot_cells *100,1),'%'), size=3.5) +
  annotate("text", x = 1.2, y = -1.8, label = paste0("Epit ",round(sum(df$x_axis > 0 & df$y_axis < 0)/ tot_cells *100,1),'%'), size=3.5) + 
  xlab('') +
  ylab('') + 
  ggtitle (x) + 
  theme_void() + 
  NoLegend()
  })
pdf ('Plots/FIGURE_S2A_neftel_diagram_on_malignant_cells_per_sample.pdf',width = 12,height = 12)
print (wrap_plots (p2l, ncol=5))
dev.off()





#cnmf_module = 'cNMF27'
genes = head (cnmf_spectra_unique_comb_full[['Cm17']],10)
genes = c(genes,'RUNX2','RUNX1','SNAI2')
#genes = rownames (srt)[grep ('KRT', rownames(srt))]
#genes = genes[!grepl ('KRTAP', genes)]
ps = as.data.frame (AverageExpression (srt, features = genes, group.by = 'crispr_calls')[[1]])
ps = ps[rowSums(ps) > 0,]
col_names = c('NTC-1','NTC-3')
ps_cntrl = rowMeans (ps[,col_names])
ps = ps - ps_cntrl
ps =  -1 * ps
ps = ps[,!colnames(ps) %in% col_names]
#ps = ps[,col_names]

ps_scaled = t(scale(t(ps)))
ps_scaled = na.omit (ps_scaled)
pdf (file.path('Plots','marker_genes_unque_difference_heatmap.pdf'),height= 7, width = 5.3)
Heatmap (ps_scaled,
#column_labels = sapply (colnames(ps3), function(x) unlist(strsplit(x, '_'))[1]), 
#column_rot = 45,
  cluster_columns = T, 
 # col = palette_cnv_fun,
  #name = 'Po-Pr',
  #top_annotation = ha,
  #left_annotation = ha1,
  row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
  column_names_gp = gpar(fontsize = 8))
  #col = brewer.pal (10,'Greys'))
dev.off()


# Show lower Cm17 and higher Cm2 in SOX9 KD compared to control
srt$SNAI2 = srt@assays$RNA@layers$data[rownames(srt) == 'SNAI2',]
ccomp = srt@meta.data[, c('Cm17','Cm2','merged_call','SNAI2')]
ccomp = ccomp[ccomp$merged_call %in% c('NTC','SOX9'),]
#ccomp[ccomp$crispr_calls == 'NTC_3'] = 'NTC_1'
gp1 = ggplot (ccomp, aes (x = merged_call, y = SNAI2, fill = merged_call)) + 
 #geom_point (shape = 21, alpha=.5) +
  # geom_violin (trim=TRUE,size=2,
  #   width=1,
  #   scale='width',
  #   linewidth = .2, alpha=0.7) + 
  geom_boxplot (
    linewidth = .2,
    width=.3,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.7
     ) + 
  scale_fill_manual (values = c(SOX9_3 = 'brown', NTC_3 = 'grey60', NTC_1 = 'grey40')) + 
  gtheme
  
gp2 = ggplot (ccomp, aes (x = merged_call, y = Cm17, fill = merged_call)) + 
  #geom_point (shape = 21, alpha=.5) +
  # geom_violin (trim=TRUE,size=2,
  #   width=1,
  #   scale='width',
  #   linewidth = .2, alpha=0.7) + 
  geom_boxplot (
    linewidth = .2,
    width=.3,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.7
     ) + 
  scale_fill_manual (values = c(SOX9_3 = 'brown', NTC_3 = 'grey60', NTC_1 = 'grey40')) + 
  gtheme
  
  library (rstatix)
  stat.test = ccomp %>% 
  wilcox_test(reformulate ('merged_call', 'SNAI2')) %>%
          adjust_pvalue (method = "none") %>%
          add_significance ()
          stat.test = stat.test %>% add_xy_position (x = metaGroupNames[1], step.increase=0.01)
          box + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
          bracket.nudge.y = 0, hide.ns = F,
          label = "p.adj.signif") + NoLegend()
        
  #stat_ellipse(aes(color = crispr_calls), type = "norm", size = 1) +

pdf (file.path ('Plots','SOX9_controls_scatter.pdf'),width = 6,3)
DotPlot (srt, group.by = 'merged_call', features= c('SNAI2','AXL','SERPINE2','SOX6','CALB2','ITLN1'), scale=F)
wrap_plots (gp1, gp2)
dev.off()


### Check relationship between Cm17 and Cm16 ####
ccomp = srt@meta.data[, c('Cm19','Cm16','crispr_calls')]
ccomp = ccomp[ccomp$crispr_calls %in% c('NTC_3','NTC_1','SOX9_3'),]
ccomp = ccomp[ccomp$Cm19 < 0.4,]
#ccomp[ccomp$crispr_calls == 'NTC_3'] = 'NTC_1'
gp1 = ggplot (ccomp, aes (x = Cm19, y = Cm16)) + 
  geom_point(shape = 21, alpha=0.5, color = "black", aes(fill = crispr_calls)) +
  geom_smooth(aes(fill = NULL),  # force single regression line
              method = "lm", se = FALSE, color = "blue") +
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "pearson") + gtheme
  #scale_fill_manual (values = c(SOX9_3 = 'brown', NTC_3 = 'grey60', NTC_1 = 'grey40')) + 
  #gtheme

pdf (file.path ('Plots','checking_cc_is_anticorrelated_EMT.pdf'), width=8, height=5)
gp1
dev.off()










### Try regressing out cc ####

# Using limma
library(limma)
srt = readRDS ('srt.rds')

# Filter out cells not following same CNV pattern ####
library (infercnv)
icnv = readRDS ('../../_cellranger_raw_Filter_400_800_25/no_harmony/infercnv/all_samples_subsampled_Inf_ref_broad_meso/infercnv.results.obj.Rds') # output folder
lexpr.data = log2(icnv@expr.data) 
lexpr.data = lexpr.data[,-icnv@reference_grouped_cell_indices[[1]]]
avg_cnv = rowMeans (lexpr.data)
cnv_cor = cor (lexpr.data, avg_cnv)
cnv_cor_filt = rownames(cnv_cor)[cnv_cor[,1] > 0.5]
srt = srt[, colnames (srt) %in% cnv_cor_filt]

# Map crispr calls ####
crispr_calls = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/cropseq/raw/Sample1Cropseq/cellranger_output/ALTS06_sample1_0_v1/per_sample_outs/ALTS06_sample1_0/count/crispr_analysis/protospacer_calls_per_cell.csv')
crispr_ref = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/cropseq/raw/Sample1Cropseq/cellranger_output/ALTS06_sample1_0_v1/per_sample_outs/ALTS06_sample1_0/count/crispr_analysis/feature_reference.csv')
srt$crispr_calls = crispr_calls$feature_call[match (colnames(srt), crispr_calls$cell_barcode)]

# Remove cells with no guides ####
srt = srt[,!is.na(srt$crispr_calls)]

res2 = as.data.frame (table (srt$crispr_calls))
res2 = res2[order(-res2$Freq),]
res2$Var1 = factor (res2$Var1, levels = unique (res2$Var1))
bp = ggplot (res2, aes (x = Var1, y = Freq)) + geom_bar (stat = 'identity') + gtheme

# Map only cells with one gRNA called or in combination with NTC
x <- srt$crispr_calls

keep <- sapply(strsplit(x, "\\|"), function(parts) {
  # Rule 1: single gene
  if (length(parts) == 1) return(TRUE)
  
  # Check if NTC is present
  ntc_present <- any(grepl("^NTC", parts))
  
  if (ntc_present) {
    # Remove NTC from list
    other_genes <- parts[!grepl("^NTC", parts)]
    # After removing NTC, get unique base gene names
    gene_names <- sub("_\\d+$", "", other_genes)
    # Keep only if there is exactly ONE unique non-NTC gene
    return(length(unique(gene_names)) == 1)
  }
  
  # Rule 3: all same gene name ignoring _number
  gene_names <- sub("_\\d+$", "", parts)
  if (length(unique(gene_names)) == 1) return(TRUE)
  
  return(FALSE)
})

srt <- srt[,keep]

res2 = as.data.frame (table (srt$crispr_calls))
res2 = res2[order(-res2$Freq),]
res2$Var1 = factor (res2$Var1, levels = unique (res2$Var1))
bp = ggplot (res2, aes (x = Var1, y = Freq)) + geom_bar (stat = 'identity') + gtheme

### Further remove guides with less than n cells ####
# Show genes left in after also filtering by lowest amount of cells 
cells_threshold = 9
gene_left = as.character(res2$Var1[res2$Freq > cells_threshold])
srt = srt[,srt$crispr_calls %in% gene_left]

# cc.genes <- readLines('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/gene_sets/regev_lab_cell_cycle_genes.txt')
# s.genes <- cc.genes[1:43]
# g2m.genes <- cc.genes[44:97]



# Merge guides
merged_guides = list (
  MEF2D = c('MEF2D-1','MEF2D-3','MEF2D−3|NTC−3','MEF2D−1|MEF2D−3','MEF2D−1|NTC−3','MEF2D−1|MEF2D−3|NTC−3','MEF2D−1|NTC−1'),
  MEF2A = c('MEF2A−3','MEF2A−1','MEF2A−3|NTC−3','MEF2A−1|NTC−3','MEF2A−3|NTC−1'),
  PITX1 = c('PITX1−2','PITX1−3','PITX1−1','PITX1−2|PITX1−3','PITX1−1|PITX1−2','PITX1−1|PITX1−3'),
  BPTF = c('BPTF−3','BPTF−2'),
  HMGA1 = 'HMGA1−3',
  TEAD2 = 'TEAD2-1',
  NTC = c('NTC-1','NTC-3'),
  TCF3 = c('TCF3-1','TCF3-3'),
  SOX9 = 'SOX9-3',
  TWIST1 = c('TWIST1-1','TWIST1-3'),
  TEAD4 = c('TEAD4-2','TEAD4-3'))
merged_guides = lapply (merged_guides, function(x) gsub("[−–-]", "_", x))
srt$merged_call = srt$crispr_calls

for (i in names(merged_guides)) srt$merged_call[srt$merged_call %in% merged_guides[[i]]] = i
table (srt$merged_call)  

srt = srt[, !grepl ('\\|', srt$merged_call)]




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

srt2 = srt

# Regress out cell cycle ####
expr_norm = srt2@assays$RNA@layers$data
expr_norm_agg = aggregate (t(expr_norm), by = list(crispr_calls = srt2$merged_call), mean)
rownames (expr_norm_agg) = expr_norm_agg[,1]
expr_norm_agg = expr_norm_agg[,-1]
expr_norm_agg = t(expr_norm_agg)
rownames (expr_norm_agg) = rownames (srt)

# Alternatively try from raw counts aggregated and logged ####
expr_norm_agg = log2(AverageExpression (srt2, group.by = 'merged_call')[[1]]+1)
colnames(expr_norm_agg) = gsub ('-','_', colnames(expr_norm_agg))
#expr_norm_agg = aggregate (t(expr_norm), by = list(crispr_calls = srt2$crispr_calls), mean)

srt = CreateSeuratObject (counts = expr_norm_agg, data = expr_norm_agg)
srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = lapply(cnmf_spectra_unique_comb_full, function(x) head(x,50)),
        cor_threshold = NULL, 
        pos_threshold = NULL, # 
        listName = 'Cm', outdir = paste0(projdir,'Plots/'))
#expr_norm_agg = srt@assays$RNA@layers$data
expr_adj = removeBatchEffect(expr_norm_agg, 
  #covariates = cbind(srt$cNMF8,srt$cNMF6,srt$cNMF1)
  #covariates = cbind(srt$S.Score, srt$G2M.Score)
  covariates = srt$Cm16)
rownames (expr_adj) = rownames (srt)

srt@assays$RNA@layers[['data']] = expr_adj
#source (file.path (scrna_pipeline_dir,'data_processing.R'))
# compute module scores ####
srt = ModScoreCor (
        seurat_obj = srt,
        geneset_list = lapply(cnmf_spectra_unique_comb_full, function(x) head(x,50)),
        cor_threshold = NULL, 
        pos_threshold = NULL, # 
        listName = 'Cm', outdir = paste0(projdir,'Plots/'))


# Extract Cm columns
df = as.data.frame (srt@meta.data) 
cm_cols <- grep("^Cm\\d+", colnames(srt@meta.data), value = TRUE)
# Get mean expression for NTC_1
ntc_df1 <- df["NTC", ]
#ntc_df3 <- df["NTC_3", ]
#ntc_means1 <- colMeans(ntc_df1[, cm_cols], na.rm = TRUE)
#ntc_means3 <- colMeans(ntc_df3[, cm_cols], na.rm = TRUE)
#ntc_means = rowMeans (cbind(ntc_means1, ntc_means3))
ntc_means = ntc_df1[, cm_cols]

df = df[,names(cnmf_spectra_unique_comb_full)]
#Cm_df = aggregate (Cm_df, by = list(crispr_calls = df$crispr_calls), FUN = mean)
#rownames (Cm_df) = Cm_df[,1]
#Cm_df = Cm_df[,-1]
#df = Cm_df
#rownames (df) = gsub ('-','_', rownames(df))
#Cm_df = as.data.frame (t(Cm_df))
#Cm_df$ntc_means = ntc_means
#Cm_df = removeBatchEffect (t(Cm_df), covariates = Cm_df$Cm16)
#df = t(Cm_df)
# Get all unique crispr_calls (including NTC_1) 
all_calls = rownames(df)
# Initialize matrix
diff_matrix <- matrix(nrow = length(all_calls), ncol = length(cm_cols))
rownames(diff_matrix) <- all_calls
colnames(diff_matrix) <- cm_cols

# Compute difference: NTC_1 mean - crispr_call mean
for (call in all_calls) {
  group_df <- df[call, ]
  # group_means <- colMeans(group_df[, cm_cols], na.rm = TRUE)
  
  diff_matrix[call, ] <- unlist(ntc_means) - unlist(group_df)
}

# # Initialize p-value matrix
# pval_matrix <- matrix(NA, nrow = length(all_calls), ncol = length(cm_cols))
# rownames(pval_matrix) <- all_calls
# colnames(pval_matrix) <- cm_cols

# # Loop over calls and Cm columns
# for (call in all_calls) {
  
#   group_df <- df[df$crispr_calls == call, ]
  
#   for (cm in cm_cols) {
#     control_vals <- df[df$crispr_calls %in% c("NTC_1", "NTC_3"), cm]
#     group_vals   <- group_df[[cm]]
    
#     # Check if both groups have >1 unique values
#     if (length(unique(control_vals)) > 1 && length(unique(group_vals)) > 1) {
#       pval_matrix[call, cm] <- wilcox.test(control_vals, group_vals, exact = FALSE)$p.value
#     } else {
#       pval_matrix[call, cm] <- NA
#     }
#   }
# }
#diff_matrix[diff_matrix > 0.2] = 0.2
#diff_matrix[diff_matrix < -0.2] = -0.2
rownames(diff_matrix) = gsub ('-','_',rownames(diff_matrix))
ha = HeatmapAnnotation (cells = anno_barplot(as.vector(unname(table (srt2$crispr_calls)[rownames(diff_matrix)]))))
pdf (file.path ('Plots','prolif_regressed_ps_from_norm_knock_downs_heatmap2.pdf'), height=5, width=6)
Heatmap (t(diff_matrix), row_names_gp=gpar(fontsize = 7), top_annotation=ha,
  column_names_rot=45,
  border=T,
  column_split = ifelse (grepl ('NTC',rownames (diff_matrix)),'control','KOs'),
  column_names_gp = gpar (fontsize=7))#,
  # cell_fun = function (j, i, x, y, width, height, fill) 
  #           {
  #          if (t(pval_matrix)[i, j] < 0.001)
  #             {
  #              grid.text("***", x, y, just='center', vjust=.8,
  #               gp = gpar(fontsize = 5, col='black'))
  #             } else {
  #             if(t(pval_matrix)[i, j] < 0.01)
  #                 {
  #                 grid.text("**", x, y, just='center', vjust=.8,
  #                 gp = gpar(fontsize = 5, col='black'))   
  #                 } else {
  #                 if(t(pval_matrix)[i, j] < 0.05)
  #                   {
  #                   grid.text("*", x, y, just='center', vjust=.8,
  #                   gp = gpar(fontsize = 5, col='black'))         
  #                   }}}
  #     })
dev.off()


# Run pathway enrichment 
expr_adj
expr_adj_diff = lapply (colnames(expr_adj), function(x) expr_adj[,x]- expr_adj[,'NTC'])
expr_adj_diff = do.call (cbind ,expr_adj_diff)
colnames (expr_adj_diff) = colnames (expr_adj)
expr_adj_diff = expr_adj_diff[,colnames(expr_adj_diff) != 'NTC']


fgseaResAll = list()
fgseaRanks = list()
for (ann in gmt_annotations)  
  {
  gmt.file = paste0 ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/',org,'/',ann)
  pathways = gmtPathways (gmt.file)
  pathways = lapply (pathways, function(x) x[!is.na(x)])
  pathways = lapply (pathways, function(x) x[!x == 'NA'])
  pathways = lapply (pathways, function(x) x[!x == ''])
  message (paste('Compute enrichment per cluster using annotation:', ann))
  fgseaResCluster = list()
  for (i in colnames (expr_adj_diff))
    {
    message (paste ('fGSEA running cluster',i))
    deg2Cluster = expr_adj_diff[,i]
    #if (rankby == 'pval_signedLFC') fgsea_ranks = -log10 (deg2Cluster$p_val + 1e-300) * sign (deg2Cluster$avg_log2FC)
    #if (rankby == 'LFC') fgsea_ranks = deg2Cluster$avg_log2FC
    #fgsea_ranks = -log10 (deg2Cluster$p_val + 1e-300) * sign (deg2Cluster$avg_log2FC)
    
    #fgsea_ranks = setNames (fgsea_ranks, deg2Cluster$gene)    
    #fgsea_ranks = fgsea_ranks[fgsea_ranks != 0]
    #degCluster = degCluster[degCluster$p_val_adj < 01,]
    #fgsea_ranks = -log10(degCluster$p_val_adj) * sign (degCluster$avg_log2FC)
    #fgsea_ranks = (-log10(degCluster$p_val_adj) + 1e-10) + abs(degCluster$avg_log2FC) * sign(degCluster$avg_log2FC)
    #fgsea_ranks = degCluster$avg_log2FC
    fgseaRes = fgseaMultilevel (pathways, 
        deg2Cluster, 
        minSize=15, 
        maxSize=500, # changed this from 1500 to 1000 cause it generated an error
        BPPARAM = NULL)
    fgseaResCol = collapsePathways (fgseaRes, stats = fgsea_ranks, pathway = pathways)
    mainPathways = fgseaRes[fgseaRes$pathway %in% fgseaResCol$mainPathways]
    fgseaResCluster[[i]] = mainPathways
    fgseaRanks[[ann]][[i]] = fgsea_ranks
    }
  fgseaResAll[[ann]] = fgseaResCluster
  }
saveRDS (fgseaResAll, file.path(paste0('cc_regressed_norm_fGSEA_annotation_',paste(gmt_annotations,collapse='_'),'_rankby_',rankby,'.rds')))
saveRDS (fgseaRanks, file.path(paste0('cc_regressed_norm_fgsea_ranks_',paste(gmt_annotations,collapse='_'),'.rds')))

pvalAdjTrheshold = 1
fgseaResAll_dp = lapply (fgseaResAll, function(y) dotGSEA (y, padj_threshold = pvalAdjTrheshold, 
    type = 'fgsea',top_pathways = top_pathways,
    cluster_rows=T,
    cluster_cols=T)
  )
  lapply (seq_along(fgseaResAll_dp), function(x) {
  pdf (file.path('Plots',paste0('fGSEA_annotation_',names (fgseaResAll_dp)[x],'_dotplots.pdf')),width=8, 3 +length(unique(fgseaResAll_dp[[x]]$data$pathway))/7)
    print(fgseaResAll_dp[[x]])
    dev.off()
    })




# ### Try with DESeq2 modelling in the cell cycle ####
# expr = AverageExpression (srt2, group.by = 'crispr_calls', layer = 'counts')[[1]]
# expr <- matrix(as.integer(expr), nrow = nrow(expr), ncol = ncol(expr),
#                    dimnames = dimnames(expr))

# colnames (expr) = gsub ('-','_',colnames(expr))

# ccomp = aggregate (srt2@meta.data[,c('Cm16')], by = list (crispr_calls = srt2$crispr_calls), FUN = mean)
# ccomp = ccomp[match(colnames(expr), ccomp$crispr_calls),]
# rownames (ccomp) = ccomp$crispr_calls
# colnames (ccomp)[2] = 'cc'
# library (DESeq2)
# #cc = setNames (ccomp$x, ccomp$crispr_calls)
# dds <- DESeqDataSetFromMatrix (countData = expr,
#                               colData = ccomp,
#                               design = ~crispr_calls)

# dds <- DESeq(dds)
# res <- results(dds, contrast=c("crispr_calls","SOX9_3","NTC_1"))

# normalized_counts <- counts(dds, normalized=TRUE)
# collapsed_counts_norm = t(sapply(setdiff(unique(matched_gene_symbols),c(NA,'')),function(g) colSums(ldm$counts[which(matched_gene_symbols==g),])))


# # get normalised data
# load ('normalized_counts.Rd')
# colnames (normalized_counts)

# meta = ldm$sample_metadata
# meta[is.na(meta)] = '-'
# meta_f = meta[colnames (normalized_counts),]

# hcc = t(normalized_counts)
# hcc_f <- hcc[,apply (hcc, 2, var, na.rm=TRUE) != 0]
# hcc.pca <- prcomp (hcc_f, scale. = TRUE)

# # Demo Style
# #meta = ldm[[4]]
# #meta[is.na(meta)] = '-'
# meta_f$timepoint_response = paste0 (meta_f$timepoint, '_', meta_f$response)
# pca_p = lapply (colnames(meta_f)[c(3,4,5,6,8,9)], function(x) ggbiplot(hcc.pca, obs.scale = 1, var.scale = 1,
#          groups = meta_f[,x], ellipse = TRUE, circle = TRUE,var.axes = FALSE, filter.drivers = rownames(hcc.pca$rotation)[1],
#          labels = meta_f[,2]) +
#   scale_color_discrete(name = '') +
#   theme(legend.direction = 'horizontal', legend.position = 'top'))

# library (patchwork)
# png ('pca_bulk_groups.png', 3000,1000)
# wrap_plots (pca_p, ncol=5)
# dev.off()









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

lapply (cnmf_spectra_unique, function(x) head(x, 20))

srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = lapply(cnmf_spectra_unique, function(x) head(x,50)),
        cor_threshold = NULL, 
        pos_threshold = NULL, # 
        listName = 'Cm', outdir = paste0(projdir,'Plots/'))

#ccomp = srt@meta.data[,c(names (cnmf_spectra_unique),names (cnmf_spectra_unique_comb_full))]
ccomp = srt@meta.data[,names (cnmf_spectra_unique)]
ccomp_cor = cor (ccomp)
pdf (file.path ('Plots','cnmf_cor_heatmap.pdf'))
Heatmap (ccomp_cor)
dev.off()

ccomp = srt@meta.data
bp = lapply (names (cnmf_spectra_unique), function(x) 
  ggplot (ccomp, aes_string (x = 'merged_call', y = x)) + geom_boxplot () + gtheme)

pdf (file.path ('Plots','cNMF_boxplots.pdf'),20,20)
wrap_plots (bp)
dev.off()

sapply (unique (srt$merged_call), function(x) sapply (cnmf_spectra_unique, function(y) x %in% y))



# Compare with cnmf from patient samples
ov = ovmat (c(cnmf_spectra_unique_comb_full, cnmf_spectra_unique), compare_lists=list(names(
  cnmf_spectra_unique_comb_full), names (cnmf_spectra_unique)))

pdf (file.path ('Plots','cellcycle_overlap_cnmf_heatmap.pdf'))
ov
dev.off()



### Check distributions of Cms after cc correction
ccomp = srt@meta.data
ccomp2 = srt2@meta.data

bp = lapply (names (cnmf_spectra_unique_comb_full), function(x) 
  ggplot (ccomp, aes_string (x = 'crispr_calls', y = x)) + geom_boxplot () + gtheme)

bp2 = lapply (names (cnmf_spectra_unique_comb_full), function(x) 
  ggplot (ccomp2, aes_string (x = 'crispr_calls', y = x)) + geom_boxplot () + gtheme)

pdf (file.path ('Plots','Cm_boxplots.pdf'),20,20)
wrap_plots (bp)
wrap_plots (bp2)
dev.off()

ccomp2 = srt2@meta.data[,names (cnmf_spectra_unique_comb_full)]
pdf (file.path ('Plots','cnmf_cor_unadjusted_heatmap.pdf'))
Heatmap (cor (ccomp2))
dev.off()



#### get variable features to check how many cell cycle genes are found ####
cc_var = VariableFeatures (srt)
cc_genes = cc_var[(cc_var %in% cnmf_spectra_unique_comb_full[['Cm16']])]
norm_mat = srt@assays$RNA@layers$data
norm_mat = norm_mat[rownames(srt) %in% cc_genes,]
norm_mat = t(scale(t(norm_mat)))
rownames (norm_mat) = rownames(srt)[match (cc_genes, rownames(srt))]
ha = HeatmapAnnotation (ko = ifelse (srt$crispr_calls %in% c('NTC_1','NTC_3'),'control','KO'))
ht = Heatmap (norm_mat, top_annotation = ha, row_names_gp=gpar (fontsize=5))
pdf (file.path ('Plots','cc_genes_heatmap.pdf'), height=8)
ht
dev.off()




### Try scaling data and regress cell cycle #####

# Using limma
library(limma)
srt = readRDS ('srt.rds')

# Filter out cells not following same CNV pattern ####
library (infercnv)
icnv = readRDS ('../../_cellranger_raw_Filter_400_800_25/no_harmony/infercnv/all_samples_subsampled_Inf_ref_broad_meso/infercnv.results.obj.Rds') # output folder
lexpr.data = log2(icnv@expr.data) 
lexpr.data = lexpr.data[,-icnv@reference_grouped_cell_indices[[1]]]
avg_cnv = rowMeans (lexpr.data)
cnv_cor = cor (lexpr.data, avg_cnv)
cnv_cor_filt = rownames(cnv_cor)[cnv_cor[,1] > 0.5]
srt = srt[, colnames (srt) %in% cnv_cor_filt]

# Map crispr calls ####
crispr_calls = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/cropseq/raw/Sample1Cropseq/cellranger_output/ALTS06_sample1_0_v1/per_sample_outs/ALTS06_sample1_0/count/crispr_analysis/protospacer_calls_per_cell.csv')
crispr_ref = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/cropseq/raw/Sample1Cropseq/cellranger_output/ALTS06_sample1_0_v1/per_sample_outs/ALTS06_sample1_0/count/crispr_analysis/feature_reference.csv')
srt$crispr_calls = crispr_calls$feature_call[match (colnames(srt), crispr_calls$cell_barcode)]

# Remove cells with no guides ####
srt = srt[,!is.na(srt$crispr_calls)]

res2 = as.data.frame (table (srt$crispr_calls))
res2 = res2[order(-res2$Freq),]
res2$Var1 = factor (res2$Var1, levels = unique (res2$Var1))
bp = ggplot (res2, aes (x = Var1, y = Freq)) + geom_bar (stat = 'identity') + gtheme

# Map only cells with one gRNA called or in combination with NTC
x <- srt$crispr_calls

keep <- sapply(strsplit(x, "\\|"), function(parts) {
  # Rule 1: single gene
  if (length(parts) == 1) return(TRUE)
  
  # Check if NTC is present
  ntc_present <- any(grepl("^NTC", parts))
  
  if (ntc_present) {
    # Remove NTC from list
    other_genes <- parts[!grepl("^NTC", parts)]
    # After removing NTC, get unique base gene names
    gene_names <- sub("_\\d+$", "", other_genes)
    # Keep only if there is exactly ONE unique non-NTC gene
    return(length(unique(gene_names)) == 1)
  }
  
  # Rule 3: all same gene name ignoring _number
  gene_names <- sub("_\\d+$", "", parts)
  if (length(unique(gene_names)) == 1) return(TRUE)
  
  return(FALSE)
})

srt <- srt[,keep]

res2 = as.data.frame (table (srt$crispr_calls))
res2 = res2[order(-res2$Freq),]
res2$Var1 = factor (res2$Var1, levels = unique (res2$Var1))
bp = ggplot (res2, aes (x = Var1, y = Freq)) + geom_bar (stat = 'identity') + gtheme

### Further remove guides with less than n cells ####
# Show genes left in after also filtering by lowest amount of cells 
cells_threshold = 9
gene_left = as.character(res2$Var1[res2$Freq > cells_threshold])
srt = srt[,srt$crispr_calls %in% gene_left]

# cc.genes <- readLines('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/gene_sets/regev_lab_cell_cycle_genes.txt')
# s.genes <- cc.genes[1:43]
# g2m.genes <- cc.genes[44:97]

# Merge guides
merged_guides = list (
  MEF2D = c('MEF2D-1','MEF2D-3','MEF2D−3|NTC−3','MEF2D−1|MEF2D−3','MEF2D−1|NTC−3','MEF2D−1|MEF2D−3|NTC−3','MEF2D−1|NTC−1'),
  MEF2A = c('MEF2A−3','MEF2A−1','MEF2A−3|NTC−3','MEF2A−1|NTC−3','MEF2A−3|NTC−1'),
  PITX1 = c('PITX1−2','PITX1−3','PITX1−1','PITX1−2|PITX1−3','PITX1−1|PITX1−2','PITX1−1|PITX1−3'),
  BPTF = c('BPTF−3','BPTF−2'),
  HMGA1 = 'HMGA1−3',
  TEAD2 = 'TEAD2-1',
  NTC = c('NTC-1','NTC-3'),
  TCF3 = c('TCF3-1','TCF3-3'),
  SOX9 = 'SOX9-3',
  TWIST1 = c('TWIST1-1','TWIST1-3'),
  TEAD4 = c('TEAD4-2','TEAD4-3'))
merged_guides = lapply (merged_guides, function(x) gsub("[−–-]", "_", x))
srt$merged_call = srt$crispr_calls

for (i in names(merged_guides)) srt$merged_call[srt$merged_call %in% merged_guides[[i]]] = i
table (srt$merged_call)  

srt = srt[, !grepl ('\\|', srt$merged_call)]


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

srt2 = srt

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

srt = ScaleData (srt, features = VariableFeatures(FindVariableFeatures (srt, nfeat=10000)), vars.to.regress = 'Cm16')
#srt = ScaleData (srt, features = VariableFeatures(FindVariableFeatures (srt, nfeat=10000)), vars.to.regress = c('S.Score','G2M.Score'))
var_feat = VariableFeatures(FindVariableFeatures (srt, nfeat=10000))

mod_genes = lapply(cnmf_spectra_unique_comb_full, function(x) 
  {
  a = head(x,50)
  a = a [a %in% var_feat]
})

# srt = AddModuleScore (srt, features = mod_genes,
#   assay = 'RNA', slot = 'scale.data')

scaled_data = srt@assays$RNA@layers$scale.data
rownames (scaled_data) = VariableFeatures(FindVariableFeatures (srt, nfeat=10000))
df = lapply (mod_genes, function(x){
  df = colMeans(scaled_data[rownames(scaled_data) %in% x,])
})
df = do.call (cbind, df)
df = as.data.frame (df)
df$merged_call = srt$merged_call


# Extract Cm columns
cm_cols <- grep("^Cm\\d+", colnames(df), value = TRUE)

# Get mean expression for NTC_1
ntc_df1 <- df[df$merged_call == "NTC", ]
#ntc_df3 <- df[df$merged_call == "NTC_3", ]
ntc_means1 <- colMeans(ntc_df1[, cm_cols], na.rm = TRUE)
#ntc_means3 <- colMeans(ntc_df3[, cm_cols], na.rm = TRUE)
#ntc_means = rowMeans (cbind(ntc_means1, ntc_means3))
ntc_means = ntc_means1
#ntc_means = ntc_means3

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
pdf (file.path ('Plots','knock_downs_prolif_regressed_scaled_heatmap2.pdf'), height=5, width=4)
Heatmap (t(diff_matrix), row_names_gp=gpar(fontsize = 7), top_annotation=ha,
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





# srt = ModScoreCor (
#         seurat_obj = srt, 
#         geneset_list = lapply(cnmf_spectra_unique, function(x) head(x,50)),
#         cor_threshold = NULL, 
#         pos_threshold = NULL, # 
#         listName = 'cNMF', outdir = paste0(projdir,'Plots/'))


scaled_matrix <- matrix(nrow = length(all_calls), ncol = length(cm_cols))
rownames(scaled_matrix) <- all_calls
colnames(scaled_matrix) <- cm_cols

for (call in all_calls) {
  group_df = df[df$merged_call == call, ]
  group_means = colMeans(group_df[, cm_cols], na.rm = TRUE)
  
  scaled_matrix[call, ] = group_means
}

pdf (file.path ('Plots','knock_downs_prolif_regressed_scaled_nodiff_heatmap.pdf'), height=5, width=5)
Heatmap (t(scaled_matrix), row_names_gp=gpar(fontsize = 7), top_annotation=ha,
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





### Pathway enrichments on scaled data with regressed cell cycle ####
#kds = kds[30:length (kds)]
library (presto)
ccomp = srt@assays$RNA@layers$scale.data
rownames (ccomp) = VariableFeatures(FindVariableFeatures(srt, nfeat= 10000))
kds = unique (srt$merged_call)
kds = kds[kds != 'NTC']
wlc_res = list()
for (kd in kds)
  {
  wlc_res[[kd]] = wilcoxauc (ccomp, srt$merged_call, c('NTC',kd))
  }
head (wlc_res[[1]])
names (wlc_res) = kds

fgseaResAll = list()
fgseaRanks = list()
for (ann in gmt_annotations)  
  {
  gmt.file = paste0 ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/',org,'/',ann)
  pathways = gmtPathways (gmt.file)
  pathways = lapply (pathways, function(x) x[!is.na(x)])
  pathways = lapply (pathways, function(x) x[!x == 'NA'])
  pathways = lapply (pathways, function(x) x[!x == ''])
  message (paste('Compute enrichment per cluster using annotation:', ann))
  fgseaResCluster = list()
  for (i in names (wlc_res))
    {
    message (paste ('fGSEA running cluster',i))
    deg2Cluster = setNames (wlc_res[[i]]$logFC, wlc_res[[i]]$feature)
    #if (rankby == 'pval_signedLFC') fgsea_ranks = -log10 (deg2Cluster$p_val + 1e-300) * sign (deg2Cluster$avg_log2FC)
    #if (rankby == 'LFC') fgsea_ranks = deg2Cluster$avg_log2FC
    #fgsea_ranks = -log10 (deg2Cluster$p_val + 1e-300) * sign (deg2Cluster$avg_log2FC)
    
    #fgsea_ranks = setNames (fgsea_ranks, deg2Cluster$gene)    
    #fgsea_ranks = fgsea_ranks[fgsea_ranks != 0]
    #degCluster = degCluster[degCluster$p_val_adj < 01,]
    #fgsea_ranks = -log10(degCluster$p_val_adj) * sign (degCluster$avg_log2FC)
    #fgsea_ranks = (-log10(degCluster$p_val_adj) + 1e-10) + abs(degCluster$avg_log2FC) * sign(degCluster$avg_log2FC)
    #fgsea_ranks = degCluster$avg_log2FC
    fgseaRes = fgseaMultilevel (pathways, 
        deg2Cluster, 
        minSize=15, 
        maxSize=500, # changed this from 1500 to 1000 cause it generated an error
        BPPARAM = NULL)
    fgseaResCol = collapsePathways (fgseaRes, stats = fgsea_ranks, pathway = pathways)
    mainPathways = fgseaRes[fgseaRes$pathway %in% fgseaResCol$mainPathways]
    fgseaResCluster[[i]] = mainPathways
    fgseaRanks[[ann]][[i]] = fgsea_ranks
    }
  fgseaResAll[[ann]] = fgseaResCluster
  }
saveRDS (fgseaResAll, file.path(paste0('cc_regressed_scaled_fGSEA_annotation_',paste(gmt_annotations,collapse='_'),'_rankby_',rankby,'.rds')))
saveRDS (fgseaRanks, file.path(paste0('cc_regressed_scaled_norm_fgsea_ranks_',paste(gmt_annotations,collapse='_'),'.rds')))

pvalAdjTrheshold = 0.05
fgseaResAll_dp = lapply (fgseaResAll, function(y) dotGSEA (y, padj_threshold = pvalAdjTrheshold, 
    type = 'fgsea',top_pathways = top_pathways,
    cluster_rows=T,
    cluster_cols=T)
  )
  lapply (seq_along(fgseaResAll_dp), function(x) {
  pdf (file.path('Plots',paste0('fGSEA_annotation_scaled_',names (fgseaResAll_dp)[x],'_dotplots.pdf')),width=8, 3 +length(unique(fgseaResAll_dp[[x]]$data$pathway))/7)
    print(fgseaResAll_dp[[x]])
    dev.off()
    })












#### Try removing cell cycle genes #####
srt = srt2
norm_mat = srt@assays$RNA@layers$data
cc_cor = cor (srt$Cm16, t(as.matrix(norm_mat)))
colnames (cc_cor) = rownames(srt)
cc_cor2 = as.data.frame (t(cc_cor))
cc_cor2 = cc_cor2[order(-cc_cor2$V1),, drop=F]

pdf (file.path ('Plots','cc_genes_cor_Cm16_distribution.pdf'))
hist (cc_cor2[,1], breaks=100)
dev.off()

cc_cor2 = na.omit(cc_cor2)
cc_genes = rownames(cc_cor2)[cc_cor2$V1 > 0.2]

expr_prolif = srt[,srt$cc > 0.2]@assays$RNA@layers$data
expr_prolif = rowMeans (expr_prolif)
expr_prolif = data.frame (UMI = expr_prolif, cc = 'cc', row.names = rownames(srt))
expr_prolif$cc_gene = rownames(expr_prolif) %in% cc_genes
expr_prolif$cor = cc_cor2$V1[match(rownames(expr_prolif), rownames(cc_cor2))]
expr_rest = srt[,srt$cc <  -0.25]@assays$RNA@layers$data
expr_rest = rowMeans (expr_rest)
expr_rest = data.frame (UMI = expr_rest, cc = 'rest', row.names = rownames(srt))
expr_rest$cc_gene = rownames(expr_rest) %in% cc_genes
expr_rest$cor = cc_cor2$V1[match(rownames(expr_rest), rownames(cc_cor2))]
expr_comb = rbind (expr_prolif, expr_rest)

gp = ggplot (expr_comb, aes (x = cc, y = UMI, fill = cc_gene)) + geom_boxplot()
sp = ggplot (expr_comb, aes (x = UMI, y = cor)) + geom_point (, alpha=0.5, size=1) + facet_wrap (~cc)
pdf (file.path ('Plots','cc_distribution_boxplots.pdf'), width=10, height=6)
gp
sp
dev.off()




















### CHeck expression of TF target genes inferred from scATAC-seq ####
target_genes = readRDS ('../../../tumor_compartment/scatac_ArchR/genes_with_tf_peaks.rds')
var_feat = VariableFeatures (srt)
target_genes = lapply (target_genes, function(x) x[x %in% var_feat])
guides = names (target_genes)
names (target_genes) = paste0(names (target_genes), '_target_genes')

srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = target_genes,
        cor_threshold = NULL, 
        pos_threshold = NULL, # 
        listName = 'tg', outdir = paste0(projdir,'Plots/'))

ccomp = srt@meta.data
bp = lapply (guides, function(x)
{
bp = ggplot (ccomp[ccomp$merged_call %in% c('NTC',x),], aes_string (
  x = 'merged_call', y = paste0(x,'_target_genes'))) + 
geom_boxplot () + gtheme
bp$data$merged_call = factor(bp$data$merged_call, levels = c('NTC',x))
bp
})


pdf (file.path ('Plots','target_gene_modules_dotplot.pdf'))
wrap_plots (bp)
dev.off()







ccomp = srt@meta.data
ccomp$SOX9 = srt@assays$RNA@layers$data[rownames(srt) == 'SOX9',]
bp = ggplot (ccomp, aes (x = merged_call, y = SOX9)) + geom_boxplot()
pdf(file.path('Plots','SOX9_expression.pdf'))
DotPlot (srt, group.by = 'merged_call', feature=c('SOX9','SNAI2','RUNX1','RUNX2','TWIST1'))
VlnPlot (srt, group.by = 'merged_call', feature = 'SOX9')
bp
dev.off()





gene = c('SOX9','SNAI2','RUNX2','SNAI3','NFYC','NFIC')

pdf (file.path ('Plots','SOX9_target_genes_dotplot.pdf'), height=4, width=6)
DotPlot (srt[,srt$merged_call %in% c('SOX9','NTC')], features = gene, group.by = 'merged_call') + gtheme
dev.off()

### Assess enrichment of Sox9 targets from Fuchs paper

srt$celltype_not_cc_filtered = 'celltype_not_cc_filtered'
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
  metaGroupName1 = 'celltype_not_cc_filtered' # Specify metaGroup for celltype / clustering
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
    minSize=15, 
    maxSize=500, # changed this from 1500 to 1000 cause it generated an error
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
           fgseaRanks) + 
           labs(title=z))
        })
  dev.off()
  pdf (file.path('Plots','SOX9_Fuchs_enrichment_plots.pdf'),5,3)
  print(ep)
  dev.off() 



pdf (file.path ('Plots','SOX9_target_genes_W2to6_dotplot.pdf'), height=4, width=206)
DotPlot (srt[,srt$merged_call %in% c('SOX9','NTC')], features = unique(pathways$W2to6), group.by = 'merged_call') + gtheme
dev.off()


