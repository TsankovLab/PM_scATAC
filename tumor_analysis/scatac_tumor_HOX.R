conda activate meso_scatac
R

set.seed(1234)

####### ANALYSIS of TUMOR compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

# Load functions for hub detection ####
source (file.path('..','..','git_repo','utils','knnGen.R'))
source (file.path('..','..','git_repo','utils','addCoax.R'))
source (file.path('..','..','git_repo','utils','Hubs_finder.R'))
source (file.path('..','..','git_repo','utils','hubs_track.R'))

# Set # of threads and genome reference ####
addArchRThreads (threads = 8)
addArchRGenome ("hg38")

if (!file.exists ('Save-ArchR-Project.rds')) 
  { source (file.path('..','..','PM_scATAC','scatac_tumor_create_ArchRobj.R'))
  } else {
 archp = loadArchRProject (projdir)   
  }

archp$Sample3 = archp$Sample2
archp$Sample3[archp$Clusters == 'C1'] = 'P11_HOX'


# Load RNA ####
srt = readRDS (file.path('..','scrna','srt.rds'))
srt$sampleID2 = srt$sampleID
srt$sampleID2[srt$sampleID2 %in% c('HU37','HU62')] = 'normal_pleura'
srt_NN = srt[,!srt$sampleID2 %in% c("normal_pleura")]
sarc_order = read.csv (file.path('..','scrna','cnmf20_sarcomatoid_sample_order.csv'), row.names=1)
sarc_order = sarc_order[! sarc_order$sampleID2 %in% c('normal_pleura'),]
sarc_order = rbind (data.frame (sampleID = 'normal_pleura', x = -1),sarc_order)
#archp$Sample2 = factor (archp$Sample2, levels = sarc_order$sampleID)

# Subset ArchR project to remove normal sample ####
archp_NN = archp[archp$Sample2 != 'normal1']

archp_NN = addClusters (input = archp_NN, resolution = .6,
  reducedDims = "IterativeLSI", maxClusters = 100,
  force = TRUE)

# Plot UMAP of samples and clusters ####
metaGroupNames = c('Sample2','Clusters')
pdf()
p <- plotEmbedding(
    ArchRProj = archp_NN, 
    colorBy = "cellColData", 
    name = metaGroupNames[1], 
    embedding = "UMAP",
    pal = palette_sample,
    labelMeans =F
#    imputeWeights = getImputeWeights(archp)
)
dev.off()

pdf (file.path('Plots','samples_umap.pdf'))
p
dev.off()

pdf()
p <- plotEmbedding(
    ArchRProj = archp_NN, 
    colorBy = "cellColData", 
    name = metaGroupNames[2], 
    labelMeans=FALSE,
    embedding = "UMAP"#,
#    pal = palette_expression,
#    imputeWeights = getImputeWeights(archp)
)
dev.off()

pdf (file.path('Plots','clusters_umap.pdf'))
p
dev.off()


# ### Run TF correlation to identify TF modules across cancers #### 
# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat = as.matrix (mMat)#[selected_TF,])
# mMat_NN = mMat[,rownames(archp_NN@cellColData)]

# # Filter by RNA expression ####
# metaGroupName = 'sampleID'
# active_TFs = exp_genes (srt_NN, rownames(mMat), min_exp = 0.5, metaGroupName)
# mMat = mMat[active_TFs, ]

# sampleN=500
# set.seed(123)
# mMat_NN = lapply (unique(archp_NN$Sample2), function(sam) {
#   if (sampleN > sum(archp_NN$Sample2 == sam)) sampled_bc = rownames(archp_NN@cellColData)[archp_NN$Sample2 == sam] else
#     sampled_bc = sample (rownames(archp_NN@cellColData)[archp_NN$Sample2 == sam],sampleN)
#   mMat_NN[,sampled_bc]
# })
# mMat_NN = do.call (cbind, mMat_NN)
# mMat_cor = cor (as.matrix(t(scale(mMat_NN))), method = 'spearman')

# # Correlate module scores with TFs ####
# set.seed(123)
# centers=5
# km = kmeans (mMat_cor, centers=centers)
# km_df = as.data.frame (km$cluster)
# km_df = km_df[order (km_df[,1]),,drop=F]
# write.csv (km_df, 'kmeans_TF_modules.csv')
# pdf ()
# cor_mMat_hm = draw (Heatmap (mMat_cor,
#   clustering_distance_rows='euclidean' ,
#   clustering_distance_columns = 'euclidean', 
#   col=palette_deviation_cor_fun, 
#   row_split = km$cluster,
#   column_split = km$cluster,
#   border=T,
#   row_names_gp = gpar(fontsize = 0), 
#   column_names_gp = gpar(fontsize = 0)))
# dev.off()

# pdf (file.path ('Plots','TF_modules_heatmap2.pdf'), width = 4,height=3)
# cor_mMat_hm
# dev.off()

# Clusters_sample = paste0(archp_NN$Clusters, '_',archp_NN$Sample2)
# remove_low_clusters = !Clusters_sample %in% names(table (Clusters_sample)[table (Clusters_sample) < 10])
# Clusters_sample = Clusters_sample[remove_low_clusters]
# archp_NN = archp_NN[remove_low_clusters]

# mMat_NN = mMat[,rownames(archp_NN@cellColData)]
# tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat_NN[names(km$cluster[km$cluster == x]),]))
# names (tf_modules) = paste0('mod_',unique(km$cluster))
# tf_modules = do.call (cbind, tf_modules)
# archp_NN@cellColData = archp_NN@cellColData[!colnames(archp_NN@cellColData) %in% paste0('mod_',unique(km$cluster))]
# all (rownames(archp_NN@cellColData) == colnames(mMat_NN))
# archp_NN@cellColData = cbind (archp_NN@cellColData, tf_modules) 

# archp_NN = addImputeWeights (archp_NN)
# pdf()
# TF_p = plotEmbedding (
#     ArchRProj = archp_NN,
#     colorBy = "cellColData",
#     name = paste0('mod_',unique(km$cluster)), 
#     pal = rev(palette_deviation),
#     #useSeqnames='z',
#     embedding = "UMAP")
# dev.off()
# pdf (file.path ('Plots','TF_modules_umap.pdf'), width = 20,height=6)
# wrap_plots (TF_p, ncol=5)
# dev.off()

# Check individual TFs

tf_markers = c('HOXA13','HOXB13','HOXC13','HOXC10')
markerMotifs = getFeatures (archp, select = paste(tf_markers, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs = grep ("z:", markerMotifs, value = TRUE)
#archp = addImputeWeights (archp)
pdf()
TF_p = plotEmbedding(
    ArchRProj = archp_NN, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    pal = rev (palette_deviation),
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path ('Plots','HOX_genes_umap.pdf'), width = 30, height=18)
wrap_plots (TF_p, ncol=4)
dev.off()

# # ridge plots of TF modules ####
# tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat_NN[names(km$cluster[km$cluster == x]),]))
# names (tf_modules) = paste0('mod_',unique(km$cluster))
# tf_modules = as.data.frame (do.call (cbind, tf_modules))
# tf_modules$Clusters_Sample = Clusters_sample
# tf_modules$Sample = sapply (tf_modules$Clusters_Sample, function(x) unlist(strsplit(x, '_'))[2])
# palette_clusters_sample = setNames (palette_sample[tf_modules$Sample], tf_modules$Clusters_Sample)
# palette_clusters_sample = palette_clusters_sample[!duplicated(palette_clusters_sample)]
# tf_modules = gather (tf_modules, module, expression,1:centers)
# tf_modules$module = factor (tf_modules$module, levels = paste0('mod_',names (row_order (cor_mMat_hm))))

# dp = ggplot (tf_modules) +
#   geom_density(aes(x=expression,fill=Clusters_Sample),color='white',
#                       alpha = 0.6) +
#   facet_wrap (~module, nrow = 7, scales = 'free',strip.position = "left") +
#   scale_fill_manual (values = palette_clusters_sample) +
#   gtheme_no_rot

# pdf (file.path ('Plots','TF_modules_ridge_plots.pdf'), width = 5,height=8)
# dp
# dev.off()

# Alternatively run DAM per cluster / sample ####
archp_NN$Clusters_sample = paste0(archp_NN$Clusters, '_',archp_NN$Sample2)
Clusters_sample = paste0(archp_NN$Clusters, '_',archp_NN$Sample2)
remove_low_clusters = !Clusters_sample %in% names(table (Clusters_sample)[table (Clusters_sample) < 10])
Clusters_sample = Clusters_sample[remove_low_clusters]
archp_NN = archp_NN[remove_low_clusters]
# archp2 = archp # store whole archr object in secondary object to run script only for the subset
# archp = archp_NN

srt_NN$Clusters_sample = srt_NN$sampleID

# Find DAM ####
metaGroupName = "Clusters_sample"
top_genes = Inf
force=TRUE

  top_genes = Inf
  DAM_df = DAM (
  ArchRProj = archp_NN,
  metaGroupName = metaGroupName,
  FDR_threshold = 1e-2,
  meandiff_threshold = 0,
  top_genes=top_genes,
  filter_by_scRNA=TRUE, # Make sure has same metaGroupName
  seurat_obj = srt_NN,
  min_exp=.1,
  force = force)

# Save table for supplementary information
write.csv (DAM_df, paste0('DAM_table_',metaGroupName, '.csv'))

top_genes = 3

DAM_df <- DAM_df %>%
group_by (comparison) %>%
slice_head(n = top_genes) %>%   # keep top 5 per celltype
ungroup() %>%
  mutate(
    # extract second term after "_" e.g. "P12"
    second_term = str_extract(comparison, "(?<=_)P\\d+"),
    # extract numeric part (12)
    second_num  = as.numeric(str_extract(second_term, "\\d+"))
  ) %>%
  arrange(second_num) %>%
  select(-second_term, -second_num)


if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames(mMat) = rowData(mSE)$name
mMat = mMat[unique (DAM_df$gene), rownames(archp_NN@cellColData)]
#mMat_mg = mMat[active_DAM, ]
mMat_mg = as.data.frame (t(mMat))
mMat_mg$metaGroup = as.character (archp_NN@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
library(stringr)
order_p <- as.numeric(str_extract(rownames(mMat_mg), "(?<=_P)\\d+"))
mMat_mg <- mMat_mg[order(order_p), ]
#mMat_mg = mMat_mg[celltype_order,]

#tf_labels = c(head(DAM_df$gene[DAM_df$comparison=='C1_P11'],top_genes))
ha = rowAnnotation (foo = anno_mark(at = match(head(DAM_df$gene[DAM_df$comparison=='C1_P11'],top_genes),colnames(mMat_mg)), 
    labels = head(DAM_df$gene[DAM_df$comparison=='C1_P11'],top_genes), labels_gp = gpar(fontsize = 8, fontface = 'italic')))
ha2 = HeatmapAnnotation (sample = sapply (rownames(mMat_mg), function(x) unlist(strsplit(x,'_'))[2]), col= list(sample = palette_sample))

library(stringr)

# extract second term: "P11", "P23", etc.
split_raw <- sapply(rownames(mMat_mg), function(x) unlist(strsplit(x,"_"))[2])

# extract numeric part
p_nums <- as.numeric(str_extract(split_raw, "\\d+"))

# get unique P groups in the correct numeric order
unique_levels <- unique(split_raw[order(p_nums)])

# make factor with ordered unique levels
split_factor <- factor(split_raw, levels = unique_levels)


DAM_hm = Heatmap (t(scale(mMat_mg)), 
          top_annotation = ha2,
          #column_split = sapply (rownames(mMat_mg), function(x) unlist(strsplit(x,'_'))[2]),
          column_split = split_factor,
          right_annotation = ha,
          row_labels = colnames (mMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 0, fontface = 'italic'),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          name = 'TF activity',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_deviation_fun(scale(mMat_mg))

          #right_annotation = motif_ha
          )

  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_heatmaps.pdf')), width = 6.5, height = 3.5)
print(DAM_hm)
dev.off()

# Compare HOX+ vs HOX- P11 clusters using also scRNA ####
library (presto)

#sample_names_rna = c('P1','P14','P13','P3','P12','P5','P11','P4','P8','P14','HU37','HU62')
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat_P11 = mMat[,archp$Sample2 == 'P11']
P11_clusters = archp$Sample3[archp$Sample3 %in% c('P11','P11_HOX')]
p11_dev_rna = wilcoxauc (mMat_P11 , y = P11_clusters)
p11_dev_rna = p11_dev_rna[p11_dev_rna$group == 'P11_HOX',]
head (p11_dev_rna[order(p11_dev_rna$logFC),],50)
p11_dev_rna[p11_dev_rna$feature == 'NFATC1',]

# Get scRNA pseudobulks
ps = log2(as.data.frame (AverageExpression (srt, features = p11_dev_rna$feature, group.by = 'sampleID3')[[1]]) +1)
colnames(ps) = gsub ('-','_',colnames(ps))
#ps = ps[, colnames(ps) %in% sample_names_rna]
ps_P11_diff = ps[, 'P11_HOX', drop=F] - ps[, 'P11', drop=F]
ps_P11_diff['HOXB13',]
p11_dev_rna$rna_diff = NA
p11_dev_rna$rna_diff = ps_P11_diff[p11_dev_rna$feature,1]
p11_dev_rna$rna_diff[is.na(p11_dev_rna$rna_diff)] = 0

logfcThreshold = .1
pvalAdjTrheshold = 0.05
p11_dev_rna$sig = ifelse (abs(p11_dev_rna$logFC) > logfcThreshold & p11_dev_rna$padj < pvalAdjTrheshold, 1,0)
p11_dev_rna$sig = p11_dev_rna$sig * sign (p11_dev_rna$logFC)
p11_dev_rna$sig[p11_dev_rna$sig == 1] = 'HOX+'
p11_dev_rna$sig[p11_dev_rna$sig == -1] = 'HOX-'

p11_dev_rna$rna_sign = ifelse (abs(p11_dev_rna$logFC) > logfcThreshold & p11_dev_rna$padj < pvalAdjTrheshold, 1,0)
p11_dev_rna$rna_sign = p11_dev_rna$rna_sign * sign (p11_dev_rna$rna_diff)
p11_dev_rna$rna_sign[p11_dev_rna$rna_sign == 1] = 'HOX+'
p11_dev_rna$rna_sign[p11_dev_rna$rna_sign == -1] = 'HOX-'

res_filtered = p11_dev_rna[abs(p11_dev_rna$logFC) > logfcThreshold & p11_dev_rna$padj < pvalAdjTrheshold,]
res_filtered = head (res_filtered$feature[order (-abs(res_filtered$logFC))],20)
p11_dev_rna$labels = ''
p11_dev_rna$labels[match (res_filtered, p11_dev_rna$feature)] = res_filtered

vp = ggplot (p11_dev_rna, aes(x=logFC, y= -log10(padj))) +
    geom_point(shape=21, aes (fill = sig, color = rna_sign, size = abs(rna_diff)), alpha=.5) +
    geom_vline(xintercept = logfcThreshold, linetype="dashed", 
                color = "grey20", size=.4) +
    geom_vline(xintercept = -logfcThreshold, linetype="dashed", 
                color = "grey20", size=.4) +
    geom_hline(yintercept = -log10(pvalAdjTrheshold), linetype="dashed", 
                color = "grey20", size=.4) + 
    geom_text_repel (size=2.2, data = p11_dev_rna, aes(label = labels),segment.size=.2) + 
    ggtitle ('Hubs differential accessibility') +
    #geom_label_repel (size=2,max.overlaps=10000, data = deg2_cl, aes(label = show_genes), color='red') + 
    scale_color_manual (values = c("0"='grey77',"HOX-"='navyblue',"HOX+"='darkred')) + 
    scale_fill_manual (values = c("0"='grey77',"HOX-"='navyblue',"HOX+"='darkred')) + gtheme_no_rot


pdf (file.path ('Plots', 'P11_TF_volcano.pdf'),height=3.5,width=6.5)
vp
dev.off()


### Check HOXB13 expression in bulk and correlate with survival ####
# load bulkRNA cohorts ####
meso_bulk_l = readRDS (file.path ('..','..','bulkRNA_meso','bulk_RNA_studies.rds'))
meso_bulk_meta_l = readRDS (file.path ('..','..','bulkRNA_meso','bulk_RNA_studies_metadata.rds'))

# Correlate HOXB13 expression with TME proportions ####
bp_bueno = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/MPM_naive_study/reproduction2/bulkRNA/bayesprism_bueno_theta.rds')
bp_tcga = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/MPM_naive_study/reproduction2/bulkRNA/bayesprism_tcga_theta.rds')

bp_bueno = bp_bueno[colnames (meso_bulk_l[['bueno']]),]
cor (t(meso_bulk_l[['bueno']]['HOXB13',]),bp_bueno, method='spearman')

bp_tcga = bp_tcga[colnames (meso_bulk_l[['tcga']]),]
bp_tcga = na.omit (bp_tcga)
meso_bulk_tcga = meso_bulk_l[['tcga']][,rownames(bp_tcga)]
cor (meso_bulk_tcga['HOXB13',],bp_tcga, method='spearman')


### Check specificity of HOXB13 for tumor cells and genes in same module ####
srt = readRDS (file.path ('..','..','main','scrna','srt.rds')) # Replace srt object with whole sc data and subset for P11
srt = srt[,srt$sampleID == 'P11']
hoxb13_sig = read.csv (file.path('..','scatac_scrna_P11','HOXB13_cnmf.csv'))
hox_sig = head(hoxb13_sig[,2],Inf)

pdf (file.path('Plots','HOX_signature_dotplot.pdf'), width=30)
DotPlot (srt, feature = hox_sig) + gtheme
dev.off()

hox_sig_specific = c('HOXB13')#,'CA8','SULT1E1','SYT1')#,#'GAS2',#'WDR72','CHST9','RPRM',
#'LY6G6D')#,'COL9A1','MPPED2','TEKT3','CLIC5','NKX2−5','S100A7','ASPG','BEX1','GABRA2','PKP2','TDRD10','GRM1','PLPPR3',
  #'PI3','ACTR3B','LY6H','TNNT2','CXCL14')

exp_mat = log2(meso_bulk_l[[2]]+1)[rownames(meso_bulk_l[[2]]) %in% hox_sig_specific,]
ha = HeatmapAnnotation (subtype = meso_bulk_meta_l[[2]][,c(
  'subtype'#,
  #'CANCER_TYPE_DETAILED',
  #'AJCC_PATHOLOGIC_TUMOR_STAGE',
  #'AGE',
  #'ANEUPLOIDY_SCORE'
  )], col = list(subtype = palette_bulk))
#ha = HeatmapAnnotation (df = meso_bulk_meta_l[[2]][,grep ('STATUS',colnames(meso_bulk_meta_l[[2]]))])
hm = Heatmap (t(scale(t(exp_mat))), top_annotation = ha, border=T, column_names_gp=gpar(fontsize = 0),
  row_names_gp=gpar(fontface = 'italic'))

pdf (file.path('Plots','HOX_signature_tcga_heatmap.pdf'), width=9, height=6)
hm
dev.off()

exp_mat = log2(meso_bulk_l[[1]]+1)[rownames(meso_bulk_l[[1]]) %in% hox_sig_specific,]
#ha = HeatmapAnnotation (df = meso_bulk_meta_l[[1]][,c('subtype','PD.L1.expression..RPKM.','Asbestos.exposure','NF2..FISH.','FISH.chrom3','Type.of.pre.op.Treatment','Fibers.gm.lung','Sex')])
ha = HeatmapAnnotation (subtype = meso_bulk_meta_l[[1]][,c(
  'subtype'#,
  #'CANCER_TYPE_DETAILED',
  #'AJCC_PATHOLOGIC_TUMOR_STAGE',
  #'AGE',
  #'ANEUPLOIDY_SCORE'
  )], col = list(subtype = palette_bulk))
hm = Heatmap (t(scale(t(exp_mat))), top_annotation = ha, border=T, column_names_gp=gpar(fontsize = 0),
  row_names_gp=gpar(fontface = 'italic'))

pdf (file.path('Plots','HOX_signature_bueno_heatmap.pdf'), width=9, height=6)
hm
dev.off()

exp_mat = log2(meso_bulk_l[[3]]+1)[rownames(meso_bulk_l[[3]]) %in% hox_sig_specific,]
#ha = HeatmapAnnotation (df = meso_bulk_meta_l[[3]][,c('subtype','Subtype2','Cytological.Variant1','Stroma1','Chimio','Radioth','Immunoth','History.Location')])
ha = HeatmapAnnotation (subtype = meso_bulk_meta_l[[3]][,c(
  'subtype'#,
  #'CANCER_TYPE_DETAILED',
  #'AJCC_PATHOLOGIC_TUMOR_STAGE',
  #'AGE',
  #'ANEUPLOIDY_SCORE'
  )], col = list(subtype = palette_bulk))
hm = Heatmap (t(scale(t(exp_mat))), top_annotation = ha, border=T, column_names_gp=gpar(fontsize = 0),
  row_names_gp=gpar(fontface = 'italic'))

pdf (file.path('Plots','HOX_signature_mesomics_heatmap.pdf'), width=9, height=6)
hm
dev.off()



###########################
### Survival Analysis #####
###########################
library ('survminer')
library ('survival')

module_l = list(HOX_sig = hox_sig_specific)
meso_bulk_l_scaled = lapply (meso_bulk_l, function(x) t(scale(t(x))))
meso_bulk_meta_l2 = lapply (seq_along(meso_bulk_l_scaled), function(x) 
  {
   tmp = lapply (module_l, function(y) colMeans (meso_bulk_l_scaled[[x]][rownames(meso_bulk_l_scaled[[x]]) %in% y,,drop=F]))
   tmp = do.call (cbind, tmp)
   colnames (tmp) = gsub('-','_',colnames(tmp))
   tmp = tmp[,apply(tmp, 2, function(t) !any(is.na(t)))]
   meso_bulk_meta_l[[x]] = as.data.frame (meso_bulk_meta_l[[x]])
   meso_bulk_meta_l[[x]] = cbind (meso_bulk_meta_l[[x]], tmp)
   if(length(module_l) == 1) colnames(meso_bulk_meta_l[[x]])[ncol(meso_bulk_meta_l[[x]])] = names(module_l)
   meso_bulk_meta_l[[x]]
  })
names (meso_bulk_meta_l2) = c('bueno','tcga','mesomics')
#survival_name = 'activeTFs'
#survival_name = 'megahub'
survival_name = 'HOXs'



# Run Cox hazard ratio regression survival analysis ####
# Set variables per dataset to use
low='1st Qu.'
high='3rd Qu.'  
studies = c('bueno','tcga','mesomics')
cfit_study = list()
cox_l = list()
sc_p = list()
plot_study= list()
for (study in studies)
    {
    cfit = list()
    module_l = lapply (module_l, function(x) gsub('-','_',x))
    names (module_l) = gsub('-','_',names(module_l))
    mods = colnames(meso_bulk_meta_l2[[study]])[colnames(meso_bulk_meta_l2[[study]]) %in% names(module_l)]
    for (mod in mods)
        {
        meta_surv = meso_bulk_meta_l2[[study]]
        meta_surv = meta_surv[!is.na(meta_surv$census),]
  
        form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,"census"])),
                            status) ~', mod, '+ strata (subtype)'))
        # form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,"census"])),
        #                     status) ~', mod, '+ sarc_score'))
        cfit[[mod]] = coxph(form , data=meta_surv) 
        CI <- round(exp(confint(cfit[[mod]])), 2)
        cox_df = data.frame (
          HR = round(exp(coef(cfit[[mod]])), 2),
          CI = paste0('(',paste (CI, collapse=' - '),')'),
          LL = CI[1],
          UL = CI[2],
          P_value_C = round(summary(cfit[[mod]])$coefficients[, 5],2),
          label = mod
          )
        
        raw.vec=meta_surv[,mod]
        classified.vec=NA
        lowExpr = as.numeric(summary(raw.vec)[low])
        classified.vec[raw.vec < lowExpr]='Low'
        highExpr = as.numeric(summary(raw.vec)[high])
        classified.vec[raw.vec > highExpr]='High'
        classified.vec[is.na (classified.vec)] = 'Med'
        meta_surv[,mod] = factor (classified.vec, levels = c('Low','Med','High'))

        form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,"census"])),
                            status) ~', mod, '+ strata (subtype)'))
        # form = as.formula (paste('Surv(as.numeric(as.character(meta_surv[,"census"])),
        #                     status) ~', mod, '+ sarc_score'))
        cfit[[mod]] = coxph(form , data=meta_surv) 
        s = summary (cfit[[mod]])
        cox_df$P_value_S = round(s$logtest[3],2)  #s$logtest[3]
        cox_l[[mod]] = cox_df
        
        sc_p[[mod]] = ggadjustedcurves (cfit[[mod]], 
                data = meta_surv, 
                method = "conditional",
                variable=mod,
                palette = 'aaas',
                size=0.4,
                surv.median.line = 'hv',
                ggtheme = theme_classic()) +
                labs (title = mod,
                subtitle = paste0('log-rank = ',round(s$logtest[3],2)),
                caption = paste("n = ", nrow(meta_surv))) +
                      ylim(0,1)+ geom_hline(yintercept = 0.5,c(0.5,0.5),linetype='dotted', col='grey22') 
        }
    cox_df = do.call (rbind, cox_l)
    cox_df$Index = rownames (cox_df)
    cox_df$Index = factor (cox_df$Index, levels = rev(cox_df$Index))
    cfit_study[[study]] = cox_df
    plot_study[[study]] = sc_p
    }

for (study in names (cfit_study))
  {
  forest <- ggplot(cfit_study[[study]], aes(y = Index, x = HR)) + 
    geom_point(shape = 18, size = 5) +  
    geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
    #scale_y_continuous(name = "", breaks=1:nrow(cfit_study[[study]]), labels = cfit_study[[study]]$label, trans = "reverse", expand = expansion(add = 0.5)) +
    #scale_x_continuous(trans = 'log10')   + 
    xlab("Hazard Ratio") + 
    ylab(" ") + 
    theme_classic()
      
  tab <- ggplot(cfit_study[[study]], aes(y = Index)) +
    geom_text(aes(x = 0, label = sprintf("%0.1f", round(HR, digits = 2))), size = 4) +
    geom_text(aes(x = 1, label = CI), size = 4) + 
    geom_text(aes(x = 2, label = P_value_C), size = 4) + 
    geom_text(aes(x = 3, label = P_value_S), size = 4) + 
    #scale_y_continuous(trans = 'reverse', expand = expansion(add = 0.5)) +
    scale_x_continuous(
      breaks = 0:3, labels = c('HR', 'CI', 'P value (C)','P value (S)'), 
      position = 'top', expand = expansion(add = 0.5)) +
    theme_void() +
    theme(axis.text.x = element_text(face = 'bold'))
    
    pdf (file.path('Plots',paste0('cox_regression_',study,'_',survival_name,'.pdf')), height=8,5)
    print (forest + tab + plot_layout(ncol = 2, widths = c(1, 3)))
    dev.off()
    pdf (file.path('Plots',paste0('cox_regression_',study,'_stratified_',survival_name,'.pdf')), height = 2.8,3)
    print (plot_study[[study]])
    dev.off()
  }
  
  # Merge forest plots
  cfit_study_df = do.call (rbind, cfit_study)
  cfit_study_df$Index = rownames(cfit_study_df)
  cfit_study_df$Index = factor (cfit_study_df$Index, levels = rev(cfit_study_df$Index))

   forest <- ggplot(cfit_study_df, aes(y = Index, x = HR)) + 
    geom_point(shape = 18, size = 5) +  
    geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
    #scale_y_continuous(name = "", breaks=1:nrow(cfit_study[[study]]), labels = cfit_study[[study]]$label, trans = "reverse", expand = expansion(add = 0.5)) +
    #scale_x_continuous(trans = 'log10')   + 
    xlab("Hazard Ratio") + 
    ylab(" ") + 
    theme_classic()

  tab <- ggplot(cfit_study_df, aes(y = Index)) +
    geom_text(aes(x = 0, label = sprintf("%0.1f", round(HR, digits = 2))), size = 4) +
    geom_text(aes(x = 1, label = CI), size = 4) + 
    geom_text(aes(x = 2, label = P_value_C), size = 4) + 
    #geom_text(aes(x = 3, label = P_value_S), size = 3) + 
    #scale_y_continuous(trans = 'reverse', expand = expansion(add = 0.5)) +
    scale_x_continuous(
      breaks = 0:3, labels = c('HR', 'CI', 'P value (C)','P value (S)'), 
      position = 'top', expand = expansion(add = 0.3)) +
    theme_void() +
    theme(axis.text.x = element_text(face = 'bold'))
    
    pdf (file.path('Plots',paste0('cox_regression_',survival_name,'.pdf')), height=4,5)
    print (forest + tab + plot_layout(ncol = 2, widths = c(3, 3)))
    dev.off()
   






### Run chrombpnet on HOX cluster to identify peaks where it binds and with which TFs ####
pdf()
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = 'Clusters', 
    embedding = "UMAP"#,
#    pal = palette_expression,
#    imputeWeights = getImputeWeights(archp)
)
dev.off()

pdf (file.path('Plots','clusters_umap.pdf'))
p
dev.off()

if (!all(file.exists (file.path ('chromBPnet',c('fragments_C15.tsv','peakset_C15.bed')))))
source (file.path ('..','..','git_repo','tumor_analysis','chromBPnet_prepare_file_P11_HOX.R'))


### Combine TF modisco and finemo outputs to build network of co-occurring TFs across peaks ####
library (httr)
library (XML)
library (igraph)

chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet'
fold_number = 0
celltypes = c('C1','C2') 


count_to_adj = function(data = NULL)
  {
  count1s <- function(x, y) colSums(x == 1 & y == 1)
  n <- 1:ncol(data)
  mat <- outer(n, n, function(x, y) count1s(data[, x], data[, y]))
  diag(mat) <- 0
  dimnames(mat) <- list(colnames(data), colnames(data))
  return (mat)
  }
    
ov_motif_peaks_adj_l = list()
finemo_hits_l2 = list()
for (celltype in celltypes)
  {
  #celltype= 'IL1B'
  modisco_motifs = as.data.frame(readHTMLTable(file.path(chromBPdir, paste0(celltype,'_model'),'fold_0','modisco','report','motifs.html')))
  modisco_motifs$motif_match0 = sapply (modisco_motifs$NULL.match0, function(x) unlist(strsplit (x, '_'))[1])
  modisco_motifs$motif_match1 = sapply (modisco_motifs$NULL.match1, function(x) unlist(strsplit (x, '_'))[1])
  modisco_motifs$motif_match2 = sapply (modisco_motifs$NULL.match2, function(x) unlist(strsplit (x, '_'))[1])
  
  finemo_hits = read.table(file.path(chromBPdir,paste0(celltype,'_model'),'finemo_out','hits.tsv'), sep='\t', header=T)
  finemo_hits$motif_name0 = modisco_motifs$motif_match0[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
  finemo_hits$motif_name1 = modisco_motifs$motif_match1[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
  finemo_hits$motif_name2 = modisco_motifs$motif_match2[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
  finemo_hits$sample = celltype
  finemo_hits_l2[[celltype]] = finemo_hits
  write.table (finemo_hits[,c(1,2,3,14,15,16)], paste0(celltype,'_finemo_to_genome_browser.tsv'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  peakset = read.table (file.path(chromBPdir,paste0('peakset_',celltype,'.bed')))
  peakset = peakset[,1:3]
  colnames (peakset) = c ('chr','start','end')
  peakset = makeGRangesFromDataFrame (peakset)
  finemo_hits = makeGRangesFromDataFrame (finemo_hits, keep.extra.columns = T)
  
  match_rank = c('motif_name0')#,'motif_name1','motif_name2')
  ov_motif_peaks_mat_combined = list()
  for (motif_match in match_rank)
    {
    finemo_hits_l = split (finemo_hits, finemo_hits@elementMetadata[,motif_match])
    
    ov_motif_peaks = lapply (finemo_hits_l, function(x) findOverlaps (peakset, x, select='first'))
    ov_motif_peaks_mat = do.call (cbind, ov_motif_peaks)
    rownames (ov_motif_peaks_mat) = as.character(peakset)
    ov_motif_peaks_mat[is.na(ov_motif_peaks_mat)] = 0
    ov_motif_peaks_mat[ov_motif_peaks_mat > 0] = 1
    ov_motif_peaks_mat_combined[[motif_match]] = ov_motif_peaks_mat
    #ov_motif_peaks_df = as.data.frame (ov_motif_peaks_mat)
    }
  ov_motif_peaks_mat_combined = do.call (cbind, ov_motif_peaks_mat_combined)  
  tf_columns = colnames(ov_motif_peaks_mat_combined)
  ov_motif_peaks_adj_l2 = as.data.frame (count_to_adj (data = ov_motif_peaks_mat_combined))
  # From https://stackoverflow.com/questions/66515117/convert-dummy-coded-matrix-to-adjacency-matrix
  ov_motif_peaks_adj_l[[celltype]] = ov_motif_peaks_adj_l2 
  }

### Make barplot of top pioneer factors between C1 and C2 ####
finemo_hits_df = do.call (rbind, finemo_hits_l2)
top_finemo_hits = c(
  head(table (finemo_hits_df$motif_name0[finemo_hits_df$sample == 'C1'])[order(-table (finemo_hits_df$motif_name0[finemo_hits_df$sample == 'C1']))],8),
  head(table (finemo_hits_df$motif_name0[finemo_hits_df$sample == 'C2'])[order(-table (finemo_hits_df$motif_name0[finemo_hits_df$sample == 'C2']))],8))
finemo_hits_df_sub = finemo_hits_df[finemo_hits_df$motif_name0 %in% names(top_finemo_hits),]
palette_finemo_HOX = paletteer::paletteer_d("ochRe::olsen_seq")
bp = cellComp(
  seurat_obj = finemo_hits_df_sub,
  metaGroups = c('sample','motif_name0'),
  plot_as = 'bar',
  pal = palette_finemo_HOX
  )

pdf (file.path('Plots',paste0('chrombpnet_finemo_calls_',paste(celltypes, collapse='_'),'.pdf')))
bp
dev.off()

# Get partners in crime of HOXB13
celltype= 'C1'
pic = unlist(as.vector(ov_motif_peaks_adj_l[[celltype]]['HXB13',]))
pic[order(-pic)]

TFs = c(head(names(pic[order(-pic)]),5),'HOXB13')
TFs = c(TFs, 'ELK4')
pdf (file.path ('Plots','expression_of_TFs_binding_HOXB13_locus.pdf'),width=20)
DotPlot (srt, TFs, group.by='seurat_clusters') + gtheme
dev.off()


### Make barplot of top pioneer factors between C1 and C2 ####
finemo_hits_df = do.call (rbind, finemo_hits_l2)
peakset = lapply (celltypes, function(x) 
  {
  peakset = read.table (file.path(chromBPdir,paste0('peakset_',x,'.bed')))
  peakset = peakset[,1:3]
  colnames (peakset) = c ('chr','start','end')
  makeGRangesFromDataFrame (peakset)})
names(peakset) = celltypes
peaks_diff = peakset[['C1']][-queryHits(findOverlaps (peakset[['C1']], peakset[['C2']]))]
finemo_hits_gr = makeGRangesFromDataFrame(finemo_hits_l2[['C1']], keep.extra.columns=T)
finemo_hits_gr = finemo_hits_gr[finemo_hits_gr$motif_name0 == 'HXB13',]
hxb13_peaks = countOverlaps (peaks_diff, finemo_hits_gr)
hxb13_peaks = ifelse(hxb13_peaks >0,1,0)
table (hxb13_peaks)

hxb13_peaks = peaks_diff[unique(queryHits (findOverlaps (peaks_diff, finemo_hits_gr)))]
job = submitGreatJob (hxb13_peaks, species = 'hg38')
tbl = getEnrichmentTables(job)


# Try with heatmaps 
palette_cooccurence2 = paletteer::paletteer_c("grDevices::Oslo",100)
all_TF = unique(unname(unlist(lapply (ov_motif_peaks_adj_l, function(x) unique(colnames(x))))))
all_TF = all_TF[all_TF != 'NaN']
ht_list = NULL 
set.seed (1234)
hm_mat = as.data.frame (ov_motif_peaks_adj_l[[celltype]][match(all_TF, colnames(ov_motif_peaks_adj_l[[celltype]])),])
hm_mat = as.data.frame (t(hm_mat[all_TF,]))[all_TF,]
rownames (hm_mat) = all_TF
colnames (hm_mat) = all_TF
hm_mat[is.na(hm_mat)] = 0
#hm_mat = scale (hm_mat)
hm_mat = log2(hm_mat+1)
if (is.null(ht_list)) {
d = as.dist(t(1-cor(hm_mat)))
d[is.na(d)] = 1
hc = hclust(d)
column_order = colnames(hm_mat)[hc$order]
d = as.dist(1-cor(hm_mat))
d[is.na(d)] = 1
hc = hclust(as.dist(d))
row_order = colnames(hm_mat)[hc$order]
ht_list = Heatmap (
  hm_mat, 
  column_names_gp= gpar(fontsize=0),
  row_names_gp= gpar(fontsize=6), column_order=column_order, row_order=row_order,
  col = palette_cooccurence2, border=T, column_title = celltype)
} else {
ht_list = ht_list + Heatmap (
hm_mat[,column_order], cluster_columns=F,
column_names_gp= gpar(fontsize=0),
row_names_gp= gpar(fontsize=6),
col = palette_cooccurence2, border=T, column_title = celltype)
}

pdf (file.path ('Plots','combined_chromBPnet_cooccurrence_heatmaps.pdf'),width=4,height=3)
ht_list
dev.off()

  
# Map peaks with HOXB13 seqlet and run GSEA enrichment using nearby genes ####
TF = 'HXB13'
tf_peaks = finemo_hits[finemo_hits$motif_name0 == TF,]
peakHits = queryHits (findOverlaps(getPeakSet(archp),tf_peaks))
peakHits_top = head (table (peakHits)[order(-table(peakHits))],200)
tf_peaks = getPeakSet(archp)[as.numeric(names(peakHits_top))]
nearby_genes = unique(tf_peaks$nearestGene)
write.csv (nearby_genes, 'HOXB13_chrombpnet_top_hits.csv')
# use RNA to select for most expressed genes
metaGroupname = 'seurat_clusters'
ps = log2(as.data.frame (AverageExpression (srt, features = nearby_genes, group.by = metaGroupname)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
hox_cluster = 'g1'
ps = ps[, hox_cluster,drop=F]
ps = ps[order(-ps$g1),,drop=F]
head(ps,150)

#### Too many peaks have HOXB13!!


# Check expression of TF binding HOXB13 locus as predicted by chrombpnet and chromVAR ####

# Motif enrichment of peaks found around HOX13 genes ####
tf_match = getMatches (archp)
C1_peaks = readRDS (file.path ('PeakCalls','C1-reproduciblePeaks.gr.rds'))
tf_match = tf_match[queryHits (findOverlaps(tf_match, C1_peaks))]
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)
region = GRanges (c(
  'chr17:48726626-48729823'))

region_peaks = bg_peakSet[queryHits(findOverlaps(bg_peakSet, region))]
#tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
region_TF =  hyperMotif (
  selected_peaks = region_peaks, 
  motifmatch = tf_match)

head (region_TF, 100)

scrna_cor = read.csv (file.path ('..','scrna','correlated_genes_p11s_HOXB13.csv'))

cor_df = data.frame (rna_cor = scrna_cor$x[match(rownames(region_TF), scrna_cor$X)],
  TF_enrich_log10 = -log10(region_TF$padj), 
  TF_enrich_padj = region_TF$padj,
  row.names = rownames(region_TF)
  )
cor_df$label = ifelse (cor_df$TF_enrich_padj < 0.05, rownames(cor_df), '')
gp = ggplot (cor_df, aes (x = TF_enrich_log10, y = rna_cor, label = label)) + 
  geom_point(size = 2, shape = 21, stroke=0.3) + 
  geom_text_repel() + 
  gtheme_no_rot

pdf (file.path ('Plots','scrna_cor_vs_TF_enrich_HOX_scatter.pdf'),width=13,3)
gp
dev.off()


TFs = c('HOXB13','NFYB')
pdf (file.path ('Plots','expression_of_TFs_binding_HOXB13_locus.pdf'),width=20)
DotPlot (srt, TFs, group.by='seurat_clusters') + gtheme
dev.off()

# Check TF expression bindings HOBX13 locus without running hyper ####
tf_match = getMatches (archp)
region = GRanges (c(
  'chr17:48726626-48729823'))
region_peaks = tf_match[queryHits(findOverlaps(rowRanges(tf_match), region))]
TF_hits = colnames(assay(region_peaks))[colSums (assay(region_peaks)) > 0]
TF_hits = sapply (TF_hits, function(x) unlist(strsplit(x, '_'))[1])
pdf (file.path ('Plots','expression_of_TFs_binding_HOXB13_locus.pdf'),width=20)
DotPlot (srt, unname(TF_hits), group.by='seurat_clusters') + gtheme
dev.off()


















