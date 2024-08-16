conda activate meso_scatac
use UGER
R

####### ANALYSIS of TUMOR compartment #######
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/stroma/scatac_ArchR'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','PM_scATAC','utils','load_packages.R'))
source (file.path('..','..','PM_scATAC','utils','useful_functions.R'))
source (file.path('..','..','PM_scATAC','utils','ggplot_aestetics.R'))
source (file.path('..','..','PM_scATAC','utils','scATAC_functions.R'))
source (file.path('..','..','PM_scATAC','utils','palettes.R'))

# Load functions for hub detection ####
source (file.path('..','PM_scATAC','utils','knnGen.R'))
source (file.path('..','PM_scATAC','utils','addCoax.R'))
source (file.path('..','PM_scATAC','utils','Hubs_finder.R'))
source (file.path('..','PM_scATAC','utils','hubs_track.R'))

# Set # of threads and genome reference ####
addArchRThreads(threads = 8) 
addArchRGenome("hg38")
  
sample_names = c(
    # Tumor  
    'P1', # p786
    'P4', # p811
    'P8', # p826
    'P3', # p846
    'P5', #'p848'
    'P10', # p10
    'P11', # p11
    'P12', # p12
    'P13', # p13
    'P14',# p14
    # Normal
    'RPL_280_neg_1',
    'RPL_280_neg_2',
    'RPL_Epi_1',
    'RPL_Epi_2'#,
    #'cf_distal'
    )

# Load RNA
#srt = readRDS ('../tumor_compartment/scrna/scRNA_meso.rds')
#sarc_order = read.csv ('../scrna/cnmf20_sarcomatoid_sample_order.csv', row.names=1)

# Load RNA
srt = readRDS ('../scrna/srt.rds')
#sarc_order = read.csv ('../scrna/cnmf20_sarcomatoid_sample_order.csv', row.names=1)

# Load last istance
if (!file.exists ('Save-ArchR-Project.rds'))
   {
    # Fix the fragment file of multiome sample by removing the header 
    #system ('zcat atac_fragments.tsv.gz | grep -v ^\# | bgzip > atac_fragments_fixed.tzv.gz')
    
    fragment_paths =c(
    '/ahg/regevdata/projects/lungCancerBueno/10x/191121/scATAC_Pt_mesothelioma_CD45_neg_cellranger_atac_v1.2/138_ATACseq_CD45_neg_Lung_ATAC/outs/fragments.tsv.gz',
    '/ahg/regevdata/projects/lungCancerBueno/10x/200128/scATAC_Pt811_mesothelioma_CD45pos_neg_cellranger_atac_v1.2/161_ATACseq_Pt811_mesothelioma_CD45pos_CD45neg_ATAC/outs/fragments.tsv.gz',
    '/ahg/regevdata/projects/lungCancerBueno/10x/200331/10X_Single_Cell_ATAC/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01218_Robby_AlexTsankov/202_ATAC_826CD45pos_826CD45neg/outs/fragments.tsv.gz',  
    '/ahg/regevdata/projects/lungCancerBueno/10x/200721/scATAC_Pt846/10X_Single_Cell_RNA/TD01729_AlexTsankov/846-MesoPool-CD45-pos-CD45-neg-nuclei/outs/fragments.tsv.gz',
    '/ahg/regevdata/projects/lungCancerBueno/10x/200911/scATAC_Pt848/10X_Single_Cell_ATAC/TD01814_AlexTsankov/848/outs/fragments.tsv.gz',
    '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_10/230510/P10_scATAC/cellranger_output/ALTS03_Zhao6ATAC_0_v1/fragments.tsv.gz',
    '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_11/230714/ZHAO8mesotheliomaATAC/cellranger_output/ALTS04_Zhao8ATAC_0_v1/fragments.tsv.gz',
    '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_12/230718/ZHAO9mesotheliomaATAC/cellranger_output/ALTS04_Zhao9ATAC_0_v1/fragments.tsv.gz',
    '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_13/231018/ZHAO12mesotheliomaATAC/cellranger_output/ALTS04_Zhao12ATAC_0_v1/fragments.tsv.gz',
    '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_14/240109/ZHAO13mesotheliomaATAC/cellranger_output/ALTS04_Zhao13ATAC_0_v1/fragments.tsv.gz',#,
    '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz',
    '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz',
    '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz',
    '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz'#,
    #"/ahg/regevdata/projects/ICA_Lung/10x/200116/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01396_AlexTsankov/aDcd45n/outs/fragments.tsv.gz"
    )
     
      #setwd (projdir)  
      ArrowFiles = createArrowFiles (inputFiles = fragment_paths,
      sampleNames = sample_names,
      minTSS = 4, #Dont set this too high because you can always increase later
      minFrags = 1000,
      maxFrags = Inf,
      addTileMat = TRUE,
      addGeneScoreMat = TRUE,
      force = FALSE,
      subThreading = T
      )
      
  archp = ArchRProject (
    ArrowFiles = ArrowFiles, 
    outputDirectory = projdir,
    copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  )
  
  ### Subset ArchR object only for cells retained in Signac analysis ####
  archp_main = loadArchRProject ('../../main/scatac_ArchR/')
  archp$celltype = 0
  archp$celltype[match (rownames(archp_main@cellColData), rownames(archp@cellColData))] = archp_main$celltype
  
  # Use Annotation from signac ####
  normal = readRDS ('../../per_sample_QC_signac/signac_normal.rds')
  normal_annotation = data.frame (barcode= colnames(normal), celltype = normal$predicted.id)
  normal_annotation$barcode = gsub (paste0(sample_names[11],'_'),paste0(sample_names[11],'#'),normal_annotation$barcode)
  normal_annotation$barcode = gsub (paste0(sample_names[12],'_'),paste0(sample_names[12],'#'),normal_annotation$barcode)
  normal_annotation$barcode = gsub (paste0(sample_names[13],'_'),paste0(sample_names[13],'#'),normal_annotation$barcode)
  normal_annotation$barcode = gsub (paste0(sample_names[14],'_'),paste0(sample_names[14],'#'),normal_annotation$barcode)
  normal_annotation = normal_annotation[normal_annotation$barcode %in% rownames(archp@cellColData),]
  archp$celltype[match (normal_annotation$barcode, rownames(archp@cellColData))] = normal_annotation$celltype
  
  archp = archp[archp$celltype != 0]
  archp = archp[archp$celltype != 'bad_quality']
  archp = archp[archp$celltype %in% names(table (archp$celltype)[table(archp$celltype) > 10])]
  archp$celltype[archp$celltype == 'Fibroblast'] = 'Fibroblasts'
  archp$celltype[archp$celltype == 'Myeloid'] = 'MonoMac'
  archp$celltype[archp$celltype == 'B.cells'] = 'B_cells'
  
  archp = archp[!archp$celltype %in% c('AT1','AT2','Alveolar','B_cells','Ciliated','Mast','MonoMac','NK','Plasma','Secretory','T.NK.cells','T_cells','cDCs','pDCs','Malignant')]

  archp$Sample2 = archp$Sample
  archp$Sample2[grep ('RPL',archp$Sample2)] = 'normal_pleura'

  ## and replace with Mo's annotation ##
  # normal_mo = read.csv ('../../tumor_compartment/all_tissues_ArchR/metadata_external_data/tsankov_refined_annotation.csv')



  
  #archp = archp[!is.na(archp$celltype)]

  # Dimensionality reduction and clustering
  varfeat = 25000
  LSI_method = 2
  archp = addIterativeLSI (ArchRProj = archp,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = FALSE, LSIMethod = LSI_method,
    varFeatures = varfeat)

  archp = addHarmony (
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample2", force=FALSE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony",
    name='Clusters_H',
    force = TRUE)

  archp = addClusters (input = archp, resolution = 3,
    reducedDims = "IterativeLSI", maxClusters = 100,
    force = TRUE)
  archp = addUMAP (ArchRProj = archp, 
    reducedDims = "IterativeLSI",
    force = TRUE)
  # archp = addTSNE (ArchRProj = archp, 
  #   reducedDims = "IterativeLSI",
  #   force = TRUE)
  
  #archp = saveArchRProject (archp)

  umap_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
   name = "celltype", embedding = "UMAP")
  umap_p2 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample2",
     embedding = "UMAP")
  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample2",
     embedding = "UMAP_H")
  umap_p4 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "celltype",
     embedding = "UMAP_H")
  umap_p5 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters_H",
     embedding = "UMAP_H")
  
  pdf (file.path('Plots','celltype_umap_signac_filtered.pdf'),5,5)
  print (umap_p1)
  print (umap_p2)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()
  
  
  archp = saveArchRProject (archp, load = T, dropCells=T)
  
  } else {
  archp = loadArchRProject (projdir)
  }


### Call peaks on celltypes ####
metaGroupName = 'Clusters_H'
force=TRUE
if(!all(file.exists(file.path('PeakCalls', unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds')))) | force) 
source (file.path('..','PM_scATAC','callPeaks.R'))

### chromVAR analysis ####
force=FALSE
if (!all(file.exists(file.path('Annotations',
  c('Motif-Matches-In-Peaks.rds',
    'Motif-Positions-In-Peaks.rds',
    'Motif-In-Peaks-Summary.rds')))))
source (file.path ('..','PM_scATAC','chromVAR.R'))
  
# Find activating and repressing TFs #### 
if (!file.exists ('TF_activators_genescore.rds')) 
  {
    source (file.path('..','PM_scATAC','activeTFs.R'))
  } else {
    corGSM_MM = readRDS ('TF_activators_genescore.rds') 
  }



### Co-expression of TFs #### 
metaGroupName = 'Sample2'
if (!any (ls() == 'mSE')) mSE = ArchR::getMatrixFromProject (archp, useMatrix = 'MotifMatrix', logFile=NULL)
mSE = mSE[, archp$cellNames]
all (colnames(mSE) == rownames(archp))

# # Get deviation matrix ####
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Subset only for positively correlated TF with genescore ####
positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0.1]
mMat = mMat[positive_TF,]

# mMat_agg = as.data.frame (t(mMat))
# mMat_agg$metaGroup = as.character (archp_meta[,metaGroupName])
# mMat_agg = aggregate (.~ metaGroup, mMat_agg, mean)
# rownames (mMat_agg) = mMat_agg[,1]
# mMat_agg = mMat_agg[,-1]
# mMat_agg = t(mMat_agg)
# rownames (mMat_agg) = active_TF

mMat_cor = cor (as.matrix(t(mMat)), method = 'pearson')
#d = as.dist (1-cor(as.matrix(t(mMat))))

#d = dist (mMat, method ='euclidean')
#hc1 <- hclust(d, method = "complete" ) # Hierarchical clustering using Complete Linkage

# row_filt = rowSums (mMat_cor) != 0
# tf_name = rownames(mMat_cor)[row_filt]
km = kmeans (mMat_cor, centers=10)
#km_ha = rowAnnotation (
#  km = anno_simple(as.character(km$cluster), width = unit(2, "mm"),
#    col = setNames (km_col, as.character(unique(km$cluster))), border=T))

# cor_mMat_hm = Heatmap (mMat_cor,
#         cluster_rows = T,
#         cluster_columns= T,
#         #left_annotation = km_ha,
#         row_split = km$cluster,
#         column_title = 'all_MPM',
#         column_gap = unit(.2, "mm"),
#         row_gap = unit(.2, "mm"),
#         name = 'r',
#         #column_km_repeats=20,
#         col = viridis::turbo(100),
#         #clustering_distance_rows = 'pearson', 
#         #clustering_distance_columns = 'pearson', 
#         column_split = km$cluster,
#         #col = pal_corr1,
#         row_names_gp = gpar(fontsize = 0),
#         column_names_gp = gpar(fontsize = 0)#,
#         #width = unit(5, "cm")
#         )
cor_mMat_hm = draw (Heatmap (mMat_cor,# row_km=15,
  #left_annotation = ha,
  #rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=palette_module_correlation_fun, 
  row_split = km$cluster,
  column_split = km$cluster,
  #row_km=2, 
  #column_km=2,
#  right_annotation = ha,
  border=F,
  row_names_gp = gpar(fontsize = 0),
  column_names_gp = gpar(fontsize = 0)))#,
  # ,
  # cell_fun = function(j, i, x, y, w, h, fill) {
  #       if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
  #           grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#        }}))

pdf (file.path ('Plots','TF_modules.pdf'), width = 4,height=3)
cor_mMat_hm
dev.off()

tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = do.call (cbind, tf_modules)
archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% paste0('mod_',unique(km$cluster))]
archp@cellColData = cbind (archp@cellColData, tf_modules)

archp = addImputeWeights (archp)
TF_p = lapply (paste0('mod_',unique(km$cluster)), function(x) plotEmbedding (
    ArchRProj = archp,
    colorBy = "cellColData",
    name = x, 
    pal = palette_deviation,
    #useSeqnames='z',
    imputeWeights = getImputeWeights(archp),
    embedding = "UMAP_H"))

pdf (file.path ('Plots','TF_modules_umap.pdf'), width = 30,height=14)
wrap_plots (TF_p, ncol=5)
dev.off()


# Try with ridge plots ####
library(ggridges)
library(ggplot2)
library(viridis)
#library(hrbrthemes)
tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = as.data.frame (do.call (cbind, tf_modules))
tf_modules$Sample = archp$Sample2
tf_modules$celltype = archp$celltype

# Plot
rp = lapply (paste0('mod_',unique(km$cluster)), function(x) 
  ggplot(tf_modules, aes_string(x = x, y = 'Sample', fill = '..x..')) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, alpha=.5) +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue-White Diverging", direction = -1) +
    facet_wrap (.~celltype) +
    theme_classic())
pdf (file.path ('Plots','TF_modules_ridge_plots.pdf'), width = 40,height=6)
wrap_plots (rp, ncol=5)
dev.off()


# Compare normal vs tumor CAF / smooth muscle and endothelial ####

# Make data.frame of deviation difference and expression difference between normal and tumors ####
celltypes = c('SmoothMuscle','Fibroblasts','Endothelial')
min_cells = 10
for (celltype in celltypes)
  {

metaGroupName = 'Sample2'
archp_meta = as.data.frame (archp@cellColData)
archp_meta = archp_meta[archp_meta$celltype == celltype,]
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat_agg = as.data.frame (t(mMat))
mMat_agg = mMat_agg[rownames(archp_meta), ]
mMat_agg$metaGroup = as.character (archp_meta[,metaGroupName])
mMat_agg = aggregate (.~ metaGroup, mMat_agg, mean)
rownames (mMat_agg) = mMat_agg[,1]
mMat_agg = mMat_agg[,-1]
mMat_agg = t(mMat_agg)
mMat_agg = mMat_agg[rownames(mMat_agg) %in% rownames(srt),]
metaGroupName = 'sampleID'
srt_celltype = srt[,srt$celltype == celltype]
DefaultAssay(srt_celltype) = 'RNA'
normal_sample = names(table (srt_celltype$sampleID[grep ('HU',srt_celltype$sampleID)])[table (srt_celltype$sampleID[grep ('HU',srt_celltype$sampleID)]) > min_cells])
sample_names_rna = c('P1','P14','P13','P3','P12','P5','P11','P4','P8','P14',normal_sample)
ps = log2(as.data.frame (AverageExpression (srt_celltype, features = rownames(mMat_agg), group.by = metaGroupName)[[1]]) +1)
ps = ps[, colnames(ps) %in% sample_names_rna]

TF_diff_rna = data.frame (
  tumor_dev = apply (mMat_agg[,unique(archp_meta$Sample2)[!unique(archp_meta$Sample2) == 'normal_pleura']], 1, mean),
  normal_dev = mMat_agg[,unique(archp_meta$Sample2)[unique(archp_meta$Sample2) == 'normal_pleura']],
  normal_rna = if (length(normal_sample) > 1) apply (ps[,normal_sample], 1, mean) else ps[,normal_sample],
  tumor_rna = apply (ps[, !colnames(ps) %in% normal_sample], 1, mean),
  genescore = corGSM_MM$cor[match(rownames(mMat_agg), corGSM_MM$GeneScoreMatrix_name)])
TF_diff_rna$dev_diff = TF_diff_rna$tumor_dev - TF_diff_rna$normal_dev
TF_diff_rna$rna_diff = TF_diff_rna$tumor_rna - TF_diff_rna$normal_rna


# Compute significance per sample vs normal in scRNA and dev ####
library (presto)
pval_threshold = 0.01
occurrence_threshold = 2

srt_celltype$sampleID2 = srt_celltype$sampleID
srt_celltype$sampleID2[srt_celltype$sampleID2 %in% normal_sample] = 'normal'
rna_comparisons = list(
  P1 = c('P1','normal'),
  P11 = c('P11','normal'),
  P12 = c('P12','normal'),
  P13 = c('P13','normal'),
  P14 = c('P14','normal'),
  P3 = c('P3','normal'),
  P4 = c('P4','normal'),
  P5 = c('P5','normal'),
  P8 = c('P8','normal'))

rna_res = lapply (rna_comparisons, function(x) 
  wilcoxauc (srt_celltype[rowData (mSE)$name,], group_by = 'sampleID2', groups_use = x))
rna_res_df = lapply (names (rna_comparisons), function(x) rna_res[[x]][rna_res[[x]]$group == x,'padj',drop=FALSE])
rna_res_df = do.call (cbind , rna_res_df)
rownames (rna_res_df) = rownames (srt_celltype[rowData (mSE)$name,])
colnames (rna_res_df) = names (rna_comparisons)

# Repeat using chromVAR deviations ####
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = mMat[,rownames(archp_meta)]
#all (colnames(mMat) == rownames(archp@cellColData))

dev_comparisons = lapply(paste(unique(archp_meta$Sample2))[paste(unique(archp_meta$Sample2)) != 'normal_pleura'], function(x) c(x, 'normal_pleura'))
names (dev_comparisons) = paste(unique(archp_meta$Sample2))[paste(unique(archp_meta$Sample2)) != 'normal_pleura']

dev_res = lapply (dev_comparisons, function(x) 
  wilcoxauc (mMat, y = archp_meta$Sample2, groups_use = x))
dev_res_df = lapply (names (dev_comparisons), function(x) dev_res[[x]][dev_res[[x]]$group == x,'padj',drop=FALSE])
dev_res_df = do.call (cbind , dev_res_df)
rownames (dev_res_df) = rowData (mSE)$name
colnames (dev_res_df) = names(dev_comparisons)
#occurence_filter = apply (dev_res_df, 1, function(x) sum (x < pval_threshold))
# dev_res_df_filtered = dev_res_df[occurence_filter > occurrence_threshold, ]
# dev_res_df_filtered = dev_res_df_filtered[rownames(dev_res_df_filtered) %in% tf_tumor_pos,]
# dev_selected_TF = rownames(dev_res_df_filtered)
# selected_TF = intersect (rna_selected_TF, dev_selected_TF)

rna_res_df_logic = rna_res_df < pval_threshold
dev_res_df_logic = dev_res_df < pval_threshold
rna_res_df_logic = rna_res_df_logic[unique (intersect (rownames(dev_res_df_logic), rownames(rna_res_df_logic))),]
dev_res_df_logic = dev_res_df_logic[unique (intersect (rownames(dev_res_df_logic), rownames(rna_res_df_logic))),]
common_samples = intersect (colnames(rna_res_df_logic), colnames(dev_res_df_logic))
comb_res_df = rna_res_df_logic[,common_samples] + dev_res_df_logic[,common_samples]
comb_res_df[comb_res_df == 1] = 0
comb_res_df[comb_res_df == 2] = 1
selected_TF = rownames(comb_res_df) [rowSums (comb_res_df) > occurrence_threshold] 
tf_tumor_pos = rownames(TF_diff_rna)[TF_diff_rna$dev_diff > 0 & TF_diff_rna$rna_diff > 0]
selected_TF = selected_TF[selected_TF %in% tf_tumor_pos]


# Order by mean logFC ####
res_df2 = lapply (names (rna_comparisons), function(x) rna_res[[x]][rna_res[[x]]$group == x,'logFC',drop=FALSE])
res_df2 = do.call (cbind , res_df2)
rownames (res_df2) = rownames(mMat_agg)
tf_order = rownames(res_df2)[order (-apply (res_df2, 1, mean))]
selected_TF_ordered = tf_order[tf_order %in% selected_TF]
corGSM_MM_filtered = as.data.frame (corGSM_MM [match (selected_TF_ordered, corGSM_MM$GeneScoreMatrix_name),])
head (corGSM_MM_filtered [order (corGSM_MM_filtered$cor),c('GeneScoreMatrix_name','cor') ],100)

# Export selected TFs ####
write.csv (selected_TF_ordered, 'Active_TFs.csv')

# Plot distribution of diff deviation and diff expression between tumor and normal ####
diff_line = 0
TF_diff_rna$label = ''
TF_diff_rna$label[match (selected_TF_ordered,rownames(TF_diff_rna))] = selected_TF_ordered
TF_diff_rna$color = TF_diff_rna$label != ''
TF_diff_rna$label_top = ''
TF_diff_rna$label_top[match (head (selected_TF_ordered,20),rownames(TF_diff_rna))] = head (selected_TF_ordered,20)
TF_diff_rna$alpha = 0.7
TF_diff_rna$alpha[match (selected_TF_ordered,rownames(TF_diff_rna))] = 1
tf_diff_p = ggplot (TF_diff_rna, aes (x= dev_diff, y = rna_diff,label = label_top)) + 
  geom_point(aes(fill=genescore, color=color, alpha = alpha), size = 2, shape = 21, stroke=0.3) + # Color points based on x value
  scale_color_manual(values = c('FALSE' = "grey", 'TRUE' = "red")) + # Customize colors
  scale_fill_gradient(low = "white", high = "black") +
  #scale_fill_manual(values = c('FALSE' = "grey", 'TRUE' = "red")) + # Customize colors
  geom_vline(xintercept = diff_line, linetype = "dashed", color = "grey44", linewidth=.3) + # Vertical dashed line
  geom_hline(yintercept = diff_line, linetype = "dashed", color = "grey44", linewidth=.3) + # Vertical dashed line
  gtheme_no_rot + # Use a minimal theme
   geom_text_repel(
    segment.size=.05,
    max.overlaps = 100,
#    point.padding = 0.2, 
    size=2#,
#   nudge_x = .25,
#    nudge_y = .2,
#    segment.curvature = -1e-20
    ) +
    xlab ('deviation difference') + 
    ylab ('RNA difference') + 
    xlim (c(-0.2, .2)) + 
    ylim (c(-0.6, .6))

pdf (paste0 ('Plots/Diff_normal_tumor_deviation_and_rna_',celltype,'_scatterplot2.pdf'),5,height = 4)
print (tf_diff_p)
dev.off()


motifs <- unique(TF_diff_rna$label_top)
motifs = motifs[motifs != '']
markerMotifs <- getFeatures (archp, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep ("z:", markerMotifs, value = TRUE)

archp = addImputeWeights (archp)
TF_p = plotEmbedding(
    ArchRProj = archp, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP_H",
    pal = palette_deviation,
    imputeWeights = getImputeWeights(archp)
)

pdf (file.path ('Plots',paste0('DEG_TF_',celltype,'_umap.pdf')), width = 30,height=14)
print (wrap_plots (TF_p, ncol=5))
dev.off()
}


