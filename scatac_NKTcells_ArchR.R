conda activate meso_scatac
use UGER
R

set.seed(1234)

packages = c(
  'Signac',
  'Seurat',
  'biovizBase',
  'ggplot2',
  'patchwork',
  'scATACutils',
  'SummarizedExperiment',
  'epiAneufinder',
  'JASPAR2020',
  'TFBSTools',
  'TxDb.Hsapiens.UCSC.hg38.knownGene',
  'EnsDb.Hsapiens.v86',
  'gplots',
  'regioneR',
  'ComplexHeatmap',
  'ArchR',
  'BSgenome.Hsapiens.UCSC.hg38',
  'tidyverse',
  'ggrepel',
  'RColorBrewer')
lapply(packages, require, character.only = TRUE)

####### ANALYSIS of TUMOR compartment #######
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

#devtools::install_github("immunogenomics/presto") #needed for DAA
source ('../PM_scATAC/useful_functions.R')
source ('../PM_scATAC/ggplot_aestetics.R')
source ('../PM_scATAC/scATAC_functions.R')
source ('../PM_scATAC/palettes.R')

set.seed (1234)
addArchRThreads (threads = 8) 
addArchRGenome ("Hg38")

archp = loadArchRProject (projdir)

archp$status = ifelse (grepl ('^P',archp$Sample2), 'tumor','normal')  
# Dimensionality reduction and clustering ####
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


## harmonize project object ####
archp$project = ifelse (grepl ('^P', archp$Sample2), 'PM','normal')

  # Dimensionality reduction and clustering
  varfeat = 25000
  LSI_method = 2
  archp = addIterativeLSI (ArchRProj = archp,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = FALSE, LSIMethod = LSI_method,
    varFeatures = varfeat)

  archp = addClusters (input = archp, resolution = 3,
    reducedDims = "IterativeLSI", maxClusters = 100,
    force = TRUE)
  archp = addUMAP (ArchRProj = archp, 
    reducedDims = "IterativeLSI",
    force = TRUE)

  archp = addHarmony (
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    name = "Harmony_project",
    groupBy = c('project', 'Sample2'), force=TRUE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony_project", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony_project",
    name='Clusters_H',
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
  
  pdf (file.path('Plots','celltype_umap_signac_filtered_harmony_on_project.pdf'),5,5)
  print (umap_p1)
  print (umap_p2)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()

# Immune markers ####
meso_markers = c('FOXP3','CD8A','CD4','GNLY')
archp = addImputeWeights (archp)
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = meso_markers, 
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)

#archp$celltype[archp$Clusters == 'C30'] = 'Fibroblasts_WT1'
#p = lapply (p, function(x) x + theme_void() + NoLegend ()) #+ ggtitle scale_fill_gradient2 (rev (viridis::plasma(100))))

pdf (file.path('Plots','TNK_markers_fplots_cleaned.pdf'), width = 25, height = 25)
print (p)
dev.off()


run_peakCall = TRUE
if (run_peakCall)
  {
  ### Call peaks on celltypes ###
  metaGroupName = 'Clusters_H'
  archp = addGroupCoverages (
    ArchRProj = archp, 
    groupBy = metaGroupName,  
    force = FALSE,
    minCells= 20, # I think this should be set corresponding to the smallest cluster in the group or lower
    maxCells = 500,
    minReplicates = 2,
    sampleRatio = 0.8,
    useLabels = TRUE)
  
  archp = addReproduciblePeakSet (
      archp,
      groupBy= metaGroupName,
      peakMethod = 'Macs2',
      reproducibility = "1",
      maxPeaks = 500000, 
      minCells=20,
      force =TRUE) # I think this should be set corresponding to the smallest cluster in the group or lower
  archp = addPeakMatrix (archp)
  
  archp = saveArchRProject (archp, load=TRUE)
  
  metaGroupNames = c('TSSEnrichment','nFrags','ReadsInTSS','FRIP')  
    umap_p12 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
     name = x, embedding = "UMAP"))
      
  pdf (paste0(projdir,'/Plots/qc_umap_after_filtering.pdf'), 15,15)
  wrap_plots (umap_p12, ncol=5)
  dev.off()
  }

### chromVAR analysis ####
run_chromVAR = TRUE

if (run_chromVAR)
  {  
  archp = addBgdPeaks (archp, force= FALSE)
  archp = addMotifAnnotations (ArchRProj = archp,
      motifSet = "cisbp",
      #motifSet = 'JASPAR2020',
      #name = "JASPAR2020_Motif",
      force=TRUE)
  archp = addDeviationsMatrix (
    ArchRProj = archp, 
    peakAnnotation = "Motif",
    force = TRUE
  )}
  


# Find activating and repressing TFs ####
run_activeTF = FALSE

devMethod = 'ArchR'
 if (devMethod == 'ArchR')
    {
    TF_db='Motif'
    mSE = ArchR::getMatrixFromProject (archp, useMatrix = paste0(TF_db,'Matrix'))
    mSE = mSE[, archp$cellNames]
    rowData(mSE)$name = gsub ('_.*','',rowData(mSE)$name)
    rowData(mSE)$name = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", rowData(mSE)$name)
    }
  
if (!file.exists ('TF_activators_genescore.rds'))
  {
  seGroupMotif <- getGroupSE(ArchRProj = archp, useMatrix = "MotifMatrix", groupBy = "Clusters")
  seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
  rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
    rowMaxs(assay(seZ) - assay(seZ)[,x])
  }) %>% Reduce("cbind", .) %>% rowMaxs
  corGSM_MM <- correlateMatrices(
      ArchRProj = archp,
      useMatrix1 = "GeneScoreMatrix",
      useMatrix2 = "MotifMatrix",
      reducedDims = "IterativeLSI"
  )
  corGSM_MM = corGSM_MM[!grepl ('-AS',corGSM_MM$GeneScoreMatrix_name),]
  corGSM_MM = corGSM_MM[!grepl ('-DT',corGSM_MM$GeneScoreMatrix_name),]
  corGSM_MM = corGSM_MM[!grepl ('-OT',corGSM_MM$GeneScoreMatrix_name),]
  corGSM_MM = corGSM_MM[!grepl ('-RAB5IF',corGSM_MM$GeneScoreMatrix_name),]
  corGSM_MM = corGSM_MM[!grepl ('-IT2',corGSM_MM$GeneScoreMatrix_name),]
  corGSM_MM = corGSM_MM[!grepl ('-C8orf76',corGSM_MM$GeneScoreMatrix_name),]
  corGSM_MM = na.omit (corGSM_MM)
  saveRDS (corGSM_MM, 'TF_activators_genescore.rds')
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


table (archp$celltype)
archp$celltype2 = archp$celltype
archp$celltype2[archp$celltype2 == 'NaturalKillerTCell'] = 'NK'
archp$celltype2[archp$celltype2 == 'NaiveTcell'] = 'T_cells'
archp$celltype2[archp$celltype2 == 'TLymphocyte1(CD8+)'] = 'T_cells'
archp$celltype2[archp$celltype2 == 'Tlymphocyte2(CD4+)'] = 'T_cells'

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
tf_modules$celltype2 = archp$celltype2

# Plot
rp = lapply (paste0('mod_',unique(km$cluster)), function(x) 
  ggplot(tf_modules, aes_string(x = x, y = 'Sample', fill = '..x..')) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, alpha=.5) +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue-White Diverging", direction = -1) +
    facet_wrap (.~celltype2) +
    theme_classic())
pdf (file.path ('Plots','TF_modules_ridge_plots.pdf'), width = 30,height=8)
wrap_plots (rp, ncol=5)
dev.off()






