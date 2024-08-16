conda activate meso_scatac
use UGER
R

set.seed(1234)

####### ANALYSIS of Myeloid compartment #######
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','PM_scATAC','utils','load_packages.R'))
source (file.path('..','..','PM_scATAC','utils','useful_functions.R'))
source (file.path('..','..','PM_scATAC','utils','ggplot_aestetics.R'))
source (file.path('..','..','PM_scATAC','utils','scATAC_functions.R'))
source (file.path('..','..','PM_scATAC','utils','palettes.R'))

# Set # of threads and genome reference ####
addArchRThreads(threads = 8) 
addArchRGenome("hg38")

# Load ArchR project ####
archp = loadArchRProject (projdir)

# Add metadata ####
archp$status = ifelse (grepl ('^P',archp$Sample2), 'tumor','normal')  

## Subset only for tumor samples ####
archp = archp[archp$Sample2 %in% c('P1','P10','P11','P12','P13','P14','P3','P4','P5','P8')]

## Reduce dimension and harmonize ####
  varfeat = 25000
  LSI_method = 2
  archp = addIterativeLSI (ArchRProj = archp,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = TRUE, LSIMethod = LSI_method,
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
    groupBy = c('Sample2'), force=TRUE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony_project", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony_project",
    name='Clusters_H',
    force = TRUE)


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
  
  pdf (file.path('Plots','celltype_umap_harmony_sample_umap.pdf'),5,5)
  print (umap_p1)
  print (umap_p2)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()

# Remove doublets ####
archp = archp[archp$Clusters_H != 'C1']

# Dimensionality reduction and clustering
varfeat = 25000
LSI_method = 2
archp = addIterativeLSI (ArchRProj = archp,
  useMatrix = "TileMatrix", name = "IterativeLSI",
  force = TRUE, LSIMethod = LSI_method,
  varFeatures = varfeat)

archp = addClusters (input = archp, resolution = 5,
  reducedDims = "IterativeLSI", maxClusters = 100,
  force = TRUE)
archp = addUMAP (ArchRProj = archp, 
  reducedDims = "IterativeLSI",
  force = TRUE)

archp = addHarmony (
  ArchRProj = archp,
  reducedDims = "IterativeLSI",
  name = "Harmony_sample",
  groupBy = c('Sample2'), force=TRUE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony_sample", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp, resolution = 5,
    reducedDims = "Harmony_sample",
    name='Clusters_H',
    force = TRUE)


umap_p3 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample2",
   embedding = "UMAP_H")
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters_H",
   embedding = "UMAP_H")

pdf (file.path('Plots','celltype_umap_harmony_on_project_sample.pdf'),5,5)
print (umap_p3)
print (umap_p5)
dev.off()


# Check for doublets ####
meso_markers = c('C1QA','APOE','IL1B','CD3D')
archp = addImputeWeights (archp)
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = meso_markers, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)

pdf (file.path('Plots','myeloid_markers_fplots.pdf'), width = 18, height = 15)
wrap_plots (p, ncol=3)
dev.off()


# Remove doublets ####
archp = archp[archp$Clusters_H != 'C1']

# Dimensionality reduction and clustering
varfeat = 25000
LSI_method = 2
archp = addIterativeLSI (ArchRProj = archp,
  useMatrix = "TileMatrix", name = "IterativeLSI",
  force = TRUE, LSIMethod = LSI_method,
  varFeatures = varfeat)

archp = addClusters (input = archp, resolution = 10,
  reducedDims = "IterativeLSI", maxClusters = 100,
  force = TRUE)
archp = addUMAP (ArchRProj = archp, 
  reducedDims = "IterativeLSI",
  force = TRUE)

archp = addHarmony (
  ArchRProj = archp,
  reducedDims = "IterativeLSI",
  name = "Harmony_sample",
  groupBy = c('Sample2'), force=TRUE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony_sample", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp, resolution = 10,
    reducedDims = "Harmony_sample",
    name='Clusters_H',
    force = TRUE)


umap_p3 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample2",
   embedding = "UMAP_H")
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters_H",
   embedding = "UMAP_H")

pdf (file.path('Plots','harmony_sample_cleaned_umap.pdf'),5,5)
print (umap_p3)
print (umap_p5)
dev.off()


# Remove doublets ####
archp = archp[!archp$Clusters_H %in% c('C20','C23')]

# Dimensionality reduction and clustering
varfeat = 25000
LSI_method = 2
archp = addIterativeLSI (ArchRProj = archp,
  useMatrix = "TileMatrix", name = "IterativeLSI",
  force = TRUE, LSIMethod = LSI_method,
  varFeatures = varfeat)

archp = addHarmony (
  ArchRProj = archp,
  reducedDims = "IterativeLSI",
  name = "Harmony_sample",
  groupBy = c('Sample2'), force=TRUE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony_sample", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp, resolution = .3,
    reducedDims = "Harmony_sample",
    name='Clusters_H',
    force = TRUE)


umap_p3 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample2",
   embedding = "UMAP_H")
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters_H",
   embedding = "UMAP_H")

pdf (file.path('Plots','harmony_sample_cleaned2_umap.pdf'),5,5)
print (umap_p3)
print (umap_p5)
dev.off()


# Get markers for gene score ####
immune_markers = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/scRNA_immune_markers_humanLUAD_Samarth_Assaf.csv')
immune_markers = immune_markers [immune_markers$group %in% c('Neutrophil','TRMac','IM','DC2','DC1','pDC','mregDC','CD14 mono',
  'CD16 mono','NK','Mast cell','Mgk','B/Plasma',' T cell','Treg','MoMac'),]
#immune_markers = immune_markers[immune_markers$group %in% c('CD14 mono','CD16 mono','DC1','DC2','MoMac'),]
immune_markers = immune_markers$gene
immune_markers = immune_markers[!immune_markers %in% c('CD14 MONO','IHBA','SEPP1','IL3RA')]
archp = addImputeWeights (archp)
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = immune_markers, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)

png (file.path('Plots','myeloid_markers_fplots.png'), width = 18000, height = 15000, res=300)
wrap_plots (p)
dev.off()

### Cell annotation ####
archp$celltype = 0
archp$celltype = ifelse (archp$Clusters_H == 'C4','Monocytes','Macs')

archp$celltype2 = 0
archp$celltype2 = ifelse (archp$Clusters_H == 'C4','Monocytes','Macs')
archp$celltype2[archp$celltype2 == 'Macs'] = paste0('Macs',archp$Clusters_H[archp$celltype2 == 'Macs'])


umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
  pal = palette_sample,
   embedding = "UMAP_H")

umap_p6 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype",
   embedding = "UMAP_H")

umap_p4 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype2",
   embedding = "UMAP_H")

pdf (file.path('Plots','harmony_celltype_umap.pdf'),5,width = 10)
wrap_plots (umap_p6, umap_p5,umap_p4)
dev.off()








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


### Co-expression of TFs to identify TF modules #### 

# # Get deviation matrix ####
if (!exists ('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Subset only for positively correlated TF with genescore ####
positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0.1]
mMat = mMat[positive_TF,]

# Differential Accessed motifs ####
metaGroupName = "celltype2"
force=FALSE
source (file.path('..','PM_scATAC/DAM.R'))

# Correlate TF against each others and set number of kmeans clusters ####
mMat_cor = cor (as.matrix(t(mMat)), method = 'pearson')
km = kmeans (mMat_cor, centers=10)

# Generate heatmap ####
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

# Generate UMAPs of modules expression ####
tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),rownames(archp@cellColData)]))
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

pdf (file.path ('Plots','TF_modules_umap.pdf'), width = 30, height=14)
wrap_plots (TF_p, ncol=5)
dev.off()

# Generate ridge plots ####
library (ggridges)
library (ggplot2)
library (viridis)
#library(hrbrthemes)
tf_modules = lapply (unique (km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0 ('mod_',unique(km$cluster))
tf_modules = as.data.frame (do.call (cbind, tf_modules))
tf_modules$Sample = archp$Sample2
tf_modules$celltype = archp$celltype
tf_modules$celltype2 = archp$celltype2

# Plot
rp = lapply (paste0 ('mod_',unique(km$cluster)), function(x) 
  ggplot(tf_modules, aes_string(x = x, y = 'Sample', fill = '..x..')) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, alpha=.5) +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue-White Diverging", direction = -1) +
    facet_wrap (.~celltype2) +
    theme_classic())
pdf (file.path ('Plots','TF_modules_ridge_plots.pdf'), width = 30,height=8)
wrap_plots (rp, ncol=5)
dev.off()
