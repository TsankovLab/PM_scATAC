conda activate meso_scatac
#use UGER
R

set.seed(1234)

####### ANALYSIS of NKT compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR'
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
#source (file.path('..','..','git_repo','utils','scATAC_functions.R'))

addArchRThreads (threads = 4) 
addArchRGenome ("Hg38")

####### ANALYSIS of stroma compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/stroma/scatac_ArchR'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)
archp = loadArchRProject (projdir)


# Load RNA
srt = readRDS ('../scrna/srt.rds')
#srt$celltype[srt$celltype == 'SmoothMuscle'] = 'Smooth Muscle'
#saveRDS (srt, '../scrna/srt.rds')
#sarc_order = read.csv ('../scrna/cnmf20_sarcomatoid_sample_order.csv', row.names=1)

# Load last istance
if (!file.exists ('Save-ArchR-Project.rds'))
   {
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
    name = "Harmony",
    groupBy = "Sample", force=TRUE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony",
    name='Clusters_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony",
    name='Clusters_H_2', resolution = 2,
    force = FALSE)
  
# Run genescore DAG ####
metaGroupName = "Clusters_H"
force = TRUE
if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source (file.path('..','..','git_repo','utils','DAG.R'))

pdf()
  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample",labelMeans =F,
     embedding = "UMAP_H", pal = palette_sample)
  umap_p5 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters_H",
     embedding = "UMAP_H")
  umap_p6 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters_H_2",
   embedding = "UMAP_H")
dev.off()

pdf (file.path('Plots','clusters_umap.pdf'),5,5)
print (umap_p3)
print (umap_p5)
print (umap_p6)
dev.off()
  
  # Check LEC markers to reannotate cluster
  genes = c('TFF3','PECAM1','VWF','PLVAP','NOTCH4','COL1A1','COL1A2','ACTA2','MYL9','HP','WT1')
  archp = addImputeWeights (archp)

  pdf (file.path ('Plots',paste0('markers_LEC_umap.pdf')), width = 30,height=14)
  TF_p = plotEmbedding(
      ArchRProj = archp, 
      colorBy = "GeneScoreMatrix", 
      name = genes, 
      embedding = "UMAP_H",
      pal = palette_expression,
      imputeWeights = getImputeWeights(archp)
  )
  TF_p
  dev.off()

  pdf (file.path ('Plots',paste0('markers_LEC_umap.pdf')), width = 30,height=14)
  wrap_plots (TF_p)
  dev.off()

  #archp$celltype[archp$Clusters == 'C7'] = 'LEC'
  archp$celltype = 0
  archp$celltype[archp$Clusters_H %in% c('C2','C3','C4')] = 'Endothelial'
  archp$celltype[archp$Clusters_H %in% c('C7','C8','C9','C10')] = 'Fibroblasts'
  archp$celltype[archp$Clusters_H %in% c('C6')] = 'SmoothMuscle'
  archp$celltype[archp$Clusters_H %in% c('C1')] = 'LEC'
  archp$celltype[archp$Clusters_H_2 %in% c('C10')] = 'Mesothelium'
  archp$celltype[archp$Clusters_H %in% c('C5')] = 'Unknown'
  archp = archp[archp$celltype != 'Unknown']

  pdf()
  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample",labelMeans =F,
     embedding = "UMAP_H", pal = palette_sample)
  umap_p5 = plotEmbedding (ArchRProj = archp, labelMeans =F,
    colorBy = "cellColData", name = "celltype",
    pal = palette_stroma,
     embedding = "UMAP_H")
  dev.off()

  pdf (file.path('Plots','celltype_umap.pdf'),5,5)
  print (umap_p3)
  print (umap_p5)
  dev.off()
  
  archp = saveArchRProject (archp, load = T, dropCells=T)
  
  } else {
  archp = loadArchRProject (projdir)
  }


### Call peaks on celltypes ####
metaGroupName = 'Clusters_H'
force = FALSE
peak_reproducibility='2'
pdf() # This is necessary cause cairo throws error and stops the script
if(!all(file.exists(file.path('PeakCalls', unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds'))) | force)
source (file.path('..','..','git_repo','utils','callPeaks.R'))
dev.off()


### chromVAR analysis ####
force=FALSE
if (!all(file.exists(file.path('Annotations',
  c('Motif-Matches-In-Peaks.rds',
    'Motif-Positions-In-Peaks.rds',
    'Motif-In-Peaks-Summary.rds')))) | force)
source (file.path ('..','..','git_repo','utils','chromVAR.R'))
  
# # Find activating and repressing TFs #### 
# if (!file.exists ('TF_activators_genescore.rds')) 
#   {
#     source (file.path('..','..','git_repo','utils','activeTFs.R'))
#   } else {
#     corGSM_MM = readRDS ('TF_activators_genescore.rds') 
#   }

# Run genescore DAG ####
metaGroupName = "celltype"
force = TRUE
if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source (file.path('..','..','git_repo','utils','DAG.R'))

# Differential Accessed motifs ####
metaGroupName = "celltype"
force=FALSE
source (file.path('..','..','git_repo','utils','DAM.R'))

mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat_mg = mMat[active_DAM, ]
mMat_mg = as.data.frame (t(mMat_mg))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
mMat_mg = mMat_mg[names (DAM_list),]

# Generate RNA pseudobulk of matching cell types ####
metaGroupName = 'celltype'
#selected_TF = c(rownames(DAM_hm@matrix), 'NR4A3','NR4A2','NR4A1')
ps = log2(as.data.frame (AverageExpression (srt, features = active_DAM, group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
ps = ps[, colnames(DAM_hm@matrix)]
ps_tf = ps[active_DAM,]

  
DAM_hm = Heatmap (t(scale(mMat_mg)), 
          row_labels = colnames (mMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          name = 'chromVAR',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = rev(palette_deviation)#,
          #width = unit(2, "cm")
          #right_annotation = motif_ha
          )

scaled_ps = t(scale(t(ps_tf)))
scaled_ps[is.na(scaled_ps)] = 0
TF_exp_selected_hm = Heatmap (scaled_ps,
        #right_annotation=tf_mark,
        #column_split = column_split_rna,
        cluster_rows = F, #km = 4, 
        name = 'expression (scaled)',
        column_gap = unit(.5, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=F, 
        col = palette_expression,
        row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
        column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
        border=T,
        width = unit(2, "cm"))

TF_exp_selected_hm2 = Heatmap (ps_tf,
        #right_annotation=tf_mark,
        #column_split = column_split_rna,
        cluster_rows = F, #km = 4, 
        name = 'expression',
        column_gap = unit(.5, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=F, 
        col = palette_expression,
        row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
        column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
        border=T,
        width = unit(2, "cm"))

pdf (file.path ('Plots','DAM_with_rna_expression_heatmaps.pdf'), width = 3,height=4)
draw (DAM_hm )#+ TF_exp_selected_hm + TF_exp_selected_hm2)
dev.off()
   

### Co-expression of TFs across cells #### 

# # Get deviation matrix ####
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Subset only for expressed TFs ####
metaGroupName = 'celltype'
ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat), group.by = metaGroupName)[[1]]) +1)
min_exp = 0.5
ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
active_TFs = rownames(ps)[rowSums(ps) > 0]
#positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0]
metaGroupName = 'celltype'
mMat = mMat[active_TFs,]
mMat = as.data.frame (t(mMat))
mMat$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat = aggregate (.~ metaGroup, mMat, mean)
rownames (mMat) = mMat[,1]
mMat = mMat[,-1]

mMat_cor = cor (as.matrix(mMat), method = 'pearson')
km = kmeans (mMat_cor, centers=3)

pdf (file.path ('Plots','TF_modules_heatmap.pdf'), width = 4,height=3)
cor_mMat_hm = draw (Heatmap (mMat_cor,# row_km=15,
  #left_annotation = ha,
  #rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=palette_deviation_cor_fun, 
  row_split = km$cluster,
  column_split = km$cluster,
  #row_km=2, 
  #column_km=2,
#  right_annotation = ha,
  border=T,
#   ,
  row_names_gp = gpar(fontsize = 0), 
  column_names_gp = gpar(fontsize = 0)
# cell_fun = function(j, i, x, y, w, h, fill) {# THIS DOESNT WORK NEED TO USE LAYER_FUN
#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#         }}
  ))
  # ,
  # cell_fun = function(j, i, x, y, w, h, fill) {
  #       if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
  #           grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#        }}))
dev.off()

pdf (file.path ('Plots','TF_modules_heatmap.pdf'), width = 4,height=3)
cor_mMat_hm
dev.off()

if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = do.call (cbind, tf_modules)
archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% paste0('mod_',unique(km$cluster))]
archp@cellColData = cbind (archp@cellColData, tf_modules) 

pdf ()
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "cellColData",
    name = paste0('mod_',unique(km$cluster)), 
    pal = rev(palette_deviation),
    #useSeqnames='z',
    embedding = "UMAP_H")
dev.off()
pdf (file.path ('Plots','TF_modules_umap.pdf'), width = 20,height=6)
wrap_plots (TF_p, ncol=5)
dev.off()

tf_markers = c("WT1",'GATA4','GATA6','HIF1A','TGIF1')
markerMotifs = getFeatures (archp, select = paste(tf_markers, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs = grep ("z:", markerMotifs, value = TRUE)
archp = addImputeWeights (archp)
pdf ()
TF_p = plotEmbedding(
    ArchRProj = archp, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP_H",
    pal = rev (palette_deviation),
#    imputeWeights = getImputeWeights(archp),
    imputeWeights = NULL
)
dev.off()
pdf (file.path ('Plots','selected_TF_umap.pdf'), width = 20,height=6)
wrap_plots (TF_p, ncol=5)
dev.off()


#colnames(TF_hm@ht_list$chromVAR@matrix)[unlist(column_order(TF_hm)[c('2','3','4','5')])]
#which (colnames(mMat) == 'NR4A2')
#sapply (column_order(TF_hm), function(x) 497 %in% x)

# Distance matrix ####
d <- as.dist(1 - cor(t(mMat), method='pearson'))

# Hierarchical clustering ####
hc <- hclust(d)

# Dendrogram ####
pdf (file.path ('Plots',paste0('TF_',metaGroupName,'_no_km_dendrogram.pdf')), width=3, height=3.6)
plot(hc)
dev.off()





# Subset Endothelial cells ####
metaGroupName = 'celltype'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% 'Endothelial'],
  outputDirectory = file.path('..','..','Endothelial','scatac_ArchR'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

# Subset Fibroblasts cells ####
metaGroupName = 'celltype'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% 'Endothelial'],
  outputDirectory = file.path('..','..','Endothelial','scatac_ArchR'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

metaGroupName='celltype'
subsetArchRProject_light (ArchRProject = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('Fibroblasts')],
  projdir_new = file.path('..','..','Fibroblasts','scatac_ArchR')
  )



