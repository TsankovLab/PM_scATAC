set.seed(1234)

####### ANALYSIS of stroma compartment #######
projdir = 'stroma'
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

addArchRThreads (threads = 8) 
addArchRGenome ("Hg38")

archp = loadArchRProject (projdir)
cell_subsets_order = c('Mesothelium','Fibroblasts','SmoothMuscle','Endothelial','LEC')

# Load RNA
srt = readRDS ('srt.rds')
#srt$celltype[srt$celltype == 'SmoothMuscle'] = 'Smooth Muscle'
#saveRDS (srt, '../scrna/srt.rds')
#sarc_order = read.csv ('../scrna/cnmf20_sarcomatoid_sample_order.csv', row.names=1)

# Load last istance

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
  archp$celltype_lv2 = 0
  archp$celltype_lv2[archp$Clusters_H %in% c('C2','C3','C4')] = 'Endothelial'
  archp$celltype_lv2[archp$Clusters_H %in% c('C7','C8','C9','C10')] = 'Fibroblasts'
  archp$celltype_lv2[archp$Clusters_H %in% c('C6')] = 'SmoothMuscle'
  archp$celltype_lv2[archp$Clusters_H %in% c('C1')] = 'LEC'
  archp$celltype_lv2[archp$Clusters_H_2 %in% c('C10')] = 'Mesothelium'
  archp$celltype_lv2[archp$Clusters_H %in% c('C5')] = 'Unknown'
  archp = archp[archp$celltype_lv2 != 'Unknown']

  pdf()
  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample",labelMeans =F,
     embedding = "UMAP_H", pal = palette_sample)
  umap_p5 = plotEmbedding (ArchRProj = archp, labelMeans =F,
    colorBy = "cellColData", name = "celltype_lv2",
    pal = palette_stroma,
     embedding = "UMAP_H")
  dev.off()

  pdf (file.path('Plots','celltype_lv2_umap2.pdf'),5,5)
  print (umap_p3)
  print (umap_p5)
  dev.off()
  

### Call peaks on celltypes ####
metaGroupName = 'Clusters_H'
metaGroupName = 'celltype_lv2'
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
  

  # Find DAM ####
  metaGroupName = "celltype_lv2"  
  force = T
  top_genes = Inf
  
  DAM_df = DAM (
  ArchRProj = archp,
  metaGroupName = metaGroupName,
  FDR_threshold = 1e-2,
  meandiff_threshold = 0,
  top_genes=top_genes,
  filter_by_scRNA=TRUE, # Make sure has same metaGroupName
  seurat_obj = srt,
  min_exp=.1,
  force = force)

# Save table for supplementary information
write.csv (DAM_df, paste0('DAM_table_',metaGroupName, '.csv'))

# Take only top five to show heatmap
DAM_df <- DAM_df %>%
  mutate(comparison = factor(comparison, levels = cell_subsets_order)) %>%
  group_by(comparison) %>%
#  arrange(desc(Log2FC), .by_group = TRUE) %>%
  slice_head(n = 5) %>%   # keep top 5 per celltype
  ungroup() %>%
  arrange(comparison)

if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames(mMat) = rowData(mSE)$name

mMat = mMat[unique (DAM_df$gene), ]

#mMat_mg = mMat[active_DAM, ]
mMat_mg = as.data.frame (t(mMat))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
mMat_mg = mMat_mg[cell_subsets_order,]

# Generate heatmap ####

 DAM_hm = Heatmap (t(scale(mMat_mg)), 
          row_labels = colnames (mMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          row_names_side = 'left',
          cluster_rows = F,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          name = 'chromVAR',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = rev(palette_deviation)#,
          #width = unit(2, "cm")
          #right_annotation = motif_ha
          )

#DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_heatmap.pdf')), width = 2.5, height = 4)
print(DAM_hm)
dev.off()


### Co-expression of TFs across cells #### 

# # Get deviation matrix ####
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name


# Subset Endothelial cells ####
metaGroupName = 'celltype_lv2'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% 'Endothelial'],
  outputDirectory = file.path('..','Endothelial','scatac_ArchR'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)
saveRDS (srt[,srt$celltype_lv2 == 'Endothelial'], file.path ('..','Endothelial','srt.rds'))

