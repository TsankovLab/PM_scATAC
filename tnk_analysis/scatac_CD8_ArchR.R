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

####### ANALYSIS of NKT compartment #######
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/CD8'
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


## Reduce dimension and harmonize ####

  # Dimensionality reduction and clustering
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
    groupBy = c('project', 'Sample2'), force=TRUE
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
  
  pdf (file.path('Plots','celltype_harmony_sample_umap.pdf'),5,5)
  print (umap_p1)
  print (umap_p2)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()

# ## Further remove outlier clusters ####
archp = archp[!archp$Clusters_H %in% c('C1','C9')]
archp = archp[archp$Sample2 != 'P8'] # Remove sample with very few cells to avoid harmony weirdness 

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
    name = "Harmony_sample",
    groupBy = c('Sample2'), force=TRUE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony_sample", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony_sample",
    name='Clusters_H',
    force = TRUE)


  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample2",
     embedding = "UMAP_H")
  umap_p4 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "celltype",
     embedding = "UMAP_H")
  umap_p5 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters_H",
     embedding = "UMAP_H")
  
  pdf (file.path('Plots','celltype_harmony_sample_umap2.pdf'),5,5)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()

# ## Further remove outlier clusters ####
archp = archp[!archp$Clusters_H %in% c('C1','C2','C3','C4')]

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
    name = "Harmony_sample",
    groupBy = c('Sample2'), force=TRUE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony_sample", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony_sample",
    name='Clusters_H',
    force = TRUE)


  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample2",
     embedding = "UMAP_H")
  umap_p4 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "celltype",
     embedding = "UMAP_H")
  umap_p5 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters_H",
     embedding = "UMAP_H")
  
  pdf (file.path('Plots','celltype_harmony_sample_umap3.pdf'),5,5)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()



run_GS_analysis = TRUE

if (run_GS_analysis)
  {
  # Find DAG ####
  metaGroupName = "Clusters_H"
  force = FALSE
  if (!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force)
    {
    DAG_list = getMarkerFeatures (
      ArchRProj = archp, 
      testMethod = "wilcoxon",
            #useGroups = "ClusterA",
            #bgdGroups = "Clusters1B",
      binarize = FALSE,
      useMatrix = "GeneScoreMatrix",
      groupBy = metaGroupName
    #  useSeqnames="z"
    )

    listnames = colnames (DAG_list)
    DAG_list = lapply (1:ncol (DAG_list), function(x) 
      {
      df = DAG_list[,x]  
      df = do.call (cbind, (assays(df)))
      colnames(df) = names (assays(DAG_list))
      df$gene = rowData (DAG_list)$name
      df
      })
    names (DAG_list) = listnames
    saveRDS (DAG_list, paste0 ('DAG_',metaGroupName,'.rds'))    
    } else {
    DAG_list = readRDS (paste0('DAG_',metaGroupName,'.rds'))
    }
  
  FDR_threshold = 1e-8
  lfc_threshold = 1
  top_genes = 5
  DAG_top_list = DAG_list[sapply (DAG_list, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$Log2FC) > lfc_threshold,]) > 0)]
  DAG_top_list = lapply (seq_along(DAG_top_list), function(x) {
    res = DAG_top_list[[x]]
    res = na.omit (res)
    res = res[res$FDR < FDR_threshold,]
    res = res[order (res$FDR), ]
    res = res[abs(res$Log2FC) > lfc_threshold,]
    res$comparison = names(DAG_top_list)[x]
    if (nrow(res) < top_genes) 
      {
      res
      } else {
      head (res,top_genes)
      }
    })
  DAG_df = Reduce (rbind ,DAG_top_list)
  
  if (!any (ls() == 'gsSE')) gsSE = ArchR::getMatrixFromProject (archp, useMatrix = 'GeneScoreMatrix')
  gsSE = gsSE[, archp$cellNames]
  gsMat = assays (gsSE)[[1]]
  rownames (gsMat) = rowData (gsSE)$name
  gsMat_mg = gsMat[rownames (gsMat) %in% DAG_df$gene, ]
  gsMat_mg = as.data.frame (t(gsMat_mg))
  gsMat_mg$metaGroup = as.character(archp@cellColData[,metaGroupName])
  gsMat_mg = aggregate (.~ metaGroup, gsMat_mg, mean)
  rownames (gsMat_mg) = gsMat_mg[,1]
  gsMat_mg = gsMat_mg[,-1]
  gsMat_mg = gsMat_mg[names(table (archp@cellColData[,metaGroupName])[table (archp@cellColData[,metaGroupName]) > 50]),]
  DAG_hm = Heatmap (t(scale(gsMat_mg)), 
          row_labels = colnames (gsMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = T,
          #col = pals_heatmap[[5]],
          cluster_columns=T,#col = pals_heatmap[[1]],
          column_names_rot = 45,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_expression,
          name = 'GS'
          #right_annotation = motif_ha
          )

  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (file.path ('Plots',paste0('DAG_clusters_',metaGroupName,'_heatmaps_only_meso.pdf')), width = 2.5, height = 6)
print(DAG_hm)
dev.off()
}



### Call peaks ####
run_peakCall = TRUE
if (run_peakCall)
  {
  metaGroupName = 'Clusters_H'
  archp = addGroupCoverages (
    ArchRProj = archp, 
    groupBy = metaGroupName,  
    force = TRUE,
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

### Refine annotation ####
archp$celltype2 = 0
archp$celltype2[archp$Clusters_H %in% c('C3')] = 'NK_FGFBP2'
archp$celltype2[archp$Clusters_H %in% c('C2')] = 'NK'
archp$celltype2[archp$Clusters_H %in% c('C11')] = 'Tregs'
archp$celltype2[archp$Clusters_H %in% c('C12','C13','C5')] = 'CD4'
archp$celltype2[archp$Clusters_H %in% c('C10','C4','C6','C7','C8','C9')] = 'CD8'

p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "cellColData", 
    name = 'celltype2', 
    embedding = "UMAP_H"#,
    #pal = palette_expression
)

pdf (file.path('Plots','celltype2_umap.pdf'), width = 5, height = 5)
print (p)
dev.off()

### chromVAR analysis ####
run_chromVAR = TRUE

if (run_chromVAR)
  {  
  archp = addBgdPeaks (archp, force= TRUE)
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
  
### Run peak calling ####
metaGroupName = "Clusters"
force=FALSE
if(!all(file.exists(file.path('PeakCalls', paste0(unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds')))) | force) source ('../../PM_scATAC/callPeaks.R')
  



### chromVAR analysis ####
force=FALSE
if (!all(file.exists(file.path('Annotations',
  c('Motif-Matches-In-Peaks.rds',
    'Motif-Positions-In-Peaks.rds',
    'Motif-In-Peaks-Summary.rds')))))
source (file.path('..','..','PM_scATAC','chromVAR.R'))


# Find activating and repressing TFs #### 
if (!file.exists ('TF_activators_genescore.rds')) 
  {
    source (file.path('..','..','PM_scATAC','activeTFs.R'))
  } else {
    corGSM_MM = readRDS ('TF_activators_genescore.rds') 
  }


### Co-expression of TFs #### 
metaGroupName = 'Clusters_H'
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
all (colnames(mSE) == rownames(archp))

# # Get deviation matrix ####
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Subset only for positively correlated TF with genescore ####
positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0.1]
mMat = mMat[positive_TF,]

mMat_cor = cor (as.matrix(t(mMat)), method = 'pearson')

km = kmeans (mMat_cor, centers=5)

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

TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "cellColData",
    name = paste0('mod_',unique(km$cluster)), 
    pal = rev(palette_deviation),
    #useSeqnames='z',
    embedding = "UMAP_H")

pdf (file.path ('Plots','TF_modules_umap.pdf'), width = 10,height=4)
wrap_plots (TF_p, ncol=5)
dev.off()




### ChromVAR based analysis ####
run_chromVAR_analysis = TRUE

if (run_chromVAR_analysis)
  {
  # Find DAM ####
  metaGroupName = "Clusters_H"
  force = FALSE
  if (!file.exists (paste0('DAM_',metaGroupName,'.rds')) | force)
    {
    DAM_list = getMarkerFeatures (
      ArchRProj = archp, 
      testMethod = "wilcoxon",
            #useGroups = "ClusterA",
            #bgdGroups = "Clusters1B",
      binarize = FALSE,
      useMatrix = "MotifMatrix",
      groupBy = metaGroupName
    #  useSeqnames="z"
    )

    listnames = colnames (DAM_list)
    DAM_list = lapply (1:ncol (DAM_list), function(x) 
      {
      df = DAM_list[,x]  
      df = do.call (cbind, (assays(df)))
      colnames(df) = names (assays(DAM_list))
      df$gene = rowData (DAM_list)$name
      df
      })
    names (DAM_list) = listnames
    saveRDS (DAM_list, paste0 ('DAM_',metaGroupName,'.rds'))    
    } else {
    DAM_list = readRDS (paste0('DAM_',metaGroupName,'.rds'))
    }
  
  active_genes = corGSM_MM$MotifMatrix_name[corGSM_MM$cor > 0.1]
  DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% active_genes,])
  names (DAM_list2) = names (DAM_list)
  FDR_threshold = 1e-3
  meandiff_threshold = 0
  top_genes = 3
  DAM_top_list = DAM_list2[sapply (DAM_list2, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$MeanDiff) > meandiff_threshold,]) > 0)]
  DAM_top_list = lapply (seq_along(DAM_top_list), function(x) {
    res = DAM_top_list[[x]]
    #res = na.omit (res)
    res = res[res$FDR < FDR_threshold,]
    res = res[order (res$FDR), ]
    res = res[res$MeanDiff > meandiff_threshold,]
    res$comparison = names(DAM_top_list)[x]
    if (nrow(res) < top_genes) 
      {
      res
      } else {
      head (res,top_genes)
      }
    })
  DAM_df = Reduce (rbind ,DAM_top_list)
  
  
  DAM_df$gene = gsub ('_.*','',DAM_df$gene)
  DAM_df$gene = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", DAM_df$gene)
  mMat = assays (mSE)[[1]]
  rownames (mMat) = rowData (mSE)$name
  mMat_mg = mMat[DAM_df$gene, ]
  mMat_mg = as.data.frame (t(mMat_mg))
  mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
  mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
  rownames (mMat_mg) = mMat_mg[,1]
  mMat_mg = mMat_mg[,-1]
  mMat_mg = mMat_mg[names (DAM_list),]
  #mMat_mg = mMat_mg[names(table (archp@cellColData[,metaGroupName])[table (archp@cellColData[,metaGroupName]) > 50]),]
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
          col = rev(palette_deviation)

          #right_annotation = motif_ha
          )

  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_heatmaps.pdf')), width = 3, height = 5)
print(DAM_hm)
dev.off()
}




# TNK markers ####
meso_markers = c('CD8A','CTLA4','PDCD1','HAVCR2','TIGIT','TOX','GZMB','IL7R')
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

pdf (file.path('Plots','CD8_markers_fplots_only_normal.pdf'), width = 15, height = 8)
wrap_plots(p, ncol=4)
dev.off()



