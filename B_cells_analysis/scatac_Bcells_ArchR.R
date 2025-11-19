conda activate meso_scatac
R

set.seed(1234)

####### ANALYSIS of B compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/B_cells/scatac_ArchR'
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

# Set # of threads and genome reference ####
addArchRThreads(threads = 1) 
addArchRGenome("hg38")

archp = loadArchRProject (projdir)
srt = readRDS ('../scrna/srt.rds')
process = FALSE

if (process) 
{

## Reduce dimension and harmonize ####
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
    name = "Harmony_project",
    groupBy ='Sample', force=TRUE
)

archp = addUMAP (ArchRProj = archp,
    reducedDims = "Harmony_project", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony_project",
    name='Clusters_H', res=2,
    force = TRUE)

pdf()
umap_p3 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
  pal = palette_sample,
   embedding = "UMAP_H")
umap_p4 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype_revised2",
   embedding = "UMAP_H")
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters_H",
   embedding = "UMAP_H")
dev.off()

  pdf (file.path('Plots','celltype_umap_harmony_on_project_sample.pdf'),5,5)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()

### Remove C2 after Harmonized clustering ####
archp = archp[archp$Clusters_H != 'C2']
    archp = addIterativeLSI (ArchRProj = archp,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = TRUE, LSIMethod = LSI_method,
    varFeatures = varfeat)

  archp = addHarmony (
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    name = "Harmony_project",
    groupBy ='Sample', force=TRUE
)

archp = addUMAP (ArchRProj = archp,
    reducedDims = "Harmony_project", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony_project",
    name='Clusters_H', res=.8,
    force = TRUE)

archp = archp[archp$Clusters_H != 'C1'] # Further remove low quality clusters 

archp = addClusters (input = archp,
    reducedDims = "Harmony_project",
    name='Clusters_H', res=.6,
    force = TRUE)

archp = saveArchRProject (archp)

pdf()
umap_p3 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
  pal = palette_sample,
   embedding = "UMAP_H",
   labelMeans=FALSE)
umap_p4 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype_revised2",
   embedding = "UMAP_H")
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters_H",
   embedding = "UMAP_H")
dev.off()

  pdf (file.path('Plots','celltype_umap_harmony_on_project_sample.pdf'),5,5)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()

pdf()
qc_p = plotEmbedding(
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = c('TSSEnrichment','nFrags','FRIP'), 
    embedding = "UMAP_H",
   # pal = palette_deviation,
    imputeWeights = getImputeWeights(archp)
)
dev.off()

  pdf (file.path('Plots','QC_featp.pdf'),12,12)
  print (qc_p)
  dev.off()
}



### Cell Annotation ####
archp$celltype_lv2 = archp$Clusters_H
archp$celltype_lv2[archp$celltype_lv2 == 'C1'] = 'Plasma'
archp$celltype_lv2[archp$celltype_lv2 == 'C2'] = 'B_cell_Memory_C2'
archp$celltype_lv2[archp$celltype_lv2 == 'C3'] = 'B_cell_Memory_C3'
archp$celltype_lv2[archp$celltype_lv2 == 'C4'] = 'B_cell_Naive_C4'
archp$celltype_lv2[archp$celltype_lv2 == 'C5'] = 'B_cell_Naive_C5'
write.csv (data.frame (barcode = archp@cellColData, celltype_lv2 = archp$celltype_lv2), 'barcode_annotation.csv')

# Run genescore DAG ####
metaGroupName = "celltype_lv2"
force = T
if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source (file.path('..','..','git_repo','utils','DAG.R'))


### Compare expression of DAG in scatac vs scrna in B cells ####
DAG_list = readRDS (paste0('DAG_',metaGroupName,'.rds'))
FDR_threshold = .05
lfc_threshold = 0.5
top_genes = 5
DAG_top_list = DAG_list[sapply (DAG_list, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$Log2FC) > lfc_threshold,]) > 0)]
DAG_top_list = lapply (seq_along (DAG_top_list), function(x) {
  res = DAG_top_list[[x]]
  #res = na.omit (res)
  res = res[res$FDR < FDR_threshold,]
  res = res[order (res$FDR), ]
  res = res[res$Log2FC > lfc_threshold,]
  res$comparison = names(DAG_top_list)[x]
  if (nrow(res) < top_genes) 
    {
    res
    } else {
     head (res,top_genes)
    #res[order (-res$Log2FC),]
    }
  })
DAG_df = Reduce (rbind ,DAG_top_list)

if (!any (ls() == 'gsSE')) gsSE = ArchR::getMatrixFromProject (archp, useMatrix = 'GeneScoreMatrix')
gsSE = gsSE[, archp$cellNames]
gsMat = assays (gsSE)[[1]]
rownames (gsMat) = rowData (gsSE)$name
gsMat_mg = gsMat[rownames (gsMat) %in% unique(DAG_df$gene), ]
gsMat_mg = as.data.frame (t(gsMat_mg))
gsMat_mg$metaGroup = as.character(archp@cellColData[,metaGroupName])
gsMat_mg = aggregate (.~ metaGroup, gsMat_mg, mean)
rownames (gsMat_mg) = gsMat_mg[,1]
gsMat_mg = gsMat_mg[,-1]
gsMat_mg = gsMat_mg[names(table (archp@cellColData[,metaGroupName])[table (archp@cellColData[,metaGroupName]) > 50]),]
gsMat_mg = gsMat_mg[,unique(DAG_df$gene)]
DAG_hm = Heatmap (t(scale(gsMat_mg)), 
        row_labels = colnames (gsMat_mg),
        column_title = paste('top',top_genes),
        #clustering_distance_columns = 'euclidean',
        #clustering_distance_rows = 'euclidean',
        cluster_rows = F,
        #col = pals_heatmap[[5]],
        cluster_columns=F,#col = pals_heatmap[[1]],
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        column_names_rot = 45,
        #rect_gp = gpar(col = "white", lwd = .5),
        border=TRUE,
        col = palette_expression
        #right_annotation = motif_ha
        )
         
  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (paste0('Plots/DAG_clusters_',metaGroupName,'_heatmaps.pdf'), width = 3, height = 4)
print (DAG_hm)
dev.off()

# Show distribution of protocadherins compared to all DAGs (in B cells only) ####
archp_B = archp[archp$celltype_lv1 == 'B_cells']
top_genes = Inf
DAG_top_list = DAG_list[sapply (DAG_list, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$Log2FC) > lfc_threshold,]) > 0)]
DAG_top_list = lapply (seq_along(DAG_top_list), function(x) {
  res = DAG_top_list[[x]]
  #res = na.omit (res)
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
DAG_df = DAG_df[DAG_df$comparison != 'Plasma',]


if (!any (ls() == 'gsSE')) gsSE = fetch_mat (archp, 'GeneScore')
gsSE = gsSE[, archp$cellNames]
gsMat = assays (gsSE)[[1]]
rownames (gsMat) = rowData (gsSE)$name
gsMat_mg = gsMat[rownames (gsMat) %in% unique(DAG_df$gene), ]
gsMat_mg = as.data.frame (t(gsMat_mg))
gsMat_mg = colMeans (gsMat_mg)
srt = srt[,srt$celltype_lv1 == 'B_cells']
expmat = as.data.frame (log2(AverageExpression (srt, feature = names (gsMat_mg), group.by = 'celltype_lv1')[[1]]+1))

exp_gs = data.frame (expression = expmat[names(gsMat_mg),], genescore = gsMat_mg)
exp_gs$protocadherins = ifelse (grepl ('PCDHG',rownames(exp_gs)), 'protocadherin','other')
sp = ggplot (exp_gs, aes (x = expression, y = genescore)) + geom_point(color = 'grey44', size=0.9, alpha=.4) + 
geom_point(data = exp_gs[exp_gs$protocadherins == 'protocadherin',], color = 'hotpink', size=0.9) + 
ylim (c(0,1)) + xlim (c(0,0.5)) +
geom_smooth(method = "lm", se = FALSE, color = "grey44", size=0.3) +
  stat_cor(method = "pearson",
           label.x.npc = "left",  # place in left side of plot
           label.y.npc = "top")  +
  gtheme


pdf (file.path ('Plots','genescore_expression_bcells_markers_scatter.pdf'), 3,height=3)
sp
dev.off()


protocaderins_region = GRanges ('chr5:141322500-141507335')
metaGroupName = 'Clusters_H'

pdf()
meso_markers = plotBrowserTrack2 (
    ArchRProj = archp, 
    sizes = c(6, 1, 1, 6,1,1),
    groupBy = metaGroupName, 
#    geneSymbol = TF,
    normMethod = "ReadsInTSS",
    scCellsMax=3000,
    plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
    hubs_regions = NULL,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    region = protocaderins_region,
    upstream = 0,
    #pal = palette_tnk_cells_ext2,
    #ylim=c(0,0.1),
    downstream = 0,
    #loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
#    loops = getCoAccessibility (archp, corCutOff = 0.25),
    #  returnLoops = TRUE),
    useGroups= NULL,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2,returnLoops = TRUE),
    #hubs = hubs_obj$peakLinks2
)
dev.off()

plotPDF (meso_markers, ArchRProj = archp, 
  width=5,height=6, 
  name =paste0(paste('protocaderins', collapse='_'),'_coveragePlots.pdf'),
  addDOC = F)


### Call peaks on celltypes ####
metaGroupName = 'Clusters_H'
force=F
peak_reproducibility=2
if(!all(file.exists(file.path('PeakCalls', unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds'))) | force) 
source (file.path('..','..','git_repo','utils','callPeaks.R'))

### chromVAR analysis ####
force=F
if (!all(file.exists(file.path('Annotations',
  c('Motif-Matches-In-Peaks.rds',
    'Motif-Positions-In-Peaks.rds',
    'Motif-In-Peaks-Summary.rds')))) | force)
source (file.path ('..','..','git_repo','utils','chromVAR.R'))
  
# Differential Accessed motifs ####
metaGroupName = "Clusters_H"
force=FALSE
source (file.path('..','..','git_repo','utils','DAM.R'))

DAM_list = readRDS (paste0('DAM_',metaGroupName,'.rds'))
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames(mMat) = rowData(mSE)$name

# Clean TF names
  DAM_list = lapply (DAM_list, function(x)
       {
       x$gene = gsub ('_.*','',x$gene)
       x$gene = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", x$gene)
       x
       })
  

  if (metaGroupName %in% colnames(srt@meta.data))
  {
  #Get active genes from RNA
  ps = log2(as.data.frame (AverageExpression (srt, 
  features = sapply (unique(unlist(lapply(DAM_list, function(x) x$gene))), function(x) unlist(strsplit (x, '_'))[1]), 
  group.by = metaGroupName)[[1]]) +1)
  min_exp = .1
  ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
  active_TFs = rownames(ps)

  #active_genes = corGSM_MM$MotifMatrix_name[corGSM_MM$cor > 0.1]
  DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% active_TFs,])    
  } else {
  DAM_list2 = DAM_list  
  }
  names (DAM_list2) = names (DAM_list)
  FDR_threshold = 1e-2
  meandiff_threshold = 0
  top_genes = 5
  DAM_top_list = DAM_list2[sapply (DAM_list2, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$MeanDiff) > meandiff_threshold,]) > 0)]
  DAM_top_list = lapply (seq_along(DAM_top_list), function(x) {
    res = DAM_top_list[[x]]
    #res = na.omit (res)
    res = res[res$FDR < FDR_threshold,]
    res = res[order (res$FDR), ]
    #res = res[order (-res$MeanDiff), ]
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
  active_DAM = unique(DAM_df$gene)
  
  # # Get deviation matrix ####
  if (!exists ('mSE') | force) mSE = fetch_mat (archp, 'Motif')
  mMat = assays (mSE)[[1]]
  rownames (mMat) = rowData (mSE)$name
  mMat_mg = mMat[active_DAM, ]
  mMat_mg = as.data.frame (t(mMat_mg))
  mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
  mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
  rownames (mMat_mg) = mMat_mg[,1]
  mMat_mg = mMat_mg[,-1]
  mMat_mg = mMat_mg[names (DAM_list),]
  mMat_mg = na.omit (mMat_mg)
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
          name = 'TF activity',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_deviation_fun(scale(mMat_mg))
          #right_annotation = motif_ha
          )

#DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_heatmaps.pdf')), width = 3, height = 4)
print(DAM_hm)
dev.off()




















### Co-expression of TFs across cells #### 

### Run TF correlation to identify TF modules across TNK cells #### 
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat)#[selected_TF,])

# Filter by RNA expression ####
metaGroupName = 'seurat_clusters'
active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)
mMat = mMat[active_TFs, ]

mMat_cor = cor (as.matrix(t(scale(mMat))), method = 'spearman')

set.seed(1234)
centers=5
km = kmeans (mMat_cor, centers=centers)

pdf (file.path ('Plots','TF_modules_heatmap2.pdf'), width = 4,height=3)
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

pdf (file.path ('Plots','TF_modules_heatmap2.pdf'), width = 4,height=3)
cor_mMat_hm
dev.off()

tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = do.call (cbind, tf_modules)
archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% paste0('mod_',unique(km$cluster))]
archp@cellColData = cbind (archp@cellColData, tf_modules) 

pdf()
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "cellColData",
    name = paste0('mod_',unique(km$cluster)), 
    pal = rev(palette_deviation_correlation),
    #useSeqnames='z',
    embedding = "UMAP_H")
dev.off()
pdf (file.path ('Plots','TF_modules_umap2.pdf'), width = 20,height=6)
wrap_plots (TF_p, ncol=5)
dev.off()


library(ggplot2)
ccomp = as.data.frame (archp@cellColData)
gp1 = ggplot(ccomp, aes(x = mod_2, fill = Sample)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
    scale_fill_manual(values= palette_sample) + gtheme

gp2 = ggplot(ccomp, aes(x = mod_1, fill = Sample)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  scale_fill_manual(values= palette_sample) + gtheme

pdf (file.path ('Plots','TF_modules_dist_sample_hist.pdf'),7,3)
wrap_plots (gp1, gp2)
dev.off()

### Correlation of TF modules with seq quality ####
library (ggpubr)


sp <- lapply (unique(km$cluster), function(x) ggplot(ccomp, aes_string(x = 'nFrags', y = paste0('mod_',x))) + #, fill = sampleID, color = sampleID)) +
  geom_point(alpha = .8, shape = 21, stroke = .1) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson",
           label.x.npc = "left",  # place in left side of plot
           label.y.npc = "top") +
  gtheme)
sp1 <- lapply (unique(km$cluster), function(x) ggplot(ccomp, aes_string(x = 'TSSEnrichment', y = paste0('mod_',x))) + #, fill = sampleID, color = sampleID)) +
  geom_point(alpha = .8, shape = 21, stroke = .1) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson",
           label.x.npc = "left",  # place in left side of plot
           label.y.npc = "top") +
  gtheme)
sp2 <- lapply (unique(km$cluster), function(x) ggplot(ccomp, aes_string(x = 'FRIP', y = paste0('mod_',x))) + #, fill = sampleID, color = sampleID)) +
  geom_point(alpha = .8, shape = 21, stroke = .1) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson",
           label.x.npc = "left",  # place in left side of plot
           label.y.npc = "top") +
  gtheme)
pdf (file.path ('Plots','qc_scatter_cor_with_TF_modules_scatter.pdf'), height=3,width=12.5)
wrap_plots (sp, ncol=5)
wrap_plots (sp1, ncol=5)
wrap_plots (sp2, ncol=5)
dev.off()


# Check NR4A2 deviations
tf_markers = c('RUNX1')
tf_markers = c('NR4A2')
tf_markers = c('JUN')
markerMotifs = getFeatures (archp, select = paste(tf_markers, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs = grep ("z:", markerMotifs, value = TRUE)
#archp = addImputeWeights (archp)
pdf()
TF_p = plotEmbedding(
    ArchRProj = archp, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP_H",
   # pal = palette_deviation,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path ('Plots',paste0(tf_markers,'_umap.pdf')), width = 5,height=6)
TF_p
dev.off()









# Export bigwig files for IGV viewer
getGroupBW (
  ArchRProj = archp,
  groupBy = "Clusters_H",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 5000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)












































