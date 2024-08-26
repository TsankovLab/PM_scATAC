conda activate meso_scatac
use UGER
R

set.seed(1234)

####### ANALYSIS of NKT compartment #######
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR'
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

archp = loadArchRProject (projdir)
srt = readRDS (file.path ('..','scrna','srt.rds'))

# Add metadata ####
archp$status = ifelse (grepl ('^P',archp$Sample2), 'tumor','normal')  
archp$project = ifelse (grepl ('^P', archp$Sample2), 'PM','normal')

trim_clusters = FALSE
if (trim_clusters) 
{

## Reduce dimension and harmonize ####

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
  
  pdf (file.path('Plots','celltype_umap_harmony_on_project_sample.pdf'),5,5)
  print (umap_p1)
  print (umap_p2)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()

## Further remove outlier clusters ####
archp = archp[!archp$Clusters_H %in% c('C10','C1')]


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
  
  pdf (file.path('Plots','celltype_umap_harmony_on_project_sample.pdf'),5,5)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()

# Remove cluster C1 containing only 6 cells ####
archp = archp[!archp$Clusters_H %in% c('C1')]
}


# Recluster only using meso samples ####
archp = archp[archp$Sample2 %in% c('P1','P10','P11','P12','P13','P14','P3','P4','P5','P8')]

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

  umap_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
   name = "celltype", embedding = "UMAP")
  umap_p2 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample2",
     embedding = "UMAP")
  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample2", labelMeans =F,
     embedding = "UMAP_H", pal = palette_sample)
  umap_p5 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters_H",
     embedding = "UMAP_H")
  umap_p3 = umap_p3 + theme_void()
  
  pdf (file.path('Plots','celltype_umap_signac_filtered_harmony_on_project3.pdf'),5,5)
  print (umap_p1)
  print (umap_p2)
  print (umap_p3)
  print (umap_p5)
  dev.off()

pdf (file.path('Plots','sample_celltype_umaps.pdf'),5,5)
  print (umap_p3)
  print (umap_p4)
  dev.off()

# TNK markers ####
tnk_markers = c('CD3D','CD8A','PDCD1','CD4', 'FOXP3','GNLY',
  'FGFBP2','KLRC1')
archp = addImputeWeights (archp)
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = tnk_markers, 
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
p = lapply (seq_along(p), function(x) p[[x]] + theme_void() + ggtitle (tnk_markers[x]) + NoLegend())
#archp$celltype[archp$Clusters == 'C30'] = 'Fibroblasts_WT1'
#p = lapply (p, function(x) x + theme_void() + NoLegend ()) #+ ggtitle scale_fill_gradient2 (rev (viridis::plasma(100))))

pdf (file.path('Plots','TNK_markers_fplots.pdf'), width = 10, height = 10)
print (wrap_plots(p, ncol=4))
dev.off()

### Annotate meso cells ####
archp$celltype2 = 0
archp$celltype2[archp$Clusters_H %in% c('C6')] = 'NK_FGFBP2'
archp$celltype2[archp$Clusters_H %in% c('C5')] = 'NK_KLRC1'
archp$celltype2[archp$Clusters_H %in% c('C2')] = 'Tregs'
archp$celltype2[archp$Clusters_H %in% c('C3','C4','C7')] = 'CD4'
archp$celltype2[archp$Clusters_H %in% c('C1','C8','C9','C10','C11')] = 'CD8'


# Subset CD8 cells ####
metaGroupName = 'celltype2'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% 'CD8'],
  outputDirectory = file.path('..','CD8'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)


cd8_ct = read.csv (file.path('..','..','CD8','scatac_ArchR','barcode_annotation.csv')) # get annotation from subclustered CD8 cells
archp$celltype2[match(cd8_ct$barcode, rownames(archp@cellColData))] = cd8_ct$celltype
archp$celltype2[archp$celltype2 %in% c('C1','C2','C4')] = 'CD8'

umap_p4 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "celltype2",
     embedding = "UMAP_H")
  
pdf (file.path('Plots','celltype_umap_signac_filtered_harmony_on_project_meso_only.pdf'),5,5)
print (umap_p4)
dev.off()

# Check expression of GZMB PRF1 and KLRC1 ####
metaGroupName = 'celltype2'
archp = addImputeWeights (archp)
markers = c('GZMA','GZMB','PRF1','PDCD1','CTLA4','TIGIT','KLRC1')
p2 <- plotGroups(
    ArchRProj = archp, 
    groupBy = metaGroupName, 
    colorBy = "GeneScoreMatrix", 
    name = markers,
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE,
    pal = palette_tnk_cells
   )
p2 = lapply (p2, function(x) x + gtheme)
pdf (file.path ('Plots','NK_cytotoxicity_markers.pdf'))
wrap_plots (p2)
dev.off()


### Call peaks on celltypes ####
metaGroupName = 'celltype2'
force=FALSE
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


# Differential Accessed motifs ####
metaGroupName = "celltype2"
force=FALSE
source (file.path('..','..','PM_scATAC','utils','DAM.R'))

mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat_mg = mMat[selected_TF, ]
mMat_mg = as.data.frame (t(mMat_mg))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
mMat_mg = mMat_mg[names (DAM_list),]

# Generate RNA pseudobulk of matching cell types ####
metaGroupName = 'celltype2'
selected_TF = c(rownames(DAM_hm@matrix), 'NR4A3','NR4A2','NR4A1')
ps = log2(as.data.frame (AverageExpression (srt, features = selected_TF, group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
ps = ps[, colnames(DAM_hm@matrix)]
ps_tf = ps[selected_TF,]

  
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
          col = rev(palette_deviation),
          width = unit(2, "cm")
          #right_annotation = motif_ha
          )

scaled_ps = t(scale(t(ps_tf)))
scaled_ps[is.na(scaled_ps)] = 0
TF_exp_selected_hm = Heatmap (scaled_ps,
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

pdf (file.path ('Plots','DAM_with_rna_expression_heatmaps.pdf'), width = 8,height=4)
draw (DAM_hm + TF_exp_selected_hm + TF_exp_selected_hm2)
dev.off()
   

### Co-expression of TFs across cells #### 

# # Get deviation matrix ####
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Subset only for expressed TFs ####
metaGroupName = 'celltype2'
ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat), group.by = metaGroupName)[[1]]) +1)
min_exp = 0.1
ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
selected_TF = rownames(ps)[rowSums(ps) > 0]
#positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0]
metaGroupName = 'celltype2'
mMat = mMat[selected_TF,]
mMat = as.data.frame (t(mMat))
mMat$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat = aggregate (.~ metaGroup, mMat, mean)
rownames (mMat) = mMat[,1]
mMat = mMat[,-1]

TF_hm = draw(Heatmap (scale(mMat), 
          #row_labels = colnames (mMat_mg),
          #column_title = paste('top',top_genes),
          clustering_distance_columns = 'pearson',
          clustering_distance_rows = 'pearson',
          cluster_rows = T,
          column_km = 20,
          #col = pals_heatmap[[5]],
          cluster_columns=T,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
          column_names_gp = gpar(fontsize = 0),
          column_names_rot = 45,
          name = 'chromVAR',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_deviation_centered#,
          #width = unit(2, "cm")
          #right_annotation = motif_ha
          ))

pdf (file.path ('Plots',paste0('TF_',metaGroupName,'heatmap.pdf')), width=8, height=1.6)
TF_hm
dev.off()

colnames(TF_hm@ht_list$chromVAR@matrix)[unlist(column_order(TF_hm)[c('2','3','4','5')])]
which (colnames(mMat) == 'NR4A2')
sapply (column_order(TF_hm), function(x) 497 %in% x)

# Distance matrix ####
d <- as.dist(1 - cor(t(mMat), method='pearson'))

# Hierarchical clustering ####
hc <- hclust(d)

# Dendrogram ####
pdf (file.path ('Plots',paste0('TF_',metaGroupName,'_no_km_dendrogram.pdf')), width=3, height=3.6)
plot(hc)
dev.off()

# # Generate RNA pseudobulk of matching cell types ####
# metaGroupName = 'celltype2'
# selected_TF = colnames(TF_hm@ht_list$chromVAR@matrix)[unlist(column_order(TF_hm)[c('2','3','4','5')])]
# ps = log2(as.data.frame (AverageExpression (srt, features = selected_TF, group.by = metaGroupName)[[1]]) +1)
# colnames (ps) = gsub ('-','_',colnames(ps))
# #ps = ps[, colnames(DAM_hm@matrix)]
# ps_tf = ps[selected_TF,]

# metaGroupName = 'celltype2'
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat_mg = mMat[selected_TF, ]
# mMat_mg = as.data.frame (t(mMat_mg))
# mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
# mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
# rownames (mMat_mg) = mMat_mg[,1]
# mMat_mg = mMat_mg[,-1]
  
#  DAM_hm = Heatmap (t(scale(mMat_mg)), 
#           row_labels = colnames (mMat_mg),
          
#           clustering_distance_columns = 'euclidean',
#           clustering_distance_rows = 'euclidean',
#           cluster_rows = F,
#           #col = pals_heatmap[[5]],
#           cluster_columns=F,#col = pals_heatmap[[1]],
#           row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
#           column_names_gp = gpar(fontsize = 8),
#           column_names_rot = 45,
#           name = 'chromVAR',
#           #rect_gp = gpar(col = "white", lwd = .5),
#           border=TRUE,
#           col = rev(palette_deviation),
#           width = unit(2, "cm")
#           #right_annotation = motif_ha
#           )

# scaled_ps = t(scale(t(ps_tf)))
# scaled_ps[is.na(scaled_ps)] = 0
# TF_exp_selected_hm = Heatmap (scaled_ps,
#         #right_annotation=tf_mark,
#         #column_split = column_split_rna,
#         cluster_rows = F, #km = 4, 
#         name = 'expression',
#         column_gap = unit(.5, "mm"),
#         row_gap = unit(.2, "mm"),
#         clustering_distance_rows = 'euclidean',
#         clustering_distance_columns = 'euclidean',
#         cluster_columns=F, 
#         col = palette_expression,
#         row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
#         column_names_gp = gpar(fontsize = 8),
#           column_names_rot = 45,
#         border=T,
#         width = unit(2, "cm"))

# TF_exp_selected_hm2 = Heatmap (ps_tf,
#         #right_annotation=tf_mark,
#         #column_split = column_split_rna,
#         cluster_rows = F, #km = 4, 
#         name = 'expression',
#         column_gap = unit(.5, "mm"),
#         row_gap = unit(.2, "mm"),
#         clustering_distance_rows = 'euclidean',
#         clustering_distance_columns = 'euclidean',
#         cluster_columns=F, 
#         col = palette_expression,
#         row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
#         column_names_gp = gpar(fontsize = 8),
#           column_names_rot = 45,
#         border=T,
#         width = unit(2, "cm"))

# pdf (file.path ('Plots','NK_KLRC1_CD8_ext_with_rna_expression_heatmaps.pdf'), width = 8,height=14)
# draw (DAM_hm + TF_exp_selected_hm + TF_exp_selected_hm2)
# dev.off()
   

# Run DAM NK KLRC1 + CD8 ext vs FGFBP2 + CD8 ####
library (presto)
metaGroupName = 'celltype2'

if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = mMat[selected_TF,]
mMat = as.data.frame (mMat)

mMat_ext = mMat[,as.character(archp@cellColData[,metaGroupName]) %in% c('CD8','CD8_exhausted','NK_KLRC1','NK_FGFBP2')]
ext_vs_ctx = as.character(archp@cellColData[,metaGroupName])[as.character(archp@cellColData[,metaGroupName]) %in% c('CD8','CD8_exhausted','NK_KLRC1','NK_FGFBP2')]
ext_vs_ctx = ifelse (ext_vs_ctx %in% c('CD8','NK_FGFBP2'),'ctx','ext')
res = wilcoxauc (mMat_ext, ext_vs_ctx)

res_l = lapply (split (res, res$group), function(x){
  tmp = x[x$logFC > 0,]
  tmp = tmp[order (tmp$pval),]
  tmp
})

tf_ext = res_l[['ext']]$feature[res_l[['ext']]$padj < 0.01]

# Try just taking highest mean of NK KLRC1 + CD8 exhausted ####
metaGroupName = 'celltype2'

if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = mMat[selected_TF,]
mMat = as.data.frame (mMat)

mMat_mg = as.data.frame (t(mMat))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
mMat_mg = t (mMat_mg)

ext_vs_ctx_mean = rowMeans (mMat_mg[,c('CD8_exhausted','NK_KLRC1')])
ext_vs_ctx_mean = ext_vs_ctx_mean[order (-ext_vs_ctx_mean)]
selected_TF = names(head(ext_vs_ctx_mean,100))
#mMat_mg = mMat_mg[names (DAM_list),]

ps = log2(as.data.frame (AverageExpression (srt, features = selected_TF, group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
#ps = ps[, colnames(DAM_hm@matrix)]
ps_tf = ps[selected_TF,c('CD8_exhausted','NK_KLRC1')]
ps_tf_mean = rowMeans(ps_tf)
top_exp_tf = head(names(ps_tf_mean)[order(-ps_tf_mean)],10)
mMat_mg = mMat_mg[selected_TF,c('CD8_exhausted','NK_KLRC1')]

# Import table of cor tfs to KLRC1 and PDCD1 in nk and cd8
nk_cd8_ext_cor = read.csv (file.path('..','scrna','nk_cd8_ext_cor.csv'), row.names=1)

cd8_bar = HeatmapAnnotation (mark = anno_mark(at = match(top_exp_tf,names(ps_tf_mean)), 
  labels_gp=gpar(fontsize = 6, fontface = 'italic'), 
    labels = top_exp_tf,side='left'),  
'  ' = anno_barplot (-ps_tf[,c('CD8_exhausted')],border=F,gp = gpar(color = "white")),
  ' ' = nk_cd8_ext_cor[rownames(ps_tf), 'cd8'],
  which='row', col=list(' ' = palette_expression_cor_fun),
   simple_anno_size = unit(.2, "cm"))
nk_bar = HeatmapAnnotation (
  ' ' = nk_cd8_ext_cor[rownames(ps_tf), 'nk'],
  '  ' = anno_barplot(ps_tf[,c('NK_KLRC1')],border=F,gp = gpar(color = "white")), which='row',
   col=list(' ' = palette_expression_cor_fun),
    simple_anno_size = unit(.2, "cm"))


TF_hm = Heatmap (mMat_mg, 
          #row_labels = colnames (mMat_mg),
          #column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          left_annotation = cd8_bar,
          right_annotation = nk_bar,
          row_names_side = 'left',
          #column_km = 20,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 0, fontface = 'italic'),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          name = 'chromVAR',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_deviation_fun (mMat_mg[selected_TF,])#,
          #width = unit(1, "cm")
          #right_annotation = motif_ha
          )

# # Generate RNA pseudobulk of matching cell types ####

# TF_exp_selected_hm2 = Heatmap (ps_tf[,colnames(mMat_mg)],
#         #right_annotation=tf_mark,
#         #column_split = column_split_rna,
#         cluster_rows = F, #km = 4, 
#         name = 'expression',
#         column_gap = unit(.5, "mm"),
#         row_gap = unit(.2, "mm"),
#         clustering_distance_rows = 'pearson',
#         clustering_distance_columns = 'pearson',
#         cluster_columns=F, 
#         col = palette_expression,
#         row_names_gp = gpar(fontsize = 5, fontface = 'italic'),
#         column_names_gp = gpar(fontsize = 8),
#           column_names_rot = 45,
#         border=T,
#         width = unit(2, "cm"))

pdf (file.path ('Plots','chromvar_rna_expression_NK_CD8_EXT_heatmaps.pdf'), width = 3,height=3)
draw (TF_hm)# + TF_exp_selected_hm2)
dev.off()
   




### Co-expression of TFs between themself across cells #### 

# # Get deviation matrix ####
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
all (colnames(mSE) == rownames(archp))
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Subset only for positively correlated TF with genescore ####
positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0]
mMat = mMat[positive_TF,]

mMat_cor = cor (as.matrix(t(mMat)), method = 'pearson')
km = kmeans (mMat_cor, centers=20)

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

pdf (file.path ('Plots','TF_modules_meso_only.pdf'), width = 4,height=3)
cor_mMat_hm
dev.off()

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

pdf (file.path ('Plots','TF_modules_umap_meso_only.pdf'), width = 30, height=14)
wrap_plots (TF_p, ncol=5)
dev.off()



# Try with ridge plots ####
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
pdf (file.path ('Plots','TF_modules_ridge_plots_meso_only.pdf'), width = 30,height=8)
wrap_plots (rp, ncol=5)
dev.off()






# Try with ridge plots ####
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

### Plot cell type markers on genome tracks ####
metaGroupName = 'celltype2'
celltype_markers = c('PDCD1','TIGIT','TOX','HAVCR2','CTLA4')
#celltype_markers = c('WT1','CALB2','GATA4','MSLN','KRT5','KRT18','ITLN1','HP','SOX9')
meso_markers <- plotBrowserTrack(
    ArchRProj = archp, 
    groupBy = metaGroupName, 
    geneSymbol = celltype_markers,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 250000,
    downstream = 250000,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
    #loops = getCoAccessibility (archp, corCutOff = 0.3,
    #  returnLoops = TRUE),
    useGroups= NULL
)
plotPDF (meso_markers, ArchRProj = archp, width=14, name ='MPM_markers_coveragePlots.pdf')
}


# ### Use P2G analysis and cNMF from RNA to identify active TF via regulons  ####
run_p2g_TF = T

if (run_p2g_TF)
  {
  run_p2g = T  
  if (run_p2g)
    {
    maxDist = 250000
    archp = addPeak2GeneLinks(
        ArchRProj = archp,
        useMatrix = 'GeneScoreMatrix',
        reducedDims = "IterativeLSI",
        maxDist = maxDist
    )
    }  
    
  p2g_corr = .2
  p2g = getPeak2GeneLinks(
      ArchRProj = archp,
      corCutOff = p2g_corr,
      resolution = 1,
      returnLoops = FALSE
  )
}  


### Plot deviation and expression of NR4A2 ####
metaGroupName = 'celltype2'
selected_TF = 'NR4A2'
selected_TF = 'IRF4'
selected_TF = 'EOMES'
selected_TF = 'CBFB'
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat_mg = mMat[selected_TF, , drop=F]
mMat_mg = as.data.frame (t(mMat_mg))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = setNames(mMat_mg[,2], rownames(mMat_mg))

# Generate RNA pseudobulk of matching cell types ####
metaGroupName = 'celltype2'
ps = log2(as.data.frame (AverageExpression (srt, features = selected_TF, group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))

ps_df = data.frame (deviation = mMat_mg, expression = unlist(ps[names(mMat_mg)]))

ps_sp = ggplot (ps_df, aes (x = deviation, y = expression, fill= rownames(ps_df), label = rownames(ps_df))) + 
  geom_point (aes (fill = rownames(ps_df)), shape=21, size=5) + 
  scale_fill_manual (values = palette_tnk_cells) + 
  gtheme_no_rot + geom_text_repel (size=3.5, data = ps_df, aes(label = rownames(ps_df))) 

pdf (file.path ('Plots',paste0(selected_TF,'_dev_exp_scatterplot.pdf')), width=4,height=2.5)
ps_sp
dev.off()






  
