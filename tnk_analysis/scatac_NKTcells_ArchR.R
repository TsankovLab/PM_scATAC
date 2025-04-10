conda activate meso_scatac
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

# Set # of threads and genome reference ####
addArchRThreads(threads = 8) 
addArchRGenome("hg38")

archp = loadArchRProject (projdir)
srt = readRDS (file.path ('..','scrna','srt.rds'))


trim_clusters = FALSE
if (trim_clusters) 
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

# Run genescore DAG ####
metaGroupName = "Clusters_H"
force = TRUE
if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source (file.path('..','..','git_repo','utils','DAG.R'))

# TNK markers ####
tnk_markers = c('CD3D','CD8A','PDCD1','HAVCR2','CD4', 'FOXP3','GNLY',
  'FGFBP2','KLRC1','XCL1','ICOS')
tnk_markers = c('CD4','CD8A',
  'PDCD1','FOXP3','GNLY', 'HAVCR2',
  'FGFBP2','KLRC1','CTLA4')#,'ICOS')
archp = addImputeWeights (archp)
pdf()
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = tnk_markers, 
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
p = lapply (seq_along(p), function(x) p[[x]] + theme_void() + ggtitle (tnk_markers[x]) + NoLegend())
#archp$celltype[archp$Clusters == 'C30'] = 'Fibroblasts_WT1'
#p = lapply (p, function(x) x + theme_void() + NoLegend ()) #+ ggtitle scale_fill_gradient2 (rev (viridis::plasma(100))))

pdf (file.path('Plots','TNK_markers_fplots.pdf'), width = 7, height = 13)
print (wrap_plots(p, ncol=3))
dev.off()

### Annotate meso cells ####
archp$celltype2 = 0
archp$celltype2[archp$Clusters_H %in% c('C7')] = 'NK_FGFBP2'
archp$celltype2[archp$Clusters_H %in% c('C8','C9')] = 'NK_KLRC1'
archp$celltype2[archp$Clusters_H %in% c('C16')] = 'Tregs'
archp$celltype2[archp$Clusters_H %in% c('C14','C20','C11','C17','C12','C13')] = 'CD4'
archp$celltype2[archp$Clusters_H %in% c('C1','C5','C6','C18','C15','C19','C10')] = 'CD8'
archp$celltype2[archp$Clusters_H %in% c('C3','C2','C4')] = 'CD8_exhausted'

write.csv (data.frame (barcode = rownames(archp@cellColData), celltype = archp$celltype2), 'barcode_annotation.csv')





# TNK markers ####
meso_markers = c('CD8A','CTLA4','PDCD1','HAVCR2','TIGIT','TOX','GZMB','IL7R','KLRC1','GNLY','ICOS')
#archp = addImputeWeights (archp)
qc_param = c('nFrags','TSSEnrichment','ReadsInTSS')

sams = unique(archp$Sample)
sams = c('P10','P13','P14','P23')
for (sam in sams)
{
archp_sam = archp[archp$Sample == sam]  
varfeat = 25000
  LSI_method = 2
  archp_sam = addIterativeLSI (ArchRProj = archp_sam,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = TRUE, LSIMethod = LSI_method,
    varFeatures = varfeat)

archp_sam = addUMAP (ArchRProj = archp_sam, 
    reducedDims = "IterativeLSI", name='UMAP',
    force = TRUE)

archp_sam = addClusters (input = archp_sam,
    reducedDims = "IterativeLSI",
    name='Clusters_H',
    force = TRUE)

pdf()
umap_p3 = plotEmbedding (ArchRProj = archp_sam, 
  colorBy = "cellColData", name = "Sample",
   embedding = "UMAP")
umap_p4 = plotEmbedding (ArchRProj = archp_sam, 
  colorBy = "cellColData", name = "celltype2",
   embedding = "UMAP")
umap_p5 = plotEmbedding (ArchRProj = archp_sam, 
  colorBy = "cellColData", name = "Clusters_H",
   embedding = "UMAP")
dev.off()

pdf (file.path('Plots',paste0(sam,'_celltype_harmony_sample_umap.pdf')),5,5)
print (umap_p3)
print (umap_p4)
print (umap_p5)
dev.off()  

archp_sam = addImputeWeights (archp_sam)

pdf()
p <- plotEmbedding(
    ArchRProj = archp_sam,
    colorBy = "GeneScoreMatrix", 
    name = meso_markers,
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights (archp_sam)
)
dev.off()

pdf (file.path('Plots',paste0('TNK_',sam,'.pdf')), width = 18, height = 12)
print(wrap_plots (p, ncol=4))
dev.off()
}

pdf()
p3 <- plotGroups(
    ArchRProj = archp, 
    groupBy = "Clusters_H", 
    colorBy = "GeneScoreMatrix", 
    name = c('ICOS','PDCD1','HAVCR2','KLRC1','LAG3','TIGIT'),
    plotAs = "violin",
    #pal = palette_sample,
    alpha = 0.4,
    addBoxPlot = TRUE
   )
dev.off()

pdf (file.path ('Plots','ext_genescore_clusters.pdf'), height=8, width=10)
wrap_plots (p3, ncol=3)
dev.off()


# Check QC in P10 to assess if CD8 exhausted KLRC1+ are doublets
qc_param = c('TSSEnrichment','nFrags','ReadsInTSS')
archp$celltype2_smaple = paste0(archp$celltype2, '_',archp$Sample)
pdf()
p3 <- lapply (unique(archp$Sample), function (x) plotGroups(
    ArchRProj = archp[archp$Sample == x], 
    groupBy = "celltype2_smaple", 
    colorBy = "cellColData", 
    name = qc_param,
    plotAs = "violin",
    pal = palette_sample,
    alpha = 0.4,
    addBoxPlot = TRUE
   ))
dev.off()
pdf (file.path ('Plots','qc_celltype_samples.pdf'), height=4, width=10)
lapply (p3, function (x) wrap_plots (x, ncol=3))
dev.off()

cp = cellComp (as.data.frame (archp@cellColData[,c('Sample','celltype_revised2')]),
   metaGroups = c('Sample','celltype_revised2'),
   subset = 'CD8_exhausted',
   prop=FALSE,
   plot_as = 'bar',
   pal = palette_tnk_cells
   ) + gtheme
pdf (file.path ('Plots','cell_composition_P23_exluded.pdf'))
cp
dev.off()


# # Subset CD8 cells ####
metaGroupName = 'celltype2'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('CD8','CD8_exhausted')],
  outputDirectory = file.path('..','..','CD8'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)


# cd8_ct = read.csv (file.path('..','..','CD8','scatac_ArchR','barcode_annotation.csv')) # get annotation from subclustered CD8 cells
# archp$celltype2[match(cd8_ct$barcode, rownames(archp@cellColData))] = cd8_ct$celltype
# archp$celltype2[archp$celltype2 %in% c('C1','C2','C4')] = 'CD8'
pdf()
umap_p1 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "celltype2",
     embedding = "UMAP_H",labelMeans =F,
     pal = palette_tnk_cells
     )
umap_p2 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample", labelMeans =F,
     embedding = "UMAP_H", pal = palette_sample)
dev.off()

pdf (file.path('Plots','celltype_umap_harmony.pdf'),5,width=8)
print (wrap_plots (umap_p1,umap_p2))
dev.off()

# Export annotation
write.csv (data.frame (barcode= rownames(archp@cellColData), celltype = archp$celltype2), 'barcode_annotation.csv')


# Check expression of GZMB PRF1 and KLRC1 ####
metaGroupName = 'celltype2'
archp = addImputeWeights (archp)
markers = c('GZMA','GZMB','PRF1','PDCD1','HAVCR2','CTLA4','TIGIT','KLRC1')
activating_markers = c('GZMA','GZMB','GZMH','PRF1','ITGAE','ITGA1')
inhibitory_markers = c('KIR3DL1','KIR2DL3')
pdf()
p2 <- plotGroups(
    ArchRProj = archp,#[archp$celltype2 %in% c('NK_KLRC1','NK_FGFBP2')], 
    groupBy = metaGroupName, 
    colorBy = "GeneScoreMatrix", 
    name = markers,
    plotAs = "violin",
    alpha = 0.4,
    imputeWeights = getImputeWeights(archp),
    addBoxPlot = TRUE,
    pal = palette_tnk_cells
   )
dev.off()
p2 = lapply (p2, function(x) x + gtheme)
p2 = lapply (p2, function(x) {x$data = x$data[x$data[,1] %in% c('NK_FGFBP2','NK_KLRC1'),]; x})
pdf (file.path ('Plots','NK_cytotoxicity_markers.pdf'), width=8,height=3)
wrap_plots (p2, ncol=6)
dev.off()


### Call peaks on celltypes ####
metaGroupName = 'celltype2'
force=TRUE
peak_reproducibility=2
if(!all(file.exists(file.path('PeakCalls', unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds'))) | force) 
source (file.path('..','..','git_repo','utils','callPeaks.R'))

### chromVAR analysis ####
force=TRUE
if (!all(file.exists(file.path('Annotations',
  c('Motif-Matches-In-Peaks.rds',
    'Motif-Positions-In-Peaks.rds',
    'Motif-In-Peaks-Summary.rds')))) | force)
source (file.path ('..','..','git_repo','utils','chromVAR.R'))
  


# Differential Accessed motifs ####
metaGroupName = "celltype2"
force=FALSE
source (file.path('..','..','git_repo','utils','DAM.R'))

if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames(mMat) = rowData(mSE)$name

# Filter by RNA expression ####
metaGroupName = 'celltype2'
active_TFs = exp_genes (srt, active_DAM, min_exp = 0.1, metaGroupName)
mMat = mMat[active_TFs, ]

#mMat_mg = mMat[active_DAM, ]
mMat_mg = as.data.frame (t(mMat))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
mMat_mg = mMat_mg[names (DAM_list),]

# Generate RNA pseudobulk of matching cell types ####
metaGroupName = 'celltype2'
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
draw (DAM_hm) # + TF_exp_selected_hm + TF_exp_selected_hm2)
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
active_TFs = rownames(ps)

mMat = mMat[active_TFs,]

set.seed (1234)
km = kmeans (t(scale(mMat)), centers=3)

ha = HeatmapAnnotation (df = data.frame (
  celltype = as.character(archp@cellColData[,metaGroupName]),
  sample = archp$Sample), col=list (celltype = palette_tnk_cells, sample=palette_sample))

TF_hm = Heatmap (scale(mMat), 
          top_annotation= ha,
          #row_labels = colnames (mMat_mg),
          #column_title = paste('top',top_genes),
          clustering_distance_columns = 'pearson',
          clustering_distance_rows = 'pearson',
          column_split = km$cluster,
          cluster_rows = T,
#          column_km = 3,
          #col = pals_heatmap[[5]],
          cluster_columns=T,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 0, fontface = 'italic'),
          column_names_gp = gpar(fontsize = 0),
          column_names_rot = 45,
          name = 'chromVAR',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_deviation_centered#,
          #width = unit(2, "cm")
          #right_annotation = motif_ha
          )

pdf (file.path ('Plots',paste0('TF_',metaGroupName,'heatmap2.pdf')), width=12, height=5)
TF_hm
dev.off()

archp$cell_kmeans = paste0('km',km$cluster)
pdf()
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "cell_kmeans",
   embedding = "UMAP_H")
dev.off()

pdf (file.path('Plots','celltype_kmeans_umap.pdf'),5,5)
print (umap_p5)
dev.off()





### Run TF correlation to identify TF modules across TNK cells #### 
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat)#[selected_TF,])

# Filter by RNA expression ####
metaGroupName = 'celltype2'
active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)
mMat = mMat[active_TFs, ]

mMat_cor = cor (as.matrix(t(scale(mMat))), method = 'spearman')

set.seed(1234)
centers=3
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




# Differential Peaks in CD8 exhausted vs CD8 ####
# Find DAP ####
#force = FALSE
metaGroupName = 'celltype2'
force=F
if (!file.exists ('DAP_CD8_CD8_ext_pairwise.rds') | force)
  {
  DAP_list = getMarkerFeatures (
    ArchRProj = archp, 
    testMethod = "wilcoxon",
          useGroups = "CD8_exhausted",
          bgdGroups = "CD8",
    k=100,
    binarize = FALSE,
    useMatrix = "PeakMatrix",
    groupBy = metaGroupName
  #  useSeqnames="z"
  )
  saveRDS (DAP_list, 'DAP_CD8_CD8_ext_pairwise.rds')
  } else {
  DAP_list = readRDS ('DAP_CD8_CD8_ext_pairwise.rds')
  }
DAP_res = do.call (cbind, (assays(DAP_list)))
colnames (DAP_res) = names(assays(DAP_list))
DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames(DAP_res) = as.character(DAP_res_regions)

pdf(file.path('Plots','Exhausted_peaks_MA_plot.pdf'), width=5,height=5)
pma <- markerPlot(seMarker = DAP_list, name = 'CD8_exhausted', cutOff = "FDR <= 0.01", plotAs = "MA")
pma
dev.off()

# Take only significant regions ####
DAP_res_sig = DAP_res[DAP_res$FDR < .01 & DAP_res$Log2FC > 0, ]
saveRDS (GRanges(rownames(DAP_res_sig)), 'T_cell_exhaustion_peaks.rds')


# Run GSEA enrichment analysis #### 
#!!! To run this analysis load only ArcHR and clusterprofiler packages !!!!
library (fgsea)    
options(warn = 0)
ps = getPeakSet (archp)
gmt.file = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_DN.v2024.1.Hs.gmt'
gmt.file = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP.v2024.1.Hs.gmt'
gmt.file = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/GSE24026_PD1_LIGATION_VS_CTRL_IN_ACT_TCELL_LINE_UP.v2024.1.Hs.gmt'
csv.file = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/Tcell_exhastuion_genes_PMID37091230.csv'
pathways = read.gmt (gmt.file)
pathways = read.csv (csv.file)
pathways = list(Tcell_exhastuion_genes_PMID37091230 = pathways[,1])
#pathways = gmtPathways (gmt.file)
DAP_list = readRDS ('DAP_CD8_CD8_ext_pairwise.rds')
DAP_res = do.call (cbind, (assays(DAP_list)))
colnames (DAP_res) = names(assays(DAP_list))
DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames(DAP_res) = as.character(DAP_res_regions)

peak_genes = unname(ps[queryHits (findOverlaps(ps, GRanges (rownames(DAP_res))))]$nearestGene)
names (peak_genes) = as.character(ps)[queryHits(findOverlaps (ps,GRanges (rownames(DAP_res))))]
peak_genes = peak_genes[rownames(DAP_res)]
peak_genes = setNames (-log10(DAP_res$Pval) * sign (DAP_res$Log2FC), peak_genes)
#names (peak_genes) = 
peak_genes = peak_genes[!duplicated(names(peak_genes))]
peak_genes = peak_genes[!is.na(names(peak_genes))]
fgseaRes = fgseaMultilevel (pathways, 
          peak_genes,#, 
          minSize=15, 
          maxSize=1500,
          nPermSimple=100000,
          BPPARAM = NULL
          )
pdf (file.path ('Plots','DAP_exhaustion_enrichment_plot.pdf'), width=5, height=3)
plotEnrichment(pathways[["Tcell_exhastuion_genes_PMID37091230"]],
               peak_genes) + labs(title="Tcell_exhastuion_genes_PMID37091230")
dev.off()


# Calculate hypergeometric test of peaks with motifs ####
mm = getMatches (archp)
mmMat = assay(mm)
colnames(mmMat) = gsub ('_.*','',colnames(mmMat))
colnames(mmMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames(mmMat))

hyper_res = lapply (unique(km$cluster), function(mod) {
  module_peaks = rowRanges(mm)[which (rowSums (mmMat[,names(km$cluster[km$cluster == mod])]) > 0)]
  q = length (queryHits (findOverlaps (module_peaks, GRanges (rownames(DAP_res_sig)))))
  n = length (getPeakSet(archp)) - length(GRanges (rownames(DAP_res_sig)))
  m = length (GRanges (rownames(DAP_res_sig)))
  k = length (module_peaks)
  #hyper_res[[mod]] = 
  print (phyper (q, m, n, k, lower.tail = FALSE, log.p = FALSE))
  })

# Barplot of module enrichment in CD8 ext peaks ####
mod_df = data.frame (module = paste0('module',unique(km$cluster)), pval = -log10(unlist(hyper_res)+1e-9))
bp = ggplot (mod_df, aes (x = module, y = pval, fill = module)) + 
geom_bar (stat = 'identity',color='grey22') + 
scale_fill_manual (values = c(module1 = 'darkred',module2='grey',module3='grey')) + 
gtheme
pdf (file.path ('Plots','TFmodule_enrichments2.pdf'), width=4,height=3)
bp
dev.off()

# Km 1 is enriched in CD8 exhausted peaks. Find which TF in it is mostly enriched ####
ext_module = 1
hyper_res_TF = lapply (names(km$cluster[km$cluster== ext_module]), function(mod) {
  module_peaks = rowRanges(mm)[which (rowSums (mmMat[,mod, drop=F]) > 0)]
  q = length (queryHits (findOverlaps (module_peaks, GRanges (rownames(DAP_res_sig)))))
  n = length(getPeakSet(archp)) - length(GRanges (rownames(DAP_res_sig)))
  m = length(GRanges (rownames(DAP_res_sig)))
  k = length(module_peaks)
  #hyper_res[[mod]] = 
  print (phyper (q, m, n, k, lower.tail = FALSE, log.p = FALSE))
  })
hyper_res_TF_df = data.frame (hyper = unlist(hyper_res_TF), TF= names(km$cluster[km$cluster== ext_module]))
hyper_res_TF_df = hyper_res_TF_df[order(hyper_res_TF_df$hyper), ]

# Barplot of enriched TFs in ext module enrichment in CD8 ext peaks ####
hyper_res_TF_df$pval= -log10(hyper_res_TF_df$hyper)#+1e-9)
hyper_res_TF_df$TF = factor (hyper_res_TF_df$TF, levels = rev(hyper_res_TF_df$TF))
top_TF = 13
hyper_res_TF_df$label = ifelse (hyper_res_TF_df$TF %in% head(hyper_res_TF_df$TF, top_TF), as.character(head(hyper_res_TF_df$TF, top_TF)), '')

bp = ggplot (head(hyper_res_TF_df,50), aes (x = TF, y = pval, fill = pval)) + 
geom_bar (stat = 'identity') + 
scale_fill_gradient(low = "white", high = "darkred") +
geom_text(aes(label = label, y = pval + .1), hjust = -0.05,na.rm = TRUE, angle=0,size=3) +
#ylim (c(0,11.5)) +
gtheme_no_rot + 
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
coord_flip()

pdf (file.path ('Plots','TF_in_module_enrichments3.pdf'), width=3,height=7)
bp
dev.off()


# Compute co-occurrence of exthausted TFs ####
TFs = names(km$cluster)[km$cluster == 1]
motifMat = getPositions (archp)
matches = getMatches (archp)
matchesMat = assay (matches)
colnames (matchesMat) = gsub ('_.*','',colnames (matchesMat))
colnames (matchesMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames (matchesMat))

matchesMat = matchesMat[,TFs]
matchesMat = matchesMat[rowSums (matchesMat) > 0,]

cooc = matrix (ncol = ncol(matchesMat), nrow= ncol(matchesMat))

for (i in 1:ncol(matchesMat))
  {
  for (z in 1:ncol(matchesMat)) 
    {
    ov = sum (rowSums (matchesMat[,c(i,z)]) == 2) / min (colSums(matchesMat[,c(i,z)]))
    cooc[i,z] = ov
    }
  }

colnames (cooc) = TFs
rownames (cooc) = TFs
diag (cooc) = 1

cooc_hm = Heatmap (
    cooc,
    column_names_rot =45, 
    row_names_gp = gpar(fontsize = 3),
    column_names_gp = gpar(fontsize = 3),
    col =palette_cooccurrence, 
    cluster_rows=F,
    cluster_columns = F,
rect_gp = gpar (col = "white", lwd = 1))

pdf (file.path ('Plots','selected_TF_cooccurence_heatmaps2.pdf'), width = 6,height=5)
cooc_hm 
dev.off()

# Compute co-occurrence of exthausted TFs in exthausted peaks ####
TFs = as.character(head(hyper_res_TF_df$TF, top_TF))
motifMat = getPositions (archp)
matches = getMatches (archp)
matches = subsetByOverlaps (matches, GRanges(rownames(DAP_res_sig)))
matchesMat = assay (matches)
colnames (matchesMat) = gsub ('_.*','',colnames (matchesMat))
colnames (matchesMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames (matchesMat))

matchesMat = matchesMat[,TFs]
matchesMat = matchesMat[rowSums (matchesMat) > 0,]
colSums (matchesMat)
cooc = matrix (ncol = ncol(matchesMat), nrow= ncol(matchesMat))

for (i in 1:ncol(matchesMat))
  {
  for (z in 1:ncol(matchesMat)) 
    {
    ov = sum (rowSums (matchesMat[,c(i,z)]) == 2) / min (colSums(matchesMat[,c(i,z)]))
    cooc[i,z] = ov
    }
  }

colnames (cooc) = TFs
rownames (cooc) = TFs
diag (cooc) = 0

cooc_hm = Heatmap (
    cooc,
    column_names_rot =45, 
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    col =palette_cooccurrence, 
    cluster_rows=F,
    cluster_columns = F,
rect_gp = gpar (col = "white", lwd = 1))

pdf (file.path ('Plots','selected_TF_cooccurence_in_ext_peaks_heatmaps2.pdf'), width = 6,height=5)
cooc_hm 
dev.off()


# Make Venn Diagram of overlap of exhauted peaks with RUNX and NR4A2 and AP1 ####
library (ggVennDiagram)
peak_ov = lapply (as.list(as.data.frame(matchesMat)), function(x) paste0('p',1:nrow(matchesMat))[x])
peak_ov = list(
  RUNX = unique(unlist(peak_ov[c('RUNX1','RUNX2')])),
  NR4A2 = peak_ov[['NR4A2']],
  AP1 = unique(unlist(peak_ov[c('FOS','JUND','JUN','BATF','JUNB','BACH1','FOSL2','FOSB','SMARCC1','BACH2')])))
vp = ggVennDiagram(peak_ov) + scale_fill_gradient(low="grey90",high = "darkred")
pdf (file.path('Plots','ext_TF_overlap_venn.pdf'),4,4)
vp
dev.off()


# # Run DAM NK KLRC1 + CD8 ext vs FGFBP2 + CD8 ####
# library (presto)
# metaGroupName = 'celltype2'

# # Subset only for expressed TFs ####
# metaGroupName = 'celltype2'
# ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat), group.by = metaGroupName)[[1]]) +1)
# min_exp = 0.1
# ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
# active_TFs = rownames(ps)[rowSums(ps) > 0]

# if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat = mMat[active_TFs,]
# mMat = as.data.frame (mMat)

# mMat_ext = mMat[,as.character(archp@cellColData[,metaGroupName]) %in% c('CD8','CD8_exhausted','NK_KLRC1','NK_FGFBP2')]
# ext_vs_ctx = as.character(archp@cellColData[,metaGroupName])[as.character(archp@cellColData[,metaGroupName]) %in% c('CD8','CD8_exhausted','NK_KLRC1','NK_FGFBP2')]
# ext_vs_ctx = ifelse (ext_vs_ctx %in% c('CD8','NK_FGFBP2'),'ctx','ext')
# res = wilcoxauc (mMat_ext, ext_vs_ctx)

# res_l = lapply (split (res, res$group), function(x){
#   tmp = x[x$logFC > 0,]
#   tmp = tmp[order (tmp$pval),]
#   tmp
# })

# tf_ext = res_l[['ext']]$feature[res_l[['ext']]$padj < 0.01]

# Use TF in module instead
tf_ext = names (km$cluster[km$cluster == 1])

# # Take highest mean of NK KLRC1 + CD8 exhausted TFs ####
metaGroupName = 'celltype2'

if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = scale(mMat[tf_ext,])
mMat = as.data.frame (mMat)

mMat_mg = as.data.frame (t(mMat))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
mMat_mg = t (mMat_mg)

ext_vs_ctx_mean = rowMeans (mMat_mg[,c('CD8_exhausted','NK_KLRC1')])
ext_vs_ctx_mean = ext_vs_ctx_mean[order (-ext_vs_ctx_mean)]

ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat_mg), group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
#ps = ps[, colnames(DAM_hm@matrix)]
ps_tf = ps[tf_ext,c('CD8_exhausted','NK_KLRC1')]
ps_tf_mean = rowMeans(ps_tf)
top_exp_tf = head(names(ps_tf_mean)[order(-ps_tf_mean)],10)
mMat_mg = mMat_mg[tf_ext,c('CD8_exhausted','NK_KLRC1')]



# Import table of cor tfs to KLRC1 and PDCD1 in nk and cd8
#nk_cd8_ext_cor = read.csv (file.path('..','scrna','nk_cd8_ext_cor.csv'), row.names=1)

cd8_bar = HeatmapAnnotation (mark = anno_mark(at = match(top_exp_tf,names(ext_vs_ctx_mean)), 
  labels_gp=gpar(fontsize = 6), 
  labels_rot = 45,
  labels = top_exp_tf,side='top'),  
  ' ' = anno_barplot (ps_tf[names(ext_vs_ctx_mean),c('CD8_exhausted')],border=F,gp = gpar(color = "white")),
#  ' ' = nk_cd8_ext_cor[rownames(ps_tf), 'cd8'],
  #which='row', 
  col=list(' ' = palette_expression_cor_fun),
  simple_anno_size = unit(.2, "cm"))
nk_bar = HeatmapAnnotation (
  #'' = nk_cd8_ext_cor[rownames(ps_tf), 'nk'],
  ' ' = anno_barplot(-ps_tf[names(ext_vs_ctx_mean),c('NK_KLRC1')],
    border=F,gp = gpar(color = "white")), 
  #which='row',
   col=list(' ' = palette_expression_cor_fun),
    simple_anno_size = unit(.2, "cm"))

TF_hm = Heatmap (t(mMat_mg[names(ext_vs_ctx_mean),]), 
          #row_labels = colnames (mMat_mg),
          #column_title = paste('top',top_genes),
          #clustering_distance_columns = 'euclidean',
          #clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          top_annotation = cd8_bar,
          bottom_annotation = nk_bar,
          row_names_side = 'left',
          #column_km = 20,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 0),
          column_names_rot = 45,
          name = 'TF activity',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_deviation_fun (mMat_mg)#,
          #width = unit(1, "cm")
          #right_annotation = motif_ha
          )

pdf (file.path ('Plots','chromvar_rna_expression_NK_CD8_EXT_heatmaps2.pdf'), width = 5,height=2.5)
draw (TF_hm)# + TF_exp_selected_hm2)
dev.off()
   


# Export table
all (rownames(mMat_mg) == rownames (ps_tf))
mat_combined = cbind(mMat_mg, ps_tf)
colnames (mat_combined) = c('CD_ext_activity','NK_KLRC1_activity','CD8_ext_RNA','NK_KLRC1_RNA')
write.csv (mat_combined, 'top_TF_CD8_NK_dual_ext_TF_activity_RNA.csv')






# Generate a scatterplot of NK KRLC1 and CD8 ext TFs ####
# Run DAM NK KLRC1 + CD8 ext vs FGFBP2 + CD8 ####

# # Differential Accessed motifs ####
# metaGroupName = "celltype2"
# force=FALSE
# source (file.path('..','..','git_repo','utils','DAM.R'))

# library (presto)
# metaGroupName = 'celltype2'

# if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat = mMat[active_TFs,]
# mMat = as.data.frame (mMat)

# mMat_ext = mMat[,as.character(archp@cellColData[,metaGroupName]) %in% c('CD8','CD8_exhausted','NK_KLRC1','NK_FGFBP2')]
# ext_vs_ctx = as.character(archp@cellColData[,metaGroupName])[as.character(archp@cellColData[,metaGroupName]) %in% c('CD8','CD8_exhausted','NK_KLRC1','NK_FGFBP2')]
# ext_vs_ctx = ifelse (ext_vs_ctx %in% c('CD8','NK_FGFBP2'),'ctx','ext')
# res = wilcoxauc (mMat_ext, ext_vs_ctx)
# res = res[res$group == 'ext',]

# ### Run same for RNA ####
# eMat = srt@assays$RNA@data
# eMat = eMat[rownames(eMat) %in% active_TFs,]
# #ps = ps[, colnames(DAM_hm@matrix)]
# eMat_ext = eMat[,srt@meta.data[,metaGroupName] %in% c('CD8','CD8_exhausted','NK_KLRC1','NK_FGFBP2')]
# ext_vs_ctx = srt@meta.data[,metaGroupName][srt@meta.data[,metaGroupName] %in% c('CD8','CD8_exhausted','NK_KLRC1','NK_FGFBP2')]
# ext_vs_ctx = ifelse (ext_vs_ctx %in% c('CD8','NK_FGFBP2'),'ctx','ext')
# eRes = wilcoxauc (eMat_ext, ext_vs_ctx)
# eres = eRes[eRes$group == 'ext',]
# eres = eres[match (res$feature, eres$feature), ]
# res_combined = data.frame (res, padj_exp = eres$padj, logFC_exp = eres$logFC, logFC_exp_sign = sign(eres$logFC))

# # volcano plot
# logfcThreshold = 0.02
# pvalAdjTrheshold = 1e-5
# res_combined$sig = 'ns'
# res_combined$sig[res_combined$padj < pvalAdjTrheshold & res_combined$logFC > logfcThreshold] = 'Ext'
# res_combined$sig[res_combined$padj < pvalAdjTrheshold & res_combined$logFC < -logfcThreshold] = 'Naive'
# res_combined$exp_sig = res_combined$sig
# res_combined$exp_sig[res_combined$exp_sig != 'ns'] = ifelse (res_combined$logFC_exp[res_combined$exp_sig != 'ns'] > 0, 'Ext','Naive')


# vp = ggplot (res_combined, aes(x=logFC, y= -log10(padj))) +
#     geom_point(shape=21, aes (fill = sig, color = exp_sig, size = -log10(padj_exp)), alpha=.5) +
#     geom_vline(xintercept = logfcThreshold, linetype="dashed", 
#                 color = "grey20", size=.5) +
#     geom_vline(xintercept = -logfcThreshold, linetype="dashed", 
#                 color = "grey20", size=.5) +
#     geom_hline(yintercept = -log10(pvalAdjTrheshold), linetype="dashed", 
#                 color = "grey20", size=.5) + 
# #    geom_text_repel (size=2, data = p11_dev_rna, aes(label = labels),segment.size=.2) + 
#     ggtitle ('Dysfunctional vs naive') +
#     #geom_label_repel (size=2,max.overlaps=10000, data = deg2_cl, aes(label = show_genes), color='red') + 
#     scale_color_manual (values = c(ns='grey77',Naive='green',Ext='red')) + 
#     scale_fill_manual (values = c(ns='grey77',Naive='green',Ext='red')) + 
#     gtheme_no_rot

# pdf (file.path ('Plots', 'dysfunctional_vs_naive_volcano2.pdf'),height=3,width=5)
# vp
# dev.off()


### Make scatterplot of sig exhausted TF comparing dev and expression in CD8 ext and NK KLRC1 ####
#ext_TF = res_combined$feature[res_combined$sig == 'Ext' & res_combined$logFC_exp > 0]

if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = mMat[tf_ext,]
mMat = as.data.frame (mMat)
mMat_mg = as.data.frame (t(mMat))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
mMat_mg = t (mMat_mg)
mMat_mg = mMat_mg[,c('CD8_exhausted','NK_KLRC1')]

ps = log2(as.data.frame (AverageExpression (srt, features = tf_ext, group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
#ps = ps[, colnames(DAM_hm@matrix)]
ps_tf = ps[tf_ext,c('CD8_exhausted','NK_KLRC1')]
colnames (ps_tf) = paste0('exp_',colnames (ps_tf))
ext_TF_df = cbind (mMat_mg, ps_tf)
ext_TF_df$exp_diff = ext_TF_df$exp_CD8_exhausted - ext_TF_df$exp_NK_KLRC1 
ext_TF_df$exp_mean = rowMeans (ext_TF_df[,c('exp_CD8_exhausted', 'exp_NK_KLRC1')])
ext_TF_df$TF_mean = rowMeans (ext_TF_df[,c('CD8_exhausted', 'NK_KLRC1')])
ext_TF_df$label = ''
ext_TF_df$label[match (head(rownames(ext_TF_df)[order(-ext_TF_df$TF_mean)],10),rownames(ext_TF_df))] = head(rownames(ext_TF_df)[order(-ext_TF_df$TF_mean)],10)
#ext_TF_df$logFC_exp = res_combined$logFC_exp[match (rownames(ext_TF_df),res_combined$feature)]
custom_palette <- c("white",'navyblue', "white")
vp = ggplot (ext_TF_df, aes(x=CD8_exhausted, y= NK_KLRC1)) +
    geom_point(shape=21, aes (fill = exp_diff, size = exp_mean)) +
    geom_text_repel (size=2, aes(label = label),segment.size=.2) + 
  scale_fill_gradientn(
    colors = custom_palette,   # Custom colors
    limits = c(-1, 1),         # Range for color scaling
    values = c(0, 0.5, 1),     # Relative positions of the colors in the palette
    oob = scales::squish       # Handle out-of-bounds values
  ) +
    # ggtitle ('Top Dysfunctional') +
    # #geom_label_repel (size=2,max.overlaps=10000, data = deg2_cl, aes(label = show_genes), color='red') + 
    # scale_color_manual (values = c(ns='grey77',Naive='green',Ext='red')) + 
    # scale_fill_manual (values = c(ns='grey77',Naive='green',Ext='red')) + 
    gtheme_no_rot

pdf (file.path ('Plots', 'top_dysfunctional_scatter3.pdf'),height=3,width=4.5)
vp
dev.off()


### Export TF high in exhausted T NK cells and supported by RNA
write.csv (ext_TF_df,'exhausted_TF_rna.csv')






# res_l = lapply (split (res, res$group), function(x){
#   tmp = x[x$logFC > 0,]
#   tmp = tmp[order (tmp$pval),]
#   tmp
# })

# tf_ext = res_l[['ext']]$feature[res_l[['ext']]$padj < 0.01]

# tf_ext = res_l[['ext']]$feature[res_l[['ext']]$padj < 0.01]

# # Take highest mean of NK KLRC1 + CD8 exhausted TFs ####
# metaGroupName = 'celltype2'

# if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat = mMat[tf_ext,]
# mMat = as.data.frame (mMat)

# mMat_mg = as.data.frame (t(mMat))
# mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
# mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
# rownames (mMat_mg) = mMat_mg[,1]
# mMat_mg = mMat_mg[,-1]
# mMat_mg = t (mMat_mg)
# mMat_mg = mMat_mg[,c('CD8_exhausted','NK_KLRC1')]
# mMat_mg = mMat_mg[tf_ext,]

# ps = log2(as.data.frame (AverageExpression (srt, features = tf_ext, group.by = metaGroupName)[[1]]) +1)
# colnames (ps) = gsub ('-','_',colnames(ps))
# #ps = ps[, colnames(DAM_hm@matrix)]
# ps_tf = ps[tf_ext,c('CD8_exhausted','NK_KLRC1')]
# ps_tf = ps_tf[tf_ext,]





### Plot deviation and expression of NR4A2 across samples ####
sample_names = c('P1','P10','P11','P12','P13','P3','P4','P5','P8') # Use only matching samples
archp$sample_celltype2 = paste0(archp$Sample,'|',archp$celltype2)
metaGroupName = 'sample_celltype2'
TF = 'IRF4'
TF = 'EOMES'
TF = 'CBFB'
TF = 'NR4A2'
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat_mg = mMat[TF, , drop=F]
mMat_mg = as.data.frame (t(mMat_mg))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = setNames(mMat_mg[,2], rownames(mMat_mg))

# Generate RNA pseudobulk of matching cell types ####
srt$sample_celltype2 = paste0(srt$sampleID,'|',srt$celltype2)
metaGroupName = 'sample_celltype2'
ps = log2(as.data.frame (AverageExpression (srt, features = TF, group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
mMat_mg = mMat_mg[sapply (names(mMat_mg), function(x) unlist(strsplit(x, '\\|'))[1] %in% sample_names)]
ps_df = data.frame (deviation = mMat_mg, expression = unlist(ps[names(mMat_mg)]))
ps_df$celltype = sapply (rownames(ps_df), function(x) unlist (strsplit (x, '\\|'))[2])
ps_df$sample = sapply (rownames(ps_df), function(x) unlist (strsplit (x, '\\|'))[1])
ps_df1 = ps_df[, c('deviation','celltype','sample')]
colnames (ps_df1)[1] = 'intensity'
ps_df1$type = 'activity'
ps_df2 = ps_df[, c('expression','celltype','sample')]
colnames (ps_df2)[1] = 'intensity'
ps_df2$type = 'expression'
ps_df3 = rbind (ps_df1, ps_df2)
ps_df3$celltype = factor (ps_df3$celltype, levels = c('NK_KLRC1','CD8_exhausted','Tregs','CD8','CD4','NK_FGFBP2'))
ps_sp = ggplot (ps_df3, aes (x = celltype, y= intensity, fill= celltype)) + 
  geom_boxplot () + 
  scale_fill_manual (values = palette_tnk_cells) + 
  facet_wrap (~type, scales = 'free_y') +
  ggtitle (TF) +
  gtheme #+ geom_text_repel (size=3.5, data = ps_df, aes(label = rownames(ps_df))) 

pdf (file.path ('Plots',paste0(TF,'_dev_exp_boxplot.pdf')), width=5,height=2.5)
ps_sp
dev.off()





# Compare NKT pseudobulks to peakset of exhausted CD8 and from human meta-analysis Riegel et al ####
if (!exists('fragments')) 
  {
  fragments = getFragmentsFromProject (archp)
  fragments = unlist(fragments)
  }

# CD8 ext DAP ####
DAP_res = readRDS ('DAP_CD8_CD8_ext_pairwise.rds')
ext_ps_gr = GRanges (rownames(DAP_res[DAP_res$FDR < .01 & DAP_res$Log2FC > 0, ]))
ext_ps_gr = extendGR (gr = ext_ps_gr, upstream = 3000, downstream = 3000)
peak_windows = slidingWindows (x = ext_ps_gr, width = 50, step = 25)

# Create peak enrichment plots with CD8 ext peaks ####
metaGroupName = 'celltype2'
force = TRUE
ext_l = list()

pb =progress::progress_bar$new(total = length (unique (as.character(archp@cellColData[,metaGroupName]))))
if (!file.exists(paste0('ext_peaks_',metaGroupName,'.rds')) | force)
  {
  for (metagroup in unique (as.character(archp@cellColData[,metaGroupName])))
    {
    pb$tick()            
    #fragments = ReadFragments(fragment_paths[sam], cutSite = FALSE)
    fragments_metagroup = fragments[fragments$RG %in% rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) == metagroup]]
    fragments_metagroup_counts = lapply (peak_windows, function(x) countOverlaps (x, fragments_metagroup))
    fragments_metagroup_counts_df = do.call (cbind, fragments_metagroup_counts)
    fragments_metagroup_counts_df = (fragments_metagroup_counts_df / sum (archp$ReadsInTSS[as.logical(archp@cellColData[,metaGroupName] == metagroup)])) * 10e6
    #colnames (fragments_metagroup_counts_df) = as.character(ext_hg38_ov)
    # write.table (fragments_metagroup_counts_df, file.path(paste0('riegel_ext_peaks_windows_',metagroup,'.tsv')), sep='\t')
    ext_l[[metagroup]] = fragments_metagroup_counts_df
    }
  saveRDS (ext_l, paste0('ext_peaks_',metaGroupName,'.rds'))
  } else {
  ext_l = readRDS (paste0('ext_peaks_',metaGroupName,'.rds'))
  }

#ext_df_long = gather (ext_df, coverage, celltype, 1:(ncol(ext_df)-1))
ext_den = lapply (ext_l, function(x) rowSums(x))
ext_den = as.data.frame (do.call (cbind, ext_den))
#ext_den = ext_den[,1,drop=F]
#ext_den = gather (as.data.frame(ext_den), celltype, coverage)
ext_den$bin = 1:260
#ext_den = split (ext_den, ext_den$celltype)
den = lapply (c('NK_FGFBP2','CD4','CD8','NK_KLRC1','CD8_exhausted','Tregs'), function(x) ggplot (ext_den[,c(x,'bin')], aes_string (x= 'bin', y = x)) + 
#geom_bar (stat= 'identity')
geom_density(color="navyblue", fill="navyblue", alpha=0.5, stat='identity') +
ylim (c(0, max(ext_den))) +
gtheme_no_rot)
#facet_wrap (~celltype, ncol=length(unique(ext_den$celltype)))

pdf (file.path ('Plots',paste0('density_coverage_',metaGroupName,'_CD8_ext.pdf')), height=4,width=12)
wrap_plots (den, ncol= length(unique(names (ext_l))))
dev.off()

median_celltype = lapply (ext_l, function(x) colMeans (x))
median_celltype_order = sapply (median_celltype, mean)
median_celltype_order = names(median_celltype_order)[order(median_celltype_order)]
ext_l2 = lapply (seq_along(ext_l), function(x) ext_l[[x]][,order(-median_celltype[[x]])])
ext_df = do.call (rbind, ext_l2)
ext_df = apply (ext_df, c(1,2), function(x) as.numeric(x))
ext_df = as.data.frame (t(log10(ext_df+1)))
colnames (ext_df) = rep (names(ext_l), each=dim(ext_l[[1]])[1])
#rownames(ext_df) = as.character(ext_hg38_sub)
#ext_df$region = rownames(ext_df)
#rownames (ext_df) = NULL
#ext_df = as.data.frame(apply(ext_df, 2, function(x) sort(x, decreasing = TRUE)))
#palette_fragments = paletteer::paletteer_c("ggthemes::Classic Area-Brown",n=40)
ha = HeatmapAnnotation (bar1 = anno_barplot(colSums(ext_df),gp = gpar(fill = "azure4",border =NA,lty='blank'),border =FALSE, baseline=200,lty='blank'))
ext_df2 = ext_df

pdf (file.path ('Plots','CD8_ext_peakset_enrichments.pdf'), height=4,width=6)
Heatmap (ext_df,
  top_annotation = ha,
#  column_split = , 
  column_split=factor(colnames(ext_df), levels=median_celltype_order),
  cluster_rows=F,
  column_title_gp = gpar(
fontsize = 8),
  column_names_gp = gpar(
fontsize = 0),
  row_names_gp = gpar(
fontsize = 0),
  col = rev(palette_fragments), 
  cluster_columns=F,
  border=T)
dev.off()


# Riegel CD8 ext peaks ####
library (liftOver)
ext_ps = read.csv ('riegel_Text_list.csv')
ext_clusters = c('C1','C2','C3','C4')
ext_ps = ext_ps[ext_ps$group_name %in% ext_clusters,]
ext_ps_gr = GRanges (ext_ps)
if(!file.exists('hg19ToHg38.over.chain'))
  {
  download.file("https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz", "hg19ToHg38.over.chain.gz")
  system("gzip -d hg19ToHg38.over.chain.gz")
  #system (paste0('wget (https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg38.over.chain.gz)', '-P', getwd()))
  }
ch = import.chain ('hg19ToHg38.over.chain')
#seqlevelsStyle(hub_hg19) = "UCSC"  # necessary
ext_hg38 = unlist (liftOver(ext_ps_gr, ch))
ext_hg38_sub = head (ext_hg38[ext_hg38$Log2FC > 3, ],1000)
ext_hg38_sub = resize (ext_hg38_sub, 1, "center")
ext_hg38_sub = extendGR (gr = ext_hg38_sub, upstream = 3000, downstream = 3000)
peak_windows_riegel = slidingWindows (x = ext_hg38_sub, width = 50, step = 25)
# peak_windows = peak_windows[sapply(peak_windows, length) == 50]
# peak_windows_l =  lapply (1:length(peak_windows[[1]]), function(x) unlist(lapply (peak_windows, function(y) y[x])))

# Create peak enrichment plots with Riegel peaks ####
metaGroupName = 'celltype2'
force = FALSE
ext_l = list()

pb =progress::progress_bar$new(total = length (unique (as.character(archp@cellColData[,metaGroupName]))))
if (!file.exists(paste0('riegel_ext_peaks_',metaGroupName,'.rds')) | force)
  {
  for (metagroup in unique (as.character(archp@cellColData[,metaGroupName])))
    {
    pb$tick()            
    #fragments = ReadFragments(fragment_paths[sam], cutSite = FALSE)
    fragments_metagroup = fragments[fragments$RG %in% rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) == metagroup]]
    fragments_metagroup_counts = lapply (peak_windows, function(x) countOverlaps (x, fragments_metagroup))
    fragments_metagroup_counts_df = do.call (cbind, fragments_metagroup_counts)
    fragments_metagroup_counts_df = (fragments_metagroup_counts_df / sum (archp$ReadsInTSS[as.logical(archp@cellColData[,metaGroupName] == metagroup)])) * 10e6
    #colnames (fragments_metagroup_counts_df) = as.character(ext_hg38_ov)
    # write.table (fragments_metagroup_counts_df, file.path(paste0('riegel_ext_peaks_windows_',metagroup,'.tsv')), sep='\t')
    ext_l[[metagroup]] = fragments_metagroup_counts_df
    }
  saveRDS (ext_l, paste0('riegel_ext_peaks_',metaGroupName,'.rds'))
  } else {
  ext_l = readRDS (paste0('riegel_ext_peaks_',metaGroupName,'.rds'))
  }

#ext_df_long = gather (ext_df, coverage, celltype, 1:(ncol(ext_df)-1))
ext_den = lapply (ext_l, function(x) rowSums(x))
ext_den = as.data.frame (do.call (cbind, ext_den))
#ext_den = ext_den[,1,drop=F]
#ext_den = gather (as.data.frame(ext_den), celltype, coverage)
ext_den$bin = 1:240
#ext_den = split (ext_den, ext_den$celltype)
den = lapply (c('NK_FGFBP2','CD4','CD8','NK_KLRC1','CD8_exhausted','Tregs'), function(x) ggplot (ext_den[,c(x,'bin')], aes_string (x= 'bin', y = x)) + 
#geom_bar (stat= 'identity')
geom_density(color="navyblue", fill="navyblue", alpha=0.5, stat='identity') +
ylim (c(0, max(ext_den))) +
gtheme_no_rot)
#facet_wrap (~celltype, ncol=length(unique(ext_den$celltype)))

pdf (file.path ('Plots',paste0('density_coverage_',metaGroupName,'2.pdf')), height=4,width=12)
wrap_plots (den, ncol= length(unique(names (ext_l))))
dev.off()

median_celltype = lapply (ext_l, function(x) colMeans (x))
median_celltype_order = sapply (median_celltype, mean)
median_celltype_order = names(median_celltype_order)[order(median_celltype_order)]
ext_l2 = lapply (seq_along(ext_l), function(x) ext_l[[x]][,order(-median_celltype[[x]])])
ext_df = do.call (rbind, ext_l2)
ext_df = apply (ext_df,c(1,2), function(x) as.numeric(x))
ext_df = as.data.frame (t(log10(ext_df+1)))
colnames (ext_df) = rep (names(ext_l), each=240)
rownames(ext_df) = as.character(ext_hg38_sub)
#ext_df$region = rownames(ext_df)
#rownames (ext_df) = NULL
#ext_df = as.data.frame(apply(ext_df, 2, function(x) sort(x, decreasing = TRUE)))
#palette_fragments = paletteer::paletteer_c("ggthemes::Classic Area-Brown",n=40)
ha = HeatmapAnnotation (bar1 = anno_barplot(colSums(ext_df),gp = gpar(fill = "azure4",border =NA,lty='blank'),border =FALSE, baseline=200,lty='blank'))
pdf (file.path ('Plots','Riegel_ext_peakset_enrichments.pdf'), height=4,width=6)
Heatmap (ext_df,
  top_annotation = ha,
#  column_split = , 
  column_split=factor(rep (names(ext_l), each=240), levels=median_celltype_order),
  cluster_rows=F,
  column_title_gp = gpar(
fontsize = 8),
  column_names_gp = gpar(
fontsize = 0),
  row_names_gp = gpar(
fontsize = 0),
  col = rev(palette_fragments), 
  cluster_columns=F,
  border=T)
dev.off()







### Find all peaks around exhausted TFs and correlate accessibility and RNA expression across celltype pseudobulks ####
metaGroupName = 'celltype2'
TF = read.csv ('exhausted_TF_rna.csv')
TF = c('NR4A2','RUNX2') # include only those that show balanced expression between NK KLRC1 and CD8 ext
# Get all peaks correlated with exhausted TF
library (org.Hs.eg.db)
gene_regions = genes (TxDb.Hsapiens.UCSC.hg38.knownGene)
eg_sym = toTable (org.Hs.egSYMBOL)
gene_regions$symbol = eg_sym$symbol[match(gene_regions$gene_id, eg_sym$gene_id)]
gene_regions = gene_regions[gene_regions$symbol %in% TF ]

extend_region = 250000
gene_regions_extended = GRangesList (lapply (seq_along(gene_regions), function(x) extendGR (gr = gene_regions[x], upstream = extend_region, downstream = extend_region)))

gene_regions_peaks = lapply (seq_along(gene_regions), function(x) queryHits (findOverlaps(getPeakSet(archp), gene_regions_extended[x])))
gene_regions_peaks = lapply (gene_regions_peaks, function(x) getPeakSet(archp)[x])

force = FALSE
if (!file.exists ('pMats.rds') | force)
  {
  pMats = lapply (gene_regions_peaks, function(x) getGroupSE(
    ArchRProj = archp,
    useMatrix = 'PeakMatrix',
    groupBy = metaGroupName,
    divideN = TRUE,
    scaleTo = NULL,
    threads = getArchRThreads(),
    verbose = TRUE,
    logFile = createLogFile("getGroupSE")
  ))
  pMats2 = lapply (seq_along(gene_regions_peaks), function(x) pMats[[x]][queryHits(findOverlaps(GRanges(rowData(pMats[[x]])), gene_regions_peaks[[x]]))])
  pMats2 = lapply (pMats2, function(x) {tmp = assay(x); rownames(tmp) = as.character(GRanges (rowData(x))); tmp})
  names (pMats2) = gene_regions$symbol
  saveRDS (pMats2, 'pMats2.rds')
  } else {
  pMats2 = readRDS ('pMats2.rds')  
  }


ps = log2(as.data.frame (AverageExpression (srt, features = gene_regions$symbol, group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))


# Compute distance between TSS and correlated peaks ####
library (rtracklayer)
library (AnnotationDbi)
gene_ids = toTable (org.Hs.egSYMBOL)$gene_id[match (TF,toTable (org.Hs.egSYMBOL)$symbol)]
mygenes.transcripts = subset (genes(TxDb.Hsapiens.UCSC.hg38.knownGene, columns=c("tx_id", "tx_name","gene_id")), gene_id %in% gene_ids)
mygenes.tss = resize (mygenes.transcripts, width=1, fix='start')

mygenes.tss$symbol = toTable (org.Hs.egSYMBOL)$symbol[match (mygenes.tss$gene_id ,toTable (org.Hs.egSYMBOL)$gene_id)]

peak_gene_cor = lapply (gene_regions$symbol, function(x) cor (t(pMats2[[x]]), t(ps[x, colnames(pMats2[[x]])])))
peak_gene_cor = lapply (peak_gene_cor, function(x) x[order(-x[,1]),,drop=F])
names (peak_gene_cor) = gene_regions$symbol

top_e = 5
peak_gene_cor = lapply (names(peak_gene_cor), function(x) {
  peak_gene_cor[[x]] = as.data.frame (peak_gene_cor[[x]])
  gr = GRanges (rownames(peak_gene_cor[[x]]))
  dis = distanceToNearest (gr,mygenes.tss[mygenes.tss$symbol == x])@elementMetadata$distance
  pr = follow (gr,mygenes.tss[mygenes.tss$symbol == x])
  pr[is.na(pr)] = 0
  dis = ifelse (pr == 1, dis * -1, dis) 
  peak_gene_cor[[x]]$distance = dis
  peak_gene_cor[[x]]$name = ''
  peak_gene_cor[[x]]$name[1:top_e] = head (paste0('E',peak_gene_cor[[x]]$distance),top_e)
  peak_gene_cor[[x]]$max = apply(pMats2[[x]],1, max)[match(rownames(peak_gene_cor[[x]]),rownames(pMats2[[x]]))]
  colnames(peak_gene_cor[[x]]) = c('cor','distance','name','max')
  peak_gene_cor[[x]]
  })

names (peak_gene_cor) = gene_regions$symbol
vp = lapply (names(peak_gene_cor), function(x)
    {
    ggplot (peak_gene_cor[[x]], aes(x=distance, y= cor)) +
    geom_bar (stat = 'identity') +
    geom_point (aes (size = max), shape=21, color='white',alpha=.5, fill = 'red') +
    xlim(c(-extend_region,extend_region)) + 
    ggtitle (x) + 
    geom_label_repel (size=3, color='grey22', aes(label = name), segment.size=.2) + gtheme_no_rot
    })

pdf (file.path ('Plots','enhancers_distance_ext_TF3.pdf'),width=7,height=3)
vp
dev.off()



# # ### Use P2G analysis and cNMF from RNA to identify active TF via regulons  ####
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
    
### Plot NR4A2 region ####
archp$celltype3 =  archp$celltype2
archp$celltype3[archp$celltype3 == 'CD4'] =  c('C1_CD4')
archp$celltype3[archp$celltype3 == 'CD8'] =  c('C2_CD8')
archp$celltype3[archp$celltype3 == 'NK_FGFBP2'] =  c('C3_NK_FGFBP2')
archp$celltype3[archp$celltype3 == 'Tregs'] =  c('C4_Tregs')
archp$celltype3[archp$celltype3 == 'CD8_exhausted'] =  c('C5_CD8_exhausted')
archp$celltype3[archp$celltype3 == 'NK_KLRC1'] =  c('C6_NK_KLRC1')
palette_tnk_cells_ext2 = palette_tnk_cells
names (palette_tnk_cells_ext2) = c(
  'C2_CD8',
  'C1_CD4',
  'C4_Tregs',
  'C6_NK_KLRC1',
  'C3_NK_FGFBP2',
  'C5_CD8_exhausted')


### Hubs analysis #####
metaGroupName = "Clusters_H"
cor_cutoff = 0.2
#max_dist = 12500
max_dist = 12500
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
hubs_obj = readRDS (file.path(hubs_dir,'global_hubs_obj.rds'))  

metaGroupName = 'celltype3'
TF = 'NR4A2'
pdf()
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp, 
    sizes = c(6, 1, 1, 1,1,1),
    groupBy = metaGroupName, 
    geneSymbol = TF,
    normMethod = "ReadsInTSS",
    scCellsMax=3000,
    plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
    hubs_regions = hubs_obj$hubsCollapsed,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 50000,
    pal = palette_tnk_cells_ext2,
    #ylim=c(0,0.1),
    downstream = 300000,
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
  width=5,height=3, 
  name =paste0(paste(TF, collapse='_'),'_coveragePlots.pdf'),
  addDOC = F)

metaGroupName = 'celltype2'
top_dah = data.frame (
gene = srt@assays$RNA@data[TF,],
group = srt@meta.data[,metaGroupName])
top_dah$group = factor (top_dah$group, levels = 
  rev (c('CD4','CD8','NK_FGFBP2','Tregs','CD8_exhausted',
  'NK_KLRC1','TFH')))
top_dah = na.omit(top_dah)
bp = ggplot (top_dah, aes (x = gene, y = group, fill = group)) + 
vlp + 
bxpv + 
scale_fill_manual (values = palette_tnk_cells) +
#geom_point (position='identity', alpha=.3, color="grey44", size=1) +
gtheme_no_rot

pdf (file.path ('Plots', paste0('scrna_',TF,'_boxplots.pdf')), height=4, width=4)
bp
dev.off()


# Run chromBPnet 
# Variables
chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet'
repodir='/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
grefdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet'
celltype='CD8'
fold_numbers = c(0,1,2,3,4)


### Check expression of RUNX3 and CTCF predicted by chromBPnet
pdf (file.path ('Plots','TF_exp_dotplots.pdf'), width=5, height=3)
DotPlot (srt[,srt$celltype2 != 'Proliferating'], 
  features= c('RUNX1','RUNX2','RUNX3','NR4A2','NR4A1','NR4A3','EGR2'), group.by = 'celltype2') + gtheme
dev.off()

deg_res = FindMarkers (srt, ident.1 = 'NK_KLRC1', ident.2 = 'NK_FGFBP2', group.by='celltype2')
deg_res[rownames(deg_res) %in% c('RUNX1','RUNX2','RUNX3'),]
head (deg_res,10)







metaGroupName = 'celltype3'
TF = c('CTLA4','PDCD1','ICOS','HAVCR2')
pdf()
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp, 
    sizes = c(6, 1, 1, 1,1,1),
    groupBy = metaGroupName, 
    geneSymbol = TF,
    normMethod = "ReadsInTSS",
    scCellsMax=3000,
    loop_size = .2,
        plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
        hubs_regions = hubs_obj$hubsCollapsed,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 100000,
    pal = palette_tnk_cells_ext2,
    #ylim=c(0,0.1),
    downstream = 100000,
    #loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
#    loops = getCoAccessibility (archp, corCutOff = 0.25),
    #  returnLoops = TRUE),
    useGroups= NULL,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2,returnLoops = TRUE),
    hubs = hubs_obj$peakLinks2
)
dev.off()
plotPDF (meso_markers, ArchRProj = archp, 
  width=5,height=3, 
  name =paste0(paste(TF, collapse='_'),'_coveragePlots.pdf'),
  addDOC = F)

  
TF='ICOS'
metaGroupName = 'celltype2'
top_dah = data.frame (
gene = srt@assays$RNA@data[TF,],
group = srt@meta.data[,metaGroupName])
top_dah$group = factor (top_dah$group, levels = 
  rev (c('CD4','CD8','NK_FGFBP2','Tregs','CD8_exhausted',
  'NK_KLRC1','TFH')))
top_dah = na.omit(top_dah)
bp = ggplot (top_dah, aes (x = gene, y = group, fill = group)) + 
vlp + 
bxpv + 
scale_fill_manual (values = palette_tnk_cells) +
#geom_point (position='identity', alpha=.3, color="grey44", size=1) +
gtheme_no_rot

pdf (file.path ('Plots', paste0('scrna_',TF,'_boxplots.pdf')), height=4, width=4)
bp
dev.off()






# Check footprint of RUNX and NR4A2 across celltypes ####
metaGroupName='celltype2'
archp <- addGroupCoverages (ArchRProj = archp, groupBy = metaGroupName)
motifPositions <- getPositions (archp)

motifs <- c('RUNX1','RUNX2','RUNX3','NR4A2')
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

seFoot <- getFootprints(
  ArchRProj = archp, 
  #positions = motifPositions_sample[markerMotifs], 
  positions = motifPositions[markerMotifs], 
  groupBy = metaGroupName
)
  
plotFootprints(
seFoot = seFoot,
ArchRProj = archp, 
normMethod = "Subtract",
plotName = "Footprints-Subtract-Bias_",
addDOC = FALSE, height=7.5, width=5,
pal = palette_tnk_cells,
smoothWindow = 25)
  





# # Compare similarity of celltypes ATAC vs RNA ####
# # # Get deviation matrix ####
# if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name

# # Subset only for expressed TFs ####
# metaGroupName = 'celltype2'
# ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat), group.by = metaGroupName)[[1]]) +1)
# min_exp = 0.1
# ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
# active_TFs = rownames(ps)[rowSums(ps) > 0]
# mMat = mMat[active_TFs,]
# mMat = as.data.frame (t(mMat))
# mMat$metaGroup = as.character (archp@cellColData[,metaGroupName])
# mMat = aggregate (.~ metaGroup, mMat, mean)
# rownames (mMat) = mMat[,1]
# mMat = mMat[,-1]
          

# mMat_cor = cor (t(mMat))
# hm_dev = draw (Heatmap (mMat_cor, 
#   rect_gp = gpar(type = "none"),
#   clustering_distance_rows='pearson' ,
#   clustering_distance_columns = 'pearson', 
#   col=palette_deviation_cor_fun, 
#   border=F,
#   ,
#   cell_fun = function(j, i, x, y, w, h, fill) {
#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#         }}))

# pdf (file.path('Plots','celltypes_similarity_triangle_deviations_heatmap.pdf'),width = 4,height=3.2)
# hm_dev = draw (Heatmap (mMat_cor, 
#   rect_gp = gpar(type = "none"),
#   clustering_distance_rows='pearson' ,
#   clustering_distance_columns = 'pearson', 
#   col=palette_deviation_cor_fun, 
#   border=F,
#   ,
#   cell_fun = function(j, i, x, y, w, h, fill) {
#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#         }}))

# print (hm_dev)
# dev.off()

# # correlate cell types using selected TF in scrna ####
# metaGroupName = 'celltype2'
# nfeats = active_TFs
# ps = log2(as.data.frame (AverageExpression (srt, features = nfeats, group.by = metaGroupName)[[1]]) +1)
# colnames (ps) = gsub ('-','_',colnames(ps))
# ps = ps[,rownames(mMat)]
# ps_cor = cor (ps)
# ps_cor_scaled = cor (t(scale(t(ps))))

# ps_cor_scaled = ps_cor_scaled[row_order(hm_dev),row_order(hm_dev)]
# ps_cor_scaled = ps_cor_scaled[rev(seq(nrow(ps_cor_scaled))),rev(seq(nrow(ps_cor_scaled)))]
# # triangle heatmap
# #ps_cor = ps_cor[rownames(ps_cor) != 'Proliferating', colnames(ps_cor) != 'Proliferating']

# pdf (file.path('Plots','celltypes_similarity_triangle_expression_heatmap.pdf'),width = 3.7,height=3)
# hm = draw (Heatmap (ps_cor_scaled, 
#   rect_gp = gpar(type = "none"),
#   cluster_columns = F,
#   cluster_rows = F,
#   #clustering_distance_rows='pearson' ,
#   #clustering_distance_columns = 'pearson', 
#   col=palette_expression_cor_fun, 
#   border=F,
#   ,
#   cell_fun = function(j, i, x, y, w, h, fill) {
#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#         }}))

# ps_cor = ps_cor[row_order(hm_dev),row_order(hm_dev)]
# ps_cor = ps_cor[rev(seq(nrow(ps_cor))),rev(seq(nrow(ps_cor)))]
# # triangle heatmap
# #ps_cor = ps_cor[rownames(ps_cor) != 'Proliferating', colnames(ps_cor) != 'Proliferating']
# hm2 = draw (Heatmap (ps_cor, 
#   rect_gp = gpar(type = "none"),
#   cluster_columns = F,
#   cluster_rows = F,
#   #clustering_distance_rows='pearson' ,
#   #clustering_distance_columns = 'pearson', 
#   col=palette_expression_correlation, 
#   border=F,
#   ,
#   cell_fun = function(j, i, x, y, w, h, fill) {
#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#         }}))

# print (hm)
# print (hm2)
# dev.off()


  
# # correlate cell types using variable genes in scrna ####
# metaGroupName = 'celltype2'
# ps = log2(as.data.frame (AverageExpression (srt, features = VariableFeatures(srt), group.by = metaGroupName)[[1]]) +1)
# colnames (ps) = gsub ('-','_',colnames(ps))
# ps = ps[,rownames(mMat)]
# ps_cor = cor (ps)
# ps_cor_scaled = cor (t(scale(t(ps))))
# #ps_cor = cor (ps)

# ps_cor_scaled = ps_cor_scaled[row_order(hm_dev),row_order(hm_dev)]
# ps_cor_scaled = ps_cor_scaled[rev(seq(nrow(ps_cor_scaled))),rev(seq(nrow(ps_cor_scaled)))]
# # triangle heatmap
# #ps_cor = ps_cor[rownames(ps_cor) != 'Proliferating', colnames(ps_cor) != 'Proliferating']
# hm = draw (Heatmap (ps_cor_scaled, 
#   rect_gp = gpar(type = "none"),
#   cluster_columns = F,
#   cluster_rows = F,
#   #clustering_distance_rows='pearson' ,
#   #clustering_distance_columns = 'pearson', 
#   col=palette_expression_cor_fun, 
#   border=F,
#   ,
#   cell_fun = function(j, i, x, y, w, h, fill) {
#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#         }}))

# ps_cor = ps_cor[row_order(hm_dev),row_order(hm_dev)]
# ps_cor = ps_cor[rev(seq(nrow(ps_cor))),rev(seq(nrow(ps_cor)))]
# # triangle heatmap
# #ps_cor = ps_cor[rownames(ps_cor) != 'Proliferating', colnames(ps_cor) != 'Proliferating']
# hm2 = draw (Heatmap (ps_cor, 
#   rect_gp = gpar(type = "none"),
#   cluster_columns = F,
#   cluster_rows = F,
#   #clustering_distance_rows='pearson' ,
#   #clustering_distance_columns = 'pearson', 
#   col=palette_expression_correlation, 
#   border=F,
#   ,
#   cell_fun = function(j, i, x, y, w, h, fill) {
#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#         }}))
# pdf (file.path('Plots','celltypes_similarity_triangle_variable_features_expression_heatmap.pdf'),width = 3.7,height=3)
# hm = draw (Heatmap (ps_cor_scaled, 
#   rect_gp = gpar(type = "none"),
#   cluster_columns = F,
#   cluster_rows = F,
#   #clustering_distance_rows='pearson' ,
#   #clustering_distance_columns = 'pearson', 
#   col=palette_expression_cor_fun, 
#   border=F,
#   ,
#   cell_fun = function(j, i, x, y, w, h, fill) {
#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#         }}))

# ps_cor = ps_cor[row_order(hm_dev),row_order(hm_dev)]
# ps_cor = ps_cor[rev(seq(nrow(ps_cor))),rev(seq(nrow(ps_cor)))]
# # triangle heatmap
# #ps_cor = ps_cor[rownames(ps_cor) != 'Proliferating', colnames(ps_cor) != 'Proliferating']
# hm2 = draw (Heatmap (ps_cor, 
#   rect_gp = gpar(type = "none"),
#   cluster_columns = F,
#   cluster_rows = F,
#   #clustering_distance_rows='pearson' ,
#   #clustering_distance_columns = 'pearson', 
#   col=palette_expression_correlation, 
#   border=F,
#   ,
#   cell_fun = function(j, i, x, y, w, h, fill) {
#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#         }}))
# print (hm)
# print (hm2)
# dev.off()



#### Import Hi-ChiP data from GSE168881 ####
library (HiCExperiment)
library (HiTC)

samples = c('SRR13961069','SRR13961070','SRR13961071','SRR13961072','SRR13961073','SRR13961074')

hic = list()
for (sam in samples)
  {
  matrix_files = paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/HiC_Texhaustion_GSE168881_hg38/hic_results/matrix/',sam,'/raw/10000/',sam,'_10000.matrix') # raw
  #matrix_files = paste0('/sc/arion/scratch/giottb01/HiCPro_matrices/hic_results/matrix/',sam,'/iced/10000/',sam,'_10000_iced.matrix') # iced normalized
  bed_files = paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/HiC_Texhaustion_GSE168881_hg38/hic_results/matrix/',sam,'/raw/10000/',sam,'_10000_abs.bed')
  bed = read.table (bed_files)
  
  hic[[sam]]<-importC(matrix_files,
               bed_files)
  }

E14exp = list()
E14norm.binned = list()

for (sam in samples)
  {
  # Focus on chr2:156312347−156492348
  
  #chr2 156292347−156642348
  E14subset = extractRegion (hic[[sam]]$chr2chr2, c(1,2),
  chr="chr2", from=156292347, to=156642348)
  ## Binning of 5C interaction map
  E14subset.binned <- binningC (E14subset, binsize=10000, method="median", step=3)
  pdf()
  E14exp[[sam]] <- getExpectedCounts (E14subset.binned, method="loess", stdev=TRUE, plot=TRUE)
  dev.off()
  E14norm.binned[[sam]] <- normPerExpected (E14subset.binned, method="loess", stdev=TRUE)
  }

projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR'

#E14norm.binned = E14exp

pdf (file.path (projdir,'Plots','HiC_map_iced_normalized2.pdf'),height=3,width=4)
lapply (E14norm.binned, function(x) mapC(x,
#mapC(E14norm.binned[[samples[1]]] - E14norm.binned[[samples[6]]],
#tracks=list(RefSeqGene=gene, CTCF=ctcf),
maxrange=10))
dev.off()



samples2 = c('SRR13961069','SRR13961074')
for (sam in samples2)
{
# Focus on chr2:156312347−156492348

#chr2 156292347−156642348
E14subset = extractRegion(hic[['SRR13961069']]$chr2chr2, c(1,2),
chr="chr2", from=156292347, to=156642348)

E14subset2 = extractRegion(hic[['SRR13961074']]$chr2chr2, c(1,2),
chr="chr2", from=156292347, to=156642348)

## Binning of 5C interaction map
E14subset.binned <- binningC (E14subset, binsize=10000, method="median", step=3)
E14subset2.binned <- binningC (E14subset2, binsize=10000, method="median", step=3)

E14norm.binned <- normPerExpected (E14subset.binned, method="loess", stdev=TRUE)
E14norm2.binned <- normPerExpected (E14subset2.binned, method="loess", stdev=TRUE)

E14norm3.binned = E14norm.binned
E14norm.binned@intdata[is.na(E14norm.binned@intdata)] = 0
E14norm2.binned@intdata[is.na(E14norm2.binned@intdata)] = 0
E14norm3.binned@intdata = E14norm2.binned@intdata - E14norm.binned@intdata

pdf (file.path (projdir,'Plots','HiC_map_diff.pdf'),height=3,width=4)
mapC (E14norm3.binned, maxrange=10)
#mapC(E14norm.binned[[samples[1]]] - E14norm.binned[[samples[6]]],
#tracks=list(RefSeqGene=gene, CTCF=ctcf),

dev.off()

