conda activate meso_scatac
#use UGER
R


set.seed(1234)

projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/Endothelial/scatac_ArchR'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

####### ANALYSIS of endothelial cells #######
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


# Load RNA
srt = readRDS (file.path('..','..','stroma','scrna','srt.rds'))
srt = srt[, srt$celltype == 'Endothelial']

#sarc_order = read.csv ('../scrna/cnmf20_sarcomatoid_sample_order.csv', row.names=1)
archp = loadArchRProject (projdir)
#devtools::load_all('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/ArchR_fixed_branch/ArchR/')

  #archp = archp[!is.na(archp$celltype)]

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

  

  pdf (file.path('Plots','celltype_umap.pdf'),5,5)
  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample",labelMeans =F,
     embedding = "UMAP_H", pal = palette_sample)
  # umap_p4 = plotEmbedding (ArchRProj = archp, 
  #   colorBy = "cellColData", name = "celltype",labelMeans =F,
  #    embedding = "UMAP_H", pal = palette_stroma)
  umap_p5 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters_H",
     embedding = "UMAP_H")
  
  print (umap_p3)
  print (umap_p5)
  dev.off()
  
  archp = saveArchRProject (archp, load = T, dropCells=T)
  
  } else {
  archp = loadArchRProject (projdir)
  }

### Call peaks on celltypes ####
metaGroupName = 'Clusters_H'
force=TRUE
peak_reproducibility='2'
pdf() # This is necessary cause cairo throws error and stops the script
if(!all(file.exists(file.path('PeakCalls', unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds'))) | force) 
source (file.path('..','..','git_repo','utils','callPeaks.R'))
dev.off()


### chromVAR analysis ####
force=TRUE
if (!all(file.exists(file.path('Annotations',
  c('Motif-Matches-In-Peaks.rds',
    'Motif-Positions-In-Peaks.rds',
    'Motif-In-Peaks-Summary.rds')))) | force)
source (file.path ('..','..','git_repo','utils','chromVAR.R'))
  


# Read in peak files from scATAC studies ####
projects = c('rawlins_fetal_lung')
projects_peaks = lapply (seq_along(projects), function(x) {
  bed_files = list.files (file.path('..','..','tumor_compartment','all_tissues_ArchR',projects[x],'PeakCalls'), pattern = '.rds')
  grlist = lapply (seq_along(bed_files), 
    function(y) readRDS (file.path('..','..','tumor_compartment','all_tissues_ArchR',projects[x],'PeakCalls',bed_files[y])))
names (grlist) = paste0(projects[x], '_', sapply (bed_files, function(z) unlist(strsplit (z, '-'))[1]))
grlist
})
projects_peaks = unlist (projects_peaks, recursive=F)

#### chromVAR analysis ####
archp = addBgdPeaks (archp, force= TRUE)
archp = addPeakAnnotations (ArchRProj = archp, 
     regions = projects_peaks, name = "scATAC_datasets", force= TRUE)

archp = addDeviationsMatrix (
  ArchRProj = archp, 
  peakAnnotation = "scATAC_datasets",
  force = TRUE
)

projects = c('Tsankov_lung')
projects_peaks = lapply (seq_along(projects), function(x) {
  bed_files = list.files (file.path('..','..','tumor_compartment','all_tissues_ArchR',projects[x],'PeakCalls'), pattern = '.rds')
  grlist = lapply (seq_along(bed_files), 
    function(y) readRDS (file.path('..','..','tumor_compartment','all_tissues_ArchR',projects[x],'PeakCalls',bed_files[y])))
names (grlist) = paste0(projects[x], '_', sapply (bed_files, function(z) unlist(strsplit (z, '-'))[1]))
grlist
})
projects_peaks = unlist (projects_peaks, recursive=F)

#### chromVAR analysis ####
archp = addBgdPeaks (archp, force= TRUE)
archp = addPeakAnnotations (ArchRProj = archp, 
     regions = projects_peaks, name = "scATAC_normal2",force=T)

archp = addDeviationsMatrix (
  ArchRProj = archp, 
  peakAnnotation = "scATAC_normal2",
  force = TRUE
)
archp=saveArchRProject(archp)


# # Differential Accessed motifs ####
# metaGroupName = "celltype"
# force=FALSE
# source (file.path('..','..','git_repo','utils','DAM.R'))

# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat_mg = mMat[active_DAM, ]
# mMat_mg = as.data.frame (t(mMat_mg))
# mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
# mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
# rownames (mMat_mg) = mMat_mg[,1]
# mMat_mg = mMat_mg[,-1]
# mMat_mg = mMat_mg[names (DAM_list),]

# # Generate RNA pseudobulk of matching cell types ####
# metaGroupName = 'celltype'
# #selected_TF = c(rownames(DAM_hm@matrix), 'NR4A3','NR4A2','NR4A1')
# ps = log2(as.data.frame (AverageExpression (srt, features = active_DAM, group.by = metaGroupName)[[1]]) +1)
# colnames (ps) = gsub ('-','_',colnames(ps))
# ps = ps[, colnames(DAM_hm@matrix)]
# ps_tf = ps[active_DAM,]

  
#  DAM_hm = Heatmap (t(scale(mMat_mg)), 
#           row_labels = colnames (mMat_mg),
#           column_title = paste('top',top_genes),
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
#         name = 'expression (scaled)',
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

# pdf (file.path ('Plots','DAM_with_rna_expression_heatmaps.pdf'), width = 8,height=4)
# draw (DAM_hm + TF_exp_selected_hm + TF_exp_selected_hm2)
# dev.off()
   

# ### Co-expression of TFs across cells #### 

# # # Get deviation matrix ####
# if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name

# # Subset only for expressed TFs ####
# metaGroupName = 'celltype'
# ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat), group.by = metaGroupName)[[1]]) +1)
# min_exp = 0.1
# ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
# active_TFs = rownames(ps)[rowSums(ps) > 0]
# #positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0]
# metaGroupName = 'celltype'
# mMat = mMat[active_TFs,]
# mMat = as.data.frame (t(mMat))
# mMat$metaGroup = as.character (archp@cellColData[,metaGroupName])
# mMat = aggregate (.~ metaGroup, mMat, mean)
# rownames (mMat) = mMat[,1]
# mMat = mMat[,-1]

# mMat_cor = cor (as.matrix(mMat), method = 'pearson')
# km = kmeans (mMat_cor, centers=5)

# pdf (file.path ('Plots','TF_modules_heatmap.pdf'), width = 4,height=3)
# cor_mMat_hm = draw (Heatmap (mMat_cor,# row_km=15,
#   #left_annotation = ha,
#   #rect_gp = gpar(type = "none"),
#   clustering_distance_rows='euclidean' ,
#   clustering_distance_columns = 'euclidean', 
#   col=palette_deviation_cor_fun, 
#   row_split = km$cluster,
#   column_split = km$cluster,
#   #row_km=2, 
#   #column_km=2,
# #  right_annotation = ha,
#   border=T,
# #   ,
#   row_names_gp = gpar(fontsize = 0), 
#   column_names_gp = gpar(fontsize = 0)
# # cell_fun = function(j, i, x, y, w, h, fill) {# THIS DOESNT WORK NEED TO USE LAYER_FUN
# #         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
# #             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
# #         }}
#   ))
#   # ,
#   # cell_fun = function(j, i, x, y, w, h, fill) {
#   #       if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#   #           grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
# #        }}))
# dev.off()

# pdf (file.path ('Plots','TF_modules_heatmap.pdf'), width = 4,height=3)
# cor_mMat_hm
# dev.off()

# if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name

# tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
# names (tf_modules) = paste0('mod_',unique(km$cluster))
# tf_modules = do.call (cbind, tf_modules)
# archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% paste0('mod_',unique(km$cluster))]
# archp@cellColData = cbind (archp@cellColData, tf_modules) 

# pdf ()
# TF_p = plotEmbedding (
#     ArchRProj = archp,
#     colorBy = "cellColData",
#     name = paste0('mod_',unique(km$cluster)), 
#     pal = rev(palette_deviation),
#     #useSeqnames='z',
#     embedding = "UMAP")
# dev.off()
# pdf (file.path ('Plots','TF_modules_umap.pdf'), width = 20,height=6)
# wrap_plots (TF_p, ncol=5)
# dev.off()

# #colnames(TF_hm@ht_list$chromVAR@matrix)[unlist(column_order(TF_hm)[c('2','3','4','5')])]
# #which (colnames(mMat) == 'NR4A2')
# #sapply (column_order(TF_hm), function(x) 497 %in% x)

# # Distance matrix ####
# d <- as.dist(1 - cor(t(mMat), method='pearson'))

# # Hierarchical clustering ####
# hc <- hclust(d)

# # Dendrogram ####
# pdf (file.path ('Plots',paste0('TF_',metaGroupName,'_no_km_dendrogram.pdf')), width=3, height=3.6)
# plot(hc)
# dev.off()

# ### 

# if (!exists('fSE')) fSE = fetch_mat(archp, 'scATAC_datasets')
# fMat = assays (fSE)[[1]]
# rownames (fMat) = gsub ('rawlins_fetal_lung_','', rownames (fMat))

# if (!exists('nSE')) nSE = fetch_mat(archp, 'scATAC_normal')
# nMat = assays (nSE)[[1]]
# rownames (nMat) = gsub ('Tsankov_lung_','', rownames (nMat))

# if (!exists('nSE2')) nSE2 = fetch_mat(archp, 'scATAC_normal2')
# nMat2 = assays (nSE2)[[1]]
# rownames (nMat2) = gsub ('Tsankov_lung_','', rownames (nMat2))
# rownames(nMat2) = paste0(rownames(nMat2),'2')

# if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name

# # Check TF activity correlation with FRIP
# res = cor (as.matrix(t (mMat)) , archp$FRIP)
# res[order(-res), ,drop=F]

# if (!exists('gsSE')) gsSE = fetch_mat(archp, 'GeneScore')
# gsMat = assays (gsSE)[[1]]
# rownames (gsMat) = rowData (gsSE)$name

# Compare fetal vs normal endothelial cells using deviations ####
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = scale (assays (mSE)[[1]])
rownames (mMat) = rowData (mSE)$name

if (!exists('gsSE')) gsSE = fetch_mat(archp, 'GeneScore')
gsMat = scale (assays (gsSE)[[1]])
rownames (gsMat) = rowData (gsSE)$name

if (!exists('fSE')) fSE = fetch_mat(archp, 'scATAC_datasets')
fMat = scale(assays (fSE)[[1]])
rownames (fMat) = gsub ('rawlins_fetal_lung_','', rownames (fMat))

if (!exists('nSE2')) nSE2 = fetch_mat(archp, 'scATAC_normal2')
nMat2 = scale(assays (nSE2)[[1]])
rownames (nMat2) = gsub ('Tsankov_lung_','', rownames (nMat2))
rownames(nMat2) = paste0(rownames(nMat2),'2')
rownames(nMat2) = gsub ('2$','',rownames(nMat2))


fetal_normal_endo = c('Arterialendo','Latecap','Midcap',#'VascularSMC',
  'Venousendo','Endothelial')
fnmMat = cbind (t(fMat), t(nMat2))
fnmMat = fnmMat[,fetal_normal_endo]
fnmMat = scale(t(fnmMat))
all (rownames(archp@cellColData) == colnames(fnmMat))
mMats = t(mMat[c('MEF2C','ETS1'),])

ha = HeatmapAnnotation (depth = archp$FRIP, which='row')
ha2 = HeatmapAnnotation (depth = archp$Sample, which='row',col=list(depth = palette_sample))
hm = Heatmap (t(fnmMat), left = ha,
  right =ha2, row_split = archp$Sample,
  column_names_gp = gpar(fontsize = 8), col= rev(palette_deviation),
  #column_split = colnames(fnmMat) %in% c('MEF2C','ETS1'),
  row_names_gp = gpar(fontsize = 0))
hm1 = Heatmap (mMats,
  column_names_gp = gpar(fontsize = 8), col= rev(palette_deviation),
  #column_split = colnames(fnmMat) %in% c('MEF2C','ETS1'),
  row_names_gp = gpar(fontsize = 0))

gsMats = gsMat[c('PLVAP','COL4A1','MEF2C','ETS1'),]
hm2 = Heatmap (as.matrix((t(gsMats))), 
  row_split = archp$Sample,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 0),
  col=palette_expression)

pdf (file.path ('Plots','fetal_TFs_heatmap.pdf'), width = 4,height=6)
draw (hm + hm1 + hm2)
dev.off()


hm4 = Heatmap (cor (as.matrix(fnmMat)),
  column_km = 3,
  row_km = 3,
  #row_split = archp$Sample,
  column_names_gp = gpar(fontsize = 0),
  row_names_gp = gpar(fontsize = 0),
  col=palette_deviation_cor_fun)

pdf (file.path ('Plots','fetal_normal_cor_heatmap.pdf'), width = 4,height=4)
hm4
dev.off()
   
fnmMatso = fnmMat[,order(-fnmMat['Endothelial',])]
#fnmMatso = fnmMatso[,!colnames(fnmMatso) %in% 'VascularSMC']

# Also correlate TFs
all (colnames(fnmMat) == colnames(mMat))
metaGroupName = 'celltype'
ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat), group.by = metaGroupName)[[1]]) +1)
min_exp = 0.5
active_TFs = rownames(ps)[rowSums(ps) > min_exp]
#positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0]
metaGroupName = 'celltype'
mMat = mMat[active_TFs,]

cor_tf = as.data.frame (t(cor (t(fnmMat['Endothelial',, drop=F]), t(as.matrix(scale(mMat))))))
head (cor_tf[order(cor_tf$Endothelial),,drop=F], 50)
#cor_tf[grep ('HOX', rownames(cor_tf)),,drop=F]

ha = HeatmapAnnotation (depth = archp$FRIP)
hm5 = Heatmap (fnmMatso,
  cluster_columns=FALSE,
  top = ha,
  cluster_rows=FALSE,
  #row_split = archp$Sample,
  column_names_gp = gpar(fontsize = 0),
  row_names_gp = gpar(fontsize = 12),
  col=palette_deviation_cor_fun,
  border=T)

pdf (file.path ('Plots','fetal_normal_ordered_heatmap.pdf'), width = 4,height=1.3)
hm5
dev.off()



## Check RNA fetal markers ####
logfcThreshold = 0.2
Idents (srt) = 'celltype2'
degClusters = FindAllMarkers (srt, max.cells.per.ident = 1000, min.pct = .1, logfc.threshold = logfcThreshold, verbose = T)
head (degClusters[degClusters$cluster == 'COL4A1',], 10)
top_genes=20
degClusters = degClusters[degClusters$avg_log2FC > 0,]
top_deg = degClusters %>% arrange(p_val_adj) %>% group_by (cluster)  %>% slice(1:top_genes)

if (!exists('gsSE')) gsSE = fetch_mat(archp, 'GeneScore')
gsMat = assays (gsSE)[[1]]
rownames (gsMat) = rowData (gsSE)$name
gsMat = gsMat [top_deg$gene,]
gsMat = as.data.frame (t(gsMat))
all (rownames(gsMat) == rownames(archp@cellColData))
gsMat$metagroup = archp$fetal2
gsMat = aggregate (.~metagroup, gsMat, mean)
rownames (gsMat) = gsMat[,1]
gsMat = gsMat[,-1]

pdf (file.path ('Plots','genescore_fetal_markers_heatmap.pdf'), height=8,width=4)
Heatmap (t(scale(gsMat)), col = palette_expression, column_names_gp = gpar(fontsize = 9), column_names_rot=45,
row_names_gp = gpar(fontsize = 9))
dev.off()

### Use gene scores from RNA markers to drive annotation of endothelial cells ####
archp = addModuleScore (
  ArchRProj = archp,
  useMatrix = 'GeneScoreMatrix',
  name = "fetal",
  features = split(top_deg$gene, top_deg$cluster),
  nBin = 25,
  nBgd = 100,
  seed = 1,
  threads = getArchRThreads(),
  logFile = createLogFile("addModuleScore")
)


# Export BigWig  ####
all (rownames(archp@cellColData) == colnames(fnmMat))
archp$fetal_group = ifelse (archp$fetal.COL4A1 > .5, 'fetal','adult')

metaGroupName = 'fetal_group'
exp_bigwig = T
if (exp_bigwig)
  {
  getGroupBW(
    ArchRProj = archp,
    groupBy = metaGroupName,
    normMethod = "ReadsInTSS",
    tileSize = 100,
    maxCells = 1000,
    ceiling = 4,
    verbose = TRUE,
    threads = getArchRThreads(),
    logFile = createLogFile("getGroupBW")
  )
  }
archp = saveArchRProject (archp)




### Make heatmap ordered by module score of fetal markers ####
sams = c('P1','P10','P11','P14','P4') # Select samples that have at least 100 endothelial cells

if (!exists('fSE')) fSE = fetch_mat(archp, 'scATAC_datasets')
fMat = assays (fSE)[[1]]
rownames (fMat) = gsub ('rawlins_fetal_lung_','', rownames (fMat))

if (!exists('nSE')) nSE = fetch_mat(archp, 'scATAC_normal2')
nMat = assays (nSE)[[1]]
rownames (nMat) = gsub ('Tsankov_lung_','', rownames (nMat))

if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = as.data.frame (t(scale(assays (mSE)[[1]])))
colnames (mMat) = rowData (mSE)$name

metaGroupName = 'celltype2'
ps = log2(as.data.frame (AverageExpression (srt, features = colnames(mMat), group.by = metaGroupName)[[1]]) +1)
min_exp = 0.5
active_TFs = rownames(ps)[rowSums(ps) > min_exp]
#positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0]
mMat = mMat[,active_TFs]

modMat = as.data.frame (t(scale(t(archp@cellColData[, c('fetal.Artery','fetal.COL4A1','fetal.Vein')]))))

fetal_normal_endo = c('Latecap','Midcap','Endothelial')
fnmMat = cbind (scale(t(fMat)), scale(t(nMat)))
#fnmMat = fnmMat[, fetal_normal_endo]
#fnmMat = as.data.frame (t(scale (t(fnmMat))))
fnmMat = fnmMat[,!duplicated(colnames(fnmMat))]

fscore_l = list()
mMats_l = list()
fnmMats_l = list()
top_TFs=20

for (sam in sams)
  {
  modMats = modMat[archp$Sample %in% sam,]
  fetal_order = order(modMats$fetal.COL4A1)
  fscore = modMats$fetal.COL4A1[fetal_order]
  fscore_l[[sam]] = fscore
  
  mMats = mMat[archp$Sample %in% sam,]
  cor_tfs = cor (mMats, modMats$fetal.COL4A1, method='spearman')
  mMats = as.data.frame(scale(mMats[fetal_order,head(order (-cor_tfs[,1]),top_TFs)]))
  mMats_l[[sam]] = mMats
  
  fnmMats = fnmMat[archp$Sample %in% sam,]
  fnmMats = fnmMats[fetal_order,]
  fnmMats_l[[sam]] = fnmMats
  }



# Generate heatmap of one sample ####
sam = 'P1'

modMats = modMat[archp$Sample %in% sam,]
fetal_order = order(modMats$fetal.COL4A1)
fscore = modMats$fetal.COL4A1[fetal_order]

mMats = mMat[archp$Sample %in% sam,]
cor_tfs = cor (mMats, modMats$fetal.COL4A1, method='spearman')
mMats = as.data.frame(scale(mMats[fetal_order,head(order (-cor_tfs[,1]),top_TFs)]))

fnmMats = fnmMat[archp$Sample %in% sam,]
fnmMats = as.data.frame(fnmMats[fetal_order,])

# Apply rolling mean for each column with overlapping bins
library (zoo)
bin_width <- 10   # Number of observations per bin
overlap <- 1

fscore = rollapply(fscore, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
fnmMats <- as.data.frame(lapply(fnmMats, function(x) {
  rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
}))
mMats <- as.data.frame(lapply(mMats, function(x) {
  rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
}))


ha = HeatmapAnnotation (
  fetal = fscore,#,which='row'#,
  col =list(fetal = colorRamp2(c(min(fscore),0,max(fscore)), paletteer::paletteer_d("beyonce::X41",3))))

# Add distribution of fetal and adult endothelial for the other 4 samples    
fnmMats_df = do.call (cbind, lapply (fnmMats_l, function(x) colMeans (tail (x, round(nrow(x)/100*30)))))
ha1 = HeatmapAnnotation (' ' = anno_boxplot(fnmMats_df,axis = FALSE, 
  width = unit(1, "cm"), outline=F, border=F,
    gp = gpar(fill = c(paletteer_d("beyonce::X41",3)[3],
      paletteer_d("beyonce::X41",3)[3],'darkgreen'))), 
  which = 'row')

hm1 = Heatmap (t(fnmMats),
  cluster_columns=FALSE,
  right = ha1,
  top = ha,
  cluster_rows=FALSE,
  #row_split = archp$Sample,
  column_names_gp = gpar(fontsize = 0),
  row_names_gp = gpar(fontsize = 8),
  col=palette_deviation2,
  border=T)



# Add recurrence of TFs in top 20 most correlated across the other 4 samples
mMats_df = table (unlist(lapply (mMats_l, function(x) head(colnames(x),top_TFs))))
mMats_df = setNames(as.vector(mMats_df[colnames(mMats)]), names(mMats_df[colnames(mMats)]))
ha3 = HeatmapAnnotation(' '= anno_barplot(mMats_df, bar_width=1, border=F, 
  gp = gpar(fill = 'white'),
  width = unit(1, "cm")), which='row')

hm2 = Heatmap (t(mMats),
  right = ha3,
  #top = ha2,
  cluster_columns=FALSE,
  cluster_rows=FALSE,
  #row_split = archp$Sample,
  column_names_gp = gpar(fontsize = 0),
  row_names_gp = gpar(fontsize = 8, fontface='italic'),
  col=palette_deviation_cor_fun,
  border=T, width=15)

pdf (file.path ('Plots',paste0('genescore_fetal_markers_heatmap.pdf')), height=4, width=7)
draw (hm1 %v% hm2)
dev.off()


# assign fetal based on scaled fetal score
modMat = archp@cellColData[, c('fetal.Artery','fetal.COL4A1','fetal.Vein')]
modMat = as.data.frame (t(scale(t(modMat))))
fscore = modMat$fetal.COL4A1
fetal_barcodes = unlist(lapply (c('P1','P10','P11','P14'), function(x) 
  {
   tmp = data.frame (fscore = fscore[archp$Sample == x], barcode = rownames(archp@cellColData)[archp$Sample == x])
   head (tmp[order(-tmp$fscore),'barcode'],round(nrow(tmp)/100*30))
  }))
archp$fetal_gss2 = 'adult'
archp$fetal_gss2[match(fetal_barcodes, rownames(archp@cellColData))] = 'fetal'





## Generate scatterplot of fetal and adult TA vs GS ####
if (!exists('fSE')) fSE = fetch_mat(archp, 'scATAC_datasets')
fMat = assays (fSE)[[1]]
rownames (fMat) = gsub ('rawlins_fetal_lung_','', rownames (fMat))

if (!exists('nSE')) nSE = fetch_mat(archp, 'scATAC_normal2')
nMat = assays (nSE)[[1]]
rownames (nMat) = gsub ('Tsankov_lung_','', rownames (nMat))

if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = as.data.frame (t(scale(assays (mSE)[[1]])))
colnames (mMat) = rowData (mSE)$name

if (!exists('gsSE')) gsSE = fetch_mat(archp, 'GeneScore')
gsMat = assays (gsSE)[[1]]
rownames (gsMat) = rowData (gsSE)$name

metaGroupName = 'celltype2'
ps = log2(as.data.frame (AverageExpression (srt, features = colnames(mMat), group.by = metaGroupName)[[1]]) +1)
min_exp = 0.5
active_TFs = rownames(ps)[rowSums(ps) > min_exp]
#positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0]
mMat = mMat[,active_TFs]

# Get endothelial subset module and scale them per cell to reduce seq-depth bias
modMat = as.data.frame (t(scale(t(archp@cellColData[, c('fetal.Artery','fetal.COL4A1','fetal.Vein')]))))

# Combine fetal and normal deviations and scale per cell to reduce seq-depth bias
fetal_normal_endo = c('Latecap','Midcap','Endothelial')
fnmMat = cbind (scale(t(fMat)), scale(t(nMat)))
#fnmMat = fnmMat[, fetal_normal_endo]
#fnmMat = as.data.frame (t(scale (t(fnmMat))))
fnmMat = fnmMat[,!duplicated(colnames(fnmMat))]

# Extract the data frame of interest
sams = c('P1','P10','P11','P14','P4') # Select samples that have at least 100 endothelial cells

ext_dev = 'Endothelial'
p_l=list()
for (sam in sams)
  {

frip = archp$FRIP
all (rownames(archp@cellColData)[archp$Sample == sam] == rownames(as.data.frame(fnmMat[archp$Sample == sam,])))
all (rownames(archp@cellColData)[archp$Sample == sam] == rownames(as.data.frame(modMat[archp$Sample == sam,])))
all (rownames(archp@cellColData)[archp$Sample == sam] == rownames(as.data.frame(mMat[archp$Sample == sam,])))

cells = modMat$fetal.COL4A1[archp$Sample == sam] > -Inf
df <- data.frame (
  fetal_ta = fnmMat[archp$Sample == sam,ext_dev][cells], 
  fetal_gs = modMat$fetal.COL4A1[archp$Sample == sam][cells], 
  MEF2C = mMat$ETS1[archp$Sample == sam][cells],
  FRIP = frip[archp$Sample == sam][cells])


# Create the scatterplot with a polished aesthetic
p_l[[sam]] <- ggplot(df, aes(x = fetal_ta, y = fetal_gs, color= FRIP)) +
  geom_point(
    alpha = 0.7, 
    size = 3
  ) + # Scatterplot points with transparency and size
  geom_smooth(
    method = "lm", 
    color = "dodgerblue", 
    fill = "lightblue", 
    se = TRUE, 
    linetype = "solid", 
    size = 1
  ) + # Regression line with confidence interval
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = " | ")), 
    method = "spearman", 
    #label.x = min(df$fetal) + 0.1 * diff(range(df$fetal)), 
    #label.y = max(df$fetal_gs) - 0.1 * diff(range(df$fetal_gs)),
    color = "darkred",
    size = 5
  ) + # Add correlation coefficient and p-value
  theme_classic() + # Clean theme
  labs(
    title = "fetal gs vs. fetal ta Relationship",
    #subtitle = "Scatterplot with Regression Line and Correlation Coefficient",
    x = "fetal ta",
    y = "fetal gs"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# Print the plot
  }
pdf (file.path ('Plots',paste0('adult_vs_fetal_scatterplot_',ext_dev,'.pdf')), width=12)
wrap_plots(p_l)
dev.off()


# pdf (file.path('Plots','celltype_umap.pdf'),5,5)
#   umap_p3 = plotEmbedding (ArchRProj = archp, 
#     colorBy = "cellColData", name = "Sample",labelMeans =F,
#      embedding = "UMAP_H", pal = palette_sample)
#   # umap_p4 = plotEmbedding (ArchRProj = archp, 
#   #   colorBy = "cellColData", name = "celltype",labelMeans =F,
#   #    embedding = "UMAP_H", pal = palette_stroma)
#   umap_p5 = plotEmbedding (ArchRProj = archp, 
#     colorBy = "cellColData", name = "Clusters_H",
#      embedding = "UMAP_H")
  
#   print (umap_p3)
#   print (umap_p5)
#   dev.off()



# Apply rolling mean for each column with overlapping bins
sams = c('P1','P10','P11','P14','P4') # Select samples that have at least 100 endothelial cells
ext_dev = c('Endothelial')
df = list()

for (sam in sams)
{

library (zoo)
bin_width <- 30   # Number of observations per bin
overlap <- 30
frip = archp$FRIP
fscore_order = modMat$fetal.COL4A1[archp$Sample == sam]
fscore_order = order (fscore_order)

df[[sam]] <- data.frame (
  fetal_ta = rollapply (fnmMat[archp$Sample == sam, ext_dev][fscore_order],width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"), 
  fetal_gs = rollapply (modMat$fetal.COL4A1[archp$Sample == sam][fscore_order], width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"),
  #MEF2C = rollapply (mMat$ELF2[archp$Sample == sam][fscore_order], width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"),
  FRIP = rollapply (frip[archp$Sample == sam][fscore_order], width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"),
  sample = sam)  
}
#df = tail (df,50)
# Create the scatterplot with a polished aesthetic
df = do.call (rbind, df)
p1 <- ggplot(df, aes(x = fetal_ta, y = fetal_gs, color = sample)) +
  geom_point(
    alpha = 0.7, 
    size = 2
  ) + # Scatterplot points with transparency and size
  geom_smooth(
    aes (group = sample),
    method = "lm", 
    color = "dodgerblue", 
    fill = "lightblue", 
    se = FALSE, 
    linetype = "solid", 
    size = 1
  ) + 
  #ylim(c(-1.2,1.2)) + 
  geom_smooth(method=lm, fullrange=FALSE, se=FALSE)+
  # stat_cor(
  # #  aes(label = paste("r = ", ..r.label.., ", p = ", ..p.label.., sep = "")), 
  #   method = "pearson",  # Pearson correlation
  #   label.x = 0, label.y = max(df$fetal_gs),  # Position of label
  #   size = 2, color = "black"
  # ) +
  facet_wrap (~sample, scales = 'free', nrow=length(sams), strip.position = "left") +
  gtheme_no_rot +
  theme(
    strip.placement = "inside",
    panel.spacing = unit(0.5, "lines"), # Adjust spacing between facets
     strip.text = element_blank(), # Customize strip text
    strip.background = element_blank(), # Optionally remove the strip background
    axis.text = element_blank(),       # Remove axis numbers
    axis.ticks = element_blank() 
  ) +
  scale_color_manual(values = palette_sample)   # Regression line with confidence interval
 

# Apply rolling mean for each column with overlapping bins
ext_dev = c('Midcap')
df = list()

for (sam in sams)
{

library (zoo)
bin_width <- 30   # Number of observations per bin
overlap <- 30

fscore_order = modMat$fetal.COL4A1[archp$Sample == sam]
fscore_order = order (fscore_order)

df[[sam]] <- data.frame (
  fetal_ta = rollapply (fnmMat[archp$Sample == sam, ext_dev][fscore_order],width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"), 
  fetal_gs = rollapply (modMat$fetal.COL4A1[archp$Sample == sam][fscore_order], width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"),
  FRIP = rollapply (frip[archp$Sample == sam][fscore_order], width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"),
  sample = sam)
}
#df = tail (df,50)
# Create the scatterplot with a polished aesthetic
df = do.call (rbind, df)
p2 <- ggplot(df, aes(x = fetal_ta, y = fetal_gs, color = sample)) +
  geom_point(
    alpha = 0.7, 
    size = 2
  ) + # Scatterplot points with transparency and size
  geom_smooth(
    aes (group = sample),
    method = "lm", 
    color = "dodgerblue", 
    fill = "lightblue", 
    se = FALSE, 
    linetype = "solid", 
    size = 1
  ) +  
  #ylim(c(-1.2,1.2)) +
  geom_smooth(method=lm, fullrange=FALSE, se=FALSE) +
  # stat_cor(
  #   #aes(label = paste("r = ", ..r.label.., ", p = ", ..p.label.., sep = "")), 
  #   method = "pearson",  # Pearson correlation
  #   size = 2, color = "black"
  # ) +
  facet_wrap (~sample, scales = 'free', nrow=length(sams), strip.position = "left") +
 gtheme_no_rot +
  theme(
    strip.placement = "inside",
    panel.spacing = unit(0.5, "lines"), # Adjust spacing between facets
     strip.text = element_blank(), # Customize strip text
    strip.background = element_blank(), # Optionally remove the strip background
    axis.text = element_blank(),       # Remove axis numbers
    axis.ticks = element_blank() 
  ) +
  scale_color_manual(values = palette_sample)  # Regression line with confidence interval
# Print the plot
pdf (file.path ('Plots',paste0('adult_vs_fetal_scatterplot_smoothed.pdf')), width=5, height=6)
wrap_plots (p1, p2, ncol = 2)
dev.off()



### Find most correlated TF to fetal endothelial score ####

# Apply rolling mean for each column with overlapping bins
sams = c('P1','P10','P11','P14','P4') # Select samples that have at least 100 endothelial cells

cor_l=list()
for (sam in sams)
{
library (zoo)
bin_width <- 30   # Number of observations per bin
overlap <- 30
frip = archp$FRIP
fscore_order = modMat$fetal.COL4A1[archp$Sample == sam]
fscore_order = order (fscore_order)

fetal_gs = rollapply (modMat$fetal.COL4A1[archp$Sample == sam][fscore_order], width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
smoothed_mMat = rollapply (mMat[archp$Sample == sam,][fscore_order,], width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
cor_l[[sam]] = as.data.frame (cor (smoothed_mMat, fetal_gs))
cor_l[[sam]]$sample = sam
cor_l[[sam]]$gene = rownames(cor_l[[sam]])
}
top_genes=30
cor_df = do.call (rbind, cor_l)
TF_order = head(cor_df %>% group_by(gene) %>% summarise(median_value = median(V1)) %>% arrange (-median_value),top_genes)
cor_df = cor_df[cor_df$gene %in% TF_order$gene,]
cor_df$gene = factor (cor_df$gene, levels = rev (unique(TF_order$gene)))
bp = ggplot (cor_df, aes (x = V1, y = gene, alpha=.3), alpha=.5) + 
geom_point (position = position_dodge(.8), alpha= 0.5, color = 'grey22', size=1) +
geom_boxplot (color = 'grey66',fill = 'darkblue',
    linewidth = .1,
    width=0.7,
    outlier.alpha = 0.2,
    outlier.shape = NA,
     size=0.5, alpha=.5
     ) + 
gtheme_no_rot #+
#scale_fill_manual (values = c(activity = 'darkred', genescore = 'navyblue')) + 
#geom_vline (xintercept = 0, color='red',  linetype='dashed')
pdf (paste0 ('Plots/fetal_endothelial_TFs_boxplots.pdf'), width = 5,height=5)
bp
dev.off()




