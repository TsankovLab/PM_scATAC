conda activate meso_scatac

R

set.seed(1234)

####### ANALYSIS of Myeloid compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR'
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

# Set # of threads and genome reference ####
addArchRThreads(threads = 1) 
addArchRGenome("hg38")

# Load ArchR project ####
archp = loadArchRProject (projdir)

# Load RNA ####
srt = readRDS (file.path('..','scrna','srt.rds'))
sample_levels = c('Monocytes','TREM2','SPP1','cDCs','IFN_CXCLs','IM')

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
    groupBy = c('Sample'), force=TRUE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony_project", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony_project",
    name='Clusters_H', res=1,
    force = TRUE)


### Annotate cell types ####
archp$celltype2 = 0
archp$celltype2[archp$Clusters_H %in% c('C6','C7','C9')] = 'Monocytes'
archp$celltype2[archp$Clusters_H %in% c('C8')] = 'cDCs'
archp$celltype2[archp$Clusters_H %in% c('C1','C2','C3','C4','C5','C10','C11')] = 'TAMs'

pdf()
umap_p3 = plotEmbedding (ArchRProj = archp, labelMeans = F,
  colorBy = "cellColData", name = "Sample",
  pal= palette_sample,
   embedding = "UMAP_H")
umap_p4 = plotEmbedding (ArchRProj = archp, labelMeans = F,
  colorBy = "cellColData", name = "celltype2",
   embedding = "UMAP_H",
   pal = palette_myeloid2)
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters_H",
   embedding = "UMAP_H")
dev.off()

pdf (file.path('Plots','celltype_umap_harmony_sample_umap.pdf'),5,5)
print (umap_p3)
print (umap_p4)
print (umap_p5)
dev.off()


# # Check for doublets ####
# meso_markers = c('CD3D','CD3E','EPCAM','KRT19','KRT5','VWF','PECAM1')
# #archp = addImputeWeights (archp)
# pdf()
# p <- plotEmbedding(
#     ArchRProj = archp,
#     colorBy = "GeneScoreMatrix", 
#     name = meso_markers, 
#     size=1,
#     embedding = "UMAP_H",
#     pal = palette_expression,
#     imputeWeights = getImputeWeights(archp)
# )
# dev.off()
# pdf (file.path('Plots','doublets_fplots.pdf'), width = 18, height = 15)
# wrap_plots (p, ncol=3)
# dev.off()


macs_markers=read.csv (file.path ('..','..','..','..','gene_sets','scRNA_immune_markers_humanLUAD_Samarth_Assaf.csv'))
mono_markers = macs_markers[macs_markers$group == 'CD14 mono','gene']
mono_markers = mono_markers[mono_markers != 'CD14 MONO']
mono_markers = c('VCAN','FCN1','CXCL8','CXCL2','IL1B','EREG','TIMP1','THBS1','CCR2','FLT3','FOXM1','CDK1','PCNA','FCGR3A','C1QA','C1QB','CD1C','CD1A','FCER1A')
dc_markers = macs_markers[macs_markers$group %in% c('DC1','DC2','mregDC'),'gene']
ap1_complex = c('JUN','FOSB','FOS','BACH1','SMARCC1','FOSL2','JUND','JDP2','BATF','CEBPB','CEBPA','CEBPZ','FOSL1','NFE2','NFE2L2','NFE2L1')
selected_markers = c('FCN1','EREG','TIMP1', 'CCR7','LAMP3','HLA-DQA1','C1QA','CD68','C5AR2')
pdf()
p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = selected_markers, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path ('Plots','TF_selected_featureplots.pdf'), width = 18,height=14)
wrap_plots (p2)
dev.off()

pdf()
p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = ap1_complex, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path ('Plots','TF_ap1_complex_featureplots.pdf'), width = 28,height=24)
wrap_plots (p2)
dev.off()

# Check sample quality
qc = c('ReadsInTSS','TSSEnrichment','nFrags')
pdf()
p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "cellColData", 
    name = qc, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = NULL
)
dev.off()
pdf (file.path ('Plots','qc_featureplots.pdf'), width = 18,height=14)
wrap_plots (p2)
dev.off()



# # Get markers for gene score ####
# immune_markers = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/scRNA_immune_markers_humanLUAD_Samarth_Assaf.csv')
# immune_markers = immune_markers [immune_markers$group %in% c('Neutrophil','TRMac','IM','DC2','DC1','pDC','mregDC','CD14 mono',
#   'CD16 mono','NK','Mast cell','Mgk','B/Plasma',' T cell','Treg','MoMac'),]
# #immune_markers = immune_markers[immune_markers$group %in% c('CD14 mono','CD16 mono','DC1','DC2','MoMac'),]
# immune_markers = immune_markers$gene
# immune_markers = immune_markers[!immune_markers %in% c('CD14 MONO','IHBA','SEPP1','IL3RA')]
# #archp = addImputeWeights (archp)
# pdf()
# p <- plotEmbedding(
#     ArchRProj = archp,
#     colorBy = "GeneScoreMatrix", 
#     name = immune_markers, 
#     size=1,
#     embedding = "UMAP_H",
#     pal = palette_expression,
#     imputeWeights = getImputeWeights(archp)
# )
# dev.off()
# png (file.path('Plots','myeloid_markers_fplots.png'), width = 18000, height = 15000, res=300)
# wrap_plots (p)
# dev.off()


### Call peaks on celltypes ####
pdf(file.path('Plots','peakcalls.pdf'))
metaGroupName = 'Clusters_H'
force=TRUE
peak_reproducibility=2
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
  


# Add cNMF modules from scRNA-seq ####
shared_cnmf = readRDS (file.path('..','scrna','shared_cnmf_myeloid.rds'))
shared_cnmf = lapply (shared_cnmf, function(x) x[x %in% getFeatures (archp)])
shared_cnmf = lapply (shared_cnmf, function(x) head (x, 50))
remove_modules = c('cDCs2','CC_S','CC_G2M','LILRA') # remove monocyres cDC and CC modules. Consider re-inculding CC 
shared_cnmf = shared_cnmf[!names(shared_cnmf) %in% remove_modules]

force = FALSE
if (!all (names (shared_cnmf) %in% colnames (archp@cellColData)) | force)
  {
  archp@cellColData = archp@cellColData[,!colnames(archp@cellColData) %in% names(shared_cnmf)]
  archp = addModuleScore (
      ArchRProj = archp,
      useMatrix = 'GeneScoreMatrix',
      name = '',
      features = shared_cnmf,
      nBin = 25,
      nBgd = 100,
      seed = 1,
      threads = getArchRThreads(),
      logFile = createLogFile("addModuleScore")
    )
  colnames (archp@cellColData) = gsub ('^\\.','',colnames(archp@cellColData))    
  }

# Assign TAMs to cNMF modules from scRNA-seq ####
archp_MAC = archp[!archp$celltype2 %in% c('Monocytes','cDCs')]
cnmf_scatac = as.data.frame (scale(t(scale(t(archp@cellColData[,names(shared_cnmf)])))))
if (!file.exists ('cnmf_scatac.rds')) saveRDS (cnmf_scatac, 'cnmf_scatac.rds')


cap = 3
cnmf_scatac_cap = cnmf_scatac
cnmf_scatac_cap[cnmf_scatac_cap > cap] = cap
cnmf_scatac_cap[cnmf_scatac_cap < -cap] = -cap
set.seed (123)
km_cnmf = kmeans (t(scale(t(cnmf_scatac_cap))), centers=6) # double scale modules and cluster using k-means
ha = HeatmapAnnotation (sample = archp$Sample, col=list(sample = palette_sample))
hm = Heatmap (t(cnmf_scatac_cap), 
  col = palette_genescore_fun(cnmf_scatac_cap), 
  top_annotation = ha,
  clustering_distance_columns = 'pearson',
  clustering_distance_rows = 'pearson',
  show_column_dend = F,
  column_split = km_cnmf$cluster,
  #column_km=3,
  row_names_gp = gpar (fontsize = 8),
  column_names_gp = gpar (fontsize = 0),
  border=T)

pdf (file.path ('Plots','cnmf_scatac_scaled_only_MAC2_heatmap4.pdf'), height=1.5)
hm
dev.off()
  

# Show TAM modules in UMAPs ####
archp = addImputeWeights (archp)
pdf()
p <- plotEmbedding (
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = names (shared_cnmf), 
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
# p2 <- plotEmbedding (
#     ArchRProj = archp, 
#     colorBy = "cellColData", 
#     name = names (shared_cnmf), 
#     embedding = "UMAP_H",
#     pal = palette_expression,
#     imputeWeights = getImputeWeights(archp)
# )
dev.off()
pdf (file.path ('Plots','shared_cnmf_TAMs_fplots2.pdf'),14,14)
wrap_plots (p, ncol=6)
dev.off()


#### Identify TF regulators correlated to each scRNA cnmf #####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
#mMat = scale(as.matrix(mMat))#[selected_TF,])

# Filter by RNA expression ####
metaGroupName = 'celltype2'
active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)
mMat = t(scale (mMat[active_TFs, ]))
if (!file.exists ('mMat_scaled_active.rds')) saveRDS (mMat, 'mMat_scaled_active.rds')


cnmf_scatac = readRDS ('cnmf_scatac.rds')
mMat = readRDS ('mMat_scaled_active.rds')
tf_cnmf_cor = cor (mMat, cnmf_scatac, method='spearman')
cell_subsets = colnames(tf_cnmf_cor)[c(4,1,3,5,6,2)]
top_5 = unlist(lapply (cell_subsets, function(x) head (rownames(tf_cnmf_cor[order(-tf_cnmf_cor[,x]),]),5)))


DAM_hm = Heatmap (tf_cnmf_cor[top_5,cell_subsets], 
          #row_labels = colnames (mMat_mg),
          #column_title = paste('top',top_genes),
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

pdf (file.path ('Plots','cnmf_clusters_DAM_heatmap4.pdf'), width = 3,height=4)
draw (DAM_hm)
dev.off()

### Try using metacells for each cnmf ####
sams = c('P1','P10','P11','P12','P13','P14','P23','P5') # Select samples that have at least 100 endothelial cells
df = list()

library (zoo)
bin_width <- 30   # Number of observations per bin
overlap <- 30
cnmf_mods = c('Mono','cDCs','TREM2','SPP1','IFN_CXCLs','IM')
#archp$Sample2 = 'sample'
#sams = 'sample'

for (sam in sams)
{
  cnmf_l = list()
  for (cnmf_mod in cnmf_mods)
  {
  metacells_order = cnmf_scatac[,cnmf_mod][archp$Sample == sam]
  metacells_order = order (-metacells_order)
  mMat_sam = mMat[archp$Sample == sam, ]
  mMat_sam = mMat_sam[metacells_order,]
  cnmf_sam = cnmf_scatac[,cnmf_mod][archp$Sample == sam][metacells_order]
  
 cnmf_l[[cnmf_mod]] <- cor (
  rollapply (mMat_sam, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"),
  rollapply (cnmf_sam, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"), 
  method = 'spearman')  
  }
df[[sam]] = as.data.frame (do.call (cbind, cnmf_l))
colnames (df[[sam]]) = cnmf_mods
#df[[sam]]$sample = sam
}

# Convert list to 3D array: rows x cols x n_matrices
array_3d <- array(unlist(df), dim = c(nrow(df[[1]]), ncol(df[[1]]), length(df)))

# Apply median across the third dimension (i.e., across matrices)
median_matrix <- apply(array_3d, c(1, 2), median)
rownames (median_matrix) = rownames (df[[1]])
colnames (median_matrix) = colnames (df[[1]])
top_5 = unlist(lapply (cnmf_mods, function(x) head (rownames(median_matrix[order(-median_matrix[,x]),]),5)))
top_5 = unlist(lapply (cnmf_mods, function(x) head (rownames(df[[1]][order(-df[[1]][,x]),]),5)))

DAM_hm = Heatmap (df[[1]][top_5,cnmf_mods], 
          #row_labels = colnames (mMat_mg),
          #column_title = paste('top',top_genes),
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

pdf (file.path ('Plots','cnmf_clusters_DAM_metacells_heatmap.pdf'), width = 3,height=4)
draw (DAM_hm)
dev.off()


#df = tail (df,50)
# Create the scatterplot with a polished aesthetic
df = do.call (rbind, df)

# Re-Annotate based on cnmf clustering ####
all (names(km_cnmf$cluster) == rownames(archp@cellColData))
archp$cnmf_cluster = paste0('cnmf_cluster_',km_cnmf$cluster)
archp$cnmf_celltypes = archp$cnmf_cluster
archp$cnmf_celltypes[archp$cnmf_cluster == 'cnmf_cluster_1'] = 'IM'
archp$cnmf_celltypes[archp$cnmf_cluster == 'cnmf_cluster_2'] = 'Mono'
archp$cnmf_celltypes[archp$cnmf_cluster == 'cnmf_cluster_3'] = 'SPP1'
archp$cnmf_celltypes[archp$cnmf_cluster == 'cnmf_cluster_4'] = 'IFN_CXCLs'
archp$cnmf_celltypes[archp$cnmf_cluster == 'cnmf_cluster_5'] = 'cDCs'
archp$cnmf_celltypes[archp$cnmf_cluster == 'cnmf_cluster_6'] = 'TREM2'

# Integrate MACs with myeloid annotation
#archp$celltype2[match(rownames(archp_MAC@cellColData), rownames(archp@cellColData))] = archp_MAC$cnmf_celltypes

pdf()
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "cnmf_celltypes",
   embedding = "UMAP_H")
dev.off()

pdf (file.path('Plots','celltype_umap_TAM_annotated_umap2.pdf'),5,5)
print (umap_p5)
dev.off()


write.csv (data.frame (barcode = rownames(archp@cellColData), celltype = archp$celltype2), 'barcode_annotation.csv')
write.csv (data.frame (barcode = rownames(archp@cellColData), celltype = archp$cnmf_celltypes), 'barcode_annotation_cnmf_celltypes.csv')


# Differential Accessed motifs ####
metaGroupName = "cnmf_celltypes"
force=F
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

### Run TF correlation to identify TF modules across TNK cells #### 
mMat = readRDS ('mMat_scaled_active.rds')
mMat_cor = cor (as.matrix(mMat), method = 'spearman')

set.seed(1234)
centers=2
km = kmeans (mMat_cor, centers=centers)
if (!file.exists ('TF_activity_modules.rds')) saveRDS (km, 'TF_activity_modules.rds')


AP1 = c('JUNB',
'FOSL2',
'JUN',
'SMARCC1',
'FOSL1',
'JUND',
'FOS',
'JDP2',
'BACH1',
'FOSB',
'BATF',
'NFE2',
'NFE2L2')
ha2 = rowAnnotation (foo = anno_mark(at = match(AP1,colnames(mMat_cor)), 
    labels = AP1, labels_gp = gpar(fontsize = 7)))

pdf (file.path ('Plots','TF_modules_heatmap3.pdf'), width = 4,height=3)
cor_mMat_hm = draw (Heatmap (mMat_cor,# row_km=15,
  right_annotation = ha2,
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

pdf (file.path ('Plots','TF_modules_heatmap3.pdf'), width = 4.6, height=3)
cor_mMat_hm
dev.off()

# Add metacolumns of average TF modules activity ####
tf_modules = lapply (unique(km$cluster), function(x) colMeans (t(mMat)[rownames(t(mMat)) %in% names(km$cluster[km$cluster == x]),]))
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
    imputeWeights=NULL,
    #useSeqnames='z',
    embedding = "UMAP_H")
dev.off()

pdf (file.path ('Plots','TF_modules_umap2.pdf'), width = 20,height=16)
wrap_plots (TF_p, ncol=5)
dev.off()

# Generate same heatmap but using scrna ####
metacells = readRDS (file.path ('..','scrna','metacells.rds'))
metacells_mat = metacells@assays$RNA$data[rownames(mMat_cor),]
metacells_mat = cor (t(metacells_mat), method = 'spearman')
all (rownames(metacells_mat) == rownames(mMat_cor))
tf_order = unname(unlist(row_order(cor_mMat_hm)))

pdf (file.path ('Plots','TF_modules_RNA_heatmap2.pdf'), width = 4,height=3)
cor_mMat_hm2 = draw (Heatmap (metacells_mat[tf_order,tf_order],# row_km=15,
  row_split = km$cluster,
  column_split = km$cluster,
  cluster_rows=F,
  cluster_columns = F,
  col=palette_expression_cor_fun, 
  border=T,
  row_names_gp = gpar(fontsize = 0), 
  column_names_gp = gpar(fontsize = 0)))
dev.off()

pdf (file.path ('Plots','TF_modules_ATAC_RNA_heatmap2.pdf'), width = 4,height=3)
cor_mMat_hm2
dev.off()


pdf (file.path ('Plots','TF_modules_RNA_heatmap3.pdf'), width = 4,height=3)
cor_mMat_hm2 = draw (Heatmap (metacells_mat[tf_order,tf_order],# row_km=15,
  col=palette_expression_cor_fun, 
  border=T,
  row_names_gp = gpar(fontsize = 0), 
  column_names_gp = gpar(fontsize = 0)))
dev.off()

pdf (file.path ('Plots','TF_modules_ATAC_RNA_heatmap3.pdf'), width = 4,height=3)
cor_mMat_hm2
dev.off()



### Plot UMAP using chromvar TF ####
mMat = readRDS ('mMat_scaled_active.rds')
library (uwot)
set.seed(42)  # for reproducibility
umap_result <- umap(mMat, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
library(ggplot2)

umap_df = as.data.frame(umap_result)
umap_df$celltypes = archp$cnmf_celltypes
umap_df$FRIP = log2 (archp$nFrags+1)
umap_df$Sample = archp$Sample
umap_df$mod_2 = archp$mod_2
umap_df = cbind (umap_df, t(scale(t(as.data.frame (archp@cellColData[,names(shared_cnmf)])))))
sp = lapply (unique(archp$cnmf_celltypes), function(x) 
  ggplot (umap_df[umap_df$celltypes == x,], aes(V1, V2, color = celltypes)) +
  geom_point(size = 2) +
  labs(title = "UMAP of Iris Data", x = "UMAP1", y = "UMAP2") +
  theme_minimal())
sp2 = lapply (names (shared_cnmf), function(x) 
  ggplot (umap_df, aes_string('V1', 'V2', color = x)) +
  geom_point(size = .2) +
  scale_color_viridis_b() +
  labs(title = "chromvar UMAP cnmf scores", x = "UMAP1", y = "UMAP2") +
  theme_minimal())
spq = ggplot (umap_df, aes(V1, V2, color = Sample)) +
  geom_point(size = 2) +
  labs(title = "UMAP of Iris Data", x = "UMAP1", y = "UMAP2") +
  theme_minimal()
spq2 = ggplot (umap_df, aes(V1, V2, color = mod_2)) +
  geom_point(size = 2) +
  scale_color_viridis_b() +
  labs(title = "UMAP of Iris Data", x = "UMAP1", y = "UMAP2") +
  theme_minimal()  
png (file.path ('Plots','chromvar_cnmf_scores_umap.png'), 3000,3000,res=300)
wrap_plots (sp2)
dev.off()

png (file.path ('Plots','chromvar_infl_scores_umap.png'), 3000,3000,res=300)
spq2
dev.off()

pdf (file.path ('Plots','chromvar_umap.pdf'))
sp
spq
dev.off()

sp1 = lapply (unique(umap_df$Sample), function(y) 
  {
  lapply (unique(archp$cnmf_celltypes), function(x) 
  ggplot (umap_df[umap_df$celltypes == x & umap_df$Sample == y,], aes(V1, V2, color = celltypes)) +
  geom_point(size = 2) +
  labs(title = "UMAP of Iris Data", x = "UMAP1", y = "UMAP2") +
  theme_minimal())
  })

pdf (file.path ('Plots','chromvar_per_sample_umap.pdf'), width=10)
lapply (sp1, function(x) wrap_plots(x))
dev.off()



# Make scatterplots of inflammation vs cnmfs ####
cnmf_mods = c('Mono','SPP1','TREM2','IFN_CXCLs','IM')
cnmf_scatac = readRDS ('cnmf_scatac.rds')


### Try using metacells for each cnmf ####
sams = c('P1','P10','P11','P12','P13','P14','P23','P5') # Select samples that have at least 100 endothelial cells
df = list()
archp$Sample2 = 'sample'
sams = 'sample'
library (zoo)
bin_width <- 30   # Number of observations per bin
overlap <- 30
for (sam in sams)
{
  cnmf_l = list()
  for (cnmf_mod in cnmf_mods)
  {
  metacells_order = cnmf_scatac[,cnmf_mod][archp$Sample2 == sam]
  metacells_order = order (-metacells_order)
  mMat_sam = archp$mod_2[archp$Sample2 == sam]
  mMat_sam = mMat_sam[metacells_order]
  cnmf_sam = cnmf_scatac[,cnmf_mod][archp$Sample2 == sam]
  cnmf_sam = cnmf_sam[metacells_order]
  
 cnmf_l[[cnmf_mod]] <- cor (
  rollapply (mMat_sam, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"),
  rollapply (cnmf_sam, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"), 
  method = 'spearman')  
  }
df[[sam]] = as.data.frame (do.call (cbind, cnmf_l))
colnames (df[[sam]]) = cnmf_mods
#df[[sam]]$sample = sam
}

# Convert list to 3D array: rows x cols x n_matrices
array_3d <- array(unlist(df), dim = c(nrow(df[[1]]), ncol(df[[1]]), length(df)))

# Apply median across the third dimension (i.e., across matrices)
median_matrix <- apply(array_3d, c(1, 2), median)
rownames (median_matrix) = rownames (df[[1]])
colnames (median_matrix) = colnames (df[[1]])
top_5 = unlist(lapply (cnmf_mods, function(x) head (rownames(median_matrix[order(-median_matrix[,x]),]),5)))




atac_mat = cbind (cnmf_scatac, mod_2 = archp$mod_2, sample = archp$Sample)
atac_mat_long = gather (as.data.frame(atac_mat), cnmf, score, 1:(ncol(atac_mat)-2))
atac_mat_long$cnmf = factor (atac_mat_long$cnmf, levels = cnmf_mods)

library (ggpointdensity)
# remove outliers
#atac_mat_longL = split (atac_mat_long, atac_mat_long$cnmf)
#atac_mat_longL = lapply (atac_mat_longL, function(x) x[x$score > quantile(x$score,.005) & x$score < quantile(x$score,.995),])
#atac_mat_long = do.call (rbind,atac_mat_longL)
sp = ggplot (atac_mat_long, aes(x = score, y = mod_2, fill = sample)) +
  geom_point(shape = 21,
    alpha = 0.6, 
    size = 2
  ) + 
  facet_wrap (~ interaction(cnmf, sample, sep = " | "), scales = 'free', ncol= length(cnmf_mods)) #+
   #geom_pointdensity (alpha=1, size=.1)# +
#   scale_color_viridis (option='F') +
# # Scatterplot points with transparency and size
#   geom_smooth(
#     method = "lm", 
#     color = "white", 
#     fill = "white", 
#     se = FALSE, 
#     linetype = "dashed", 
#     size = .4
#   ) + # Regression line with confidence interval
#   theme_void() + # Clean theme
#   labs(
#     title = "mod_2 vs cnmfs",
#     #subtitle = "Scatterplot with Regression Line and Correlation Coefficient",
#     x = "cnmf",
#     y = "mod_2")
#   # ) +
  # theme(
  #   plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
  #   plot.subtitle = element_text(size = 14, hjust = 0.5),
  #   axis.title = element_text(size = 14),
  #   axis.text = element_text(size = 12),
  #   panel.grid.major = element_line(color = "grey90"),
  #   panel.grid.minor = element_blank()
  # ) 

png (file.path ('Plots','cnmf_km_sample_scatterplots.png'),width=3000, height=5500, res=300)
sp
dev.off()

sp = sp + stat_cor (
    aes(label = paste(..rr.label.., ..p.label.., sep = " | ")), 
    method = "spearman", 
    #label.x = min(df$fetal) + 0.1 * diff(range(df$fetal)), 
    #label.y = max(df$fetal_gs) - 0.1 * diff(range(df$fetal_gs)),
    color = "grey11",
    size = 1
  )   
pdf (file.path ('Plots','cnmf_AP1_scatterplots3.pdf'),width=7, height=1.4)
sp
dev.off()

# # Try with ridge plots ####
library (ggridges)
library (ggplot2)
library (viridis)
library (tidyr)
#library(hrbrthemes)

# Plot
ccomp = as.data.frame (archp@cellColData)
#ccomp = ccomp[ccomp$cnmf_celltypes %in% c('cDCs'),]
ccomp$cnmf_celltypes = factor (ccomp$cnmf_celltypes, levels = rev(c('Mono','TREM2','SPP1','IFN_CXCLs','cDCs','IM')))
ccomp$module = archp$mod_2
rp <- ggplot(ccomp, aes(x = module, y = cnmf_celltypes, fill = ..x..)) +
  geom_density_ridges_gradient(
  scale = 3,
  rel_min_height = 0.01,
  linewidth = 0.4,
  color='white',
  alpha = 0.3
) +

  scale_fill_viridis_c(option = "C") +  # Optional: nice color gradient
  theme_ridges() +                      # Optional: clean ridge plot theme
  theme(legend.position = "right")     # Adjust legend position
#   theme_classic() + facet_wrap (~sample, ncol=5)
pdf (file.path ('Plots','cnmf_inflammation_module_ridge_plots.pdf'), width = 5,height=3)
rp
dev.off()


# Make heatmap of correlation of inflammation module vs cnmf modules ####
# TF_modules = split(names(km$cluster), km$cluster)
# srt = ModScoreCor (
#     seurat_obj = srt, 
#     geneset_list = TF_modules, 
#     cor_threshold = NULL, 
#     pos_threshold = NULL, # threshold for fetal_pval2
#     listName = 'TF_module', outdir = NULL)



if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = scale(as.matrix (mMat[active_TFs,]))#[selected_TF,])
mMat = mMat[names(km$cluster[km$cluster == 2]),]

mMat = colMeans (mMat)
all (colnames (mMat) == rownames(archp@cellColData))
archp$AP1 = mMat
rna_mod_cor = cor (srt@meta.data[,cnmf_mods],srt@meta.data[,'2'])
atac_mod_cor = cor (scale(archp@cellColData[,cnmf_mods]),mMat, method='spearman')
#cor.test (archp@cellColData[,names(shared_cnmf)],archp@cellColData[,'mod_2'])

hm2 = Heatmap (atac_mod_cor, border=T, col = palette_deviation_fun(atac_mod_cor))
hm1 = Heatmap (rna_mod_cor, cluster_rows=F, border=T, col = rev (palette_deviation_correlation))
pdf (file.path ('Plots','atac_module_cor3.pdf'), width=3.5, height=3)
hm2 + hm1
dev.off()





# # Show all TFs included in inflammation module ####
# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = scale(assays (mSE)[[1]])
# rownames (mMat) = rowData (mSE)$name
# mMat = as.matrix(mMat)#[selected_TF,])

# metaGroupName = 'celltype2'
# mMat_mg = mMat[names (km$cluster)[km$cluster==2], ]
# mMat_mg = as.data.frame (t(mMat_mg))
# mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
# mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
# rownames (mMat_mg) = mMat_mg[,1]
# mMat_mg = mMat_mg[,-1]

# hm = Heatmap (
#     t(mMat_mg[,AP1]),
# #    right_annotation = ha2,
#     column_names_rot =45, 
#     row_names_gp = gpar(fontsize = 5),
#     column_names_gp = gpar(fontsize = 6),
#     col = rev(as.character(palette_deviation)), 
#     cluster_rows=T,
#     cluster_columns = T,
#     border=T
# #rect_gp = gpar (col = "white", lwd = 1)
# )


# ### Check AP1 on the RNA side ####
# metaGroupName = 'shared_cnmf2_r_max'
# ps = log2(as.data.frame (AverageExpression (srt, features = AP1, group.by = metaGroupName)[[1]]) +1)

# hm2 = Heatmap (
#     t(scale(t(ps))),
#     #right_annotation = ha2,
#     column_names_rot =45, 
#     row_names_gp = gpar(fontsize = 5),
#     column_names_gp = gpar(fontsize = 6),
#     col = rev(as.character(palette_expression_correlation)), 
#     cluster_rows=T,
#     cluster_columns = T,
#     border=T
# #rect_gp = gpar (col = "white", lwd = 1)
# )


# pdf (file.path ('Plots','inflammation_module_atac_rna_TFs_heatmap2.pdf'), width = 4.2,height=2)
# hm + hm2
# dev.off()




# ### Plot regulon score of TFs found in km2 along with km2 average score from atac #####
# auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)

# regulon_TFs_in_modules = list(
#   km1 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 1])],
#   km2 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 2])]
#   )
# rownames (auc_mtx) = auc_mtx[,1]
# auc_mtx = auc_mtx[,-1]
# colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))
# auc_mtx = auc_mtx[, colnames(auc_mtx) %in% regulon_TFs_in_modules$km2]
# auc_mtx_avg = aggregate (auc_mtx, by=as.list(srt@meta.data[,'shared_cnmf2_r_max',drop=F]), mean)
# rownames (auc_mtx_avg) = auc_mtx_avg[,1]
# auc_mtx_avg = auc_mtx_avg[,-1]
# auc_mtx_avg = auc_mtx_avg[!rownames(auc_mtx_avg) %in% c('IL1B','cDCs'),]
# auc_mtx_avg_scaled = as.data.frame (scale (auc_mtx_avg))
# auc_mtx_avg_scaled$celltype = rownames(auc_mtx_avg_scaled)
# auc_mtx_avg_scaled_l = gather (auc_mtx_avg_scaled, TF, score,1:(ncol(auc_mtx_avg_scaled)-1))
# auc_mtx_avg_scaled_l$celltype = factor (auc_mtx_avg_scaled_l$celltype, levels = c('Monocytes','SPP1','TREM2','IFN','C1Q','IM'))


# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat = as.matrix(mMat)#[selected_TF,])
# atac_mod_2 = as.data.frame(t(mMat[regulon_TFs_in_modules$km2,]))
# atac_mod_2 = aggregate (atac_mod_2, by=list(celltype = archp$celltype2), mean)
# atac_mod_2 = atac_mod_2[atac_mod_2$celltype != 'cDCs',]
# atac_mod_2 = atac_mod_2 %>%
#     mutate_if(is.numeric, scale)
# atac_mod_2 = gather (atac_mod_2, TF, score,2:(ncol(atac_mod_2)))

# ### Plot  TF activity
# atac_mod_2_summary <- atac_mod_2 %>%
#   group_by(celltype) %>%
#   dplyr::  summarize(
#     mean_score = mean(score, na.rm = TRUE),
#     sd_score = sd(score, na.rm = TRUE),
#     n = n(),
#     se = sd_score / sqrt(n)  # standard error
#   )
# atac_mod_2_summary$celltype = factor (atac_mod_2_summary$celltype, levels = c('Monocytes','SPP1','TREM2','C1Q','IFN','IM'))

# gp = ggplot(atac_mod_2_summary, aes(x = celltype, y = mean_score, group = 1)) +
#   geom_line(color = "darkred", size = 1.5) +
#   geom_ribbon(aes(ymin = mean_score - se, ymax = mean_score + se),
#               fill = "darkred", alpha = 0.2) +
#   theme_minimal() +
#   labs(
#     x = "Celltype",
#     y = "Mean Score",
#     title = "Mean Score per Celltype with Standard Error Shading"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ### Plot SCENIC regulon scores

# rna_mod_2_summary <- auc_mtx_avg_scaled_l %>%
#   group_by(celltype) %>%
#   dplyr::  summarize(
#     mean_score = mean(score, na.rm = TRUE),
#     sd_score = sd(score, na.rm = TRUE),
#     n = n(),
#     se = sd_score / sqrt(n)  # standard error
#   )
# rna_mod_2_summary$celltype = factor (rna_mod_2_summary$celltype, levels = c('Monocytes','SPP1','TREM2','IFN','C1Q','IM'))


# gp = ggplot() +
#   # First dataset (normal scale, left axis)
#   geom_line(data = atac_mod_2_summary, aes(x = celltype, y = mean_score, group=1), color = "darkred", size = .5) +
#   geom_ribbon(data = atac_mod_2_summary, aes(x = celltype, ymin = mean_score - se, ymax = mean_score + se, group=1), 
#               fill = "darkred", alpha = 0.2) +

#   # Second dataset (scaled, right axis)
#   geom_line(data = rna_mod_2_summary, aes(x = celltype, y = mean_score, group=1), color = "navyblue", size = .5) +
#   geom_ribbon(data = rna_mod_2_summary, aes(x = celltype, ymin = mean_score - se, ymax = mean_score + se, group=1),
#               fill = "navyblue", alpha = 0.2) +

#   theme_minimal() +
#   labs(
#     x = "Celltype",
#     title = "Overlayed Mean Score Lines with Separate Y Axes"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


# pdf (file.path ('Plots','inflamed_module_atac_rna_lineplot2.pdf'),width=5,height=5)
# gp
# dev.off()






### Plot correlation of regulon score of TFs found in km2 along with correlation of km2 average score from atac #####
#auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))
regulon_TFs_in_modules = list(
  km1 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 1])],
  km2 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 2])]
  )
auc_mtx = auc_mtx[, colnames(auc_mtx) %in% regulon_TFs_in_modules$km2]

#colnames(srt@meta.data)[colnames(srt@meta.data) == 'Mono'] = 'Monocytes'
cnmf_mods = c('Mono','SPP1','TREM2','IFN_CXCLs','IM')
rownames(auc_mtx) = gsub ('\\.','-',rownames(auc_mtx))
all (rownames(auc_mtx) == colnames(srt))

auc_mtx_cor = as.data.frame (cor (auc_mtx, srt@meta.data[,cnmf_mods]))
auc_mtx_cor$TF = rownames(auc_mtx_cor)

auc_mtx_avg_scaled_l = gather (auc_mtx_cor, celltype, score,1:(ncol(auc_mtx_cor)-1))
auc_mtx_avg_scaled_l$celltype = factor (auc_mtx_avg_scaled_l$celltype, levels = c('Mono','SPP1','TREM2','IFN_CXCLs','IM'))


if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat)#[selected_TF,])
atac_mod_2 = as.data.frame(t(mMat[regulon_TFs_in_modules$km2,]))

atac_mod_2 = as.data.frame (cor (atac_mod_2, t(scale(t(as.data.frame(archp@cellColData[,cnmf_mods]))))))
atac_mod_2$TF = rownames(atac_mod_2)
atac_mod_2 = gather (atac_mod_2, celltype, score,1:(ncol(atac_mod_2)-1))
atac_mod_2$celltype = factor (atac_mod_2$celltype, levels = cnmf_mods)

### Plot  TF activity
atac_mod_2_summary <- atac_mod_2 %>%
  group_by(celltype) %>%
  dplyr::  summarize(
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    n = n(),
    se = sd_score / sqrt(n)  # standard error
  )
atac_mod_2_summary$celltype = factor (atac_mod_2_summary$celltype, levels = cnmf_mods)

gp = ggplot(atac_mod_2_summary, aes(x = celltype, y = mean_score, group = 1)) +
  geom_line(color = "darkred", size = 1.5) +
  geom_ribbon(aes(ymin = mean_score - se, ymax = mean_score + se),
              fill = "darkred", alpha = 0.2) +
  theme_minimal() +
  labs(
    x = "Celltype",
    y = "Mean Score",
    title = "Mean Score per Celltype with Standard Error Shading"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Plot SCENIC regulon scores
rna_mod_2_summary <- auc_mtx_avg_scaled_l %>%
  group_by(celltype) %>%
  dplyr::  summarize(
    mean_score = mean(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    n = n(),
    se = sd_score / sqrt(n)  # standard error
  )
rna_mod_2_summary$celltype = factor (rna_mod_2_summary$celltype, levels = cnmf_mods)


gp = ggplot() +
  # First dataset (normal scale, left axis)
  geom_line(data = atac_mod_2_summary, aes(x = celltype, y = mean_score, group=1), color = "darkred", size = .5) +
  geom_ribbon(data = atac_mod_2_summary, aes(x = celltype, ymin = mean_score - se, ymax = mean_score + se, group=1), 
              fill = "darkred", alpha = 0.2) +

  # Second dataset (scaled, right axis)
  geom_line(data = rna_mod_2_summary, aes(x = celltype, y = mean_score, group=1), color = "navyblue", size = .5) +
  geom_ribbon(data = rna_mod_2_summary, aes(x = celltype, ymin = mean_score - se, ymax = mean_score + se, group=1),
              fill = "navyblue", alpha = 0.2) +

  theme_minimal() +
  labs(
    x = "Celltype",
    title = "Overlayed Mean Score Lines with Separate Y Axes"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey44", size = .5) 


pdf (file.path ('Plots','inflamed_module_atac_rna_lineplot3.pdf'),width=5,height=5)
gp
dev.off()








# ### Plot correlation of regulon score of TFs found in km1 along with correlation of km1 average score from atac #####
# auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
# rownames (auc_mtx) = auc_mtx[,1]
# auc_mtx = auc_mtx[,-1]
# colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))
# auc_mtx = auc_mtx[, colnames(auc_mtx) %in% regulon_TFs_in_modules$km1]

# colnames(srt@meta.data)[colnames(srt@meta.data) == 'Mono'] = 'Monocytes'
# cnmf_mods = c('IL1B','Monocytes','IFN','C1Q','IM','SPP1','TREM2')
# rownames(auc_mtx) = gsub ('\\.','-',rownames(auc_mtx))
# all (rownames(auc_mtx) == colnames(srt))

# auc_mtx_cor = as.data.frame (cor (auc_mtx, srt@meta.data[,cnmf_mods]))
# auc_mtx_cor$TF = rownames(auc_mtx_cor)

# auc_mtx_avg_scaled_l = gather (auc_mtx_cor, celltype, score,1:(ncol(auc_mtx_cor)-1))
# auc_mtx_avg_scaled_l$celltype = factor (auc_mtx_avg_scaled_l$celltype, levels = c('IL1B','Monocytes','SPP1','TREM2','IFN','C1Q','IM'))


# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat = as.matrix(mMat)#[selected_TF,])
# atac_mod_2 = as.data.frame(t(mMat[regulon_TFs_in_modules$km1,]))

# new_cnmf_names = c(
#   cnmf.1 = 'TREM2', 
#   cnmf.3 = 'Monocytes',
#   cnmf.4 ='SPP1',
#   cnmf.5 = 'IL1B',
#   cnmf.6 = 'cDCs',
#   cnmf.8 = 'IFN',
#   cnmf.9 = 'IM',
#   cnmf.10 = 'C1Q'
#   )

# if (all (names(new_cnmf_names) %in% colnames(archp@cellColData))) colnames(archp@cellColData)[match(names(new_cnmf_names), colnames(archp@cellColData))] = new_cnmf_names

# atac_mod_2 = as.data.frame (cor (atac_mod_2, t(scale(t(as.data.frame(archp@cellColData[,cnmf_mods]))))))
# atac_mod_2$TF = rownames(atac_mod_2)
# atac_mod_2 = gather (atac_mod_2, celltype, score,1:(ncol(atac_mod_2)-1))
# atac_mod_2$celltype = factor (atac_mod_2$celltype, levels = c('IL1B','Monocytes','SPP1','TREM2','IFN','C1Q','IM'))

# ### Plot  TF activity
# atac_mod_2_summary <- atac_mod_2 %>%
#   group_by(celltype) %>%
#   dplyr::  summarize(
#     mean_score = mean(score, na.rm = TRUE),
#     sd_score = sd(score, na.rm = TRUE),
#     n = n(),
#     se = sd_score / sqrt(n)  # standard error
#   )
# atac_mod_2_summary$celltype = factor (atac_mod_2_summary$celltype, levels = c('IL1B','Monocytes','SPP1','TREM2','C1Q','IFN','IM'))

# gp = ggplot(atac_mod_2_summary, aes(x = celltype, y = mean_score, group = 1)) +
#   geom_line(color = "darkred", size = 1.5) +
#   geom_ribbon(aes(ymin = mean_score - se, ymax = mean_score + se),
#               fill = "darkred", alpha = 0.2) +
#   theme_minimal() +
#   labs(
#     x = "Celltype",
#     y = "Mean Score",
#     title = "Mean Score per Celltype with Standard Error Shading"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ### Plot SCENIC regulon scores
# rna_mod_2_summary <- auc_mtx_avg_scaled_l %>%
#   group_by(celltype) %>%
#   dplyr::  summarize(
#     mean_score = mean(score, na.rm = TRUE),
#     sd_score = sd(score, na.rm = TRUE),
#     n = n(),
#     se = sd_score / sqrt(n)  # standard error
#   )
# rna_mod_2_summary$celltype = factor (rna_mod_2_summary$celltype, levels = c('IL1B','Monocytes','SPP1','TREM2','IFN','C1Q','IM'))


# gp = ggplot() +
#   # First dataset (normal scale, left axis)
#   geom_line(data = atac_mod_2_summary, aes(x = celltype, y = mean_score, group=1), color = "darkred", size = .5) +
#   geom_ribbon(data = atac_mod_2_summary, aes(x = celltype, ymin = mean_score - se, ymax = mean_score + se, group=1), 
#               fill = "darkred", alpha = 0.2) +

#   # Second dataset (scaled, right axis)
#   geom_line(data = rna_mod_2_summary, aes(x = celltype, y = mean_score, group=1), color = "navyblue", size = .5) +
#   geom_ribbon(data = rna_mod_2_summary, aes(x = celltype, ymin = mean_score - se, ymax = mean_score + se, group=1),
#               fill = "navyblue", alpha = 0.2) +

#   theme_minimal() +
#   labs(
#     x = "Celltype",
#     title = "Overlayed Mean Score Lines with Separate Y Axes"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey44", size = .5) 


# pdf (file.path ('Plots','km1_module_atac_rna_lineplot2.pdf'),width=5,height=5)
# gp
# dev.off()




# Show all TFs included in inflammation module ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat)#[selected_TF,])
atac_mod_2 = as.data.frame(t(mMat[regulon_TFs_in_modules$km2,]))

atac_mod_2 = as.data.frame (cor (atac_mod_2, t(scale(t(as.data.frame(archp@cellColData[,cnmf_mods]))))))
#atac_mod_2$TF = rownames(atac_mod_2)

hm = Heatmap (
    atac_mod_2[,cnmf_mods],
#    right_annotation = ha2,
    column_names_rot =45, 
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 6),
    col = rev(as.character(palette_deviation[-1])), 
    cluster_rows=T,
    cluster_columns = F,
    border=T
#rect_gp = gpar (col = "white", lwd = 1)
)


### Check AP1 on the RNA side ####
auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))
auc_mtx = auc_mtx[, colnames(auc_mtx) %in% regulon_TFs_in_modules$km2]

rownames(auc_mtx) = gsub ('\\.','-',rownames(auc_mtx))
all (rownames(auc_mtx) == colnames(srt))

auc_mtx_cor = as.data.frame (cor (auc_mtx, srt@meta.data[,cnmf_mods]))
#auc_mtx_cor$TF = rownames(auc_mtx_cor)

hm2 = Heatmap (
    auc_mtx_cor[,cnmf_mods],
    #right_annotation = ha2,
    column_names_rot =45, 
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 6),
    col = rev(as.character(palette_expression_correlation[-1])), 
    cluster_rows=T,
    cluster_columns = F,
    border=T
#rect_gp = gpar (col = "white", lwd = 1)
)


pdf (file.path ('Plots','inflammation_module_atac_rna_TFs_cor_heatmap2.pdf'), width = 3.6,height=3)
hm 
dev.off()


















### Check footprint of AP1-complex and NFKB1 ####
metaGroupName='inflamed'
archp <- addGroupCoverages (ArchRProj = archp, groupBy = metaGroupName)
motifPositions <- getPositions (archp)

motifs <- c('NFKB1','JUNB','FOS','JUND')

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
plotName = "Footprints-Subtract-Bias_inflamed_",
addDOC = FALSE, height=4, width=2,
#pal = palette_tnk_cells,
smoothWindow = 25)



### Check inflammation score across TAMs ####
mod_df = data.frame (
  #celltype = archp$cnmf_celltypes,
  #sample = archp$Sample,
  #Infl_module = archp$mod_2,
  Infl_module = archp$mod_2)
mod_df = aggregate(archp$mod_2, by=list(
  celltype = archp$cnmf_celltypes,
  sample = archp$Sample), FUN=mean)
head (mod_df)
df_order = mod_df %>% 
group_by (celltype) %>% 
summarize (avg_module = median(x)) %>% 
arrange(avg_module)
mod_df$celltype = factor (mod_df$celltype, levels = rev(df_order$celltype))

bp = ggplot (mod_df, aes (x = celltype, y = x, fill=celltype)) +
vlp + 
bxpv + 
scale_fill_manual (values = palette_myeloid) +
#geom_point (position='identity', alpha=.3, color="grey44", size=1) +
gtheme

pdf (file.path ('Plots','celltype_infl_module_sample_boxplots.pdf'),2.3,width=4)
bp
dev.off()


# # Compute co-occurrence of TFs ####
km = readRDS ('TF_activity_modules.rds')
infl_TF = names(km$cluster[km$cluster == 2])
metaGroupName = 'inflamed'
celltypes = unique (archp@cellColData[,metaGroupName])
motifMat = getPositions (archp)
matches = getMatches (archp)

# Find DAP ####
#force = FALSE
force=F
if (!file.exists (paste0('DAP_',metaGroupName,'.rds')) | force)
  {
  DAP_list = getMarkerFeatures (
    ArchRProj = archp, 
    testMethod = "wilcoxon",
          useGroups = celltypes[1],
          bgdGroups = celltypes[2],
    k=100,
    binarize = FALSE,
    useMatrix = "PeakMatrix",
    groupBy = metaGroupName
  #  useSeqnames="z"
  )
  saveRDS (DAP_list, paste0('DAP_',metaGroupName,'.rds'))
  } else {
  DAP_list = readRDS (paste0('DAP_',metaGroupName,'.rds'))
  }
DAP_res = do.call (cbind, (assays(DAP_list)))
colnames (DAP_res) = names(assays(DAP_list))
DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames(DAP_res) = as.character(DAP_res_regions)
DAP_res = DAP_res[!is.na(DAP_res$FDR),]
DAP_res = DAP_res[DAP_res$FDR < 0.05,]
DAP_res_l = list(up = GRanges(rownames(DAP_res)[DAP_res$Log2FC > 0]),
                down = GRanges(rownames(DAP_res)[DAP_res$Log2FC < 0]))
sapply (DAP_res_l, length)
names (DAP_res_l) = celltypes

cooc_l = list()
ov_size_l=list()
for (celltype in celltypes)
  {
  matches_ct = matches[queryHits(findOverlaps (matches, DAP_res_l[[celltype]]))]
  matchesMat = assay (matches_ct)
  colnames (matchesMat) = gsub ('_.*','',colnames (matchesMat))
  colnames (matchesMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames (matchesMat))
  matchesMat = matchesMat[,infl_TF]
  matchesMat = matchesMat[rowSums (matchesMat) > 0,]
  
  cooc = matrix (ncol = ncol(matchesMat), nrow= ncol(matchesMat))
  ov_size = matrix (ncol = ncol(matchesMat), nrow= ncol(matchesMat))
  
  for (i in 1:ncol(matchesMat))
    {
    for (z in 1:ncol(matchesMat)) 
      {
      #if (sum(rowSums (matchesMat[,c(i,z)]) == 2) < 10)
     #   {
      #  cooc[i,z] = 0
      #  } else {
        ov = sum (rowSums (matchesMat[,c(i,z)]) == 2) / min (colSums(matchesMat[,c(i,z)]))
        cooc[i,z] = ov
        ov_size[i,z] = sum(rowSums (matchesMat[,c(i,z)]) == 2)
      #  }      
      }
    }
  
  colnames (cooc) = infl_TF
  rownames (cooc) = infl_TF
  diag (cooc) = 0
  cooc_l[[celltype]] = cooc#[rowSums(cooc) >0,rowSums(cooc) >0]
  colnames (ov_size) = infl_TF
  rownames (ov_size) = infl_TF
  ov_size_l[[celltype]] = ov_size
  }
#cooc_l[[celltype]][lower.tri (cooc_l[[celltype]])]
cooc_diff = cooc_l[[1]] - cooc_l[[2]]
cooc_diff[is.na(cooc_diff)] = 0
ov_size_max = pmin(ov_size_l[[1]], ov_size_l[[2]])
diag(ov_size_max) = 0

cooc_hm = Heatmap (
    cooc_diff,
    column_names_rot =45, 
    row_names_gp = gpar(fontsize = 3),
    column_names_gp = gpar(fontsize = 3),
    col =palette_cooccurrence_cor, 
    cluster_rows=T,
    cluster_columns = T#,
#rect_gp = gpar (col = "white", lwd = 1)
)

pdf (file.path ('Plots','selected_TF_cooccurence_heatmaps_diff.pdf'), width = 8,height=6)
cooc_hm 
dev.off()

## Try with scatterplot 
cooc_diff_df = as.data.frame (cooc_diff)
cooc_diff_df$TF2 = rownames (cooc_diff_df)
cooc_diff_df = gather (cooc_diff_df, TF, overlap, 1:(ncol(cooc_diff_df)- 1))

# overlap size
ov_size_max_df = as.data.frame(ov_size_max)
ov_size_max_df$TF2 = rownames(ov_size_max_df)
ov_size_max_df = gather (ov_size_max_df, TF, size, 1:(ncol(ov_size_max_df)- 1))

ov_diff_size_df = cbind (cooc_diff_df, ov_size_max_df)
ov_diff_size_df$TF_pair = paste0(ov_diff_size_df$TF2,'_', ov_diff_size_df$TF)
ov_diff_size_df = ov_diff_size_df[, c(3,6,7)]
ov_diff_size_df$topTFs = ifelse (ov_diff_size_df$overlap > 0.5 & ov_diff_size_df$size > 1000, ov_diff_size_df$TF_pair, '') 
sp = ggplot (ov_diff_size_df, aes (x = overlap, y = size, label = topTFs)) + 
geom_point() + 
geom_text_repel (size=2) +
gtheme_no_rot

pdf (file.path ('Plots','peak_overlap_inflamed_scatterplot.pdf'))
sp
dev.off()


# Compare cNMF modules with inflammatory program in scATAC-seq and scRNA-seq ####
shared_cnmf = readRDS (file.path('..','scrna','shared_cnmf_myeloid.rds'))
shared_cnmf = lapply (shared_cnmf, function(x) x[x %in% getFeatures (archp)])
#remove_modules = c('cnmf.3','cnmf.6','cnmf.7','cnmf.5') # remove monocyres cDC and CC modules. Consider re-inculding CC 

pdf (file.path ('Plots','scrna_celltype_dimplot.pdf'))
DimPlot (srt, group.by = 'celltype', reduction = 'umap')
dev.off()

srt = ModScoreCor (
    seurat_obj = srt, 
    geneset_list = shared_cnmf, 
    cor_threshold = NULL, 
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'shared_cnmf', outdir = NULL)




# TF_modules = c(
# 'JUNB
# FOSL2
# JUN
# SMARCC1
# FOSL1
# JUND
# FOS
# JDP2
# BACH1
# FOSB')
# TF_modules = strsplit(TF_modules, '\n')
# names (TF_modules) = 'AP1'



# Differential Peaks in cells positive for km2 vs rest ####
# Find DAP ####
#force = FALSE
archp$inflamed = ifelse (archp$mod_2 > 0, 'inflamed','non_inflamed')
metaGroupName = 'inflamed'
force=F
if (!file.exists ('DAP_inflamed_pairwise.rds') | force)
  {
  DAP_list = getMarkerFeatures (
    ArchRProj = archp, 
    testMethod = "wilcoxon",
          useGroups = "inflamed",
          bgdGroups = "non_inflamed",
    k=100,
    binarize = FALSE,
    useMatrix = "PeakMatrix",
    groupBy = metaGroupName
  #  useSeqnames="z"
  )
  saveRDS (DAP_list, 'DAP_inflamed_pairwise.rds')
  } else {
  DAP_list = readRDS ('DAP_inflamed_pairwise.rds')
  }
DAP_res = do.call (cbind, (assays(DAP_list)))
colnames (DAP_res) = names(assays(DAP_list))
DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames (DAP_res) = as.character(DAP_res_regions)

pdf(file.path('Plots','inflamed_MA_plot.pdf'), width=5,height=5)
pma <- markerPlot (seMarker = DAP_list, name = 'inflamed', cutOff = "FDR <= 0.01", plotAs = "MA")
pma
dev.off()
  
# Take only significant regions ####
DAP_res_sig = DAP_res[DAP_res$FDR < .01 & DAP_res$Log2FC > 0, ]
saveRDS (GRanges(rownames(DAP_res_sig)), 'inflamed_peaks.rds')

### Perform enrichment on DAP ####
archp = addBgdPeaks (archp, force= T)
archp = addMotifAnnotations (ArchRProj = archp, 
      motifSet = "cisbp", 
      #motifSet = 'JASPAR2020',
      #name = "JASPAR2020_Motif",
      force=T)
enrichMotifs <- peakAnnoEnrichment(
    seMarker = DAP_list,
    ArchRProj = archp,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
pdf (file.path ('Plots','enrich_inflamation_heatmap.pdf'))
heatmapEM
dev.off()

# Run GSEA enrichment analysis #### 
#!!! To run this analysis load only ArcHR and clusterprofiler packages !!!!
library (fgsea)    
options(warn = 0)
ps = getPeakSet (archp)

gmt_annotations = c(
'h.all.v7.4.symbol.gmt',#,
'c5.bp.v7.1.symbol.gmt',
'c3.tft.v7.1.symbol.gmt'
)

gmt.file = paste0 ('../../git_repo/files/h.all.v7.4.symbols.gmt')
gmt.file = paste0 ('../../git_repo/files/c5.bp.v7.1.symbol.gmt')
# gmt.file = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_DN.v2024.1.Hs.gmt'
# gmt.file = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP.v2024.1.Hs.gmt'
# gmt.file = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/GSE24026_PD1_LIGATION_VS_CTRL_IN_ACT_TCELL_LINE_UP.v2024.1.Hs.gmt'
# csv.file = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/Tcell_exhastuion_genes_PMID37091230.csv'
pathways = clusterProfiler::read.gmt (gmt.file)
pathways = pathways[grep('inflammatory', pathways$term,ignore.case=T),]
#pathways = read.csv (csv.file)
pathways = split(pathways$gene, pathways$term)
#pathways = gmtPathways (gmt.file)
DAP_list = readRDS ('DAP_inflamed_pairwise.rds')
DAP_res = do.call (cbind, (assays(DAP_list)))
colnames (DAP_res) = names(assays(DAP_list))
DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames(DAP_res) = as.character(DAP_res_regions)

peak_genes = unname(ps[queryHits (findOverlaps(ps, GRanges (rownames(DAP_res))))]$nearestGene)
names (peak_genes) = as.character(ps)[queryHits(findOverlaps (ps,GRanges (rownames(DAP_res))))]
peak_genes = peak_genes[rownames(DAP_res)]
#peak_genes = setNames (-log10(DAP_res$Pval) * sign (DAP_res$Log2FC), peak_genes)
peak_genes = setNames (DAP_res$Log2FC, peak_genes)
#peak_genes = peak_genes[!duplicated(names(peak_genes))]
#names (peak_genes) = 
peak_genes = peak_genes[!duplicated(names(peak_genes))]
peak_genes = peak_genes[!is.na(names(peak_genes))]
library (BiocParallel)
BiocParallel::register(BiocParallel::SerialParam())

### fgsea throws a BiocParallel error when I load all packages including clusterProfiler...try avoiding loading packages except ArchR and fgsea
#peak_genes2 = setNames(order(peak_genes), names(peak_genes))
fgseaRes = fgseaMultilevel (pathways, 
          peak_genes,#, 
          minSize=15, 
        #  scoreType='pos',
          maxSize=1500,
          nproc=1,
          nPermSimple=100000,
          BPPARAM = NULL
          )
pvalAdjTrheshold = 0.05
top_pathways=5
fgseaRes$padj = fgseaRes$pval
fgseaResAll_dp = dotGSEA (
  list(fgseaRes), 
  padj_threshold = pvalAdjTrheshold, 
  type = 'fgsea',
  top_pathways = top_pathways,
  cluster_rows=F,
  cluster_cols=F)

pdf (file.path ('Plots','fgsea_dotplot.pdf'), width=7, height=3)
fgseaResAll_dp
dev.off()


pdf (file.path ('Plots','GO_Inflammatory_enrichment_plot.pdf'), width=5, height=3)
plotEnrichment(pathways[["GO_INFLAMMATORY_RESPONSE"]],
               peak_genes) + labs(title="GO_INFLAMMATORY_RESPONSE")
dev.off()


















## Run peak2genes results with hubs links ####
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
    

# Import hubs from myeloid analysis ####
metaGroupName = "Clusters_H"
cor_cutoff = 0.3
#max_dist = 12500
max_dist = 12500
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
hubs_obj = readRDS (file.path(hubs_dir,'global_hubs_obj.rds'))


# Generate matrix of fragment counts of hubs x barcodes ####
if (!file.exists(file.path (hubs_dir, paste0('hubs_cells_mat.rds'))))
  {
  if (!exists ('fragments')) fragments = unlist (getFragmentsFromProject (
    ArchRProj = archp))    
  hubsCell_mat = matrix (ncol = length(rownames(archp@cellColData)), nrow = length(hubs_obj$hubsCollapsed))
  colnames (hubsCell_mat) = rownames(archp@cellColData)
  rownames (hubsCell_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (rownames(archp@cellColData)))
  for (cell in rownames(archp@cellColData)) 
    {
    pb$tick()  
    fragments_in_cell = fragments[fragments$RG %in% cell]  
    fragments_in_cell_in_hubs = countOverlaps (hubs_obj$hubsCollapsed, fragments_in_cell)
    hubsCell_mat[,cell] = fragments_in_cell_in_hubs
    }
  all (colnames (hubsCell_mat) == rownames(archp@cellColData))  
  hubsCell_mat = t(t(hubsCell_mat) * (10^6 / archp$nFrags)) # scale
  saveRDS (hubsCell_mat, file.path (hubs_dir,paste0('hubs_cells_mat.rds')))
  } else {
  hubsCell_mat = readRDS (file.path (hubs_dir,paste0('hubs_cells_mat.rds')))  
  }
hubsCell_mat = as.data.frame (hubsCell_mat)

all (colnames(hubsCell_mat) == rownames(archp@cellColData))
#ha = HeatmapAnnotation (fetal = archp$fetal, which='row')
# hm = Heatmap (
#   scale (t(hubsCell_mat)), 
#  # left_annotation = ha, 
#   column_names_gp = gpar(fontsize = 3),
#   row_names_gp = gpar(fontsize = 0),
#   show_column_dend = T,
#   #column_km = 5,
#   #row_dend_width = unit(5,'mm'),
#   row_dend_side = 'left',
#   col = rev(palette_hubs_accessibility),
#   border=T,
#   name = 'Hubs')
# pdf (file.path (hubs_dir,'Plots',paste0('hubs_cells_',metaGroupName,'_heatmap.pdf')), height=2.2, width = 5)
# hm
# dev.off()



### Compute differential hub analysis 
# Compute differential hub accessibility DHA ####
library (presto)
metaGroupName = 'celltype2'
all (colnames(hubsCell_mat) == rownames(archp@cellColData))
metagroup = as.character (archp@cellColData[,metaGroupName])
res = wilcoxauc (log2(hubsCell_mat+1), metagroup)
res = res[res$logFC > 0,]

res_l = lapply (split (res, res$group), function(x){
  tmp = x[order (x$padj),]
  tmp
})

res_df = do.call (rbind, res_l)
res_df$gene = hubs_obj$hubsCollapsed$gene[match(res_df$feature, hubs_obj$hubs_id)]
head (res_df[res_df$group == 'TREM2',],20)

top_hubs = 5
res_df_top = res_df %>% group_by (group) %>%
  slice_head(n = top_hubs)

levels_order = unique(res_df_top$group)
metagroup = factor (metagroup, levels = levels_order, ordered=T)


# Generate matrix of fragment counts of hubs x metagroup ####
metaGroupName = 'celltype2'
if (!file.exists(file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds'))))
  {
  if (!exists ('fragments')) fragments = unlist (getFragmentsFromProject (
    ArchRProj = archp))   
  hubsSample_mat = matrix (ncol = length(unique(archp@cellColData[,metaGroupName])), nrow = length(hubs_obj$hubsCollapsed))
  colnames (hubsSample_mat) = unique(archp@cellColData[,metaGroupName])
  rownames (hubsSample_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (unique(archp@cellColData[,metaGroupName])))
  for (sam in unique(archp@cellColData[,metaGroupName]))
    {
    pb$tick()  
    fragments_in_sample = fragments[fragments$RG %in% rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) == sam]]  
    fragments_in_sample_in_hubs = countOverlaps (hubs_obj$hubsCollapsed, fragments_in_sample)
    hubsSample_mat[,sam] = fragments_in_sample_in_hubs
    }
  frags_in_sample = sapply (unique(archp@cellColData[,metaGroupName]), function(x) sum (archp$ReadsInTSS[as.character(archp@cellColData[,metaGroupName]) == x]))
  hubsSample_mat = t(t(hubsSample_mat) * (10^6 / frags_in_sample)) # scale
  saveRDS (hubsSample_mat, file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))
  } else {
  hubsSample_mat = readRDS (file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))  
  }
hubsSample_mat = as.data.frame (hubsSample_mat)

#as.character(archp@cellColData[,metaGroupName])[order(metagroup, na.last=T)]
#ha = HeatmapAnnotation (celltype = metagroup[order (metagroup, na.last=T)])
#hubs_mat
dah_hm = t(scale(t(log2(hubsSample_mat[res_df_top$feature,unique(res_df_top$group)]))))
#dah_hm[dah_hm > 2] = 2
#dah_hm[dah_hm < -2] = -2
hm = Heatmap (dah_hm,
  cluster_rows=F,
#  top_annotation = ha,
  cluster_columns=F,
  row_labels = res_df_top$gene,
  col = rev(palette_hubs_accessibility),
  column_names_gp= gpar (fontsize=11),
  row_names_gp= gpar (fontsize=9),
  column_names_rot=45,
  border=T)

pdf (file.path ('Plots','top_DAH_cnmf_celltypes_heatmap.pdf'), width=5)
hm
dev.off()


#TF = 'SNAI1'
TF = sapply (unique(res_df_top$gene), function(x) unlist(strsplit(x, '-'))[1])
metaGroupName = 'celltype2'

# Compare inflamed vs non-inflamed Momacs ####
library (presto)
archp$inflamed = ifelse (archp$mod_2 > 0, 'inflamed','non_inflamed')
metaGroupName = 'inflamed'
all (colnames(hubsCell_mat) == rownames(archp@cellColData))
metagroup = as.character (archp@cellColData[,metaGroupName])
res = wilcoxauc (log2(hubsCell_mat+1), metagroup)
res = res[res$logFC > 0,]

res_l = lapply (split (res, res$group), function(x){
  tmp = x[order (x$padj),]
  #tmp = x[order (-x$logFC),]
  tmp
})

res_df = do.call (rbind, res_l)
res_df$gene = hubs_obj$hubsCollapsed$gene[match(res_df$feature, hubs_obj$hubs_id)]

top_hubs = 10
res_df_top = res_df %>% group_by (group) %>%
  slice_head(n = top_hubs)

hub = res_df_top$feature
#hub = hubs_obj$hubs_id[grep ('NFKB1', hubs_obj$hubsCollapsed$gene)]
#sample_levels = c('Monocytes','cDCs','SPP1','TREM2','C1Q','IFN','IM')

sample_levels = c('Monocytes','cDCs','SPP1','TREM2','C1Q','IFN','IM')
metaGroupName = 'celltype2'
pdf()
#archp$fetal_sample = paste0(archp$Sample, archp$fetal_group)
#metaGroupName = 'fetal_group'
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp,#[!archp$Sample3 %in% c('P11_HOX')], 
    sample_levels = sample_levels, 
    hubs_regions = hubs_obj$hubsCollapsed,
    #ylim = c(0,0.30),
    groupBy = metaGroupName, 
    #sample_levels = sample_sarc_order,
    minCells = 10,
    #geneSymbol = TF,
    plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
    #pal = palette_sample,
    #pal = palette_fetal,
    threads=1,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    region = ext_range (GRanges (hubs_obj$hubsCollapsed[match(hub[1], hubs_obj$hubs_id)]),10000,10000),
    #upstream = 100000,
    #downstream = 100000,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    pal = palette_myeloid,
    #loops = getCoAccessibility (archp, corCutOff = 0.3,
    #  returnLoops = TRUE),
    useGroups= NULL
)
dev.off()
plotPDF (meso_markers, ArchRProj = archp,height=3.5, width=8, name =paste0('MPM_markers_inflamed_coveragePlots.pdf'),addDOC=F)


metaGroupName = 'shared_cnmf2_r_max'
hub_gr = ext_range (GRanges (hubs_obj$hubsCollapsed[match(hub, hubs_obj$hubs_id)]),100000,100000)
genes_in_region = unique(getPeakSet (archp)[subjectHits (findOverlaps (hub_gr, getPeakSet (archp)))]$nearestGene)
genes_in_region = c(genes_in_region)
genes_in_region = 'CD44'
top_dah = data.frame (
  gene = colMeans (srt@assays$RNA@data[rownames(srt) %in% genes_in_region,,drop=F]),
  group = srt@meta.data[,metaGroupName])
top_dah$group = factor (top_dah$group, levels = rev(sample_levels))
top_dah = na.omit(top_dah)
bp = ggplot (top_dah, aes (x = gene, y = group, fill = group)) + 
vlp + 
bxpv + 
scale_fill_manual (values = palette_myeloid) +
#geom_point (position='identity', alpha=.3, color="grey44", size=1) +
gtheme_no_rot

pdf (file.path ('Plots', paste0('scrna_region_boxplots.pdf')), height=4, width=4)
bp
dev.off()


# Check footprint across celltypes ####
metaGroupName='cnmf_celltypes'
archp <- addGroupCoverages (ArchRProj = archp, groupBy = metaGroupName)
motifPositions <- getPositions (archp)

motifs <- c('NFKB1','JUNB','FOS','JUND','SPI1','SPIB','MITF','RUNX1','CEBPA','SRF','NFAT5','STAT2','MEF2C','PRDM1','IRF3','IRF8','IRF1','IRF2','IRF9')
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
pal = palette_myeloid,
smoothWindow = 25)
  


# Check TF deviations
TF = 'E2F3'
getFeatures (archp, 'MotifMatrix')[grep (TF, getFeatures (archp, 'MotifMatrix'))]
TF1 = c('z:JUN_143','z:FOS_137','z:NFKB1_719','z:FOXM1_352','z:TFDP1_310','z:E2F3_313')

pdf ()
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "MotifMatrix",
    name = TF1, 
    useSeqnames='z',
    pal = rev (palette_deviation),    
    embedding = "UMAP_H",
    imputeWeights = NULL
    )
dev.off()

pdf(file.path('Plots','avg_AP1_deviation_fplot.pdf'),15,15)
wrap_plots(TF_p)
wrap_plots(umap_p1)
wrap_plots (umap_p0,umap_p2)#,umap_p3)
dev.off()

pdf (file.path ('Plots','FRIP_umap.pdf'))
plotEmbedding (
    ArchRProj = archp,
    colorBy = "cellColData",
    name = 'FRIP', 
    useSeqnames='z',
    pal = rev (palette_deviation),    
    embedding = "UMAP_H",
    imputeWeights = NULL
    )
dev.off()




# Export bigiwg files ####
metaGroupName = 'cnmf_celltypes'
exp_bigwig = TRUE
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



### Compare expression of genes in inflammation peaks vs rest ####
tf_match = getMatches (archp)
bg_peaks = getPeakSet (archp)
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
#ap1_complex = c('JUN','FOSB','FOS','BACH1','SMARCC1','FOSL2','JUND','JDP2','BATF')
km1 = names (km$cluster[km$cluster == 1])
km2 = names (km$cluster[km$cluster == 2])
tf_match2 = tf_match[,km2] 
tf_match2 = tf_match2[rowSums(assay(tf_match2)) > 0,]
tf_match1 = tf_match[,km1]
tf_match1 = tf_match1[rowSums(assay(tf_match1)) > 0,]

peakSet1 = rowRanges(tf_match1)[queryHits(findOverlaps(tf_match1, bg_peaks))] 
peakSet2 = rowRanges(tf_match2)[queryHits(findOverlaps(tf_match2, bg_peaks))] 
identical (peakSet1, peakSet2)
#hub_regions = hubs_obj$hubsCollapsed[which(hubs_obj$hubs_id %in% x)]
hub_regions_peaks = bg_peakSet[queryHits(findOverlaps(bg_peakSet, hub_regions))]

metaGroupName = 'celltype2'
pMats = getGroupSE(
  ArchRProj = archp,
  useMatrix = 'PeakMatrix',
  groupBy = metaGroupName,
  divideN = TRUE,
  scaleTo = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)

is_sequential <- function(x) {
  length(x) > 1 && all(diff(x) == 1)
}
peakset = getPeakSet(archp)
pmat_peakset = makeGRangesFromDataFrame(rowData(pMats))
pmat_peakset$nearestGene = peakset$nearestGene
#is_sequential(queryHits (findOverlaps (pmat_peakset, peakset)))

metaGroupName = 'shared_cnmf2_r_max'

for (tf in km2)
  {
  tf_peaks = rowRanges(tf_match[,tf][rowSums(assay(tf_match[,tf]))>0,])
  tf_peaks_fragments = pMats[queryHits (findOverlaps (pmat_peakset, tf_peaks)),]
  pmat_peakset_sub = pmat_peakset[queryHits (findOverlaps (pmat_peakset, tf_peaks)),]
  #tf_peaks = tf_peaks[queryHits (findOverlaps (tf_peaks, makeGRangesFromDataFrame(rowData(pMats))))]
  tf_peaks_fragments = as.data.frame(assay(tf_peaks_fragments))
  
  ps = log2(as.data.frame (AverageExpression (srt, 
  features = unique(pmat_peakset_sub$nearestGene), 
  group.by = metaGroupName)[[1]]) +1)

  sapply (rownames(ps), function(x) cor(ps[x,], 
    na.omit(tf_peaks_fragments[unname(pmat_peakset_sub$nearestGene) == x,]))
  }
min_exp = .1




















# ## Add column on DAM heatmap showing if TF is pioneer or not from chrombpnet ####
# ## Show barplots of top TF occurrence using finemo chrombpnet outputs ####

# ### Compare TF expression from scRNA and inferred by chrombpnet per cell type ####
# library (httr)
# library (XML)
# library (igraph)
#BiocManager::install("universalmotif")
library ('universalmotif')

metaGroupName = 'inflamed'
if (!any (ls() == 'mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData(mSE)$name
#mMat_mg = mMat[DAM_df$gene, ]
# mMat_mg = as.data.frame (t(mMat))
# mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
# mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
# rownames (mMat_mg) = mMat_mg[,1]
# mMat_mg = mMat_mg[,-1]


# #Get active genes from RNA
# metaGroupName = 'celltype_simplified2'
# ps = log2(as.data.frame (AverageExpression (srt, 
# features = colnames(mMat_mg),
# group.by = metaGroupName)[[1]]) +1)
# min_exp = .1
# #ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
# #active_TFs = rownames(ps)[rowSums(ps) > 0]

# #active_genes = corGSM_MM$MotifMatrix_name[corGSM_MM$cor > 0.1]
# #DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% active_TFs,])    
# mMat_l = as.list (as.data.frame (t(mMat_mg)))
# mMat_l = lapply (mMat_l, function(x) data.frame (dev = x, row.names = colnames(mMat_mg)))
# #mMat_l = lapply (mMat_l, function(x) x[rownames(x) %in% active_TFs,,drop=F])

# metaGroupName = 'celltype_lv1'
chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/chromBPnet'
# metaGroupName = 'celltype_lv1'
# celltypes = unique (archp@cellColData[,metaGroupName])

# tf_database = read_meme('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/HOCOMOCO_db/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme', skip = 0, readsites = FALSE, readsites.meta = FALSE)
# tf_database = unique(unlist(lapply(tf_database, function(x) unlist(strsplit(x@name,'_'))[1])))

# list.files (file.path(chromBPdir, celltypes[3],'no_bias_model'))
chrombpnet_counts = list()
metaGroupName = 'inflamed'
celltypes = unique (archp@cellColData[,metaGroupName])
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_counts[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
  }

# Check overlap NFKB with AP-1 complex
motif_pairs_l = list(
  c('FOS','NFKB1'),
  c('JUN','NFKB1'),
  c('SNAI1','NFKB1'),
  c('ATF1','NFKB1'),
  c('CTCF','NFKB1'),
  c('SPI1','NFKB1'),
  c('KLF12','NFKB1'),
  c('CEBPA','NFKB1'),
  c('RUNX1','NFKB1'),
  c('GABPA','NFKB1'),
  c('ZBT7A','NFKB1'))

peak_overlap = NULL
motif_pairs = list()
for (i in seq_along(motif_pairs_l))
  {
  for (celltype in celltypes)
    {
    peakset = read.table (file.path(chromBPdir,paste0('MACS2_',celltype),paste0(celltype, '_peaks.narrowPeak')))
    colnames(peakset) = c('chr','start','end')
    peakset = makeGRangesFromDataFrame (peakset)
    chrombpnet_counts_gr1 = chrombpnet_counts[[celltype]][grep (motif_pairs_l[[i]][1], chrombpnet_counts[[celltype]]$V4),]
    chrombpnet_counts_gr2 = chrombpnet_counts[[celltype]][grep (motif_pairs_l[[i]][2], chrombpnet_counts[[celltype]]$V4),]
    colnames(chrombpnet_counts_gr1) = c('chr','start','end')
    colnames(chrombpnet_counts_gr2) = c('chr','start','end')
    chrombpnet_counts_gr1 = makeGRangesFromDataFrame (chrombpnet_counts_gr1)
    chrombpnet_peaks_gr1 = peakset[unique(queryHits(findOverlaps(peakset, chrombpnet_counts_gr1)))]
    chrombpnet_counts_gr2 = makeGRangesFromDataFrame (chrombpnet_counts_gr2)
    chrombpnet_peaks_gr2 = peakset[unique(queryHits(findOverlaps(peakset, chrombpnet_counts_gr2)))]
  
    peak_overlap[[celltype]] = sum (countOverlaps (chrombpnet_peaks_gr1,chrombpnet_peaks_gr2) > 0) / 
    min (c(length(chrombpnet_peaks_gr1),length(chrombpnet_peaks_gr2)))
    }
  motif_pairs[[i]] = peak_overlap
  }
 motif_pairs_l2 = unlist(lapply (motif_pairs_l, function(x) paste (x, collapse='_')))
 names (motif_pairs) = motif_pairs_l2
#motif_pairs = unlist (motif_pairs_l, recursive=F)


# Intersect DAP with MACS2 peaks and look at chrombpnet predicted TFs inflamed vs not ####
DAP_list = readRDS (paste0('DAP_',metaGroupName,'.rds'))
DAP_res = do.call (cbind, (assays(DAP_list)))
colnames (DAP_res) = names(assays(DAP_list))
DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames(DAP_res) = as.character(DAP_res_regions)
DAP_res = DAP_res[!is.na(DAP_res$FDR),]
DAP_res = DAP_res[DAP_res$FDR < 0.05,]
DAP_res_l = list(up = GRanges(rownames(DAP_res)[DAP_res$Log2FC > 0]),
                down = GRanges(rownames(DAP_res)[DAP_res$Log2FC < 0]))
celltypes = c('inflamed','non_inflamed')
names (DAP_res_l) = celltypes

chrombpnet_tfs_DAP_l = list()
ov_mat_l = list()
for (celltype in celltypes)
  {
  peakset = read.table (file.path(chromBPdir,paste0('MACS2_',celltype),paste0(celltype, '_peaks.narrowPeak')))
  colnames(peakset) = c('chr','start','end')
  peakset = makeGRangesFromDataFrame (peakset)
  peakset = peakset[unique(queryHits(findOverlaps(peakset, DAP_res_l[[celltype]])))]
  chrombpnet_counts_gr = chrombpnet_counts[[celltype]]
  colnames(chrombpnet_counts_gr) = c('chr','start','end','TF')
  chrombpnet_counts_gr = makeGRangesFromDataFrame (chrombpnet_counts_gr, keep.extra.columns=T)
  chrombpnet_counts_gr = chrombpnet_counts_gr[queryHits(findOverlaps(chrombpnet_counts_gr, peakset))]
  
  peakset_l = split (peakset, 1:length(peakset))
  ov_mat = sapply (unique(chrombpnet_counts_gr$TF), function(y) countOverlaps (peakset, chrombpnet_counts_gr[chrombpnet_counts_gr$TF == y]))
  rownames(ov_mat) = as.character(peakset)
  ov_mat_l[[celltype]] = ov_mat
  chrombpnet_tfs_DAP_l[[celltype]] = chrombpnet_counts_gr
  }
table (chrombpnet_tfs_DAP_l[[1]]$TF)[order(-table (chrombpnet_tfs_DAP_l[[1]]$TF))]
table (chrombpnet_tfs_DAP_l[[2]]$TF)[order(-table (chrombpnet_tfs_DAP_l[[2]]$TF))]

ov_mat_cor = lapply (ov_mat_l, function(x) cor (x))

# Filter using only TF from inflammation module
km = readRDS ('TF_activity_modules.rds')
infl_TF = names(km$cluster[km$cluster == 2])
ov_mat_cor_fl = lapply (ov_mat_cor, function(x) {
  tfmatch = unique(unlist(sapply (infl_TF, function(y) (grep(y, rownames(x))))))
  x[tfmatch, tfmatch]
})
pdf (file.path ('Plots','chrombpnet_TF_cor_in_peaks.pdf'),height=4,width=4.5)
lapply (ov_mat_cor_fl, function (x) Heatmap (
  x, col = palette_cooccurrence_cor_fun,
  row_names_gp = gpar (fontsize = 6),
  column_names_gp = gpar (fontsize = 6), 
border=T))
dev.off()

# Check expression of FOS JUNB and JUND
genes = c('FOS','JUNB','JUND','NFKB1','NFKB2','SPI1','SPIB','CEBPA','SRF','MITF')
pdf (file.path ('Plots','expression_FOS_JUNB_JUND.pdf'))
DotPlot (srt, features = genes, group.by = 'shared_cnmf2_r_max') + gtheme
dev.off()

# Assume your data is in ov_mat (rows = regions, columns = TFs)
hm_l = list()
for (celltype in celltypes)
  {
  # Ensure it's a binary matrix (0s and 1s)
  ov_mat_bin <- ov_mat_l[[celltype]] > 0
  
  # Convert to a matrix if it's a data.frame or tibble
  ov_mat_bin <- as.matrix(ov_mat_bin)
  
  # Initialize an empty matrix to store Jaccard indices
  n <- ncol(ov_mat_bin)
  jaccard_matrix <- matrix(0, nrow = n, ncol = n)
  colnames(jaccard_matrix) <- colnames(ov_mat_bin)
  rownames(jaccard_matrix) <- colnames(ov_mat_bin)
  
  # Compute Jaccard index for each pair of TFs
  for (i in 1:n) {
    for (j in i:n) {
      a <- ov_mat_bin[, i]
      b <- ov_mat_bin[, j]
      intersection <- sum(a & b)
      union <- sum(a | b)
      jaccard <- ifelse(union == 0, NA, intersection / union)
      jaccard_matrix[i, j] <- jaccard
      jaccard_matrix[j, i] <- jaccard  # symmetry
    }
  }
  hm_l[[celltype]] = Heatmap (jaccard_matrix, column_names_gp = gpar(fontsize = 4), row_names_gp = gpar(fontsize = 4))
  }


pdf (file.path ('Plots','inflamed_overlap_heatmap2.pdf'), width=5, height=5)
hm_l
dev.off()


### Make barplots of most abundant TFs identified in inflamed and non-inflamed cells
bp_df = data.frame (
  Freq = c(proportions(head(table (chrombpnet_tfs_DAP_l[[1]]$TF)[order (-table (chrombpnet_tfs_DAP_l[[1]]$TF))],5)),
proportions(head (table (chrombpnet_tfs_DAP_l[[2]]$TF)[order (-table (chrombpnet_tfs_DAP_l[[2]]$TF))],5))),
  TF = names (c(head(table (chrombpnet_tfs_DAP_l[[1]]$TF)[order (-table (chrombpnet_tfs_DAP_l[[1]]$TF))],5),
head (table (chrombpnet_tfs_DAP_l[[2]]$TF)[order (-table (chrombpnet_tfs_DAP_l[[2]]$TF))],5))),
  type = c(rep('inflamed',5), rep('noninflamed',5)))


library(dplyr)

df_ordered <- bp_df %>%
  group_by(type) %>%
  arrange(desc(Freq), .by_group = TRUE) %>%
  ungroup()
  
inflTF_palette = c(
  FOS_JUNB_JUNDinflamed = 'red', 
  SPI1_SPIB_IRF8inflamed='#007FFFFF',
  FIGLA_MESP1_SNAI1inflamed = '#FFEFB2FF',
  KLF12_SP1_SP1inflamed = '#001933FF',
  CEBPB_CEBPD_CEBPAinflamed = '#A89797FF',
  SPI1_SPIB_ELK1noninflamed = '#007FFFFF', 
  SNAI1_SNAI2_FIGLAnoninflamed = '#FFEFB2FF',
  KLF12_SP1_SP1noninflamed = '#001933FF',
  ATF1_CREB1_CREB5noninflamed = '#7FBFFFFF', 
  SPI1_SPIB_IRF4noninflamed = '#007FFFFF')


df_ordered$TF_type = paste0(df_ordered$TF, df_ordered$type)
df_ordered$TF_type = factor (df_ordered$TF_type, levels = unique (df_ordered$TF_type))
bp = ggplot (df_ordered, aes (x = type, y = Freq, fill = TF_type)) +
  geom_bar (stat = 'identity', position = 'stack') + 
  scale_fill_manual (values = inflTF_palette) + gtheme

pdf (file.path ('Plots', 'TF_abundance_inflamed_noninflamed_barplot.pdf'),4,width=4.5)
bp
dev.off()




chrombnet_counts_2 = list()
for (celltype in celltypes)
  {
  chrombpnet_counts_tmp = chrombpnet_counts[[celltype]]
  chrombnet_counts_2[[celltype]] = table (chrombpnet_counts_tmp$V4)[order(-table (  chrombpnet_counts_tmp$V4))]
  }


chrombpnet_profile = list()
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  }
  #chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  #chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  #chrombpnet_tf = rbind (chrombpnet_count_tf, chrombpnet_profile_tf)

#ap1_complex = c('JUN','FOSB','FOS','BACH1','SMARCC1','FOSL2','JUND','JDP2','BATF','CEBPB','CEBPA','CEBPZ','FOSL1','NFE2','NFE2L2','NFE2L1')

chrombpnet_profile_2 = list()
for (celltype in celltypes)
  {
  chrombpnet_profile_tmp = chrombpnet_profile[[celltype]]
  chrombpnet_profile_2[[celltype]] = table (chrombpnet_profile_tmp$V4)[order(-table (  chrombpnet_profile_tmp$V4))]
  }  
  # #assign_max_exp = unlist(sapply (names(chrombpnet_counts_tmp), function(x) unlist(strsplit(x, '_'))[which.max(ps[unlist(strsplit(x, '_')), celltype])]))
  # #tf_dev = mMat_l[[celltype]][assign_max_exp,]
  # chrombpnet_counts_tmp2 = data.frame (occurrence =   chrombpnet_counts_tmp[names(assign_max_exp)], TF_max_exp = assign_max_exp, TF_max_dev = tf_dev)
  # #chrombpnet_counts_tmp2 = chrombpnet_counts_tmp2[!chrombpnet_counts_tmp2$TF_max_exp %in% ap1_complex,]
  # chrombpnet_counts_tmp2 = chrombpnet_counts_tmp2[!duplicated(chrombpnet_counts_tmp2$TF_max_exp),]
  # chrombpnet_counts_tmp2 = chrombpnet_counts_tmp2[chrombpnet_counts_tmp2$TF_max_exp %in% head(chrombpnet_counts_tmp2$TF_max_exp[order(-chrombpnet_counts_tmp2$TF_max_dev)],10),]
  # chrombpnet_counts_tmp2$celltype = celltype
  # chrombpnet_counts_tmp2$order = seq(nrow(chrombpnet_counts_tmp2))
  # chrombpnet_counts2[[celltype]] = chrombpnet_counts_tmp2
  # }

chrombpnet_counts_df = do.call (rbind, chrombpnet_counts2)
chrombpnet_counts_df = chrombpnet_counts_df %>% group_by (celltype) %>% mutate(Proportion = occurrence.Freq / sum(occurrence.Freq))
chrombpnet_counts_df$TF_max_exp2 = chrombpnet_counts_df$TF_max_exp
chrombpnet_counts_df$TF_max_exp[chrombpnet_counts_df$Proportion < 0.05] = ''
chrombpnet_counts_df$TF_max_exp = factor (chrombpnet_counts_df$TF_max_exp, levels =unique(chrombpnet_counts_df$TF_max_exp))
chrombpnet_counts_df$order = factor (chrombpnet_counts_df$order, levels =unique(chrombpnet_counts_df$order))
# Create stacked bar plot with text beside each band
bp = ggplot (chrombpnet_counts_df, aes(x = celltype, y = Proportion, fill = order)) +
  geom_bar (stat = "identity", color = 'white') +
  geom_text (aes(label = TF_max_exp), 
            position = position_stack (vjust = 0.5), 
            hjust = 0.5,  # Move text outside the bar
            size = 3) + 
  #coord_flip() +  # Flip to make text more readable
  gtheme
pdf (file.path ('Plots','chrombpnet_counts_TF_barplot.pdf'), width=7, height=4)
bp
dev.off()


chrombpnet_profile = list()
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  }
  #chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  #chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  #chrombpnet_tf = rbind (chrombpnet_count_tf, chrombpnet_profile_tf)
chrombpnet_profile2 = list()
for (celltype in celltypes)
  {
  chrombpnet_profile_tmp = chrombpnet_profile[[celltype]]
  chrombpnet_profile_tmp = table (  chrombpnet_profile_tmp$V4)[order(-table (  chrombpnet_profile_tmp$V4))]
  assign_max_exp = unlist(sapply (names(chrombpnet_profile_tmp), function(x) unlist(strsplit(x, '_'))[which.max(ps[unlist(strsplit(x, '_')), celltype])]))
  tf_dev = mMat_l[[celltype]][assign_max_exp,]
  chrombpnet_profile_tmp2 = data.frame (occurrence =   chrombpnet_profile_tmp[names(assign_max_exp)], TF_max_exp = assign_max_exp, TF_max_dev = tf_dev)
  #chrombpnet_profile_tmp2 = chrombpnet_profile_tmp2[!chrombpnet_profile_tmp2$TF_max_exp %in% ap1_complex,]
  chrombpnet_profile_tmp2 = chrombpnet_profile_tmp2[!duplicated(chrombpnet_profile_tmp2$TF_max_exp),]
  chrombpnet_profile_tmp2 = chrombpnet_profile_tmp2[chrombpnet_profile_tmp2$TF_max_exp %in% head(chrombpnet_profile_tmp2$TF_max_exp[order(-chrombpnet_profile_tmp2$TF_max_dev)],10),]
  chrombpnet_profile_tmp2$celltype = celltype
  chrombpnet_profile_tmp2$order = seq(nrow(chrombpnet_profile_tmp2))
  chrombpnet_profile2[[celltype]] = chrombpnet_profile_tmp2
  }

chrombpnet_profile_df = do.call (rbind, chrombpnet_profile2)
chrombpnet_profile_df = chrombpnet_profile_df %>% group_by (celltype) %>% mutate(Proportion = occurrence.Freq / sum(occurrence.Freq))
chrombpnet_profile_df$TF_max_exp2 = chrombpnet_profile_df$TF_max_exp
chrombpnet_profile_df$TF_max_exp[chrombpnet_profile_df$Proportion < 0.05] = ''
chrombpnet_profile_df$TF_max_exp = factor (chrombpnet_profile_df$TF_max_exp, levels =unique(chrombpnet_profile_df$TF_max_exp))
chrombpnet_profile_df$order = factor (chrombpnet_profile_df$order, levels =unique(chrombpnet_profile_df$order))
# Create stacked bar plot with text beside each band
bp = ggplot (chrombpnet_profile_df, aes(x = celltype, y = Proportion, fill = order)) +
  geom_bar (stat = "identity", color = 'white') +
  geom_text (aes(label = TF_max_exp), 
            position = position_stack (vjust = 0.5), 
            hjust = 0.5,  # Move text outside the bar
            size = 3) + 
  #coord_flip() +  # Flip to make text more readable
  gtheme
pdf (file.path ('Plots','chrombpnet_profile_TF_barplot.pdf'), width=7, height=4)
bp
dev.off()


### Check if AP1 peaks have less correlated gene expression nearby ####
###-- Annotate hubs using p2g links ---###
# maxDist = 500000
# force=F
# if (!file.exists (paste0('p2g_links_',maxDist,'.rds')) | force)
#   {
#   archp = addPeak2GeneLinks (
#     ArchRProj = archp,
#     maxDist = maxDist,
#     reducedDims = "IterativeLSI",
#     overlapCutoff = 0.5,
#     #cellsToUse = metaGroup_df$barcode
#     )
#   saveRDS (p2g_links, paste0('p2g_links_',maxDist,'.rds'))
#   } else {
#   p2g_links = readRDS (paste0('p2g_links_',maxDist,'.rds'))
#   }

# # Get p2g data.frames
# p2g_cor_threshold = 0.3
# p2g_links = getPeak2GeneLinks (
#     ArchRProj = archp,
#     corCutOff = p2g_cor_threshold,
#     resolution = 1,
#     returnLoops = FALSE
#   )
# p2g_mean_cor = sapply (unique(p2g_links$idxATAC), function(x) mean(abs(p2g_links[p2g_links$idxATAC == x,'Correlation'])))
# p2g_idx = unique(p2g_links$idxATAC)

# motifmatch = getMatches (archp)
# pset = getPeakSet(archp)
# pset$mean_cor = 0
# pset$mean_cor[p2g_idx] = p2g_mean_cor
# tf_peaks_meancor = list()
# for (i in 1:ncol(motifmatch))
#   {
#    tf_peaks_meancor[[i]] = mean (pset$mean_cor[as.logical(assay(motifmatch[,i]))])
#   }
# names (tf_peaks_meancor) = colnames (motifmatch)
# tf_peaks_meancor = unlist(tf_peaks_meancor)
# names (tf_peaks_meancor) = gsub ('_.*','',names (tf_peaks_meancor))
# names (tf_peaks_meancor) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", names (tf_peaks_meancor))
# head (tf_peaks_meancor[names(km$cluster)][order(unlist(tf_peaks_meancor[names(km$cluster)]))],150)




















