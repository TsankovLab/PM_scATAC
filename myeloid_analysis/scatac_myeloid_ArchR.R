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

cell_subsets_order = c('Mono_CD14','Mono_CD16','TAM_CXCLs','TAM_TREM2','TAM_MARCO','cDCs','TAM_interstitial')

# Load RNA ####
srt = readRDS (file.path('..','scrna','srt.rds'))
#sample_levels = c('Monocytes','TREM2','SPP1','cDCs','IM')

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
archp$celltype_lv2 = 0
archp$celltype_lv2[archp$Clusters_H %in% c('C6','C7','C9')] = 'Monocytes'
archp$celltype_lv2[archp$Clusters_H %in% c('C8')] = 'cDCs'
archp$celltype_lv2[archp$Clusters_H %in% c('C1','C2','C3','C4','C5','C10','C11')] = 'TAMs'

pdf()
umap_p3 = plotEmbedding (ArchRProj = archp, labelMeans = F,
  colorBy = "cellColData", name = "Sample",
  pal= palette_sample,
   embedding = "UMAP_H")
umap_p4 = plotEmbedding (ArchRProj = archp, labelMeans = F,
  colorBy = "cellColData", name = "celltype_lv2",
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
write.csv (patchvecs (shared_cnmf), 'cnmf_consensus_myeloid_modules.csv')

remove_modules = c('Stress','TAM_unknown','CC')
shared_cnmf = shared_cnmf[!names(shared_cnmf) %in% remove_modules]

force = T
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
#archp_MAC = archp[!archp$celltype_lv2 %in% c('Monocytes','cDCs')]
cnmf_scatac = as.data.frame (t(scale(t(scale(archp@cellColData[,names(shared_cnmf)])))))
#cnmf_scatac = as.data.frame (scale(archp@cellColData[,names(shared_cnmf)]))
if (!file.exists ('cnmf_scatac.rds')) saveRDS (cnmf_scatac, 'cnmf_scatac.rds')

cap = 3
cnmf_scatac_cap = cnmf_scatac
cnmf_scatac_cap[cnmf_scatac_cap > cap] = cap
cnmf_scatac_cap[cnmf_scatac_cap < -cap] = -cap
set.seed (123)
km_cnmf = kmeans (cnmf_scatac_cap, centers=7) # double scale modules and cluster using k-means
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

pdf (file.path ('Plots','cnmf_scatac_scaled_heatmap.pdf'), height=1.5)
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
#     colorBy = "GeneScoreMatrix", 
#     name = c('C3','FOSB','CSF1R','A2M','NAIP','ADAP2','FOLR2','MALAT1' ),
#     embedding = "UMAP_H",
#     pal = palette_expression,
#     imputeWeights = getImputeWeights(archp)
# )
dev.off()
pdf (file.path ('Plots','shared_cnmf_TAMs_fplots.pdf'),14,14)
wrap_plots (p, ncol=6)
# p2
dev.off()


# Re-Annotate based on cnmf clustering ####
all (names(km_cnmf$cluster) == rownames(archp@cellColData))
archp$cnmf_cluster = paste0('cnmf_cluster_',km_cnmf$cluster)
archp$celltype_lv3 = archp$cnmf_cluster
archp$celltype_lv3[archp$celltype_lv3 == 'cnmf_cluster_1'] = 'Mono_CD16'
archp$celltype_lv3[archp$celltype_lv3 == 'cnmf_cluster_2'] = 'TAM_MARCO'
archp$celltype_lv3[archp$celltype_lv3 == 'cnmf_cluster_3'] = 'Mono_CD14'
archp$celltype_lv3[archp$celltype_lv3 == 'cnmf_cluster_4'] = 'TAM_CXCLs'
archp$celltype_lv3[archp$celltype_lv3 == 'cnmf_cluster_5'] = 'TAM_TREM2'
archp$celltype_lv3[archp$celltype_lv3 == 'cnmf_cluster_6'] = 'TAM_interstitial'
archp$celltype_lv3[archp$celltype_lv3 == 'cnmf_cluster_7'] = 'cDCs'
archp$celltype_lv3[archp$celltype_lv2 == 'cDCs'] = 'cDCs'
#archp$celltype_lv2[archp$celltype_lv2 == 'cnmf_cluster_6'] = 'TREM2'

write.csv (data.frame (barcode = rownames(archp@cellColData), celltype_lv2 = archp$celltype_lv3), 'barcode_annotation.csv')
# Integrate MACs with myeloid annotation
#archp$celltype2[match(rownames(archp_MAC@cellColData), rownames(archp@cellColData))] = archp_MAC$celltype_lv2

pdf()
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype_lv2",
   embedding = "UMAP_H")
dev.off()

pdf (file.path('Plots','celltype_umap_TAM_annotated_umap2.pdf'),5,5)
print (umap_p5)
dev.off()

  # Find DAM ####
  metaGroupName = "celltype_lv3"
  force = F
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
pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_heatmap.pdf')), width = 3, height = 4)
print(DAM_hm)
dev.off()

### Show IRF TFs activity ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames(mMat) = rowData(mSE)$name
mMat = as.data.frame(mMat)
mMat = mMat[rownames(mMat)[grep('IRF',rownames(mMat))], ]

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
pdf (file.path ('Plots',paste0('IRF_TFs_',metaGroupName,'_heatmap.pdf')), width = 3, height = 4)
print(DAM_hm)
dev.off()


# #### Identify TF regulators correlated to each scRNA cnmf #####
# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# #mMat = scale(as.matrix(mMat))#[selected_TF,])

# # Filter by RNA expression ####
# metaGroupName = 'celltype_lv2'
# min_exp = .1
# active_TFs = exp_genes (srt, rownames(mMat), min_exp = min_exp, metaGroupName)
# mMat = t(scale (mMat[active_TFs, ]))
# if (!file.exists ('mMat_scaled_active.rds')) saveRDS (mMat, 'mMat_scaled_active.rds')


# cnmf_scatac = readRDS ('cnmf_scatac.rds')
# mMat = readRDS ('mMat_scaled_active.rds')
# tf_cnmf_cor = cor (mMat, cnmf_scatac, method='spearman')
# cell_subsets = colnames(tf_cnmf_cor)[c(1, 5, 2, 4, 3)]
# top_5 = unique (unlist(lapply (cell_subsets, function(x) head (rownames(tf_cnmf_cor[order(-tf_cnmf_cor[,x]),]),5))))


# DAM_hm = Heatmap (tf_cnmf_cor[top_5, cell_subsets], 
#           #row_labels = colnames (mMat_mg),
#           #column_title = paste('top',top_genes),
#           clustering_distance_columns = 'euclidean',
#           clustering_distance_rows = 'euclidean',
#           cluster_rows = F,
#           #col = pals_heatmap[[5]],
#           cluster_columns=F,#col = pals_heatmap[[1]],
#           row_names_gp = gpar(fontsize = 8),
#           column_names_gp = gpar(fontsize = 8),
#           column_names_rot = 45,
#           name = 'chromVAR',
#           #rect_gp = gpar(col = "white", lwd = .5),
#           border=TRUE,
#           col = rev(palette_deviation)#,
#           #width = unit(2, "cm")
#           #right_annotation = motif_ha
#           )

# pdf (file.path ('Plots','cnmf_clusters_DAM_heatmap.pdf'), width = 3,height=4)
# draw (DAM_hm)
# dev.off()

# ### Try using metacells for each cnmf ####
# sams = c('P1','P10','P11','P12','P13','P14','P23','P5') # Select samples that have at least 100 endothelial cells
# df = list()

# library (zoo)
# bin_width <- 30   # Number of observations per bin
# overlap <- 30
# cnmf_mods = c('Mono','cDCs','TREM2','SPP1','IFN_CXCLs','IM')
# #archp$Sample2 = 'sample'
# #sams = 'sample'

# for (sam in sams)
# {
#   cnmf_l = list()
#   for (cnmf_mod in cnmf_mods)
#   {
#   metacells_order = cnmf_scatac[,cnmf_mod][archp$Sample == sam]
#   metacells_order = order (-metacells_order)
#   mMat_sam = mMat[archp$Sample == sam, ]
#   mMat_sam = mMat_sam[metacells_order,]
#   cnmf_sam = cnmf_scatac[,cnmf_mod][archp$Sample == sam][metacells_order]
  
#  cnmf_l[[cnmf_mod]] <- cor (
#   rollapply (mMat_sam, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"),
#   rollapply (cnmf_sam, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"), 
#   method = 'spearman')  
#   }
# df[[sam]] = as.data.frame (do.call (cbind, cnmf_l))
# colnames (df[[sam]]) = cnmf_mods
# #df[[sam]]$sample = sam
# }

# # Convert list to 3D array: rows x cols x n_matrices
# array_3d <- array(unlist(df), dim = c(nrow(df[[1]]), ncol(df[[1]]), length(df)))

# # Apply median across the third dimension (i.e., across matrices)
# median_matrix <- apply(array_3d, c(1, 2), median)
# rownames (median_matrix) = rownames (df[[1]])
# colnames (median_matrix) = colnames (df[[1]])
# top_5 = unlist(lapply (cnmf_mods, function(x) head (rownames(median_matrix[order(-median_matrix[,x]),]),5)))
# top_5 = unlist(lapply (cnmf_mods, function(x) head (rownames(df[[1]][order(-df[[1]][,x]),]),5)))

# DAM_hm = Heatmap (df[[1]][top_5,cnmf_mods], 
#           #row_labels = colnames (mMat_mg),
#           #column_title = paste('top',top_genes),
#           clustering_distance_columns = 'euclidean',
#           clustering_distance_rows = 'euclidean',
#           cluster_rows = F,
#           #col = pals_heatmap[[5]],
#           cluster_columns=F,#col = pals_heatmap[[1]],
#           row_names_gp = gpar(fontsize = 8),
#           column_names_gp = gpar(fontsize = 8),
#           column_names_rot = 45,
#           name = 'chromVAR',
#           #rect_gp = gpar(col = "white", lwd = .5),
#           border=TRUE,
#           col = rev(palette_deviation)#,
#           #width = unit(2, "cm")
#           #right_annotation = motif_ha
#           )

# pdf (file.path ('Plots','cnmf_clusters_DAM_metacells_heatmap.pdf'), width = 3,height=4)
# draw (DAM_hm)
# dev.off()


#df = tail (df,50)
# Create the scatterplot with a polished aesthetic
# df = do.call (rbind, df)

### Co-expression of TFs across cells #### 

### Run TF correlation to identify TF modules across TNK cells ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
#mMat = scale(as.matrix(mMat))#[selected_TF,])

# # Filter by RNA expression ####
metaGroupName = 'celltype_lv2'
min_exp = .1
active_TFs = exp_genes (srt, rownames(mMat), min_exp = min_exp, metaGroupName)
mMat = t(scale (mMat[active_TFs, ]))
mMat_cor = cor (as.matrix(mMat), method = 'spearman')

set.seed(123)
centers=2
km = kmeans (mMat_cor, centers=centers)
if (!file.exists ('TF_activity_modules.rds')) saveRDS (km, 'TF_activity_modules.rds')
write.csv (patchvecs (split (names(km$cluster),km$cluster)), 'regMye_modules.csv')

genes_highlight = c(
'IKZF1
HIVEP3
HIVEP1
NFKB2
RELB
RELA
NFKB1
REL
NFE2L2
NFE2
BATF
JDP2
FOSB
SMARCC1
JUND
FOSL1
BACH1
JUNB
FOSL2
FOS
JUN')

genes_highlight2 = c(
'NFE2L2
NFE2
BATF
JDP2
FOSB
SMARCC1
JUND
FOSL1
BACH1
JUNB
FOSL2
FOS
JUN')
# saveRDS (genes_highlight,'inflammation_TFs.rds')

genes_highlight = unlist(strsplit(genes_highlight,'\n'))
genes_highlight2 = unlist(strsplit(genes_highlight2,'\n'))

ha2 = rowAnnotation (foo = anno_mark(at = match(genes_highlight,colnames(mMat_cor)), 
    labels = genes_highlight, labels_gp = gpar(fontsize = 7)))

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

pdf (file.path ('Plots','TF_modules_heatmap.pdf'), width = 4.6, height=3)
cor_mMat_hm
dev.off()

# Add metacolumns of average TF modules activity ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
#mMat = scale(as.matrix(mMat))#[selected_TF,])

# # Filter by RNA expression ####
metaGroupName = 'celltype_lv2'
min_exp = .1
active_TFs = exp_genes (srt, rownames(mMat), min_exp = min_exp, metaGroupName)
mMat = t(scale (mMat[active_TFs, ]))

tf_modules = lapply (unique(km$cluster), function(x) colMeans (t((mMat))[rownames(t(mMat)) %in% names(km$cluster[km$cluster == x]),]))
tf_modules = c(tf_modules, 
  list(AP1NFKB1 = colMeans (t((mMat))[rownames(t(mMat)) %in% genes_highlight,]),
  AP1 = colMeans (t((mMat))[rownames(t(mMat)) %in% genes_highlight2,])))
# tf_modules = c(tf_modules, list(AP1 = colMeans (t(scale(mMat))[rownames(t(mMat)) %in% c('JUN','FOSB','FOS','BACH1','SMARCC1','FOSL2','JUND','JDP2','BATF'),])))
#tf_module_infl = colMeans (t)


names (tf_modules)[1:2] = paste0('mod_',unique(km$cluster))
tf_modules = do.call (cbind, tf_modules)
archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% colnames(tf_modules)]
archp@cellColData = cbind (archp@cellColData, tf_modules) 
archp$inflamed = ifelse (archp$mod_2 > 0, 'inflamed','non_inflamed')

pdf()
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "cellColData",
    name = paste0('mod_',c(1:2)),
    pal = rev(palette_deviation_correlation),
    imputeWeights=NULL,
    #useSeqnames='z',
    embedding = "UMAP_H")
dev.off()

pdf (file.path ('Plots','TF_modules_fplots.pdf'), width = 20,height=16)
wrap_plots (TF_p, ncol=5)
dev.off()

# Generate same heatmap but using scrna ####
metacells = readRDS (file.path ('..','scrna','metacells.rds'))
metacells_mat = as.matrix(metacells@assays$RNA$data)[rownames(mMat_cor),]
metacells_mat = cor (t(metacells_mat), method = 'spearman')
all (rownames(metacells_mat) == rownames(mMat_cor))
tf_order = unname(unlist(row_order(cor_mMat_hm)))

pdf (file.path ('Plots','TF_modules_RNA_heatmap.pdf'), width = 4,height=3)
cor_mMat_hm2 = draw (Heatmap (metacells_mat[tf_order,tf_order],# row_km=15,
  row_split = km$cluster,
  column_split = km$cluster,
  cluster_rows=F,
  cluster_columns = F,
  col=palette_expression_correlation, 
  border=T,
  row_names_gp = gpar(fontsize = 0), 
  column_names_gp = gpar(fontsize = 0)))
dev.off()


# ### Plot UMAP using chromvar TF ####
# mMat = readRDS ('mMat_scaled_active.rds')
# library (uwot)
# set.seed(42)  # for reproducibility
# umap_result <- umap(mMat, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
# library(ggplot2)

# umap_df = as.data.frame(umap_result)
# umap_df$celltypes = archp$celltype_lv2
# umap_df$FRIP = log2 (archp$nFrags+1)
# umap_df$Sample = archp$Sample
# umap_df$mod_2 = archp$mod_2
# umap_df = cbind (umap_df, t(scale(t(as.data.frame (archp@cellColData[,names(shared_cnmf)])))))
# sp = lapply (unique(archp$celltype_lv2), function(x) 
#   ggplot (umap_df[umap_df$celltypes == x,], aes(V1, V2, color = celltypes)) +
#   geom_point(size = 2) +
#   labs(title = "UMAP of Iris Data", x = "UMAP1", y = "UMAP2") +
#   theme_minimal())
# sp2 = lapply (names (shared_cnmf), function(x) 
#   ggplot (umap_df, aes_string('V1', 'V2', color = x)) +
#   geom_point(size = .2) +
#   scale_color_viridis_b() +
#   labs(title = "chromvar UMAP cnmf scores", x = "UMAP1", y = "UMAP2") +
#   theme_minimal())
# spq = ggplot (umap_df, aes(V1, V2, color = Sample)) +
#   geom_point(size = 2) +
#   labs(title = "UMAP of Iris Data", x = "UMAP1", y = "UMAP2") +
#   theme_minimal()
# spq2 = ggplot (umap_df, aes(V1, V2, color = mod_2)) +
#   geom_point(size = 2) +
#   scale_color_viridis_b() +
#   labs(title = "UMAP of Iris Data", x = "UMAP1", y = "UMAP2") +
#   theme_minimal()  
# png (file.path ('Plots','chromvar_cnmf_scores_umap.png'), 3000,3000,res=300)
# wrap_plots (sp2)
# dev.off()

# png (file.path ('Plots','chromvar_infl_scores_umap.png'), 3000,3000,res=300)
# spq2
# dev.off()

# pdf (file.path ('Plots','chromvar_umap.pdf'))
# sp
# spq
# dev.off()

# sp1 = lapply (unique(umap_df$Sample), function(y) 
#   {
#   lapply (unique(archp$celltype_lv2), function(x) 
#   ggplot (umap_df[umap_df$celltypes == x & umap_df$Sample == y,], aes(V1, V2, color = celltypes)) +
#   geom_point(size = 2) +
#   labs(title = "UMAP of Iris Data", x = "UMAP1", y = "UMAP2") +
#   theme_minimal())
#   })

# pdf (file.path ('Plots','chromvar_per_sample_umap.pdf'), width=10)
# lapply (sp1, function(x) wrap_plots(x))
# dev.off()



# # Make scatterplots of inflammation vs cnmfs ####
# cnmf_mods = c('Mono','SPP1','TREM2','IFN_CXCLs','IM')
# cnmf_scatac = readRDS ('cnmf_scatac.rds')


# ### Try using metacells for each cnmf ####
# sams = c('P1','P10','P11','P12','P13','P14','P23','P5') # Select samples that have at least 100 myeloid cells
# df = list()
# archp$Sample2 = 'sample'
# sams = 'sample'
# library (zoo)
# bin_width <- 30   # Number of observations per bin
# overlap <- 30
# for (sam in sams)
# {
#   cnmf_l = list()
#   for (cnmf_mod in cnmf_mods)
#   {
#   metacells_order = cnmf_scatac[,cnmf_mod][archp$Sample2 == sam]
#   metacells_order = order (-metacells_order)
#   mMat_sam = archp$mod_2[archp$Sample2 == sam]
#   mMat_sam = mMat_sam[metacells_order]
#   cnmf_sam = cnmf_scatac[,cnmf_mod][archp$Sample2 == sam]
#   cnmf_sam = cnmf_sam[metacells_order]
  
#  cnmf_l[[cnmf_mod]] <- cor (
#   rollapply (mMat_sam, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"),
#   rollapply (cnmf_sam, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"), 
#   method = 'spearman')  
#   }
# df[[sam]] = as.data.frame (do.call (cbind, cnmf_l))
# colnames (df[[sam]]) = cnmf_mods
# #df[[sam]]$sample = sam
# }

# # Convert list to 3D array: rows x cols x n_matrices
# array_3d <- array(unlist(df), dim = c(nrow(df[[1]]), ncol(df[[1]]), length(df)))

# # Apply median across the third dimension (i.e., across matrices)
# median_matrix <- apply(array_3d, c(1, 2), median)
# rownames (median_matrix) = rownames (df[[1]])
# colnames (median_matrix) = colnames (df[[1]])
# top_5 = unlist(lapply (cnmf_mods, function(x) head (rownames(median_matrix[order(-median_matrix[,x]),]),5)))




# atac_mat = cbind (cnmf_scatac, mod_2 = archp$mod_2, sample = archp$Sample)
# atac_mat_long = gather (as.data.frame(atac_mat), cnmf, score, 1:(ncol(atac_mat)-2))
# atac_mat_long$cnmf = factor (atac_mat_long$cnmf, levels = cnmf_mods)

# library (ggpointdensity)
# # remove outliers
# #atac_mat_longL = split (atac_mat_long, atac_mat_long$cnmf)
# #atac_mat_longL = lapply (atac_mat_longL, function(x) x[x$score > quantile(x$score,.005) & x$score < quantile(x$score,.995),])
# #atac_mat_long = do.call (rbind,atac_mat_longL)
# sp = ggplot (atac_mat_long, aes(x = score, y = mod_2, fill = sample)) +
#   geom_point(shape = 21,
#     alpha = 0.6, 
#     size = 2
#   ) + 
#   facet_wrap (~ interaction(cnmf, sample, sep = " | "), scales = 'free', ncol= length(cnmf_mods)) #+
#    #geom_pointdensity (alpha=1, size=.1)# +
# #   scale_color_viridis (option='F') +
# # # Scatterplot points with transparency and size
# #   geom_smooth(
# #     method = "lm", 
# #     color = "white", 
# #     fill = "white", 
# #     se = FALSE, 
# #     linetype = "dashed", 
# #     size = .4
# #   ) + # Regression line with confidence interval
# #   theme_void() + # Clean theme
# #   labs(
# #     title = "mod_2 vs cnmfs",
# #     #subtitle = "Scatterplot with Regression Line and Correlation Coefficient",
# #     x = "cnmf",
# #     y = "mod_2")
# #   # ) +
#   # theme(
#   #   plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
#   #   plot.subtitle = element_text(size = 14, hjust = 0.5),
#   #   axis.title = element_text(size = 14),
#   #   axis.text = element_text(size = 12),
#   #   panel.grid.major = element_line(color = "grey90"),
#   #   panel.grid.minor = element_blank()
#   # ) 

# png (file.path ('Plots','cnmf_km_sample_scatterplots.png'),width=3000, height=5500, res=300)
# sp
# dev.off()

# sp = sp + stat_cor (
#     aes(label = paste(..rr.label.., ..p.label.., sep = " | ")), 
#     method = "spearman", 
#     #label.x = min(df$fetal) + 0.1 * diff(range(df$fetal)), 
#     #label.y = max(df$fetal_gs) - 0.1 * diff(range(df$fetal_gs)),
#     color = "grey11",
#     size = 1
#   )   
# pdf (file.path ('Plots','cnmf_AP1_scatterplots3.pdf'),width=7, height=1.4)
# sp
# dev.off()

# # Try with PCA on mMat ####
library (ggridges)
library (ggplot2)
library (viridis)
library (tidyr)
#library(hrbrthemes)


# Add metacolumns of average TF modules activity ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
#mMat = scale(as.matrix(mMat))#[selected_TF,])

# # Filter by RNA expression ####
metaGroupName = 'celltype_lv2'
min_exp = .1
active_TFs = exp_genes (srt, rownames(mMat), min_exp = min_exp, metaGroupName)
mMat = t(scale (mMat[active_TFs, ]))
library(uwot)
library(ggplot2)
library(patchwork)

#-------------------------
# 0. RUN PCA
#-------------------------
set.seed(123)  # for reproducibility
p <- prcomp(mMat, center = TRUE, scale. = TRUE)
plot_df <- data.frame(p$x, celltype = archp$celltype_lv3)
df_sub <- plot_df[plot_df$celltype %in% 
                    c("Mono_CD14","TAM_interstitial","TAM_TREM2",
                      "TAM_MARCO","TAM_CXCLs"), ]

plot_df <- data.frame(
  PC1 = p$x[,1],
  PC2 = p$x[,2],
  celltype = archp$celltype_lv3
)

df_sub <- plot_df[plot_df$celltype %in%
                    c("Mono_CD14","TAM_interstitial","TAM_TREM2",
                      "TAM_MARCO","TAM_CXCLs"), ]


#-------------------------
# 1. MAIN SCATTER
#-------------------------
p_scatter <- ggplot(df_sub, aes(PC1, PC2, color = celltype)) +
  geom_point(alpha = 0.4, size = 0.1) +
  geom_density_2d(size = 0.4) +
  scale_color_manual(values = palette_myeloid) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )


#-------------------------
# 2. TOP DENSITY (PC1)
#-------------------------
p_density_x <- ggplot(df_sub, aes(PC1, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  gtheme_no_rot


#-------------------------
# 3. RIGHT DENSITY (PC2)
#-------------------------
p_density_y <- ggplot(df_sub, aes(PC2, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  gtheme_no_rot


#-------------------------
# 4. SAVE PDF
#-------------------------
pdf(file.path("Plots","TF_activity_pca.pdf"),
    width = 5, height = 5)

p_density_x + plot_spacer() + p_scatter + p_density_y +
  plot_layout(
    ncol = 2, nrow = 2,
    widths  = c(4, 1),
    heights = c(1, 4),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

dev.off()



### Ridge plot of PC1 ####
library(ggridges)
df_sub$celltype = factor (df_sub$celltype, levels = rev(cell_subsets_order))
p_ridge <- ggplot(df_sub, aes(x = PC1, y = celltype, fill = celltype)) +
  geom_density_ridges(scale = 3, alpha = 0.7, color = NA) +
  scale_fill_manual (values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y  = element_text(size = 8),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  gtheme_no_rot

pdf(file.path("Plots","momac_PC1_ridge_plots.pdf"),
    width = 7, height = 2)
p_ridge
dev.off()

### Correlate PCs of TF activity with scRNA-seq derived PCs of TF regulons ####
regulons_PCs = readRDS ('../scrna/PCA_regulons.rds') 
momac_axis_rot = p$rotation[,1]
regulons_rotations = regulons_PCs$rotation
pc_var <- regulons_PCs$sdev^2
names(pc_var) = colnames(regulons_rotations)
#rownames(regulons_rotations) = regulons_rotations$X
#regulons_rotations = regulons_rotations[,-1]
regulons_rotations = as.data.frame(regulons_rotations)[names(momac_axis_rot), ]
cor_df = as.data.frame (cor (momac_axis_rot, regulons_rotations, use = 'pairwise.complete.obs', method='pearson'))


df <- data.frame(
  PC = names(cor_df),
  loading = as.numeric(cor_df),
  var = pc_var
) %>%
  mutate(abs_loading = abs(loading)) %>%
  arrange(desc(abs_loading)) %>%
  mutate(PC = factor(PC, levels = PC))


top_n <- 5   # number of PCs to label per metric

df_lab <- bind_rows(
  df %>% top_n(top_n, abs_loading),
  df %>% top_n(top_n, var)
) %>% distinct(PC, .keep_all = TRUE)

sp = ggplot(df, aes(x = var, y = abs_loading)) +
  geom_point(alpha = 0.7, color = "steelblue", size = 2) +
  geom_text(
    data = df_lab,
    aes(label = PC),
    vjust = -0.6,
    size = 3.2,
    color = "black"
  ) +
  gtheme_no_rot +
  labs(
    x = "PC Variance",
    y = "Absolute Loading",
    title = "PC Absolute Loading vs PC Variance")


pdf (file.path ('Plots','PC_correlation_scrna_barplot.pdf'),height=3,4)
sp
dev.off()

### Filter TFs by genescore correlation ####
seGroupMotif <- getGroupSE (ArchRProj = archp, useMatrix = "MotifMatrix", groupBy = "Clusters_H")
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

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
#corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.3 & corGSM_MM$padj < 0.05 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
positive_TF = sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
#corGSM_MM = corGSM_MM[corGSM_MM$GeneScoreMatrix_name %in% active_TFs,]
#corGSM_MM = data.frame(corGSM_MM)

pdf ()
pg = ggplot(corGSM_MM, aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )+
  # Add labels for TFRegulator == YES
  geom_text_repel(
    data = subset(corGSM_MM, TFRegulator == "YES"),
    aes(label = MotifMatrix_name),
    size = 3,
    max.overlaps = 20
  )
dev.off()

pdf (file.path ('Plots','positive_expression_TFs.pdf'), width = 4,height=4)
pg
dev.off()



### Make momac module taking all positive TF contributing to PC1 ####
momac_TF = p$rotation[,1][order(p$rotation[,1])]
momac_TF = momac_TF[momac_TF < 0]
momac_TF = momac_TF[names(momac_TF) %in% positive_TF]
momac_module = rowMeans (mMat[,names(momac_TF)])
archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% c('momac_module','momac')] 
archp$momac_module = momac_module
archp$momac = ifelse (momac_module > 0, 'momac','resident')


table (archp$momac, archp$inflamed)

pdf (file.path ('Plots','momac_module_umap.pdf'))
plotEmbedding (
    ArchRProj = archp,
    colorBy = "cellColData",
    name = 'momac_module', 
    #useSeqnames='z',
    pal = rev (palette_deviation),    
    embedding = "UMAP_H",
    imputeWeights = NULL
    )
dev.off()

# # Ridge plot of PCA comp / km cluster ####
# ccomp = data.frame (module = plot_df[,2], celltype = plot_df$celltype)
# #ccomp = ccomp[ccomp$cnmf_celltypes %in% c('cDCs'),]
# cell_subsets_order = c("TAM_interstitial","TAM_MARCO","cDCs","TAM_TREM2","TAM_CXCLs","Mono_CD16","Mono_CD14")
# ccomp$celltype = factor (ccomp$celltype, levels = cell_subsets_order)
# #ccomp$module = srt$mod_2
# #ccomp = ccomp[!is.na(ccomp$cnmf_celltypes),]
# rp <- ggplot(ccomp, aes(x = module, y = celltype, fill = ..x..)) +
#   geom_density_ridges_gradient(
#   scale = 3,
#   rel_min_height = 0.01,
#   linewidth = 0.4,
#   color='white',
#   alpha = 0.3
# ) +

#   scale_fill_gradientn (colors = viridis::rocket(100)) +  # Optional: nice color gradient
#   theme_ridges() +                      # Optional: clean ridge plot theme
#   theme(legend.position = "right")     # Adjust legend position
# #   theme_classic() + facet_wrap (~sample, ncol=5)
# pdf (file.path ('Plots','inflammation_module_ridge_plots.pdf'), width = 5,height=3)
# rp
# dev.off()


# Make heatmap of correlation of inflammation module vs cnmf modules ####
# TF_modules = split(names(km$cluster), km$cluster)
# srt = ModScoreCor (
#     seurat_obj = srt, 
#     geneset_list = TF_modules, 
#     cor_threshold = NULL, 
#     pos_threshold = NULL, # threshold for fetal_pval2
#     listName = 'TF_module', outdir = NULL)



# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat = scale(as.matrix (mMat[active_TFs,]))#[selected_TF,])
# mMat = mMat[names(km$cluster[km$cluster == 2]),]

# mMat = colMeans (mMat)
# all (colnames (mMat) == rownames(archp@cellColData))
# archp$AP1 = mMat
# rna_mod_cor = cor (srt@meta.data[,cnmf_mods],srt@meta.data[,'2'])
# atac_mod_cor = cor (scale(archp@cellColData[,cnmf_mods]),mMat, method='spearman')
# #cor.test (archp@cellColData[,names(shared_cnmf)],archp@cellColData[,'mod_2'])

# hm2 = Heatmap (atac_mod_cor, border=T, col = palette_deviation_fun(atac_mod_cor))
# hm1 = Heatmap (rna_mod_cor, cluster_rows=F, border=T, col = rev (palette_deviation_correlation))
# pdf (file.path ('Plots','atac_module_cor3.pdf'), width=3.5, height=3)
# hm2 + hm1
# dev.off()

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

### Plot correlation of regulon score of TFs found in km1 along with correlation of km2 average score from atac #####
# #auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
# auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
# rownames (auc_mtx) = auc_mtx[,1]
# auc_mtx = auc_mtx[,-1]
# colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))
# regulon_TFs_in_modules = list(
#   km1 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 1])],
#   km2 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 2])]
#   )
# auc_mtx = auc_mtx[, colnames(auc_mtx) %in% regulon_TFs_in_modules$km1]

# #colnames(srt@meta.data)[colnames(srt@meta.data) == 'Mono'] = 'Monocytes'
# cnmf_mods = cell_subsets_order
# rownames(auc_mtx) = gsub ('\\.','-',rownames(auc_mtx))
# all (rownames(auc_mtx) == colnames(srt))

# auc_mtx_cor = as.data.frame (cor (auc_mtx, srt@meta.data[,cnmf_mods]))
# auc_mtx_cor$TF = rownames(auc_mtx_cor)

# auc_mtx_avg_scaled_l = gather (auc_mtx_cor, celltype, score,1:(ncol(auc_mtx_cor)-1))
# auc_mtx_avg_scaled_l$celltype = factor (auc_mtx_avg_scaled_l$celltype, levels = c('Mono','SPP1','TREM2','IFN_CXCLs','IM'))


# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat = as.matrix(mMat)#[selected_TF,])
# atac_mod_2 = as.data.frame(t(mMat[regulon_TFs_in_modules$km2,]))

# atac_mod_2 = as.data.frame (cor (atac_mod_2, t(scale(t(as.data.frame(archp@cellColData[,cnmf_mods]))))))
# atac_mod_2$TF = rownames(atac_mod_2)
# atac_mod_2 = gather (atac_mod_2, celltype, score,1:(ncol(atac_mod_2)-1))
# atac_mod_2$celltype = factor (atac_mod_2$celltype, levels = cnmf_mods)

# ### Plot  TF activity
# atac_mod_2_summary <- atac_mod_2 %>%
#   group_by(celltype) %>%
#   dplyr::  summarize(
#     mean_score = mean(score, na.rm = TRUE),
#     sd_score = sd(score, na.rm = TRUE),
#     n = n(),
#     se = sd_score / sqrt(n)  # standard error
#   )
# atac_mod_2_summary$celltype = factor (atac_mod_2_summary$celltype, levels = cnmf_mods)

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
# rna_mod_2_summary$celltype = factor (rna_mod_2_summary$celltype, levels = cnmf_mods)


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


# pdf (file.path ('Plots','inflamed_module_atac_rna_lineplot3.pdf'),width=5,height=5)
# gp
# dev.off()








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




# # Show all TFs included in inflammation module ####
# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat = as.matrix(mMat)#[selected_TF,])
# atac_mod_2 = as.data.frame(t(mMat[regulon_TFs_in_modules$km2,]))

# atac_mod_2 = as.data.frame (cor (atac_mod_2, t(scale(t(as.data.frame(archp@cellColData[,cnmf_mods]))))))
# #atac_mod_2$TF = rownames(atac_mod_2)

# hm = Heatmap (
#     atac_mod_2[,cnmf_mods],
# #    right_annotation = ha2,
#     column_names_rot =45, 
#     row_names_gp = gpar(fontsize = 5),
#     column_names_gp = gpar(fontsize = 6),
#     col = rev(as.character(palette_deviation[-1])), 
#     cluster_rows=T,
#     cluster_columns = F,
#     border=T
# #rect_gp = gpar (col = "white", lwd = 1)
# )


# ### Check AP1 on the RNA side ####
# auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
# rownames (auc_mtx) = auc_mtx[,1]
# auc_mtx = auc_mtx[,-1]
# colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))
# auc_mtx = auc_mtx[, colnames(auc_mtx) %in% regulon_TFs_in_modules$km2]

# rownames(auc_mtx) = gsub ('\\.','-',rownames(auc_mtx))
# all (rownames(auc_mtx) == colnames(srt))

# auc_mtx_cor = as.data.frame (cor (auc_mtx, srt@meta.data[,cnmf_mods]))
# #auc_mtx_cor$TF = rownames(auc_mtx_cor)

# hm2 = Heatmap (
#     auc_mtx_cor[,cnmf_mods],
#     #right_annotation = ha2,
#     column_names_rot =45, 
#     row_names_gp = gpar(fontsize = 5),
#     column_names_gp = gpar(fontsize = 6),
#     col = rev(as.character(palette_expression_correlation[-1])), 
#     cluster_rows=T,
#     cluster_columns = F,
#     border=T
# #rect_gp = gpar (col = "white", lwd = 1)
# )


# pdf (file.path ('Plots','inflammation_module_atac_rna_TFs_cor_heatmap2.pdf'), width = 3.6,height=3)
# hm 
# dev.off()

# ### Check footprint of AP1-complex and NFKB1 ####
# metaGroupName='inflamed'
# archp <- addGroupCoverages (ArchRProj = archp, groupBy = metaGroupName)
# motifPositions <- getPositions (archp)

# motifs <- c('NFKB2','JUNB','FOS','JUND')

# markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

# seFoot <- getFootprints(
#   ArchRProj = archp, 
#   #positions = motifPositions_sample[markerMotifs], 
#   positions = motifPositions[markerMotifs], 
#   groupBy = metaGroupName
# )
  
# plotFootprints(
# seFoot = seFoot,
# ArchRProj = archp, 
# normMethod = "Subtract",
# plotName = "Footprints-Subtract-Bias_inflamed_",
# addDOC = FALSE, height=4, width=2,
# #pal = palette_tnk_cells,
# smoothWindow = 25)



# ### Check inflammation score across TAMs ####
# mod_df = data.frame (
#   #celltype = archp$celltype_lv2,
#   #sample = archp$Sample,
#   #Infl_module = archp$mod_2,
#   Infl_module = archp$mod_2)
# mod_df = aggregate(archp$mod_2, by=list(
#   celltype = archp$celltype_lv2,
#   sample = archp$Sample), FUN=mean)
# head (mod_df)
# df_order = mod_df %>% 
# group_by (celltype) %>% 
# summarize (avg_module = median(x)) %>% 
# arrange(avg_module)
# mod_df$celltype = factor (mod_df$celltype, levels = rev(df_order$celltype))

# bp = ggplot (mod_df, aes (x = celltype, y = x, fill=celltype)) +
# vlp + 
# bxpv + 
# scale_fill_manual (values = palette_myeloid) +
# #geom_point (position='identity', alpha=.3, color="grey44", size=1) +
# gtheme

# pdf (file.path ('Plots','celltype_infl_module_sample_boxplots.pdf'),2.3,width=4)
# bp
# dev.off()




# # Compare cNMF modules with inflammatory program in scATAC-seq and scRNA-seq ####
# shared_cnmf = readRDS (file.path('..','scrna','shared_cnmf_myeloid.rds'))
# shared_cnmf = lapply (shared_cnmf, function(x) x[x %in% getFeatures (archp)])
# #remove_modules = c('cnmf.3','cnmf.6','cnmf.7','cnmf.5') # remove monocyres cDC and CC modules. Consider re-inculding CC 

# pdf (file.path ('Plots','scrna_celltype_dimplot.pdf'))
# DimPlot (srt, group.by = 'celltype', reduction = 'umap')
# dev.off()

# srt = ModScoreCor (
#     seurat_obj = srt, 
#     geneset_list = shared_cnmf, 
#     cor_threshold = NULL, 
#     pos_threshold = NULL, # threshold for fetal_pval2
#     listName = 'shared_cnmf', outdir = NULL)




# # TF_modules = c(
# # 'JUNB
# # FOSL2
# # JUN
# # SMARCC1
# # FOSL1
# # JUND
# # FOS
# # JDP2
# # BACH1
# # FOSB')
# # TF_modules = strsplit(TF_modules, '\n')
# # names (TF_modules) = 'AP1'



# Differential Peaks in cells positive for km2 vs rest ####
# Find DAG ####
#force = FALSE
#archp$inflamed = ifelse (archp$momac > 0, 'inflamed','non_inflamed')
metaGroupName = 'momac'
force=T
if (!file.exists ('DAG_inflamed_pairwise.rds') | force)
  {
  DAG_list = getMarkerFeatures (
    ArchRProj = archp, 
    testMethod = "wilcoxon",
          useGroups = "momac",
          bgdGroups = "resident",
    k=100,
    binarize = FALSE,
    useMatrix = "GeneScoreMatrix",
    groupBy = metaGroupName
  #  useSeqnames="z"
  )
  saveRDS (DAG_list, 'DAG_inflamed_pairwise.rds')
  } else {
  DAG_list = readRDS ('DAG_inflamed_pairwise.rds')
  }
DAG_res = do.call (cbind, (assays(DAG_list)))
colnames (DAG_res) = names(assays(DAG_list))
#DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames (DAG_res) = rowData(DAG_list)$name

pdf(file.path('Plots','inflamed_MA_plot.pdf'), width=5,height=5)
gma <- markerPlot (seMarker = DAG_list, name = 'momac', cutOff = "FDR <= 0.01", plotAs = "MA")
gma
dev.off()
  
# Take only significant regions ####
DAG_res_sig = DAG_res[DAG_res$FDR < .01 & DAG_res$Log2FC > 0, ]
saveRDS (DAG_res_sig, 'momac_genescore.rds')

# ### Perform enrichment on DAP ####
# archp = addBgdPeaks (archp, force= T)
# archp = addMotifAnnotations (ArchRProj = archp, 
#       motifSet = "cisbp", 
#       #motifSet = 'JASPAR2020',
#       #name = "JASPAR2020_Motif",
#       force=T)
# enrichMotifs <- peakAnnoEnrichment(
#     seMarker = DAP_list,
#     ArchRProj = archp,
#     peakAnnotation = "Motif",
#     cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
#   )
# heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
# pdf (file.path ('Plots','enrich_inflamation_heatmap.pdf'))
# heatmapEM
# dev.off()

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
#pathways = pathways[grep('inflammatory', pathways$term,ignore.case=T),]
#pathways = read.csv (csv.file)
pathways = split(pathways$gene, pathways$term)
#pathways = gmtPathways (gmt.file)
DAG_list = readRDS ('DAG_inflamed_pairwise.rds')
DAG_res = do.call (cbind, (assays(DAG_list)))
colnames (DAG_res) = names(assays(DAG_list))
#DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames(DAG_res) = rowData (DAG_list)$name

# peak_genes = unname(ps[queryHits (findOverlaps(ps, GRanges (rownames(DAP_res))))]$nearestGene)
# names (peak_genes) = as.character(ps)[queryHits(findOverlaps (ps,GRanges (rownames(DAP_res))))]
# peak_genes = peak_genes[rownames(DAP_res)]
# #peak_genes = setNames (-log10(DAP_res$Pval) * sign (DAP_res$Log2FC), peak_genes)
# peak_genes = setNames (DAP_res$Log2FC, peak_genes)
# #peak_genes = peak_genes[!duplicated(names(peak_genes))]
# #names (peak_genes) = 
# peak_genes = peak_genes[!duplicated(names(peak_genes))]
# peak_genes = peak_genes[!is.na(names(peak_genes))]
ranked_genes = setNames(DAG_res$Log2FC,rownames(DAG_res))
library (BiocParallel)
BiocParallel::register(BiocParallel::SerialParam())

### fgsea throws a BiocParallel error when I load all packages including clusterProfiler...try avoiding loading packages except ArchR and fgsea
#peak_genes2 = setNames(order(peak_genes), names(peak_genes))
fgseaRes = fgseaMultilevel (pathways, 
          ranked_genes,#, 
          minSize=15, 
          scoreType='std',
          maxSize=5000,
          nproc=1,
          nPermSimple=10000,
          BPPARAM = NULL
          )
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                              pathways, ranked_genes)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
                                 order(-NES), pathway]
fgseaRes_df = as.data.frame (fgseaRes[pathway %in% mainPathways,])
#head (fgseaRes_df[order (fgseaRes_df$pval), ],100)
pvalAdjTrheshold = 0.05
top_pathways=10
#fgseaRes_df = as.data.frame (fgseaRes)
#fgseaRes_df$padj = fgseaRes_df$pval
#fgseaRes$padj = fgseaRes$pval
fgseaResAll_dp = dotGSEA (
  list(fgseaRes_df), 
  padj_threshold = pvalAdjTrheshold, 
  type = 'fgsea',
  top_pathways = top_pathways,
  cluster_rows=F,
  cluster_cols=F)

pdf (file.path ('Plots','fgsea_dotplot2.pdf'), width=8, height=4)
fgseaResAll_dp
dev.off()


pdf (file.path ('Plots','GO_Inflammatory_enrichment_plot.pdf'), width=5, height=3)
plotEnrichment(pathways[["GO_INFLAMMATORY_RESPONSE"]],
               ranked_genes) + labs(title="GO_INFLAMMATORY_RESPONSE")
dev.off()






## Run peak2genes results with hubs links ####
run_p2g = F
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
force=F
if (!file.exists(file.path (hubs_dir, paste0('hubs_cells_mat.rds'))) | force)
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
metaGroupName = 'celltype_lv3'
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
res_df_top$group = factor (res_df_top$group, levels = cell_subsets_order)
res_df_top = res_df_top[order(res_df_top$group), ]


# Generate matrix of fragment counts of hubs x metagroup ####
metaGroupName = 'celltype_lv3'
force = F
if (!file.exists(file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds'))) | force)
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
  frags_in_sample = sapply (unique(archp@cellColData[,metaGroupName]), function(x) sum (archp$nFrags[as.character(archp@cellColData[,metaGroupName]) == x]))
  hubsSample_mat = t(t(hubsSample_mat) * (10^6 / frags_in_sample)) # scale
  saveRDS (hubsSample_mat, file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))
  } else {
  hubsSample_mat = readRDS (file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))  
  }
hubsSample_mat = as.data.frame (hubsSample_mat)

#as.character(archp@cellColData[,metaGroupName])[order(metagroup, na.last=T)]
#ha = HeatmapAnnotation (celltype = metagroup[order (metagroup, na.last=T)])
#hubs_mat
hubsSample_mat = hubsSample_mat[,levels(res_df_top$group)]
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

pdf (file.path ('Plots','top_DAH_celltype_lv2_heatmap2.pdf'), width=5, height=5)
hm
dev.off()


#TF = 'SNAI1'
# TF = sapply (unique(res_df_top$gene), function(x) unlist(strsplit(x, '-'))[1])
# metaGroupName = 'celltype_lv2'

# Compare momacs vs resident ####
library (presto)
#archp$inflamed = ifelse (archp$mod_2 > 0, 'inflamed','non_inflamed')
metaGroupName = 'momac'
all (colnames(hubsCell_mat) == rownames(archp@cellColData))
metagroup = as.character (archp@cellColData[,metaGroupName])
res = wilcoxauc (log2(hubsCell_mat+1), metagroup)
res = res[res$group == 'momac',]
res = res[res$logFC > 0,]
res_df = res[order(res$padj), ]
#res_df = res_df[res_df$logFC > 0.2,]
res_df = head (res_df, 10)
res_df$gene = hubs_obj$hubsCollapsed$gene[match(res_df$feature, hubs_obj$hubs_id)]

hubsSample_mat = hubsSample_mat[,cell_subsets_order]
dah_hm = t(scale(t(log2(hubsSample_mat[res_df$feature,]))))

hm = Heatmap (dah_hm,
  cluster_rows=F,
#  top_annotation = ha,
  cluster_columns=F,
  row_labels = res_df$gene,
  col = rev(palette_hubs_accessibility),
  column_names_gp= gpar (fontsize=11),
  row_names_gp= gpar (fontsize=4),
  column_names_rot=45,
  border=T)

pdf (file.path ('Plots','DAH_inflamed_celltype_lv2_heatmap.pdf'), width=5, height=10)
hm
dev.off()



top_hubs = 20
res_df_top = res_df %>% group_by (group) %>%
  slice_head(n = top_hubs)
res_df_top$gene = hubs_obj$hubsCollapsed$gene[match(res_df_top$feature, hubs_obj$hubs_id)]
hub = res_df_top$feature[10]
TF = 'REL'
hub = res_df[grep (TF, res_df$gene),][1,]$feature
hub = 'HUB417'
#hub = hubs_obj$hubs_id[grep ('NFKB1', hubs_obj$hubsCollapsed$gene)]
#sample_levels = c('Monocytes','cDCs','SPP1','TREM2','C1Q','IFN','IM')

metaGroupName = 'celltype_lv3'
metaGroupName = 'inflamed'
ccomp = as.data.frame (archp@cellColData)
median_order = sort (unlist(lapply (split (ccomp$mod_2, ccomp$celltype_lv2), function(x) median(x))))
cell_subsets_order2 = rev(names (median_order))
palette_inflamed = c(inflamed = 'darkred',non_inflamed='grey22')
pdf()
#archp$fetal_sample = paste0(archp$Sample, archp$fetal_group)
#metaGroupName = 'fetal_group'
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp,#[!archp$Sample3 %in% c('P11_HOX')], 
    sample_levels = cell_subsets_order, 
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
    region = ext_range (GRanges (hubs_obj$hubsCollapsed[match(hub[1], hubs_obj$hubs_id)]),50000,50000),
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    pal = palette_myeloid,
    #pal = palette_inflamed,
    #loops = getCoAccessibility (archp, corCutOff = 0.3,
    #  returnLoops = TRUE),
    useGroups= NULL
)
dev.off()
plotPDF (meso_markers, ArchRProj = archp,height=3.5, width=4, name =paste0(metaGroupName,'_MPM_markers_inflamed_coveragePlots.pdf'),addDOC=F)


# metaGroupName = 'shared_cnmf2_r_max'
# hub_gr = ext_range (GRanges (hubs_obj$hubsCollapsed[match(hub, hubs_obj$hubs_id)]),10000,10000)
# genes_in_region = unique(getPeakSet (archp)[subjectHits (findOverlaps (hub_gr, getPeakSet (archp)))]$nearestGene)
# genes_in_region = c(genes_in_region)
# genes_in_region = 'CD44'
# top_dah = data.frame (
#   gene = colMeans (srt@assays$RNA@data[rownames(srt) %in% genes_in_region,,drop=F]),
#   group = srt@meta.data[,metaGroupName])
# top_dah$group = factor (top_dah$group, levels = rev(sample_levels))
# top_dah = na.omit(top_dah)
# bp = ggplot (top_dah, aes (x = gene, y = group, fill = group)) + 
# vlp + 
# bxpv + 
# scale_fill_manual (values = palette_myeloid) +
# #geom_point (position='identity', alpha=.3, color="grey44", size=1) +
# gtheme_no_rot

# pdf (file.path ('Plots', paste0('scrna_region_boxplots.pdf')), height=4, width=4)
# bp
# dev.off()


# # Check footprint across celltypes ####
# metaGroupName='celltype_lv2'
# archp <- addGroupCoverages (ArchRProj = archp, groupBy = metaGroupName)
# motifPositions <- getPositions (archp)

# motifs <- c('NFKB1','JUNB','FOS','JUND','SPI1','SPIB','MITF','RUNX1','CEBPA','SRF','NFAT5','STAT2','MEF2C','PRDM1','IRF3','IRF8','IRF1','IRF2','IRF9')
# markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

# seFoot <- getFootprints(
#   ArchRProj = archp, 
#   #positions = motifPositions_sample[markerMotifs], 
#   positions = motifPositions[markerMotifs], 
#   groupBy = metaGroupName
# )

# plotFootprints(
# seFoot = seFoot,
# ArchRProj = archp, 
# normMethod = "Subtract",
# plotName = "Footprints-Subtract-Bias_",
# addDOC = FALSE, height=7.5, width=5,
# pal = palette_myeloid,
# smoothWindow = 25)
  


# # Check TF deviations
# TF = 'E2F3'
# getFeatures (archp, 'MotifMatrix')[grep (TF, getFeatures (archp, 'MotifMatrix'))]
# TF1 = c('z:JUN_143','z:FOS_137','z:NFKB1_719','z:FOXM1_352','z:TFDP1_310','z:E2F3_313')

# pdf ()
# TF_p = plotEmbedding (
#     ArchRProj = archp,
#     colorBy = "MotifMatrix",
#     name = TF1, 
#     useSeqnames='z',
#     pal = rev (palette_deviation),    
#     embedding = "UMAP_H",
#     imputeWeights = NULL
#     )
# dev.off()

# pdf(file.path('Plots','avg_AP1_deviation_fplot.pdf'),15,15)
# wrap_plots(TF_p)
# wrap_plots(umap_p1)
# wrap_plots (umap_p0,umap_p2)#,umap_p3)
# dev.off()

# pdf (file.path ('Plots','FRIP_umap.pdf'))
# plotEmbedding (
#     ArchRProj = archp,
#     colorBy = "cellColData",
#     name = 'FRIP', 
#     useSeqnames='z',
#     pal = rev (palette_deviation),    
#     embedding = "UMAP_H",
#     imputeWeights = NULL
#     )
# dev.off()




# Export bigiwg files ####
metaGroupName = 'inflamed'
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

# is_sequential = function(x) {
#   length(x) > 1 && all(diff(x) == 1)
# }
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



# Import chromBPnet finemo motifs ####
#archp_P1 = archp[archp$Clusters %in% c('C2') & archp$Sample == 'P1']
library ('universalmotif')

ps = getPeakSet (archp)

chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/chromBPnet'

chrombpnet_counts = list()
celltypes = c('inflamed','non_inflamed')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_counts[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
  gr = makeGRangesFromDataFrame (chrombpnet_counts[[celltype]], keep.extra.columns=T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  chrombpnet_counts[[celltype]]$peak_type = ps$peakType[findOverlaps(gr, ps, select='first')]
  chrombpnet_counts[[celltype]] = chrombpnet_counts[[celltype]][chrombpnet_counts[[celltype]]$V4 != 'NaN_NaN_NaN',] # Some seqlets have NA motif match and qvalues...removing those
  nonsig_motifs = chrombpnet_counts[[celltype]]$V7 > 0.05
  chrombpnet_counts[[celltype]]$V4[nonsig_motifs] = chrombpnet_counts[[celltype]]$V6[nonsig_motifs]
  }


chrombpnet_profile = list()
#celltypes = c('Mesothelium','Malignant')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  gr = makeGRangesFromDataFrame (chrombpnet_profile[[celltype]], keep.extra.columns=T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  chrombpnet_profile[[celltype]]$peak_type = ps$peakType[findOverlaps(gr, ps, select='first')]
  chrombpnet_profile[[celltype]] = chrombpnet_profile[[celltype]][chrombpnet_profile[[celltype]]$V4 != 'NaN_NaN_NaN',] # Some seqlets have NA motif match and qvalues...removing those
  nonsig_motifs = chrombpnet_profile[[celltype]]$V7 > 0.05
  chrombpnet_profile[[celltype]]$V4[nonsig_motifs] = chrombpnet_profile[[celltype]]$V6[nonsig_motifs]
  }

#chrombpnet_profile = lapply (chrombpnet_profile, function(x) x[x$V5 != 'NaN_NaN_NaN',])


top_n <- 5
n <- length(chrombpnet_counts)

bp_list <- lapply(seq_len(n), function(i) {
  tbl <- table(chrombpnet_counts [[i]]$V4)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_counts[[i]]$V5[chrombpnet_counts [[i]]$V4 == tf][1]
  })
  
  data.frame(
    Freq = proportions(top_tbl),
    TF   = tf_names,
    direction = directions,
    type = rep(celltypes[[i]], length(top_tbl))
  )
})

bp_df <- do.call(rbind, bp_list)

# Make neg values negative
bp_df <- bp_df %>%
  mutate(Freq = ifelse(direction == "neg", -Freq.Freq, Freq.Freq))

# Create custom ordering per type
bp_df <- bp_df %>%
  group_by(type, direction) %>%
  mutate(
    TF_order = ifelse(direction == "pos",
                      rank(-Freq, ties.method = "first"),  # descending
                      rank(Freq, ties.method = "first"))   # ascending for neg (opposite)
  ) %>%
  ungroup()

# Build a combined factor: ensures pos stack from bottom up, neg from top down
bp_df <- bp_df %>%
  arrange(type, direction, TF_order)

bp_df$TF_id <- paste(bp_df$TF, bp_df$type, sep = "_")
bp_df$TF_id <- factor(bp_df$TF_id, levels = unique(bp_df$TF_id))
bp_df$type = factor (bp_df$type, levels = celltypes)
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF)) ) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

  

pdf (file.path ('Plots', 'TF_abundance_counts_barplot.pdf'),6,width=6.5)
bp
dev.off()


n <- length(chrombpnet_profile)

bp_list <- lapply(seq_len(n), function(i) {
  tbl <- table(chrombpnet_profile[[i]]$V4)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_profile[[i]]$V5[chrombpnet_profile[[i]]$V4 == tf][1]
  })
  
  data.frame(
    Freq = proportions(top_tbl),
    TF   = tf_names,
    direction = directions,
    type = rep(celltypes[[i]], length(top_tbl))
  )
})

bp_df <- do.call(rbind, bp_list)

# Make neg values negative
bp_df <- bp_df %>%
  mutate(Freq = ifelse(direction == "neg", -Freq.Freq, Freq.Freq))

# Create custom ordering per type
bp_df <- bp_df %>%
  group_by(type, direction) %>%
  mutate(
    TF_order = ifelse(direction == "pos",
                      rank(-Freq, ties.method = "first"),  # descending
                      rank(Freq, ties.method = "first"))   # ascending for neg (opposite)
  ) %>%
  ungroup()

# Build a combined factor: ensures pos stack from bottom up, neg from top down
bp_df <- bp_df %>%
  arrange(type, direction, TF_order)

bp_df$TF_id <- paste(bp_df$TF, bp_df$type, sep = "_")
bp_df$TF_id <- factor(bp_df$TF_id, levels = unique(bp_df$TF_id))
bp_df$type = factor (bp_df$type, levels = celltypes)
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF)) ) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

  

pdf (file.path ('Plots', 'TF_abundance_profile_barplot.pdf'),6,width=6.5)
bp
dev.off()




### Compare TF regulon expression with inflammation regulatory module to identify which TF might be driving ####
auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

km = readRDS (file.path ('..','scatac_ArchR','TF_activity_modules.rds'))
#genes_highlight = readRDS (file.path ('..','scatac_ArchR','inflammation_TFs.rds'))
regulon_TFs_in_modules = list(
  km1 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 1])],
  km2 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 2])]
  )
auc_mtx = auc_mtx[, colnames(auc_mtx) %in% regulon_TFs_in_modules$km2]

srt$mod_2 = unname(rowMeans (auc_mtx))

# # Try with ridge plots ####
library (ggridges)
library (ggplot2)
library (viridis)
library (tidyr)
#library(hrbrthemes)

# Plot
ccomp = as.data.frame (srt@meta.data)
#ccomp = ccomp[ccomp$cnmf_celltypes %in% c('cDCs'),]
cell_subsets_order = c("TAM_interstitial","TAM_MARCO","cDCs","TAM_TREM2","TAM_CXCLs","Mono_CD16","Mono_CD14")
ccomp$cnmf_celltypes = factor (ccomp$shared_cnmf_r_max, levels = cell_subsets_order)
ccomp$module = srt$mod_2
#ccomp = ccomp[!is.na(ccomp$cnmf_celltypes),]
rp <- ggplot(ccomp, aes(x = module, y = cnmf_celltypes, fill = ..x..)) +
  geom_density_ridges_gradient(
  scale = 3,
  rel_min_height = 0.01,
  linewidth = 0.4,
  color='white',
  alpha = 0.3
) +

  scale_fill_gradientn (colors = palette_expression) +  # Optional: nice color gradient
  theme_ridges() +                      # Optional: clean ridge plot theme
  theme(legend.position = "right")     # Adjust legend position
#   theme_classic() + facet_wrap (~sample, ncol=5)
pdf (file.path ('Plots','cnmf_inflammation_module_ridge_plots.pdf'), width = 5,height=3)
rp
dev.off()


# #### Show regulon of TF in inflammation module ####
# samples_to_use = c('P1','P10','P11','P12','P13','P14','P3','P4','P5','P8')
# all (rownames (auc_mtx) == colnames (srt)) srt$celltype_lv3[match(rownames (auc_mtx), colnames (srt))]
# auc_mtx_agg = aggregate (auc_mtx, by = list(sample= srt$sampleID[match(rownames (auc_mtx), colnames (srt))], 
#   celltype = srt$celltype_lv3[match(rownames (auc_mtx), colnames (srt))]), mean)


# ### bind module score of inflammatory module ####
# infl_mod = 'mod_2'
# metaGroup= 'celltype_lv3'
# ccomp = as.data.frame (archp@cellColData[,infl_mod])
# colnames (ccomp) = infl_mod
# ccomp_mg = aggregate (ccomp, by = list(sample= archp$Sample, celltype = archp$celltype_lv3), mean)
# ccomp_mg = ccomp_mg[ccomp_mg$sample %in% samples_to_use,]


# # 1. Merge mod_2 with TF AUCs
# df = auc_mtx_agg %>%
#   left_join(ccomp_mg, by = c("sample", "celltype"))
# tf_cols <- setdiff(names(auc_mtx_agg), c("sample", "celltype"))
# cor_per_sample <- df %>%
#   group_by(sample) %>%
#   group_modify(~ {
#     tf_mat <- .x %>% select(all_of(tf_cols))
#     mod_vec <- .x$mod_2

#     tibble(
#       gene = tf_cols,
#       correlation = map_dbl(tf_mat, ~ cor(.x, mod_vec, use = "pairwise.complete.obs"))
#     )
#   }) %>%
#   ungroup()

# # Do the same but with RNA expression of TFs ####
# ps = log2(as.data.frame (AverageExpression (srt, 
#         features = cor_per_sample$gene,
#         group.by = c('sampleID','celltype_lv3'))[[1]]) +1)
# ps = as.data.frame (t(ps))
# ps2 = data.frame (sample = sapply (rownames(ps), function(x) unlist(strsplit(x,'_'))[1]), 
#   celltype = sapply (rownames(ps), function(x) unlist(strsplit(x,'_'))[2]))
# ps2$celltype = gsub ('-','_',ps2$celltype)
# ps = cbind (ps2, ps)
# ps = ps[, colnames(auc_mtx_agg)]

# df2 = ps %>%
#   left_join(ccomp_mg, by = c("sample", "celltype"))
# cor_per_sample2 <- df2 %>%
#   group_by(sample) %>%
#   group_modify(~ {
#     tf_mat <- .x %>% select(all_of(tf_cols))
#     mod_vec <- .x$mod_2

#     tibble(
#       gene = tf_cols,
#       correlation = map_dbl(tf_mat, ~ cor(.x, mod_vec, use = "pairwise.complete.obs"))
#     )
#   }) %>%
#   ungroup()


# # ---- 1. Select top 10 TFs with highest mean correlation ----
# top10_TFs <- cor_per_sample %>%
#   group_by(gene) %>%
#   summarise(mean_cor = median(correlation, na.rm = TRUE)) %>%
#   arrange(desc(mean_cor)) %>%
#   slice_head(n = 10) %>%
#   pull(gene)

# # ---- 2. Subset for the top 10 ----
# cor_top10 <- cor_per_sample %>%
#   filter(gene %in% top10_TFs)
# cor_top10$type = 'regulon'  

# cor_top10_2 <- cor_per_sample2 %>%
#   filter(gene %in% top10_TFs)
# cor_top10_2$type = 'TF_expression'  

# cor_top10_3 = rbind (cor_top10, cor_top10_2)

# # 1. Compute median correlation per gene for type 'regulon'
# gene_order <- cor_top10_3 %>%
#   filter(type == "regulon") %>%
#   group_by(gene) %>%
#   summarise(median_cor = median(correlation, na.rm = TRUE)) %>%
#   arrange(desc(median_cor)) %>%
#   pull(gene)

# # 2. Convert gene to factor with levels ordered by 'regulon' median
# cor_top10_3 <- cor_top10_3 %>%
#   mutate(gene = factor(gene, levels = gene_order))

# # 3. Plot with fill by type
# gp = ggplot(cor_top10_3[cor_top10_3$type == 'regulon',], aes(x = gene, y = correlation, fill = type)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(aes(color = type),  alpha = 0.7, size = 0.4) +
#   gtheme
# pdf (file.path ('Plots','corr_regulons_to_inflammatory_mod_boxplots.pdf'), height=3, width=5)
# gp
# dev.off()


#### Correlate gneescore with TF inflammatory module ####
#  selected_TF = names(km$cluster[km$cluster == 2])
#  sams = unique (archp$Sample)
#  # Get deviations ####
#   if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
#   mMat = scale (assays (mSE)[[1]])
#   rownames (mMat) = rowData (mSE)$name
#   mMat = mMat[selected_TF,]

#   # Get genescore ####
#   if (!exists('gsSE')) gsSE = fetch_mat (archp, 'GeneScore')
#   gsMat = assays (gsSE)[[1]]
#   rownames (gsMat) = rowData (gsSE)$name
#   gsMat = scale (gsMat[rownames (gsMat) %in% selected_TF,])
  
# #  mMat = scale(t(mMat))
#   mMat = lapply (sams, function(x) mMat[,archp$Sample == x])
#   names (mMat) = sams

#   gsMat = lapply (sams, function(x) gsMat[,archp$Sample == x])
#   names (gsMat) = sams

#   # Average mats along sarc module score ####
#   library(zoo)

#   bin_width <- 40   # Number of observations per bin
#   overlap <- 40    
#   mMat_ordered = lapply (sams, function(sam) mMat[[sam]][order(cnmf_mat[[sam]][,sarc_module]),])
#   names(mMat_ordered) = sams
#   mMat_ordered_avg = lapply (sams, function (sam) as.data.frame (lapply (as.data.frame (mMat_ordered[[sam]]), function(x) {
#     zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
#   })))
#   names (mMat_ordered_avg) = sams

#   gsMat_ordered = lapply (sams, function(sam) gsMat[[sam]][order(cnmf_mat[[sam]][,sarc_module]),])
#   names(gsMat_ordered) = sams
#   gsMat_ordered_avg = lapply (sams, function (sam) 
#     {
#     tmp_df = as.data.frame (lapply (as.data.frame (gsMat_ordered[[sam]]), function(x) {
#     zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
#     }))
#     tmp_df[is.na(tmp_df)] = 0 # HIC2 returns NaN for one or some samples
#     tmp_df
#     })

#   names (gsMat_ordered_avg) = sams

#   cnmfMat_ordered = lapply (sams, function(sam) cnmf_mat[[sam]][order(cnmf_mat[[sam]][,sarc_module]),])
#   names(cnmfMat_ordered) = sams
#   cnmfMat_ordered_avg = lapply (sams, function (sam) as.data.frame (lapply (as.data.frame (cnmfMat_ordered[[sam]]), function(x) {
#     zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
#   })))
#   names (cnmfMat_ordered_avg) = sams

#   m_cor = lapply (sams, function(sam) 
#     {
#     cor_tmp = as.data.frame (cor (mMat_ordered_avg[[sam]], cnmfMat_ordered_avg[[sam]][,sarc_module, drop=F], method = 'spearman'))
#     colnames(cor_tmp)[1] = 'score'
#     cor_tmp$sample = sam
#     cor_tmp$TF = rownames(cor_tmp)
#     cor_tmp
#     })
#   m_cor_df = do.call (rbind, m_cor)
#   m_cor_df$type = 'activity'
#   m_cor_levels = m_cor_df %>% group_by (TF) %>% summarise (median_value = median(score)) %>% arrange(-median_value)

#   gs_cor = lapply (sams, function(sam) 
#     {
#     cor_tmp = as.data.frame(cor (gsMat_ordered_avg[[sam]], cnmfMat_ordered_avg[[sam]][,sarc_module, drop=F], method = 'spearman'))
#     colnames(cor_tmp)[1] = 'score'
#     cor_tmp$sample = sam
#     cor_tmp$TF = rownames(cor_tmp)
#     cor_tmp
#     })
#   gs_cor_df = do.call (rbind, gs_cor)
#   gs_cor_df$type = 'genescore'

#   # Combine and plot ####
#   combined_df = rbind (m_cor_df, gs_cor_df)
  
#   # ranked_data <- combined_df %>%
#   # group_by(TF, type) %>%
#   # summarise(median_score = median(score), .groups = "drop") %>%
#   # group_by(TF) %>%
#   # summarise(mean_of_medians = mean(median_score), .groups = "drop") %>%
#   # arrange(desc(mean_of_medians)) %>%
#   # mutate(rank = row_number())

#   #ps = ps[, colnames(ps) %in% sample_names_rna]

#   #lapply (tc_cor, function(x) {x = x['cNMF19',]; head(x[order(-x)],10)})

# top_sarc_TF = head(m_cor_levels$TF,20)
# saveRDS (top_sarc_TF, 'top_sarc_TF.rds')
# combined_df$TF = factor (combined_df$TF, levels = top_sarc_TF)
# combined_df = combined_df[combined_df$TF %in% top_sarc_TF, ]
# palette_expression_disc = paletteer::paletteer_c("grDevices::Purple-Blue", length(top_sarc_TF))

# bp = ggplot (combined_df, aes(x = TF, y = score, fill = TF, color = type)) + 
#   geom_boxplot(
#     aes(group = interaction(TF, type)),  # Split by both TF and type
#     position = position_dodge(0.8),
#     #color = 'grey20',
#     linewidth = 0.4,
#     width = 0.7,
#     outlier.alpha = 0.2,
#     outlier.shape = NA,
#     size = 0.5,
#     alpha = 0.6
#   ) +
#   geom_point(
#     aes(group = interaction(TF, type), color = type),  # Align points with boxplots
#     position = position_dodge(0.8),  # Ensure same dodge width
#     alpha = 0.5,
#     size = 1
#   ) +
#   gtheme +
#   scale_fill_manual(values = palette_expression_disc) +
#   scale_color_manual(values = c(activity = '#AE123AFF', genescore = '#001260FF')) +
#   geom_hline(yintercept = 0, color = 'red', linetype = 'dashed')

#   pdf (paste0 ('Plots/sarcomatoid_score_TF_activity_boxplots2.pdf'), width = 8,height=3)
#   bp
#   dev.off()



### Filter TFs by genescore correlation ####
seGroupMotif <- getGroupSE(ArchRProj = archp, useMatrix = "MotifMatrix", groupBy = "Clusters_H")
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

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
#corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.3 & corGSM_MM$padj < 0.05 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
active_TF = sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
#corGSM_MM = corGSM_MM[corGSM_MM$GeneScoreMatrix_name %in% names(km$cluster[km$cluster == 2]),]
#corGSM_MM = data.frame(corGSM_MM)

pdf ()
p = ggplot(corGSM_MM, aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )+
  # Add labels for TFRegulator == YES
  geom_text_repel(
    data = subset(corGSM_MM, TFRegulator == "YES"),
    aes(label = MotifMatrix_name),
    size = 3,
    max.overlaps = 20
  )
dev.off()

pdf (file.path ('Plots','positive_expression_TFs.pdf'), width = 4,height=4)
p
dev.off()


# # # # Compute co-occurrence of TFs ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
#mMat = scale(as.matrix(mMat))#[selected_TF,])

# # Filter by RNA expression ####
metaGroupName = 'celltype_lv2'
min_exp = .1
active_TFs = exp_genes (srt, rownames(mMat), min_exp = min_exp, metaGroupName)

metaGroupName = 'momac'
celltypes = unique (archp@cellColData[,metaGroupName])
motifMat = getPositions (archp)
matches = getMatches (archp)
colnames(matches) = gsub ('_.*','',colnames(matches) )
colnames(matches) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames(matches))
matches = matches[,active_TFs]

# # Find DAP ####
#force = FALSE
force=F
if (!file.exists (paste0('DAP_',metaGroupName,'.rds')) | force)
  {
  DAP_list = getMarkerFeatures (
    ArchRProj = archp, 
    testMethod = "wilcoxon",
          useGroups = 'momac',
          bgdGroups = 'resident',
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

cooc_l = list()
ov_size_l=list()
celltype = 'up'
ps = getPeakSet (archp)
  matches_ct = matches[queryHits(findOverlaps (matches, ps))]
  matchesMat = assay (matches_ct)
  # colnames (matchesMat) = gsub ('_.*','',colnames (matchesMat))
  # colnames (matchesMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames (matchesMat))
  # matchesMat = matchesMat[,infl_TF]
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
        ov = sum (rowSums (matchesMat[,c(i,z)]) == 2) / sum (colSums(matchesMat[,c(i,z)]))
        cooc[i,z] = ov
        ov_size[i,z] = sum(rowSums (matchesMat[,c(i,z)]) == 2)
      #  }      
      }
    }
  
  colnames (cooc) = colnames(matchesMat)
  rownames (cooc) = colnames(matchesMat)
  diag (cooc) = 0
  cooc#[rowSums(cooc) >0,rowSums(cooc) >0]
  colnames (ov_size) = colnames(matchesMat)
  rownames (ov_size) = colnames(matchesMat)
  ov_size
  
#cooc_l[[celltype]][lower.tri (cooc_l[[celltype]])]
# cooc_diff = cooc[[1]] - cooc_l[[2]]
# cooc_diff[is.na(cooc_diff)] = 0
# ov_size_max = pmin(ov_size_l[[1]], ov_size_l[[2]])
# diag(ov_size_max) = 0

cooc_hm = Heatmap (
    cooc,
    column_names_rot =45, 
    row_names_gp = gpar(fontsize = 3),
    column_names_gp = gpar(fontsize = 3),
    col =palette_cooccurrence_cor, 
    cluster_rows=T,
    cluster_columns = T#,
#rect_gp = gpar (col = "white", lwd = 1)
)

pdf (file.path ('Plots','active_TF_cooccurence_momac_DAP_peaks_heatmaps.pdf'), width = 16,height=16)
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




















### Compare TF regulon expression with inflammation regulatory module to identify which TF might be driving ####
auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))
auc_mtx = auc_mtx[, colnames(auc_mtx) %in% active_TF]

srt$mod_2 = unname(rowMeans (auc_mtx))

# # Try with ridge plots regulons ####
library (ggridges)
library (ggplot2)
library (viridis)
library (tidyr)
#library(hrbrthemes)

# Plot
ccomp = as.data.frame (srt@meta.data)
#ccomp = ccomp[ccomp$cnmf_celltypes %in% c('cDCs'),]
cell_subsets_order = c("TAM_interstitial","TAM_MARCO","cDCs","TAM_TREM2","TAM_CXCLs","Mono_CD16","Mono_CD14")
ccomp$cnmf_celltypes = factor (ccomp$shared_cnmf_r_max, levels = cell_subsets_order)
ccomp$module = srt$mod_2
#ccomp = ccomp[!is.na(ccomp$cnmf_celltypes),]
rp <- ggplot(ccomp, aes(x = module, y = cnmf_celltypes, fill = ..x..)) +
  geom_density_ridges_gradient(
  scale = 3,
  rel_min_height = 0.01,
  linewidth = 0.4,
  color='white',
  alpha = 0.3
) +

  scale_fill_gradientn (colors = palette_expression) +  # Optional: nice color gradient
  theme_ridges() +                      # Optional: clean ridge plot theme
  theme(legend.position = "right")     # Adjust legend position
#   theme_classic() + facet_wrap (~sample, ncol=5)
pdf (file.path ('Plots','cnmf_inflammation_module_filtererd_ridge_plots.pdf'), width = 5,height=3)
rp
dev.off()



# # Try with ridge plots deviations ####
# Add metacolumns of average TF modules activity ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
#mMat = scale(as.matrix(mMat))#[selected_TF,])

# # Filter by RNA expression ####
metaGroupName = 'celltype_lv2'
min_exp = .1
active_TFs = exp_genes (srt, rownames(mMat), min_exp = min_exp, metaGroupName)
mMat = t(scale (mMat[active_TFs, ]))

infl_active = colMeans (t(scale(mMat))[rownames(t(mMat)) %in% active_TF,])
#tf_module_infl = colMeans (t)

archp@cellColData$infl_active = infl_active

library (ggridges)
library (ggplot2)
library (viridis)
library (tidyr)
#library(hrbrthemes)

# Plot
ccomp = as.data.frame (archp@cellColData)
#ccomp = ccomp[ccomp$celltype_lv2 %in% c('cDCs'),]
#ccomp$celltype_lv2 = factor (ccomp$celltype_lv2, levels = rev(cell_subsets_order))
ccomp$module = ccomp$infl_active
median_order = sort (unlist(lapply (split (ccomp$module, ccomp$celltype_lv3), function(x) median(x))))
ccomp$celltype_lv3 = factor (ccomp$celltype_lv3, levels = names (median_order))
rp <- ggplot(ccomp, aes(x = module, y = celltype_lv3, fill = ..x..)) +
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
pdf (file.path ('Plots','scatac_cnmf_inflammation_module_filtered_ridge_plots.pdf'), width = 7,height=5)
rp
dev.off()


### scatterplot ####
library(ggplot2)

# assuming your data frame is called df
# df <- your dataframe

gp = ggplot(ccomp, aes(x = mod_1, y = mod_2, color = celltype_lv3)) +
  geom_point(size = .8, alpha = 0.4) +
  geom_density_2d(size = 0.7, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = palette_myeloid) +
  labs(
    x = "Module 1 score",
    y = "Module 2 score",
    color = "Cell type",
    title = "mod_1 vs mod_2 with Cell-Type Density Contours"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

pdf (file.path ('Plots','regulon_inflammation_module_scatter_plots.pdf'), width = 9,height=9)
gp
dev.off()


infl_mod = c(
  "ATF3",   
  "BCL11A", 
  "BCL3",   
  "CEBPB",  
  "CREM",   
  "FOSB",   
  "FOSL1",  
  "FOSL2",
  "GTF2B",
  "IRF1",
  "IRF4",   
  "JUN",    
  "NFE2L3", 
  "NFKB1",  
  "NFKB2",  
  "PRDM1",
  "REL",    
  "RUNX3",  
  "SATB1",  
  "SPIB",   
  "TCF7L2", 
  "ZBTB7A")

ccomp = mMat[,colnames(mMat) %in% infl_mod]
ccomp = data.frame (
  ccomp, 
  celltype_lv3 = archp$celltype_lv3[match(rownames(ccomp),rownames(archp@cellColData))])
#ccomp = ccomp[ccomp$celltype_lv2 %in% c('cDCs'),]
#ccomp$celltype_lv2 = factor (ccomp$celltype_lv2, levels = rev(cell_subsets_order))
#median_order = sort (unlist(lapply (split (ccomp$module, ccomp$celltype_lv3), function(x) median(x))))
#ccomp$celltype_lv3 = factor (ccomp$celltype_lv3, levels = names (median_order))
rp <- lapply (infl_mod[infl_mod %in% colnames(ccomp)], function(x) ggplot(ccomp, aes_string(x = x, y = 'celltype_lv3', fill = '..x..')) +
  geom_density_ridges_gradient(
  scale = 3,
  rel_min_height = 0.01,
  linewidth = 0.4,
  color='white',
  alpha = 0.3
) +

  scale_fill_viridis_c(option = "C") +  # Optional: nice color gradient
  theme_ridges() +                      # Optional: clean ridge plot theme
  theme(legend.position = "right"))     # Adjust legend position
#   theme_classic() + facet_wrap (~sample, ncol=5)
pdf (file.path ('Plots','regulon_inflammation_module_ridge_plots.pdf'), width = 20,height=20)
wrap_plots(rp)
dev.off()

### Check genes around AP1:NFKB1 composite motif  ####
motif = 'pos_patterns.pattern_42'
finemo_hits = read.table(file.path(chromBPdir,celltype,'no_bias_model','finemo_out_counts','hits.tsv'), sep='\t', header=T)
ps = getPeakSet(archp)
composite_hits = finemo_hits[finemo_hits$motif_name == motif,]
composite_hits_gr = makeGRangesFromDataFrame (composite_hits)
composite_hits$gene = unname(ps$nearestGene[findOverlaps(composite_hits_gr, ps, select = 'first')])

composite_hits[queryHits (findOverlaps (composite_hits_gr, hubs_obj$hubsCollapsed)),]
as.character(composite_hits_gr[which(composite_hits$gene == 'LMO1')])
hubs_obj$hubs_id[queryHits (findOverlaps (hubs_obj$hubsCollapsed, makeGRangesFromDataFrame (composite_hits[which(composite_hits$gene == 'CD274'),])))]




### Get genes for regulons associated to inflammatory module ####
scenic_genes = read.csv (file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs','motifs.csv'), header=T, skip=1)

genes <- sapply(seq(nrow(scenic_genes)), function(x) unlist(regmatches(scenic_genes[x,9], gregexpr("(?<=\\(')[A-Za-z0-9_-]+", scenic_genes[x,9], perl = TRUE))))
names (genes) = scenic_genes[,1]

gene_regulon = 'FOSL2'
ps = as.data.frame (log2(AverageExpression (srt, group.by = 'celltype_lv3', features = unique(unlist(genes[names(genes) %in% gene_regulon])))[[1]]+1))

srt = ModScoreCor (
    seurat_obj = srt, 
    geneset_list = list(gene_regulon = rownames(ps)), 
    cor_threshold = NULL, 
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'TF_module', outdir = NULL)

ccomp = srt@meta.data[,c('gene_regulon','celltype_lv3')]
gp = ggplot(ccomp, aes(x = celltype_lv3, y = gene_regulon, fill = celltype_lv3)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = "", y = "Expression", title = "Expression distribution per cell type")

pdf (file.path ('Plots','NFKB1_gene_regulons_score.pdf'))
gp
reductionName='sampleID_harmony_umap'
DimPlot (srt, group.by = 'celltype_lv3', reduction=reductionName)
fp (srt, 'gene_regulon')
dev.off()


### Check Genescore of relative regulon ####
gs = getFeatures (archp, 'GeneScoreMatrix')
genes_gs = unique(unlist(genes[names(genes) %in% gene_regulon]))
archp@cellColData = archp@cellColData[,!colnames(archp@cellColData) %in% 'regulon']
archp = addModuleScore (
    ArchRProj = archp,
    useMatrix = 'GeneScoreMatrix',
    name = '',
    features = list(regulon = genes_gs[genes_gs %in% gs]),
    nBin = 25,
    nBgd = 100,
    seed = 1,
    threads = getArchRThreads(),
    logFile = createLogFile("addModuleScore")
  )
colnames (archp@cellColData) = gsub ('^\\.','',colnames(archp@cellColData))

df_long = as.data.frame (archp@cellColData[,c('regulon','celltype_lv3')])
# Plot boxplots for each column
gp = ggplot(df_long, aes(x = celltype_lv3, y = regulon, fill = celltype_lv3)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = "", y = "Expression", title = "Expression distribution per cell type")

pdf (file.path ('Plots','NFKB1_genescore_regulons.pdf'))
gp
dev.off()

pdf()
p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = 'cellColData', 
    name = 'regulon', 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path ('Plots','regulon_featureplots.pdf'), width = 18,height=14)
wrap_plots (p2)
dev.off()


### Check accessibilty at peaks associated with genes 
ps = getPeakSet (archp)
ps_sub = ps[ps$nearestGene %in% unique(unlist(genes[names(genes) %in% 'NFKB1']))]

mat = 'Peak'
metaGroupName='celltype_lv3'
pmat = ArchR::getMatrixFromProject (archp, useMatrix = paste0(mat,'Matrix'))
pmat = pmat[queryHits (findOverlaps (pmat, ps_sub))]
pmat = pmat[,rownames(archp)]
pmat = as.data.frame (t(assays (pmat)[[1]]))
pmat$metaGroup = paste0(as.character(archp@cellColData[,'Sample']),'__',as.character(archp@cellColData[,metaGroupName]))
pmat = aggregate (.~ metaGroup, pmat, mean)
rownames (pmat) = pmat[,1]
pmat = pmat[,-1]
pmat = rowMeans (pmat)
pmat = as.data.frame (pmat)
pmat$sample = sapply (rownames (pmat), function(x) unlist(strsplit(x,'__'))[1])
pmat$celltype = sapply (rownames (pmat), function(x) unlist(strsplit(x,'__'))[2])

#pmat$celltype = factor (pmat$celltype, levels = c('CD8_exhausted','NK_KLRC1','Tregs','NK_FGFBP2','CD8','CD4'))
bp = ggplot (pmat, aes (x = celltype, y = pmat, fill = celltype)) + 
      #geom_violin (trim=TRUE, aes (fill = treatment), alpha=.6) +
      #geom_violin (aes_string(fill = metaGroupNames[3])) +
      geom_point (aes (x = celltype, y = pmat), position='identity', alpha=.7, color="blue", size=1.2) +
      geom_boxplot(width=0.5, aes (fill = celltype, color = pmat), alpha=.6) +
      scale_fill_manual (values= palette_myeloid) + 
      scale_color_manual (values= palette_myeloid) + 
      geom_line (data = pmat, aes(x = celltype, y = pmat, group = sample), color='blue',linewidth=.2, alpha=.7) +
      gtheme + 
      NoLegend()
      #facet_wrap (~module, scales = 'free_y',strip.position = "left") +
      #theme(strip.background = element_rect(fill = "white", color = "black")) 

pdf (file.path ('Plots','NFKB1_genes_regulon_per_sample_celltype_boxplot.pdf'),height=3, width =3)
bp
dev.off()

### Get all peaks with NFKB1 motif ####


mat = 'Peak'
metaGroupName='celltype_lv3'
#pmat = ArchR::getMatrixFromProject (archp, useMatrix = paste0(mat,'Matrix'))
#motifMat = getPositions (archp)
matches = getMatches (archp)
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

all (as.character(rowRanges(matches)) == as.character (makeGRangesFromDataFrame (rowData(pMats))))
idx <- match(as.character(rowRanges(matches)), as.character(makeGRangesFromDataFrame (rowData(pMats))))
#ov = findOverlaps (makeGRangesFromDataFrame (rowData(pMats)), rowRanges(matches))
#reorder_pmat = subjectHits (ov)
pMatso = pMats[idx,]
all (as.character(rowRanges(matches)) == as.character (makeGRangesFromDataFrame (rowData(pMatso))))

gene_regulon = c('NFKB1')
ps = getPeakSet (archp)
#ps = ps[ps$peakType == 'Distal']
ps_sub = ps[ps$nearestGene %in% unique(unlist(genes[names(genes) %in% gene_regulon]))]

gene_TF = sapply (gene_regulon, function(x) colnames (matches)[grep(x,colnames(matches))])
if (length(gene_TF) > 1) {
  idxs = rowSums(assays(matches)[[1]][,gene_TF]) == length(gene_TF)
} else {
  idxs = assays(matches)[[1]][,gene_TF]
}
pmat_sub = pMatso[idxs,]

# # Find DAP ####
#force = FALSE
force=T
metaGroupName='inflamed'
if (!file.exists (paste0('DAP_',metaGroupName,'.rds')) | force)
  {
  DAP_list = getMarkerFeatures (
    ArchRProj = archp, 
    testMethod = "wilcoxon",
          useGroups = 'inflamed',
          bgdGroups = 'non_inflamed',
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
#names (DAP_res_l) = celltypes

pmat_gr = makeGRangesFromDataFrame (rowData(pmat_sub))
pmat_sub = pmat_sub[queryHits(findOverlaps(pmat_gr, DAP_res_l[[1]]))]
ps_regulon = queryHits (findOverlaps (makeGRangesFromDataFrame (rowData(pmat_sub)), ps_sub))
ps_regulon_idx = seq(nrow (pmat_sub)) %in% ps_regulon
pmat_sub_mat = assay(pmat_sub)

pmat_sub_mat2 = as.data.frame (pmat_sub_mat)
pmat_sub_mat2$regulon = ps_regulon_idx
df_long <- pmat_sub_mat2 %>%
  pivot_longer(cols = -regulon,
               names_to = "celltype",
               values_to = "expression")
df_long$celltype = factor(df_long$celltype, levels = rev(cell_subsets_order))
# Plot boxplots for each column
gp = ggplot(df_long, aes(x = celltype, y = expression, fill = regulon)) +
  geom_boxplot() + gtheme + facet_wrap (~regulon)

pdf (file.path ('Plots','NFKB1_peaks_gene_regulons_boxplot.pdf'))
gp
dev.off()

pdf (file.path ('Plots','NFKB1_regulon_peaks_vs_rest_heatmap.pdf'))
Heatmap (t(scale(t(pmat_sub_mat[ps_regulon_idx,]))), col=palette_fragments, 
  row_names_gp = gpar(fontsize = 0))
Heatmap (t(scale(t(pmat_sub_mat[!ps_regulon_idx,]))), col=palette_fragments, 
  row_names_gp = gpar(fontsize = 0))
dev.off()



## Find most abundant composites 
chromBPdir='/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/chromBPnet'
celltype='inflamed'
finemo_res = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
table (finemo_res$V6)[order(table (finemo_res$V6))]
head (finemo_res$V4[finemo_res$V6 == 'pos_patterns.pattern_33'],1)




