conda activate meso_scatac

R

set.seed(1234)

####### ANALYSIS of Myeloid compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

# Set # of threads and genome reference ####
addArchRThreads(threads = 1) 
addArchRGenome("hg38")

# Load ArchR project ####
archp = loadArchRProject (projdir)

# Add metadata ####
archp$status = ifelse (grepl ('^P',archp$Sample2), 'tumor','normal')  

## Subset only for tumor samples ####
archp = archp[archp$Sample2 %in% c('P1','P10','P11','P12','P13','P14','P3','P4','P5','P8')]

# Load RNA ####
srt = readRDS (file.path('..','scrna','srt.rds'))

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


pdf (file.path('Plots','harmony_celltype_umap.pdf'),5,width = 10)
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

wrap_plots (umap_p6, umap_p5,umap_p4)
dev.off()


### Call peaks on celltypes ####
pdf(file.path('Plots','peakcalls.pdf'))
metaGroupName = 'celltype2'
force=TRUE
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
  

# Differential Accessed motifs ####
metaGroupName = "celltype2"
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
metaGroupName = 'Clusters_H'
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

pdf (file.path ('Plots','DAM_with_rna_expression_heatmaps.pdf'), width = 8,height=4)
draw (DAM_hm + TF_exp_selected_hm + TF_exp_selected_hm2)
dev.off()


# Import shared cNMF from scRNA and make module scores ####
shared_cnmf = readRDS (file.path('..','scrna','shared_cnmf_myeloid.rds'))
archp = addModuleScore (
    ArchRProj = archp,
    useMatrix = 'GeneScoreMatrix',
    name = "shared_cnmf",
    features = shared_cnmf,
    nBin = 25,
    nBgd = 100,
    seed = 1,
    threads = getArchRThreads(),
    logFile = createLogFile("addModuleScore")
  )
pdf()
archp = addImputeWeights (archp)
p <- plotEmbedding (
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = colnames(archp@cellColData)[grep('shared_cnmf',colnames(archp@cellColData))], 
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)

pdf (file.path ('Plots','shared_cnmf.pdf'),10,10)
wrap_plots (p, ncol=3)
dev.off()


# Correlate module scores with TFs ####
metaGroupName = "celltype2"
force=FALSE
source (file.path('..','..','git_repo','utils','DAM.R'))
if (!exists('mSE')) gsSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.data.frame (t(as.matrix(scale(mMat[active_TFs,]))))
cnmfs = as.data.frame (archp@cellColData[,grep ('shared_cnmf', colnames(archp@cellColData))])
cnmfs = t(scale (t(cnmfs)))
tf_cnmf_df = cor (cnmfs, mMat)
top_TFs = colnames(tf_cnmf_df)[order(-tf_cnmf_df['shared_cnmf.8',])]
top_TFs = c(head (top_TFs,10), tail (top_TFs,10))
ha = HeatmapAnnotation(' ' = anno_mark(at = match(top_TFs, colnames(tf_cnmf_df)), 
    labels = top_TFs, labels_gp = gpar(fontsize = 6, fontface='italic')))
cnmf_tf_df = Heatmap (tf_cnmf_df,
  column_dend_side = 'bottom',
  top_annotation = ha,
  row_names_gp = gpar (fontsize = 8),
  col = rev(palette_deviation),
  column_km = 2,
  #rect_gp = gpar(type = "none"),
  border=T,
  column_names_gp = gpar(fontsize = 0, fontface='italic'))

pdf (paste0 ('Plots/active_TF_cnmf_cor_heatmap.pdf'), width = 6,height=3)
cnmf_tf_df
dev.off()


# Check distribution of cnmf modules across samples ####
cnmfs = as.data.frame (archp@cellColData[,grep ('shared_cnmf', colnames(archp@cellColData))])
cnmfs$sample = archp$Sample
# cnmfs_agg = aggregate (cnmfs, by = list(sampleID= archp$Sample), mean)
# rownames (cnmfs_agg) = cnmfs_agg$sampleID
# cnmfs_agg = cnmfs_agg[,-1]

# FRIP = as.data.frame (archp$FRIP)
# frip_agg = aggregate (FRIP, by = list(sampleID= archp$Sample), mean)   


# cnmfs_agg_norm = cnmfs_agg / frip_agg[,2]
# # Check if they correlate to seq-deph
# p.adjust(sapply (seq(ncol(cnmfs_agg)), function(x) cor.test (cnmfs_agg[[x]], frip_agg[[2]])$p.value))

# cnmfs_agg$sample = rownames(cnmfs_agg)
cnmfs_df = gather (cnmfs, cnmf, score, 1:(ncol (cnmfs)-1))

cp = ggplot (cnmfs_df, aes (x= cnmf, y = score, fill = sample)) + geom_boxplot () + gtheme + 
scale_fill_manual (values = palette_sample)
pdf (file.path ('Plots','cnmf_score_sample_boxplots.pdf'), width=12)
cp
dev.off()








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



### Compute correlation of myeloid modules with hubs
archp_meta = as.data.frame (archp@cellColData)
cnmfs = as.data.frame (archp@cellColData[,grep ('shared_cnmf', colnames(archp@cellColData))])
#cnmfs = scale (t(cnmfs))

module_driver = 'shared_cnmf.6'
top_hubs = 30
traj_sample = list()
sams = unique(archp$Sample2)
#for (sam in unique(archp$Sample2))
#    {
    sam = sams  
    library(zoo)
    bin_width <- 100   # Number of observations per bin
    overlap <- 10    
    cnmf_ordered_sample = as.data.frame(scale(t(cnmfs[archp_meta$Sample %in% sam,])))
    cnmf_ordered_sample = cnmf_ordered_sample[, order(unlist(cnmf_ordered_sample[module_driver,]))]
    hubs_sample = scale(hubsCell_mat[,archp_meta$Sample %in% sam])[,colnames(cnmf_ordered_sample)]

    cnmf_ordered_sample = as.data.frame(lapply(as.data.frame(t(cnmf_ordered_sample)), function(x) {
      zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
    }))
    hub_ordered_sample = as.data.frame(lapply(as.data.frame(t(hubs_sample)), function(x) {
      zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
    }))
    colnames (hub_ordered_sample) = rownames(hubs_sample)

    cor_cnmf_hub = cor (cnmf_ordered_sample, hub_ordered_sample, method = 'spearman')
    cor_order = order (-cor_cnmf_hub[module_driver,])
    hub_ordered_sample = hub_ordered_sample[, cor_order]
    hub_ordered_sample = as.data.frame (t(hub_ordered_sample))
    hub_ordered_sample = rbind (head(hub_ordered_sample,top_hubs), tail (hub_ordered_sample, top_hubs))
    row_names = hubs_obj$hubsCollapsed$gene[match(rownames(hub_ordered_sample), hubs_obj$hubs_id)]
    traj_sample = Heatmap (
      t(scale(t(hub_ordered_sample))), 
      col = rev(palette_hubs_accessibility), 
      cluster_columns=F,
    #  name = sam,
      cluster_rows=F,
      column_names_gp = gpar(fontsize = 0),
      row_names_gp = gpar(fontsize = 7, fontface='italic'),
      border=T)
    #}

pdf (file.path('Plots',paste0('hubs_cor_to_cnmf_',module_driver,'_all.pdf')), height=6, width=3)
traj_sample
dev.off()

pdf (file.path('Plots',paste0('hubs_cor_to_cnmf_',module_driver,'_sample.pdf')), height=8, width=3)
traj_sample
dev.off()


### Check TF enrichment in HUBS correlated with TREM2 vs anticorrelated ####
### TF Enrichment in peaks in large hubs ####
tf_match = getMatches (archp)
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)
hubs_to_enrich = hubs_obj$hubsCollapsed[match (rownames (hub_ordered_sample), hubs_obj$hubs_id)]
hubs_to_enrich_peaks = bg_peakSet[queryHits(findOverlaps(bg_peakSet, hubs_to_enrich))]
#tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
hubs_TF_res =  hyperMotif (
  selected_peaks = hubs_to_enrich_peaks, 
  motifmatch = tf_match)


### Group macrophages by TREM2 module score ####
#module_driver = 'shared_cnmf.6'
cnmfs = as.data.frame (archp@cellColData[,grep ('shared_cnmf', colnames(archp@cellColData))])
#cnmfs = as.data.frame (do.call (rbind, lapply (unique(archp$Sample), function(x) t(scale(t(cnmfs[archp$Sample == x,]))))))
cnmfs = as.data.frame (t(scale(t(archp@cellColData[,grep ('shared_cnmf', colnames(archp@cellColData))]))))

TREM2_state = cnmfs$shared_cnmf.6
names (TREM2_state) = rownames(archp@cellColData)
TREM2_state = TREM2_state[order(-TREM2_state)]

n <- length(TREM2_state) # Get the length of the vector
part_size <- ceiling(n / 5) # Calculate the size of each part

# Create indices for splitting
splits <- split(TREM2_state, ceiling(seq_along(TREM2_state) / part_size))
labeled_splits <- setNames(splits, c("very_high","high", "mid", "low","very_low"))

df <- data.frame(
  value = unlist(labeled_splits),
  label = rep(names(labeled_splits), times = sapply(labeled_splits, length))
)

archp$TREM2_state = df$label[match(rownames(archp@cellColData), sapply(rownames(df), function(x) unlist(strsplit(x,'\\.'))[2]))] 
#table (archp$TREM2_state2, archp$TREM2_state)

# Export bigiwg files ####
metaGroupName = 'TREM2_state'
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











