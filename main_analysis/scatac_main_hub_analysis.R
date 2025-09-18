# Load functions for hub detection
source (file.path('..','..','git_repo','utils','knnGen.R'))
source (file.path('..','..','git_repo','utils','addCoax.R'))
source (file.path('..','..','git_repo','utils','Hubs_finder.R'))
source (file.path('..','..','git_repo','utils','hubs_track.R'))


# Export bigiwg files ####
#archp$celltype_status = paste0(archp$celltype2, '_', archp$status)
metaGroupName = 'celltype_lv1'
archp$celltype_lv1a = archp$celltype_lv1
archp$celltype_lv1a[archp$celltype_lv1a == 'Fibroblasts' & archp$Sample == 'P1'] = 'Fibroblasts_P1'

metaGroupName = 'celltype_lv1a'
exp_bigwig = F
if (exp_bigwig)
  {
  getGroupBW(
    ArchRProj = archp[archp$celltype_lv1 %in% c('Mesothelium','Fibroblasts','Malignant')],
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


if (!file.exists ('peak_regions.bed'))
  {
  peak_regions = as.data.frame (getPeakSet (archp), row.names=NULL)
  peak_regions = peak_regions[,c(1:3)]
  write.table (peak_regions, file.path('peak_regions.bed'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  }


### Hubs analysis #####
metaGroupName = "celltype_lv1"
cor_cutoff = 0.3
#max_dist = 12500
max_dist = 12500
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
dir.create(file.path (hubs_dir, 'Plots'), recursive=T)


# Generate cluster-aware knn groups ####
k= 30
metaGroupName = 'celltype_lv1'

force = FALSE
if (!file.exists(file.path (hubs_dir, paste0 ('KNNs_',metaGroupName,'k_',k,'.rds'))) | force)
  {
  KNNs = knnGen (
    ArchRProj = archp, 
    k = k,
    reducedDims = 'IterativeLSI', 
    group = metaGroupName,
    overlapCutoff = 0.6,
    #cellsToUse = metaGroup_df$barcode,
    #min.cells_in_group = min_cells,
    min_knn_cluster = 2
    )  
  saveRDS (KNNs, file.path (hubs_dir, paste0 ('KNNs_',metaGroupName,'k_',k,'.rds')))  
  } else {
  KNNs = readRDS (file.path (hubs_dir, paste0 ('KNNs_',metaGroupName,'k_',k,'.rds')))
  }

# Make KNN data.frame ####
KNNs_df = lapply (seq_along(KNNs), function(x) data.frame (
  cell = KNNs[[x]], 
  group=paste0('KNN_',x), 
  group2 = unlist(strsplit (names(KNNs)[x],'KNN'))[1]))

KNNs_df = do.call (rbind, KNNs_df)
archp$knn_groups = KNNs_df$group[match (archp$cellNames, KNNs_df$cell)]

# Plot KNNs on UMAP ####
umap_knn = plotEmbedding (ArchRProj = archp, embedding = 'UMAP_H',
  colorBy = "cellColData", name = 'knn_groups',plotAs ='hex',
    baseSize=0, labelMeans=FALSE) + NoLegend() 
pdf (file.path(hubs_dir, 'Plots',paste0('knn_',k,'.pdf')), height=5, width=5)
umap_knn
dev.off()



### Run peak calling to define peak background ####
metaGroupName = "celltype_lv1"
force=FALSE
peak_reproducibility=2
if(!all(file.exists(file.path('PeakCalls', paste0(unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds')))) | force) source ('../../git_repo/utils/callPeaks.R')
  

# Run Co-accessibility ####
run_coax = TRUE
if (run_coax)
  {
  archp = addCoAx (
    archp, 
    KNNs,
    maxDist = max_dist)
  }

### Run hub finder ####
force=F
if (!file.exists (file.path(hubs_dir,'global_hubs_obj.rds')) | force)
  {
  hubs_obj = hubs_finder (
    ArchRProj = archp, 
    group_by = NULL,
    cor_cutoff = cor_cutoff,
    #select_group = metaGroup_df$metaGroup,
    cor_FDR = 1, 
    cor_var = 0,
    min_peaks = min_peaks,
    macs_score = 1,
    dgs = dgs,
    cores=1,
    remove_chr = c('chrX')
    ) 
  saveRDS (hubs_obj, file.path(hubs_dir,'global_hubs_obj.rds'))
  hubs_regions = as.data.frame (hubs_obj$hubsCollapsed)
  hubs_regions$width = hubs_obj$hubs_id
  write.table (hubs_regions, file.path(hubs_dir, 'hub_regions.bed'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  } else {
  hubs_obj = readRDS (file.path(hubs_dir,'global_hubs_obj.rds'))  
  }



# Generate matrix of fragment counts of hubs x metagroup ####
metaGroupName = 'celltype_lv1'
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
  frags_in_sample = sapply (unique(archp@cellColData[,metaGroupName]), function(x) sum (archp$nFrags[as.character(archp@cellColData[,metaGroupName]) == x]))
  hubsSample_mat = t(t(hubsSample_mat) * (10^6 / frags_in_sample)) # scale
  saveRDS (hubsSample_mat, file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))
  } else {
  hubsSample_mat = readRDS (file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))  
  }
hubsSample_mat = as.data.frame (hubsSample_mat)

ha = HeatmapAnnotation (size = anno_barplot(width (hubs_obj$hubsCollapsed), gp = gpar(color = "red"), height =  unit(8, "mm")))

TF = 'WT1'
int_genes = which (rownames(hubsSample_mat) %in% hubs_obj$hubs_id[grepl (TF, hubs_obj$hubsCollapsed$gene)])
int_genes_label = paste0(hubs_obj$hubs_id, '_', hubs_obj$hubsCollapsed$gene)[grepl (TF, hubs_obj$hubsCollapsed$gene)]
ha2 = HeatmapAnnotation(foo = anno_mark(at = int_genes, side = 'column',
  labels_rot = 45,
    labels = int_genes_label, labels_gp = gpar(fontsize = 7, fontface='italic')))
hm = Heatmap (
  scale (t(hubsSample_mat)), 
  bottom_annotation = ha2, 
  column_names_gp = gpar(fontsize = 0),
  show_column_dend = F,
  #row_dend_width = unit(5,'mm'),
  row_dend_side = 'left',
  col = rev(palette_hubs_accessibility),
  border=T,
  name = 'Hubs')
pdf (file.path (hubs_dir,'Plots',paste0('hubs_',metaGroupName,'_heatmap.pdf')), height=4)
hm
dev.off()


# Generate matrix of fragment counts of hubs x barcodes ####
metaGroupName = 'celltype_lv1'
if (!file.exists(file.path (hubs_dir, paste0('hubs_cells_',metaGroupName,'_mat.rds'))))
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
  saveRDS (hubsCell_mat, file.path (hubs_dir,paste0('hubs_cells_',metaGroupName,'_mat.rds')))
  } else {
  hubsCell_mat = readRDS (file.path (hubs_dir,paste0('hubs_cells_',metaGroupName,'_mat.rds')))  
  }
hubsCell_mat = as.data.frame (hubsCell_mat)


# Compute differential hub accessibility DHA ####
library (presto)
metaGroupName = 'celltype_lv1'
all (colnames(hubsCell_mat) == rownames(archp@cellColData))
res = wilcoxauc (log2(hubsCell_mat+1), as.character (archp@cellColData[,metaGroupName]))

res_l = lapply (split (res, res$group), function(x){
  tmp = x[order (x$padj),]
  tmp
})


# celltype specific hubs to input to IGV ####
c(HUB277 HUB77 HUB68 HUB32 HUB522 HUB5713 HUB26 HUB48 HUB6929 HUB1459 HUB16 HUB2371)

HUB522 HUB5713 HUB277 HUB32 HUB68 HUB1459 HUB26 HUB16 HUB48 HUB77 HUB6929 HUB2371


# Generate matrix of fragment counts of hubs x metagroup ####
table (archp$celltype_lv1)
archp$compartment = archp$celltype_lv1
archp$compartment[archp$celltype_lv1 %in% c('Alveolar','Mesothelium','Malignant')] = 'epithelial'
archp$compartment[archp$celltype_lv1 %in% c('Endothelial','Fibroblasts','SmoothMuscle')] = 'stroma'
archp$compartment[archp$celltype_lv1 %in% c('Myeloid','pDCs')] = 'myeloid'
archp$compartment[archp$celltype_lv1 %in% c('B_cells','NK','Plasma','T_cells')] = 'lymphoid'

metaGroupName = 'compartment'
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
  frags_in_sample = sapply (unique(archp@cellColData[,metaGroupName]), function(x) sum (archp$nFrags[as.character(archp@cellColData[,metaGroupName]) == x]))
  hubsSample_mat = t(t(hubsSample_mat) * (10^6 / frags_in_sample)) # scale
  saveRDS (hubsSample_mat, file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))
  } else {
  hubsSample_mat = readRDS (file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))  
  }
hubsSample_mat = as.data.frame (hubsSample_mat)

# Find shared hubs across cell types ####
hubsSample_rank = hubsSample_mat / width(hubs_obj$hubsCollapsed)
hubsSample_rank = apply (hubsSample_rank, 2, function(x) order(-x))
rownames (hubsSample_rank) = rownames(hubsSample_mat)
shared_hubs = apply (hubsSample_rank, 1, mean)
shared_hubs = shared_hubs[order(shared_hubs)]
which (names (shared_hubs) == 'HUB322')

# Plot some of shared hubs ####
hubsCell_mat_t = as.data.frame (t(hubsCell_mat))
shared_hubs = names(head (shared_hubs[order(shared_hubs)]))
hubsCell_mat_t = hubsCell_mat_t[rownames(archp@cellColData),]
archp@cellColData = cbind(archp@cellColData, hubsCell_mat_t[,shared_hubs])
umap_p1 =  plotEmbedding (
  ArchRProj = archp, 
  colorBy = "cellColData",
 name = shared_hubs, 
 embedding = "UMAP",
 imputeWeights = getImputeWeights(archp)
 )
  
pdf (file.path(hubs_dir,'Plots','shared_hubs_umap.pdf'), 15,15)
wrap_plots (umap_p1, ncol=4)
dev.off()

# Check if VIM is also expressed across cell types in scRNA ####
# Load RNA 
srt = readRDS ('../scrna/srt.rds')
#sarc_order = read.csv ('../scrna/cnmf20_sarcomatoid_sample_order.csv', row.names=1)

vp = VlnPlot (srt, feature = c('ERG','VIM'), group.by = 'celltype_simplified', col=palette_celltype_lv1) + NoLegend()

hubsCell_mat = t(hubsCell_mat)[rownames(archp@cellColData),]
archp@cellColData = cbind(archp@cellColData, log2(hubsCell_mat[,c('HUB322'), drop=F]+1))
#archp@cellColData = archp@cellColData[, -ncol(archp@cellColData)]
p <- plotGroups(
  ArchRProj = archp, 
  groupBy = "celltype_lv1", 
  colorBy = "cellColData", 
  name = "HUB322",
  plotAs = "violin",
  pal = palette_celltype_lv1,
  alpha = 0.4,
  addBoxPlot = TRUE
 )
pdf (file.path ('Plots','VIM_expression_accessibility.pdf'), height=5, width=12)
wrap_plots (vp, p, ncol=2)
dev.off()

### it is....check correlation of HUB322 with TFs ####

devMethod = 'ArchR'
 if (devMethod == 'ArchR')
    {
    TF_db='Motif'
    if (!exists ('mSE')) mSE = ArchR::getMatrixFromProject (archp, useMatrix = paste0(TF_db,'Matrix'))
    mSE = mSE[, archp$cellNames]
    rowData(mSE)$name = gsub ('_.*','',rowData(mSE)$name)
    rowData(mSE)$name = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", rowData(mSE)$name)
    }
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
all (colnames (mMat) == rownames(hubsCell_mat))
hub_exp = log2(hubsCell_mat[,'HUB322']+1)
zero_val = hub_exp != 0

hub_nz = hub_exp[which(zero_val)]
mMat_nz = t(as.matrix(mMat))[which(zero_val),]

archp_meta = archp@cellColData
archp_meta = archp_meta[which(zero_val),]
hub_tf_cor = lapply (unique(archp@cellColData[,'celltype_lv1']), function(x)
  as.data.frame (t(cor (hub_nz[match(rownames(archp_meta)[as.character(archp_meta$celltype_lv1) == x],
    names(hub_nz))], mMat_nz[match(rownames(archp_meta)[as.character(archp_meta$celltype_lv1) == x],
    rownames(mMat_nz)),], method='spearman')))
  )
hub_tf_cor = lapply (hub_tf_cor, function(x) x[order(-x[,1]),, drop=F])
hub_tf_cor = as.data.frame (t(hub_tf_cor))
head (hub_tf_cor[order (-hub_tf_cor[,1]),,drop=F],10)


cor_df = data.frame (hub = hub_exp[which(zero_val)], tf = t(as.matrix(mMat))[which(zero_val),'ETS1'])
cor_df$celltype = archp$celltype_lv1[match(rownames(cor_df), rownames(archp@cellColData))]
cor_df$sample = archp$Sample[match(rownames(cor_df), rownames(archp@cellColData))]

cor (cor_df[cor_df$celltype == 'Malignant','hub'], cor_df[cor_df$celltype == 'Malignant','tf'])

cor_p = ggplot (cor_df[cor_df$celltype == 'T_cells',], aes (x = hub, y = tf, color = celltype)) + geom_point(alpha=.4)
cor_p2 = ggplot (cor_df[cor_df$celltype == 'T_cells',], aes (x = hub, y = tf, color = sample)) + geom_point(alpha=.4)

pdf (file.path('Plots','HUB322_VIM_correlation_ETS1.pdf'), width=9)
wrap_plots (cor_p,cor_p2)
dev.off()









# Compare DHA to list of SE to validate approach ####
hubs_obj = readRDS (file.path(hubs_dir,'global_hubs_obj.rds'))  
library (presto)
metaGroupName = 'celltype_lv1'
hubsCell_mat = readRDS (file.path (hubs_dir,paste0('hubs_cells_',metaGroupName,'_mat.rds')))  
all (colnames(hubsCell_mat) == rownames(archp@cellColData))
res = wilcoxauc (log2(hubsCell_mat+1), as.character (archp@cellColData[,metaGroupName]))


res_l = lapply (split (res, res$group), function(x){
  tmp = x[x$padj < 0.05,]
  tmp[tmp$logFC > 1,]
})



# Run hypergeometric test on SE database 2
path = system.file (package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain (path)
SE_path = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/SEdb_ALL/'
SE_regions = lapply (list.files(SE_path, pattern = '.bed'), function(x) {
  y = read.table (paste0(SE_path,x), header=T, sep='\t')
  y = y[,c(3:6)]
  colnames (y) = c('chr','start','end','rank')
  y = makeGRangesFromDataFrame (y)
  })
names(SE_regions) = list.files(SE_path, pattern = '.bed')
SE_regions = SE_regions[grepl ('ENCODE', names(SE_regions))]

# Read in also SE celltypes from Wooseung
SE_path2 = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/SE_db2_from_Wooseung'
SE_regions2 = lapply (list.files(SE_path2, pattern = '.bed'), function(x) {
  y = read.table (file.path(SE_path2,x), header=T, sep='\t')
  y = y[,c(1:3,6)]
  colnames (y) = c('chr','start','end','rank')
  y = makeGRangesFromDataFrame (y)
  })
names(SE_regions2) = list.files(SE_path2, pattern = '.bed')
# SE_regions2 = lapply (SE_regions2, function(x) 
#   {
#   seqlevelsStyle(x) = "UCSC"  # necessary
#   x = unlist (liftOver (SE_regions2, ch))
#   })

SE_regions = c (SE_regions, SE_regions2)
names (SE_regions) = sub ('.bed','',names(SE_regions))

hyper_SE_cluster2 = list()
for (celltype in unique(archp$celltype_lv1))
  {
  bg = readRDS (file.path('PeakCalls',paste0(celltype,'-reproduciblePeaks.gr.rds')))
  hubs_h38 = hubs_obj$hubsCollapsed[match(res_l[[celltype]]$feature, hubs_obj$hubs_id)]
  seqlevelsStyle(hubs_h38) = "UCSC"  # necessary
  hubs_h19 = unlist (liftOver (hubs_h38, ch))


  q = sapply (SE_regions, function(x) 
    {
    q = countOverlaps (hubs_h19, x) 
    q[q > 0] = 1
    sum (q)
    })
  k = length (hubs_h19)
  m = sapply (SE_regions, function(x) 
    {
    m = countOverlaps (bg, x) 
    m[m > 0] = 1
    sum(m)
    })
  n = length(bg) - m
message (paste ('Run hyper tests for',celltype,' hubs'))
hyper_SE = NULL
for (i in 1:length(SE_regions))
  {
  hyper_SE[i] = phyper (q = q[i],m = m[i],n = n[i],k = k, lower.tail = FALSE, log.p = FALSE)
  }
hyper_SE_cluster2 [[celltype]] = -log10(p.adjust(hyper_SE))
#hyper_SE_cluster2 [[celltype]][hyper_SE_cluster2 [[celltype]] < 0.05] = NA
}

hyper_SE_cluster_df = do.call (cbind, hyper_SE_cluster2)
rownames (hyper_SE_cluster_df) = names(SE_regions)
hyper_SE_cluster_df[is.infinite(hyper_SE_cluster_df)] = -99
hyper_SE_cluster_df[hyper_SE_cluster_df == -99] = max (hyper_SE_cluster_df)
hyper_SE_cluster_df[hyper_SE_cluster_df < -log10(0.05)] = 0
#hyper_SE_cluster_df[hyper_SE_cluster_df > 100] = 100
rownames(hyper_SE_cluster_df) = sub ('ENCODE_','',rownames(hyper_SE_cluster_df))
rownames(hyper_SE_cluster_df) = sub ('.bed','',rownames(hyper_SE_cluster_df))
SE_hm2 = Heatmap (hyper_SE_cluster_df,
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        cluster_rows = T,
        #col = pals_heatmap[[5]],
        cluster_columns=T,#col = pals_heatmap[[1]],
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 10),
        #rect_gp = gpar(col = "white", lwd = .5),
        border=TRUE, 
        col=viridis::viridis(100),
        column_names_rot = 45)
        #right_annotation = motif_ha
        
ct1 = c('Malignant','Mesothelium','Alveolar',
  'Fibroblasts','SmoothMuscle','Endothelial',
  'Myeloid','T_cells',#'NK',
  'B_cells','Plasma','pDCs')
ct2 = c('Mammary_gland_primary_cell_mammary-epithelial-cell',
  #'Skin_in_vitro_differentiated_cells_bipolar-neuron',
  'Skin_primary_cell_keratinocyte',
  'AT1','AT2',
  'Lung_primary_cell_fibroblast-of-lung',
  'Smoothmuscle',
  'Endothelial',
  'Blood_primary_cell_CD14-positive-monocyte',
  'Tcell',
  #'NK',
  'Blood_primary_cell_Bcell')        
rownames(hyper_SE_cluster_df)
SE_hm3 = Heatmap (as.data.frame(hyper_SE_cluster_df)[ct2,ct1],
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        cluster_rows = F,
        #col = pals_heatmap[[5]],
        cluster_columns=F,#col = pals_heatmap[[1]],
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 10),
        #rect_gp = gpar(col = "white", lwd = .5),
        border=TRUE, 
        col=viridis::viridis(100),
        column_names_rot = 45)
        #right_annotation = motif_ha
        

# pdf (file.path('Plots', 'hyper_hub_SE2.pdf'),width=5.5, height=5)
# SE_hm2
# dev.off()
pdf (file.path('Plots', 'hyper_hub_SE3.pdf'),width=5.5, height=3)
SE_hm3
dev.off()




#### Compute P2G ####
run_p2g_TF = TRUE

    maxDist = 250000
    archp = addPeak2GeneLinks(
        ArchRProj = archp,
        useMatrix = 'GeneScoreMatrix',
        reducedDims = "IterativeLSI",
        maxDist = maxDist
    )
  

### Show example of hub ####
TF = 'WT1'
pdf()
#archp$fetal_sample = paste0(archp$Sample, archp$fetal_group)
#metaGroupName = 'fetal_group'
metaGroupName = 'celltype_lv1'
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp,#[!archp$Sample3 %in% c('P11_HOX')], 
    sample_levels = celltype_order, 
    #ylim = c(0,0.01),
    groupBy = metaGroupName, 
    hubs_regions = hubs_obj$hubsCollapsed,
    #sample_levels = sample_sarc_order,
    minCells = 10,
    geneSymbol = TF,
    plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
    #pal = palette_sample,
    #pal = palette_fetal,
    threads=1,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (hubs_obj$hubsCollapsed[match(hub, hubs_obj$hubs_id)]),100000,100000),
    upstream = 150000,
    downstream = 100000,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    pal = palette_celltype_lv1,
    #loops = getCoAccessibility (archp, corCutOff = 0.3,
    #  returnLoops = TRUE),
    useGroups= NULL
)
dev.off()
plotPDF (meso_markers, ArchRProj = archp,height=3.5, width=6, name =paste0('MPM_markers_coveragePlots.pdf'),addDOC=F)
  

metaGroupName = 'celltype_simplified2'
top_dah = data.frame (
gene = srt@assays$RNA@data[TF,],
group = srt@meta.data[,metaGroupName])
top_dah$group = factor (top_dah$group, levels = rev(celltype_order))
top_dah = na.omit(top_dah)
bp = ggplot (top_dah, aes (x = gene, y = group, fill = group)) + 
vlp + 
bxpv + 
scale_fill_manual (values = palette_celltype_lv1) +
#geom_point (position='identity', alpha=.3, color="grey44", size=1) +
gtheme_no_rot

pdf (file.path ('Plots', paste0('scrna_',TF,'_boxplots.pdf')), height=4, width=4)
bp
dev.off()




### Find WT1 hub to check proprortion of chrombBPnet TFs between mesothelium and fibroblasts and maybe malignant ####
wt1_hubs = c(81,286,5951)
wt1_hub_peaks = do.call (rbind, hubs_obj$hubsClusters[[1]][wt1_hubs])

# Import chromBPnet finemo motifs for P1 too ####
#archp_P1 = archp[archp$Clusters %in% c('C2') & archp$Sample == 'P1']
library ('universalmotif')

# metaGroupName = 'Clusters2'
# if (!any (ls() == 'mSE')) mSE = fetch_mat (archp_P1, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData(mSE)$name
chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet'
# metaGroupName = 'celltype_lv1'


chrombpnet_counts = list()
celltypes = c('Mesothelium','Fibroblasts_P1','Malignant')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))
  chrombpnet_counts[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
  }

# Subset for WT1 hub peaks 
wt1_hub_peaks_gr = makeGRangesFromDataFrame (wt1_hub_peaks)
chrombpnet_counts_sub = lapply (chrombpnet_counts, function(x) 
{
 colnames (x) = c('chr','start','end') 
 x = makeGRangesFromDataFrame (x, keep.extra.columns=T)
 x = x[queryHits (findOverlaps (x, wt1_hub_peaks_gr))]
 x = as.data.frame (x)
 colnames (x) = paste0('V',1:ncol(x))
 x
})


chrombpnet_profile = list()
celltypes = c('Mesothelium','Fibroblasts_P1','Malignant')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  }

chrombpnet_profile = lapply (chrombpnet_profile, function(x) x[x$V5 != 'NaN_NaN_NaN',])

# Subset for WT1 hub peaks 
wt1_hub_peaks_gr = makeGRangesFromDataFrame (wt1_hub_peaks)
chrombpnet_profile_sub = lapply (chrombpnet_profile, function(x) 
{
 colnames (x) = c('chr','start','end') 
 x = makeGRangesFromDataFrame (x, keep.extra.columns=T)
 x = x[queryHits (findOverlaps (x, wt1_hub_peaks_gr))]
 x = as.data.frame (x)
 colnames (x) = paste0('V',1:ncol(x))
 x
})


library(dplyr)
library(ggplot2)


top_n <- 5
n <- length(chrombpnet_counts_sub)

bp_list <- lapply(seq_len(n), function(i) {
  tbl <- table(chrombpnet_counts_sub[[i]]$V7)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_counts_sub[[i]]$V6[chrombpnet_counts_sub[[i]]$V7 == tf][1]
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
bp_df$type = factor (bp_df$type, levels = c('Mesothelium','Malignant','Fibroblasts_P1'))
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF)) ) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

  

pdf (file.path ('Plots', 'TF_abundance_WT1_cHub_counts_barplot.pdf'),6,width=6.5)
bp
dev.off()


library(dplyr)
library(ggplot2)

top_n <- 5
n <- length(chrombpnet_profile_sub)

bp_list <- lapply(seq_len(n), function(i) {
  tbl <- table(chrombpnet_profile_sub[[i]]$V7)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_profile_sub[[i]]$V6[chrombpnet_profile_sub[[i]]$V7 == tf][1]
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
bp_df$type = factor (bp_df$type, levels = c('Mesothelium','Malignant','Fibroblasts_P1'))
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF)) ) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

  

pdf (file.path ('Plots', 'TF_abundance_WT1_cHub_profile_barplot.pdf'),6,width=6.5)
bp
dev.off()



















