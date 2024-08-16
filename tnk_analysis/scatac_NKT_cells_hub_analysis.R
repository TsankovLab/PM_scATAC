# Load functions for hub detection ####
source (file.path('..','PM_scATAC','utils','knnGen.R'))
source (file.path('..','PM_scATAC','utils','addCoax.R'))
source (file.path('..','PM_scATAC','utils','Hubs_finder.R'))
source (file.path('..','PM_scATAC','utils','hubs_track.R'))


if (makeBW)
{

# Export bigiwg files ####
archp$celltype_status = paste0(archp$celltype2, '_', archp$status)
metaGroupName = 'celltype2'
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
}


if (!file.exists ('peak_regions.bed'))
  {
  peak_regions = as.data.frame (getPeakSet (archp), row.names=NULL)
  peak_regions = peak_regions[,c(1:3)]
  write.table (peak_regions, file.path('peak_regions.bed'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  }


### Hubs analysis #####
metaGroupName = "Clusters_H"
cor_cutoff = 0.3
#max_dist = 12500
max_dist = 12500
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
dir.create(file.path (hubs_dir, 'Plots'), recursive=T)


# Generate cluster-aware knn groups ####
k= 30
metaGroupName = 'Clusters_H'

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


# Run Co-accessibility ####
run_coax = FALSE
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
metaGroupName = 'celltype3'
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
hm = Heatmap (
  scale (t(hubsSample_mat)), 
#  top_annotation = ha, 
  column_names_gp = gpar(fontsize = 0),
  show_column_dend = F,
  #row_dend_width = unit(5,'mm'),
  row_dend_side = 'left',
  col = rev(palette_hubs_accessibility),
  border=T,
  name = 'Hubs')
pdf (file.path (hubs_dir,'Plots',paste0('hubs_',metaGroupName,'_heatmap.pdf')), height=2.2)
hm
dev.off()

# Annotate also exhausted CD8 cells and generate heatmap of hubs x metagroup ####
cd8_ct = read.csv (file.path('..','..','CD8','scatac_ArchR','barcode_annotation.csv')) # get annotation from subclustered CD8 cells
archp$celltype3_ext = archp$celltype3
archp$celltype3_ext[match(cd8_ct$barcode, rownames(archp@cellColData))] = cd8_ct$celltype
archp$celltype3_ext[archp$celltype3_ext %in% c('C1','C2','C4')] = 'CD8'
metaGroupName = 'celltype3_ext'

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
hm = Heatmap (
  scale (t(hubsSample_mat)), 
#  top_annotation = ha, 
  column_names_gp = gpar(fontsize = 0),
  show_column_dend = F,
  #row_dend_width = unit(5,'mm'),
  row_dend_side = 'left',
  col = rev(palette_hubs_accessibility),
  border=T,
  name = 'Hubs')
pdf (file.path (hubs_dir,'Plots',paste0('hubs_',metaGroupName,'_heatmap.pdf')), height=2.2)
hm
dev.off()

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



# Compute differential hub accessibility to identify DAH in Tregs ####
library (presto)
metaGroupName = 'celltype3'
all (colnames(hubsCell_mat) == rownames(archp@cellColData))
res = wilcoxauc (log2(hubsCell_mat+1), as.character (archp@cellColData[,metaGroupName]))

res_l = lapply (split (res, res$group), function(x){
  tmp = x[x$logFC > 0,]
  tmp = tmp[order (tmp$pval),]
  tmp
})

head (res_l[['Tregs']],100)

## Hubs to show in IGV ####
HUB84 HUB499 HUB1324 HUB575 HUB178 HUB733 HUB429 HUB369 HUB242 HUB602



### Call peaks on celltypes ####
metaGroupName = 'celltype3'
force=FALSE
if(!all(file.exists(file.path('PeakCalls', paste0(unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds')))) | force) 
source (file.path('..','PM_scATAC','callPeaks.R'))
  
### TF Enrichment in hubs in each cell type ####
tf_enr_l = list()
for (ct in unique(archp@cellColData[,metaGroupName]))
{
tf_match = getMatches (archp)
bg_peaks = readRDS (file.path('PeakCalls',paste0(ct,'-reproduciblePeaks.gr.rds')))
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)[queryHits(findOverlaps(tf_match,bg_peaks))]
hubs_regions = hubs_obj$hubsCollapsed[which(hubs_obj$hubs_id %in% head(res_l[[ct]]$feature,100))]
hubs_peaks = bg_peakSet[queryHits(findOverlaps(bg_peakSet, hubs_regions))]
#tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
hubs_regions_TF =  hyperMotif (
  selected_peaks = hubs_peaks, 
  motifmatch = tf_match)

tf_enr_l[[ct]] = hubs_regions_TF
}
names (tf_enr_l)




tf_markers = c('FOXP3','MAFF','JDP2','FOSB','FOS','BACH1','NFEL2L2','NFE2')
markerMotifs = getFeatures (archp, select = paste(tf_markers, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs = grep ("z:", markerMotifs, value = TRUE)
#archp = addImputeWeights (archp)
TF_p = plotEmbedding(
    ArchRProj = archp, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP_H",
    pal = rev (palette_deviation),
    imputeWeights = getImputeWeights(archp)
)

pdf (file.path ('Plots','top_TF_markers_Tregs_meso_only.pdf'), width = 30, height=18)
wrap_plots (TF_p, ncol=4)
dev.off()


# Compute differential hub accessibility between NK cells subsets ####
library (presto)
metaGroupName = 'celltype3'
all (colnames(hubsCell_mat) == rownames(archp@cellColData))

hubsCell_mat_nk = hubsCell_mat[,as.character (archp@cellColData[,metaGroupName]) %in% c('NK_FGFBP2','NK_KLRC1')]
nk_groups = as.character (archp@cellColData[,metaGroupName])[as.character (archp@cellColData[,metaGroupName]) %in% c('NK_FGFBP2','NK_KLRC1')]
res = wilcoxauc (log2(hubsCell_mat_nk+1), nk_groups)

res_l = lapply (split (res, res$group), function(x){
  tmp = x[x$logFC > 0,]
  tmp = tmp[order (tmp$pval),]
  tmp
})

head (res_l[['NK_KLRC1']],10)



### Order cells using trajectory from CD4 to Tregs and find hubs that correlate with it ####
archp = addImputeWeights (archp)
seGS <- getMatrixFromProject (archp)
#celltype_markers = c('FOXP3',rowData (seGS)$name[grep ('FOXP3', rowData (seGS)$name)])
celltype_markers = c('FOXP3','ILRA2','CTLA4')
seGS = seGS[rowData (seGS)$name %in% celltype_markers,]
matGS <- imputeMatrix (assay(seGS), getImputeWeights(archp))
rownames(matGS) = rowData (seGS)$name
#rownames (matGS) = celltype_markers


archp = addTrajectory(
    ArchRProj = archp, 
    name = "Treg_maturation", 
    groupBy = "Clusters_H",
    trajectory = c('C7','C3','C2'), 
    embedding = "UMAP_H", 
    force = FALSE
)

p = plotTrajectory (archp, trajectory = "Treg_maturation", colorBy = "cellColData", name = "Treg_maturation", embedding = 'UMAP_H')
pdf (file.path ('Plots','Treg_maturation_trajectory.pdf'))
p
dev.off()

trajMM  <- getTrajectory(ArchRProj = archp, name = "Treg_maturation", useMatrix = "MotifMatrix", log2Norm = FALSE)
trajGS  <- getTrajectory(ArchRProj = archp, name = "Treg_maturation", useMatrix = "GeneScoreMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = palette_deviation)
p2 <- plotTrajectoryHeatmap(trajGS, pal = palette_expression)

pdf (file.path ('Plots','Treg_TF_trajectory.pdf'))
p1
p2
dev.off()

# Get pseudotime scores and bin them ####
pseudotime_scores = archp$Treg_maturation
names (pseudotime_scores) = rownames(archp@cellColData)
pseudotime_scores = na.omit (pseudotime_scores)
pseudotime_scores = pseudotime_scores[order(pseudotime_scores)]

pseudotime_scores_binned = paste0('PT',ceiling(seq_along(pseudotime_scores)/400))
names (pseudotime_scores_binned) = names (pseudotime_scores)

archp$pseudotime_binned = pseudotime_scores_binned[match(rownames(archp@cellColData), names(pseudotime_scores))]
archp2 = archp[!is.na(archp$pseudotime_binned)]

# Generate matrix of fragment counts using pseudotime bins ####
metaGroupName = 'pseudotime_binned'
force = TRUE
if (!file.exists(file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds'))) | force)
  {
  if (!exists ('fragments')) fragments = unlist (getFragmentsFromProject (
    ArchRProj = archp))   
  hubsSample_mat = matrix (ncol = length(unique(archp2@cellColData[,metaGroupName])), nrow = length(hubs_obj$hubsCollapsed))
  colnames (hubsSample_mat) = unique(archp2@cellColData[,metaGroupName])
  rownames (hubsSample_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (unique(archp2@cellColData[,metaGroupName])))
  for (sam in unique(archp2@cellColData[,metaGroupName]))
    {
    pb$tick()  
    fragments_in_sample = fragments[fragments$RG %in% rownames(archp2@cellColData)[as.character(archp2@cellColData[,metaGroupName]) == sam]]  
    fragments_in_sample_in_hubs = countOverlaps (hubs_obj$hubsCollapsed, fragments_in_sample)
    hubsSample_mat[,sam] = fragments_in_sample_in_hubs
    }
  frags_in_sample = sapply (unique(archp2@cellColData[,metaGroupName]), function(x) sum (archp2$nFrags[as.character(archp2@cellColData[,metaGroupName]) == x]))
  hubsSample_mat = t(t(hubsSample_mat) * (10^6 / frags_in_sample)) # scale
  saveRDS (hubsSample_mat, file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))
  } else {
  hubsSample_mat = readRDS (file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))  
  }
hubsSample_mat = as.data.frame (hubsSample_mat)


# Correlate hubs vs binned pseudotime score ####
pseudotime_binned_avg = unlist(lapply (split (archp2$Treg_maturation, archp2$pseudotime_binned), mean))

hubsSample_mat_cor = as.data.frame (t(cor (pseudotime_binned_avg, t(hubsSample_mat[,names(pseudotime_binned_avg)]))), method='spearman')
top_cor_hubs = rownames(head (hubsSample_mat_cor[order(-hubsSample_mat_cor[,1]), ,drop=F],100))


hubsSample_mat_top = hubsSample_mat[top_cor_hubs, names(pseudotime_binned_avg[order (pseudotime_binned_avg)])]

top_cor_hubs_labels = head (paste0(hubs_obj$hubs_id,':',hubs_obj$hubsCollapsed$gene)[match (top_cor_hubs, hubs_obj$hubs_id)],10)



ha = HeatmapAnnotation (psuedotime = pseudotime_binned_avg[order (pseudotime_binned_avg)])
ha2 = rowAnnotation(foo = anno_mark(at = seq_along(top_cor_hubs_labels), 
    labels = top_cor_hubs_labels, labels_gp = gpar(fontsize = 6)))
hm = Heatmap (
  t(scale (t(hubsSample_mat_top))), 
  top_annotation = ha, 
  right_annotation = ha2,
  column_names_gp = gpar(fontsize = 0),
  row_names_gp = gpar(fontsize = 0),
  show_column_dend = F,
  #row_dend_width = unit(5,'mm'),
  row_dend_side = 'left',
  cluster_rows = F,
  cluster_columns = F,
  col = palette_pseudotime,
  border=T,
  name = 'pseudotime_Hubs')
pdf (file.path (hubs_dir,'Plots',paste0('hubs_',metaGroupName,'_pseudotime_heatmap.pdf')), width=4.2, height=3.2)
hm
dev.off()

# Export bigiwg files ####
metaGroupName = 'pseudotime_binned'
exp_bigwig = T
if (exp_bigwig)
  {
  getGroupBW(
    ArchRProj = archp2,
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
}

# Select HUBs to show on IGV along CD4 > Treg pseudotime ####
HUB135 HUB219 HUB531 HUB396 HUB1198 HUB620 HUB532 HUB90 HUB1223 HUB257 


### TF Enrichment in hubs in each cell type ####
metaGroupName = 'celltype3'
tf_enr_l = list()
# for (ct in unique(archp@cellColData[,metaGroupName]))
# {
ct = 'Tregs'
tf_match = getMatches (archp)
bg_peaks = readRDS (file.path('PeakCalls',paste0(ct,'-reproduciblePeaks.gr.rds')))
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)[queryHits(findOverlaps(tf_match,bg_peaks))]
mega_hubs = hubs_obj$hubsCollapsed[which(hubs_obj$hubs_id %in% top_cor_hubs[2])]
mega_hubs_peaks = bg_peakSet[queryHits(findOverlaps(bg_peakSet, mega_hubs))]
#tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
mega_hubs_TF =  hyperMotif (
  selected_peaks = mega_hubs_peaks, 
  motifmatch = tf_match)

tf_enr_l[[ct]] = mega_hubs_TF
#}
ct = 'Tregs'
head (tf_enr_l[[ct]],20)






# Compute differential hub accessibility between CD8 exhausted and rest of CD8 ####
library (presto)
cd8_ct = read.csv (file.path('..','..','CD8','scatac_ArchR','barcode_annotation.csv')) # get annotation from subclustered CD8 cells
archp$celltype3_ext = archp$celltype3
archp$celltype3_ext[match(cd8_ct$barcode, rownames(archp@cellColData))] = cd8_ct$celltype
archp$celltype3_ext[archp$celltype3_ext %in% c('C1','C2','C4')] = 'CD8'
metaGroupName = 'celltype3_ext'
all (colnames(hubsCell_mat) == rownames(archp@cellColData))

hubsCell_mat_cd8 = hubsCell_mat[,as.character (archp@cellColData[,metaGroupName]) %in% c('CD8','CD8_exhausted')]
cd8_groups = as.character (archp@cellColData[,metaGroupName])[as.character (archp@cellColData[,metaGroupName]) %in% c('CD8','CD8_exhausted')]
res = wilcoxauc (log2(hubsCell_mat_cd8+1), cd8_groups)

res_l = lapply (split (res, res$group), function(x){
  tmp = x[x$logFC > 0,]
  tmp = tmp[order (tmp$pval),]
  tmp
})

head (res_l[['CD8_exhausted']],10)

### TF Enrichment in hubs in each cell type ####
tf_enr_l = list()

tf_match = getMatches (archp)
bg_peaks = readRDS (file.path('PeakCalls',paste0(ct,'-reproduciblePeaks.gr.rds')))
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)[queryHits(findOverlaps(tf_match,bg_peaks))]
mega_hubs = hubs_obj$hubsCollapsed[which(hubs_obj$hubs_id %in% head(res_l[[ct]]$feature,100))]
mega_hubs_peaks = bg_peakSet[queryHits(findOverlaps(bg_peakSet, mega_hubs))]
#tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
mega_hubs_TF =  hyperMotif (
  selected_peaks = mega_hubs_peaks, 
  motifmatch = tf_match)

tf_enr_l[[ct]] = mega_hubs_TF
}
names (tf_enr_l)


# Check DAH between NK KLRC1+ and CD8 exhausted to find commonalities ####
cd8_ct = read.csv (file.path('..','..','CD8','scatac_ArchR','barcode_annotation.csv')) # get annotation from subclustered CD8 cells
archp$celltype3_ext = archp$celltype3
archp$celltype3_ext[match(cd8_ct$barcode, rownames(archp@cellColData))] = cd8_ct$celltype
archp$celltype3_ext[archp$celltype3_ext %in% c('C1','C2','C4')] = 'CD8'
