# Load functions for hub detection
source ('../PM_scATAC/knnGen.R')
source ('../PM_scATAC/addCoax.R')
source ('../PM_scATAC/Hubs_finder.R')
source ('../PM_scATAC/hubs_track.R')


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
    ArchRProj = archp_sub))   
  hubsSample_mat = matrix (ncol = length(unique(archp_sub@cellColData[,metaGroupName])), nrow = length(hubs_obj$hubsCollapsed))
  colnames (hubsSample_mat) = unique(archp_sub@cellColData[,metaGroupName])
  rownames (hubsSample_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (unique(archp_sub@cellColData[,metaGroupName])))
  for (sam in unique(archp_sub@cellColData[,metaGroupName]))
    {
    pb$tick()  
    fragments_in_sample = fragments[fragments$RG %in% rownames(archp_sub@cellColData)[as.character(archp_sub@cellColData[,metaGroupName]) == sam]]  
    fragments_in_sample_in_hubs = countOverlaps (hubs_obj$hubsCollapsed, fragments_in_sample)
    hubsSample_mat[,sam] = fragments_in_sample_in_hubs
    }
  frags_in_sample = sapply (unique(archp_sub@cellColData[,metaGroupName]), function(x) sum (archp_sub$nFrags[as.character(archp_sub@cellColData[,metaGroupName]) == x]))
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
    ArchRProj = archp_sub))    
  hubsCell_mat = matrix (ncol = length(rownames(archp_sub@cellColData)), nrow = length(hubs_obj$hubsCollapsed))
  colnames (hubsCell_mat) = rownames(archp_sub@cellColData)
  rownames (hubsCell_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (rownames(archp_sub@cellColData)))
  for (cell in rownames(archp_sub@cellColData)) 
    {
    pb$tick()  
    fragments_in_cell = fragments[fragments$RG %in% cell]  
    fragments_in_cell_in_hubs = countOverlaps (hubs_obj$hubsCollapsed, fragments_in_cell)
    hubsCell_mat[,cell] = fragments_in_cell_in_hubs
    }
  all (colnames (hubsCell_mat) == rownames(archp_sub@cellColData))  
  hubsCell_mat = t(t(hubsCell_mat) * (10^6 / archp_sub$nFrags)) # scale
  saveRDS (hubsCell_mat, file.path (hubs_dir,paste0('hubs_cells_mat.rds')))
  } else {
  hubsCell_mat = readRDS (file.path (hubs_dir,paste0('hubs_cells_mat.rds')))  
  }
hubsCell_mat = as.data.frame (hubsCell_mat)



# Compute differential hub accessibility DHA ####
library (presto)
metaGroupName = 'celltype3'
all (colnames(hubsCell_mat) == rownames(archp_sub@cellColData))
res = wilcoxauc (log2(hubsCell_mat+1), as.character (archp_sub@cellColData[,metaGroupName]))

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
archp_sub = addGroupCoverages (
  ArchRProj = archp_sub, 
  groupBy = metaGroupName,  
  force = FALSE,
  minCells= 20, # I think this should be set corresponding to the smallest cluster in the group or lower
  maxCells = 500,
  minReplicates = 2,
  sampleRatio = 0.8,
  useLabels = TRUE)

archp_sub = addReproduciblePeakSet (
    archp_sub,
    groupBy= metaGroupName,
    peakMethod = 'Macs2',
    reproducibility = "1",
    maxPeaks = 500000, 
    minCells=20,
    force =TRUE) # I think this should be set corresponding to the smallest cluster in the group or lower
archp_sub = addPeakMatrix (archp_sub)
  
### TF Enrichment in hubs in each cell type ####
tf_enr_l = list()
for (ct in unique(archp_sub@cellColData[,metaGroupName]))
{
tf_match = getMatches (archp_sub)
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


# Check hub size distribution between normal and tumor ####
#res$width = width(hubs_obj$hubsCollapsed)[match (res$feature, hubs_obj$hubs_id)]
vp = ggplot (size_comp_df, aes(x=group, y=logFC)) +
    geom_boxplot() + 
    facet_wrap (~ celltype) +
    ggtitle ('Hub size') + gtheme

pdf (file.path (hubs_dir, 'Plots','DAH_size.pdf'),4,4)
vp
dev.off()






tf_markers = c('FOXP3','MAFF','JDP2','FOSB','FOS','BACH1','NFEL2L2','NFE2')
markerMotifs = getFeatures (archp_sub, select = paste(tf_markers, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs = grep ("z:", markerMotifs, value = TRUE)
#archp = addImputeWeights (archp)
TF_p = plotEmbedding(
    ArchRProj = archp_sub, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP_H",
    pal = rev (palette_deviation),
    imputeWeights = getImputeWeights(archp_sub)
)

pdf (file.path ('Plots','top_TF_markers_Tregs_meso_only.pdf'), width = 30, height=18)
wrap_plots (TF_p, ncol=4)
dev.off()




### Order cells by FOXP3 expression and find hubs that correlate with it ####
archp_sub = addImputeWeights (archp_sub)
seGS <- getMatrixFromProject (archp_sub)
#celltype_markers = c('FOXP3',rowData (seGS)$name[grep ('FOXP3', rowData (seGS)$name)])
celltype_markers = c('FOXP3','ILRA2','CTLA4')
seGS = seGS[rowData (seGS)$name %in% celltype_markers,]
matGS <- imputeMatrix (assay(seGS), getImputeWeights(archp))
rownames(matGS) = rowData (seGS)$name
#rownames (matGS) = celltype_markers





archp_sub = addTrajectory(
    ArchRProj = archp_sub, 
    name = "Treg_maturation", 
    groupBy = "Clusters_H",
    trajectory = c('C7','C3','C2'), 
    embedding = "UMAP_H", 
    force = FALSE
)

p = plotTrajectory(archp_sub, trajectory = "Treg_maturation", colorBy = "cellColData", name = "Treg_maturation", embedding = 'UMAP_H')
pdf (file.path ('Plots','Treg_maturation_trajectory.pdf'))
p
dev.off()

trajMM  <- getTrajectory(ArchRProj = archp_sub, name = "Treg_maturation", useMatrix = "MotifMatrix", log2Norm = FALSE)
trajGS  <- getTrajectory(ArchRProj = archp_sub, name = "Treg_maturation", useMatrix = "GeneScoreMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = palette_deviation)
p2 <- plotTrajectoryHeatmap(trajGS, pal = palette_expression)

pdf (file.path ('Plots','Treg_TF_trajectory.pdf'))
p1
p2
dev.off()

# Get pseudotime scores and bin them ####
pseudotime_scores = archp_sub$Treg_maturation
names (pseudotime_scores) = rownames(archp_sub@cellColData)
pseudotime_scores = na.omit (pseudotime_scores)
pseudotime_scores = pseudotime_scores[order(pseudotime_scores)]

pseudotime_scores_binned = paste0('PT',ceiling(seq_along(pseudotime_scores)/400))
names (pseudotime_scores_binned) = names (pseudotime_scores)

archp_sub$pseudotime_binned = pseudotime_scores_binned[match(rownames(archp_sub@cellColData), names(pseudotime_scores))]
archp_sub2 = archp_sub[!is.na(archp_sub$pseudotime_binned)]

# Generate matrix of fragment counts using pseudotime bins ####
metaGroupName = 'pseudotime_binned'
force = TRUE
if (!file.exists(file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds'))) | force)
  {
  if (!exists ('fragments')) fragments = unlist (getFragmentsFromProject (
    ArchRProj = archp_sub))   
  hubsSample_mat = matrix (ncol = length(unique(archp_sub2@cellColData[,metaGroupName])), nrow = length(hubs_obj$hubsCollapsed))
  colnames (hubsSample_mat) = unique(archp_sub2@cellColData[,metaGroupName])
  rownames (hubsSample_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (unique(archp_sub2@cellColData[,metaGroupName])))
  for (sam in unique(archp_sub2@cellColData[,metaGroupName]))
    {
    pb$tick()  
    fragments_in_sample = fragments[fragments$RG %in% rownames(archp_sub2@cellColData)[as.character(archp_sub2@cellColData[,metaGroupName]) == sam]]  
    fragments_in_sample_in_hubs = countOverlaps (hubs_obj$hubsCollapsed, fragments_in_sample)
    hubsSample_mat[,sam] = fragments_in_sample_in_hubs
    }
  frags_in_sample = sapply (unique(archp_sub2@cellColData[,metaGroupName]), function(x) sum (archp_sub2$nFrags[as.character(archp_sub2@cellColData[,metaGroupName]) == x]))
  hubsSample_mat = t(t(hubsSample_mat) * (10^6 / frags_in_sample)) # scale
  saveRDS (hubsSample_mat, file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))
  } else {
  hubsSample_mat = readRDS (file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))  
  }
hubsSample_mat = as.data.frame (hubsSample_mat)


# Correlate hubs vs binned pseudotime score ####
pseudotime_binned_avg = unlist(lapply (split (archp_sub2$Treg_maturation, archp_sub2$pseudotime_binned), mean))

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
    ArchRProj = archp_sub2,
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

# Select HUBs to show on IGV ####
HUB135 HUB219 HUB531 HUB396 HUB1198 HUB620 HUB532 HUB90 HUB1223 HUB257 


### TF Enrichment in hubs in each cell type ####
metaGroupName = 'celltype3'
tf_enr_l = list()
# for (ct in unique(archp_sub@cellColData[,metaGroupName]))
# {
ct = 'Tregs'
tf_match = getMatches (archp_sub)
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









# differential hub t-test on sample normal vs tumor ####
# Generate matrix of fragment counts of hubs x metagroup ####
archp_sub$celltype_status_sample = paste0(archp_sub$celltype_status,'_',archp_sub$Sample2)
keep_samples = names(table (archp_sub$celltype_status_sample)[table (archp_sub$celltype_status_sample) > 10])

metaGroupName = 'celltype_status_sample'
if (!file.exists(file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds'))))
  {
  if (!exists ('fragments')) fragments = unlist (getFragmentsFromProject (
    ArchRProj = archp_sub))   
  hubsSample_mat = matrix (ncol = length(unique(archp_sub@cellColData[,metaGroupName])), nrow = length(hubs_obj$hubsCollapsed))
  colnames (hubsSample_mat) = unique(archp_sub@cellColData[,metaGroupName])
  rownames (hubsSample_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (unique(archp_sub@cellColData[,metaGroupName])))
  for (sam in unique(archp_sub@cellColData[,metaGroupName]))
    {
    pb$tick()  
    fragments_in_sample = fragments[fragments$RG %in% rownames(archp_sub@cellColData)[as.character(archp_sub@cellColData[,metaGroupName]) == sam]]  
    fragments_in_sample_in_hubs = countOverlaps (hubs_obj$hubsCollapsed, fragments_in_sample)
    hubsSample_mat[,sam] = fragments_in_sample_in_hubs
    }
  frags_in_sample = sapply (unique(archp_sub@cellColData[,metaGroupName]), function(x) sum (archp_sub$nFrags[as.character(archp_sub@cellColData[,metaGroupName]) == x]))
  hubsSample_mat = t(t(hubsSample_mat) * (10^6 / frags_in_sample)) # scale
  saveRDS (hubsSample_mat, file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))
  } else {
  hubsSample_mat = readRDS (file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))  
  }
hubsSample_mat = as.data.frame (hubsSample_mat)

hubsSample_mat = hubsSample_mat[, keep_samples]

meta_group = c('NK','T_cells')
res_l = list()
for (i in meta_group)
  {
  hubsSample_mat_df = log2(hubsSample_mat[, grepl (i, colnames(hubsSample_mat))]+1)
  res = sapply (seq(nrow(hubsSample_mat_df)), function(x) t.test (
    hubsSample_mat_df[x,grepl ('tumor',colnames(hubsSample_mat_df))],
    hubsSample_mat_df[x,grepl ('normal',colnames(hubsSample_mat_df))])$p.value)
  lfc = sapply (seq(nrow(hubsSample_mat_df)), function(x) rowMeans(hubsSample_mat_df[x,grepl ('tumor',colnames(hubsSample_mat_df))]) - rowMeans(hubsSample_mat_df[x,grepl ('normal',colnames(hubsSample_mat_df))]))
  res_l[[i]] = data.frame(
    hub = rownames(hubsSample_mat_df), 
    pval = res,
    padj = p.adjust (res), 
    logFC = lfc, 
    celltype = i)
  }
head (res_l[[2]][order(res_l[[2]]$padj),],30)

## Plot volcano plots ####
res = do.call (rbind, res_l)
res = res[order(res$padj), ]
logfcThreshold = 1
pvalAdjTrheshold = 0.05
res$sig = ifelse (abs(res$logFC) > logfcThreshold & res$padj < pvalAdjTrheshold, 1,0)
res$sig = res$sig * sign (res$logFC)
res$sig = as.character (res$sig)
res_filtered = res[abs(res$logFC) > logfcThreshold & res$padj < pvalAdjTrheshold,]
res_filtered = head (res_filtered$hub[order (-abs(res_filtered$logFC))],50)
res$labels = ''
res$labels[match (res_filtered, res$hub)] = res_filtered
vp = ggplot (res, aes(x=logFC, y=-log10(padj))) +
     geom_point(size=1, shape=19, aes (color = sig), alpha=.5) +
     geom_vline(xintercept = logfcThreshold, linetype="dotted", 
                 color = "grey20", size=1) +
     geom_vline(xintercept = -logfcThreshold, linetype="dotted", 
                 color = "grey20", size=1) +
     geom_hline(yintercept = -log10(pvalAdjTrheshold), linetype="dotted", 
                 color = "grey20", size=1) + 
     geom_text_repel (size=2, data = res, aes(label = labels)) + 
     ggtitle ('Hubs differential accessibility') +
     #geom_label_repel (size=2,max.overlaps=10000, data = deg2_cl, aes(label = show_genes), color='red') + 
     scale_color_manual (values = c("0"='grey77',"-1"='#666666FF',"1"='#F8A02EFF')) + theme_light() +
     facet_wrap (.~celltype, ncol=4)

pdf (file.path (hubs_dir, 'Plots', 'DAH_volcano_ttest.pdf'),width = 6,height = 3)
print (vp)
dev.off()  


# Generate matrix of fragment counts of hubs x barcodes ####
metaGroupName = 'celltype_status'
if (!file.exists(file.path (hubs_dir, paste0('hubs_cells_',metaGroupName,'_mat.rds'))))
  {
  fragments = unlist (getFragmentsFromProject (
    ArchRProj = archp_sub))    
  hubsCell_mat = matrix (ncol = length(rownames(archp_sub@cellColData)), nrow = length(hubs_obj$hubsCollapsed))
  colnames (hubsCell_mat) = rownames(archp_sub@cellColData)
  rownames (hubsCell_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (rownames(archp_sub@cellColData)))
  for (cell in rownames(archp_sub@cellColData)) 
    {
    pb$tick()  
    fragments_in_cell = fragments[fragments$RG %in% cell]  
    fragments_in_cell_in_hubs = countOverlaps (hubs_obj$hubsCollapsed, fragments_in_cell)
    hubsCell_mat[,cell] = fragments_in_cell_in_hubs
    }
  all (colnames (hubsCell_mat) == rownames(archp_sub@cellColData))  
  hubsCell_mat = t(t(hubsCell_mat) * (10^6 / archp_sub$nFrags)) # scale
  saveRDS (hubsCell_mat, file.path (hubs_dir,paste0('hubs_cells_',metaGroupName,'_mat.rds')))
  } else {
  hubsCell_mat = readRDS (file.path (hubs_dir,paste0('hubs_cells_',metaGroupName,'_mat.rds')))  
  }
hubsCell_mat = as.data.frame (hubsCell_mat)




## Plot volcano plots ####
res = size_comp_df
logfcThreshold = 1
pvalAdjTrheshold = 0.05
res$sig = ifelse (abs(res$logFC) > logfcThreshold & res$padj < pvalAdjTrheshold, 1,0)
res$sig = res$sig * sign (res$logFC)
res$sig = as.character (res$sig)
res_filtered = res[abs(res$logFC) > logfcThreshold & res$padj < pvalAdjTrheshold,]
res_filtered = head (res_filtered$feature[order (-abs(res_filtered$logFC))],50)
res$labels = ''
res$labels[match (res_filtered, res$feature)] = res_filtered
vp = ggplot (res, aes(x=logFC, y=-log10(padj))) +
     geom_point(size=1, shape=19, aes (color = sig), alpha=.5) +
     geom_vline(xintercept = logfcThreshold, linetype="dotted", 
                 color = "grey20", size=1) +
     geom_vline(xintercept = -logfcThreshold, linetype="dotted", 
                 color = "grey20", size=1) +
     geom_hline(yintercept = -log10(pvalAdjTrheshold), linetype="dotted", 
                 color = "grey20", size=1) + 
     geom_text_repel (size=2, data = res, aes(label = labels)) + 
     ggtitle ('Hubs differential accessibility') +
     #geom_label_repel (size=2,max.overlaps=10000, data = deg2_cl, aes(label = show_genes), color='red') + 
     scale_color_manual (values = c("0"='grey77',"-1"='#666666FF',"1"='#F8A02EFF')) + theme_light() +
     facet_wrap (.~celltype, ncol=4)

pdf (file.path (hubs_dir, 'Plots', 'DAH_volcano.pdf'),width = 14,height = 4)
print (vp)
dev.off()  



# check hubs that are up across all normal vs tumor ####
wlc_res_pos = lapply (wlc_res, function(x) x[x$logFC > 0,])
wlc_res_pos = do.call (rbind, wlc_res_pos)
head (table (wlc_res_pos$feature)[order(-table (wlc_res_pos$feature))])

wlc_res_neg = lapply (wlc_res, function(x) x[x$logFC < 0,])
wlc_res_neg = do.call (rbind, wlc_res_neg)
head (table (wlc_res_neg$feature)[order(-table (wlc_res_neg$feature))])

wlc_res_t = wlc_res[['T_cells']]
head (wlc_res_t[order(wlc_res_t$logFC),],20)





p2 <- plotGroups(
    ArchRProj = archp_sub, 
    groupBy = metaGroupName2, 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p1 <- plotGroups(
    ArchRProj = archp_sub, 
    groupBy = metaGroupName2, 
    colorBy = "cellColData", 
    name = "nFrags",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
pdf (file.path ('Plots','QC-Sample-Statistics.pdf'), height=5)
wrap_plots (p1, p2)
dev.off()




