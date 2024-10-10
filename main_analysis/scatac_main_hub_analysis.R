# Load functions for hub detection
source (file.path('..','..','git_repo','utils','knnGen.R'))
source (file.path('..','..','git_repo','utils','addCoax.R'))
source (file.path('..','..','git_repo','utils','Hubs_finder.R'))
source (file.path('..','..','git_repo','utils','hubs_track.R'))


# Export bigiwg files ####
#archp$celltype_status = paste0(archp$celltype2, '_', archp$status)
metaGroupName = 'celltype_revised'
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


if (!file.exists ('peak_regions.bed'))
  {
  peak_regions = as.data.frame (getPeakSet (archp), row.names=NULL)
  peak_regions = peak_regions[,c(1:3)]
  write.table (peak_regions, file.path('peak_regions.bed'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  }


### Hubs analysis #####
metaGroupName = "celltype_revised"
cor_cutoff = 0.3
#max_dist = 12500
max_dist = 12500
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
dir.create(file.path (hubs_dir, 'Plots'), recursive=T)


# Generate cluster-aware knn groups ####
k= 30
metaGroupName = 'celltype_revised'

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
metaGroupName = 'celltype_revised'
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
pdf (file.path (hubs_dir,'Plots',paste0('hubs_',metaGroupName,'_heatmap.pdf')), height=3)
hm
dev.off()


# Generate matrix of fragment counts of hubs x barcodes ####
metaGroupName = 'celltype_revised'
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
metaGroupName = 'celltype_revised'
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
table (archp$celltype_revised)
archp$compartment = archp$celltype_revised
archp$compartment[archp$celltype_revised %in% c('Alveolar','Mesothelium','Malignant')] = 'epithelial'
archp$compartment[archp$celltype_revised %in% c('Endothelial','Fibroblasts','SmoothMuscle')] = 'stroma'
archp$compartment[archp$celltype_revised %in% c('Myeloid','pDCs')] = 'myeloid'
archp$compartment[archp$celltype_revised %in% c('B_cells','NK','Plasma','T_cells')] = 'lymphoid'

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

vp = VlnPlot (srt, feature = c('ERG','VIM'), group.by = 'celltype_simplified', col=palette_celltype_simplified) + NoLegend()

hubsCell_mat = t(hubsCell_mat)[rownames(archp@cellColData),]
archp@cellColData = cbind(archp@cellColData, log2(hubsCell_mat[,c('HUB322'), drop=F]+1))
#archp@cellColData = archp@cellColData[, -ncol(archp@cellColData)]
p <- plotGroups(
  ArchRProj = archp, 
  groupBy = "celltype_revised", 
  colorBy = "cellColData", 
  name = "HUB322",
  plotAs = "violin",
  pal = palette_celltype_simplified,
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
hub_tf_cor = lapply (unique(archp@cellColData[,'celltype_revised']), function(x)
  as.data.frame (t(cor (hub_nz[match(rownames(archp_meta)[as.character(archp_meta$celltype_revised) == x],
    names(hub_nz))], mMat_nz[match(rownames(archp_meta)[as.character(archp_meta$celltype_revised) == x],
    rownames(mMat_nz)),], method='spearman')))
  )
hub_tf_cor = lapply (hub_tf_cor, function(x) x[order(-x[,1]),, drop=F])
hub_tf_cor = as.data.frame (t(hub_tf_cor))
head (hub_tf_cor[order (-hub_tf_cor[,1]),,drop=F],10)


cor_df = data.frame (hub = hub_exp[which(zero_val)], tf = t(as.matrix(mMat))[which(zero_val),'ETS1'])
cor_df$celltype = archp$celltype_revised[match(rownames(cor_df), rownames(archp@cellColData))]
cor_df$sample = archp$Sample[match(rownames(cor_df), rownames(archp@cellColData))]

cor (cor_df[cor_df$celltype == 'Malignant','hub'], cor_df[cor_df$celltype == 'Malignant','tf'])

cor_p = ggplot (cor_df[cor_df$celltype == 'T_cells',], aes (x = hub, y = tf, color = celltype)) + geom_point(alpha=.4)
cor_p2 = ggplot (cor_df[cor_df$celltype == 'T_cells',], aes (x = hub, y = tf, color = sample)) + geom_point(alpha=.4)

pdf (file.path('Plots','HUB322_VIM_correlation_ETS1.pdf'), width=9)
wrap_plots (cor_p,cor_p2)
dev.off()



















### TF Enrichment in hubs ####
tf_match = getMatches (archp)
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)
mega_hubs = range (hubs_obj$hubsCollapsed[which (hubs_obj$hubs_id %in% c('HUB322'))])
mega_hubs_peaks = bg_peakSet[queryHits(findOverlaps(bg_peakSet, mega_hubs))]
#tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
mega_hubs_TF =  hyperMotif (
  selected_peaks = mega_hubs_peaks, 
  motifmatch = tf_match)

head (mega_hubs_TF, 20)





# Compute differential hub accessibility DHA to check hub size normal vs tumor cells ####
library (presto)
metaGroupName = 'celltype_revised'
all (colnames(hubsCell_mat) == rownames(archp@cellColData))
res = wilcoxauc (log2(hubsCell_mat+1), as.character (archp@cellColData[,metaGroupName]))


res_l = lapply (split (res, res$group), function(x){
  tmp = x[order (x$padj),]
  tmp
})

# Generate boxplot 
res2 = res[res$group == 'Malignant',]
res2$group2 = as.character (sign (res2$logFC))
res2$width = width (hubs_obj$hubsCollapsed)
rp = ggplot (res2, aes(x = group2, y= width)) + geom_boxplot() + gtheme

pdf (file.path (hubs_dir,'Plots','difference_width_tumor_normal_boxplot.pdf'))
rp
dev.off()





# Look into mega hubs identified in P11 HOX- cluster ####
# Map large hubs to UMAP ####
hubsCell_mat = t(hubsCell_mat)[rownames(archp@cellColData),]
archp@cellColData = cbind(archp@cellColData, hubsCell_mat[,c('HUB1','HUB2','HUB3','HUB4','HUB5','HUB6','HUB7','HUB8','HUB9','HUB10')])
umap_p1 =  plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
 name = c('HUB1','HUB2','HUB3','HUB4','HUB5','HUB6','HUB7','HUB8','HUB9','HUB10'), embedding = "UMAP")
  
pdf (file.path(hubs_dir,'Plots','megahubs_umap.pdf'), 15,15)
wrap_plots (umap_p1, ncol=4)
dev.off()













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
    ArchRProj = archp, 
    groupBy = metaGroupName2, 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p1 <- plotGroups(
    ArchRProj = archp, 
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


# Identify hubs with no annotated genes within 1MB ####
hubs = hubs_obj$hubsCollapsed
hubs_ext = extendGR (hubs, upstream = 1000000, downstream = 1000000)
geneless_hubs = findOverlaps (hubs_ext, archp@geneAnnotation[[1]])
hubs_ext[!1:length(hubs) %in% unique(queryHits(geneless_hubs))]









