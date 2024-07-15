
# Load functions for hub detection
source ('../PM_scATAC/knnGen.R')
source ('../PM_scATAC/addCoax.R')
source ('../PM_scATAC/Hubs_finder.R')
source ('../PM_scATAC/hubs_track.R')

# Export bigiwg files ####
archp$celltype_status = paste0(archp$celltype2, '_', archp$status)
metaGroupName = 'celltype_status'
exp_bigwig = FALSE
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
metaGroupName = "Clusters_H"
cor_cutoff = 0.2
#max_dist = 12500
max_dist = 5000
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
dir.create(file.path (hubs_dir, 'Plots'), recursive=T)

# Generate cluster-aware knn groups ####
k= 30
metaGroupName = 'Clusters_H'

force = FALSE
if (!file.exists(paste0 ('KNNs_',metaGroupName,'k_',k,'.rds')) | force)
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
  saveRDS (KNNs, paste0 ('KNNs_',metaGroupName,'k_',k,'.rds'))  
  } else {
  KNNs = readRDS (paste0 ('KNNs_',metaGroupName,'k_',k,'.rds'))
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
    cores=1
    ) 
  saveRDS (hubs_obj, file.path(hubs_dir,'global_hubs_obj.rds'))
  hubs_regions = as.data.frame (hubs_obj$hubsCollapsed)
  hubs_regions$width = hubs_obj$hubs_id
  write.table (hubs_regions, file.path(hubs_dir, 'hub_regions.bed'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  } else {
  hubs_obj = readRDS (file.path(hubs_dir,'global_hubs_obj.rds'))  
  }



# Generate matrix of fragment counts of hubs x metagroup ####
metaGroupName = 'celltype_status'
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
  frags_in_sample = sapply (unique(archp@cellColData[,metaGroupName]), function(x) sum (archp$nfrags[as.character(archp@cellColData[,metaGroupName]) == x]))
  hubsSample_mat = t(t(hubsSample_mat) * (10^6 / frags_in_sample)) # scale
  saveRDS (hubsSample_mat, file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))
  } else {
  hubsSample_mat = readRDS (file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))  
  }
hubsSample_mat = as.data.frame (hubsSample_mat)

ha = HeatmapAnnotation (size = anno_barplot(width (hubs_obj$hubsCollapsed), gp = gpar(color = "red"), height =  unit(8, "mm")))
hm = Heatmap (
  scale (t(hubsSample_mat)), 
  top_annotation = ha, 
  column_names_gp = gpar(fontsize = 0),
  show_column_dend = F,
  row_dend_width = unit(3,'mm'),
  row_dend_side = 'left')
pdf (file.path (hubs_dir,'Plots',paste0('hubs_',metaGroupName,'_heatmap.pdf')), height=3)
hm
dev.off()

# differential hub t-test on sample normal vs tumor ####
# Generate matrix of fragment counts of hubs x metagroup ####
archp$celltype_status_sample = paste0(archp$celltype_status,'_',archp$Sample2)
keep_samples = names(table (archp$celltype_status_sample)[table (archp$celltype_status_sample) > 10])

metaGroupName = 'celltype_status_sample'
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
  hubsCell_mat = t(t(hubsCell_mat) * (10^6 / archp$nfrags)) # scale
  saveRDS (hubsCell_mat, file.path (hubs_dir,paste0('hubs_cells_',metaGroupName,'_mat.rds')))
  } else {
  hubsCell_mat = readRDS (file.path (hubs_dir,paste0('hubs_cells_',metaGroupName,'_mat.rds')))  
  }
hubsCell_mat = as.data.frame (hubsCell_mat)

ha = HeatmapAnnotation (size = anno_barplot (width (hubs_obj$hubsCollapsed), gp = gpar(color = "red"), height =  unit(8, "mm")))
ha2 = HeatmapAnnotation (celltype = as.character(archp@cellColData[,metaGroupName]), height =  unit(8, "mm"), which = 'row')
hm = Heatmap (
  scale (t(hubsCell_mat)), 
  left_annotation = ha2,
  top_annotation = ha, 
  column_names_gp = gpar(fontsize = 0),
  show_column_dend = F,
  row_dend_width = unit(3,'mm'),
  row_dend_side = 'left')
pdf (file.path (hubs_dir,'Plots',paste0('hubs_cells_',metaGroupName,'_heatmap.pdf')), height=3)
hm
dev.off()


# Compute differential hub accessibility DHA ####
library (presto)
metaGroupName = 'celltype2'
metaGroupName2 = 'status'
meta_group = c('NK','T_cells')
comparison = unique(archp@cellColData[,metaGroupName2])
wlc_res = lapply (meta_group, function(x) 
  {
  hubsCell_mat_comp = hubsCell_mat[,as.character(archp@cellColData[,metaGroupName]) == x]
  comparison = as.character(archp@cellColData[,metaGroupName2])[match(colnames (hubsCell_mat_comp), rownames(archp@cellColData))]
  res = wilcoxauc (log2(hubsCell_mat_comp+1), comparison)  
  res[res$group == 'normal',]
  })


names (wlc_res) = meta_group
size_comp = lapply (seq_along(wlc_res), function(x)
   {
    tmp = wlc_res[[x]]
    tmp$width = width (hubs_obj$hubsCollapsed)
    tmp$group = ifelse (tmp$logFC < 0, 'tumor','normal')
    tmp$celltype = names(wlc_res)[x]
    tmp
   })
size_comp_df = do.call (rbind, size_comp)
size_comp_df_t = size_comp_df[size_comp_df$celltype == 'T_cells',]
head (size_comp_df_t[order(-size_comp_df_t$logFC),])

# Check hub size distribution between normal and tumor ####
#res$width = width(hubs_obj$hubsCollapsed)[match (res$feature, hubs_obj$hubs_id)]
vp = ggplot (size_comp_df, aes(x=group, y=logFC)) +
    geom_boxplot() + 
    facet_wrap (~ celltype) +
    ggtitle ('Hub size') + gtheme

pdf (file.path (hubs_dir, 'Plots','DAH_size.pdf'),4,4)
vp
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




