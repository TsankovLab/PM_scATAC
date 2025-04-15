# Load functions for hub detection ####
source (file.path('..','..','git_repo','utils','knnGen.R'))
source (file.path('..','..','git_repo','utils','addCoax.R'))
source (file.path('..','..','git_repo','utils','Hubs_finder.R'))
source (file.path('..','..','git_repo','utils','hubs_track.R'))

# Export bigiwg files ####
metaGroupName = 'celltype'
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
metaGroupName = "celltype"
cor_cutoff = 0.3
#max_dist = 12500
max_dist = 12500
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
dir.create(file.path (hubs_dir, 'Plots'), recursive=T)

# Generate cluster-aware knn groups ####
k= 100
#metaGroupName = 'celltype_status'

force = F
if (!file.exists(paste0 ('KNNs_',metaGroupName,'k_',k,'.rds')) | force)
  {
  KNNs = knnGen (
    ArchRProj = archp, 
    k = k,
    reducedDims = 'Harmony', 
    group = metaGroupName,
    overlapCutoff = 0.7,
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
pdf()
umap_knn = plotEmbedding (ArchRProj = archp, embedding = 'UMAP_H',
  colorBy = "cellColData", name = 'knn_groups',plotAs ='hex',
    baseSize=0, labelMeans=FALSE) + NoLegend() 
dev.off()
pdf (file.path(hubs_dir, 'Plots',paste0('knn_',k,'.pdf')), height=20, width=20)
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
    cores=2
    ) 
  saveRDS (hubs_obj, file.path(hubs_dir,'global_hubs_obj.rds'))
  hubs_regions = as.data.frame (hubs_obj$hubsCollapsed)
  hubs_regions$width = hubs_obj$hubs_id
  write.table (hubs_regions, file.path(hubs_dir, 'hub_regions.bed'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  } else {
  hubs_obj = readRDS (file.path(hubs_dir,'global_hubs_obj.rds'))  
  }



# Generate matrix of fragment counts of hubs x sample ####
metaGroupName = 'celltype'
if (!file.exists(file.path (hubs_dir,paste0('hubs_',metaGroupName,'_mat.rds'))))
  {
  fragments = unlist (getFragmentsFromProject (
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
  saveRDS (hubsSample_mat, file.path (hubs_dir,paste0('hubs_',metaGroupName,'_mat.rds')))
  } else {
  hubsSample_mat = readRDS (file.path (hubs_dir,paste0('hubs_',metaGroupName,'_mat.rds')))
  }
  hubsSample_mat = as.data.frame (hubsSample_mat)

#hubsSample_mat = hubsSample_mat[, !colnames (hubsSample_mat) %in% c('normal_Pericytes','LEC'),]
#ha = HeatmapAnnotation (size = anno_barplot(width (hubs_obj$hubsCollapsed), gp = gpar(color = "red"), height =  unit(8, "mm")))
#ha2 = HeatmapAnnotation (status = ifelse (grepl('tumor', colnames(hubsSample_mat)),'tumor','normal'), which='row')
hm = Heatmap (
  t(scale (t(hubsSample_mat))), 
 # top_annotation = ha, 
 # left_annotation = ha2,
  column_names_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 0),
  column_names_rot = 45,
  show_column_dend = F,
  row_dend_width = unit(3,'mm'),
  row_dend_side = 'left',
  col=rev (palette_fragments))

pdf (file.path (hubs_dir,'Plots',paste0('hubs_',metaGroupName,'_heatmap.pdf')), height=4, width=3)
hm
dev.off()

# Find hubs high in mesothelium and fibroblasts ####
hubsSample_mat_scaled = t(scale (t(hubsSample_mat)))
fib_meso_avg = rowMeans (hubsSample_mat_scaled[, c('Fibroblasts','Mesothelium')])
fib_meso_avg = fib_meso_avg[order(-fib_meso_avg)]


# Check top correlated hubs on browser track ####
top_hubs = as.character(head(hub_cor_levels$hub_id, 50))
metaGroupName='celltype'

TF = 'RUNX2'

pdf()
meso_markers = plotBrowserTrack2 (
    ArchRProj = archp, 
    sizes = c(6, 1, 1, 1,1,1),
    groupBy = metaGroupName, 
    #region = ext_range(hubs_obj$hubsCollapsed[match(top_hubs, hubs_obj$hubs_id)],50000,50000),
    region = ext_range (hubs_obj$hubsCollapsed[hubs_obj$hubs_id %in% head(names (fib_meso_avg))],10000,50000),
    sample_levels = NULL,
    #geneSymbol = TF,
    genelabelsize=2,
    #geneSymbol = TF,
    normMethod = "ReadsInTSS",
    scCellsMax=100,
    hubs_regions = hubs_obj$hubsCollapsed,
    plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    #upstream = 1000,
    pal = palette_stroma,
    minCells = 25,
    #ylim=c(0,0.4),
    #downstream = 100000,
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
  name =paste0('cor_hubs_sarc_TFs_coveragePlots.pdf'),
  addDOC = F) 

metaGroupName = 'sampleID3'
srt_NN = srt[,srt@meta.data[,metaGroupName] != 'normal_pleura']
hub_id = 'HUB7701'
hub_gr = hubs_obj$hubsCollapsed[match (hub_id, hubs_obj$hubs_id)]
genes_in_region = unique(getPeakSet (archp)[subjectHits (findOverlaps (hub_gr, getPeakSet (archp)))]$nearestGene)
genes_in_region = c(genes_in_region)
top_dah = data.frame (
gene = colMeans (srt_NN@assays$RNA@data[rownames(srt_NN) %in% genes_in_region,,drop=F]),
group = srt_NN@meta.data[,metaGroupName])
top_dah$group = factor (top_dah$group, levels = rev(sample_sarc_order_levels[!sample_sarc_order_levels %in% c('P3','P13')]))
top_dah = na.omit(top_dah)
bp = ggplot (top_dah, aes (x = gene, y = group, fill = group)) + 
vlp + 
bxpv + 
scale_fill_manual (values = palette_sample) +
#geom_point (position='identity', alpha=.3, color="grey44", size=1) +
gtheme_no_rot

pdf (file.path ('Plots', paste0('scrna_region_boxplots.pdf')), height=4, width=4)
bp
dev.off()


# Plot hubs fibroblasts vs mesothelium in scatterplot ####

# Import peaks mesothelium and fibroblasts ####
fib_peaks = readRDS ('PeakCalls/Fibroblasts-reproduciblePeaks.gr.rds')
fib_peaks = getPeakSet (archp)[queryHits(findOverlaps(getPeakSet(archp),fib_peaks))]
meso_peaks = readRDS ('PeakCalls/Mesothelium-reproduciblePeaks.gr.rds')
meso_peaks = getPeakSet (archp)[queryHits(findOverlaps(getPeakSet(archp),meso_peaks))]


  if(!exists('fragments')) fragments = unlist (getFragmentsFromProject (
    ArchRProj = archp))
  
  fib_mat = matrix (ncol = 1, nrow = length(fib_peaks))
  colnames (fib_mat) = 'Fibroblasts'
  rownames (fib_mat) = as.character(fib_peaks)
    fragments_in_sample = fragments[fragments$RG %in% rownames(archp@cellColData)[archp$celltype == 'Fibroblasts']]  
    fragments_in_sample_in_peaks = countOverlaps (fib_peaks, fragments_in_sample)
    fib_mat[,1] = fragments_in_sample_in_peaks
  frags_in_sample = sum (archp$nFrags[as.character(archp$celltype) == 'Fibroblasts'])
  fib_mat = t(t(fib_mat) * (10^6 / frags_in_sample)) # scale

  meso_mat = matrix (ncol = 1, nrow = length(meso_peaks))
  colnames (meso_mat) = 'Mesothelium'
  rownames (meso_mat) = as.character(meso_peaks)
    fragments_in_sample = fragments[fragments$RG %in% rownames(archp@cellColData)[archp$celltype == 'Mesothelium']]  
    fragments_in_sample_in_peaks = countOverlaps (meso_peaks, fragments_in_sample)
    meso_mat[,1] = fragments_in_sample_in_peaks
  frags_in_sample = sum (archp$nFrags[as.character(archp$celltype) == 'Mesothelium'])
  meso_mat = t(t(meso_mat) * (10^6 / frags_in_sample)) # scale

# Combine peak mats
meso_fib_mat = cbind (as.data.frame(meso_mat),as.data.frame (fib_mat)[rownames(as.data.frame(meso_mat)),])
colnames (meso_fib_mat) = c('Mesothelium','Fibroblasts')

library (ggpointdensity)
library (smplot2)
mf_p = ggplot (log2(hubsSample_mat[,c('Fibroblasts','Mesothelium')]+1),
  aes (x = Fibroblasts, y = Mesothelium)) + 
geom_point (color = 'grey22',alpha=0.4) +  
geom_pointdensity (alpha=1, size=.1) +
scale_color_viridis (option='F') +
# geom_smooth(
#   method = "lm", 
#   color = "white", 
#   fill = "white", 
#   se = T, 
#   linetype = "dashed", 
#   size = .4
# ) + # Regression line with confidence interval
labs(
  title = "cHubs",
  #subtitle = "Scatterplot with Regression Line and Correlation Coefficient",
  x = "Fibroblasts",
  y = "Mesothelium") +
  sm_statCorr(linetype = 'dashed',color = "white")+ 
  gtheme_no_rot
  # ) +
png (file.path ('Plots','Fibroblasts_mesothelium_hubs_scatter.png'),1000,1000,res=300)
mf_p
dev.off()

mf_p = ggplot (log2(meso_fib_mat[,c('Fibroblasts','Mesothelium')]+1),
  aes (x = Fibroblasts, y = Mesothelium)) + 
geom_point (color = 'grey22',alpha=0.4) +  
geom_pointdensity (alpha=1, size=.1) +
scale_color_viridis (option='F') +
# geom_smooth(
#   method = "lm", 
#   color = "white", 
#   fill = "white", 
#   se = T, 
#   linetype = "dashed", 
#   size = .4
# ) + # Regression line with confidence interval
labs(
  title = "Peaks",
  #subtitle = "Scatterplot with Regression Line and Correlation Coefficient",
  x = "Fibroblasts",
  y = "Mesothelium") +
  sm_statCorr(linetype = 'dashed',color = "white") + 
  gtheme_no_rot
  # ) +
png (file.path ('Plots','Fibroblasts_mesothelium_peaks_scatter.png'),1000,1000,res=300)
mf_p
dev.off()






# Generate matrix of fragment counts of hubs x barcodes ####
if (!file.exists(file.path (hubs_dir,'hubs_cells_mat.rds')))
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
  hubsCell_mat = t(t(hubsCell_mat) * (10^6 / archp$nFrags)) # scale
  saveRDS (hubsCell_mat, file.path (hubs_dir,'hubs_cells_mat.rds'))
  } else {
  hubsCell_mat = readRDS (file.path (hubs_dir,'hubs_cells_mat.rds'))  
  }
hubsCell_mat = as.data.frame (hubsCell_mat)




# Compute differential hub accessibility DHA ####
library (presto)
metaGroupName = 'celltype'

res = wilcoxauc (log2(hubsCell_mat+1), archp$celltype)
res = res[res$logFC > 0,]
res_df = lapply (split (res, res$group), function(x) x[order(x$padj),])
res_df = do.call (rbind, res_df)
res_df$genes = hubs_obj$hubsCollapsed$gene[match(res_df$feature, hubs_obj$hubs_id)]
head (res_df[res_df$group == 'LEC',])



TF = 'ETS1'
TF = 'MIR4777'
metaGroupName = 'celltype'

pdf()
#archp$fetal_sample = paste0(archp$Sample, archp$fetal_group)
#metaGroupName = 'fetal_group'
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp,#[!archp$Sample3 %in% c('P11_HOX')], 
    #group_order = sample_levels, 
#    ylim = c(0,0.30),
    groupBy = metaGroupName, 
    hubs_regions = hubs_obj$hubsCollapsed,
    #sample_levels = sample_sarc_order,
    minCells = 10,
    geneSymbol = TF,
    plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
    #pal = palette_sample,
    pal = palette_stroma,
    threads=1,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 100000,
    downstream = 100000,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
    #loops = getCoAccessibility (archp, corCutOff = 0.3,
    #  returnLoops = TRUE),
    useGroups= NULL
)
dev.off()
plotPDF (meso_markers, ArchRProj = archp,height=3.5, width=6, name =paste0('MPM_markers_',TF,'_coveragePlots.pdf'),addDOC=F)
  