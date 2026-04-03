# Load functions for hub detection ####
source (file.path('..','..','git_repo','utils','knnGen.R'))
source (file.path('..','..','git_repo','utils','addCoax.R'))
source (file.path('..','..','git_repo','utils','Hubs_finder.R'))
source (file.path('..','..','git_repo','utils','hubs_track.R'))
#source (file.path('..','..','git_repo','utils','scATAC_functions.R'))

# Export bigiwg files ####
archp$celltype_status = paste0(archp$celltype2, '_', archp$status)
archp$celltype_sample = paste0(archp$celltype2, '_', archp$Sample)
metaGroupName = 'celltype_sample'
metaGroupName = 'Clusters_H'
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
metaGroupName = "Clusters_H"
cor_cutoff = 0.2
#max_dist = 12500
max_dist = 12500
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
dir.create(file.path (hubs_dir, 'Plots'), recursive=T)


# Generate cluster-aware knn groups ####
k= 50
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
run_coax = T
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

## Bind peak2genes results with hubs links ####
hubs_obj$peakLinks
p2gl = getPeak2GeneLinks (archp, corCutOff = 0.2)[[1]]
elementMetadata(p2gl) = elementMetadata(p2gl)[,-ncol(elementMetadata(p2gl))]
colnames(elementMetadata(p2gl)) = 'value'
hubs_obj$peakLinks2 = unlist(GRangesList(c(hubs_obj$peakLinks, p2gl)))

# Annotate also exhausted CD8 cells and generate heatmap of hubs x metagroup ####
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
  frags_in_sample = sapply (unique(archp@cellColData[,metaGroupName]), function(x) sum (archp$nFrags[as.character(archp@cellColData[,metaGroupName]) == x]))
  hubsSample_mat = t(t(hubsSample_mat) * (10^6 / frags_in_sample)) # scale
  saveRDS (hubsSample_mat, file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))
  } else {
  hubsSample_mat = readRDS (file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))  
  }
hubsSample_mat = as.data.frame (hubsSample_mat)

#ha = HeatmapAnnotation (size = anno_barplot(width (hubs_obj$hubsCollapsed), gp = gpar(color = "red"), height =  unit(8, "mm")))
int_genes = c('NR4A2','ICOS','CTLA4','PDCD1','HAVCR2')
int_genes_pos = unlist(sapply (int_genes, function(x) grep (x, hubs_obj$hubsCollapsed$gene)))
int_genes_label = paste0(hubs_obj$hubs_id, '_', hubs_obj$hubsCollapsed$gene)[unlist(sapply (int_genes, function(x) grep (x, hubs_obj$hubsCollapsed$gene)))]

ha2 = rowAnnotation(foo = anno_mark(at = int_genes_pos, 
    labels = int_genes_label, labels_gp = gpar(fontsize = 7, fontface='italic')))
hm = Heatmap (
  t(scale (t(hubsSample_mat))), 
  right_annotation = ha2, 
  #column_names_gp = gpar(fontsize = 3),
  show_column_dend = F,
  show_row_dend = F,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 0),
  column_names_rot = 45,
  #column_km = 5,
  #row_dend_width = unit(5,'mm'),
  row_dend_side = 'left',
  col = rev(palette_hubs_accessibility),
  border=T,
  name = 'Hubs')
pdf (file.path (hubs_dir,'Plots',paste0('hubs_',metaGroupName,'_heatmap.pdf')), height=4, width = 3.5)
hm
dev.off()

# # Compare with similarity using all peaks called ####
# if (!exists('pSE')) pSE = ArchR::getMatrixFromProject (archp, 'PeakMatrix')
# pSE = pSE[,rownames(archp@cellColData)]
# metaGroupName = 'celltype2'
# pmat = assays(pSE)[[1]]
# pmat_l = list()
# for (i in unique(archp@cellColData[,metaGroupName])) 
#   {
#   pmat_l[[i]] = data.frame (i = rowMeans (pmat[,as.character(archp@cellColData[,metaGroupName]) == i]))
#   colnames(pmat_l[[i]]) = i
#   }
# pmat_agr = do.call (cbind, pmat_l)
# rownames(pmat_agr) = as.character(rowRanges(pSE))
# pmat_agr_cor = cor (pmat_agr)

# hm = Heatmap (
#   scale (t(pmat_agr)), 
# #  top_annotation = ha, 
#   column_names_gp = gpar(fontsize = 0),
#   show_column_dend = T,
#   #column_km = 2,
#   #row_dend_width = unit(5,'mm'),
#   row_dend_side = 'left',
#   col = rev(palette_hubs_accessibility),
#   border=T,
#   name = 'peaks')
# pdf (file.path (hubs_dir,'Plots',paste0('peaks_',metaGroupName,'_heatmap.pdf')), height=2.2, width = 100)
# hm
# dev.off()





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
peak_reproducibility=2
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






# Compute differential hub accessibility between CD8 exhausted and rest of CD8 ####
library (presto)
all (colnames(hubsCell_mat) == rownames(archp@cellColData))
metaGroupName = 'celltype2'
hubsCell_mat_cd8 = hubsCell_mat[,as.character (archp@cellColData[,metaGroupName]) %in% c('CD8','CD8_exhausted')]
cd8_groups = as.character (archp@cellColData[,metaGroupName])[as.character (archp@cellColData[,metaGroupName]) %in% c('CD8','CD8_exhausted')]
res = wilcoxauc (log2(hubsCell_mat_cd8+1), cd8_groups)

res_l = lapply (split (res, res$group), function(x){
  tmp = x[x$logFC > 0,]
  tmp = tmp[order (tmp$pval),]
  tmp
})

head (res_l[['CD8_exhausted']],200)

# Plot only hubs up in exhaustion across TNK ####
hubsSample_mat_ext = hubsSample_mat[res_l[['CD8_exhausted']]$feature[res_l[['CD8_exhausted']]$padj < 0.00001],]
hubsSample_mat_ext = hubsSample_mat[head(res_l[['CD8_exhausted']]$feature,500),]
#ha = HeatmapAnnotation (size = anno_barplot(width (hubs_obj$hubsCollapsed), gp = gpar(color = "red"), height =  unit(8, "mm")))
hm = draw(Heatmap (
  scale (t(hubsSample_mat_ext)),
#  top_annotation = ha, 
  column_names_gp = gpar(fontsize = 3),
  column_names_rot = 45,
  show_column_dend = F,
  column_km = 5,
  #row_dend_width = unit(5,'mm'),
  row_dend_side = 'left',
  col = rev(palette_hubs_accessibility),
  border=T,
  name = 'Hubs'))
pdf (file.path (hubs_dir,'Plots',paste0('hubs_',metaGroupName,'_Ext_up_heatmap.pdf')), height=2.2,width=5)
hm
dev.off()

hm@cellColData  
### TF Enrichment in hubs high in KLRC1+ and CD8 exhausted ####
tf_enr_l = list()
selected_hubs = hubs_obj$hubsCollapsed[unlist(c(column_order(hm)[[2]],column_order(hm)[[3]]))]
selected_hubs = hubs_obj$hubsCollapsed[grep ('SPRY1', hubs_obj$hubsCollapsed$gene)]
tf_match = getMatches (archp)
ct = 'CD8'
bg_peaks = readRDS (file.path('PeakCalls',paste0(ct,'-reproduciblePeaks.gr.rds')))
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)[queryHits(findOverlaps(tf_match,bg_peaks))]

selected_hubs_peaks = bg_peakSet[queryHits(findOverlaps(bg_peakSet, selected_hubs))]
#tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
selected_hubs_TF =  hyperMotif (
  selected_peaks = selected_hubs_peaks, 
  motifmatch = tf_match)
head (selected_hubs_TF,30)

# Plot NR4A3 deviation has showing HUB31 higher in KLRC1 and CD8 exhausted cells ####
tf_markers = c('NR4A3')
tf_markers = c('NR4A2')
markerMotifs = getFeatures (archp, select = paste(tf_markers, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs = grep ("z:", markerMotifs, value = TRUE)
#archp = addImputeWeights (archp)
TF_p = plotEmbedding(
    ArchRProj = archp, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP_H",
    pal = palette_deviation,
    imputeWeights = getImputeWeights(archp)
)

# Check expression of NR4A3 ####
metaGroupName = 'celltype2'
#metaGroupName = 'Sample'
tf = c('NR4A2')
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
ext_tf = data.frame (
TF = mMat[tf,],
celltype =as.character(archp@cellColData[,metaGroupName]),
sample = archp$Sample)

ct_order = unlist(lapply(split (ext_tf, ext_tf$celltype), function(x) mean(x$TF)))
ext_tf$celltype = factor (ext_tf$celltype, levels = names(ct_order)[order(-ct_order)])
metaGroup = 'Sample'
ext_tf_agr = aggregate (ext_tf[,1,drop=F], by = list(celltype= ext_tf$celltype, sample = ext_tf$sample), mean)

bp = ggplot (ext_tf_agr, aes (x = celltype, y = TF)) + 
geom_boxplot (aes (fill=celltype),
    linewidth = .2,
    width=.8,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.3, alpha=0.7
     ) + 
# geom_violin (aes (fill=celltype), trim=TRUE,size=2,
#     width=1,
#     scale='width',
#     linewidth = .2, alpha=0.7) +
scale_fill_manual (values= palette_tnk_cells) + 
#geom_point (position='identity', alpha=.2, color="grey44", size=.1) +
gtheme

stat.test = bp$data %>%
t_test(reformulate ('celltype', 'TF')) %>%
adjust_pvalue (method = "fdr") %>%
add_significance ()
stat.test = stat.test %>% add_xy_position (x = 'celltype', step.increase=.4)
bp = bp + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
bracket.nudge.y = -0.1, hide.ns = T,
label = "p.adj.signif") + NoLegend()

pdf (file.path ('Plots',paste0(tf,'_boxplots.pdf')),height=3.3,width=2.5)
bp
dev.off()

 
  metaGroupName = 'celltype3'
  markers = 'SPRY1'
  markers = 'NR4A2'
  markers = 'TOX2'
  #celltype_markers = c('WT1','CALB2','GATA4','MSLN','KRT5','KRT18','ITLN1','HP','SOX9')
archp$celltype3 =  archp$celltype2
archp$celltype3[archp$celltype3 == 'CD4'] =  c('C1_CD4')
archp$celltype3[archp$celltype3 == 'CD8'] =  c('C2_CD8')
archp$celltype3[archp$celltype3 == 'NK_FGFBP2'] =  c('C3_NK_FGFBP2')
archp$celltype3[archp$celltype3 == 'Tregs'] =  c('C4_Tregs')
archp$celltype3[archp$celltype3 == 'CD8_exhausted'] =  c('C5_CD8_exhausted')
archp$celltype3[archp$celltype3 == 'NK_KLRC1'] =  c('C6_NK_KLRC1')
palette_tnk_cells_ext2 = palette_tnk_cells
names (palette_tnk_cells_ext2) = c(
  'C2_CD8',
  'C1_CD4',
  'C4_Tregs',
  'C6_NK_KLRC1',
  'C3_NK_FGFBP2',
  'C5_CD8_exhausted')

meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp, 
    sizes = c(6, 1, 1, 1,1),
    groupBy = metaGroupName, 
    geneSymbol = markers,
    normMethod = "ReadsInTSS",
    plotSummary = c("bulkTrack", "featureTrack", 
        #"loopTrack", 
      "geneTrack", "hubTrack"),
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 50000,
    pal = palette_tnk_cells_ext2,
    #ylim=c(0,0.1),
    downstream = 200000,
    #loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
#    loops = getCoAccessibility (archp, corCutOff = 0.25),
    #  returnLoops = TRUE),
    useGroups= NULL,
    hubs = hubs_obj$peakLinks2
)
plotPDF (meso_markers, ArchRProj = archp, 
  width=5,height=3, 
  name =paste0(paste(markers, collapse='_'),'_coveragePlots.pdf'),
  addDOC = F)

metaGroupName = 'celltype2'
top_dah = data.frame (
gene = srt@assays$RNA@data[markers,],
group = srt@meta.data[,metaGroupName])
top_dah$group = factor (top_dah$group, levels = 
  rev (c('CD4','CD8','NK_FGFBP2','Tregs','CD8_exhausted',
  'NK_KLRC1','TFH')))
top_dah = na.omit(top_dah)
bp = ggplot (top_dah, aes (x = gene, y = group, fill = group)) + 
vlp + 
bxpv + 
scale_fill_manual (values = palette_tnk_cells) +
#geom_point (position='identity', alpha=.3, color="grey44", size=1) +
gtheme_no_rot

pdf (file.path ('Plots', paste0('scrna_',markers,'_boxplots.pdf')), height=4, width=4)
bp
dev.off()


# Check DAH between NK KLRC1+ and CD8 exhausted to find commonalities ####
cd8_ct = read.csv (file.path('..','..','CD8','scatac_ArchR','barcode_annotation.csv')) # get annotation from subclustered CD8 cells
archp$celltype2 = archp$celltype3
archp$celltype2[match(cd8_ct$barcode, rownames(archp@cellColData))] = cd8_ct$celltype
archp$celltype2[archp$celltype2 %in% c('C1','C2','C4')] = 'CD8'






