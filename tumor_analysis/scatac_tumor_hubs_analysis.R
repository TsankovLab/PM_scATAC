conda activate meso_scatac
R

set.seed(1234)

####### ANALYSIS of TUMOR compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR'
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
addArchRThreads (threads = 8)
addArchRGenome ("hg38")

if (!file.exists ('Save-ArchR-Project.rds')) 
  { source (file.path('..','..','PM_scATAC','scatac_tumor_create_ArchRobj.R'))
  } else {
 archp = loadArchRProject (projdir)   
  }

archp$Sample3 = archp$Sample2
archp$Sample3[archp$Clusters == 'C1'] = 'P11_HOX'


# Load RNA ####
srt = readRDS (file.path('..','scrna','srt.rds'))
srt$sampleID3[srt$sampleID3 %in% c('HU37','HU62')] = 'normal_pleura'
srt_NN = srt[,!srt$sampleID %in% c("HU37','HU62")]
sarc_order = read.csv (file.path('..','scrna','cnmf20_sarcomatoid_sample_order.csv'), row.names=1)
sarc_order = sarc_order[! sarc_order$sampleID %in% c('HU37','HU62'),]
sarc_order = rbind (data.frame (sampleID = 'normal_pleura', x = -1),sarc_order)
#archp$Sample2 = factor (archp$Sample2, levels = sarc_order$sampleID)

# Export bigiwg files ####
metaGroupName = 'SampleP11'
metaGroupName = 'Sample3'
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
metaGroupName = "Sample2"
cor_cutoff = 0.3
#max_dist = 12500
max_dist = 12500
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
dir.create(file.path (hubs_dir, 'Plots'), recursive=T)

# Generate cluster-aware knn groups ####
k= 100
metaGroupName = 'Sample3'
archp_NN = archp[archp$Sample3 != 'normal_pleura']

force = F
if (!file.exists(paste0 ('KNNs_',metaGroupName,'k_',k,'.rds')) | force)
  {
  KNNs = knnGen (
    ArchRProj = archp_NN, 
    k = k,
    reducedDims = 'IterativeLSI', 
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
archp_NN$knn_groups = KNNs_df$group[match (archp_NN$cellNames, KNNs_df$cell)]

# Plot KNNs on UMAP ####
pdf()
umap_knn = plotEmbedding (ArchRProj = archp_NN, embedding = 'UMAP',
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
  archp_NN = addCoAx (
    archp_NN, 
    KNNs,
    maxDist = max_dist)
  }

### Run hub finder ####
force=F
if (!file.exists (file.path(hubs_dir,'global_hubs_obj.rds')) | force)
  {
  hubs_obj = hubs_finder (
    ArchRProj = archp_NN, 
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
metaGroupName = 'Sample3'
if (!file.exists(file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds'))))
  {
  if (!exists ('fragments')) fragments = unlist (getFragmentsFromProject (
    ArchRProj = archp_NN))   
  hubsSample_mat = matrix (ncol = length(unique(archp_NN@cellColData[,metaGroupName])), nrow = length(hubs_obj$hubsCollapsed))
  colnames (hubsSample_mat) = unique(archp_NN@cellColData[,metaGroupName])
  rownames (hubsSample_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (unique(archp_NN@cellColData[,metaGroupName])))
  for (sam in unique(archp_NN@cellColData[,metaGroupName]))
    {
    pb$tick()  
    fragments_in_sample = fragments[fragments$RG %in% rownames(archp_NN@cellColData)[as.character(archp_NN@cellColData[,metaGroupName]) == sam]]  
    fragments_in_sample_in_hubs = countOverlaps (hubs_obj$hubsCollapsed, fragments_in_sample)
    hubsSample_mat[,sam] = fragments_in_sample_in_hubs
    }
  frags_in_sample = sapply (unique(archp_NN@cellColData[,metaGroupName]), function(x) sum (archp_NN$nFrags[as.character(archp_NN@cellColData[,metaGroupName]) == x]))
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
  #row_dend_width = unit(5,'mm'),
  row_dend_side = 'left',
  col = rev(palette_hubs_accessibility),
  border=T,
  name = 'Hubs')
pdf (file.path (hubs_dir,'Plots',paste0('hubs_',metaGroupName,'_heatmap.pdf')), height=3.2,width=6)
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



# # Check if there are hub size differences between normal and tumor sample ####
# library (presto)
# tumor_sample = unique (archp$Sample2)[unique(archp$Sample2) != 'normal_pleura']
# wlc_res = lapply (tumor_sample, function(x) 
#   {
#   hubsCell_mat_comp = hubsCell_mat[,archp$Sample2 %in% c(x,'normal_pleura')]
#   comparison = archp$Sample2[archp$Sample2 %in% c(x,'normal_pleura')]
#   res = wilcoxauc (log2(hubsCell_mat_comp+1), comparison)
#   res[res$group == 'normal_pleura',]
#   })
# names (wlc_res) = tumor_sample
# size_comp = lapply (wlc_res, function(x)
#    {
#     x$width = width (hubs_obj$hubsCollapsed)
#     x$group = ifelse (x$logFC < 0, 'normal_pleura','tumor')
#     x
#    })

# size_comp_df = do.call (rbind, size_comp)



### Try correlating BRD4 with Hubs load per cell ####
archp = addImputeWeights (archp)
seGS <- getMatrixFromProject (archp)
celltype_markers = c('BRD4',rowData (seGS)$name[grep ('MED', rowData (seGS)$name)])
seGS = seGS[rowData (seGS)$name %in% celltype_markers,]
matGS <- imputeMatrix (assay(seGS), getImputeWeights(archp))
rownames(matGS) = rowData (seGS)$name
#rownames (matGS) = celltype_markers

hubsCell_mean = hubsCell_mat[,colnames(matGS)]
hubsCell_mean = colMeans (log2(hubsCell_mean+1))
matGS = t (matGS)
# sox9_cor = as.data.frame (t(cor (matGS[,'SOX9'], hubsCell_mat)))
# sox9_cor$gene = hubs_obj$hubsCollapsed$gene[match (rownames(sox9_cor), hubs_obj$hubs_id)]
# sox9_cor = sox9_cor[order (-sox9_cor$V1),, drop=F]
# sox9_cor$region = as.character (hubs_obj$hubsCollapsed)[match (rownames(sox9_cor), hubs_obj$hubs_id)]
# head (sox9_cor,20)

TF_hub_cor = lapply (unique(archp$Sample2)[unique(archp$Sample2) != 'normal_pleura'], function(x) 
  cor (as.matrix(matGS[grepl (paste0(x,'#'), rownames(matGS)),]), 
    hubsCell_mean[grepl (paste0(x,'#'), names(hubsCell_mean))], method = 'pearson'))

TF_hub_cor_o = lapply (TF_hub_cor, function(x) x[order(-x[,1]),])
hubs_genes_cor = sapply (TF_hub_cor, function(x) x[grepl ('BRD4', rownames(x))])
names (hubs_genes_cor) = unique(archp$Sample2)[unique(archp$Sample2) != 'normal_pleura']
sapply (TF_hub_cor, function(x) x[grepl ('^MED1$', rownames(x))])



# # Compute differential hub accessibility (DHA) between normal and tumor samples ####
library (presto)
archp$tumor_vs_normal = ifelse (archp$Sample2 == 'normal_pleura', 'normal','tumor')
#hubsCell_mat_comp = hubsCell_mat[,archp$Sample2 %in% c('P1','normal_pleura')]
#comparison = archp$Sample2[archp$Sample2 %in% c('P1','normal_pleura')]
#hubsCell_mat_comp2 = hubsCell_mat[,archp$Sample2 %in% c('P1','P8','P5','P4')]
#comparison2 = archp$Sample2[archp$Sample2 %in% c('P1','P8','P5','P4')] 
#comparison2 = ifelse(comparison2 == 'P1','P1','epit')  
res = wilcoxauc (log2(hubsCell_mat+1), archp$tumor_vs_normal)
res = res[res$group == 'tumor',]
res = res[order (res$padj),]
head (res[order(-res$logFC),])

logfcThreshold = 1
pvalAdjTrheshold = 0.05
res$sig = ifelse (abs(res$logFC) > logfcThreshold & res$padj < pvalAdjTrheshold, 1,0)
res$sig = res$sig * sign (res$logFC)
res$sig = as.character(res$sig)
res_filtered = res[abs(res$logFC) > logfcThreshold & res$padj < pvalAdjTrheshold,]
res_filtered = head (res_filtered$feature[order (-abs(res_filtered$logFC))],20)
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
    scale_color_manual (values = c("0"='grey77',"-1"='#666666FF',"1"='#F8A02EFF')) + theme_light()

pdf (file.path (hubs_dir,'Plots', 'DAH_volcano.pdf'),4,4)
vp
dev.off()


# # Compute differential hub accessibility (DHA) between cancers ####
library (presto)
res = wilcoxauc (log2(hubsCell_mat+1), archp$Sample3)
res_sub = res[res$group == 'P4',]
res_sub = res_sub[order (res_sub$padj),]
head (res_sub[order(-res_sub$logFC),])
as.character(seqnames(hubs_obj$hubsCollapsed[match(head (res_sub$feature[order(-res_sub$logFC)],20), hubs_obj$hubs_id)]))
hubs_obj$hubsCollapsed[match(head (res_sub$feature[order(-res_sub$logFC)],20), hubs_obj$hubs_id)]


logfcThreshold = 1
pvalAdjTrheshold = 0.05
res$sig = ifelse (abs(res$logFC) > logfcThreshold & res$padj < pvalAdjTrheshold, 1,0)
res$sig = res$sig * sign (res$logFC)
res$sig = as.character(res$sig)
res_filtered = res[abs(res$logFC) > logfcThreshold & res$padj < pvalAdjTrheshold,]
res_filtered = head (res_filtered$feature[order (-abs(res_filtered$logFC))],20)
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
    scale_color_manual (values = c("0"='grey77',"-1"='#666666FF',"1"='#F8A02EFF')) + theme_light()

pdf (file.path (hubs_dir,'Plots', 'DAH_volcano.pdf'),4,4)
vp
dev.off()


# Check hub size distribution between normal and tumor ####
res$width = width(hubs_obj$hubsCollapsed)[match (res$feature, hubs_obj$hubs_id)]
vp = ggplot (res, aes(x=log10(width), y=logFC)) +
    geom_point(size=1, shape=19, aes (color = sig), alpha=.5) + 
    ggtitle ('Hub size')

pdf (file.path (hubs_dir, 'Plots','DAH_size.pdf'),4,4)
vp
dev.off()

mega_hubs = range (hubs_obj$hubsCollapsed[1:4])

# Check mega hubs found in P11 ####
metaGroupName='Sample3'
matching_samples=c('normal_pleura','P1','P4','P5','P8','P11','P11_HOX','P12','P14')
pdf()
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp[archp$Sample3 %in% matching_samples], 
    sizes = c(6, 1, 1, 1,1,1),
    groupBy = metaGroupName, 
    region = mega_hubs,
    genelabelsize=0,
    #geneSymbol = TF,
    normMethod = "ReadsInTSS",
    scCellsMax=3000,
    plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 50000,
    pal = palette_sample,
    #ylim=c(0,0.1),
    downstream = 300000,
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
  name =paste0('region_coveragePlots.pdf'),
  addDOC = F)

metaGroupName = 'sampleID3'
genes_in_region = unique(getPeakSet (archp)[subjectHits (findOverlaps (mega_hubs, getPeakSet (archp)))]$nearestGene)
top_dah = data.frame (
gene = colMeans (srt@assays$RNA@data[rownames(srt) %in% genes_in_region,]),
group = srt@meta.data[,metaGroupName])
top_dah$group = factor (top_dah$group, levels =rev(c('normal_pleura','P1','P4','P5','P8','P11','P11_HOX','P12','P14')))
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


# Export ranges of top HUBs
saveRDS (mega_hubs, 'P11_chr18_region.rds')


# Look into mega hubs identified in P11 HOX- cluster ####
# Map large hubs to UMAP ####
hubsCell_mat = hubsCell_mat[rownames(archp@cellColData),]
archp@cellColData = cbind(archp@cellColData, hubsCell_mat[,c('HUB1','HUB2','HUB3','HUB4','HUB5','HUB6','HUB7','HUB8')])
umap_p1 =  plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
 name = c('HUB1','HUB2','HUB3','HUB4','HUB5','HUB6','HUB7','HUB8'), embedding = "UMAP")
  
pdf (file.path(hubs_dir,'Plots','HUB1_2_3_umap.pdf'), 15,15)
wrap_plots (umap_p1, ncol=4)
dev.off()

### TF Enrichment in peaks in large hubs ####
tf_match = getMatches (archp)
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)
mega_hubs = range (hubs_obj$hubsCollapsed[c(1:8)])
mega_hubs_peaks = bg_peakSet[queryHits(findOverlaps(bg_peakSet, mega_hubs))]
#tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
mega_hubs_TF =  hyperMotif (
  selected_peaks = mega_hubs_peaks, 
  motifmatch = tf_match)

head (mega_hubs_TF, 20)

### Correlate mega hubs on chr18 with other genes / TF ####
archp = addImputeWeights (archp)
#celltype_markers = c('WT1', 'CALB2','GATA4','HP','SOX9','MESP1','SOX6','TWIST1','SNAI2')
celltype_markers = rownames (head (mega_hubs_TF, 20))
seGS = getMatrixFromProject (archp)
seGS = seGS[rowData (seGS)$name %in% celltype_markers,]
matGS = imputeMatrix (assay(seGS), getImputeWeights(archp))
rownames (matGS) = rowData (seGS)$name [rowData (seGS)$name %in% celltype_markers]
matGS = t(matGS)

hubsCell_mat2 = hubsCell_mat[,colnames(matGS)]
hubsCell_mat2 = log2(t(hubsCell_mat2)+1)
  cor (as.matrix(hubsCell_mat2[archp$Clusters %in% c('C15'),c('HUB1','HUB2','HUB3','HUB4','HUB5')]),
  as.matrix(matGS)[archp$Clusters %in% c('C15'),])
# sox9_cor = as.data.frame (t(cor (matGS[


# res2 = wilcoxauc (log2(hubsCell_mat_comp2+1), comparison2)
# res2 = res2[res2$logFC >0,]
# res2 = res2[order(res2$padj),]
# head (res2, 20)
# comparison2 = archp$Sample2 == 'P1'
# res_P1_rest = wilcoxauc(log2(hubsCell_mat[,comparison2 %in% c('P1','normal_pleura')]+1), comparison2[comparison2 %in% c('P1','normal_pleura')])

#res_p1 = res[res$group == 'P1',]
#res_p1$region = as.character(hubs_obj$hubsCollapsed)
#res_p1 = res_p1[res_p1$padj < 0.05 & res_p1$logFC > 0 & res_p1$avgExpr > 3,]

# # Intersect with hubs higher in P1 vs epithelioid samples ####
# res2_p1 = res2[res2$group == 'P1',]
# res2_p1$region = as.character(hubs_obj$hubsCollapsed)
# res2_p1 = res2_p1[res2_p1$padj < 0.05 & res2_p1$logFC > 0 & res2_p1$avgExpr > 3,]
# res_p1 = res_p1[res_p1$feature %in% res2_p1$feature,]


# # intersect with correlated hubs to SOX9 ####
# res_p1$SOX9_cor = TF_hub_cor2[res_p1$feature,]$V1
# res_p1 = res_p1[order (-res_p1$SOX9_cor),]
# res_p1 = res_p1[res_p1$SOX9_cor > 0.1, ]

# tf_match = getMatches (archp)
# colnames(tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
# tf_match = tf_match[queryHits(findOverlaps(rowRanges(tf_match),hubs_obj$hubsCollapsed)),]
# top_cor_hubs = res_p1$feature
# top_cor_hubs_peaks = lapply(top_cor_hubs, function(x) getPeakSet(archp) [queryHits(findOverlaps(getPeakSet(archp), hubs_obj$hubsCollapsed[which(hubs_obj$hubs_id %in% x)]))])
# top_cor_hubs_TF = lapply (top_cor_hubs_peaks, function(x) hyperMotif (
#   selected_peaks = x,
#   motifmatch = tf_match))
# names (top_cor_hubs_TF) = top_cor_hubs

# top_cor_hubs_TF = lapply (top_cor_hubs_TF, function(x) x[rownames(top_cor_hubs_TF[[1]]), ])
# top_cor_hubs_TF_df = do.call (cbind, top_cor_hubs_TF)
# top_cor_hubs_TF_df = top_cor_hubs_TF_df[, grep ('padj', colnames(top_cor_hubs_TF_df))]
# top_cor_hubs_TF_df[top_cor_hubs_TF_df > 0.05] = 1
# top_cor_hubs_TF_df = -log10(top_cor_hubs_TF_df)
# top_cor_hubs_TF_df = top_cor_hubs_TF_df[rowSums (top_cor_hubs_TF_df) != 0, ]
# top_cor_hubs_TF_df[sapply(top_cor_hubs_TF_df, is.infinite)] <- 300

# TF_ht = Heatmap (top_cor_hubs_TF_df, row_names_gp = gpar (fontsize=3), column_names_gp = gpar (fontsize=5))

# pdf (paste0('Plots/TF_top_DHA_P1_hubs_heatmap.pdf'),width = 40,height=25)
# print (TF_ht)
# dev.off()

# top_cor_hubs_TF_df['SOX9',]
# res_p1[res_p1$feature == 'HUB22404',]
# # all peaks of top hubs
# top_cor_hubs = res_p1$feature[1:1000]
# top_cor_hubs_peaks = getPeakSet(archp)[queryHits(findOverlaps(getPeakSet(archp), hubs_obj$hubsCollapsed[which(hubs_obj$hubs_id %in% top_cor_hubs)]))]
# top_cor_hubs_TF = hyperMotif (
#   selected_peaks = top_cor_hubs_peaks,
#   motifmatch = tf_match)
# head (top_cor_hubs_TF[order (top_cor_hubs_TF$padj),],100)






# library (presto)
# pSE = getMatrixFromProject (archp, useMatrix = 'PeakMatrix')
# pMat = assays (pSE)[[1]]
# hubs_peaks_idx1 = subjectHits (findOverlaps (hubs_obj$peaksMerged,rowRanges(pSE)))
# pMat = pMat[hubs_peaks_idx1,]

# pMat_comp = pMat[,archp$Sample2 %in% c('P1','normal_pleura')]
# comparison = archp$Sample2[archp$Sample2 %in% c('P1','normal_pleura')]
# pMat_comp2 = pMat[,archp$Sample2 %in% c('P1','P8','P5','P4')]
# comparison2 = archp$Sample2[archp$Sample2 %in% c('P1','P8','P5','P4')] 
# comparison2 = ifelse(comparison2 == 'P1','P1','epit')  
# res = wilcoxauc (pMat_comp, comparison)
# res2 = wilcoxauc (log2(pMat_comp2), comparison2)

# # comparison2 = archp$Sample2 == 'P1'
# # res_P1_rest = wilcoxauc(log2(hubsCell_mat[,comparison2 %in% c('P1','normal_pleura')]+1), comparison2[comparison2 %in% c('P1','normal_pleura')])

# res_p1 = res[res$group == 'P1',]
# res_p1 = res_p1[res_p1$avgExpr != 0,]
# #res_p1$region = as.character(hubs_obj$hubsCollapsed)

# # Try using fgsea ####
# TF = 'MESP1'
# tf_matches = getMatches (archp)
# TF = colnames(tf_matches)[grep (TF, colnames(tf_matches))]

# hubs_peaks_idx2 = subjectHits (findOverlaps (hubs_obj$peaksMerged, rowRanges(tf_matches)))
# tf_matches = tf_matches[hubs_peaks_idx2, colnames(tf_matches) == TF]
# tf_matches_hit = which (assays (tf_matches)[[1]][,1])
# hub_peaks_rank = setNames (-log10(res_p1$pval) * res_p1$logFC, paste0('h',1:nrow(res_p1)))
# #tf_matches_peaks = unique(queryHits (findOverlaps (hubs_obj$hubsCollapsed[match(TF_hubs_id, hubs_obj$hubs_id)], rowRanges(tf_matches)[tf_matches_hit])))
# pathways = list(SOX9 = paste0('h',tf_matches_hit))

# library (fgsea)
# #hub_peaks_rank = na.omit (hub_peaks_rank)
# fgseaRes = fgseaMultilevel (pathways, 
#           hub_peaks_rank#, 
#           #minSize=15, 
#           #maxSize=1500,
#           #BPPARAM = NULL
#           )


# tf_match = getMatches (archp)
# colnames(tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
# bg_peakSet = rowRanges (tf_match)
# top_cor_hubs = head (sapply (rownames(TF_hub_cor2), function(x) unlist(strsplit (x, '_'))[1]),100)
# top_cor_hubs_peaks = lapply(top_cor_hubs, function(x) getPeakSet(archp) [queryHits(findOverlaps(getPeakSet(archp), hubs_obj$hubsCollapsed[which(hubs_obj$hubs_id %in% x)]))])
# top_cor_hubs_TF = lapply (top_cor_hubs_peaks, function(x) hyperMotif (
#   selected_peaks = x,
#   motifmatch = tf_match))
  
# top_cor_hubs_TF = lapply (top_cor_hubs_TF, function(x) x[rownames(top_cor_hubs_TF[[1]]), ])
#   top_cor_hubs_TF_df = do.call (cbind, top_cor_hubs_TF)
#   top_cor_hubs_TF_df = top_cor_hubs_TF_df[, grep ('padj', colnames(top_cor_hubs_TF_df))]
#   top_cor_hubs_TF_df[top_cor_hubs_TF_df > 0.05] = 1
#   top_cor_hubs_TF_df = -log10(top_cor_hubs_TF_df)
#   top_cor_hubs_TF_df = top_cor_hubs_TF_df[rowSums (top_cor_hubs_TF_df) != 0, ]
#   top_cor_hubs_TF_df[sapply(top_cor_hubs_TF_df, is.infinite)] <- 300
  
#   TF_ht = Heatmap (top_cor_hubs_TF_df, row_names_gp = gpar (fontsize=3), column_names_gp = gpar (fontsize=5))
  
#   pdf (paste0('Plots/TF_top_cor_hubs_heatmap.pdf'),width = 15,height=15)
#   print (TF_ht)
#   dev.off()






# # Try using clusters as KNN to improve accuracy of correlations ####
# # DOESNT WORK! ####
# fragments = unlist(getFragmentsFromProject (
#   ArchRProj = archp))  
  
# ###--- Generate matrces of hubs x cells / knn from collapsed non-redundant hubs to identify differential hubs  ---###
# metaGroupName = 'Clusters'
# metaGroup = as.character (archp@cellColData[,metaGroupName])
# metaGroup_df = data.frame (barcode = rownames(archp@cellColData), metaGroup = metaGroup)

# addArchRThreads (threads = 1) 
# peaks_sb = getGroupSE (archp, useMatrix = 'PeakMatrix', groupBy = 'Clusters')
# peaks_sb_mat = assays (peaks_sb)[[1]]
# gene_sb = getGroupSE (archp, useMatrix = 'GeneScoreMatrix', groupBy = 'Clusters')
# gene_sb_mat = assays (gene_sb)[[1]]
# rownames (gene_sb_mat) = rowData(gene_sb)$name

# cor (gene_sb_mat['SOX9',],gene_sb_mat['VIM',], method='spearman')

# # Correlate with SOX9 genescore
# if (!any (ls() == 'gsSE')) gsSE = ArchR::getMatrixFromProject (archp, useMatrix = 'GeneScoreMatrix')
# gsSE = gsSE[, archp$cellNames]
# gsMat = assays (gsSE)[[1]]
# rownames (gsMat) = rowData (gsSE)$name

# head (KNNs_df)
# gsMat = gsMat[c(tf_name_selected, 'AXL','VIM','CALB2','ITLN1'),]
# gsMat = as.data.frame (t(gsMat))
# gsMat = gsMat[unlist(KNNs),]
# gsMat$cellGroup = KNNs_df$group
# gsMat = aggregate (.~ cellGroup, data = gsMat, FUN = mean)
# rownames (gsMat) = gsMat[,1]
# gsMat = gsMat[,-1]

# hubsKnn_mat = t(hubsKnn_mat)
# hubsKnn_mat = as.data.frame (hubsKnn_mat)
# rownames (hubsKnn_mat) = rownames (gsMat)
# TF_hub_cor = lapply (unique(KNNs_df$group2), function(x) cor (gsMat[unique(KNNs_df$group[KNNs_df$group2 == x]),'SOX9'], hubsKnn_mat[unique(KNNs_df$group[KNNs_df$group2 == x]),], method = 'spearman'))
# names (TF_hub_cor) = unique(KNNs_df$group2)

# TF_hub_cor2 = as.data.frame (t(TF_hub_cor[['P1']]))
# TF_hub_cor2 = TF_hub_cor2[order (-TF_hub_cor2$V1),, drop=F]
# head (TF_hub_cor2, 10)
# rownames (TF_hub_cor2) = sapply (rownames(TF_hub_cor2), function(x) unlist(strsplit(x, '_'))[1])


# KNNs_df2 = KNNs_df[!duplicated(KNNs_df$group),]
# table (KNNs_df2$group2)
# a = lapply (unique(KNNs_df$group2), function(x) cor (gsMat[unique(KNNs_df$group[KNNs_df$group2 == x]),'MESP1'], gsMat[unique(KNNs_df$group[KNNs_df$group2 == x]),'SOX9'], method='spearman'))
# names(a) = unique(KNNs_df$group2)

# pdf ('Plots/cor_genescore_genes.pdf')
# plot (gsMat[,'SOX9'], gsMat[,'VIM'])
# dev.off 

# # Check enrichment of TF in hubs using fgsea
# TF = 'MESP1'
# hubs_id = sapply (rownames(TF_hub_cor2), function(x) unlist(strsplit (x, '_'))[1])
# tf_matches = getMatches (archp)
# TF = colnames(tf_matches)[grep (TF, colnames(tf_matches))]


# hub_peaks = rep (hubs_obj$hubs_id, sapply (hubs_obj$hubsClusters[[1]], nrow))
# hubs_obj$peaksMerged$hub_id = hub_peaks
# hubs_peaks_idx = subjectHits (findOverlaps (hubs_obj$peaksMerged,rowRanges(tf_matches)))
# tf_matches = tf_matches[hubs_peaks_idx, colnames(tf_matches) == TF]
# tf_matches_hit = which (assays (tf_matches)[[1]][,1])
# hub_peaks_cor = setNames (TF_hub_cor2[hubs_obj$peaksMerged$hub_id,]$V1, paste0('h',1:nrow(tf_matches))) 
# #tf_matches_peaks = unique(queryHits (findOverlaps (hubs_obj$hubsCollapsed[match(TF_hubs_id, hubs_obj$hubs_id)], rowRanges(tf_matches)[tf_matches_hit])))
# pathways = list(SOX9 = paste0('h',tf_matches_hit))

# library (fgsea)
# hub_peaks_cor = na.omit (hub_peaks_cor)
# fgseaRes = fgseaMultilevel (pathways, 
#           hub_peaks_cor#, 
#           #minSize=15, 
          #maxSize=1500,
          #BPPARAM = NULL
          )



# # Compute differential hub accessibility DHA ####
# library (presto)
# hubsCell_mat_comp = hubsCell_mat[,archp$Sample2 %in% c('P1','normal_pleura')]
# comparison = archp$Sample2[archp$Sample2 %in% c('P1','normal_pleura')]
# hubsCell_mat_comp2 = hubsCell_mat[,archp$Sample2 %in% c('P1','P8','P5','P4')]
# comparison2 = archp$Sample2[archp$Sample2 %in% c('P1','P8','P5','P4')] 
# comparison2 = ifelse(comparison2 == 'P1','P1','epit')  
# res = wilcoxauc (log2(hubsCell_mat_comp+1), comparison)
# res2 = wilcoxauc (log2(hubsCell_mat_comp2+1), comparison2)

# # comparison2 = archp$Sample2 == 'P1'
# # res_P1_rest = wilcoxauc(log2(hubsCell_mat[,comparison2 %in% c('P1','normal_pleura')]+1), comparison2[comparison2 %in% c('P1','normal_pleura')])

# res_p1 = res[res$group == 'P1',]
# res_p1$region = as.character(hubs_obj$hubsCollapsed)
# #res_p1 = res_p1[res_p1$padj < 0.05 & res_p1$logFC > 0 & res_p1$avgExpr > 3,]

# # # Intersect with hubs higher in P1 vs epithelioid samples ####
# # res2_p1 = res2[res2$group == 'P1',]
# # #res2_p1$region = as.character(hubs_obj$hubsCollapsed)
# # #res2_p1 = res2_p1[res2_p1$padj < 0.05 & res2_p1$logFC > 0 & res2_p1$avgExpr > 3,]
# # res_p1 = res_p1[res_p1$feature %in% res2_p1$feature,]


# # # intersect with correlated hubs to SOX9 ####
# # res_p1$SOX9_cor = TF_hub_cor2[res_p1$feature,]$V1
# # res_p1 = res_p1[order (-res_p1$SOX9_cor),]
# # res_p1 = res_p1[res_p1$SOX9_cor > 0.1, ]

# TF = 'SOX9'
# tf_matches = getMatches (archp)
# TF = colnames(tf_matches)[grep (TF, colnames(tf_matches))]

# hub_peaks = rep (hubs_obj$hubs_id, sapply (hubs_obj$hubsClusters[[1]], nrow))
# hubs_obj$peaksMerged$hub_id = hub_peaks

# tf_matches = tf_matches[subjectHits(findOverlaps(hubs_obj$peaksMerged, rowRanges(tf_matches))),]
# tf_matches = tf_matches[,colnames(tf_matches) == TF]
# tf_matches_mat = as.data.frame(assays(tf_matches)[[1]])
# tf_matches_mat$hubs_id = hubs_obj$peaksMerged$hub_id
# tf_matches_mat = aggregate (. ~ hubs_id, data = tf_matches_mat, FUN = sum)
# tf_matches_mat = tf_matches_mat[match(res_p1$feature, tf_matches_mat$hubs_id),]
# tf_matches_mat[tf_matches_mat >0] = 1
# hub_rank = setNames (res_p1$logFC, res_p1$feature)
# pathways = list(res_p1$feature[tf_matches_mat[,2] > 0])
# names (pathways) = TF
# library (fgsea)

# table (hub_rank[pathways[[1]]]>0)
# table (res_p1$logFC>0)
# #hub_peaks_rank = na.omit (hub_peaks_rank)
# fgseaRes = fgseaMultilevel (pathways, 
#           hub_rank#, 
#           #minSize=15, 
#           #maxSize=1500,
#           #BPPARAM = NULL
#           )

# pdf (paste0('Plots/hubs_enrichmentplot_',TF,'.pdf'))
# plotEnrichment(pathways[[TF]],
#                hub_rank) + labs(title=TF)
# dev.off()



### Try correlating hubs to imputed gene score of TFs ####
archp = addImputeWeights (archp)
celltype_markers = c('WT1','CALB2','GATA4','HP','SOX9','MESP1','SOX6','TWIST1','SNAI2')
seGS <- getMatrixFromProject (archp)
seGS = seGS[rowData (seGS)$name %in% celltype_markers,]
matGS <- imputeMatrix (assay(seGS), getImputeWeights(archp))
rownames (matGS) = celltype_markers

hubsCell_mat = hubsCell_mat[,colnames(matGS)]
hubsCell_mat = log2(t(hubsCell_mat)+1)
matGS = t(matGS)
# sox9_cor = as.data.frame (t(cor (matGS[,'SOX9'], hubsCell_mat)))
# sox9_cor$gene = hubs_obj$hubsCollapsed$gene[match (rownames(sox9_cor), hubs_obj$hubs_id)]
# sox9_cor = sox9_cor[order (-sox9_cor$V1),, drop=F]
# sox9_cor$region = as.character (hubs_obj$hubsCollapsed)[match (rownames(sox9_cor), hubs_obj$hubs_id)]
# head (sox9_cor,20)

TF_hub_cor = lapply (unique(KNNs_df$group2)[unique(KNNs_df$group2) != 'normal_pleura'], function(x) cor (matGS[grepl (paste0(x,'#'), rownames(matGS)),'SOX9'], hubsCell_mat[grepl (paste0(x,'#'), rownames(matGS)),], method = 'spearman'))
TF_hub_cor = lapply (TF_hub_cor, function(x) {
 as.data.frame (t(x))
  #x = x[order (-x$V1),drop=F,]
  })
TF_hub_cor = do.call (cbind, TF_hub_cor)
TF_hub_cor_means = as.data.frame (rowMeans (TF_hub_cor))
TF_hub_cor_means = na.omit (TF_hub_cor_means)
TF_hub_cor_means = TF_hub_cor_means[order (-TF_hub_cor_means[,1]),, drop=F]
colnames (TF_hub_cor_means) = 'cor'

TF_hub_cor_means$gene = hubs_obj$hubsCollapsed$gene[match (rownames(TF_hub_cor_means), hubs_obj$hubs_id)]
TF_hub_cor_means$region = as.character (hubs_obj$hubsCollapsed)[match (rownames(TF_hub_cor_means), hubs_obj$hubs_id)]

TF_hub_cor = TF_hub_cor[] 

hubs_obj$hubsCollapsed[5118]


pdf ('Plots/hubs_barcodes.pdf')
plot (matGS[,'SOX9'], hubsCell_mat[,'HUB95'])
dev.off()


### Correlate hubs to imputed TF gene score sorting cells by TF gene score ####
archp = addImputeWeights (archp)
celltype_markers = c('WT1','CALB2','GATA4','HP','SOX9','MESP1','SOX6','TWIST1','SNAI2')
seGS <- getMatrixFromProject (archp)
seGS = seGS[rowData (seGS)$name %in% celltype_markers,]
matGS <- imputeMatrix (assay(seGS), getImputeWeights(archp))
rownames (matGS) = celltype_markers

hubsCell_mat = hubsCell_mat[,colnames(matGS)]
hubsCell_mat = log2(t(hubsCell_mat)+1)

TF = 'SOX9'


hubs_mat = hubsCell_mat
hub_to_TF = function(TF= NULL, genescore_mat = NULL, hubs_mat = NULL, slice_size = 30)
  {
   genescore_TF = genescore_mat[TF, ] 
   genescore_TF = genescore_TF[order(genescore_TF)]
   hubs_mat_sorted = as.data.frame (hubs_mat[names(genescore_TF),])
   hubs_mat_sorted$TF = genescore_TF
   slice_number = floor (nrow (hubs_mat_sorted) / slice_size)
   slices = rep (seq(slice_number), each = slice_size)
   hubs_mat_sorted_adj = hubs_mat_sorted[1:length(slices),]
   hubs_mat_sorted_adj$slice = slices
   hubs_mat_agg = aggregate (.~ slice, data = hubs_mat_sorted_adj, FUN = 'mean')
   hubs_mat_agg = hubs_mat_agg[,-1]
   return (cor (hubs_mat_agg[,1:(ncol(hubs_mat_agg)-1)], hubs_mat_agg[,ncol(hubs_mat_agg)], method='spearman'))
  }

pdf (file.path (hubs_dir, 'Plots','sanity_check_corr_hub_TF.pdf'))
plot (hubs_mat_agg[,'HUB8564'], hubs_mat_agg[,'TF'])
dev.off()

print (all (colnames(matGS) == rownames (hubsCell_mat)))

# Run TF vs hub correlation ####
archp_meta = archp@cellColData
archp_meta = archp_meta[archp_meta$Clusters != 'C14',]
archp_meta = archp_meta[archp_meta$Sample2 != 'normal_pleura',]
hubs_to_TF_list = lapply (unique (archp_meta$Sample2), function(y)
    hub_to_TF (TF = 'SNAI2', 
      genescore_mat = matGS[,rownames(archp_meta)[as.logical(archp_meta$Sample2 == y)]], 
      hubs_mat = hubsCell_mat[rownames(archp_meta)[as.logical(archp_meta$Sample2 == y)],], 
      slice_size=30))

hubs_to_TF_list = lapply (hubs_to_TF_list, function(x) {
  x = as.data.frame (x)  
  x$gene = hubs_obj$hubsCollapsed$gene[match (rownames(x), hubs_obj$hubs_id)]
  x$region = as.character (hubs_obj$hubsCollapsed)[match (rownames(x), hubs_obj$hubs_id)]
  x
})
names (hubs_to_TF_list) = unique(archp_meta$Sample2)
hubs_to_TF_list_ordered = lapply (hubs_to_TF_list, function (x) {x = x[order(-x[,1]),]; x})
names (hubs_to_TF_list_ordered) = unique(archp_meta$Sample2)
head (hubs_to_TF_list_ordered[['P5']],10)

high_sample = c('P4','P5','P11','P12','P10')
hubs_to_TF_list_fl = hubs_to_TF_list[high_sample]
hubs_to_TF_median = as.data.frame (apply (do.call(cbind, lapply(hubs_to_TF_list_fl, function(x) x[,'V1',drop=F])), 1, median))
hubs_to_TF_median = hubs_to_TF_median[order(-hubs_to_TF_median[,1]),,drop=F]

hubs_obj$hubsCollapsed[match (head (rownames(hubs_to_TF_median),50), hubs_obj$hubs_id)]$gene

write.csv (hubs_obj$hubsCollapsed[match (head (rownames(hubs_to_TF_median),200), hubs_obj$hubs_id)]$gene, 
  file.path (hubs_dir, paste0('hubs_correlated_',TF,'.csv')))


### Pull hubs from genes correlated to TF from the scRNA ####
tf_cor_genes = readRDS ('../scrna/TF_cor_genes_per_sample.rds')

# check hubs in top 300 correlated genes per TF
TF = 'SOX9'
cor_genes = tf_cor_genes[[TF]]
cor_genes = apply (cor_genes,1,median)
cor_genes = cor_genes[order(-cor_genes)]
hubs_with_genes = sapply(head(names(cor_genes),50), 
    function(y) which(grepl(y, hubs_obj$hubsCollapsed$gene)))

#hubs_P1 = res2$feature[res2$logFC > 0.3 & res2$avgExpr > 0]
#hub_genes = sapply (names (head (cor_genes, 100)), function(x) as.character(hubs_obj[hubs_obj$hubs_id %in% hubs_P1]$hubsCollapsed[grepl (x, hubs_obj[hubs_obj$hubs_id %in% hubs_P1]$hubsCollapsed$gene)]))
hub_genes = hub_genes[sapply(hub_genes, function(x) length(x) > 0)]


  metaGroupName = 'Sample2'
  celltype_markers = c('HIC1','CD44','BHLHE40', 'SERPINE1')
  #celltype_markers = c('WT1','CALB2','GATA4','MSLN','KRT5','KRT18','ITLN1','HP','SOX9')
  meso_markers <- plotBrowserTrack(
      ArchRProj = archp, 
      groupBy = metaGroupName, 
      geneSymbol = celltype_markers,
      #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
      #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
      upstream = 150000,
      downstream = 150000,
      loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
      #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
      #loops = getCoAccessibility (archp, corCutOff = 0.3,
      #  returnLoops = TRUE),
      useGroups= NULL
  )
  plotPDF (meso_markers, ArchRProj = archp, width=14, name ='MPM_markers_coveragePlots.pdf')


### Check peaks around HOXB13 ####
tf_match = getMatches (archp)
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)
mega_hubs = GRanges ('chr17:48727159-48731662')
mega_hubs_peaks = bg_peakSet[queryHits(findOverlaps(bg_peakSet, mega_hubs))]
#tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
mega_hubs_TF =  hyperMotif (
  selected_peaks = mega_hubs_peaks, 
  motifmatch = tf_match)

head (mega_hubs_TF, 40)







# Compare peaks overlap between scatac celltypes and meso tumors / normal ####

# # Compute differential hub accessibility DHA for each tumor ####
library (presto)
res = wilcoxauc (log2(hubsCell_mat+1), archp$SampleP11)
res_flt = lapply (split (res, res$group), function(x) x[x$logFC > 1, ])
ps = getPeakSet (archp)
res_peaks = lapply (res_flt, function(x) ps[queryHits (findOverlaps(ps, hubs_obj$hubsCollapsed[which(hubs_obj$hubs_id %in% x$feature)]))])

# Compare peaks overlap between scatac celltypes and meso tumors / normal ####
projects = c('yang_kidney','Tsankov_lung','rawlins_fetal_lung','JShendure','greenleaf_colon','greenleaf_brain','bingren_pan')
projects_peaks = lapply (seq_along(projects), function(x) {
  bed_files = list.files (file.path('..','all_tissues_ArchR',projects[x],'PeakCalls'), pattern = '.rds')
  grlist = lapply (seq_along(bed_files), 
    function(y) readRDS (file.path('..','all_tissues_ArchR',projects[x],'PeakCalls',bed_files[y])))
names (grlist) = paste0(projects[x], '_', sapply (bed_files, function(z) unlist(strsplit (z, '-'))[1]))
grlist
})
projects_peaks = unlist (projects_peaks, recursive=F)
projects_peaks2 = GRangesList (projects_peaks)

meso_peaks = res_peaks
peaks_ov_mat = sapply (meso_peaks, function(x) sapply (projects_peaks2, function(y) sum(countOverlaps (x, y)) / min (c(length(x), length(y)))))

# prop_sample_df[is.na(prop_sample_df)] = 0
ht = Heatmap (
  # prop_sample_df, 
  peaks_ov_mat,
  col = rev(palette_deviation), 
  row_names_gp= gpar (fontsize=6), 
  column_names_gp= gpar (fontsize=6), 
  column_names_rot = 45)
pdf (file.path ('Plots', 'scatacDatasets_overlap_DAH_sampleP11_heatmap.pdf'),width = 6,height=44)
ht
dev.off()

# Plot only fetal mesothelium ####
ht = Heatmap (
  # prop_sample_df, 
  t(scale(t(peaks_ov_mat[c('rawlins_fetal_lung_Earlymeso',
    'rawlins_fetal_lung_Mid_latemeso',
    'bingren_pan_Mesothelial.Cell'),order(-peaks_ov_mat['rawlins_fetal_lung_Earlymeso',])]))),
  cluster_rows=F,
  cluster_columns=F,
  col = rev(palette_deviation), 
  row_names_gp= gpar (fontsize=6), 
  column_names_gp= gpar (fontsize=6),
  name = 'overlap'#, 
  #column_names_rot = 45
  )
pdf (file.path ('Plots', 'scatacDatasets_overlap_DAH_sampleP11_only_meso_heatmap.pdf'),width = 3.6,height=1.2)
ht
dev.off()






# Overlap hubs with CNV data from TCGA to prioritize oncogenic hubs ####
source ('../../git_repo/tumor_analysis/compile_TCGA_CNV_data.R')

cnv_hubs = matrix (ncol = length(meso_CNV_gr_hg38), nrow = length(hubs_obj$hubsCollapsed))
rownames (cnv_hubs) = as.character(hubs_obj$hubsCollapsed)

gain = .25
loss = -0.25
for (pat in seq_along(meso_CNV_gr_hg38))
  {
  cnv_hubs_ov = findOverlaps (hubs_obj$hubsCollapsed, meso_CNV_gr_hg38[[pat]], select='first')
  cnv_hubs[,pat] = meso_CNV_gr_hg38[[pat]]$Segment_Mean[cnv_hubs_ov]
  cnv_hubs[,pat][cnv_hubs[,pat] >= gain] = 1
  cnv_hubs[,pat][cnv_hubs[,pat] <= loss] = -1
  cnv_hubs[,pat][cnv_hubs[,pat] > loss & cnv_hubs[,pat] < gain] = 0
  }
cnv_hubs_sum = rowMeans (cnv_hubs)
cnv_hubs_df = data.frame (hub = factor (names(cnv_hubs_sum), levels = names(cnv_hubs_sum)[order (-cnv_hubs_sum)]), cnv = cnv_hubs_sum)
rownames (cnv_hubs_df) = hubs_obj$hubs_id[match(cnv_hubs_df$hub, as.character(hubs_obj$hubsCollapsed))]

bp = ggplot (cnv_hubs_df, aes (x = hub, y = cnv)) + 
geom_bar (stat= 'identity')

pdf (file.path('Plots','hubs_overlapping_CNV.pdf'))
bp
dev.off()

# Plot CNV vs DAH ####
cnv_hubs_df = cbind (cnv_hubs_df, res[match(rownames(cnv_hubs_df), res$feature),])
cnv_hubs_df$avgExpr_sign = cnv_hubs_df$avgExpr * sign (cnv_hubs_df$logFC)
sp = ggplot (cnv_hubs_df, aes (x = cnv, y = logFC)) + 
geom_point ()

pdf (file.path('Plots','hubs_overlapping_CNV_DAH.pdf'))
sp
dev.off()

cnv_hubs_df$cnv_sign = ifelse (cnv_hubs_df$cnv > 0, 'Amp','Del')
cnv_hubs_df2 = cnv_hubs_df[!is.na (cnv_hubs_df$cnv),]

library (rstatix)
library (ggpubr)
bp = ggplot (cnv_hubs_df2, aes (x = cnv_sign, y = logFC)) + 
geom_boxplot ()
stat.test = cnv_hubs_df2 %>%
t_test(reformulate ('cnv_sign', 'logFC')) %>%
adjust_pvalue (method = "fdr") %>%
add_significance ()
stat.test = stat.test %>% add_xy_position (x = 'cnv_sign', step.increase=.4)
bp = bp + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
bracket.nudge.y = 1, hide.ns = T,
label = "p.adj.signif") + NoLegend() + 
theme_classic()


pdf (file.path('Plots','hubs_overlapping_CNV_DAH_boxplot.pdf'), width=2, height=4)
bp
dev.off()


logfcThreshold = 1
pvalAdjTrheshold = 0.05
cnv_hubs_df2$sig = ifelse (abs(cnv_hubs_df2$logFC) > logfcThreshold & cnv_hubs_df2$padj < pvalAdjTrheshold, 1,0)
cnv_hubs_df2$sig = cnv_hubs_df2$sig * sign (cnv_hubs_df2$logFC)
cnv_hubs_df2$sig[cnv_hubs_df2$sig == -1] = 'N'
cnv_hubs_df2$sig[cnv_hubs_df2$sig == 1] = 'T'

cnv_hubs_df2$cnv_sign = ifelse (abs(cnv_hubs_df2$logFC) > logfcThreshold & cnv_hubs_df2$padj < pvalAdjTrheshold, 1,0)
cnv_hubs_df2$cnv_sign = cnv_hubs_df2$cnv_sign * sign (cnv_hubs_df2$cnv)
cnv_hubs_df2$cnv_sign[cnv_hubs_df2$cnv_sign == -1] = 'Del'
cnv_hubs_df2$cnv_sign[cnv_hubs_df2$cnv_sign == 1] = 'Amp'

res_filtered = cnv_hubs_df2[abs(cnv_hubs_df2$logFC) > logfcThreshold & cnv_hubs_df2$padj < pvalAdjTrheshold,]
res_filtered = head (res_filtered$feature[order (-abs(res_filtered$cnv))], 5)
cnv_hubs_df2$labels = ''
cnv_hubs_df2$labels[match (res_filtered, cnv_hubs_df2$feature)] = res_filtered
vp = ggplot (cnv_hubs_df2, aes(x=logFC, y= -log10(padj))) +
    geom_point(shape=19, aes (color = sig), alpha=.5) +
    geom_vline(xintercept = logfcThreshold, linetype="dotted", 
                color = "grey20", size=1) +
    geom_vline(xintercept = -logfcThreshold, linetype="dotted", 
                color = "grey20", size=1) +
    geom_hline(yintercept = -log10(pvalAdjTrheshold), linetype="dotted", 
                color = "grey20", size=1) + 
    geom_text_repel (
      size=2, 
      data = cnv_hubs_df2, 
      aes(label = labels),
      segment.size=.2,
      max.overlaps = 10000) + 
    ggtitle ('Hubs differential accessibility') +
    #geom_label_repel (size=2,max.overlaps=10000, data = deg2_cl, aes(label = show_genes), color='red') + 
    scale_color_manual (values = c("0"='grey77',"N"='green',"T"='red')) + 
    scale_fill_manual (values = c("0"='grey77',"Del"='green',"Amp"='red')) + 
    theme_light()

pdf (file.path ('Plots', 'HUB_CNV_volcano.pdf'),height=3,width=3)
vp
dev.off()

res_filtered = cnv_hubs_df2[cnv_hubs_df2$logFC > logfcThreshold & cnv_hubs_df2$padj < pvalAdjTrheshold,]
res_filtered = head (res_filtered$feature[order (-res_filtered$cnv)], 30)


## Generate circos plot for showing hubs on recurrent CNVs ####
# source script to load TCGA_CNV data
source (file.path('..','..','git_repo','tumor_analysis','compile_TCGA_CNV_data.R'))

# Add labels of hubs in top amplified regions and differentially accessed in tumors ####
cnv_hubs_df2 = cnv_hubs_df2[cnv_hubs_df2$avgExpr_sign > 1, ]
cnv_hubs_df2 = cnv_hubs_df2[order(-cnv_hubs_df2$cnv),]
head (cnv_hubs_df2)
cnv_hubs = rbind (head (cnv_hubs_df2,20), tail (cnv_hubs_df2,20))
cnv_hubs$seqnames = sapply (as.character(cnv_hubs$hub), function(x) unlist(strsplit(x, '\\:'))[1])
cnv_hubs$start = sapply (as.character(cnv_hubs$hub), function(x) unlist(strsplit(x, '\\-'))[1])
cnv_hubs$start = as.numeric(sapply (as.character(cnv_hubs$start), function(x) unlist(strsplit(x, '\\:'))[2]))
cnv_hubs$end = as.numeric(sapply (as.character(cnv_hubs$hub), function(x) unlist(strsplit(x, '\\-'))[2]))
cnv_hubs = cnv_hubs[,c('seqnames','start','end','feature')]

# hubs track
cnv_hubs_df_track = cnv_hubs_df
cnv_hubs_df_track$seqnames = sapply (as.character(cnv_hubs_df_track$hub), function(x) unlist(strsplit(x, '\\:'))[1])
cnv_hubs_df_track$start = sapply (as.character(cnv_hubs_df_track$hub), function(x) unlist(strsplit(x, '\\-'))[1])
cnv_hubs_df_track$start = as.numeric(sapply (as.character(cnv_hubs_df_track$start), function(x) unlist(strsplit(x, '\\:'))[2]))
cnv_hubs_df_track$end = as.numeric(sapply (as.character(cnv_hubs_df_track$hub), function(x) unlist(strsplit(x, '\\-'))[2]))
cnv_hubs_df_track = cnv_hubs_df_track[,c('seqnames','start','end','logFC')]

col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))


pdf (file.path ('Plots','CNV_hubs_circos.pdf'))
circos.initializeWithIdeogram (species = "hg38", chromosome.index= paste0('chr',1:22))
circos.genomicTrack (cnv_mat_avg,
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, bg.border = rep(0, 22),
            col = ifelse(value[[1]] > 0, "red", "blue"), border=NA, ...)
        #circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
})
# circos.genomicTrack (cnv_hubs_df_track,
#     panel.fun = function(region, value, ...) {
#         circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
#             col = ifelse(value[[1]] > 0, "red", "blue"), border=NA, ...)
#         #circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
# })
# circos.genomicHeatmap(cnv_hubs_df_track, col = col_fun, side = "inside", border = "white")
circos.genomicLabels(cnv_hubs, labels.column = 4, side = "inside",
    col = 'grey22', line_col = 'grey22',padding = 0.6, cex=0.5)
dev.off()

pdf (file.path ('Plots','CNV_megahubs_circos.pdf'))
circos.initializeWithIdeogram (species = "hg38", chromosome.index= paste0('chr',18))
circos.genomicTrack(cnv_mat_avg,
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
            col = ifelse(value[[1]] > 0, "red", "blue"), border=NA, ...)
        #circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
})
circos.genomicTrack (cnv_hubs_df_track,
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
            col = ifelse(value[[1]] > 0, "red", "blue"), border=NA, ...)
        #circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
})
# circos.genomicHeatmap(cnv_hubs_df_track, col = col_fun, side = "inside", border = "white")
# circos.genomicLabels(cnv_hubs, labels.column = 4, side = "inside",
#     col = 'grey22', line_col = 'grey22',padding = 0.6, cex=0.5)
dev.off()


# Map point mutations from Waddel and MESOMICS WGS to identify hyper-mutated hubs ####
snp_msm = read.table ('/ahg/regevdata/projects/ICA_Lung/Wooseung/Mesothelioma/Data/var_annovar_maf_corr_allvariants.txt', sep='\t',header=T, quote='')
snp_wdl = read.table ('/ahg/regevdata/projects/ICA_Lung/Wooseung/SuperEnhancer/Data/Mutation/Waddellhg38.bed', sep='\t',header=F, quote='')
colnames (snp_wdl) = c('chr','start','end','gene','gene2')
snp_wdl$chr = paste0('chr',snp_wdl$chr)

snp_msm_gr = makeGRangesFromDataFrame (snp_msm)
snp_wdl_gr = makeGRangesFromDataFrame (snp_wdl)

snp_msm_df = as.data.frame (snp_msm_gr)
snp_msm_df = snp_msm_df[,1:3]
snp_wdl_df = as.data.frame (snp_wdl_gr)
snp_wdl_df = snp_wdl_df[,1:3]
write.table (snp_msm_df, 'msm_mutation_snp.bed', sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table (snp_wdl_df, 'wdl_mutation_snp.bed', sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
mutated_hubs = countOverlaps  (hubs_obj$hubsCollapsed, snp_wdl_gr)
mut_density_hubs = mutated_hubs / width (hubs_obj$hubsCollapsed)
mut_density_hubs = mut_density_hubs * mutated_hubs

head (hubs_obj$hubs_id[order (-mut_density_hubs)],10)



### Find common oncogenic hubs  
# # Compute differential hub accessibility DHA of subsampled tumors vs normal ####
library (presto)
archp$status = ifelse (archp$Sample2 == 'normal_pleura','normal','tumor')
# Subsample
tumor_barcodes_sub = unlist(lapply (unique(archp$Sample2), function(x) sample (rownames(archp@cellColData)[archp$Sample2 == x], 100 , replace=TRUE)))
tumor_barcodes_sub = unique (tumor_barcodes_sub)
res = wilcoxauc (log2(hubsCell_mat[,tumor_barcodes_sub]+1), archp$status[match(tumor_barcodes_sub, rownames(archp@cellColData))])
res = res[res$group == 'tumor',]
head (res [order(-res$logFC),],10)

# Run DAH each tumor vs normal and then take intersection
dah_comp = as.list(unique (archp$Sample2))
dah_comp = lapply (dah_comp, function(x)  c(x, 'normal_pleura'))
dah_comp = dah_comp[-11]
res = lapply (dah_comp, function(x) wilcoxauc (log2(hubsCell_mat[,archp$Sample2 %in% x]+1), archp$status[archp$Sample2 %in% x]))

head (res[[1]][order(res[[1]]$padj),])

# Plot CNV vs DAH ####
cnv_hubs_df = cbind (cnv_hubs_df, res[match(rownames(cnv_hubs_df), res$feature),])
cnv_hubs_df$avgExpr_sign = cnv_hubs_df$avgExpr * sign (cnv_hubs_df$logFC)
cnv_hubs_df$label = ''
cnv_hubs_df = na.omit (cnv_hubs_df)
cnv_hubs_df$label[cnv_hubs_df$cnv > 0 & cnv_hubs_df$logFC > 2] = rownames(cnv_hubs_df)[cnv_hubs_df$cnv > 0 & cnv_hubs_df$logFC > 2]
sp = ggplot (cnv_hubs_df, aes (x = cnv, y = logFC, label=label)) + 
geom_point () + theme_light() + geom_text_repel()

pdf (file.path('Plots','hubs_overlapping_CNV_DAH.pdf'))
sp
dev.off()



# Add labels of hubs in top amplified regions and differentially accessed in tumors ####
cnv_hubs_df2 = cnv_hubs_df#[cnv_hubs_df$avgExpr_sign > 1, ]
cnv_hubs_df2 = cnv_hubs_df2[order(-cnv_hubs_df2$logFC),]
head (cnv_hubs_df2)
cnv_hubs = rbind (head (cnv_hubs_df2,20), tail (cnv_hubs_df2,20))
cnv_hubs$seqnames = sapply (as.character(cnv_hubs$hub), function(x) unlist(strsplit(x, '\\:'))[1])
cnv_hubs$start = sapply (as.character(cnv_hubs$hub), function(x) unlist(strsplit(x, '\\-'))[1])
cnv_hubs$start = as.numeric(sapply (as.character(cnv_hubs$start), function(x) unlist(strsplit(x, '\\:'))[2]))
cnv_hubs$end = as.numeric(sapply (as.character(cnv_hubs$hub), function(x) unlist(strsplit(x, '\\-'))[2]))
cnv_hubs = cnv_hubs[,c('seqnames','start','end','feature')]
cnv_hubs$logFC = res$logFC[match(rownames(cnv_hubs), res$feature)]
cnv_hubs = cnv_hubs[cnv_hubs$seqnames != 'chrX',]

# hubs track
cnv_hubs_df_track = cnv_hubs_df
cnv_hubs_df_track$seqnames = sapply (as.character(cnv_hubs_df_track$hub), function(x) unlist(strsplit(x, '\\:'))[1])
cnv_hubs_df_track$start = sapply (as.character(cnv_hubs_df_track$hub), function(x) unlist(strsplit(x, '\\-'))[1])
cnv_hubs_df_track$start = as.numeric(sapply (as.character(cnv_hubs_df_track$start), function(x) unlist(strsplit(x, '\\:'))[2]))
cnv_hubs_df_track$end = as.numeric(sapply (as.character(cnv_hubs_df_track$hub), function(x) unlist(strsplit(x, '\\-'))[2]))
cnv_hubs_df_track = cnv_hubs_df_track[,c('seqnames','start','end','logFC')]

col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
dah_pal = c('TRUE' = 'purple','FALSE' = 'darkgreen')

pdf (file.path ('Plots','CNV_hubs_circos.pdf'), 5,5)
circos.initializeWithIdeogram (species = "hg38", chromosome.index= paste0('chr',1:22))
circos.genomicTrack (cnv_mat_avg,
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, bg.border = rep(0, 22),
            col = ifelse(value[[1]] > 0, "red", "blue"), border=NA, ...)
        #circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
})
# circos.genomicTrack (cnv_hubs_df_track,
#     panel.fun = function(region, value, ...) {
#         circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
#             col = ifelse(value[[1]] > 0, "red", "blue"), border=NA, ...)
#         #circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
# })
# circos.genomicHeatmap(cnv_hubs_df_track, col = col_fun, side = "inside", border = "white")
circos.genomicLabels(cnv_hubs, labels.column = 4, side = "inside",
    col = dah_pal[as.character(cnv_hubs[[5]] > 0)], line_col = 'grey22',padding = 0.6, cex=0.5)

dev.off()

