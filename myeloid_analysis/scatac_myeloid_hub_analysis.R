# Load functions for hub detection ####
source (file.path('..','git_repo','utils','knnGen.R'))
source (file.path('..','git_repo','utils','addCoax.R'))
source (file.path('..','git_repo','utils','Hubs_finder.R'))
source (file.path('..','git_repo','utils','hubs_track.R'))



# Export bigiwg files ####
metaGroupName = 'cnmf_celltypes'
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



if (!file.exists ('myeloid_peak_regions.bed'))
  {
  peak_regions = as.data.frame (getPeakSet (archp), row.names=NULL)
  peak_regions = peak_regions[,c(1:3)]
  write.table (peak_regions, file.path('myeloid_peak_regions.bed'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
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
k= 100
metaGroupName = 'Clusters_H'

force = F
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
force=FALSE
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
metaGroupName = 'cnmf_celltypes'
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

ha = HeatmapAnnotation (size = anno_barplot(width (hubs_obj$hubsCollapsed), gp = gpar(color = "red"), height =  unit(8, "mm")))
hm = Heatmap (
  scale (t(hubsSample_mat)), 
#  top_annotation = ha, 
  column_names_gp = gpar(fontsize = 0),
  #column_km = 2,
  #row_dend_width = unit(5,'mm'),
  row_dend_side = 'left',
  clustering_distance_columns = 'pearson',
  clustering_distance_rows = 'pearson',
  col = rev(palette_hubs_accessibility),
  border=T,
  name = 'Hubs')
pdf (file.path (hubs_dir,'Plots',paste0('hubs_',metaGroupName,'_heatmap.pdf')), height=2.2)
hm
dev.off()





# Generate matrix of fragment counts of hubs x barcodes ####
force=FALSE
if (!file.exists(file.path (hubs_dir, paste0('hubs_cells_mat.rds')))| force)
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



# Compute differential hub accessibility DHA ####
library (presto)
# archp$monomac_vs_resident = ifelse (archp$cnmf_celltypes %in% 'IM', 'resident','monomac')
# metaGroupName = 'monomac_vs_resident'
metaGroupName='momac'
all (colnames(hubsCell_mat) == rownames(archp@cellColData))
compare_groups = c('momac','resident')
res = wilcoxauc (log2(hubsCell_mat[,as.character (archp@cellColData[,metaGroupName]) %in% compare_groups]+1), 
  as.character (archp@cellColData[,metaGroupName][as.character (archp@cellColData[,metaGroupName]) %in% compare_groups]))

# res_l = lapply (split (res, res$group), function(x){
#   tmp = x[x$logFC > 0,]
#   tmp = tmp[order (tmp$pval),]
#   tmp
# })
# res = res_l[[1]]
# res$gene = hubs_obj$hubsCollapsed$gene[match (res$feature, hubs_obj$hubs_id)]
# head (res,50)

res = res[res$group == 'momac',]
res$gene = hubs_obj$hubsCollapsed$gene[match (res$feature, hubs_obj$hubs_id)]

DAH_df = data.frame (Log2FC = res$logFC, 
                    #Mean = res$Mean[,1],
                    FDR = res$padj,
                    #Pval = res$Pval[,1],
                    #MeanDiff = res$MeanDiff[,1],
                    #AUC = res$AUC[,1],
                    #MeanBGD = res$MeanBGD[,1],
                    feature = res$feature)

logFCthreshold = 0.5
pValThreshold = 1e-5
DAH_df$color = 'ns'
DAH_df$label = ''
#DAH_df$label[abs (DAH_df$Log2FC) > logFCthreshold & DAH_df$FDR < pValThreshold] = rownames(DAH_df)[abs (DAH_df$Log2FC) > logFCthreshold & DAH_df$FDR < pValThreshold]
DAH_df$label[DAH_df$feature == 'HUB427'] = 'HUB427'

# Make volcano plots of DAM per each comparison
jitter = position_jitter (width = 0.01, height = 0.01)
vol_p = ggplot(DAH_df, aes(x = Log2FC, y = -log10(FDR), label = label)) +
  geom_point(aes(color = Log2FC), size = 0.2, alpha = 0.5, position = jitter) +
  geom_text_repel(
    min.segment.length = 0.2,
    box.padding = 0.2,
    size = 2,
    max.overlaps = 10000
  ) +
  ggtitle("DAH moMac vs IM") +
  xlab("LFC") + 
  ylab("-log10 adjusted p-value") +
  scale_color_gradientn(colors = rev(palette_hubs_accessibility)) +
  geom_hline(yintercept = -log10(pValThreshold), linetype = "dashed", color = "black", size = 0.2) +
  geom_vline(xintercept = -logFCthreshold, linetype = "dashed", color = "black", size = 0.2) +
  geom_vline(xintercept =  logFCthreshold, linetype = "dashed", color = "black", size = 0.2) +
  gtheme_no_rot

pdf (file.path('Plots','DAH_volcano.pdf'),width=4,height=3)
print (vol_p)
dev.off()
#palette_hubs_accessibility


## Hubs to show in IGV ####
HUB84 HUB499 HUB1324 HUB575 HUB178 HUB733 HUB429 HUB369 HUB242 HUB602


# Find top DAH ####
res = wilcoxauc (hubsCell_mat, archp$cnmf_celltypes)

res_l = lapply (split (res, res$group), function(x){
  tmp = x[x$logFC > 0,]
  tmp = tmp[order (tmp$pval),]
  head (tmp,5)
})
res_df = do.call (rbind, res_l)

DAH_df = hubsSample_mat[res_df$feature, unique(res_df$group)]
rownames (DAH_df) = paste0(rownames(DAH_df), ':',hubs_obj$hubsCollapsed$gene[match(res_df$feature, hubs_obj$hubs_id)])

hm = Heatmap (
  t(scale (t(DAH_df))), 
#  top_annotation = ha, 
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8),
  column_names_rot=45,
  #column_km = 2,
  #row_dend_width = unit(5,'mm'),
  row_dend_side = 'left',
  clustering_distance_columns = 'pearson',
  clustering_distance_rows = 'pearson',
  col = rev(palette_hubs_accessibility),
  cluster_columns = F,
  cluster_rows = F,
  border=T,
  name = 'Hubs')
pdf (file.path (hubs_dir,'Plots',paste0('hubs_DAH_heatmap.pdf')), height=6.2, width=4)
hm
dev.off()





metaGroupName = 'cnmf_celltypes'
TF = 'IL1A'
pdf()
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp,
    sample_levels = c('Mono','TREM2','SPP1','cDCs','IFN_CXCLs','IM'), 
    sizes = c(6, 1, 1, 1,1,1),
    groupBy = metaGroupName, 
    geneSymbol = TF,
    normMethod = "ReadsInTSS",
    scCellsMax=3000,
    plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
    hubs_regions = hubs_obj$hubsCollapsed,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 50000,
    pal = palette_myeloid,
    #ylim=c(0,0.1),
    downstream = 100000,
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
  name =paste0(paste(TF, collapse='_'),'_coveragePlots.pdf'),
  addDOC = F)

metaGroupName = 'celltype2'
top_dah = data.frame (
gene = srt@assays$RNA@data[TF,],
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

pdf (file.path ('Plots', paste0('scrna_',TF,'_boxplots.pdf')), height=4, width=4)
bp
dev.off()





### Find all peaks around exhausted TFs and correlate accessibility and RNA expression across celltype pseudobulks ####
TF = c('EREG',' IL1B')
metaGroupName = 'cnmf_celltypes'
# Get all peaks correlated with exhausted TF
library (org.Hs.eg.db)
gene_regions = genes (TxDb.Hsapiens.UCSC.hg38.knownGene)
eg_sym = toTable (org.Hs.egSYMBOL)
gene_regions$symbol = eg_sym$symbol[match(gene_regions$gene_id, eg_sym$gene_id)]
gene_regions = gene_regions[gene_regions$symbol %in% TF ]

extend_region = 250000
gene_regions_extended = GRangesList (lapply (seq_along(gene_regions), function(x) extendGR (gr = gene_regions[x], upstream = extend_region, downstream = extend_region)))

gene_regions_peaks = lapply (seq_along(gene_regions), function(x) queryHits (findOverlaps(getPeakSet(archp), gene_regions_extended[x])))
gene_regions_peaks = lapply (gene_regions_peaks, function(x) getPeakSet(archp)[x])

force = FALSE
if (!file.exists ('pMats.rds') | force)
  {
  pMats = lapply (gene_regions_peaks, function(x) getGroupSE(
    ArchRProj = archp,
    useMatrix = 'PeakMatrix',
    groupBy = metaGroupName,
    divideN = TRUE,
    scaleTo = NULL,
    threads = getArchRThreads(),
    verbose = TRUE,
    logFile = createLogFile("getGroupSE")
  ))
  pMats2 = lapply (seq_along(gene_regions_peaks), function(x) pMats[[x]][queryHits(findOverlaps(GRanges(rowData(pMats[[x]])), gene_regions_peaks[[x]]))])
  pMats2 = lapply (pMats2, function(x) {tmp = assay(x); rownames(tmp) = as.character(GRanges (rowData(x))); tmp})
  names (pMats2) = gene_regions$symbol
  saveRDS (pMats2, 'pMats2.rds')
  } else {
  pMats2 = readRDS ('pMats2.rds')  
  }


ps = log2(as.data.frame (AverageExpression (srt, features = gene_regions$symbol, group.by = 'shared_cnmf2_r_max')[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
colnames (ps)[colnames (ps) == 'Monocytes'] = 'Mono'

# Compute distance between TSS and correlated peaks ####
library (rtracklayer)
library (AnnotationDbi)
gene_ids = toTable (org.Hs.egSYMBOL)$gene_id[match (TF,toTable (org.Hs.egSYMBOL)$symbol)]
mygenes.transcripts = subset (genes(TxDb.Hsapiens.UCSC.hg38.knownGene, columns=c("tx_id", "tx_name","gene_id")), gene_id %in% gene_ids)
mygenes.tss = resize (mygenes.transcripts, width=1, fix='start')

mygenes.tss$symbol = toTable (org.Hs.egSYMBOL)$symbol[match (mygenes.tss$gene_id ,toTable (org.Hs.egSYMBOL)$gene_id)]

pMats2[[1]] = pMats2[[1]][,!colnames(pMats2[[1]]) %in% 'CC']
peak_gene_cor = lapply (gene_regions$symbol, function(x) cor (t(pMats2[[x]]), t(ps[x, colnames(pMats2[[x]])])))
peak_gene_cor = lapply (peak_gene_cor, function(x) x[order(-x[,1]),,drop=F])
names (peak_gene_cor) = gene_regions$symbol

top_e = 5
peak_gene_cor = lapply (names(peak_gene_cor), function(x) {
  peak_gene_cor[[x]] = as.data.frame (peak_gene_cor[[x]])
  gr = GRanges (rownames(peak_gene_cor[[x]]))
  dis = distanceToNearest (gr,mygenes.tss[mygenes.tss$symbol == x])@elementMetadata$distance
  pr = follow (gr,mygenes.tss[mygenes.tss$symbol == x])
  pr[is.na(pr)] = 0
  dis = ifelse (pr == 1, dis * -1, dis) 
  peak_gene_cor[[x]]$distance = dis
  peak_gene_cor[[x]]$name = ''
  peak_gene_cor[[x]]$name[1:top_e] = head (paste0('E',peak_gene_cor[[x]]$distance),top_e)
  peak_gene_cor[[x]]$max = apply(pMats2[[x]],1, max)[match(rownames(peak_gene_cor[[x]]),rownames(pMats2[[x]]))]
  colnames(peak_gene_cor[[x]]) = c('cor','distance','name','max')
  peak_gene_cor[[x]]
  })

names (peak_gene_cor) = gene_regions$symbol
vp = lapply (names(peak_gene_cor), function(x)
    {
    ggplot (peak_gene_cor[[x]], aes(x=distance, y= cor)) +
    geom_bar (stat = 'identity') +
    geom_point (aes (size = max), shape=21, color='white',alpha=.5, fill = 'red') +
    xlim(c(-extend_region,extend_region)) + 
    ggtitle (x) + 
    geom_label_repel (size=3, color='grey22', aes(label = name), segment.size=.2) + gtheme_no_rot
    })

pdf (file.path ('Plots','enhancers_distance_ext_TF3.pdf'),width=7,height=3)
vp
dev.off()


















### Call peaks on celltypes ####
metaGroupName = 'cnmf_cluster'
archp = addGroupCoverages (
  ArchRProj = archp, 
  groupBy = metaGroupName,  
  force = FALSE,
  minCells= 20, # I think this should be set corresponding to the smallest cluster in the group or lower
  maxCells = 500,
  minReplicates = 2,
  sampleRatio = 0.8,
  useLabels = TRUE)

archp = addReproduciblePeakSet (
    archp,
    groupBy= metaGroupName,
    peakMethod = 'Macs2',
    reproducibility = "2",
    maxPeaks = 500000, 
    minCells=20,
    force =TRUE) # I think this should be set corresponding to the smallest cluster in the group or lower
archp = addPeakMatrix (archp)
  
### TF Enrichment in hubs in each cell type ####
tf_enr_l = list()
for (ct in unique(archp@cellColData[,metaGroupName]))
{
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




### Order cells by FOXP3 expression and find hubs that correlate with it ####
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

p = plotTrajectory(archp, trajectory = "Treg_maturation", colorBy = "cellColData", name = "Treg_maturation", embedding = 'UMAP_H')
pdf (file.path ('Plots','Treg_maturation_trajectory.pdf'))
p
dev.off()

trajMM  <- getTrajectory (ArchRProj = archp, name = "Treg_maturation", useMatrix = "MotifMatrix", log2Norm = FALSE)
trajGS  <- getTrajectory (ArchRProj = archp, name = "Treg_maturation", useMatrix = "GeneScoreMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap (trajMM, pal = palette_deviation)
p2 <- plotTrajectoryHeatmap (trajGS, pal = palette_expression)

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
ha2 = rowAnnotation (foo = anno_mark(at = seq_along(top_cor_hubs_labels), 
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

# Select HUBs to show on IGV ####
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
  hubsCell_mat = t(t(hubsCell_mat) * (10^6 / archp$nFrags)) # scale
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




