# Differential Peaks in cells positive for km2 vs rest ####
# Find DAP ####
#force = FALSE

archp = addClusters (input = archp, resolution = 1.5,
  reducedDims = "IterativeLSI", name = 'Clusters2',
  maxClusters = 100,
  force = TRUE)


force=T
metaGroupName = 'Clusters2'
if (!file.exists (paste0('DAP_',metaGroupName,'.rds')) | force)
  {
  DAP_list = getMarkerFeatures (
    ArchRProj = archp, 
    testMethod = "wilcoxon",
          useGroups = "C9",
          bgdGroups = "C4",
    k=100,
    binarize = FALSE,
    useMatrix = "PeakMatrix",
    groupBy = metaGroupName
  #  useSeqnames="z"
  )
  saveRDS (DAP_list, paste0('DAP_',metaGroupName,'.rds'))
  } else {
  DAP_list = readRDS (paste0('DAP_',metaGroupName,'.rds'))
  }
DAP_res = do.call (cbind, (assays(DAP_list)))
colnames (DAP_res) = names(assays(DAP_list))
DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames (DAP_res) = as.character(DAP_res_regions)

pdf(file.path('Plots','P23_Clusters2_MA_plot.pdf'), width=5,height=5)
pma <- markerPlot (seMarker = DAP_list, name = 'C9', cutOff = "FDR <= 0.01", plotAs = "MA")
pma
dev.off()

# Take only significant regions ####
DAP_res_sig = DAP_res[DAP_res$FDR < .01 & DAP_res$Log2FC > 0, ]
saveRDS (GRanges(rownames(DAP_res_sig)), 'P23_Clusters2_C9_peaks.rds')

### Perform enrichment on DAP ####
archp = addBgdPeaks (archp, force= F)
archp = addMotifAnnotations (ArchRProj = archp, 
      motifSet = "cisbp", 
      #motifSet = 'JASPAR2020',
      #name = "JASPAR2020_Motif",
      force=F)
enrichMotifs <- peakAnnoEnrichment(
    seMarker = DAP_list,
    ArchRProj = archp,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1"
  )
heatmapEM = plotEnrichHeatmap (enrichMotifs, n = 7, transpose = TRUE)
pdf (file.path ('Plots','enrich_P23_Clusters2_heatmap.pdf'))
heatmapEM
dev.off()



archp_P23 = archp[archp$Clusters2 %in% c('C9','C4')]

# ## Add column on DAM heatmap showing if TF is pioneer or not from chrombpnet ####
# ## Show barplots of top TF occurrence using finemo chrombpnet outputs ####

# ### Compare TF expression from scRNA and inferred by chrombpnet per cell type ####
# library (httr)
# library (XML)
# library (igraph)
#BiocManager::install("universalmotif")
library ('universalmotif')

metaGroupName = 'Clusters2'
if (!any (ls() == 'mSE')) mSE = fetch_mat (archp_P23, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData(mSE)$name

chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet'
# metaGroupName = 'celltype_lv1'
# celltypes = unique (archp_P23@cellColData[,metaGroupName])

# tf_database = read_meme('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/HOCOMOCO_db/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme', skip = 0, readsites = FALSE, readsites.meta = FALSE)
# tf_database = unique(unlist(lapply(tf_database, function(x) unlist(strsplit(x@name,'_'))[1])))

# list.files (file.path(chromBPdir, celltypes[3],'no_bias_model'))
chrombpnet_counts = list()
celltypes = c('SOX9_low_P23','SOX9_high_P23')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_counts[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
  }


### Make barplots of most abundant TFs identified in inflamed and non-inflamed cells
bp_df = data.frame (
  Freq = c(proportions(head(table (chrombpnet_counts[[1]]$V4)[order (-table (chrombpnet_counts[[1]]$V4))],5)),
proportions(head (table (chrombpnet_counts[[2]]$V4)[order (-table (chrombpnet_counts[[2]]$V4))],5))),
  TF = names (c(head(table (chrombpnet_counts[[1]]$V4)[order (-table (chrombpnet_counts[[1]]$V4))],5),
head (table (chrombpnet_counts[[2]]$V4)[order (-table (chrombpnet_counts[[2]]$V4))],5))),
  type = c(rep(celltypes[[1]],5), rep(celltypes[[2]],5)))


library(dplyr)

df_ordered <- bp_df %>%
  group_by(type) %>%
  arrange(desc(Freq), .by_group = TRUE) %>%
  ungroup()
  
P23_motif_palette = rev(paletteer::paletteer_d("palettetown::seadra") )


df_ordered$TF_type = paste0(df_ordered$TF, df_ordered$type)
df_ordered$TF_type = factor (df_ordered$TF_type, levels = unique (df_ordered$TF_type))
bp = ggplot (df_ordered, aes (x = type, y = Freq, fill = TF_type)) +
  geom_bar (stat = 'identity', position = 'stack') + 
  scale_fill_manual (values = P23_motif_palette) + gtheme

pdf (file.path ('Plots', 'TF_abundance_P23_barplot.pdf'),4,width=4.5)
bp
dev.off()


chrombpnet_profile = list()
celltypes = c('SOX9_low_P23','SOX9_high_P23')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  }


### Make barplots of most abundant TFs identified in inflamed and non-inflamed cells
bp_df = data.frame (
  Freq = c(proportions(head(table (chrombpnet_profile[[1]]$V4)[order (-table (chrombpnet_profile[[1]]$V4))],5)),
proportions(head (table (chrombpnet_profile[[2]]$V4)[order (-table (chrombpnet_profile[[2]]$V4))],5))),
  TF = names (c(head(table (chrombpnet_profile[[1]]$V4)[order (-table (chrombpnet_profile[[1]]$V4))],5),
head (table (chrombpnet_profile[[2]]$V4)[order (-table (chrombpnet_profile[[2]]$V4))],5))),
  type = c(rep(celltypes[[1]],5), rep(celltypes[[2]],5)))


library(dplyr)

df_ordered <- bp_df %>%
  group_by(type) %>%
  arrange(desc(Freq), .by_group = TRUE) %>%
  ungroup()
  
P23_motif_palette = rev(paletteer::paletteer_d("palettetown::seadra") )


df_ordered$TF_type = paste0(df_ordered$TF, df_ordered$type)
df_ordered$TF_type = factor (df_ordered$TF_type, levels = unique (df_ordered$TF_type))
bp = ggplot (df_ordered, aes (x = type, y = Freq, fill = TF_type)) +
  geom_bar (stat = 'identity', position = 'stack') + 
  scale_fill_manual (values = P23_motif_palette) + gtheme

pdf (file.path ('Plots', 'TF_abundance_P23_profile_barplot.pdf'),4,width=4.5)
bp
dev.off()





### Check SOX9 in peaks in P23 ####
matches = getMatches (archp)
TFpositions = getPositions (archp)
matchesMat = assays(matches)[[1]]
colnames (matchesMat) = gsub ('_.*','', colnames (matchesMat))
colnames (matchesMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames (matchesMat))

names (TFpositions) = gsub ('_.*','', names (TFpositions))
names (TFpositions) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", names (TFpositions))
TFpositions_sox9 = TFpositions['SOX9'][[1]]

P23_high_peaks = read.table ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet/MACS2_SOX9_high_P23/SOX9_high_P23_peaks_capped.narrowPeak')
colnames (P23_high_peaks) = c('chr','start','end')
p23_high_gr = makeGRangesFromDataFrame (P23_high_peaks[,1:3])
#sox9_hits = rowRanges (matches)[matchesMat[,'SOX9']]
sox9_hits = TFpositions_sox9[queryHits(findOverlaps (TFpositions_sox9, p23_high_gr))]
elementMetadata(sox9_hits)$gene = rowRanges (matches)$nearestGene[subjectHits(findOverlaps(sox9_hits,rowRanges(matches)))]

write.table (as.data.frame (sox9_hits, row.names = NULL), file = 'SOX9_motifs_positions_in_P23_high_peaks.bed', row.names=FALSE, col.names=FALSE, sep='\t', quote=F)
write.table (P23_high_peaks, file = 'P23_high_peaks.bed', row.names=FALSE, col.names=FALSE, sep='\t', quote=F)


write.table (as.data.frame (rowRanges (matches)[matchesMat[,'NFYB']], row.names = NULL), file = 'NFYB_motifs.bed', row.names=FALSE, col.names=FALSE, sep='\t', quote=F)


# check SOX9 deviation again
archp_P23 = archp[archp$Sample %in% c('P23')]
archp_P23 = archp_P23[archp_P23@embeddings$UMAP[[1]][[2]] < 0]
tf_markers = c('SOX9','SNAI2','RUNX2','RUNX1','NFIC')
markerMotifs = getFeatures (archp_P23, select = paste(tf_markers, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs = grep ("z:", markerMotifs, value = TRUE)
archp_P23 = addImputeWeights (archp_P23)

pdf()
TF_p = plotEmbedding(
    ArchRProj = archp_P23, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    pal = rev(palette_deviation),
    imputeWeights = getImputeWeights(archp_P23)
)

TF_p2 = plotEmbedding(
    ArchRProj = archp_P23, 
    colorBy = "GeneScoreMatrix", 
    name = c('SOX9','SNAI2','RUNX2','RUNX1','NFIC'), 
    embedding = "UMAP",
    pal = rev(palette_genescore),
    imputeWeights = getImputeWeights(archp_P23)
)

umap_p1 = plotEmbedding (ArchRProj = archp_P23, labelMeans = T, 
  colorBy = "cellColData", name = "Clusters2", 
  #pal = palette_sample,
   embedding = "UMAP")

dev.off()

pdf (file.path('Plots','SOX9_dev_P23_fplot.pdf'), width = 14, height = 12)
wrap_plots (TF_p)
wrap_plots(TF_p2)
dev.off()





# Import chromBPnet finemo motifs for P1 too ####
archp_P1 = archp[archp$Clusters %in% c('C2') & archp$Sample == 'P1']
library ('universalmotif')

# metaGroupName = 'Clusters2'
# if (!any (ls() == 'mSE')) mSE = fetch_mat (archp_P1, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData(mSE)$name
chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet'
# metaGroupName = 'celltype_lv1'
# celltypes = unique (archp_P23@cellColData[,metaGroupName])

# tf_database = read_meme('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/HOCOMOCO_db/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme', skip = 0, readsites = FALSE, readsites.meta = FALSE)
# tf_database = unique(unlist(lapply(tf_database, function(x) unlist(strsplit(x@name,'_'))[1])))

# list.files (file.path(chromBPdir, celltypes[3],'no_bias_model'))
chrombpnet_counts = list()
celltypes = c('SOX9_high_P1')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_counts[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
  }


### Make barplots of most abundant TFs identified in inflamed and non-inflamed cells
bp_df = data.frame (
  Freq = c(proportions(head(table (chrombpnet_counts[[1]]$V4)[order (-table (chrombpnet_counts[[1]]$V4))],5))),
  TF = names (c(head(table (chrombpnet_counts[[1]]$V4)[order (-table (chrombpnet_counts[[1]]$V4))],5))),
  type = c(rep(celltypes[[1]],5)))


library(dplyr)

df_ordered <- bp_df %>%
  group_by(type) %>%
  arrange(desc(Freq), .by_group = TRUE) %>%
  ungroup()
  
P1_motif_palette = rev(paletteer::paletteer_d("palettetown::seadra") )


df_ordered$TF_type = paste0(df_ordered$TF, df_ordered$type)
df_ordered$TF_type = factor (df_ordered$TF_type, levels = unique (df_ordered$TF_type))
bp = ggplot (df_ordered, aes (x = type, y = Freq, fill = TF_type)) +
  geom_bar (stat = 'identity', position = 'stack') + 
  scale_fill_manual (values = P1_motif_palette) + gtheme

pdf (file.path ('Plots', 'TF_abundance_P1_barplot.pdf'),4,width=4.5)
bp
dev.off()


chrombpnet_profile = list()
celltypes = c('SOX9_high_P1')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  }


### Make barplots of most abundant TFs identified in inflamed and non-inflamed cells
bp_df = data.frame (
  Freq = c(proportions(head(table (chrombpnet_profile[[1]]$V4)[order (-table (chrombpnet_profile[[1]]$V4))],5))),
  TF = names (c(head(table (chrombpnet_profile[[1]]$V4)[order (-table (chrombpnet_profile[[1]]$V4))],5))),
  type = c(rep(celltypes[[1]],5)))


library(dplyr)

df_ordered <- bp_df %>%
  group_by(type) %>%
  arrange(desc(Freq), .by_group = TRUE) %>%
  ungroup()
  
P1_motif_palette = rev(paletteer::paletteer_d("palettetown::seadra") )


df_ordered$TF_type = paste0(df_ordered$TF, df_ordered$type)
df_ordered$TF_type = factor (df_ordered$TF_type, levels = unique (df_ordered$TF_type))
bp = ggplot (df_ordered, aes (x = type, y = Freq, fill = TF_type)) +
  geom_bar (stat = 'identity', position = 'stack') + 
  scale_fill_manual (values = P1_motif_palette) + gtheme

pdf (file.path ('Plots', 'TF_abundance_P1_profile_barplot.pdf'),4,width=4.5)
bp
dev.off()
