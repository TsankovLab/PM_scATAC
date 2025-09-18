# Differential Peaks in cells positive for km2 vs rest ####
# Find DAP ####
#force = FALSE

archp = addClusters (input = archp, resolution = 1.5,
  reducedDims = "IterativeLSI", name = 'Clusters2',
  maxClusters = 100,
  force = TRUE)

archp = addClusters (input = archp, resolution = 3,
  reducedDims = "IterativeLSI", name = 'Clusters3',
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

# # ## Add column on DAM heatmap showing if TF is pioneer or not from chrombpnet ####
# # ## Show barplots of top TF occurrence using finemo chrombpnet outputs ####

# # ### Compare TF expression from scRNA and inferred by chrombpnet per cell type ####
# # library (httr)
# # library (XML)
# # library (igraph)
# #BiocManager::install("universalmotif")
# library ('universalmotif')

# metaGroupName = 'Clusters2'
# if (!any (ls() == 'mSE')) mSE = fetch_mat (archp_P23, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData(mSE)$name

# chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet'
# # metaGroupName = 'celltype_lv1'
# # celltypes = unique (archp_P23@cellColData[,metaGroupName])

# # tf_database = read_meme('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/HOCOMOCO_db/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme', skip = 0, readsites = FALSE, readsites.meta = FALSE)
# # tf_database = unique(unlist(lapply(tf_database, function(x) unlist(strsplit(x@name,'_'))[1])))

# # list.files (file.path(chromBPdir, celltypes[3],'no_bias_model'))
# chrombpnet_counts = list()
# celltypes = c('SOX9_low_P23','SOX9_high_P23')
# for (celltype in celltypes)
#   {
#   message (paste0('reading finemo output for ', celltype))  
#   chrombpnet_counts[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
#   }


# ### Make barplots of most abundant TFs identified in inflamed and non-inflamed cells
# bp_df = data.frame (
#   Freq = c(proportions(head(table (chrombpnet_counts[[1]]$V4)[order (-table (chrombpnet_counts[[1]]$V4))],5)),
# proportions(head (table (chrombpnet_counts[[2]]$V4)[order (-table (chrombpnet_counts[[2]]$V4))],5))),
#   TF = names (c(head(table (chrombpnet_counts[[1]]$V4)[order (-table (chrombpnet_counts[[1]]$V4))],5),
# head (table (chrombpnet_counts[[2]]$V4)[order (-table (chrombpnet_counts[[2]]$V4))],5))),
#   type = c(rep(celltypes[[1]],5), rep(celltypes[[2]],5)))


# library(dplyr)

# df_ordered <- bp_df %>%
#   group_by(type) %>%
#   arrange(desc(Freq), .by_group = TRUE) %>%
#   ungroup()
  
# P23_motif_palette = rev(paletteer::paletteer_d("palettetown::seadra") )


# df_ordered$TF_type = paste0(df_ordered$TF, df_ordered$type)
# df_ordered$TF_type = factor (df_ordered$TF_type, levels = unique (df_ordered$TF_type))
# bp = ggplot (df_ordered, aes (x = type, y = Freq, fill = TF_type)) +
#   geom_bar (stat = 'identity', position = 'stack') + 
#   scale_fill_manual (values = P23_motif_palette) + gtheme

# pdf (file.path ('Plots', 'TF_abundance_P23_barplot.pdf'),4,width=4.5)
# bp
# dev.off()


# chrombpnet_profile = list()
# celltypes = c('SOX9_low_P23','SOX9_high_P23')
# for (celltype in celltypes)
#   {
#   message (paste0('reading finemo output for ', celltype))  
#   chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
#   }


# ### Make barplots of most abundant TFs identified in inflamed and non-inflamed cells
# bp_df = data.frame (
#   Freq = c(proportions(head(table (chrombpnet_profile[[1]]$V4)[order (-table (chrombpnet_profile[[1]]$V4))],5)),
# proportions(head (table (chrombpnet_profile[[2]]$V4)[order (-table (chrombpnet_profile[[2]]$V4))],5))),
#   TF = names (c(head(table (chrombpnet_profile[[1]]$V4)[order (-table (chrombpnet_profile[[1]]$V4))],5),
# head (table (chrombpnet_profile[[2]]$V4)[order (-table (chrombpnet_profile[[2]]$V4))],5))),
#   type = c(rep(celltypes[[1]],5), rep(celltypes[[2]],5)))


# library(dplyr)

# df_ordered <- bp_df %>%
#   group_by(type) %>%
#   arrange(desc(Freq), .by_group = TRUE) %>%
#   ungroup()
  
# P23_motif_palette = rev(paletteer::paletteer_d("palettetown::seadra") )


# df_ordered$TF_type = paste0(df_ordered$TF, df_ordered$type)
# df_ordered$TF_type = factor (df_ordered$TF_type, levels = unique (df_ordered$TF_type))
# bp = ggplot (df_ordered, aes (x = type, y = Freq, fill = TF_type)) +
#   geom_bar (stat = 'identity', position = 'stack') + 
#   scale_fill_manual (values = P23_motif_palette) + gtheme

# pdf (file.path ('Plots', 'TF_abundance_P23_profile_barplot.pdf'),4,width=4.5)
# bp
# dev.off()





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
tf_markers = c('SOX9','SNAI2','RUNX2','SOX6','JUND','FOS')
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
    name = c('SOX9','SNAI2','RUNX2','RUNX1','NFIC','SOX6'), 
    embedding = "UMAP",
    pal = rev(palette_genescore),
    imputeWeights = getImputeWeights(archp_P23)
)

TF_p3 = plotEmbedding(
    ArchRProj = archp_P23, 
    colorBy = "GeneScoreMatrix", 
    name = c('SOX9','SNAI2','RUNX2','RUNX1','NFIC','SOX6'), 
    embedding = "UMAP",
    pal = rev(palette_genescore),
    imputeWeights = NULL
)

umap_p1 = plotEmbedding (ArchRProj = archp_P23, labelMeans = T, 
  colorBy = "cellColData", name = "Clusters2", 
  #pal = palette_sample,
   embedding = "UMAP")
umap_p2 = plotEmbedding (ArchRProj = archp_P23, labelMeans = T, 
  colorBy = "cellColData", name = "Clusters3", 
  #pal = palette_sample,
   embedding = "UMAP")

dev.off()

pdf (file.path('Plots','SOX9_dev_P23_fplot.pdf'), width = 14, height = 12)
wrap_plots (TF_p)
wrap_plots(TF_p2)
wrap_plots (TF_p3)
wrap_plots(umap_p1, umap_p2)
dev.off()








# Import chromBPnet finemo motifs ####
#archp_P1 = archp[archp$Clusters %in% c('C2') & archp$Sample == 'P1']
library ('universalmotif')

ps = getPeakSet (archp)

chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet'

chrombpnet_counts = list()
celltypes = c('SOX9_low_P23','SOX9_high_P23')#,'SOX9_high_P1')
#celltypes = c('SOX9_high_P1')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_counts[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
  gr = makeGRangesFromDataFrame (chrombpnet_counts[[celltype]], keep.extra.columns=T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  chrombpnet_counts[[celltype]]$peak_type = ps$peakType[findOverlaps(gr, ps, select='first')]
  }


chrombpnet_profile = list()
#celltypes = c('Mesothelium','Malignant')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  gr = makeGRangesFromDataFrame (chrombpnet_profile[[celltype]], keep.extra.columns=T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  chrombpnet_profile[[celltype]]$peak_type = ps$peakType[findOverlaps(gr, ps, select='first')]
  }

chrombpnet_profile = lapply (chrombpnet_profile, function(x) x[x$V5 != 'NaN_NaN_NaN',])


top_n <- 5
n <- length(chrombpnet_counts )

bp_list <- lapply(seq_len(n), function(i) {
  tbl <- table(chrombpnet_counts [[i]]$V4)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_counts[[i]]$V5[chrombpnet_counts [[i]]$V4 == tf][1]
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
bp_df$type = factor (bp_df$type, levels = celltypes)
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF)) ) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle(paste0("Top ",top_n," TFs (counts head)")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

  

pdf (file.path ('Plots', 'TF_abundance_counts_barplot.pdf'),6,width=6.5)
bp
dev.off()


top_n <- 5
n <- length(chrombpnet_profile)

bp_list <- lapply(seq_len(n), function(i) {
  tbl <- table(chrombpnet_profile[[i]]$V4)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_profile[[i]]$V5[chrombpnet_profile[[i]]$V4 == tf][1]
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
bp_df$type = factor (bp_df$type, levels = celltypes)
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF))) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle(paste0("Top ",top_n," TFs (profile head)")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

  

pdf (file.path ('Plots', 'TF_abundance_profile_barplot.pdf'),6,width=6.5)
bp
dev.off()


### Check how many genes are associated with SNA2 in the counts head in SOX9 high ####
snai_hits = chrombpnet_counts[['SOX9_high_P23']]
snai_hits = snai_hits[grep ('SNAI2', snai_hits$V4), ]
colnames (snai_hits) = c('chr','start','end','motif','type','peak_type')
ps = getPeakSet (archp)
snai_hits = makeGRangesFromDataFrame (snai_hits)
ps_snai = ps[unique (queryHits (findOverlaps (ps, snai_hits)))]
length (unique(ps_snai$nearestGene))









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
  gr = makeGRangesFromDataFrame (chrombpnet_counts[[celltype]], keep.extra.columns=T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  chrombpnet_counts[[celltype]]$peak_type = ps$peakType[findOverlaps(gr, ps, select='first')]
  }


chrombpnet_profile = list()
#celltypes = c('Mesothelium','Malignant')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  gr = makeGRangesFromDataFrame (chrombpnet_profile[[celltype]], keep.extra.columns=T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  chrombpnet_profile[[celltype]]$peak_type = ps$peakType[findOverlaps(gr, ps, select='first')]
  }

chrombpnet_profile = lapply (chrombpnet_profile, function(x) x[x$V5 != 'NaN_NaN_NaN',])


top_n <- 5
n <- length(chrombpnet_counts )

bp_list <- lapply(seq_len(n), function(i) {
  tbl <- table(chrombpnet_counts [[i]]$V4)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_counts[[i]]$V5[chrombpnet_counts [[i]]$V4 == tf][1]
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
bp_df$type = factor (bp_df$type, levels = celltypes)
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF)) ) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle(paste0("Top ",top_n," TFs (counts head)")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

  

pdf (file.path ('Plots', 'TF_abundance_counts_P1_barplot.pdf'),6,width=5.5)
bp
dev.off()


top_n <- 5
n <- length(chrombpnet_profile)

bp_list <- lapply(seq_len(n), function(i) {
  tbl <- table(chrombpnet_profile[[i]]$V4)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_profile[[i]]$V5[chrombpnet_profile[[i]]$V4 == tf][1]
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
bp_df$type = factor (bp_df$type, levels = celltypes)
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF))) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle(paste0("Top ",top_n," TFs (profile head)")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

  

pdf (file.path ('Plots', 'TF_abundance_profile_P1_barplot.pdf'),6,width=5.5)
bp
dev.off()




### Check footprints #####
clusters_to_keep = c('C10','C4','C5','C6','C7','C8','C9')
archp_P23 = archp_P23[archp_P23$Clusters2 %in% clusters_to_keep]
metaGroupName='Clusters2'
# archp2 = archp
# archp = archp[archp$sarc_score != 'mid']
archp <- addGroupCoverages (ArchRProj = archp_P23, groupBy = metaGroupName)
motifPositions <- getPositions (archp_P23)

motifs <- c("SOX9", "SOX6",'RUNX2','SNAI2','FOS')
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

#markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
#markerMotifs


archp_P23 = addGroupCoverages (
  ArchRProj = archp_P23, 
  groupBy = metaGroupName,  
  force = TRUE,
  minCells= 20, # I think this should be set corresponding to the smallest cluster in the group or lower
  maxCells = 500,
  minReplicates = 2,
  sampleRatio = 0.8,
  useLabels = TRUE)

#sams = c('P10','P12',
#sams = 'P23'#),'P4','P5')

#for (sam in sams)
#  {
  #peaks_sample = readRDS (file.path ('PeakCalls',paste0(sam,'-reproduciblePeaks.gr.rds')))
  #motifPositions_sample = lapply (motifPositions, function(x) x[queryHits(findOverlaps(x, peaks_sample))])  
  seFoot <- getFootprints(
    ArchRProj = archp_P23, 
    flank = 1000,
    #positions = motifPositions_sample[markerMotifs], 
    positions = motifPositions[markerMotifs], 
    groupBy = metaGroupName
  )
  
plotFootprints(
seFoot = seFoot[,c('C4._.Rep1','C9._.Rep1')],
ArchRProj = archp_P23, 
flank = 1000,
normMethod = "Subtract",
plotName = 'P23_SOX9_clusters_footprintes',
addDOC = FALSE, height=4.5, width=3,
smoothWindow = 25)
  #}





### Cross regions of SOX9 motif from chromvar and finemo to find out if there is a different motif discovered by modisco ####
### Check SOX9 in peaks in P23 ####
matches = getMatches (archp)
TFpositions = getPositions (archp)
matchesMat = assays(matches)[[1]]
colnames (matchesMat) = gsub ('_.*','', colnames (matchesMat))
colnames (matchesMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames (matchesMat))

names (TFpositions) = gsub ('_.*','', names (TFpositions))
names (TFpositions) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", names (TFpositions))
TFpositions_sox9 = TFpositions['SOX9'][[1]]

chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet'
celltype = 'SOX9_high_P23'
fine = read.table (file.path(chromBPdir, celltype, 'no_bias_model', paste0(celltype,'_finemo_counts_to_genome_browser.tsv')), sep= '\t')
fine = makeGRangesFromDataFrame (fine, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3', keep.extra.columns = T)
fine_ov = fine[queryHits (findOverlaps (fine, TFpositions_sox9))]
head (as.character(fine_ov[grep ('SNAI', fine_ov$V4)]),100)






# Correlate module scores with TFs ####
force=F
metaGroupName = 'Clusters2'
if (!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force)
  {
  DAG_list = getMarkerFeatures (
    ArchRProj = archp, 
    testMethod = "wilcoxon",
          useGroups = "C9",
          bgdGroups = "C4",
    k=100,
    binarize = FALSE,
    useMatrix = "GeneScoreMatrix",
    groupBy = metaGroupName
  #  useSeqnames="z"
  )
  saveRDS (DAG_list, paste0('DAG_',metaGroupName,'.rds'))
  } else {
  DAG_list = readRDS (paste0('DAG_',metaGroupName,'.rds'))
  }
DAG_res = do.call (cbind, (assays(DAG_list)))
colnames (DAG_res) = names(assays(DAG_list))
#DAG_res_regions = makeGRangesFromDataFrame(rowData(DAG_list)[,c(1,3,4)])
#rownames (DAG_res) = as.character(DAG_res_regions)

pdf(file.path('Plots','P23_Clusters2_MA_plot.pdf'), width=5,height=5)
pma <- markerPlot (seMarker = DAG_list, name = 'C9', cutOff = "FDR <= 0.01", plotAs = "MA")
pma
dev.off()

# Take only significant regions ####
DAG_res_sig = DAG_res[DAG_res$FDR < .01 & DAG_res$Log2FC > 0, ]
saveRDS (GRanges(rownames(DAG_res_sig)), 'SOX9_high_genes.rds')

# Run GSEA enrichment analysis #### 
#!!! To run this analysis load only ArcHR and clusterprofiler packages !!!!
library (fgsea)    
options(warn = 0)
ps = getPeakSet (archp)

gmt_annotations = c(
'h.all.v7.4.symbol.gmt',#,
'c5.bp.v7.1.symbol.gmt',
'c3.tft.v7.1.symbol.gmt'
)

gmt.file = paste0 ('../../git_repo/files/c5.bp.v7.1.symbol.gmt')
gmt.file = paste0 ('../../git_repo/files/h.all.v7.4.symbols.gmt')
# gmt.file = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_DN.v2024.1.Hs.gmt'
# gmt.file = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP.v2024.1.Hs.gmt'
# gmt.file = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/GSE24026_PD1_LIGATION_VS_CTRL_IN_ACT_TCELL_LINE_UP.v2024.1.Hs.gmt'
# csv.file = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/human/Tcell_exhastuion_genes_PMID37091230.csv'
pathways = clusterProfiler::read.gmt (gmt.file)
#pathways = pathways[grep('inflammatory', pathways$term,ignore.case=T),]
#pathways = read.csv (csv.file)
pathways = split(pathways$gene, pathways$term)
#pathways = gmtPathways (gmt.file)
DAG_list = readRDS (paste0('DAG_',metaGroupName,'.rds'))
DAG_res = do.call (cbind, (assays(DAG_list)))
colnames (DAG_res) = names(assays(DAG_list))
rownames(DAG_res) = rowData(DAG_list)$name
#DAG_res_regions = makeGRangesFromDataFrame(rowData(DAG_list)[,c(1,3,4)])
#rownames(DAG_res) = as.character(DAG_res_regions)

# peak_genes = unname(ps[queryHits (findOverlaps(ps, GRanges (rownames(DAP_res))))]$nearestGene)
# names (peak_genes) = as.character(ps)[queryHits(findOverlaps (ps,GRanges (rownames(DAP_res))))]
# peak_genes = peak_genes[rownames(DAP_res)]
#peak_genes = setNames (-log10(DAP_res$Pval) * sign (DAP_res$Log2FC), peak_genes)
ranked_genes = setNames (-log10(DAG_res$FDR) * sign(DAG_res$Log2FC), rownames(DAG_res))
ranked_genes = na.omit (ranked_genes)
#peak_genes = peak_genes[!duplicated(names(peak_genes))]
#names (peak_genes) = 
ranked_genes = ranked_genes[!duplicated(names(ranked_genes))]
ranked_genes = ranked_genes[!is.na(names(ranked_genes))]
library (BiocParallel)
BiocParallel::register(BiocParallel::SerialParam())

### fgsea throws a BiocParallel error when I load all packages including clusterProfiler...try avoiding loading packages except ArchR and fgsea
#peak_genes2 = setNames(order(peak_genes), names(peak_genes))
fgseaRes = fgsea (pathways, 
          ranked_genes,#, 
          minSize=15, 
        #  scoreType='pos',
          maxSize=1500,
          nproc=1,
          nPermSimple=100000,
          BPPARAM = NULL
          )
pvalAdjTrheshold = 0.05
top_pathways=5
#fgseaRes$padj = fgseaRes$pval
fgseaResAll_dp = dotGSEA (
  list(fgseaRes), 
  padj_threshold = pvalAdjTrheshold, 
  type = 'fgsea',
  top_pathways = top_pathways,
  cluster_rows=F,
  cluster_cols=F)

pdf (file.path ('Plots','fgsea_dotplot.pdf'), width=7, height=3)
fgseaResAll_dp
dev.off()

### Compare DEG SOX9 high vs SOX9 low to list of genes from Elaine Fuchs ####
rankby = 'pval_signedLFC'

sox9_gs = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/gene_sets/Sox9_targets_mouse.csv')
sox9_gs_W1 = sox9_gs$symbol[sox9_gs$W1 > 0]
sox9_gs_W2 = sox9_gs$symbol[sox9_gs$W2 > 0]
sox9_gs_W6 = sox9_gs$symbol[sox9_gs$W6 > 0]
sox9_gs_W12 = sox9_gs$symbol[sox9_gs$W12 > 0]
sox9_gs_class = split (sox9_gs$symbol, sox9_gs$class)
pathways = sox9_gs_class
pathways = lapply (pathways, function(x) mouseMan (x))
pathways = lapply (pathways, function(x) x$HGNC.symbol)

#fgsea_ranks = -log10 (deg2Cluster$p_val + 1e-300) * sign (deg2Cluster$avg_log2FC)
#fgsea_ranks = degCluster$avg_log2FC
fgseaRes_fuchs = fgsea (pathways, 
          ranked_genes,#, 
          minSize=15, 
        #  scoreType='pos',
          maxSize=5000,
          nproc=1,
          nPermSimple=100000,
          BPPARAM = NULL
          )
#fgseaResCol = collapsePathways (fgseaRes, stats = fgsea_ranks, pathway = pathways)
pvalAdjTrheshold = 0.05
top_pathways=5
fgseaRes_fuchs$padj = fgseaRes_fuchs$padj
fgseaResAll_dp = dotGSEA (
  list(fgseaRes_fuchs), 
  padj_threshold = pvalAdjTrheshold, 
  type = 'fgsea',
  top_pathways = top_pathways,
  cluster_rows=F,
  cluster_cols=F)

pdf (file.path ('Plots','fgsea_fuchs_dotplot.pdf'), width=4, height=3)
fgseaResAll_dp
dev.off()

#fgsea_ranks

message ('Print enrichment plots for each signficant pathway and cell type')
  
        sig_pathways = fgseaRes_fuchs$pathway[fgseaRes_fuchs$padj < pvalAdjTrheshold]
        sig_pathways = 'W12'
        sig_pathways = 'W2to6'
        #gmt.file = paste0 ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/GSEA_gs/',org,'/',names(fgseaResAll)[x])
        #pathways = gmtPathways (gmt.file)
  pdf()
  ep = lapply (sig_pathways, function(z) {
           print (plotEnrichment(pathways[[z]],
           ranked_genes) + 
           labs(title=z))
        })
  dev.off()
  pdf (file.path('Plots','SOX9_Fuchs_enrichment_plots.pdf'),width= 5, height= 3)
  print(ep)
  dev.off() 

