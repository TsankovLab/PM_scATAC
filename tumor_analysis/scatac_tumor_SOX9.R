conda activate meso_scatac
R

set.seed(1234)

####### ANALYSIS of TUMOR compartment #######
projdir = 'tumor_compartment'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

# Set # of threads and genome reference ####
addArchRThreads (threads = 1)
addArchRGenome ("hg38")

archp = loadArchRProject (projdir)   

#### Create ArchR object archp_norm with scATAC-seq data from normal distal lung samples at NIH SRA database: PRJNA1044083. ####
sample_names_distal_normal_lung = c( # Include only samples from distal normal lung 
  # Normal
  'RPL_280_neg_1',
  'RPL_280_neg_2',
  'RPL_Epi_1',
  'RPL_Epi_2'
)
inputFiles = 'path/to/fragment_files_of_distal_normal_lung'
#setwd (ArrowDir)
ArrowFiles = createArrowFiles (inputFiles = inputFiles,
  sampleNames = sample_names_distal_normal_lung,
  #validBarcodes = c(meso_im, normal_im),
  minTSS = 0, #Dont set this too high because you can always increase later
  minFrags = 0, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = T,
  subThreading= FALSE
)

#projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/combined_data/normal_all_meso/'
#archp = ArchRProject (

projdir = 'normal_distal_lung'
archp_norm = ArchRProject (
    ArrowFiles = ArrowFiles, 
    outputDirectory = projdir,
    copyArrows = FALSE 
  )

normal_mesothelium_scatac_barcodes = read.csv (file.path ('..','..','git_repo','files','normal_mesothelium_scatac_barcodes.csv')) 
archp_norm = archp_norm[colnames(archp_norm@cellColData) %in% normal_mesothelium_barcodes$x] # Subset for only normal mesothelial cells
saveArchRProject (archp_norm) # Save to update ArrowFiles 
  
### Make ArchR object including PM tumor cells and normal mesothelium from Shah et al, PRJNA1044083 ####
ArrowFiles_tumor_dir = 'tumor_compartment/ArrowFiles'
ArrowFiles_tumor_dir = file.path(ArrowFiles_tumor_dir, paste0(unique(archp$Sample),'.arrow'))

ArrowFiles_normal_dir = '/path/to/normal_distal_lung_ArchR'
ArrowFiles_normal_dir = file.path(ArrowFiles_normal_dir, paste0(unique(archp_norm$Sample),'.arrow'))

ArrowFiles_combined_dir = c(ArrowFiles_tumor_dir, ArrowFiles_normal_dir)

projdir = 'tumor_compartment_with_normal'
archp = ArchRProject (
ArrowFiles = ArrowFiles_combined_dir,
outputDirectory = projdir,
copyArrows = FALSE) 
saveArchRProject (archp)  
setwd (projdir)

### Similarly, create srt_norm seurat object from Shah et al, PRJNA1044083 (use only HU37 and HU62 samples) to merge with PM tumor scRNA-seq ####
normal_mesothelium_scrna_barcodes = read.csv (file.path ('..','..','git_repo','files','normal_mesothelium_scrna_barcodes.csv'))
srt_NN = readRDS ('GSM9326197_scRNAseq_srt.rds')
srt_NN$sampleID3 = srt_NN$sampleID

srt_norm = readRDS ('path/to/seurat_distal_normal_lung.rds')
srt_norm = srt_norm[, colnmaes(srt_norm) %in% normal_mesothelium_scrna_barcodes$x]  # Subset for only normal mesothelial cells
srt_norm$sampleID3 = 'normal_pleura'
srt = merge (srt_NN, srt_norm)

srt$sampleID3[srt$sampleID3 %in% c('HU37','HU62')] = 'normal_pleura'
sarc_order = read.csv (file.path('..','scrna','cnmf20_sarcomatoid_sample_order.csv'), row.names=1)
sarc_order = sarc_order[! sarc_order$sampleID %in% c('HU37','HU62'),]
sarc_order = rbind (data.frame (sampleID = 'normal_pleura', x = -1),sarc_order) 


archp = addClusters (input = archp, resolution = 1.5,
  reducedDims = "IterativeLSI", name = 'Clusters2',
  maxClusters = 100,
  force = TRUE)
archp_NN = archp[!archp$Sample3 %in% c('normal1')]

pdf ()
umap_p1 = plotEmbedding (ArchRProj = archp_NN, labelMeans = F, 
  colorBy = "cellColData", name = "Sample3", 
  pal = palette_sample,
   embedding = "UMAP")
umap_p2 = plotEmbedding (ArchRProj = archp_NN, labelMeans = T, 
  colorBy = "cellColData", name = "Clusters",
   embedding = "UMAP")
umap_p3 = plotEmbedding (ArchRProj = archp_NN, labelMeans = T, 
  colorBy = "cellColData", name = "Clusters2",
   embedding = "UMAP")
dev.off()

pdf (file.path ('Plots','sample_clusters_umap.pdf'))
umap_p1
umap_p2
umap_p3
dev.off()


# Rank tumor samples by sarcomatoid score ####
# Add cNMF identified in scRNA to archr object ####
if (!exists('gsSE')) gsSE = fetch_mat (archp, 'GeneScore')
gsMat = assays (gsSE)[[1]]
rownames (gsMat) = rowData (gsSE)$name

cnmf_spectra_unique = readRDS (file.path('..','..','git_repo','files','malignant_cnmf_genelist_25_nfeat_5000.rds'))
cnmf_spectra_unique = lapply (cnmf_spectra_unique, function(x) head(x, top_genes)[head(x, top_genes) %in% rownames (gsMat)])

force = FALSE
if (!all (names (cnmf_spectra_unique) %in% colnames (archp@cellColData)) | force)
  {
  archp@cellColData = archp@cellColData[,!grepl ('mod',colnames(archp@cellColData), ignore.case=T)]
  archp = addModuleScore (
      ArchRProj = archp,
      useMatrix = 'GeneScoreMatrix',
      name = '',
      features = cnmf_spectra_unique,
      nBin = 25,
      nBgd = 100,
      seed = 1,
      threads = getArchRThreads(),
      logFile = createLogFile("addModuleScore")
    )
  colnames (archp@cellColData) = gsub ('^\\.','',colnames(archp@cellColData))
  }

sams = as.character(unique(archp_meta$Sample3))
tumor_sams = sams[!sams %in% c('normal1','P11_HOX')] # remove normal,low cell numbers and outlier samples

cnmf_mat = t(scale(t(as.data.frame (archp@cellColData[,names(cnmf_spectra_unique)]))))
archp_meta = as.data.frame (archp@cellColData)
archp_meta = archp_meta[archp_meta$Sample3 %in% tumor_sams,]
average_by_group <- as.data.frame(archp_meta) %>%
  group_by (Sample3) %>%
  summarise(Average = median(cNMF20))
average_by_group = average_by_group[order(-average_by_group$Average),]
sample_sarc_order = factor (archp_meta$Sample3, levels = average_by_group$Sample3)

# ----------------------------
# Define sample groups
# ----------------------------
archp_meta = as.data.frame (archp@cellColData)

sarc_module = 'cNMF20'
p_l2=list()
for (sam in tumor_sams)
  {
  archp_sub = archp[archp$Sample3 == sam]
  varfeat = 25000
  LSI_method=2
  archp_sub = addIterativeLSI (ArchRProj = archp_sub, 
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force=TRUE, LSIMethod=LSI_method,
    varFeatures = varfeat)
  archp_sub = addUMAP (ArchRProj = archp_sub, 
    reducedDims = "IterativeLSI", seed = 2,
    force = TRUE)
  archp_sub = addImputeWeights (archp_sub)
  # Dimensionality reduction and clustering
  varfeat = 25000
  LSI_method=2
  pdf()
  p_l2[[sam]] <- plotEmbedding(
      ArchRProj = archp_sub, 
      colorBy = "cellColData", 
      name = sarc_module, rastr=F,
      embedding = "UMAP",
      pal = palette_expression,
      imputeWeights = getImputeWeights(archp_sub),
      plotAs='points'
  )
  dev.off()
  }

pdf (file.path('Plots','sarcomatoid_score_feature_plots.pdf'), width = 12, height = 15)
print (wrap_plots (p_l2), ncol = 5)
dev.off()

# Make violin + boxplot aof sarcomatoid scores and barplot of cell numbers in samples 
ccomp_df = as.data.frame (archp@cellColData)
ccomp_df = ccomp_df[ccomp_df$Sample3 %in% tumor_sams,]
ccomp_df$Sample3 = factor (ccomp_df$Sample3, levels = levels (sample_sarc_order))
box = ggplot (ccomp_df, aes_string (x= 'Sample3', y= 'cNMF20')) +
  geom_violin (trim=TRUE, aes_string (fill = 'Sample3'),size=2,
    width=1,
    scale='width',
    linewidth = .2, alpha=0.7) +
  geom_boxplot (aes_string(fill = 'Sample3'),
    linewidth = .2,
    width=0.2,
    outlier.alpha = 0.2,
    outlier.size = 1,
    outlier.shape = NA,
     size=0.3, alpha=0.7
     ) + ylim (c(-0.4,.8))+
  gtheme +
  scale_fill_manual (values= palette_sample) +
  NoLegend()
  
pdf(paste0('Plots/Sarcomatoid_signatures_boxplot.pdf'),width=3,2) #width = 10, height = 11,
print (box)
dev.off()

# barplot cell abundances


### Run peak calling ####
metaGroupName = "Clusters"
force=FALSE
peak_reproducibility='1' # Set to 1 to better identify inter-tumor heterogeneity
if(!all(file.exists(file.path('PeakCalls', paste0(unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds')))) | force) source (file.path('..','..','git_repo','utils','callPeaks.R'))
  
### chromVAR analysis ####
force=TRUE
if (!all(file.exists(file.path('Annotations',
  c('Motif-Matches-In-Peaks.rds',
    'Motif-Positions-In-Peaks.rds',
    'Motif-In-Peaks-Summary.rds'))))| force)
source (file.path('..','..','git_repo','utils','chromVAR.R'))


#### Compute P2G ####
run_p2g_TF = FALSE

if (run_p2g_TF)
  {
  run_p2g = FALSE  
  if (run_p2g)
    {
    maxDist = 250000
    archp = addPeak2GeneLinks(
        ArchRProj = archp,
        useMatrix = 'GeneScoreMatrix',
        reducedDims = "IterativeLSI",
        maxDist = maxDist
    )
    }  
    
  p2g_corr = .2
  p2g = getPeak2GeneLinks(
      ArchRProj = archp,
      corCutOff = p2g_corr,
      resolution = 1,
      returnLoops = FALSE
  )



# Discover oncogenic drivers ####
if (!file.exists ('selected_TF.rds')) 
  source (file.path('..','..','git_repo','tumor_analysis','discover_oncogenic_TFs.R'))
selected_TF = readRDS ('selected_TF.rds')


# Make coexpression network for each sample using top TFs deviations ####
selected_TF = readRDS ('selected_TF.rds')
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
archp_meta = as.data.frame (archp@cellColData)

# ----------------------------
# Define sample groups
# ----------------------------

sams = as.character(unique(archp_meta$Sample3))
tumor_sams = sams[!sams %in% c('normal1','P3','P13','P11_HOX')] # remove normal,low cell numbers and outlier samples

#archp_meta = archp_meta[archp_meta$Sample3 %in% sams,]
mMat = assays (mSE)[[1]]
  rownames (mMat) = rowData (mSE)$name
mMat = t(as.matrix(scale(mMat[selected_TF,])))

all (rownames(mMat) == rownames(archp_meta))
cor_TF_l = list()
for (sam in tumor_sams)
  {
  cor_TF_l[[sam]] = cor (mMat[archp_meta$Sample3 == sam,], method = 'spearman')
  }

corTF_array <- simplify2array (cor_TF_l)
median_matrix <- apply (corTF_array, c(1, 2), median)

set.seed (123)
hr = hclust(as.dist(1-median_matrix))#, method = "average")
clusters = dendextend::cutree(hr, k = 2)

pdf()
set.seed (1234)
cor_TF_df = draw (Heatmap (median_matrix,
  row_names_gp = gpar(fontsize = 6),
  clustering_distance_rows='pearson',
  clustering_distance_columns='pearson',
  column_names_gp = gpar(fontsize = 6),
  col = palette_deviation_cor_fun, border=T))
dev.off()

pdf (file.path ('Plots','selected_TF_dev_corr_heatmaps.pdf'), width = 8,height=7)
cor_TF_df 
dev.off()

# Add pathways of TF correlated genes from scRNA ####
if (!file.exists('enrichment_pathways_TFs.rds')) source (file.path('..','..','git_repo','tumor_analysis','enrichment_cnmfs.R'))
TF_cor_sum = readRDS (file.path('enrichment_pathways_TFs.rds'))
TFrow_order = unname(unlist(row_order (cor_TF_df)))

rownames (TF_cor_sum) = gsub ('HALLMARK_','',rownames(TF_cor_sum))
rownames (TF_cor_sum) = gsub ('_', ' ',rownames (TF_cor_sum))
TF_cor_sum = TF_cor_sum[apply (TF_cor_sum, 1, function(x) any(x > 1)),]
pdf()
hm = draw (Heatmap (
    TF_cor_sum[,TFrow_order],
 #   column_split = ifelse(TFrow_order_split =='1',2,1),
    column_names_rot =45, 
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 5),
    col =palette_enrichment, 
    cluster_rows=T,
    cluster_columns = F,
#     cell_fun = function(j, i, x, y, width, height, fill) {
#         grid.text(sprintf("%.0f", t(TF_cor_sum)[i, j]), x, y, gp = gpar(fontsize = 10, col='white'))
# },
    border=T))
dev.off()

pdf (file.path ('Plots','selected_TF_dev_corr_pathways_heatmaps.pdf'), width = 9.5,height=2.4)
hm
dev.off()

# Alternatively run DAM per cluster / sample ####
archp_NN$Clusters_sample = paste0(archp_NN$Clusters, '_',archp_NN$Sample2)
Clusters_sample = paste0(archp_NN$Clusters, '_',archp_NN$Sample2)
remove_low_clusters = !Clusters_sample %in% names(table (Clusters_sample)[table (Clusters_sample) < 10])
Clusters_sample = Clusters_sample[remove_low_clusters]
archp_NN = archp_NN[remove_low_clusters]
# archp2 = archp # store whole archr object in secondary object to run script only for the subset
# archp = archp_NN

srt_NN$Clusters_sample = srt_NN$sampleID

# Find DAM ####
metaGroupName = "Clusters_sample"
top_genes = Inf
force=TRUE

  top_genes = Inf
  DAM_df = DAM (
  ArchRProj = archp_NN,
  metaGroupName = metaGroupName,
  FDR_threshold = 1e-2,
  meandiff_threshold = 0,
  top_genes=top_genes,
  filter_by_scRNA=TRUE, # Make sure has same metaGroupName
  seurat_obj = srt_NN,
  min_exp=.1,
  force = force)

# Save table for supplementary information
write.csv (DAM_df, paste0('DAM_table_',metaGroupName, '.csv'))

top_genes = 3

DAM_df <- DAM_df %>%
group_by (comparison) %>%
slice_head(n = top_genes) %>%   # keep top 5 per celltype
ungroup() %>%
  mutate(
    # extract second term after "_" e.g. "P12"
    second_term = str_extract(comparison, "(?<=_)P\\d+"),
    # extract numeric part (12)
    second_num  = as.numeric(str_extract(second_term, "\\d+"))
  ) %>%
  arrange(second_num) %>%
  select(-second_term, -second_num)


if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames(mMat) = rowData(mSE)$name
mMat = mMat[unique (DAM_df$gene), rownames(archp_NN@cellColData)]
#mMat_mg = mMat[active_DAM, ]
mMat_mg = as.data.frame (t(mMat))
mMat_mg$metaGroup = as.character (archp_NN@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
library(stringr)
order_p <- as.numeric(str_extract(rownames(mMat_mg), "(?<=_P)\\d+"))
mMat_mg <- mMat_mg[order(order_p), ]
#mMat_mg = mMat_mg[celltype_order,]

#tf_labels = c(head(DAM_df$gene[DAM_df$comparison=='C1_P11'],top_genes))
ha = rowAnnotation (foo = anno_mark(at = match(head(DAM_df$gene[DAM_df$comparison=='C1_P11'],top_genes),colnames(mMat_mg)), 
    labels = head(DAM_df$gene[DAM_df$comparison=='C1_P11'],top_genes), labels_gp = gpar(fontsize = 8, fontface = 'italic')))
ha2 = HeatmapAnnotation (sample = sapply (rownames(mMat_mg), function(x) unlist(strsplit(x,'_'))[2]), col= list(sample = palette_sample))

library(stringr)

# extract second term: "P11", "P23", etc.
split_raw <- sapply(rownames(mMat_mg), function(x) unlist(strsplit(x,"_"))[2])

# extract numeric part
p_nums <- as.numeric(str_extract(split_raw, "\\d+"))

# get unique P groups in the correct numeric order
unique_levels <- unique(split_raw[order(p_nums)])

# make factor with ordered unique levels
split_factor <- factor(split_raw, levels = unique_levels)

DAM_hm = Heatmap (t(scale(mMat_mg)), 
          top_annotation = ha2,
          #column_split = sapply (rownames(mMat_mg), function(x) unlist(strsplit(x,'_'))[2]),
          column_split = split_factor,
          right_annotation = ha,
          row_labels = colnames (mMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 0, fontface = 'italic'),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          name = 'TF activity',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_deviation_fun(scale(mMat_mg))

          #right_annotation = motif_ha
          )

  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_heatmaps.pdf')), width = 6.5, height = 3.5)
print(DAM_hm)
dev.off()


### Discover epigenomic features correlated with sarcomatoid score ####
source (file.path('..','..','git_repo','tumor_analysis','epigenomic_features_correlated_to_scS_score.R'))
top_sarc_TF

# Compute co-occurrence of sarcomatoid TFs ####
motifMat = getPositions (archp)
matches = getMatches (archp)
matchesMat = assay (matches)
colnames (matchesMat) = gsub ('_.*','',colnames (matchesMat))
colnames (matchesMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames (matchesMat))
names (motifMat) = gsub ('_.*','',names (motifMat))
names (motifMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", names (motifMat))

matchesMat = matchesMat[,top_sarc_TF]
matchesMat = matchesMat[rowSums (matchesMat) > 0,]

cooc = matrix (ncol = ncol(matchesMat), nrow= ncol(matchesMat))

for (i in 1:ncol(matchesMat))
  {
  for (z in 1:ncol(matchesMat)) 
    {
    ov = sum (rowSums (matchesMat[,c(i,z)]) == 2) /  sum(colSums(matchesMat[,c(i,z)]))
    cooc[i,z] = ov
    }
  }

colnames (cooc) = top_sarc_TF
rownames (cooc) = top_sarc_TF
diag (cooc) = 0

cooc_hm = Heatmap (
    cooc,
    column_names_rot =45, 
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 6),
    col = rev(palette_cooccurrence), 
    cluster_rows=F,
    cluster_columns = F,
rect_gp = gpar (col = "white", lwd = 1)
)

pdf (file.path ('Plots','selected_TF_cooccurence_heatmaps.pdf'), width = 6,height=5)
cooc_hm 
dev.off()

# Distribution of SOX9 SOX6 motif distances in overlapping peaks ####
tf_pair = c('SOX9','SOX6')
top_sarc_TF = readRDS ('top_sarc_TF.rds')
motifMat = getPositions (archp)
matches = getMatches (archp)
matchesMat = assay (matches)
colnames (matchesMat) = gsub ('_.*','',colnames (matchesMat))
colnames (matchesMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames (matchesMat))
names (motifMat) = gsub ('_.*','',names (motifMat))
names (motifMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", names (motifMat))

ov_peaks_soxs = rowRanges (matches)[rowSums(matchesMat[,colnames (matchesMat) %in% tf_pair]) == 2]
#ranges_soxs =  [rowSums(matchesMat[, colnames (matchesMat) %in% tf_pair]) == 2]
#tmp = findOverlaps (motifMat[[tf_pair[2]]], ov_peaks_soxs, select='all')
#any (duplicated (queryHits (tmp)))

pos_in_peaks_tf1 =  motifMat[[tf_pair[1]]][subjectHits (findOverlaps (ov_peaks_soxs, motifMat[[tf_pair[1]]], select='all'))]
pos_in_peaks_tf1$peak = queryHits (findOverlaps (ov_peaks_soxs, motifMat[[tf_pair[1]]], select='all'))
pos_in_peaks_tf1 = pos_in_peaks_tf1[!duplicated (pos_in_peaks_tf1$peak)]
pos_in_peaks_tf2 =  motifMat[[tf_pair[2]]][subjectHits (findOverlaps (ov_peaks_soxs, motifMat[[tf_pair[2]]], select='all'))]
pos_in_peaks_tf2$peak = queryHits (findOverlaps (ov_peaks_soxs, motifMat[[tf_pair[2]]], select='all'))
pos_in_peaks_tf2 = pos_in_peaks_tf2[!duplicated (pos_in_peaks_tf2$peak)]
tf_peaks = unique (c(pos_in_peaks_tf1$peak,as.character(pos_in_peaks_tf2$peak)))
pos_in_peaks_tf1 = pos_in_peaks_tf1[match (tf_peaks, pos_in_peaks_tf1$peak)]
pos_in_peaks_tf2 = pos_in_peaks_tf2[match (tf_peaks, pos_in_peaks_tf2$peak)]
mid_tf1 = start(pos_in_peaks_tf1) + end(pos_in_peaks_tf1) / 2
mid_tf2 = start(pos_in_peaks_tf2) + end(pos_in_peaks_tf2) / 2
mi_diff = mid_tf1 - mid_tf2
#mi_diff[mi_diff > 100] = 100
#mi_diff[mi_diff < -100] = -100
mi_diff_df = as.data.frame (mi_diff)
hp = ggplot(mi_diff_df) +
  geom_histogram(
    aes(x = mi_diff),
    bins = 100
  ) +
  geom_vline(
    xintercept = 0,
    color = "red",
    linewidth = .1
  ) +
  scale_y_log10() +
  xlim(c(-100, 100)) +
  gtheme_no_rot

pdf (file.path ('Plots',paste0('distribution_distance_TF_pair_', paste0(tf_pair, collapse='_'),'hist.pdf')), height=2.5, width=3)
hp
dev.off()

### Blum et al barplot #### 
# Import Blum meta-analysis to compare with top TF correlated with scS-score
top_sarc_TF = readRDS ('top_sarc_TF.rds')
blum_df = read.csv (file.path('..','Blum_et_al_SE_score.csv')) 
blum_dfE = data.frame (gene = blum_df$Gene.Name, 
  score = blum_df$correlation.E.score, 
  SE_score = 'epithelioid',
  pval = blum_df$adjpvalue.E.score)
blum_dfS = data.frame (
  gene = blum_df$Gene.Name.1, 
  score = blum_df$correlation.S.score, 
  SE_score = 'sarcomatoid',
  pval = blum_df$adjpvalue.S.score)
blum_df = rbind (blum_dfE, blum_dfS)
blum_df = na.omit (blum_df)
blum_df$sig = blum_df$pval
blum_df$sig[blum_df$pval < 0.05 & blum_df$pval > 0.01] = '*'
blum_df$sig[blum_df$pval < 0.01 & blum_df$pval  > 0.001] = '**'
blum_df$sig[blum_df$pval < 0.001] = '***'
rownames (blum_df) = blum_df$gene
blum_df = blum_df[top_sarc_TF, ]
blum_df$gene = rownames (blum_df)
blum_df$score[is.na(blum_df$score)] = 0
blum_df$SE_score [is.na(blum_df$SE_score == 'epithelioid')] = 0
blum_df$score[blum_df$SE_score == 'epithelioid'] = -1 * blum_df$score[blum_df$SE_score == 'epithelioid']
blum_df$gene = factor (blum_df$gene, levels = rownames (blum_df)[order (-blum_df$score)])
blum_df = na.omit (blum_df)
# Create the dot plot
dp = ggplot(blum_df, aes(x = gene, y = score)) +
  geom_bar (position="stack", stat="identity",alpha=.7, aes(fill = SE_score)) + # Adds the points
  labs(title = "Correlation to Blum et al") + # Labels
  scale_fill_manual (values = c(epithelioid = 'blue',sarcomatoid = 'purple')) + 
  gtheme +
  geom_text(aes(label=sig), vjust=0.5, hjust = 0.5) #+
  #geom_text(aes(label=sig), vjust=1, hjust=0.5)

pdf (file.path ('Plots','Blum_top_sarc_TF.pdf'), height=3, width = 4)
dp
dev.off()


### SOX6/9 sc-Score_scatterplot ####
cnmf_spectra_unique = readRDS (file.path('..','..','git_repo','files','malignant_cnmf_genelist_25_nfeat_5000.rds'))
cnmf_spectra_unique = lapply (cnmf_spectra_unique, function(x) head(x, top_genes)[head(x, top_genes) %in% rownames (gsMat)])

if (!all (names(cnmf_spectra_unique) %in% colnames (srt@meta.data))) # Add module score of cNMF modules
  {
  srt = ModScoreCor (
          seurat_obj = srt, 
          geneset_list = cnmf_spectra_unique, 
          cor_threshold = NULL, 
          pos_threshold = NULL, # threshold for fetal_pval2
          listName = 'cNMF_', outdir = paste0(projdir,'Plots/'))

  }

avg_exp = log2(AverageExpression (srt, group.by = 'sampleID', feature = c('SOX6','SOX9'))[[1]]+1)
avg_exp = t(avg_exp)
ccomp_df = srt@meta.data
ccomp_df = aggregate (ccomp_df[,c('cNMF20')], by=as.list(srt@meta.data[,'sampleID',drop=F]), 'mean')
ccomp_df = cbind (ccomp_df, avg_exp[ccomp_df$sampleID,])
ccomp_df = ccomp_df[!rownames (ccomp_df) %in% c('HU37','HU62'),]
sp <- ggplot(ccomp_df, aes(x = x, y = SOX6)) + #, fill = sampleID, color = sampleID)) +
  geom_point(alpha = .8, shape = 21, stroke = 1, aes (fill = sampleID)) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson",
           label.x.npc = "left",  # place in left side of plot
           label.y.npc = "top") + # place near top
  scale_fill_manual(values = rev(palette_sample)) +
  scale_color_manual(values = rev(palette_sample)) +
  gtheme
sp2 <- ggplot(ccomp_df, aes(x = x, y = SOX9)) + #, fill = sampleID, color = sampleID)) +
  geom_point(alpha = .8, shape = 21, stroke = 1, aes (fill = sampleID)) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson",
           label.x.npc = "left",  # place in left side of plot
           label.y.npc = "top") + # place near top
  scale_fill_manual(values = rev(palette_sample)) +
  scale_color_manual(values = rev(palette_sample)) +
  gtheme  
pdf (file.path ('Plots','SOX6_9_scSscore_scatterplot.pdf'), height=3,width=3.5)
sp
sp2
dev.off()




#### Run inferCNV ####
library (infercnv)
library (gtools)

projdir_out = 'infercnv' # output folder 
dir.create (projdir_out, showWarnings = FALSE)

srt_cnv = srt[,!srt$sampleID %in% c('HU37','HU62')]
srt_cnv$cnv_type = srt_cnv$sampleID
ref = srt[,srt$sampleID %in% c('HU37','HU62')]
ref$cnv_type = 'reference'
srt_cnv = merge (srt_cnv, ref)

# load TxDbi object for the species and get genomic positions of genes
require (org.Hs.eg.db)
require (TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(txdb) = paste0('chr',1:22) # select chromosomes to include in the analysis
gene_regions = as.data.frame (genes (txdb))

# map entrez id to symbol ids
symbol = toTable (org.Hs.egSYMBOL)
gene_regions$symbol = symbol$symbol[match(gene_regions$gene_id, symbol$gene_id)]
gene_regions = gene_regions[gene_regions$symbol %in% rownames(srt_cnv),]

#Prepare files for inferCNV input
message ('Generate expression and annotation files')
#srt_cnv$orig.ident = srt_cnv@meta.data [,metaGroupName2]
malig_samples = unique(srt_cnv$cnv_type[srt_cnv$cnv_type != 'reference'])
        
exprMat = srt_cnv@assays$RNA@counts        
gene_regions2 = gene_regions[mixedorder(gene_regions$seqnames), ]
exprMat = exprMat[gene_regions2$symbol, ]

all (rownames(exprMat) == gene_regions2$symbol)
rownames (gene_regions2) = gene_regions2$symbol
colnames (gene_regions2) = NULL
message ('save gene_regions file')    
write.table (gene_regions2, file.path(projdir_out,'gene_regions.txt'), sep='\t')
#exprMat = exprMat[,rownames(annot_df)]
annot_df = srt_cnv@meta.data[,'cnv_type', drop=F]
annot_df[,1] = srt_cnv$cnv_type

colnames (annot_df) = NULL
message ('save annotation file')    
write.table (annot_df, file.path(projdir_out,'annot_df.txt'), sep='\t', col.names=NA)

chr_exclude = c('chrX','chry')
infercnv_obj = CreateInfercnvObject (raw_counts_matrix=exprMat,
                                    annotations_file= file.path(projdir_out,'annot_df.txt'),
                                    delim="\t",
                                    gene_order_file=file.path(projdir_out,'gene_regions.txt'),
                                    ref_group_names=c('reference'),
                                    chr_exclude= chr_exclude
                                    )

infercnv_result = infercnv::run(infercnv_obj,
                           cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                           out_dir=projdir_out,  # dir is auto-created for storing outputs
                           cluster_by_groups=T,  # cluster
                           denoise=T,
                           HMM=F,
                           save_rds = F,
                           no_prelim_plot = T,
                           no_plot = T,
                           plot_probabilities = FALSE
                           )

saveRDS (infercnv_result, file.path(projdir_out, 'infercnv.results.obj.Rds'))


#### Import inferCNV results to generate inferCNV plot but averaging by sample ####
library (biomaRt)
library (GenomicRanges)
library (BSgenome.Hsapiens.UCSC.hg38)
# cluster infercnv heatmap not by cluster
library (infercnv)
infercnv_result = readRDS (file.path(projdir_out, 'infercnv.results.obj.Rds'))
icnf_exp = infercnv_result@expr.data

makeWindows <- function(genome, blacklist, windowSize = 10e6, slidingSize = 2e6,
  gene_level_annotation = TxDb.Hsapiens.UCSC.hg38.knownGene){
  chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
  chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
  windows <- slidingWindows(x = chromSizes, width = windowSize, step = slidingSize) %>% unlist %>% .[which(width(.)==windowSize),]
  mcols(windows)$wSeq <- as.character(seqnames(windows))
    mcols(windows)$wStart <- start(windows)
    mcols(windows)$wEnd <- end(windows)
  message("Subtracting Blacklist...")
  windowsBL <- lapply(seq_along(windows), function(x){
      if(x %% 100 == 0){
        message(sprintf("%s of %s", x, length(windows)))
      }
      gr <- GenomicRanges::setdiff(windows[x,], blacklist)
      mcols(gr) <- mcols(windows[x,])
      return(gr)
    })
  names(windowsBL) <- paste0("w",seq_along(windowsBL))
  windowsBL <- unlist(GRangesList(windowsBL), use.names = TRUE)
  mcols(windowsBL)$name <- names(windowsBL)
  message("Adding Nucleotide Information...")
  windowSplit <- split (windowsBL, as.character(seqnames(windowsBL)))
  windowNuc <- lapply(seq_along(windowSplit), function(x){
    message(sprintf("%s of %s", x, length(windowSplit)))
      chrSeq <- Biostrings::getSeq(genome,chromSizes[which(as.character(seqnames(chromSizes))==names(windowSplit)[x])])
      grx <- windowSplit[[x]]
      aFreq <- alphabetFrequency(Biostrings::Views(chrSeq[[1]], ranges(grx)))
      mcols(grx)$GC <- rowSums(aFreq[, c("G","C")]) / rowSums(aFreq)
      mcols(grx)$AT <- rowSums(aFreq[, c("A","T")]) / rowSums(aFreq)
      return(grx)
    }) %>% GRangesList %>% unlist %>% sortSeqlevels %>% sort
  windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
  # get gene density
  gene_density = genes (gene_level_annotation)
  mcols(windowNuc)$gene_density = countOverlaps (windowNuc, gene_density)
  windowNuc
}

# Get GRanges of bins excluding black list regions
ws = 100000
ss = ws
windowSize = ws
slidingSize = ss
genome = BSgenome.Hsapiens.UCSC.hg38
chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
windows <- slidingWindows(x = chromSizes, width = windowSize, step = slidingSize) %>% unlist %>% .[which(width(.)==windowSize),]
  
#region_df = do.call (rbind, region)
# specify the database
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# loop through rows, get genes, then paste with collapse,
# and finally bind back with data d.
gene_info <- getBM (attributes = c("chromosome_name", "start_position", "end_position","external_gene_name"),
                   filters = "external_gene_name",
                   values = rownames(icnf_exp),
                   mart = ensembl)

gene_info_filtered = gene_info[gene_info$chromosome_name %in% c(1:22),]
colnames (gene_info_filtered) = c('chr','start','end')
gene_info_filtered$chr = paste0('chr',gene_info_filtered$chr)
gene_info_gr = makeGRangesFromDataFrame (gene_info_filtered, keep.extra.columns=T)

ov = findOverlaps (windows, gene_info_gr)
qhits = queryHits (ov)
shits = subjectHits (ov)
icnv_regions = lapply (unique(qhits), function(x) 
tmp = colMeans(icnf_exp[gene_info_gr@elementMetadata[,1][shits[which (qhits == x)]],,drop=F]))

region_names = as.character(windows)[unique(qhits)]
chr_names = sapply (region_names, function(x) unlist(strsplit (x,'\\:'))[1])
icnv_regions_df = do.call (cbind, icnv_regions)
sample_row = srt_cnv$cnv_type
all (names (sample_row) == colnames(icnf_exp))
icnv_regions_df_sample = lapply (unique(sample_row), function(x) colMeans (icnv_regions_df[sample_row == x,, drop=F]))
icnv_regions_df_sample = do.call (cbind, icnv_regions_df_sample)
colnames (icnv_regions_df_sample) = unique(sample_row)
sample_to_keep = unique (sample_row)[unique (sample_row) != 'reference']
icnv_regions_df_sample = icnv_regions_df_sample[,sample_to_keep]
rownames (icnv_regions_df_sample) = region_names
icnv_regions_df_sample = t(icnv_regions_df_sample)
colnames (icnv_regions_df_sample) = chr_names
chr_names= factor (chr_names, levels = unique(chr_names))
#rownames(icnv_regions_df_sample) = c('p786','p811','p826','p846','p848','p4','p8','p7','p9','p11','p12','p13')

icnv_regions_df_sample_c = icnv_regions_df_sample
icnv_regions_df_sample_c[icnv_regions_df_sample_c > 2] = 2
icnv_regions_df_sample_c[icnv_regions_df_sample_c < 0] = 0
palette_cnv = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging",3))
icnv_regions_df_sample_c = icnv_regions_df_sample_c[c('P1','P3','P4','P5','P8','P11','P12','P13','P14'),]
ht = Heatmap (
  icnv_regions_df_sample_c,
  #row_order = ,
  column_gap = unit(1,'mm'),
  column_names_gp = gpar(fontsize = 0),
  cluster_column_slices = FALSE,
  column_split = chr_names, 
  border=T,
  cluster_columns=F,
  cluster_rows=F)
  #colorRamp2(c(0.8, 1, 1.2), c(palette_cnv[1], palette_cnv[2], palette_cnv[3])))
pdf (file.path ('Plots','cnv_sample_mean_heatmap.pdf'),height=5, width=9)
ht
dev.off()

