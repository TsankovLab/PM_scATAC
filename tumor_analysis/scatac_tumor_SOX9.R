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
addArchRThreads (threads = 1)
addArchRGenome ("hg38")

if (!file.exists ('Save-ArchR-Project.rds')) 
  { source (file.path('..','..','PM_scATAC','scatac_tumor_create_ArchRobj.R'))
  } else {
 archp = loadArchRProject (projdir)   
  }

#archp$Sample3 = archp$Sample2
#archp$Sample3[archp$Clusters == 'C12'] = 'P11_HOX'
#archp$Sample3[grep ('normal',archp$Sample3)] = 'normal'

# Load RNA ####
srt = readRDS (file.path('..','scrna','srt.rds'))
srt$sampleID3[srt$sampleID3 %in% c('HU37','HU62')] = 'normal_pleura'
sarc_order = read.csv (file.path('..','scrna','cnmf20_sarcomatoid_sample_order.csv'), row.names=1)
sarc_order = sarc_order[! sarc_order$sampleID %in% c('HU37','HU62'),]
sarc_order = rbind (data.frame (sampleID = 'normal_pleura', x = -1),sarc_order) 
#archp$Sample2 = factor (archp$Sample2, levels = sarc_order$sampleID)

# Check SOX15 and SOX9 expression
pdf (file.path('Plots','SOX15_SOX9_expression.pdf'),width=5)
DotPlot (srt, features =c('SOX9','SOX15','SOX6','RUNX2','RUNX1','RUNX3','SNAI2','TWIST1','TEAD1','ZN784'), group.by = 'sampleID') + gtheme
dev.off()

markers = c('RUNX1','RUNX2','THAP1','MZF1', 'GLIS3', 'PLAGL1', 'SOX9','VIM','SERPINE2','SNAI2')
metaGroupName='Sample'
#archp = addImputeWeights (archp)
pdf()
p2_l <- plotGroups(
  ArchRProj = archp, 
  groupBy = metaGroupName, 
  colorBy = "GeneScoreMatrix", 
  name = markers,
  plotAs = "violin",
  alpha = 0.4,
  getImputeWeights=NULL,
  addBoxPlot = TRUE#,
  #pal = palette_tnk_cells
 )
dev.off()
pdf (file.path ('Plots','_RUNXs_genescores.pdf'),14,14)
print (wrap_plots (p2_l))
dev.off()


archp = addClusters (input = archp, resolution = 1.5,
  reducedDims = "IterativeLSI", name = 'Clusters2',
  maxClusters = 100,
  force = TRUE)
archp_NN = archp[!archp$Sample3 %in% c('normal1','normal2','normal3')]

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

# Run genescore DAG ####
archp$status = ifelse (archp$Sample3 == 'normal_pleura','normal','tumor')
metaGroupName = "status"
force = FALSE
if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source (file.path('..','..','git_repo','utils','DAG.R'))

archp = addImputeWeights (archp)
celltype_markers = c('WT1','CALB2','RUNX2','TCF3','SOX9','VIM','AXL','SOX6','MESP1','HMGA1','TWIST1','SNAI2')
uncommitted_markers = c('AR','CDCA7','DNMT3A','HHIP','HMGA2','KCNQ1OT1','MEST','NOS1','PAPLN','PRDM6')
pdf()
p2 <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "GeneScoreMatrix", 
    name = celltype_markers, 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()

pdf()
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "GeneScoreMatrix", 
    name = uncommitted_markers, 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path('Plots','marker_genes_feature_plots.pdf'), width = 20, height = 20)
print (wrap_plots (p2, ncol = 4))
print (wrap_plots (p, ncol= 4))
dev.off()



# Rank tumor samples by sarcomatoid score ####
# Add cNMF identified in scRNA to archr object ####
if (!exists('gsSE')) gsSE = fetch_mat (archp, 'GeneScore')
gsMat = assays (gsSE)[[1]]
rownames (gsMat) = rowData (gsSE)$name

nfeat=5000
k=25
top_genes=50
cnmf_spectra_unique = readRDS (paste0('../scrna/cnmf_genelist_',k,'_nfeat_',nfeat,'.rds'))
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
sams = as.character(unique(archp_meta$Sample3))
tumor_sams = sams[!sams %in% c('normal1','P11_HOX')] # remove normal,low cell numbers and outlier samples

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
peak_reproducibility='1' # Set to 1 to better identify tumor heterogeneity
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
#any(lapply(corTF_array, function(x) any(is.na(x))))
# Take element-wise median
median_matrix <- apply (corTF_array, c(1, 2), median)

# set.seed(123)
# centers=3
# km = kmeans (median_matrix, centers=centers)
# km_df = as.data.frame (km$cluster)
# km_df = km_df[order (km_df[,1]),,drop=F]
set.seed (123)
# km = kmeans (median_matrix, centers=2)
# km$cluster[km$cluster == 1]
hr = hclust(as.dist(1-median_matrix))#, method = "average")
clusters = dendextend::cutree(hr, k = 2)

pdf()
set.seed (1234)
cor_TF_df = draw (Heatmap (median_matrix,
  #left= ha1,
#  row_split = clusters,
#  column_split = clusters,
  row_names_gp = gpar(fontsize = 6),
  clustering_distance_rows='pearson',
  clustering_distance_columns='pearson',
  column_names_gp = gpar(fontsize = 6),
  col = palette_deviation_cor_fun, border=T))

  # rect_gp = gpar(type = "none"),
  # cell_fun = function(j, i, x, y, w, h, fill) {
  #       if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
  #           grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
  #       }}))
dev.off()

pdf (file.path ('Plots','selected_TF_dev_corr_heatmaps.pdf'), width = 8,height=7)
cor_TF_df 
dev.off()

# Add pathways of TF correlated genes from scRNA ####
if (!file.exists('enrichment_pathways_TFs.rds')) source (file.path('..','..','git_repo','tumor_analysis','enrichment_cnmfs.R'))
TF_cor_sum = readRDS (file.path('enrichment_pathways_TFs.rds'))
TFrow_order = unname(unlist(row_order (cor_TF_df)))
#TFrow_order_split = rep (names(row_order (cor_TF_df)), lapply (row_order (cor_TF_df),length))
#TFrow_order_split = TFrow_order_split[TFrow_order]
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


# #### Generate heatmap of significant Cox proportion hazards from bulk-RNA ####
# coxht = readRDS (file.path('..','..','bulkRNA_meso','oncoTF_sig_cox_bulkRNA.rds'))



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

# Check distribution of SOX9 SOX6 motif distances in overlapping peaks ####
tf_pair = c('ZEB1','MESP1')
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

pdf (file.path ('Plots',paste0('distribution_distance_TF_pair_', paste0(tf_pair, collapse='_'),'hist.pdf')), height=3, width=3)
x = hist (mi_diff, breaks=200)
dev.off()

x$mid[which.max (x$counts)]


### Import Blum meta-analysis to compare with top TF correlated with scS-score ####
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



# # Check overall genescore and RNA expression of TFs correlated with scS-score ####
# metaGroupName = 'Sample'
# gsgSE = getGroupSE(
#   ArchRProj = archp,
#   useMatrix = 'GeneScoreMatrix',
#   groupBy = metaGroupName,
#   divideN = TRUE,
#   scaleTo = NULL,
#   threads = getArchRThreads(),
#   verbose = TRUE,
#   logFile = createLogFile("getGroupSE")
# )
# gsmat = assays (gsgSE)[[1]]
# rownames (gsmat) = rowData (gsgSE)$name
# gsmat = gsmat[top_sarc_TF,sams]
# gsmat = as.data.frame (gsmat)
# gsmat$TF = rownames (gsmat)
# gs_sarc_TF_df = gather (gsmat, sample, score, 1:(ncol(gsmat)-1))
# gs_sarc_TF_df$type = 'genescore'

# rna_sarc_tf = as.data.frame (log2(AverageExpression (srt, features = top_sarc_TF, group.by = 'sampleID')[[1]]+1))
# rna_sarc_tf = rna_sarc_tf[,colnames (rna_sarc_tf) %in% sams]
# rna_sarc_tf$TF = rownames (rna_sarc_tf)
# rna_sarc_TF_df = gather (rna_sarc_tf, sample, score, 1:(ncol(rna_sarc_tf)-1))
# rna_sarc_TF_df$type = 'expression'

# sarc_TF_df = rbind (gs_sarc_TF_df, rna_sarc_TF_df)
# sarc_TF_df$TF = factor (sarc_TF_df$TF, levels = top_sarc_TF)

# bp1 = ggplot (sarc_TF_df[sarc_TF_df$type == 'genescore',], aes(x = TF, y = score, fill = TF, color=type)) + 
#   geom_boxplot(
#     #aes(group = interaction(TF, type)),  # Split by both TF and type
#     position = position_dodge(0.8),
#     #color = 'grey20',
#     linewidth = 0.4,
#     width = 0.7,
#     outlier.alpha = 0.2,
#     outlier.shape = NA,
#     size = 0.5,
#     alpha = 0.6
#   ) +
#   geom_point(
#     #aes(group = interaction(TF, type), color = type),  # Align points with boxplots
#     position = position_dodge(0.8),  # Ensure same dodge width
#     alpha = 0.5,
#     size = 1
#   ) +
#   gtheme +
#   scale_fill_manual(values = palette_expression_disc) +
#   scale_color_manual (values = c(genescore = 'blue', expression = 'darkgreen')) +
#   ylim (c(0,2.5))
#   #geom_hline(yintercept = 0, color = 'blue', linetype = 'dashed')
# bp2 = ggplot (sarc_TF_df[sarc_TF_df$type == 'expression',], aes(x = TF, y = score, fill = TF, color=type)) + 
#   geom_boxplot(
#     #aes(group = interaction(TF, type)),  # Split by both TF and type
#     position = position_dodge(0.8),
#     #color = 'grey20',
#     linewidth = 0.4,
#     width = 0.7,
#     outlier.alpha = 0.2,
#     outlier.shape = NA,
#     size = 0.5,
#     alpha = 0.6
#   ) +
#   geom_point(
#     #aes(group = interaction(TF, type), color = type),  # Align points with boxplots
#     position = position_dodge(0.8),  # Ensure same dodge width
#     alpha = 0.5,
#     size = 1
#   ) +
#   gtheme +
#   scale_fill_manual(values = palette_expression_disc) +
#   scale_color_manual (values = c(genescore = 'blue', expression = 'darkgreen'))# +
#   #geom_hline(yintercept = 0, color = 'blue', linetype = 'dashed')

# pdf (file.path ('Plots','sarcomatoid_score_TF_genescore_and_expression_boxplots.pdf'), width = 5,height=3)
# bp1
# bp2
# dev.off()




### Identify hubs / peaks correlated with sarcomatoid score ####

### Hubs analysis #####
metaGroupName = "Sample3"
cor_cutoff = 0.3
#max_dist = 12500
max_dist = 12500
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
dir.create(file.path (hubs_dir, 'Plots'), recursive=T)
hubs_obj = readRDS (file.path(hubs_dir,'global_hubs_obj.rds'))  
hubsCell_mat = readRDS (file.path (hubs_dir,paste0('hubs_cells_mat.rds')))  

# Select samples ####
sams = as.character(unique(archp$Sample3))
sams = sams[!sams %in% c('normal2','normal3','normal4','P3','P13','P11_HOX')] # remove normal,low cell numbers and outlier samples

# Identify hubs higher in malignant ####
library (presto)
pval_threshold = 0.01
occurrence_threshold = 4

test_comparisons = as.list(
  paste0(sams, ' normal1'))
names(test_comparisons) = sams
test_comparisons = test_comparisons[-10]
test_comparisons = lapply (test_comparisons, function(x) unlist(strsplit(x, ' ')))
#archp = archp[rownames(archp@cellColData) %in% colnames(hubsCell_mat)]
all (colnames(hubsCell_mat) == rownames(archp@cellColData))
test_res = lapply (test_comparisons, function(sam_pair) 
  wilcoxauc (log2(hubsCell_mat+1), y= archp$Sample3, groups_use = sam_pair))
test_res = lapply (test_res, function(x) x[x$logFC > 0,])
test_res = do.call (rbind, test_res)
test_res$gene = hubs_obj$hubsCollapsed$gene[match(test_res$feature, hubs_obj$hubs_id)]
test_res = test_res[test_res$padj < pval_threshold,]
test_res = test_res[test_res$group != 'normal1',]
test_res = test_res[test_res$avgExpr > .5,]
test_res = test_res[test_res$logFC > .5,]
malig_hubs = names(table (test_res$feature)[table (test_res$feature) > occurrence_threshold])

hubsCell_mat_sub = scale (log2(hubsCell_mat[malig_hubs,]+1))

bin_width <- 50   # Number of observations per bin
overlap <- 50    
sarc_module = 'cNMF20'


# Get cnmf modules ####
archp_meta = as.data.frame (archp@cellColData)
cnmf_mat = archp@cellColData[,grep ('cNMF',colnames(archp@cellColData))]
cnmf_mat = as.data.frame (t(scale (t(cnmf_mat))))
cnmf_mat = lapply (sams, function(x) cnmf_mat[archp_meta$Sample3 == x, ])
names (cnmf_mat) = sams

cnmfMat_ordered = lapply (sams, function(sam) cnmf_mat[[sam]][order(cnmf_mat[[sam]][,sarc_module]),])
names(cnmfMat_ordered) = sams
cnmfMat_ordered_avg = lapply (sams, function (sam) as.data.frame (lapply (as.data.frame (cnmfMat_ordered[[sam]]), function(x) {
  zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
})))
names (cnmfMat_ordered_avg) = sams

all (colnames(hubsCell_mat_sub) == rownames(archp@cellColData))
hubsCell_mat_sub = lapply (sams, function(x) hubsCell_mat_sub[,archp_meta$Sample3 == x])
names (hubsCell_mat_sub) = sams
hubsMat_ordered = lapply (sams, function(sam) hubsCell_mat_sub[[sam]][,order(cnmf_mat[[sam]][,sarc_module])])
names(hubsMat_ordered) = sams
hubsMat_ordered_avg = lapply (sams, function (sam) as.data.frame (lapply (as.data.frame (t(hubsMat_ordered[[sam]])), function(x) {
  zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
})))
names (hubsMat_ordered_avg) = sams
  
# Correlate hubs with sarcomatoid module per sample ####
hub_cor_l = list()
for (sam in sams)
  {
  tmp = cor (hubsMat_ordered_avg[[sam]], cnmfMat_ordered_avg[[sam]][,sarc_module], method='spearman')
  tmp = as.data.frame (tmp)
  tmp$sample = sam
  tmp$hub_id = rownames(tmp)
  colnames(tmp) = c('score','sample','hub_id')
  hub_cor_l[[sam]] = tmp
  }

# load cnmf modules to exclude labelling any genes in the sarcomatoid module ####
nfeat=5000
k=25
top_genes=50
sarc_module = 'cNMF20'
cnmf_spectra_unique = readRDS (paste0('../scrna/cnmf_genelist_',k,'_nfeat_',nfeat,'.rds'))
cnmf_spectra_unique = lapply (cnmf_spectra_unique, function(x) head(x,top_genes))
sarc_module_genes = cnmf_spectra_unique[[sarc_module]]
#cnmf_spectra_unique = lapply (cnmf_spectra_unique, function(x) head(x, 50)[head(x, 50) %in% rownames (gsMat)])

hubs_cor_df = do.call (rbind, hub_cor_l)
hubs_cor_df$score[is.na(hubs_cor_df$score)] = 0
hub_cor_levels = hubs_cor_df %>% group_by(hub_id) %>% summarise(median_value = median(score)) %>% arrange (-median_value)
hub_cor_levels$hub_id = factor (hub_cor_levels$hub_id, levels = hub_cor_levels$hub_id)
hub_cor_levels$gene = paste0(hubs_obj$hubs_id,'-',hubs_obj$hubsCollapsed$gene)[match(hub_cor_levels$hub_id, hubs_obj$hubs_id)]

selected_TF = readRDS ('selected_TF.rds')
hub_cor_levels$in_module = as.character(rowSums (sapply (sarc_module_genes, function(x) grepl (x, hub_cor_levels$gene))))
hub_cor_levels$sarc_TF = ''
hub_cor_levels$sarc_TF[unlist(sapply (selected_TF, function(x) grep(x, hub_cor_levels$gene)))] = hub_cor_levels$gene[unlist(sapply (selected_TF, function(x) grep(x, hub_cor_levels$gene)))]
hub_cor_levels$hit = ifelse (hub_cor_levels$sarc_TF != '','hit','nohit')

bp = ggplot (hub_cor_levels[hub_cor_levels$median_value > 0,], 
  aes (x= hub_id, y = median_value, label = sarc_TF, fill = hit), color='white') + 
geom_bar (stat = 'identity') +
scale_fill_manual (values = c(`hit` = 'darkred',`nohit` = 'grey66')) + 
#scale_color_manual (values = c(`TRUE` = 'darkred',`FALSE` = 'grey66')) + 
geom_text_repel (size=1.8, 
  max.overlaps=10000, 
  nudge_y = .1, 
  nudge_x = .2, 
  segment.curvature = 0.1,
  segment.size=.2) +
#geom_text(aes(label = labels, y = median_value + .001), hjust = -0.05,na.rm = TRUE, angle=90,size=1) +
gtheme_no_rot +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) 


pdf (file.path ('Plots','hubs_cor_ranked.pdf'), width=4, height=3)
bp
dev.off()





# Check top correlated hubs on browser track ####
top_hubs = as.character(head(hub_cor_levels$hub_id, 50))
metaGroupName='Sample3'
matching_samples=c('normal_pleura','P1','P3','P4','P5','P8','P11','P11_HOX','P12','P14','P23')
TF = 'RUNX2'
TF = 'SOX6'
TF = 'SOX9'
TF = 'SNAI2'
sample_sarc_order_levels = levels(sample_sarc_order)
#sample_sarc_order_levels = sample_sarc_order_levels[c(7,1:6,8:13)]
sample_sarc_order_levels = sample_sarc_order_levels[! sample_sarc_order_levels %in% 'normal1']
palette_sample_scs = as.character(paletteer::paletteer_c("grDevices::Purple-Blue",length (unique(archp_NN@cellColData[,metaGroupName]))))
names (palette_sample_scs) = sample_sarc_order_levels
pdf()
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp_NN, 
    sizes = c(6, 1, 1, 1,1,1),
    groupBy = metaGroupName, 
    #region = ext_range(hubs_obj$hubsCollapsed[match(top_hubs, hubs_obj$hubs_id)],50000,50000),
    region = ext_range (hubs_obj$hubsCollapsed[grep('SOX9', hubs_obj$hubsCollapsed$gene)],10000,50000),
    sample_levels = sample_sarc_order_levels,
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
    pal = palette_sample_scs,
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
  name =paste0(TF,'_sarc_cor_hubs_sarc_TFs_coveragePlots.pdf'),
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




# TF motifs in hubs correlating with sarcomatoid score ####
TFmatch = getMatches (archp)
TFpositions = getPositions (archp)
names(TFpositions) = gsub ('_.*','',names(TFpositions))
names(TFpositions) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", names(TFpositions))
TFs = c('SOX9','RUNX2','HIC2','SOX5','SOX6','ZEB1','TWIST1')
TFpositions = TFpositions[TFs]
top_hubs_gr = lapply (top_hubs, function(x) unlist(lapply (TFpositions, function(y)
  sum(countOverlaps(y, hubs_obj$hubsCollapsed[match(x,hubs_obj$hubs_id)])))))
top_hubs_mat = do.call (cbind, top_hubs_gr)
top_hubs_df = as.data.frame (top_hubs_mat)
colnames(top_hubs_df) = top_hubs
top_hubs_df$TF = rownames(top_hubs_df)
top_hubs_df = gather (top_hubs_df, hub, hit, 1:(ncol(top_hubs_df) - 1))
top_hubs_df$hub = factor (top_hubs_df$hub, levels = top_hubs_df$hub)
#colnames(top_hubs_df)[1] = 'hit'
top_hubs_df$hit[top_hubs_df$hit == 0] = NA
dp = ggplot(top_hubs_df, aes(x = hub, y = TF)) +
  geom_point(aes (size = hit), color = 'darkred') + # Adds the points
  labs(title = "motif hit") +#+ # Labels
  gtheme
  #scale_color_manual (values = c(epithelioid = 'darkgreen',sarcomatoid = 'firebrick2')) + gtheme

pdf (file.path ('Plots','Motif_hit_sarc_hubs.pdf'), height=3)
dp
dev.off()


# plot top correlated hubs ####
hubs_cor_df$hub_id = factor (hubs_cor_df$hub_id, levels = hub_cor_levels$hub_id)
#lapply (tc_cor, function(x) {x = x['cNMF19',]; head(x[order(-x)],10)})

top_hubs = head(hub_cor_levels$hub_id, 20)
hubs_cor_df_top = hubs_cor_df[hubs_cor_df$hub_id %in% top_hubs, ]
bp = ggplot (hubs_cor_df_top, aes (x = hub_id, y = score), alpha=.5) + 
geom_boxplot (color = 'grey66',
    fill = 'brown',
    linewidth = .1,
    width=.5,
    outlier.alpha = 0.2,
    outlier.shape = NA,
     size=0.5, alpha=0.7
     ) + 
geom_point (position = 'jitter', alpha= 0.2, color = 'grey33', size=1) +
gtheme +
#scale_fill_manual (values = c(activity = '#B2183BFF', expression = '#7663A3FF')) + 
geom_hline (yintercept = 0, color='red',  linetype='dashed')

pdf (paste0 ('Plots/sarcomatoid_score_hubs_cor_boxplots.pdf'), width = 7,height=3)
bp
dev.off()








# Show sarcomatoid module ####
#metaGroupName = "Clusters"
#force=FALSE
#if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source ('../../PM_scATAC/DAG.R')

sarc_module = 'cNMF20'
tf_markers = c('SOX9','SNAI2','TWIST1','SOX6','HIC2','RUNX2')

archp = addImputeWeights (archp)
pdf()
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = sarc_module, 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
markerMotifs = getFeatures (archp, select = paste(tf_markers, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs = grep ("z:", markerMotifs, value = TRUE)
TF_p = plotEmbedding(
    ArchRProj = archp, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    pal = rev(palette_deviation),
    imputeWeights = getImputeWeights(archp)
)

TF_p2 = plotEmbedding(
    ArchRProj = archp, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    pal = rev(palette_deviation),
    imputeWeights = NULL
)

dev.off()

pdf (file.path('Plots','sarcomatoid_score_feature_plots2.pdf'), width = 20, height = 20)
print (wrap_plots (c(list(p), TF_p), ncol = 3))
print (wrap_plots (TF_p2), ncol = 3)
dev.off()

# Plot the same but scaled ####
scaled_cnmfs = archp@cellColData[grep ('cNMF',colnames(archp@cellColData))]
scaled_cnmfs = as.data.frame (t(scale(t(scaled_cnmfs))))
umap_df = data.frame (archp@embeddings$UMAP[[1]], scaled_cnmfs)
umap_p1 = ggplot(data = umap_df) + 
geom_point (aes_string (x = 'IterativeLSI.UMAP_Dimension_1', y= 'IterativeLSI.UMAP_Dimension_2', color = cnmf_module), size = .1) + 
scale_colour_gradientn (colours = rev(brewer.pal (n = 11, name = "RdBu")),limits=c(-max (abs(umap_df[,cnmf_module])), max (abs(umap_df[,cnmf_module])))) +
ggtitle ('sarcomatoid_score') + 
#facet_wrap (as.formula(paste("~", metaGroupNames[3]))) + 
theme_classic() +
theme_void()

pdf (file.path('Plots','sarcomatoid_score_scaled_feature_plots.pdf'), width = 6, height = 6)
umap_p1
dev.off()
  
# # Order cells per samples along SOX9 deviation and correlated TFs ####
# selected_TF = readRDS ('selected_TF.rds')
# library (scales)
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name

# # Filter by RNA expression ####
# metaGroupName = 'sampleID'
# active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)
# mMat = mMat[active_TFs, ]

# archp_meta = as.data.frame (archp@cellColData)
# archp_meta$Sample3 = archp_meta$Sample2
# #archp_meta$Sample3[archp_meta$Clusters == 'C15'] = 'P11_HOX'

# all (colnames(mMat) == rownames(archp_meta))

# TF_driver = 'SOX9'
# top_TFs = 50
# traj_sample = list()
# for (sam in unique(archp_meta$Sample3))
#     {
#     library(zoo)
#     bin_width <- 100   # Number of observations per bin
#     overlap <- 1    
#     mMat_ordered_sample = as.data.frame(scale(mMat[,archp_meta$Sample3 == sam]))
#     mMat_ordered_sample = mMat_ordered_sample[, order(unlist(mMat_ordered_sample[TF_driver,]))]
#     cor_to_tf = order(-as.vector(cor (t(mMat_ordered_sample)[,TF_driver],t(mMat_ordered_sample), method='pearson')))
#     mMat_ordered_sample = mMat_ordered_sample[c(head(cor_to_tf,top_TFs),tail(cor_to_tf,top_TFs))  ,]
#     mMat_ordered_sample = as.data.frame(lapply(as.data.frame(t(mMat_ordered_sample)), function(x) {
#       zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
#     }))
#     traj_sample[[sam]] = Heatmap (
#       t(as.data.frame(lapply(mMat_ordered_sample, rescale, to = c(-10,10)))), 
#       col = rev(palette_deviation), 
#       cluster_columns=F,
#       name = sam,
#       cluster_rows=F,
#       column_names_gp = gpar(fontsize = 0),
#       row_names_gp = gpar(fontsize = 7, fontface='italic'))
#     }

# pdf (file.path('Plots',paste0('sarc_trajectory_',TF_driver,'_sample.pdf')), height=12, width=3)
# traj_sample
# dev.off()



# # Plot UMAP per sample showing sarcomatoid score and correlated TFs e.g. SOX9 HIC2.. ####
# sams = unique(archp$Sample2)
# sams = sams[!sams %in% c('normal_pleura','P3','P13')]
# cnmf_mat = as.matrix(archp@cellColData[,grep ('sarcomatoid', colnames(archp@cellColData))])
# cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp$Sample2 == x, ])))
# names(cnmf_mat) = sams
# #cnmf_mat = do.call (rbind, lapply (cnmf_mat, function(x) t(x)))
# #umap_df = data.frame(archp@embeddings$UMAP[[1]], archp$Sample2)

# varfeat = 1000
# LSI_method = 2
# pdf()
# sample_LSI = lapply (sams, function(x) {
#   tmp = addIterativeLSI (ArchRProj = archp[archp$Sample2 == x],
#     useMatrix = "MotifMatrix", name = "IterativeLSI",
#     force = TRUE, LSIMethod = LSI_method,
#     varFeatures = varfeat)
#     addUMAP (ArchRProj = tmp, 
#     reducedDims = "IterativeLSI",
#     force = TRUE)@embeddings$UMAP[[1]]
#   })
# dev.off()


# names (sample_LSI) = sams
# umap_df = lapply (sams, function(x) data.frame (sample_LSI[[x]][rownames(t(cnmf_mat[[x]])),], t(cnmf_mat[[x]])))
# umap_df = do.call (rbind, umap_df)
# umap_df$Sample2 = sapply (rownames(umap_df), function(x) unlist(strsplit(x, '\\#'))[1])

# umap_p1 = ggplot (data = umap_df) + 
# geom_point (aes (x = IterativeLSI.UMAP_Dimension_1, y= IterativeLSI.UMAP_Dimension_2, color = sarcomatoid.cNMF20), size = .01) + 
# scale_colour_gradientn (colours = rev(palette_deviation)) +
# ggtitle ('sarcomatoid_score') +
# scale_size (range = c(0.1, .4)) + 
# facet_wrap (~Sample2) + 
# theme_void ()

# pdf (file.path('Plots','sarcomatoid_score_scaled_per_sample_feature_plots.pdf'), width = 10, height = 5)
# umap_p1
# dev.off()

  












# Check correlation of SOX9 expression and deviation with sarcomatoid module to pick sample for chrombpnet analysis ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat[selected_TF,])

archp_meta = as.data.frame (archp@cellColData)
sams = unique(archp_meta$Sample3)
sams = sams[!sams %in% c('normal_pleura','P3','P13','P11_HOX')]

# Compute correlation of sarcomatoid cNMF with TFs ####
cnmf_mat = as.matrix(archp@cellColData[,grep ('sarcomatoid', colnames(archp@cellColData))])
cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp_meta$Sample3 == x, ])))
names (cnmf_mat) = sams
tf_mat = lapply (sams, function(x) scale(mMat[,archp_meta$Sample3 == x]))
names (tf_mat) = sams

tf = 'SOX9'
mod = 'cNMF20'
sarc_tf = lapply (sams, function(sam)  data.frame (barcode = colnames(tf_mat[[sam]]), TF = as.vector(tf_mat[[sam]][tf,]), module = as.vector(cnmf_mat[[sam]][mod,,drop=F]), sample=sam))
sarc_tf_df = do.call (rbind,sarc_tf)

# Compute TF correlation to sarcomatoid in scRNA ####
metacells = readRDS (file.path('..','scrna','metacells.rds'))
metacells$sampleID = metacells$sampleID3
nfeat=5000
k=25
top_genes=50
cnmf_spectra_unique = readRDS (paste0('../scrna/cnmf_genelist_',k,'_nfeat_',nfeat,'.rds'))


if (!all (names(cnmf_spectra_unique) %in% colnames(metacells@meta.data)))
{
metacells = ModScoreCor (
        seurat_obj = metacells, 
        geneset_list = cnmf_spectra_unique, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'cNMF_', outdir = paste0(projdir,'Plots/'))  
}
metacells_assay = metacells@assays$RNA@layers$data
rownames (metacells_assay) = rownames (srt)
metacells_assay = metacells_assay[tf,,drop=F]

sams = c('P1','P4','P5','P8','P11','P12')
tc_cor = lapply (sams, function(sam)
  {
  tmp = data.frame(TF = metacells_assay[tf,metacells$sampleID == sam], module = metacells$cNMF20[metacells$sampleID == sam])
  tmp$sample = sam
  tmp
  })
scrna_tf_cor_df = do.call (rbind, tc_cor)

# plot correlation TF vs module score in atac and RNA space
sp = ggplot (sarc_tf_df, aes (x= TF, y = module, color = sample)) + 
geom_point(alpha = .4) + facet_wrap (~sample, ncol = length(unique(sarc_tf_df$sample)), scales = 'free') +
scale_color_manual (values = palette_sample) + gtheme_no_rot +                                      
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") +
  stat_cor(method = "spearman", color ='black')
sp2 = ggplot (scrna_tf_cor_df, aes (x= TF, y = module, color = sample)) + geom_point(alpha = .4) + 
facet_wrap (~sample, ncol = length(unique(sarc_tf_df$sample)), scales = 'free') + 
scale_color_manual (values = palette_sample) + gtheme_no_rot +                                      
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") +
  stat_cor(method = "spearman", color ='black')
pdf (file.path ('Plots',paste0('correlation',tf,'_module',mod,'.pdf')), width=20, height = 3)
sp
sp2
dev.off()

# use P23 as having most cells and define sarc score high and low cells ####
cnmf_module = 'cNMF20'
cnmfs = archp@cellColData[,grep ('cNMF',colnames(archp@cellColData))]
cnmfs = as.data.frame(cnmfs)
cnmfs_scaled = as.data.frame (t(scale (t(cnmfs))))
top_cells = 2000
cnmf_labelled = lapply (split(cnmfs_scaled, archp$Sample3), function(x)
  {
  df1 = data.frame (barcode = head (rownames(x)[order(x[,cnmf_module])], top_cells),
    subtype = 'epithelioid')
  df2 = data.frame (barcode = tail (rownames(x)[order(x[,cnmf_module])], top_cells),
    subtype = 'sarcomatoid')
  rbind (df1, df2)
  })
cnmf_labelled_df = do.call (rbind,cnmf_labelled)

archp$epit_sarc = 0
archp$epit_sarc = cnmf_labelled_df$subtype[match (rownames(archp@cellColData), cnmf_labelled_df$barcode)]
archp$epit_sarc = paste0(archp$Sample3, '__',archp$epit_sarc)
table (archp$epit_sarc) #, archp$Clusters)

# Check groups have been assigned correctly ####
pdf()
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = cnmf_module, 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
p1 <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "GeneScoreMatrix", 
    name = 'SOX9', 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
p2 <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "GeneScoreMatrix", 
    name = 'SOX6', 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
groups_p = plotEmbedding(
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = 'epit_sarc', 
    embedding = "UMAP"#,
#    pal = rev (palette_deviation),
 #   imputeWeights = getImputeWeights(archp)
)

dev.off()

pdf (file.path('Plots','sarcomatoid_score_groups_plots2.pdf'), width = 20, height = 20)
print (wrap_plots (p,p1,p2, groups_p, ncol = 4))
dev.off()








#### Check regions with SOX9 and SOX6 ####
tf_name = 'SOX6'
TFmatch = getMatches (archp)
TFpositions = getPositions (archp)
x= 'P5__epithelioid'
tmp = readRDS (file.path('PeakCalls',paste0(x, '-reproduciblePeaks.gr.rds')))
TFmatch_sub = subsetByOverlaps (TFmatch, tmp)
TFmatch_assay = assay (TFmatch_sub)
motif_peaks = rowRanges (TFmatch_sub)[which(TFmatch_assay[,which(grepl(tf_name, colnames(TFmatch_sub)))])]
#write.csv (as.character(motif_peaks), 'SOX9_peaks_P5_sarcomatoid.csv')
#head(table (motif_peaks$nearestGene)[order(-table (motif_peaks$nearestGene))],50)
motif_peaks=as.data.frame(motif_peaks, row.names=NULL)
write.table (motif_peaks[,c(1:4)], file.path('SOX9_peaks_P5_sarcomatoid.bed'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
#'RUNX2' %in% unique(motif_peaks$nearestGene)
#table (TFMatch_assay)
head (motif_peaks,100)
motif_peaks = rowRanges (TFmatch_sub)[which(TFmatch_assay[,which(grepl(tf_name, colnames(TFmatch_sub)))])]
TFpositions_sub = TFpositions[[grep (tf_name, names(TFpositions))]]
TFpositions_sub = TFpositions_sub[unique(queryHits(findOverlaps (TFpositions_sub, motif_peaks)))]
head (as.character(TFpositions_sub[order(-TFpositions_sub$score)]))

### Check most co-occurring TFs with SOX9 / SOX6 ####
# # Get deviation matrix ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]] 
rownames (mMat) = rowData (mSE)$name

# Filter by RNA expression ####
metaGroupName = 'sampleID'
active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)

tf_name = 'SOX9'
peaksets = c('P5__sarcomatoid','P5__epithelioid')
cooc_TF_l = list()
for (peakset in peaksets)
{
tmp = readRDS (file.path('PeakCalls',paste0(peakset, '-reproduciblePeaks.gr.rds')))
tmp= tmp[tmp$score > 200,]
TFmatch_sub = subsetByOverlaps (TFmatch, tmp)
colnames(TFmatch_sub) = gsub ('_.*','',colnames(TFmatch_sub))
colnames(TFmatch_sub) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames(TFmatch_sub))
TFmatch_assay = assay (TFmatch_sub[,active_TFs])
motif_peaks = rowRanges (TFmatch_sub)[which(TFmatch_assay[,tf_name])]
TFmatch_sub = subsetByOverlaps (TFmatch_sub, motif_peaks)
cooc_TFs = colSums (TFmatch_assay)
cooc_TFs = proportions(cooc_TFs)
cooc_TFs = data.frame (cooc = cooc_TFs, TF = names(cooc_TFs))
cooc_TFs$TF = factor (cooc_TFs$TF, levels = cooc_TFs$TF[order (-cooc_TFs$cooc)])
cooc_TFs$peakset = peakset
cooc_TF_l[[peakset]] = cooc_TFs
}
cooc_TF_df = do.call (rbind, cooc_TF_l)
bp = ggplot (cooc_TF_df, aes (x =TF, y = cooc, fill = peakset)) + geom_bar (stat = 'identity', position = 'dodge') + gtheme
pdf (file.path ('Plots','cooc_TF.pdf'), width=70)
bp
dev.off()





# Assign high vs low sc-S cells per sample ####
metaGroupName = 'sarc_score'
mods = as.data.frame (t(scale(t(archp@cellColData[,grep('cNMF',colnames(archp@cellColData))]))))
#mods = as.data.frame (archp@cellColData[,grep('cNMF',colnames(archp@cellColData))])
mods$Sample = archp$Sample2
mods = mods[order(-mods$cNMF20),]

sams = unique (archp$Sample2)
table (archp$Sample2)
sams = c('P10','P12','P23','P4','P5','P8','P14')
sarc_score = lapply (sams, function(sam) 
  {
  x = mods[mods$Sample == sam,'cNMF20', drop=F]
  barcodes_high = data.frame(
    barcode = rownames(x)[x[[1]] > quantile(x[[1]],0.6)],
    score='high')
  barcodes_low = data.frame(
    barcode= rownames(x)[x[[1]] < quantile(x[[1]],0.4)],
    score = 'low')
  rbind (barcodes_high, barcodes_low)
  })
sarc_score_df = do.call (rbind, sarc_score)
archp$sarc_score = 'mid'
archp$sarc_score[match(sarc_score_df$barcode, rownames(archp@cellColData))] = sarc_score_df$score

pdf (file.path ('Plots','sarc_score_groups_scaled.pdf'))
 umap_p1 = plotEmbedding (
  ArchRProj = archp, 
  colorBy = "cellColData",
   name = 'sarc_score', embedding = "UMAP")
umap_p1 
dev.off()

# Assign sarcomatoid score using UMAP clustering for P1 and P23
archp$sarc_score_cluster = ''
archp$sarc_score_cluster[archp$Clusters %in% c('C5')] = 'SOX9_high'
archp$sarc_score_cluster[archp$Clusters %in% c('C3')] = 'SOX9_low'
archp$sarc_score_cluster[archp$Sample == 'P1'] = 'SOX9_high'
archp$sarc_score_sample = paste0(archp$sarc_score_cluster, archp$Sample)

# NEW CLUSTERING
archp$sarc_score_cluster2 = ''
archp$sarc_score_cluster2[archp$Clusters2 %in% c('C4')] = 'SOX9_low'
archp$sarc_score_cluster2[archp$Clusters2 %in% c('C9')] = 'SOX9_high'
archp$sarc_score_cluster2[archp$Clusters2 %in% c('C22')] = 'SOX9_low'
archp$sarc_score_cluster2[archp$Clusters2 %in% c('C23')] = 'SOX9_high'
archp$sarc_score_sample2 = paste0(archp$sarc_score_cluster2, archp$Sample)


# Check SOX9 footprinting ####
archp = archp[archp$sarc_score_sample2 %in% c('SOX9_lowP23','SOX9_highP23','SOX9_lowP11','SOX9_highP11')]
metaGroupName='sarc_score_sample2'
# archp2 = archp
# archp = archp[archp$sarc_score != 'mid']
archp <- addGroupCoverages (ArchRProj = archp, groupBy = metaGroupName)
motifPositions <- getPositions (archp)

motifs <- c("SOX9", "SOX6",'RUNX2','RUNX1','SNAI2','GLIS3')
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

#markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
#markerMotifs

# ### Run peak calling ####
# metaGroupName = "Sample"
# force=FALSE
# peak_reproducibility='1' # Set to 1 to better identify tumor heterogeneity
# if(!all(file.exists(file.path('PeakCalls', paste0(unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds')))) | force) source (file.path('..','..','git_repo','utils','callPeaks.R'))

#sams = c('P10','P12',
sams = 'P23'#),'P4','P5')

#for (sam in sams)
#  {
  #peaks_sample = readRDS (file.path ('PeakCalls',paste0(sam,'-reproduciblePeaks.gr.rds')))
  #motifPositions_sample = lapply (motifPositions, function(x) x[queryHits(findOverlaps(x, peaks_sample))])  
  seFoot <- getFootprints(
    ArchRProj = archp, 
    flank = 1000,
    #positions = motifPositions_sample[markerMotifs], 
    positions = motifPositions[markerMotifs], 
    groupBy = metaGroupName
  )
  
plotFootprints(
seFoot = seFoot[,grepl('P11', colnames(seFoot))],
ArchRProj = archp, 
flank = 1000,
normMethod = "Subtract",
plotName = 'Sarcomatoid TFs',
addDOC = FALSE, height=4.5, width=3,
smoothWindow = 25)
  #}

















markers = c('RUNX1','RUNX2','RUNX3','SOX9','VIM','SERPINE2','SNAI2')
metaGroupName='sarc_score'
#archp = addImputeWeights (archp)
sams = c('P10','P12','P23','P4','P5','P8','P14')
p2_l=list()
for (sam in sams)
  {
  p2_l[[sam]] <- plotGroups(
    ArchRProj = archp[archp$Sample == sam & archp$sarc_score != 'mid'], 
    groupBy = metaGroupName, 
    colorBy = "GeneScoreMatrix", 
    name = markers,
    plotAs = "violin",
    alpha = 0.4,
    getImputeWeights=NULL,
    addBoxPlot = TRUE#,
    #pal = palette_tnk_cells
   )
  }
p2_l2 = unlist (p2_l, recursive=F)
pdf (file.path ('Plots',paste0(sam,'_RUNXs_genescores.pdf')),14,14)
print (wrap_plots (p2_l2, ncol=length(markers)))
dev.off()



markers = c('RUNX1','RUNX2','RUNX3','SOX9','VIM','SERPINE2','SNAI2')
metaGroupName='Clusters'
#archp = addImputeWeights (archp)
sams = c('P10','P12','P23','P4','P5','P8','P14')
p2_l=list()
for (sam in sams)
  {
  p2_l[[sam]] <- plotGroups(
    ArchRProj = archp[archp$Sample == sam & archp$sarc_score != 'mid'], 
    groupBy = metaGroupName, 
    colorBy = "GeneScoreMatrix", 
    name = markers,
    plotAs = "violin",
    alpha = 0.4,
    getImputeWeights=NULL,
    addBoxPlot = TRUE#,
    #pal = palette_tnk_cells
   )
  }
p2_l2 = unlist (p2_l, recursive=F)
pdf (file.path ('Plots',paste0(sam,'_RUNXs_genescores_clusters.pdf')),14,14)
print (wrap_plots (p2_l2, ncol=length(markers)))
dev.off()


# Correlate module scores with TFs ####
metaGroupName = "sarc_score"
#archp$celltype2 = archp$Clusters_H
force=FALSE
# source (file.path('..','..','git_repo','utils','DAG.R'))






  ### Generate feature plots of uncommited state from Bueno ####
force = FALSE
uncommitted = list (uncommitted = c('MEST','HHIP','DNMT3A','TET1'))
archp@cellColData = archp@cellColData[,!grepl ('uncommitted',colnames(archp@cellColData), ignore.case=T)]
archp = addModuleScore (
    ArchRProj = archp,
    useMatrix = 'GeneScoreMatrix',
    name = '',
    features = uncommitted,
    nBin = 25,
    nBgd = 100,
    seed = 1,
    threads = getArchRThreads(),
    logFile = createLogFile("addModuleScore")
  )
colnames (archp@cellColData) = gsub ('^\\.','',colnames(archp@cellColData))

archp = addImputeWeights (archp)
#celltype_markers = c('WT1','CALB2','RUNX2','TCF3','SOX9','SOX6','MESP1','HMGA1','TWIST1','SNAI2')
pdf()
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = 'uncommitted', 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path('Plots','uncommitted_feature_plots.pdf'), width = 10, height = 10)
p
dev.off()






# Compute co-occurrence of sarcomatoid TFs ####
motifMat = getPositions (archp)
matches = getMatches (archp)
matchesMat = assay (matches)
colnames (matchesMat) = gsub ('_.*','',colnames (matchesMat))
colnames (matchesMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames (matchesMat))

ko_tfs = c('BPTF','TEAD1','TEAD4','SOX9','TWIST1','TCF3','HMGA1','MEF2A','MEF2D','PITX1')
#matchesMat = matchesMat[,]
#matchesMat = matchesMat[rowSums (matchesMat) > 0,]

nfeat=5000
k=25
top_genes=50
library(readxl)
cnmf_spectra_unique_comb_full = as.list (read_excel( "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/MPM_naive_study/reproduction2/PM_scRNA_atlas/data/cnmf_per_compartment.xlsx", sheet = "Cms_full"))
cnmf_spectra_unique_comb_full = lapply (cnmf_spectra_unique_comb_full, function(x) na.omit (x[x != 'NA']))


p2g_corr = .2
  p2g = getPeak2GeneLinks(
      ArchRProj = archp,
      corCutOff = p2g_corr,
      resolution = 1,
      returnLoops = FALSE
  )
p2g$gene = getPeakSet(archp)$nearestGene[match(p2g$idxATAC, getPeakSet(archp)$idx)]

p2g_tf = lapply (lapply (cnmf_spectra_unique_comb_full, function(y) head(y,50)), function(genes_in_module) p2g$idxATAC[p2g$gene %in% genes_in_module])
tf_ov = list()
for (ko_tf in ko_tfs)
  {
  message (paste0('compute overlap peaks in modules for TF: ', ko_tf))
  tf_ov[[ko_tf]] = unlist(lapply (p2g_tf, function(peaks_in_module)
    sum(matchesMat[rowRanges(matches)$idx %in% peaks_in_module, ko_tf]) /
     length(matchesMat[rowRanges(matches)$idx %in% peaks_in_module,ko_tf])))
  }

tf_df = do.call (cbind, tf_ov)

tf_scaled = scale(tf_df)
tf_scaled[is.na(tf_scaled)] = 0
pdf (file.path ('Plots','tf_target_genes_peaks_overlap.pdf'))
Heatmap (tf_scaled)
dev.off()

# Find genes with correlated peaks enriched in crop-seq TFs ####
tf_target_genes = list()

p2g_corr = -Inf
  p2g = getPeak2GeneLinks(
      ArchRProj = archp,
      corCutOff = p2g_corr,
      resolution = 1,
      returnLoops = FALSE
  )
p2g$gene = getPeakSet(archp)$nearestGene[match(p2g$idxATAC, getPeakSet(archp)$idx)]
p2g = p2g[!is.na(p2g$gene),]

for (ko_tf in ko_tfs)
  {
  message (paste0('find target genes for TF: ', ko_tf))
  tf_peaks = matchesMat[,ko_tf]
  gene_peaks = table (p2g$gene[match (rowRanges(matches)$idx[tf_peaks], p2g$idxATAC)])
  gene_peaks = names (head (gene_peaks[order(-gene_peaks)],100))
  tf_target_genes[[ko_tf]] = gene_peaks
  }

### p2g genes only have a little more than 2,000 genes. Using full peakset instead... #####
#ps = getPeakSet (archp)
#ps2 = ps[!is.na(ps$nearestGene), ]
for (ko_tf in ko_tfs)
  {
  message (paste0('find target genes for TF: ', ko_tf))
  tf_peaks = table (rowRanges(matches)$nearestGene[matchesMat[,ko_tf]])
  #gene_peaks = table (rowRanges(matches)$nearestGene[match(ps$idx, tf_peaks)])
  tf_peaks = names (head (tf_peaks[order(-tf_peaks)],Inf))
  tf_target_genes[[ko_tf]] = tf_peaks
  }

saveRDS (tf_target_genes, 'genes_with_tf_peaks.rds')


# Export positions of SOX9 motifs ####
matches = getMatches (archp)
matchesMat = assays(matches)[[1]]
colnames (matchesMat) = gsub ('_.*','', colnames (matchesMat))
colnames (matchesMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames (matchesMat))
write.table (as.data.frame (rowRanges (matches)[matchesMat[,'SOX9']], row.names = NULL), file = 'SOX9_motifs.bed', row.names=FALSE, col.names=FALSE, sep='\t', quote=F)


write.table (as.data.frame (rowRanges (matches)[matchesMat[,'NFYB']], row.names = NULL), file = 'NFYB_motifs.bed', row.names=FALSE, col.names=FALSE, sep='\t', quote=F)

## Find SOX9 in peaks of Cm17 module genes
ps1 = rowRanges (matches)[matchesMat[,'SOX9']]
ps2 = ps1[ps1$nearestGene %in% cnmf_spectra_unique_comb_full[['Cm17']]]
write.csv (ps2, 'SOX9_Cm17_hits.csv')









markers = c('RUNX1','RUNX2','SOX9','VIM','SERPINE2','SNAI2',
  'AXL','HAPLN1','SOX6','LAMC2','FBN1','LOXL2','CAV1')
metaGroupName='Clusters2'
#archp = addImputeWeights (archp)
  p2_l <- plotGroups(
    ArchRProj = archp, 
    groupBy = metaGroupName, 
    colorBy = "GeneScoreMatrix", 
    name = markers,
    plotAs = "violin",
    alpha = 0.4,
    getImputeWeights=NULL,
    addBoxPlot = TRUE#,
    #pal = palette_tnk_cells
   )

pdf (file.path ('Plots','markers_genescores_boxplots.pdf'),14,14)
wrap_plots (p2_l)
dev.off()













### use TF activity to cluster P23 ####


# # Try with PCA on mMat ####
library (ggridges)
library (ggplot2)
library (viridis)
library (tidyr)
#library(hrbrthemes)


# Add metacolumns of average TF modules activity ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = scale(mMat[,archp$Sample == 'P23'])
#mMat = scale(as.matrix(mMat))#[selected_TF,])

# # # Filter by RNA expression ####
metaGroupName = 'sampleID'
min_exp = .1
active_TFs = exp_genes (srt, rownames(mMat), min_exp = min_exp, metaGroupName)
mMat = t(scale (mMat[active_TFs, ]))
library(uwot)
library(ggplot2)
library(patchwork)

#-------------------------
# 0. RUN PCA
#-------------------------
set.seed(123)  # for reproducibility
p <- prcomp(t(mMat), center = TRUE, scale. = TRUE)
plot_df <- data.frame(p$x, celltype = t(mMat[c('SOX9','SNAI2','RUNX2'),]))
# df_sub <- plot_df[plot_df$celltype %in% 
#                     c("Mono_CD14","TAM_interstitial","TAM_TREM2",
#                       "TAM_MARCO","TAM_CXCLs"), ]

# plot_df <- data.frame(
#   PC1 = p$x[,1],
#   PC2 = p$x[,2],
#   SOX9 = plot_df[,'celltype.SOX9'],
#   SNAI2 = plot_df[,'celltype.SNAI2'],
#   RUNX2 = plot_df[,'celltype.RUNX2']
# )
df_sub= plot_df
# df_sub <- plot_df[plot_df$celltype %in%
#                     c("Mono_CD14","TAM_interstitial","TAM_TREM2",
#                       "TAM_MARCO","TAM_CXCLs"), ]


#-------------------------
# 1. MAIN SCATTER
#-------------------------
library(viridis)
which.max(abs (cor (df_sub$celltype.SOX9, p$x,method = 'pearson')))
pca_cor[order(pca_cor),drop=F]
p_scatter <- ggplot(df_sub, aes(PC1, PC2, color = celltype.SNAI2)) +
  geom_point(alpha = 0.4, size = 0.1) +
  #geom_density_2d(size = 0.4) +
  scale_color_viridis_c(option = "viridis", direction = 1) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )

p_scatter2 <- ggplot(df_sub, aes(PC1, PC2, color = celltype.RUNX2)) +
  geom_point(alpha = 0.4, size = 0.1) +
  #geom_density_2d(size = 0.4) +
  scale_color_viridis_c(option = "viridis", direction = 1) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )

  
# Compute percentiles
#q10 <- quantile(plot_df$celltype.SOX9, 0.20, na.rm = TRUE)
#q90 <- quantile(plot_df$celltype.SOX9, 0.80, na.rm = TRUE)

# Subset the dataframe
#df_sub2 <- subset(plot_df, celltype.SOX9 >= q10 & celltype.SOX9 <= q90)

p_scatter3 <- ggplot(df_sub, aes(PC1, PC2, color = celltype.SOX9)) +
  geom_point(alpha = 0.4, size = 0.1) +
  #geom_density_2d(size = 0.4) +
  scale_color_viridis_c(option = "viridis", direction = 1) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )

#-------------------------
# 4. SAVE PDF
#-------------------------
pdf(file.path("Plots","TF_activity_pca.pdf"),
    width = 8, height = 2)

wrap_plots (p_scatter, p_scatter2, p_scatter3, ncol=3)
dev.off()

library(uwot)

set.seed(123)

# Transpose so cells are rows (required for UMAP)
mat_t <- t(scale(mMat))

# Run UMAP
um <- umap(as.matrix(mat_t), n_neighbors = 30, min_dist = 0.3, metric = "euclidean",n_components=10)
cor (um, mat_t[,'SOX9'], method='spearman')
plot_df <- data.frame(um, celltype = t(mMat[c('SOX9','SNAI2','RUNX2'),]))


p_scatter <- ggplot(plot_df, aes(X1, X2, color = celltype.SNAI2)) +
  geom_point(alpha = 0.4, size = 0.1) +
  #geom_density_2d(size = 0.4) +
  scale_color_viridis_c(option = "viridis", direction = -1) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )

p_scatter2 <- ggplot(plot_df, aes(X1, X2, color = celltype.RUNX2)) +
  geom_point(alpha = 0.4, size = 0.1) +
  #geom_density_2d(size = 0.4) +
  scale_color_viridis_c(option = "viridis", direction = -1) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )

# Compute percentiles
q10 <- quantile(plot_df$celltype.SOX9, 0.10, na.rm = TRUE)
q90 <- quantile(plot_df$celltype.SOX9, 0.90, na.rm = TRUE)

# Subset the dataframe
df_sub <- subset (plot_df, celltype.SOX9 >= q10 & celltype.SOX9 <= q90)

# Plot
p_scatter3 <- ggplot(df_sub, aes(X1, X2, color = celltype.SOX9)) +
  geom_point(alpha = 0.4, size = 0.1) +
  #geom_density_2d(size = 0.4) +
  scale_color_viridis_c(option = "viridis", direction = -1) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )

pdf(file.path ("Plots","TF_activity_umap.pdf"),
    width = 10, height = 3)

wrap_plots (p_scatter, p_scatter2, p_scatter3, ncol=3)
dev.off()





### Try to order cells by scS-score and then look for SOX9 expression / deviation and maybe use to re-run chromBPnet ####


  # Get deviations ####
  if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
  mMat = assays (mSE)[[1]]
  rownames (mMat) = rowData (mSE)$name
  mMat = scale (mMat[selected_TF,])

  # Get genescore ####
  if (!exists('gsSE')) gsSE = fetch_mat (archp, 'GeneScore')
  gsMat = assays (gsSE)[[1]]
  rownames (gsMat) = rowData (gsSE)$name
  gsMat = scale(gsMat[selected_TF,])
  
  mMat = t(mMat)
  mMat = lapply (sams, function(x) mMat[archp_meta$Sample3 == x,])
  names (mMat) = sams

  # Get cnmf modules ####
  cnmf_mat = archp@cellColData[,grep ('cNMF',colnames(archp@cellColData))]
  cnmf_mat = as.data.frame (t(scale (t(cnmf_mat))))
  cnmf_mat = lapply (sams, function(x) cnmf_mat[archp_meta$Sample3 == x, ])
  names (cnmf_mat) = sams

  gsMat = scale(t(gsMat))
  gsMat = lapply (sams, function(x) gsMat[archp_meta$Sample3 == x,])
  names (gsMat) = sams

  # Average mats along sarc module score ####
  library(zoo)

  bin_width <- 40   # Number of observations per bin
  overlap <- 40    
  mMat_ordered = lapply (sams, function(sam) mMat[[sam]][order(cnmf_mat[[sam]][,sarc_module]),])
  names(mMat_ordered) = sams
  mMat_ordered_avg = lapply (sams, function (sam) as.data.frame (lapply (as.data.frame (mMat_ordered[[sam]]), function(x) {
    zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
  })))
  names (mMat_ordered_avg) = sams

  gsMat_ordered = lapply (sams, function(sam) gsMat[[sam]][order(cnmf_mat[[sam]][,sarc_module]),])
  names(gsMat_ordered) = sams
  gsMat_ordered_avg = lapply (sams, function (sam) 
    {
    tmp_df = as.data.frame (lapply (as.data.frame (gsMat_ordered[[sam]]), function(x) {
    zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
    }))
    tmp_df[is.na(tmp_df)] = 0 # HIC2 returns NaN for one or some samples
    tmp_df
    })

  names (gsMat_ordered_avg) = sams

  cnmfMat_ordered = lapply (sams, function(sam) cnmf_mat[[sam]][order(cnmf_mat[[sam]][,sarc_module]),])
  names(cnmfMat_ordered) = sams
  cnmfMat_ordered_avg = lapply (sams, function (sam) as.data.frame (lapply (as.data.frame (cnmfMat_ordered[[sam]]), function(x) {
    zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
  })))
  names (cnmfMat_ordered_avg) = sams

  m_cor = lapply (sams, function(sam) 
    {
    cor_tmp = as.data.frame (cor (mMat_ordered_avg[[sam]], cnmfMat_ordered_avg[[sam]][,sarc_module, drop=F], method = 'spearman'))
    colnames(cor_tmp)[1] = 'score'
    cor_tmp$sample = sam
    cor_tmp$TF = rownames(cor_tmp)
    cor_tmp
    })
  m_cor_df = do.call (rbind, m_cor)
  m_cor_df$type = 'activity'
  m_cor_levels = m_cor_df %>% group_by (TF) %>% summarise (median_value = median(score)) %>% arrange(-median_value)

  gs_cor = lapply (sams, function(sam) 
    {
    cor_tmp = as.data.frame(cor (gsMat_ordered_avg[[sam]], cnmfMat_ordered_avg[[sam]][,sarc_module, drop=F], method = 'spearman'))
    colnames(cor_tmp)[1] = 'score'
    cor_tmp$sample = sam
    cor_tmp$TF = rownames(cor_tmp)
    cor_tmp
    })
  gs_cor_df = do.call (rbind, gs_cor)
  gs_cor_df$type = 'genescore'





### Make module score of SOX9 co-expressed genes from bulk ####
sox9_module = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/bulkRNA_meso/SOX9_gene_correlations.csv')
markerMotifs = getFeatures (archp, useMatrix = "GeneScoreMatrix")
sox9_module = sox9_module[order(-sox9_module$bueno),]  
sox9_module = sox9_module[sox9_module[,1] %in% markerMotifs,]  
sox9_module = list (SOX9 = head(sox9_module$X, 100))
  archp = addModuleScore (
      ArchRProj = archp,
      useMatrix = 'GeneScoreMatrix',
      name = '',
      features = sox9_module,
      nBin = 25,
      nBgd = 100,
      seed = 1,
      threads = getArchRThreads(),
      logFile = createLogFile("addModuleScore")
    )
colnames (archp@cellColData) = gsub ('^\\.','',colnames(archp@cellColData))
  
pdf()
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = 'SOX9', 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path('Plots','bulkRNA_SOX9_module_feature_plots.pdf'), width = 20, height = 20)
#print (wrap_plots (p2, ncol = 4))
print (wrap_plots (p, ncol= 4))
dev.off()

# Try with SOX9 regulon from scRNA ####
scenic_genes = read.csv (file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scrna/SCENIC/vg_5000_mw_tss500bp/tumor_programs','motifs.csv'), header=T, skip=1)
genes <- sapply(seq(nrow(scenic_genes)), function(x) unlist(regmatches(scenic_genes[x,9], gregexpr("(?<=\\(')[A-Za-z0-9_-]+", scenic_genes[x,9], perl = TRUE))))
names (genes) = scenic_genes[,1]
names (genes) = paste0('SCENIC_',names(genes))
archp = addModuleScore (
      ArchRProj = archp,
      useMatrix = 'GeneScoreMatrix',
      name = '',
      features = genes['SCENIC_SOX9'],
      nBin = 25,
      nBgd = 100,
      seed = 1,
      threads = getArchRThreads(),
      logFile = createLogFile("addModuleScore")
    )
colnames (archp@cellColData) = gsub ('^\\.','',colnames(archp@cellColData))
  
pdf()
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = 'SCENIC_SOX9', 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()

pdf (file.path('Plots','SCENIC_SOX9_module_feature_plots.pdf'), width = 20, height = 20)
#print (wrap_plots (p2, ncol = 4))
print (wrap_plots (p, ncol= 4))
dev.off()

### Ridge plot of SOX9 regulon ####
library(ggridges)
ccomp = archp@cellColData[,c('SCENIC_SOX9','Sample3')]
p_ridge <- ggplot(ccomp, aes(x = SCENIC_SOX9, y = Sample3, fill = Sample3)) +
  geom_density_ridges(scale = 3, alpha = 0.7, color = NA) +
  scale_fill_manual (values = palette_sample) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y  = element_text(size = 8),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  gtheme_no_rot

pdf(file.path("Plots","SCENIC_SOX9_ridge_plots.pdf"),
    width = 7, height = 2)
p_ridge
dev.off()

ccomp_P23 = ccomp[ccomp$Sample3 == 'P23',]
hp = ggplot (ccomp_P23, aes (x = SCENIC_SOX9)) + geom_histogram()
pdf(file.path("Plots","SCENIC_SOX9_histogram.pdf"),
    width = 7, height = 2)
hp
dev.off()

alpha <- 0.05  # total tail probability
lower <- quantile(ccomp_P23$SCENIC_SOX9, alpha / 2, na.rm = TRUE)
upper <- quantile(ccomp_P23$SCENIC_SOX9, 1 - alpha / 2, na.rm = TRUE)
ccomp_P23$breaks = 'not_selected'
ccomp_P23$breaks[ccomp_P23$SCENIC_SOX9 <= lower] = 'SOX9_regulon_low'
ccomp_P23$breaks[ccomp_P23$SCENIC_SOX9 >= upper] = 'SOX9_regulon_high'

bp = ggplot(ccomp_P23, aes(x = SCENIC_SOX9, fill = breaks)) +
  geom_histogram(bins = 30, color = "black")
pdf(file.path("Plots","SCENIC_SOX9_histogram.pdf"),
    width = 7, height = 2)
bp
dev.off()





cor (archp$SCENIC_SOX9, archp$SOX9)


