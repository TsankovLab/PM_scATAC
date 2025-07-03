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
#archp$Sample3[grep ('normal',archp$Sample3)] = 'normal'

# Load RNA ####
srt = readRDS (file.path('..','scrna','srt.rds'))
srt$sampleID3[srt$sampleID3 %in% c('HU37','HU62')] = 'normal_pleura'
sarc_order = read.csv (file.path('..','scrna','cnmf20_sarcomatoid_sample_order.csv'), row.names=1)
sarc_order = sarc_order[! sarc_order$sampleID %in% c('HU37','HU62'),]
sarc_order = rbind (data.frame (sampleID = 'normal_pleura', x = -1),sarc_order) 
#archp$Sample2 = factor (archp$Sample2, levels = sarc_order$sampleID)


archp_NN = archp[!archp$Sample3 %in% c('normal1','normal2','normal3')]

pdf ()
umap_p1 = plotEmbedding (ArchRProj = archp_NN, labelMeans = F, 
  colorBy = "cellColData", name = "Sample2", 
  pal = palette_sample,
   embedding = "UMAP")
umap_p2 = plotEmbedding (ArchRProj = archp_NN, labelMeans = T, 
  colorBy = "cellColData", name = "Clusters",
   embedding = "UMAP")
dev.off()

pdf (file.path ('Plots','sample_clusters_umap.pdf'))
umap_p1
umap_p2
dev.off()

# Run genescore DAG ####
archp$status = ifelse (archp$Sample3 == 'normal_pleura','normal','tumor')
metaGroupName = "status"
force = FALSE
if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source (file.path('..','..','git_repo','utils','DAG.R'))

archp = addImputeWeights (archp)
celltype_markers = c('WT1','CALB2','RUNX2','TCF3','SOX9','SOX6','MESP1','HMGA1','TWIST1','SNAI2')
pdf()
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "GeneScoreMatrix", 
    name = celltype_markers, 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path('Plots','marker_genes_feature_plots.pdf'), width = 20, height = 20)
print (wrap_plots (p, ncol = 4))
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
average_by_group <- as.data.frame(archp@cellColData) %>%
  group_by(Sample3) %>%
  summarise(Average = mean(cNMF20))
average_by_group = average_by_group[order(-average_by_group$Average),]
sample_sarc_order = factor (archp$Sample3, levels = average_by_group$Sample3)



# ### Compare bins malignants normal ####
# run_bin_analysis = FALSE

# if (run_bin_analysis)
#   {
#   # Load fragments
#   fragments = unlist (getFragmentsFromProject (
#        ArchRProj = archp))
  
#   ws = 1e6
#   ss = 2e5
#   if (!file.exists (paste0('bins_',ws,'_ss_',ss,'.rds')))
#     {
#     blacklist = toGRanges(paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/blacklisted_regions/ENCODE_blacklist/',"hg38-blacklist.v2.bed")) # taken from https://github.com/Boyle-Lab/Blacklist/tree/master/lists
#     windows = makeWindows(genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist,
#       windowSize = ws, slidingSize = ss)
#     saveRDS (windows, paste0('bins_',ws,'_ss_',ss,'.rds'))
#     } else {
#     windows = readRDS (paste0('bins_',ws,'_ss_',ss,'.rds'))  
#     }
  
#   flt_celltypes = 100
  
#   #archp$celltype_sample = paste0(archp$celltype, '-', archp$Sample2)
#   metaGroupName = 'Clusters'
  
#   #fragments$RG = as.character(fragments$RG)  
#   barcode_metaGroup = as.data.frame (archp@cellColData[,c(metaGroupName, 'TSSEnrichment','nFrags')])
#   colnames (barcode_metaGroup)[colnames (barcode_metaGroup) == metaGroupName] = 'metaGroup'
#   barcode_metaGroup$barcode = rownames(barcode_metaGroup)
#   celltype_bins = lapply (unique(archp@cellColData[,metaGroupName]), function(mg) 
#     {
#     high_quality_cells = barcode_metaGroup[barcode_metaGroup$metaGroup == mg,]
#     high_quality_cells = high_quality_cells[order (-high_quality_cells$TSSEnrichment),] # ordering by TSSEnrichment seem to work the best
#     high_quality_cells = head (high_quality_cells$barcode, flt_celltypes)
#     fragments_group = fragments[fragments$RG %in% high_quality_cells]
#     #names (fragments_cell) = NULL 
#     #fragments_cts = sample (unique(fragments_ct$RG), 200)
#     fr_ov = countOverlaps (windows, fragments_group)
#   #  fr_ov = fr_ov / sum (archp$nFrags[archp$celltype_sample == x])
#     #(fr_ov / sum (fr_ov)) * 1000
#     })
#   celltype_bins = do.call (cbind, celltype_bins)
#   colnames (celltype_bins) = sapply (unique(archp@cellColData[,metaGroupName]), function(x) unlist(strsplit (x,'-'))[1])
#   head (celltype_bins)
  
#   celltype_bins_cor = cor (celltype_bins, method = 'spearman')
#   ha = HeatmapAnnotation (sample = sapply (unique(archp@cellColData[,metaGroupName]), function(x) unlist(strsplit (x,'-'))[2]), which='row')
#   binH = Heatmap (celltype_bins_cor, col = viridis::rocket(100),# row_names_gp= gpar (fontsize=6), 
#     #column_names_gp= gpar (fontsize=6), 
#     right_annotation = ha,
#     cluster_rows = T, #row_km = 2, column_km = 2,
#     clustering_distance_rows = 'pearson', 
#     clustering_distance_columns = 'pearson',
#     column_title = paste('Binned',ws))
  
#   png (paste0('Plots/binned_fragments_binned_',ws,'_celltype_heatmap.png'),width=600, height=500)
#   binH
#   dev.off()
#   }


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



# Make data.frame of deviation difference and expression difference between normal and tumors ####
force = F
tumor_sams = c('P1','P10','P11','P12','P13','P14','P23','P3','P4','P5','P8')
normal_sams = c('normal1')
tumor_sams_rna = c('P1','P11','P12','P13','P14','P3','P4','P5','P8')
normal_sams_rna = 'normal_pleura'

if(!file.exists('selected_TF.rds') | force)
  {
  # Get activity matrix and format
  metaGroupName = 'Sample3'
  if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
  archp_meta = as.data.frame (archp@cellColData)
  mMat = assays (mSE)[[1]]
  rownames (mMat) = rowData (mSE)$name
  mMat_agg = as.data.frame (t(mMat))
  mMat_agg$metaGroup = as.character (archp_meta[,metaGroupName])
  mMat_agg = aggregate (.~ metaGroup, mMat_agg, mean)
  rownames (mMat_agg) = mMat_agg[,1]
  mMat_agg = mMat_agg[,-1]
  mMat_agg = t(mMat_agg)
  mMat_agg = mMat_agg[rownames(mMat_agg) %in% rownames(srt),]
  
  # Get genescore matrix and format
  metaGroupName = 'Sample3'
  if (!exists('gsSE')) gsSE = fetch_mat (archp, 'GeneScore')
  archp_meta = as.data.frame (archp@cellColData)
  gsMat = assays (gsSE)[[1]]
  rownames (gsMat) = rowData (gsSE)$name
  gsMat = scale(gsMat[rownames(gsMat) %in% rownames(mMat_agg),])
  gsMat_agg = as.data.frame (t(gsMat))
  gsMat_agg$metaGroup = as.character (archp_meta[,metaGroupName])
  gsMat_agg = aggregate (.~ metaGroup, gsMat_agg, mean)
  rownames (gsMat_agg) = gsMat_agg[,1]
  gsMat_agg = gsMat_agg[,-1]
  gsMat_agg = t(gsMat_agg)
  
  # Get expression matrix and format
  metaGroupName = 'sampleID3'
  DefaultAssay (srt) = 'RNA'
  sample_names_rna = c('P1','P3','P4','P5','P8','P11','P12','P13','P14','normal_pleura')
  ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat_agg), group.by = metaGroupName)[[1]]) +1)
  colnames(ps) = gsub ('-','_',colnames(ps))
  ps = ps[, colnames(ps) %in% sample_names_rna]
  ps = ps[rownames(mMat_agg),]
  
# Combine matrices
  TF_diff_rna = data.frame (
    tumor_dev = apply (mMat_agg[,tumor_sams], 1, mean),
    normal_dev = mMat_agg[,normal_sams],
    tumor_gs = apply (gsMat_agg[,tumor_sams], 1, mean)[rownames(mMat_agg)],
    normal_gs =  gsMat_agg[,normal_sams][rownames(mMat_agg)],
    tumor_rna = apply (ps[,tumor_sams_rna],1,mean),
    normal_rna = ps[, normal_sams_rna])
    
  # Compute differences between tumor and normal per type    
  TF_diff_rna$dev_diff = TF_diff_rna$tumor_dev - TF_diff_rna$normal_dev
  TF_diff_rna$gs_diff = TF_diff_rna$tumor_gs - TF_diff_rna$normal_gs
  TF_diff_rna$rna_diff = TF_diff_rna$tumor_rna - TF_diff_rna$normal_rna
  
  
# Compare chromVAR deviations ####
library (presto)
mMat = scale (assays (mSE)[[1]])
rownames (mMat) = rowData (mSE)$name

all (colnames(mMat) == rownames(archp@cellColData))

# Set comparisons and run test for activity
dev_comparisons <- setNames(lapply(tumor_sams, function(x) c(x, normal_sams)), tumor_sams)

dev_res = lapply (dev_comparisons, function(x) 
  {
  tmp = wilcoxauc (mMat, y = archp$Sample3, groups_use = x)
  tmp = tmp[tmp$group == x[1],]
  tmp[tmp$logFC > 0,]
  })
pos_features = unique(unlist(lapply (dev_res, function(x) x$feature)))

dev_res_df = lapply (names (dev_comparisons), function(x) 
  {
  tmp = dev_res[[x]]
  rownames(tmp) = dev_res[[x]]$feature
  tmp = tmp[,c('padj'),drop=FALSE]
  tmp[pos_features,,drop=F]
  })
dev_res_df = do.call (cbind , dev_res_df)
# dev_res_df = dev_res_df[, !grepl ('feature',colnames(dev_res_df))]
#rownames (dev_res_df) = rowData (mSE)$name
colnames (dev_res_df) = names(dev_comparisons)
rownames (dev_res_df) = pos_features

# # Repeat using genescore
# gsMat = assays (gsSE)[[1]]
# rownames (gsMat) = rowData (gsSE)$name
# gsMat = scale (gsMat[rownames(gsMat) %in% rownames(mMat),])

# # Set comparisons and run test for genescore
# gs_res = lapply (dev_comparisons, function(x) 
#   {
#   tmp = wilcoxauc (gsMat, y = archp$Sample3, groups_use = x)
#   tmp = tmp[tmp$group == x[1],]
#   tmp[tmp$logFC > 0,]
#   })
# # b_df = data.frame(activity = gsMat['SOX9',], sample = archp$Sample3)
# # pdf (file.path('Plots','bxp.pdf'))
# # ggplot(b_df, aes (x = sample, y = activity)) + geom_boxplot()
# # dev.off()

# #pos_features = unique(unlist(lapply (gs_res, function(x) x$feature)))

# gs_res_df = lapply (names (dev_comparisons), function(x) 
#   {
#   tmp = gs_res[[x]]
#   rownames(tmp) = gs_res[[x]]$feature
#   tmp = tmp[,c('padj'),drop=FALSE]
#   tmp[pos_features,,drop=F]
#   })
# gs_res_df = do.call (cbind , gs_res_df)
# colnames (gs_res_df) = names(dev_comparisons)
# rownames (gs_res_df) = pos_features

# # rownames(gs_res_df) = gs_res_df[,1]
# # gs_res_df = gs_res_df[, !grepl ('feature',colnames(gs_res_df))]
# #rownames (dev_res_df) = rowData (mSE)$name

# # Repeat for scRNA  ##
# rna_comparisons <- setNames(lapply(tumor_sams_rna, function(x) c(x, normal_sams_rna)), tumor_sams_rna)

# # Set comparisons and run test for genescore
# rna_res = lapply (rna_comparisons, function(x) 
#   {
#   tmp = wilcoxauc (srt[rowData(mSE)$name,], group_by = 'sampleID3', groups_use = x)
#   tmp = tmp[tmp$group == x[1],]
#   tmp[tmp$logFC > 0,]
#   })
# #pos_features = unique(unlist(lapply (rna_res, function(x) x$feature)))

# rna_res_df = lapply (names (rna_comparisons), function(x) 
#   {
#   tmp = rna_res[[x]]
#   rownames(tmp) = rna_res[[x]]$feature
#   tmp = tmp[,c('padj'),drop=FALSE]
#   tmp[pos_features,,drop=F]
#   })
# rna_res_df = do.call (cbind , rna_res_df)
# colnames (rna_res_df) = names(rna_comparisons)
# rownames (rna_res_df) = pos_features

  
# Also test genetic and chromatin regulators ####
# rna_res_2 = lapply (rna_comparisons, function(x) 
#   wilcoxauc (srt, group_by = 'sampleID2', groups_use = x))
# chromatin_regulators = c('SETD5, ASH1L, CREBBP, PRDM2, KDM2B, KMT2D, EZH2, SETDB1')
# chromatin_regulators = unlist(strsplit (chromatin_regulators, ', '))
# genetic_drivers = c('BAP1, NF2, CDKN2B, CDKN2A, TP53, LATS2, SETD2, FAT4, PTCH1')
# chromatin_regulators2 = c('LATS1, DDX3X, ULK2, RYR2, CFAP45, SETDB1, DDX51, SF3B1, TRAF7, PTEN, RBFOX1, CSMD1, MTAP, TTC28, PCDH15, USH2A, CNTNAP2, DNAH1, KCNH7, PTK2, ROBO2, DLG2, PBRM1, PTPRD, ANTXR2, CTNNA3, LINGO2, LRP1B, PLCB1, UNC79, WWOX')
# chromatin_regulators2 = unlist (strsplit (chromatin_regulators2, ', '))
# genetic_drivers = unlist(strsplit (genetic_drivers, ', '))
# check_genes = unique (c(chromatin_regulators, genetic_drivers, chromatin_regulators2))

# rna_res_sum_p2 = apply (do.call (cbind, lapply (rna_res_2, function(x) x[match(check_genes, x$feature),c('padj'), drop=F])), 1, median)
# rna_res_sum_lfc2 = apply (do.call (cbind, lapply (rna_res_2, function(x) x[match(check_genes, x$feature),c('logFC'), drop=F])),1, median)
# rna_res_sum_ae2 = apply (do.call (cbind, lapply (rna_res_2, function(x) x[match(check_genes, x$feature),c('avgExpr'), drop=F])),1, median)
# rna_res_sum2 = data.frame (row.names = check_genes, RNA_pvalue_median = rna_res_sum_p2, RNA_lfc_median = rna_res_sum_lfc2, RNA_avExpr_median = rna_res_sum_ae2)
# write.csv (rna_res_sum2, 'normal_tumor_DEG_chromatin_reg_genetic_drivers.csv')

#occurence_filter = apply (rna_res_df, 1, function(x) sum (x < pval_threshold))
# rna_res_df_filtered = rna_res_df[occurence_filter > occurrence_threshold, ]
# tf_tumor_pos = rownames(TF_diff_rna)[TF_diff_rna$dev_diff > 0 & TF_diff_rna$rna_diff > 0]
# rna_res_df_filtered = rna_res_df_filtered[rownames(rna_res_df_filtered) %in% tf_tumor_pos,]
# rna_selected_TF = rownames(rna_res_df_filtered)
  
# Filter non-significant and low occurence 
pval_threshold = 0.05
occurrence_threshold = 5

# DEV + GS
dev_res_df_logic = dev_res_df < pval_threshold
dev_res_df_logic[is.na(dev_res_df_logic)] = FALSE
# gs_res_df_logic = gs_res_df < pval_threshold
# gs_res_df_logic[is.na(gs_res_df_logic)] = FALSE
# comb_res_df = dev_res_df_logic + gs_res_df_logic

# # DEV + RNA
# dev_res_df_logic = dev_res_df < pval_threshold
# dev_res_df_logic[is.na(dev_res_df_logic)] = FALSE
# rna_res_df_logic = rna_res_df < pval_threshold
# rna_res_df_logic[is.na(rna_res_df_logic)] = FALSE
# shared_features = intersect (rownames(dev_res_df_logic), rownames(rna_res_df_logic))
# rna_res_df_logic = rna_res_df_logic[shared_features,]
# comb_res_df = dev_res_df_logic[,colnames(rna_res_df_logic)] + rna_res_df_logic

comb_res_df <- matrix(as.numeric(dev_res_df_logic), nrow = nrow(dev_res_df_logic), ncol = ncol(dev_res_df_logic))
rownames(comb_res_df) = rownames(dev_res_df_logic)
#comb_res_df[comb_res_df == 1] = 0
#comb_res_df[comb_res_df == 2] = 1
selected_TF = rownames(comb_res_df) [rowSums (comb_res_df) >= occurrence_threshold]

metaGroupName = 'sampleID3'
ps = log2(as.data.frame (AverageExpression (
  srt[,srt$sampleID3 %in% c(tumor_sams_rna)], 
  features = selected_TF, 
  group.by = metaGroupName)[[1]]) +1)
min_exp = .5
active_TFs = rownames(ps)[apply (ps,1,function(x) any(x > min_exp))]
selected_TF = selected_TF[selected_TF %in% active_TFs]
tf_tumor_pos = rownames(TF_diff_rna)[TF_diff_rna$rna_diff > 0]
selected_TF = selected_TF[selected_TF %in% tf_tumor_pos]

# Order by mean logFC
tf_order = rownames(TF_diff_rna)[order (-rowMeans(TF_diff_rna[,c('dev_diff','rna_diff')]))]
selected_TF_ordered = tf_order[tf_order %in% selected_TF]


# Plot distribution of diff deviation and diff expression between tumor and normal
diff_line = 0
TF_diff_rna$label = ''
TF_diff_rna$label[match (selected_TF,rownames(TF_diff_rna))] = selected_TF
TF_diff_rna$color = TF_diff_rna$label != ''
TF_diff_rna$label_top = ''
TF_diff_rna$label_top[match (head (selected_TF_ordered,30),rownames(TF_diff_rna))] = head (selected_TF_ordered,30)
TF_diff_rna$alpha = 0.7
TF_diff_rna$alpha[match (selected_TF,rownames(TF_diff_rna))] = .9
TF_diff_rna$direction_rna = ifelse (apply (TF_diff_rna[,c('tumor_rna','normal_rna')], 1, function(x) which.max (x)) == 1, 'tumor','normal')
TF_diff_rna$highest_rna =apply (TF_diff_rna[,c('tumor_rna','normal_rna')], 1, function(x) max (x))
tf_diff_p = ggplot (TF_diff_rna, aes (x= dev_diff, y = rna_diff, size = highest_rna, fill = direction_rna, label = label_top)) + 
geom_point(aes(color=color, alpha = alpha, size=(highest_rna+1)^2),shape = 21, stroke=0.1) + # Color points based on x value
scale_color_manual(values = c('FALSE' = "grey", 'TRUE' = "black")) + # Customize colors
scale_fill_manual(values = c(normal = "grey20", tumor = "darkred")) +
#scale_fill_manual(values = c('FALSE' = "grey", 'TRUE' = "red")) + # Customize colors
geom_vline(xintercept = diff_line, linetype = "dashed", color = "grey44", linewidth=.3) + # Vertical dashed line
geom_hline(yintercept = diff_line, linetype = "dashed", color = "grey44", linewidth=.3) + # Vertical dashed line
gtheme_no_rot + # Use a minimal theme
 geom_text_repel(
  segment.size=.05,
  max.overlaps = 100,
# point.padding = 0.2, 
  size=1.8#,
# nudge_x = .25,
# nudge_y = .2,
# segment.curvature = -1e-20
  ) +
xlab ('activity difference') + 
ylab ('RNA difference') + 
xlim (c(-0.2, .2)) + 
ylim (c(-0.6, .6))

pdf (paste0 ('Plots/Diff_normal_tumor_deviation_and_rna_scatterplot3.pdf'),5,height = 4)
print (tf_diff_p)
dev.off()
saveRDS(selected_TF, 'selected_TF.rds')
} else {
selected_TF = readRDS ('selected_TF.rds')
}







# ### Try running wilcox test for each normal and tumor sample comb ####
#   # Lists of TFs are quite different using ENCODE normal samples (they have very low # cells) so consider discard them
# metaGroupName = 'Sample3'
# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# if (!exists('gsSE')) gsSE = fetch_mat (archp, 'GeneScore')
# archp_meta = as.data.frame (archp@cellColData)

# mMat = scale(assays (mSE)[[1]])
# rownames (mMat) = rowData (mSE)$name
# # mMat_agg = as.data.frame (t(mMat))
# # mMat_agg$metaGroup = as.character (archp_meta[,metaGroupName])
# # mMat_agg = aggregate (.~ metaGroup, mMat_agg, mean)
# # rownames (mMat_agg) = mMat_agg[,1]
# # mMat_agg = as.data.frame(t(mMat_agg[,-1]))

# gsMat = assays (gsSE)[[1]]
# rownames (gsMat) = rowData (gsSE)$name
# gsMat = scale(gsMat[rownames(gsMat) %in% rownames(mMat_agg),])
# # gsMat_agg = as.data.frame (t(gsMat))
# # gsMat_agg$metaGroup = as.character (archp_meta[,metaGroupName])
# # gsMat_agg = aggregate (.~ metaGroup, gsMat_agg, mean)
# # rownames (gsMat_agg) = gsMat_agg[,1]
# # gsMat_agg = as.data.frame(t(gsMat_agg[,-1]))

# # Set up comparisons ####
# normal_sams = unique (archp_meta[,metaGroupName][grep ('normal',archp_meta[,metaGroupName])])
# tumor_sams = unique (archp_meta[,metaGroupName][grep ('P',archp_meta[,metaGroupName])])
# outlier_sams = 'P11_HOX'
# tumor_sams = tumor_sams[!tumor_sams %in% outlier_sams]

# pval_threshold = 0.05
# occurrence_threshold_normal = length(normal_sams) - 1
# occurrence_threshold_tumor = 4

# dev_res_l=list()
# for (tumor_sam in tumor_sams)
#   {
#   dev_res = lapply (normal_sams, function(normal_sam)
#   {
#   res = wilcoxauc (mMat, y = as.character (archp_meta[,metaGroupName]), groups_use = c(tumor_sam, normal_sam))
#   res = res[res$group == tumor_sam,]
#   res[res$logFC > 0 & res$padj < pval_threshold,'feature']
#   })
#   sig_features = table (unlist(dev_res)[table (unlist(dev_res)) >= occurrence_threshold_normal]) 
#   dev_res_l[[tumor_sam]] = names(sig_features)
#   }
# sig_features_dev_final = names(table(unlist(unname(dev_res_l)))[table(unlist(unname(dev_res_l))) > occurrence_threshold_tumor])
# sig_features_dev_final

# gs_res_l=list()
# for (tumor_sam in tumor_sams)
#   {
#   gs_res = lapply (normal_sams, function(normal_sam)
#   {
#   res = wilcoxauc (gsMat, y = as.character (archp_meta[,metaGroupName]), groups_use = c(tumor_sam, normal_sam))
#   res[res$logFC > 0 & res$padj < pval_threshold,'feature']
#   })
#   sig_features = table (unlist(gs_res)[table (unlist(gs_res)) >= occurrence_threshold_normal]) 
#   gs_res_l[[tumor_sam]] = names(sig_features)
#   }
# sig_features_gs_final = names(table(unlist(unname(gs_res_l)))[table(unlist(unname(gs_res_l))) > occurrence_threshold_tumor])
# intersect (sig_features_dev_final, sig_features_gs_final)


# head (mMat_agg['RUNX1',])
# rownames(mMat_agg)[p.adjust(sapply (seq(nrow(mMat_agg)), function(x) t.test (mMat_agg[x,grep ('normal', colnames(mMat_agg))], mMat_agg[x,!grepl ('normal', colnames(mMat_agg))])$p.value)) < 0.05]

# # dev_res = lapply (dev_comparisons, function(x) 
# #     wilcoxauc (mMat, y = archp$Sample2, groups_use = x))


# Make coexpression network for each sample using top TFs deviations ####
selected_TF = readRDS ('selected_TF.rds')
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
archp_meta = as.data.frame (archp@cellColData)
sams = as.character(unique(archp_meta$Sample3))
sams = sams[!sams %in% c('normal1','P3','P13','P11_HOX')] # remove normal,low cell numbers and outlier samples

#archp_meta = archp_meta[archp_meta$Sample3 %in% sams,]
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = t(as.matrix(scale(mMat[selected_TF,])))

all (rownames(mMat) == rownames(archp_meta))
cor_TF_l = list()
for (sam in sams)
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
  row_split = clusters,
  column_split = clusters,
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

pdf (file.path ('Plots','selected_TF_dev_corr_heatmaps3.pdf'), width = 8,height=7)
cor_TF_df 
dev.off()

# Add pathways of TF correlated genes from scRNA ####
if (!file.exists('enrichment_pathways_TFs.rds')) source (file.path('..','..','git_repo','tumor_analysis','enrichment_cnmfs.R'))
TF_cor_sum = readRDS (file.path('enrichment_pathways_TFs.rds'))
TFrow_order = unname(unlist(row_order (cor_TF_df)))
TFrow_order_split = rep (names(row_order (cor_TF_df)), lapply (row_order (cor_TF_df),length))
#TFrow_order_split = TFrow_order_split[TFrow_order]
rownames (TF_cor_sum) = gsub ('HALLMARK_','',rownames(TF_cor_sum))
rownames (TF_cor_sum) = gsub ('_', ' ',rownames (TF_cor_sum))
TF_cor_sum = TF_cor_sum[apply (TF_cor_sum, 1, function(x) any(x > 1)),]
pdf()
hm = draw (Heatmap (
    TF_cor_sum[,TFrow_order],
    column_split = ifelse(TFrow_order_split =='1',2,1),
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

pdf (file.path ('Plots','selected_TF_dev_corr_pathways_heatmaps2.pdf'), width = 9.5,height=2.6)
hm
dev.off()


#### Generate heatmap of significant Cox proportion hazards from bulk-RNA ####
coxht = readRDS (file.path('..','..','bulkRNA_meso','oncoTF_sig_cox_bulkRNA.rds'))

ht = Heatmap (coxht[[1]], 
    clustering_distance_columns = 'euclidean',
    clustering_distance_rows = 'euclidean',
    #column_title = paste('* = < ',pvalAdjTrheshold,'top',topGenes,'genes'),
    col = plotcol,
    #heatmap_legend_param = list(
    #     at = c(col_limit, 0, -col_limit),#,
    #    direction = "horizontal",#,
    #    #legend_height = unit(5, "cm"),
    #    just = c('top')
    #    #title_position = "leftcenter-rot"
   #    ),
    bottom_annotation = labels_annotation,
    column_names_gp = gpar(fontsize = 0),
    row_names_gp = gpar(fontsize = 5),
    cell_fun = function (j, i, x, y, width, height, fill) 
            {
           if (t(pval_mat)[i, j] < pvalAdjTrheshold3)
              {
               grid.text("***", x, y, just='center', vjust=.8,
                gp = gpar(fontsize = 5, col='black'))
              } else {
              if(t(pval_mat)[i, j] < pvalAdjTrheshold2)
                  {
                  grid.text("**", x, y, just='center', vjust=.8,
                  gp = gpar(fontsize = 5, col='black'))   
                  } else {
                  if(t(pval_mat)[i, j] < pvalAdjTrheshold)
                    {
                    grid.text("*", x, y, just='center', vjust=.8,
                    gp = gpar(fontsize = 5, col='black'))         
                    }}}
      }, ...)




# ### Show sarcomatoid score per sample regenerating UMAPs ####
sarc_module = 'cNMF20'
p_l2=list()
for (sam in sams)
  {
  archp_sub = archp[archp$Sample3 == sam]
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

pdf (file.path('Plots','sarcomatoid_score_feature_plots3.pdf'), width = 12, height = 15)
print (wrap_plots (p_l2), ncol = 5)
dev.off()


# Plot boxplots ordered by correlation to sarcomatoid score using TF activity and scRNA ####
  
  selected_TF = readRDS ('selected_TF.rds')
  sarc_module = 'cNMF20'
  archp_meta = as.data.frame (archp@cellColData)
  sams = as.character(unique(archp_meta$Sample3))
  sams = sams[!sams %in% c('normal1','normal2','normal3','normal4','P3','P13','P11_HOX')] # Remove normal lung, low cell numbers and outlier samples

  # # Take only TF expressed in RNA tumors ####
  # ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat_agg), group.by = 'sampleID3')[[1]]) +1)
  # colnames(ps) = gsub ('-','_',colnames(ps))
  # ps = ps[,!colnames(ps) %in% c('HU37','HU62')]
  # min_exp = .5
  # active_TFs = rownames(ps)[apply (ps, 1, function(x) any (x > min_exp))]

  # Get deviations ####
  if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
  mMat = assays (mSE)[[1]]
  rownames (mMat) = rowData (mSE)$name
  mMat = mMat[selected_TF,]

  # Get genescore ####
  if (!exists('gsSE')) gsSE = fetch_mat (archp, 'GeneScore')
  gsMat = assays (gsSE)[[1]]
  rownames (gsMat) = rowData (gsSE)$name
  gsMat = gsMat[selected_TF,]
  
  mMat = scale(t(mMat))
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

  # Combine and plot ####
  combined_df = rbind (m_cor_df, gs_cor_df)
  
  # ranked_data <- combined_df %>%
  # group_by(TF, type) %>%
  # summarise(median_score = median(score), .groups = "drop") %>%
  # group_by(TF) %>%
  # summarise(mean_of_medians = mean(median_score), .groups = "drop") %>%
  # arrange(desc(mean_of_medians)) %>%
  # mutate(rank = row_number())

  #ps = ps[, colnames(ps) %in% sample_names_rna]

  #lapply (tc_cor, function(x) {x = x['cNMF19',]; head(x[order(-x)],10)})

top_sarc_TF = head(m_cor_levels$TF,20)
combined_df$TF = factor (combined_df$TF, levels = top_sarc_TF)
combined_df = combined_df[combined_df$TF %in% top_sarc_TF, ]
palette_expression_disc = paletteer::paletteer_c("grDevices::Purple-Blue", length(top_sarc_TF))

bp = ggplot(combined_df, aes(x = TF, y = score, fill = TF, color = type)) + 
  geom_boxplot(
    aes(group = interaction(TF, type)),  # Split by both TF and type
    position = position_dodge(0.8),
    #color = 'grey20',
    linewidth = 0.4,
    width = 0.7,
    outlier.alpha = 0.2,
    outlier.shape = NA,
    size = 0.5,
    alpha = 0.6
  ) +
  geom_point(
    aes(group = interaction(TF, type), color = type),  # Align points with boxplots
    position = position_dodge(0.8),  # Ensure same dodge width
    alpha = 0.5,
    size = 1
  ) +
  gtheme +
  scale_fill_manual(values = palette_expression_disc) +
  scale_color_manual(values = c(activity = '#AE123AFF', genescore = '#001260FF')) +
  geom_hline(yintercept = 0, color = 'red', linetype = 'dashed')

  pdf (paste0 ('Plots/sarcomatoid_score_TF_activity_boxplots.pdf'), width = 8,height=3)
  bp
  dev.off()



# Compare with TF correlation to sarcomatoid in scRNA ####
metacells = readRDS (file.path('..','scrna','metacells.rds'))
metacells$sampleID = metacells$sampleID3
nfeat=5000
k=25
cnmf_spectra_unique = readRDS (paste0('../scrna/cnmf_genelist_',k,'_nfeat_',nfeat,'.rds'))

sams_rna = c('P1','P4','P5','P8','P11','P12')

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
metacells_assay = metacells_assay[selected_TF,]

rna_cor = lapply (sams_rna, function(sam)
  {
  tmp = cor(as.matrix(t(metacells_assay[,metacells$sampleID == sam])), as.matrix(metacells@meta.data[,sarc_module])[metacells$sampleID == sam,,drop=F], method='spearman')
  tmp = as.data.frame (tmp)
  tmp$sample = sam
  tmp$TF = rownames(tmp)
  colnames (tmp) = c('score','sample','TF')
  tmp
  })

rna_tf_cor_df = do.call (rbind, rna_cor)
rna_tf_cor_df$type = 'expression'


# Combine and plot ####
combined_df_rna = rbind (m_cor_df, rna_tf_cor_df)
combined_df_rna$TF = factor (combined_df_rna$TF, levels = m_cor_levels$TF)
#lapply (tc_cor, function(x) {x = x['cNMF19',]; head(x[order(-x)],10)})

top_sarc_TF = head(m_cor_levels$TF,20)
combined_df_rna = combined_df_rna[combined_df_rna$TF %in% top_sarc_TF, ]
palette_expression_disc = paletteer::paletteer_c("grDevices::Purple-Blue", length(top_sarc_TF))

bp = ggplot(combined_df_rna, aes(x = TF, y = score, fill = TF, color = type)) + 
  geom_boxplot(
    aes(group = interaction(TF, type)),  # Split by both TF and type
    position = position_dodge(0.8),
    #color = 'grey20',
    linewidth = 0.4,
    width = 0.7,
    outlier.alpha = 0.2,
    outlier.shape = NA,
    size = 0.5,
    alpha = 0.6
  ) +
  geom_point(
    aes(group = interaction(TF, type), color = type),  # Align points with boxplots
    position = position_dodge(0.8),  # Ensure same dodge width
    alpha = 0.5,
    size = 1
  ) +
  gtheme +
  scale_fill_manual(values = palette_expression_disc) +
  scale_color_manual (values = c(activity = '#AE123AFF', expression = 'darkgreen')) +
  geom_hline(yintercept = 0, color = 'red', linetype = 'dashed')

  pdf (paste0 ('Plots/sarcomatoid_score_TF_expression_boxplots.pdf'), width = 8,height=3)
  bp
  dev.off()


# Compute co-occurrence of sarcomatoid TFs ####
motifMat = getPositions (archp)
matches = getMatches (archp)
matchesMat = assay (matches)
colnames (matchesMat) = gsub ('_.*','',colnames (matchesMat))
colnames (matchesMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames (matchesMat))

matchesMat = matchesMat[,top_sarc_TF]
matchesMat = matchesMat[rowSums (matchesMat) > 0,]

cooc = matrix (ncol = ncol(matchesMat), nrow= ncol(matchesMat))

for (i in 1:ncol(matchesMat))
  {
  for (z in 1:ncol(matchesMat)) 
    {
    ov = sum (rowSums (matchesMat[,c(i,z)]) == 2) / min (colSums(matchesMat[,c(i,z)]))
    cooc[i,z] = ov
    }
  }

colnames (cooc) = top_sarc_TF
rownames (cooc) = top_sarc_TF
diag (cooc) = 1

cooc_hm = Heatmap (
    cooc,
    column_names_rot =45, 
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 6),
    col =palette_cooccurrence, 
    cluster_rows=F,
    cluster_columns = F,
rect_gp = gpar (col = "white", lwd = 1))

pdf (file.path ('Plots','selected_TF_cooccurence_heatmaps2.pdf'), width = 6,height=5)
cooc_hm 
dev.off()


### Import Blum meta-analysis to compare with top TF correlated with scS-score ####
blum_df = read.csv (file.path('..','Blum_et_al_SE_score.csv'))
blum_dfE = data.frame (gene = blum_df$Gene.Name, score = blum_df$correlation.E.score, SE_score = 'epithelioid')
blum_dfS = data.frame (gene = blum_df$Gene.Name.1, score = blum_df$correlation.S.score, SE_score = 'sarcomatoid')
blum_df = rbind (blum_dfE, blum_dfS)
blum_df = na.omit (blum_df)
rownames (blum_df) = blum_df$gene
blum_df = blum_df[top_sarc_TF, ]
blum_df$gene = rownames (blum_df)
blum_df$gene = factor (blum_df$gene, levels = blum_df$gene)
# Create the dot plot
dp = ggplot(blum_df, aes(x = gene, y = 1)) +
  geom_point(aes (size = score, color = SE_score)) + # Adds the points
  labs(title = "Correlation to Blum et al") + # Labels
  scale_color_manual (values = c(epithelioid = 'darkgreen',sarcomatoid = 'firebrick2')) + gtheme

pdf (file.path ('Plots','Blum_top_sarc_TF2.pdf'))
dp
dev.off()



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

dev.off()

pdf (file.path('Plots','sarcomatoid_score_feature_plots2.pdf'), width = 20, height = 20)
print (wrap_plots (c(list(p), TF_p), ncol = 3))
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

# Prepare files for chromBPnet analysis ####
source (file.path ('..','..','git_repo','tumor_analysis','chromBPnet_prepare_files_sarc_score.R'))



### Combine TF modisco and finemo outputs to build network of co-occurring TFs across peaks ####
library (httr)
library (XML)
library (igraph)

chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet'
fold_number = 0
celltypes = c('P23__epithelioid','P23__sarcomatoid') # in new clustering is C1


count_to_adj = function(data = NULL)
  {
  count1s <- function(x, y) colSums(x == 1 & y == 1)
  n <- 1:ncol(data)
  mat <- outer(n, n, function(x, y) count1s(data[, x], data[, y]))
  diag(mat) <- 0
  dimnames(mat) <- list(colnames(data), colnames(data))
  return (mat)
  }
    
ov_motif_peaks_adj_l = list()
for (celltype in celltypes)
  {
  #celltype= 'IL1B'
  modisco_motifs = as.data.frame(readHTMLTable(file.path(chromBPdir, paste0(celltype,'_model'),'fold_0','modisco','report','motifs.html')))
  modisco_motifs$motif_match0 = sapply (modisco_motifs$NULL.match0, function(x) unlist(strsplit (x, '_'))[1])
  modisco_motifs$motif_match1 = sapply (modisco_motifs$NULL.match1, function(x) unlist(strsplit (x, '_'))[1])
  modisco_motifs$motif_match2 = sapply (modisco_motifs$NULL.match2, function(x) unlist(strsplit (x, '_'))[1])
  
  finemo_hits = read.table(file.path(chromBPdir,paste0(celltype,'_model'),'finemo_out','hits.tsv'), sep='\t', header=T)
  finemo_hits$motif_name0 = modisco_motifs$motif_match0[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
  finemo_hits$motif_name1 = modisco_motifs$motif_match1[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
  finemo_hits$motif_name2 = modisco_motifs$motif_match2[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
  write.table (finemo_hits[,c(1,2,3,14,15,16)], paste0(celltype,'_finemo_to_genome_browser.tsv'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  peakset = read.table (file.path(chromBPdir,paste0('peakset_',celltype,'.bed')))
  peakset = peakset[,1:3]
  colnames (peakset) = c ('chr','start','end')
  peakset = makeGRangesFromDataFrame (peakset)
  finemo_hits = makeGRangesFromDataFrame (finemo_hits, keep.extra.columns = T)
  
  match_rank = c('motif_name0')#,'motif_name1','motif_name2')
  ov_motif_peaks_mat_combined = list()
  for (motif_match in match_rank)
    {
    finemo_hits_l = split (finemo_hits, finemo_hits@elementMetadata[,motif_match])
    
    ov_motif_peaks = lapply (finemo_hits_l, function(x) findOverlaps (peakset, x, select='first'))
    ov_motif_peaks_mat = do.call (cbind, ov_motif_peaks)
    rownames (ov_motif_peaks_mat) = as.character(peakset)
    ov_motif_peaks_mat[is.na(ov_motif_peaks_mat)] = 0
    ov_motif_peaks_mat[ov_motif_peaks_mat > 0] = 1
    ov_motif_peaks_mat_combined[[motif_match]] = ov_motif_peaks_mat
    #ov_motif_peaks_df = as.data.frame (ov_motif_peaks_mat)
    }
  ov_motif_peaks_mat_combined = do.call (cbind, ov_motif_peaks_mat_combined)  
  tf_columns = colnames(ov_motif_peaks_mat_combined)
  ov_motif_peaks_adj_l2 = as.data.frame (count_to_adj (data = ov_motif_peaks_mat_combined))
  # From https://stackoverflow.com/questions/66515117/convert-dummy-coded-matrix-to-adjacency-matrix
  ov_motif_peaks_adj_l[[celltype]] = ov_motif_peaks_adj_l2
  }

# Try with heatmaps 

# Try with heatmaps 
palette_cooccurence2 = paletteer::paletteer_c("grDevices::Oslo",100)
all_TF = unique(unname(unlist(lapply (ov_motif_peaks_adj_l, function(x) unique(colnames(x))))))
all_TF = all_TF[all_TF != 'NaN']
ht_list = NULL 
set.seed (1234)
for (celltype in celltypes)
{
  hm_mat = as.data.frame (ov_motif_peaks_adj_l[[celltype]][match(all_TF, colnames(ov_motif_peaks_adj_l[[celltype]])),])
  hm_mat = as.data.frame (t(hm_mat[all_TF,]))[all_TF,]
  rownames (hm_mat) = all_TF
  colnames (hm_mat) = all_TF
  hm_mat[is.na(hm_mat)] = 0
  #hm_mat = scale (hm_mat)
  hm_mat = log2(hm_mat+1)
  if (is.null(ht_list)) {
  d = as.dist(t(1-cor(hm_mat)))
  d[is.na(d)] = 1
  hc = hclust(d)
  column_order = colnames(hm_mat)[hc$order]
  d = as.dist(1-cor(hm_mat))
  d[is.na(d)] = 1
  hc = hclust(as.dist(d))
  row_order = colnames(hm_mat)[hc$order]
  ht_list = Heatmap (
    hm_mat, 
    column_names_gp= gpar(fontsize=0),
    row_names_gp= gpar(fontsize=3), column_order=column_order, row_order=row_order,
    col = palette_cooccurence2, border=T, column_title = celltype)
  } else {
  ht_list = ht_list + Heatmap (
  hm_mat[,column_order], cluster_columns=F,
  column_names_gp= gpar(fontsize=0),
  row_names_gp= gpar(fontsize=3),
  col = palette_cooccurence2, border=T, column_title = celltype)
  }
}
pdf (file.path ('Plots','combined_chromBPnet_cooccurrence_heatmaps.pdf'),width=10,height=4)
ht_list
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


# Check SOX9 footprinting ####
archp$sarc_score_sample = paste0(archp$sarc_score, archp$Sample)
# archp2 = archp
# archp = archp[archp$sarc_score != 'mid']
archp <- addGroupCoverages (ArchRProj = archp, groupBy = "Clusters")
motifPositions <- getPositions (archp)

motifs <- c("SOX9", "SOX6",'RUNX2','TCF3','SNAI2','SNAI1','TWIST2','TWIST1','HIC2','HIC1')
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
metaGroupName='Clusters'
#for (sam in sams)
#  {
  #peaks_sample = readRDS (file.path ('PeakCalls',paste0(sam,'-reproduciblePeaks.gr.rds')))
  #motifPositions_sample = lapply (motifPositions, function(x) x[queryHits(findOverlaps(x, peaks_sample))])  
  seFoot <- getFootprints(
    ArchRProj = archp, 
    #positions = motifPositions_sample[markerMotifs], 
    positions = motifPositions[markerMotifs], 
    groupBy = metaGroupName
  )
  
    plotFootprints(
    seFoot = seFoot[,grepl('C10', colnames(seFoot)) | grepl('C7', colnames(seFoot))],
    ArchRProj = archp, 
    normMethod = "Subtract",
    plotName = paste0("Footprints-Subtract-Bias_",sam),
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
