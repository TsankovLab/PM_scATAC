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
addArchRThreads(threads = 1) 
addArchRGenome("hg38")

if (!file.exists ('Save-ArchR-Project.rds')) 
  { source (file.path('..','..','PM_scATAC','scatac_tumor_create_ArchRobj.R'))
  } else {
 archp = loadArchRProject (projdir)   
  }

# Load RNA ####
srt = readRDS (file.path('..','scrna','srt.rds'))
sarc_order = read.csv (file.path('..','scrna','cnmf20_sarcomatoid_sample_order.csv'), row.names=1)
sarc_order = sarc_order[! sarc_order$sampleID %in% c('HU37','HU62'),]
sarc_order = rbind (data.frame (sampleID = 'normal_pleura', x = -1),sarc_order)
#archp$Sample2 = factor (archp$Sample2, levels = sarc_order$sampleID)

# Run genescore DAG ####
metaGroupName = "Clusters"
force=FALSE
if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source ('../../PM_scATAC/DAG.R')

celltype_markers = c('WT1','CALB2','RUNX2','TCF3','SOX9','MESP1','HMGA1','TWIST1','SNAI2')
archp = addImputeWeights (archp)
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "GeneScoreMatrix", 
    name = celltype_markers, 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)

pdf (file.path('Plots','marker_genes_feature_plots.pdf'), width = 20, height = 20)
print (wrap_plots (p, ncol = 4))
dev.off()

### Plot cell type markers on genome tracks ####
metaGroupName = 'Sample2'
celltype_markers = c('WT1','CALB2','RUNX2','TCF3','SOX9','MESP1','HMGA1','TWIST1','SNAI2')
#celltype_markers = c('WT1','CALB2','GATA4','MSLN','KRT5','KRT18','ITLN1','HP','SOX9')
meso_markers <- plotBrowserTrack(
    ArchRProj = archp, 
    groupBy = metaGroupName, 
    geneSymbol = celltype_markers,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 250000,
    downstream = 250000,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
    #loops = getCoAccessibility (archp, corCutOff = 0.3,
    #  returnLoops = TRUE),
    useGroups= NULL
)
plotPDF (meso_markers, ArchRProj = archp, width=14, name ='MPM_markers_coveragePlots.pdf')
}


### Compare bins malignants normal ####
run_bin_analysis = FALSE

if (run_bin_analysis)
  {
  # Load fragments
  fragments = unlist (getFragmentsFromProject (
       ArchRProj = archp))
  
  ws = 1e6
  ss = 2e5
  if (!file.exists (paste0('bins_',ws,'_ss_',ss,'.rds')))
    {
    blacklist = toGRanges(paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/blacklisted_regions/ENCODE_blacklist/',"hg38-blacklist.v2.bed")) # taken from https://github.com/Boyle-Lab/Blacklist/tree/master/lists
    windows = makeWindows(genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist,
      windowSize = ws, slidingSize = ss)
    saveRDS (windows, paste0('bins_',ws,'_ss_',ss,'.rds'))
    } else {
    windows = readRDS (paste0('bins_',ws,'_ss_',ss,'.rds'))  
    }
  
  flt_celltypes = 100
  
  #archp$celltype_sample = paste0(archp$celltype, '-', archp$Sample2)
  metaGroupName = 'Clusters'
  
  #fragments$RG = as.character(fragments$RG)  
  barcode_metaGroup = as.data.frame (archp@cellColData[,c(metaGroupName, 'TSSEnrichment','nFrags')])
  colnames (barcode_metaGroup)[colnames (barcode_metaGroup) == metaGroupName] = 'metaGroup'
  barcode_metaGroup$barcode = rownames(barcode_metaGroup)
  celltype_bins = lapply (unique(archp@cellColData[,metaGroupName]), function(mg) 
    {
    high_quality_cells = barcode_metaGroup[barcode_metaGroup$metaGroup == mg,]
    high_quality_cells = high_quality_cells[order (-high_quality_cells$TSSEnrichment),] # ordering by TSSEnrichment seem to work the best
    high_quality_cells = head (high_quality_cells$barcode, flt_celltypes)
    fragments_group = fragments[fragments$RG %in% high_quality_cells]
    #names (fragments_cell) = NULL 
    #fragments_cts = sample (unique(fragments_ct$RG), 200)
    fr_ov = countOverlaps (windows, fragments_group)
  #  fr_ov = fr_ov / sum (archp$nFrags[archp$celltype_sample == x])
    #(fr_ov / sum (fr_ov)) * 1000
    })
  celltype_bins = do.call (cbind, celltype_bins)
  colnames (celltype_bins) = sapply (unique(archp@cellColData[,metaGroupName]), function(x) unlist(strsplit (x,'-'))[1])
  head (celltype_bins)
  
  celltype_bins_cor = cor (celltype_bins, method = 'spearman')
  ha = HeatmapAnnotation (sample = sapply (unique(archp@cellColData[,metaGroupName]), function(x) unlist(strsplit (x,'-'))[2]), which='row')
  binH = Heatmap (celltype_bins_cor, col = viridis::rocket(100),# row_names_gp= gpar (fontsize=6), 
    #column_names_gp= gpar (fontsize=6), 
    right_annotation = ha,
    cluster_rows = T, #row_km = 2, column_km = 2,
    clustering_distance_rows = 'pearson', 
    clustering_distance_columns = 'pearson',
    column_title = paste('Binned',ws))
  
  png (paste0('Plots/binned_fragments_binned_',ws,'_celltype_heatmap.png'),width=600, height=500)
  binH
  dev.off()
  }


### Run peak calling ####
metaGroupName = "Clusters"
force=FALSE
if(!all(file.exists(file.path('PeakCalls', paste0(unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds')))) | force) source ('../../PM_scATAC/callPeaks.R')
  



### chromVAR analysis ####
force=FALSE
if (!all(file.exists(file.path('Annotations',
  c('Motif-Matches-In-Peaks.rds',
    'Motif-Positions-In-Peaks.rds',
    'Motif-In-Peaks-Summary.rds')))))
source (file.path('..','..','git_repo','utils','chromVAR.R'))



### Co-expression of TFs #### 
metaGroupName = 'Sample2'
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
all (colnames(mSE) == rownames(archp))

# # Get deviation matrix ####
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Filter by RNA expression ####
metaGroupName = 'sampleID'
active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)
mMat = mMat[active_TFs, ]


mMat_cor = cor (as.matrix(t(scale(mMat))), method = 'spearman')

km = kmeans (mMat_cor, centers=5)

pdf (file.path ('Plots','TF_modules_heatmap.pdf'), width = 4,height=3)
cor_mMat_hm = draw (Heatmap (mMat_cor,# row_km=15,
  #left_annotation = ha,
  #rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=palette_deviation_cor_fun, 
  row_split = km$cluster,
  column_split = km$cluster,
  #row_km=2, 
  #column_km=2,
#  right_annotation = ha,
  border=T,
#   ,
  row_names_gp = gpar(fontsize = 0), 
  column_names_gp = gpar(fontsize = 0)
# cell_fun = function(j, i, x, y, w, h, fill) {# THIS DOESNT WORK NEED TO USE LAYER_FUN
#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#         }}
  ))
  # ,
  # cell_fun = function(j, i, x, y, w, h, fill) {
  #       if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
  #           grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#        }}))
dev.off()

pdf (file.path ('Plots','TF_modules_heatmap.pdf'), width = 4,height=3)
cor_mMat_hm
dev.off()

tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = do.call (cbind, tf_modules)
archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% paste0('mod_',unique(km$cluster))]
archp@cellColData = cbind (archp@cellColData, tf_modules) 

pdf()
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "cellColData",
    name = paste0('mod_',unique(km$cluster)), 
    pal = rev(palette_deviation),
    #useSeqnames='z',
    embedding = "UMAP")
dev.off()
pdf (file.path ('Plots','TF_modules_umap.pdf'), width = 20,height=6)
wrap_plots (TF_p, ncol=5)
dev.off()


# ridge plots of TF modules ####
library(ggridges)
library(ggplot2)
library(viridis)
#library(hrbrthemes)
tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = as.data.frame (do.call (cbind, tf_modules))
tf_modules$Sample = archp$Sample2
# sapply (rownames(tf_modules), function(x) unlist(strsplit (x, '\\#'))[1]) == archp$Sample2
tf_modules = gather (tf_modules, module, expression,1:5)
tf_modules$module = factor (tf_modules$module, levels = paste0('mod_',names (row_order (cor_mMat_hm))))

# Plot
# rp = lapply (paste0('mod_',unique(km$cluster)), function(x) 
#   ggplot(tf_modules,  aes_string(x='Sample', y=x, fill='..x..')) +

#   geom_vridgeline(stat="ydensity", trim=FALSE, alpha = 0.85, scale = 2, width=.10) +
#   palette_deviation_ggplot_fill +
#     theme_classic())
# pdf (file.path ('Plots','TF_modules_ridge_plots.pdf'), width = 20,height=3)
# wrap_plots (rp, ncol=5)
# dev.off()


dp = ggplot (tf_modules) +
  geom_density(aes(x=expression,fill=Sample),
                      alpha = 0.4) +
  # geom_vline(aes(xintercept = mean, group = tf_modules, linetype = Sample),
  #            data = combined_sla_means) +
  facet_wrap (~module, nrow = 5, scales = 'free',strip.position = "left") +
  scale_fill_manual (values = palette_sample) +
  gtheme_no_rot

pdf (file.path ('Plots','TF_modules_ridge_plots2.pdf'), width = 7,height=8)
dp
dev.off()

# violin plots of TF modules ####
# tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
# names (tf_modules) = paste0('mod_',unique(km$cluster))
# tf_modules = as.data.frame (do.call (cbind, tf_modules))
# tf_modules$Sample = archp$Sample2
# tf_modules = gather (tf_modules, module, expression,1:5)

# box = ggplot (tf_modules, aes_string (x= 'Sample', y= 'expression')) +
#   geom_violin (trim=TRUE, aes_string (fill = 'Sample'),size=2,
#     width=1,
#     scale='width',
#     linewidth = .2, alpha=0.7) +
#   geom_boxplot (aes_string(fill = 'expression'),
#     linewidth = .2,
#     width=0.2,
#     outlier.alpha = 0.2,
#     outlier.size = 1,
#      size=0.3, alpha=0.7
#      ) +
#   gtheme + 
#   facet_wrap (~module, nrow = 5) + 
# #  palette_deviation_ggplot_fill + 
#   NoLegend()
  
# pdf(file.path('Plots','TF_modules_boxplots.pdf'),width=3,height = 8) #width = 10, height = 11,
# print (box)
# dev.off()


# rp = lapply (paste0('mod_',unique(km$cluster)), function(x) 
#   ggplot(tf_modules, aes_string(x = x, y = 'Sample', fill = '..x..')) +
#   geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, alpha=.5) +
#   palette_deviation_ggplot_fill +
#     theme_classic())
# pdf (file.path ('Plots','TF_modules_ridge_plots.pdf'), width = 20,height=3)
# wrap_plots (rp, ncol=5)
# dev.off()


# Check expression of SOX9 and others - SOX9 should cluster with module 5 but it doesnt probably cause of low accuracy of deviation ####
tf_name2 = 'SOX9_756'

tf_name2 = unlist(sapply (c('SOX9','TWIST1','MESP1','SNAI2','TEAD1'), function(x) rownames(assay(mSE))[grepl (x, rownames(assay(mSE)))]))
tf_name2 = paste0('z:',tf_name2)
archp = addImputeWeights (archp)
pdf ()
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "MotifMatrix",
    name = tf_name2, 
    useSeqnames='z',
    pal = rev(as.character(palette_deviation)),    
    embedding = "UMAP",
    imputeWeights = getImputeWeights(archp)
    )
dev.off()
pdf (file.path ('Plots','TF_umap.pdf'), width = 20,height=6)
wrap_plots (TF_p, ncol=5)
dev.off()


# ### Use P2G analysis and cNMF from RNA to identify active TF via regulons  ####
# run_p2g_TF = FALSE

# if (run_p2g_TF)
#   {
#   run_p2g = FALSE  
#   if (run_p2g)
#     {
#     maxDist = 250000
#     archp = addPeak2GeneLinks(
#         ArchRProj = archp,
#         useMatrix = 'GeneScoreMatrix',
#         reducedDims = "IterativeLSI",
#         maxDist = maxDist
#     )
#     }  
    
#   p2g_corr = .2
#   p2g = getPeak2GeneLinks(
#       ArchRProj = archp,
#       corCutOff = p2g_corr,
#       resolution = 1,
#       returnLoops = FALSE
#   )
  
#   # # Convert df in Granges add gene Name and correlation ####
#   # gene = 'WT1'
#   p2g$geneName = mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
#   p2g = p2g[!is.na (p2g$FDR),] # remove NaN correlations (not sure why there are some)
#   p2g = p2g[p2g$Correlation > p2g_corr, ]
#   p2gGR = metadata (p2g)$peakSet[p2g$idxATAC]
#   p2gGR$geneName = p2g$geneName
#   p2gGR$correlation = p2g$Correlation
  


#   cnmf_list = lapply (cnmf_list, function(x) head (x,200))
#   p2g_cnmf = lapply (cnmf_list, function(x) p2gGR[p2gGR$geneName %in% x])
  
#   tf_match = getMatches (archp)
#   colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
#   bg_peakSet = rowRanges (tf_match)
#   #tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
#   nmf_TF = lapply (p2g_cnmf, function(x) hyperMotif (
#     selected_peaks = x, 
#     motifmatch = tf_match))
  
#   nmf_TF = lapply (nmf_TF, function(x) x[rownames(nmf_TF[[1]]), ])
#   nmf_TF_df = do.call (cbind, nmf_TF)
#   nmf_TF_df = nmf_TF_df[, grep ('padj', colnames(nmf_TF_df))]
#   nmf_TF_df[nmf_TF_df > 0.05] = 1
#   nmf_TF_df = -log10(nmf_TF_df)
#   nmf_TF_df = nmf_TF_df[rowSums (nmf_TF_df) != 0, ]
#   nmf_TF_df[sapply(nmf_TF_df, is.infinite)] <- 300
#   saveRDS (nmf_TF_df, 'nmf_TF_p2g_enrichments.rds')
#   TF_ht = Heatmap (nmf_TF_df, row_names_gp = gpar (fontsize=3), column_names_gp = gpar (fontsize=5))
  
#   pdf (paste0('Plots/TF_nmf_',k,'_nfeat_',nfeat,'_heatmap.pdf'),width = 3,height=25)
#   print (TF_ht)
#   dev.off()
#   } else {
#   nmf_TF_df = readRDS ('nmf_TF_p2g_enrichments.rds')
#   }

  # Make data.frame of deviation difference and expression difference between normal and tumors ####
if(!file.exists('selected_TF.rds'))
  {
  metaGroupName = 'Sample2'
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
  
  metaGroupName = 'sampleID'
  DefaultAssay (srt) = 'RNA'
  sample_names_rna = c('P1','P14','P13','P3','P12','P5','P11','P4','P8','P14','HU37','HU62')
  ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat_agg), group.by = metaGroupName)[[1]]) +1)
  ps = ps[, colnames(ps) %in% sample_names_rna]
  
  TF_diff_rna = data.frame (
    tumor_dev = apply (mMat_agg[,unique(archp$Sample2)[!unique(archp$Sample2) == 'normal_pleura']], 1, mean),
    normal_dev = mMat_agg[,unique(archp$Sample2)[unique(archp$Sample2) == 'normal_pleura']],
    normal_rna = apply (ps[,c('HU37','HU62')], 1, mean),
    tumor_rna = apply (ps[, !colnames(ps) %in% c('HU37','HU62')], 1, mean))
  TF_diff_rna$dev_diff = TF_diff_rna$tumor_dev - TF_diff_rna$normal_dev
  TF_diff_rna$rna_diff = TF_diff_rna$tumor_rna - TF_diff_rna$normal_rna
  
  
  # Compute significance per sample vs normal in scRNA and dev ####
  library (presto)
  pval_threshold = 0.01
  occurrence_threshold = 4
  
  srt$sampleID2 = srt$sampleID
  srt$sampleID2[srt$sampleID2 %in% c('HU37','HU62')] = 'normal'
  rna_comparisons = list(
    P1 = c('P1','normal'),
    P11 = c('P11','normal'),
    P12 = c('P12','normal'),
    P13 = c('P13','normal'),
    P14 = c('P14','normal'),
    P3 = c('P3','normal'),
    P4 = c('P4','normal'),
    P5 = c('P5','normal'),
    P8 = c('P8','normal'))
  
  rna_res = lapply (rna_comparisons, function(x) 
    wilcoxauc (srt[rowData (mSE)$name,], group_by = 'sampleID2', groups_use = x))
  rna_res_df = lapply (names (rna_comparisons), function(x) rna_res[[x]][rna_res[[x]]$group == x,'padj',drop=FALSE])
  rna_res_df = do.call (cbind , rna_res_df)
  rownames (rna_res_df) = rownames (srt[rowData (mSE)$name,])
  colnames (rna_res_df) = names (rna_comparisons)
  
  # Also test genetic and chromatin regulators ####
  rna_res_2 = lapply (rna_comparisons, function(x) 
    wilcoxauc (srt, group_by = 'sampleID2', groups_use = x))
  chromatin_regulators = c('SETD5, ASH1L, CREBBP, PRDM2, KDM2B, KMT2D, EZH2, SETDB1')
  chromatin_regulators = unlist(strsplit (chromatin_regulators, ', '))
  genetic_drivers = c('BAP1, NF2, CDKN2B, CDKN2A, TP53, LATS2, SETD2, FAT4, PTCH1')
  chromatin_regulators2 = c('LATS1, DDX3X, ULK2, RYR2, CFAP45, SETDB1, DDX51, SF3B1, TRAF7, PTEN, RBFOX1, CSMD1, MTAP, TTC28, PCDH15, USH2A, CNTNAP2, DNAH1, KCNH7, PTK2, ROBO2, DLG2, PBRM1, PTPRD, ANTXR2, CTNNA3, LINGO2, LRP1B, PLCB1, UNC79, WWOX')
  chromatin_regulators2 = unlist (strsplit (chromatin_regulators2, ', '))
  genetic_drivers = unlist(strsplit (genetic_drivers, ', '))
  check_genes = unique (c(chromatin_regulators, genetic_drivers, chromatin_regulators2))
  
  rna_res_sum_p2 = apply (do.call (cbind, lapply (rna_res_2, function(x) x[match(check_genes, x$feature),c('padj'), drop=F])), 1, median)
  rna_res_sum_lfc2 = apply (do.call (cbind, lapply (rna_res_2, function(x) x[match(check_genes, x$feature),c('logFC'), drop=F])),1, median)
  rna_res_sum_ae2 = apply (do.call (cbind, lapply (rna_res_2, function(x) x[match(check_genes, x$feature),c('avgExpr'), drop=F])),1, median)
  rna_res_sum2 = data.frame (row.names = check_genes, RNA_pvalue_median = rna_res_sum_p2, RNA_lfc_median = rna_res_sum_lfc2, RNA_avExpr_median = rna_res_sum_ae2)
  write.csv (rna_res_sum2, 'normal_tumor_DEG_chromatin_reg_genetic_drivers.csv')
  
  #occurence_filter = apply (rna_res_df, 1, function(x) sum (x < pval_threshold))
  # rna_res_df_filtered = rna_res_df[occurence_filter > occurrence_threshold, ]
  # tf_tumor_pos = rownames(TF_diff_rna)[TF_diff_rna$dev_diff > 0 & TF_diff_rna$rna_diff > 0]
  # rna_res_df_filtered = rna_res_df_filtered[rownames(rna_res_df_filtered) %in% tf_tumor_pos,]
  # rna_selected_TF = rownames(rna_res_df_filtered)
  
  # Repeat using chromVAR deviations ####
  mMat = assays (mSE)[[1]]
  rownames (mMat) = rowData (mSE)$name
  
  all (colnames(mMat) == rownames(archp@cellColData))
  
  dev_comparisons = list(
    P1 = c('P1','normal_pleura'),
    P11 = c('P11','normal_pleura'),
    P12 = c('P12','normal_pleura'),
    P13 = c('P13','normal_pleura'),
    P14 = c('P14','normal_pleura'),
    P3 = c('P3','normal_pleura'),
    P4 = c('P4','normal_pleura'),
    P5 = c('P5','normal_pleura'),
    P8 = c('P8','normal_pleura'))
  
  dev_res = lapply (dev_comparisons, function(x) 
    wilcoxauc (mMat, y = archp$Sample2, groups_use = x))
  dev_res_df = lapply (names (dev_comparisons), function(x) dev_res[[x]][dev_res[[x]]$group == x,'padj',drop=FALSE])
  dev_res_df = do.call (cbind , dev_res_df)
  rownames (dev_res_df) = rowData (mSE)$name
  colnames (dev_res_df) = names(dev_comparisons)
  #occurence_filter = apply (dev_res_df, 1, function(x) sum (x < pval_threshold))
  # dev_res_df_filtered = dev_res_df[occurence_filter > occurrence_threshold, ]
  # dev_res_df_filtered = dev_res_df_filtered[rownames(dev_res_df_filtered) %in% tf_tumor_pos,]
  # dev_selected_TF = rownames(dev_res_df_filtered)
  # selected_TF = intersect (rna_selected_TF, dev_selected_TF)
  
  rna_res_df_logic = rna_res_df < pval_threshold
  dev_res_df_logic = dev_res_df < pval_threshold
  rna_res_df_logic = rna_res_df_logic[unique (intersect (rownames(dev_res_df_logic), rownames(rna_res_df_logic))),]
  dev_res_df_logic = dev_res_df_logic[unique (intersect (rownames(dev_res_df_logic), rownames(rna_res_df_logic))),]
  comb_res_df = rna_res_df_logic + dev_res_df_logic
  comb_res_df[comb_res_df == 1] = 0
  comb_res_df[comb_res_df == 2] = 1
  selected_TF = rownames(comb_res_df) [rowSums (comb_res_df) > occurrence_threshold] 
  tf_tumor_pos = rownames(TF_diff_rna)[TF_diff_rna$dev_diff > 0 & TF_diff_rna$rna_diff > 0]
  selected_TF = selected_TF[selected_TF %in% tf_tumor_pos]
  
  
  # Order by mean logFC ####
  res_df2 = lapply (names (rna_comparisons), function(x) rna_res[[x]][rna_res[[x]]$group == x,'logFC',drop=FALSE])
  res_df2 = do.call (cbind , res_df2)
  rownames (res_df2) = rownames(mMat_agg)
  tf_order = rownames(res_df2)[order (-apply (res_df2, 1, mean))]
  selected_TF_ordered = tf_order[tf_order %in% selected_TF]
  
  # Export selected TFs ####
  write.csv (selected_TF_ordered, 'Active_TFs.csv')
  
  # Plot distribution of diff deviation and diff expression between tumor and normal ####
  diff_line = 0
  TF_diff_rna$label = ''
  TF_diff_rna$label[match (selected_TF_ordered,rownames(TF_diff_rna))] = selected_TF_ordered
  TF_diff_rna$color = TF_diff_rna$label != ''
  TF_diff_rna$label_top = ''
  TF_diff_rna$label_top[match (head (selected_TF_ordered,20),rownames(TF_diff_rna))] = head (selected_TF_ordered,20)
  TF_diff_rna$alpha = 0.7
  TF_diff_rna$alpha[match (selected_TF_ordered,rownames(TF_diff_rna))] = 1
  tf_diff_p = ggplot (TF_diff_rna, aes (x= dev_diff, y = rna_diff,label = label_top)) + 
  geom_point(aes(color=color, alpha = alpha), size = 2, shape = 21, stroke=0.3) + # Color points based on x value
  scale_color_manual(values = c('FALSE' = "grey", 'TRUE' = "red")) + # Customize colors
  scale_fill_gradient(low = "white", high = "black") +
  #scale_fill_manual(values = c('FALSE' = "grey", 'TRUE' = "red")) + # Customize colors
  geom_vline(xintercept = diff_line, linetype = "dashed", color = "grey44", linewidth=.3) + # Vertical dashed line
  geom_hline(yintercept = diff_line, linetype = "dashed", color = "grey44", linewidth=.3) + # Vertical dashed line
  gtheme_no_rot + # Use a minimal theme
   geom_text_repel(
    segment.size=.05,
    max.overlaps = 100,
#    point.padding = 0.2, 
    size=2#,
#   nudge_x = .25,
#    nudge_y = .2,
#    segment.curvature = -1e-20
    ) +
    xlab ('deviation difference') + 
    ylab ('RNA difference') + 
    xlim (c(-0.2, .2)) + 
    ylim (c(-0.6, .6))

  pdf (paste0 ('Plots/Diff_normal_tumor_deviation_and_rna_scatterplot2.pdf'),5,height = 4)
  print (tf_diff_p)
  dev.off()
  saveRDS(selected_TF, 'selected_TF.rds')
  } else {
  selected_TF = readRDS ('selected_TF.rds')
  }

### Look at overlap of TF with CNV from TCGA ####
## Generate circos plot for showing hubs on recurrent CNVs ####
# source script to load TCGA_CNV data
source (file.path('..','..','git_repo','tumor_analysis','compile_TCGA_CNV_data.R'))

# Find genomic loci of active TFs ####
all_genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
all_genes$gene = as.data.frame(org.Hs.egSYMBOL)[match (all_genes$gene_id, as.data.frame(org.Hs.egSYMBOL)[,1]),'symbol']
all_genes = as.data.frame (all_genes[all_genes$gene %in% selected_TF,])
all_genes = all_genes[, c('seqnames','start','end','gene')]


pdf (file.path ('Plots','CNV_TFs_circos.pdf'))
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
circos.genomicLabels(all_genes, labels.column = 4, side = "inside",
    col = 'grey22', line_col = 'grey22',padding = 0.6, cex=0.5)
dev.off()


# Also get CNV value for each TF ####
all_genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
all_genes$gene = as.data.frame(org.Hs.egSYMBOL)[match (all_genes$gene_id, as.data.frame(org.Hs.egSYMBOL)[,1]),'symbol']
all_genes = as.data.frame (all_genes[all_genes$gene %in% rownames(dev_res_df),])
all_genes = all_genes[, c('seqnames','start','end','gene')]

all_genes_gr = makeGRangesFromDataFrame (all_genes, keep.extra.columns = TRUE)
cnv_mat_avg_gr = makeGRangesFromDataFrame (cnv_mat_avg, keep.extra.columns = TRUE)
tf_idx = subjectHits (findOverlaps(all_genes_gr, cnv_mat_avg_gr))
tf_idx2 = queryHits (findOverlaps(all_genes_gr, cnv_mat_avg_gr))
cnv_mat_avg_gr = cnv_mat_avg_gr[tf_idx]
cnv_mat_avg_gr$gene = all_genes_gr$gene [tf_idx2]
cnv_mat_avg_df = as.data.frame (cnv_mat_avg_gr, row.names=NULL)
cnv_mat_avg_df = do.call (rbind, lapply (split (cnv_mat_avg_df, cnv_mat_avg_df$gene), function(x) data.frame (gene = x$gene[1], CNV_avg_log = mean (x$CNV_avg_log))))
write.csv (cnv_mat_avg_df, 'TF_CNV_TCGA.csv')

# Show differences in TF CNV score between up and down regulated in normal ####
metaGroupName = 'sampleID'
active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)

cnv_mat_avg_df$selected_TFs = cnv_mat_avg_df$gene %in% selected_TF
#cnv_mat_avg_df = cnv_mat_avg_df[cnv_mat_avg_df$gene %in% active_TFs,]

bp = ggplot (cnv_mat_avg_df, aes (x = selected_TFs, y= CNV_avg_log)) + 
bxp + gtheme

stat.test = bp$data %>%
          wilcox_test (reformulate ('selected_TFs', 'CNV_avg_log')) %>%
          adjust_pvalue (method = "none") %>%
          add_significance ()
stat.test = stat.test %>% add_xy_position (x = 'selected_TFs', step.increase=0.5)

    bp = bp + stat_pvalue_manual (stat.test, 
      remove.bracket=FALSE,
      bracket.nudge.y = 0, 
      hide.ns = F,
      label = "p.adj") + 
      NoLegend()
    
pdf (file.path ('Plots','CNV_score_in_selected_TF.pdf'))
bp
dev.off()

### Check selected TF in TCGA ATAC ####
# Add deviations from TCGA ATAC bulk ####
tcga_dev = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/TCGA_atac/tcga_meso_atac_deviations_clinical_info.rds')
tcga_mat = assays (tcga_dev)$z
tcga_mat = aggregate (t(tcga_mat), by = list(sampleID= colnames(tcga_mat)), mean)
rownames (tcga_mat) = tcga_mat[,1]
tcga_mat = tcga_mat[,-1]
tcga_mat = as.data.frame (t(tcga_mat))
rownames (tcga_mat) = unname (sapply (rownames(tcga_mat), function(x) unlist (strsplit (x, '_'))[3]))



# Plot heatmaps of candidate TFs ####
# Heatmap of TF deviations ####
column_split = ifelse (grepl ('normal_pleura', colnames(mMat_agg)), 'Normal','Tumor')
mMat_agg_tf = mMat_agg[selected_TF_ordered,]

TF_cluster_selected_hm = Heatmap (mMat_agg_tf,
        #right_annotation=tf_mark,
        column_split = column_split,
        cluster_rows = F, #km = 4, 
        name = 'z-score\ndeviations',
        column_gap = unit(.8, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T, 
        col = palette_deviation,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 5), 
        border=T,
        width = unit(2, "cm"))

#ps_order = ps[tf_name,sarc_order$sampleID]
# Add scRNA data ####
column_split_rna = ifelse (grepl ('HU', colnames(ps)), 'Normal','Tumor')
ps_tf = ps[selected_TF_ordered,]
TF_exp_selected_hm = Heatmap (ps_tf,
        #right_annotation=tf_mark,
        column_split = column_split_rna,
        cluster_rows = T, #km = 4, 
        name = 'expression',
        column_gap = unit(.5, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T, 
        col = palette_expression,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 5), 
        border=T,
        width = unit(2, "cm"))

# Add regulons enriched for each TF ####
TF_regulons_hm = Heatmap (t(scale(t(nmf_TF_df[selected_TF_ordered,]))),
        #right_annotation=tf_mark,
#        column_split = column_split,
        cluster_rows = F, #km = 4, 
        name = 'expression',
        column_gap = unit(.2, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T,
        column_names_rot = 45, 
        col = viridis::inferno (100),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5), 
        border=F,
        width = unit(5, "cm"))

# Add correlation to bulk RNA sarcomatoid score ####
tf_sarc_cor = read.csv ('../../bulkRNA_meso/activeTF_sarcomatoid_correlation.csv', row.names = 1)
tf_sarc_cor_tf = tf_sarc_cor[selected_TF_ordered,]
TF_bulk_sarc_cor_hm = Heatmap (tf_sarc_cor_tf,
        #right_annotation=tf_mark,
        column_split = colnames (tf_sarc_cor),
        cluster_rows = F, #km = 4, 
        name = 'sarc_cor_bulk',
        column_gap = unit(.2, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T,
        column_names_rot = 45, 
        #col = palette_expression (100),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5), 
        border=F,
        width = unit(.6, "cm"))

# Add cNMF vs TF correlation from scRNA ####


cnmf_tf_cor_hm = Heatmap (t(cnmf_tf_cor)[selected_TF_ordered,],
        #right_annotation=tf_mark,
#        column_split = colnames (cnmf_tf_cor),
        cluster_rows = F, #km = 4, 
        name = 'sarc_cor_bulk',
        column_gap = unit(.2, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T,
        column_names_rot = 45, 
        #col = palette_expression (100),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5), 
        border=F,
        width = unit(6, "cm"))


pdf (paste0 ('Plots/selected_TF_dev_exp_significant_bulkcor_heatmaps2.pdf'), width = 8,height=5)
draw (TF_cluster_selected_hm + TF_exp_selected_hm + cnmf_tf_cor_hm + TF_bulk_sarc_cor_hm)
dev.off()

# Import Viestra archetypes to add to the table below ####
motif_clusters = read.csv ('../../Viestra_motif_clustered.csv')
colnames (motif_clusters)[1] = 'Cluster'
head (motif_clusters)
motif_clusters$Cluster2 = sapply (motif_clusters$Motif, function(x) paste(unique(motif_clusters[motif_clusters$Motif == x,'Cluster']),collapse=','))

# Add correlation to bulk RNA sarcomatoid score ####
tf_sarc_cor = read.csv ('../../bulkRNA_meso/activeTF_sarcomatoid_correlation.csv', row.names = 1)
tf_sarc_cor_tf = tf_sarc_cor[selected_TF_ordered,]

# Export table ####
tf_cnv_tcga = read.csv ('TF_CNV_TCGA.csv')
ccle_exp = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/CCLE/meso_scATAC_active_TFS.csv', row.names=1)
ccle_exp = ccle_exp[selected_TF,]
ccle_exp = ccle_exp[colnames(ccle_exp) != 'gene']
enr_path = t (readRDS ('../scrna/enrichment_pathways_TFs.rds'))
enr_path = enr_path[selected_TF,]
rna_res_sum_p = apply (do.call (cbind, lapply (rna_res, function(x) x[match(selected_TF, x$feature),c('padj'), drop=F])), 1, mean)
rna_res_sum_lfc = apply (do.call (cbind, lapply (rna_res, function(x) x[match(selected_TF, x$feature),c('logFC'), drop=F])),1, mean)
rna_res_sum_ae = apply (do.call (cbind, lapply (rna_res, function(x) x[match(selected_TF, x$feature),c('avgExpr'), drop=F])),1, mean)
rna_res_sum = data.frame (RNA_pvalue_mean = rna_res_sum_p, RNA_lfc_mean = rna_res_sum_lfc, RNA_avExpr_mean = rna_res_sum_ae)

dev_res_sum_p = apply (do.call (cbind, lapply (dev_res, function(x) x[match(selected_TF, x$feature),c('padj'), drop=F])), 1, mean)
dev_res_sum_lfc = apply (do.call (cbind, lapply (dev_res, function(x) x[match(selected_TF, x$feature),c('logFC'), drop=F])), 1, mean)
dev_res_sum_ae = apply (do.call (cbind, lapply (dev_res, function(x) x[match(selected_TF, x$feature),c('avgExpr'), drop=F])), 1, mean)
dev_res_sum = data.frame (DEV_pvalue_mean = dev_res_sum_p, DEV_lfc_mean = dev_res_sum_lfc, DEV_avExpr_mean = dev_res_sum_ae)



tf_table = cbind (
  dev_res_sum, 
  rna_res_sum, 
  tf_sarc_cor_tf[selected_TF,], 
  tf_cnv_tcga[match(selected_TF,tf_cnv_tcga$gene),'CNV_avg_log', drop=F], 
  motif_clusters[match (selected_TF, motif_clusters$Motif),'Cluster2',drop=F],
  ccle_exp,   
  enr_path)
rownames (tf_table) = selected_TF
write.csv (tf_table, 'candidate_TFs_table.csv')






























# Add cNMF identified in scRNA to archr object ####
nfeat=5000
k=25
cnmf_spectra_unique = readRDS (paste0('../scrna/cnmf_genelist_',k,'_nfeat_',nfeat,'.rds'))
cnmf_spectra_unique = cnmf_spectra_unique[!names(cnmf_spectra_unique) %in% c( 'cNMF16','cNMF24')]
cnmf_spectra_unique = lapply (cnmf_spectra_unique, function(x) head(x, 50)[head(x, 50) %in% rownames (gsMat)])

  archp = addModuleScore (
    ArchRProj = archp,
    useMatrix = 'GeneScoreMatrix',
    name = "sarcomatoid",
    features = cnmf_spectra_unique,
    nBin = 25,
    nBgd = 100,
    seed = 1,
    threads = getArchRThreads(),
    logFile = createLogFile("addModuleScore")
  )

cnmfs = archp@cellColData[,grep ('sarcomatoid',colnames(archp@cellColData))]
cnmfs = as.data.frame (t(scale (t(cnmfs))))




# Make coexpression network for each sample using top TFs deviations ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
archp_meta = as.data.frame (archp@cellColData)
archp_meta$Sample3 = archp_meta$Sample2
archp_meta$Sample3[archp_meta$Clusters == 'C14'] = 'P11_HOX'
sams = unique(archp_meta$Sample3)
sams = sams[!sams %in% c('normal_pleura','P3','P13')]

archp_meta = archp_meta[archp_meta$Sample3 %in% sams,]
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat[selected_TF,])

#all (colnames(mMat) == rownames(archp_meta))
cor_TF_l = list()
for (sam in unique(archp_meta$Sample3))
  {
  cor_TF_l[[sam]] = cor (t(as.matrix(scale(mMat[,archp_meta$Sample3 == sam]))), method = 'spearman')
  }

corTF_array <- simplify2array(cor_TF_l)
#any(lapply(corTF_array, function(x) any(is.na(x))))
# Take element-wise median
median_matrix <- apply(corTF_array, c(1, 2), mean)

# Compute correlation of sarcomatoid cNMF with TFs ####
cnmf_mat = as.matrix(archp@cellColData[,grep ('sarcomatoid', colnames(archp@cellColData))])
cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp_meta$Sample3 == x, ])))
names (cnmf_mat) = sams
tf_mat = lapply (sams, function(x) scale(mMat[,archp_meta$Sample3 == x]))
names (tf_mat) = sams

sarc_tf = lapply (sams, function(sam) cor (t(tf_mat[[sam]]), t(cnmf_mat[[sam]]['sarcomatoid.cNMF20',,drop=F])))
sarc_tf_df = do.call (cbind,sarc_tf)

ha1 = HeatmapAnnotation (' ' = anno_boxplot(-sarc_tf_df,axis = T, 
  width = unit(1.2, "cm"), outline=F, border=F,box_width = .8,
    gp = gpar(fill = 'white')), 
  which = 'row')

pdf()
cor_TF_df = draw (Heatmap (median_matrix,
  left= ha1,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 7, fontface='italic'),
  col = palette_deviation_cor_fun,
  rect_gp = gpar(type = "none"),
  cell_fun = function(j, i, x, y, w, h, fill) {
        if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }}))
dev.off()

# Add pathways of TF correlated genes from scRNA
TFrow_order = row_order (cor_TF_df)
TF_cor_sum = readRDS (file.path('..','scrna','enrichment_pathways_TFs.rds'))
rownames (TF_cor_sum) = gsub ('HALLMARK_','',rownames(TF_cor_sum))
TF_cor_sum = TF_cor_sum[apply (TF_cor_sum, 1, function(x) any(x > 3)),]
pdf()
hm = draw (Heatmap (
    t(TF_cor_sum[,TFrow_order]),
    column_names_rot =45, 
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 6),
    col =palette_enrichment, 
    cluster_rows=F,
#     cell_fun = function(j, i, x, y, width, height, fill) {
#         grid.text(sprintf("%.0f", t(TF_cor_sum)[i, j]), x, y, gp = gpar(fontsize = 10, col='white'))
# },
    border=T))
dev.off()

pdf (paste0 ('Plots/selected_TF_dev_corr_heatmaps.pdf'), width = 6,height=5)
cor_TF_df 
dev.off()

pdf (paste0 ('Plots/selected_TF_dev_corr_pathways_heatmaps.pdf'), width = 3,height=6)
hm
dev.off()

# # Plot also scRNA 
# tf_order = row_order (cor_TF_df)

# diag(cor_TF) = 0
# cor_TF2 = cor_TF
# cor_TF2[cor_TF2 > 0.5] = 0.5
# cor_TF2[cor_TF2 < -0.5] = -0.5
# pdf (file.path('Plots','cor_selected_TF_metacells_heatmap.pdf'),width=5, height=4)
# Heatmap (cor_TF[tf_order,tf_order], 
#   column_names_gp = gpar(fontsize = 5),
#   row_names_gp = gpar(fontsize = 5),
#   cluster_rows=T,
#   cluster_columns=T,
#     rect_gp = gpar(type = "none"),
#   # clustering_distance_rows = 'pearson',
#   # clustering_distance_columns='pearson',
#   col = palette_expression_cor_fun(cor_TF2),
#   cell_fun = function(j, i, x, y, w, h, fill) {
#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#                     }})

# dev.off()




# Plot boxplots ordered by correlation to sarcomatoid score using TF activity and scRNA ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat[selected_TF,])

archp_meta = as.data.frame (archp@cellColData)
archp_meta$Sample3 = archp_meta$Sample2
archp_meta$Sample3[archp_meta$Clusters == 'C14'] = 'P11_HOX'
sams = unique(archp_meta$Sample3)
sams = sams[!sams %in% c('normal_pleura','P3','P13')]

# Compute correlation of sarcomatoid cNMF with TFs ####
cnmf_mat = as.matrix(archp@cellColData[,grep ('sarcomatoid', colnames(archp@cellColData))])
cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp_meta$Sample3 == x, ])))
names (cnmf_mat) = sams
tf_mat = lapply (sams, function(x) scale(mMat[,archp_meta$Sample3 == x]))
names (tf_mat) = sams

sarc_tf = lapply (sams, function(sam) cor (t(tf_mat[[sam]]), t(cnmf_mat[[sam]]['sarcomatoid.cNMF20',,drop=F])))
sarc_tf_df = do.call (cbind,sarc_tf)

sarc_tf_med = apply (sarc_tf_df,1, median)
sarc_tf_med = order(-sarc_tf_med)
sarc_tf_df2 = as.data.frame (sarc_tf_df)
colnames(sarc_tf_df2) = names (tf_mat) 
sarc_tf_df2$TF = rownames(sarc_tf_df2)
sarc_tf_df2 = gather (sarc_tf_df2, scS, score, 1:(ncol(sarc_tf_df2) - 1))
sarc_tf_df2$TF = factor (sarc_tf_df2$TF, levels = unique(sarc_tf_df2$TF)[sarc_tf_med])
sarc_tf_df2$type = 'scATAC'

# Compute TF correlation to sarcomatoid in scRNA ####
metacells = readRDS (file.path('..','scrna','metacells.rds'))
nfeat=5000
k=25
cnmf_spectra_unique = readRDS (paste0('../scrna/cnmf_genelist_',k,'_nfeat_',nfeat,'.rds'))

sams = c('P1','P11','P12','P13','P4','P5','P8')

metacells = ModScoreCor (
        seurat_obj = metacells, 
        geneset_list = cnmf_spectra_unique, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'cNMF_', outdir = paste0(projdir,'Plots/'))

metacells_assay = metacells@assays$RNA@layers$data
rownames (metacells_assay) = rownames (srt)
metacells_assay = metacells_assay[selected_TF,]

tc_cor = lapply (sams, function(sam)
  {
  tmp = cor(as.matrix(t(metacells_assay[,metacells$sampleID == sam])), as.matrix(metacells$cNMF20)[metacells$sampleID == sam,,drop=F], method='spearman')
  tmp = as.data.frame (tmp)
  tmp$scS = sam
  tmp$TF = rownames(tmp)
  colnames (tmp) = c('score','scS','TF')
  tmp = tmp[,c('TF','scS','score')]
  })
scrna_tf_cor_df = do.call (rbind, tc_cor)
scrna_tf_cor_df$type = 'scRNA'
combined_df = rbind (sarc_tf_df2, scrna_tf_cor_df)
#lapply (tc_cor, function(x) {x = x['cNMF19',]; head(x[order(-x)],10)})

top_sarc_TF = head(levels (combined_df$TF),20)
combined_df = combined_df[combined_df$TF %in% top_sarc_TF, ]
bp = ggplot (combined_df, aes (x = TF, y = score, fill = type), alpha=.5) + 
geom_boxplot (
    linewidth = .1,
    width=1,
    outlier.alpha = 0.2,
    outlier.shape = NA,
     size=0.5, alpha=0.7
     ) + 
geom_point (position = 'jitter', alpha= 0.2, color = 'grey22', size=1) +
gtheme +
scale_fill_manual (values = c(scATAC = '#B2183BFF', scRNA = '#7663A3FF')) + 
geom_hline(yintercept = 0, color='red',  linetype='dashed')

pdf (paste0 ('Plots/sarcomatoid_score_TF_boxplots.pdf'), width = 7,height=3)
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

pdf (paste0 ('Plots/selected_TF_cooccurence_heatmaps.pdf'), width = 6,height=5)
cooc_hm 
dev.off()

# Check overlap of co-occurring motifs #### WORK IN PROGRESS
sox9 = motifMat[[grep ('SOX9',names(motifMat))]]
sox9p = getPeakSet (archp)[unique(queryHits (findOverlaps (getPeakSet(archp), sox9)))]
sox6 = motifMat[[grep ('SOX6',names(motifMat))]]
sox6p = getPeakSet (archp)[unique(queryHits (findOverlaps (getPeakSet(archp), sox6)))]
peak_ovs = sox9p[queryHits(findOverlaps(sox9p, sox6p))]
#'SERPINE1' %in% unique(peak_ovs$nearestGene)
distance (motifMat[[grep('SOX9', names(motifMat))]][queryHits (findOverlaps (motifMat[[grep('SOX9', names(motifMat))]], peak_ovs[i]))],
motifMat[[grep('SOX6', names(motifMat))]][queryHits (findOverlaps (motifMat[[grep('SOX6', names(motifMat))]], peak_ovs))])







  


# Try with PCA ####
library (uwot)
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat)#[selected_TF,])
tf_mat = lapply (sams, function(x) t(scale(mMat[,archp$Sample3 == x ])))
names (tf_mat) = sams
#tf_mat = do.call (rbind,tf_mat)
cnmf_mat = as.matrix(archp@cellColData[,grep ('sarcomatoid', colnames(archp@cellColData))])
cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp$Sample3 == x, ])))
names(cnmf_mat) = sams

TF_driver='SOX9'

TFclass = list(
stripeTF = c('FOSL2','BACH2','FOSB','JUND','JDP2','BACH1','FOS','JUNB','FOSL1','JUN','SMARCC1'),
lineageTF = c(
  'SOX5','SOX6','SOX9',
  #rownames(mMat)[grep ('HOX',rownames(mMat))],
  #rownames(mMat)[grep ('DLX',rownames(mMat))],
  'GATA4',
  #rownames(mMat)[grep ('NKX',rownames(mMat))],
  'RUNX2'),
emtTF = c('SNAI2','TWIST1','ZEB1','TCF3','TWIST2','HIC2','MESP1','HIC1'),
inflTF = c('IRF2','IRF9','E2F3','E2F7'))

TFclass = setNames (rep (names(TFclass), sapply (TFclass, length)),unlist(TFclass))
TFclass = TFclass[rownames(mMat)]

pca_result <- lapply (tf_mat, function(x) prcomp(t(x), center = TRUE, scale = TRUE))
names(pca_result) = sams

pca_data = lapply (sams, function (sam) 
  {
  tmp = as.data.frame(pca_result[[sam]]$x)[,c('PC1','PC2')]
  umap_result = tmp
  # umap_result <- umap(tmp, n_neighbors = 15, min_dist = 0.1, n_components = 2)
  # umap_result = as.data.frame (umap_result)
  # colnames(umap_result) <- c("UMAP1", "UMAP2")
  umap_result$tf = ''
  umap_result$tf[match (c('SOX9','SOX5','SOX6','RUNX2','GATA4'), rownames(umap_result))] = c('SOX9','SOX5','SOX6','RUNX2','GATA4')
  #umap_result$TF = tf_mat[[sam]][,TF_driver] 
  #umap_result$sarc = cnmf_mat[[sam]]['sarcomatoid.cNMF20',]
  umap_result$class = TFclass
  umap_result$class[is.na(umap_result$class)] = 'ND'
  umap_result$class2 = umap_result$class
  umap_result$class2[umap_result$class2 %in% c('emtTF','inflTF')] = 'nonlineageTF'
  umap_result = umap_result[umap_result$class != 'stripeTF', ]
  
  ggplot(umap_result, aes(x = PC1, y = PC2, label=tf)) +#, colors=class))+#, color = sarc)) +
  geom_point(color='grey',alpha=0.5) + gtheme +
  geom_point(data = umap_result[umap_result$class != 'ND',], aes(color=class),alpha=0.5) +
  geom_smooth(data = umap_result[umap_result$class != 'ND',],
    method = "lm", se = TRUE, aes (color = class2)) + 
  geom_point(data = umap_result[umap_result$tf != '',],size=4) +
  geom_text_repel(data= umap_result[umap_result$tf != '',]) + 
  labs(title = paste("UMAP - ",sam)) + 
  scale_color_manual (values = c(ND = 'grey', emtTF = 'red', inflTF = 'blue', lineageTF= 'green',nonlineageTF = 'purple'))
  })

pdf (file.path('Plots','TF_umap.pdf'), width = 44, height = 6)
wrap_plots (pca_data, ncol=length(sams))
dev.off()


# Correlate module scores with TFs ####
if (!exists('gsSE')) gsSE = fetch_mat(archp, 'GeneScore')
if (!exists('mSE')) gsSE = fetch_mat(archp, 'Motif')
gsMat = assays (gsSE)[[1]]
rownames (gsMat) = rowData (gsSE)$name
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.data.frame (t(as.matrix(scale(mMat[selected_TF,]))))

cnmf_tf_l = list()
for (sam in unique(archp_meta$Sample3))
  {
  cnmf_tf_l[[sam]] = cor (cnmfs[archp_meta$Sample3 == sam,], mMat[archp_meta$Sample3 == sam,], method = 'spearman')
  }

cnmf_tf_array = simplify2array (cnmf_tf_l)
#any(lapply(corTF_array, function(x) any(is.na(x))))
# Take element-wise median
mean_matrix = apply(cnmf_tf_array, c(1, 2), mean)

cnmf_tf_df = Heatmap (mean_matrix,
  row_names_gp = gpar (fontsize = 8),
  #rect_gp = gpar(type = "none"),
  column_names_gp = gpar(fontsize = 8, fontface='italic'))

pdf (paste0 ('Plots/selected_TF_cnmf_corr_heatmaps.pdf'), width = 6,height=5)
cnmf_tf_df
dev.off()






### check correlation of sarcomatoid score to external celltypes ####

# Compute correlation of sarcomatoid cNMF external celltypes ####
archp$Sample3 = archp$Sample2
archp$Sample3[archp$Clusters == 'C14'] = 'P11_HOX'

sams = unique(archp$Sample3)
sams = sams[!sams %in% c('normal_pleura','P3','P13')]
cnmf_mat = as.matrix(archp@cellColData[,grep ('sarcomatoid', colnames(archp@cellColData))])
cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp$Sample3 == x, ])))
names (cnmf_mat) = sams


# # Get external matrix ####
if (!exists('fSE')) fSE = fetch_mat (archp, 'scATAC_datasets')
all (colnames(fSE) == rownames(archp))
fetal_Mat = assays (fSE)[[1]]
rownames (fetal_Mat) = rownames(rowData (fSE))

fetal_Mat = lapply (sams, function(x) scale(fetal_Mat[,archp$Sample3 == x]))
names (fetal_Mat) = sams

sarc_tf = lapply (sams, function(sam) cor (t(fetal_Mat[[sam]]), t(cnmf_mat[[sam]]['sarcomatoid.cNMF8',,drop=F])))
sarc_tf = lapply (sarc_tf, function(x) head(x[order(-x[,1]),],10))

sarc_tf_df = do.call (cbind,sarc_tf)

# ha1 = HeatmapAnnotation (' ' = anno_boxplot(-sarc_tf_df,axis = T, 
#   width = unit(1.2, "cm"), outline=F, border=F,box_width = .8,
#     gp = gpar(fill = 'white')), 
#   which = 'row')

cnmf_fetal_l = list()
for (sam in sams)
  {
  cnmf_fetal_l[[sam]] = cor (t(cnmf_mat[[sam]]), t(fetal_Mat[[sam]]), method = 'spearman')
  }

fetal_array <- simplify2array(cnmf_fetal_l)
#any(lapply(corTF_array, function(x) any(is.na(x))))
# Take element-wise median
median_matrix <- apply(fetal_array, c(1, 2), mean)

cor_fetal_df = Heatmap (median_matrix[rownames(median_matrix) != 'sarcomatoid.cNMF9',],
  row_names_gp = gpar(fontsize = 8),
  #rect_gp = gpar(type = "none"),
  column_names_gp = gpar(fontsize = 7, fontface='italic'),
  #col = palette_deviation_cor_fun
  )

pdf (paste0 ('Plots/selected_fetal_cnmf_corr_heatmaps.pdf'), width = 36,height=5)
cor_fetal_df
dev.off()


# # Get ENCODE matrix ####
if (!exists('enSE')) enSE = fetch_mat (archp, 'ENCODE_H3K4me3')
all (colnames(enSE) == rownames(archp))
encode_Mat = assays (enSE)[[1]]
rownames (encode_Mat) = rownames(rowData (enSE))

cnmf_mat = as.matrix(archp@cellColData[,grep ('sarcomatoid', colnames(archp@cellColData))])
cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp$Sample3 == x, ])))
names (cnmf_mat) = sams

sams = unique(archp$Sample3)
sams = sams[!sams %in% c('normal_pleura','P3','P13')]
encode_Mat = lapply (sams, function(x) scale(encode_Mat[,archp$Sample3 == x]))
names (encode_Mat) = sams

sarc_tf = lapply (sams, function(sam) cor (t(encode_Mat[[sam]]), t(cnmf_mat[[sam]]['sarcomatoid.cNMF20',,drop=F])))
sarc_tf = lapply (sarc_tf, function(x) head(x[order(-x[,1]),],10))
sarc_tf_df = do.call (cbind,sarc_tf)

# ha1 = HeatmapAnnotation (' ' = anno_boxplot(-sarc_tf_df,axis = T, 
#   width = unit(1.2, "cm"), outline=F, border=F,box_width = .8,
#     gp = gpar(fill = 'white')), 
#   which = 'row')

cnmf_encode_l = list()
for (sam in sams)
  {
  cnmf_encode_l[[sam]] = cor (t(cnmf_mat[[sam]]), t(encode_Mat[[sam]]), method = 'spearman')
  }

cnmf_encode_array <- simplify2array(cnmf_encode_l)
#any(lapply(corTF_array, function(x) any(is.na(x))))
# Take element-wise median
median_matrix <- apply(cnmf_encode_array, c(1, 2), mean)

cor_cnmf_encode_df = Heatmap (median_matrix,
  row_names_gp = gpar(fontsize = 8),
  #rect_gp = gpar(type = "none"),
  column_names_gp = gpar(fontsize = 7, fontface='italic'),
  #col = palette_deviation_cor_fun
  )

pdf (paste0 ('Plots/ENCODE_cnmf_corr_heatmaps.pdf'), width = 10,height=7)
cor_cnmf_encode_df
dev.off()


# Also plot boxplots ordered by correlation to sarcomatoid score ####
sarc_tf_med = apply (sarc_tf_df,1, mean)
sarc_tf_order = head (order(-sarc_tf_med),20)
sarc_tf_med = sarc_tf_med[sarc_tf_order]
sarc_tf_df2 = sarc_tf_df[sarc_tf_order,]
sarc_tf_df2 = as.data.frame (sarc_tf_df2)
colnames(sarc_tf_df2) = names (encode_Mat) 
#sarc_tf_df2$TF = rownames(sarc_tf_df2)
sarc_tf_df2$ENCODE_celltype = factor (names(sarc_tf_med), levels = names (sarc_tf_med))
sarc_tf_df2 = gather (sarc_tf_df2, scS, score, 1:(ncol(sarc_tf_df2)-1))
bp = ggplot (sarc_tf_df2, aes (x = ENCODE_celltype, y = score)) + geom_boxplot () + gtheme

pdf (paste0 ('Plots/sarcomatoid_score_ENCODE_boxplots.pdf'), width = 6,height=5)
bp
dev.off()








# Show sarcomatoid module ####
#metaGroupName = "Clusters"
#force=FALSE
#if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source ('../../PM_scATAC/DAG.R')

celltype_markers = 'sarcomatoid.cNMF20'

pdf()
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = celltype_markers, 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()

pdf (file.path('Plots','sarcomatoid_score_feature_plots.pdf'), width = 20, height = 20)
print (wrap_plots (p, ncol = 4))
dev.off()

# Plot the same but scaled ####
scaled_cnmfs = archp@cellColData[grep ('sarcomatoid',colnames(archp@cellColData))]
scaled_cnmfs = as.data.frame (t(scale(t(scaled_cnmfs))))
umap_df = data.frame (archp@embeddings$UMAP[[1]], scaled_cnmfs)
umap_p1 = ggplot(data = umap_df) + 
geom_point (aes (x = IterativeLSI.UMAP_Dimension_1, y= IterativeLSI.UMAP_Dimension_2, color = sarcomatoid.cNMF20), size = .1) + 
scale_colour_gradientn (colours = rev(brewer.pal (n = 11, name = "RdBu")),limits=c(-max (abs(umap_df$sarcomatoid.cNMF20)), max (abs(umap_df$sarcomatoid.cNMF20)))) +
ggtitle ('sarcomatoid_score') + 
#facet_wrap (as.formula(paste("~", metaGroupNames[3]))) + 
theme_classic() +
theme_void()

pdf (file.path('Plots','sarcomatoid_score_scaled_feature_plots.pdf'), width = 6, height = 6)
umap_p1
dev.off()
  
# Order cells per samples along SOX9 deviation and correlated TFs ####
library (scales)
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
# Filter by RNA expression ####
metaGroupName = 'sampleID'
active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)
mMat = mMat[active_TFs, ]

archp_meta = as.data.frame (archp@cellColData)
archp_meta$Sample3 = archp_meta$Sample2
archp_meta$Sample3[archp_meta$Clusters == 'C14'] = 'P11_HOX'

all (colnames(mMat) == rownames(archp_meta))

TF_driver = 'SOX9'
top_TFs = 50
traj_sample = list()
for (sam in unique(archp_meta$Sample3))
    {
    library(zoo)
    bin_width <- 100   # Number of observations per bin
    overlap <- 1    
    mMat_ordered_sample = as.data.frame(scale(mMat[,archp_meta$Sample3 == sam]))
    mMat_ordered_sample = mMat_ordered_sample[, order(unlist(mMat_ordered_sample[TF_driver,]))]
    cor_to_tf = order(-as.vector(cor (t(mMat_ordered_sample)[,TF_driver],t(mMat_ordered_sample), method='pearson')))
    mMat_ordered_sample = mMat_ordered_sample[c(head(cor_to_tf,top_TFs),tail(cor_to_tf,top_TFs))  ,]
    mMat_ordered_sample = as.data.frame(lapply(as.data.frame(t(mMat_ordered_sample)), function(x) {
      zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
    }))
    traj_sample[[sam]] = Heatmap (
      t(as.data.frame(lapply(mMat_ordered_sample, rescale, to = c(-10,10)))), 
      col = rev(palette_deviation), 
      cluster_columns=F,
      name = sam,
      cluster_rows=F,
      column_names_gp = gpar(fontsize = 0),
      row_names_gp = gpar(fontsize = 7, fontface='italic'))
    }

pdf (file.path('Plots',paste0('sarc_trajectory_',TF_driver,'_sample.pdf')), height=12, width=3)
traj_sample
dev.off()



# Plot UMAP per sample showing sarcomatoid score and correlated TFs e.g. SOX9 HIC2.. ####
sams = unique(archp$Sample2)
sams = sams[!sams %in% c('normal_pleura','P3','P13')]
cnmf_mat = as.matrix(archp@cellColData[,grep ('sarcomatoid', colnames(archp@cellColData))])
cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp$Sample2 == x, ])))
names(cnmf_mat) = sams
#cnmf_mat = do.call (rbind, lapply (cnmf_mat, function(x) t(x)))
#umap_df = data.frame(archp@embeddings$UMAP[[1]], archp$Sample2)

varfeat = 1000
LSI_method = 2
pdf()
sample_LSI = lapply (sams, function(x) {
  tmp = addIterativeLSI (ArchRProj = archp[archp$Sample2 == x],
    useMatrix = "MotifMatrix", name = "IterativeLSI",
    force = TRUE, LSIMethod = LSI_method,
    varFeatures = varfeat)
    addUMAP (ArchRProj = tmp, 
    reducedDims = "IterativeLSI",
    force = TRUE)@embeddings$UMAP[[1]]
  })
dev.off()


names (sample_LSI) = sams
umap_df = lapply (sams, function(x) data.frame (sample_LSI[[x]][rownames(t(cnmf_mat[[x]])),], t(cnmf_mat[[x]])))
umap_df = do.call (rbind, umap_df)
umap_df$Sample2 = sapply (rownames(umap_df), function(x) unlist(strsplit(x, '\\#'))[1])

umap_p1 = ggplot (data = umap_df) + 
geom_point (aes (x = IterativeLSI.UMAP_Dimension_1, y= IterativeLSI.UMAP_Dimension_2, color = sarcomatoid.cNMF20), size = .01) + 
scale_colour_gradientn (colours = rev(palette_deviation)) +
ggtitle ('sarcomatoid_score') +
scale_size (range = c(0.1, .4)) + 
facet_wrap (~Sample2) + 
theme_void ()

pdf (file.path('Plots','sarcomatoid_score_scaled_per_sample_feature_plots.pdf'), width = 10, height = 5)
umap_p1
dev.off()

  









### Correlate HOX genes with CNV in P11 HOX+ cluster ####
p11cnv = readRDS (file.path('..','..','per_sample_QC_signac','CNV_analysis','CNV_LFC_GC_P11_ws_1e+07_ss_5e+06.rds'))
p11cnv_mat = assays(p11cnv)$counts
p11cnv_mat = scale (p11cnv_mat)

if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
# # Get deviation matrix ####
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Get active TFs ####
selected_TF = read.csv ('Active_TFs.csv', row.names=1)
p11mMat = mMat[, archp$Clusters == 'C14']
#p11mMat = mMat[, archp$Sample2 == 'P11']
colnames(p11mMat) = sapply (colnames(p11mMat), function(x) unlist(strsplit(x, '#'))[2])
p11cnv_mat = p11cnv_mat[,colnames (p11mMat)]

#p11mMat = p11mMat[grepl ('^', rownames(p11mMat)), ]

# correlate
p11mMat = t(p11mMat)
p11cnv_mat = t(p11cnv_mat)

cnv_hox_cor = cor (as.matrix(p11mMat), as.matrix(p11cnv_mat), method='spearman')

hm = Heatmap  (cnv_hox_cor, cluster_columns = F, 
  column_names_gp = gpar(fontsize = 4),
  row_names_gp = gpar(fontsize = 4))
pdf (file.path ('Plots','cnv_HOX_P11_cor_heatmap.pdf'), width=20, height=30)
hm
dev.off()

min_cor = cnv_hox_cor['HOXB13',]
min_cor = min_cor[order(min_cor)]

pdf (file.path ('Plots','cnv_TF_cor.pdf'))
plot (p11mMat[,'HOXB13'], p11cnv_mat[,names(min_cor)[1]])
dev.off()
cor (p11mMat[,'HOXB13'], p11cnv_mat[,names(min_cor)[1]], method = 'spearman')

cnv_tf_cor = p11cnv_mat[,names(min_cor)[1]]
names(cnv_tf_cor) = paste0('P11#',names(cnv_tf_cor))
archp$p11_hox_cnv = cnv_tf_cor[match(rownames(archp@cellColData), names(cnv_tf_cor))]
archp$p11_hox_cnv[is.na(archp$p11_hox_cnv)] = 0
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = 'p11_hox_cnv', 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)

pdf (file.path('Plots','cnv_TF_cor_umap.pdf'))
p
dev.off()



### Try with genescore ####
p11cnv = readRDS (file.path('..','..','per_sample_QC_signac','CNV_analysis','CNV_LFC_GC_P11_ws_1e+07_ss_5e+06.rds'))
p11cnv_mat = assays(p11cnv)$counts
p11cnv_mat = scale (p11cnv_mat)

if (!exists('gsSE')) gsSE = fetch_mat(archp, 'GeneScore')
# # Get deviation matrix ####
gsMat = assays (gsSE)[[1]]
rownames (gsMat) = rowData (gsSE)$name

# Get active TFs ####
selected_TF = read.csv ('Active_TFs.csv', row.names=1)
p11mMat = gsMat[, archp$Clusters == 'C14']
#p11mMat = mMat[, archp$Sample2 == 'P11']
colnames(p11mMat) = sapply (colnames(p11mMat), function(x) unlist(strsplit(x, '#'))[2])
p11cnv_mat = p11cnv_mat[,colnames (p11mMat)]

p11mMat = p11mMat[grepl ('^HOX', rownames(p11mMat)), ]

# correlate
p11mMat = t(p11mMat)
p11cnv_mat = t(p11cnv_mat)

cnv_hox_cor = cor (as.matrix(p11mMat), as.matrix(p11cnv_mat), method='spearman')

hm = Heatmap  (cnv_hox_cor, cluster_columns = F, 
  column_names_gp = gpar(fontsize = 4),
  row_names_gp = gpar(fontsize = 4))
pdf (file.path ('Plots','cnv_HOX_P11_cor_gs_heatmap.pdf'), width=7, height=5)
hm
dev.off()

min_cor = cnv_hox_cor['HOXB13',]
min_cor = min_cor[order(min_cor)]

pdf (file.path ('Plots','cnv_TF_cor_gs.pdf'))
plot (p11mMat[,'HOXB13'], p11cnv_mat[,names(min_cor)[1]])
dev.off()
cor (p11mMat[,'HOXB13'], p11cnv_mat[,names(min_cor)[1]], method = 'spearman')


# Motif enrichment of peaks found around HOX13 genes ####
tf_match = getMatches (archp)
C14_peaks = readRDS (file.path ('PeakCalls','C14-reproduciblePeaks.gr.rds'))
tf_match = tf_match[queryHits (findOverlaps(tf_match, C14_peaks))]
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)
region = GRanges (c(
  'chr17:48722763-48730749',
  'chr2:176092139-176093775',
  'chr7:27192364-27202091',
  'chr12:53936831-53948544'))

region_peaks = bg_peakSet[queryHits(findOverlaps(bg_peakSet, region))]
#tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
region_TF =  hyperMotif (
  selected_peaks = region_peaks, 
  motifmatch = tf_match)

head (region_TF, 40)

scrna_cor = read.csv (file.path('..','scrna','correlated_genes_p11s_HOXB13.csv'))

cor_df = data.frame (rna_cor = scrna_cor$x[match(rownames(region_TF), scrna_cor$X)],
  TF_enrich_log10 = -log10(region_TF$padj), 
  TF_enrich_padj = region_TF$padj,
  row.names = rownames(region_TF)
  )
cor_df$label = ifelse (cor_df$TF_enrich_padj < 0.05, rownames(cor_df), '')
gp = ggplot (cor_df, aes (x = TF_enrich_log10, y = rna_cor, label = label)) + 
  geom_point(size = 2, shape = 21, stroke=0.3) + 
  geom_text_repel() + 
  gtheme_no_rot

pdf (file.path ('Plots','scrna_cor_vs_TF_enrich_HOX_scatter.pdf'),3,3)
gp
dev.off()






### Label tumor subtypes for chrombpnet ####
archp$Sample3 = archp$Sample2
archp$Sample3[archp$Clusters == 'C14'] = 'P11_HOX'

cnmfs = archp@cellColData[,grep ('sarcomatoid',colnames(archp@cellColData))]
cnmfs_scaled = lapply (unique(archp$Sample2), function(x) as.data.frame (t(scale (t(cnmfs)))))
cnmfs_scaled = do.call (rbind,cnmfs_scaled)

cnmfs_scaled = cnmfs_scaled[rownames(archp@cellColData),]
all (rownames (cnmfs_scaled) == rownames(archp@cellColData))
archp$epit_sarc2 = ifelse (cnmfs$sarcomatoid.cNMF20 > mean (cnmfs$sarcomatoid.cNMF20), 'sarc','epit')
archp$epit_sarc2 = paste0(archp$Sample3, '__',archp$epit_sarc2)
table (archp$epit_sarc2, archp$Clusters)

archp$epit_sarc = 'NS'
archp$epit_sarc[archp$Clusters == 'C15'] = 'sarc'
archp$epit_sarc[archp$Clusters == 'C16'] = 'epit'

pdf (file.path ('Plots','sarc_epit_umap.pdf'),width=10)
p1 <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = 'epit_sarc', 
    embedding = "UMAP",
    #pal = palette_expression,
    #imputeWeights = getImputeWeights(archp)
)
p2 <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = 'epit_sarc2', 
    embedding = "UMAP",
    #pal = palette_expression,
    #imputeWeights = getImputeWeights(archp)
)
wrap_plots (p1,p2)
dev.off()


### Check correlation with seq-depth for module scores ####
cnmf = archp@cellColData[,grep ('sarcomatoid',colnames(archp@cellColData))]
cnmf_scaled = as.data.frame (t(scale (t(cnmf))))

metaGroupNames = c('FRIP','nFrags','TSSEnrichment','nMonoFrags','nDiFrags','nMultiFrags','ReadsInTSS','ReadsInPromoter','PromoterRatio','ReadsInPeaks','NucleosomeRatio')

# hp = list()
# for (metagroupname in metaGroupNames)
# {
# cnmf_cor = data.frame (cor = cor (as.matrix(cnmf), as.matrix(archp@cellColData[,metagroupname])), tr = 'not scaled')
# cnmf_cor_scaled = data.frame (cor = cor (as.matrix(cnmf_scaled), as.matrix(archp@cellColData[,metagroupname])), tr = 'scaled')
# cnmf_df = rbind (cnmf_cor, cnmf_cor_scaled)
# hp[[metagroupname]] = ggplot (cnmf_df, aes (y = cor, x = tr)) +
#   #geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
#   geom_boxplot(alpha=.2)+
#   stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "red") +
#   geom_jitter(width = 0.2, alpha = 1, color = "blue") +
#   geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 1) + 
#   labs(title = paste0('cNMF correlation to ',metagroupname,' Distributions'), x = "Value", y = "Frequency") +
#   scale_fill_manual(values = c("blue", "red")) +
#   theme_minimal()
# } 

# pdf (file.path ('Plots','histogram_cnmf_cor_FRIP.pdf'), width=10, height=5)
# wrap_plots(hp)
# dev.off()





hps_gs = list()
for (metagroupname in metaGroupNames)
{
cnmf_cor_sample_scaled = lapply (unique(archp$Sample2), function(x) data.frame (
  cor = cor (as.matrix(t(scale(t(cnmf[archp$Sample2 == x,])))), 
  as.matrix(archp@cellColData[,metagroupname][archp$Sample2 == x])),
  tr = 'scaled', sample = x))

cnmf_cor_sample_scaled_df = do.call (rbind, cnmf_cor_sample_scaled)
cnmf_sample_cor = lapply (unique(archp$Sample2), function(x) data.frame (
  cor = cor (as.matrix(cnmf[archp$Sample2 == x,]), 
  as.matrix(archp@cellColData[,metagroupname][archp$Sample2 == x])), 
tr = 'not scaled', sample = x))
cnmf_sample_cor_df = do.call (rbind, cnmf_sample_cor)
cnmf_df = rbind (cnmf_cor_sample_scaled_df, cnmf_sample_cor_df)
hps_gs[[metagroupname]] = ggplot (cnmf_df, aes (y = cor, x = tr, fill = sample)) +
  #geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  geom_boxplot(alpha=.2,outlier.shape = NA)+
  #stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "red") +
  #geom_jitter(position=position_jitterdodge(), width = 0.2, alpha = 0.4, color = "blue", size=1) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 1) + 
  labs(title = paste0('cNMF correlation to ',metagroupname,' Distributions'), x = "Value", y = "Frequency") +
#  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal()
} 

pdf (file.path ('Plots','histogram_cnmf_cor_FRIP_per_sample.pdf'), width=25, height=15)
wrap_plots(hps_gs, ncol = 5)
dev.off()

### Try with Deviations ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
all (colnames(mSE) == rownames(archp))
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name


hps = list()

for (metagroupname in metaGroupNames)
{
cnmf_cor_sample_scaled = lapply (unique(archp$Sample2), function(x) data.frame (
  cor = cor (as.matrix(t(scale(t(t(mMat)[archp$Sample2 == x,])))), 
  as.matrix(archp@cellColData[,metagroupname][archp$Sample2 == x])),
  tr = 'scaled', sample = x))

cnmf_cor_sample_scaled_df = do.call (rbind, cnmf_cor_sample_scaled)
cnmf_sample_cor = lapply (unique(archp$Sample2), function(x) data.frame (
  cor = cor (as.matrix(t(mMat)[archp$Sample2 == x,]), 
  as.matrix(archp@cellColData[,metagroupname][archp$Sample2 == x])), 
tr = 'not scaled', sample = x))
cnmf_sample_cor_df = do.call (rbind, cnmf_sample_cor)
cnmf_df = rbind (cnmf_cor_sample_scaled_df, cnmf_sample_cor_df)
hps[[metagroupname]] = ggplot (cnmf_df, aes (y = cor, x = tr, fill = sample)) +
  #geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  geom_boxplot(alpha=.2,outlier.shape = NA)+
  #stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "red") +
  #geom_jitter(position=position_jitterdodge(), width = 0.2, alpha = 0.4, color = "blue", size=1) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 1) + 
  labs(title = paste0('Deviation correlation to ',metagroupname,' Distributions'), x = "Value", y = "Frequency") +
#  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal()
} 

pdf (file.path ('Plots','histogram_deviations_cor_FRIP_per_sample.pdf'), width=25, height=15)
wrap_plots(hps, ncol=5)
dev.off()


head (hps$ReadsInTSS$data[hps$ReadsInTSS$data$tr == 'not scaled',][order (-hps$ReadsInTSS$data$cor[hps$ReadsInTSS$data$tr == 'not scaled']),],50)
hps$FRIP$data[grep ('JUN',rownames(hps$FRIP$data)),]


