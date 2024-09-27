conda activate meso_scatac
R

set.seed(1234)

####### ANALYSIS of TUMOR compartment #######
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

# Load functions for hub detection ####
source (file.path('..','git_repo','utils','knnGen.R'))
source (file.path('..','git_repo','utils','addCoax.R'))
source (file.path('..','git_repo','utils','Hubs_finder.R'))
source (file.path('..','git_repo','utils','hubs_track.R'))

# Set # of threads and genome reference ####
addArchRThreads(threads = 8) 
addArchRGenome("hg38")

if (!file.exists ('Save-ArchR-Project.rds')) 
  { source ('../../PM_scATAC/scatac_tumor_create_ArchRobj.R')
  } else {
 archp = loadArchRProject (projdir)   
  }

# Load RNA ####
srt = readRDS ('../scrna/srt.rds')
sarc_order = read.csv ('../scrna/cnmf20_sarcomatoid_sample_order.csv', row.names=1)



# Run genescore DAG ####
metaGroupName = "Clusters"
force=FALSE
if(!file.exists (paste0('DAG_',metaGroupName,'.rds') | force) source ('../../PM_scATAC/DAG.R')

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



# Find activating and repressing TFs #### 
if (!file.exists ('TF_activators_genescore.rds')) 
  {
    source (file.path('..','..','git_repo','utils','activeTFs.R'))
  } else {
    corGSM_MM = readRDS ('TF_activators_genescore.rds') 
  }



### Co-expression of TFs #### 
metaGroupName = 'Sample2'
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
all (colnames(mSE) == rownames(archp))

# # Get deviation matrix ####
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Subset only for positively correlated TF with genescore ####
positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0.1]
mMat = mMat[positive_TF,]

mMat_cor = cor (as.matrix(t(mMat)), method = 'pearson')

km = kmeans (mMat_cor, centers=5)

cor_mMat_hm = draw (Heatmap (mMat_cor,# row_km=15,
  #left_annotation = ha,
  #rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=palette_module_correlation_fun, 
  row_split = km$cluster,
  column_split = km$cluster,
  #row_km=2, 
  #column_km=2,
#  right_annotation = ha,
  border=F,
  row_names_gp = gpar(fontsize = 0),
  column_names_gp = gpar(fontsize = 0)))#,
  # ,
  # cell_fun = function(j, i, x, y, w, h, fill) {
  #       if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
  #           grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#        }}))

pdf (file.path ('Plots','TF_modules.pdf'), width = 4,height=3)
cor_mMat_hm
dev.off()

tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = do.call (cbind, tf_modules)
archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% paste0('mod_',unique(km$cluster))]
archp@cellColData = cbind (archp@cellColData, tf_modules) 

TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "cellColData",
    name = paste0('mod_',unique(km$cluster)), 
    pal = rev(palette_deviation),
    #useSeqnames='z',
    embedding = "UMAP")

pdf (file.path ('Plots','TF_modules_umap.pdf'), width = 20,height=6)
wrap_plots (TF_p, ncol=5)
dev.off()

# # Make radar plot to show oncogenic signatures per sample ####
# tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
# names (tf_modules) = paste0('mod_',unique(km$cluster))
# tf_modules = do.call (cbind, tf_modules)
# tf_modules_2 = apply (tf_modules, 1, function(x) colnames(tf_modules)[which.max (x)])
# tf_modules_2 = split (tf_modules_2, archp$Sample2)
# tf_modules_2 = lapply (tf_modules_2, function(x) 
#   {
#   x = as.data.frame(table(x))
#   rownames (x) = x$x
#   x = x[paste0('mod_',unique(km$cluster)),]
#   rownames (x) = paste0('mod_',unique(km$cluster))
#   as.data.frame (t(x[,2]))
#   })
# tf_modules_2 = do.call (rbind, tf_modules_2)
# colnames (tf_modules_2) = paste0('mod_',unique(km$cluster))
# tf_modules_2[is.na(tf_modules_2)] = 0

# tf_modules_2 = as.data.frame(t(apply (tf_modules_2, 1, proportions)))

# pdf (file.path ('Plots','TF_modules_radar_plots.pdf'), width = 10,height=10)
# for (i in rownames (tf_modules_2))
#   {
#   maxmin = as.data.frame (t(data.frame (rep(max (tf_modules_2),length(colnames (tf_modules_2))), 0)  ))
#   colnames (maxmin) = colnames (tf_modules_2)
#   tf_modules_3 = rbind (maxmin, tf_modules_2[i,, drop=F])

#   radarchart (tf_modules_3,
#       axistype=0 , 
#       maxmin=T,
#       cglcol="grey", axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
#       #custom polygon
#       pcol=palette_sample , pfcol=alpha(palette_sample,0.3) , plwd=1 , plty=1,
#       #custom the grid
#       cglty=1,
#       #custom labels
#       vlcex=0.8 
#       )
#   legend(x=0.7, y=1.3, legend = rownames(tf_modules_2), bty = "n", pch=20 , col=palette_sample , text.col = "grey", cex=.6, pt.cex=2)
#   }

# dev.off()



# Try with ridge plots ####
library(ggridges)
library(ggplot2)
library(viridis)
#library(hrbrthemes)
tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = as.data.frame (do.call (cbind, tf_modules))
tf_modules$Sample = archp$Sample2

# Plot
rp = lapply (paste0('mod_',unique(km$cluster)), function(x) 
  ggplot(tf_modules, aes_string(x = x, y = 'Sample', fill = '..x..')) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, alpha=.5) +
  paletteer::scale_fill_paletteer_c("grDevices::PuRd", direction=-1) +
    theme_classic())
pdf (file.path ('Plots','TF_modules_ridge_plots.pdf'), width = 20,height=3)
wrap_plots (rp, ncol=5)
dev.off()


# Check expression of SOX9 and others - SOX9 should cluster with module 5 but it doesnt probably cause of low accuracy of deviation ####
tf_name2 = 'SOX9_756'

tf_name2 = unlist(sapply (c('SOX9','TWIST1','MESP1','NKX2-5'), function(x) rownames(assay(mSE))[grepl (x, rownames(assay(mSE)))]))
tf_name2 = paste0('z:',tf_name2)
archp = addImputeWeights (archp)
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "MotifMatrix",
    name = tf_name2, 
    useSeqnames='z',
    col = palette_deviation,    
    embedding = "UMAP",
    imputeWeights = getImputeWeights(archp)
    )

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
  

#   # Import cNMF results and intersect with p2g ####
#   nfeat=5000
#   k=25
#   cnmf_list = readRDS (paste0('../scrna/cnmf_genelist_',k,'_nfeat_',nfeat,'.rds'))
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
metaGroupName = 'Sample2'
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
DefaultAssay(srt) = 'RNA'
sample_names_rna = c('P1','P14','P13','P3','P12','P5','P11','P4','P8','P14','HU37','HU62')
ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat_agg), group.by = metaGroupName)[[1]]) +1)
ps = ps[, colnames(ps) %in% sample_names_rna]

TF_diff_rna = data.frame (
  tumor_dev = apply (mMat_agg[,unique(archp$Sample2)[!unique(archp$Sample2) == 'normal_pleura']], 1, mean),
  normal_dev = mMat_agg[,unique(archp$Sample2)[unique(archp$Sample2) == 'normal_pleura']],
  normal_rna = apply (ps[,c('HU37','HU62')], 1, mean),
  tumor_rna = apply (ps[, !colnames(ps) %in% c('HU37','HU62')], 1, mean),
  genescore = corGSM_MM$cor[match(rownames(mMat_agg), corGSM_MM$GeneScoreMatrix_name)])
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
corGSM_MM_filtered = as.data.frame (corGSM_MM [match (selected_TF_ordered, corGSM_MM$GeneScoreMatrix_name),])
head (corGSM_MM_filtered [order (corGSM_MM_filtered$cor),c('GeneScoreMatrix_name','cor') ],100)

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
  geom_point(aes(fill=genescore, color=color, alpha = alpha), size = 2, shape = 21, stroke=0.3) + # Color points based on x value
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
tf_diff_p
dev.off()


### Look at overlap of TF with CNV from TCGA ####
## Generate circos plot for showing hubs on recurrent CNVs ####
# source script to load TCGA_CNV data
source ('../../PM_scATAC/compile_TCGA_CNV_data.R')

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
all_genes_gr = makeGRangesFromDataFrame (all_genes, keep.extra.columns = TRUE)
cnv_mat_avg_gr = makeGRangesFromDataFrame (cnv_mat_avg, keep.extra.columns = TRUE)
tf_idx = subjectHits (findOverlaps(all_genes_gr, cnv_mat_avg_gr))
tf_idx2 = queryHits (findOverlaps(all_genes_gr, cnv_mat_avg_gr))
cnv_mat_avg_gr = cnv_mat_avg_gr[tf_idx]
cnv_mat_avg_gr$gene = all_genes_gr$gene [tf_idx2]
cnv_mat_avg_df = as.data.frame (cnv_mat_avg_gr, row.names=NULL)
cnv_mat_avg_df = do.call (rbind, lapply (split (cnv_mat_avg_df, cnv_mat_avg_df$gene), function(x) data.frame (gene = x$gene[1], CNV_avg_log = mean (x$CNV_avg_log))))
write.csv (cnv_mat_avg_df, 'TF_CNV_TCGA.csv')


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
cnmf_tf_cor = readRDS ('../scrna/cNMF_TF_correlation_subsampled.rds')

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











































# Make coexpression network for each sample using top TFs deviations
archp_meta = as.data.frame (archp@cellColData)
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

mMat = mMat[selected_TF_ordered,]
all (colnames(mMat) == rownames(archp_meta))
cor_TF_l = list()
for (sam in unique(archp_meta$Sample2))
  {
  cor_TF_l[[sam]] = cor (t(as.matrix(mMat[,archp_meta$Sample2 == sam])))
  cor_TF_l[[sam]] = Heatmap (cor_TF_l[[sam]], name = sam,
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 5))
  }

pdf (paste0 ('Plots/selected_TF_dev_corr_heatmaps.pdf'), width = 8,height=9)
cor_TF_l
dev.off()

# Take median of all correlations
# Sample data: list of matrices
all (colnames(mMat) == rownames(archp_meta))
cor_TF_l = list()
for (sam in unique(archp_meta$Sample2)) cor_TF_l[[sam]] = cor (t(as.matrix(mMat[,archp_meta$Sample2 == sam])))

cor_TF <- cor_TF_l[[1]]
cor_TF[] <- tapply(unlist(cor_TF_l), rep(seq(length(cor_TF_l[[1]])),length(cor_TF_l)), FUN=median)

cor_TF = Heatmap (cor_TF,
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 5))

pdf (paste0 ('Plots/selected_TF_dev_corr_heatmap.pdf'), width = 8,height=9)
cor_TF
dev.off()
 





# Order cells per samples along SOX9 deviation and plot the rest of TF deviations together
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = mMat[tf_name_selected,]
archp_meta = as.data.frame (archp@cellColData)
all (colnames(mMat) == rownames(archp_meta))

traj_sample = list()
for (sam in unique(archp_meta$Sample2))
    {
    mMat_ordered_sample = mMat[,archp_meta$Sample2 == sam]
    mMat_ordered_sample = mMat_ordered_sample[, order(mMat_ordered_sample['MESP1',])]
    traj_sample[[sam]] = Heatmap (
      t(scale(t(mMat_ordered_sample))), 
      col = viridis::plasma(100), 
      cluster_columns=F,
      row_names_gp = gpar(fontsize = 4))
    }

pdf ('Plots/sarc_trajectory_per_sample.pdf', height=10)
traj_sample
dev.off()



### Check for TF which co-occur in same peaks
selected_TF_ordered
tf_match = getMatches (archp)
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
colnames (tf_match) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames (tf_match))

tf_match = tf_match[, selected_TF_ordered]
bg_peakSet = rowRanges (tf_match)
  #tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
  




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



