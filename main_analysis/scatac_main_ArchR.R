conda activate meso_scatac
R

set.seed(1234)

packages = c(
  'Signac',
  'Seurat',
  'biovizBase',
  'ggplot2',
  'patchwork',
  'scATACutils',
  'SummarizedExperiment',
  'epiAneufinder',
  'JASPAR2020',
  'TFBSTools',
  'TxDb.Hsapiens.UCSC.hg38.knownGene',
  'EnsDb.Hsapiens.v86',
  'gplots',
  'regioneR',
  'ComplexHeatmap',
  'ArchR',
  'BSgenome.Hsapiens.UCSC.hg38',
  'tidyverse',
  'ggrepel',
  'RColorBrewer')
lapply(packages, require, character.only = TRUE)

####### ANALYSIS of main data #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)


#devtools::install_github("immunogenomics/presto") #needed for DAA
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))
source (file.path('..','..','git_repo','utils','hubs_track.R'))

set.seed (1234)
addArchRThreads (threads = 1) 
addArchRGenome ("Hg38")

sample_names = c(
    # Tumor  
    'P1', # p786
    'P3', # p846
    'P4', # p811
    'P5', #'p848'
    'P8', # p826
    'P10', # p10
    'P11', # p11
    'P12', # p12
    'P13', # p13
    'P14',#,# p14
    'P23'
    )

# Load RNA
srt = readRDS ('../scrna/srt.rds')
srt$celltype_simplified2[srt$celltype_simplified2 == 'pDC'] = 'pDCs'
#sarc_order = read.csv ('../scrna/cnmf20_sarcomatoid_sample_order.csv', row.names=1)

archp = loadArchRProject (projdir)

# Set order of celltype for displaying purposes ####
celltype_order = c('Malignant','Mesothelium','Alveolar','Fibroblasts','SmoothMuscle','Endothelial','Myeloid','T_cells','NK','B_cells','Plasma','pDCs')

#sarc_order = c('P1','P13','P3','P12','P5','P11','P4','P8','P14','P10')
#archp$Sample2 = archp$Sample
#archp$Sample2 = factor (archp$Sample2, levels = sarc_order)


  ### QC plots ####
  pdf()
  p1 = plotFragmentSizes(ArchRProj = archp, groupBy = 'Sample', pal = palette_sample)
  p2 = plotTSSEnrichment(ArchRProj = archp, groupBy = 'Sample', pal = palette_sample)

  p3 <- plotGroups(
    ArchRProj = archp, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    pal = palette_sample,
    alpha = 0.4,
    addBoxPlot = TRUE
   )

  p4 <- plotGroups(
    ArchRProj = archp, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "nFrags",
    plotAs = "violin",
    pal = palette_sample,
    alpha = 0.4,
    addBoxPlot = TRUE
   )
dev.off()
  pdf (file.path ('Plots', 'QC_plots.pdf'))
  wrap_plots (p1, p2 ,p3, p4)
  dev.off()


  varfeat = 25000
  LSI_method = 2
  archp = addIterativeLSI (ArchRProj = archp,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = TRUE, LSIMethod = LSI_method,
    varFeatures = varfeat)

  archp = addClusters (input = archp, resolution = 3,
    reducedDims = "IterativeLSI", maxClusters = 100,
    force = TRUE)
  archp = addUMAP (ArchRProj = archp, 
    reducedDims = "IterativeLSI",
    force = TRUE)

pdf ()  
umap_p0 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters",
   embedding = "UMAP",
   #pal = palette_celltype_lv1,
   labelMeans = FALSE)

umap_p1 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype_revised",
   embedding = "UMAP",
   pal = palette_celltype_lv1,
   labelMeans = FALSE)

umap_p2 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
   embedding = "UMAP",
   pal = palette_sample,
   labelMeans = FALSE)
dev.off()

pdf (file.path ('Plots','celltype_revised_umap.pdf'))
umap_p0
umap_p1
umap_p2
dev.off()


# Plot gene score of cell type markers ####
meso_markers = c(
  'KRT19','AMOTL2','EGFR',
  'HP','BDKRB1','ZBTB7C',
  'SFTA3','SFTPB','LINC00261',
  'COL1A1','FBN1','COL6A2',
  'MYH11','COL4A2','COL4A1',
  'CLDN5','VWF','CDH5',
  'LYZ','IL1B','CD83',
  'CD3E','BCL11B','RUNX3',
  'GNLY', 'KLRD1','ITGAE',
  'CD79A','PAX5','BLK',
  'IGLL5','ELL2','ELL2',
  'VASH2','MAD1L1','ZFAT')
#meso_markers = rev (meso_markers)

if (!any (ls() == 'gsSE')) gsSE = fetch_mat (archp, mat = 'GeneScore')
gsMat = assays (gsSE)[[1]]

rownames (gsMat) = rowData (gsSE)$name
gsMat_mg = gsMat[meso_markers, ]
gsMat_mg = as.data.frame (t(gsMat_mg))
gsMat_mg$metaGroup = as.character(archp@cellColData[,metaGroupName])
gsMat_mg = aggregate (.~ metaGroup, gsMat_mg, mean)
rownames (gsMat_mg) = gsMat_mg[,1]
gsMat_mg = gsMat_mg[,-1]
gsMat_mg = gsMat_mg[names(table (archp@cellColData[,metaGroupName])[table (archp@cellColData[,metaGroupName]) > 50]),]
gsMat_mg = gsMat_mg[celltype_order,]

DAG_hm = Heatmap (t(t(scale(gsMat_mg))), 
        column_labels = colnames (gsMat_mg),
        column_title = paste('top',top_genes),
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        cluster_rows = F,
        row_names_side = 'left',
        #col = palette_genescore_fun(gsMat_mg),
        col = palette_expression,
        #col = pals_heatmap[[5]],
        cluster_columns=F,#col = pals_heatmap[[1]],
        row_names_gp = gpar(fontsize = 8),
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 6),
        rect_gp = gpar(col = "white", lwd = .1),
        border=TRUE
        #right_annotation = motif_ha
        )

  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (file.path('Plots',paste0('DAG_clusters_',metaGroupName,'_heatmaps2.pdf')), width = 6, height = 3)
print(DAG_hm)
dev.off()





# ### Gene score based analysis ####

#   # Find DAG ####
#   #metaGroupName = "Clusters"
#   metaGroupName = 'celltype_lv1'
#   force = FALSE
#   if (!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force)
#     {
#     DAG_list = getMarkerFeatures (
#       ArchRProj = archp, 
#       testMethod = "wilcoxon",
#             #useGroups = "ClusterA",
#             #bgdGroups = "Clusters1B",
#       binarize = FALSE,
#       useMatrix = "GeneScoreMatrix",
#       groupBy = metaGroupName
#     #  useSeqnames="z"
#     )

#     listnames = colnames (DAG_list)
#     DAG_list = lapply (1:ncol (DAG_list), function(x) 
#       {
#       df = DAG_list[,x]  
#       df = do.call (cbind, (assays(df)))
#       colnames(df) = names (assays(DAG_list))
#       df$gene = rowData (DAG_list)$name
#       df
#       })
#     names (DAG_list) = listnames
#     saveRDS (DAG_list, paste0 ('DAG_',metaGroupName,'.rds'))    
#     } else {
#     DAG_list = readRDS (paste0('DAG_',metaGroupName,'.rds'))
#     }
  
#   FDR_threshold = 1e-2
#   lfc_threshold = 0
#   top_genes = 3
#   DAG_top_list = DAG_list[sapply (DAG_list, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$Log2FC) > lfc_threshold,]) > 0)]
#   DAG_top_list = lapply (seq_along(DAG_top_list), function(x) {
#     res = DAG_top_list[[x]]
#     res = na.omit (res)
#     res = res[res$FDR < FDR_threshold,]
#     res = res[order (res$FDR), ]
#     res = res[abs(res$Log2FC) > lfc_threshold,]
#     res$comparison = names(DAG_top_list)[x]
#     if (nrow(res) < top_genes) 
#       {
#       res
#       } else {
#       head (res,top_genes)
#       }
#     })
#   DAG_df = Reduce (rbind ,DAG_top_list)
  
#   if (!any (ls() == 'gsSE')) gsSE = ArchR::getMatrixFromProject (archp, useMatrix = 'GeneScoreMatrix')
#   gsSE = gsSE[, archp$cellNames]
#   gsMat = assays (gsSE)[[1]]
#   rownames (gsMat) = rowData (gsSE)$name
#   gsMat_mg = gsMat[rownames (gsMat) %in% DAG_df$gene, ]
#   gsMat_mg = as.data.frame (t(gsMat_mg))
#   gsMat_mg$metaGroup = as.character(archp@cellColData[,metaGroupName])
#   gsMat_mg = aggregate (.~ metaGroup, gsMat_mg, mean)
#   rownames (gsMat_mg) = gsMat_mg[,1]
#   gsMat_mg = gsMat_mg[,-1]
#   gsMat_mg = gsMat_mg[names(table (archp@cellColData[,metaGroupName])[table (archp@cellColData[,metaGroupName]) > 50]),]
#   DAG_hm = Heatmap (t(scale(gsMat_mg)), 
#           row_labels = colnames (gsMat_mg),
#           column_title = paste('top',top_genes),
#           clustering_distance_columns = 'euclidean',
#           clustering_distance_rows = 'euclidean',
#           cluster_rows = F,
#           #col = pals_heatmap[[5]],
#           cluster_columns=F,#col = pals_heatmap[[1]],
#           row_names_gp = gpar(fontsize = 6),
#           column_names_gp = gpar(fontsize = 4),
#           rect_gp = gpar(col = "white", lwd = .5),
#           border=TRUE
#           #right_annotation = motif_ha
#           )

#   #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
# pdf (paste0('Plots/DAG_clusters_',metaGroupName,'_heatmaps2.pdf'), width = 8, height = 50)
# print(DAG_hm)
# dev.off()


# # Plot gene score of cell type markers ####
# meso_markers = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/gene_sets/highlevel_MPM_markers.csv')[[1]]
# meso_markers = c(meso_markers, 'IGLL5')
# meso_markers = meso_markers[meso_markers != 'IGHM']
# #meso_markers = c(meso_markers, 'KRT5','LILRA4','MS4A1')
# meso_markers = c()
# meso_markers = c('TOX','PCDHGA6','PCDHGC3')
# archp = addImputeWeights (archp)

# pdf()
# p <- plotEmbedding(
#     ArchRProj = archp,
#     colorBy = "GeneScoreMatrix", 
#     name = meso_markers, 
#     embedding = "UMAP",
#     pal = palette_expression,
#     imputeWeights = getImputeWeights(archp)
# )
# dev.off()
# #p = lapply (p, function(x) x + theme_void() + NoLegend ()) #+ ggtitle scale_fill_gradient2 (rev (viridis::plasma(100))))

# pdf (file.path('Plots','marker_genes_feature_plots_3.pdf'), width = 8, height = 5)
# print (wrap_plots (p, ncol = 3))
# dev.off()

# Generate CNV map for each sample ####
# Get Granges of blacklist regions 
blacklist = toGRanges (paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/blacklisted_regions/ENCODE_blacklist/',"hg38-blacklist.v2.bed")) # taken from https://github.com/Boyle-Lab/Blacklist/tree/master/lists
projdir_cnv = file.path('..','..','per_sample_QC_signac','CNV_analysis')
dir.create (projdir_cnv)

# Get GRanges of bins excluding black list regions
ws = 10e5
ws = 10e6
ss = 10e5
ss = 5e6
if (!file.exists (file.path (projdir_cnv, paste0 ('windows_',ws,'_',ss,'.rds'))))
  {
  windows = makeWindows (genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist,
    windowSize = ws, slidingSize = ss)
  saveRDS (windows, file.path (projdir_cnv, paste0 ('windows_',ws,'_',ss,'.rds')))
  } else {
  windows = readRDS (file.path (projdir_cnv, paste0 ('windows_',ws,'_',ss,'.rds')))
  }

deleted_chr = c('chr4','chr22','chr13')
if (!exists ('fragments_l')) fragments_l = getFragmentsFromProject (archp)
print_mat = F
force=F

# Load also normal mesothelium to compare CNV against ####
archp_n = loadArchRProject ('../../tumor_compartment/scatac_ArchR')
archp_n = archp_n[archp_n$Sample2 == 'normal1']
if (!exists ('fragments_n')) fragments_n = getFragmentsFromProject (archp_n)
fragments_n = unlist (fragments_n)
fragments_l[['normal1']] = fragments_n

# Loop per sample and run CNV analysis
sample_names = c(sample_names, 'normal1')
cnaObj_l = list()
for (sam in sample_names)
  {
  # Run CNV analysis
  #force=FALSE
  if (!file.exists (file.path(projdir_cnv, paste0('CNV_LFC_GC_',sam,'_ws_',ws,'_ss_',ss,'.rds'))) | force)
    {
    cnaObj_l[[sam]] = scCNA (windows, fragments_l[[sam]], neighbors = 100, LFC = 1.5, FDR = 0.1, force = FALSE, remove = c("chrX","chrM","chrY"))
    saveRDS(cnaObj_l[[sam]], file.path(projdir_cnv, paste0('CNV_LFC_GC_',sam,'_ws_',ws,'_ss_',ss,'.rds')))
    } else {
    message ('cnaObj found! loading...')
    cnaObj_l[[sam]] = readRDS (file.path(projdir_cnv, paste0('CNV_LFC_GC_',sam,'_ws_',ws,'_ss_',ss,'.rds')))
    }
  }

  #rownames(meta.data_df2) = sapply (rownames(meta.data_df2), function(x) unlist(strsplit (x,'\\.'))[2])
  # sgn_l[[sam]]@meta.data = sgn_l[[sam]]@meta.data[,!colnames(sgn_l[[sam]]@meta.data) %in% colnames(meta.data)]
  # sgn_l[[sam]]@meta.data = cbind (sgn_l[[sam]]@meta.data, meta.data[match(colnames(sgn_l[[sam]]), rownames(meta.data)),])

  # Print CNV heatmap
  colnames(cnaObj_l[['normal1']]) = sapply (colnames(cnaObj_l[['normal1']]), function(x) unlist(strsplit(x,'\\#'))[2])
  colnames(cnaObj_l[['normal1']]) = paste0('normal1#',colnames(cnaObj_l[['normal1']]))
  malignant_cells = c(rownames(archp@cellColData)[archp$celltype_lv1 == 'Malignant'], colnames(cnaObj_l[['normal1']]))
  for (sam in sample_names)
    {
    mat_type = 'z'
    cnaObj_mat = t(cnaObj_l[[sam]]@assays@data@listData[[mat_type]])
    colnames (cnaObj_mat) = paste0(seqnames(rowRanges(cnaObj_l[[sam]])),':', ranges(rowRanges(cnaObj_l[[sam]])))
    rownames (cnaObj_mat) = paste0(sam,'#',colnames (cnaObj_l[[sam]]))
    cnaObj_mat = cnaObj_mat[rownames(cnaObj_mat) %in% malignant_cells,]
    cnaObj_mat[is.na(cnaObj_mat)] = 0
    cnaObj_mat[is.infinite(cnaObj_mat)]= 0
    # ha = HeatmapAnnotation (bar=row_ann,
    #   which='row')
    # #col_fun = colorRamp2(c(min(cnaObj_z), 2,max(cnaObj_z)), c("blue", "white", "red"))  
    row_chr = gsub ('\\:.*','',colnames(cnaObj_mat))
    row_chr = factor (row_chr, levels = unique (row_chr))
    png (file.path('Plots', paste0('GL_method_',sam,'_',mat_type,'_gdn_heatmap.png')), width=1500,height=800)
    print (Heatmap (cnaObj_mat, 
    cluster_rows=T,#col=col_fun,
    cluster_columns=F, 
     #left_annotation = ha, 
     #row_km = 2,
    row_names_gp = gpar(fontsize = 0.01),column_names_gp = gpar(fontsize = 0.1),
    column_split = row_chr,
    column_gap = unit(2, "mm"),
    cluster_column_slices = F,row_title_gp = gpar(fontsize = 1),column_title_rot=90))
    dev.off()
    }

cnaMat_obj_l = list()
for (sam in sample_names)
    {
    mat_type = 'z'
    cnaObj_mat = t(cnaObj_l[[sam]]@assays@data@listData[[mat_type]])
    colnames (cnaObj_mat) = paste0(seqnames(rowRanges(cnaObj_l[[sam]])),':', ranges(rowRanges(cnaObj_l[[sam]])))
    bc = colnames(cnaObj_l[[sam]])
    if (any (grep ('\\#', colnames(cnaObj_l[[sam]])))) bc = sapply (colnames(cnaObj_l[[sam]]), function(x) unlist(strsplit(x, '\\#'))[2])
    rownames (cnaObj_mat) = paste0(sam,'#',bc)
    cnaObj_mat = cnaObj_mat[rownames(cnaObj_mat) %in% malignant_cells,]
    cnaObj_mat[is.na(cnaObj_mat)] = 0
    cnaObj_mat[is.infinite(cnaObj_mat)]= 0
    cnaMat_obj_l[[sam]] = cnaObj_mat
    }
cnaMat_obj_avg = lapply (cnaMat_obj_l, function(x) colMeans (x))
cnaMat_obj_avg = do.call (cbind, cnaMat_obj_avg)
cnaMat_obj_avg = cnaMat_obj_avg[,sample_names]
ha = rowAnnotation (sample = colnames(cnaMat_obj_avg), col = list(sample = palette_sample))
pdf (file.path('Plots', paste0('GL_method_',mat_type,'_CNV_heatmap2.pdf')), width=6,height=3)
    print (Heatmap (t(cnaMat_obj_avg), 
    cluster_rows=F,#col=col_fun,
    cluster_columns=F, 
    left_annotation = ha, 
    row_names_side = 'left',
    column_title_side = 'bottom',
    border=T,
     #row_km = 2,
    row_names_gp = gpar (fontsize = 8),
    column_names_gp = gpar(fontsize = 0),
    column_title_gp = gpar(fontsize = 7),
    column_split = row_chr,
    row_split = ifelse ('normal1' == colnames(cnaMat_obj_avg),'normal','tumor'),
    column_gap = unit(.5, "mm"),
    cluster_column_slices = F,
    row_title_gp = gpar(fontsize = 1),
    column_title_rot=45))
dev.off()




cnv_cols = c('cnvload_z','cnvload_lg','cnvload_counts','cCNV_score')
cnv_cols = c('chr22','chr13','chr4')
fp = list()
for (sam in sample_names)
  {
  fp[[sam]] = FeaturePlot (sgn_l[[sam]], feature = cnv_cols, combine=FALSE)
  for (ft in 1:length(fp[[sam]])) {fp[[sam]][[ft]] = fp[[sam]][[ft]] + scale_colour_gradientn (colours = viridis::turbo(100),na.value="white")}
  }
pdf (file.path(projdir,'Plots','CNV_score_per_sample_umap.pdf'),width=10,height=5)
lapply (sample_names, function(x) print (wrap_plots (fp[[x]])))
dev.off()



### Run peak calling on celltype annotation ####
# Add tumor sample info in celltype metagroup
archp$celltype_revised_sample = archp$celltype_revised
archp$celltype_revised_sample[archp$celltype_revised_sample == 'Malignant'] = paste0('Malignant_',archp$Sample[archp$celltype_revised_sample == 'Malignant'])

### Run peak calling ####
metaGroupName = "Clusters"
force=TRUE
peak_reproducibility=2
if(!all(file.exists(file.path('PeakCalls', paste0(unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds')))) | force) source ('../../git_repo/utils/callPeaks.R')
  

# Make barplots of number of peaks per cell type ####
ct_peaks = lapply (file.path('PeakCalls',paste0(names(palette_celltype_lv1)[names(palette_celltype_lv1) %in% unique(archp$celltype_lv1)],'-reproduciblePeaks.gr.rds')), function(x) readRDS (x))

ct_peaks = do.call (cbind,  lapply(ct_peaks, function(x) table (x$peakType)))
colnames (ct_peaks) = names(palette_celltype_lv1)[names(palette_celltype_lv1) %in% unique(archp$celltype_lv1)]

ct_peaks = as.data.frame (ct_peaks)
ct_peaks$peaktype = rownames(ct_peaks)
ct_peaks = gather (ct_peaks, celltype, peaks, 1:(ncol(ct_peaks)-1))

df = ct_peaks %>% 
group_by (celltype) %>% 
summarize (sum = sum(peaks)) %>% 
arrange(sum)

ct_peaks$celltype = factor (ct_peaks$celltype, levels = rev(df$celltype[order(df$sum)]))

bp = ggplot (ct_peaks, aes(x = celltype, y = peaks, fill = peaktype)) + 
geom_bar (stat= 'identity') + 
paletteer::scale_fill_paletteer_d("nbapalettes::knicks_city") +
gtheme

pdf (file.path ('Plots','peakcalls.pdf'))
bp
dev.off()

archp = saveArchRProject (archp, load=TRUE)
  
  metaGroupNames = c('TSSEnrichment','nFrags','ReadsInTSS','FRIP')  
    umap_p12 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
     name = x, embedding = "UMAP"))
      
  pdf (paste0(projdir,'/Plots/qc_umap_after_filtering.pdf'), 15,15)
  wrap_plots (umap_p12, ncol=5)
  dev.off()


### chromVAR analysis ####
run_chromVAR = TRUE

if (run_chromVAR)
  {  
  archp = addBgdPeaks (archp, force= TRUE)
  archp = addMotifAnnotations (ArchRProj = archp,
      motifSet = "cisbp",
      #motifSet = 'JASPAR2020',
      #name = "JASPAR2020_Motif",
      force=TRUE)
  archp = addDeviationsMatrix (
    ArchRProj = archp, 
    peakAnnotation = "Motif",
    force = TRUE
  )
  
  archp = saveArchRProject (ArchRProj = archp,  
      load = TRUE)
  }



### ChromVAR based analysis ####


  # Find DAM ####
  metaGroupName = "celltype_lv1"
  force = TRUE
  if (!file.exists (paste0('DAM_',metaGroupName,'.rds')) | force)
    {
    DAM_list = getMarkerFeatures (
      ArchRProj = archp, 
      testMethod = "wilcoxon",
            #useGroups = "ClusterA",
            #bgdGroups = "Clusters1B",
      binarize = FALSE,
      useMatrix = "MotifMatrix",
      groupBy = metaGroupName
    #  useSeqnames="z"
    )

    listnames = colnames (DAM_list)
    DAM_list = lapply (1:ncol (DAM_list), function(x) 
      {
      df = DAM_list[,x]  
      df = do.call (cbind, (assays(df)))
      colnames(df) = names (assays(DAM_list))
      df$gene = rowData (DAM_list)$name
      df
      })
    names (DAM_list) = listnames
    saveRDS (DAM_list, paste0 ('DAM_',metaGroupName,'.rds'))    
    } else {
    DAM_list = readRDS (paste0('DAM_',metaGroupName,'.rds'))
    DAM_list = lapply (DAM_list, function(x) {x$gene = gsub ('_.*','',x$gene); x})
    DAM_list = lapply (DAM_list, function(x) {x$gene = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", x$gene); x})    
    }
  DAM_list = DAM_list[celltype_order]
  # Filter by active genes as computed by correlation with gene score
  # active_genes = corGSM_MM$MotifMatrix_name[corGSM_MM$cor > -Inf]
  # DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% active_genes,])

  # Filter by genes minimally expressed from scRNA (in at least 1 celltype) ####
  ps = log2(as.data.frame (AverageExpression (srt, 
    features = sapply (unique(unlist(lapply(DAM_list, function(x) x$gene))), function(x) unlist(strsplit (x, '_'))[1]), 
    group.by = 'celltype_lv1')[[1]]) +1)
  colnames (ps) = gsub ('-','_',colnames(ps))
  min_exp = 0.1
  DAM_list2 = lapply (names(DAM_list), function(x) DAM_list[[x]][DAM_list[[x]]$gene %in% rownames(ps)[ps[,x] > min_exp],])
  #active_TFs = rownames(ps)[rowSums(ps) > 0]

  #DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% active_TFs,])

  names (DAM_list2) = names (DAM_list)
  FDR_threshold = 1e-3
  meandiff_threshold = 0
  top_genes = 5
  DAM_top_list = DAM_list2[sapply (DAM_list2, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$MeanDiff) > meandiff_threshold,]) > 0)]
  DAM_top_list = lapply (seq_along(DAM_top_list), function(x) {
    res = DAM_top_list[[x]]
    #res = na.omit (res)
    res = res[res$FDR < FDR_threshold,]
    res = res[order (res$FDR), ]
    res = res[res$MeanDiff > meandiff_threshold,]
    res$comparison = names(DAM_top_list)[x]
    if (nrow(res) < top_genes) 
      {
      res
      } else {
      head (res,top_genes)
      }
    })
  DAM_df = Reduce (rbind ,DAM_top_list)
  

  # # Get deviation matrix ####
  force = F
  if (!exists ('mSE') | force) mSE = fetch_mat (archp, 'Motif')
  mMat = assays (mSE)[[1]]
  rownames (mMat) = rowData (mSE)$name
  mMat_mg = mMat[DAM_df$gene, ]
  mMat_mg = as.data.frame (t(mMat_mg))
  mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
  mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
  rownames (mMat_mg) = mMat_mg[,1]
  mMat_mg = mMat_mg[,-1]
  mMat_mg = mMat_mg[names (DAM_list),]
  mMat_mg = na.omit (mMat_mg)
  #mMat_mg = mMat_mg[names(table (archp@cellColData[,metaGroupName])[table (archp@cellColData[,metaGroupName]) > 50]),]
  DAM_hm = Heatmap (scale(mMat_mg), 
          column_labels = colnames (mMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 7),
          column_names_rot = 45,
          name = 'TF activity',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_deviation_fun(scale(mMat_mg))
          #right_annotation = motif_ha
          )

#DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_heatmaps2.pdf')), width = 8, height = 3)
print(DAM_hm)
dev.off()
# ### Run TF correlation to identify TF modules across all cells #### 
# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat = as.matrix(mMat)#[selected_TF,])
# mMat = mMat[names(km$cluster)[km$cluster ==3],]

# mMat_cor = cor (as.matrix(t(scale(mMat))), method = 'spearman')

# set.seed(1234)
# centers=4
# km = kmeans (mMat_cor, centers=centers)

# pdf (file.path ('Plots','TME_TF_modules_heatmap.pdf'), width = 4,height=3)
# cor_mMat_hm = draw (Heatmap (mMat_cor,
#   clustering_distance_rows='euclidean' ,
#   clustering_distance_columns = 'euclidean', 
#   col=palette_deviation_cor_fun, 
#   row_split = km$cluster,
#   column_split = km$cluster,
#   border=T,
#   row_names_gp = gpar(fontsize = 0), 
#   column_names_gp = gpar(fontsize = 0)
#   ))
# dev.off()

# pdf (file.path ('Plots','TF_modules_heatmap2.pdf'), width = 4,height=3)
# cor_mMat_hm
# dev.off()

# tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
# names (tf_modules) = paste0('mod_',unique(km$cluster))
# tf_modules = do.call (cbind, tf_modules)
# archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% paste0('mod_',unique(km$cluster))]
# archp@cellColData = cbind (archp@cellColData, tf_modules) 

# pdf()
# TF_p = plotEmbedding (
#     ArchRProj = archp,
#     colorBy = "cellColData",
#     name = paste0('mod_',unique(km$cluster)), 
#     pal = rev(palette_deviation_correlation),
#     #useSeqnames='z',
#     embedding = "UMAP")
# dev.off()
# pdf (file.path ('Plots','TF_modules_umap2.pdf'), width = 20,height=6)
# wrap_plots (TF_p, ncol=5)
# dev.off()

# # Find shared TF across celltypes ####
# tf_name2 = unlist(sapply (c('TGIF1','TWIST2','NFKB2','HMGA2'), function(x) rownames(assay(mSE))[grepl (x, rownames(assay(mSE)))]))
# tf_name2 = paste0('z:',tf_name2)
# archp = addImputeWeights (archp)
# pdf()
# TF_p = plotEmbedding (
#     ArchRProj = archp,
#     colorBy = "MotifMatrix",
#     name = tf_name2, 
#     useSeqnames='z',
#     pal = rev (palette_deviation),    
#     embedding = "UMAP",
#     imputeWeights = getImputeWeights(archp)
#     )
# dev.off()
# pdf (file.path ('Plots','pan_TF_fplots.pdf'), width = 30,height=16)
# wrap_plots (TF_p, ncol=5)
# dev.off()


# tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
# names (tf_modules) = paste0('mod_',unique(km$cluster))
# tf_modules = do.call (cbind, tf_modules)
# archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% paste0('mod_',unique(km$cluster))]
# archp@cellColData = cbind (archp@cellColData, tf_modules)

# archp = addImputeWeights (archp)
# TF_p = lapply (paste0('mod_',unique(km$cluster)), function(x) plotEmbedding (
#     ArchRProj = archp,
#     colorBy = "cellColData",
#     name = x, 
#     pal = palette_deviation,
#     #useSeqnames='z',
#     imputeWeights = getImputeWeights(archp),
#     embedding = "UMAP"))

# pdf (file.path ('Plots','TF_modules_umap.pdf'), width = 30,height=14)
# wrap_plots (TF_p, ncol=5)
# dev.off()


### Subset ArchR project ####
run_dropcells = FALSE
if (run_dropcells) archp = saveArchRProject (archp, dropCells = T) # Make sure to run this first before subsetting with the custom subset function

# Subset T cells ####
metaGroupName = 'celltype_revised'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('T_cells','NK')],
  outputDirectory = file.path('..','..','NKT_cells','scatac_ArchR'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

# Subset Myeloid ####
metaGroupName = 'celltype_revised'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('Myeloid')],
  outputDirectory = file.path('..','..','myeloid_cells','scatac_ArchR'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

# Subset Stroma ####
metaGroupName = 'celltype_lv1'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('Endothelial','Fibroblasts','Mesothelium','SmoothMuscle')],
  outputDirectory = file.path('..','..','stroma','scatac_ArchR'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

# Subset Malignant ####
metaGroupName = 'celltype_lv1'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('Malignant')],
  outputDirectory = file.path('..','..','tumor_compartment','scatac_ArchR'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

metaGroupName='celltype_lv1'
subsetArchRProject_light (ArchRProject = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('Malignant')],
  projdir_new = file.path('..','..','tumor_compartment','scatac_ArchR')
  )

# Subset B cells ####
metaGroupName='celltype_lv1'
subsetArchRProject_light (ArchRProject = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('B_cells','Plasma')],
  projdir_new = file.path('..','..','B_cells','scatac_ArchR')
  )


#### Annotate barcode with level2  annotations in scATAC-seq ####
tnk_ann = read.csv ('../../NKT_cells/scatac_ArchR/barcode_annotation.csv')
mye_ann = read.csv ('../../myeloid_cells/scatac_ArchR/barcode_annotation.csv')
stro_ann = read.csv ('../../stroma/scatac_ArchR/barcode_annotation.csv')
bcells_ann = read.csv ('../../B_cells/scatac_ArchR/barcode_annotation.csv')
celltype_lv2_df = rbind (tnk_ann, mye_ann, stro_ann, bcells_ann)
archp$celltype_lv2 = celltype_lv2_df$celltype_lv2[match(rownames(archp@cellColData), celltype_lv2_df$barcode)]
archp$celltype_lv2[archp$celltype_lv1 == 'Malignant'] = 'Malignant'
archp$celltype_lv2[archp$celltype_lv1 == 'pDCs'] = 'pDCs'
archp$celltype_lv2[archp$celltype_lv1 == 'Alveolar'] = 'Alveolar'

# Save annotation in repo
annotation_df = data.frame (barcode = rownames(archp@cellColData), celltype_lv1 = archp$celltype_lv1, celltype_lv2 = archp$celltype_lv2)
write.csv (annotation_df, '../../git_repo/barcode_anntation.csv')



#### Annotate barcode with level2 annotations in scRNA-seq ####
srt_tnk = readRDS (file.path('..','..','NKT_cells','scrna','srt.rds'))
tnk_ann = srt_tnk$celltype_lv2
tnk_ann = tnk_ann[names(tnk_ann) %in% colnames(srt)]
srt$celltype_lv2 = srt$celltype_lv1
srt$celltype_lv2[match (names(tnk_ann), colnames(srt))] = tnk_ann

mye_ann = read.csv ('../../myeloid_cells/scatac_ArchR/barcode_annotation.csv')
stro_ann = read.csv ('../../stroma/scatac_ArchR/barcode_annotation.csv')
bcells_ann = read.csv ('../../B_cells/scatac_ArchR/barcode_annotation.csv')
celltype_lv2_df = rbind (tnk_ann, mye_ann, stro_ann, bcells_ann)
archp$celltype_lv2 = celltype_lv2_df$celltype_lv2[match(rownames(archp@cellColData), celltype_lv2_df$barcode)]
archp$celltype_lv2[archp$celltype_lv1 == 'Malignant'] = 'Malignant'
archp$celltype_lv2[archp$celltype_lv1 == 'pDCs'] = 'pDCs'
archp$celltype_lv2[archp$celltype_lv1 == 'Alveolar'] = 'Alveolar'

# Save annotation in repo
annotation_df = data.frame (barcode = rownames(archp@cellColData), celltype_lv1 = archp$celltype_lv1, celltype_lv2 = archp$celltype_lv2)
write.csv (annotation_df, '../../git_repo/barcode_anntation.csv')


# Show that enhancer linked to NR4A2 is only up in NK KLRC1 and CD8 exhausted across all cells ####
#nkt_ann = read.csv (file.path('..','..','NKT_cells','scatac_ArchR','barcode_annotation_nkt.csv'))
#mye_ann = read.csv (file.path('..','..','myeloid_cells','scatac_ArchR','barcode_annotation.csv'))
# archp$celltype2 = archp$celltype_revised
# archp$celltype2[match(nkt_ann$barcode, rownames(archp@cellColData))] = nkt_ann$celltype
# archp$celltype2[match(mye_ann$barcode, rownames(archp@cellColData))] = mye_ann$celltype
metaGroupName = 'celltype_lv2'
pMats = getGroupSE(
  ArchRProj = archp,
  useMatrix = 'PeakMatrix',
  groupBy = metaGroupName,
  divideN = TRUE,
  scaleTo = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)
enhancer_region = GRanges ('chr2:156480366-156480866')
promoter_region = getPeakSet(archp)
promoter_region =  promoter_region[which(promoter_region$peakType == 'Promoter' & promoter_region$nearestGene == 'NR4A2')]
promoter_region = GRanges (seqnames = as.character(seqnames(promoter_region))[1], ranges= IRanges(start = min(start(promoter_region)), end = max(end(promoter_region))))

pmat_enhancer = unlist(as.data.frame (assay (pMats[queryHits(findOverlaps (GRanges(rowData(pMats)), enhancer_region)),])))
pmat_enhancer_df = data.frame (region = pmat_enhancer, celltype =names(pmat_enhancer))
pmat_enhancer_df$celltype = factor (pmat_enhancer_df$celltype, pmat_enhancer_df$celltype[order(-pmat_enhancer_df$region)])
pmat_enhancer_df$type = 'enhancer'
pmat_promoter = unlist(colSums (as.data.frame (assay (pMats[queryHits(findOverlaps (GRanges(rowData(pMats)), promoter_region)),]))))
pmat_promoter_df = data.frame (region = pmat_promoter, celltype =names(pmat_promoter))

pmat_promoter_df$type = 'promoter'
pmat_df = rbind (pmat_enhancer_df, pmat_promoter_df)
pmat_df$type = factor (pmat_df$type, levels =c ('promoter','enhancer'))

pdf (file.path ('Plots','eNR4F2_accessibility_barplot.pdf'), height=3, width=6)
ep = ggplot (pmat_df, aes (x = celltype, y = region, fill = type)) + 
geom_bar(stat = 'identity',position ='stack', color='grey22') + gtheme +
scale_fill_manual (values = c(enhancer = '#C1D32FFF', promoter = 'grey'))
#+ scale_fill_manual (values = c(palette_tnk_cells, palette_myeloid, palette_celltype_lv1))
ep
dev.off()


# Check expression of NR4A2 across cell types ####
ps = log2(as.data.frame (AverageExpression (srt, 
    features = 'NR4A2',
    group.by = 'celltype_lv1')[[1]]) +1)
ps = as.data.frame (t(ps))
ee = ggplot (ps, aes (x = rownames(ps), y = V1)) + 
geom_bar(stat = 'identity',position ='stack', color='grey22') + gtheme +
scale_fill_manual (values = c(enhancer = '#C1D32FFF', promoter = 'grey'))

pdf (file.path ('Plots','eNR4F2_expression_barplot.pdf'), height=3, width=6)
ee
dev.off()

pdf(file.path ('Plots','NR4A2_dotplot.pdf'), width = 9)
VlnPlot (srt, feature = 'NR4A2', group.by = 'celltype3', cols = c(palette_tnk_cells, palette_celltype_lv1))
VlnPlot (srt[,!srt$celltype_lv1 %in% c('NK','T_cells')], feature = 'NR4A2', group.by = 'celltype_lv2', cols = c(palette_tnk_cells, palette_celltype_lv1))
dev.off()

### Show as UMAP #####

# # first need to add features to peakMatrix ####
# ps <- archp@peakSet
# ps$name <- paste0(seqnames(ps),"_peak",ps$idx)
# archp@peakSet <- ps
# archp =  addPeakMatrix (ArchRProj = archp, force = TRUE, threads = 16)

# ps[queryHits (findOverlaps (ps,enhancer_region))]
# ps[queryHits (findOverlaps (ps,promoter_region))]

# archp = addImputeWeights (archp)

# pdf ()
# e_p = plotEmbedding(
#   ArchRProj = archp,
#   colorBy = "PeakMatrix",
#   name = c("chr2_peak20773"),
#   embedding = "UMAP",
#   quantCut = c(0.01, 0.95),
#   imputeWeights = getImputeWeights(ArchRProj = archp),
#   rastr = TRUE,
#   plotAs = "points"
# )

# p_p = plotEmbedding(
#   ArchRProj = archp,
#   colorBy = "PeakMatrix",
#   name = c("chr2_peak20739"),
#   embedding = "UMAP",
#   quantCut = c(0.01, 0.95),
#   imputeWeights = getImputeWeights(ArchRProj = archp),
#   rastr = TRUE,
#   plotAs = "points"
# )

# ct = plotEmbedding(
#   ArchRProj = archp,
#   colorBy = "cellColData",
#   name = c("celltype_lv2"),
#   embedding = "UMAP",
# #  quantCut = c(0.01, 0.95),
# #  imputeWeights = getImputeWeights(ArchRProj = archp),
#   rastr = TRUE,
#   plotAs = "points"
# )
# dev.off()

# pdf(file.path ('Plots','NR4A2_enhancer_accessibility_umap.pdf'), width = 12)
# wrap_plots (e_p, p_p, ct)
# dev.off()

# Show regions linked to IL1B across all cells ####
nkt_ann = read.csv (file.path('..','..','NKT_cells','scatac_ArchR','barcode_annotation_nkt.csv'))
mye_ann = read.csv (file.path('..','..','myeloid_cells','scatac_ArchR','barcode_annotation_cnmf_celltypes.csv'))
archp$celltype2 = archp$celltype_revised
archp$celltype2[match(nkt_ann$barcode, rownames(archp@cellColData))] = nkt_ann$celltype
archp$celltype2[match(mye_ann$barcode, rownames(archp@cellColData))] = mye_ann$celltype
metaGroupName = 'celltype2'  
pMats = getGroupSE(
  ArchRProj = archp,
  useMatrix = 'PeakMatrix',
  groupBy = metaGroupName,
  divideN = TRUE,
  scaleTo = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)
hub_region = GRanges ('chr2:112774866-112842997')
hub_region = GRanges (c(
  'chr2:112828807-112833348',
  'chr2:112835783-112838053',
  'chr2:112838685-112840955'))
hub_region_names = c('momac_specific','unspecific_(promoter)','momac_IM')
pset = getPeakSet (archp)
hub_region_peaks = lapply (seq_along(hub_region), function(x) pset[queryHits (findOverlaps(pset,hub_region[x]))])

pmat_enhancer = do.call(cbind, lapply (hub_region_peaks, function(x) as.data.frame (colMeans(assay (pMats[queryHits(findOverlaps (GRanges(rowData(pMats)), x)),])))))
colnames (pmat_enhancer) = hub_region_names
#pmat_enhancer_df$celltype = factor (pmat_enhancer_df$celltype, pmat_enhancer_df$celltype[order(-pmat_enhancer_df$region)])
#pmat_enhancer_df$type = 'enhancer'
#pmat_promoter = unlist(colSums (as.data.frame (assay (pMats[queryHits(findOverlaps (GRanges(rowData(pMats)), promoter_region)),]))))
#pmat_promoter_df = data.frame (region = pmat_promoter, celltype =names(pmat_promoter))
# panno = rowData (pMats[queryHits(findOverlaps (GRanges(rowData(pMats)), hub_region_peaks)),])
# panno = makeGRangesFromDataFrame (panno)
# panno = pset[queryHits (findOverlaps(pset, panno))]

# panno = HeatmapAnnotation (ptype = panno$peakType, which='row')
pmat_enhancer_scaled = scale(pmat_enhancer)
#pmat_enhancer_scaled[pmat_enhancer_scaled < -1] = -1
#pmat_enhancer_scaled[pmat_enhancer_scaled > 1] = 1
pmat_enhancer_scaled = pmat_enhancer_scaled[,c('momac_specific','momac_IM','unspecific_(promoter)')]
celltype_order = c('Mono','TREM2','SPP1','IFN_CXCLs','IM','Mesothelium','Malignant','Fibroblasts','B_cells','pDCs','Plasma','SmoothMuscle','Endothelial','Tregs','NK_KLRC1','Alveolar','CD4','CD8','NK_FGFBP2','CD8_exhausted')
pmat_enhancer_scaled = pmat_enhancer_scaled[celltype_order,]
pdf (file.path ('Plots','hub_peaks_heatmap.pdf'), height=2)
Heatmap (
  t(scale(pmat_enhancer_scaled)), 
  cluster_rows=F, 
  cluster_columns=F, 
  col = palette_fragments_fun, 
  column_names_rot=45)
dev.off()



pmat_promoter_df$type = 'promoter'
pmat_df = rbind (pmat_enhancer_df, pmat_promoter_df)
pmat_df$type = factor (pmat_df$type, levels =c ('promoter','enhancer'))

pdf (file.path ('Plots','eNR4F2_accessibility_barplot.pdf'), height=3, width=6)
ep = ggplot (pmat_df, aes (x = celltype, y = region, fill = type)) + 
geom_bar(stat = 'identity',position ='stack', color='grey22') + gtheme +
scale_fill_manual (values = c(enhancer = '#C1D32FFF', promoter = 'grey'))
#+ scale_fill_manual (values = c(palette_tnk_cells, palette_myeloid, palette_celltype_lv1))
ep
dev.off()




# # Input exhausted peaks as chromvar score ####
# T_ext = list(
#   T_ext = readRDS (file.path('..','..','NKT_cells','scatac_ArchR','T_cell_exhaustion_peaks.rds')),
#   T_ext2 = readRDS (file.path('..','..','NKT_cells','scatac_ArchR','T_cell_exhaustion_peaks.rds')))

# # Import chrombpnet output and cross-reference with scRNA-seq TF expression
# archp = addBgdPeaks (archp, force= TRUE)
# archp = addPeakAnnotations (ArchRProj = archp, 
#      regions = T_ext, name = "T_exhausted",force=T)

# archp = addDeviationsMatrix (
#   ArchRProj = archp, 
#   peakAnnotation = "T_exhausted",
#   force = TRUE
# )

# pdf()
# if (!any (ls() == 'tSE')) tSE = fetch_mat (archp, 'T_exhausted')
# TF_p = plotEmbedding (
#     ArchRProj = archp,
#     colorBy = "T_exhaustedMatrix",
#     name = 'z:T_ext', 
#     #useSeqnames='z',
#     pal = rev (palette_deviation),    
#     embedding = "UMAP",
#     imputeWeights = getImputeWeights(archp)
#     )
# dev.off()
# pdf (file.path ('Plots','T_exhausted_fplots.pdf'), width = 30,height=16)
# TF_p
# dev.off()



### Show cell type proportions per sample ####
archp_meta = as.data.frame (archp@cellColData)
archp_meta$celltype_lv1 = factor (archp_meta$celltype_lv1, levels = names(palette_celltype_lv1))
archp_meta$Sample2 = archp_meta$Sample
archp_meta$Sample2 = factor (archp_meta$Sample2, levels= rev(c('P1','P3','P4','P5','P8','P10','P11','P12','P13','P14','P23')))
bp = cellComp(
  seurat_obj = archp_meta,
  metaGroups = c('Sample2','celltype_lv1'),
  plot_as = 'bar',
  pal = palette_celltype_lv1
  ) + gtheme + coord_flip()

pdf (file.path('Plots',paste0('celltype_barplots.pdf')), height=6, width=5)
bp
dev.off()






# Import chromBPnet finemo motifs ####
#archp_P1 = archp[archp$Clusters %in% c('C2') & archp$Sample == 'P1']
library ('universalmotif')

ps = getPeakSet (archp)

chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet'
chromBPdir = '/sc/arion/scratch/giottb01/chromBPnet'

chrombpnet_counts = list()
celltypes = unique (archp$celltype_lv1)
celltypes = c('Malignant','Mesothelium','Alveolar','Fibroblasts','SmoothMuscle','Endothelial','Myeloid','T_cells','NK','B_cells','Plasma','pDCs') 
#celltypes = celltypes[celltypes != 'pDCs']
annotated_motifs = read.table (file.path(chromBPdir,'modisco_merged_counts','compiled','modisco_compiled.tsv'), sep='\t', header=T)
#annotated_motifs$pattern2 = sapply (annotated_motifs$pattern, function(x) unlist(strsplit(x,'__'))[2])
#rownames (annotated_motifs) = annotated_motifs$pattern2
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_counts[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
  chrombpnet_counts[[celltype]]$V8 = sub ('pos_patterns.','',chrombpnet_counts[[celltype]]$V6)
  chrombpnet_counts[[celltype]]$V8 = sub ('neg_patterns.','',chrombpnet_counts[[celltype]]$V8)
  chrombpnet_counts[[celltype]]$V9 = paste(annotated_motifs$match0,annotated_motifs$match1,annotated_motifs$match2)[match(chrombpnet_counts[[celltype]]$V8, annotated_motifs$pattern)]
  chrombpnet_counts[[celltype]]$V9 <- gsub('_HUMAN\\.H11MO\\.(0|1|2)\\.[A-D]', '', chrombpnet_counts[[celltype]]$V9)

  gr = makeGRangesFromDataFrame (chrombpnet_counts[[celltype]], keep.extra.columns=T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  chrombpnet_counts[[celltype]]$peak_type = ps$peakType[findOverlaps(gr, ps, select='first')]
  chrombpnet_counts[[celltype]]$direction = ifelse (grepl ('pos_pattern', chrombpnet_counts[[celltype]]$V6), 'positive','negative')
  #chrombpnet_counts[[celltype]] = chrombpnet_counts[[celltype]][chrombpnet_counts[[celltype]]$V4 != 'NaN_NaN_NaN',] # Some seqlets have NA motif match and qvalues...removing those
  #nonsig_motifs = chrombpnet_counts[[celltype]]$V7 > 0.05
  #chrombpnet_counts[[celltype]]$V4[nonsig_motifs] = chrombpnet_counts[[celltype]]$V6[nonsig_motifs]
  }

# Export finemo hits as tsv files 
lapply (names(chrombpnet_counts), function(x) write.table (chrombpnet_counts[[x]][,c(1:3,9)],file.path ('chromBPnet', x,'no_bias_model',paste0(x, '_finemo_counts_to_genome_browser_harmonized.tsv')), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE))

annotated_motifs_profile = read.table (file.path(chromBPdir,'modisco_merged_profile','compiled','modisco_compiled.tsv'), sep='\t', header=T)

chrombpnet_profile = list()
#celltypes = c('Mesothelium','Malignant')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  chrombpnet_profile[[celltype]]$V8 = sub ('pos_patterns.','',chrombpnet_profile[[celltype]]$V6)
  chrombpnet_profile[[celltype]]$V8 = sub ('neg_patterns.','',chrombpnet_profile[[celltype]]$V8)
  chrombpnet_profile[[celltype]]$V9 = paste(annotated_motifs_profile$match0,annotated_motifs_profile$match1, annotated_motifs_profile$match2)[match(chrombpnet_profile[[celltype]]$V8, annotated_motifs_profile$pattern)]
  chrombpnet_profile[[celltype]]$V9 <- gsub('_HUMAN\\.H11MO\\.(0|1|2)\\.[A-D]', '', chrombpnet_profile[[celltype]]$V9)
  gr = makeGRangesFromDataFrame (chrombpnet_profile[[celltype]], keep.extra.columns=T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  chrombpnet_profile[[celltype]]$peak_type = ps$peakType[findOverlaps(gr, ps, select='first')]
  chrombpnet_profile[[celltype]]$direction = ifelse (grepl ('pos_pattern', chrombpnet_profile[[celltype]]$V6), 'positive','negative')
  }

lapply (names(chrombpnet_profile), function(x) write.table (chrombpnet_profile[[x]][,c(1:3,9)],file.path (chromBPdir, x,'no_bias_model',paste0(x, '_finemo_profile_to_genome_browser_harmonized.tsv')), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE))

### Check active motifs in selected peaks from hubs ####
hubs_obj = readRDS (file.path('hubs_obj_cor_0.3_md_12500_dgs_0_min_peaks_5','global_hubs_obj.rds'))  
selected_hubs = names(hubs_obj$hubsCollapsed[grep ('WT1', hubs_obj$hubsCollapsed$gene)])
selected_peaks = do.call (rbind, hubs_obj$hubsClusters[[1]][as.numeric(selected_hubs)])
selected_peaks = makeGRangesFromDataFrame (selected_peaks)
chrombpnet_counts_selected = lapply (chrombpnet_counts, function(x) {
  gr = makeGRangesFromDataFrame (x, seqnames.field= 'V1', start.field = 'V2', end.field = 'V3')
  x[queryHits(findOverlaps(gr, selected_peaks)),]
})

chrombpnet_profile_selected = lapply (chrombpnet_profile, function(x) {
  gr = makeGRangesFromDataFrame (x, seqnames.field= 'V1', start.field = 'V2', end.field = 'V3')
  x[queryHits(findOverlaps(gr, selected_peaks)),]
})


selected_celltypes = c('Mesothelium','Fibroblasts')
chrombpnet_counts_selected = chrombpnet_counts_selected[selected_celltypes]
top_n <- 5

bp_list <- lapply(names (chrombpnet_counts_selected), function(i) {
  tbl <- table(chrombpnet_counts_selected [[i]]$V9)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_counts_selected[[i]]$direction[chrombpnet_counts_selected [[i]]$V9 == tf][1]
  })
  
  data.frame(
    Freq = proportions(top_tbl),
    TF   = tf_names,
    direction = directions,
    type = rep(i, length(top_tbl)),
    pattern = chrombpnet_counts_selected[[i]]$V6[match(tf_names, chrombpnet_counts_selected[[i]]$V9)]
  )
})

bp_df <- do.call(rbind, bp_list)

# Make neg values negative
# bp_df <- bp_df %>%
#   mutate(Freq = ifelse(direction == "negative", -Freq.Freq, Freq.Freq))

# Create custom ordering per type
bp_df <- bp_df %>%
  group_by(type, direction) %>%
  mutate(
    TF_order = ifelse(direction == "positive",
                      rank(-Freq, ties.method = "first"),  # descending
                      rank(Freq, ties.method = "first"))   # ascending for neg (opposite)
  ) %>%
  ungroup()

# Build a combined factor: ensures pos stack from bottom up, neg from top down
bp_df <- bp_df %>%
  arrange(type, direction, TF_order)

bp_df$TF_id <- bp_df$TF
bp_df$TF_id <- factor(bp_df$TF_id, levels = unique(bp_df$TF_id))
#bp_df$type = factor (bp_df$type, levels = celltypes)
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq.Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("trekcolors::lcars_2379",8)) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

  

pdf (file.path ('Plots', 'TF_abundance_counts_selected_peaks_barplot.pdf'),6,width=5.5)
bp
dev.off()

selected_celltypes = c('Mesothelium','Fibroblasts')
chrombpnet_profile_selected = chrombpnet_profile_selected[selected_celltypes]
top_n <- 10
#n <- length(chrombpnet_profile_selected)

bp_list <- lapply(names(chrombpnet_profile_selected), function(i) {
  tbl <- table(chrombpnet_profile_selected [[i]]$V9)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_profile_selected[[i]]$direction[chrombpnet_profile_selected [[i]]$V9 == tf][1]
  })
  
  data.frame(
    Freq = proportions(top_tbl),
    TF   = tf_names,
    direction = directions,
    type = rep(i, length(top_tbl)),
    pattern = chrombpnet_profile_selected[[i]]$V6[match(tf_names, chrombpnet_profile_selected[[i]]$V9)]
  )
})

bp_df <- do.call(rbind, bp_list)

# Make neg values negative
bp_df <- bp_df %>%
  mutate(Freq = ifelse(direction == "negative", -Freq.Freq, Freq.Freq))

# Create custom ordering per type
bp_df <- bp_df %>%
  group_by(type, direction) %>%
  mutate(
    TF_order = ifelse(direction == "positive",
                      rank(-Freq, ties.method = "first"),  # descending
                      rank(Freq, ties.method = "first"))   # ascending for neg (opposite)
  ) %>%
  ungroup()

# Build a combined factor: ensures pos stack from bottom up, neg from top down
bp_df <- bp_df %>%
  arrange(type, direction, TF_order)

bp_df$TF_id <- bp_df$TF
bp_df$TF_id <- factor(bp_df$TF_id, levels = unique(bp_df$TF_id))
#bp_df$type = factor (bp_df$type, levels = celltypes)
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(unique(bp_df$TF)))) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

pdf (file.path ('Plots', 'TF_abundance_profile_selected_peaks_barplot.pdf'),height = 6, width = 5.5)
bp
dev.off()


# # Subset for peak type (promoter) ####
# chrombpnet_counts_p1 = lapply (chrombpnet_counts, function(x) {
#   x = x[x$peak_type == 'Promoter',]
#   x = x[!is.na(x$peak_type),]
#   })
# chrombpnet_profile_p1 = lapply (chrombpnet_profile, function(x) 
#   {
#   x = x[x$peak_type == 'Promoter',]
#   x = x[!is.na(x$peak_type),]
#   })

# top_n <- 5
# n <- length(chrombpnet_counts_p1 )

# bp_list <- lapply(seq_len(n), function(i) {
#   tbl <- table(chrombpnet_counts_p1 [[i]]$V4)
#   tbl_sorted <- sort(tbl, decreasing = TRUE)
#   top_tbl <- head(tbl_sorted, top_n)
  
#   tf_names <- names(top_tbl)
#   directions <- sapply(tf_names, function(tf) {
#     chrombpnet_counts_p1[[i]]$V5[chrombpnet_counts_p1 [[i]]$V4 == tf][1]
#   })
  
#   data.frame(
#     Freq = proportions(top_tbl),
#     TF   = tf_names,
#     direction = directions,
#     type = rep(celltypes[[i]], length(top_tbl)),
#     ptype = 'promoter'
#   )
# })

# bp_df <- do.call(rbind, bp_list)

# # Make neg values negative
# bp_df <- bp_df %>%
#   mutate(Freq = ifelse(direction == "neg", -Freq.Freq, Freq.Freq))

# # Create custom ordering per type
# bp_df <- bp_df %>%
#   group_by(type, direction) %>%
#   mutate(
#     TF_order = ifelse(direction == "pos",
#                       rank(-Freq, ties.method = "first"),  # descending
#                       rank(Freq, ties.method = "first"))   # ascending for neg (opposite)
#   ) %>%
#   ungroup()

# # Build a combined factor: ensures pos stack from bottom up, neg from top down
# bp_df <- bp_df %>%
#   arrange(type, direction, TF_order)

# bp_df$TF_id <- paste(bp_df$TF, bp_df$type, sep = "_")
# bp_df$TF_id <- factor(bp_df$TF_id, levels = unique(bp_df$TF_id))
# bp_df$type = factor (bp_df$type, levels = celltypes)

# bp_df_cp = bp_df
# # Plot stacked bars
# bp <- ggplot(bp_df_cp, aes(x = type, y = Freq, fill = TF_id)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF))) +
#   theme_minimal(base_size = 14) +
#   ylab("Proportion of counts") +
#   xlab("Cell type") +
#   ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

# pdf (file.path ('Plots', 'TF_abundance_counts_promoter_barplot.pdf'),6,width=12.5)
# bp
# dev.off()


# top_n <- 5
# n <- length(chrombpnet_profile_p1)

# bp_list <- lapply(seq_len(n), function(i) {
#   tbl <- table(chrombpnet_profile_p1[[i]]$V4)
#   tbl_sorted <- sort(tbl, decreasing = TRUE)
#   top_tbl <- head(tbl_sorted, top_n)
  
#   tf_names <- names(top_tbl)
#   directions <- sapply(tf_names, function(tf) {
#     chrombpnet_profile_p1[[i]]$V5[chrombpnet_profile_p1[[i]]$V4 == tf][1]
#   })
  
#   data.frame(
#     Freq = proportions(top_tbl),
#     TF   = tf_names,
#     direction = directions,
#     type = rep(celltypes[[i]], length(top_tbl)),
#     ptype = 'promoter'
#   )
# })

# bp_df <- do.call(rbind, bp_list)

# # Make neg values negative
# bp_df <- bp_df %>%
#   mutate(Freq = ifelse(direction == "neg", -Freq.Freq, Freq.Freq))

# # Create custom ordering per type
# bp_df <- bp_df %>%
#   group_by(type, direction) %>%
#   mutate(
#     TF_order = ifelse(direction == "pos",
#                       rank(-Freq, ties.method = "first"),  # descending
#                       rank(Freq, ties.method = "first"))   # ascending for neg (opposite)
#   ) %>%
#   ungroup()

# # Build a combined factor: ensures pos stack from bottom up, neg from top down
# bp_df <- bp_df %>%
#   arrange(type, direction, TF_order)

# bp_df$TF_id <- paste (bp_df$TF, bp_df$type, sep = "_")
# bp_df$TF_id <- factor (bp_df$TF_id, levels = unique(bp_df$TF_id))
# bp_df$type = factor (bp_df$type, levels = celltypes)
# bp_df_pp = bp_df
# # Plot stacked bars
# bp <- ggplot(bp_df_pp, aes(x = type, y = Freq, fill = TF_id)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF)) ) +
#   theme_minimal(base_size = 14) +
#   ylab("Proportion of counts") +
#   xlab("Cell type") +
#   ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

# pdf (file.path ('Plots', 'TF_abundance_profile_promoter_barplot.pdf'),6,width=12.5)
# bp
# dev.off()





# # Subset for peak type (NOT promoter) ####
# chrombpnet_counts_p2 =  lapply (chrombpnet_counts, function(x) {
#   x = x[x$peak_type == 'Distal',]
#   x = x[!is.na(x$peak_type),]
#   })
# chrombpnet_profile_p2 =  lapply (chrombpnet_profile, function(x) {
#   x = x[x$peak_type == 'Distal',]
#   x = x[!is.na(x$peak_type),]
#   })

# top_n <- 5
# n <- length(chrombpnet_counts_p2)

# bp_list <- lapply(seq_len(n), function(i) {
#   tbl <- table(chrombpnet_counts_p2[[i]]$V4)
#   tbl_sorted <- sort(tbl, decreasing = TRUE)
#   top_tbl <- head(tbl_sorted, top_n)
  
#   tf_names <- names(top_tbl)
#   directions <- sapply(tf_names, function(tf) {
#     chrombpnet_counts_p2[[i]]$V5[chrombpnet_counts_p2[[i]]$V4 == tf][1]
#   })
  
#   data.frame(
#     Freq = proportions(top_tbl),
#     TF   = tf_names,
#     direction = directions,
#     type = rep(celltypes[[i]], length(top_tbl)),
#     ptype = 'distal'
#   )
# })

# bp_df <- do.call (rbind, bp_list)

# # Make neg values negative
# bp_df <- bp_df %>%
#   mutate(Freq = ifelse(direction == "neg", -Freq.Freq, Freq.Freq))

# # Create custom ordering per type
# bp_df <- bp_df %>%
#   group_by(type, direction) %>%
#   mutate(
#     TF_order = ifelse(direction == "pos",
#                       rank(-Freq, ties.method = "first"),  # descending
#                       rank(Freq, ties.method = "first"))   # ascending for neg (opposite)
#   ) %>%
#   ungroup()

# # Build a combined factor: ensures pos stack from bottom up, neg from top down
# bp_df <- bp_df %>%
#   arrange(type, direction, TF_order)

# bp_df$TF_id <- paste(bp_df$TF, bp_df$type, sep = "_")
# bp_df$TF_id <- factor(bp_df$TF_id, levels = unique(bp_df$TF_id))
# bp_df$type = factor (bp_df$type, levels = celltypes)
# bp_df_cd = bp_df

# # Plot stacked bars
# bp <- ggplot(bp_df_cd, aes(x = type, y = Freq, fill = TF_id)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF)) ) +
#   theme_minimal(base_size = 14) +
#   ylab("Proportion of counts") +
#   xlab("Cell type") +
#   ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

  

# pdf (file.path ('Plots', 'TF_abundance_counts_NOT_promoter_barplot.pdf'),6,width=12.5)
# bp
# dev.off()


# top_n <- 5
# n <- length(chrombpnet_profile_p2)

# bp_list <- lapply(seq_len(n), function(i) {
#   tbl <- table(chrombpnet_profile_p2[[i]]$V4)
#   tbl_sorted <- sort(tbl, decreasing = TRUE)
#   top_tbl <- head(tbl_sorted, top_n)
  
#   tf_names <- names(top_tbl)
#   directions <- sapply(tf_names, function(tf) {
#     chrombpnet_profile_p2[[i]]$V5[chrombpnet_profile_p2[[i]]$V4 == tf][1]
#   })
  
#   data.frame(
#     Freq = proportions(top_tbl),
#     TF   = tf_names,
#     direction = directions,
#     type = rep(celltypes[[i]], length(top_tbl),
#     ptype = 'distal'
#   ))
#   )
# })

# bp_df <- do.call(rbind, bp_list)

# # Make neg values negative
# bp_df <- bp_df %>%
#   mutate(Freq = ifelse(direction == "neg", -Freq.Freq, Freq.Freq))

# # Create custom ordering per type
# bp_df <- bp_df %>%
#   group_by(type, direction) %>%
#   mutate(
#     TF_order = ifelse(direction == "pos",
#                       rank(-Freq, ties.method = "first"),  # descending
#                       rank(Freq, ties.method = "first"))   # ascending for neg (opposite)
#   ) %>%
#   ungroup()

# # Build a combined factor: ensures pos stack from bottom up, neg from top down
# bp_df <- bp_df %>%
#   arrange(type, direction, TF_order)

# bp_df$TF_id <- paste(bp_df$TF, bp_df$type, sep = "_")
# bp_df$TF_id <- factor(bp_df$TF_id, levels = unique(bp_df$TF_id))
# bp_df$type = factor (bp_df$type, levels = celltypes)
# bp_df_pd = bp_df

# # Plot stacked bars
# bp <- ggplot(bp_df_pd, aes(x = type, y = Freq, fill = TF_id)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF)) ) +
#   theme_minimal(base_size = 14) +
#   ylab("Proportion of counts") +
#   xlab("Cell type") +
#   ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#   geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

  

# pdf (file.path ('Plots', 'TF_abundance_profile_NOT_promoter_barplot.pdf'),6,width=12.5)
# bp
# dev.off()


# ### Combine barplots ####
# bp_df_comb = rbind (bp_df_cd, bp_df_cp)
# bp_df_comb$TF_id = as.character(bp_df_comb$TF_id)
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[1],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[2],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[3],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[4],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[5],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[6],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[7],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[8],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[9],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[10],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[11],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[12],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[grepl ('pattern',bp_df_comb$TF_id)] = 'pattern_not_matching'

# # Plot stacked bars
# bp <- ggplot(bp_df_comb, aes(x = type, y = Freq, fill = TF_id)) +
#   geom_bar(stat = "identity") +
# #  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF)) ) +
#   theme_minimal(base_size = 14) +
#   ylab("Proportion of counts") +
#   xlab("Cell type") +
#   ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) + 
#   facet_wrap (~ptype)

  

# pdf (file.path ('Plots', 'TF_abundance_combined_counts_barplot.pdf'),6,width=32.5)
# bp
# dev.off()

# bp_df_comb = rbind (bp_df_pd, bp_df_pp)
# bp_df_comb$TF_id = as.character(bp_df_comb$TF_id)
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[1],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[2],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[3],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[4],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[5],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[6],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[7],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[8],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[9],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[10],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[11],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] = gsub (celltypes[12],'',bp_df_comb$TF_id[!grepl ('pattern', bp_df_comb$TF_id)] )
# bp_df_comb$TF_id[grepl ('pattern',bp_df_comb$TF_id)] = 'pattern_not_matching'

# # Plot stacked bars
# bp <- ggplot(bp_df_comb, aes(x = type, y = Freq, fill = TF_id)) +
#   geom_bar(stat = "identity") +
# #  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF)) ) +
#   theme_minimal(base_size = 14) +
#   ylab("Proportion of counts") +
#   xlab("Cell type") +
#   ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) + 
#   facet_wrap (~ptype)

  

# pdf (file.path ('Plots', 'TF_abundance_combined_profile_barplot.pdf'),6,width=32.5)
# bp
# dev.off()



# ps = getPeakSet(archp)
# chrombpnet_counts = list()
# celltypes = c('Fibroblasts_P1')
# for (celltype in celltypes)
#   {
#   message (paste0('reading finemo output for ', celltype))  
#   chrombpnet_counts[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
#   gr = makeGRangesFromDataFrame (chrombpnet_counts[[celltype]], keep.extra.columns=T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
#   chrombpnet_counts[[celltype]]$peak_type = ps$peakType[findOverlaps(gr, ps, select='first')]
#   }







# # Import chromBPnet tfmodisco results to check for composite motifs ####
# #archp_P1 = archp[archp$Clusters %in% c('C2') & archp$Sample == 'P1']
# library ('universalmotif')

# ps = getPeakSet (archp)

# chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet'
# celltype = 'Mesothelium'
# library (httr)
# library (XML)


# # Load tfmodisco calls of contribution counts ####
# modisco_motifs = as.data.frame(readHTMLTable(file.path(chromBPdir,celltype, 'no_bias_model','modisco_profile',paste0(celltype,'_report'),'motifs.html')))

# # Check expression of profile motifs ####
# table (srt$celltype_lv1)

# genes = c('SP3','SP4','KLF12','FOS','FOSB','JUNB','ELK4','ELK1','GABPA','NFIC','NFIB','WT1','SALL4','EGR2','PITX2','NKX2-5','ZSCAN22','GLI1','GLI2','GLI3')
# pdf (file.path('Plots','chrombpnet_motifs_expression_dotplot.pdf'))
# DotPlot (srt, genes, group.by = 'celltype_lv1') + gtheme
# dev.off()



### Export fragment files of annotated cells for GEO submission ####
# Get fragments from each sample
library (rtracklayer)
fragments_l = getFragmentsFromProject (archp)
tmpdir = '/sc/arion/scratch/giottb01/meso_scatac_fragments'
for (sam in names (fragments_l))
  {
    if (!file.exists(file.path(tmpdir,paste0(sam,'_fragments.bed'))))
      export(fragments_l[[sam]], file.path(tmpdir, paste0(sam,'_fragments.bed')))
  }
    


