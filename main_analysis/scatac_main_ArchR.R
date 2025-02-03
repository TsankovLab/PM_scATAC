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

####### ANALYSIS of TUMOR compartment #######
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
addArchRThreads (threads = 8) 
addArchRGenome ("Hg38")

sample_names = c(
    # Tumor  
    'P1', # p786
    'P4', # p811
    'P8', # p826
    'P3', # p846
    'P5', #'p848'
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

#sarc_order = c('P1','P13','P3','P12','P5','P11','P4','P8','P14','P10')
#archp$Sample2 = archp$Sample
#archp$Sample2 = factor (archp$Sample2, levels = sarc_order)

### Gene score based analysis ####
run_GS_analysis = FALSE

if (run_GS_analysis)
  {
  # Find DAG ####
  metaGroupName = "Clusters"
  force = FALSE
  if (!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force)
    {
    DAG_list = getMarkerFeatures (
      ArchRProj = archp, 
      testMethod = "wilcoxon",
            #useGroups = "ClusterA",
            #bgdGroups = "Clusters1B",
      binarize = FALSE,
      useMatrix = "GeneScoreMatrix",
      groupBy = metaGroupName
    #  useSeqnames="z"
    )

    listnames = colnames (DAG_list)
    DAG_list = lapply (1:ncol (DAG_list), function(x) 
      {
      df = DAG_list[,x]  
      df = do.call (cbind, (assays(df)))
      colnames(df) = names (assays(DAG_list))
      df$gene = rowData (DAG_list)$name
      df
      })
    names (DAG_list) = listnames
    saveRDS (DAG_list, paste0 ('DAG_',metaGroupName,'.rds'))    
    } else {
    DAG_list = readRDS (paste0('DAG_',metaGroupName,'.rds'))
    }
  
  FDR_threshold = 1e-8
  lfc_threshold = 1
  top_genes = 20
  DAG_top_list = DAG_list[sapply (DAG_list, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$Log2FC) > lfc_threshold,]) > 0)]
  DAG_top_list = lapply (seq_along(DAG_top_list), function(x) {
    res = DAG_top_list[[x]]
    res = na.omit (res)
    res = res[res$FDR < FDR_threshold,]
    res = res[order (res$FDR), ]
    res = res[abs(res$Log2FC) > lfc_threshold,]
    res$comparison = names(DAG_top_list)[x]
    if (nrow(res) < top_genes) 
      {
      res
      } else {
      head (res,top_genes)
      }
    })
  DAG_df = Reduce (rbind ,DAG_top_list)
  
  if (!any (ls() == 'gsSE')) gsSE = ArchR::getMatrixFromProject (archp, useMatrix = 'GeneScoreMatrix')
  gsSE = gsSE[, archp$cellNames]
  gsMat = assays (gsSE)[[1]]
  rownames (gsMat) = rowData (gsSE)$name
  gsMat_mg = gsMat[rownames (gsMat) %in% DAG_df$gene, ]
  gsMat_mg = as.data.frame (t(gsMat_mg))
  gsMat_mg$metaGroup = as.character(archp@cellColData[,metaGroupName])
  gsMat_mg = aggregate (.~ metaGroup, gsMat_mg, mean)
  rownames (gsMat_mg) = gsMat_mg[,1]
  gsMat_mg = gsMat_mg[,-1]
  gsMat_mg = gsMat_mg[names(table (archp@cellColData[,metaGroupName])[table (archp@cellColData[,metaGroupName]) > 50]),]
  DAG_hm = Heatmap (t(scale(gsMat_mg)), 
          row_labels = colnames (gsMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 6),
          column_names_gp = gpar(fontsize = 4),
          rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE
          #right_annotation = motif_ha
          )

  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (paste0('Plots/DAG_clusters_',metaGroupName,'_heatmaps.pdf'), width = 8, height = 50)
print(DAG_hm)
dev.off()




# Plot gene score of cell type markers ####
meso_markers = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/gene_sets/highlevel_MPM_markers.csv')[[1]]
meso_markers = c(meso_markers, 'IGLL5')
meso_markers = meso_markers[meso_markers != 'IGHM']
#meso_markers = c(meso_markers, 'KRT5','LILRA4','MS4A1')
archp = addImputeWeights (archp)

pdf()
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = meso_markers, 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
#p = lapply (p, function(x) x + theme_void() + NoLegend ()) #+ ggtitle scale_fill_gradient2 (rev (viridis::plasma(100))))

pdf (file.path('Plots','marker_genes_feature_plots_3.pdf'), width = 25, height = 25)
print (wrap_plots (p, ncol = 8))
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
  

archp = saveArchRProject (archp, load=TRUE)
  
  metaGroupNames = c('TSSEnrichment','nFrags','ReadsInTSS','FRIP')  
    umap_p12 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
     name = x, embedding = "UMAP"))
      
  pdf (paste0(projdir,'/Plots/qc_umap_after_filtering.pdf'), 15,15)
  wrap_plots (umap_p12, ncol=5)
  dev.off()
  }

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
run_chromVAR_analysis = FALSE

if (run_chromVAR_analysis)
  {
  # Find DAM ####
  metaGroupName = "celltype_revised"
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
  
  # Filter by active genes as computed by correlation with gene score
  # active_genes = corGSM_MM$MotifMatrix_name[corGSM_MM$cor > -Inf]
  # DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% active_genes,])

  # Filter by genes minimally expressed from scRNA (in at least 1 celltype) ####
  ps = log2(as.data.frame (AverageExpression (srt, 
    features = sapply (unique(unlist(lapply(DAM_list, function(x) x$gene))), function(x) unlist(strsplit (x, '_'))[1]), 
    group.by = 'celltype_simplified2')[[1]]) +1)
  min_exp = 0.1
  ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
  selected_TF = rownames(ps)[rowSums(ps) > 0]

  DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% selected_TF,])

  names (DAM_list2) = names (DAM_list)
  FDR_threshold = 1e-3
  meandiff_threshold = 0
  top_genes = 3
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
  
devMethod = 'ArchR'
 if (devMethod == 'ArchR')
    {
    TF_db='Motif'
    if (!exists ('mSE')) mSE = ArchR::getMatrixFromProject (archp, useMatrix = paste0(TF_db,'Matrix'))
    mSE = mSE[, archp$cellNames]
    rowData(mSE)$name = gsub ('_.*','',rowData(mSE)$name)
    rowData(mSE)$name = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", rowData(mSE)$name)
    }
  DAM_df$gene = gsub ('_.*','',DAM_df$gene)
  DAM_df$gene = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", DAM_df$gene)
  mMat = assays (mSE)[[1]]
  rownames (mMat) = rowData (mSE)$name
  mMat_mg = mMat[DAM_df$gene, ]
  mMat_mg = as.data.frame (t(mMat_mg))
  mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
  mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
  rownames (mMat_mg) = mMat_mg[,1]
  mMat_mg = mMat_mg[,-1]
  mMat_mg = mMat_mg[names (DAM_list),]
  #mMat_mg = mMat_mg[names(table (archp@cellColData[,metaGroupName])[table (archp@cellColData[,metaGroupName]) > 50]),]
  DAM_hm = Heatmap (t(scale(mMat_mg)), 
          row_labels = colnames (mMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          name = 'chromVAR',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = rev(palette_deviation)

          #right_annotation = motif_ha
          )

  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_heatmaps.pdf')), width = 3, height = 5)
print(DAM_hm)
dev.off()
}


### Co-expression of TFs #### 
metaGroupName = 'celltype_revised'
if (!any (ls() == 'mSE')) mSE = fetch_mat (archp, 'Motif')

all (colnames(mSE) == rownames(archp@cellColData))

# # Get deviation matrix ####
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Filter by genes minimally expressed from scRNA (in at least 1 celltype) ####
ps = log2(as.data.frame (AverageExpression (srt, 
  features = sapply (unique(unlist(lapply(DAM_list, function(x) x$gene))), function(x) unlist(strsplit (x, '_'))[1]), 
  group.by = 'celltype_simplified2')[[1]]) +1)
min_exp = 0.1
ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
selected_TF = rownames(ps)[rowSums(ps) > 0]
# positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0.1]
mMat = mMat[selected_TF,]

# mMat_agg = as.data.frame (t(mMat))
# mMat_agg$metaGroup = as.character (archp_meta[,metaGroupName])
# mMat_agg = aggregate (.~ metaGroup, mMat_agg, mean)
# rownames (mMat_agg) = mMat_agg[,1]
# mMat_agg = mMat_agg[,-1]
# mMat_agg = t(mMat_agg)
# rownames (mMat_agg) = active_TF

mMat_cor = cor (as.matrix(t(mMat)), method = 'pearson')
#d = as.dist (1-cor(as.matrix(t(mMat))))

#d = dist (mMat, method ='euclidean')
#hc1 <- hclust(d, method = "complete" ) # Hierarchical clustering using Complete Linkage

# row_filt = rowSums (mMat_cor) != 0
# tf_name = rownames(mMat_cor)[row_filt]
km = kmeans (mMat_cor, centers=20)

pdf (file.path ('Plots','TF_modules.pdf'), width = 4,height=3)
cor_mMat_hm = draw (Heatmap (mMat_cor,# row_km=15,
  #left_annotation = ha,
  #rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=palette_deviation_correlation, 
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
dev.off()

pdf (file.path ('Plots','TF_modules.pdf'), width = 4,height=3)
cor_mMat_hm
dev.off()

# Make heatmap of TF by cell types ####
pan_TF = c('TGIF2','TGIF1','NFKB2','SNAI2','TWIST2','HMGA2')
ha = HeatmapAnnotation (celltype = archp$celltype_revised, col = list(celltype = palette_celltype_simplified))
ha2 = rowAnnotation (foo = anno_mark(at = match(pan_TF,rownames(mMat2)), 
    labels = pan_TF, labels_gp = gpar(fontsize = 6, fontface = 'italic')))
mMat = as.matrix (mMat)
mMat2 = mMat
mMat2[mMat2 > 1] = 1
mMat2[mMat2 < -1] = -1
pdf (file.path ('Plots','TF_modules_celltypes_heatmap.pdf'), width = 4,height=3)
cor_mMat_hm = draw (Heatmap (t(scale(t(mMat2))),# row_km=15,
  top_annotation = ha,
  right_annotation = ha2,
  #left_annotation = ha,
  #rect_gp = gpar(type = "none"),
  clustering_distance_rows='pearson' ,
  clustering_distance_columns = 'pearson', 
  col=palette_deviation_cor_fun, 
  #column_split = km$cluster,
  row_km=3, 
  #column_km=2,
#  right_annotation = ha,
  border=T,
  row_names_gp = gpar(fontsize = 0),
  column_names_gp = gpar(fontsize = 0)))#,
  # ,
  # cell_fun = function(j, i, x, y, w, h, fill) {
  #       if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
  #           grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#        }}))
dev.off()

pdf (file.path ('Plots','TF_modules_celltypes_heatmap.pdf'), width = 8,height=4)
cor_mMat_hm
dev.off()


# Find shared TF across celltypes ####
tf_name2 = unlist(sapply (c('TGIF2','TGIF1','TWIST2','NFKB2','HMGA2'), function(x) rownames(assay(mSE))[grepl (x, rownames(assay(mSE)))]))
tf_name2 = paste0('z:',tf_name2)
archp = addImputeWeights (archp)
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "MotifMatrix",
    name = tf_name2, 
    useSeqnames='z',
    pal = rev (palette_deviation),    
    embedding = "UMAP",
    imputeWeights = getImputeWeights(archp)
    )

pdf (file.path ('Plots','pan_TF_fplots.pdf'), width = 30,height=16)
wrap_plots (TF_p, ncol=5)
dev.off()


tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = do.call (cbind, tf_modules)
archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% paste0('mod_',unique(km$cluster))]
archp@cellColData = cbind (archp@cellColData, tf_modules)

archp = addImputeWeights (archp)
TF_p = lapply (paste0('mod_',unique(km$cluster)), function(x) plotEmbedding (
    ArchRProj = archp,
    colorBy = "cellColData",
    name = x, 
    pal = palette_deviation,
    #useSeqnames='z',
    imputeWeights = getImputeWeights(archp),
    embedding = "UMAP"))

pdf (file.path ('Plots','TF_modules_umap.pdf'), width = 30,height=14)
wrap_plots (TF_p, ncol=5)
dev.off()

#### make heatmap of all positively regulated TFs across cell types to find shared programs ####
# # Get deviation matrix ####
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = scale(assays (mSE)[[1]])
rownames (mMat) = rowData (mSE)$name

# Filter for only expressed TFs ####
metaGroupName = 'celltype_simplified2'
ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat), group.by = metaGroupName)[[1]]) +1)
min_exp = .1
active_TFs = rownames(ps)[rowSums(ps) > min_exp]
mMat = mMat[active_TFs,]
mMat = as.data.frame (t(mMat))
#mMat$metaGroup = as.character (archp@cellColData[,metaGroupName])
#mMat = aggregate (.~ metaGroup, mMat, mean)
#rownames(mMat) = mMat[,1]
#mMat = mMat[,-1]
metaGroupName = 'celltype_lv1'
ha = HeatmapAnnotation (celltype = as.character (archp@cellColData[,metaGroupName]), col=list(celltype = palette_celltype_simplified))
DAM_hm = Heatmap (t(mMat), 
          row_labels = colnames (mMat),
          top_annotation = ha,
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = T,
          #col = pals_heatmap[[5]],
          cluster_columns=T,#col = pals_heatmap[[1]],
          row_names_gp = gpar (fontsize = 0),
          column_names_gp = gpar(fontsize = 0),
          column_names_rot = 45,
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_deviation_centered
          #right_annotation = motif_ha
          )

pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_all_TF_heatmaps.pdf')), width = 12.6, height = 5)
print(DAM_hm)
dev.off()











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


# Show that enhancer linked to NR4A2 is only up in NK KLRC1 and CD8 exhausted across all cells ####
nkt_ann = read.csv (file.path('..','..','NKT_cells','scatac_ArchR','barcode_annotation_nkt.csv'))
mye_ann = read.csv (file.path('..','..','myeloid_cells','scatac_ArchR','barcode_annotation_main.csv'))
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
#+ scale_fill_manual (values = c(palette_tnk_cells, palette_myeloid, palette_celltype_simplified))
ep
dev.off()


# Check expression of NR4A2 across cell types ####
srt_tnk = readRDS (file.path('..','..','NKT_cells','scrna','srt.rds'))
tnk_ann = srt_tnk$celltype2
tnk_ann = tnk_ann[names(tnk_ann) %in% colnames(srt)]
srt$celltype3 = srt$celltype
srt$celltype3[match (names(tnk_ann), colnames(srt))] = tnk_ann

srt$celltype_simplified3 = srt$celltype_simplified
srt$celltype_simplified3[match (names(tnk_ann), colnames(srt))] = tnk_ann
pdf(file.path ('Plots','NR4A2_dotplot.pdf'), width = 9)
VlnPlot (srt, feature = 'NR4A2', group.by = 'celltype3', cols = c(palette_tnk_cells, palette_celltype_simplified))
VlnPlot (srt[,!srt$celltype_simplified3 %in% c('NK','T_cells')], feature = 'NR4A2', group.by = 'celltype_simplified3', cols = c(palette_tnk_cells, palette_celltype_simplified))
dev.off()


# Input exhausted peaks as chromvar score ####
T_ext = list(
  T_ext = readRDS (file.path('..','..','NKT_cells','scatac_ArchR','T_cell_exhaustion_peaks.rds')),
  T_ext2 = readRDS (file.path('..','..','NKT_cells','scatac_ArchR','T_cell_exhaustion_peaks.rds')))

# Import chrombpnet output and cross-reference with scRNA-seq TF expression
archp = addBgdPeaks (archp, force= TRUE)
archp = addPeakAnnotations (ArchRProj = archp, 
     regions = T_ext, name = "T_exhausted",force=T)

archp = addDeviationsMatrix (
  ArchRProj = archp, 
  peakAnnotation = "T_exhausted",
  force = TRUE
)

pdf()
if (!any (ls() == 'tSE')) tSE = fetch_mat (archp, 'T_exhausted')
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "T_exhaustedMatrix",
    name = 'z:T_ext', 
    #useSeqnames='z',
    pal = rev (palette_deviation),    
    embedding = "UMAP",
    imputeWeights = getImputeWeights(archp)
    )
dev.off()
pdf (file.path ('Plots','T_exhausted_fplots.pdf'), width = 30,height=16)
TF_p
dev.off()





### Compare TF expression from scRNA and inferred by chrombpnet per cell type ####
# library (httr)
# library (XML)
# library (igraph)
#BiocManager::install("universalmotif")
library ('universalmotif')


# # Differential Accessed motifs ####
metaGroupName = "celltype_lv1"
force=FALSE
# source (file.path('..','..','git_repo','utils','DAM.R'))
DAM_list = readRDS (paste0('DAM_',metaGroupName,'.rds'))
DAM_list = lapply (DAM_list, function(x)
       {
       x$gene = gsub ('_.*','',x$gene)
       x$gene = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", x$gene)
       x
       })

metaGroupName = 'celltype_lv1'
if (!any (ls() == 'mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
#mMat_mg = mMat[DAM_df$gene, ]
mMat_mg = as.data.frame (t(mMat))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]


#Get active genes from RNA
metaGroupName = 'celltype_simplified2'
ps = log2(as.data.frame (AverageExpression (srt, 
features = colnames(mMat_mg),
group.by = metaGroupName)[[1]]) +1)
#min_exp = .1
#ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
#active_TFs = rownames(ps)[rowSums(ps) > 0]

#active_genes = corGSM_MM$MotifMatrix_name[corGSM_MM$cor > 0.1]
#DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% active_TFs,])    
mMat_l = as.list (as.data.frame (t(mMat_mg)))
mMat_l = lapply (mMat_l, function(x) data.frame (dev = x, row.names = colnames(mMat_mg)))
#mMat_l = lapply (mMat_l, function(x) x[rownames(x) %in% active_TFs,,drop=F])

metaGroupName = 'celltype_lv1'
chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet'
metaGroupName = 'celltype_lv1'
celltypes = unique (archp@cellColData[,metaGroupName])

tf_database = read_meme('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/HOCOMOCO_db/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme', skip = 0, readsites = FALSE, readsites.meta = FALSE)
tf_database = unique(unlist(lapply(tf_database, function(x) unlist(strsplit(x@name,'_'))[1])))

list.files (file.path(chromBPdir, celltypes[3],'no_bias_model'))

DAM_tf_l = list()
sp_l = list()
celltypes=c('SmoothMuscle','Mesothelium','Alveolar','Plasma','Fibroblasts','NK','Myeloid')
ap1_complex = c('FOS','FOSL2','FOSL1','JDP2','JUN','JUND','JUNB','FOSB')

for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_count_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))

  chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  chrombpnet_tf = rbind (chrombpnet_count_tf, chrombpnet_profile_tf)
  chrombpnet_tf = chrombpnet_tf[!duplicated(chrombpnet_tf$V4),]
  #DAM_tf = DAM_list2[[celltype]]
  DAM_tf = mMat_l[[celltype]]
  DAM_tf$expression = ps[rownames(DAM_tf), celltype]
  DAM_tf = DAM_tf[order (-DAM_tf$dev),,drop=F]
  DAM_tf = DAM_tf[!rownames(DAM_tf) %in% ap1_complex,]
  DAM_tf = DAM_tf[DAM_tf$expression > min_exp,]
  DAM_tf$pioneer = sapply (seq(nrow(DAM_tf)), function(x) rownames(DAM_tf)[x] %in% unlist(strsplit(chrombpnet_tf[, 'V4'],'_')))
  DAM_tf$pioneer_TF = ifelse (DAM_tf$pioneer, rownames(DAM_tf),'')
  DAM_tf$celltype = celltype
  DAM_tf_l[[celltype]] = DAM_tf
  }

#DAM_tf_l = lapply (DAM_tf_l, function(x) x[order(-x$dev),,drop=F])
DAM_tf2 = do.call (rbind, DAM_tf_l)
DAM_tf2 = DAM_tf2[!DAM_tf2$pioneer_TF %in% ap1_complex,]
DAM_tf3 = do.call (rbind, lapply (split(DAM_tf2, DAM_tf2$celltype), function(x){head(x[x$pioneer_TF !='',],10)}))
DAM_tf2$alpha = ifelse (DAM_tf2$pioneer, 0.4,.2)
sp_l = ggplot(DAM_tf2, aes(x = celltype, y = dev), linewidth=0.3) +  # x=1 makes it similar to a boxplot
geom_jitter(shape=21, aes (fill = pioneer, alpha =alpha, size=expression),width = 0.2, color='white') +  # Jitter prevents overlap
geom_text_repel(data = DAM_tf3, aes(label = pioneer_TF), 
                size = 2, nudge_x = 0.2, na.rm = TRUE) +
scale_fill_manual(values = c("TRUE" = "brown", "FALSE" = "grey")) +  # Custom colors
gtheme_no_rot

pdf (file.path ('Plots','Pioneer_TF_celltype_scatter.pdf'), height=4,width=10)
sp_l
dev.off()


## Show barplots of top TF occurrence using finemo chrombpnet outputs ####
chrombpnet_cl = list()
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_count_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
  #chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  #chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  #chrombpnet_tf = rbind (chrombpnet_count_tf, chrombpnet_profile_tf)
  chrombpnet_tf = table (chrombpnet_count_tf$V4)[order(-table (chrombpnet_count_tf$V4))]
  assign_max_exp = unlist(sapply (names(chrombpnet_tf), function(x) unlist(strsplit(x, '_'))[which.max(ps[unlist(strsplit(x, '_')), celltype])]))
  tf_dev = mMat_l[[celltype]][assign_max_exp,]
  chrombpnet_df = data.frame (occurrence = chrombpnet_tf[names(assign_max_exp)], TF_max_exp = assign_max_exp, TF_max_dev = tf_dev)
  chrombpnet_df = chrombpnet_df[!chrombpnet_df$TF_max_exp %in% ap1_complex,]
  chrombpnet_df = chrombpnet_df[!duplicated(chrombpnet_df$TF_max_exp),]
  chrombpnet_df = chrombpnet_df[chrombpnet_df$TF_max_exp %in% head(chrombpnet_df$TF_max_exp[order(-chrombpnet_df$TF_max_dev)],5),]
  chrombpnet_df$celltype = celltype
  chrombpnet_df$order = seq(nrow(chrombpnet_df))
  chrombpnet_cl[[celltype]] = chrombpnet_df
  }

chrombpnet_cl_df = do.call (rbind, chrombpnet_cl)
chrombpnet_cl_df$TF_max_exp = factor (chrombpnet_cl_df$TF_max_exp, levels =unique(chrombpnet_cl_df$TF_max_exp))
chrombpnet_cl_df$order = factor (chrombpnet_cl_df$order, levels =unique(chrombpnet_cl_df$order))
# Create stacked bar plot with text beside each band
bp = ggplot(chrombpnet_cl_df, aes(x = celltype, y = occurrence.Freq, fill = order)) +
  geom_bar(stat = "identity", color = 'black') +
  geom_text(aes(label = TF_max_exp), 
            position = position_stack(vjust = 0.5), 
            hjust = 0.4,  # Move text outside the bar
            size = 4) + 
  #coord_flip() +  # Flip to make text more readable
  gtheme
pdf (file.path ('Plots','chrombpnet_count_TF_barplot.pdf'))
bp
dev.off()

chrombpnet_cl = list()
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_count_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  #chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  #chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  #chrombpnet_tf = rbind (chrombpnet_count_tf, chrombpnet_profile_tf)
  chrombpnet_tf = table (chrombpnet_count_tf$V4)[order(-table (chrombpnet_count_tf$V4))]
  assign_max_exp = unlist(sapply (names(chrombpnet_tf), function(x) unlist(strsplit(x, '_'))[which.max(ps[unlist(strsplit(x, '_')), celltype])]))
  tf_dev = mMat_l[[celltype]][assign_max_exp,]
  chrombpnet_df = data.frame (occurrence = chrombpnet_tf[names(assign_max_exp)], TF_max_exp = assign_max_exp, TF_max_dev = tf_dev)
  chrombpnet_df = chrombpnet_df[!chrombpnet_df$TF_max_exp %in% ap1_complex,]
  chrombpnet_df = chrombpnet_df[!duplicated(chrombpnet_df$TF_max_exp),]
  chrombpnet_df = chrombpnet_df[chrombpnet_df$TF_max_exp %in% head(chrombpnet_df$TF_max_exp[order(-chrombpnet_df$TF_max_dev)],5),]
  chrombpnet_df$celltype = celltype
  chrombpnet_df$order = seq(nrow(chrombpnet_df))
  chrombpnet_cl[[celltype]] = chrombpnet_df
  }

chrombpnet_cl_df = do.call (rbind, chrombpnet_cl)
chrombpnet_cl_df$TF_max_exp = factor (chrombpnet_cl_df$TF_max_exp, levels =unique(chrombpnet_cl_df$TF_max_exp))
chrombpnet_cl_df$order = factor (chrombpnet_cl_df$order, levels =unique(chrombpnet_cl_df$order))
# Create stacked bar plot with text beside each band
bp = ggplot(chrombpnet_cl_df, aes(x = celltype, y = occurrence.Freq, fill = order)) +
  geom_bar(stat = "identity", color = 'black') +
  geom_text(aes(label = TF_max_exp), 
            position = position_stack(vjust = 0.5), 
            hjust = 0.4,  # Move text outside the bar
            size = 4) + 
  #coord_flip() +  # Flip to make text more readable
  gtheme
pdf (file.path ('Plots','chrombpnet_profile_TF_barplot.pdf'))
bp
dev.off()

