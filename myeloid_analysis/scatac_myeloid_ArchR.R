conda activate meso_scatac

R

set.seed(1234)

####### ANALYSIS of Myeloid compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/'
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

# Load ArchR project ####
archp = loadArchRProject (projdir)

# Load RNA ####
srt = readRDS (file.path('..','scrna','srt.rds'))
sample_levels = c('Monocytes','cDCs','SPP1','TREM2','C1Q','IFN','IM')


## Reduce dimension and harmonize ####
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

  archp = addHarmony (
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    name = "Harmony_project",
    groupBy = c('Sample'), force=TRUE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony_project", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony_project",
    name='Clusters_H', res=1,
    force = TRUE)


### Annotate cell types ####
archp$celltype2 = 0
archp$celltype2[archp$Clusters_H %in% c('C6','C7','C9')] = 'Monocytes'
archp$celltype2[archp$Clusters_H %in% c('C8')] = 'cDCs'
archp$celltype2[archp$Clusters_H %in% c('C1','C2','C3','C4','C5','C10','C11')] = 'TAMs'

pdf()
umap_p3 = plotEmbedding (ArchRProj = archp, labelMeans = F,
  colorBy = "cellColData", name = "Sample",
  pal= palette_sample,
   embedding = "UMAP_H")
umap_p4 = plotEmbedding (ArchRProj = archp, labelMeans = F,
  colorBy = "cellColData", name = "celltype2",
   embedding = "UMAP_H",
   pal = palette_myeloid2)
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters_H",
   embedding = "UMAP_H")
dev.off()

pdf (file.path('Plots','celltype_umap_harmony_sample_umap.pdf'),5,5)
print (umap_p3)
print (umap_p4)
print (umap_p5)
dev.off()


# # Check for doublets ####
# meso_markers = c('CD3D','CD3E','EPCAM','KRT19','KRT5','VWF','PECAM1')
# #archp = addImputeWeights (archp)
# pdf()
# p <- plotEmbedding(
#     ArchRProj = archp,
#     colorBy = "GeneScoreMatrix", 
#     name = meso_markers, 
#     size=1,
#     embedding = "UMAP_H",
#     pal = palette_expression,
#     imputeWeights = getImputeWeights(archp)
# )
# dev.off()
# pdf (file.path('Plots','doublets_fplots.pdf'), width = 18, height = 15)
# wrap_plots (p, ncol=3)
# dev.off()


macs_markers=read.csv (file.path ('..','..','..','..','gene_sets','scRNA_immune_markers_humanLUAD_Samarth_Assaf.csv'))
mono_markers = macs_markers[macs_markers$group == 'CD14 mono','gene']
mono_markers = mono_markers[mono_markers != 'CD14 MONO']
mono_markers = c('VCAN','FCN1','CXCL8','CXCL2','IL1B','EREG','TIMP1','THBS1','CCR2','FLT3','FOXM1','CDK1','PCNA','FCGR3A','C1QA','C1QB','CD1C','CD1A','FCER1A')
dc_markers = macs_markers[macs_markers$group %in% c('DC1','DC2','mregDC'),'gene']
ap1_complex = c('JUN','FOSB','FOS','BACH1','SMARCC1','FOSL2','JUND','JDP2','BATF','CEBPB','CEBPA','CEBPZ','FOSL1','NFE2','NFE2L2','NFE2L1')
selected_markers = c('FCN1','EREG','TIMP1', 'CCR7','LAMP3','HLA-DQA1','C1QA','CD68','C5AR2')
pdf()
p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = selected_markers, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path ('Plots','TF_selected_featureplots.pdf'), width = 18,height=14)
wrap_plots (p2)
dev.off()

pdf()
p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = ap1_complex, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path ('Plots','TF_ap1_complex_featureplots.pdf'), width = 28,height=24)
wrap_plots (p2)
dev.off()

# Check sample quality
qc = c('ReadsInTSS','TSSEnrichment','nFrags')
pdf()
p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "cellColData", 
    name = qc, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = NULL
)
dev.off()
pdf (file.path ('Plots','qc_featureplots.pdf'), width = 18,height=14)
wrap_plots (p2)
dev.off()



# # Get markers for gene score ####
# immune_markers = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/scRNA_immune_markers_humanLUAD_Samarth_Assaf.csv')
# immune_markers = immune_markers [immune_markers$group %in% c('Neutrophil','TRMac','IM','DC2','DC1','pDC','mregDC','CD14 mono',
#   'CD16 mono','NK','Mast cell','Mgk','B/Plasma',' T cell','Treg','MoMac'),]
# #immune_markers = immune_markers[immune_markers$group %in% c('CD14 mono','CD16 mono','DC1','DC2','MoMac'),]
# immune_markers = immune_markers$gene
# immune_markers = immune_markers[!immune_markers %in% c('CD14 MONO','IHBA','SEPP1','IL3RA')]
# #archp = addImputeWeights (archp)
# pdf()
# p <- plotEmbedding(
#     ArchRProj = archp,
#     colorBy = "GeneScoreMatrix", 
#     name = immune_markers, 
#     size=1,
#     embedding = "UMAP_H",
#     pal = palette_expression,
#     imputeWeights = getImputeWeights(archp)
# )
# dev.off()
# png (file.path('Plots','myeloid_markers_fplots.png'), width = 18000, height = 15000, res=300)
# wrap_plots (p)
# dev.off()


### Call peaks on celltypes ####
pdf(file.path('Plots','peakcalls.pdf'))
metaGroupName = 'Clusters_H'
force=TRUE
peak_reproducibility=2
if(!all(file.exists(file.path('PeakCalls', unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds'))) | force) 
source (file.path('..','..','git_repo','utils','callPeaks.R'))
dev.off()

### chromVAR analysis ####
force=TRUE
if (!all(file.exists(file.path('Annotations',
  c('Motif-Matches-In-Peaks.rds',
    'Motif-Positions-In-Peaks.rds',
    'Motif-In-Peaks-Summary.rds')))) | force)
source (file.path ('..','..','git_repo','utils','chromVAR.R'))
  


# Add cNMF modules from scRNA-seq ####
shared_cnmf = readRDS (file.path('..','scrna','shared_cnmf_myeloid.rds'))
shared_cnmf = lapply (shared_cnmf, function(x) x[x %in% getFeatures (archp)])
remove_modules = c('cnmf.3','cnmf.6','cnmf.7','cnmf.5') # remove monocyres cDC and CC modules. Consider re-inculding CC 
shared_cnmf_MAC = shared_cnmf[!names(shared_cnmf) %in% remove_modules]

force = FALSE
if (!all (names (shared_cnmf) %in% colnames (archp@cellColData)) | force)
  {
  archp@cellColData = archp@cellColData[,!grepl ('cnmf',colnames(archp@cellColData))]
  archp = addModuleScore (
      ArchRProj = archp,
      useMatrix = 'GeneScoreMatrix',
      name = '',
      features = shared_cnmf,
      nBin = 25,
      nBgd = 100,
      seed = 1,
      threads = getArchRThreads(),
      logFile = createLogFile("addModuleScore")
    )
  colnames (archp@cellColData) = gsub ('^\\.','',colnames(archp@cellColData))    
  }

# Assign TAMs to cNMF modules from scRNA-seq ####
archp_MAC = archp[!archp$celltype2 %in% c('Monocytes','cDCs')]
cnmf_scatac = as.data.frame (t(scale(t(archp_MAC@cellColData[,names(shared_cnmf_MAC)]))))

set.seed (123)
km = kmeans (scale(cnmf_scatac), centers=5) # double scale modules and cluster using k-means

ha = HeatmapAnnotation (sample = archp_MAC$Sample, col=list(sample = palette_sample))
hm = Heatmap (t(scale(cnmf_scatac)), 
  col = palette_genescore_fun(scale(cnmf_scatac)), 
  top_annotation = ha,
  show_column_dend = F,
  column_split = km$cluster,
#  column_km=9,
  row_names_gp = gpar (fontsize = 8),
  column_names_gp = gpar (fontsize = 0),
  border=T)

pdf (file.path ('Plots','cnmf_scatac_scaled_only_MAC2_heatmap2.pdf'), height=2.5)
hm
dev.off()


# Re-Annotate based on cnmf clustering ####
all (names(km$cluster) == rownames(archp_MAC@cellColData))
archp_MAC$cnmf_cluster = paste0('cnmf_cluster_',km$cluster)
archp_MAC$cnmf_celltypes = archp_MAC$cnmf_cluster
archp_MAC$cnmf_celltypes[archp_MAC$cnmf_cluster == 'cnmf_cluster_1'] = 'C1Q'
archp_MAC$cnmf_celltypes[archp_MAC$cnmf_cluster == 'cnmf_cluster_2'] = 'SPP1'
archp_MAC$cnmf_celltypes[archp_MAC$cnmf_cluster == 'cnmf_cluster_3'] = 'IFN'
archp_MAC$cnmf_celltypes[archp_MAC$cnmf_cluster == 'cnmf_cluster_4'] = 'TREM2'
archp_MAC$cnmf_celltypes[archp_MAC$cnmf_cluster == 'cnmf_cluster_5'] = 'IM'

# Integrate MACs with myeloid annotation
archp$celltype2[match(rownames(archp_MAC@cellColData), rownames(archp@cellColData))] = archp_MAC$cnmf_celltypes

pdf()
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype2",
   embedding = "UMAP_H")
dev.off()

pdf (file.path('Plots','celltype_umap_TAM_annotated_umap.pdf'),5,5)
print (umap_p5)
dev.off()


write.csv (data.frame (barcode = rownames(archp@cellColData), celltype = archp$celltype2), 'barcode_annotation.csv')

# Find DAM in cnmf_clusters ####
sample_levels = c('Monocytes','cDCs','SPP1','TREM2','C1Q','IFN','IM')
#metaGroupName = "celltype2"
metaGroupName = "celltype2"
#archp$celltype2 = archp$cnmf_cluster
force=F
source (file.path('..','..','git_repo','utils','DAM.R'))

DAM_list = readRDS (paste0('DAM_',metaGroupName,'.rds'))
# Clean TF names
  DAM_list = lapply (DAM_list, function(x)
       {
       x$gene = gsub ('_.*','',x$gene)
       x$gene = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", x$gene)
       x
       })

if (metaGroupName %in% colnames(srt@meta.data))
{
#Get active genes from RNA
ps = log2(as.data.frame (AverageExpression (srt, 
features = sapply (unique(unlist(lapply(DAM_list, function(x) x$gene))), function(x) unlist(strsplit (x, '_'))[1]), 
group.by = metaGroupName)[[1]]) +1)
min_exp = .1
ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
active_TFs = rownames(ps)
#active_genes = corGSM_MM$MotifMatrix_name[corGSM_MM$cor > 0.1]
DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% active_TFs,])    
} else {
DAM_list2 = DAM_list  
}

names (DAM_list2) = names (DAM_list)
DAM_list2 = DAM_list2[sample_levels]
FDR_threshold = 0.05
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
active_DAM = unique(DAM_df$gene)
  
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

metaGroupName = 'celltype2'
mMat_mg = mMat[active_DAM, ]
mMat_mg = as.data.frame (t(mMat_mg))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
mMat_mg = mMat_mg[sample_levels,]
#mMat_mg = mMat_mg[names (DAM_list),]

#selected_TF = c(rownames(DAM_hm@matrix), 'NR4A3','NR4A2','NR4A1')

DAM_hm = Heatmap (t(scale(mMat_mg)), 
          row_labels = colnames (mMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          name = 'chromVAR',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = rev(palette_deviation)#,
          #width = unit(2, "cm")
          #right_annotation = motif_ha
          )

pdf (file.path ('Plots','cnmf_clusters_DAM_heatmap2.pdf'), width = 3,height=4)
draw (DAM_hm)
dev.off()

# Find DAG in cnmf_clusters ####
metaGroupName = "celltype2"
#metaGroupName = "celltype2"
#archp$celltype2 = archp$cnmf_cluster
force=FALSE
source (file.path('..','..','git_repo','utils','DAG.R'))

top_genes = 5
DAG_top_list = DAG_list[sapply (DAG_list, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$Log2FC) > lfc_threshold,]) > 0)]
DAG_top_list = lapply (seq_along(DAG_top_list), function(x) {
  res = DAG_top_list[[x]]
  #res = na.omit (res)
  res = res[res$FDR < FDR_threshold,]
  res = res[order (res$FDR), ]
  res = res[res$Log2FC > lfc_threshold,]
  res$comparison = names(DAG_top_list)[x]
  if (nrow(res) < top_genes) 
    {
    res
    } else {
    head (res,top_genes)
    }
  })
DAG_df = Reduce (rbind ,DAG_top_list)
head (DAG_df[DAG_df$comparison == 'IM',],20)

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
DAG_hm = Heatmap (t(scale(gsMat_mg[,DAG_df$gene])), 
        row_labels = DAG_df$gene,
        column_title = paste('top',top_genes),
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        cluster_rows = F,
        cluster_columns=F,
        col = palette_expression,
        #cluster_columns=T,#col = pals_heatmap[[1]],
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        rect_gp = gpar(col = "white", lwd = .5),
        border=TRUE,
        column_names_rot=45
        #right_annotation = motif_ha
        )
         
  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (paste0('Plots/DAG_clusters_',metaGroupName,'_heatmaps2.pdf'), width = 2.5, height = 5)
print (DAG_hm)
dev.off()

  

# Show TAM modules in UMAPs
archp = addImputeWeights (archp)
pdf()
p <- plotEmbedding (
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = names (shared_cnmf_MAC), 
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path ('Plots','shared_cnmf_TAMs_fplots.pdf'),16,16)
wrap_plots (p, ncol=3)
dev.off()

  






### Co-expression of TFs across cells #### 

### Run TF correlation to identify TF modules across TNK cells #### 
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat)#[selected_TF,])

# Filter by RNA expression ####
metaGroupName = 'celltype2'
active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)
mMat = mMat[active_TFs, ]

mMat_cor = cor (as.matrix(t(scale(mMat))), method = 'spearman')

set.seed(1234)
centers=2
km = kmeans (mMat_cor, centers=centers)

pdf (file.path ('Plots','TF_modules_heatmap2.pdf'), width = 4,height=3)
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

pdf (file.path ('Plots','TF_modules_heatmap2.pdf'), width = 4,height=3)
cor_mMat_hm
dev.off()


# Generate same heatmap but using scrna ####
metacells = readRDS (file.path ('..','scrna','metacells.rds'))
metacells_mat = metacells@assays$RNA$data[rownames(mMat_cor),]
metacells_mat = cor (t(metacells_mat), method = 'spearman')
all (rownames(metacells_mat) == rownames(mMat_cor))
tf_order = unname(unlist(row_order(cor_mMat_hm)))

pdf (file.path ('Plots','TF_modules_RNA_heatmap2.pdf'), width = 4,height=3)
cor_mMat_hm2 = draw (Heatmap (metacells_mat[tf_order,tf_order],# row_km=15,
  row_split = km$cluster,
  column_split = km$cluster,
  cluster_rows=F,
  cluster_columns = F,
  col=palette_expression_cor_fun, 
  border=T,
  row_names_gp = gpar(fontsize = 0), 
  column_names_gp = gpar(fontsize = 0)))
dev.off()

pdf (file.path ('Plots','TF_modules_ATAC_RNA_heatmap2.pdf'), width = 4,height=3)
cor_mMat_hm2
dev.off()


pdf (file.path ('Plots','TF_modules_RNA_heatmap3.pdf'), width = 4,height=3)
cor_mMat_hm2 = draw (Heatmap (metacells_mat[tf_order,tf_order],# row_km=15,
  col=palette_expression_cor_fun, 
  border=T,
  row_names_gp = gpar(fontsize = 0), 
  column_names_gp = gpar(fontsize = 0)))
dev.off()

pdf (file.path ('Plots','TF_modules_ATAC_RNA_heatmap3.pdf'), width = 4,height=3)
cor_mMat_hm2
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
    pal = rev(palette_deviation_correlation),
    imputeWeights=NULL,
    #useSeqnames='z',
    embedding = "UMAP_H")
dev.off()

pdf (file.path ('Plots','TF_modules_umap2.pdf'), width = 20,height=16)
wrap_plots (TF_p, ncol=5)
dev.off()

# Show all TFs included in inflammation module ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat)#[selected_TF,])

metaGroupName = 'celltype2'
mMat_mg = mMat[names (km$cluster)[km$cluster==2], ]
mMat_mg = as.data.frame (t(mMat_mg))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
hm = Heatmap (
    t(mMat_mg),
    column_names_rot =45, 
    row_names_gp = gpar(fontsize = 4),
    column_names_gp = gpar(fontsize = 8),
    col = rev(as.character(palette_deviation)), 
    cluster_rows=T,
    cluster_columns = T#,
#rect_gp = gpar (col = "white", lwd = 1)
)

pdf (file.path ('Plots','inflammation_module_TFs_heatmap.pdf'), width = 2.3,height=8)
hm
dev.off()


### Check inflammation score across TAMs ####
mod_df = data.frame (
  celltype = archp$celltype2,
  #Infl_module = archp$mod_2,
  Infl_module = archp$AP1)
head (mod_df)
df_order = mod_df %>% 
group_by (celltype) %>% 
summarize (avg_module = median(Infl_module)) %>% 
arrange(avg_module)
mod_df$celltype = factor (mod_df$celltype, levels = rev(df_order$celltype))

bp = ggplot (mod_df, aes (x = celltype, y = Infl_module, fill=celltype)) +
vlp + 
bxpv + 
scale_fill_manual (values = palette_myeloid) +
#geom_point (position='identity', alpha=.3, color="grey44", size=1) +
gtheme

pdf (file.path ('Plots','celltype_infl_module_boxplots2.pdf'),2.3,width=4)
bp
dev.off()


# # Compute co-occurrence of TFs ####
infl_TF = names(km$cluster[km$cluster == 2])
metaGroupName = 'inflamed'
celltypes = unique (archp@cellColData[,metaGroupName])
motifMat = getPositions (archp)
matches = getMatches (archp)

# Find DAP ####
#force = FALSE
force=F
if (!file.exists (paste0('DAP_',metaGroupName,'.rds')) | force)
  {
  DAP_list = getMarkerFeatures (
    ArchRProj = archp, 
    testMethod = "wilcoxon",
          useGroups = celltypes[1],
          bgdGroups = celltypes[2],
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
rownames(DAP_res) = as.character(DAP_res_regions)
DAP_res = DAP_res[!is.na(DAP_res$FDR),]
DAP_res = DAP_res[DAP_res$FDR < 0.05,]
DAP_res_l = list(up = GRanges(rownames(DAP_res)[DAP_res$Log2FC > 0]),
                down = GRanges(rownames(DAP_res)[DAP_res$Log2FC < 0]))
sapply (DAP_res_l, length)
names (DAP_res_l) = celltypes

cooc_l = list()
ov_size_l=list()
for (celltype in celltypes)
  {
  matches_ct = matches[queryHits(findOverlaps (matches, DAP_res_l[[celltype]]))]
  matchesMat = assay (matches_ct)
  colnames (matchesMat) = gsub ('_.*','',colnames (matchesMat))
  colnames (matchesMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames (matchesMat))
  matchesMat = matchesMat[,infl_TF]
  matchesMat = matchesMat[rowSums (matchesMat) > 0,]
  
  cooc = matrix (ncol = ncol(matchesMat), nrow= ncol(matchesMat))
  ov_size = matrix (ncol = ncol(matchesMat), nrow= ncol(matchesMat))
  
  for (i in 1:ncol(matchesMat))
    {
    for (z in 1:ncol(matchesMat)) 
      {
      #if (sum(rowSums (matchesMat[,c(i,z)]) == 2) < 10)
     #   {
      #  cooc[i,z] = 0
      #  } else {
        ov = sum (rowSums (matchesMat[,c(i,z)]) == 2) / min (colSums(matchesMat[,c(i,z)]))
        cooc[i,z] = ov
        ov_size[i,z] = sum(rowSums (matchesMat[,c(i,z)]) == 2)
      #  }      
      }
    }
  
  colnames (cooc) = infl_TF
  rownames (cooc) = infl_TF
  diag (cooc) = 0
  cooc_l[[celltype]] = cooc#[rowSums(cooc) >0,rowSums(cooc) >0]
  colnames (ov_size) = infl_TF
  rownames (ov_size) = infl_TF
  ov_size_l[[celltype]] = ov_size
  }
#cooc_l[[celltype]][lower.tri (cooc_l[[celltype]])]
cooc_diff = cooc_l[[1]] - cooc_l[[2]]
cooc_diff[is.na(cooc_diff)] = 0
ov_size_max = pmin(ov_size_l[[1]], ov_size_l[[2]])
diag(ov_size_max) = 0

cooc_hm = Heatmap (
    cooc_diff,
    column_names_rot =45, 
    row_names_gp = gpar(fontsize = 3),
    column_names_gp = gpar(fontsize = 3),
    col =palette_cooccurrence_cor, 
    cluster_rows=T,
    cluster_columns = T#,
#rect_gp = gpar (col = "white", lwd = 1)
)

pdf (file.path ('Plots','selected_TF_cooccurence_heatmaps_diff.pdf'), width = 8,height=6)
cooc_hm 
dev.off()

## Try with scatterplot 
cooc_diff_df = as.data.frame (cooc_diff)
cooc_diff_df$TF2 = rownames (cooc_diff_df)
cooc_diff_df = gather (cooc_diff_df, TF, overlap, 1:(ncol(cooc_diff_df)- 1))

# overlap size
ov_size_max_df = as.data.frame(ov_size_max)
ov_size_max_df$TF2 = rownames(ov_size_max_df)
ov_size_max_df = gather (ov_size_max_df, TF, size, 1:(ncol(ov_size_max_df)- 1))

ov_diff_size_df = cbind (cooc_diff_df, ov_size_max_df)
ov_diff_size_df$TF_pair = paste0(ov_diff_size_df$TF2,'_', ov_diff_size_df$TF)
ov_diff_size_df = ov_diff_size_df[, c(3,6,7)]
ov_diff_size_df$topTFs = ifelse (ov_diff_size_df$overlap > 0.5 & ov_diff_size_df$size > 1000, ov_diff_size_df$TF_pair, '') 
sp = ggplot (ov_diff_size_df, aes (x = overlap, y = size, label = topTFs)) + 
geom_point() + 
geom_text_repel (size=2) +
gtheme_no_rot

pdf (file.path ('Plots','peak_overlap_inflamed_scatterplot.pdf'))
sp
dev.off()


# Compare inflamed vs non-inflamed Momacs ####
library (presto)
archp$inflamed = ifelse (archp$mod_2 > 0, 'inflamed','non_inflamed')
metaGroupName = 'inflamed'
all (colnames(hubsCell_mat) == rownames(archp@cellColData))
metagroup = as.character (archp@cellColData[,metaGroupName])
res = wilcoxauc (log2(hubsCell_mat+1), metagroup)
res = res[res$logFC > 0,]

res_l = lapply (split (res, res$group), function(x){
  tmp = x[order (x$padj),]
  tmp
})

res_df = do.call (rbind, res_l)
res_df$gene = hubs_obj$hubsCollapsed$gene[match(res_df$feature, hubs_obj$hubs_id)]

top_hubs = 20
res_df_top = res_df %>% group_by (group) %>%
  slice_head(n = top_hubs)

hub = res_df_top$feature
hub = hubs_obj$hubs_id[grep ('NFKB1', hubs_obj$hubsCollapsed$gene)]
#sample_levels = c('Monocytes','cDCs','SPP1','TREM2','C1Q','IFN','IM')

sample_levels = c('Monocytes','cDCs','SPP1','TREM2','C1Q','IFN','IM')

pdf()
#archp$fetal_sample = paste0(archp$Sample, archp$fetal_group)
#metaGroupName = 'fetal_group'
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp,#[!archp$Sample3 %in% c('P11_HOX')], 
    sample_levels = sample_levels, 
    hubs_regions = hubs_obj$hubsCollapsed,
    #ylim = c(0,0.30),
    groupBy = metaGroupName, 
    #sample_levels = sample_sarc_order,
    minCells = 10,
    #geneSymbol = TF,
    plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
    #pal = palette_sample,
    #pal = palette_fetal,
    threads=1,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    region = ext_range (GRanges (hubs_obj$hubsCollapsed[match(hub, hubs_obj$hubs_id)]),100000,100000),
    #upstream = 100000,
    #downstream = 100000,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    pal = palette_myeloid,
    #loops = getCoAccessibility (archp, corCutOff = 0.3,
    #  returnLoops = TRUE),
    useGroups= NULL
)
dev.off()
plotPDF (meso_markers, ArchRProj = archp,height=3.5, width=6, name =paste0('MPM_markers_inflamed_coveragePlots.pdf'),addDOC=F)
  


# Compare cNMF modules with inflammatory program in scATAC-seq and scRNA-seq ####
shared_cnmf = readRDS (file.path('..','scrna','shared_cnmf_myeloid.rds'))
shared_cnmf = lapply (shared_cnmf, function(x) x[x %in% getFeatures (archp)])
#remove_modules = c('cnmf.3','cnmf.6','cnmf.7','cnmf.5') # remove monocyres cDC and CC modules. Consider re-inculding CC 

pdf (file.path ('Plots','scrna_celltype_dimplot.pdf'))
DimPlot (srt, group.by = 'celltype', reduction = 'umap')
dev.off()

srt = ModScoreCor (
    seurat_obj = srt, 
    geneset_list = shared_cnmf, 
    cor_threshold = NULL, 
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'shared_cnmf', outdir = NULL)

TF_modules = split(names(km$cluster), km$cluster)
TF_modules = c('JUNB
FOSL2
JUN
SMARCC1
FOSL1
JUND
FOS
JDP2
BACH1
FOSB')
TF_modules = strsplit(TF_modules, '\n')
names (TF_modules) = 'AP1'


srt = ModScoreCor (
    seurat_obj = srt, 
    geneset_list = TF_modules, 
    cor_threshold = NULL, 
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'AP1', outdir = NULL)

if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat)#[selected_TF,])
mMat = mMat[TF_modules[[1]],]
mMat = colMeans (mMat)
all (colnames (mMat) == rownames(archp@cellColData))
archp$AP1 = mMat
rna_mod_cor = cor (srt@meta.data[,names(shared_cnmf)],srt@meta.data[,'AP1'])
atac_mod_cor = cor (archp@cellColData[,names(shared_cnmf)],mMat)
#cor.test (archp@cellColData[,names(shared_cnmf)],archp@cellColData[,'mod_2'])

hm2 = Heatmap (atac_mod_cor, border=T, col = rev (palette_deviation_correlation))
hm1 = Heatmap (rna_mod_cor, cluster_rows=F, border=T, col = rev (palette_deviation_correlation))
pdf (file.path ('Plots','rna_atac_module_cor.pdf'), width=2.5, height=3)
hm2 + hm1
dev.off()

# Make scatterplots of inflammation vs cnmfs ####
atac_mat = cbind(archp@cellColData[,names(shared_cnmf)], AP1 = archp$AP1)
atac_mat_long = gather (as.data.frame(atac_mat), cnmf, score, 1:(ncol(atac_mat)-1))
atac_mat_long$cnmf = factor (atac_mat_long$cnmf, levels = rownames(atac_mod_cor)[order(-atac_mod_cor[,1])])
library (ggpointdensity)
# remove outliers
atac_mat_longL = split (atac_mat_long, atac_mat_long$cnmf)
atac_mat_longL = lapply (atac_mat_longL, function(x) x[x$score > quantile(x$score,.005) & x$score < quantile(x$score,.995),])
atac_mat_long = do.call (rbind,atac_mat_longL)
sp = ggplot(atac_mat_long, aes(x = score, y = AP1)) +
  # geom_point(
  #   alpha = 0.3, 
  #   size = .1
  # ) + 
  facet_wrap (~cnmf, scales = 'free', ncol= nrow(atac_mod_cor)) +
   geom_pointdensity (alpha=1, size=.1) +
  scale_color_viridis (option='F') +
# Scatterplot points with transparency and size
  geom_smooth(
    method = "lm", 
    color = "white", 
    fill = "white", 
    se = FALSE, 
    linetype = "dashed", 
    size = .4
  ) + # Regression line with confidence interval
  theme_void() + # Clean theme
  labs(
    title = "AP1 vs cnmfs",
    #subtitle = "Scatterplot with Regression Line and Correlation Coefficient",
    x = "cnmf",
    y = "AP1")
  # ) +
  # theme(
  #   plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
  #   plot.subtitle = element_text(size = 14, hjust = 0.5),
  #   axis.title = element_text(size = 14),
  #   axis.text = element_text(size = 12),
  #   panel.grid.major = element_line(color = "grey90"),
  #   panel.grid.minor = element_blank()
  # ) 

png (file.path ('Plots','cnmf_AP1_scatterplots2.png'),width=3900, height=500, res=300)
sp
dev.off()

sp = sp + stat_cor (
    aes(label = paste(..rr.label.., ..p.label.., sep = " | ")), 
    method = "spearman", 
    #label.x = min(df$fetal) + 0.1 * diff(range(df$fetal)), 
    #label.y = max(df$fetal_gs) - 0.1 * diff(range(df$fetal_gs)),
    color = "grey11",
    size = 1
  )   
pdf (file.path ('Plots','cnmf_AP1_scatterplots2.pdf'),width=9, height=1.4)
sp
dev.off()


## Run peak2genes results with hubs links ####
run_p2g = T
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
    

# Import hubs from myeloid analysis ####
metaGroupName = "Clusters_H"
cor_cutoff = 0.3
#max_dist = 12500
max_dist = 12500
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
hubs_obj = readRDS (file.path(hubs_dir,'global_hubs_obj.rds'))


# Generate matrix of fragment counts of hubs x barcodes ####
if (!file.exists(file.path (hubs_dir, paste0('hubs_cells_mat.rds'))))
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

all (colnames(hubsCell_mat) == rownames(archp@cellColData))
#ha = HeatmapAnnotation (fetal = archp$fetal, which='row')
# hm = Heatmap (
#   scale (t(hubsCell_mat)), 
#  # left_annotation = ha, 
#   column_names_gp = gpar(fontsize = 3),
#   row_names_gp = gpar(fontsize = 0),
#   show_column_dend = T,
#   #column_km = 5,
#   #row_dend_width = unit(5,'mm'),
#   row_dend_side = 'left',
#   col = rev(palette_hubs_accessibility),
#   border=T,
#   name = 'Hubs')
# pdf (file.path (hubs_dir,'Plots',paste0('hubs_cells_',metaGroupName,'_heatmap.pdf')), height=2.2, width = 5)
# hm
# dev.off()



### Compute differential hub analysis 
# Compute differential hub accessibility DHA ####
library (presto)
metaGroupName = 'celltype2'
all (colnames(hubsCell_mat) == rownames(archp@cellColData))
metagroup = as.character (archp@cellColData[,metaGroupName])
res = wilcoxauc (log2(hubsCell_mat+1), metagroup)
res = res[res$logFC > 0,]

res_l = lapply (split (res, res$group), function(x){
  tmp = x[order (x$padj),]
  tmp
})

res_df = do.call (rbind, res_l)
res_df$gene = hubs_obj$hubsCollapsed$gene[match(res_df$feature, hubs_obj$hubs_id)]
head (res_df[res_df$group == 'TREM2',],20)

top_hubs = 5
res_df_top = res_df %>% group_by (group) %>%
  slice_head(n = top_hubs)

levels_order = unique(res_df_top$group)
metagroup = factor (metagroup, levels = levels_order, ordered=T)


# Generate matrix of fragment counts of hubs x metagroup ####
metaGroupName = 'celltype2'
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
  frags_in_sample = sapply (unique(archp@cellColData[,metaGroupName]), function(x) sum (archp$ReadsInTSS[as.character(archp@cellColData[,metaGroupName]) == x]))
  hubsSample_mat = t(t(hubsSample_mat) * (10^6 / frags_in_sample)) # scale
  saveRDS (hubsSample_mat, file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))
  } else {
  hubsSample_mat = readRDS (file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))  
  }
hubsSample_mat = as.data.frame (hubsSample_mat)

#as.character(archp@cellColData[,metaGroupName])[order(metagroup, na.last=T)]
#ha = HeatmapAnnotation (celltype = metagroup[order (metagroup, na.last=T)])
#hubs_mat
dah_hm = t(scale(t(log2(hubsSample_mat[res_df_top$feature,unique(res_df_top$group)]))))
#dah_hm[dah_hm > 2] = 2
#dah_hm[dah_hm < -2] = -2
hm = Heatmap (dah_hm,
  cluster_rows=F,
#  top_annotation = ha,
  cluster_columns=F,
  row_labels = res_df_top$gene,
  col = rev(palette_hubs_accessibility),
  column_names_gp= gpar (fontsize=11),
  row_names_gp= gpar (fontsize=9),
  column_names_rot=45,
  border=T)

pdf (file.path ('Plots','top_DAH_cnmf_celltypes_heatmap.pdf'), width=5)
hm
dev.off()


#TF = 'SNAI1'
TF = sapply (unique(res_df_top$gene), function(x) unlist(strsplit(x, '-'))[1])
metaGroupName = 'celltype2'

pdf()
#archp$fetal_sample = paste0(archp$Sample, archp$fetal_group)
#metaGroupName = 'fetal_group'
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp,#[!archp$Sample3 %in% c('P11_HOX')], 
    #group_order = sample_levels, 
    #ylim = c(0,0.30),
    groupBy = metaGroupName, 
    hubs_regions = hubs_obj$hubsCollapsed,
    #sample_levels = sample_sarc_order,
    minCells = 10,
    geneSymbol = TF,
    plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
    #pal = palette_sample,
    #pal = palette_fetal,
    threads=1,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 100000,
    downstream = 100000,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
    #loops = getCoAccessibility (archp, corCutOff = 0.3,
    #  returnLoops = TRUE),
    useGroups= NULL
)
dev.off()
plotPDF (meso_markers, ArchRProj = archp,height=3.5, width=6, name =paste0('MPM_markers_',TF,'_coveragePlots.pdf'),addDOC=F)
  


# Plot largest chubs ####  
TF = head (hubs_obj$hubsCollapsed,10)
metaGroupName = 'celltype2'
TF = 'NFKB1'
pdf()
#archp$fetal_sample = paste0(archp$Sample, archp$fetal_group)
#metaGroupName = 'fetal_group'
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp,#[!archp$Sample3 %in% c('P11_HOX')], 
    #group_order = sample_levels, 
    #ylim = c(0,0.30),
    groupBy = metaGroupName, 
    hubs_regions = hubs_obj$hubsCollapsed,
    #sample_levels = sample_sarc_order,
    minCells = 10,
    #geneSymbol = TF,
    region = ext_range(TF, 100000,100000),
    plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
    #pal = palette_sample,
    #pal = palette_fetal,
    threads=1,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 100000,
    downstream = 100000,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
    #loops = getCoAccessibility (archp, corCutOff = 0.3,
    #  returnLoops = TRUE),
    useGroups= NULL
)
dev.off()
plotPDF (meso_markers, ArchRProj = archp,height=3.5, width=6, name =paste0('MPM_markers_largest_cHubs_coveragePlots.pdf'),addDOC=F)



# Check footprint of RUNX and NR4A2 across celltypes ####
metaGroupName='celltype2'
archp <- addGroupCoverages (ArchRProj = archp, groupBy = metaGroupName)
motifPositions <- getPositions (archp)

motifs <- c('NFKB1','JUNB','FOS','JUND','SPI1','SPIB','MITF','RUNX1','CEBPA','SRF','NFAT5','STAT2','MEF2C','PRDM1','IRF3','IRF8','IRF1','IRF2','IRF9')
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

seFoot <- getFootprints(
  ArchRProj = archp, 
  #positions = motifPositions_sample[markerMotifs], 
  positions = motifPositions[markerMotifs], 
  groupBy = metaGroupName
)
  
plotFootprints(
seFoot = seFoot,
ArchRProj = archp, 
normMethod = "Subtract",
plotName = "Footprints-Subtract-Bias_",
addDOC = FALSE, height=7.5, width=5,
pal = palette_myeloid,
smoothWindow = 25)
  


# Check TF deviations
TF = 'E2F3'
getFeatures (archp, 'MotifMatrix')[grep (TF, getFeatures (archp, 'MotifMatrix'))]
TF1 = c('z:JUN_143','z:FOS_137','z:NFKB1_719','z:FOXM1_352','z:TFDP1_310','z:E2F3_313')

pdf ()
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "MotifMatrix",
    name = TF1, 
    useSeqnames='z',
    pal = rev (palette_deviation),    
    embedding = "UMAP_H",
    imputeWeights = NULL
    )
dev.off()
pdf(file.path('Plots','avg_AP1_deviation_fplot.pdf'),15,15)
wrap_plots(TF_p)
wrap_plots(umap_p1)
wrap_plots (umap_p0,umap_p2)#,umap_p3)
dev.off()





# ############ INTEGRATION WITH scRNA ######
# archp <- addGeneIntegrationMatrix (
#     ArchRProj = archp, 
#     useMatrix = "GeneScoreMatrix",
#     matrixName = "GeneIntegrationMatrix",
#     reducedDims = "IterativeLSI",
#     seRNA = srt,
#     addToArrow = TRUE,
#     groupRNA = "celltype2",
#     nameCell = "predictedCell_Un",
#     nameGroup = "predictedGroup_Un",
#     nameScore = "predictedScore_Un",
#     force = TRUE
# )

# #rna_p = DimPlot (archp, group.by = 'celltype')
# pdf()
# int_p <- plotEmbedding(
#     archp, 
#     colorBy = "cellColData", 
#     name = "predictedGroup_Un", 
#     embedding = 'UMAP_H'
# #    pal = palD
# )
# int_p2 <- plotEmbedding(
#     archp, 
#     colorBy = "cellColData", 
#     name = "predictedGroup_Un", 
#     highlightCells = rownames(archp@cellColData)[archp$predictedGroup_Un == 'Mono_CD14'],
#     embedding = 'UMAP_H'
# #    pal = palD
# )
# dev.off()
# pdf (paste0 (projdir,'/Plots/RNA_integration_UMAPs.pdf'), width=15)
# int_p
# int_p2
# dev.off()



# ## Add RNA integration annotation for monocytes 
# archp_MAC = archp[!archp$predictedGroup_Un %in% c('DC2','Mono_CD14','Mono_CD16')]

# # Check if monocytes are better defined with label transfer
# TF = c('S100A9','EREG','CCL2')
# metaGroupName='predictedGroup_Un'
# pdf()
# #archp$fetal_sample = paste0(archp$Sample, archp$fetal_group)
# #metaGroupName = 'fetal_group'
# meso_markers <- plotBrowserTrack2 (
#     ArchRProj = archp,#[!archp$Sample3 %in% c('P11_HOX')], 
#     #group_order = sample_levels, 
#     #ylim = c(0,0.30),
#     groupBy = metaGroupName, 
#     #sample_levels = sample_sarc_order,
#     minCells = 10,
#     geneSymbol = TF,
#     plotSummary = c("bulkTrack", "featureTrack", 
#         "loopTrack","geneTrack", 
#         "hubTrack",'hubregiontrack'),
#     #pal = palette_sample,
#     #pal = palette_fetal,
#     threads=1,
#     #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
#     #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
#     upstream = 100000,
#     downstream = 100000,
#     loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
#     #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
#     #loops = getCoAccessibility (archp, corCutOff = 0.3,
#     #  returnLoops = TRUE),
#     useGroups= NULL
# )
# dev.off()
# plotPDF (meso_markers, ArchRProj = archp,height=3.5, width=6, name =paste0('Monocytes_markers_coveragePlots.pdf'),addDOC=F)


# # Generate heatmap of cnmf modules across cells and use kmeans to cluster ####
# shared_cnmf = readRDS (file.path('..','scrna','shared_cnmf_myeloid.rds'))
# shared_cnmf = lapply (shared_cnmf, function(x) x[x %in% getFeatures (archp)])

# remove_modules = c('cnmf.3','cnmf.6') # remove monocyres and cDC modules
# shared_cnmf_MAC = shared_cnmf[!names(shared_cnmf) %in% remove_modules]
# cnmf_scatac = as.data.frame (t(scale(t(archp_MAC@cellColData[,names(shared_cnmf_MAC)]))))

# set.seed (123)
# cnmf_remove = c('cnmf.7') # remove redundant cell cycle module 
# km = kmeans (scale(cnmf_scatac[,!colnames(cnmf_scatac) %in% cnmf_remove]), centers=6)  

# ha = HeatmapAnnotation (sample = archp_MAC$Sample, col=list(sample = palette_sample))
# hm = Heatmap (t(scale(cnmf_scatac[,!colnames(cnmf_scatac) %in% cnmf_remove])), 
#   col = palette_genescore_fun(scale(cnmf_scatac)), 
#   top_annotation = ha,
#   show_column_dend = F,
#   column_split = km$cluster,
# #  column_km=9,
#   row_names_gp = gpar (fontsize = 8),
#   column_names_gp = gpar (fontsize = 0),
#   border=T)

# pdf (file.path ('Plots','cnmf_scatac_barcodes_heatmap_scaled_only_MAC2.pdf'), height=2.5)
# hm
# dev.off()

# archp$cnmf_celltypes2 = archp$predictedGroup_Un
# archp$cnmf_celltypes2[match(names(km$cluster), rownames(archp@cellColData))] = paste0('cnmf_',km$cluster)
# archp$cnmf_celltypes2[archp$cnmf_celltypes2 =='cnmf_1'] = 'IM'
# archp$cnmf_celltypes2[archp$cnmf_celltypes2 =='cnmf_2'] = 'IL1B'
# archp$cnmf_celltypes2[archp$cnmf_celltypes2 =='cnmf_3'] = 'TREM2'
# archp$cnmf_celltypes2[archp$cnmf_celltypes2 =='cnmf_4'] = 'SPP1'
# archp$cnmf_celltypes2[archp$cnmf_celltypes2 =='cnmf_5'] = 'IFN'
# archp$cnmf_celltypes2[archp$cnmf_celltypes2 =='cnmf_6'] = 'C1QA'







# Export bigiwg files ####
metaGroupName = 'cnmf_celltypes'
exp_bigwig = TRUE
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



























## Add column on DAM heatmap showing if TF is pioneer or not from chrombpnet ####
## Show barplots of top TF occurrence using finemo chrombpnet outputs ####

### Compare TF expression from scRNA and inferred by chrombpnet per cell type ####
# library (httr)
# library (XML)
# library (igraph)
#BiocManager::install("universalmotif")
library ('universalmotif')

metaGroupName = 'inflamed'
if (!any (ls() == 'mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
#mMat_mg = mMat[DAM_df$gene, ]
# mMat_mg = as.data.frame (t(mMat))
# mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
# mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
# rownames (mMat_mg) = mMat_mg[,1]
# mMat_mg = mMat_mg[,-1]


# #Get active genes from RNA
# metaGroupName = 'celltype_simplified2'
# ps = log2(as.data.frame (AverageExpression (srt, 
# features = colnames(mMat_mg),
# group.by = metaGroupName)[[1]]) +1)
# min_exp = .1
# #ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
# #active_TFs = rownames(ps)[rowSums(ps) > 0]

# #active_genes = corGSM_MM$MotifMatrix_name[corGSM_MM$cor > 0.1]
# #DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% active_TFs,])    
# mMat_l = as.list (as.data.frame (t(mMat_mg)))
# mMat_l = lapply (mMat_l, function(x) data.frame (dev = x, row.names = colnames(mMat_mg)))
# #mMat_l = lapply (mMat_l, function(x) x[rownames(x) %in% active_TFs,,drop=F])

# metaGroupName = 'celltype_lv1'
chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/chromBPnet'
# metaGroupName = 'celltype_lv1'
# celltypes = unique (archp@cellColData[,metaGroupName])

# tf_database = read_meme('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/HOCOMOCO_db/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme', skip = 0, readsites = FALSE, readsites.meta = FALSE)
# tf_database = unique(unlist(lapply(tf_database, function(x) unlist(strsplit(x@name,'_'))[1])))

# list.files (file.path(chromBPdir, celltypes[3],'no_bias_model'))


chrombpnet_counts = list()
metaGroupName = 'inflamed'
celltypes = unique (archp@cellColData[,metaGroupName])
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_counts[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
  }

# Check overlap NFKB with AP-1 complex
motif_pairs_l = list(
  c('FOS','NFKB1'),
  c('JUN','NFKB1'),
  c('SNAI1','NFKB1'),
  c('ATF1','NFKB1'),
  c('CTCF','NFKB1'),
  c('SPI1','NFKB1'),
  c('KLF12','NFKB1'),
  c('CEBPA','NFKB1'),
  c('RUNX1','NFKB1'),
  c('GABPA','NFKB1'),
  c('ZBT7A','NFKB1'))
peak_overlap = NULL
motif_pairs=list()
for (i in seq_along(motif_pairs_l))
  {
  for (celltype in celltypes)
    {
    peakset = read.table (file.path(chromBPdir,paste0('MACS2_',celltype),paste0(celltype, '_peaks.narrowPeak')))
    colnames(peakset) = c('chr','start','end')
    peakset = makeGRangesFromDataFrame (peakset)
    chrombpnet_counts_gr1 = chrombpnet_counts[[celltype]][grep (motif_pairs_l[[i]][1], chrombpnet_counts[[celltype]]$V4),]
    chrombpnet_counts_gr2 = chrombpnet_counts[[celltype]][grep (motif_pairs_l[[i]][2], chrombpnet_counts[[celltype]]$V4),]
    colnames(chrombpnet_counts_gr1) = c('chr','start','end')
    colnames(chrombpnet_counts_gr2) = c('chr','start','end')
    chrombpnet_counts_gr1 = makeGRangesFromDataFrame (chrombpnet_counts_gr1)
    chrombpnet_peaks_gr1 = peakset[unique(queryHits(findOverlaps(peakset, chrombpnet_counts_gr1)))]
    chrombpnet_counts_gr2 = makeGRangesFromDataFrame (chrombpnet_counts_gr2)
    chrombpnet_peaks_gr2 = peakset[unique(queryHits(findOverlaps(peakset, chrombpnet_counts_gr2)))]
  
    peak_overlap[[celltype]] = sum (countOverlaps (chrombpnet_peaks_gr1,chrombpnet_peaks_gr2) > 0) / 
    min (c(length(chrombpnet_peaks_gr1),length(chrombpnet_peaks_gr2)))
    }
  motif_pairs[[i]] = peak_overlap
  }
 motif_pairs_l2 = unlist(lapply (motif_pairs_l, function(x) paste (x, collapse='_')))
 names (motif_pairs) = motif_pairs_l2
#motif_pairs = unlist (motif_pairs_l, recursive=F)




  

# Intersect DAP with MACS2 peaks and look at chrombpnet predicted TFs inflamed vs not ####
DAP_list = readRDS (paste0('DAP_',metaGroupName,'.rds'))
DAP_res = do.call (cbind, (assays(DAP_list)))
colnames (DAP_res) = names(assays(DAP_list))
DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames(DAP_res) = as.character(DAP_res_regions)
DAP_res = DAP_res[!is.na(DAP_res$FDR),]
DAP_res = DAP_res[DAP_res$FDR < 0.05,]
DAP_res_l = list(up = GRanges(rownames(DAP_res)[DAP_res$Log2FC > 0]),
                down = GRanges(rownames(DAP_res)[DAP_res$Log2FC < 0]))
celltypes = c('inflamed','non_inflamed')
names (DAP_res_l) = celltypes

chrombpnet_tfs_DAP_l = list()
ov_mat_l = list()
for (celltype in celltypes)
  {
  peakset = read.table (file.path(chromBPdir,paste0('MACS2_',celltype),paste0(celltype, '_peaks.narrowPeak')))
  colnames(peakset) = c('chr','start','end')
  peakset = makeGRangesFromDataFrame (peakset)
  peakset = peakset[unique(queryHits(findOverlaps(peakset, DAP_res_l[[celltype]])))]
  chrombpnet_counts_gr = chrombpnet_counts[[celltype]]
  colnames(chrombpnet_counts_gr) = c('chr','start','end','TF')
  chrombpnet_counts_gr = makeGRangesFromDataFrame (chrombpnet_counts_gr, keep.extra.columns=T)
  chrombpnet_counts_gr = chrombpnet_counts_gr[queryHits(findOverlaps(chrombpnet_counts_gr, peakset))]
  
  peakset_l = split (peakset, 1:length(peakset))
  ov_mat = sapply (unique(chrombpnet_counts_gr$TF), function(y) countOverlaps (peakset, chrombpnet_counts_gr[chrombpnet_counts_gr$TF == y]))
  rownames(ov_mat) = as.character(peakset)
  ov_mat_l[[celltype]] = ov_mat
  chrombpnet_tfs_DAP_l[[celltype]] = chrombpnet_counts_gr

  }
table (chrombpnet_tfs_DAP_l[[1]]$TF)[order(-table (chrombpnet_tfs_DAP_l[[1]]$TF))]
table (chrombpnet_tfs_DAP_l[[2]]$TF)[order(-table (chrombpnet_tfs_DAP_l[[2]]$TF))]

ov_mat_cor = lapply (ov_mat_l, function(x) cor (x))

# Filter using only TF from inflammation module
infl_TF
ov_mat_cor_fl = lapply (ov_mat_cor, function(x) {
  tfmatch = unique(unlist(sapply (infl_TF, function(y) (grep(y, rownames(x))))))
  x[tfmatch, tfmatch]
})
pdf (file.path ('Plots','chrombpnet_TF_cor_in_peaks.pdf'),height=4,width=4.5)
lapply (ov_mat_cor_fl, function (x) Heatmap (
  x, col = palette_cooccurrence_cor_fun,
  row_names_gp = gpar (fontsize = 6),
  column_names_gp = gpar (fontsize = 6), 
border=T))
dev.off()

# Check expression of FOS JUNB and JUND
genes = c('FOS','JUNB','JUND','NFKB1','NFKB2','SPI1','SPIB','CEBPA','SRF','MITF')
pdf (file.path ('Plots','expression_FOS_JUNB_JUND.pdf'))
DotPlot (srt, features = genes, group.by = 'celltype2') + gtheme
dev.off()






chrombnet_counts_2 = list()
for (celltype in celltypes)
  {
  chrombpnet_counts_tmp = chrombpnet_counts[[celltype]]
  chrombnet_counts_2[[celltype]] = table (chrombpnet_counts_tmp$V4)[order(-table (  chrombpnet_counts_tmp$V4))]
  }


chrombpnet_profile = list()
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  }
  #chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  #chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  #chrombpnet_tf = rbind (chrombpnet_count_tf, chrombpnet_profile_tf)

#ap1_complex = c('JUN','FOSB','FOS','BACH1','SMARCC1','FOSL2','JUND','JDP2','BATF','CEBPB','CEBPA','CEBPZ','FOSL1','NFE2','NFE2L2','NFE2L1')

chrombpnet_profile_2 = list()
for (celltype in celltypes)
  {
  chrombpnet_profile_tmp = chrombpnet_profile[[celltype]]
  chrombpnet_profile_2[[celltype]] = table (chrombpnet_profile_tmp$V4)[order(-table (  chrombpnet_profile_tmp$V4))]
  }  
  # #assign_max_exp = unlist(sapply (names(chrombpnet_counts_tmp), function(x) unlist(strsplit(x, '_'))[which.max(ps[unlist(strsplit(x, '_')), celltype])]))
  # #tf_dev = mMat_l[[celltype]][assign_max_exp,]
  # chrombpnet_counts_tmp2 = data.frame (occurrence =   chrombpnet_counts_tmp[names(assign_max_exp)], TF_max_exp = assign_max_exp, TF_max_dev = tf_dev)
  # #chrombpnet_counts_tmp2 = chrombpnet_counts_tmp2[!chrombpnet_counts_tmp2$TF_max_exp %in% ap1_complex,]
  # chrombpnet_counts_tmp2 = chrombpnet_counts_tmp2[!duplicated(chrombpnet_counts_tmp2$TF_max_exp),]
  # chrombpnet_counts_tmp2 = chrombpnet_counts_tmp2[chrombpnet_counts_tmp2$TF_max_exp %in% head(chrombpnet_counts_tmp2$TF_max_exp[order(-chrombpnet_counts_tmp2$TF_max_dev)],10),]
  # chrombpnet_counts_tmp2$celltype = celltype
  # chrombpnet_counts_tmp2$order = seq(nrow(chrombpnet_counts_tmp2))
  # chrombpnet_counts2[[celltype]] = chrombpnet_counts_tmp2
  # }

chrombpnet_counts_df = do.call (rbind, chrombpnet_counts2)
chrombpnet_counts_df = chrombpnet_counts_df %>% group_by (celltype) %>% mutate(Proportion = occurrence.Freq / sum(occurrence.Freq))
chrombpnet_counts_df$TF_max_exp2 = chrombpnet_counts_df$TF_max_exp
chrombpnet_counts_df$TF_max_exp[chrombpnet_counts_df$Proportion < 0.05] = ''
chrombpnet_counts_df$TF_max_exp = factor (chrombpnet_counts_df$TF_max_exp, levels =unique(chrombpnet_counts_df$TF_max_exp))
chrombpnet_counts_df$order = factor (chrombpnet_counts_df$order, levels =unique(chrombpnet_counts_df$order))
# Create stacked bar plot with text beside each band
bp = ggplot (chrombpnet_counts_df, aes(x = celltype, y = Proportion, fill = order)) +
  geom_bar (stat = "identity", color = 'white') +
  geom_text (aes(label = TF_max_exp), 
            position = position_stack (vjust = 0.5), 
            hjust = 0.5,  # Move text outside the bar
            size = 3) + 
  #coord_flip() +  # Flip to make text more readable
  gtheme
pdf (file.path ('Plots','chrombpnet_counts_TF_barplot.pdf'), width=7, height=4)
bp
dev.off()


chrombpnet_profile = list()
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  }
  #chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  #chrombpnet_profile_tf = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  #chrombpnet_tf = rbind (chrombpnet_count_tf, chrombpnet_profile_tf)
chrombpnet_profile2 = list()
for (celltype in celltypes)
  {
  chrombpnet_profile_tmp = chrombpnet_profile[[celltype]]
  chrombpnet_profile_tmp = table (  chrombpnet_profile_tmp$V4)[order(-table (  chrombpnet_profile_tmp$V4))]
  assign_max_exp = unlist(sapply (names(chrombpnet_profile_tmp), function(x) unlist(strsplit(x, '_'))[which.max(ps[unlist(strsplit(x, '_')), celltype])]))
  tf_dev = mMat_l[[celltype]][assign_max_exp,]
  chrombpnet_profile_tmp2 = data.frame (occurrence =   chrombpnet_profile_tmp[names(assign_max_exp)], TF_max_exp = assign_max_exp, TF_max_dev = tf_dev)
  #chrombpnet_profile_tmp2 = chrombpnet_profile_tmp2[!chrombpnet_profile_tmp2$TF_max_exp %in% ap1_complex,]
  chrombpnet_profile_tmp2 = chrombpnet_profile_tmp2[!duplicated(chrombpnet_profile_tmp2$TF_max_exp),]
  chrombpnet_profile_tmp2 = chrombpnet_profile_tmp2[chrombpnet_profile_tmp2$TF_max_exp %in% head(chrombpnet_profile_tmp2$TF_max_exp[order(-chrombpnet_profile_tmp2$TF_max_dev)],10),]
  chrombpnet_profile_tmp2$celltype = celltype
  chrombpnet_profile_tmp2$order = seq(nrow(chrombpnet_profile_tmp2))
  chrombpnet_profile2[[celltype]] = chrombpnet_profile_tmp2
  }

chrombpnet_profile_df = do.call (rbind, chrombpnet_profile2)
chrombpnet_profile_df = chrombpnet_profile_df %>% group_by (celltype) %>% mutate(Proportion = occurrence.Freq / sum(occurrence.Freq))
chrombpnet_profile_df$TF_max_exp2 = chrombpnet_profile_df$TF_max_exp
chrombpnet_profile_df$TF_max_exp[chrombpnet_profile_df$Proportion < 0.05] = ''
chrombpnet_profile_df$TF_max_exp = factor (chrombpnet_profile_df$TF_max_exp, levels =unique(chrombpnet_profile_df$TF_max_exp))
chrombpnet_profile_df$order = factor (chrombpnet_profile_df$order, levels =unique(chrombpnet_profile_df$order))
# Create stacked bar plot with text beside each band
bp = ggplot (chrombpnet_profile_df, aes(x = celltype, y = Proportion, fill = order)) +
  geom_bar (stat = "identity", color = 'white') +
  geom_text (aes(label = TF_max_exp), 
            position = position_stack (vjust = 0.5), 
            hjust = 0.5,  # Move text outside the bar
            size = 3) + 
  #coord_flip() +  # Flip to make text more readable
  gtheme
pdf (file.path ('Plots','chrombpnet_profile_TF_barplot.pdf'), width=7, height=4)
bp
dev.off()


