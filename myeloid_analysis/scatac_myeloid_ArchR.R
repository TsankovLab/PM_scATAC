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
archp$celltype2[archp$Clusters_H == 'C6'] = 'Monocytes'
archp$celltype2[archp$Clusters_H %in% c('C7','C9')] = 'Inflammatory'
archp$celltype2[archp$Clusters_H %in% c('C8')] = 'DCs'
archp$celltype2[archp$Clusters_H %in% c('C10','C11','C3')] = 'TREM2_SPP1'
archp$celltype2[archp$Clusters_H %in% c('C4','C5','C2','C1')] = 'Interstitial'


pdf()
umap_p3 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
  pal= palette_sample,
   embedding = "UMAP_H")
umap_p4 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype2",
   embedding = "UMAP_H")
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters_H",
   embedding = "UMAP_H")
dev.off()

  pdf (file.path('Plots','celltype_umap_harmony_sample_umap.pdf'),5,5)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()

# Check for doublets ####
meso_markers = c('C1QA','SPP1','APOE','IL1B','CD3D','TREM2','C1QB','C3')
#archp = addImputeWeights (archp)
pdf()
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = meso_markers, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = NULL
)
dev.off()
pdf (file.path('Plots','myeloid_markers_fplots.pdf'), width = 18, height = 15)
wrap_plots (p, ncol=3)
dev.off()



# Get markers for gene score ####
immune_markers = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/scRNA_immune_markers_humanLUAD_Samarth_Assaf.csv')
immune_markers = immune_markers [immune_markers$group %in% c('Neutrophil','TRMac','IM','DC2','DC1','pDC','mregDC','CD14 mono',
  'CD16 mono','NK','Mast cell','Mgk','B/Plasma',' T cell','Treg','MoMac'),]
#immune_markers = immune_markers[immune_markers$group %in% c('CD14 mono','CD16 mono','DC1','DC2','MoMac'),]
immune_markers = immune_markers$gene
immune_markers = immune_markers[!immune_markers %in% c('CD14 MONO','IHBA','SEPP1','IL3RA')]
#archp = addImputeWeights (archp)
pdf()
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = immune_markers, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
png (file.path('Plots','myeloid_markers_fplots.png'), width = 18000, height = 15000, res=300)
wrap_plots (p)
dev.off()


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
  


# Check correlation of cnmf modules in scrna and scatac-seq ####
shared_cnmf = readRDS (file.path('..','scrna','shared_cnmf_myeloid.rds'))
shared_cnmf = lapply (shared_cnmf, function(x) x[x %in% getFeatures (archp)])
#srt_tam = srt[,srt$celltype2 == 'TAMs']
remove_samples = c('P3','P8', 'P4')
sample_names = unique(archp$Sample)#[unique(archp$Sample) %in% unique(srt$sampleID)]
sample_names_rna = unique(archp$Sample)[unique(archp$Sample) %in% unique(srt$sampleID)]
sample_names = sample_names[! sample_names %in% remove_samples]
sample_names_rna = sample_names_rna[! sample_names_rna %in% remove_samples]

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


# Generate heatmap of cnmf modules across cells and use kmeans to cluster ####
cnmf_scatac = as.data.frame (t(scale(t(archp@cellColData[,names(shared_cnmf)]))))

set.seed (123)
cnmf_remove = c('cnmf.7') # remove redundant cell cycle module 
km = kmeans (scale(cnmf_scatac[,!colnames(cnmf_scatac) %in% cnmf_remove]), centers=9)  

ha = HeatmapAnnotation (sample = archp$Sample, col=list(sample = palette_sample))
hm = Heatmap (t(scale(cnmf_scatac[,!colnames(cnmf_scatac) %in% cnmf_remove])), 
  col = palette_genescore_fun(scale(cnmf_scatac)), 
  top_annotation = ha,
  show_column_dend = F,
  column_split = km$cluster,
#  column_km=9,
  row_names_gp = gpar (fontsize = 8),
  column_names_gp = gpar (fontsize = 0),
  border=T)

pdf (file.path ('Plots','cnmf_scatac_barcodes_heatmap_scaled3.pdf'), height=3)
hm
dev.off()

# Re-Annotate based on cnmf clustering ####
all (names(km$cluster) == rownames(archp@cellColData))
archp$cnmf_cluster2 = paste0('cnmf_cluster_',km$cluster)
archp$cnmf_celltypes = archp$cnmf_cluster2
archp$cnmf_celltypes[archp$cnmf_cluster2 == 'cnmf_cluster_1'] = 'cDCs'
archp$cnmf_celltypes[archp$cnmf_cluster2 == 'cnmf_cluster_2'] = 'IM'
archp$cnmf_celltypes[archp$cnmf_cluster2 == 'cnmf_cluster_3'] = 'Mono'
archp$cnmf_celltypes[archp$cnmf_cluster2 == 'cnmf_cluster_4'] = 'C1Q'
archp$cnmf_celltypes[archp$cnmf_cluster2 == 'cnmf_cluster_5'] = 'SPP1'
archp$cnmf_celltypes[archp$cnmf_cluster2 == 'cnmf_cluster_6'] = 'TREM2'
archp$cnmf_celltypes[archp$cnmf_cluster2 == 'cnmf_cluster_7'] = 'IL1B'
archp$cnmf_celltypes[archp$cnmf_cluster2 == 'cnmf_cluster_8'] = 'IFN'
archp$cnmf_celltypes[archp$cnmf_cluster2 == 'cnmf_cluster_9'] = 'CC'

write.csv (data.frame (barcode = rownames(archp@cellColData), celltype = archp$cnmf_celltypes), 'barcode_annotation.csv')

# Find DAM in cnmf_clusters ####
#metaGroupName = "celltype2"
metaGroupName = "cnmf_celltypes"
#archp$celltype2 = archp$cnmf_cluster
force=FALSE
source (file.path('..','..','git_repo','utils','DAM.R'))

mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

metaGroupName = 'celltype2'
ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat), group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
#ps = ps[, colnames(DAM_hm@matrix)]
min_exp = .1
ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
active_TFs = rownames(ps)[rowSums(ps) > min_exp]

metaGroupName = 'cnmf_celltypes'
mMat_mg = mMat[active_DAM[active_DAM%in%active_TFs], ]
mMat_mg = as.data.frame (t(mMat_mg))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
#mMat_mg = mMat_mg[names (DAM_list),]

#selected_TF = c(rownames(DAM_hm@matrix), 'NR4A3','NR4A2','NR4A1')

DAM_hm = Heatmap (t(scale(mMat_mg)), 
          row_labels = colnames (mMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = T,
          #col = pals_heatmap[[5]],
          cluster_columns=T,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          name = 'chromVAR',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = rev(palette_deviation)#,
          #width = unit(2, "cm")
          #right_annotation = motif_ha
          )

scaled_ps = t(scale(t(ps[colnames(mMat_mg),])))
scaled_ps[is.na(scaled_ps)] = 0
TF_exp_selected_hm = Heatmap (scaled_ps[colnames(mMat_mg),],
        #right_annotation=tf_mark,
        #column_split = column_split_rna,
        cluster_rows = F, #km = 4, 
        name = 'expression (scaled)',
        column_gap = unit(.5, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=F, 
        col = palette_expression,
        row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
        column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
        border=T,
        width = unit(2, "cm"))

TF_exp_selected_hm2 = Heatmap (ps[colnames(mMat_mg),],
        #right_annotation=tf_mark,
        #column_split = column_split_rna,
        cluster_rows = F, #km = 4, 
        name = 'expression',
        column_gap = unit(.5, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=F, 
        col = palette_expression,
        row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
        column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
        border=T,
        width = unit(2, "cm"))

pdf (file.path ('Plots','cnmf_clusters_DAM_heatmap2.pdf'), width = 3,height=4)
draw (DAM_hm)
dev.off()

# Find DAG in cnmf_clusters ####
metaGroupName = "cnmf_celltypes"
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

  


archp = addImputeWeights (archp)
pdf()
p <- plotEmbedding (
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = names (shared_cnmf), 
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path ('Plots','shared_cnmf_fplots.pdf'),16,16)
wrap_plots (p, ncol=3)
dev.off()
  
srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = shared_cnmf, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'cNMF_', outdir = paste0(projdir,'Plots/'))
cnmf_scrna = srt@meta.data[,c(names(shared_cnmf))]
#cnmf_scrna = t(scale (t(cnmf_scrna)))
#colnames (cnmf_scrna) = sapply (colnames(cnmf_scrna), function(x) unlist(strsplit(x, '\\.'))[2])
cnmf_scrna = lapply (sample_names_rna, function(x) cnmf_scrna[srt$sampleID == x, ])
names (cnmf_scrna) = sample_names_rna
cnmf_scrna_celltype = lapply (sample_names_rna, function(x) srt$celltype2[srt$sampleID == x])
names (cnmf_scrna_celltype) = sample_names_rna
#cnmf_scatac = archp@cellColData[,names(shared_cnmf)]
cnmf_scatac = as.data.frame (t(scale(t(archp@cellColData[,names(shared_cnmf)]))))
cnmf_scatac = lapply (sample_names, function(x) cnmf_scatac[archp$Sample == x,])
names (cnmf_scatac) = sample_names

# Take median of module module correlation across samples ####
#sample_names_rna = sample_names[sample_names %in% unique(srt$sampleID)]
scrna_tf_cor = lapply (sample_names_rna, function(x) cor  (cnmf_scrna[[x]], method = 'spearman'))

sample_array <- simplify2array (scrna_tf_cor)
#any(lapply(corTF_array, function(x) any(is.na(x))))
# Take element-wise median
scrna_tf_cor <- apply (sample_array, c(1, 2), median)

pdf(file.path ('Plots','scrna_module_correlation_median_heatmap2.pdf'),width=5,height=4)
Heatmap (scrna_tf_cor, col = palette_expression_cor_fun(scrna_tf_cor), border=T)
dev.off()

scatac_tf_cor = lapply (sample_names, function(x) cor  (cnmf_scatac[[x]], method = 'spearman'))

sample_array <- simplify2array (scatac_tf_cor)
#any(lapply(corTF_array, function(x) any(is.na(x))))
# Take element-wise median
scatac_tf_cor <- apply (sample_array, c(1, 2), median)

pdf(file.path ('Plots','scatac_module_correlation_median_heatmap2.pdf'),width=5,height=4)
Heatmap (scatac_tf_cor, col = palette_deviation_cor_fun, border=T)
dev.off()


# Correlate module scores with TFs ####
metaGroupName = "Clusters_H"
archp$celltype2 = archp$Clusters_H
force=FALSE
source (file.path('..','..','git_repo','utils','DAM.R'))

metaGroupName = "sampleID"
#Get active genes from RNA
ps = log2(as.data.frame (AverageExpression (srt,
features = sapply (unique(unlist(lapply(DAM_list, function(x) x$gene))), function(x) unlist(strsplit (x, '_'))[1]), 
group.by = metaGroupName)[[1]]) +1)
min_exp = 0.5
ps = ps[apply(ps, 1, function(x) all (x > min_exp)),]
active_TFs = rownames(ps)[rowSums(ps) > 0]
write.csv (active_TFs, 'active_TFs.csv')
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.data.frame (t(as.matrix(scale(mMat[active_TFs,]))))
cnmf_scatac = as.data.frame (archp@cellColData[,names (shared_cnmf)])
cnmf_scatac = t(scale (t(cnmf_scatac)))


#colnames (cnmfs) = sapply (colnames(cnmfs), function(x) unlist(strsplit(x, '\\.'))[2])
cnmf_scatac_cor = lapply (sample_names, function(x) cor (cnmf_scatac[archp$Sample ==x ,], mMat[archp$Sample == x,], method = 'spearman'))


sample_array <- simplify2array(cnmf_scatac_cor)
#any(lapply(corTF_array, function(x) any(is.na(x))))
# Take element-wise median
cnmf_scatac_cor <- apply(sample_array, c(1, 2), median)

top_TFs = colnames(cnmf_scatac_cor)[order(-cnmf_scatac_cor['cnmf.8',])]
top_TFs = c(head (top_TFs,50), tail (top_TFs,50))

# Run TF nmf module correlation in scrna metacells #### 
metacells = readRDS (file.path ('..','scrna','metacells.rds'))
shared_cnmf = readRDS (file.path('..','scrna','shared_cnmf_myeloid.rds'))
metacells = ModScoreCor (
        seurat_obj = metacells, 
        geneset_list = shared_cnmf, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'cNMF_', outdir = paste0(projdir,'Plots/'))

cnmf_metacells = metacells@meta.data[,c(names(shared_cnmf))]
metacells_assay = metacells@assays$RNA@layers$data
rownames (metacells_assay) = rownames (metacells)
metacells_assay = t(metacells_assay[active_TFs,])
#sample_names_rna = sample_names[sample_names %in% unique(metacells$sampleID)]
#sample_names_rna = sample_names_rna[! sample_names_rna %in% remove_samples]
cnmf_scrna_cor = lapply (sample_names_rna, function(x) cor  (t(scale(t(metacells_assay[metacells$sampleID == x,]))), t(scale(t(cnmf_metacells[metacells$sampleID == x,]))), method = 'spearman'))

sample_array <- simplify2array(cnmf_scrna_cor)
#any(lapply(corTF_array, function(x) any(is.na(x))))
# Take element-wise median
cnmf_scrna_cor <- apply(sample_array, c(1, 2), median)

ha = rowAnnotation(' ' = anno_mark(at = match(top_TFs, colnames(cnmf_scatac_cor)), 
    labels = top_TFs, labels_gp = gpar(fontsize = 6, fontface='italic')))
hm1 = Heatmap (t(cnmf_scatac_cor),
  column_dend_side = 'bottom',
  #top_annotation = ha,
  row_names_gp = gpar (fontsize = 0),
  col = rev(palette_deviation),
#  row_km = 2,
  #rect_gp = gpar(type = "none"),
  border=T,
  column_names_gp = gpar (fontsize = 6, fontface='italic'))
hm2 = Heatmap (cnmf_scrna_cor,
  column_dend_side = 'bottom',
  #right_annotation = ha,
  row_names_gp = gpar (fontsize = 5),
  col = rev(palette_expression_correlation),
  #column_km = 2,
  #rect_gp = gpar(type = "none"),
  border=T,
  column_names_gp = gpar(fontsize = 8, fontface='italic'))

pdf (paste0 ('Plots/active_TF_cnmf_cor_heatmap.pdf'), width = 4,height=6)
hm1 + hm2
hm2
dev.off()

hm3 = Heatmap (cor (cnmf_scrna_cor, t(cnmf_scatac_cor)), cluster_rows=F, cluster_columns=F, border=T)
pdf (paste0 ('Plots/scrna_scatac_cor_cor_heatmap.pdf'), width = 5,height=4)
hm3
dev.off()


trem2_mod = cnmf_scrna_cor[,'cnmf.3']
trem2_mod = trem2_mod[order (-trem2_mod)]
trem2_mod [grep ('NFKB',names(trem2_mod))]
head (trem2_mod, 30)

# Show TF correlation to TREM2 module in scatac and rna 
module = 'cnmf.9'
TF = 'NR1H3'
scrna_module = cnmf_scrna_cor[,module]
scatac_module = cnmf_scatac_cor[module,]
cor (cnmf_scrna_cor, t(cnmf_scatac_cor))

mod_tf = data.frame (expression = scrna_module, activity = scatac_module)
scrna_tf = data.frame (module = cnmf_metacells[,module], TF = metacells_assay[,TF], sample = metacells$sampleID)
dev_tf = data.frame (module = cnmf_scatac[,module], TF = mMat[,TF], sample = archp$Sample)

library (smplot2)
pdf (file.path ('Plots',paste0('expression_vs_activity_TFs_module',module,'.pdf')),4,4)
ggplot (scrna_tf, aes (x = module, y = TF, color=sample)) + geom_point() + gtheme +# + facet_wrap (~sample) +
  sm_statCorr(
    color = "#0f993d", corr_method = "spearman",
    linetype = "dashed"
  )
ggplot (mod_tf, aes (x = expression, y = activity)) + geom_point() + gtheme +
  sm_statCorr(
    color = "#0f993d", corr_method = "spearman",
    linetype = "dashed"
  )
ggplot (dev_tf, aes (x = module, y = TF, color = sample)) + geom_point() + gtheme + facet_wrap (~sample) +
  sm_statCorr(
    color = "#0f993d", corr_method = "spearman",
    linetype = "dashed"
  )
dev.off()


# Highlight TF with high correlation in scatac and scrna space ####

head (mod_tf)
mod_tf = data.frame (expression = scrna_module, activity = scatac_module)
top_TFs = 30
mod_tf$mean = rowMeans (mod_tf)
mod_tf$high_TF = ''
mod_tf$high_TF[head(order(-mod_tf$mean),top_TFs)] = rownames(mod_tf)[head(order(-mod_tf$mean),top_TFs)]
mod_tf$high_TF[mod_tf$expression < 0 | mod_tf$activity < 0 ] = ''
mod_tf$hits = ifelse (mod_tf$high_TF != '','hit','nohit')
mod_tf$hits2 = ifelse (mod_tf$high_TF != '',1,0.5)
sp = ggplot (mod_tf, aes (x = expression, y = activity, label = high_TF, color = hits, alpha = hits2), color='grey22') + geom_point() + gtheme_no_rot + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred", size = 1) + # Dashed horizontal line
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkred", size = 1) + 
  scale_color_manual (values = c(hit = 'darkblue',nohit = 'grey22')) +
  geom_text_repel (size=3) 

pdf (file.path ('Plots',paste0('expression_activity_',module,'_scatter.pdf')),4,width=5)
sp
dev.off()




# Make barplot of module TFs ####
module = 'cnmf.4'
top_activity = cnmf_scatac_cor[module,][order(-cnmf_scatac_cor[module,])]
top_exp = scrna_tf_cor[,module][names(top_activity)]
module_df = data.frame (
  TF = factor(names(top_activity), levels =names(top_activity)),
 expression = top_exp,
 activity = top_activity)

bp = ggplot (module_df, aes (x = TF, y = expression)) + 
geom_bar (stat='identity',fill = 'darkgreen') +
gtheme

pdf (file.path ('Plots',paste0('activity_ranked_module_',module,'_barplot.pdf')),width=50)
bp
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
metaGroupName = 'cnmf_celltypes'
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
metaGroupName = 'cnmf_celltypes'
if (!file.exists(file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds'))))
  {
  if (!exists ('fragments')) fragments = unlist (getFragmentsFromProject (
    ArchRProj = archp_NN))   
  hubsSample_mat = matrix (ncol = length(unique(archp_NN@cellColData[,metaGroupName])), nrow = length(hubs_obj$hubsCollapsed))
  colnames (hubsSample_mat) = unique(archp_NN@cellColData[,metaGroupName])
  rownames (hubsSample_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (unique(archp_NN@cellColData[,metaGroupName])))
  for (sam in unique(archp_NN@cellColData[,metaGroupName]))
    {
    pb$tick()  
    fragments_in_sample = fragments[fragments$RG %in% rownames(archp_NN@cellColData)[as.character(archp_NN@cellColData[,metaGroupName]) == sam]]  
    fragments_in_sample_in_hubs = countOverlaps (hubs_obj$hubsCollapsed, fragments_in_sample)
    hubsSample_mat[,sam] = fragments_in_sample_in_hubs
    }
  frags_in_sample = sapply (unique(archp_NN@cellColData[,metaGroupName]), function(x) sum (archp_NN$nFrags[as.character(archp_NN@cellColData[,metaGroupName]) == x]))
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
metaGroupName = 'cnmf_celltypes'

pdf()
#archp$fetal_sample = paste0(archp$Sample, archp$fetal_group)
#metaGroupName = 'fetal_group'
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp,#[!archp$Sample3 %in% c('P11_HOX')], 
    #group_order = sample_levels, 
    #ylim = c(0,0.30),
    groupBy = metaGroupName, 
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
metaGroupName = 'cnmf_celltypes'

pdf()
#archp$fetal_sample = paste0(archp$Sample, archp$fetal_group)
#metaGroupName = 'fetal_group'
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp,#[!archp$Sample3 %in% c('P11_HOX')], 
    #group_order = sample_levels, 
    #ylim = c(0,0.30),
    groupBy = metaGroupName, 
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
  

# Separate myeloid cells in AP1-high and AP1-low ####
ap1_complex = c('FOS','FOSL2','FOSL1','JDP2','JUN','JUND','JUNB','FOSB')

mMat = scale (assays (mSE)[[1]])
rownames (mMat) = rowData (mSE)$name

srt
metaGroupName='celltype2'
ps = log2 (as.data.frame (AverageExpression (srt, features = rownames(mMat), group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
#ps = ps[, colnames(DAM_hm@matrix)]
min_exp = .1
ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
active_TFs = rownames(ps)[rowSums(ps) > min_exp]

# metaGroupName = 'cnmf_celltypes'
mMat_mg = mMat[active_TFs, ]
#mMat_mg = as.data.frame (t(mMat_mg))
# mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
# mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
# rownames (mMat_mg) = mMat_mg[,1]
# mMat_mg = mMat_mg[,-1]

colnames(mMat)

avg_AP1 = colMeans (mMat[ap1_complex,])

head (avg_AP1)
ap1_df = as.data.frame (avg_AP1)
# Create histogram
dp = ggplot(ap1_df, aes(x = avg_AP1, fill = avg_AP1 > 0)) + 
  geom_histogram(binwidth = 0.5, color = "black", alpha = 0.7) +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("Negative", "Positive")) +
  labs(title = "Histogram of Values", x = "Value", y = "Count", fill = "Sign") +
  theme_minimal()

pdf (file.path('Plots','avg_AP1_deviation_density.pdf'))
dp
dev.off()

all (rownames(ap1_df) == rownames(archp@cellColData))
archp$avg_AP1_dev = ifelse (ap1_df > 0, 'high','low')

pdf()
umap_p0 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
  #pal= palette_sample,
   embedding = "UMAP_H")
umap_p1 = lapply (unique(archp$cnmf_celltypes), function(x) plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "cnmf_celltypes",
  highlightCells = rownames(archp@cellColData)[archp$cnmf_celltypes == x],
  #pal= palette_sample,
   embedding = "UMAP_H"))

umap_p2 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "ReadsInTSS",
  #pal= palette_sample,
   embedding = "UMAP_H")
umap_p3 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "avg_AP1_dev",
  #pal= palette_sample,
   embedding = "UMAP_H")
dev.off()

archp = addImputeWeights (archp)

TF = 'NFKB1'
getFeatures (archp, 'MotifMatrix')[grep (TF, getFeatures (archp, 'MotifMatrix'))]
TF1 = c('z:JUN_143','z:FOS_137','z:NFKB1_719')

pdf ()
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "MotifMatrix",
    name = TF1, 
    useSeqnames='z',
    pal = rev (palette_deviation),    
    embedding = "UMAP_H",
    imputeWeights = getImputeWeights(archp)
    )
dev.off()
pdf(file.path('Plots','avg_AP1_deviation_fplot.pdf'),15,15)
wrap_plots(TF_p)
wrap_plots(umap_p1)
wrap_plots (umap_p0,umap_p2,umap_p3)
dev.off()





############ INTEGRATION WITH scRNA ######
archp <- addGeneIntegrationMatrix (
    ArchRProj = archp, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = srt,
    addToArrow = TRUE,
    groupRNA = "celltype2",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    force = TRUE
)

#rna_p = DimPlot (archp, group.by = 'celltype')
pdf()
int_p <- plotEmbedding(
    archp, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    embedding = 'UMAP_H'
#    pal = palD
)
int_p2 <- plotEmbedding(
    archp, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    highlightCells = rownames(archp@cellColData)[archp$predictedGroup_Un == 'Mono_CD14'],
    embedding = 'UMAP_H'
#    pal = palD
)
dev.off()
pdf (paste0 (projdir,'/Plots/RNA_integration_UMAPs.pdf'), width=15)
int_p
int_p2
dev.off()



## Add RNA integration annotation for monocytes 
archp_MAC = archp[!archp$predictedGroup_Un %in% c('DC2','Mono_CD14','Mono_CD16')]

# Check if monocytes are better defined with label transfer
TF = c('S100A9','EREG','CCL2')
metaGroupName='predictedGroup_Un'
pdf()
#archp$fetal_sample = paste0(archp$Sample, archp$fetal_group)
#metaGroupName = 'fetal_group'
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp,#[!archp$Sample3 %in% c('P11_HOX')], 
    #group_order = sample_levels, 
    #ylim = c(0,0.30),
    groupBy = metaGroupName, 
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
plotPDF (meso_markers, ArchRProj = archp,height=3.5, width=6, name =paste0('Monocytes_markers_coveragePlots.pdf'),addDOC=F)


# Generate heatmap of cnmf modules across cells and use kmeans to cluster ####
shared_cnmf = readRDS (file.path('..','scrna','shared_cnmf_myeloid.rds'))
shared_cnmf = lapply (shared_cnmf, function(x) x[x %in% getFeatures (archp)])

remove_modules = c('cnmf.3','cnmf.6') # remove monocyres and cDC modules
shared_cnmf_MAC = shared_cnmf[!names(shared_cnmf) %in% remove_modules]
cnmf_scatac = as.data.frame (t(scale(t(archp_MAC@cellColData[,names(shared_cnmf_MAC)]))))

set.seed (123)
cnmf_remove = c('cnmf.7') # remove redundant cell cycle module 
km = kmeans (scale(cnmf_scatac[,!colnames(cnmf_scatac) %in% cnmf_remove]), centers=6)  

ha = HeatmapAnnotation (sample = archp_MAC$Sample, col=list(sample = palette_sample))
hm = Heatmap (t(scale(cnmf_scatac[,!colnames(cnmf_scatac) %in% cnmf_remove])), 
  col = palette_genescore_fun(scale(cnmf_scatac)), 
  top_annotation = ha,
  show_column_dend = F,
  column_split = km$cluster,
#  column_km=9,
  row_names_gp = gpar (fontsize = 8),
  column_names_gp = gpar (fontsize = 0),
  border=T)

pdf (file.path ('Plots','cnmf_scatac_barcodes_heatmap_scaled_only_MAC2.pdf'), height=2.5)
hm
dev.off()

archp$cnmf_celltypes2 = archp$predictedGroup_Un
archp$cnmf_celltypes2[match(names(km$cluster), rownames(archp@cellColData))] = paste0('cnmf_',km$cluster)
archp$cnmf_celltypes2[archp$cnmf_celltypes2 =='cnmf_1'] = 'IM'
archp$cnmf_celltypes2[archp$cnmf_celltypes2 =='cnmf_2'] = 'IL1B'
archp$cnmf_celltypes2[archp$cnmf_celltypes2 =='cnmf_3'] = 'TREM2'
archp$cnmf_celltypes2[archp$cnmf_celltypes2 =='cnmf_4'] = 'SPP1'
archp$cnmf_celltypes2[archp$cnmf_celltypes2 =='cnmf_5'] = 'IFN'
archp$cnmf_celltypes2[archp$cnmf_celltypes2 =='cnmf_6'] = 'C1QA'

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


### Check inflammation score across TAMs ####
mod_df = data.frame (
  celltype = archp$cnmf_celltypes2,
  Infl_module = archp$mod_2)

bp = ggplot (mod_df, aes (x = celltype, y = Infl_module, fill=celltype)) +
vlp + 
bxpv + 
#scale_fill_manual (values = palette_tnk_cells) +
#geom_point (position='identity', alpha=.3, color="grey44", size=1) +
gtheme

pdf (file.path ('Plots','celltype_infl_module_boxplots.pdf'),2.3,width=4)
bp
dev.off()








### Compute correlation of myeloid modules with hubs ####
archp_meta = as.data.frame (archp@cellColData)
cnmfs = as.data.frame (archp@cellColData[,names (shared_cnmf)])
#cnmfs = scale (t(cnmfs))

module_driver = 'cnmf.8'
top_hubs = 30
traj_sample = list()
sams = unique(archp$Sample2)
#for (sam in unique(archp$Sample2))
#    {
    sam = sams  
    library(zoo)
    bin_width <- 100   # Number of observations per bin
    overlap <- 10    
    cnmf_ordered_sample = as.data.frame(scale(t(cnmfs[archp_meta$Sample %in% sam,])))
    cnmf_ordered_sample = cnmf_ordered_sample[, order(unlist(cnmf_ordered_sample[module_driver,]))]
    hubs_sample = scale(hubsCell_mat[,archp_meta$Sample %in% sam])[,colnames(cnmf_ordered_sample)]

    cnmf_ordered_sample = as.data.frame(lapply(as.data.frame(t(cnmf_ordered_sample)), function(x) {
      zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
    }))
    hub_ordered_sample = as.data.frame(lapply(as.data.frame(t(hubs_sample)), function(x) {
      zoo::rollapply(x, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left")
    }))
    colnames (hub_ordered_sample) = rownames(hubs_sample)

    cor_cnmf_hub = cor (cnmf_ordered_sample, hub_ordered_sample, method = 'spearman')
    cor_order = order (-cor_cnmf_hub[module_driver,])
    hub_ordered_sample = hub_ordered_sample[, cor_order]
    hub_ordered_sample = as.data.frame (t(hub_ordered_sample))
    hub_ordered_sample = head(hub_ordered_sample,top_hubs)#, tail (hub_ordered_sample, top_hubs))
    row_names = hubs_obj$hubsCollapsed$gene[match(rownames(hub_ordered_sample), hubs_obj$hubs_id)]
    traj_sample = Heatmap (
      t(scale(t(log2(hub_ordered_sample+1)))),
      row_labels = hubs_obj$hubsCollapsed$gene[match(rownames(hub_ordered_sample), hubs_obj$hubs_id)], 
      col = rev(palette_hubs_accessibility), 
      cluster_columns=F,
    #  name = sam,
      cluster_rows=F,
      column_names_gp = gpar(fontsize = 0),
      row_names_gp = gpar(fontsize = 7, fontface='italic'),
      border=T)
    #}

pdf (file.path('Plots',paste0('hubs_cor_to_cnmf_',module_driver,'_all.pdf')), height=4, width=5)
traj_sample
dev.off()
  
# pdf (file.path('Plots',paste0('hubs_cor_to_cnmf_',module_driver,'_sample.pdf')), height=8, width=3)
# traj_sample
# dev.off()


### Check TF enrichment in HUBS correlated with TREM2 vs anticorrelated ####
### TF Enrichment in peaks in large hubs ####
tf_match = getMatches (archp)
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)
hubs_to_enrich = hubs_obj$hubsCollapsed[match (rownames (hub_ordered_sample), hubs_obj$hubs_id)]
hubs_to_enrich_peaks = bg_peakSet[queryHits(findOverlaps(bg_peakSet, hubs_to_enrich))]
#tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
hubs_TF_res =  hyperMotif (
  selected_peaks = hubs_to_enrich_peaks, 
  motifmatch = tf_match)


### Group macrophages by TREM2 module score ####
#module_driver = 'shared_cnmf.6'
cnmfs = as.data.frame (t(scale(t(archp@cellColData[,names(shared_cnmf)]))))
#cnmfs = as.data.frame (do.call (rbind, lapply (unique(archp$Sample), function(x) t(scale(t(cnmfs[archp$Sample == x,]))))))
archp$TREM2_state3 = ifelse (cnmfs$cnmf.3 > -0.08,'high','low')
module= 'cnmf.4'
TREM2_state = cnmfs[,module]

# Try sorting based on TREM2 genescore alone ####
if (!exists('gsSE')) gsSE = fetch_mat (archp, 'GeneScoreMatrix')
gsMat = assay (gsSE)
rownames (gsMat) = rowData(gsSE)$name
gsMat = log2(gsMat['TREM2',]+1)


pdf (file.path ('Plots','TREM2_genescore_distribution.pdf'))
hist (gsMat)
dev.off()

all (names(gsMat) == rownames(archp@cellColData))
archp$TREM2_state3 = ifelse (gsMat>1, 'high','low')

names (TREM2_state) = rownames(archp@cellColData)
TREM2_state = TREM2_state[order(-TREM2_state)]

n <- length(TREM2_state) # Get the length of the vector
part_size <- ceiling(n / 6) # Calculate the size of each part

# Create indices for splitting
splits <- split(TREM2_state, ceiling(seq_along(TREM2_state) / part_size))
labeled_splits <- setNames(splits, c("very_high","high", "low","very_low"))

df <- data.frame(
  value = unlist(labeled_splits),
  label = rep(names(labeled_splits), times = sapply(labeled_splits, length))
)

archp$TREM2_state2 = df$label[match(rownames(archp@cellColData), sapply(rownames(df), function(x) unlist(strsplit(x,'\\.'))[2]))] 
#table (archp$TREM2_state2, archp$TREM2_state)

pdf()
metaGroupName = 'TREM2_state3'
celltype_markers = c('TREM2','NFATC2','APOE','NR1H3','FTL','LIPA','VIM','CD9','CD63','APOC1','C1QA')
meso_markers <- plotBrowserTrack(
    ArchRProj = archp, 
    groupBy = metaGroupName,
    scCellsMax = 2000, 
    geneSymbol = celltype_markers,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 200000,
    downstream = 200000,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
    #loops = getCoAccessibility (archp, corCutOff = 0.3,
    #  returnLoops = TRUE),
    useGroups= NULL
)
dev.off()
plotPDF (meso_markers, ArchRProj = archp, width=14, name ='MPM_markers_coveragePlots.pdf')

pdf (file.path ('Plots','TREM2_state_dimplot.pdf'))
umap_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
 name = "TREM2_state3", embedding = "UMAP_H")
umap_p1
dev.off()

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




TF = 'JUN'
TF = 'VCAN'
getFeatures (archp, 'MotifMatrix')[grep (TF, getFeatures (archp, 'MotifMatrix'))]
TF1 = 'z:JUN_143'
archp = addImputeWeights (archp)
pdf ()
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "MotifMatrix",
    name = TF1, 
    useSeqnames='z',
    pal = rev (palette_deviation),    
    embedding = "UMAP_H",
    imputeWeights = getImputeWeights(archp)
    )
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = TF, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path ('Plots','TF_featureplots.pdf'), width = 8,height=4)
wrap_plots (TF_p, p)
dev.off()

mono_markers = c('VCAN','FCN1','CXCL8','CXCL2','IL1B','EREG','TIMP1','THBS1','CCR2','FLT3','FOXM1','CDK1','PCNA','FCGR3A')
pdf()
p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = mono_markers, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path ('Plots','TF_mono_featureplots.pdf'), width = 18,height=14)
wrap_plots (p2)
dev.off()

# Cluster each sample individually to better find monocytes ####
  varfeat = 25000
  LSI_method = 2
sams = unique(archp$Sample)
remove_samples = c('P8','P3','P4') # remove low number samples
sams = sams[! sams %in% remove_samples]
archp_sub_l = list()
for (sam in sams)
  {
  archp_sub = archp[archp$Sample == sam]
  archp_sub = addIterativeLSI (ArchRProj = archp_sub,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = TRUE, LSIMethod = LSI_method,
    varFeatures = varfeat)

  archp_sub = addClusters (input = archp_sub, resolution = 3,
    reducedDims = "IterativeLSI", maxClusters = 100,
    force = TRUE)
  archp_sub = addUMAP (ArchRProj = archp_sub, 
    reducedDims = "IterativeLSI",
    force = TRUE)
  archp_sub_l[[sam]] = archp_sub
  }

mono_markers = c('VCAN','FCN1','CXCL8','CXCL2','IL1B','EREG','TIMP1','THBS1','CCR2','FLT3','FOXM1','CDK1','PCNA','FCGR3A')
pdf()
p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = mono_markers, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path ('Plots',paste0(sam,'monocyte_markers_umap.pdf')))
p2
dev.off()


























### Combine TF modisco and finemo outputs to build network of co-occurring TFs across peaks ####
library (httr)
library (XML)
library (igraph)

chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/chromBPnet'
fold_number = 0
celltypes = c('IL1B','Mono','TREM2','SPP1','C1Q','IM')


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
  modisco_motifs = as.data.frame(readHTMLTable(paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/chromBPnet/',celltype,'_model/fold_0/modisco/report/motifs.html')))
  modisco_motifs$motif_match0 = sapply (modisco_motifs$NULL.match0, function(x) unlist(strsplit (x, '_'))[1])
  modisco_motifs$motif_match1 = sapply (modisco_motifs$NULL.match1, function(x) unlist(strsplit (x, '_'))[1])
  modisco_motifs$motif_match2 = sapply (modisco_motifs$NULL.match2, function(x) unlist(strsplit (x, '_'))[1])

  finemo_hits = read.table(paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/chromBPnet/',celltype,'_model/finemo_out/hits.tsv'), sep='\t', header=T)
  finemo_hits$motif_name0 = modisco_motifs$motif_match0[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
  finemo_hits$motif_name1 = modisco_motifs$motif_match1[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
  finemo_hits$motif_name2 = modisco_motifs$motif_match2[match(finemo_hits$motif_name, modisco_motifs$NULL.pattern)]
  peakset = read.table (file.path(chromBPdir,paste0('peakset_',celltype,'.bed')))
  peakset = peakset[,1:3]
  colnames (peakset) = c ('chr','start','end')
  peakset = makeGRangesFromDataFrame (peakset)
  finemo_hits = makeGRangesFromDataFrame (finemo_hits, keep.extra.columns = T)
  
  match_rank = c('motif_name0','motif_name1','motif_name2')
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




