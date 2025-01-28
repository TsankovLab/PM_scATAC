conda activate meso_scatac
R

set.seed(1234)

packages = c(
  'ArchR',
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
  'BSgenome.Hsapiens.UCSC.hg38')
lapply(packages, require, character.only = TRUE)

set.seed (1234)
addArchRThreads (threads = 8) 
addArchRGenome ("Hg38")

####### START ANALYSIS #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/ENCODE_normal'
dir.create (file.path(projdir,'Plots'), recursive =T)
setwd (projdir)

# Load function (JGranja)
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))


# combine fragments files from different channels of the normal RPL 
#system ('cat /ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz \
#/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz \
#/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz \
#/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz > RPL_normal_fragments.tsv.gz')

# Fix the fragment file of multiome sample by removing the header 
#system ('zcat atac_fragments.tsv.gz | grep -v ^\# | bgzip > atac_fragments_fixed.tzv.gz')
sample_names = list.files ('/sc/arion/scratch/giottb01/ENCODE_scATAC-seq_lung/tarfiles/compressed_fragments', pattern='tsv.gz')
fragment_paths = file.path ('/sc/arion/scratch/giottb01/ENCODE_scATAC-seq_lung/tarfiles/compressed_fragments',sample_names)
sample_names = sub ('.tsv.gz','',sample_names)

names (fragment_paths) = sample_names


ArrowFiles = createArrowFiles (inputFiles = fragment_paths,
      sampleNames = sample_names,
      minTSS = 4, #Dont set this too high because you can always increase later
      minFrags = 1000,
      maxFrags = Inf,
      addTileMat = TRUE,
      addGeneScoreMat = TRUE,
      force = TRUE,
      subThreading = T
      )
 

  archp = ArchRProject (
    ArrowFiles = ArrowFiles, 
    outputDirectory = projdir,
    copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  )
  

  ### Subset ArchR object only for cells retained in Signac analysis ####
  ### QC plots ####
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
pdf()
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
  p2 <- plotGroups(
    ArchRProj = archp, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    pal = palette_sample,
    alpha = 0.4,
    addBoxPlot = TRUE
   )  
dev.off()
  pdf (file.path ('Plots', 'QC_plots.pdf'))
  wrap_plots (p1, p2 ,p3, p4)
  dev.off()

archp = loadArchRProject (projdir)
archp = archp[archp$TSSEnrichment > 6 & archp$ReadsInTSS > 500]

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

# Harmony ####
archp = addHarmony (
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample", force=TRUE
)
archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony", name='UMAP_H',
    force = TRUE)
archp = addClusters (input = archp,
    reducedDims = "Harmony", resolution = 0.2,
    name='Clusters_H',
    force = TRUE)

pdf ()  
umap_p0 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters",
   embedding = "UMAP",
   #pal = palette_celltype_simplified,
   labelMeans = FALSE)

umap_p2 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
   embedding = "UMAP",
  # pal = palette_sample,
   labelMeans = FALSE)

umap_p1 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters_H",
   embedding = "UMAP_H",
   #pal = palette_celltype_simplified,
   labelMeans = FALSE)

umap_p3 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
   embedding = "UMAP_H",
  # pal = palette_sample,
   labelMeans = FALSE)
umap_p4 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
   embedding = "UMAP",
  # pal = palette_sample,
   labelMeans = FALSE)

dev.off()

pdf (file.path('Plots','clusters_umap.pdf'), width = 7, height = 7)
umap_p0
umap_p2
umap_p1
umap_p3
umap_p4
dev.off()



genes = c('HP','WT1',  
  'CD3D', 'LYZ','PECAM1','COL1A1',
  'VWF','KLRC1','GNLY','SFTPA1',
  'COL1A2','EPCAM','CALB2','CCR7','GATA2',
  'ITLN1','TEAD1')
genes_meso = c('HAS1','PHYHIP','CALB2','WT1','HP', 'GATA4')

pdf()
p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = genes_meso, 
    embedding = "UMAP_H",
    pal = viridis::plasma(100),
    imputeWeights = NULL
)

pdf (file.path('Plots','marker_genes_feature_plots_3.pdf'), width = 25, height = 25)
wrap_plots (p2)
dev.off()

#archp = archp[archp$Clusters != 'C1']

#   varfeat = 25000
#   LSI_method = 2
#   archp = addIterativeLSI (ArchRProj = archp,
#     useMatrix = "TileMatrix", name = "IterativeLSI",
#     force = TRUE, LSIMethod = LSI_method,
#     varFeatures = varfeat)

#   archp = addClusters (input = archp, resolution = 3,
#     reducedDims = "IterativeLSI", maxClusters = 100,
#     force = TRUE)
#   archp = addUMAP (ArchRProj = archp, 
#     reducedDims = "IterativeLSI",
#     force = TRUE)

# pdf ()  
# umap_p0 = plotEmbedding (ArchRProj = archp, 
#   colorBy = "cellColData", name = "Clusters",
#    embedding = "UMAP",
#    #pal = palette_celltype_simplified,
#    labelMeans = FALSE)

# umap_p2 = plotEmbedding (ArchRProj = archp, 
#   colorBy = "cellColData", name = "Sample",
#    embedding = "UMAP",
#   # pal = palette_sample,
#    labelMeans = FALSE)
# dev.off()

genes = c('HP','WT1',  
  'CD3D', 'LYZ','PECAM1','COL1A1',
  'VWF','KLRC1','GNLY','SFTPA1',
  'COL1A2','EPCAM','CALB2','CCR7','GATA2',
  'ITLN1','TEAD1')
genes_meso = c('HAS1','PHYHIP','CALB2','WT1','HP')


# archp = addImputeWeights (archp)
pdf()
# p <- plotEmbedding(
#     ArchRProj = archp,
#     colorBy = "GeneScoreMatrix", 
#     name = genes_meso, 
#     embedding = "UMAP",
#     pal = rev(viridis::plasma(100)),
#     imputeWeights = getImputeWeights(archp)
# )

p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = genes, 
    embedding = "UMAP",
    pal = rev(viridis::plasma(100)),
    imputeWeights = NULL
)

umap_p2 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
   embedding = "UMAP",
  # pal = palette_sample,
   labelMeans = FALSE)
dev.off()

pdf (file.path('Plots','marker_genes_feature_plots_3.pdf'), width = 25, height = 25)
wrap_plots (p2)
umap_p2
dev.off()


samples = unique(archp$Sample)
archp_l = list()
for (sample in samples)
{
archp_sub = archp[archp$Sample == sample]
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

archp_sub = addImputeWeights(archp_sub)
pdf ()
p <- plotEmbedding(
    ArchRProj = archp_sub,
    colorBy = "GeneScoreMatrix", 
    name = genes, 
    embedding = "UMAP",
    #pal = palette_expression,
    imputeWeights = getImputeWeights(archp_sub)
)

umap_p2 = plotEmbedding (ArchRProj = archp_sub, 
  colorBy = "cellColData", name = "Clusters",
   embedding = "UMAP",
  # pal = palette_sample,
   labelMeans = FALSE)
dev.off()

archp_l[[sample]] = archp_sub
pdf (file.path('Plots',paste0('marker_genes_feature_plots_3_',sample,'.pdf')), width = 25, height = 25)
print (wrap_plots (p, ncol = 8))
dev.off()
pdf (file.path ('Plots',paste0('celltype_revised_umap_',sample,'.pdf')))
umap_p2
dev.off()

}

#p = lapply (p, function(x) x + theme_void() + NoLegend ()) #+ ggtitle scale_fill_gradient2 (rev (viridis::plasma(100))))

pdf (file.path ('Plots','celltype_revised_umap.pdf'))
umap_p0
#umap_p1
umap_p2
dev.off()

pdf (file.path('Plots','marker_genes_feature_plots_3.pdf'), width = 25, height = 25)
print (wrap_plots (p, ncol = 8))
dev.off()


pdf()
p1 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "cellColData", 
    name = 'nFrags', 
    embedding = "UMAP",
    #pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)

p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "cellColData", 
    name = 'TSSEnrichment', 
    embedding = "UMAP",
    #pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()

pdf (file.path('Plots','QC_umap.pdf'), width = 25, height = 25)
print (wrap_plots (list(p1, p2)))
dev.off()

# Run genescore DAG ####
metaGroupName = "Clusters"
force=TRUE
if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source (file.path('..','..','git_repo','utils','DAG.R'))

pdf (paste0('Plots/DAG_clusters_',metaGroupName,'_heatmaps.pdf'), width = 8, height = 15)
print (DAG_hm)
dev.off()

 

# Make module score of mesothelial markers
meso_markers = list(mesothelium = c('KRT6B','CDCP1','BNC1','MIR31','SLC28A3','OR1M1','WT1','WT1-AS','HAS1','UPK1B','MIR4725','ALDH8A1'))
archp = addModuleScore (
  ArchRProj = archp,
  useMatrix = 'GeneScoreMatrix',
  name = "normal",
  features = meso_markers,
  nBin = 25,
  nBgd = 100,
  seed = 1,
  threads = getArchRThreads(),
  logFile = createLogFile("addModuleScore")
)

pdf()
p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "cellColData", 
    name = 'normal.mesothelium', 
    embedding = "UMAP_H",
    #pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()

pdf (file.path('Plots','meso_module_umap.pdf'), width = 7, height = 7)
p2
dev.off()

pdf()
p2 <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "cellColData", 
    name = 'normal.mesothelium', 
    embedding = "UMAP",
    #pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()

archp = addClusters (input = archp, resolution = 10,
    reducedDims = "IterativeLSI", maxClusters = 100,
    force = TRUE)
umap_p0 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters",
   embedding = "UMAP", 
   #pal = palette_celltype_simplified,
   labelMeans = TRUE) + NoLegend()

pdf (file.path('Plots','meso_module_umap2.pdf'), width = 12, height = 7)
wrap_plots (p2, umap_p0)
dev.off()

### Select mesothelial cells cluster
write.csv (rownames(archp@cellColData)[archp$Clusters == 'C80'], 'mesothelium_barcodes.csv')