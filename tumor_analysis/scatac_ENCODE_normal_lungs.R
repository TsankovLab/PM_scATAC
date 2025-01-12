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
sample_names = list.files ('/sc/arion/scratch/giottb01/ENCODE_scATAC-seq_lung', pattern='bam$')
fragment_paths = file.path ('/sc/arion/scratch/giottb01/ENCODE_scATAC-seq_lung',sample_names)
sample_names = sub ('.bam','',sample_names)

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
  cell_annotation = read.csv ('../../per_sample_QC_signac/cell_annotation.csv')
  colnames (cell_annotation) = c('barcode','celltype')
  cell_annotation$barcode = sub ('\\.','#', cell_annotation$barcode)
  
  archp$celltype = cell_annotation$celltype[match(rownames(archp@cellColData), cell_annotation$barcode)]
  archp$celltype[archp$celltype == 'Myeloid'] = 'MonoMac'
  archp$celltype[archp$celltype == 'pDC'] = 'pDCs'
  archp = archp[!is.na(archp$celltype)]
  archp = archp[archp$celltype != 'bad_quality']
  
  # Import revised cell annotation after multiple iterations of Umaps (NOT SAVED)
  ann = read.csv ('../../git_repo/files/barcode_annotation.csv')
  archp$celltype_revised = ann$celltype[match(rownames(archp@cellColData), ann$barcode)]
  archp = archp[!is.na(archp$celltype_revised)]

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
   #pal = palette_celltype_simplified,
   labelMeans = FALSE)

umap_p1 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype_revised",
   embedding = "UMAP",
   pal = palette_celltype_simplified,
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

