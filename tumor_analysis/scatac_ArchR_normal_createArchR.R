conda activate meso_scatac
R

projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/normal_lung/scatac_ArchR'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

set.seed(1234)


# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

# Set # of threads and genome reference ####
addArchRThreads (threads = 8) 
addArchRGenome ("hg38")

sample_names = c(
  'RPL_280_neg_1',
  'RPL_280_neg_2',
  'RPL_Epi_1',
  'RPL_Epi_2'
  )

  # Fix the fragment file of multiome sample by removing the header 
  #system ('zcat atac_fragments.tsv.gz | grep -v ^\# | bgzip > atac_fragments_fixed.tzv.gz')

fragment_paths =c('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/normal_lung')
fragment_paths = file.path(fragment_paths, paste0(sample_names,'_fragments.tsv.gz'))
   
ArrowFiles = createArrowFiles (
inputFiles = fragment_paths,
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
  copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

# Dimensionality reduction and clustering
varfeat = 25000
LSI_method=2
archp = addIterativeLSI (ArchRProj = archp, 
  useMatrix = "TileMatrix", name = "IterativeLSI",
  force=TRUE, LSIMethod=LSI_method,
  varFeatures = varfeat)
archp = addClusters (input = archp, resolution = .8,
  reducedDims = "IterativeLSI", maxClusters = 100,
  force = TRUE)
archp = addUMAP (ArchRProj = archp, 
  reducedDims = "IterativeLSI",
  force = TRUE)

meso_markers = c('C1QA','WT1','HP','GATA4','CD3D','COL1A1','PECAM1')
#archp = addImputeWeights (archp)
pdf()
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = meso_markers, 
    size=1,
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = NULL
)
dev.off()
pdf (file.path('Plots','markers_fplots.pdf'), width = 18, height = 15)
wrap_plots (p, ncol=3)
dev.off()

pdf()
umap_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
 name = "Sample", embedding = "UMAP", pal = palette_sample)
umap_p2 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters",
   embedding = "UMAP")
umap_p3 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "nFrags",
   embedding = "UMAP")
umap_p4 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "TSSEnrichment",
   embedding = "UMAP")
dev.off()

pdf (file.path('Plots','celltype_umap.pdf'))
print (umap_p1)
print (umap_p2)
dev.off()

# Save ArchR object ####
archp = saveArchRProject (archp[archp$Clusters == 'C12'], dropCells = T)
