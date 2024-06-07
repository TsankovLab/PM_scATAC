ish -q interactive -l h_vmem=8g -pe smp 24 -binding linear:24
use UGER
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
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/all_tissues_ArchR'
datadir = '/broad/hptmp/mbairakd/data/bed_files'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

#devtools::install_github("immunogenomics/presto") #needed for DAA
source ('../../PM_scATAC/useful_functions.R')
source ('../../PM_scATAC/ggplot_aestetics.R')
source ('../../PM_scATAC/scATAC_functions.R')
source ('../../PM_scATAC/palettes.R')

set.seed (1234)
addArchRThreads (threads = 24) 
addArchRGenome ("Hg38")

   
fragment_paths = c(
  file.path(datadir,'JShendure_scATAC', 'migrated_to_hg38', list.files (file.path (datadir, 'JShendure_scATAC','migrated_to_hg38'), pattern = 'bgz$')),
  #file.path(datadir,'bingren_adult_brain', list.files (file.path (datadir, 'bingren_adult_brain'), pattern = 'bgz')),
  file.path(datadir,'greenleaf_brain_scATAC', list.files (file.path (datadir, 'greenleaf_brain_scATAC'), pattern = 'bgz$')),
  file.path(datadir,'yang_kidney_scATAC','migrated_to_hg38', list.files (file.path (datadir, 'yang_kidney_scATAC','migrated_to_hg38'), pattern = 'bgz$')),
  file.path(datadir,'Tsankov_scATAC', list.files (file.path (datadir, 'Tsankov_scATAC'), pattern = 'bgz$')),
  file.path(datadir,'bingren_scATAC', list.files (file.path (datadir, 'bingren_scATAC'), pattern = 'bgz$')),
  file.path(datadir,'greenleaf_colon_scATAC', list.files (file.path (datadir, 'greenleaf_colon_scATAC'), pattern = 'bgz$')),
  file.path(datadir,'rawlins_fetal_lung_scATAC', list.files (file.path (datadir, 'rawlins_fetal_lung_scATAC'), pattern = 'bgz$'))
  )
  sample_names = c(
   paste0 ('JShendure', 1:length(list.files (file.path (datadir, 'JShendure_scATAC','migrated_to_hg38'), pattern = 'bgz$'))),
    #paste0 ('bingren_adult_brain', 1:length(list.files (file.path (datadir, 'bingren_adult_brain'), pattern = 'bgz'))),
    paste0 ('greenleaf_brain', 1:length(list.files (file.path (datadir, 'greenleaf_brain_scATAC'), pattern = 'bgz$'))),
    paste0 ('yang_kidney', 1:length(list.files (file.path (datadir, 'yang_kidney_scATAC','migrated_to_hg38'), pattern = 'bgz$'))),
    paste0 ('Tsankov_lung', 1:length(list.files (file.path (datadir, 'Tsankov_scATAC'), pattern = 'bgz$'))),
    paste0 ('bingren_pan', 1:length(list.files (file.path (datadir, 'bingren_scATAC'), pattern = 'bgz$'))),
    paste0 ('greenleaf_colon', 1:length(list.files (file.path (datadir, 'greenleaf_colon_scATAC'), pattern = 'bgz$'))),
    paste0 ('rawlins_fetal_lung', 1:length(list.files (file.path (datadir, 'rawlins_fetal_lung_scATAC'), pattern = 'bgz$')))
    )

  #setwd (projdir)  
  ArrowFiles = createArrowFiles (inputFiles = fragment_paths,
  sampleNames = sample_names,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000,
  maxFrags = Inf,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = FALSE,
  subThreading = T
  )
  
  archp = ArchRProject (
    ArrowFiles = ArrowFiles, 
    outputDirectory = projdir,
    copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  )
  

  ### Run peak calling on celltype annotation ####

run_peakCall = FALSE
if (run_peakCall)
  {
  ### Call peaks on celltypes ###
  metaGroupName = 'celltype'
  archp = addGroupCoverages (
    ArchRProj = archp, 
    groupBy = metaGroupName,  
    force = FALSE,
    minCells= 20, # I think this should be set corresponding to the smallest cluster in the group or lower
    maxCells = 500,
    minReplicates = 2,
    sampleRatio = 0.8,
    useLabels = TRUE)
  
  archp = addReproduciblePeakSet (
      archp,
      groupBy= metaGroupName,
      peakMethod = 'Macs2',
      reproducibility = "1",
      maxPeaks = 500000, 
      minCells=20,
      force =TRUE) # I think this should be set corresponding to the smallest cluster in the group or lower
  archp = addPeakMatrix (archp)
  
  archp = saveArchRProject (archp, load=TRUE)
  
  metaGroupNames = c('TSSEnrichment','nFrags','ReadsInTSS','FRIP')  
    umap_p12 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
     name = x, embedding = "UMAP"))
      
  pdf (paste0(projdir,'/Plots/qc_umap_after_filtering.pdf'), 15,15)
  wrap_plots (umap_p12, ncol=5)
  dev.off()
  }
