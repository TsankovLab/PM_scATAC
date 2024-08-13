#ish -q interactive -l h_vmem=64g -pe smp 4 -binding linear:4
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
addArchRThreads (threads = 8) 
addArchRGenome ("Hg38")

   
fragment_paths = c(
  #file.path(datadir,'JShendure_scATAC', 'migrated_to_hg38', list.files (file.path (datadir, 'JShendure_scATAC','migrated_to_hg38'), pattern = 'bgz$')),
  #file.path(datadir,'bingren_adult_brain', list.files (file.path (datadir, 'bingren_adult_brain'), pattern = 'bgz$'))#,
  #file.path(datadir,'greenleaf_brain_scATAC', list.files (file.path (datadir, 'greenleaf_brain_scATAC'), pattern = 'bgz$')),
  #file.path(datadir,'yang_kidney_scATAC','migrated_to_hg38', list.files (file.path (datadir, 'yang_kidney_scATAC','migrated_to_hg38'), pattern = 'bgz$')),
  #file.path(datadir,'Tsankov_scATAC', list.files (file.path (datadir, 'Tsankov_scATAC'), pattern = 'bgz$')),
  #file.path(datadir,'bingren_scATAC', list.files (file.path (datadir, 'bingren_scATAC'), pattern = 'bgz$')),
  #file.path(datadir,'greenleaf_colon_scATAC', list.files (file.path (datadir, 'greenleaf_colon_scATAC'), pattern = 'bgz$'))#,
  #file.path(datadir,'rawlins_fetal_lung_scATAC', list.files (file.path (datadir, 'rawlins_fetal_lung_scATAC'), pattern = 'bgz$'))
  )
#fragment_paths = file.path(datadir,'bingren_adult_brain','sorted.bgzip')
# fragment_paths = c(
#   file.path(datadir,'JShendure_scATAC', 'migrated_to_hg38', 'sample_10_muscle.fragments.txt.bgz'),
#   file.path(datadir,'bingren_adult_brain', 'GSM7822133_MM_439.bed.bgz')
#   )
  sample_names = c(
   #paste0 ('JShendure', 1:length(list.files (file.path (datadir, 'JShendure_scATAC','migrated_to_hg38'), pattern = 'bgz$'))),
    #paste0 ('bingren_adult_brain', 1:length(list.files (file.path (datadir, 'bingren_adult_brain'), pattern = 'bgz$')))#,
    #paste0 ('greenleaf_brain', 1:length(list.files (file.path (datadir, 'greenleaf_brain_scATAC'), pattern = 'bgz$'))),
    #paste0 ('yang_kidney','migrated_to_hg38', 1:length(list.files (file.path (datadir, 'yang_kidney_scATAC','migrated_to_hg38'), pattern = 'bgz$'))),
    #paste0 ('Tsankov_lung', 1:length(list.files (file.path (datadir, 'Tsankov_scATAC'), pattern = 'bgz$'))),
    #paste0 ('bingren_pan', 1:length(list.files (file.path (datadir, 'bingren_scATAC'), pattern = 'bgz$'))),
    #paste0 ('greenleaf_colon', 1:length(list.files (file.path (datadir, 'greenleaf_colon_scATAC'), pattern = 'bgz$')))#,
    #paste0 ('rawlins_fetal_lung', 1:length(list.files (file.path (datadir, 'rawlins_fetal_lung_scATAC'), pattern = 'bgz$')))
    )
  

  archp = ArchRProject (
    ArrowFiles = list.files (projdir, pattern = '.arrow'), 
    outputDirectory = projdir,
    copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  )
  
table (archp$Sample)
archp$project =  gsub("^\\d+|\\d+$", "", archp$Sample)    
table (archp$project)

### Import metadata ####
meta_dir = '/broad/hptmp/bgiotti/BingRen_scATAC_atlas/data/metadata'

meta_files = c(
  'GSE184462_metadata.csv',
  #bingren_adult_brain_metadata.csv
  'GSE149683_File_S2.Metadata_of_high_quality_cells.csv',
  'tsankov_refined_annotation.csv',
  'GSE162170_atac_cell_metadata.csv',
  'yang_kidney.csv',
  'greenleaf_colon_metadata.csv',
  'rawlins_fetal_lung_metadata.csv'
  )
#gse184 = read.table (file.path (meta_dir2, 'GSE184462_metadata.tsv'), sep='\t', header=T)
meta_dir2 = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/all_tissues_ArchR/metadata_external_data'
meta = lapply (meta_files, function(x) read.csv (file.path (meta_dir2, x), row.names=1))


arrows = list.files (projdir, pattern = '.arrow')

source ('../../PM_scATAC/all_tissues_scatac_annotate.R')

arrows_controls = arrows[grep ('rawlins', arrows)]
arrows_controls = append (arrows_controls, arrows[grep ('Tsankov', arrows)])
arrows_controls = append (arrows_controls, arrows[grep ('bingren', arrows)])

metaGroupName = 'celltype'
archp_prj = ArchRProject (
ArrowFiles = arrows_controls, 
outputDirectory = file.path(projdir,'normal_controls'),
copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

dir.create (file.path(projdir,'normal_controls','Plots')) 
archp_prj$celltype = as.character(archp$celltype[match(rownames(archp_prj@cellColData), rownames(archp@cellColData))])
archp_prj = archp_prj[!is.na(archp_prj$celltype)]
archp_prj = archp_prj[archp_prj$celltype != 0]
archp_prj$celltype = gsub (' ','',archp_prj$celltype)
archp_prj$celltype = gsub ('/','_',archp_prj$celltype)
archp_prj$celltype = gsub ('-','_',archp_prj$celltype)
archp_prj$project = gsub("^\\d+|\\d+$", "", archp_prj$Sample)    
archp_prj$celltype_project = paste0(archp_prj$project,'_', archp_prj$celltype)
table (archp_prj$celltype_project[grep ('meso', archp_prj$celltype_project, ignore.case = T)])
table (archp_prj$celltype_project[grep ('bingren_pan_Adult_lung', archp_prj$celltype_project, ignore.case = T)])

archp_prj = addIterativeLSI (
ArchRProj = archp_prj, 
useMatrix = "TileMatrix", 
name = "IterativeLSI",
force=FALSE)
archp_prj = addUMAP (
ArchRProj = archp_prj, 
reducedDims = "IterativeLSI",
force = FALSE)
  
  umap_p1 = lapply (metaGroupName, function(x) plotEmbedding (ArchRProj = archp_prj, colorBy = "cellColData",
   name = x, embedding = "UMAP"))
    
  pdf (file.path(projdir,prj,'Plots',paste0('celltype_umap.pdf')), 10,10)
  print (umap_p1)
  dev.off()
  
  archp_prj = addGroupCoverages (
    ArchRProj = archp_prj, 
    groupBy = metaGroupName,  
    force = TRUE,
    minCells= 20, # I think this should be set corresponding to the smallest cluster in the group or lower
    maxCells = 500,
    minReplicates = 2,
    sampleRatio = 0.8,
    useLabels = TRUE)  
  archp_prj = addReproduciblePeakSet (
      archp_prj,
      groupBy= metaGroupName,
      peakMethod = 'Macs2',
      reproducibility = "1",
      maxPeaks = 500000, 
      minCells=20,
      force =FALSE) 
  }


