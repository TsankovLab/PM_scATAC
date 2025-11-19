conda activate meso_scatac
R

set.seed(1234)

####### ANALYSIS of TUMOR compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

# Set # of threads and genome reference ####
addArchRThreads (threads = 8) 
addArchRGenome ("hg38")

sample_names_tumor = c(
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
  'P14',# p14
  'P23'
  )

sample_names_normal = c(
  # Normal
  'RPL_280_neg_1',
  'RPL_280_neg_2',
  'RPL_Epi_1',
  'RPL_Epi_2',#,
  'ENCSR592WYM-1_fragments',
  'ENCSR606VJQ-1_fragments',
  'ENCSR548BXU-1_fragments'
  )
  #'cf_distal')
  # Fix the fragment file of multiome sample by removing the header 
  #system ('zcat atac_fragments.tsv.gz | grep -v ^\# | bgzip > atac_fragments_fixed.tzv.gz')
ArrowFiles_tumor_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR'
ArrowFiles_tumor_dir = file.path(ArrowFiles_tumor_dir, paste0(sample_names_tumor,'.arrow'))
ArrowFiles_normal_dir = c(
  '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/normal_lung/scatac_ArchR/ArrowFiles/RPL_280_neg_1.arrow',
  '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/normal_lung/scatac_ArchR/ArrowFiles/RPL_280_neg_2.arrow',
  '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/normal_lung/scatac_ArchR/ArrowFiles/RPL_Epi_1.arrow',
  '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/normal_lung/scatac_ArchR/ArrowFiles/RPL_Epi_2.arrow',
  '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/ENCODE_normal/ArrowFiles/ENCSR592WYM-1_fragments.arrow',
  '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/ENCODE_normal/ArrowFiles/ENCSR606VJQ-1_fragments.arrow',
  '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/ENCODE_normal/ArrowFiles/ENCSR548BXU-1_fragments.arrow'
  )
ArrowFiles_dir = c(ArrowFiles_tumor_dir, ArrowFiles_normal_dir)

# if (!all (file.exists(paste0('ArrowFiles/',c(paste0(sample_names_tumor,'.arrow'),paste0(sample_names_normal,'.arrow')))))) 
#   {
system ('mv ArrowFiles/*.arrow ../')
archp = ArchRProject (
ArrowFiles = ArrowFiles_dir,
outputDirectory = projdir,
copyArrows = FALSE) 
  # }
  # } else {
  # archp = ArchRProject (
  # ArrowFiles = ArrowFiles_dir,
  # outputDirectory = projdir,
  # copyArrows = FALSE)  
  # # }

# Further subset malignant cells as there is a bug in ArchR that doesnt update arrow metadata when re-reading arrow files from subset object
ann_tumor = read.csv (file.path ('..','..','git_repo','files','barcode_annotation.csv'))
ann_tumor = ann_tumor[ann_tumor$celltype_lv1 == 'Malignant','barcode']
ann_normal = read.csv (file.path ('..','..','normal_lung','scatac_ArchR','mesothelium_annotation.csv'))
ann_normal = ann_normal$x
ann_normal2 = read.csv (file.path ('..','ENCODE_normal','mesothelium_barcodes.csv'))
ann_normal2 = ann_normal2$x

ann_combined = c(ann_tumor, ann_normal, ann_normal2)
archp = archp[rownames(archp@cellColData) %in% ann_combined]

#if (!all (file.exists(paste0('ArrowFiles/',c(paste0(sample_names_tumor,'.arrow'),paste0(sample_names_normal,'.arrow'))))))
archp = saveArchRProject (archp)

archp$Sample2 = archp$Sample
archp$Sample2[grep ('RPL',archp$Sample2)] = 'normal1'
archp$Sample2[grep ('ENCSR548BXU-1_fragments',archp$Sample2)] = 'normal2'
archp$Sample2[grep ('ENCSR592WYM-1_fragments',archp$Sample2)] = 'normal3'
archp$Sample2[grep ('ENCSR606VJQ-1_fragments',archp$Sample2)] = 'normal4'

# Dimensionality reduction and clustering
varfeat = 25000
LSI_method=2
archp = addIterativeLSI (ArchRProj = archp, 
  useMatrix = "TileMatrix", name = "IterativeLSI",
  force=TRUE, LSIMethod=LSI_method,
  varFeatures = varfeat)
archp = addClusters (input = archp, resolution = .6,
  reducedDims = "IterativeLSI", maxClusters = 100,
  force = TRUE)
archp = addUMAP (ArchRProj = archp, 
  reducedDims = "IterativeLSI", seed = 2,
  force = TRUE)

pdf()
umap_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
 name = "Sample2", embedding = "UMAP", pal = palette_sample)
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

pdf (file.path('Plots','celltype_umap.pdf'),10,10)
print (umap_p1)
print (umap_p2)
dev.off()
