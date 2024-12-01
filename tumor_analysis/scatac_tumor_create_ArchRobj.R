library (ArchR)

projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

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
addArchRThreads (threads = 1) 
addArchRGenome ("hg38")

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
  'P14',# p14
  'P23',
  # Normal
  'RPL_280_neg_1',
  'RPL_280_neg_2',
  'RPL_Epi_1',
  'RPL_Epi_2'#,
  #'cf_distal'
  )

  # Fix the fragment file of multiome sample by removing the header 
  #system ('zcat atac_fragments.tsv.gz | grep -v ^\# | bgzip > atac_fragments_fixed.tzv.gz')

fragment_paths =c(
  '/ahg/regevdata/projects/lungCancerBueno/10x/191121/scATAC_Pt_mesothelioma_CD45_neg_cellranger_atac_v1.2/138_ATACseq_CD45_neg_Lung_ATAC/outs/fragments.tsv.gz',
  '/ahg/regevdata/projects/lungCancerBueno/10x/200128/scATAC_Pt811_mesothelioma_CD45pos_neg_cellranger_atac_v1.2/161_ATACseq_Pt811_mesothelioma_CD45pos_CD45neg_ATAC/outs/fragments.tsv.gz',
  '/ahg/regevdata/projects/lungCancerBueno/10x/200331/10X_Single_Cell_ATAC/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01218_Robby_AlexTsankov/202_ATAC_826CD45pos_826CD45neg/outs/fragments.tsv.gz',  
  '/ahg/regevdata/projects/lungCancerBueno/10x/200721/scATAC_Pt846/10X_Single_Cell_RNA/TD01729_AlexTsankov/846-MesoPool-CD45-pos-CD45-neg-nuclei/outs/fragments.tsv.gz',
  '/ahg/regevdata/projects/lungCancerBueno/10x/200911/scATAC_Pt848/10X_Single_Cell_ATAC/TD01814_AlexTsankov/848/outs/fragments.tsv.gz',
  '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_10/230510/P10_scATAC/cellranger_output/ALTS03_Zhao6ATAC_0_v1/fragments.tsv.gz',
  '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_11/230714/ZHAO8mesotheliomaATAC/cellranger_output/ALTS04_Zhao8ATAC_0_v1/fragments.tsv.gz',
  '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_12/230718/ZHAO9mesotheliomaATAC/cellranger_output/ALTS04_Zhao9ATAC_0_v1/fragments.tsv.gz',
  '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_13/231018/ZHAO12mesotheliomaATAC/cellranger_output/ALTS04_Zhao12ATAC_0_v1/fragments.tsv.gz',
  '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_14/240109/ZHAO13mesotheliomaATAC/cellranger_output/ALTS04_Zhao13ATAC_0_v1/fragments.tsv.gz',
  '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz',
  '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz',
  '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz',
  '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz'#,
  #"/ahg/regevdata/projects/ICA_Lung/10x/200116/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01396_AlexTsankov/aDcd45n/outs/fragments.tsv.gz"
  )
   
  
ArrowFiles_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/ArrowFiles'
ArrowFiles_dir = file.path(ArrowFiles_dir,paste0(sample_names,'.arrow'))
if (!all (paste0(sample_names,'.arrow') %in% list.files(ArrowFiles_dir)))
  {
  #setwd (projdir)  
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
  } else {
  ArrowFiles = paste0(ArrowFiles_dir, paste0(sample_names,'.arrow'))
  }


archp = ArchRProject (
  ArrowFiles = ArrowFiles_dir, 
  outputDirectory = projdir,
  copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

archp$Sample2 = archp$Sample
archp$Sample2[grep ('RPL',archp$Sample2)] = 'normal_pleura'

  
### Subset ArchR object only for cells retained in Signac analysis ####
#tumor_l = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/signac_list.rds')
tumor_annotation = read.csv ('../../git_repo/files/barcode_annotation.csv', row.names=1)
normal = readRDS ('../../per_sample_QC_signac/signac_normal.rds')
normal_annotation = data.frame (barcode= colnames(normal), celltype = normal$predicted.id)
normal_annotation$barcode = gsub (paste0(sample_names[11],'_'),paste0(sample_names[11],'#'),normal_annotation$barcode)
normal_annotation$barcode = gsub (paste0(sample_names[12],'_'),paste0(sample_names[12],'#'),normal_annotation$barcode)
normal_annotation$barcode = gsub (paste0(sample_names[13],'_'),paste0(sample_names[13],'#'),normal_annotation$barcode)
normal_annotation$barcode = gsub (paste0(sample_names[14],'_'),paste0(sample_names[14],'#'),normal_annotation$barcode)

normal_annotation = normal_annotation[normal_annotation$celltype == 'Mesothelium',]
tumor_annotation = tumor_annotation[tumor_annotation$celltype == 'Malignant',]
tumor_annotation$sample = sapply (tumor_annotation$barcode, function(x) unlist(strsplit(x, '#'))[1])
keep_barcodes = c(normal_annotation$barcode, tumor_annotation$barcode)

archp = archp[rownames(archp@cellColData) %in% keep_barcodes]
  
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
archp = addTSNE (ArchRProj = archp, 
  reducedDims = "IterativeLSI",
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

pdf (file.path('Plots','celltype_umap.pdf'))
print (umap_p1)
print (umap_p2)
dev.off()

# Remove outlier cells ####
archp = archp[archp$Clusters != 'C1']

# Save ArchR object ####
archp = saveArchRProject (archp, dropCells = T)
