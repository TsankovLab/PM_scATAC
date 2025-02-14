conda activate meso_scatac
# use UGER # Add this before running R to be able to run cNMF scripts using UGER 

R

set.seed(1234)

# Set project directory
proj_name = '10x_normal_pbmc'
projdir_init = file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/',proj_name)

# Set project dir
# Load utils functions palettes and packages ####
source (file.path(projdir_init,'..','..','git_repo','utils','load_packages.R'))
source (file.path(projdir_init,'..','..','git_repo','utils','useful_functions.R'))
source (file.path(projdir_init,'..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path(projdir_init,'..','..','git_repo','utils','scATAC_functions.R'))
source (file.path(projdir_init,'..','..','git_repo','utils','palettes.R'))

sample_names = c(
	'normal_pbmc')
samples_path = paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/10x_normal_pbmc/filtered_gene_bc_matrices/hg19')


meta = data.frame (
  sampleID = sample_names,
  sample_path = samples_path
  )

#meta = cbind (meta, sample_path = samples_path)
#meta = meta[!duplicated(meta$assignment),]
#samples_path = '/broad/hptmp/bgiotti/JiaMoPool/cellranger_output/RAPA05_JiaPool1_1_v1/raw_feature_bc_matrix/'

# Set cellranger_to_seurat parameters
cr_to_seurat = list(
  run_cellbender = FALSE, # if this is set to TRUE then input raw (and not filtered) cellranger count matrices
  cellbender_dir = NULL,
  cellbender_samples=NULL,
  cellbender_parameters=NULL,
  org = 'human',
  datatype = 'RNA',
  cr_output = 'raw',
  is.hashed = FALSE
  )

# Set QC parameters
qc_params = list(
  filtering = 'hard', # 'emptyDrops' or 'hard' filtering
  nFeat = 400,#400, # Number of features per cells. default 400
  nCounts = 1000, #1000, # Number of UMI per cell. Default 800
  pchM = 25, # Percent mitochondrial genes. Default 25 
  remove.samples = NULL, # Remove bad samples. Takes vector of sampleIDs 
  processInd = FALSE # run preprocessing per sample before running it on the merged data
  )

### Data processing and clustering variables ###
harmony_params = list(
  batch = c('no')
  #batch = 'no'
  )

data_processing_param = list(
  variablefeatures = 'seurat', # options are 'scran', 'seurat','manual'
  nfeat = 3000, # number of variable genes to consider for dimentionality reduction
  sigPCs = 15,
  vars_to_regress=NULL,
  metaGroupNames = c('sampleID'),
  res = c(0.2, 0.8, 2, 5, 10) # denovo cluster resolutions 
  )
    
# Initiate pipeline ####
force=FALSE # re run pipeline from the beginning no matter if objects are found
subclustername = NULL
scrna_pipeline_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scrna_pipeline'
source (file.path(scrna_pipeline_dir,'master_scrna.R'))

# De novo marker discovery ####
enricher_universe = 'all'
logfcThreshold = .25
pvalAdjTrheshold = 0.01
metaGroupName = c('RNA_snn_res.0.8')
top_pathways = 10
top_genes = 5
force = F
source (file.path(scrna_pipeline_dir,'DEG_standard.R'))

srt$celltype = as.character(srt$RNA_snn_res.0.8)
srt$celltype[srt$celltype == '6'] = 'NK'
pdf (file.path('Plots','cell_annotation.pdf'))
DimPlot(srt, reduction = "umap", group.by = 'celltype', label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

saveRDS (srt[, srt$celltype == 'NK'], 'srt_nk.rds')

