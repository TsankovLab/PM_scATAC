conda activate meso_scatac
# use UGER # Add this before running R to be able to run cNMF scripts using UGER 

R

set.seed(1234)

# Set project directory
proj_name = 'NK_cell_line_AH'
projdir_init = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/NK_cell_line_scrna_AH'

# Set project dir
# Load utils functions palettes and packages ####
source (file.path(projdir_init,'..','..','git_repo','utils','load_packages.R'))
source (file.path(projdir_init,'..','..','git_repo','utils','useful_functions.R'))
source (file.path(projdir_init,'..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path(projdir_init,'..','..','git_repo','utils','scATAC_functions.R'))
source (file.path(projdir_init,'..','..','git_repo','utils','palettes.R'))

sample_names = c(
	'AMHO18_199untreated_0_v1',
	'AMHO18_247AZA_0_v1',
	'AMHO18_247untreated_0_v1',
	'AMHO18_199AZA_0_v1',
	'AMHO18_245untreated_0_v1',
	'AMHO18_245AZA_0_v1')
samples_path = paste0('/sc/arion/projects/nmibc_bcg/NK_azacitidine/5GEX/RAIscSeq/cellranger_output/', sample_names,'/filtered_feature_bc_matrix.h5')


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


# # Get highlevel markers for doublets identification ####
# markers = read.csv (paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/gene_sets/highlevel_MPM_markers.csv'))
# colnames (markers)[1] = 'gene'
# markers_name = 'doublets_markers'
# metaGroupName = paste0(reductionGraphSnn, '_res.10')
# source (file.path(scrna_pipeline_dir, 'markers&clAssigner.R'))

# # De novo marker discovery ####
# enricher_universe = 'all'
# logfcThreshold = .25
# pvalAdjTrheshold = 0.01
# metaGroupName = c('RNA_snn_res.2')
# top_pathways = 10
# top_genes = 5
# force = F
# source (file.path(scrna_pipeline_dir,'DEG_standard.R'))

# Plot NR4A2 and other markers of dysfunction
genes = c('NR4A2','TIGIT','HAVCR2','EGR2','TOX','TOX2','GZMA','GMZH','GZMB','PRF1','GNLY','KLRC1','CD8A','CD8B')
pdf (file.path ('Plots','NK_markers.pdf'))
wrap_plots (fp (srt, genes))
dev.off()

celline = setNames (c('AMHO18_199','AMHO18_247','AMHO18_247','AMHO18_199','AMHO18_245','AMHO18_245'),unique(srt$orig.ident))
condition = setNames (c('untreated','AZA','untreated','AZA','untreated','AZA'),unique(srt$orig.ident))
srt$line = celline[srt$orig.ident]
srt$condition = condition[srt$orig.ident]
srt$line_condition = paste0(srt$line,'_', srt$condition)
srt$sampleID = paste0('NKCL_',srt$line_condition)

# Import NK cells in meso tumor
srt_meso = readRDS ('../../../scrna/srt.rds')
srt_meso = srt_meso[,srt_meso$celltype2 %in% c('NK_KLRC1','NK_FGFBP2')]
srt_meso$sampleID = paste0('MESO_',srt_meso$sampleID)

# #srt$line_condition[is.na(srt$line_condition)] = 0
# srt_merged$celltype2[is.na(srt_merged$celltype2)] = paste0('NKCL_',srt_merged$RNA_snn_res.2[is.na(srt_merged$celltype2)])
# srt_merged$celltype2[srt_merged$RNA_snn_res.2 =='17'] = 'NKCL_NR4A2+'

genes = c('NR4A2','TIGIT','HAVCR2','EGR2','GZMA','GMZH','GZMB','PRF1','GNLY','NKG7','KLRC1','CD8A','CD8B')
# pdf (file.path ('Plots','NK_markers_dotplots.pdf'), width=7, height=6)
# #VlnPlot (srt_merged, feature = genes, group.by = 'celltype2')
# DotPlot (srt_merged, feature = genes, group.by = 'celltype2') + gtheme
# DotPlot (srt_merged[,srt_merged$RNA_snn_res.2 == '17'], feature = genes, group.by = 'orig.ident') + gtheme
# dev.off()

# Split by cell line and condition 
# srt_merged$line_condition[is.na(srt_merged$line_condition)] = srt_merged$sampleID[is.na(srt_merged$line_condition)]
# pdf (file.path ('Plots','NK_markers_sample_split_dotplots.pdf'), width=7, height=6)
# #VlnPlot (srt_merged, feature = genes, group.by = 'celltype2')
# DotPlot (srt_merged, feature = genes, group.by = 'line_condition') + gtheme
# dev.off()

# write.csv (as.data.frame(srt$RNA_snn_res.2), 'barcode_annotation.csv')

# Import PBMC from meso patient
srt_pbmc = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/MPM_naive_study/MPM_naive_PBMC_CITEseq2_analysis/_cellranger_filtered_Filter_400_1000_10/sampleID_harmony_cc_nCount_RNA_regressed/srt.rds')
srt_pbmc = srt_pbmc[,srt_pbmc$celltype_simplified == 'NK']
srt_pbmc$sampleID = paste0('pbmc_MESO_', srt_pbmc$sampleID)

# pdf (file.path ('Plots','NK_pbmc_meso_markers_dotplots.pdf'), width=7, height=6)
# DotPlot (srt_pbmc, feature = genes, group.by = 'celltype_simplified') + gtheme
# #DotPlot (srt_pbmc[,srt_merged$RNA_snn_res.2 == '17'], feature = genes, group.by = 'orig.ident') + gtheme
# dev.off()

# pdf (file.path ('Plots','NK_pbmc_meso_markers.pdf'))
# reductionName = 'sampleID_harmony_umap'
# wrap_plots (fp (srt_pbmc, genes))
# dev.off()

# Import NK cells from normal PBMC
srt_pbmc_normal = readRDS('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/10x_normal_pbmc/_cellranger_raw_Filter_400_1000_25/no_harmony/srt_nk.rds')
srt_pbmc_normal$sampleID = 'normal_pbmc'

# Combine datasets
srt_merged = merge (srt, srt_meso)
srt_merged = merge (srt_merged, srt_pbmc)
srt_merged = merge (srt_merged , srt_pbmc_normal)

genes = c('NR4A2','TIGIT','HAVCR2','EGR2','GZMA','GMZH','GZMB','PRF1','GNLY','NKG7','KLRC1','CD8A','CD8B')
pdf (file.path ('Plots','NK_markers_all_merged_dotplots.pdf'), width=12, height=6)
#VlnPlot (srt_merged, feature = genes, group.by = 'celltype2')
DotPlot (srt_merged, feature = genes, group.by = 'sampleID') + gtheme
dev.off()

genes = c('NR4A2')#,'TIGIT','HAVCR2','EGR2','GZMA','GMZH','GZMB','PRF1','GNLY','NKG7','KLRC1','CD8A','CD8B')
pdf (file.path ('Plots','NR4A2_all_merged_dotplots.pdf'), width=6, height=6)
#VlnPlot (srt_merged, feature = genes, group.by = 'celltype2')
DotPlot (srt_merged, feature = genes, group.by = 'sampleID') + gtheme
dev.off()





# #### Run cNMF ####
# nfeat = 5000
# force=F
# k_list = c(5:30)
# k_selections = c(5:30)
# cores= 100
# cnmf_name = 'TNK'
# cnmf_out = paste0('cNMF/cNMF_',cnmf_name,'_',paste0(k_list[1],'_',k_list[length(k_list)]),'_vf',nfeat)
# dir.create (file.path(cnmf_out,'Plots'), recursive=T)
# repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'

# ### RUN consensus NMF ####
# source (file.path ('..','..','git_repo','utils','cnmf_prepare_inputs.R')) 


# ### Import and format spectra files ####
# k_selection = 19
# source (file.path ('..','..','git_repo','utils','cnmf_format_spectra_files.R')) 
# sapply (cnmf_spectra_unique, function(x) 'KLRC1' %in% x)
# cnmf_spectra_unique[[6]]

