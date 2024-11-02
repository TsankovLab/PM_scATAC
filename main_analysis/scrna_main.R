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
  'BSgenome.Hsapiens.UCSC.hg38',
  'tidyverse',
  'ggrepel',
  'RColorBrewer')
lapply(packages, require, character.only = TRUE)

####### ANALYSIS of TUMOR compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scrna'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)


#devtools::install_github("immunogenomics/presto") #needed for DAA
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))
source (file.path('..','..','git_repo','utils','hubs_track.R'))

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
    'P14'#,# p14
    # # Normal
    # 'RPL_280_neg_1',
    # 'RPL_280_neg_2',
    # 'RPL_Epi_1',
    # 'RPL_Epi_2'#,
    # #'cf_distal'
    )

# Load RNA
srt = readRDS ('../scrna/srt.rds')
srt$celltype_simplified2[srt$celltype_simplified2 == 'pDC'] = 'pDCs'
sarc_order = read.csv ('../scrna/cnmf20_sarcomatoid_sample_order.csv', row.names=1)

pdf (file.path ('Plots','pan_TFs_dotplot.pdf'))
DotPlot (srt, features = c('RUNX2','NR4A2','TGFBI','TGIF1','TGIF2'))
dev.off()

pdf (file.path ('Plots','pan_TFs_vlnplot.pdf'), width=20)
VlnPlot (srt, features = c('RUNX2','NR4A2','TGFBI','TGIF1','TGIF2'), group.by = 'celltype2')
dev.off()




### Run cNMF per compartment ####

celltypes = names(table (srt$celltype_simplified))[table (srt$celltype_simplified) > 300]

srt_main = srt
for (ct in celltypes)
  {
  #### Run cNMF ####
  nfeat = 5000
  force=F
  k_list = c(20:30)
  k_selections = c(20:30)
  cores= 100
  cnmf_name = ct
  srt = srt_main[,srt_main$celltype_simplified == ct]

  repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
  
  ### RUN consensus NMF ####
  source (file.path ('..','..','git_repo','utils','cnmf_prepare_inputs.R')) 
  }


### Import and format spectra files ####
k_selection = 20
cnmf_name = 'Malignant'
source (file.path ('..','..','git_repo','utils','cnmf_format_spectra_files.R')) 
sapply (cnmf_spectra_unique, function(x) 'KLRC1' %in% x)
cnmf_spectra_unique[[6]]








