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
# srt_p14 = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/P14_analysis/_cellranger_raw_Filter_400_1000_25/no_harmony/srt.rds')
# srt_p14 = subset(x = srt_p14, features = rownames(srt))
# srt2 = SplitObject (srt, split.by = 'sampleID')
# srt_p14$celltype_lv1 = srt_p14$celltype
# srt_p14$celltype_lv1[srt_p14$celltype_lv1 == 'LEC'] = 'Endothelial'
# srt_p14$celltype_lv1[srt_p14$celltype_lv1 == 'Smooth Muscle'] = 'SmoothMuscle'
# srt_p14 = srt_p14[,srt_p14$celltype_lv1 != 'NA']
# srt_p14 = srt_p14[,srt_p14$celltype_lv1 != 'doublets']
# srt = merge (x = srt_p14, y= srt2)
# srt[["RNA"]] = JoinLayers(srt[["RNA"]])



# Initiate pipeline
scrna_pipeline_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scrna_pipeline'

org = 'human'

### Data processing and clustering variables ###
batch = 'no'
variablefeatures = 'scran'
nfeat = 2000 # number of variable genes to consider for dimentionality reduction
sigPCs = 15
ccRegress = T # Regress cell cycle gene expression 
vars_to_regress=NULL
metaGroupNames = c('sampleID','celltype_lv1')
res = c(2) # denovo cluster resolutions 

source (file.path(scrna_pipeline_dir,'data_processing.R'))
source (file.path(scrna_pipeline_dir,'harmony_and_clustering.R'))

srt_tnk = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scrna/srt.rds')
srt_stroma = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/stroma/scrna/srt.rds')
srt_mye = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/srt.rds')
srt_tumor = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scrna/srt.rds')
srt_B = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/B_cells/scrna/srt.rds')


srt$celltype_lv2 = srt$celltype_lv1
srt$celltype_lv2[match(colnames(srt_tnk), colnames(srt))] = srt_tnk$celltype_lv2
srt$celltype_lv2[match(colnames(srt_stroma), colnames(srt))] = srt_stroma$celltype_lv2
srt$celltype_lv2[match(colnames(srt_mye), colnames(srt))] = srt_mye$celltype_lv2
srt$celltype_lv2[match(colnames(srt_B), colnames(srt))] = srt_B$celltype_lv2
srt$celltype_lv2[match(colnames(srt_tumor), colnames(srt))] = 'Malignant'

table (srt$celltype_lv1, useNA = 'always')
table (srt$celltype_lv2, useNA = 'always')

# Add normal samples but add tag to avoid using them in compartment subsets except for tumor analysis
srt_normal = srt_tumor[,srt_tumor$sampleID %in% c('HU37','HU62')]
srt_light = merge (x = srt_normal, y= srt)
srt_light[["RNA"]] = JoinLayers(srt_light[["RNA"]])

srt_light$celltype_lv1[is.na(srt_light$celltype_lv1)] = 'normal_mesothelium'
srt_light$celltype_lv2[is.na(srt_light$celltype_lv2)] = 'normal_mesothelium'
table (srt_light$celltype_lv1, useNA = 'always')
table (srt_light$celltype_lv2, useNA = 'always')


coldata_to_keep = c('sampleID','celltype_lv1','celltype_lv2')
srt_light@meta.data = srt_light@meta.data[,coldata_to_keep]

saveRDS (srt_light,'srt_light.rds')


table (colnames (srt_tnk) %in% colnames (srt),srt_tnk$sampleID)
table (colnames (srt_tumor[,!srt_tumor$sampleID %in% c('HU37','HU62')]) %in% colnames (srt))
table (colnames (srt_stroma) %in% colnames (srt))

saveRDS (srt, 'srt.rds')

head (srt)

pdf (file.path ('Plots','celltypes_lv1_umap.pdf'))
DimPlot (srt, group.by = 'celltype_lv1')
DimPlot (srt, group.by = 'sampleID')
dev.off()



srt$celltype_simplified2[srt$celltype_simplified2 == 'pDC'] = 'pDCs'
srt = srt[,srt$sampleID %in% sample_names]
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









scatac_markers = c('ADAM29
MAST4
ALPL
ST14
FBLN7
RBM47
MIR4802
TNFSF9
PCDH8
TRPV3
ZBTB20−AS5
COL4A3
ZBTB20
ATXN1
CHRNA1
LINC00944
TRAM2
ENO4
OXCT2
CDYL2
ZBTB17
LRRK1
KLHL6
ZNF532
SLC23A2
LINC01814
ITFG1−AS1
LPCAT1
TNFRSF13B
PALD1
PAQR8
LOC285626
CD80
LINC01845
PCDHGB2
PCDHGB1
PCDHGA4
PCDHGA2
PCDHGA3
PCDHGA1
PTPN3
PCDHGA7
PCDHGA9
PCDHGB5
PCDHGA8
PCDHGB4
PCDHGB3
PCDHGA5
PCDHGA6
PCDHGB6
PCDHGA10
PCDHGB7
PCDHGA11
UBE2QL1
TOX
SHANK1
OSBPL1A
INSR
PCDH9
SNTG2
MYADM
CERS6
SATB1−AS1
MIR4265
THRB−AS1
GAS7
SPRY1
YBX3
MAML3
SH3RF3
MIR4266
SPRY4−AS1
ZNF667
ZNF667−AS1
TCL1A
TCL1B
IL4R
TCL6
MIR3139
GAB1
ZBTB16
IL21R−AS1
IL21R
MSI2
FOXP1
ZNF516
C18orf65')
scatac_markers = unlist (strsplit(scatac_markers, '\\\n'))
pdf (file.path ('Plots','Bcells_protocaderins_markers.pdf'), width=15)
DotPlot (srt, feature = scatac_markers) + gtheme
dev.off()
