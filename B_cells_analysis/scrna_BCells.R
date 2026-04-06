conda activate scrnatools
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
  'RColorBrewer',
  'SingleCellExperiment',
  'scran',
  'harmony')
lapply(packages, require, character.only = TRUE)

####### ANALYSIS of TUMOR compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/B_cells/scrna'
scrna_pipeline_dir = file.path('..','..','git_repo','utils')
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
if (!file.exists ('srt.rds'))
{
srt = readRDS ('../../main/scrna/srt.rds')
srt = srt[,srt$celltype_simplified %in% c('B_cells','Plasma')]
srt = srt[,srt$sampleID %in% sample_names]
srt = srt[, srt$sampleID %in% names (table (srt$sampleID)[table (srt$sampleID) > 20])]

srt = NormalizeData (srt)
sce = SingleCellExperiment (list(counts=srt@assays$RNA@layers$counts, logcounts = srt@assays$RNA@layers$data),
rowData=rownames(srt)) 
rownames(sce) = rownames(srt)
sce = modelGeneVar (sce)
# remove batchy genes
batchy_genes = c('RPL','RPS','MT-')
sce = sce[!apply(sapply(batchy_genes, function(x) grepl (x, rownames(sce))),1,any),]
nfeat = 5000
vf = getTopHVGs (sce, n=nfeat)

# remove ambient RNA genes
vf = vf[!grepl ('SFTP',vf)]
vf = vf[!grepl ('HB',vf)]
VariableFeatures (srt) = vf
srt <- ScaleData (srt)
srt <- RunPCA (srt)

batch = 'sampleID'
reductionSave = 'harmony'
srt = srt %>% 
RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
RunUMAP (reduction = reductionSave, dims = 1:15)

srt = FindNeighbors (object = srt, reduction = reductionSave, dims = 1:15, k.param = 30,
                            verbose = TRUE)
srt = FindClusters (srt, resolution = 1, verbose = T, n.start = 100)

dp = DimPlot (srt, group.by = 'seurat_clusters', reduction='umap')
dp2 = DimPlot (srt, group.by = 'sampleID', reduction='umap')
dp3 = DimPlot (srt, group.by = 'celltype_simplified', reduction='umap')

pdf (file.path('Plots','seurat_clusters_umap.pdf'), width=12)
print (wrap_plots (dp, dp2,dp3))
dev.off()
saveRDS (srt, 'srt.rds')
} else {
srt = readRDS ('srt.rds')  
}

# # De novo marker discovery ####
org='human'
enricher_universe = 'all'
logfcThreshold = .25
pvalAdjTrheshold = 0.01
metaGroupName = c('seurat_clusters')
top_pathways = 10
top_genes = 5
force = T
source (file.path(scrna_pipeline_dir,'DEG_standard.R'))

reductionName = 'umap'
fps = fp (srt, gene = c('CD79B','CD79A','CD19','CD37','RGS13','MEF2B','LMO2','IGHG3','IGHGP','IGHG4',
  ))
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
reductionName = 'umap'
fps2 = fp (srt, gene = scatac_markers)

pdf (file.path('Plots','Bcell_markers2.pdf'),20,20)
wrap_plots (fps2)
dev.off()
  



