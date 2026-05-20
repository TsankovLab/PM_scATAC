ml anaconda3
conda activate meso_scatac
R

# Load packages
packages = c('ArchR','rGREAT','igraph','karyoploteR','Seurat','ggpubr','org.Hs.eg.db',
  'TxDb.Hsapiens.UCSC.hg38.knownGene','EnsDb.Hsapiens.v86','gplots','ggrepel','patchwork','ComplexHeatmap',
  'RColorBrewer','ggrepel','gplots', 'rstatix','dplyr','tidyr')
lapply(packages, require, character.only = TRUE)

#devtools::install_github("immunogenomics/presto") #needed for DAA
#source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/useful_functions.R')
meso_repo = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
source (file.path (meso_repo, 'utils','useful_functions.R'))
source (file.path (meso_repo, 'utils','scATAC_functions.R'))

# Load functions for hub detection
# source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/cooltools/hubs_tools/knnGen.R')
# source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/cooltools/hubs_tools/addCoAx.R')
# source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/cooltools/hubs_tools/Hubs_finder.R')
# source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/cooltools/hubs_tools/hubs_track.R')
# source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/scrna_pipeline/useful_functions.R')
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/pbmc_myeloid'
setwd (projdir)


set.seed(1234)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

# Load utils functions palettes and packages ####
source (file.path('..','..','..','mesothelioma','scATAC_PM','git_repo','utils','load_packages.R'))
source (file.path('..','..','..','mesothelioma','scATAC_PM','git_repo','utils','useful_functions.R'))
source (file.path('..','..','..','mesothelioma','scATAC_PM','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','..','mesothelioma','scATAC_PM','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','..','mesothelioma','scATAC_PM','git_repo','utils','palettes.R'))


#meta = meta[meta$ASSAY.Type == 'ATAC',]
#meta = meta[meta$tissue %in% c('Normal','Tumor'),]
#meta = meta[meta$Project == 'lung',]
#fragmentsDir = '/broad/hptmp/bgiotti/AS_human_lung_scATAC/fragment_files/'
#fragment_filenames = list.files (fragmentsDir, pattern = '*gz$')
#sampleID = sapply (fragment_filenames, function(x) unlist(strsplit(x, '_'))[1])

#meta = meta[meta$sample_ID %in% sampleID, ]
#meta = meta[meta$tissue %in% c('Normal','Tumor'),]
if (!file.exists('Save-ArchR-Project.rds'))
  {

#fragment_filenames = fragment_filenames[!grepl('.tbi',fragment_filenames)]
#fragment_filenames = fragment_filenames[sapply (meta$sample_ID, function(x) grep (x, fragment_filenames))]
#fragment_filepaths = paste0(fragmentsDir, fragment_filenames)
meta = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/sample_annots.csv')
meta = meta[meta$library_chemistry == 'V0_multiome',]
meta = meta[meta$tissue %in% c('Normal','Tumor','PBMC'), ]
patient_ID = unique (meta$patient_ID)

# Copy arrows of multiome from NT_myeloid3 and PBMC_myeloid projects
nt_myeloid_path = '/ahg/regevdata/projects/ICA_Lung/Bruno/AS_human_lung_scatac/analysis/NT_myeloid3/ArrowFiles/'
pbmc_myeloid_path = '/ahg/regevdata/projects/ICA_Lung/Bruno/AS_human_lung_scatac/analysis/pbmc_myeloid/ArrowFiles/'
nt_myeloid_arrows = c(paste0(patient_ID, c('N.arrow')), paste0(patient_ID, c('T.arrow')))
pbmc_myeloid_arrows = paste0(patient_ID, 'P.arrow')
if (!all(file.exists(c(paste0(projdir,'/ArrowFiles/',nt_myeloid_arrows),paste0(projdir,'ArrowFiles/',pbmc_myeloid_arrows))))
list.files (pbmc_myeloid_path)
lapply (paste0(nt_myeloid_path, nt_myeloid_arrows), function(x) system (paste0('cp ', x, ' ',projdir,'ArrowFiles/')))
lapply (paste0(pbmc_myeloid_path, pbmc_myeloid_arrows), function(x) system (paste0('cp ', x, ' ',projdir,'ArrowFiles/')))

####### START ANALYSIS #######
system (paste('mkdir -p',projdir))
setwd (projdir)

  fragment_filepaths = paste0('/broad/hptmp/bgiotti/samarth_fragment_files/',list.files('/broad/hptmp/bgiotti/samarth_fragment_files/',pattern='gz$'))
  sample_names = list.files('/broad/hptmp/bgiotti/samarth_fragment_files/',pattern='gz$')
  ArrowFiles = createArrowFiles (inputFiles = fragment_filepaths,
    sampleNames = sample_names,
    minTSS = 6, #Dont set this too high because you can always increase later
    minFrags = 1000,
    maxFrags = 100000,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    force = TRUE
  )  

archp = ArchRProject (
  ArrowFiles = ArrowFiles, 
  outputDirectory = projdir,
  copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

archp = saveArchRProject (archp)
} else {
archp = loadArchRProject (projdir)
}

### Run peak calling ####
metaGroupName = "celltype_simplified"
force=F
peak_reproducibility=2
if(!all(file.exists(file.path('PeakCalls', paste0(unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds')))) | force) source (file.path(meso_repo,'utils','callPeaks.R'))

# pbmc_ann = read.csv ('../pbmc_myeloid/barcode_celltype_simplified.csv')
# tn_ann = read.csv ('../NT_myeloid3/barcode_celltype_simplified.csv')
# ann_df = rbind (pbmc_ann, tn_ann)
# barcodes = rownames (archp@cellColData)
# barcodes = sub ('_fragments.tsv.gz','',barcodes)
# archp@cellColData = cbind (archp@cellColData, ann_df[match(barcodes, ann_df$barcode),])
# table (barcodes %in% tn_ann$barcode)
# barcode_myeloid = which(barcodes %in% ann_df$barcode)
# archp = archp[barcode_myeloid]

# archp$compartment = 0
# archp$compartment[grep ('P',archp$Sample)] = 'PBMC'
# archp$compartment[grep ('T',archp$Sample)] = 'Tumor'
# archp$compartment[grep ('N',archp$Sample)] = 'Normal'
metaGroupNames = c('celltype_simplified','Sample')
pdf(file.path('Plots','celltype_umap.pdf'))
umap_p1 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
name = x, embedding = "UMAP_H"))
umap_p1
dev.off()

# # Check TF deviations
TF = 'FOS'
TF = getFeatures (archp, 'MotifMatrix')[grep (TF, getFeatures (archp, 'MotifMatrix'))]
TF1 = c('z:JUN_143','z:FOS_137','z:NFKB1_719')

pdf ()
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "MotifMatrix",
    name = TF1, 
    useSeqnames='z',
    pal = rev (palette_deviation),    
    embedding = "UMAP_H",
    imputeWeights = NULL
    )
dev.off()

pdf(file.path('Plots','AP1_deviation_fplot.pdf'),7,7)
wrap_plots (TF_p)
#wrap_plots(umap_p1)
#wrap_plots (umap_p0,umap_p2)#,umap_p3)
dev.off()


# Make TF activity average per cell type 
metaGroupName = 'celltype'
mMats = getGroupSE(
  ArchRProj = archp[archp$compartment == 'PBMC'],
  useMatrix = 'MotifMatrix',
  groupBy = metaGroupName,
  divideN = TRUE,
  scaleTo = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)

mmat = assays(mMats)[[1]]
rownames (mmat) = rowData (mMats)$name
rownames(mmat) = gsub ('_.*','',rownames(mmat))
rownames(mmat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", rownames(mmat))

write.csv (mmat, 'average_TF_activity_PBMC2.csv')



### Make groups of peakmatrix 
metaGroupName = 'celltype_simplified'
pMats = getGroupSE(
  ArchRProj = archp,
  useMatrix = 'PeakMatrix',
  groupBy = metaGroupName,
  divideN = TRUE,
  scaleTo = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)

mmat = assays(mMats)[[1]]

MDM_TFs = c('FOSL1','FOSL2','BACH1','PPARG','NFKB2','KLF12','HIVEP3','SMAD1','NFKB1','REL','RUNX1','SNAI1','RUNX2','NFAT5')


# # Check TF deviations
TF = 'FOS'
TF = getFeatures (archp, 'MotifMatrix')[unlist(sapply (MDM_TFs, function(x) grep (x, getFeatures (archp, 'MotifMatrix'))))]
TF = TF[grep ('z:',TF)]
#TF1 = c('z:JUN_143','z:FOS_137','z:NFKB1_719')

pdf ()
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "MotifMatrix",
    name = TF, 
    useSeqnames='z',
    pal = rev (palette_deviation),    
    embedding = "UMAP_H",
    imputeWeights = NULL
    )
dev.off()

pdf(file.path('Plots','AP1_deviation_fplot.pdf'),14,14)
wrap_plots (TF_p)
dev.off()


# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames(mMat) = rowData(mSE)$name


### Make groups of peakmatrix 
metaGroupName = 'celltype_simplified'
mMats = getGroupSE(
  ArchRProj = archp,
  useMatrix = 'MotifMatrix',
  groupBy = metaGroupName,
  divideN = TRUE,
  scaleTo = NULL,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)
mMat = assays(mMats)[[1]]
rownames(mMat) = rowData (mMats)$name
rownames(mMat) = gsub ('_.*','',rownames(mMat))
rownames(mMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", rownames(mMat))

mMat = mMat[MDM_TFs, ]
#mMat = mMat[,c(7:10,12,13,15)]

pdf(file.path('Plots','MDM_TFs_heatmap.pdf'),width = 3,5)
Heatmap (t(scale(t(mMat))), column_names_rot = 45, row_names_side = 'left',
  col = rev(palette_deviation), border = T)
dev.off()











