conda activate meso_scatac
use UGER
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
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/immune/scatac_ArchR'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

#devtools::install_github("immunogenomics/presto") #needed for DAA
source ('../../PM_scATAC/useful_functions.R')
source ('../../PM_scATAC/ggplot_aestetics.R')
source ('../../PM_scATAC/scATAC_functions.R')
source ('../../PM_scATAC/palettes.R')

set.seed (1234)
addArchRThreads (threads = 1) 
addArchRGenome ("Hg38")

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
    # Normal
    # 'RPL_280_neg_1',
    # 'RPL_280_neg_2',
    # 'RPL_Epi_1',
    # 'RPL_Epi_2'#,
    #'cf_distal'
    )

# Load RNA
#srt = readRDS ('../tumor_compartment/scrna/scRNA_meso.rds')
#sarc_order = read.csv ('../scrna/cnmf20_sarcomatoid_sample_order.csv', row.names=1)

# Load last istance
if (!file.exists ('Save-ArchR-Project.rds'))
   {
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
    '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_14/240109/ZHAO13mesotheliomaATAC/cellranger_output/ALTS04_Zhao13ATAC_0_v1/fragments.tsv.gz'#,#,
    # '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz',
    # '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz',
    # '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz',
    # '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz'#,
    #"/ahg/regevdata/projects/ICA_Lung/10x/200116/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01396_AlexTsankov/aDcd45n/outs/fragments.tsv.gz"
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
  
  bingren_immune_arrows = list.files (file.path ('..','..','tumor_compartment','all_tissues_ArchR','normal_controls_immune','ArrowFiles'))
  ArrowFiles = append (ArrowFiles, (file.path ('..','..','tumor_compartment','all_tissues_ArchR','normal_controls_immune','ArrowFiles',bingren_immune_arrows)))
  archp = ArchRProject (
    ArrowFiles = ArrowFiles, 
    outputDirectory = projdir,
    copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  )
  
  ### Subset ArchR object only for immune PM cells retained in Signac analysis ####
  archp_main = loadArchRProject ('../../main/scatac_ArchR/')
  archp$celltype = 0
  archp$celltype[match (rownames(archp_main@cellColData), rownames(archp@cellColData))] = archp_main$celltype
  archp_tumor = loadArchRProject ('../../tumor_compartment/scatac_ArchR/')
  archp_tumor = archp_tumor[rownames(archp_tumor@cellColData) %in% rownames(archp@cellColData)]
  archp$celltype[match (rownames(archp_tumor@cellColData), rownames(archp@cellColData))] = archp_tumor$Sample2
  archp = archp[archp$celltype != 'Malignant']

  # Pass annotation from bingren normal ####
  archp_imm = read.csv ('../../tumor_compartment/all_tissues_ArchR/normal_controls_immune/archp_meta.csv')
  archp$celltype[match (archp_imm$X, rownames(archp@cellColData))] = archp_imm$celltype_project
  archp$celltype[match (archp_imm$X, rownames(archp@cellColData))] = archp_imm$celltype
  archp$Sample2 = archp$Sample
  archp$Sample2[match (archp_imm$X, rownames(archp@cellColData))] = sapply (archp_imm$celltype_project, function(x) unlist(strsplit(x,'_'))[6])
  archp = archp [archp$celltype != 0]
  archp = archp [archp$celltype %in% c(
    'AlverolarType2,Immune',
    'B_cells',
    'Macrophage(General)',
    'Macrophage(General,Alveolar)',
    'MonoMac','NK','NaiveTcell','NaturalKillerTCell',
    'Plasma','PlasmaCell','TLymphocyte1(CD8+)','T_cells','Tlymphocyte2(CD4+)')]

  #archp = archp[!is.na(archp$celltype)]

  # Dimensionality reduction and clustering
  varfeat = 25000
  LSI_method = 2
  archp = addIterativeLSI (ArchRProj = archp,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = FALSE, LSIMethod = LSI_method,
    varFeatures = varfeat)

  archp = addHarmony (
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample2", force=FALSE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony",
    name='Clusters_H',
    force = TRUE)

  archp = addClusters (input = archp, resolution = 3,
    reducedDims = "IterativeLSI", maxClusters = 100,
    force = TRUE)
  archp = addUMAP (ArchRProj = archp, 
    reducedDims = "IterativeLSI",
    force = TRUE)
  # archp = addTSNE (ArchRProj = archp, 
  #   reducedDims = "IterativeLSI",
  #   force = TRUE)
  
  #archp = saveArchRProject (archp)

  umap_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
   name = "celltype", embedding = "UMAP")
  umap_p2 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample2",
     embedding = "UMAP")
  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample2",
     embedding = "UMAP_H")
  umap_p4 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "celltype",
     embedding = "UMAP_H")
  umap_p5 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters_H",
     embedding = "UMAP_H")
  
  pdf (file.path('Plots','celltype_umap_signac_filtered.pdf'),5,5)
  print (umap_p1)
  print (umap_p2)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()


## harmonize project object ####
archp$project = ifelse (grepl ('^P', archp$Sample2), 'PM','normal')

  # Dimensionality reduction and clustering
  varfeat = 25000
  LSI_method = 2
  archp = addIterativeLSI (ArchRProj = archp,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = FALSE, LSIMethod = LSI_method,
    varFeatures = varfeat)

  archp = addClusters (input = archp, resolution = 3,
    reducedDims = "IterativeLSI", maxClusters = 100,
    force = TRUE)
  archp = addUMAP (ArchRProj = archp, 
    reducedDims = "IterativeLSI",
    force = TRUE)

  archp = addHarmony (
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    name = "Harmony_project",
    groupBy = c('project'), force=TRUE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony_project", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony_project",
    name='Clusters_H',
    force = TRUE)

  # archp = addTSNE (ArchRProj = archp, 
  #   reducedDims = "IterativeLSI",
  #   force = TRUE)
  
  #archp = saveArchRProject (archp)


  umap_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
   name = "celltype", embedding = "UMAP")
  umap_p2 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample2",
     embedding = "UMAP")
  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample2",
     embedding = "UMAP_H")
  umap_p4 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "celltype",
     embedding = "UMAP_H")
  umap_p5 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters_H",
     embedding = "UMAP_H")
  
  pdf (file.path('Plots','celltype_umap_signac_filtered_harmony_on_project.pdf'),5,5)
  print (umap_p1)
  print (umap_p2)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()

## Clean object ####
archp = archp[!archp$Clusters_H %in% c('C5','C1')]
archp = archp[!archp$celltype %in% c('AlverolarType2,Immune')]

  archp = addHarmony (
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    name = "Harmony_project",
    groupBy = c('project'), force=TRUE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony_project", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony_project",
    name='Clusters_H',
    force = TRUE)

  # archp = addTSNE (ArchRProj = archp, 
  #   reducedDims = "IterativeLSI",
  #   force = TRUE)
  
  #archp = saveArchRProject (archp)

  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Sample2",
     embedding = "UMAP_H")
  umap_p4 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "celltype",
     embedding = "UMAP_H")
  umap_p5 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters_H",
     embedding = "UMAP_H")
  
  pdf (file.path('Plots','celltype_umap_signac_filtered_harmony_on_project_cleaned.pdf'),5,5)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()



# Immune markers ####
meso_markers = c('PTPRC','CD3E','CD3D','CD3G','FOXP3','CD8A','CD4','GNLY','LYZ','C1QA')
archp = addImputeWeights (archp)
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = meso_markers, 
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)

#archp$celltype[archp$Clusters == 'C30'] = 'Fibroblasts_WT1'
#p = lapply (p, function(x) x + theme_void() + NoLegend ()) #+ ggtitle scale_fill_gradient2 (rev (viridis::plasma(100))))

pdf (file.path('Plots','TNK_markers_fplots_cleaned.pdf'), width = 25, height = 25)
print (p)
dev.off()


# Subset T cells ####
metaGroupName = 'celltype'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c(
    'T_cells','Tlymphocyte2(CD4+)','TLymphocyte1(CD8+)','NaturalKillerTCell','NaiveTcell','NK')],
  outputDirectory = file.path('..','..','NKT_cells'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)


# Subset Myeloid cells ####
metaGroupName = 'celltype'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c(
    'Macrophage(General)','Macrophage(General,Alveolar)','MonoMac')],
  outputDirectory = file.path('..','..','myeloid_cells'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)





