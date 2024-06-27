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
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/stroma/scatac_ArchR'
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
    # Normal
    'RPL_280_neg_1',
    'RPL_280_neg_2',
    'RPL_Epi_1',
    'RPL_Epi_2'#,
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
    '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_14/240109/ZHAO13mesotheliomaATAC/cellranger_output/ALTS04_Zhao13ATAC_0_v1/fragments.tsv.gz',#,
    '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz',
    '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz',
    '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz',
    '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz'#,
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
      
  archp = ArchRProject (
    ArrowFiles = ArrowFiles, 
    outputDirectory = projdir,
    copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  )
  
  ### Subset ArchR object only for cells retained in Signac analysis ####
  archp_main = loadArchRProject ('../../main/scatac_ArchR/')
  archp$celltype = 0
  archp$celltype[match (rownames(archp_main@cellColData), rownames(archp@cellColData))] = archp_main$celltype
  archp_tumor = loadArchRProject ('../../tumor_compartment/scatac_ArchR/')
  archp$celltype[match (rownames(archp_tumor@cellColData), rownames(archp@cellColData))] = archp_tumor$Sample2
  archp = archp[archp$celltype != 'Malignant']

  # Use Annotation from signac
  normal = readRDS ('../../per_sample_QC_signac/signac_normal.rds')
  normal_annotation = data.frame (barcode= colnames(normal), celltype = normal$predicted.id)
  normal_annotation$barcode = gsub (paste0(sample_names[11],'_'),paste0(sample_names[11],'#'),normal_annotation$barcode)
  normal_annotation$barcode = gsub (paste0(sample_names[12],'_'),paste0(sample_names[12],'#'),normal_annotation$barcode)
  normal_annotation$barcode = gsub (paste0(sample_names[13],'_'),paste0(sample_names[13],'#'),normal_annotation$barcode)
  normal_annotation$barcode = gsub (paste0(sample_names[14],'_'),paste0(sample_names[14],'#'),normal_annotation$barcode)
  normal_annotation = normal_annotation[normal_annotation$barcode %in% rownames(archp@cellColData),]
  archp$celltype[match (normal_annotation$barcode, rownames(archp@cellColData))] = normal_annotation$celltype
  
  archp = archp[archp$celltype != 0]
  archp = archp[archp$celltype != 'bad_quality']
  archp = archp[archp$celltype %in% names(table (archp$celltype)[table(archp$celltype) > 10])]
  archp$celltype[archp$celltype == 'Fibroblast'] = 'Fibroblasts'
  archp$celltype[archp$celltype == 'Myeloid'] = 'MonoMac'
  archp$celltype[archp$celltype == 'B.cells'] = 'B_cells'
  
  archp = archp[!archp$celltype %in% c('AT1','AT2','Alveolar','B_cells','Ciliated','Mast','MonoMac','NK','Plasma','Secretory','T.NK.cells','T_cells','cDCs','pDCs')]

  archp$Sample2 = archp$Sample
  archp$Sample2[grep ('RPL',archp$Sample2)] = 'normal_pleura'

  ## and replace with Mo's annotation ##
  # normal_mo = read.csv ('../../tumor_compartment/all_tissues_ArchR/metadata_external_data/tsankov_refined_annotation.csv')



  
  #archp = archp[!is.na(archp$celltype)]

  # Dimensionality reduction and clustering
  varfeat = 25000
  LSI_method = 2
  archp = addIterativeLSI (ArchRProj = archp,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = FALSE, LSIMethod = LSI_method,
    varFeatures = varfeat)

#   archp = addHarmony (
#     ArchRProj = archp,
#     reducedDims = "IterativeLSI",
#     name = "Harmony",
#     groupBy = "Sample", force=FALSE
# )

# archp = addUMAP (ArchRProj = archp, 
#     reducedDims = "Harmony", name='UMAP_H',
#     force = TRUE)

# archp = addClusters (input = archp,
#     reducedDims = "Harmony",
#     name='Clusters_H',
#     force = TRUE)

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
   name = "Sample", embedding = "UMAP")
  umap_p2 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "celltype",
     embedding = "UMAP")
  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "nFrags",
     embedding = "UMAP")
  umap_p4 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "TSSEnrichment",
     embedding = "UMAP")
  umap_p5 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters",
     embedding = "UMAP")
  
  pdf (file.path('Plots','celltype_umap_signac_filtered.pdf'),12,12)
  print (umap_p1)
  print (umap_p2)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()
  
  plotPDF (umap_p1, umap_p2, umap_p3, umap_p4,
   name = paste0('Plot-UMAP-Sample-Clusters_',LSI_method,'_',length(rownames(archp)),'_varfeat_',varfeat,'.pdf'),
          ArchRProj = archp, addDOC = FALSE, width = 5, height = 5,logFile=NULL)
  
  archp = saveArchRProject (archp, load = T)
  
  } else {
  archp = loadArchRProject (projdir)
  }

metaGroupName = 'Clusters'
  archp = addGroupCoverages (
    ArchRProj = archp, 
    groupBy = metaGroupName,  
    force = FALSE,
    minCells= 20, # I think this should be set corresponding to the smallest cluster in the group or lower
    maxCells = 500,
    minReplicates = 2,
    sampleRatio = 0.8,
    useLabels = TRUE)
  
  archp = addReproduciblePeakSet (
      archp,
      groupBy= metaGroupName,
      peakMethod = 'Macs2',
      reproducibility = "1",
      maxPeaks = 500000, 
      minCells=20,
      force =TRUE) # I think this should be set corresponding to the smallest cluster in the group or lower
  archp = addPeakMatrix (archp)
  
  archp = saveArchRProject (archp, load=TRUE)
  
  archp = addBgdPeaks (archp, force= TRUE)
  archp = addMotifAnnotations (ArchRProj = archp,
      motifSet = "cisbp",
      #motifSet = 'JASPAR2020',
      #name = "JASPAR2020_Motif",
      force=TRUE)
  archp = addDeviationsMatrix (
    ArchRProj = archp, 
    peakAnnotation = "Motif",
    force = TRUE
  )
  
  archp = saveArchRProject (ArchRProj = archp,  
      load = TRUE)


# Find activating and repressing TFs ####
run_activeTF = FALSE

devMethod = 'ArchR'
 if (devMethod == 'ArchR')
    {
    TF_db='Motif'
    mSE = ArchR::getMatrixFromProject (archp, useMatrix = paste0(TF_db,'Matrix'))
    mSE = mSE[, archp$cellNames]
    rowData(mSE)$name = gsub ('_.*','',rowData(mSE)$name)
    rowData(mSE)$name = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", rowData(mSE)$name)
    }
  
if (if (!file.exists ('TF_activators_genescore.rds')))
  {
  seGroupMotif <- getGroupSE(ArchRProj = archp, useMatrix = "MotifMatrix", groupBy = "Clusters")
  seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
  rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
    rowMaxs(assay(seZ) - assay(seZ)[,x])
  }) %>% Reduce("cbind", .) %>% rowMaxs
  corGSM_MM <- correlateMatrices(
      ArchRProj = archp,
      useMatrix1 = "GeneScoreMatrix",
      useMatrix2 = "MotifMatrix",
      reducedDims = "IterativeLSI"
  )
  corGSM_MM = corGSM_MM[!grepl ('-AS',corGSM_MM$GeneScoreMatrix_name),]
  corGSM_MM = corGSM_MM[!grepl ('-DT',corGSM_MM$GeneScoreMatrix_name),]
  corGSM_MM = corGSM_MM[!grepl ('-OT',corGSM_MM$GeneScoreMatrix_name),]
  corGSM_MM = corGSM_MM[!grepl ('-RAB5IF',corGSM_MM$GeneScoreMatrix_name),]
  corGSM_MM = corGSM_MM[!grepl ('-IT2',corGSM_MM$GeneScoreMatrix_name),]
  corGSM_MM = corGSM_MM[!grepl ('-C8orf76',corGSM_MM$GeneScoreMatrix_name),]
  corGSM_MM = na.omit (corGSM_MM)
  saveRDS (corGSM_MM, 'TF_activators_genescore.rds')
  } else {
  corGSM_MM = readRDS ('TF_activators_genescore.rds') 
  }


### Co-expression of TFs #### 
metaGroupName = 'Sample2'
if (!any (ls() == 'mSE')) mSE = ArchR::getMatrixFromProject (archp, useMatrix = 'MotifMatrix', logFile=NULL)
mSE = mSE[, archp$cellNames]
all (colnames(mSE) == rownames(archp))

# # Get deviation matrix ####
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Subset only for positively correlated TF with genescore ####
positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0.1]
mMat = mMat[positive_TF,]

# mMat_agg = as.data.frame (t(mMat))
# mMat_agg$metaGroup = as.character (archp_meta[,metaGroupName])
# mMat_agg = aggregate (.~ metaGroup, mMat_agg, mean)
# rownames (mMat_agg) = mMat_agg[,1]
# mMat_agg = mMat_agg[,-1]
# mMat_agg = t(mMat_agg)
# rownames (mMat_agg) = active_TF

mMat_cor = cor (as.matrix(t(mMat)), method = 'pearson')
#d = as.dist (1-cor(as.matrix(t(mMat))))

#d = dist (mMat, method ='euclidean')
#hc1 <- hclust(d, method = "complete" ) # Hierarchical clustering using Complete Linkage

# row_filt = rowSums (mMat_cor) != 0
# tf_name = rownames(mMat_cor)[row_filt]
km = kmeans (mMat_cor, centers=20)
#km_ha = rowAnnotation (
#  km = anno_simple(as.character(km$cluster), width = unit(2, "mm"),
#    col = setNames (km_col, as.character(unique(km$cluster))), border=T))

# cor_mMat_hm = Heatmap (mMat_cor,
#         cluster_rows = T,
#         cluster_columns= T,
#         #left_annotation = km_ha,
#         row_split = km$cluster,
#         column_title = 'all_MPM',
#         column_gap = unit(.2, "mm"),
#         row_gap = unit(.2, "mm"),
#         name = 'r',
#         #column_km_repeats=20,
#         col = viridis::turbo(100),
#         #clustering_distance_rows = 'pearson', 
#         #clustering_distance_columns = 'pearson', 
#         column_split = km$cluster,
#         #col = pal_corr1,
#         row_names_gp = gpar(fontsize = 0),
#         column_names_gp = gpar(fontsize = 0)#,
#         #width = unit(5, "cm")
#         )
cor_mMat_hm = draw (Heatmap (mMat_cor,# row_km=15,
  #left_annotation = ha,
  #rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=palette_module_correlation_fun, 
  row_split = km$cluster,
  column_split = km$cluster,
  #row_km=2, 
  #column_km=2,
#  right_annotation = ha,
  border=F,
  row_names_gp = gpar(fontsize = 0),
  column_names_gp = gpar(fontsize = 0)))#,
  # ,
  # cell_fun = function(j, i, x, y, w, h, fill) {
  #       if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
  #           grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#        }}))

pdf (file.path ('Plots','TF_modules.pdf'), width = 4,height=3)
cor_mMat_hm
dev.off()

tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = do.call (cbind, tf_modules)
archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% paste0('mod_',unique(km$cluster))]
archp@cellColData = cbind (archp@cellColData, tf_modules)

TF_p = lapply (paste0('mod_',unique(km$cluster)), function(x) plotEmbedding (
    ArchRProj = archp,
    colorBy = "cellColData",
    name = x, 
    pal = palette_deviation,
    #useSeqnames='z',
    embedding = "UMAP"))

pdf (file.path ('Plots','TF_modules_umap.pdf'), width = 40,height=20)
wrap_plots (TF_p, ncol=5)
dev.off()




