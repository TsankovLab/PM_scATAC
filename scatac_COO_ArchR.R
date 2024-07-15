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
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/COO/scatac_ArchR'
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

# Annotate WT1+ fibroblasts ####
meso_markers = 'WT1'
archp = addImputeWeights (archp)
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = meso_markers, 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)

archp$celltype[archp$Clusters == 'C30'] = 'Fibroblasts_WT1'
#p = lapply (p, function(x) x + theme_void() + NoLegend ()) #+ ggtitle scale_fill_gradient2 (rev (viridis::plasma(100))))

pdf (file.path('Plots','WT1_expression.pdf'), width = 25, height = 25)
print (p)
dev.off()


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


### Compare bins malignants normal ####
run_bin_analysis = TRUE

if (run_bin_analysis)
  {
  # Load fragments
  fragments = unlist (getFragmentsFromProject (
       ArchRProj = archp))
  
  ws = 1e7
  ss = 1e7
  if (!file.exists (paste0('bins_',ws,'_ss_',ss,'.rds')))
    {
    blacklist = toGRanges(paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/blacklisted_regions/ENCODE_blacklist/',"hg38-blacklist.v2.bed")) # taken from https://github.com/Boyle-Lab/Blacklist/tree/master/lists
    windows = makeWindows(genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist,
      windowSize = ws, slidingSize = ss)
    saveRDS (windows, paste0('bins_',ws,'_ss_',ss,'.rds'))
    } else {
    windows = readRDS (paste0('bins_',ws,'_ss_',ss,'.rds'))  
    }
  
  flt_celltypes = 100
  
  #archp$celltype_sample = paste0(archp$celltype, '-', archp$Sample2)

  metaGroupName = 'celltype'
  #fragments$RG = as.character(fragments$RG)  
  barcode_metaGroup = as.data.frame (archp@cellColData[,c(metaGroupName, 'TSSEnrichment','nFrags','Sample2')])
  colnames (barcode_metaGroup)[colnames (barcode_metaGroup) == metaGroupName] = 'metaGroup'
  barcode_metaGroup_tumor = barcode_metaGroup[barcode_metaGroup$metaGroup %in% c('P1','P10','P11','P12','P13','P14','P3','P4','P5','P8'),]
  barcode_metaGroup_normal = barcode_metaGroup[barcode_metaGroup$metaGroup %in% c('Endothelial','Fibroblasts','Fibroblasts_WT1','Pericytes','Mesothelium','SmoothMuscle') & barcode_metaGroup$Sample2 == 'normal_pleura',]
  barcode_metaGroup = rbind (barcode_metaGroup_tumor, barcode_metaGroup_normal)
  barcode_metaGroup$barcode = rownames(barcode_metaGroup)
  celltype_bins = lapply (unique(barcode_metaGroup$metaGroup), function(mg) 
    {
    high_quality_cells = barcode_metaGroup[barcode_metaGroup$metaGroup == mg,]
    high_quality_cells = high_quality_cells[order (-high_quality_cells$TSSEnrichment),] # ordering by TSSEnrichment seem to work the best
    high_quality_cells = head (high_quality_cells$barcode, flt_celltypes)
    fragments_group = fragments[fragments$RG %in% high_quality_cells]
    #names (fragments_cell) = NULL 
    #fragments_cts = sample (unique(fragments_ct$RG), 200)
    fr_ov = countOverlaps (windows, fragments_group)
  #  fr_ov = fr_ov / sum (archp$nFrags[archp$celltype_sample == x])
    #(fr_ov / sum (fr_ov)) * 1000
    })
  celltype_bins = do.call (cbind, celltype_bins)
  colnames (celltype_bins) = sapply (unique(barcode_metaGroup$metaGroup), function(x) unlist(strsplit (x,'-'))[1])
  head (celltype_bins)
  
  celltype_bins_cor = cor (celltype_bins, method = 'pearson')
  celltype_bins_cor = celltype_bins_cor[colnames(celltype_bins_cor) != 'P3',colnames(celltype_bins_cor) != 'P3']
  ha = HeatmapAnnotation (sample = sapply (colnames(celltype_bins_cor), function(x) unlist(strsplit (x,'-'))[2]), which='row')
  binH = Heatmap (celltype_bins_cor, col = viridis::rocket(100),# row_names_gp= gpar (fontsize=6), 
    #column_names_gp= gpar (fontsize=6), 
    right_annotation = ha,
    cluster_rows = T, #row_km = 2, column_km = 2,
    clustering_distance_rows = 'pearson', 
    clustering_distance_columns = 'pearson',
    column_title = paste('Binned',ws))
  
  png (file.path('Plots',paste0('binned_fragments_binned_',ws,'_celltype_heatmap.png')),width=600, height=500)
  binH
  dev.off()
  }

