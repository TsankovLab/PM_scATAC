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
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR'
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
sarc_order = read.csv ('../scrna/cnmf20_sarcomatoid_sample_order.csv', row.names=1)

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
    #'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz',
    #'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz',
    #'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz',
    #'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz'#,
    #"/ahg/regevdata/projects/ICA_Lung/10x/200116/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01396_AlexTsankov/aDcd45n/outs/fragments.tsv.gz"
    )
     
    
    ArrowFiles_dir = '/ahg/regevdata/projects/ICA_Lung/Bruno/ArrowFiles/'
    if (!all (paste0(sample_names,'.arrow') %in% list.files(ArrowFiles_dir)))
      {
      #setwd (projdir)  
      ArrowFiles = createArrowFiles (inputFiles = fragment_paths,
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
    ArrowFiles = ArrowFiles, 
    outputDirectory = projdir,
    copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  )
  
  ### Subset ArchR object only for cells retained in Signac analysis ####
  cell_annotation = read.csv ('../../per_sample_QC_signac/cell_annotation.csv')
  colnames (cell_annotation) = c('barcode','celltype')
  cell_annotation$barcode = sub ('\\.','#', cell_annotation$barcode)
  
  archp$celltype = cell_annotation$celltype[match(rownames(archp@cellColData), cell_annotation$barcode)]
  archp$celltype[archp$celltype == 'Myeloid'] = 'MonoMac'
  archp$celltype[archp$celltype == 'pDC'] = 'pDCs'
  archp = archp[!is.na(archp$celltype)]
  archp = archp[archp$celltype != 'bad_quality']
  
  #archp = archp[rownames(archp) %in% keep_barcodes]
  
  ### QC plots ####
  p1 = plotFragmentSizes(ArchRProj = archp, groupBy = 'Sample', pal = palette_sample)
  p2 = plotTSSEnrichment(ArchRProj = archp, groupBy = 'Sample', pal = palette_sample)

  p3 <- plotGroups(
    ArchRProj = archp, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    pal = palette_sample,
    alpha = 0.4,
    addBoxPlot = TRUE
   )

  p4 <- plotGroups(
    ArchRProj = archp, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "nFrags",
    plotAs = "violin",
    pal = palette_sample,
    alpha = 0.4,
    addBoxPlot = TRUE
   )

  pdf (file.path ('Plots', 'QC_plots.pdf'))
  wrap_plots (p1, p2 ,p3, p4)
  dev.off()


  # Dimensionality reduction and clustering
  varfeat = 25000
  LSI_method = 2
  archp = addIterativeLSI (ArchRProj = archp,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = TRUE, LSIMethod = LSI_method,
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
  archp = addTSNE (ArchRProj = archp, 
    reducedDims = "IterativeLSI",
    force = TRUE)
  
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
  
  tsne_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
   name = "Sample", embedding = "TSNE")
  tsne_p2 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "celltype",
     embedding = "TSNE")
  tsne_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "nFrags",
     embedding = "TSNE")
  tsne_p4 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "TSSEnrichment",
     embedding = "TSNE")
  
  pdf (file.path('Plots','celltype_tsne_signac_filtered2.pdf'))
  print (tsne_p1)
  print (tsne_p2)
  print (tsne_p3)
  print (tsne_p4)
  dev.off()
  
  archp = saveArchRProject (archp, load = T)
  
  } else {
  archp = loadArchRProject (projdir)
  }

sarc_order = c('P1','P13','P3','P12','P5','P11','P4','P8','P14','P10')
archp$Sample2 = archp$Sample
archp$Sample2 = factor (archp$Sample2, levels = sarc_order)

# Plot gene score of cell type markers ####
meso_markers = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/highlevel_MPM_markers.csv')[[1]]
meso_markers = meso_markers[meso_markers != 'IGHM']
meso_markers = c(meso_markers, 'KRT5','LILRA4','MS4A1')
archp = addImputeWeights (archp)
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "GeneScoreMatrix", 
    name = meso_markers, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(archp)
)

pdf (file.path('Plots','marker_genes_feature_plots.pdf'), width = 25, height = 25)
print (wrap_plots (p, ncol = 8))
dev.off()


### Further remove low quality / doublets clusters ####
#low_quality_clusters = c('C63','C62','C50','C61','C56','C2')
low_quality_clusters = c('C36','C43','C53','C54','C49','C55','C52')
table (archp$celltype[!archp$Clusters %in% low_quality_clusters], archp$Sample[!archp$Clusters %in% low_quality_clusters])

archp = archp[!archp$Clusters %in% low_quality_clusters]
  
  #archp = archp[rownames(archp) %in% keep_barcodes]
  
  # Dimensionality reduction and clustering
  varfeat = 25000
  LSI_method = 2
  archp = addIterativeLSI (ArchRProj = archp,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = TRUE, LSIMethod = LSI_method,
    varFeatures = varfeat)

  archp = addClusters (input = archp, resolution = 3,
    reducedDims = "IterativeLSI", maxClusters = 100,
    force = TRUE)
  archp = addClusters (input = archp, resolution = 20,
    name = 'Clusters_20',
    reducedDims = "IterativeLSI", maxClusters = 100,
    force = TRUE)
  archp = addUMAP (ArchRProj = archp, 
    reducedDims = "IterativeLSI",
    force = TRUE)
  archp = addTSNE (ArchRProj = archp, 
    reducedDims = "IterativeLSI",
    force = TRUE)
  
  archp = saveArchRProject (archp)

  umap_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
   name = "Sample", embedding = "UMAP")
  umap_p2 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "celltype",
     embedding = "UMAP")
  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters",
     embedding = "UMAP")
  umap_p4 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters_20",
     embedding = "UMAP")
  
  pdf (file.path('Plots','celltype_umap_signac_filtered_2.pdf'),22,22)
  print (umap_p1)
  print (umap_p2)
  print (umap_p3)
  print (umap_p4)
  dev.off()
  
### Gene score based analysis ####
run_GS_analysis = TRUE

if (run_GS_analysis)
  {
  # Find DAG ####
  metaGroupName = "Clusters"
  force = TRUE
  if (!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force)
    {
    DAG_list = getMarkerFeatures (
      ArchRProj = archp, 
      testMethod = "wilcoxon",
            #useGroups = "ClusterA",
            #bgdGroups = "Clusters1B",
      binarize = FALSE,
      useMatrix = "GeneScoreMatrix",
      groupBy = metaGroupName
    #  useSeqnames="z"
    )

    listnames = colnames (DAG_list)
    DAG_list = lapply (1:ncol (DAG_list), function(x) 
      {
      df = DAG_list[,x]  
      df = do.call (cbind, (assays(df)))
      colnames(df) = names (assays(DAG_list))
      df$gene = rowData (DAG_list)$name
      df
      })
    names (DAG_list) = listnames
    saveRDS (DAG_list, paste0 ('DAG_',metaGroupName,'.rds'))    
    } else {
    DAG_list = readRDS (paste0('DAG_',metaGroupName,'.rds'))
    }
  
  FDR_threshold = 1e-8
  lfc_threshold = 1
  top_genes = 20
  DAG_top_list = DAG_list[sapply (DAG_list, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$Log2FC) > lfc_threshold,]) > 0)]
  DAG_top_list = lapply (seq_along(DAG_top_list), function(x) {
    res = DAG_top_list[[x]]
    res = na.omit (res)
    res = res[res$FDR < FDR_threshold,]
    res = res[order (res$FDR), ]
    res = res[abs(res$Log2FC) > lfc_threshold,]
    res$comparison = names(DAG_top_list)[x]
    if (nrow(res) < top_genes) 
      {
      res
      } else {
      head (res,top_genes)
      }
    })
  DAG_df = Reduce (rbind ,DAG_top_list)
  
  if (!any (ls() == 'gsSE')) gsSE = ArchR::getMatrixFromProject (archp, useMatrix = 'GeneScoreMatrix')
  gsSE = gsSE[, archp$cellNames]
  gsMat = assays (gsSE)[[1]]
  rownames (gsMat) = rowData (gsSE)$name
  gsMat_mg = gsMat[rownames (gsMat) %in% DAG_df$gene, ]
  gsMat_mg = as.data.frame (t(gsMat_mg))
  gsMat_mg$metaGroup = as.character(archp@cellColData[,metaGroupName])
  gsMat_mg = aggregate (.~ metaGroup, gsMat_mg, mean)
  rownames (gsMat_mg) = gsMat_mg[,1]
  gsMat_mg = gsMat_mg[,-1]
  gsMat_mg = gsMat_mg[names(table (archp@cellColData[,metaGroupName])[table (archp@cellColData[,metaGroupName]) > 50]),]
  DAG_hm = Heatmap (t(scale(gsMat_mg)), 
          row_labels = colnames (gsMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 6),
          column_names_gp = gpar(fontsize = 4),
          rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE
          #right_annotation = motif_ha
          )

  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (paste0('Plots/DAG_clusters_',metaGroupName,'_heatmaps.pdf'), width = 8, height = 50)
print(DAG_hm)
dev.off()



# Re-annotate clusters ####
archp$celltype_revised = archp$celltype
archp$celltype_revised[archp$Clusters %in% c('C11')] = 'Mesothelium'
archp$celltype_revised[archp$Clusters %in% c('C20')] = 'SmoothMuscle'
archp$celltype_revised[archp$Clusters %in% c('C21','C23','C24')] = 'Fibroblasts'
archp$celltype_revised[archp$Clusters %in% c('C28')] = 'Alveolar'
archp$celltype_revised[archp$Clusters %in% c('C15','C16','C7','C8','C9','C25','C12','C13','C14','C4','C5','C6','C1','C2','C3','C27')] = 'Malignant'
archp$celltype_revised[archp$Clusters %in% c('C17','C18','C19')] = 'Endothelial'
archp$celltype_revised[archp$Clusters %in% c('C26')] = 'bad_quality' # Further remove suspected low quality cluster C26 ####
archp$celltype_revised[archp$Clusters_20 %in% c('C28','C29','C32','C33','C31','C26','C24','C23','C25','C30','C27','C41','C44','C42','C43','C45','C39','C38','C22','C37')] = 'T_cells'
archp$celltype_revised[archp$Clusters_20 %in% c('C34','C35','C36','C40')] = 'NK'
archp$celltype_revised[archp$Clusters %in% c('C29','C30','C33','C32','C36','C34','C35','C31')] = 'Myeloid'
archp$celltype_revised[archp$Clusters %in% c('C53')] = 'pDCs'
archp$celltype_revised[archp$Clusters %in% c('C54')] = 'Plasma'
archp$celltype_revised[archp$Clusters %in% c('C58','C57','C56','C55')] = 'B_cells'



# Further remove low quality clusters ####
archp = archp[archp$celltype_revised != 'bad_quality']

umap_p1 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype_revised",
   embedding = "UMAP",
   pal = palette_celltype_simplified,
   labelMeans = FALSE)

umap_p2 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample2",
   embedding = "UMAP",
   pal = palette_sample,
   labelMeans = FALSE)

pdf (file.path ('Plots','celltype_revised_umap.pdf'))
umap_p1
umap_p2
dev.off()


# Plot gene score of cell type markers ####
meso_markers = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/highlevel_MPM_markers.csv')[[1]]
meso_markers = meso_markers[meso_markers != 'IGHM']
meso_markers = c(meso_markers, 'KRT5','LILRA4','MS4A1')
archp = addImputeWeights (archp)
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = meso_markers, 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)

#p = lapply (p, function(x) x + theme_void() + NoLegend ()) #+ ggtitle scale_fill_gradient2 (rev (viridis::plasma(100))))

pdf (file.path('Plots','marker_genes_feature_plots_3.pdf'), width = 25, height = 25)
print (wrap_plots (p, ncol = 8))
dev.off()

archp = saveArchRProject (archp)

# Export cell annotation ####
write.csv (data.frame (barcode = rownames(archp@cellColData), celltype = archp$celltype_revised), 'barcode_annotation.csv')

### Run peak calling on celltype annotation ####
# Add tumor sample info in celltype metagroup
archp$celltype_revised_sample = archp$celltype_revised
archp$celltype_revised_sample[archp$celltype_revised_sample == 'Malignant'] = paste0('Malignant_',archp$Sample[archp$celltype_revised_sample == 'Malignant'])

run_peakCall = FALSE
if (run_peakCall)
  {
  ### Call peaks on celltypes ###
  metaGroupName = 'celltype_revised_sample'
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
  
  metaGroupNames = c('TSSEnrichment','nFrags','ReadsInTSS','FRIP')  
    umap_p12 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
     name = x, embedding = "UMAP"))
      
  pdf (paste0(projdir,'/Plots/qc_umap_after_filtering.pdf'), 15,15)
  wrap_plots (umap_p12, ncol=5)
  dev.off()
  }

### chromVAR analysis ####
run_chromVAR = TRUE

if (run_chromVAR)
  {  
  archp = addBgdPeaks (archp, force= FALSE)
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
  }



# Find activating and repressing TFs ####
run_activeTF = FALSE

devMethod = 'ArchR'
 if (devMethod == 'ArchR')
    {
    TF_db='Motif'
    if (!exists('mSE')) mSE = ArchR::getMatrixFromProject (archp, useMatrix = paste0(TF_db,'Matrix'))
    mSE = mSE[, archp$cellNames]
    rowData(mSE)$name = gsub ('_.*','',rowData(mSE)$name)
    rowData(mSE)$name = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", rowData(mSE)$name)
    }
  
if (!file.exists ('TF_activators_genescore.rds'))
  {
  seGroupMotif <- getGroupSE(ArchRProj = archp, useMatrix = "MotifMatrix", groupBy = "Clusters")
  seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
  rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
    rowMaxs(assay(seZ) - assay(seZ)[,x])
  }) %>% Reduce("cbind", .) %>% rowMaxs
  corGSM_MM <- correlateMatrices (
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


### ChromVAR based analysis ####
run_chromVAR_analysis = TRUE

if (run_chromVAR_analysis)
  {
  # Find DAM ####
  metaGroupName = "celltype_revised"
  force = FALSE
  if (!file.exists (paste0('DAM_',metaGroupName,'.rds')) | force)
    {
    DAM_list = getMarkerFeatures (
      ArchRProj = archp, 
      testMethod = "wilcoxon",
            #useGroups = "ClusterA",
            #bgdGroups = "Clusters1B",
      binarize = FALSE,
      useMatrix = "MotifMatrix",
      groupBy = metaGroupName
    #  useSeqnames="z"
    )

    listnames = colnames (DAM_list)
    DAM_list = lapply (1:ncol (DAM_list), function(x) 
      {
      df = DAM_list[,x]  
      df = do.call (cbind, (assays(df)))
      colnames(df) = names (assays(DAM_list))
      df$gene = rowData (DAM_list)$name
      df
      })
    names (DAM_list) = listnames
    saveRDS (DAM_list, paste0 ('DAM_',metaGroupName,'.rds'))    
    } else {
    DAM_list = readRDS (paste0('DAM_',metaGroupName,'.rds'))
    }
  
  active_genes = corGSM_MM$MotifMatrix_name[corGSM_MM$cor > 0.1]
  DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% active_genes,])
  names (DAM_list2) = names (DAM_list)
  FDR_threshold = 1e-3
  meandiff_threshold = 0
  top_genes = 3
  DAM_top_list = DAM_list2[sapply (DAM_list2, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$MeanDiff) > meandiff_threshold,]) > 0)]
  DAM_top_list = lapply (seq_along(DAM_top_list), function(x) {
    res = DAM_top_list[[x]]
    #res = na.omit (res)
    res = res[res$FDR < FDR_threshold,]
    res = res[order (res$FDR), ]
    res = res[res$MeanDiff > meandiff_threshold,]
    res$comparison = names(DAM_top_list)[x]
    if (nrow(res) < top_genes) 
      {
      res
      } else {
      head (res,top_genes)
      }
    })
  DAM_df = Reduce (rbind ,DAM_top_list)
  
  devMethod = 'ArchR'
 if (devMethod == 'ArchR')
    {
    TF_db='Motif'
    if (!exists ('mSE')) mSE = ArchR::getMatrixFromProject (archp, useMatrix = paste0(TF_db,'Matrix'))
    mSE = mSE[, archp$cellNames]
    rowData(mSE)$name = gsub ('_.*','',rowData(mSE)$name)
    rowData(mSE)$name = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", rowData(mSE)$name)
    }
  DAM_df$gene = gsub ('_.*','',DAM_df$gene)
  DAM_df$gene = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", DAM_df$gene)
  mMat = assays (mSE)[[1]]
  rownames (mMat) = rowData (mSE)$name
  mMat_mg = mMat[DAM_df$gene, ]
  mMat_mg = as.data.frame (t(mMat_mg))
  mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
  mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
  rownames (mMat_mg) = mMat_mg[,1]
  mMat_mg = mMat_mg[,-1]
  mMat_mg = mMat_mg[names (DAM_list),]
  #mMat_mg = mMat_mg[names(table (archp@cellColData[,metaGroupName])[table (archp@cellColData[,metaGroupName]) > 50]),]
  DAM_hm = Heatmap (t(scale(mMat_mg)), 
          row_labels = colnames (mMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          name = 'chromVAR',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = rev(palette_deviation)

          #right_annotation = motif_ha
          )

  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_heatmaps.pdf')), width = 3, height = 5)
print(DAM_hm)
dev.off()
}

# Find shared TF across celltypes ####
tf_name2 = unlist(sapply (c('SOX9','TWIST1','MESP1','NKX2-5'), function(x) rownames(assay(mSE))[grepl (x, rownames(assay(mSE)))]))
tf_name2 = paste0('z:',tf_name2)
archp = addImputeWeights (archp)
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "MotifMatrix",
    name = tf_name2, 
    useSeqnames='z',
    col = palette_deviation,    
    embedding = "UMAP",
    imputeWeights = getImputeWeights(archp)
    )

pdf (file.path ('Plots','TF_umap.pdf'), width = 20,height=6)
wrap_plots (TF_p, ncol=5)
dev.off()


### Co-expression of TFs #### 
metaGroupName = 'celltype_revised'
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

archp = addImputeWeights (archp)
TF_p = lapply (paste0('mod_',unique(km$cluster)), function(x) plotEmbedding (
    ArchRProj = archp,
    colorBy = "cellColData",
    name = x, 
    pal = palette_deviation,
    #useSeqnames='z',
    imputeWeights = getImputeWeights(archp),
    embedding = "UMAP"))

pdf (file.path ('Plots','TF_modules_umap.pdf'), width = 30,height=14)
wrap_plots (TF_p, ncol=5)
dev.off()

#### make heatmap of all positively regulated TFs across cell types to find shared programs ####
# # Get deviation matrix ####
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Subset only for positively correlated TF with genescore ####
metaGroupName = 'celltype_revised'
positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0.1]
mMat = mMat[positive_TF,]
mMat = as.data.frame (t(mMat))
#mMat$metaGroup = as.character (archp@cellColData[,metaGroupName])
#mMat = aggregate (.~ metaGroup, mMat, mean)
#rownames(mMat) = mMat[,1]
#mMat = mMat[,-1]
ha = HeatmapAnnotation (celltype = as.character (archp@cellColData[,metaGroupName]), col=list(celltype = palette_celltype_simplified))
DAM_hm = Heatmap (t(mMat), 
          row_labels = colnames (mMat),
          top_annotation = ha,
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = T,
          #col = pals_heatmap[[5]],
          cluster_columns=T,#col = pals_heatmap[[1]],
          row_names_gp = gpar (fontsize = 0),
          column_names_gp = gpar(fontsize = 0),
          column_names_rot = 45,
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_deviation
          #right_annotation = motif_ha
          )

pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_all_TF_heatmaps.pdf')), width = 12.6, height = 5)
print(DAM_hm)
dev.off()





