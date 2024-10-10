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

  archp = saveArchRProject (archp, load = T)
  
  # Export cell annotation ####
  write.csv (data.frame (barcode = rownames(archp@cellColData), celltype = archp$celltype_revised), 'barcode_annotation.csv')


  } else {
  archp = loadArchRProject (projdir)
  }
