set.seed (1234)
addArchRThreads (threads = 8) 
addArchRGenome ("Hg38")

projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR'
setwd (projdir)
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
    'P14',#,# p14
    'P23'
    )


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
    '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_14/240109/ZHAO13mesotheliomaATAC/cellranger_output/ALTS04_Zhao13ATAC_0_v1/fragments.tsv.gz',
    '/sc/arion/projects/Tsankov_Normal_Lung/data/meso_polyICLC/23Mesothelioma/cellranger_output/ALTS04_P22_0_v1/fragments.tsv.gz'#,#,
    #'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz',
    #'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz',
    #'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz',
    #'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz'#,
    #"/ahg/regevdata/projects/ICA_Lung/10x/200116/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01396_AlexTsankov/aDcd45n/outs/fragments.tsv.gz"
    )
     
    #fragment_paths = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/ArrowFiles'
    #fragment_paths = file.path(fragment_paths,sample_names)
    ArrowFiles_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/ArrowFiles'
    ArrowFiles = file.path(ArrowFiles_dir, paste0(sample_names,'.arrow'))
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
  
  # Import revised cell annotation after multiple iterations of Umaps (NOT SAVED)
  ann = read.csv ('../../git_repo/files/barcode_annotation.csv')
  archp$celltype_revised = ann$celltype[match(rownames(archp@cellColData), ann$barcode)]
  archp = archp[!is.na(archp$celltype_revised)]

  ### QC plots ####
  pdf()
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
dev.off()
  pdf (file.path ('Plots', 'QC_plots.pdf'))
  wrap_plots (p1, p2 ,p3, p4)
  dev.off()


  varfeat = 25000
  LSI_method = 2
  archp = addIterativeLSI (ArchRProj = archp,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = TRUE, LSIMethod = LSI_method,
    varFeatures = varfeat)

  archp = addClusters (input = archp, resolution = 3,
    reducedDims = "IterativeLSI", maxClusters = 100,
    force = TRUE)
  archp = addUMAP (ArchRProj = archp, 
    reducedDims = "IterativeLSI",
    force = TRUE)

pdf ()  
umap_p0 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters",
   embedding = "UMAP",
   #pal = palette_celltype_lv1,
   labelMeans = FALSE)

umap_p1 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype_revised",
   embedding = "UMAP",
   pal = palette_celltype_lv1,
   labelMeans = FALSE)

umap_p2 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
   embedding = "UMAP",
   pal = palette_sample,
   labelMeans = FALSE)
dev.off()
pdf (file.path ('Plots','celltype_revised_umap.pdf'))
umap_p0
umap_p1
umap_p2
dev.off()


# Plot gene score of cell type markers ####
metaGroupName = 'celltype_revised'
  DAG_list = readRDS (paste0('DAG_',metaGroupName,'.rds'))
  
  FDR_threshold = 1e-2
  lfc_threshold = 0
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
  
pdf()
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "GeneScoreMatrix", 
    name = head(DAG_df$gene[DAG_df$comparison == 'Mesothelium'],20), 
    embedding = "UMAP",
    pal = palette_expression,
    #imputeWeights = getImputeWeights(archp)
    imputeWeights = NULL
)
dev.off()

p = lapply (seq_along(p), function(x) p[[x]] + theme_void()  + theme(
    axis.title = element_blank(),       # Remove axis titles
    axis.text = element_blank(),        # Remove axis tick labels       # Remove plot title
    plot.subtitle = element_blank(),    # Remove subtitle
    plot.caption = element_blank(),     # Remove caption
    legend.position = "none"            # Remove legend
  ) + ggtitle (head(DAG_df$gene[DAG_df$comparison == 'Mesothelium'],20)[x])
)
pdf (file.path('Plots','marker_genes_feature_plots2.pdf'), width = 10, height = 10)
print (wrap_plots (p, ncol = 4))
dev.off()

meso_markers = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/gene_sets/highlevel_MPM_markers.csv')[[1]]
meso_markers = meso_markers[meso_markers != 'IGHM']
meso_markers = c(meso_markers, 'KRT5','LILRA4','MS4A1')
meso_markers = c(
  'KRT19','AMOTL2','EGFR',
  'HP','BDKRB1','ZBTB7C',
  'SFTA3','SFTPB','LINC00261',
  'COL1A1','FBN1','COL6A2',
  'MYH11','COL4A2','COL4A1',
  'CLDN5','VWF','CDH5',
  'LYZ','IL1B','CD83',
  'CD3E','BCL11B','RUNX3',
  'GNLY', 'KLRD1','ITGAE',
  'CD79A','PAX5','BLK',
  'IGLL5','ELL2','ELL2',
  'VASH2','MAD1L1','ZFAT')
#meso_markers = rev (meso_markers)
archp = addImputeWeights (archp)
pdf()
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "GeneScoreMatrix", 
    name = meso_markers, 
    embedding = "UMAP",
    pal = palette_expression,
    #imputeWeights = getImputeWeights(archp)
    imputeWeights = NULL
)
dev.off()

p = lapply (seq_along(p), function(x) p[[x]] + theme_void()  + theme(
    axis.title = element_blank(),       # Remove axis titles
    axis.text = element_blank(),        # Remove axis tick labels       # Remove plot title
    plot.subtitle = element_blank(),    # Remove subtitle
    plot.caption = element_blank(),     # Remove caption
    legend.position = "none"            # Remove legend
  ) + ggtitle (meso_markers[x])
)
pdf (file.path('Plots','marker_genes_feature_plots.pdf'), width = 10, height = 10)
print (wrap_plots (p, ncol = 4))
dev.off()


  if (!any (ls() == 'gsSE')) gsSE = ArchR::getMatrixFromProject (archp, useMatrix = 'GeneScoreMatrix')
  gsSE = gsSE[, archp$cellNames]
  gsMat = assays (gsSE)[[1]]
  celltype_order = c('Malignant','Mesothelium','Alveolar','Fibroblast',
    'SmoothMuscle','Endothelial','Myeloid','T_cells','NK','B_cells','Plasma','pDCs')
  rownames (gsMat) = rowData (gsSE)$name
  gsMat_mg = gsMat[meso_markers, ]
  gsMat_mg = as.data.frame (t(gsMat_mg))
  gsMat_mg$metaGroup = as.character(archp@cellColData[,metaGroupName])
  gsMat_mg = aggregate (.~ metaGroup, gsMat_mg, mean)
  rownames (gsMat_mg) = gsMat_mg[,1]
  gsMat_mg = gsMat_mg[,-1]
  gsMat_mg = gsMat_mg[names(table (archp@cellColData[,metaGroupName])[table (archp@cellColData[,metaGroupName]) > 50]),]
  gsMat_mg = gsMat_mg[celltype_order,]
  DAG_hm = Heatmap (t(t(scale(gsMat_mg))), 
          column_labels = colnames (gsMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          row_names_side = 'left',
          #col = palette_genescore_fun(gsMat_mg),
          col = palette_expression,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          column_names_gp = gpar(fontsize = 6),
          rect_gp = gpar(col = "white", lwd = .1),
          border=TRUE
          #right_annotation = motif_ha
          )

  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (file.path('Plots',paste0('DAG_clusters_',metaGroupName,'_heatmaps.pdf')), width = 6, height = 3)
print(DAG_hm)
dev.off()





# Run genescore DAG ####
metaGroupName = "celltype_lv1"
force = FALSE
if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source (file.path('..','..','git_repo','utils','DAG.R'))

pdf()
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "GeneScoreMatrix", 
    name = DAG_df$gene[DAG_df$comparison %in% c('Alveolar','pDCs','Plasma')], 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
pdf (file.path('Plots','marker_genes_pDC_plasma_Alveolar.pdf'), width = 40, height = 40)
print (wrap_plots (p, ncol = 8))
dev.off()

archp = saveArchRProject (archp, load = T)

# Export cell annotation ####
write.csv (data.frame (barcode = rownames(archp@cellColData), celltype = archp$celltype_revised), '../../git_repo/files/barcode_annotation.csv')

