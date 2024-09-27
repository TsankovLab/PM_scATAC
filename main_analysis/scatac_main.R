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
dir.create (paste0(projdir,'/Plots/'), recursive =T)
setwd (projdir)


#devtools::install_github("immunogenomics/presto") #needed for DAA
source ('../../git_repo/utils/useful_functions.R')
source ('../../git_repo/utils/ggplot_aestetics.R')
source ('../../git_repo/utils/scATAC_functions.R')
source ('../../git_repo/utils/palettes.R')

set.seed(1234)
addArchRThreads (threads = 8) 
addArchRGenome("Hg38")

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
  
  # Dimensionality reduction and clustering
  varfeat = 25000
  LSI_method = 2
  archp = addIterativeLSI (ArchRProj = archp,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = TRUE, LSIMethod = LSI_method,
    varFeatures = varfeat)

  archp = addHarmony (
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample", force=FALSE
)

archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony",
    name='Clusters_H',
    force = TRUE)

  archp = addClusters (input = archp, resolution = .8,
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
    colorBy = "cellColData", name = "nFrags",
     embedding = "UMAP")
  umap_p4 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "TSSEnrichment",
     embedding = "UMAP")
  
  
  pdf (file.path('Plots','celltype_umap_signac_filtered2.pdf'))
  print (umap_p1)
  print (umap_p2)
  print (umap_p3)
  print (umap_p4)
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


### Gene score based analysis ####
run_GS_analysis = FALSE

if (run_GS_analysis)
  {
  # Find DAG ####
  metaGroupName = "Clusters"
  force = FALSE
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
          cluster_rows = T,
          #col = pals_heatmap[[5]],
          cluster_columns=T,#col = pals_heatmap[[1]],
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
  
# Coverage plots of meso developmental TFs
meso_markers = c('WT1','HP','KRT5','COL1A1','ACTA2','VWF','CD3D','GNLY','LYZ','LILRA4','MS4A1')
archp = addImputeWeights (archp)
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "GeneScoreMatrix", 
    name = meso_markers, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(archp)
)

pdf (file.path('Plots',paste0('marker_genes_',metaGroupName,'_feature_plots.pdf')), width = 15, height = 15)
print (wrap_plots (p, ncol = 4))
dev.off()

  ### Plot cell type markers using GeneScores ####
  metaGroupName = 'Sample2'
  celltype_markers = c('WT1','CALB2','GATA4','HP','SOX9','MESP1','SOX6','TWIST1','SNAI2')
  #celltype_markers = c('WT1','CALB2','GATA4','MSLN','KRT5','KRT18','ITLN1','HP','SOX9')
  meso_markers <- plotBrowserTrack(
      ArchRProj = archp, 
      groupBy = metaGroupName, 
      geneSymbol = celltype_markers,
      #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
      #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
      upstream = 250000,
      downstream = 250000,
      loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
      #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
      #loops = getCoAccessibility (archp, corCutOff = 0.3,
      #  returnLoops = TRUE),
      useGroups= NULL
  )
  plotPDF (meso_markers, ArchRProj = archp, width=14, name ='MPM_markers_coveragePlots.pdf')
}


### Run peak calling on celltype annotation ####
run_peakCall = FALSE

if (run_peakCall)
  {
  ### Call peaks on celltypes ###
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
      force =FALSE) # I think this should be set corresponding to the smallest cluster in the group or lower
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




