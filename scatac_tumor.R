conda activate signac
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
  'ggrepel')
lapply(packages, require, character.only = TRUE)

####### ANALYSIS of TUMOR compartment #######
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR'
dir.create (paste0(projdir,'/Plots/'), recursive =T)
setwd (projdir)


#devtools::install_github("immunogenomics/presto") #needed for DAA
source ('../../PM_scATAC/useful_functions.R')
source ('../../PM_scATAC/ggplot_aestetics.R')
source ('../../PM_scATAC/scATAC_functions.R')

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
    'P14',# p14
    # Normal
    'RPL_280_neg_1',
    'RPL_280_neg_2',
    'RPL_Epi_1',
    'RPL_Epi_2'#,
    #'cf_distal'
    )

# Load RNA
srt = readRDS ('../scrna/scRNA_meso.rds')
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
    '/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_14/090124/ZHAO13mesotheliomaATAC/cellranger_output/ALTS04_Zhao13ATAC_0_v1/fragments.tsv.gz',#,
    '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz',
    '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz',
    '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz',
    '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz'#,
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
  
  archp$Sample2 = archp$Sample
  archp$Sample2[grep ('RPL',archp$Sample2)] = 'normal_pleura'
  

  ### Subset ArchR object only for cells retained in Signac analysis ####
  tumor_l = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/signac_list.rds')
  normal = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/signac_normal.rds')
  normal_annotation = data.frame (barcode= colnames(normal), celltype = normal$predicted.id)
  normal_annotation$barcode = gsub (paste0(sample_names[11],'_'),paste0(sample_names[11],'#'),normal_annotation$barcode)
  normal_annotation$barcode = gsub (paste0(sample_names[12],'_'),paste0(sample_names[12],'#'),normal_annotation$barcode)
  normal_annotation$barcode = gsub (paste0(sample_names[13],'_'),paste0(sample_names[13],'#'),normal_annotation$barcode)
  normal_annotation$barcode = gsub (paste0(sample_names[14],'_'),paste0(sample_names[14],'#'),normal_annotation$barcode)
  tumor_annotation = do.call (rbind, lapply (sample_names[1:10], function(x) data.frame (barcode= paste0(x,'#',colnames(tumor_l[[x]])), celltype = tumor_l[[x]]$celltype)))
  
  normal_annotation = normal_annotation[normal_annotation$celltype == 'Mesothelium',]
  tumor_annotation = tumor_annotation[tumor_annotation$celltype == 'Malignant',]
  tumor_annotation$sample = sapply (tumor_annotation$barcode, function(x) unlist(strsplit(x, '#'))[1])
  keep_barcodes = c(normal_annotation$barcode, tumor_annotation$barcode)
  
  archp = archp[rownames(archp) %in% keep_barcodes]
  
  # Dimensionality reduction and clustering
  varfeat = 25000
  LSI_method=2
  archp = addIterativeLSI (ArchRProj = archp, 
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force=TRUE, LSIMethod=LSI_method,
    varFeatures = varfeat)
  archp = addClusters (input = archp, resolution = .8,
    reducedDims = "IterativeLSI", maxClusters = 100,
    force = TRUE)
  archp = addUMAP (ArchRProj = archp, 
    reducedDims = "IterativeLSI",
    force = TRUE)
  archp = addTSNE (ArchRProj = archp, 
    reducedDims = "IterativeLSI",
    force = TRUE)
  
  archp = saveArchRProject (archp, dropCells = T)


  umap_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
   name = "Sample2", embedding = "UMAP")
  umap_p2 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters",
     embedding = "UMAP")
  umap_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "nFrags",
     embedding = "UMAP")
  umap_p4 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "TSSEnrichment",
     embedding = "UMAP")
  
  
  pdf ('Plots/celltype_umap.pdf')
  print (umap_p1)
  dev.off()
  
  plotPDF (umap_p1, umap_p2, umap_p3, umap_p4,
   name = paste0('Plot-UMAP-Sample-Clusters_',LSI_method,'_',length(rownames(archp)),'_varfeat_',varfeat,'.pdf'),
          ArchRProj = archp, addDOC = FALSE, width = 5, height = 5,logFile=NULL)
  
  tsne_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
   name = "Sample2", embedding = "TSNE")
  tsne_p2 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "Clusters",
     embedding = "TSNE")
  tsne_p3 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "nFrags",
     embedding = "TSNE")
  tsne_p4 = plotEmbedding (ArchRProj = archp, 
    colorBy = "cellColData", name = "TSSEnrichment",
     embedding = "TSNE")
  
   plotPDF (tsne_p1, tsne_p2, tsne_p3, tsne_p4,
   name = paste0('Plot-TSNE-Sample-Clusters_',LSI_method,'_',length(rownames(archp)),'_varfeat_',varfeat,'.pdf'),
          ArchRProj = archp, addDOC = FALSE, width = 5, height = 5,logFile=NULL)
  
  } else {
  #RITfilter = 1000
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
  DAG_hm
  dev.off()
  
  ### Plot cell type markers using GeneScores ####
  metaGroupName = 'Sample2'
  celltype_markers = c('WT1','CALB2','GATA4','MSLN','KRT5','KRT18','ITLN1','HP','SOX9')
  #celltype_markers = c('WT1','CALB2','GATA4','MSLN','KRT5','KRT18','ITLN1','HP','SOX9')
  meso_markers <- plotBrowserTrack(
      ArchRProj = archp, 
      groupBy = metaGroupName, 
      geneSymbol = celltype_markers,
      #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
      #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
      upstream = 150000,
      downstream = 150000,
      loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
      #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
      #loops = getCoAccessibility (archp, corCutOff = 0.3,
      #  returnLoops = TRUE),
      useGroups= NULL
  )
  plotPDF (meso_markers, ArchRProj = archp, width=14, name ='MPM_markers_coveragePlots.pdf')
  }


### Compare bins malignants normal ####
run_bin_analysis = FALSE

if (run_bin_analysis)
  {
  # Load fragments
  fragments = unlist (getFragmentsFromProject (
       ArchRProj = archp))
  
  ws = 1e6
  ss = 2e5
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
  metaGroupName = 'Clusters'
  
  #fragments$RG = as.character(fragments$RG)  
  barcode_metaGroup = as.data.frame (archp@cellColData[,c(metaGroupName, 'TSSEnrichment','nFrags')])
  colnames (barcode_metaGroup)[colnames (barcode_metaGroup) == metaGroupName] = 'metaGroup'
  barcode_metaGroup$barcode = rownames(barcode_metaGroup)
  celltype_bins = lapply (unique(archp@cellColData[,metaGroupName]), function(mg) 
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
  colnames (celltype_bins) = sapply (unique(archp@cellColData[,metaGroupName]), function(x) unlist(strsplit (x,'-'))[1])
  head (celltype_bins)
  
  celltype_bins_cor = cor (celltype_bins, method = 'spearman')
  ha = HeatmapAnnotation (sample = sapply (unique(archp@cellColData[,metaGroupName]), function(x) unlist(strsplit (x,'-'))[2]), which='row')
  binH = Heatmap (celltype_bins_cor, col = viridis::rocket(100),# row_names_gp= gpar (fontsize=6), 
    #column_names_gp= gpar (fontsize=6), 
    right_annotation = ha,
    cluster_rows = T, #row_km = 2, column_km = 2,
    clustering_distance_rows = 'pearson', 
    clustering_distance_columns = 'pearson',
    column_title = paste('Binned',ws))
  
  png (paste0('Plots/binned_fragments_binned_',ws,'_celltype_heatmap.png'),width=600, height=500)
  binH
  dev.off()
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

run_chromVAR = FALSE

if (run_chromVAR)
  {  
  archp = addBgdPeaks (archp, force= FALSE)
  archp = addMotifAnnotations (ArchRProj = archp, 
      motifSet = "cisbp", 
      #motifSet = 'JASPAR2020',
      #name = "JASPAR2020_Motif",
      force=FALSE)
  archp = addDeviationsMatrix (
    ArchRProj = archp, 
    peakAnnotation = "Motif",
    force = FALSE
  )
  
  archp = saveArchRProject (ArchRProj = archp,  
      load = TRUE)
  }


##########################################
### Find activating and repressing TFs ###
##########################################
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
  
if (run_activeTF)
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
  activators = readRDS ('TF_activators_genescore.rds') 
  }


### Use P2G analysis and cNMF from RNA to identify active TF via regulons  ####
run_p2g_TF = FALSE

if (run_p2g_TF)
  {
  run_p2g = FALSE  
  if (run_p2g)
    {
    maxDist = 250000
    archp = addPeak2GeneLinks(
        ArchRProj = archp,
        useMatrix = 'GeneScoreMatrix',
        reducedDims = "IterativeLSI",
        maxDist = maxDist
    )
    }  
    
  p2g_corr = .2
  p2g = getPeak2GeneLinks(
      ArchRProj = archp,
      corCutOff = p2g_corr,
      resolution = 1,
      returnLoops = FALSE
  )
  
  # # Convert df in Granges add gene Name and correlation
  # gene = 'WT1'
  p2g$geneName = mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
  p2g = p2g[!is.na (p2g$FDR),] # remove NaN correlations (not sure why there are some)
  p2g = p2g[p2g$Correlation > p2g_corr, ]
  p2gGR = metadata (p2g)$peakSet[p2g$idxATAC]
  p2gGR$geneName = p2g$geneName
  p2gGR$correlation = p2g$Correlation
  

  # Import cNMF results and intersect with p2g
  nfeat=5000
  k=25
  cnmf_list = readRDS (paste0('../scrna/cnmf_genelist_',k,'_nfeat_',nfeat,'.rds'))
  cnmf_list = lapply (cnmf_list, function(x) head (x,200))
  p2g_cnmf = lapply (cnmf_list, function(x) p2gGR[p2gGR$geneName %in% x])
  
  tf_match = getMatches (archp)
  colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
  bg_peakSet = rowRanges (tf_match)
  #tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
  nmf_TF = lapply (p2g_cnmf, function(x) hyperMotif (
    selected_peaks = x, 
    motifmatch = tf_match))
  
  nmf_TF = lapply (nmf_TF, function(x) x[rownames(nmf_TF[[1]]), ])
  nmf_TF_df = do.call (cbind, nmf_TF)
  nmf_TF_df = nmf_TF_df[, grep ('padj', colnames(nmf_TF_df))]
  nmf_TF_df[nmf_TF_df > 0.05] = 1
  nmf_TF_df = -log10(nmf_TF_df)
  nmf_TF_df = nmf_TF_df[rowSums (nmf_TF_df) != 0, ]
  nmf_TF_df[sapply(nmf_TF_df, is.infinite)] <- 300
  saveRDS (nmf_TF_df, paste0('nmf_TF_p2g_enrichments.rds'))
  TF_ht = Heatmap (nmf_TF_df, row_names_gp = gpar (fontsize=3), column_names_gp = gpar (fontsize=5))
  
  pdf (paste0('Plots/TF_nmf_',k,'_nfeat_',nfeat,'_heatmap.pdf'),width = 3,height=25)
  print (TF_ht)
  dev.off()
  } else {
  nmf_TF_df = readRDS (paste0('nmf_TF_p2g_enrichments.rds'))
  }

# ### Co-expression of TFs ###
metaGroupName = 'Sample2'
if (!any (ls() == 'mSE')) mSE = ArchR::getMatrixFromProject (archp, useMatrix = 'MotifMatrix', logFile=NULL)
mSE = mSE[, archp$cellNames]
all (colnames(mSE) == rownames(archp))

# # Get deviation matrix and subset for relevant TF and aggregate
archp_meta = as.data.frame (archp@cellColData)
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

mMat = mMat[activators,]

mMat_agg = as.data.frame (t(mMat))
mMat_agg$metaGroup = as.character (archp_meta[,metaGroupName])
mMat_agg = aggregate (.~ metaGroup, mMat_agg, mean)
rownames (mMat_agg) = mMat_agg[,1]
mMat_agg = mMat_agg[,-1]
mMat_agg = t(mMat_agg)
#rownames (mMat_agg) = active_TF

d = as.dist (1-cor(as.matrix(tf_gs_cor$motifAvKnn_mean[,activators])))
d = as.dist (1-cor(t(as.matrix(mMat_agg))))
#d = dist (mMat, method ='euclidean')
hc1 <- hclust(d, method = "complete" ) # Hierarchical clustering using Complete Linkage
cor_TF = cor (as.matrix(tf_gs_cor$motifAvKnn_mean[,activators]))
cor_TF = cor (t(as.matrix(mMat_agg)))

row_filt = rowSums (cor_TF) != 0
tf_name = rownames(cor_TF)[row_filt]
km = kmeans (cor_TF[tf_name,tf_name], centers=9)
#km_ha = rowAnnotation (
#  km = anno_simple(as.character(km$cluster), width = unit(2, "mm"),
#    col = setNames (km_col, as.character(unique(km$cluster))), border=T))


metaGroupName = 'sampleID'
DefaultAssay(srt) = 'RNA'
sample_names_rna = c('P1','P4','P8','P3','P5','P10','P11','P12','P13','P14','HU37','HU62')
ps = log2(as.data.frame (AverageExpression (srt, features = activators, group.by = metaGroupName)[[1]]) +1)
ps = ps[, colnames(ps) %in% sample_names_rna]

# cor_TF_hm = Heatmap (cor_TF[tf_name,tf_name],
#         cluster_rows = T,
#         cluster_columns= T,
#         #left_annotation = km_ha,
#         #row_split = km$cluster,
#         column_title = 'all_MPM',
#         column_gap = unit(.2, "mm"),
#         row_gap = unit(.2, "mm"),
#         name = 'r',
#         #column_km_repeats=20,
#         col = viridis::viridis(100),
#         #clustering_distance_rows = 'pearson', 
#         #clustering_distance_columns = 'pearson', 
#         #column_split = km$cluster,
#         #col = pal_corr1,
#         row_names_gp = gpar(fontsize = 0),
#         column_names_gp = gpar(fontsize = 0),
#         width = unit(2, "cm"))

column_split = ifelse (grepl ('normal_pleura', colnames(mMat_agg)), 'Normal','Tumor')
TF_cluster_hm = Heatmap (t(scale(t(mMat_agg[tf_name,]))),
        #right_annotation=tf_mark,
        column_split = column_split,
        cluster_rows = T, #km = 4, 
        name = 'z-score\ndeviations',
        column_gap = unit(.2, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T, 
        col = viridis::magma (100),
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 5), 
        border=F,
        width = unit(2, "cm"))

#ps_order = ps[tf_name,sarc_order$sampleID]
column_split_rna = ifelse (grepl ('HU', colnames(ps)), 'Normal','Tumor')
TF_exp_hm = Heatmap (t(scale(t(ps))),
        #right_annotation=tf_mark,
        column_split = column_split_rna,
        cluster_rows = F, #km = 4, 
        name = 'expression',
        column_gap = unit(.2, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T, 
        col = viridis::mako (100),
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 5), 
        border=F,
        width = unit(2, "cm"))

TF_regulons_hm = Heatmap (t(scale(t(nmf_TF_df[tf_name,]))),
        #right_annotation=tf_mark,
#        column_split = column_split,
        cluster_rows = F, #km = 4, 
        name = 'expression',
        column_gap = unit(.2, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T,
        column_names_rot = 45, 
        col = viridis::inferno (100),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5), 
        border=F,
        width = unit(5, "cm"))

# # Add deviations from TCGA ATAC bulk
# tcga_dev = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/TCGA_atac/tcga_meso_atac_deviations_clinical_info.rds')
# tcga_mat = assays (tcga_dev)$z
# tcga_mat = aggregate (t(tcga_mat), by = list(sampleID= colnames(tcga_mat)), mean)
# rownames (tcga_mat) = tcga_mat[,1]
# tcga_mat = tcga_mat[,-1]
# tcga_mat = as.data.frame (t(tcga_mat))
# rownames (tcga_mat) = unname (sapply (rownames(tcga_mat), function(x) unlist (strsplit (x, '_'))[3]))
# #tcga_mat = rowMeans (tcga_mat)
# tcga_mat = tcga_mat[tf_name,]
# colnames (tcga_mat) = unname (sapply (colnames(tcga_mat), function(x) unlist (strsplit (x, '_'))[2]))
# heat_bulkATAC = Heatmap (t(scale(t(tcga_mat[tf_name,]))), 
#        clustering_distance_columns = 'pearson', 
#        cluster_columns=T, 
#        #col = pals_heatmap[[2]],
#        row_names_gp = gpar(fontsize = 5), 
#        name = 'z-score\ndeviations',
#        column_names_gp = gpar(fontsize = 6),
#        #rect_gp = gpar(col = "white", lwd = .5),
#        col = viridis::mako(10),
#        border=TRUE, width = 2)#, left_annotation=bulk_box_ha)

pdf (paste0 ('Plots/TF_cancer_modules_heatmaps2.pdf'), width = 18,height=25)
draw (TF_cluster_hm + TF_exp_hm + TF_regulons_hm)
dev.off()


# # make oncogenic signatures and plot on UMAP
# # for (i in tf_name)
# # {
# #   # archp = addCellColData (archp, name = colnames(malig_modules_dev)[i], 
# #   #   data = malig_modules_dev[,i], cells =rownames(malig_modules_dev), force=TRUE)
# # }

# tf_name2 = sapply (tf_name, function(x) rowData (mSE)$name[grepl (paste0('^',x,'_'), rowData (mSE)$name)])
# tf_name2 = paste0('z:',tf_name2)
# TF_p = plotEmbedding (
#     ArchRProj = archp,
#     colorBy = "MotifMatrix",
#     name = tf_name2, 
#     useSeqnames='z',
#     embedding = "UMAP")

# png (paste0(projdir,'Plots/coexp_TF_UMAPs.png'),
#    width=16800,height=16200, res=300)
# wrap_plots (TF_p)
# dev.off()

# tf_name2 = unlist(sapply (c('SOX9','TWIST1','MESP1','NKX2-5'), function(x) rownames(assay(mSE))[grepl (x, rownames(assay(mSE)))]))
# tf_name2 = paste0('z:',tf_name2)
# TF_p = plotEmbedding (
#     ArchRProj = archp,
#     colorBy = "MotifMatrix",
#     name = tf_name2, 
#     useSeqnames='z',
#     embedding = "UMAP")

# tf_name2 = c('SOX9','TWIST1','MESP1','NKX2-5')
# TF_p2 = plotEmbedding (
#     ArchRProj = archp,
#     colorBy = "GeneScoreMatrix",
#     name = tf_name2, 
#     #useSeqnames='z',
#     embedding = "UMAP")

# png (paste0('Plots/coexp_selected_TF_UMAPs.png'),
#    width=16800,height=16200, res=300)
# wrap_plots (c(TF_p, TF_p2))
# dev.off()

# # Show distribution per sample
# ccomp_df = as.data.frame (mMat_agg)
# ccomp_df$TF = rownames (ccomp_df)
# ccomp_df = gather (ccomp_df, 'sample','deviation', c(1:ncol(ccomp_df)-1))
# ccomp_df$status = ifelse (grepl ('N',ccomp_df$sample), 'Normal','Tumor')
# library (dplyr)
# tf_avg = ccomp_df %>% group_by (TF) %>% summarise_at(vars(deviation),
#                 ~summary(.)[[3]])
# ccomp_df$TF = factor (ccomp_df$TF, levels = tf_avg$TF[order(-tf_avg$deviation)])
# ccomp_df = ccomp_df[ccomp_df$sample != 'N3',]
# bp = ggplot (ccomp_df, aes (x= TF, y= deviation), fill = status) +
#             geom_boxplot (aes (fill = status),trim=TRUE) +
#             #geom_boxplot (aes (fill = status)) +
#             #geom_jitter (color="black", size=0.4, alpha=0.9) +
#             theme_classic() + 
#             #scale_fill_manual (values= module_pal) + 
#             ggtitle ('Mean TF activity per sample') + 
#             theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
#     #if (length(unique(metaGroupNames)) == 3) box_p = lapply (seq_along(meta_modules_names), function(x) box_p[[x]] + facet_wrap (as.formula(paste("~", metaGroupNames[3]))) + theme_classic()+ 
#     #        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))) 

# library (rstatix)
# stat.test2 = ccomp_df %>% group_by (TF) %>%
#   t_test(reformulate('status', 'deviation'), paired=F) %>%
#   adjust_pvalue(method = "fdr") %>%
#   add_significance()
# stat.test2 = stat.test2 %>% add_xy_position (x = "TF", step.increase=0.05)
# bp = bp + stat_pvalue_manual (stat.test2, remove.bracket=FALSE,
#    bracket.nudge.y = 0, hide.ns = TRUE,
#     label = "p.adj.signif") 

# pdf (paste0 ('Plots/Sample_distribution_coexp_TF.pdf'), height = 6,width=24)
# bp
# dev.off()


# # Make scatterplot of normal deviation vs median tumor deviations of active TF and the same for RNA
# #mMat_agg_scaled = t(scale(t(mMat_agg)))
# TF_diff = data.frame (
#   tumor = apply (mMat_agg[,unique(archp$Sample2)[!unique(archp$Sample2) == 'normal_pleura']], 1, mean),
#   normal = mMat_agg[,unique(archp$Sample2)[unique(archp$Sample2) == 'normal_pleura']])
# diff_line = 0
# TF_diff$label = ''
# TF_diff$label[head (order(-(TF_diff$tumor - TF_diff$normal)),10)] =  head (rownames(TF_diff)[rev(order(TF_diff$tumor - TF_diff$normal))],10)
# TF_diff$color = TF_diff$tumor > diff_line & TF_diff$normal < diff_line
# tf_diff_p = ggplot (TF_diff, aes (x= tumor, color = color, y = normal, label = label)) + 
#   geom_point( size = .5) + # Color points based on x value
#   scale_color_manual(values = c('FALSE' = "grey", 'TRUE' = "red")) + # Customize colors
#   geom_vline(xintercept = diff_line, linetype = "dashed", color = "grey44") + # Vertical dashed line
#   geom_hline(yintercept = diff_line, linetype = "dashed", color = "grey44") + # Vertical dashed line
#   gtheme_no_rot + # Use a minimal theme
#    geom_text_repel(
#     segment.size=.2,
# #    point.padding = 0.2, 
#     size=1#,
# #   nudge_x = .25,
# #    nudge_y = .2,
# #    segment.curvature = -1e-20
#     ) +
#     xlab ('median tumor dev') + 
#     ylab ('normal dev')

# pdf (paste0 ('Plots/Diff_normal_tumor_deviation_var_scatterplot.pdf'),4,height = 3)
# tf_diff_p
# dev.off()

# Make scatterplot of deviation difference and expression difference between normal and tumors ####
metaGroupName = 'Sample2'
archp_meta = as.data.frame (archp@cellColData)
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat_agg = as.data.frame (t(mMat))
mMat_agg$metaGroup = as.character (archp_meta[,metaGroupName])
mMat_agg = aggregate (.~ metaGroup, mMat_agg, mean)
rownames (mMat_agg) = mMat_agg[,1]
mMat_agg = mMat_agg[,-1]
mMat_agg = t(mMat_agg)
mMat_agg = mMat_agg[rownames(mMat_agg) %in% rownames(srt),]
metaGroupName = 'sampleID'
DefaultAssay(srt) = 'RNA'
sample_names_rna = c('P1','P4','P8','P3','P5','P10','P11','P12','P13','P14','HU37','HU62')
ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat_agg), group.by = metaGroupName)[[1]]) +1)
ps = ps[, colnames(ps) %in% sample_names_rna]

TF_diff_rna = data.frame (
  tumor_dev = apply (mMat_agg[,unique(archp$Sample2)[!unique(archp$Sample2) == 'normal_pleura']], 1, mean),
  normal_dev = mMat_agg[,unique(archp$Sample2)[unique(archp$Sample2) == 'normal_pleura']],
  normal_rna = apply (ps[,c('HU37','HU62')], 1, median),
  tumor_rna = apply (ps[, !colnames(ps) %in% c('HU37','HU62')], 1, mean),
  genescore = corGSM_MM$cor[match(rownames(mMat_agg), corGSM_MM$GeneScoreMatrix_name)])
TF_diff_rna$dev_diff = TF_diff_rna$tumor_dev - TF_diff_rna$normal_dev
TF_diff_rna$rna_diff = TF_diff_rna$tumor_rna - TF_diff_rna$normal_rna
diff_line = 0
TF_diff_rna$label = ifelse (TF_diff_rna$dev_diff > 0 & TF_diff_rna$rna_diff > 0, rownames(TF_diff_rna),'') 
TF_diff_rna$label_color = TF_diff_rna$label != ''
TF_diff_rna$color = TF_diff_rna$dev_diff > diff_line & TF_diff_rna$rna_diff > diff_line
tf_diff_p = ggplot (TF_diff_rna, aes (x= dev_diff, y = rna_diff)) + 
  geom_point(aes(fill=genescore, color=color), size = .3, shape = 21, stroke=0.1) + # Color points based on x value
  scale_color_manual(values = c('FALSE' = "grey", 'TRUE' = "red")) + # Customize colors
  scale_fill_gradient(low = "white", high = "black") +
  #scale_fill_manual(values = c('FALSE' = "grey", 'TRUE' = "red")) + # Customize colors
  geom_vline(xintercept = diff_line, linetype = "dashed", color = "grey44", linewidth=.3) + # Vertical dashed line
  geom_hline(yintercept = diff_line, linetype = "dashed", color = "grey44", linewidth=.3) + # Vertical dashed line
  gtheme_no_rot + # Use a minimal theme
#    geom_text_repel(
#     segment.size=.2,
#     max.overlaps = 10000,
# #    point.padding = 0.2, 
#     size=1#,
# #   nudge_x = .25,
# #    nudge_y = .2,
# #    segment.curvature = -1e-20
#     ) +
    xlab ('deviation difference') + 
    ylab ('RNA difference') + 
    xlim (c(-0.2, .2)) + 
    ylim (c(-0.6, .6))

pdf (paste0 ('Plots/Diff_normal_tumor_deviation_and_rna_scatterplot.pdf'),3,height = 2)
tf_diff_p
dev.off()

# Compute significance per sample vs normal in RNA space ####
library (presto)
srt$sampleID2 = srt$sampleID
srt$sampleID2[srt$sampleID2 %in% c('HU37','HU62')] = 'normal'
comparisons = list(
  P1 = c('P1','normal'),
  P11 = c('P11','normal'),
  P12 = c('P12','normal'),
  P13 = c('P13','normal'),
  P14 = c('P14','normal'),
  P3 = c('P3','normal'),
  P4 = c('P4','normal'),
  P5 = c('P5','normal'),
  P8 = c('P8','normal'))

res = lapply (comparisons, function(x) 
  wilcoxauc (srt[rownames(mMat_agg),], group_by = 'sampleID2', groups_use = x))
res_df = lapply (names (comparisons), function(x) res[[x]][res[[x]]$group == x,'padj',drop=FALSE])
res_df = do.call (cbind , res_df)
occurence_filter = apply (res_df, 1, function(x) sum (x < 0.01))
rownames (res_df) = rownames(mMat_agg)
colnames (res_df) = names(comparisons)
res_df_filtered = res_df[occurence_filter > 4, ]
tf_tumor_pos = TF_diff_rna$label[grepl('.',TF_diff_rna$label)]
res_df_filtered = res_df_filtered[rownames(res_df_filtered) %in% tf_tumor_pos,]
selected_TF = rownames(res_df_filtered)

palette_deviation = rev (as.character (paletteer::paletteer_c("grDevices::Rocket",100)))
palette_expression = rev (as.character(paletteer::paletteer_c("grDevices::Purple-Blue",100)))
palette_deviation = colorRamp2(c(-.2,-.1,0,.1,.2), c("white",'white','white', "brown",'black'))  

# Order by median logFC
res_df2 = lapply (names (comparisons), function(x) res[[x]][res[[x]]$group == x,'logFC',drop=FALSE])
res_df2 = do.call (cbind , res_df2)
rownames (res_df2) = rownames(mMat_agg)
tf_order = rownames(res_df2)[order (-apply (res_df2, 1, median))]
selected_TF_ordered = tf_order[tf_order%in% selected_TF]

column_split = ifelse (grepl ('normal_pleura', colnames(mMat_agg)), 'Normal','Tumor')
TF_cluster_selected_hm = Heatmap (mMat_agg[selected_TF_ordered,],
        #right_annotation=tf_mark,
        column_split = column_split,
        cluster_rows = F, #km = 4, 
        name = 'z-score\ndeviations',
        column_gap = unit(.8, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T, 
        col = palette_deviation,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 5), 
        border=T,
        width = unit(2, "cm"))

#ps_order = ps[tf_name,sarc_order$sampleID]
column_split_rna = ifelse (grepl ('HU', colnames(ps)), 'Normal','Tumor')
TF_exp_selected_hm = Heatmap (ps[selected_TF_ordered,],
        #right_annotation=tf_mark,
        column_split = column_split_rna,
        cluster_rows = T, #km = 4, 
        name = 'expression',
        column_gap = unit(.5, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T, 
        col = palette_expression,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 5), 
        border=T,
        width = unit(2, "cm"))


TF_regulons_hm = Heatmap (t(scale(t(nmf_TF_df[selected_TF_ordered,]))),
        #right_annotation=tf_mark,
#        column_split = column_split,
        cluster_rows = F, #km = 4, 
        name = 'expression',
        column_gap = unit(.2, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T,
        column_names_rot = 45, 
        col = viridis::inferno (100),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5), 
        border=F,
        width = unit(5, "cm"))

# Add deviations from TCGA ATAC bulk
tcga_dev = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/TCGA_atac/tcga_meso_atac_deviations_clinical_info.rds')
tcga_mat = assays (tcga_dev)$z
tcga_mat = aggregate (t(tcga_mat), by = list(sampleID= colnames(tcga_mat)), mean)
rownames (tcga_mat) = tcga_mat[,1]
tcga_mat = tcga_mat[,-1]
tcga_mat = as.data.frame (t(tcga_mat))
rownames (tcga_mat) = unname (sapply (rownames(tcga_mat), function(x) unlist (strsplit (x, '_'))[3]))
#tcga_mat = rowMeans (tcga_mat)
tcga_mat = tcga_mat[selected_TF,]
colnames (tcga_mat) = unname (sapply (colnames(tcga_mat), function(x) unlist (strsplit (x, '_'))[2]))
heat_bulkATAC = Heatmap (t(scale(t(tcga_mat[selected_TF,]))), 
       clustering_distance_columns = 'pearson', 
       cluster_columns=T, 
       #col = pals_heatmap[[2]],
       row_names_gp = gpar(fontsize = 5), 
       name = 'z-score\ndeviations',
       column_names_gp = gpar(fontsize = 6),
       column_names_rot = 45,
       #rect_gp = gpar(col = "white", lwd = .5),
       col = viridis::mako(10),
       border=TRUE, width = unit(2, "cm"))#, left_annotation=bulk_box_ha)

pdf (paste0 ('Plots/selected_TF_dev_exp_significant_heatmaps.pdf'), width = 7,height=9)
draw (TF_regulons_hm + TF_cluster_selected_hm + TF_exp_selected_hm)
dev.off()



# Make coexpression network for each sample using top TFs deviations
archp_meta = as.data.frame (archp@cellColData)
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

mMat = mMat[tf_name_selected,]
all (colnames(mMat) == rownames(archp_meta))
cor_TF_l = list()
for (sam in unique(archp_meta$Sample2))
  {
  cor_TF_l[[sam]] = cor (t(as.matrix(mMat[,archp_meta$Sample2 == sam])))
  cor_TF_l[[sam]] = Heatmap (cor_TF_l[[sam]], name = sam,
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 5))
  }

pdf (paste0 ('Plots/selected_TF_dev_corr_heatmaps.pdf'), width = 8,height=9)
cor_TF_l
dev.off()

# Take median of all correlations
# Sample data: list of matrices
all (colnames(mMat) == rownames(archp_meta))
cor_TF_l = list()
for (sam in unique(archp_meta$Sample2)) cor_TF_l[[sam]] = cor (t(as.matrix(mMat[,archp_meta$Sample2 == sam])))

cor_TF <- cor_TF_l[[1]]
cor_TF[] <- tapply(unlist(cor_TF_l), rep(seq(length(cor_TF_l[[1]])),length(cor_TF_l)), FUN=median)

cor_TF = Heatmap (cor_TF,
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 5))

pdf (paste0 ('Plots/selected_TF_dev_corr_heatmap.pdf'), width = 8,height=9)
cor_TF
dev.off()
 

# Export selected TFs
write.csv (tf_name_selected, 'Active_TFs.csv')





# Order cells per samples along SOX9 deviation and plot the rest of TF deviations together
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = mMat[tf_name_selected,]
archp_meta = as.data.frame (archp@cellColData)
all (colnames(mMat) == rownames(archp_meta))

traj_sample = list()
for (sam in unique(archp_meta$Sample2))
    {
    mMat_ordered_sample = mMat[,archp_meta$Sample2 == sam]
    mMat_ordered_sample = mMat_ordered_sample[, order(mMat_ordered_sample['MESP1',])]
    traj_sample[[sam]] = Heatmap (
      t(scale(t(mMat_ordered_sample))), 
      col = viridis::plasma(100), 
      cluster_columns=F,
      row_names_gp = gpar(fontsize = 4))
    }

pdf ('Plots/sarc_trajectory_per_sample.pdf', height=10)
traj_sample
dev.off()








# Export bigiwg files ####
metaGroupName = 'Sample2'
getGroupBW(
  ArchRProj = archp,
  groupBy = metaGroupName,
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)




















metaGroupName = 'Clusters'
DAM = getMarkerFeatures (
    ArchRProj = archp, 
    testMethod = "wilcoxon",
    useGroups = unique(archp@cellColData[,metaGroupName]),
    #bgdGroups = normal_meso,
    binarize = FALSE,
    useMatrix = paste0(TF_db,"Matrix"),
    groupBy = metaGroupName,
    useSeqnames="z")

top_tfs = 20
pValThreshold = 0.01
DAM_top = lapply (colnames (DAM), function(x) 
  {
  df = data.frame (fdr = unlist(assays(DAM[,x])$FDR), 
                  meandiff = unlist(assays(DAM[,x])$MeanDiff),
                  name = rowData(DAM)$name,
                  celltype = x)
  df = df[df$fdr < pValThreshold,]
  df = df[df$meandiff > 0,]
  df = head (df[order (df$fdr),], top_tfs)
  df$name = gsub ('_.*','',df$name) 
  df$name = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", df$name)
  df = df [df$name %in% activators, ]
  df 
  })
DAM_top = do.call (rbind, DAM_top)  
order_celltypes = c('TNKcells','Bcells','Myeloid','Endothelial','Fibroblasts','SmoothMuscle','mesothelium','malignant','sarcomatoid')
DAM_top$celltype = factor (DAM_top$celltype, levels = order_celltypes) 
DAM_top = DAM_top[order (DAM_top$celltype),]
metaGroupName = 'celltype2'
#metaGroupNameElements = unique(archp@cellColData[,metaGroupName])[!unique(archp@cellColData[,metaGroupName]) %in% c('N1','not_assigned')]
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMatDam = mMat[DAM_top$name, ]
mMatDam = as.data.frame (t(mMatDam))
mMatDam$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMatDam = aggregate (.~ metaGroup, mMatDam, mean)
rownames (mMatDam) = mMatDam[,1]
mMatDam = mMatDam[ levels (DAM_top$celltype),-1]
mMatDam = t(mMatDam)

# Add geneScore of TFs
gsMat =  assays (gsSE)[[1]]
rownames (gsMat) = rowData (gsSE)$name
gsMat = gsMat[rownames(gsMat) %in% rownames(mMatDam),]
gsMat = as.data.frame (t (gsMat))
gsMat$cluster = as.factor(getCellColData (archp)[,metaGroupName])
gsMat = aggregate (formula= .~ cluster, data= gsMat, FUN=mean)
rownames (gsMat) = gsMat[,1]
gsMat = gsMat[,-1]
gsMat = as.data.frame (t(gsMat))
gsMat = gsMat[rownames (mMatDam),]
gsMat = gsMat[,colnames(mMatDam)]

# Add deviations from TCGA ATAC bulk
# tcga_dev = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/TCGA_atac/tcga_meso_atac_deviations_clinical_info.rds')
# tcga_mat = assays (tcga_dev)$z
# tcga_mat = aggregate (t(tcga_mat), by = list(sampleID= colnames(tcga_mat)), mean)
# rownames (tcga_mat) = tcga_mat[,1]
# tcga_mat = tcga_mat[,-1]
# tcga_mat = as.data.frame (t(tcga_mat))
# rownames (tcga_mat) = unname (sapply (rownames(tcga_mat), function(x) unlist (strsplit (x, '_'))[3]))
# #tcga_mat = rowMeans (tcga_mat)
# tcga_mat = tcga_mat[rownames(mMatDam),]
# colnames (tcga_mat) = unname (sapply (colnames(tcga_mat), function(x) unlist (strsplit (x, '_'))[2]))

# Add scRNA 
#rna_mat = gcdata@assays$RNA@data
#rna_ct = as.character (gcdata$celltype_2)
#rna_ct[rna_ct == 'Malig'] = paste0 (rna_ct[rna_ct == 'Malig'],'_', gcdata$sampleid[gcdata$celltype_2 == 'Malig']) 
#rna_mat = rna_mat[rownames(rna_mat) %in% active_DAMs,]
#rna_mat = as.data.frame (t(as.matrix(rna_mat)))
#rna_mat$cluster = rna_ct
#rna_mat = aggregate (formula= .~ cluster, data= rna_mat, FUN= mean)
#rownames (rna_mat) = rna_mat$cluster
#rna_mat = rna_mat[,-1]
#rna_mat2 = apply (rna_mat, 1, function(x) x[match(rownames(mMatDam), colnames(rna_mat))])
#rna_mat2 = rna_mat2[,grepl('Malig', colnames (rna_mat2)) | grepl('Sarco', colnames (rna_mat2))]

# Add Viestra et al. motif clusters as heatmap annotation
# mclust = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/DBs/Viestra_motif_clusters/viestra_motif_clusters.csv')
# mclust$tf = sapply (mclust$Motif, function(x) unlist(strsplit(x, '_'))[1])
# mclust = mclust[!duplicated(mclust$tf),]
# mclustS = paste0('cl',(mclust[match(rownames (mMatDam), mclust$tf) ,c(1)]))

# mclust_ha = HeatmapAnnotation (Viestra_et_al.motif_cluster = anno_simple (mclustS, pch = mclustS, 
#   col=setNames(rep('white',length(unique(mclustS))),unique(mclustS)) ,
# pt_size = unit(.15, "cm")),  annotation_name_gp = gpar(fontsize=10),
# which='row', show_annotation_name =T, show_legend=F, simple_anno_size = unit(.1, "cm"), gap=unit(.1, "cm"))


tf_name_ha = rowAnnotation (max = anno_barplot(rowMax(mMatDam),location = 2,axis_param = list(
        facing = "outside"),baseline="min"),TF = anno_text(rownames(mMatDam), location = 1, rot = 30, 
    just = "right", gp = gpar(fontsize = 6)))

#km_cols = kmeans(t(mMatDam), 2)
#km_cols = grepl ('^N', colnames (mMatDam))

heat_motifs = Heatmap (t(scale(t(mMatDam))),
        #column_order = c()
        #left_annotation = tf_name_ha, 
        clustering_distance_rows = 'pearson',
        cluster_columns=F,
        cluster_rows=F,
        name = 'chromVAR',# row_km=4,
        #col = pals_heatmap[[2]],
        col = viridis::magma (10), #column_split = km_cols,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 8),
        #rect_gp = gpar(col = "white", lwd = .5),
        border=TRUE, row_dend_side  = 'left')
heat_gs = Heatmap (t(scale(t(gsMat))), 
        #column_split = km_cols,
        #clustering_distance_columns = 'Pearson', 
        cluster_columns=F,#col = pals_heatmap[[1]],
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 0),
        name = 'GeneScore',
        #col = pals_heatmap[[1]],
        col = viridis::viridis(10),
        #rect_gp = gpar(col = "white", lwd = .5),
        border=TRUE
        #right_annotation = motif_ha
        )
# heat_bulkATAC = Heatmap (t(scale(t(tcga_mat))), 
#         clustering_distance_columns = 'pearson', 
#         cluster_columns=T, 
#         #col = pals_heatmap[[2]],
#         row_names_gp = gpar(fontsize = 5), 
#         name = 'bulk_ATAC',
#         column_names_gp = gpar(fontsize = 6),
#         #rect_gp = gpar(col = "white", lwd = .5),
#         col = viridis::mako(10),
#         border=TRUE)#, left_annotation=bulk_box_ha)
#heat_scRNA = Heatmap (t(scale(t(rna_mat2))), 
#        clustering_distance_columns = 'pearson',
#        right_annotation = mclust_ha,
#        cluster_columns=T, 
#        #col = pals_heatmap[[2]],
#        row_names_gp = gpar (fontsize = 0), 
#        name = 'scRNA',
#        column_names_gp = gpar (fontsize = 6),
#        #rect_gp = gpar(col = "white", lwd = .5),
#        col = viridis::rocket(10),
#        border=TRUE)#, left_annotation=bulk_box_ha)

pdf (paste0(projdir, 'Plots/ChromVAR_noDoublets_motifs_gene_score_heatmaps_filtered_activators_cancer_vs_normal_tfdb_',TF_db,'_',devMethod,'.pdf'), width=4, height=8)
draw (heat_motifs + heat_gs , padding = unit(c(40, 2, 2, 2), "mm"),
  column_title = paste('DAM padj ',pValThreshold,' activTreshold',actTF_threshold), column_title_gp = gpar(fontsize = 10))
dev.off()

### Correlate mutational profiles of Sarcomatoid and epithelioid with normal cells ###
wgs = readRDS ("/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/WGS/Waddell_cohort/mutation_x_patient_grange.rds")

# ws = 1e6
# ss = 2e5
# if (!file.exists (paste0(projdir, 'bins_',ws,'_ss_',ss,'.rds')))
#   {
#   blacklist = toGRanges (paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/blacklisted_regions/ENCODE_blacklist/',"hg38-blacklist.v2.bed")) # taken from https://github.com/Boyle-Lab/Blacklist/tree/master/lists
#   windows = makeWindows (genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist,
#     windowSize = ws, slidingSize = ss)
#   saveRDS (windows, paste0(projdir, 'bins_',ws,'_ss_',ss,'.rds'))
#   } else {
#   windows = readRDS (paste0(projdir, 'bins_',ws,'_ss_',ss,'.rds'))  
#   }
# genome = getSeq(BSgenome.Hsapiens.UCSC.hg38, c(paste0('chr',1:22)))
# chromSizes = GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
# chromSizes = GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")  
# windows = subsetByOverlaps(windows, chromSizes)

# windows = unlist (slidingWindows (chromSizes, width = ws, step = ss))
projdir_normal = '/ahg/regevdata/projects/ICA_Lung/Bruno/Normal_lung/archr_RPL/'
archp_normal = loadArchRProject (projdir_normal)
archp_normal$celltype = archp_normal$predictedGroup_Un

# Load fragments
fragments_normal = unlist (getFragmentsFromProject (
      ArchRProj = archp_normal))
fragments_normal_flt = unlist(GRangesList(lapply (unique (archp_normal$celltype), function(x){
  tssranked_cells = archp_normal$cellNames[archp_normal$celltype == x][order(-archp_normal$TSSEnrichment[archp_normal$celltype == x])]
  ncells = 1
  nfrags = 0
  while (nfrags < 1e6 & ncells <= length(archp_normal$cellNames[archp_normal$celltype == x]))
    {
     nfrags = length (fragments_normal[fragments_normal$RG %in% c(tssranked_cells[1:ncells])])
     ncells = ncells + 1 
    }
  message (paste(x,'ncells reaching 1e6 frag:', ncells))
  fragments_normal[fragments_normal$RG %in% tssranked_cells[1:ncells]]
  })))

library ("liftOver")
path = system.file (package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
seqlevelsStyle(fragments_normal_flt) = "UCSC"  # necessary
fragments_normal_flt19 = unlist (liftOver(fragments_normal_flt, ch))

#fragments_bin = lapply (unique (archp_normal$celltype), function(x) countOverlaps (windows,fragments_normal_flt[fragments_normal_flt$RG %in% archp_normal$cellNames[archp_normal$celltype == x]]))
# wgs_epit = countOverlaps(windows, wgs[grep('Epithelioid', names(wgs))])
# wgs_sarc = countOverlaps(windows, wgs[grepl('Sarc', names(wgs)) | grepl('Biphasic', names(wgs))])


# fr_mut_cor_epit = lapply (fragments_bin, function(x) cor (x, wgs_epit, method='spearman',use='complete.obs'))
# names (fr_mut_cor_epit) = unique (archp_normal$celltype)
# fr_mut_cor_epit
# fr_mut_cor_sarc = lapply (fragments_bin, function(x) cor (x, wgs_sarc, method='spearman',use='complete.obs'))
# names (fr_mut_cor_sarc) = unique (archp_normal$celltype)
# fr_mut_cor_sarc

projdir2 = '/ahg/regevdata/projects/ICA_Lung/Bruno/celloforigin_prj/'
load (paste0(projdir2,'hg19.1Mb.ranges.Polak.Nature2015.RData'))
interval.ranges = interval.ranges[!seqnames(interval.ranges) %in% c('chrX','chrY')]
wgs_epit = countOverlaps(interval.ranges, wgs[grep('Epithelioid', names(wgs))])
wgs_sarc = countOverlaps(interval.ranges, wgs[grepl('Sarc', names(wgs)) | grepl('Biphasic', names(wgs))])

fragments_bin = lapply (unique (archp_normal$celltype), function(x) countOverlaps (interval.ranges,fragments_normal_flt19[fragments_normal_flt19$RG %in% archp_normal$cellNames[archp_normal$celltype == x]]))

fr_mut_cor_epit = lapply (fragments_bin, function(x) cor (x, wgs_epit, method='spearman',use='complete.obs'))
names (fr_mut_cor_epit) = unique (archp_normal$celltype)
fr_mut_cor_epit
fr_mut_cor_sarc = lapply (fragments_bin, function(x) cor (x, wgs_sarc, method='spearman',use='complete.obs'))
names (fr_mut_cor_sarc) = unique (archp_normal$celltype)
fr_mut_cor_sarc

# filter out low-variablitiy bins #
require(DESeq2)
fr_mat = do.call (cbind, fragments_bin)
lib.size <- estimateSizeFactorsForMatrix(fr_mat)
ed <- t(t(fr_mat)/lib.size)
#Calculate estimates of variance, coefficient of variation
means <- rowMeans(ed)
vars <- apply(ed,1,var)
cv2 <- vars/means^2
pdf (paste0(projdir, 'COO_variable_bins.pdf'))
par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9)
smoothScatter(log(means),log(cv2))
dev.off()

flt_bin = 12600
lmeans_threshold = 4
cv2_threshold = -2
#binvar = apply (fr_mat,1, var)
fr_mat = do.call (cbind, fragments_bin)
fr_mat = fr_mat[log(means) > lmeans_threshold & abs (cv2) > cv2_threshold, ]
fragments_bin_flt = lapply (1:ncol(fr_mat), function(x) fr_mat[,x])
fr_mut_cor_epit = lapply (fragments_bin_flt, function(x) cor (x, wgs_epit[which(log(means) > lmeans_threshold), method='spearman',use='complete.obs'))
names (fr_mut_cor_epit) = unique (archp_normal$celltype)
fr_mut_cor_epit
fr_mut_cor_sarc = lapply (fragments_bin_flt, function(x) cor (x, wgs_sarc[which(log(means) > lmeans_threshold), method='spearman',use='complete.obs'))
names (fr_mut_cor_sarc) = unique (archp_normal$celltype)
fr_mut_cor_sarc

# Repeat correlation analysis with WGS from p786 and p848
wgs2 = read.csv ("/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/WGS/p786_p848/mesothelioma_p786_p848.csv")
fr_mut_cor2 = lapply (fragments_bin, function(x) cor (x[wgs2$X], wgs2[,3], method='spearman',use='complete.obs'))
names (fr_mut_cor2) = unique (archp_normal$celltype)

### Compute deviations on mutations + coaccessibility ###






# archp$DoubletScore[is.na (archp$DoubletScore)] = 0
# archp$DoubletEnrichment[is.na (archp$DoubletEnrichment)] = 0
# archp = filterDoublets (archp, filterRatio = 3) # filter doublets
# archp = archp[archp$nFrags < 1e5] # remove doublets manually
# archp = saveArchRProject (ArchRProj = archp, dropCells=TRUE, load = TRUE)

# # Dimensionality reduction and clustering
# set.seed (1234)
# archp = addIterativeLSI (
#   ArchRProj = archp, 
#   useMatrix = "TileMatrix", 
#   name = "IterativeLSI",
#   force=TRUE)
# archp = addUMAP (
#   ArchRProj = archp, 
#   reducedDims = "IterativeLSI",
#   force = TRUE)

# res = 1.2
# archp = addClusters (
#   input = archp, 
#   resolution = res,
#   reducedDims = "IterativeLSI", 
#   maxClusters = Inf,
#   force = TRUE)

# metaGroupName = 'Clusters'
# metaGroupNames = c('TSSEnrichment','nFrags','ReadsInTSS','DoubletScore','DoubletEnrichment',metaGroupName,'Sample2')  
# umap_p1 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
#  name = x, embedding = "UMAP"))
  
# pdf (paste0(projdir,'Plots/sample_clusters_doublets_removed_ReadsInTSS_',RITfilter,'normal_cf_sub_umap.pdf'), 15,15)
# wrap_plots (umap_p1, ncol=4)
# dev.off()


### SELECT ONLY TUMORS ###
# Dimensionality reduction and clustering
# set.seed (1234)
# archp2 = archp[grepl ('pt',archp$Sample2)]
# archp2$celltype3 = archp2$celltype
# archp2$celltype3[grepl ('cancer',archp2$celltype3)] = 'malignant'
# archp2 = archp2[archp2$celltype3 != 'cf_mesothelium']
# archp2 = archp2[archp2$celltype3 != 'p786_mesothelium']

# archp2 = addIterativeLSI (
#   ArchRProj = archp2, 
#   useMatrix = "TileMatrix", 
#   name = "IterativeLSI",
#   force=TRUE)
# archp2 = addUMAP (
#   ArchRProj = archp2, 
#   reducedDims = "IterativeLSI",
#   force = TRUE)

# res = 1.2
# archp2 = addClusters (
#   input = archp2, 
#   resolution = res,
#   reducedDims = "IterativeLSI", 
#   maxClusters = Inf,
#   force = TRUE)

# metaGroupName = 'celltype3'
# metaGroupNames = c('TSSEnrichment','nFrags','ReadsInTSS','DoubletScore','DoubletEnrichment',metaGroupName,'Sample2')  
# umap_p1 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp2, colorBy = "cellColData",
#  name = x, embedding = "UMAP", labelMeans=F))
  
# pdf (paste0(projdir,'Plots/MPM_sample_clusters_doublets_removed_ReadsInTSS_',RITfilter,'normal_cf_sub_umap.pdf'), 15,15)
# wrap_plots (umap_p1, ncol=4)
# dev.off()

# ### Gene score based analysis ####
# # Check genescore of celltype markers
# gsSE2 = gsSE
# gsMat2 = assays (gsSE2)[[1]]
# gsMat2 = gsMat2[, archp2$cellNames]
# #celltype_markers = c("PTPRC","LYZ","CD74","CD3E","TPSB2","CD79A","TAGLN","COL1A2","PECAM1","SFTPB","AGER","TPPP3","EPCAM","KRT19","WFDC2","TFF3",'CALB2','MSLN','WT1','GATA4','GATA6')
# celltype_markers = c('CALB2',"LYZ","CD3E","TAGLN","COL1A2","PECAM1")
# rownames (gsMat2) = rowData (gsSE2)$name
# gsMat2_mg = gsMat2[rownames (gsMat2) %in% celltype_markers, ]
# gsMat2_mg = as.data.frame (t(gsMat2_mg))
# gsMat2_mg$metaGroup = as.character(archp2@cellColData[,metaGroupName])
# gsMat2_mg = aggregate (.~ metaGroup, gsMat2_mg, mean)
# rownames (gsMat2_mg) = gsMat2_mg[,1]
# gsMat2_mg = gsMat2_mg[,-1]
# #mat_GeneScore2 = mat_GeneScore[, colnames(mat_GeneScore) %in% unique(markerGenes_top$name)]

# GeneScore2_hm = Heatmap (t(scale(gsMat2_mg)), 
#         row_labels = colnames (gsMat2_mg),
#         clustering_distance_columns = 'euclidean',
#         clustering_distance_rows = 'euclidean',
#         cluster_rows = T, name = 'GeneScore',
#         #col = pals_heatmap[[5]],
#         cluster_columns=T,#col = pals_heatmap[[1]],
#         row_names_gp = gpar(fontsize = 15),
#         column_names_gp = gpar(fontsize = 15),
#         rect_gp = gpar(col = "white", lwd = .5),
#         border=TRUE, col = viridis::viridis(100)
#         #right_annotation = motif_ha
#         )

# gs2_markers_grob = grid.grabExpr(draw(GeneScore2_hm, column_title = 'Celltype markers GeneScore', column_title_gp = gpar(fontsize = 16)))
# pdf (paste0(projdir, 'Plots/MPM_Genescore_celltype_markers_clusters_',metaGroupName,'heatmaps.pdf'), 15, 6)
# wrap_plots (umap_p1[[7]],umap_p1[[6]], gs2_markers_grob)
# dev.off()






# #### Integrate RNA annotation ####
# load ('/ahg/regevdata/projects/lungCancerBueno/Results/10x_MesoSCInt_201004/All/15PCs_gcdata.Rda')
# saveRDS (gcdata, paste0(projdir, 'scRNA_all_meso_maggie.rds'))
# gcdata = readRDS (paste0(projdir, 'scRNA_all_meso_maggie.rds'))

# archp = addGeneIntegrationMatrix (
#     ArchRProj = archp, 
#     useMatrix = "GeneScoreMatrix",
#     matrixName = "GeneIntegrationMatrix",
#     reducedDims = "IterativeLSI",
#     seRNA = gcdata,
#     force = TRUE,
#     addToArrow = TRUE,
#     groupRNA = "celltype_2",
#     nameCell = "predictedCell_Un",
#     nameGroup = "predictedGroup_Un",
#     nameScore = "predictedScore_Un"
# )    

# metaGroupNames = c(metaGroupName,'Sample2','predictedGroup_Un')  
# umap_p1 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
#  name = x, embedding = "UMAP"))
  
# pdf (paste0(projdir,'Plots/rna_integration_umap.pdf'), 15,15)
# wrap_plots (umap_p1, ncol=4)
# dev.off()

# # Annotate clusters
# archp$celltype = 0
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C1','C2','C3')] = 'myeloid'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C4','C5','C6','C7','C8','C9')] = 'lymphoid'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C10')] = 'normal_mesothelium'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C11','C12')] = 'fibroblasts'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C13')] = 'endothelial'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C14')] = 'endothelial'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C15','C16')] = 'smoothMuscle'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C17')] = 'p786_cancer'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C18','C19')] = 'fibroblasts'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C20','C21') & archp$Sample2 == 'pt_848'] = 'p848_cancer'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C20','C21') & archp$Sample2 == 'pt_846'] = 'p846_cancer'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C22')] = 'p786_mesothelium'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C23')] = 'cf_mesothelium'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C24','C26','C27')] = 'p811_cancer'
# archp$celltype[as.character(archp@cellColData[,metaGroupName]) %in% c('C25')] = 'p826_cancer'

# # Annotate clusters
# archp$celltype2 = 0
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C1','C2','C3')] = 'myeloid'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C4','C5','C6','C7','C8','C9')] = 'lymphoid'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C10')] = 'normal_mesothelium'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C11','C12')] = 'fibroblasts'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C13')] = 'endothelial'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C14')] = 'endothelial'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C15','C16')] = 'smoothMuscle'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C17')] = 'p786_cancer'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C18','C19')] = 'fibroblasts'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C20','C21') & archp$Sample2 == 'pt_848'] = 'p848_cancer'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C20','C21') & archp$Sample2 == 'pt_846'] = 'p846_cancer'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C22')] = 'p786_mesothelium'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C23')] = 'cf_mesothelium'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C24')] = 'p811_C24_cancer'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C26')] = 'p811_C26_cancer'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C27')] = 'p811_C27_cancer'
# archp$celltype2[as.character(archp@cellColData[,metaGroupName]) %in% c('C25')] = 'p826_cancer'

# metaGroupName2 = 'celltype'
# metaGroupName3 = 'celltype2'
# metaGroupNames = c('Sample2','predictedGroup_Un',metaGroupName2,metaGroupName3)  
# umap_p1 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
#  name = x, embedding = "UMAP"))

# pdf (paste0(projdir,'Plots/celltype_doublet_removed_annotation_umap.pdf'), 18,15)
# wrap_plots (umap_p1, ncol=4)
# dev.off()

# celltypes_pal_I1 = paletteDiscrete(unique(archp$celltype), set='rushmore', reverse=T)
# cc_bar1 = cellComp (
#   seurat_obj = as.data.frame (archp@cellColData[,c('Sample2',metaGroupName2)]), 
#   metaGroups = c('Sample2',metaGroupName2),
#   plot_as = 'bar',
#   stat_test = 't.test',
#   stat_comparisons = NULL, # list(c("GROUP1","GROUP2"))
#   pal = celltypes_pal_I1,
#   scale.facet.box = TRUE
#   )
# umap_p2 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
#  name = metaGroupName2, embedding = "UMAP", labelMeans=FALSE,
#  pal = celltypes_pal_I1)

# pdf (paste0(projdir,'Plots/celltype_proportions_barplot.pdf'), 10,5)
# wrap_plots (umap_p2, cc_bar1, ncol=4)
# dev.off()

# # Check genescore of celltype markers
# gsMat = assays (gsSE)[[1]]

# celltype_markers = c("PTPRC","LYZ","CD74","CD3E","TPSB2","CD79A","TAGLN","COL1A2","PECAM1","SFTPB","AGER","TPPP3","EPCAM","KRT19","WFDC2","TFF3",'CALB2','MSLN','WT1','GATA4','GATA6')
# rownames (gsMat) = rowData (gsSE)$name
# gsMat_mg = gsMat[rownames (gsMat) %in% celltype_markers, ]
# gsMat_mg = as.data.frame (t(gsMat_mg))
# gsMat_mg$metaGroup = as.character(archp@cellColData[,metaGroupName2])
# gsMat_mg = aggregate (.~ metaGroup, gsMat_mg, mean)
# rownames (gsMat_mg) = gsMat_mg[,1]
# gsMat_mg = gsMat_mg[,-1]
# #mat_GeneScore2 = mat_GeneScore[, colnames(mat_GeneScore) %in% unique(markerGenes_top$name)]

# GeneScore_hm = Heatmap (t(scale(gsMat_mg)), 
#         row_labels = colnames (gsMat_mg),
#         clustering_distance_columns = 'euclidean',
#         clustering_distance_rows = 'euclidean',
#         cluster_rows = T,
#         #col = pals_heatmap[[5]],
#         cluster_columns=T,#col = pals_heatmap[[1]],
#         row_names_gp = gpar(fontsize = 10),
#         column_names_gp = gpar(fontsize = 10),
#         rect_gp = gpar(col = "white", lwd = .5),
#         border=TRUE, name = " "
#         #right_annotation = motif_ha
#         )

# pdf (paste0(projdir, 'Plots/Genescore_celltype_markers_clusters_',metaGroupName2,'heatmaps.pdf'), 5, 6)
# draw(GeneScore_hm, column_title = 'Celltype markers GeneScore', column_title_gp = gpar(fontsize = 16))
# dev.off()

# # Export bigiwg files
# getGroupBW(
#   ArchRProj = archp,
#   groupBy = metaGroupName2,
#   normMethod = "ReadsInTSS",
#   tileSize = 100,
#   maxCells = 1000,
#   ceiling = 4,
#   verbose = TRUE,
#   threads = getArchRThreads(),
#   logFile = createLogFile("getGroupBW")
# )

# # Coverage plots of meso developmental TFs
# meso_markers = c('WT1','GATA4','TEAD1','TEAD3','TEAD4','GATA6')
# p <- plotEmbedding(
#     ArchRProj = archp, 
#     colorBy = "GeneScoreMatrix", 
#     name = meso_markers, 
#     embedding = "UMAP",
#     imputeWeights = getImputeWeights(archp)
# )

# plotPDF(plotList = p, 
#     name = "Genescore_meso_MarkerGenes-W-Imputation.pdf", 
#     ArchRProj = archp, 
#     addDOC = FALSE, width = 5, height = 5)

# meso_markers <- plotBrowserTrack(
#     ArchRProj = archp, 
#     groupBy = metaGroupName2, 
#     geneSymbol = meso_markers, 
#     upstream = 200000,
#     downstream = 200000,
#     loops = getCoAccessibility (archp, corCutOff = 0.5,
#       returnLoops = TRUE),
#     useGroups= NULL
# )

# plotPDF (meso_markers,
#     name = "Meso_markers_coverageTracks.pdf", 
#     ArchRProj = archp, 
#     addDOC = FALSE, width = 8, height = 5)

# # DAG - Find differential genes cancer vs normal meso (normal pleura + cf)
# force=T
# if (file.exists (paste0 (projdir,'DAG_',metaGroupName2,'.rds')) & !force) 
#   {
#   DAG_list = readRDS (paste0 (projdir,'DAG_',metaGroupName2,'.rds'))
#   } else {
#   DAG_list = list()
#   for (i in unique(archp@cellColData[,metaGroupName2])[grep('cancer',unique(archp@cellColData[,metaGroupName2]))])
#     {
#     DAG_list[[i]] = getMarkerFeatures (
#       ArchRProj = archp, 
#       testMethod = "wilcoxon",
#       useGroups = i,
#       maxCells = 500,
#       bgdGroups = c('normal_mesothelium','cf_mesothelium'),
#       binarize = FALSE,
#       useMatrix = "GeneScoreMatrix",
#       groupBy = metaGroupName2,
#     #  useSeqnames="z"
#     )}
#     listnames = names(DAG_list)
#     DAG_list = lapply (seq_along (DAG_list), function(x) 
#       {
#       df = do.call (cbind, (assays(DAG_list[[x]])))
#       colnames(df) = names (assays(DAG_list[[x]]))
#       df$gene = rowData (DAG_list[[x]])$name
#       df
#       })
#     names (DAG_list) = listnames
#     saveRDS (DAG_list, paste0 (projdir,'DAG_',metaGroupName2,'.rds'))
#     } 

# FDR_threshold = 0.05
# lfc_threshold = .5
# top_genes = 20
# DAG_top_list = lapply (seq_along(DAG_list), function(x) {
#   res = DAG_list[[x]]
#   res = na.omit (res)
#   res = res[res$FDR < FDR_threshold,]
#   res = res[order (res$FDR), ]
#   res = res[abs(res$Log2FC) > lfc_threshold,]
#   res$comparison = names(DAG_list)[x]
#   if (nrow(res) < top_genes) 
#     {
#     res
#     } else {
#     head (res,top_genes)
#     }
#   })
# DAG_df = Reduce (rbind ,DAG_top_list)

# gsMat = assays (gsSE)[[1]]
# rownames (gsMat) = rowData (gsSE)$name
# gsMat_mg = gsMat[rownames (gsMat) %in% DAG_df$gene, ]
# gsMat_mg = as.data.frame (t(gsMat_mg))
# gsMat_mg$metaGroup = as.character(archp@cellColData[,metaGroupName2])
# gsMat_mg = aggregate (.~ metaGroup, gsMat_mg, mean)
# rownames (gsMat_mg) = gsMat_mg[,1]
# gsMat_mg = gsMat_mg[,-1]
# gsMat_mg = gsMat_mg[grepl ('mesothelium',rownames(gsMat_mg)) | grepl('cancer',rownames(gsMat_mg)),]
# #mat_GeneScore2 = mat_GeneScore[, colnames(mat_GeneScore) %in% unique(markerGenes_top$name)]

# DAG_hm = Heatmap (t(scale(gsMat_mg)), 
#         row_labels = colnames (gsMat_mg),
#         column_title = paste('top',top_genes),
#         clustering_distance_columns = 'euclidean',
#         clustering_distance_rows = 'euclidean',
#         cluster_rows = T,
#         #col = pals_heatmap[[5]],
#         cluster_columns=T,#col = pals_heatmap[[1]],
#         row_names_gp = gpar(fontsize = 2),
#         column_names_gp = gpar(fontsize = 4),
#         rect_gp = gpar(col = "white", lwd = .5),
#         border=TRUE
#         #right_annotation = motif_ha
#         )
       
# DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
# pdf (paste0(projdir, 'Plots/DAG_clusters_',metaGroupName2,'_heatmaps.pdf'), 6, 6)
# umap_p1[[3]] | DAG_grob
# dev.off()

# ### Call peaks ###
# metaGroupName = 'Clusters'
# archp = addGroupCoverages (
#   ArchRProj = archp, 
#   groupBy = metaGroupName,  
#   force = TRUE,
#   minCells= 20, # I think this should be set corresponding to the smallest cluster in the group or lower
#   maxCells = 500,
#   minReplicates = 2,
#   sampleRatio = 0.8,
#   useLabels = TRUE)

# archp = addReproduciblePeakSet (
#     archp,
#     groupBy= metaGroupName,
#     peakMethod = 'Macs2',
#     reproducibility = "1",
#     maxPeaks = 500000, 
#     minCells=20) # I think this should be set corresponding to the smallest cluster in the group or lower
# archp = addPeakMatrix (archp)

# archp = saveArchRProject (archp, load=TRUE)

# # DAG - Find differential genes cancer vs normal meso (normal pleura + cf)
# # Object is saved as standard list cause it takes too much space as SE (> 20GB, crazy)
# force = FALSE
# if (file.exists (paste0 (projdir,'DAP_',metaGroupName2,'.rds')) & !force) 
#   {
#   DAP_list = readRDS (paste0 (projdir,'DAP_',metaGroupName2,'.rds'))
#   } else {
#   DAP_list = list()
#   for (i in unique(archp@cellColData[,metaGroupName2])[grep('cancer',unique(archp@cellColData[,metaGroupName2]))])
#     {
#     DAP_list[[i]] = getMarkerFeatures (
#       ArchRProj = archp, 
#       testMethod = "wilcoxon",
#       useGroups = i,
#       maxCells = 500,
#       normBy = 'ReadsInTSS',
#       bgdGroups = c('normal_mesothelium','cf_mesothelium'),
#       binarize = FALSE,
#       useMatrix = "PeakMatrix",
#       groupBy = metaGroupName2,
#     #  useSeqnames="z"
#     )}
#     listnames = names(DAP_list)
#     DAP_list = lapply (seq_along (DAP_list), function(x) 
#       {
#       df = do.call (cbind, (assays(DAP_list[[x]])))
#       colnames(df) = names (assays(DAP_list[[x]]))
#       df$peak = paste0(rowData (DAP_list[[x]])[,1],':',rowData (DAP_list[[x]])[,3],'-',rowData (DAP_list[[x]])[,4])
#       df
#       })
#     names (DAP_list) = listnames
#     saveRDS (DAP_list, paste0 (projdir,'DAP_',metaGroupName2,'.rds'))
#     } 

# FDR_threshold = 1e-2
# lfc_threshold = 0.2
# top_peaks = Inf
# DAP_l_all = lapply (seq_along (DAP_list), function(x) {
#   res = DAP_list[[x]]
#   res = na.omit (res)
#   res = res[res$FDR < FDR_threshold,]
#   res = res[order (res$FDR), ]
#   if (nrow(res) < 1)
#     {
#     return (NULL)
#     } else {
#   res = res[abs(res$Log2FC) > lfc_threshold,]
#   res$comparison = names(DAP_list)[x]
#   if (nrow(res) < top_peaks) 
#     {
#     res
#     } else {
#     head (res,top_peaks)
#     }}
#   })
# DAP_df_all = Reduce (rbind ,DAP_l_all)

# DAP_l_neg = lapply (seq_along (DAP_list), function(x) {
#   res = DAP_list[[x]]
#   res = na.omit (res)
#   res = res[res$FDR < FDR_threshold,]
#   res = res[order (res$FDR), ]
#   res = res[res$Log2FC < -lfc_threshold,]
#   if (nrow(res) < 1)
#     {
#     return (NULL)
#     } else {
#   res$comparison = names(DAP_list)[x]
#   if (nrow(res) < top_peaks) 
#     {
#     res
#     } else {
#     head (res,top_peaks)
#     }
#   }})
# DAP_df_neg = Reduce (rbind ,DAP_l_neg)

# DAP_l_pos = lapply (seq_along (DAP_list), function(x) {
#   res = DAP_list[[x]]
#   res = na.omit (res)
#   res = res[res$FDR < FDR_threshold,]
#   res = res[order (res$FDR), ]
#   res = res[res$Log2FC > lfc_threshold,]
#   if (nrow(res) < 1)
#     {
#     return (NULL)
#     } else {
#   res$comparison = names(DAP_list)[x]
#   if (nrow(res) < top_peaks) 
#     {
#     res
#     } else {
#     head (res,top_peaks)
#     }
#   }})
# DAP_df_pos = Reduce (rbind ,DAP_l_pos)

# # Make UpSet plot of DAP
# DAP_pos_peaks = lapply (DAP_l_pos, function(x) x$peak)
# DAP_pos_ov = UpSetR::fromList (DAP_pos_peaks)
# colnames (DAP_pos_ov) = names (DAP_list)
# DAP_neg_peaks = lapply (DAP_l_neg, function(x) x$peak)
# DAP_neg_ov = UpSetR::fromList (DAP_neg_peaks)
# colnames (DAP_neg_ov) = names (DAP_list)

# pdf (paste0(projdir, 'Plots/DAP_overlap_cancers_vs_normal_upsetPlot.pdf'))
# upset(DAP_pos_ov, order.by = "freq", mainbar.y.label = "Intersection size in MPM high peaks'")
# upset(DAP_neg_ov, order.by = "freq", mainbar.y.label = "Intersection size in MPM low peaks'")
# dev.off()

# ### GREAT analysis of DAP
# # get background peaks (filter by MACS2 score to get less than 1M peaks)
# #ps = getPeakSet(archp)
# #ps = ps[ps$score > 2.5]
# #names (ps) = NULL

# #DAP_df_neg = DAP_df_neg[queryHits (findOverlaps (GRanges (DAP_df_neg$peak), ps)),]
# write.csv (DAP_df_neg, paste0 (projdir, 'DAP_neg_cancer_vs_normal.csv'))
# DAP_l2_neg = split (DAP_df_neg, DAP_df_neg$comparison)
# job_neg = lapply (DAP_l2_neg, function(x) submitGreatJob(GRanges (x$peak),# bg=ps, 
#   species='hg38'))
# tb_neg = lapply (job_neg, function(x) getEnrichmentTables(x))
# tbf_neg = lapply (tb_neg, function(x) x[[2]][x[[2]]$Hyper_Adjp_BH < 0.01,])

# # Make UpSet plot of DAP
# tbf_neg_terms = lapply (tbf_neg, function(x) x$name)
# tbf_neg_terms_ov = UpSetR::fromList (tbf_neg_terms)
# colnames (tbf_neg_terms_ov) = names (tbf_neg)

# pdf (paste0(projdir, 'Plots/DAP_GREAT_neg_cancers_vs_normal_upsetPlot.pdf'))
# upset(tbf_neg_terms_ov, order.by = "freq")
# dev.off()

# # Show terms shared in more than 2 cancers
# shared_neg_terms = names(table (unlist(tbf_neg_terms))[table(unlist(tbf_neg_terms)) > 2])
# shared_neg_terms2 = lapply (tbf_neg, function(x) 
#   {
#   y = x[,'Hyper_Adjp_BH', drop=F]
#   rownames(y) = x$name
#   y[shared_neg_terms, ,drop=F]
#   })
# shared_neg_terms2 = do.call (cbind, shared_neg_terms2)
# rownames (shared_neg_terms2) = shared_neg_terms
# shared_neg_terms2[is.na(shared_neg_terms2)] = 1
# colnames (shared_neg_terms2) = names (tbf_neg)
# shared_neg_terms_hm = Heatmap (-log10(shared_neg_terms2), col = viridis::rocket(100),
#   row_names_gp = gpar(fontsize=4), column_names_gp = gpar(fontsize=6))
# pdf (paste0(projdir, 'Plots/DAP_GREAT_neg_shared_terms_heatmap.pdf'), width=4, height=10)
# shared_neg_terms_hm
# dev.off()


# #DAP_df_pos = DAP_df_pos[queryHits (findOverlaps (GRanges (DAP_df_pos$peak), ps)),]
# write.csv (DAP_df_pos, paste0(projdir, 'DAP_pos_cancer_vs_normal.csv'))
# DAP_l2_pos = split (DAP_df_pos, DAP_df_pos$comparison)
# job_pos = lapply (DAP_l2_pos, function(x) submitGreatJob(GRanges (x$peak),# bg=ps, 
#   species='hg38'))
# tb_pos = lapply (job_pos, function(x) getEnrichmentTables(x))
# tbf_pos = lapply (tb_pos, function(x) x[[2]][x[[2]]$Hyper_Adjp_BH < 0.01,])

# # Make UpSet plot of DAP
# tbf_pos_terms = lapply (tbf_pos, function(x) x$name)
# tbf_pos_terms_ov = UpSetR::fromList (tbf_pos_terms)
# colnames (tbf_pos_terms_ov) = names (tbf_pos)

# pdf (paste0(projdir, 'Plots/DAP_GREAT_pos_cancers_vs_normal_upsetPlot.pdf'))
# upset(tbf_pos_terms_ov, order.by = "freq")
# dev.off()

# # Show terms shared in more than 2 cancers
# shared_pos_terms = names(table (unlist(tbf_pos_terms))[table(unlist(tbf_pos_terms)) > 2])
# shared_pos_terms2 = lapply (tbf_pos, function(x) 
#   {
#   y = x[,'Hyper_Adjp_BH', drop=F]
#   rownames(y) = x$name
#   y[shared_pos_terms, ,drop=F]
#   })
# shared_pos_terms2 = do.call (cbind, shared_pos_terms2)
# rownames (shared_pos_terms2) = shared_pos_terms
# shared_pos_terms2[is.na(shared_pos_terms2)] = 1
# colnames (shared_pos_terms2) = names (tbf_pos)
# shared_pos_terms_hm = Heatmap (-log10(shared_pos_terms2), col = viridis::rocket(100),
#   row_names_gp = gpar(fontsize=4), column_names_gp = gpar(fontsize=6))
# pdf (paste0(projdir, 'Plots/DAP_GREAT_pos_shared_terms_heatmap.pdf'), width=3.5, height=4)
# shared_pos_terms_hm
# dev.off()

# # Compare pos and neg go terms
# DAP_terms_l = list(high_normal = unique (unlist(shared_neg_terms)), high_cancer = unique (unlist(shared_pos_terms)))

# DAP_terms_venn = ggvenn(DAP_terms_l,
#   fill_color = viridis::magma (4),
#   stroke_size = 0.5, set_name_size = 4
#   )
# pdf (paste0(projdir,'Plots/DAP_goTerms_cancer_normal_venn.pdf'))
# DAP_terms_venn
# dev.off()

# ###-- Hypergeometric test using Motif x peaks match matrix --###
# archp = addMotifAnnotations (ArchRProj = archp, motifSet = "cisbp", name = "Motif", force=T)

# selected_peaks = GRanges (unique(DAP_df_pos$peak))
# peakSet = getPeakSet (archp)
# motifs_in_peaks = readRDS (paste0(projdir,'Annotations/Motif-In-Peaks-Summary.rds'))

# # read in matrix of peaks x TF matches
# ps = getPeakSet (archp)
# all (rowRanges(motifs_in_peaks[[2]]) == ps) # check rowRanges between motif matches mat and main peakSet is the same
# mm_mat = assay (motifs_in_peaks[[2]])

# DAP_hyperMotif_pos_l = lapply (DAP_l_pos, function (x) {
#   x = GRanges (x$peak)
#   hyperMotif (unique(x), peakSet, mm_mat)
#   })

# pval_threshold = 1e-8
# names(DAP_hyperMotif_pos_l) = names (DAP_list)
# DAP_hyperMotif_pos_shared_df = do.call (cbind, DAP_hyperMotif_pos_l)  
# DAP_hyperMotif_pos_shared_df = DAP_hyperMotif_pos_shared_df[, grep ('padj',colnames (DAP_hyperMotif_pos_shared_df))]
# shared_tf = colSums (apply (DAP_hyperMotif_pos_shared_df, 1, function(x) x < pval_threshold))
# DAP_hyperMotif_pos_shared_df = DAP_hyperMotif_pos_shared_df[shared_tf > 3, ]

# hyper_hm = Heatmap (-log10(DAP_hyperMotif_pos_shared_df + 1e-300), 
#         #column_split = heatmap_kmeans$cluster,
#         clustering_distance_columns = 'euclidean',
#         clustering_distance_rows = 'euclidean',
#         cluster_columns=T,
#         #row_order = hc1$order, 
#         cluster_rows=T,
#         col = viridis::mako(100),
#         row_names_gp = gpar(fontsize = 8),
#         column_names_gp = gpar(fontsize = 11))
# pdf (paste0(projdir,'Plots/DAP_HyperMotif_pos_shared_cancer_normal.pdf'), height=13, width=3.2)
# hyper_hm
# dev.off()

# top_tf = 20
# top_tf_namesL = lapply (DAP_hyperMotif_pos_l, function(x) {
#   x = x[x$padj < pval_threshold,]
#   x = x[order (x$padj),]
#   if (nrow(x) > top_tf) head (rownames(x), top_tf)
#   })
# top_tf_names = unique (unlist(top_tf_namesL))
# DAP_hyperMotif_pos_top_df = do.call (cbind, lapply (DAP_hyperMotif_pos_l, function(x) x[top_tf_names,]))
# DAP_hyperMotif_pos_top_df = DAP_hyperMotif_pos_top_df[, grep ('padj',colnames (DAP_hyperMotif_pos_top_df))]
# rownames (DAP_hyperMotif_pos_top_df) = gsub ('_.*','',rownames (DAP_hyperMotif_pos_top_df)) 
# rownames (DAP_hyperMotif_pos_top_df) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", rownames (DAP_hyperMotif_pos_top_df))
# DAP_hyperMotif_pos_top_df = DAP_hyperMotif_pos_top_df[,colSums (DAP_hyperMotif_pos_top_df) != 0]
# hyper_hm = Heatmap (t(scale(t(scale(-log10(DAP_hyperMotif_pos_top_df + 1e-300))))), 
#         #column_split = heatmap_kmeans$cluster,
#         clustering_distance_columns = 'euclidean',
#         clustering_distance_rows = 'euclidean',
#         cluster_columns=T,
#         #row_order = hc1$order, 
#         cluster_rows=T,
#         col = viridis::mako(100),
#         row_names_gp = gpar(fontsize = 8),
#         column_names_gp = gpar(fontsize = 11),
#         column_title = paste('double-scaled top',top_tf, 'TFs')) 
# pdf (paste0(projdir,'Plots/DAP_HyperMotif_pos_',top_tf,'_cancer_normal.pdf'), height=8, width=3)
# hyper_hm
# dev.off()

# # Make rank plots of top TF from individual cancers (vs normal)
# top_tf = 5
# ggUp = list()
# for (i in seq_along (DAP_hyperMotif_pos_l))
#   {
#   df = DAP_hyperMotif_pos_l[[i]][, 'padj',drop=F]
#   df$tf = rownames (df)
#   df$lpadj = -log10 (df$padj + 1e-300)  
#   df = df[order (-df$lpadj),]
#   df$rank = seq_len(nrow(df))
#   ggUp[[i]] <- ggplot(df, aes(rank, lpadj, color = lpadj)) + 
#     geom_point(size = .3) +
#     ggrepel::geom_text_repel(
#           data = df[seq_len(top_tf), ], aes(x = rank, y = lpadj, label = tf), 
#           size = 2,
#           nudge_x = 20,
#           color = "black"
#     ) + scale_color_gradientn (colors = paletteContinuous(set = "comet")) + ggtitle (paste('Top Tfs ',names(DAP_hyperMotif_pos_l)[i]))
#   }
# pdf (paste0(projdir, 'DAP_HyperMotif_pos_rank_tfs_cancer_normal.pdf'), width=8, height=5)
# wrap_plots (ggUp)
# dev.off()


# # DAP higher in normal
# DAP_hyperMotif_neg_l = lapply (DAP_l_neg, function (x) {
#   x = GRanges (x$peak)
#   hyperMotif (unique(x), peakSet, mm_mat)
#   })

# pval_threshold = 1e-8
# names(DAP_hyperMotif_neg_l) = names (DAP_list)
# DAP_hyperMotif_neg_shared_df = do.call (cbind, DAP_hyperMotif_neg_l)  
# DAP_hyperMotif_neg_shared_df = DAP_hyperMotif_neg_shared_df[, grep ('padj',colnames (DAP_hyperMotif_neg_shared_df))]
# shared_tf = colSums (apply (DAP_hyperMotif_neg_shared_df, 1, function(x) x < pval_threshold))
# DAP_hyperMotif_neg_shared_df = DAP_hyperMotif_neg_shared_df[shared_tf > 3, ]

# hyper_hm = Heatmap (-log10(DAP_hyperMotif_neg_shared_df + 1e-300), 
#         #column_split = heatmap_kmeans$cluster,
#         clustering_distance_columns = 'euclidean',
#         clustering_distance_rows = 'euclidean',
#         cluster_columns=T,
#         #row_order = hc1$order, 
#         cluster_rows=T,
#         col = viridis::mako(100),
#         row_names_gp = gpar(fontsize = 8),
#         column_names_gp = gpar(fontsize = 11))
# pdf (paste0(projdir,'Plots/DAP_HyperMotif_neg_shared_cancer_normal.pdf'), height=13, width=3.2)
# hyper_hm
# dev.off()

# top_tf = 20
# pval_threshold = 1e-2
# top_tf_names = lapply (DAP_hyperMotif_neg_l, function(x) {
#   x = x[x$padj < pval_threshold,]
#   x = x[order (x$padj),]
#   if (nrow(x) > top_tf) head (rownames(x), top_tf)
#   })
# top_tf_names = unique (unlist(top_tf_names))
# DAP_hyperMotif_neg_top_df = do.call (cbind, lapply (DAP_hyperMotif_neg_l, function(x) x[top_tf_names,]))
# DAP_hyperMotif_neg_top_df = DAP_hyperMotif_neg_top_df[, grep ('padj',colnames (DAP_hyperMotif_neg_top_df))]
# rownames (DAP_hyperMotif_neg_top_df) = gsub ('_.*','',rownames (DAP_hyperMotif_neg_top_df)) 
# rownames (DAP_hyperMotif_neg_top_df) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", rownames (DAP_hyperMotif_neg_top_df))
# DAP_hyperMotif_neg_top_df = DAP_hyperMotif_neg_top_df[,colSums (DAP_hyperMotif_neg_top_df) != 0]
# hyper_hm = Heatmap (t(scale(t(scale(-log10(DAP_hyperMotif_neg_top_df + 1e-300))))), 
#         #column_split = heatmap_kmeans$cluster,
#         clustering_distance_columns = 'euclidean',
#         clustering_distance_rows = 'euclidean',
#         cluster_columns=T,
#         #row_order = hc1$order, 
#         cluster_rows=T,
#         col = viridis::mako(100),
#         row_names_gp = gpar(fontsize = 8),
#         column_names_gp = gpar(fontsize = 11),
#         column_title = paste('double-scaled top',top_tf, 'TFs')) 
# pdf (paste0(projdir,'Plots/DAP_HyperMotif_neg_',top_tf,'_cancer_normal.pdf'), height=8, width=3)
# hyper_hm
# dev.off()

# # Make rank plots of top TF from individual cancers (vs normal)
# top_tf = 5
# ggDown = list()
# for (i in seq_along (DAP_hyperMotif_neg_l))
#   {
#   df = DAP_hyperMotif_neg_l[[i]][, 'padj',drop=F]
#   df$tf = rownames (df)
#   df$lpadj = -log10 (df$padj + 1e-300)  
#   df = df[order (-df$lpadj),]
#   df$rank = seq_len(nrow(df))
#   ggDown[[i]] <- ggplot(df, aes(rank, lpadj, color = lpadj)) + 
#     geom_point(size = .3) +
#     ggrepel::geom_text_repel(
#           data = df[seq_len(top_tf), ], aes(x = rank, y = lpadj, label = tf), 
#           size = 2,
#           nudge_x = 20,
#           color = "black"
#     ) + scale_color_gradientn (colors = paletteContinuous(set = "comet")) + ggtitle (paste('Top Tfs ',names(DAP_hyperMotif_neg_l)[i]))
#   }
# pdf (paste0(projdir, 'Plots/DAP_HyperMotif_neg_rank_tfs_cancer_normal_down.pdf'), width=8, height=5)
# wrap_plots (ggDown)
# dev.off()

# # Compare pos and neg go terms
# DAP_hyperMotif_shared_all = list(
#   high_normal = unique (unlist(rownames(DAP_hyperMotif_neg_shared_df))), 
#   high_cancer = unique (unlist(rownames(DAP_hyperMotif_pos_shared_df))))

# DAP_hyperMotif_shared_all_venn = ggvenn(DAP_hyperMotif_shared_all,
#   fill_color = viridis::magma (4),
#   stroke_size = 0.5, set_name_size = 4
#   )
# pdf (paste0(projdir,'Plots/DAP_hyperMotif_shared_all_cancer_normal_venn.pdf'))
# DAP_hyperMotif_shared_all_venn
# dev.off()


# #-- Format bed files for HOMER enrichment analysis of DAP --#
# # map hubs back to peak level
# DAH_epit_df = read.csv (paste0(projdir, 'ttest_malignant_vs_normal_hubs.csv'))
# hub_ttest_sig = DAH_epit_df [DAH_epit_df$padj < 0.01,]
# hubs_id_pos = as.numeric(sapply (DAH_epit_df[DAH_epit_df$log2FC > 0, 'hub'], function(x) unlist (strsplit (x, '_'))[2]))
# hubs_id_neg = as.numeric(sapply (DAH_epit_df[DAH_epit_df$log2FC < 0, 'hub'], function(x) unlist (strsplit (x, '_'))[2]))
# hub_peaks_pos = makeGRangesFromDataFrame (do.call (rbind, hubs_obj[[1]][hubs_id_pos]))
# names (hub_peaks_pos) = NULL
# hub_peaks_neg = makeGRangesFromDataFrame (do.call (rbind, hubs_obj[[1]][hubs_id_neg]))
# names (hub_peaks_neg) = NULL

# ps = getPeakSet (archp)

# homer_dir = '/ahg/regevdata/projects/ICA_Lung/Bruno/HOMER/malignant_vs_normal_mesothelium/'
# system (paste ('mkdir -p',homer_dir))
# hubs_pos_df = ps[unique(queryHits(findOverlaps (ps, hub_peaks_pos)))]
# hubs_pos_df$strand = "+/-"
# hubs_pos_df = hubs_pos_df [, c(1,2,3,5)]
# write.table (hubs_pos_df, paste0(homer_dir,'/tohomer_hubs_malignant.txt'),
# col.names=FALSE, row.names=TRUE, quote=FALSE, sep='\t')

# hubs_neg_df = ps[unique(queryHits(findOverlaps (ps, hub_peaks_neg)))]
# hubs_neg_df$strand = "+/-"
# hubs_neg_df = hubs_neg_df [, c(1,2,3,5)]
# write.table (hubs_neg_df, paste0(homer_dir,'/tohomer_hubs_normal_mesothelium.txt'),
# col.names=FALSE, row.names=TRUE, quote=FALSE, sep='\t')

# # make background for HOMER
# background_df = as.data.frame (ps, row.names=NULL)
# background_df$strand = "+/-"
# background_df = background_df [, c(1,2,3,5)]
# write.table (background_df, paste0(homer_dir,'/background.txt'),
# col.names=FALSE, row.names=TRUE, quote=FALSE, sep='\t')

# # Read in motifs with marge package
# homer_res = list()
# top_motifs = 30
# output_homer_files = list.files(homer_dir)[grep('OUTPUT',list.files(homer_dir))]
# for (i in output_homer_files) {
# homer_res[[i]] = as.data.frame(read_known_results (paste0(homer_dir,i)))
# homer_res[[i]]$cluster = i
# }
# homer_res_top = lapply (homer_res, function(x) x[1:top_motifs,])
# homer_res_df = Reduce (rbind, homer_res)
# homer_res_df$cluster = gsub ('OUTPUTtohomer_','',homer_res_df$cluster)
# homer_res_df$cluster = gsub ('.txt','',homer_res_df$cluster)

# homer_heatmap = matrix (nrow = nrow(homer_res_df), 
#   ncol= length(unique(homer_res_df$cluster)))
# colnames (homer_heatmap) = unique(homer_res_df$cluster)
# rownames (homer_heatmap) = homer_res_df$motif_name
# for (i in unique(homer_res_df$cluster))
# {
#  homer_heatmap[,i] =  homer_res_df$log_p_value [match(
#   homer_res_df$motif_name, homer_res_df[homer_res_df$cluster == i,'motif_name'])]
# }
# homer_heatmap = homer_heatmap[rownames(homer_heatmap) %in% unlist(lapply(homer_res_top, function(x) x$motif_name)),]
# #homer_heatmap = homer_heatmap[!duplicated(rownames(homer_heatmap)),]

# homer_heatmap_p = Heatmap (homer_heatmap, cluster_rows=F, 
#         cluster_columns=T, col = pals_heatmap[[4]],
#         row_names_gp = gpar(fontsize = 4),
#         column_names_gp = gpar(fontsize = 5))
# plotPDF (homer_heatmap_p, name = "homer_motifs_heatmap.pdf",
#         ArchRProj = archp, addDOC = FALSE, width = 3, height = 9)
# pdf (paste0(projdir, 'Plots/homer_motifs_barplot.pdf'))
# par (mar=c(4,6,4,4))
# barplot (homer_res_df$log_p_value[1:top_motifs], horiz = TRUE, xlab = '-log p-value',
#   names.arg = homer_res_df$motif_name[1:top_motifs], las=2, col='black')
# dev.off()



# ### chromVAR analysis
# archp = addBgdPeaks (archp, force= TRUE)
# archp = addMotifAnnotations (ArchRProj = archp, 
#     motifSet = "cisbp", name = "Motif",
#     force=TRUE)
# archp = addDeviationsMatrix (
#   ArchRProj = archp, 
#   peakAnnotation = "Motif",
#   force = TRUE
# )

# archp = saveArchRProject (ArchRProj = archp,  
#     load = TRUE)


# # DAM - Find differential TFs (chromVAR) cancer vs normal meso (normal pleura + cf)
# force = FALSE
# metaGroupName2 = 'celltype2'
# if (file.exists (paste0 (projdir,'DAM_',metaGroupName2,'.rds')) & !force) 
#   {
#   DAM_list = readRDS(paste0 (projdir,'DAM_',metaGroupName2,'.rds'))
#   } else {
#   DAM_list = list() 
#   for (i in unique(archp@cellColData[,metaGroupName2])[grep('cancer',unique(archp@cellColData[,metaGroupName2]))])
#     {
#     DAM_list[[i]] = getMarkerFeatures (
#       ArchRProj = archp, 
#       testMethod = "wilcoxon",
#       useGroups = i,
#       bgdGroups = c("normal_mesothelium",'cf_mesothelium'),
#       binarize = FALSE,
#       useMatrix = "MotifMatrix",
#       groupBy = metaGroupName2,
#       useSeqnames="z")
#     }
#   listnames = names(DAM_list)
#   DAM_list = lapply (seq_along (DAM_list), function(x) 
#     {
#     df = do.call (cbind, (assays(DAM_list[[x]])))
#     colnames(df) = names (assays(DAM_list[[x]]))
#     df$name = rowData (DAM_list[[x]])$name
#     df$name = gsub ('_.*','',df$name) 
#     df$name = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", df$name)
#     df
#     })
#   names (DAM_list) = listnames
#   saveRDS (DAM_list, paste0 (projdir,'DAM_',metaGroupName2,'.rds'))
#   } 

# # Extract DAM per comparison and attach genescore data
# # Get genescore of the TFs
# gsMat = assays (gsSE)[[1]]
# rownames (gsMat) = rowData (gsSE)$name
# DAM_list2 = list()
# for (i in names(DAM_list))
#   {
#   DAM_df = DAM_list[[i]]   
#   gsMat = gsMat[rownames(gsMat) %in% DAM_df$name, ]
#   gsMat_normal = gsMat[, archp$cellNames[as.character(archp@cellColData[,metaGroupName2]) == 'normal_mesothelium']]
#   gsMat_cancer = gsMat[, archp$cellNames[as.character(archp@cellColData[,metaGroupName2]) == i]]
#   gsMat_diff = rowMeans(gsMat_cancer) - rowMeans(gsMat_normal)
  
#   mMat = assays (mSE)[[1]]
#   rownames (mMat) = rowData (mSE)$name
#   mMatDam = mMat[rownames (mMat) %in% unique (DAM_df$name), ]
#   mMatDam_normal = mMatDam[,as.character(archp@cellColData[,metaGroupName2]) == 'normal_mesothelium']
#   mMatDam_cancer = mMatDam[,as.character(archp@cellColData[,metaGroupName2]) == i]
#   mMatDam_normal = mMatDam_normal + abs(min (mMatDam))
#   mMatDam_cancer = mMatDam_cancer + abs(min (mMatDam))
#   mMat_diff = rowMeans(mMatDam_cancer) - rowMeans(mMatDam_normal)
  
  
#   DAM_df$geneScore = gsMat_diff[match (DAM_df$name, names(gsMat_diff)) ]
#   DAM_df$dev = mMat_diff[match (DAM_df$name, names(mMat_diff))]
#   DAM_df$GS_direction = as.factor (ifelse (sign (DAM_df$geneScore) > 0, 'cancer','normal'))
  
#   DAM_list2[[i]] = DAM_df
#   }

# # Make volcano plots of DAM per each comparison
# pValThreshold = 0.00001
# dev_threshold = 0.1
# vol_p = list()
# for (i in names(DAM_list))
#   {
#   DAM_df = DAM_list2[[i]]  
#   jitter <- position_jitter (width = 0.001, height = 0.001)
#   vol_p[[i]] = ggplot(DAM_df) +
#     geom_point(aes(x = dev, y = -log10(FDR), colour = GS_direction, size=abs (geneScore)), alpha = 0.5 , position = jitter) +
#     geom_text_repel(aes(
#       x = dev, 
#       min.segment.length = 0.2,
#       box.padding = 0.2,
#       y = -log10(FDR), 
#       label = ifelse(DAM_df$FDR < pValThreshold & abs(DAM_df$dev) > dev_threshold, DAM_df$name,"")),
#       size=2) +
#     ggtitle ("chromVAR cancer vs normal") +
#     xlab("Mean difference") + 
#     ylab("-log10 adjusted p-value") +
#     theme(legend.position = "right",
#           plot.title = element_text(size = rel(1.5), hjust = 0.5),
#           axis.title = element_text(size = rel(1.25)))  + 
#     scale_color_manual(breaks = c('normal', 'cancer',NA),
#                           values=c("brown", "darkgreen","grey")) + 
#     geom_hline(yintercept=-log10(pValThreshold), linetype="dashed", 
#                   color = "black", size=.2) +
#     geom_vline(xintercept=-dev_threshold, linetype="dashed", 
#                   color = "black", size=.2) +
#     geom_vline(xintercept= dev_threshold, linetype="dashed", 
#                   color = "black", size=.2) 
#   }
# pdf (paste0(projdir, 'Plots/chromVAR_cancer_normal_TFs_volcano.pdf'),5,5)
# print (vol_p)
# dev.off()

# # Make elbow plots showing top TFs per comparison
# ggUp = list()
# for (i in names(DAM_list))
#   {
#   DAM_df = DAM_list2[[i]]
#   DAM_df$FDR_signed = sign(DAM_df$MeanDiff) * (-log10(DAM_df$FDR))
#   DAM_df = DAM_df[order (DAM_df$FDR_signed, decreasing =T),]
#   DAM_df$Rank = 1:nrow(DAM_df)
#   #DAM_df$logFDR = -log10(DAM_df$FDR)   
#   ggUp[[i]] <- ggplot(DAM_df, aes(Rank, FDR_signed, color = FDR_signed)) + 
#     geom_point(size = 1) +
#     ggrepel::geom_label_repel(
#           data = DAM_df[rev(seq_len(10)), ], 
#           aes(x = Rank, y = FDR_signed, label = name), 
#           size = 1.5,
#           nudge_x = 2,
#           max.overlaps = 100,
#           color = "black"
#     ) + theme_bw() + 
#     ylab("-log10(P-adj)") + 
#     xlab("Rank Sorted -log10(FDR)") +
#     scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
#     ggtitle (paste('Cancer',i,'vs normal'))
#     }
# pdf (paste0(projdir, 'Plots/chromVAR_cancer_normal_TFs_elbow_plots.pdf'),4,4)
# print (ggUp)
# dev.off()



# # Select for shared TFs across comparisons
# pValThreshold = 1e-2
# meanDiffThreshold = .5
# shared = 2
# top_tf = 200
# DAM_list = DAM_list[!names(DAM_list) %in% c('p811_C24_cancer','p811_C26_cancer')] # Take only cluster C27 from 811 cancer
# top_tfs_posL = lapply (DAM_list, function(x) {
#     y = x[x$FDR < pValThreshold & x$MeanDiff > meanDiffThreshold,]
#     y = y[order (y$FDR),]
#     head (y$name, top_tf)
#     })
# top_tfs_pos = unlist (top_tfs_posL)
# top_tfs_pos = names (table (top_tfs_pos)[table (top_tfs_pos) > shared])

# top_tfs_negL = lapply (DAM_list, function(x) {
#     y = x[x$FDR < pValThreshold & x$MeanDiff < -meanDiffThreshold,]
#     y = y[order (y$FDR),]
#     head (y$name, top_tf)
#     })
# top_tfs_neg = unlist (top_tfs_negL)
# top_tfs_neg = names (table (top_tfs_neg)[table (top_tfs_neg) > shared])
# top_tfs = unique (append (top_tfs_pos, top_tfs_neg))

# # Make UpSet plot of DAP
# pdf (paste0(projdir, 'Plots/DAM_overlap_cancers_vs_normal_upsetPlot.pdf'), width=10, height=4)
# upset(UpSetR::fromList (top_tfs_posL), nsets= length(top_tfs_posL), order.by = "freq", mainbar.y.label = "Intersection size in MPM high DAM'")
# upset(UpSetR::fromList (top_tfs_negL), nsets= length(top_tfs_negL), order.by = "freq", mainbar.y.label = "Intersection size in MPM low DAM'")
# dev.off()

# # Find top TF only in Sarcomatoid vs normal
# top_tfs_sarc_posL = lapply (DAM_list[names(DAM_list) == 'p786_cancer'], function(x) {
#     y = x[x$FDR < pValThreshold & x$MeanDiff > meanDiffThreshold,]
#     y = y[order (y$FDR),]
#     head (y$name, top_tf)
#     })
# top_tfs_sarc_pos = unlist (top_tfs_sarc_posL)
# top_tfs_sarc_negL = lapply (DAM_list[names(DAM_list) == 'p786_cancer'], function(x) {
#     y = x[x$FDR < pValThreshold & x$MeanDiff < -meanDiffThreshold,]
#     y = y[order (y$FDR),]
#     head (y$name, top_tf)
#     })
# top_tfs_sarc_neg = unlist (top_tfs_sarc_negL)
# top_tfs_sarc = unique(c(top_tfs_sarc_pos, top_tfs_sarc_neg))
# top_tfs_sarc = top_tfs_sarc[!top_tfs_sarc %in% top_tfs]
# write.csv (top_tfs, paste0(projdir, 'DAM_padj_',pValThreshold,'_meanDiff_',meanDiffThreshold,'_shared_',shared,'.csv'))

# ### Find activating and repressing TFs ###
# force = FALSE
# if (!file.exists (paste0(projdir, 'activeTF.rds')) | force)
#   {
#   # Subset archr for clusters to consider for TF - genescore correlations
#   metaGroupName_sub = unique (archp@cellColData[,metaGroupName2])
#   metaGroupName_sub = metaGroupName_sub[grepl ('meso',metaGroupName_sub) | grepl('cancer',metaGroupName_sub)]
#   archp_sub = archp[as.character (archp@cellColData[,metaGroupName2]) %in% metaGroupName_sub,]
  
#   #Get Reduced Dims
#   min_cells = 90
#   knns = knnGen (ArchRProj = archp_sub, 
#     reducedDims = "IterativeLSI", 
#     corCutOff = 0.75,
#     dimsToUse = 1:30, 
#     knnIteration = 500, # number of cells sampled from which k neighbors are searched NOTE: ignored if gourp is specified
#     k = 100, # number of neighbors per knnIteration NOTE: this get overwritten if group is specified
#     overlapCutoff = 0.8,
#     seed = 1,
#     cellsToUse = NULL,
#     group = metaGroupName2,
#     min.cells_in_group = min_cells, # ignored if group is NULL
#     min_knn_cluster = 1 # ignored if group is NULL
#     )
  
#   tf_gs_corL = activeTF (archp_sub, knns, metaGroupName2)
#   saveRDS (tf_gs_corL, paste0(projdir, 'activeTF.rds'))
#   } else {
#   tf_gs_corL = readRDS(paste0(projdir, 'activeTF.rds'))
#   }

# # Plot TF deviation vs genescore and maxDelta  
# actTF_threshold = .3
# tf_gs_cor = tf_gs_corL[[1]]
# tf_gs_cor$label = ""
# tf_gs_cor$label[abs(tf_gs_cor$tf_gs_cor) > actTF_threshold] = rownames (tf_gs_cor)[abs(tf_gs_cor$tf_gs_cor) > actTF_threshold] 
# tf_gs_cor$type = ""
# tf_gs_cor$type [tf_gs_cor$tf_gs_cor > actTF_threshold] = 'activator'
# tf_gs_cor$type [tf_gs_cor$tf_gs_cor < -actTF_threshold] = 'repressor'
# tf_gs_cor = tf_gs_cor[order(-tf_gs_cor$tf_gs_cor),]

# p = ggplot (tf_gs_cor, aes(tf_gs_cor, delta, color = type, label = label)) +
#   geom_point(alpha=.5) + 
#   theme_ArchR() +
#   geom_vline(xintercept = 0, lty = "dashed") + 
#   scale_color_manual (values = c("activator"="navy", "repressor"="firebrick3")) +
#   xlab("Correlation To Gene Score") +
#   ylab("Max TF Motif Delta") + 
#   geom_label_repel(max.overlaps = 20) + 
#   scale_y_continuous(
#     expand = c(0,0), 
#     limits = c(0, max(tf_gs_cor$delta)*1.05)
#   ) + geom_vline (xintercept=actTF_threshold, linetype="dashed", color = "blue") +
#    geom_vline (xintercept=actTF_threshold, linetype="dashed", color = "blue") +
#    geom_vline (xintercept=-actTF_threshold, linetype="dashed", color = "blue")
# pdf (paste0(projdir,'Plots/TF_activation_ArchR.pdf'), 7,7)
# print (p)
# dev.off()

# # 
# activators = rownames(tf_gs_cor)[tf_gs_cor$tf_gs_cor > actTF_threshold]
  
# # Save Df of TFs
# gsAvknn_df = as.data.frame (tf_gs_corL[[3]])
# gsAvknn_df$tf  = rownames (gsAvknn_df)
# gsAvknn_dfl = gather (gsAvknn_df, knn_gs, mean_gs, 1:(ncol(gsAvknn_df)-1))
# mAvknn_df = as.data.frame (tf_gs_corL[[2]])
# mAvknn_df$tf  = rownames (mAvknn_df)
# mAvknn_dfl = gather (mAvknn_df, knn_gs, mean_deviation, 1:(ncol(mAvknn_df)-1))
# gsmAvknn_df = cbind (gsAvknn_dfl, mAvknn_dfl)
# gsmAvknn_df = gsmAvknn_df[,-2]
# gsmAvknn_df$cluster = gsub ('KNN\\d*','',gsmAvknn_df$knn_gs)
# gsmAvknn_df$cor = tf_gs_cor$tf_gs_cor[match (gsmAvknn_df$tf, rownames(tf_gs_cor))]
# activator_p = list()
# for (activator in activators[activators %in% top_tfs])
#   {
#   activator.df = gsmAvknn_df[gsmAvknn_df$tf == activator,]
#   activator_p[[activator]]= ggplot(activator.df, aes(
#     x=mean_deviation, y= mean_gs)) + 
#   geom_point(size = 2, shape = 21, alpha=0.8,
#     aes(fill=cluster, colour = "white")) + 
#   ggtitle (paste('motif vs Gene score',activator), subtitle = round(gsmAvknn_df$cor[gsmAvknn_df$tf == activator],2)) +
#     xlab('Motif deviation') + ylab('Gene Score') +
#     scale_fill_manual(values= unname(paletteDiscrete(unique(gsmAvknn_df$cluster)))) +
#     theme_bw()
#   }

# pdf (paste0(projdir,'Plots/Activators_KNN_TF_scatterplot.pdf'), width=20, height=50)
# print (wrap_plots (activator_p, ncol=4))
# dev.off()

# # Take intersection of DAM and activators  
# top_active_tfs = top_tfs[top_tfs %in% activators]
# top_active_tfs_sarc = top_tfs_sarc[top_tfs_sarc %in% activators]

# # make heatmap of TF differential in cancer vs normal across all cancers vs normal
# # Get MotifMatrix and subset for DAM+activators
# metaGroupName2 = 'celltype'
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMatDam = mMat[unique(c(top_active_tfs, top_active_tfs_sarc)), ]
# mMatDam = as.data.frame (t(mMatDam))
# mMatDam$metaGroup = as.character (archp@cellColData[,metaGroupName2])
# mMatDam = aggregate (.~ metaGroup, mMatDam, mean)
# rownames (mMatDam) = mMatDam[,1]
# mMatDam = mMatDam[,-1]
# mMatDam = t(mMatDam)
# mMatDam = mMatDam[,grepl('meso',colnames(mMatDam)) | grepl('cancer',colnames(mMatDam))] # Subset for only cancer and normal mesothelium
# mMatDam = mMatDam[,-grep('p786_mesothelium',colnames(mMatDam))]
# rownames (mMatDam) = unique(c(top_active_tfs, top_active_tfs_sarc)) # subset for only significant TFs

# # Add geneScore of TFs
# gsMat =  assays (gsSE)[[1]]
# rownames (gsMat) = rowData (gsSE)$name
# gsMat = gsMat[rownames(gsMat) %in% rownames(mMatDam),]
# gsMat = as.data.frame (t (gsMat))
# gsMat$cluster = as.factor(getCellColData (archp)[,metaGroupName2])
# gsMat = aggregate (formula= .~ cluster, data= gsMat, FUN=mean)
# rownames (gsMat) = gsMat[,1]
# gsMat = gsMat[,-1]
# gsMat = as.data.frame (t(gsMat))
# gsMat = gsMat[rownames (mMatDam),]
# gsMat = gsMat[,colnames(mMatDam)]

# # Add gene integration of TFs
# giMat =  assays (giSE)[[1]]
# rownames (giMat) = rowData (giSE)$name
# giMat = giMat[rownames(giMat) %in% rownames(mMatDam),]
# giMat = as.data.frame (t (giMat))
# giMat$cluster = as.factor(getCellColData (archp)[,metaGroupName2])
# giMat = aggregate (formula= .~ cluster, data= giMat, FUN=mean)
# rownames (giMat) = giMat[,1]
# giMat = giMat[,-1]
# giMat = as.data.frame (t(giMat))
# giMat = giMat[rownames (mMatDam),]
# giMat = giMat[,colnames(mMatDam)]

# # Add deviations from TCGA ATAC bulk
# tcga_dev = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/TCGA_atac/tcga_meso_atac_deviations_clinical_info.rds')
# tcga_mat = assays (tcga_dev)$z
# tcga_mat = aggregate (t(tcga_mat), by = list(sampleID= colnames(tcga_mat)), mean)
# rownames (tcga_mat) = tcga_mat[,1]
# tcga_mat = tcga_mat[,-1]
# tcga_mat = as.data.frame (t(tcga_mat))
# rownames (tcga_mat) = unname (sapply (rownames(tcga_mat), function(x) unlist (strsplit (x, '_'))[3]))
# #tcga_mat = rowMeans (tcga_mat)
# tcga_mat = tcga_mat[rownames(mMatDam),]
# colnames (tcga_mat) = unname (sapply (colnames(tcga_mat), function(x) unlist (strsplit (x, '_'))[2]))

# # Add scRNA 
# rna_mat = gcdata@assays$RNA@data
# rna_ct = as.character (gcdata$celltype_2)
# rna_ct[rna_ct == 'Malig'] = paste0 (rna_ct[rna_ct == 'Malig'],'_', gcdata$sampleid[gcdata$celltype_2 == 'Malig']) 
# rna_mat = rna_mat[rownames(rna_mat) %in% top_active_tfs,]
# rna_mat = as.data.frame (t(as.matrix(rna_mat)))
# rna_mat$cluster = rna_ct
# rna_mat = aggregate (formula= .~ cluster, data= rna_mat, FUN= mean)
# rownames (rna_mat) = rna_mat$cluster
# rna_mat = rna_mat[,-1]
# rna_mat2 = apply (rna_mat, 1, function(x) x[match(rownames(mMatDam), colnames(rna_mat))])
# rna_mat2 = rna_mat2[,grepl('Malig', colnames (rna_mat2)) | grepl('Sarco', colnames (rna_mat2))]

# # Add Viestra et al. motif clusters as heatmap annotation
# mclust = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/DBs/Viestra_motif_clusters/viestra_motif_clusters.csv')
# mclust$tf = sapply (mclust$Motif, function(x) unlist(strsplit(x, '_'))[1])
# mclust = mclust[!duplicated(mclust$tf),]
# mclustS = paste0('cl',(mclust[match(rownames (mMatDam), mclust$tf) ,c(1)]))

# mclust_ha = HeatmapAnnotation (Viestra_et_al.motif_cluster = anno_simple (mclustS, pch = mclustS, 
#   col=setNames(rep('white',length(unique(mclustS))),unique(mclustS)) ,
# pt_size = unit(.15, "cm")),  annotation_name_gp = gpar(fontsize=7),
# which='row', show_annotation_name =T, show_legend=F, simple_anno_size = unit(.1, "cm"), gap=unit(.1, "cm"))

# tf_name_ha = rowAnnotation (max = anno_barplot(rowMax(mMatDam),location = 2,axis_param = list(
#         facing = "outside"),baseline="min"),TF = anno_text(rownames(mMatDam), location = 1, rot = 30, 
#     just = "right", gp = gpar(fontsize = 6)))

# km_cols = kmeans(t(mMatDam), 2)
# heat_motifs = Heatmap (t(scale(t(mMatDam))),
#         left_annotation = tf_name_ha, 
#         clustering_distance_columns = 'pearson',
#         clustering_distance_rows = 'pearson',
#         cluster_columns=F,
#         cluster_rows=T,
#         name = 'chromVAR', row_km=4,
#         #col = pals_heatmap[[2]],
#         col = viridis::magma (10), column_split = km_cols$cluster,
#         row_names_gp = gpar(fontsize = 4),
#         column_names_gp = gpar(fontsize = 6),
#         #rect_gp = gpar(col = "white", lwd = .5),
#         border=TRUE, row_dend_side  = 'left')
# heat_gs = Heatmap (t(scale(t(gsMat))), 
#         column_split = km_cols$cluster,
#         #clustering_distance_columns = 'Pearson', 
#         cluster_columns=F,#col = pals_heatmap[[1]],
#         row_names_gp = gpar(fontsize = 5),
#         column_names_gp = gpar(fontsize = 6),
#         name = 'GeneScore',
#         #col = pals_heatmap[[1]],
#         col = viridis::viridis(10),
#         #rect_gp = gpar(col = "white", lwd = .5),
#         border=TRUE
#         #right_annotation = motif_ha
#         )

# #heat_gi = Heatmap (t(scale(t(giMat))), column_order = heat_motifs@column_order,
# #        #clustering_distance_columns = 'Pearson', 
# #        cluster_columns=F,#col = pals_heatmap[[1]],
# #        row_names_gp = gpar(fontsize = 5),
# #        column_names_gp = gpar(fontsize = 6),
# #        name = 'GeneIntegration',
# #        #col = pals_heatmap[[1]],
# #        col = viridis::inferno(10),
# #        #rect_gp = gpar(col = "white", lwd = .5),
# #        border=TRUE
# #        #right_annotation = motif_ha
# #        )
# heat_bulkATAC = Heatmap (t(scale(t(tcga_mat))), 
#         clustering_distance_columns = 'pearson', 
#         cluster_columns=T, 
#         #col = pals_heatmap[[2]],
#         row_names_gp = gpar(fontsize = 5), 
#         name = 'bulk_ATAC',
#         column_names_gp = gpar(fontsize = 6),
#         #rect_gp = gpar(col = "white", lwd = .5),
#         col = viridis::mako(10),
#         border=TRUE)#, left_annotation=bulk_box_ha)
# heat_scRNA = Heatmap (t(scale(t(rna_mat2))), 
#         clustering_distance_columns = 'pearson',
#         right_annotation = mclust_ha,
#         cluster_columns=T, 
#         #col = pals_heatmap[[2]],
#         row_names_gp = gpar (fontsize = 0), 
#         name = 'scRNA',
#         column_names_gp = gpar (fontsize = 6),
#         #rect_gp = gpar(col = "white", lwd = .5),
#         col = viridis::rocket(10),
#         border=TRUE)#, left_annotation=bulk_box_ha)

# pdf (paste0(projdir, 'Plots/ChromVAR_noDoublets_motifs_gene_score_heatmaps_filtered_activators_cancer_vs_normal.pdf'), width=6, height=15)
# draw (heat_motifs + heat_bulkATAC + heat_gs + heat_scRNA , padding = unit(c(40, 2, 2, 2), "mm"),
#   column_title = paste("Diff Motifs (chromVAR) FDR <",pValThreshold,'MeanDiff >', meanDiffThreshold,'shared', shared,'top',top_tf,'activTreshold',actTF_threshold), column_title_gp = gpar(fontsize = 10))
# dev.off()


# ####-- Coexpression of TFs --####
# metaGroupName2 = 'celltype2'
# select.groups = c('p848_cancer','p811_C26_cancer','p811_C24_cancer','p811_C27_cancer',
#   'p846_cancer','p786_cancer','p826_cancer')

# actTF_threshold = .3
# activators = rownames(tf_gs_cor)[tf_gs_cor$tf_gs_cor > actTF_threshold]

# mMat = assays (mSE)$deviations
# mMat = mMat[,archp$cellNames[as.character(archp@cellColData[,metaGroupName2]) %in% select.groups]]
# rownames(mMat) = rowData(mSE)$name
# mMat = mMat[activators,] # subset only for functional TFs

# # plot variance of TFs and select top variable TFs
# mMat = as.data.frame (t(mMat))
# mMat$metaGroup = as.character(archp@cellColData[,metaGroupName2])[match (rownames(mMat), rownames(archp))] 
# mMat = aggregate (.~metaGroup, mMat, mean)
# rownames (mMat) = mMat[,1]
# mMat = mMat [,-1]

# TF_cluster_hm = Heatmap (mMat,
#         cluster_rows = T,# km = 5,
#         clustering_distance_rows = 'euclidean',
#         clustering_distance_columns = 'euclidean',
#         cluster_columns=T, col = viridis::viridis (100),
#         row_names_gp = gpar(fontsize = 5),
#         column_names_gp = gpar(fontsize = 5))

# mMat = assays (mSE)$deviations
# mMat = mMat[,archp$cellNames[as.character(archp@cellColData[,metaGroupName2]) %in% select.groups]]
# rownames(mMat) = rowData(mSE)$name
# mMat = mMat[activators,]
# cor_TF = cor (as.matrix(t(mMat)))
# cor_TF_hm = Heatmap (cor_TF,
#         cluster_rows = T,# km = 5,
#         clustering_distance_rows = 'pearson',
#         clustering_distance_columns = 'pearson',
#         cluster_columns=T, col = viridis::rocket (100),
#         row_names_gp = gpar(fontsize = 5),
#         column_names_gp = gpar(fontsize = 5))

# pdf (paste0 (projdir, 'TF_cancer_modules_heatmaps.pdf'), 10,10)
# TF_cluster_hm
# cor_TF_hm
# dev.off()

# d = as.dist (1-cor(t(mMat)))
# # Hierarchical clustering using Complete Linkage
# hc1 <- hclust(d, method = "complete" )

# cor_TF = cor (t(mMat[hc1$order,]))

# km = kmeans (cor_TF, centers=5)
# set.seed (123)
# #pal_corr1 = colorRamp2 (c(-1, 0, 1), c("hotpink", "black", "cyan"))
# #ha = HeatmapAnnotation (cluster= as.character(km$cluster), which='row')
# #ha1 = HeatmapAnnotation (cluster= as.character(km$cluster), which='column')
# cor_TF_hm = Heatmap (cor_TF,
#         cluster_rows = F,
#         cluster_columns=F,
#         row_split = km$cluster,
#         column_title = 'all_MPM',
#         column_km_repeats=20,
#         col = viridis::viridis(100),
#         #clustering_distance_rows = 'pearson', 
#         #clustering_distance_columns = 'pearson', 
#         column_split = km$cluster,
#         #col = pal_corr1,
#         row_names_gp = gpar(fontsize = 12),
#         column_names_gp = gpar(fontsize = 12))
# png (paste0(projdir,'Plots/TF_coexp_activators_cancer_heatmap.png'),
#    width=2500,height=2500, res=80)
# cor_TF_hm
# dev.off()

#  # make oncogenic signatures and plot on UMAP
# malig_modules = lapply (row_order (cor_TF_hm), function (x) rownames (cor_TF)[x])
# malig_modules_dev = t(Reduce(rbind,lapply (malig_modules, function(x) colMeans(as.data.frame(mMat)[x,]))))
# colnames(malig_modules_dev) = paste0 ("module_",1:ncol(malig_modules_dev)) 

# for (i in 1:ncol(malig_modules_dev))
# {
#   archp = addCellColData (archp, name = colnames(malig_modules_dev)[i], 
#     data = malig_modules_dev[,i], cells =rownames(malig_modules_dev), force=TRUE)
# }

# malig_mod_p = lapply (colnames (malig_modules_dev),
#   function(x) plotEmbedding (ArchRProj = archp[archp$cellNames %in% rownames (malig_modules_dev)],
#    colorBy = "cellColData",
#  name = x, embedding = "UMAP"))

# png (paste0(projdir,'Plots/TF_modules_UMAPs.png'),
#    width=2200,height=2200, res=300)
# wrap_plots(malig_mod_p)
# dev.off()

 ###--- Hubs analysis --###
metaGroupName = "celltype"
maxDist = 12500
dgs=3000
cor_cutoff=0.4
min_peaks=3
force = T
projdir_hubs = paste0 (projdir, 'hubs_obj_cor_',cor_cutoff,'_md_',maxDist,'_dgs_',dgs,'_min_peaks_',min_peaks,'/')
dir.create (paste0(projdir_hubs,'Plots/'), recursive = T)
archp <- addCoAccessibility(
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    corCutOff = 0.75,
    maxDist = maxDist, #stitching distance for super-enhancers (ROSE)
    overlapCutoff = 0.8
    )

if (!file.exists (paste0(projdir_hubs,'hubs_obj.rds')) | force)
  {
  hubs_obj = hubs_finder (
    ArchRProj = archp, 
    group_by = NULL,
    cor_cutoff = cor_cutoff,
    #select_group = metaGroup_df$metaGroup,
    cor_FDR = 1, 
    cor_var = 0,
    min_peaks = min_peaks,
    macs_score = 1,
    dgs = dgs,
    cores=1
    ) 
  # Export bed file of hub regions
  hubs_bed = as.data.frame (hubs_obj$hubsCollapsed)
  hubs_bed['width'] = hubs_obj$hubs_id
  write.table (hubs_bed, paste0(projdir_hubs, 'hub_regions.bed'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  saveRDS (hubs_obj, paste0(projdir_hubs,'hubs_obj.rds'))
  } else {
  hubs_obj = readRDS (paste0(projdir_hubs,'hubs_obj.rds'))  
  }

################################################################
# Check enrichment of motifs in hubs vs rest of peaks in peakset
################################################################

metaGroupName = 'celltype2'
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
    reproducibility = "2",
    maxPeaks = 500000, 
    minCells=20,
    force =TRUE) # I think this should be set corresponding to the smallest cluster in the group or lower
archp = addPeakMatrix (archp)
#archp = saveArchRProject (archp, load=TRUE)

###--- Hubs analysis --###
metaGroupName = "celltype2"
maxDist = 12500
dgs=3000
cor_cutoff=0.4
min_peaks=5
force = T
projdir_hubs = paste0 (projdir, 'hubs_obj_cor_',cor_cutoff,'_md_',maxDist,'_dgs_',dgs,'_min_peaks_',min_peaks,'/')
dir.create (paste0(projdir_hubs,'Plots/'), recursive = T)
archp <- addCoAccessibility(
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    corCutOff = 0.75,
    maxDist = maxDist, #stitching distance for super-enhancers (ROSE)
    overlapCutoff = 0.8
    )
force=T
if (!file.exists (paste0(projdir_hubs,'hubs_obj.rds')) | force)
  {
  hubs_obj = hubs_finder (
    ArchRProj = archp, 
    group_by = NULL,
    cor_cutoff = cor_cutoff,
    #select_group = metaGroup_df$metaGroup,
    cor_FDR = 1, 
    cor_var = 0,
    min_peaks = min_peaks,
    macs_score = 1,
    dgs = dgs,
    cores=1
    ) 
  # Export bed file of hub regions
  hubs_bed = as.data.frame (hubs_obj$hubsCollapsed)
  hubs_bed['width'] = hubs_obj$hubs_id
  write.table (hubs_bed, paste0(projdir_hubs, 'hub_regions.bed'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  saveRDS (hubs_obj, paste0(projdir_hubs,'hubs_obj.rds'))
  } else {
  hubs_obj = readRDS (paste0(projdir_hubs,'hubs_obj.rds'))  
  }


tf_match = getMatches (archp, name = 'Motif')
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
colnames (tf_match) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames(tf_match))
#tf_match = tf_match[,top_all_active_tfs]

hyper_res = list()
for (gr in unique(archp@cellColData[,metaGroupName]))
  {
  bg_ps = readRDS (paste0(projdir,'PeakCalls/',gr,'-reproduciblePeaks.gr.rds'))
  idx_hubs = unlist(lapply (hubs_obj$hubsClusters[[1]], function(x) {
  hubgr = makeGRangesFromDataFrame(x) 
  hub_idx = queryHits (findOverlaps (bg_ps, hubgr))
  if (length(hub_idx) == length (hubgr)) hub_idx else NA
  }))
  idx_hubs = na.omit(idx_hubs)
  ps_hubs = bg_ps[idx_hubs]
  #ps_hubs = bg_ps[sample (length(bg_ps), length(ps_hubs))]
  tf_match2 = tf_match[queryHits (findOverlaps (tf_match, bg_ps))]
  hyper_res[[gr]] = hyperMotif (selected_peaks = ps_hubs, bg_peak_motif_RSE = tf_match2)  
  saveRDS (hyper_res, paste0(projdir_hubs, 'Enrichment_test_on_hubs.rds'))
  }

#hyper_res_flt = lapply (hyper_res, function(x) rownames(head(x[order(x$padj),],100))) # show top 50 TF per cell type
hyper_res_df = do.call (cbind, hyper_res)
#hyper_res_df = hyper_res_df[unique(unlist(hyper_res_flt)),]
hyper_res_df = hyper_res_df[, grep('padj', colnames(hyper_res_df))]
hyper_res_df = hyper_res_df[apply (apply (hyper_res_df, 2, function(x) x < 0.01), 1, any),]
hyper_res_df = -log10(hyper_res_df)
rown = rownames (hyper_res_df)
#hyper_res_df[hyper_res_df > 300] = 300
hyper_res_df = do.call(data.frame,lapply(hyper_res_df, function(x) replace(x, is.infinite(x),-300)))
hyper_res_df[hyper_res_df == -300] = max (hyper_res_df)
rownames (hyper_res_df) = rown
hyper_res_df = hyper_res_df[rownames(hyper_res_df) != 'ENSG00000250542',]
#hyper_res_df = scale(hyper_res_df)
#hyper_res_df[is.na (hyper_res_df)] = 0
hyper_res_df = hyper_res_df[,colnames(hyper_res_df) != 'not_assigned.padj']
hyper_ht = Heatmap (t(scale(t(hyper_res_df))), col = viridis(100),  row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5))
pdf (paste0(projdir, 'Plots/hyper_hubs_vs_rest_per_celltype_heatmap.pdf'), height=45.5, width = 3)
hyper_ht
dev.off()

#ps = getPeakSet (archp)


  #tf_match = tf_match[queryHits(findOverlaps (tf_match, ps_hubs))]


write.csv (hyper_res, paste0(projdir_hubs,'Hubs_TF_enrichment.csv'))

#flt_rows = apply (t(hyper_res_df), 1, function(x) x < 0.05)
#flt_rows = rowSums (flt_rows) > 0
#hyper_res_df = hyper_res_df[flt_rows,]
colnames (hyper_res_df) = names (hyper_res)
hyper_res_df = -log10(hyper_res_df)
is.na(hyper_res_df) <- do.call(cbind,lapply(hyper_res_df, is.infinite))
hyper_res_df[is.na (hyper_res_df)] = -1
hyper_res_df[hyper_res_df == -1] = max (hyper_res_df)

# Check specific TF motifs in hubs 
tf_match = getMatches (archp, name = 'cisbp_Motif')
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
colnames(tf_match) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames(tf_match))
motif_regions = assays (tf_match)[[1]][, 'TCF3']
motif_regions = which (motif_regions)
hub_ov = hubs_obj$hubsCollapsed[unique (queryHits (findOverlaps (hubs_obj$hubsCollapsed, tf_match[motif_regions,])))]

ps = getPeakSet (archp)
ps_ov = ps[unique (queryHits (findOverlaps (ps, tf_match[motif_regions,])))]


# ### Compute similarity of a genomic region across clusters around a gene ###
# maxDist = 500000
# archp = addPeak2GeneLinks (
#   ArchRProj = archp,
#   maxDist = maxDist,
#   reducedDims = "IterativeLSI",
#   overlapCutoff = 0.4
#   )
# p2g_corr = 0.3
# p2g = getPeak2GeneLinks (
#     ArchRProj = archp,
#     corCutOff = p2g_corr,
#     resolution = 1,
#     returnLoops = FALSE
# )

# # Convert df in Granges add gene Name and correlation
# gene = 'GATA4'
# p2g$geneName = mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
# p2g = p2g[!is.na (p2g$FDR),] # remove NaN correlations (not sure why there are some)
# p2g = p2g[p2g$Correlation > p2g_corr, ]
# p2gGR = metadata (p2g)$peakSet[p2g$idxATAC]
# p2gGR$geneName = p2g$geneName
# p2gGR$correlation = p2g$Correlation
# p2gGR_gene = p2gGR[p2gGR$geneName == gene]
# # Subset peakMatrix by peaks in genomic region 
# if(!any (ls() == 'pSE')) pSE = getMatrixFromProject (archp, useMatrix = 'PeakMatrix')
# pSE = pSE[, rownames(archp)]
# pSEs = subsetByOverlaps (pSE, p2gGR_gene)

# # Annotate peaks with peakSet
# metaGroupName2 = 'celltype'
# ps = getPeakSet (archp)
# ann_peaks = ps[queryHits(findOverlaps (ps, pSEs))]
# peakType = paste(ann_peaks$peakType, paste(ranges(ann_peaks))) 
# pmat = t(assays (pSEs)[[1]])
# all (rownames(archp) == rownames (pmat))
# pmat = as.data.frame (pmat)
# pmat = pmat / archp$ReadsInTSS * 1000
# pmat$metaGroup = as.character(archp@cellColData[,metaGroupName2])
# pmat = aggregate (.~ metaGroup, pmat, sum)
# cell_norm = sapply (pmat$metaGroup, function(x) sum(archp@cellColData[,metaGroupName2] == x))
# pmat = pmat[,-1]
# pmat = pmat / cell_norm
# colnames (pmat) = peakType
# rownames (pmat) = names(cell_norm)
# ha = HeatmapAnnotation (peakType = sapply(peakType, function(x) unlist(strsplit(x,' '))[1]),simple_anno_size = unit(.2, "cm"))
# pdf (paste0(projdir,'Plots/Peaks_to_',gene,'.pdf'), height=2, width=6)
# Heatmap (pmat, column_title = paste (gene, 'local epigenome -',maxDist),
#   clustering_distance_rows = "euclidean", name = 'coverage', bottom_annotation = ha,
#   cluster_columns=F, row_names_gp = gpar (fontsize = 5),
#   cluster_rows = T, column_names_gp = gpar (fontsize = 3), border=T,
#   col = rev(viridis::magma(100)))
# dev.off()

# LE_region = p2gGR_gene[c(1, length(p2gGR_gene))]
# LE_region = GRanges (seqnames = seqnames (LE_region)[1], ranges = IRanges(start = start(LE_region)[1], end=end(LE_region)[2]))
# LE_region_p <- plotBrowserTrack(
#     ArchRProj = archp, 
#     groupBy = metaGroupName2, 
#     region = LE_region, 
#     upstream = 0,
#     downstream = 0,
#     pal = paletteDiscrete(unique(archp@cellColData[,metaGroupName2]), set='rushmore', reverse=T),
#     loops = getPeak2GeneLinks (archp, corCutOff = 0.2,
#       returnLoops = TRUE)[[1]],
#     useGroups= NULL
# )
# plotPDF (LE_region_p, ArchRProj = archp, width=10, height=4)

# ann_peaks = ps[queryHits(findOverlaps (ps, pSEs))]
# peakType = paste(ann_peaks$peakType, paste(ranges(ann_peaks))) 
# pmat = t(assays (pSEs)[[1]])
# all (rownames(archp) == rownames (pmat))
# pmat = as.data.frame (pmat)
# pmat = pmat / archp$ReadsInTSS * 1000

# LE_score = data.frame (LE_score = rowSums (pmat))
# LE_score$metaGroup = as.character(archp@cellColData[, metaGroupName2])[match(rownames(LE_score), archp$cellNames)]
# #LE_score = LE_score[order(-LE_score$LE_score),,drop=F]
# #LE_score$cluster = rownames(LE_score)
# #LE_score$cluster = factor (LE_score$cluster, levels = LE_score$cluster)
# p<-ggplot(data=LE_score, aes(x=metaGroup, y=LE_score)) +
#   geom_boxplot() + theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# pdf (paste0(projdir,'Plots/local_epigenome_score_',gene,'_barplot.pdf'), 3,5)
# p
# dev.off()

# # Compare to genescore
# p1 <- plotGroups(
#     ArchRProj = archp, 
#     groupBy = metaGroupName2, 
#     colorBy = "GeneScoreMatrix", 
#     name = gene,
#     plotAs = "violin",
#     pal = paletteDiscrete(unique(archp@cellColData[,metaGroupName2]), set='rushmore', reverse=T)
#    )
# pdf (paste0(projdir,'Plots/genescore_',gene,'.pdf'),2.5,3)
# p1
# dev.off()

# ### Subset only mesothelioma and normal mesothelium ###
# projdir_meso = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/all_meso_cf_mesoSub_ATAC/'
# system (paste('mkdir -p', paste0(projdir_meso,'Plots/')))
# archp_meso = subsetArchRProject(
#   ArchRProj = archp,
#   cells = rownames(archp)[as.character(archp@cellColData[,metaGroupName2]) %in% c(
#     'cf_mesothelium','normal_mesothelium','p848_cancer','p826_cancer',
#     'p786_cancer','p811_cancer','p846_cancer','p786_mesothelium')],
#   outputDirectory = projdir_meso,
#   dropCells = TRUE,
#   logFile = NULL,
#   threads = getArchRThreads(),
#   force = TRUE
# )







