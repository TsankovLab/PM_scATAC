conda activate meso_scatac
#use UGER
R

####### ANALYSIS of TUMOR compartment #######
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
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/stroma/scatac_ArchR'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)


#devtools::install_github("immunogenomics/presto") #needed for DAA
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

set.seed (1234)
addArchRThreads (threads = 1) 
addArchRGenome ("Hg38")


# Load RNA
srt = readRDS ('../scrna/srt.rds')
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
      
  archp = ArchRProject (
    ArrowFiles = ArrowFiles, 
    outputDirectory = projdir,
    copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  )
  
  ### Subset ArchR object only for cells retained in Signac analysis ####
  archp_main = loadArchRProject ('../../main/scatac_ArchR/')
  archp$celltype = 0
  archp$celltype[match (rownames(archp_main@cellColData), rownames(archp@cellColData))] = archp_main$celltype
  
  # Use Annotation from signac ####
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
  
  archp = archp[!archp$celltype %in% c('AT1','AT2','Alveolar','B_cells','Ciliated','Mast','MonoMac','NK','Plasma','Secretory','T.NK.cells','T_cells','cDCs','pDCs','Malignant')]

  archp$Sample2 = archp$Sample
#  archp$Sample2[grep ('RPL',archp$Sample2)] = 'normal_pleura'

  ## and replace with Mo's annotation ##
  # normal_mo = read.csv ('../../tumor_compartment/all_tissues_ArchR/metadata_external_data/tsankov_refined_annotation.csv')



  
  #archp = archp[!is.na(archp$celltype)]

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
    groupBy = "Sample2", force=TRUE
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
  
  # Check LEC markers to reannotate cluster
  genes = c('TFF3')
  archp = addImputeWeights (archp)

  pdf (file.path ('Plots',paste0('markers_LEC_umap.pdf')), width = 30,height=14)
  TF_p = plotEmbedding(
      ArchRProj = archp, 
      colorBy = "GeneScoreMatrix", 
      name = genes, 
      embedding = "UMAP_H",
      pal = palette_expression,
      imputeWeights = getImputeWeights(archp)
  )
  TF_p
  dev.off()
  archp$celltype[archp$Clusters == 'C7'] = 'LEC'

  pdf (file.path('Plots','celltype_umap_signac_filtered.pdf'),5,5)
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

  print (umap_p1)
  print (umap_p2)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()
  

  archp = saveArchRProject (archp, load = T, dropCells=T)
  
  } else {
  archp = loadArchRProject (projdir)
  }


### Call peaks on celltypes ####
metaGroupName = 'Clusters_H'
force=TRUE

pdf() # This is necessary cause cairo throws error and stops the script
if(!all(file.exists(file.path('PeakCalls', unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds')))) | force) 
source (file.path('..','..','git_repo','utils','callPeaks.R'))
dev.off()


### chromVAR analysis ####
force=TRUE
if (!all(file.exists(file.path('Annotations',
  c('Motif-Matches-In-Peaks.rds',
    'Motif-Positions-In-Peaks.rds',
    'Motif-In-Peaks-Summary.rds')))) | force)
source (file.path ('..','..','git_repo','utils','chromVAR.R'))
  
# # Find activating and repressing TFs #### 
# if (!file.exists ('TF_activators_genescore.rds')) 
#   {
#     source (file.path('..','..','git_repo','utils','activeTFs.R'))
#   } else {
#     corGSM_MM = readRDS ('TF_activators_genescore.rds') 
#   }


# Differential Accessed motifs ####
metaGroupName = "celltype2"
force=FALSE
source (file.path('..','..','git_repo','utils','DAM.R'))

mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat_mg = mMat[active_DAM, ]
mMat_mg = as.data.frame (t(mMat_mg))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
mMat_mg = mMat_mg[names (DAM_list),]

# Generate RNA pseudobulk of matching cell types ####
metaGroupName = 'celltype2'
#selected_TF = c(rownames(DAM_hm@matrix), 'NR4A3','NR4A2','NR4A1')
ps = log2(as.data.frame (AverageExpression (srt, features = active_DAM, group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
ps = ps[, colnames(DAM_hm@matrix)]
ps_tf = ps[active_DAM,]

  
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
          col = rev(palette_deviation),
          width = unit(2, "cm")
          #right_annotation = motif_ha
          )

scaled_ps = t(scale(t(ps_tf)))
scaled_ps[is.na(scaled_ps)] = 0
TF_exp_selected_hm = Heatmap (scaled_ps,
        #right_annotation=tf_mark,
        #column_split = column_split_rna,
        cluster_rows = F, #km = 4, 
        name = 'expression (scaled)',
        column_gap = unit(.5, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=F, 
        col = palette_expression,
        row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
        column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
        border=T,
        width = unit(2, "cm"))

TF_exp_selected_hm2 = Heatmap (ps_tf,
        #right_annotation=tf_mark,
        #column_split = column_split_rna,
        cluster_rows = F, #km = 4, 
        name = 'expression',
        column_gap = unit(.5, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=F, 
        col = palette_expression,
        row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
        column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
        border=T,
        width = unit(2, "cm"))

pdf (file.path ('Plots','DAM_with_rna_expression_heatmaps.pdf'), width = 8,height=4)
draw (DAM_hm + TF_exp_selected_hm + TF_exp_selected_hm2)
dev.off()
   

### Co-expression of TFs across cells #### 

# # Get deviation matrix ####
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Subset only for expressed TFs ####
metaGroupName = 'celltype2'
ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat), group.by = metaGroupName)[[1]]) +1)
min_exp = 0.1
ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
active_TFs = rownames(ps)[rowSums(ps) > 0]
#positive_TF = corGSM_MM[,1][corGSM_MM[,3] > 0]
metaGroupName = 'celltype2'
mMat = mMat[active_TFs,]
mMat = as.data.frame (t(mMat))
mMat$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat = aggregate (.~ metaGroup, mMat, mean)
rownames (mMat) = mMat[,1]
mMat = mMat[,-1]

TF_hm = draw (Heatmap (scale(mMat), 
          #row_labels = colnames (mMat_mg),
          #column_title = paste('top',top_genes),
          clustering_distance_columns = 'pearson',
          clustering_distance_rows = 'pearson',
          cluster_rows = T,
          column_km = 20,
          #col = pals_heatmap[[5]],
          cluster_columns=T,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
          column_names_gp = gpar(fontsize = 0),
          column_names_rot = 45,
          name = 'chromVAR',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_deviation_centered#,
          #width = unit(2, "cm")
          #right_annotation = motif_ha
          ))

pdf (file.path ('Plots',paste0('TF_',metaGroupName,'heatmap.pdf')), width=8, height=1.6)
TF_hm
dev.off()

colnames(TF_hm@ht_list$chromVAR@matrix)[unlist(column_order(TF_hm)[c('2','3','4','5')])]
which (colnames(mMat) == 'NR4A2')
sapply (column_order(TF_hm), function(x) 497 %in% x)

# Distance matrix ####
d <- as.dist(1 - cor(t(mMat), method='pearson'))

# Hierarchical clustering ####
hc <- hclust(d)

# Dendrogram ####
pdf (file.path ('Plots',paste0('TF_',metaGroupName,'_no_km_dendrogram.pdf')), width=3, height=3.6)
plot(hc)
dev.off()
