conda activate meso_scatac
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

####### ANALYSIS of main data #######
projdir = 'main'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)


#devtools::install_github("immunogenomics/presto") #needed for DAA
source (file.path('..','git_repo','utils','useful_functions.R'))
source (file.path('..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','git_repo','utils','palettes.R'))
source (file.path('..','git_repo','utils','hubs_track.R'))

set.seed (1234)
addArchRThreads (threads = 1) 
addArchRGenome ("Hg38")

# Load RNA
srt = readRDS ('srt.rds')
srt$celltype_simplified2[srt$celltype_simplified2 == 'pDC'] = 'pDCs'

archp = loadArchRProject (projdir)

# Set order of celltype for displaying purposes ####
celltype_order = c('Malignant','Mesothelium','Alveolar','Fibroblasts','SmoothMuscle','Endothelial','Myeloid','T_cells','NK','B_cells','Plasma','pDCs')

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
umap_p1 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype_lv1",
   embedding = "UMAP",
   pal = palette_celltype_lv1,
   labelMeans = FALSE)

umap_p2 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
   embedding = "UMAP",
   pal = palette_sample,
   labelMeans = FALSE)
dev.off()

pdf (file.path ('Plots','celltype_lv1_umap.pdf'))
umap_p1
umap_p2
dev.off()


# Plot gene score of cell type markers ####
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

if (!any (ls() == 'gsSE')) gsSE = fetch_mat (archp, mat = 'GeneScore')
gsMat = assays (gsSE)[[1]]

metaGroupName = 'celltype_lv1'
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
pdf (file.path('Plots',paste0('DAG_clusters_',metaGroupName,'_heatmaps2.pdf')), width = 6, height = 3)
print(DAG_hm)
dev.off()



# Generate CNV map for each sample ####
# Get Granges of blacklist regions 
blacklist = toGRanges (file.path('..','..','git_repo','files',"hg38-blacklist.v2.bed")) # taken from https://github.com/Boyle-Lab/Blacklist/tree/master/lists
projdir_cnv = file.path('..','..','per_sample_QC_signac','CNV_analysis')
dir.create (projdir_cnv)

# Get GRanges of bins excluding black list regions
ws = 10e6
ss = 5e6
if (!file.exists (file.path (projdir_cnv, paste0 ('windows_',ws,'_',ss,'.rds'))))
  {
  windows = makeWindows (genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist,
    windowSize = ws, slidingSize = ss)
  saveRDS (windows, file.path (projdir_cnv, paste0 ('windows_',ws,'_',ss,'.rds')))
  } else {
  windows = readRDS (file.path (projdir_cnv, paste0 ('windows_',ws,'_',ss,'.rds')))
  }

deleted_chr = c('chr4','chr22','chr13')
if (!exists ('fragments_l')) fragments_l = getFragmentsFromProject (archp)
print_mat = F
force=F

# Loop per sample and run CNV analysis
sample_names = unique (archp$Sample)

cnaObj_l = list()
for (sam in sample_names)
  {
  # Run CNV analysis
  #force=FALSE
  if (!file.exists (file.path(projdir_cnv, paste0('CNV_LFC_GC_',sam,'_ws_',ws,'_ss_',ss,'.rds'))) | force)
    {
    cnaObj_l[[sam]] = scCNA (windows, fragments_l[[sam]], neighbors = 100, LFC = 1.5, FDR = 0.1, force = FALSE, remove = c("chrX","chrM","chrY"))
    saveRDS(cnaObj_l[[sam]], file.path(projdir_cnv, paste0('CNV_LFC_GC_',sam,'_ws_',ws,'_ss_',ss,'.rds')))
    } else {
    message ('cnaObj found! loading...')
    cnaObj_l[[sam]] = readRDS (file.path(projdir_cnv, paste0('CNV_LFC_GC_',sam,'_ws_',ws,'_ss_',ss,'.rds')))
    }
  }

# Print CNV heatmap
malignant_cells = c(rownames(archp@cellColData)[archp$celltype_lv1 == 'Malignant'])
  
cnaMat_obj_l = list()
for (sam in sample_names)
    {
    mat_type = 'z'
    cnaObj_mat = t(cnaObj_l[[sam]]@assays@data@listData[[mat_type]])
    colnames (cnaObj_mat) = paste0(seqnames(rowRanges(cnaObj_l[[sam]])),':', ranges(rowRanges(cnaObj_l[[sam]])))
    bc = colnames(cnaObj_l[[sam]])
    if (any (grep ('\\#', colnames(cnaObj_l[[sam]])))) bc = sapply (colnames(cnaObj_l[[sam]]), function(x) unlist(strsplit(x, '\\#'))[2])
    rownames (cnaObj_mat) = paste0(sam,'#',bc)
    cnaObj_mat = cnaObj_mat[rownames(cnaObj_mat) %in% malignant_cells,]
    cnaObj_mat[is.na(cnaObj_mat)] = 0
    cnaObj_mat[is.infinite(cnaObj_mat)]= 0
    cnaMat_obj_l[[sam]] = cnaObj_mat
    }
cnaMat_obj_avg = lapply (cnaMat_obj_l, function(x) colMeans (x))
cnaMat_obj_avg = do.call (cbind, cnaMat_obj_avg)
cnaMat_obj_avg = cnaMat_obj_avg[,sample_names]
ha = rowAnnotation (sample = colnames(cnaMat_obj_avg), col = list(sample = palette_sample))
pdf (file.path('Plots', paste0('GL_method_',mat_type,'_CNV_heatmap2.pdf')), width=6,height=3)
    print (Heatmap (t(cnaMat_obj_avg), 
    cluster_rows=F,#col=col_fun,
    cluster_columns=F, 
    left_annotation = ha, 
    row_names_side = 'left',
    column_title_side = 'bottom',
    border=T,
     #row_km = 2,
    row_names_gp = gpar (fontsize = 8),
    column_names_gp = gpar(fontsize = 0),
    column_title_gp = gpar(fontsize = 7),
    column_split = row_chr,
    row_split = ifelse ('normal1' == colnames(cnaMat_obj_avg),'normal','tumor'),
    column_gap = unit(.5, "mm"),
    cluster_column_slices = F,
    row_title_gp = gpar(fontsize = 1),
    column_title_rot=45))
dev.off()


### Run peak calling ####
metaGroupName = "Clusters"
force=TRUE
peak_reproducibility=2
if(!all(file.exists(file.path('PeakCalls', paste0(unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds')))) | force) source ('../../git_repo/utils/callPeaks.R')
  

# Make barplots of number of peaks per cell type ####
ct_peaks = lapply (file.path('PeakCalls',paste0(names(palette_celltype_lv1)[names(palette_celltype_lv1) %in% unique(archp$celltype_lv1)],'-reproduciblePeaks.gr.rds')), function(x) readRDS (x))

ct_peaks = do.call (cbind,  lapply(ct_peaks, function(x) table (x$peakType)))
colnames (ct_peaks) = names(palette_celltype_lv1)[names(palette_celltype_lv1) %in% unique(archp$celltype_lv1)]

ct_peaks = as.data.frame (ct_peaks)
ct_peaks$peaktype = rownames(ct_peaks)
ct_peaks = gather (ct_peaks, celltype, peaks, 1:(ncol(ct_peaks)-1))

df = ct_peaks %>% 
group_by (celltype) %>% 
summarize (sum = sum(peaks)) %>% 
arrange(sum)

ct_peaks$celltype = factor (ct_peaks$celltype, levels = rev(df$celltype[order(df$sum)]))

bp = ggplot (ct_peaks, aes(x = celltype, y = peaks, fill = peaktype)) + 
geom_bar (stat= 'identity') + 
paletteer::scale_fill_paletteer_d("nbapalettes::knicks_city") +
gtheme

pdf (file.path ('Plots','peakcalls.pdf'))
bp
dev.off()
  
### chromVAR analysis ####
run_chromVAR = TRUE

if (run_chromVAR)
  {  
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
  }



### ChromVAR based analysis ####


  # Find DAM ####
  metaGroupName = "celltype_lv1"
  force = F

  top_genes = Inf
  DAM_df = DAM (
  ArchRProj = archp,
  metaGroupName = metaGroupName,
  FDR_threshold = 1e-2,
  meandiff_threshold = 0,
  top_genes=top_genes,
  filter_by_scRNA=TRUE, # Make sure has same metaGroupName
  seurat_obj = srt,
  min_exp=.1,
  force = force)

# Save table for supplementary information
write.csv (DAM_df, paste0('DAM_table_',metaGroupName, '.csv'))

# Take only top five to show heatmap
DAM_df <- DAM_df %>%
  mutate(comparison = factor(comparison, levels = celltype_order)) %>%
  group_by(comparison) %>%
#  arrange(desc(Log2FC), .by_group = TRUE) %>%
  slice_head(n = 5) %>%   # keep top 5 per celltype
  ungroup() %>%
  arrange(comparison)

if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames(mMat) = rowData(mSE)$name

mMat = mMat[unique (DAM_df$gene), ]

#mMat_mg = mMat[active_DAM, ]
mMat_mg = as.data.frame (t(mMat))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
mMat_mg = mMat_mg[celltype_order,]

# Generate heatmap ####

 DAM_hm = Heatmap (t(t(scale(mMat_mg))), 
          row_labels = rownames (mMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          row_names_side = 'left',
          cluster_rows = F,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          name = 'chromVAR',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = rev(palette_deviation)#,
          #width = unit(2, "cm")
          #right_annotation = motif_ha
          )

#DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_heatmaps2.pdf')), width = 10, height = 3)
print(DAM_hm)
dev.off()

### Subset ArchR project ####
run_dropcells = FALSE
if (run_dropcells) archp = saveArchRProject (archp, dropCells = T) # Make sure to run this first before subsetting with the custom subset function

# Subset T cells ####
metaGroupName = 'celltype_lv1'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('T_cells','NK')],
  outputDirectory = file.path('..','..','NKT_cells'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)
saveRDS (srt[,srt$celltype_lv1 %in% c('T_cells','NK')], file.path ('..','..','NKT_cells','srt.rds'))

# Subset Myeloid ####
metaGroupName = 'celltype_lv1'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('Myeloid')],
  outputDirectory = file.path('..','..','myeloid_cells'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)
saveRDS (srt[,srt$celltype_lv1 == 'Myeloid'], file.path ('..','..','Myeloid','srt.rds'))

# Subset Stroma ####
metaGroupName = 'celltype_lv1'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('Endothelial','Fibroblasts','Mesothelium','SmoothMuscle')],
  outputDirectory = file.path('..','..','stroma'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)
saveRDS (srt[,srt$celltype_lv1 %in% c('Endothelial','Fibroblasts','Mesothelium','SmoothMuscle')], file.path ('..','..','stroma','srt.rds'))

# Subset Malignant ####
metaGroupName = 'celltype_lv1'
subsetArchRProject(
  ArchRProj = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('Malignant')],
  outputDirectory = file.path('..','..','tumor_compartment'),
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)
saveRDS (srt[,srt$celltype_lv1 == 'Malignant'], file.path ('..','..','tumor_compartment','srt.rds'))

metaGroupName='celltype_lv1'
subsetArchRProject_light (ArchRProject = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('Malignant')],
  projdir_new = file.path('..','..','tumor_compartment')
  )

# Subset B cells ####
metaGroupName='celltype_lv1'
subsetArchRProject_light (ArchRProject = archp,
  cells = rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% c('B_cells','Plasma')],
  projdir_new = file.path('..','..','B_cells')
  )
saveRDS (srt[,srt$celltype_lv1 %in% c('B_cells','Plasma')], file.path ('..','..','B_cells','srt.rds'))


# Show that enhancer linked to NR4A2 is only up in NK KLRC1 and CD8 exhausted across all cells ####
metaGroupName = 'celltype_lv2'
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
enhancer_region = GRanges ('chr2:156480366-156480866')
promoter_region = getPeakSet(archp)
promoter_region =  promoter_region[which(promoter_region$peakType == 'Promoter' & promoter_region$nearestGene == 'NR4A2')]
promoter_region = GRanges (seqnames = as.character(seqnames(promoter_region))[1], ranges= IRanges(start = min(start(promoter_region)), end = max(end(promoter_region))))

pmat_enhancer = unlist(as.data.frame (assay (pMats[queryHits(findOverlaps (GRanges(rowData(pMats)), enhancer_region)),])))
pmat_enhancer_df = data.frame (region = pmat_enhancer, celltype =names(pmat_enhancer))
pmat_enhancer_df$celltype = factor (pmat_enhancer_df$celltype, pmat_enhancer_df$celltype[order(-pmat_enhancer_df$region)])
pmat_enhancer_df$type = 'enhancer'
pmat_promoter = unlist(colSums (as.data.frame (assay (pMats[queryHits(findOverlaps (GRanges(rowData(pMats)), promoter_region)),]))))
pmat_promoter_df = data.frame (region = pmat_promoter, celltype =names(pmat_promoter))

pmat_promoter_df$type = 'promoter'
pmat_df = rbind (pmat_enhancer_df, pmat_promoter_df)
pmat_df$type = factor (pmat_df$type, levels =c ('promoter','enhancer'))

pdf (file.path ('Plots','eNR4A2_accessibility_barplot.pdf'), height=3, width=6)
ep = ggplot (pmat_df, aes (x = celltype, y = region, fill = type)) + 
geom_bar(stat = 'identity',position ='stack', color='grey22') + gtheme +
scale_fill_manual (values = c(enhancer = 'darkred', promoter = 'grey'))
#+ scale_fill_manual (values = c(palette_tnk_cells, palette_myeloid, palette_celltype_lv1))
ep
dev.off()




### Show cell type proportions per sample ####
archp_meta = as.data.frame (archp@cellColData)
archp_meta$celltype_lv1 = factor (archp_meta$celltype_lv1, levels = names(palette_celltype_lv1))
archp_meta$Sample2 = archp_meta$Sample
archp_meta$Sample2 = factor (archp_meta$Sample2, levels= rev(c('P1','P3','P4','P5','P8','P10','P11','P12','P13','P14','P23')))
bp = cellComp(
  seurat_obj = archp_meta,
  metaGroups = c('Sample2','celltype_lv1'),
  plot_as = 'bar',
  pal = palette_celltype_lv1
  ) + gtheme + coord_flip()

pdf (file.path('Plots',paste0('celltype_barplots.pdf')), height=6, width=5)
bp
dev.off()






# Import chromBPnet finemo calls ####
#archp_P1 = archp[archp$Clusters %in% c('C2') & archp$Sample == 'P1']
library ('universalmotif')

ps = getPeakSet (archp)

chromBPdir = 'chromBPnet'

chrombpnet_counts = list()
celltypes = unique (archp$celltype_lv1)
celltypes = c('Malignant','Mesothelium','Alveolar','Fibroblasts','SmoothMuscle','Endothelial','Myeloid','T_cells','NK','B_cells','Plasma','pDCs') 
#celltypes = celltypes[celltypes != 'pDCs']
annotated_motifs = read.table (file.path(chromBPdir,'modisco_merged_counts','compiled','modisco_compiled.tsv'), sep='\t', header=T)
#annotated_motifs$pattern2 = sapply (annotated_motifs$pattern, function(x) unlist(strsplit(x,'__'))[2])
#rownames (annotated_motifs) = annotated_motifs$pattern2
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_counts[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
  chrombpnet_counts[[celltype]]$V8 = sub ('pos_patterns.','',chrombpnet_counts[[celltype]]$V6)
  chrombpnet_counts[[celltype]]$V8 = sub ('neg_patterns.','',chrombpnet_counts[[celltype]]$V8)
  chrombpnet_counts[[celltype]]$V9 = paste(annotated_motifs$match0,annotated_motifs$match1,annotated_motifs$match2)[match(chrombpnet_counts[[celltype]]$V8, annotated_motifs$pattern)]
  chrombpnet_counts[[celltype]]$V9 <- gsub('_HUMAN\\.H11MO\\.(0|1|2)\\.[A-D]', '', chrombpnet_counts[[celltype]]$V9)

  gr = makeGRangesFromDataFrame (chrombpnet_counts[[celltype]], keep.extra.columns=T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  chrombpnet_counts[[celltype]]$peak_type = ps$peakType[findOverlaps(gr, ps, select='first')]
  chrombpnet_counts[[celltype]]$direction = ifelse (grepl ('pos_pattern', chrombpnet_counts[[celltype]]$V6), 'positive','negative')
  #chrombpnet_counts[[celltype]] = chrombpnet_counts[[celltype]][chrombpnet_counts[[celltype]]$V4 != 'NaN_NaN_NaN',] # Some seqlets have NA motif match and qvalues...removing those
  #nonsig_motifs = chrombpnet_counts[[celltype]]$V7 > 0.05
  #chrombpnet_counts[[celltype]]$V4[nonsig_motifs] = chrombpnet_counts[[celltype]]$V6[nonsig_motifs]
  }

# Export finemo hits as tsv files 
lapply (names(chrombpnet_counts), function(x) write.table (chrombpnet_counts[[x]][,c(1:3,9)],file.path ('chromBPnet', x,'no_bias_model',paste0(x, '_finemo_counts_to_genome_browser_harmonized.tsv')), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE))

annotated_motifs_profile = read.table (file.path(chromBPdir,'modisco_merged_profile','compiled','modisco_compiled.tsv'), sep='\t', header=T)

chrombpnet_profile = list()
#celltypes = c('Mesothelium','Malignant')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  chrombpnet_profile[[celltype]]$V8 = sub ('pos_patterns.','',chrombpnet_profile[[celltype]]$V6)
  chrombpnet_profile[[celltype]]$V8 = sub ('neg_patterns.','',chrombpnet_profile[[celltype]]$V8)
  chrombpnet_profile[[celltype]]$V9 = paste(annotated_motifs_profile$match0,annotated_motifs_profile$match1, annotated_motifs_profile$match2)[match(chrombpnet_profile[[celltype]]$V8, annotated_motifs_profile$pattern)]
  chrombpnet_profile[[celltype]]$V9 <- gsub('_HUMAN\\.H11MO\\.(0|1|2)\\.[A-D]', '', chrombpnet_profile[[celltype]]$V9)
  gr = makeGRangesFromDataFrame (chrombpnet_profile[[celltype]], keep.extra.columns=T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  chrombpnet_profile[[celltype]]$peak_type = ps$peakType[findOverlaps(gr, ps, select='first')]
  chrombpnet_profile[[celltype]]$direction = ifelse (grepl ('pos_pattern', chrombpnet_profile[[celltype]]$V6), 'positive','negative')
  }

lapply (names(chrombpnet_profile), function(x) write.table (chrombpnet_profile[[x]][,c(1:3,9)],file.path (chromBPdir, x,'no_bias_model',paste0(x, '_finemo_profile_to_genome_browser_harmonized.tsv')), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE))

### Check active motifs in selected peaks from hubs ####
hubs_obj = readRDS (file.path('hubs_obj_cor_0.3_md_12500_dgs_0_min_peaks_5','global_hubs_obj.rds'))  
selected_hubs = names(hubs_obj$hubsCollapsed[grep ('WT1', hubs_obj$hubsCollapsed$gene)])
selected_peaks = do.call (rbind, hubs_obj$hubsClusters[[1]][as.numeric(selected_hubs)])
selected_peaks = makeGRangesFromDataFrame (selected_peaks)
chrombpnet_counts_selected = lapply (chrombpnet_counts, function(x) {
  gr = makeGRangesFromDataFrame (x, seqnames.field= 'V1', start.field = 'V2', end.field = 'V3')
  x[queryHits(findOverlaps(gr, selected_peaks)),]
})

chrombpnet_profile_selected = lapply (chrombpnet_profile, function(x) {
  gr = makeGRangesFromDataFrame (x, seqnames.field= 'V1', start.field = 'V2', end.field = 'V3')
  x[queryHits(findOverlaps(gr, selected_peaks)),]
})


selected_celltypes = c('Mesothelium','Fibroblasts')
chrombpnet_counts_selected = chrombpnet_counts_selected[selected_celltypes]
top_n <- 5

bp_list <- lapply(names (chrombpnet_counts_selected), function(i) {
  tbl <- table(chrombpnet_counts_selected [[i]]$V9)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_counts_selected[[i]]$direction[chrombpnet_counts_selected [[i]]$V9 == tf][1]
  })
  
  data.frame(
    Freq = proportions(top_tbl),
    TF   = tf_names,
    direction = directions,
    type = rep(i, length(top_tbl)),
    pattern = chrombpnet_counts_selected[[i]]$V6[match(tf_names, chrombpnet_counts_selected[[i]]$V9)]
  )
})

bp_df <- do.call(rbind, bp_list)

# Make neg values negative
# bp_df <- bp_df %>%
#   mutate(Freq = ifelse(direction == "negative", -Freq.Freq, Freq.Freq))

# Create custom ordering per type
bp_df <- bp_df %>%
  group_by(type, direction) %>%
  mutate(
    TF_order = ifelse(direction == "positive",
                      rank(-Freq, ties.method = "first"),  # descending
                      rank(Freq, ties.method = "first"))   # ascending for neg (opposite)
  ) %>%
  ungroup()

# Build a combined factor: ensures pos stack from bottom up, neg from top down
bp_df <- bp_df %>%
  arrange(type, direction, TF_order)

bp_df$TF_id <- bp_df$TF
bp_df$TF_id <- factor(bp_df$TF_id, levels = unique(bp_df$TF_id))
#bp_df$type = factor (bp_df$type, levels = celltypes)
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq.Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("trekcolors::lcars_2379",8)) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

  

pdf (file.path ('Plots', 'TF_abundance_counts_selected_peaks_barplot.pdf'),6,width=5.5)
bp
dev.off()

selected_celltypes = c('Mesothelium','Fibroblasts')
chrombpnet_profile_selected = chrombpnet_profile_selected[selected_celltypes]
top_n <- 10
#n <- length(chrombpnet_profile_selected)

bp_list <- lapply(names(chrombpnet_profile_selected), function(i) {
  tbl <- table(chrombpnet_profile_selected [[i]]$V9)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_profile_selected[[i]]$direction[chrombpnet_profile_selected [[i]]$V9 == tf][1]
  })
  
  data.frame(
    Freq = proportions(top_tbl),
    TF   = tf_names,
    direction = directions,
    type = rep(i, length(top_tbl)),
    pattern = chrombpnet_profile_selected[[i]]$V6[match(tf_names, chrombpnet_profile_selected[[i]]$V9)]
  )
})

bp_df <- do.call(rbind, bp_list)

# Make neg values negative
bp_df <- bp_df %>%
  mutate(Freq = ifelse(direction == "negative", -Freq.Freq, Freq.Freq))

# Create custom ordering per type
bp_df <- bp_df %>%
  group_by(type, direction) %>%
  mutate(
    TF_order = ifelse(direction == "positive",
                      rank(-Freq, ties.method = "first"),  # descending
                      rank(Freq, ties.method = "first"))   # ascending for neg (opposite)
  ) %>%
  ungroup()

# Build a combined factor: ensures pos stack from bottom up, neg from top down
bp_df <- bp_df %>%
  arrange(type, direction, TF_order)

bp_df$TF_id <- bp_df$TF
bp_df$TF_id <- factor(bp_df$TF_id, levels = unique(bp_df$TF_id))
#bp_df$type = factor (bp_df$type, levels = celltypes)
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(unique(bp_df$TF)))) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle("Top 10 TFs (pos vs neg, ordered stacks)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

pdf (file.path ('Plots', 'TF_abundance_profile_selected_peaks_barplot.pdf'),height = 6, width = 5.5)
bp
dev.off()


