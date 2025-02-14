conda activate meso_scatac
R

set.seed(1234)

####### ANALYSIS of NKT compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/NK_cell_line_scatac_AH'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

# Load functions for hub detection ####
source (file.path('..','..','git_repo','utils','knnGen.R'))
source (file.path('..','..','git_repo','utils','addCoax.R'))
source (file.path('..','..','git_repo','utils','Hubs_finder.R'))
source (file.path('..','..','git_repo','utils','hubs_track.R'))
#source (file.path('..','..','git_repo','utils','scATAC_functions.R'))

# Set # of threads and genome reference ####
addArchRThreads(threads = 8) 
addArchRGenome("hg38")


sample_names = c(
  'AMHO18_247AZA_0_v1',
  'AMHO18_247untreated_0_v1',
  'AMHO18_199AZA_0_v1',
  'AMHO18_245untreated_0_v1',
  'AMHO18_245AZA_0_v1')
fragment_paths = paste0('/sc/arion/projects/nmibc_bcg/NK_azacitidine/ATACseq/RAIatac/cellranger_output/',sample_names,'/fragments.tsv.gz')

      #setwd (projdir)  
ArrowFiles = createArrowFiles (inputFiles = fragment_paths,
sampleNames = sample_names,
minTSS = 6, #Dont set this too high because you can always increase later
minFrags = 1000,
maxFrags = Inf,
addTileMat = TRUE,
addGeneScoreMat = TRUE,
force = TRUE,
subThreading = T
)

archp = ArchRProject (
  ArrowFiles = ArrowFiles, 
  outputDirectory = projdir,
  copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

## Reduce dimension and harmonize ####
#archp = archp[!archp$Clusters_H %in% c('C4','C3','C20','C2')]
#archp = archp[!archp$Clusters_H %in% c('C1','C9')]
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
    name = "Harmony_project",
    groupBy ='Sample', force=TRUE
)

archp = addUMAP (ArchRProj = archp,
    reducedDims = "Harmony_project", name='UMAP_H',
    force = TRUE)

archp = addClusters (input = archp,
    reducedDims = "Harmony_project",
    name='Clusters_H', res=1.5,
    force = TRUE)

pdf()
umap_p3 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
   embedding = "UMAP_H")
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters_H",
   embedding = "UMAP_H")
dev.off()

  pdf (file.path('Plots','celltype_umap_harmony_on_project_sample.pdf'),5,5)
  print (umap_p3)
  print (umap_p5)
  dev.off()

# Run genescore DAG ####
metaGroupName = "Clusters_H"
force = TRUE
if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source (file.path('..','..','git_repo','utils','DAG.R'))

# TNK markers ####
tnk_markers = c('CD3D','CD8A','PDCD1','HAVCR2','CD4', 'FOXP3','GNLY',
  'FGFBP2','KLRC1','XCL1','ICOS','NR4A2')
archp = addImputeWeights (archp)
pdf()
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = tnk_markers, 
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
p = lapply (seq_along(p), function(x) p[[x]] + theme_void() + ggtitle (tnk_markers[x]) + NoLegend())
#archp$celltype[archp$Clusters == 'C30'] = 'Fibroblasts_WT1'
#p = lapply (p, function(x) x + theme_void() + NoLegend ()) #+ ggtitle scale_fill_gradient2 (rev (viridis::plasma(100))))

pdf (file.path('Plots','TNK_markers_fplots.pdf'), width = 7, height = 13)
print (wrap_plots(p, ncol=3))
dev.off()


### Call peaks on celltypes ####
metaGroupName = 'Clusters_H'
force=TRUE
peak_reproducibility=2
if(!all(file.exists(file.path('PeakCalls', unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds'))) | force) 
source (file.path('..','..','git_repo','utils','callPeaks.R'))

### chromVAR analysis ####
force=TRUE
if (!all(file.exists(file.path('Annotations',
  c('Motif-Matches-In-Peaks.rds',
    'Motif-Positions-In-Peaks.rds',
    'Motif-In-Peaks-Summary.rds')))) | force)
source (file.path ('..','..','git_repo','utils','chromVAR.R'))
  


# Differential Accessed motifs ####
metaGroupName = "Clusters_H"
force=FALSE
source (file.path('..','..','git_repo','utils','DAM.R'))

if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames(mMat) = rowData(mSE)$name

# Filter by RNA expression ####
metaGroupName = 'celltype2'
active_TFs = exp_genes (srt, active_DAM, min_exp = 0.1, metaGroupName)
mMat = mMat[active_TFs, ]

#mMat_mg = mMat[active_DAM, ]
mMat_mg = as.data.frame (t(mMat))
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

pdf (file.path ('Plots','DAM_with_rna_expression_heatmaps.pdf'), width = 3,height=4)
draw (DAM_hm) # + TF_exp_selected_hm + TF_exp_selected_hm2)
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
active_TFs = rownames(ps)

mMat = mMat[active_TFs,]

set.seed (1234)
km = kmeans (t(scale(mMat)), centers=3)

ha = HeatmapAnnotation (df = data.frame (
  celltype = as.character(archp@cellColData[,metaGroupName]),
  sample = archp$Sample), col=list (celltype = palette_tnk_cells, sample=palette_sample))

TF_hm = Heatmap (scale(mMat), 
          top_annotation= ha,
          #row_labels = colnames (mMat_mg),
          #column_title = paste('top',top_genes),
          clustering_distance_columns = 'pearson',
          clustering_distance_rows = 'pearson',
          column_split = km$cluster,
          cluster_rows = T,
#          column_km = 3,
          #col = pals_heatmap[[5]],
          cluster_columns=T,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 0, fontface = 'italic'),
          column_names_gp = gpar(fontsize = 0),
          column_names_rot = 45,
          name = 'chromVAR',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_deviation_centered#,
          #width = unit(2, "cm")
          #right_annotation = motif_ha
          )

pdf (file.path ('Plots',paste0('TF_',metaGroupName,'heatmap2.pdf')), width=12, height=5)
TF_hm
dev.off()

archp$cell_kmeans = paste0('km',km$cluster)
pdf()
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "cell_kmeans",
   embedding = "UMAP_H")
dev.off()

pdf (file.path('Plots','celltype_kmeans_umap.pdf'),5,5)
print (umap_p5)
dev.off()





# # ### Use P2G analysis and cNMF from RNA to identify active TF via regulons  ####
run_p2g = T
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
    

### Hubs analysis #####
metaGroupName = "Clusters_H"
cor_cutoff = 0.2
#max_dist = 12500
max_dist = 12500
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
hubs_obj = readRDS (file.path(hubs_dir,'global_hubs_obj.rds'))  

# Import barcodes from RNA
rna_ann = read.csv ('../NK_cell_line_scrna_AH/_cellranger_raw_Filter_400_1000_25/no_harmony/barcode_annotation.csv')
head (rna_ann)
scatac_barcodes = rownames(archp@cellColData)
scatac_barcodes = sub ('#','_',scatac_barcodes)
table (scatac_barcodes %in% rna_ann$X)

anergic_sig = c('CCL4L2','CLL4','DUSP2','EGR2','GBP2','ZFP36L1','NR4A2','TNFRSF9','NFBID','CCL3','LYST',
  .'ID2','NR4A3','FASLG','SLA','CBLB','KLRD1','CCL3L1','FYN','PTPN7')

# Import DEG from scRNA
deg_list = readRDS ('../NK_cell_line_scrna_AH/_cellranger_raw_Filter_400_1000_25/no_harmony/DEG_RNA_snn_res.2_logFC_0.25_padj_0.01/DEG.rds')
deg_list= deg_list[deg_list$gene %in% getFeatures(archp),]
deg_list = deg_list[deg_list$avg_log2FC >0,]
deg_list = deg_list[deg_list$p_val_adj < 0.05,]
deg_list = lapply (split(deg_list, deg_list$cluster), function(x) head (x$gene, 20))

### Use gene scores from RNA markers to drive annotation of endothelial cells ####
archp = addModuleScore (
  ArchRProj = archp,
  useMatrix = 'GeneScoreMatrix',
  name = "scrna",
  features = deg_list,
  nBin = 25,
  nBgd = 100,
  seed = 1,
  threads = getArchRThreads(),
  logFile = createLogFile("addModuleScore")
)


# TNK markers ####
archp = addImputeWeights (archp)
pdf()
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "cellColData", 
    name = paste0('scrna.',names (deg_list)), 
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
p = lapply (seq_along(p), function(x) p[[x]] + theme_void() + ggtitle (tnk_markers[x]) + NoLegend())

pdf (file.path('Plots','TNK_modules_fplots.pdf'), width = 17, height = 18)
print (wrap_plots(p))
dev.off()

pdf()
qc_param = c('TSSEnrichment','ReadsInTSS','FRIP')
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "cellColData", 
    name = qc_param, 
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)
dev.off()
p = lapply (seq_along(p), function(x) p[[x]] + theme_void() + ggtitle (tnk_markers[x]) + NoLegend())

pdf (file.path('Plots','QC_fplots.pdf'), width = 17, height = 18)
print (wrap_plots(p))
dev.off()




archp_meta = as.data.frame (archp@cellColData[, paste0('scrna.',names(deg_list))])
archp_meta = scale(t(scale(t(archp_meta))))
archp$mod17 = ifelse (archp_meta[,'scrna.17'] > summary(archp_meta[,'scrna.17'])[4],'NR4A2_high','NR4A2_low')

metaGroupName = 'mod17'
TF = 'NR4A2'
pdf()
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp, 
    sizes = c(6, 1, 1, 1,1,1),
    groupBy = metaGroupName, 
    geneSymbol = TF,
    normMethod = "ReadsInTSS",
    scCellsMax=3000,
    plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack"#, 
        #"hubTrack",
        #'hubregiontrack'
        ),
    hubs_regions = NULL,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 50000,
    #pal = palette_tnk_cells_ext2,
    #ylim=c(0,0.1),
    downstream = 300000,
    #loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
#    loops = getCoAccessibility (archp, corCutOff = 0.25),
    #  returnLoops = TRUE),
    useGroups= NULL,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2,returnLoops = TRUE),
    #hubs = hubs_obj$peakLinks2
)
dev.off()
plotPDF (meso_markers, ArchRProj = archp, 
  width=5,height=3, 
  name =paste0(paste(TF, collapse='_'),'_coveragePlots.pdf'),
  addDOC = F)

