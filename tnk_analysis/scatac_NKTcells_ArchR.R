conda activate meso_scatac
R

set.seed(1234)

####### ANALYSIS of NKT compartment #######
projdir = 'NKT_cells'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','git_repo','utils','load_packages.R'))
source (file.path('..','git_repo','utils','useful_functions.R'))
source (file.path('..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','git_repo','utils','palettes.R'))

# Load functions for hub detection ####
source (file.path('..','git_repo','utils','knnGen.R'))
source (file.path('..','git_repo','utils','addCoax.R'))
source (file.path('..','git_repo','utils','Hubs_finder.R'))
source (file.path('..','git_repo','utils','hubs_track.R'))

# Set # of threads and genome reference ####
addArchRThreads(threads = 1) 
addArchRGenome("hg38")

archp = loadArchRProject (projdir)
srt = readRDS ('srt.rds')

cell_subsets_order = c('CD4','Tregs','CD8','CD8_exhausted','FGFBP2_NK','KLRC1_NK')

## Reduce dimension and harmonize ####

  # Dimensionality reduction and clustering
  varfeat = 25000
  LSI_method = 2
  archp = addIterativeLSI (ArchRProj = archp,
    useMatrix = "TileMatrix", name = "IterativeLSI",
    force = FALSE, LSIMethod = LSI_method,
    varFeatures = varfeat)

  archp = addHarmony (
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    name = "Harmony_project",
    groupBy ='Sample', force=FALSE
)

archp = addUMAP (ArchRProj = archp,
    reducedDims = "Harmony_project", name='UMAP_H',
    force = FALSE)

archp = addClusters (input = archp,
    reducedDims = "Harmony_project",
    name='Clusters_H', res=2,
    force = FALSE)


# Run genescore DAG ####
metaGroupName = "Clusters_H"
force = FALSE
if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source (file.path('..','..','git_repo','utils','DAG.R'))

# TNK markers ####
tnk_markers = c('CD4','CD8A',
  'PDCD1','FOXP3','GNLY', 'HAVCR2',
  'FGFBP2','KLRC1','CTLA4','GZMK',
  'CXCL13','CCR7','TCF7','TOX','IFNG',
  'CX3CR1','TBX21','ASCL2','IL7R',
  'CXCR4','NR4A2','TIGIT','ICOS','ENTPD1','LAG3','TNFRSF9','CD28')
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

pdf (file.path('Plots','TNK_markers_fplots.pdf'), width = 17, height = 13)
print (wrap_plots(p, ncol=3))
dev.off()

### Annotate meso cells ####
archp$celltype_lv2 = archp$Clusters_H
archp$celltype_lv2[archp$celltype_lv2 %in% c('C7')] = 'FGFBP2_NK'
archp$celltype_lv2[archp$celltype_lv2 %in% c('C8','C9')] = 'KLRC1_NK'
archp$celltype_lv2[archp$celltype_lv2 %in% c('C16')] = 'Tregs'
archp$celltype_lv2[archp$celltype_lv2 %in% c('C14','C20','C11','C17','C12','C13')] = 'CD4'
archp$celltype_lv2[archp$celltype_lv2 %in% c('C1','C5','C6','C18','C15','C19','C10')] = 'CD8'
archp$celltype_lv2[archp$celltype_lv2 %in% c('C3','C2','C4')] = 'CD8_exhausted'

write.csv (data.frame (barcode = rownames(archp@cellColData), celltype_lv2 = archp$celltype_lv2), 'barcode_annotation.csv')


pdf()
umap_p3 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Sample",
  pal = palette_sample,
   embedding = "UMAP_H")
umap_p4 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype_lv2",
   embedding = "UMAP_H", pal = palette_tnk_cells)
umap_p5 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters_H",
   embedding = "UMAP_H")
dev.off()

  pdf (file.path('Plots','celltype_umap_harmony_on_project_sample.pdf'),5,5)
  print (umap_p3)
  print (umap_p4)
  print (umap_p5)
  dev.off()

# Check expression of GZMB PRF1 and KLRC1 ####
metaGroupName = 'celltype_lv2'
archp = addImputeWeights (archp)
markers = c('GZMA','GZMB','GZMK','GZMH','IFNG','TNF','PRF1',
  'TOX','CD96','HAVCR2','CTLA4','ENTPD1','KLRC1','LAYN')
activating_markers = c('GZMA','GZMB','GZMH','PRF1','ITGAE','ITGA1')
inhibitory_markers = c('KIR3DL1','KIR2DL3')

pdf()
p2 <- plotGroups(
    ArchRProj = archp,#[archp$celltype_lv2 %in% c('NK_KLRC1','NK_FGFBP2')], 
    groupBy = metaGroupName, 
    colorBy = "GeneScoreMatrix", 
    name = markers,
    plotAs = "violin",
    alpha = 0.4,
    imputeWeights = getImputeWeights(archp),
    addBoxPlot = TRUE,
    pal = palette_tnk_cells
   )
dev.off()
p2 = lapply (p2, function(x) x + gtheme)
p2 = lapply (p2, function(x) {x$data = x$data[x$data[,1] %in% c('FGFBP2_NK','KLRC1_NK'),]; x})


pdf (file.path ('Plots','NK_cytotoxicity_markers.pdf'), width=12,height=6)
wrap_plots (p2, ncol=7)
dev.off()


### Call peaks on celltypes ####
metaGroupName = 'celltype_lv2'
force=FALSE
peak_reproducibility=2
if(!all(file.exists(file.path('PeakCalls', unique(archp@cellColData[,metaGroupName]), '-reproduciblePeaks.gr.rds'))) | force) 
source (file.path('..','git_repo','utils','callPeaks.R'))

### chromVAR analysis ####
force=TRUE
if (!all(file.exists(file.path('Annotations',
  c('Motif-Matches-In-Peaks.rds',
    'Motif-Positions-In-Peaks.rds',
    'Motif-In-Peaks-Summary.rds')))) | force)
source (file.path ('..','..','git_repo','utils','chromVAR.R'))
  


# Differential Accessed motifs ####

  # Find DAM ####
  metaGroupName = "celltype_lv2"  
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
  mutate(comparison = factor(comparison, levels = cell_subsets_order)) %>%
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
mMat_mg = mMat_mg[cell_subsets_order,]

# Generate heatmap ####

 DAM_hm = Heatmap (t(scale(mMat_mg)), 
          row_labels = colnames (mMat_mg),
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
pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_heatmap.pdf')), width = 3, height = 4)
print(DAM_hm)
dev.off()



### Run TF correlation to identify TF modules across TNK cells #### 
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat)#[selected_TF,])

# Filter by RNA expression ####
metaGroupName = 'celltype_lv2'
active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)
mMat = mMat[active_TFs, ]

mMat_cor = cor (as.matrix(t(scale(mMat))), method = 'spearman')

set.seed(123)
centers=3
km = kmeans (mMat_cor, centers=centers)
write.csv (patchvecs (split (names(km$cluster),km$cluster)), 'regTNK_modules.csv')

pdf (file.path ('Plots','TF_modules_heatmap2.pdf'), width = 4,height=3)
cor_mMat_hm = draw (Heatmap (mMat_cor,# row_km=15,
  #left_annotation = ha,
  #rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=palette_deviation_cor_fun, 
  row_split = km$cluster,
  column_split = km$cluster,
  #row_km=2, 
  #column_km=2,
#  right_annotation = ha,
  border=T,
#   ,
  row_names_gp = gpar(fontsize = 0), 
  column_names_gp = gpar(fontsize = 0)
  ))
dev.off()

pdf (file.path ('Plots','TF_modules_heatmap3.pdf'), width = 4,height=3)
cor_mMat_hm
dev.off()

tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = do.call (cbind, tf_modules)
archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% paste0('mod_',unique(km$cluster))]
archp@cellColData = cbind (archp@cellColData, tf_modules) 

pdf()
TF_p = plotEmbedding (
    ArchRProj = archp,
    colorBy = "cellColData",
    name = paste0('mod_',unique(km$cluster)), 
    pal = rev(palette_deviation_correlation),
    #useSeqnames='z',
    embedding = "UMAP_H")
dev.off()
pdf (file.path ('Plots','TF_modules_umap3.pdf'), width = 20,height=6)
wrap_plots (TF_p, ncol=5)
dev.off()


# Differential Peaks in CD8 exhausted vs CD8 ####
# Find DAP ####
#force = FALSE
metaGroupName = 'celltype_lv2'
force=F
if (!file.exists ('DAP_CD8_CD8_ext_pairwise.rds') | force)
  {
  DAP_list = getMarkerFeatures (
    ArchRProj = archp, 
    testMethod = "wilcoxon",
          useGroups = "CD8_exhausted",
          bgdGroups = "CD8",
    k=100,
    binarize = FALSE,
    useMatrix = "PeakMatrix",
    groupBy = metaGroupName
  #  useSeqnames="z"
  )
  saveRDS (DAP_list, 'DAP_CD8_CD8_ext_pairwise.rds')
  } else {
  DAP_list = readRDS ('DAP_CD8_CD8_ext_pairwise.rds')
  }
DAP_res = do.call (cbind, (assays(DAP_list)))
colnames (DAP_res) = names(assays(DAP_list))
DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames(DAP_res) = as.character(DAP_res_regions)

pdf(file.path('Plots','Exhausted_peaks_MA_plot.pdf'), width=5,height=5)
pma <- markerPlot(seMarker = DAP_list, name = 'CD8_exhausted', cutOff = "FDR <= 0.01", plotAs = "MA")
pma
dev.off()

# Take only significant regions ####
DAP_res_sig = DAP_res[DAP_res$FDR < .01 & DAP_res$Log2FC > 0, ]
write.csv (DAP_res_sig, 'DAP_CD8_Ex.csv')
saveRDS (GRanges(rownames(DAP_res_sig)), 'T_cell_exhaustion_peaks.rds')


# Run GSEA enrichment analysis #### 
#!!! To run this analysis load only ArcHR and clusterprofiler packages !!!!
library (fgsea)    
options(warn = 0)
ps = getPeakSet (archp)
csv.file = file.path('..','git_repo','files','Tcell_exhaustion_genes_PMID37091230')
pathways = read.csv (csv.file)
pathways = list(Tcell_exhaustion_genes_PMID37091230 = pathways[,1])
#pathways = gmtPathways (gmt.file)
DAP_list = readRDS ('DAP_CD8_CD8_ext_pairwise.rds')
DAP_res = do.call (cbind, (assays(DAP_list)))
colnames (DAP_res) = names(assays(DAP_list))
DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames(DAP_res) = as.character(DAP_res_regions)

peak_genes = unname(ps[queryHits (findOverlaps(ps, GRanges (rownames(DAP_res))))]$nearestGene)
names (peak_genes) = as.character(ps)[queryHits(findOverlaps (ps,GRanges (rownames(DAP_res))))]
peak_genes = peak_genes[rownames(DAP_res)]
peak_genes = setNames (-log10(DAP_res$Pval) * sign (DAP_res$Log2FC), peak_genes)
#names (peak_genes) = 
peak_genes = peak_genes[!duplicated(names(peak_genes))]
peak_genes = peak_genes[!is.na(names(peak_genes))]
fgseaRes = fgseaMultilevel (pathways, 
          peak_genes,#, 
          minSize=15, 
          maxSize=1500,
          nPermSimple=100000,
          BPPARAM = NULL
          )
pdf (file.path ('Plots','DAP_exhaustion_enrichment_plot.pdf'), width=5, height=3)
plotEnrichment(pathways[["Tcell_exhaustion_genes_PMID37091230"]],
               peak_genes) + labs(title="Tcell_exhaustion_genes_PMID37091230")
dev.off()

# Calculate hypergeometric test of peaks with motifs ####
mm = getMatches (archp)
mmMat = assay(mm)
colnames(mmMat) = gsub ('_.*','',colnames(mmMat))
colnames(mmMat) = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", colnames(mmMat))

hyper_res = lapply (unique(km$cluster), function(mod) {
  module_peaks = rowRanges(mm)[which (rowSums (mmMat[,names(km$cluster[km$cluster == mod])]) > 0)]
  q = length (queryHits (findOverlaps (module_peaks, GRanges (rownames(DAP_res_sig)))))
  n = length (getPeakSet(archp)) - length(GRanges (rownames(DAP_res_sig)))
  m = length (GRanges (rownames(DAP_res_sig)))
  k = length (module_peaks)
  #hyper_res[[mod]] = 
  print (phyper (q, m, n, k, lower.tail = FALSE, log.p = FALSE))
  })

# Barplot of module enrichment in CD8 ext peaks ####
mod_df = data.frame (module = paste0('module',unique(km$cluster)), pval = -log10(unlist(hyper_res)+1e-9))
bp = ggplot (mod_df, aes (x = module, y = pval, fill = module)) + 
geom_bar (stat = 'identity',color='grey22') + 
scale_fill_manual (values = c(module1 = 'darkred',module2='grey',module3='grey')) + 
gtheme
pdf (file.path ('Plots','TFmodule_enrichments.pdf'), width=4,height=3)
bp
dev.off()

# Km 1 is enriched in CD8 exhausted peaks. Find which TF in it is mostly enriched ####
ext_module = 1
hyper_res_TF = lapply (names(km$cluster[km$cluster== ext_module]), function(mod) {
  module_peaks = rowRanges(mm)[which (rowSums (mmMat[,mod, drop=F]) > 0)]
  q = length (queryHits (findOverlaps (module_peaks, GRanges (rownames(DAP_res_sig)))))
  n = length(getPeakSet(archp)) - length(GRanges (rownames(DAP_res_sig)))
  m = length(GRanges (rownames(DAP_res_sig)))
  k = length(module_peaks)
  #hyper_res[[mod]] = 
  print (phyper (q, m, n, k, lower.tail = FALSE, log.p = FALSE))
  })
hyper_res_TF_df = data.frame (hyper = p.adjust (unlist(hyper_res_TF)), TF= names(km$cluster[km$cluster== ext_module]))
hyper_res_TF_df = hyper_res_TF_df[order(hyper_res_TF_df$hyper), ]

# Barplot of enriched TFs in ext module enrichment in CD8 ext peaks ####
hyper_res_TF_df$pval= -log10(hyper_res_TF_df$hyper)#+1e-9)
hyper_res_TF_df$TF = factor (hyper_res_TF_df$TF, levels = rev(hyper_res_TF_df$TF))
top_TF = 13
hyper_res_TF_df$label = ifelse (hyper_res_TF_df$TF %in% head(hyper_res_TF_df$TF, top_TF), as.character(head(hyper_res_TF_df$TF, top_TF)), '')

bp = ggplot (hyper_res_TF_df[hyper_res_TF_df$pval > -log10(0.05),], aes (x = TF, y = pval, fill = pval)) + 
geom_bar (stat = 'identity') + 
scale_fill_gradient(low = "white", high = "darkred") +
geom_text(aes(label = label, y = pval + .1), hjust = -0.05,na.rm = TRUE, angle=0,size=3) +
#ylim (c(0,11.5)) +
gtheme_no_rot + 
theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
coord_flip()

pdf (file.path ('Plots','TF_in_module_enrichments.pdf'), width=10,height=3)
bp
dev.off()





# # Take highest mean of NK KLRC1 + CD8 exhausted TFs ####
ext_module = 1 
tf_ext = names (km$cluster[km$cluster == ext_module])
metaGroupName = 'celltype_lv2'

if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = scale(mMat[tf_ext,])
mMat = as.data.frame (mMat)

mMat_mg = as.data.frame (t(mMat))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
mMat_mg = t (mMat_mg)

ext_vs_ctx_mean = rowMeans (mMat_mg[,c('CD8_exhausted','KLRC1_NK')])
ext_vs_ctx_mean = ext_vs_ctx_mean[order (-ext_vs_ctx_mean)]
write.csv (ext_vs_ctx_mean, 'regTNK1_module.csv')


ps = log2(as.data.frame (AverageExpression (srt, features = rownames(mMat_mg), group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
#ps = ps[, colnames(DAM_hm@matrix)]
ps_tf = ps[tf_ext,c('CD8_exhausted','KLRC1_NK')]
ps_tf_mean = rowMeans(ps_tf)
top_exp_tf = head(names(ps_tf_mean)[order(-ps_tf_mean)],10)
mMat_mg = mMat_mg[tf_ext,c('CD8_exhausted','KLRC1_NK')]


cd8_bar = HeatmapAnnotation (mark = anno_mark(at = match(top_exp_tf,names(ext_vs_ctx_mean)), 
  labels_gp=gpar(fontsize = 6), 
  labels_rot = 45,
  labels = top_exp_tf,side='top'),  
  ' ' = anno_barplot (ps_tf[names(ext_vs_ctx_mean),c('CD8_exhausted')],border=F,gp = gpar(color = "white")),
#  ' ' = nk_cd8_ext_cor[rownames(ps_tf), 'cd8'],
  #which='row', 
  col=list(' ' = palette_expression_cor_fun),
  simple_anno_size = unit(.2, "cm"))
nk_bar = HeatmapAnnotation (
  #'' = nk_cd8_ext_cor[rownames(ps_tf), 'nk'],
  ' ' = anno_barplot(-ps_tf[names(ext_vs_ctx_mean),c('KLRC1_NK')],
    border=F,gp = gpar(color = "white")), 
  #which='row',
   col=list(' ' = palette_expression_cor_fun),
    simple_anno_size = unit(.2, "cm"))

TF_hm = Heatmap (t(mMat_mg[names(ext_vs_ctx_mean),]), 
          #row_labels = colnames (mMat_mg),
          #column_title = paste('top',top_genes),
          #clustering_distance_columns = 'euclidean',
          #clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          top_annotation = cd8_bar,
          bottom_annotation = nk_bar,
          row_names_side = 'left',
          #column_km = 20,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 0),
          column_names_rot = 45,
          name = 'TF activity',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = palette_deviation_fun (mMat_mg)#,
          #width = unit(1, "cm")
          #right_annotation = motif_ha
          )

pdf (file.path ('Plots','chromvar_rna_expression_KLRC1_NK_CD8_EX_heatmap.pdf'), width = 5,height=2.5)
draw (TF_hm)# + TF_exp_selected_hm2)
dev.off()
   
# Export table
all (rownames(mMat_mg) == rownames (ps_tf))
mat_combined = cbind(mMat_mg, ps_tf)
colnames (mat_combined) = c('CD_Ex_activity','KLRC1_NK_activity','CD8_Ex_RNA','KLRC1_NK_RNA')
write.csv (mat_combined, 'top_TF_CD8_NK_dual_dysfunc_TF_activity_RNA.csv')



### Plot deviation and expression of NR4A2 across samples ####
archp$sample_celltype_lv2 = paste0(archp$Sample,'|',archp$celltype_lv2)
metaGroupName = 'sample_celltype_lv2'
TF = 'IRF4'
TF = 'EOMES'
TF = 'CBFB'
TF = 'NR4A2'
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat_mg = mMat[TF, , drop=F]
mMat_mg = as.data.frame (t(mMat_mg))
mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = setNames(mMat_mg[,2], rownames(mMat_mg))

# Generate RNA pseudobulk of matching cell types ####
srt$sample_celltype_lv2 = paste0(srt$sampleID,'|',srt$celltype_lv2)
metaGroupName = 'sample_celltype_lv2'
ps = log2(as.data.frame (AverageExpression (srt, features = TF, group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
mMat_mg = mMat_mg[sapply (names(mMat_mg), function(x) unlist(strsplit(x, '\\|'))[1] %in% unique (srt$sampleID))]
ps_df = data.frame (deviation = mMat_mg, expression = unlist(ps[names(mMat_mg)]))
ps_df$celltype = sapply (rownames(ps_df), function(x) unlist (strsplit (x, '\\|'))[2])
ps_df$sample = sapply (rownames(ps_df), function(x) unlist (strsplit (x, '\\|'))[1])
ps_df1 = ps_df[, c('deviation','celltype','sample')]
colnames (ps_df1)[1] = 'intensity'
ps_df1$type = 'activity'
ps_df2 = ps_df[, c('expression','celltype','sample')]
colnames (ps_df2)[1] = 'intensity'
ps_df2$type = 'expression'
ps_df3 = rbind (ps_df1, ps_df2)
ps_df3$celltype = factor (ps_df3$celltype, levels = c('KLRC1_NK','CD8_exhausted','Tregs','CD8','CD4','FGFBP2_NK'))
ps_sp = ggplot (ps_df3, aes (x = celltype, y= intensity, fill= celltype)) + 
  geom_point (aes(x = celltype, y = intensity), 
    color = 'grey44', size=.5,
    pch = 19, alpha = 0.5 , position = 'jitter') +
  bxp + 
  scale_fill_manual (values = palette_tnk_cells) + 
  facet_wrap (~type, scales = 'free_y') +
  ggtitle (TF) +
  gtheme #+ geom_text_repel (size=3.5, data = ps_df, aes(label = rownames(ps_df))) 

pdf (file.path ('Plots',paste0(TF,'_dev_exp_boxplot2.pdf')), width=5,height=2.5)
ps_sp
dev.off()


# Compare NKT pseudobulks to peakset of exhausted CD8 and from human meta-analysis Riegel et al ####
if (!exists('fragments')) 
  {
  fragments = getFragmentsFromProject (archp)
  fragments = unlist(fragments)
  }

# CD8 ext DAP ####
DAP_res = readRDS ('DAP_CD8_CD8_ext_pairwise.rds')
DAP_res = do.call (cbind, (assays(DAP_list)))
colnames (DAP_res) = names(assays(DAP_list))
DAP_res_regions = makeGRangesFromDataFrame(rowData(DAP_list)[,c(1,3,4)])
rownames(DAP_res) = as.character(DAP_res_regions)
ext_ps_gr = GRanges (rownames(DAP_res[DAP_res$FDR < .01 & DAP_res$Log2FC > 0, ]))
ext_ps_gr = extendGR (gr = ext_ps_gr, upstream = 3000, downstream = 3000)
peak_windows = slidingWindows (x = ext_ps_gr, width = 50, step = 25)

# Create peak enrichment plots with CD8 ext peaks ####
metaGroupName = 'celltype_lv2'
force = FALSE
ext_l = list()

pb =progress::progress_bar$new(total = length (unique (as.character(archp@cellColData[,metaGroupName]))))
if (!file.exists(paste0('ext_peaks_',metaGroupName,'.rds')) | force)
  {
  for (metagroup in unique (as.character(archp@cellColData[,metaGroupName])))
    {
    pb$tick()            
    #fragments = ReadFragments(fragment_paths[sam], cutSite = FALSE)
    fragments_metagroup = fragments[fragments$RG %in% rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) == metagroup]]
    fragments_metagroup_counts = lapply (peak_windows, function(x) countOverlaps (x, fragments_metagroup))
    fragments_metagroup_counts_df = do.call (cbind, fragments_metagroup_counts)
    fragments_metagroup_counts_df = (fragments_metagroup_counts_df / sum (archp$ReadsInTSS[as.logical(archp@cellColData[,metaGroupName] == metagroup)])) * 10e6
    #colnames (fragments_metagroup_counts_df) = as.character(ext_hg38_ov)
    # write.table (fragments_metagroup_counts_df, file.path(paste0('riegel_ext_peaks_windows_',metagroup,'.tsv')), sep='\t')
    ext_l[[metagroup]] = fragments_metagroup_counts_df
    }
  saveRDS (ext_l, paste0('ext_peaks_',metaGroupName,'.rds'))
  } else {
  ext_l = readRDS (paste0('ext_peaks_',metaGroupName,'.rds'))
  }

#ext_df_long = gather (ext_df, coverage, celltype, 1:(ncol(ext_df)-1))
ext_den = lapply (ext_l, function(x) rowSums(x))
ext_den = as.data.frame (do.call (cbind, ext_den))
#ext_den = ext_den[,1,drop=F]
#ext_den = gather (as.data.frame(ext_den), celltype, coverage)
ext_den$bin = 1:260
#ext_den = split (ext_den, ext_den$celltype)
den = lapply (c('FGFBP2_NK','CD4','CD8','KLRC1_NK','CD8_exhausted','Tregs'), function(x) ggplot (ext_den[,c(x,'bin')], aes_string (x= 'bin', y = x)) + 
#geom_bar (stat= 'identity')
geom_density(color="navyblue", fill="navyblue", alpha=0.5, stat='identity') +
ylim (c(0, max(ext_den))) +
gtheme_no_rot)
#facet_wrap (~celltype, ncol=length(unique(ext_den$celltype)))

pdf (file.path ('Plots',paste0('density_coverage_',metaGroupName,'_CD8_ext.pdf')), height=4,width=12)
wrap_plots (den, ncol= length(unique(names (ext_l))))
dev.off()

median_celltype = lapply (ext_l, function(x) colMeans (x))
median_celltype_order = sapply (median_celltype, mean)
median_celltype_order = names(median_celltype_order)[order(median_celltype_order)]
ext_l2 = lapply (seq_along(ext_l), function(x) ext_l[[x]][,order(-median_celltype[[x]])])
ext_df = do.call (rbind, ext_l2)
ext_df = apply (ext_df, c(1,2), function(x) as.numeric(x))
ext_df = as.data.frame (t(log10(ext_df+1)))
colnames (ext_df) = rep (names(ext_l), each=dim(ext_l[[1]])[1])
#rownames(ext_df) = as.character(ext_hg38_sub)
#ext_df$region = rownames(ext_df)
#rownames (ext_df) = NULL
#ext_df = as.data.frame(apply(ext_df, 2, function(x) sort(x, decreasing = TRUE)))
#palette_fragments = paletteer::paletteer_c("ggthemes::Classic Area-Brown",n=40)
ha = HeatmapAnnotation (bar1 = anno_barplot(colSums(ext_df),gp = gpar(fill = "azure4",border =NA,lty='blank'),border =FALSE, baseline=200,lty='blank'))
ext_df2 = ext_df

pdf (file.path ('Plots','CD8_ext_peakset_enrichments.pdf'), height=4,width=6)
Heatmap (ext_df,
  top_annotation = ha,
#  column_split = , 
  column_split=factor(colnames(ext_df), levels=median_celltype_order),
  cluster_rows=F,
  column_title_gp = gpar(
fontsize = 8),
  column_names_gp = gpar(
fontsize = 0),
  row_names_gp = gpar(
fontsize = 0),
  col = rev(palette_fragments), 
  cluster_columns=F,
  border=T)
dev.off()


# Riegel CD8 Ex peaks ####
library (liftOver)
ext_ps = read.csv ('Riegel_CD8ex_DAP.csv')
ext_clusters = c('C1','C2','C3','C4')
ext_ps = ext_ps[ext_ps$group_name %in% ext_clusters,]
ext_ps_gr = GRanges (ext_ps)
if(!file.exists('hg19ToHg38.over.chain'))
  {
  download.file("https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz", "hg19ToHg38.over.chain.gz")
  system("gzip -d hg19ToHg38.over.chain.gz")
  #system (paste0('wget (https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg38.over.chain.gz)', '-P', getwd()))
  }
ch = import.chain ('hg19ToHg38.over.chain')
#seqlevelsStyle(hub_hg19) = "UCSC"  # necessary
ext_hg38 = unlist (liftOver(ext_ps_gr, ch))
ext_hg38_sub = head (ext_hg38[ext_hg38$Log2FC > 3, ],1000)
ext_hg38_sub = resize (ext_hg38_sub, 1, "center")
ext_hg38_sub = extendGR (gr = ext_hg38_sub, upstream = 3000, downstream = 3000)
peak_windows_riegel = slidingWindows (x = ext_hg38_sub, width = 50, step = 25)
# peak_windows = peak_windows[sapply(peak_windows, length) == 50]
# peak_windows_l =  lapply (1:length(peak_windows[[1]]), function(x) unlist(lapply (peak_windows, function(y) y[x])))

# Create peak enrichment plots with Riegel peaks ####
metaGroupName = 'celltype_lv2'
force = FALSE
ext_l = list()

pb =progress::progress_bar$new(total = length (unique (as.character(archp@cellColData[,metaGroupName]))))
if (!file.exists(paste0('Riegel_CD8ex_DAP_',metaGroupName,'.rds')) | force)
  {
  for (metagroup in unique (as.character(archp@cellColData[,metaGroupName])))
    {
    pb$tick()            
    #fragments = ReadFragments(fragment_paths[sam], cutSite = FALSE)
    fragments_metagroup = fragments[fragments$RG %in% rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) == metagroup]]
    fragments_metagroup_counts = lapply (peak_windows, function(x) countOverlaps (x, fragments_metagroup))
    fragments_metagroup_counts_df = do.call (cbind, fragments_metagroup_counts)
    fragments_metagroup_counts_df = (fragments_metagroup_counts_df / sum (archp$ReadsInTSS[as.logical(archp@cellColData[,metaGroupName] == metagroup)])) * 10e6
    #colnames (fragments_metagroup_counts_df) = as.character(ext_hg38_ov)
    # write.table (fragments_metagroup_counts_df, file.path(paste0('riegel_ext_peaks_windows_',metagroup,'.tsv')), sep='\t')
    ext_l[[metagroup]] = fragments_metagroup_counts_df
    }
  saveRDS (ext_l, paste0('Riegel_CD8ex_DAP_',metaGroupName,'.rds'))
  } else {
  ext_l = readRDS (paste0('Riegel_CD8ex_DAP_',metaGroupName,'.rds'))
  }

#ext_df_long = gather (ext_df, coverage, celltype, 1:(ncol(ext_df)-1))
ext_den = lapply (ext_l, function(x) rowSums(x))
ext_den = as.data.frame (do.call (cbind, ext_den))
#ext_den = ext_den[,1,drop=F]
#ext_den = gather (as.data.frame(ext_den), celltype, coverage)
ext_den$bin = 1:240
#ext_den = split (ext_den, ext_den$celltype)
den = lapply (c('FGFBP2_NK','CD4','CD8','KLRC1_NK','CD8_exhausted','Tregs'), function(x) ggplot (ext_den[,c(x,'bin')], aes_string (x= 'bin', y = x)) + 
#geom_bar (stat= 'identity')
geom_density(color="navyblue", fill="navyblue", alpha=0.5, stat='identity') +
ylim (c(0, max(ext_den))) +
gtheme_no_rot)
#facet_wrap (~celltype, ncol=length(unique(ext_den$celltype)))

pdf (file.path ('Plots',paste0('density_coverage_',metaGroupName,'2.pdf')), height=4,width=12)
wrap_plots (den, ncol= length(unique(names (ext_l))))
dev.off()

median_celltype = lapply (ext_l, function(x) colMeans (x))
median_celltype_order = sapply (median_celltype, mean)
median_celltype_order = names(median_celltype_order)[order(median_celltype_order)]
ext_l2 = lapply (seq_along(ext_l), function(x) ext_l[[x]][,order(-median_celltype[[x]])])
ext_df = do.call (rbind, ext_l2)
ext_df = apply (ext_df,c(1,2), function(x) as.numeric(x))
ext_df = as.data.frame (t(log10(ext_df+1)))
colnames (ext_df) = rep (names(ext_l), each=240)
rownames(ext_df) = as.character(ext_hg38_sub)
#ext_df$region = rownames(ext_df)
#rownames (ext_df) = NULL
#ext_df = as.data.frame(apply(ext_df, 2, function(x) sort(x, decreasing = TRUE)))
#palette_fragments = paletteer::paletteer_c("ggthemes::Classic Area-Brown",n=40)
ha = HeatmapAnnotation (bar1 = anno_barplot(colSums(ext_df),gp = gpar(fill = "azure4",border =NA,lty='blank'),border =FALSE, baseline=200,lty='blank'))
pdf (file.path ('Plots','Riegel_CD8ex_peakset_enrichments.pdf'), height=4,width=6)
Heatmap (ext_df,
  top_annotation = ha,
#  column_split = , 
  column_split=factor(rep (names(ext_l), each=240), levels=median_celltype_order),
  cluster_rows=F,
  column_title_gp = gpar(
fontsize = 8),
  column_names_gp = gpar(
fontsize = 0),
  row_names_gp = gpar(
fontsize = 0),
  col = rev(palette_fragments), 
  cluster_columns=F,
  border=T)
dev.off()







### Find all peaks around exhausted TFs and correlate accessibility and RNA expression across celltype pseudobulks ####
metaGroupName = 'celltype_lv2'
#TF = read.csv ('exhausted_TF_rna.csv')
TF = c('NR4A2','RUNX2') # include only those that show balanced expression between NK KLRC1 and CD8 ext
# Get all peaks correlated with exhausted TF
library (org.Hs.eg.db)
gene_regions = genes (TxDb.Hsapiens.UCSC.hg38.knownGene)
eg_sym = toTable (org.Hs.egSYMBOL)
gene_regions$symbol = eg_sym$symbol[match(gene_regions$gene_id, eg_sym$gene_id)]
gene_regions = gene_regions[gene_regions$symbol %in% TF ]

extend_region = 250000
gene_regions_extended = GRangesList (lapply (seq_along(gene_regions), function(x) extendGR (gr = gene_regions[x], upstream = extend_region, downstream = extend_region)))

gene_regions_peaks = lapply (seq_along(gene_regions), function(x) queryHits (findOverlaps(getPeakSet(archp), gene_regions_extended[x])))
gene_regions_peaks = lapply (gene_regions_peaks, function(x) getPeakSet(archp)[x])

force = FALSE
if (!file.exists ('pMats.rds') | force)
  {
  pMats = lapply (gene_regions_peaks, function(x) getGroupSE(
    ArchRProj = archp,
    useMatrix = 'PeakMatrix',
    groupBy = metaGroupName,
    divideN = TRUE,
    scaleTo = NULL,
    threads = getArchRThreads(),
    verbose = TRUE,
    logFile = createLogFile("getGroupSE")
  ))
  pMats2 = lapply (seq_along(gene_regions_peaks), function(x) pMats[[x]][queryHits(findOverlaps(GRanges(rowData(pMats[[x]])), gene_regions_peaks[[x]]))])
  pMats2 = lapply (pMats2, function(x) {tmp = assay(x); rownames(tmp) = as.character(GRanges (rowData(x))); tmp})
  names (pMats2) = gene_regions$symbol
  saveRDS (pMats2, 'pMats2.rds')
  } else {
  pMats2 = readRDS ('pMats2.rds')  
  }


ps = log2(as.data.frame (AverageExpression (srt, features = gene_regions$symbol, group.by = metaGroupName)[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))


# Compute distance between TSS and correlated peaks ####
library (rtracklayer)
library (AnnotationDbi)
gene_ids = toTable (org.Hs.egSYMBOL)$gene_id[match (TF,toTable (org.Hs.egSYMBOL)$symbol)]
mygenes.transcripts = subset (genes(TxDb.Hsapiens.UCSC.hg38.knownGene, columns=c("tx_id", "tx_name","gene_id")), gene_id %in% gene_ids)
mygenes.tss = resize (mygenes.transcripts, width=1, fix='start')

mygenes.tss$symbol = toTable (org.Hs.egSYMBOL)$symbol[match (mygenes.tss$gene_id ,toTable (org.Hs.egSYMBOL)$gene_id)]

peak_gene_cor = lapply (gene_regions$symbol, function(x) cor (t(pMats2[[x]]), t(ps[x, colnames(pMats2[[x]])])))
peak_gene_cor = lapply (peak_gene_cor, function(x) x[order(-x[,1]),,drop=F])
names (peak_gene_cor) = gene_regions$symbol

top_e = 5
peak_gene_cor = lapply (names(peak_gene_cor), function(x) {
  peak_gene_cor[[x]] = as.data.frame (peak_gene_cor[[x]])
  gr = GRanges (rownames(peak_gene_cor[[x]]))
  dis = distanceToNearest (gr,mygenes.tss[mygenes.tss$symbol == x])@elementMetadata$distance
  pr = follow (gr,mygenes.tss[mygenes.tss$symbol == x])
  pr[is.na(pr)] = 0
  dis = ifelse (pr == 1, dis * -1, dis) 
  peak_gene_cor[[x]]$distance = dis
  peak_gene_cor[[x]]$name = ''
  peak_gene_cor[[x]]$name[1:top_e] = head (paste0('E',peak_gene_cor[[x]]$distance),top_e)
  peak_gene_cor[[x]]$max = apply(pMats2[[x]],1, max)[match(rownames(peak_gene_cor[[x]]),rownames(pMats2[[x]]))]
  colnames(peak_gene_cor[[x]]) = c('cor','distance','name','max')
  peak_gene_cor[[x]]
  })

names (peak_gene_cor) = gene_regions$symbol
vp = lapply (names(peak_gene_cor), function(x)
    {
    ggplot (peak_gene_cor[[x]], aes(x=distance, y= cor)) +
    geom_bar (stat = 'identity') +
    geom_point (aes (size = max), shape=21, color='white',alpha=.5, fill = 'red') +
    xlim(c(-extend_region,extend_region)) + 
    ggtitle (x) + 
    geom_label_repel (size=3, color='grey22', aes(label = name), segment.size=.2) + gtheme_no_rot
    })

pdf (file.path ('Plots','enhancers_distance_ext_TF.pdf'),width=7,height=3)
vp
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
    
### Plot NR4A2 region ####
archp$celltype3 =  archp$celltype_lv2
archp$celltype3[archp$celltype3 == 'CD4'] =  c('C1_CD4')
archp$celltype3[archp$celltype3 == 'CD8'] =  c('C2_CD8')
archp$celltype3[archp$celltype3 == 'NK_FGFBP2'] =  c('C3_FGFBP2_NK')
archp$celltype3[archp$celltype3 == 'Tregs'] =  c('C4_Tregs')
archp$celltype3[archp$celltype3 == 'CD8_exhausted'] =  c('C5_CD8_exhausted')
archp$celltype3[archp$celltype3 == 'NK_KLRC1'] =  c('C6_KLRC1_NK')
palette_tnk_cells_ext2 = palette_tnk_cells
names (palette_tnk_cells_ext2) = c(
  'C2_CD8',
  'C1_CD4',
  'C4_Tregs',
  'C6_KLRC1_NK',
  'C3_FGFBP2_NK',
  'C5_CD8_exhausted')


### Hubs analysis #####
metaGroupName = "Clusters_H"
cor_cutoff = 0.2
#max_dist = 12500
max_dist = 12500
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
hubs_obj = readRDS (file.path(hubs_dir,'global_hubs_obj.rds'))  

metaGroupName = 'celltype3'
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
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
    hubs_regions = hubs_obj$hubsCollapsed,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 50000,
    pal = palette_tnk_cells_ext2,
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

metaGroupName = 'celltype_lv2'
top_dah = data.frame (
gene = srt@assays$RNA$data[TF,],
group = srt@meta.data[,metaGroupName])
top_dah$group = factor (top_dah$group, levels = 
  rev (c('CD4','CD8','FGFBP2_NK','Tregs','CD8_exhausted',
  'KLRC1_NK','TFH')))
top_dah = na.omit(top_dah)
bp = ggplot (top_dah, aes (x = gene, y = group, fill = group)) + 
vlp + 
bxpv + 
scale_fill_manual (values = palette_tnk_cells) +
#geom_point (position='identity', alpha=.3, color="grey44", size=1) +
gtheme_no_rot

pdf (file.path ('Plots', paste0('scrna_',TF,'_boxplots.pdf')), height=4, width=4)
bp
dev.off()



metaGroupName = 'celltype3'
TF = c('CTLA4','PDCD1','ICOS','HAVCR2')
pdf()
meso_markers <- plotBrowserTrack2 (
    ArchRProj = archp, 
    sizes = c(6, 1, 1, 1,1,1),
    groupBy = metaGroupName, 
    geneSymbol = TF,
    normMethod = "ReadsInTSS",
    scCellsMax=3000,
    loop_size = .2,
        plotSummary = c("bulkTrack", "featureTrack", 
        "loopTrack","geneTrack", 
        "hubTrack",'hubregiontrack'),
        hubs_regions = hubs_obj$hubsCollapsed,
    #pal = DiscretePalette (length (unique(sgn2@meta.data[,metaGroupName])), palette = 'stepped'), 
    #region = ext_range (GRanges (DAH_df$region[22]),1000,1000),
    upstream = 100000,
    pal = palette_tnk_cells_ext2,
    #ylim=c(0,0.1),
    downstream = 100000,
    #loops = getPeak2GeneLinks (archp, corCutOff = 0.2),
    #pal = ifelse(grepl('T',unique (archp2@cellColData[,metaGroupName])),'yellowgreen','midnightblue'),
#    loops = getCoAccessibility (archp, corCutOff = 0.25),
    #  returnLoops = TRUE),
    useGroups= NULL,
    loops = getPeak2GeneLinks (archp, corCutOff = 0.2,returnLoops = TRUE),
    hubs = hubs_obj$peakLinks2
)
dev.off()
plotPDF (meso_markers, ArchRProj = archp, 
  width=5,height=3, 
  name =paste0(paste(TF, collapse='_'),'_coveragePlots.pdf'),
  addDOC = F)

  
TF='ICOS'
metaGroupName = 'celltype_lv2'
top_dah = data.frame (
gene = srt@assays$RNA@data[TF,],
group = srt@meta.data[,metaGroupName])
top_dah$group = factor (top_dah$group, levels = 
  rev (c('CD4','CD8','FGFBP2_NK','Tregs','CD8_exhausted',
  'KLRC1_NK','TFH')))
top_dah = na.omit(top_dah)
bp = ggplot (top_dah, aes (x = gene, y = group, fill = group)) + 
vlp + 
bxpv + 
scale_fill_manual (values = palette_tnk_cells) +
#geom_point (position='identity', alpha=.3, color="grey44", size=1) +
gtheme_no_rot

pdf (file.path ('Plots', paste0('scrna_',TF,'_boxplots.pdf')), height=4, width=4)
bp
dev.off()






# Check footprint of RUNX and NR4A2 across celltypes ####
metaGroupName='celltype_lv2'
archp <- addGroupCoverages (ArchRProj = archp, groupBy = metaGroupName)
motifPositions <- getPositions (archp)

motifs <- c('RUNX1','RUNX2','RUNX3','NR4A2')
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

seFoot <- getFootprints(
  ArchRProj = archp, 
  #positions = motifPositions_sample[markerMotifs], 
  positions = motifPositions[markerMotifs], 
  groupBy = metaGroupName
)
  
plotFootprints(
seFoot = seFoot,
ArchRProj = archp, 
normMethod = "Subtract",
plotName = "Footprints-Subtract-Bias_",
addDOC = FALSE, height=7.5, width=5,
pal = palette_tnk_cells,
smoothWindow = 25)
  














### Show sample distribution of average chromatin accessibility of exhausted peaks 
#### Show sample distribution of NR4A2 enhancer ####
archp$celltype_sample = paste0(archp$celltype_lv2, '_', archp$Sample)
metaGroupName = 'celltype_sample'  
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
#promoter_region = getPeakSet (archp)
#promoter_region =  promoter_region[which(promoter_region$peakType == 'Promoter' & promoter_region$nearestGene == 'NR4A2')]
#promoter_region = GRanges (seqnames = as.character(seqnames(promoter_region))[1], ranges= IRanges(start = min(start(promoter_region)), end = max(end(promoter_region))))

pmat_enhancer = unlist(as.data.frame (assay (pMats[queryHits(findOverlaps (GRanges(rowData(pMats)), enhancer_region)),])))
pmat_enhancer_df = data.frame (region = pmat_enhancer, celltype =names(pmat_enhancer))
pmat_enhancer_df$celltype = gsub ('_P.*','',pmat_enhancer_df$celltype)
pmat_enhancer_df$sample = sub('.*(P.*)', '\\1', rownames(pmat_enhancer_df))
pmat_enhancer_df = pmat_enhancer_df[pmat_enhancer_df$sample %in% names (table (pmat_enhancer_df$sample)[table (pmat_enhancer_df$sample) ==6]), ]
pmat_enhancer_df$celltype = factor(pmat_enhancer_df$celltype, levels = c('KLRC1_NK','CD8_exhausted','CD8','FGFBP2_NK','CD4','Tregs'))

bp = ggplot (pmat_enhancer_df, aes (x= celltype, y= region)) +
      #geom_violin (trim=TRUE, aes (fill = treatment), alpha=.6) +
      #geom_violin (aes_string(fill = metaGroupNames[3])) +
      geom_point (aes (x = celltype, y = region), position='identity', alpha=.7, color="blue", size=1.2) +
      geom_boxplot(width=0.5, aes (fill = celltype), color='grey22', alpha=.6) +
      scale_fill_manual (values= palette_tnk_cells) + 
      #scale_color_manual (values= palette_tnk_cells) + 
      geom_line (data = pmat_enhancer_df, aes(x = celltype, y = region, group = sample), color='grey22',linewidth=.1, alpha=.5) +
      gtheme

stat.test2 <- bp$data %>%
  filter(celltype %in% c("KLRC1_NK", "CD8_exhausted") | TRUE) %>%   # keep all
  rstatix::pairwise_t_test(
    region ~ celltype, 
    paired = TRUE, 
    p.adjust.method = "fdr"
  ) %>%
  filter(group1 %in% c("KLRC1_NK", "CD8_exhausted") | 
         group2 %in% c("KLRC1_NK", "CD8_exhausted"))

stat.test2 = stat.test2 [grepl ('KLRC1_NK',stat.test2$group1) | grepl ('CD8_exhausted', stat.test2$group1),]
stat.test2 = stat.test2 %>% add_xy_position (x = "celltype", step.increase=0.05)
bp = bp + stat_pvalue_manual (stat.test2, remove.bracket=FALSE,
   bracket.nudge.y = 0, hide.ns = TRUE, position = position_nudge (y = -0.01),
    label = "p.adj.signif") 

pdf (file.path ('Plots','NR4A2_enhancer_sample_boxplots.pdf'),width=5, height=4)
wrap_plots (bp)
dev.off()


#pmat_df$type = factor (pmat_df$type, levels =c ('promoter','enhancer'))
pdf (file.path ('Plots','eNR4F2_accessibility_sample_enhancer_boxplot.pdf'), height=3, width=6)
ep = ggplot (pmat_enhancer_df, aes (x = celltype, y = region, fill = celltype)) + 
geom_boxplot () + gtheme +
scale_fill_manual (values = c(enhancer = '#C1D32FFF', promoter = 'grey'))
#+ scale_fill_manual (values = c(palette_tnk_cells, palette_myeloid, palette_celltype_lv1))
ep
dev.off()


