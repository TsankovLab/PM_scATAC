
# # Generate heatmap of clusters x samples across active TFs ####
# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# all (colnames(mSE) == rownames(archp))
# mMat = assays (mSE)[[1]] 
# rownames (mMat) = rowData (mSE)$name

# # Filter by RNA expression ####
# metaGroupName = 'sampleID'
# active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)
# mMat = mMat[active_TFs, ]

# mMat_agg = aggregate (as.matrix(t(scale(mMat))), by = list(archp$Clusters), FUN='mean')
# rownames (mMat_agg) = mMat_agg[,1]
# mMat_agg = mMat_agg[,-1]

# remove_clusters = 'C9' # remove normal sample
# ha = HeatmapAnnotation (sample = rownames(mMat_agg), col = list(sample = palette_sample))
# mMat_hm = Heatmap (t(mMat_agg[!rownames(mMat_agg) %in% remove_clusters,]),# row_km=15,
#   #left_annotation = ha,
#   #rect_gp = gpar(type = "none"),
#   clustering_distance_rows='euclidean' ,
#   clustering_distance_columns = 'euclidean', 
#   col=palette_deviation_cor_fun, 
#   #row_split = km$cluster,
#   #column_split = km$cluster,
#   #row_km=2, 
#   #column_km=2,
# #  top_annotation = ha,
#   border=T,
# #   ,
#   row_names_gp = gpar(fontsize = 0), 
#   column_names_gp = gpar(fontsize = 8)
# )
# pdf (file.path ('Plots','TF_samples_heatmap.pdf'), width = 4,height=8)
# mMat_hm
# dev.off()



# # Try with PCA ####
# library (uwot)
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat = as.matrix(mMat)#[selected_TF,])
# tf_mat = lapply (sams, function(x) t(scale(mMat[,archp$Sample3 == x ])))
# names (tf_mat) = sams
# #tf_mat = do.call (rbind,tf_mat)
# cnmf_mat = as.matrix(archp@cellColData[,grep ('sarcomatoid', colnames(archp@cellColData))])
# cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp$Sample3 == x, ])))
# names(cnmf_mat) = sams

# # Filter by RNA expression ####
# metaGroupName = 'sampleID'
# active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)
# mMat = mMat[active_TFs, ]

# archp$Clusters_sample = paste0(archp$Clusters, '_',archp$Sample2)
# remove_clusters = names(table (archp$Clusters_sample)[table (archp$Clusters_sample) < 10])

# mMat_agg = aggregate (as.matrix(t(scale(mMat))), by = list(archp$Clusters_sample), FUN='mean')
# rownames (mMat_agg) = mMat_agg[,1]
# mMat_agg = mMat_agg[,-1]
# mMat_agg = mMat_agg[!rownames(mMat_agg) %in% remove_clusters,]


# pca_result <- prcomp(mMat_agg, center = TRUE, scale = TRUE)

# tmp = as.data.frame(pca_result$x)[,c('PC1','PC2')]
# umap_result = tmp
# umap_result$Sample = sapply (rownames(umap_result), function(x) unlist(strsplit(x,'_'))[2])
# umap_result = umap_result[umap_result$Sample != 'normal', ]

# pca_data = ggplot(umap_result, aes(x = PC1, y = PC2, fill = Sample)) +#, colors=class))+#, color = sarc)) +
# geom_point(alpha=0.8,size=4, shape=21, color='grey22') + gtheme_no_rot +
# #geom_point(data = umap_result[umap_result$class != 'ND',], aes(color=class),alpha=0.5) +
# #geom_smooth(data = umap_result[umap_result$class != 'ND',],
# #  method = "lm", se = TRUE, aes (color = class2)) + 
# #geom_point(data = umap_result[umap_result$tf != '',],size=4) +
# #geom_text_repel(data= umap_result[umap_result$tf != '',]) + 
# #labs(title = paste("UMAP - ",sam)) #+ 
# scale_fill_manual (values = palette_sample)

# pdf (file.path('Plots','clusters_TF_umap.pdf'), width = 3.5, height = 3)
# pca_data
# dev.off()












### Run TF correlation to identify TF modules across cancers #### 
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat)#[selected_TF,])

# Filter by RNA expression ####
metaGroupName = 'sampleID'
active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)
mMat = mMat[active_TFs, ]

archp_NN = archp[archp$Sample2 != 'normal_pleura']
mMat_NN = mMat[,archp$Sample2 != 'normal_pleura']
Clusters_sample = paste0(archp_NN$Clusters, '_',archp_NN$Sample2)
remove_low_clusters = !Clusters_sample %in% names(table (Clusters_sample)[table (Clusters_sample) < 10])
Clusters_sample = Clusters_sample[remove_low_clusters]
archp_NN = archp_NN[remove_low_clusters]
mMat_NN = mMat_NN[,remove_low_clusters]
mMat_cor = cor (as.matrix(t(scale(mMat_NN))), method = 'spearman')

set.seed(1234)
centers=5
km = kmeans (mMat_cor, centers=centers)

pdf (file.path ('Plots','TF_modules_heatmap.pdf'), width = 4,height=3)
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
# cell_fun = function(j, i, x, y, w, h, fill) {# THIS DOESNT WORK NEED TO USE LAYER_FUN
#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#         }}
  ))
  # ,
  # cell_fun = function(j, i, x, y, w, h, fill) {
  #       if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
  #           grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#        }}))
dev.off()

pdf (file.path ('Plots','TF_modules_heatmap.pdf'), width = 4,height=3)
cor_mMat_hm
dev.off()

tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat_NN[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = do.call (cbind, tf_modules)
archp_NN@cellColData = archp_NN@cellColData[!colnames(archp_NN@cellColData) %in% paste0('mod_',unique(km$cluster))]
archp_NN@cellColData = cbind (archp_NN@cellColData, tf_modules) 

pdf()
TF_p = plotEmbedding (
    ArchRProj = archp_NN,
    colorBy = "cellColData",
    name = paste0('mod_',unique(km$cluster)), 
    pal = rev(palette_deviation),
    #useSeqnames='z',
    embedding = "UMAP")
dev.off()
pdf (file.path ('Plots','TF_modules_umap.pdf'), width = 20,height=6)
wrap_plots (TF_p, ncol=5)
dev.off()


# ridge plots of TF modules ####
library (ggridges)
library (ggplot2)
library (viridis)
#library(hrbrthemes)
tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat_NN[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = as.data.frame (do.call (cbind, tf_modules))
tf_modules$Clusters_Sample = Clusters_sample
tf_modules$Sample = sapply (tf_modules$Clusters_Sample, function(x) unlist(strsplit(x, '_'))[2])
palette_clusters_sample = setNames (palette_sample[tf_modules$Sample], tf_modules$Clusters_Sample)
palette_clusters_sample = palette_clusters_sample[!duplicated(palette_clusters_sample)]
# sapply (rownames(tf_modules), function(x) unlist(strsplit (x, '\\#'))[1]) == archp$Sample2
tf_modules = gather (tf_modules, module, expression,1:centers)
tf_modules$module = factor (tf_modules$module, levels = paste0('mod_',names (row_order (cor_mMat_hm))))

# Plot
# rp = lapply (paste0('mod_',unique(km$cluster)), function(x) 
#   ggplot(tf_modules,  aes_string(x='Sample', y=x, fill='..x..')) +

#   geom_vridgeline(stat="ydensity", trim=FALSE, alpha = 0.85, scale = 2, width=.10) +
#   palette_deviation_ggplot_fill +
#     theme_classic())
# pdf (file.path ('Plots','TF_modules_ridge_plots.pdf'), width = 20,height=3)
# wrap_plots (rp, ncol=5)
# dev.off()

dp = ggplot (tf_modules) +
  geom_density(aes(x=expression,fill=Clusters_Sample),color='white',
                      alpha = 0.6) +
  # geom_vline(aes(xintercept = mean, group = tf_modules, linetype = Sample),
  #            data = combined_sla_means) +
  facet_wrap (~module, nrow = 5, scales = 'free',strip.position = "left") +
  scale_fill_manual (values = palette_clusters_sample) +
  gtheme_no_rot

pdf (file.path ('Plots','TF_modules_ridge_plots2.pdf'), width = 5,height=8)
dp
dev.off()





# Read in peak files from scATAC studies ####
projects = 'rawlins_fetal_lung'
projects_peaks = lapply (seq_along(projects), function(x) {
  bed_files = list.files (file.path('..','all_tissues_ArchR',projects[x],'PeakCalls'), pattern = '.rds')
  grlist = lapply (seq_along(bed_files), 
    function(y) readRDS (file.path('..','all_tissues_ArchR',projects[x],'PeakCalls',bed_files[y])))
names (grlist) = paste0(projects[x], '_', sapply (bed_files, function(z) unlist(strsplit (z, '-'))[1]))
grlist
})
projects_peaks = unlist (projects_peaks, recursive=F)

#### chromVAR analysis ####
archp = addBgdPeaks (archp, force= TRUE)
archp = addPeakAnnotations (ArchRProj = archp, 
     regions = projects_peaks, name = "rawlins_fetal_lung")

archp = addDeviationsMatrix (
  ArchRProj = archp, 
  peakAnnotation = "scATAC_datasets",
  force = TRUE
)
