
# Read in peak files from scATAC studies ####
projects = c('Tsankov_lung','rawlins_fetal_lung')
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
archp = addPeakAnnotations (ArchRProj = archp, force=T,
     regions = projects_peaks, name = "scatac_studies")

archp = addDeviationsMatrix (
  ArchRProj = archp, 
  peakAnnotation = "scatac_studies",
  force = TRUE
)

# Correlate scatac celltype dev with motif dev only HOXB13 vs fetal meso ####
if(!exists ('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = scale(as.matrix (assays (mSE)[[1]]))
rownames(mMat) = rowData (mSE)$name
if(!exists ('feSE')) feSE = fetch_mat (archp, 'scatac_studies')
feMat = scale(as.matrix (assays (feSE)[[1]]))

all (colnames (feMat) == colnames (mMat))
all (colnames(feMat) == rownames(archp@cellColData))
#scatac_celltype = 'bingren_pan_Astrocyte.2'
fetal_celltype = c('rawlins_fetal_lung_Earlymeso',
    'rawlins_fetal_lung_Mid_latemeso',
    'Tsankov_lung_Mesothelium'
    )
fetal_celltype = rownames (feMat)
#scatac_celltype = 'rawlins_fetal_lung_Mid_latemeso'
metaGroupName = 'Sample3'
scatac_tf_cor = lapply (unique (as.character(archp@cellColData[,metaGroupName])), 
  function(x) 
  cor (t(feMat)[as.character(archp@cellColData[,metaGroupName]) == x,fetal_celltype], 
    t(mMat)[as.character(archp@cellColData[,metaGroupName]) == x,grep ('HOXB13',rownames(mMat)), drop=F]))
names (scatac_tf_cor) = unique (as.character(archp@cellColData[,metaGroupName]))
scatac_tf_cor = do.call (cbind, scatac_tf_cor)
colnames (scatac_tf_cor) = unique (as.character(archp@cellColData[,metaGroupName]))
  
# Plot only fetal mesothelium ####
ht = Heatmap (
  # prop_sample_df, 
  t(scatac_tf_cor),
  cluster_rows=T,
  cluster_columns=T,
  col = rev(palette_deviation), 
  row_names_gp= gpar (fontsize=7), 
  column_names_gp= gpar (fontsize=7),
  name = 'cor',
  border=T#, 
  #column_names_rot = 45
  )
pdf (file.path ('Plots', 'fetal_HOXB13_chromvar_cor_heatmap.pdf'),width = 12,height=5.2)
ht
dev.off()

# Show as scatterplot for P11 HOX+ ####
library (ggpubr)
library (ggpointdensity)
p11_hox_df = as.data.frame(cbind(
  t(feMat[fetal_celltype, archp$Sample3 == 'P11_HOX']),
  t(mMat['HOXB13', archp$Sample3 == 'P11_HOX', drop=F])
  ))
p11_hox_df = gather (as.data.frame(p11_hox_df), fetal_stage, fetal_score, 1:3)
p11_hox_df$fetal_stage = factor (p11_hox_df$fetal_stage, levels = fetal_celltype)
sp = ggplot (p11_hox_df, aes (HOXB13, fetal_score)) + 
  geom_point (color = 'grey22', alpha = 0.5) + 
  #geom_pointdensity () +
  #scale_color_viridis () +
  gtheme_no_rot + 
  facet_wrap (~fetal_stage, ncol=3) +
  geom_smooth(method=lm) + 
  stat_cor(method = "pearson", label.x = .1, label.y = -.35)

pdf (file.path ('Plots','fetal_HOXB13_chromvar_cor_scatter.pdf'), width=7, height=3)
sp
dev.off()



### Check fetal and adult lung pops deviation across tumor cells ####
metaGroupName = 'Sample3_Clusters'
archp$Sample3_Clusters = paste0(archp$Sample3,'_',archp$Clusters)
feMat_agg = as.data.frame (t(feMat))
feMat_agg$metaGroup = as.character(archp@cellColData[,metaGroupName])
feMat_agg = aggregate (.~ metaGroup, feMat_agg, mean)
rownames (feMat_agg) = feMat_agg[,1]
feMat_agg = feMat_agg[,-1]
filter_clusters = names (table (archp$Sample3_Clusters)[table (archp$Sample3_Clusters) > 10])
feMat_agg = feMat_agg[filter_clusters,]
#feMat_agg = t(feMat_agg)
#feMat_agg = feMat_agg[rownames(feMat_agg) %in% rownames(srt),]
  
ht = Heatmap (
  # prop_sample_df, 
  t(feMat_agg),
  cluster_rows=T,
  cluster_columns=T,
  col = palette_deviation_fun(feMat_agg), 
  row_names_gp= gpar (fontsize=7), 
  column_names_gp= gpar (fontsize=7),
  name = 'cor',
  border=T#, 
  #column_names_rot = 45
  )
pdf (file.path ('Plots', 'fetal_adult_dev_clusters_heatmap.pdf'),width = 5,height=12.2)
ht
dev.off()



# Add cancer deviations from bulk ATAC-seq TCGA data ####
archp <- addArchRAnnotations (ArchRProj = archp, collection = "ATAC")
archp <- addDeviationsMatrix(
  ArchRProj = archp, 
  peakAnnotation = "ATAC",
  force = TRUE
)

if(!exists ('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = scale(as.matrix (assays (mSE)[[1]]))
rownames(mMat) = rowData (mSE)$name
if(!exists ('atacSE')) atacSE = fetch_mat (archp, 'ATAC')
atacMat = scale(as.matrix (assays (atacSE)[[1]]))

metaGroupName = 'Sample3_Clusters'
archp$Sample3_Clusters = paste0(archp$Sample3,'_',archp$Clusters)
metaGroupName = 'Sample3'
atacMat_agg = as.data.frame (t(atacMat))
atacMat_agg$metaGroup = as.character(archp@cellColData[,metaGroupName])
atacMat_agg = aggregate (.~ metaGroup, atacMat_agg, mean)
rownames (atacMat_agg) = atacMat_agg[,1]
atacMat_agg = atacMat_agg[,-1]
filter_clusters = names (table (as.character(archp@cellColData[,metaGroupName]))[table (as.character(archp@cellColData[,metaGroupName])) > 10])
atacMat_agg = atacMat_agg[filter_clusters,]
#feMat_agg = t(feMat_agg)
#feMat_agg = feMat_agg[rownames(feMat_agg) %in% rownames(srt),]
  
ht = Heatmap (
  # prop_sample_df, 
  t(atacMat_agg),
  cluster_rows=T,
  cluster_columns=T,
  col = palette_deviation_fun(atacMat_agg), 
  row_names_gp= gpar (fontsize=7), 
  column_names_gp= gpar (fontsize=7),
  name = 'cor',
  border=T#, 
  #column_names_rot = 45
  )
pdf (file.path ('Plots', 'atac_TCGA_dev_clusters_heatmap.pdf'),width = 5,height=12.2)
ht
dev.off()

