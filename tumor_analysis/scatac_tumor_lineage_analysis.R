
# cd /ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/
# wget -r http://catlas.org/catlas_downloads/humanbrain/cCREs/
# wget -r http://catlas.org/catlas_downloads/humantissues/Peaks/

# Read ENCODE collection for peak annotation
# H3K4me3_primary_filepath = '/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/ENCODE/H3K4me3_ChIP-seq/primary_cells/'
# encode_H3K4me3_primary = readRDS (paste0(H3K4me3_primary_filepath, 'H3K4me3_primary_cells_bed_list.rds'))

# Read in peak files from scATAC studies ####
projects = c('yang_kidney','Tsankov_lung','rawlins_fetal_lung','JShendure','greenleaf_colon','greenleaf_brain','bingren_pan')
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
archp = addPeakAnnotations (ArchRProj = archp, 
     regions = projects_peaks, name = "scATAC_fetal_adult")

archp = addDeviationsMatrix (
  ArchRProj = archp, 
  peakAnnotation = "scATAC_fetal_adult",
  force = TRUE
)

if (!any (ls() == 'afSE')) afSE = fetch_mat (archp, 'scATAC_fetal_adult')
afSE = afSE[, archp$cellNames]
afMat = as.matrix (scale (assays (afSE)[[1]]))


# Generate heatmap ####
afMat = as.data.frame (t(afMat))
afMat$Sample = archp$Sample3
afMat = aggregate (.~Sample, data = afMat, FUN = 'mean')
rownames (afMat) = afMat[,1]
afMat = afMat[,-1]
# row_names = unique (as.character(unlist (lapply (prop_sample, function(x) x$Var1))))
# prop_sample = lapply (prop_sample, function(x) {rownames(x) = x$Var1; x})
# prop_sample_df = do.call (cbind, lapply (prop_sample, function(x) x[row_names,]))
# rownames (prop_sample_df) = row_names
# prop_sample_df = prop_sample_df[, grepl ('Freq', colnames(prop_sample_df))]
# colnames (prop_sample_df) = unique(archp$Sample2)

# prop_sample_df[is.na(prop_sample_df)] = 0
ht = Heatmap (
  # prop_sample_df, 
  afMat,
  col = palette_deviation_fun (afMat), 
  row_names_gp= gpar (fontsize=6), 
  column_names_gp= gpar (fontsize=6), 
  column_names_rot = 45)
pdf (file.path ('Plots', 'scatacDatasets_deviations_sample2_heatmap.pdf'),width = 10.8,height=4)
ht
dev.off()




# Correlate scatac celltype dev with motif dev only HOXB13 vs fetal meso ####
if(!exists ('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = as.matrix (scale(assays (mSE)[[1]]))
rownames(mMat) = rowData (mSE)$name
if(!exists ('afSE')) afSE = fetch_mat (archp, 'scATAC_fetal_adult')
afMat = scale(assays (afSE)[[1]])

all (colnames (afMat) == colnames (mMat))
all (colnames(afMat) == rownames(archp@cellColData))
#scatac_celltype = 'bingren_pan_Astrocyte.2'
scatac_celltype = c('rawlins_fetal_lung_Earlymeso',
    'rawlins_fetal_lung_Mid_latemeso',
    'Tsankov_lung_Mesothelium'
    )
#scatac_celltype = 'rawlins_fetal_lung_Mid_latemeso'
metaGroupName = 'Sample3'
scatac_tf_cor = lapply (unique (as.character(archp@cellColData[,metaGroupName])), 
  function(x) 
  cor (t(afMat)[as.character(archp@cellColData[,metaGroupName]) == x,], 
    t(mMat)[as.character(archp@cellColData[,metaGroupName]) == x,grep ('HOXB13',rownames(mMat)), drop=F]))
names (scatac_tf_cor) = unique (as.character(archp@cellColData[,metaGroupName]))
scatac_tf_cor = do.call (cbind, scatac_tf_cor)
colnames (scatac_tf_cor) = unique (as.character(archp@cellColData[,metaGroupName]))

# Plot only fetal mesothelium ####
ht = Heatmap (
  # prop_sample_df, 
  t(scale(scatac_tf_cor)),
  cluster_rows=T,
  cluster_columns=T,
  col = rev(palette_deviation), 
  row_names_gp= gpar (fontsize=7), 
  column_names_gp= gpar (fontsize=7),
  name = 'cor',
  border=T#, 
  #column_names_rot = 45
  )
pdf (file.path ('Plots', 'fetal_HOXB13_chromvar_cor_heatmap2.pdf'),width = 12,height=5.2)
ht
dev.off()







scatac_celltype = 'bingren_pan_Astrocyte.2'
scatac_celltype = 'rawlins_fetal_lung_Earlymeso'
scatac_celltype = 'rawlins_fetal_lung_Mid_latemeso'
markerMotifs = getFeatures (archp, select = paste(scatac_celltype, collapse="|"), useMatrix = "scATAC_datasetsMatrix")
markerMotifs = grep ("z:", markerMotifs, value = TRUE)

archp = addImputeWeights (archp)
TF_p = plotEmbedding (
    ArchRProj = archp, 
    colorBy = "scATAC_datasetsMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    pal = palette_deviation,
    imputeWeights = getImputeWeights(archp)
)

pdf (file.path ('Plots',paste0('scatac_datasets_',scatac_celltype,'_umap.pdf')), width = 4,height=4)
print (wrap_plots (TF_p))
dev.off()


# Try with ridge plots ####
library(ggridges)
library(ggplot2)
library(viridis)
scaMat2 = as.data.frame (t(scaMat))
scaMat2$Sample2 = archp$Sample2
#library(hrbrthemes)
# Plot
scatac_celltype = 'bingren_pan_Astrocyte.2'
# scatac_celltype = 'rawlins_fetal_lung_Earlymeso'
# scatac_celltype = 'rawlins_fetal_lung_Mid_latemeso'
rp = ggplot(scaMat2, aes_string(x = scatac_celltype, y = 'Sample2', fill = '..x..')) +
  geom_density_ridges_gradient (scale = 3, rel_min_height = 0.01, alpha=.5) +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue-White Diverging", direction = -1) +
    theme_classic()
pdf (file.path ('Plots','scatac_celltypes_modules_ridge_plots.pdf'), width = 20,height=3)
wrap_plots (rp, ncol=5)
dev.off()


# Show as scatterplot for P11 HOX+ ####
library (ggpubr)
library (ggpointdensity)
p11_hox_df = as.data.frame(cbind(
  t(scaMat[scatac_celltype, archp$SampleP11 == 'P11_HOX+']),
  t(mMat['HOXB13', archp$SampleP11 == 'P11_HOX+', drop=F])
  ))
p11_hox_df = gather (as.data.frame(p11_hox_df), fetal_stage, fetal_score, 1:3)
p11_hox_df$fetal_stage = factor (p11_hox_df$fetal_stage, levels = scatac_celltype)
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






# Correlate all scatac_celltypes deviations against each other per sample ####
pdf(file.path ('Plots','celltypes_scatac_cor_x_sample_heatmap.pdf'),40,40)
lapply (unique (archp$Sample2), function(x) Heatmap(cor(t(scaMat[, archp$Sample2 == x]))))
dev.off()


### Use hypergeometric enrichment instead of chromvar ####
markersPeaks <- getMarkerFeatures(
    ArchRProj = archp, 
    useMatrix = "PeakMatrix", 
    groupBy = "Sample2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

enrichRegions <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = archp,
    peakAnnotation = "scATAC_datasets",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

heatmapRegions <- plotEnrichHeatmap (enrichRegions, n = 7, transpose = TRUE)


pdf (file.path ('Plots',paste0('scatac_datasets_enrichment_heatmap.pdf')), width = 8,height=7)
heatmapRegions
dev.off()


# Compare peaks overlap between scatac celltypes and meso tumors / normal ####
projects = c('yang_kidney','Tsankov_lung','rawlins_fetal_lung','JShendure','greenleaf_colon','greenleaf_brain','bingren_pan')
projects_peaks = lapply (seq_along(projects), function(x) {
  bed_files = list.files (file.path('..','all_tissues_ArchR',projects[x],'PeakCalls'), pattern = '.rds')
  grlist = lapply (seq_along(bed_files), 
    function(y) readRDS (file.path('..','all_tissues_ArchR',projects[x],'PeakCalls',bed_files[y])))
names (grlist) = paste0(projects[x], '_', sapply (bed_files, function(z) unlist(strsplit (z, '-'))[1]))
grlist
})
projects_peaks = unlist (projects_peaks, recursive=F)
projects_peaks2 = GRangesList (projects_peaks)

meso_peaks = lapply (unique (archp$Clusters), function(x) readRDS (file.path('PeakCalls',paste0(x,'-reproduciblePeaks.gr.rds'))))
names (meso_peaks) = unique (archp$Clusters)
peaks_ov_mat = sapply (meso_peaks, function(x) sapply (projects_peaks2, function(y) sum(countOverlaps (x, y)) / min (c(length(x), length(y)))))

# prop_sample_df[is.na(prop_sample_df)] = 0
ht = Heatmap (
  # prop_sample_df, 
  peaks_ov_mat,
  col = palette_deviation, 
  row_names_gp= gpar (fontsize=6), 
  column_names_gp= gpar (fontsize=6), 
  column_names_rot = 45)
pdf (file.path ('Plots', 'scatacDatasets_overlap_sample_heatmap.pdf'),width = 6,height=44)
ht
dev.off()

ct_mat = t(scale(t(peaks_ov_mat)))
which.max (ct_mat[, 'C1'])

umap_p12 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
     name = 'Clusters', embedding = "UMAP")
      
pdf (file.path('Plots','clusters_umap.pdf'), 10,10)
umap_p12
dev.off()



# Correlate scatac celltype dev with motif dev ####
if(!exists ('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = as.matrix (assays (mSE)[[1]])
if(!exists ('scaSE')) scaSE = fetch_mat (archp, 'scATAC_datasetsMatrix')
scaMat = as.matrix (assays (scaSE)[[1]])
scaMat = scaMat[!grepl('Tsankov', rownames (scaMat)),]

all (colnames (scaMat) == colnames (mMat))
all (colnames(scaMat) == rownames(archp@cellColData))
scatac_celltype = 'bingren_pan_Astrocyte.2'
scatac_celltype = 'rawlins_fetal_lung_Earlymeso'
scatac_celltype = 'rawlins_fetal_lung_Mid_latemeso'
scatac_tf_cor = lapply (unique (archp$Sample2), function(x) cor (scaMat[scatac_celltype,archp$Sample2 == x], t(mMat)[archp$Sample2 == x,]))






### Compare TF enrichment (HOXB13) between adult and fetal mesothelium ####
projects_peaks_meso = projects_peaks[grep ('meso', names(projects_peaks), ignore.case= T)]
projects_peaks_meso = GRangesList (projects_peaks_meso)
bg_peaks = unlist(projects_peaks_meso)

### TF Enrichment in peaks in large hubs ####
tf_match = getMatches (archp)
colnames (tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)
#regions = range (hubs_obj$hubsCollapsed[c(1:8)])
region_peaks = bg_peakSet[queryHits(findOverlaps(bg_peakSet, projects_peaks_meso[['rawlins_fetal_lung_Mid_latemeso']]))]
#tf_match = tf_match[unique(queryHits (findOverlaps (bg_peakSet, p2gGR)))]
hyper_res =  hyperMotif (
  selected_peaks = region_peaks, 
  motifmatch = tf_match)

#head (hyper_res, 120)
hyper_res['HOXB13',]













