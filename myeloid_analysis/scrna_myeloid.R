conda activate meso_scatac
# use UGER # Add this before running R to be able to run cNMF scripts using UGER 

R
library (Seurat)
library (scran)
library (ggplot2)
library (RColorBrewer)
library (patchwork)
library (clusterProfiler)
library (hdWGCNA)
library (ComplexHeatmap)
library (fgsea)
library (circlize)
library (tidyverse)
library (harmony)
library (ArchR)

set.seed(1234)

# Set project dir
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna'
dir.create (file.path(projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))
scrna_pipeline_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scrna_pipeline'


sample_names = c(
    'P1', # p786
    'P4', # p811
    'P8', # p826
    'P3', # p846
    'P5', #'p848'
    'P10', # p10
    'P11', # p11
    'P12', # p12
    'P13', # p13
    'P14')

sarc_order = c('P1','P13','P3','P12','P5','P11','P4','P8','P14')
#srt = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/srt.rds')
#srt = srt[,srt$celltype_simplified %in% c('cDCs','MonoMac')]

#saveRDS (srt, file.path ('srt.rds'))
srt = readRDS ('srt.rds')
table (srt$celltype2)

cell_subsets_order = c('Mono_CD14','Mono_CD16','TAM_CXCLs','TAM_TREM2','TAM_MARCO','cDCs','TAM_interstitial')

# Initiate pipeline
scrna_pipeline_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scrna_pipeline'
org = 'human'

### Data processing and clustering variables ###
batch = 'sampleID'
variablefeatures = 'scran'
nfeat = 2000 # number of variable genes to consider for dimentionality reduction
sigPCs = 15
ccRegress = T # Regress cell cycle gene expression 
vars_to_regress='nCount_RNA'
metaGroupNames = c('sampleID','celltype2')
res = c(2) # denovo cluster resolutions 

source (file.path(scrna_pipeline_dir,'data_processing.R'))
source (file.path(scrna_pipeline_dir,'harmony_and_clustering.R'))

# De novo marker discovery
enricher_universe = 'all'
logfcThreshold = .25
pvalAdjTrheshold = 0.01
metaGroupName = 'sampleID_harmony_snn_res.2'
top_pathways = 10
top_genes = 5
force = T
source (file.path(scrna_pipeline_dir,'DEG_standard.R'))

srt$celltype_lv2 = srt$celltype2
srt$celltype_lv2[srt$sampleID_harmony_snn_res.2 %in% c(2,8,3,15)] = 'Mono_CD14'
srt$celltype_lv2[srt$sampleID_harmony_snn_res.2 %in% c(7,18)] = 'cCDs'
srt$celltype_lv2[srt$sampleID_harmony_snn_res.2 %in% c(10)] = 'Mono_CD16'
srt$celltype_lv2[srt$sampleID_harmony_snn_res.2 %in% c(0,1,4,16,11,17,13,5,6,12,14)] = 'TAM'
srt$celltype_lv2[srt$sampleID_harmony_snn_res.2 %in% c(9)] = 'Proliferating'
srt$celltype_lv2[srt$sampleID_harmony_snn_res.2 %in% c(19)] = 'Mast'


pdf (file.path ('Plots','celltype_lv2_umap.pdf'))
DimPlot (srt,group.by = 'celltype_lv2', reduction = 'sampleID_harmony_umap')
dev.off()


### Find cNMF shared between myeloid in each sample ####
#### Run cNMF ####
nfeat = 3000
force=F
k_list = c(5:10)
k_selections = c(5:10)
cores= 100

# Run cNMF only in samples matching scATAC-seq samples 
sams = unique (srt$sampleID)
remove_samples = names(table (srt$sampleID)[table(srt$sampleID) < 500])
sams = sams[!sams %in% remove_samples]
srt2 = srt
srt2 = srt2[, srt2$celltype_lv2 != 'Mast'] # remove Mast since not present in scATAC-seq

force = F
for (sam in sams)
	{
	print(sam)	
	srt = srt2[,srt2$sampleID == sam]	
	cnmf_name = paste0('myeloid_',sam)
#	repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
	
	### RUN consensus NMF ####
	source (file.path (scrna_pipeline_dir,'cnmf_prepare_inputs.R'))
	}

#srt = srt2
### Import and format spectra files ####
# sams = unique (archp$Sample)
# sams = sams[sams %in% unique(srt$sampleID)]

k_selection = 10
#sams = sams[sams != 'P10'] # double check why P10 was not included
cnmf_spectra_unique_l = list()
for (sam in sams)
	{
	cnmf_name = paste0('myeloid_',sam)	
	source (file.path (scrna_pipeline_dir,'cnmf_format_spectra_files.R')) 
	cnmf_spectra_unique_l[[sam]] = cnmf_spectra_unique
	}

cnmf_spectra_unique_l = lapply (cnmf_spectra_unique_l, function(x) lapply (x, function(y) head(y, 100)))
saveRDS (cnmf_spectra_unique_l, 'cnmf_myeloid_per_sample.rds')
    
# cnmf_overlap = do.call (cbind, lapply (names(cnmf_spectra_unique_l), function(x)
# 				rowSums (sapply (names(cnmf_spectra_unique_l), function(y)
# 					unlist(lapply (cnmf_spectra_unique_l[[x]], function(z)
# 						sum(unlist(lapply (cnmf_spectra_unique_l[[y]], function(d) sum (z %in% d) > 10)))))))))
# cnmf_overlap = cnmf_overlap - 1

# Select shared programs and plot their gene overlap ####
# overlap_cutoff = 1
# result <- which(cnmf_overlap > overlap_cutoff, arr.ind=TRUE)
# row_names <- rownames(cnmf_overlap)[result[, 1]]
# col_names <- colnames(cnmf_overlap)[result[, 2]]
# # Combine row and column names into a data frame
# output <- data.frame(Row = row_names, Column = col_names)
# print(output)

# ov_mat = do.call (cbind, (lapply (seq(nrow(output)), function(x) cnmf_spectra_unique_l[[output$Column[[x]]]][[output$Row[x]]])))
# colnames (ov_mat) = paste0(output$Column, '_',output$Row)
# ov_mat = as.list(as.data.frame(ov_mat))
force=FALSE
if (!file.exists('shared_cnmf_myeloid.rds') | force)
	{
	cnmf_spectra_unique_l = unlist(cnmf_spectra_unique_l, recursive=F)
	set.seed (123)
	ov_mat = ovmat (cnmf_spectra_unique_l, ov_threshold = 0.2, df=T) 
	km = kmeans (ov_mat, centers=10)
	ho = Heatmap (ov_mat,
	column_split = km$cluster , 
	row_split = km$cluster, 
	row_names_gp = gpar(fontsize = 2),
	column_names_gp = gpar(fontsize = 2),
		col = RColorBrewer::brewer.pal(9,'Blues'))
	pdf (file.path ('Plots','shared_cNMF_overlap4.pdf'),width = 5, height=4)
	print (ho)
	dev.off()
	
	shared_cnmf = split (rownames(ov_mat), km$cluster)
	#remove_metamodule = 1 # remove metamodule include non-overlapping modules
	#shared_cnmf = shared_cnmf[names (shared_cnmf) != as.character(remove_metamodule)]
	overlap_cutoff = 2
	shared_cnmf_genes = lapply (shared_cnmf, function(x) 
		{
		ranked_genes = table (unlist(cnmf_spectra_unique_l[x]))
		ranked_genes = ranked_genes[order(-ranked_genes)]
		names (ranked_genes[ranked_genes > overlap_cutoff])
		})
	shared_cnmf_genes = shared_cnmf_genes[sapply(shared_cnmf_genes, function(x) length(x) > 0)]
	names (shared_cnmf_genes) = paste0('cnmf.',names(shared_cnmf_genes))
	

# Add modules to srt object and plot them on UMAPs ####
srt = ModScoreCor (
    seurat_obj = srt,
    geneset_list = shared_cnmf_genes,
    cor_threshold = NULL,
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'shared_cnmf', outdir = NULL)

reductionName = 'sampleID_harmony_umap'
pdf (file.path ('Plots','shared_cnmf_featplots.pdf'),10,10)
wrap_plots (fp (srt, names (shared_cnmf_genes)))
wrap_plots (fp (srt, c('TREM2','SPP1','C3','VCAN','NFKB1')))
dev.off()


	# Remove cellcycle modules 
	#cellcycle_modules = c('cnmf.2','cnmf.7')
	#shared_cnmf_genes = shared_cnmf_genes[!names (shared_cnmf_genes) %in% cellcycle_modules]
new_cnmf_names = c(
	cnmf.1 = 'TAM_MARCO', 
	cnmf.2 = 'Mono_CD16', 
	cnmf.3 = 'TAM_CXCLs',
	cnmf.4 ='CC',
	cnmf.5 = 'cDCs',
	cnmf.6 = 'Mono_CD14',
	cnmf.7 = 'TAM_TREM2',
	cnmf.8 = 'TAM_unknown',
	cnmf.9 = 'Stress',
	cnmf.10 = 'TAM_interstitial'
	)
names (shared_cnmf_genes) = new_cnmf_names[names (shared_cnmf_genes)]
	saveRDS (shared_cnmf_genes2, paste0 ('shared_cnmf_myeloid.rds'))	
	} else {	
	shared_cnmf_genes = readRDS ('shared_cnmf_myeloid.rds')
	}

shared_cnmf_genes = lapply (shared_cnmf_genes, function(x) head (x,50))
remove_modules = c('Stress','TAM_unknown','CC')
shared_cnmf_genes = shared_cnmf_genes[!names(shared_cnmf_genes) %in% remove_modules]

# Add modules to srt object and plot them on UMAPs ####
srt = ModScoreCor (
    seurat_obj = srt,
    geneset_list = shared_cnmf_genes,
    cor_threshold = NULL,
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'shared_cnmf', outdir = NULL)
srt$celltype_lv3 = srt$shared_cnmf_r_max
#remove_modules = c('CC','Stress','TAM_unknown') # remove monocyres cDC and CC modules. Consider re-inculding CC 
#shared_cnmf_genes = shared_cnmf_genes[!names(shared_cnmf_genes) %in% remove_modules]
# Plot cnmfs 
reductionName = 'sampleID_harmony_umap'

pdf (file.path ('Plots','cnmf_umap.pdf'),10,10)
DimPlot (srt, group.by = 'celltype_lv3', reduction=reductionName)
DimPlot (srt, group.by = 'seurat_clusters', reduction=reductionName, label=T)
dev.off()
srt_cleaned = srt[,!srt$seurat_clusters %in% c('9','19','2','8','14')]

pdf (file.path ('Plots','cnmf_featplots.pdf'),10,10)
DimPlot (srt, group.by = 'sampleID', reduction=reductionName)
wrap_plots (fp (srt, names (shared_cnmf_genes)))
wrap_plots (fp (srt, c('TREM2','SPP1','C3','CSF1R','C1QA','MAF','VCAN','NFKB1','MKI67')))
dev.off()


pdf (file.path ('Plots',''))


#### Run SCENIC ####
force = FALSE
org = 'human'
motif_window = 'tss500bp'#'10kbp'
scenic_name = 'monomac_programs'
genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=5000))
source (file.path(scrna_pipeline_dir, 'SCENIC.R'))

# Run SCENIC plots ####
srt$myeloid = 'myeloid'
motif_window = 'tss500bp'#'10kbp'
genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=5000))
metaGroupNames = c('sampleID','shared_cnmf2_r_max','myeloid')
reductionName = 'sampleID_harmony_umap'
source (file.path(scrna_pipeline_dir, 'SCENIC_plots.R'))

auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

km = readRDS (file.path ('..','scatac_ArchR','TF_activity_modules.rds'))
#genes_highlight = readRDS (file.path ('..','scatac_ArchR','inflammation_TFs.rds'))
regulon_TFs_in_modules = list(
  km1 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 1])],
  km2 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 2])]
  )
auc_mtx = auc_mtx[, colnames(auc_mtx) %in% regulon_TFs_in_modules$km2]

srt$mod_2 = unname(rowMeans (auc_mtx))

# # Try with ridge plots ####
library (ggridges)
library (ggplot2)
library (viridis)
library (tidyr)
#library(hrbrthemes)

# Plot
ccomp = as.data.frame (srt@meta.data)
#ccomp = ccomp[ccomp$cnmf_celltypes %in% c('cDCs'),]
cell_subsets_order = c("TAM_interstitial","TAM_MARCO","cDCs","TAM_TREM2","TAM_CXCLs","Mono_CD16","Mono_CD14")
ccomp$cnmf_celltypes = factor (ccomp$shared_cnmf_r_max, levels = cell_subsets_order)
ccomp$module = srt$mod_2
#ccomp = ccomp[!is.na(ccomp$cnmf_celltypes),]
rp <- ggplot(ccomp, aes(x = module, y = cnmf_celltypes, fill = ..x..)) +
  geom_density_ridges_gradient(
  scale = 3,
  rel_min_height = 0.01,
  linewidth = 0.4,
  color='white',
  alpha = 0.3
) +

  scale_fill_gradientn (colors = palette_expression) +  # Optional: nice color gradient
  theme_ridges() +                      # Optional: clean ridge plot theme
  theme(legend.position = "right")     # Adjust legend position
#   theme_classic() + facet_wrap (~sample, ncol=5)
pdf (file.path ('Plots','cnmf_inflammation_module_ridge_plots.pdf'), width = 5,height=3)
rp
dev.off()


#### Show regulon of TF in inflammation module ####

auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))
auc_mtx = auc_mtx[rownames(auc_mtx) %in% colnames(srt),]
all (rownames (auc_mtx) == colnames (srt))
auc_mtx$celltype = srt$celltype_lv3[match(rownames (auc_mtx), colnames (srt))]
auc_mtx_agg = aggregate (.~ celltype, data = auc_mtx, FUN='mean')
rownames (auc_mtx_agg) = auc_mtx_agg[,1]
auc_mtx_agg = auc_mtx_agg[,-1]

library(tidyverse)

# Assuming your data frame is called df
df <- auc_mtx_agg   # replace with actual name

# Add rownames as a column
df$celltype <- rownames(df)

# Convert to long format
df_long <- df %>%
  pivot_longer(
    cols = -celltype,
    names_to = "gene",
    values_to = "value"
  )

# Plot
pl = ggplot(df_long, aes(x = celltype, y = value, color = gene, group = gene)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "Cell Type",
    y = "Expression",
    color = "Gene"
  )

pdf (file.path ('Plots','regulon_TFs_lineplot.pdf'))
pl
dev.off()



# Check whether regulons recapitulate TF activity from scatac ####
auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

auc_mtx_cor = cor (auc_mtx, method = 'spearman')
set.seed(123)
centers=2
km = kmeans (auc_mtx_cor, centers=centers)

genes_highlight = c(
'IKZF1
HIVEP3
HIVEP1
NFKB2
RELB
RELA
NFKB1
REL
NFE2L2
NFE2
BATF
JDP2
FOSB
SMARCC1
JUND
FOSL1
BACH1
JUNB
FOSL2
FOS
JUN')

genes_highlight2 = c(
'NFE2L2
NFE2
BATF
JDP2
FOSB
SMARCC1
JUND
FOSL1
BACH1
JUNB
FOSL2
FOS
JUN')
# saveRDS (genes_highlight,'inflammation_TFs.rds')

genes_highlight = unlist(strsplit(genes_highlight,'\n'))
genes_highlight2 = unlist(strsplit(genes_highlight2,'\n'))

ha2 = rowAnnotation (foo = anno_mark(at = match(genes_highlight,colnames(auc_mtx_cor)), 
    labels = genes_highlight, labels_gp = gpar(fontsize = 7)))

pdf (file.path ('Plots','TF_modules_heatmap3.pdf'), width = 4,height=3)
cor_mMat_hm = draw (Heatmap (auc_mtx_cor,# row_km=15,
  right_annotation = ha2,
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

pdf (file.path ('Plots','TF_regulon_heatmap.pdf'), width = 4.6, height=3)
cor_mMat_hm
dev.off()


auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))


srt$mod_1 = rowMeans (auc_mtx[,names (km$cluster[km$cluster == 1])])
srt$mod_2 = rowMeans (auc_mtx[,names (km$cluster[km$cluster == 2])])

pdf (file.path ('Plots','gene_regulons_score_fplots.pdf'))
reductionName = 'sampleID_harmony_umap'
fp (srt, 'mod_1')
fp (srt, 'mod_2')
dev.off()


# # Try with ridge plots ####
library (ggridges)
library (ggplot2)
library (viridis)
library (tidyr)
#library(hrbrthemes)

# Plot
ccomp = srt@meta.data[,c('mod_1', 'mod_2','celltype_lv3')]
#ccomp = ccomp[ccomp$celltype_lv2 %in% c('cDCs'),]
#ccomp$celltype_lv2 = factor (ccomp$celltype_lv2, levels = rev(cell_subsets_order))
#median_order = sort (unlist(lapply (split (ccomp$module, ccomp$celltype_lv3), function(x) median(x))))
#ccomp$celltype_lv3 = factor (ccomp$celltype_lv3, levels = names (median_order))
rp <- lapply (c('mod_1', 'mod_2'), function(x) ggplot(ccomp, aes_string(x = x, y = 'celltype_lv3', fill = '..x..')) +
  geom_density_ridges_gradient(
  scale = 3,
  rel_min_height = 0.01,
  linewidth = 0.4,
  color='white',
  alpha = 0.3
) +

  scale_fill_viridis_c(option = "C") +  # Optional: nice color gradient
  theme_ridges() +                      # Optional: clean ridge plot theme
  theme(legend.position = "right"))     # Adjust legend position
#   theme_classic() + facet_wrap (~sample, ncol=5)
pdf (file.path ('Plots','regulon_inflammation_module_ridge_plots.pdf'), width = 9,height=4)
wrap_plots(rp)
dev.off()


# Try using PCA on regulons ####
auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

p <- prcomp(auc_mtx, center = TRUE, scale. = TRUE)
plot_df <- data.frame(p$x, celltype = srt$celltype_lv3[match(rownames (p$x), colnames(srt))])
df_sub <- plot_df[plot_df$celltype %in% 
                    c("Mono_CD14","TAM_interstitial","TAM_TREM2",
                      "TAM_MARCO","TAM_CXCLs"), ]

#-------------------------
# 1. MAIN SCATTER
#-------------------------
p_scatter <- ggplot(df_sub, aes(PC1, PC3, color = celltype)) +
  geom_point(alpha = 0.4, size = 0.1) +
  geom_density_2d(size = 0.4) +
  scale_color_manual(values = palette_myeloid) +
  #scale_x_continuous(limits = x_range, expand = c(0,0)) +
  #scale_y_continuous(limits = y_range, expand = c(0,0)) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )


p_density_x <- ggplot(df_sub, aes(PC1, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) + gtheme_no_rot

p_density_y <- ggplot(df_sub, aes(PC3, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) + gtheme_no_rot

pdf (file.path ('Plots','scrna_momac_PC1_PC3_regulon_pca.pdf'), width = 5,height=5)
p_density_x + plot_spacer() + p_scatter + p_density_y + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4),
    guides = "collect"      # <- collect legends into one
  ) &
  theme(
    legend.position = "bottom",  # or "top"
    legend.box = "horizontal"
  )

dev.off()


### Ridge plot of PC1 ####
df_sub$celltype = factor (df_sub$celltype, levels = rev(cell_subsets_order))
p_ridge <- ggplot(df_sub, aes(x = PC1, y = celltype, fill = celltype)) +
  geom_density_ridges(scale = 3, alpha = 0.7, color = NA) +
  scale_fill_manual (values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y  = element_text(size = 8),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  gtheme_no_rot

pdf(file.path("Plots","momac_PC1_ridge_plots.pdf"),
    width = 5, height = 3)
p_ridge
dev.off()

saveRDS (p, 'PCA_regulons.rds')










# Try using PCA on cells ####
auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

p <- prcomp(t(auc_mtx), center = TRUE, scale. = TRUE)
plot_df <- data.frame(p$x)#, celltype = srt$celltype_lv3[match(rownames (p$x), colnames(srt))])
df_sub <- plot_df

#-------------------------
# 1. MAIN SCATTER
#-------------------------
p_scatter <- ggplot(df_sub, aes(PC1, PC3)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_density_2d(size = 0.4) +
  #scale_color_manual(values = palette_myeloid) +
  #scale_x_continuous(limits = x_range, expand = c(0,0)) +
  #scale_y_continuous(limits = y_range, expand = c(0,0)) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )


p_density_x <- ggplot(df_sub, aes(PC1)) +
  geom_density(alpha = 0.3, color = NA) +
  #scale_fill_manual(values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) + gtheme_no_rot

p_density_y <- ggplot(df_sub, aes(PC3))+
  geom_density(alpha = 0.3, color = NA) +
  #scale_fill_manual(values = palette_myeloid) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) + gtheme_no_rot

pdf (file.path ('Plots','scrna_momac_cells_pca.pdf'), width = 5,height=5)
p_density_x + plot_spacer() + p_scatter + p_density_y + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4),
    guides = "collect"      # <- collect legends into one
  ) &
  theme(
    legend.position = "bottom",  # or "top"
    legend.box = "horizontal"
  )

dev.off()

### Try with UMAPs ####
set.seed(123)  # for reproducibility
library(uwot)
umap_res <- umap(
  auc_mtx,
  n_neighbors = 30,
  min_dist = 0.3,
  metric = "cosine",
  scale = TRUE
)

plot_df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  celltype = srt$celltype_lv3[match(rownames (umap_res), colnames(srt))]
)

df_sub <- plot_df[plot_df$celltype %in%
                    c("Mono_CD14","TAM_interstitial","TAM_TREM2",
                      "TAM_MARCO","TAM_CXCLs"), ]


#-------------------------
# 1. MAIN SCATTER
#-------------------------
p_scatter <- ggplot(df_sub, aes(UMAP1, UMAP2, color = celltype)) +
  geom_point(alpha = 0.4, size = 0.1) +
  geom_density_2d(size = 0.4) +
  scale_color_manual(values = palette_myeloid) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )


#-------------------------
# 2. TOP DENSITY (UMAP1)
#-------------------------
p_density_x <- ggplot(df_sub, aes(UMAP1, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  gtheme_no_rot


#-------------------------
# 3. RIGHT DENSITY (UMAP2)
#-------------------------
p_density_y <- ggplot(df_sub, aes(UMAP2, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  gtheme_no_rot


#-------------------------
# 4. SAVE PDF
#-------------------------
pdf(file.path("Plots","scrna_momac_regulon_umap.pdf"),
    width = 5, height = 5)

p_density_x + plot_spacer() + p_scatter + p_density_y +
  plot_layout(
    ncol = 2, nrow = 2,
    widths  = c(4, 1),
    heights = c(1, 4),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

dev.off()





# Check activity of NFKB1 regulon ####
auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

auc_mtx = aggregate (auc_mtx, by = list(
	celltype = srt$celltype_lv3[match(rownames(auc_mtx), colnames(srt))],
	sample = srt$sampleID[match(rownames(auc_mtx), colnames(srt))]), FUN = mean)
auc_mtx = gather (auc_mtx, TF, score, 3:(ncol(auc_mtx)))

cell_subsets_order = c("TAM_interstitial","TAM_MARCO","TAM_TREM2","TAM_CXCLs","Mono_CD14")
auc_mtx_sub = auc_mtx[auc_mtx$TF == 'NFKB1',]
auc_mtx_sub = auc_mtx_sub[auc_mtx_sub$celltype %in% cell_subsets_order, ]
auc_mtx_sub$celltype = factor (auc_mtx_sub$celltype, levels = rev(cell_subsets_order))

bp = ggplot (auc_mtx_sub, 
	aes(x = celltype, y = score, fill = celltype)) + 
	#vlp +
	geom_boxplot (
    linewidth = .2,
    width=.5,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.9
     ) +
	scale_fill_manual (values = palette_myeloid) + 
	gtheme

pdf (file.path ('Plots','NFKB1_regulon_boxplots_celltypes.pdf'), width = 4, height=4)
bp
dev.off()

# Check activity of NFKB1 expression ####
pb = as.data.frame (log2 (AverageExpression (srt, feature = 'NFKB1', group.by = c('celltype_lv3','sampleID'))[[1]]+1))
pb$TF = 'NFKB1'
pb = gather (pb, celltype, score, 1:((ncol(pb)-1)))
pb$sample = sapply (pb$celltype, function(x) unlist(strsplit(x, '_'))[2])
pb$celltype = sapply (pb$celltype, function(x) unlist(strsplit(x, '_'))[1])
pb$celltype = gsub ('-','_',pb$celltype)
pb = pb[pb$celltype %in% cell_subsets_order,]
pb$celltype = factor (pb$celltype, levels = rev(cell_subsets_order))

bp = ggplot (pb, 
	aes(x = celltype, y = score, fill = celltype)) + 
	#vlp +
	geom_boxplot (
    linewidth = .2,
    width=.5,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.9
     ) +
	scale_fill_manual (values = palette_myeloid) + 
	gtheme

pdf (file.path ('Plots','NFKB1_expression_boxplots_celltypes.pdf'), width = 4, height=4)
bp
dev.off()








### Plot meta cnmf coexpression ####
#mono_mac_cnmfs = c('TREM2','Monocytes','SPP1','IL1B','IFN','IM','C1Q')
### Take also regulon scores 
scenic_name = 'monomac_programs'
auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

library (hdWGCNA)
force = F
# table (srt$sampleID)
# srt$sampleID
#srt_tam = srt[,srt$celltype2 == 'TAMs']
select_cells = c('Mono_CD14','TAM_CXCLs','TAM_TREM2','TAM_MARCO','TAM_interstitial')
if (!file.exists ('metacells.rds') | force)
  {
  srt_sub = srt[,srt$celltype2 %in% c('Mono_CD14','TAMs')]
  srt_sub = srt_sub[,srt_sub$celltype_lv3 %in% select_cells]
  srt_sub = srt_sub[,srt_sub$sampleID %in% c('P1','P10','P11','P12','P13','P14','P4','P5')]

reductionName = 'umap'
reductionSave = 'pca'
reductionGraphKnn = 'RNA_knn'
reductionGraphSnn = 'RNA_snn' 
sigPCs=30
srt_sub  = RunUMAP (object = srt_sub , reduction = reductionSave, dims = 1:sigPCs)
  
pdf (file.path ('Plots','selected_samples_umap.pdf'))
DimPlot (srt_sub, reduction = 'sampleID_harmony_umap', group.by = 'celltype_simplified')
dev.off()

  srt_sub
  srt_sub <- SetupForWGCNA(
    srt_sub,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "metacells" # the name of the hdWGCNA experiment
  )
  # construct metacells  in each group
  srt_sub <- MetacellsByGroups(
    seurat_obj = srt_sub,
    group.by = c("sampleID"), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'sampleID_harmony_umap', # select the dimensionality reduction to perform KNN on
    k = 20, # nearest-neighbors parameter
    max_shared = 12, # maximum number of shared cells between two metacells
    ident.group = 'sampleID' # set the Idents of the metacell seurat object
  )
  
  # normalize metacell expression matrix:
  srt_sub <- NormalizeMetacells (srt_sub)
  metacells = GetMetacellObject (srt_sub)
  saveRDS (metacells, 'metacells.rds')
  } else {
  metacells = readRDS ('metacells.rds') 
  }

auc_mtx_meta = lapply (metacells$cells_merged, function(x) {
  tmp = unlist(strsplit(x, ','))
  colMeans (auc_mtx[tmp,])})
auc_mtx_meta = as.data.frame(do.call (cbind, auc_mtx_meta))


selected_TF = c('FOSL1','FOSL2','NFKB1','NFKB2','REL','TREM2','SPP1','C1QA','C3','MAF','CCL2','EREG','VCAN','IL1B')
metacells_mat = metacells@assays$RNA@layers$data
rownames (metacells_mat) = rownames(metacells)
colnames (metacells_mat) = colnames (metacells)
auc_mtx_meta2 = auc_mtx_meta[rownames(auc_mtx_meta) %in% selected_TF,]
#mMat_meta2 = mMat_meta[selected_TF, ]
metacells_mat = metacells_mat[selected_TF,]
#colnames(metacells_mat) == colnames (mMat_meta2)
hm = Heatmap (t(scale(t(auc_mtx_meta2))), name = 'regulon')
#hm2 = Heatmap (t(scale(t(mMat_meta2))), name = 'motif')
hm3 = Heatmap ((t(scale(t(metacells_mat)))), name = 'expression')

pdf (file.path ('Plots','metacells_regulons_heatmap.pdf'), width=12,height=7)
hm %v% hm3
dev.off()

metacells_mat2 = as.data.frame (metacells_mat)
metacells_mat2$type = 'expression'
metacells_mat2$TF = rownames (metacells_mat2)
auc_mtx_meta2 = as.data.frame (auc_mtx_meta2)
auc_mtx_meta2$type = 'regulon'
auc_mtx_meta2$TF = rownames (auc_mtx_meta2)
metacells_comb = rbind (auc_mtx_meta2, metacells_mat2)
#metacells_comb = gather (metacells_comb, metacell, score, 1:(ncol(metacells_comb)-2))

library(dplyr)
library(tidyr)

gene = 'FOSL2'
tf_long <- metacells_comb %>%
  #rownames_to_column("TF") %>%
  pivot_longer(
    cols = colnames (metacells_mat),
    names_to = "cell",
    values_to = "value"
  )
nfkb1_regulon <- tf_long %>%
  filter(TF == gene, type == "regulon") %>%
  select(cell, nfkb1_regulon = value)

# nfkb1_deviation <- tf_long %>%
#   filter(TF == gene, type == "deviation") %>%
#   select(cell, nfkb1_deviation = value)

trem2_df <- as.data.frame (t(metacells_mat))
trem2_df$cell = rownames(trem2_df)
plot_df = trem2_df
plot_df$sample = sapply (trem2_df$cell, function(x) unlist(strsplit(x,'_'))[1])

# plot_df <- trem2_df %>%
#   left_join(nfkb1_regulon,  by = "cell") %>%
#   left_join(nfkb1_deviation, by = "cell")

library(dplyr)

df_scaled <- plot_df %>%
  group_by(sample) %>%
  mutate(
    across(
      where(is.numeric),
      ~ as.numeric(scale(.))
    )
  ) %>%
  ungroup()

grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "IL1B"
)
rp <- ggplot(
  df_scaled,
  aes(
    x = FOSL2,
    y = NFKB1
  )
) +
  geom_point(aes(color = IL1B), alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  #) +
  #labs(
    #x = paste0(gene," deviation score"),
    #y = paste0(gene," regulon"),
    #title = paste0(gene," deviation vs regulon on IL1B (Mono)")
  ) + facet_wrap (~sample) + gtheme

pdf (file.path ('Plots','regulon_deviation_scatter.pdf'), width=7, height=5.5)
rp
dev.off()





grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "IL1B"
)
rp <- ggplot(
  plot_df,
  aes_string(
    x = 'nfkb1_deviation',
    y = gene
  )
) +
  geom_point(aes(color = IL1B), alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste0(gene," deviation score"),
    y = paste0(gene," expression"),
    title = paste0(gene," deviation vs regulon on IL1B (Mono)")
  ) + gtheme

grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "MAF"
)
dp <- ggplot(
  plot_df,
  aes_string(
    x = 'nfkb1_deviation',
    y = gene
  )
) +
  geom_point(aes(color = MAF), alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste0(gene," deviation score"),
    y = paste0(gene," expression"),
    title = paste0(gene," deviation vs regulon on MAF (IM)")
  ) + gtheme



grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "TREM2"
)

p_main <- ggplot(
  plot_df,
  aes_string(
    x = 'nfkb1_deviation',
    y = gene,
    color = 'TREM2'
  )
) +
  geom_point(alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste0(gene," deviation score"),
    y = paste0(gene," expression")#,
    #title = "NFKB1 deviation vs regulon colored by TREM2"
  ) + gtheme_no_rot

pdf (file.path ('Plots',paste0(gene,'_expression_deviation_scatter.pdf')), width=7, height=3.5)
wrap_plots (rp, dp, p_main)
dev.off()


