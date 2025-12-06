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
pdf (file.path ('Plots','cnmf_featplots.pdf'),10,10)
wrap_plots (fp (srt, names (shared_cnmf_genes)))
wrap_plots (fp (srt, c('TREM2','SPP1','C3','VCAN','NFKB1')))
dev.off()

# ### Annotate based on cnmfs ####
# set.seed (123)
# cnmf_scrna = srt@meta.data[, names(shared_cnmf_genes)]

# km_cnmf = kmeans (scale(cnmf_scrna), centers=6) # double scale modules and cluster using k-means
# #ha = HeatmapAnnotation (sample = srt$sampleID[, col=list(sample = palette_sample[sams]))
# hm = Heatmap (t(scale(cnmf_scrna)), 
#   col = palette_genescore_fun(cnmf_scrna), 
# #  top_annotation = ha,
#   clustering_distance_columns = 'pearson',
#   clustering_distance_rows = 'pearson',
#   show_column_dend = F,
#   column_split = km_cnmf$cluster,
#   #column_km=3,
#   row_names_gp = gpar (fontsize = 8),
#   column_names_gp = gpar (fontsize = 0),
#   border=T)

# pdf (file.path ('Plots','cnmf_scrna_scaled_heatmap.pdf'), height=1.5)
# hm
# dev.off()

# ### Annotate based on cnmfs ####
# shared_cnmf_genes2 = shared_cnmf_genes

# srt = ModScoreCor (
#     seurat_obj = srt, 
#     geneset_list = shared_cnmf_genes2, 
#     cor_threshold = NULL, 
#     pos_threshold = NULL, # threshold for fetal_pval2
#     listName = 'shared_cnmf2', outdir = NULL)

# pdf (file.path ('Plots','cnmf_named_umap.pdf'))
# DimPlot (srt, group.by = 'shared_cnmf2_r_max', reduction='umap')
# dev.off()

# sv()

# ### Plot meta cnmf coexpression ####
# #mono_mac_cnmfs = c('TREM2','Monocytes','SPP1','IL1B','IFN','IM','C1Q')
# ccomp_df = srt@meta.data[,names (shared_cnmf_genes)]
# ccomp_df_cor = cor (ccomp_df, method = 'spearman')

# pdf (file.path ('Plots','scRNA_meta_cnmf_coexpression_heatmap2.pdf'),4,width=5)
# Heatmap (ccomp_df_cor, col = palette_expression_cor_fun (ccomp_df_cor), border=T)
# dev.off()


# library (hdWGCNA)
# force = F
# # table (srt$sampleID)
# # srt$sampleID
# #srt_tam = srt[,srt$celltype2 == 'TAMs']
# if (!file.exists ('metacells.rds') | force)
# 	{
# 	srt <- SetupForWGCNA(
# 	  srt,
# 	  gene_select = "fraction", # the gene selection approach
# 	  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
# 	  wgcna_name = "metacells" # the name of the hdWGCNA experiment
# 	)
# 	# construct metacells  in each group
# 	srt <- MetacellsByGroups(
# 	  seurat_obj = srt,
# 	  group.by = c("sampleID"), # specify the columns in seurat_obj@meta.data to group by
# 	  reduction = 'umap', # select the dimensionality reduction to perform KNN on
# 	  k = 100, # nearest-neighbors parameter
# 	  max_shared = 50, # maximum number of shared cells between two metacells
# 	  ident.group = 'sampleID' # set the Idents of the metacell seurat object
# 	)
	
# 	# normalize metacell expression matrix:
# 	srt <- NormalizeMetacells (srt)
# 	metacells = GetMetacellObject (srt)
# 	saveRDS (metacells, 'metacells.rds')
# 	} else {
# 	metacells = readRDS ('metacells.rds')	
# 	}

# # install.packages
# # Correlate genes with cNMFs on metacells ####
# vf = VariableFeatures (FindVariableFeatures (metacells, nfeat=20000))
# metacells_assay = metacells@assays$RNA@layers$data
# rownames (metacells_assay) = rownames (srt)
# metacells_assay = metacells_assay[unique(c(vf)),]

# metacells = ModScoreCor (
#         seurat_obj = metacells, 
#         geneset_list = shared_cnmf_genes, 
#         cor_threshold = NULL, 
#         pos_threshold = NULL, # threshold for fetal_pval2
#         listName = 'cNMF_', outdir = paste0(projdir,'Plots/'))

# ccomp_df = metacells@meta.data[,c(names(shared_cnmf_genes))]

# tc_cor = lapply (unique(metacells$sampleID), function(x)
# 	{
# 	metacells_assay_sample = t(metacells_assay[vf,metacells$sampleID == x])
# 	ccomp_df_sample = ccomp_df[metacells$sampleID == x,]
# 	res = cor (metacells_assay_sample, ccomp_df_sample, method = 'spearman')
# 	#rownames (res) = colnames (ccomp_df_sample)
# 	res
# 	})	
# names (tc_cor) = unique(metacells$sampleID)
# trem2_mod = tc_cor[[1]][,6]
# trem2_mod = trem2_mod[order(-trem2_mod)]
# genes = c('TREM2','SPP1','BACH2',
# 	'JUN','FOS','C1QA',
# 	'C1QB','STAT2','C3',
# 	'PRDM1','IRF1','CEBPZ','NFYC','SP2','IRF3',
# 	'PBX3','SP2','SMARCC1','JUNB','JDP2','NFEL2L2')
# trem2_mod[genes]
# #lapply (tc_cor, function(x) {x = x['cNMF19',]; head(x[order(-x)],10)})

# # Check expression of ELK1 and other TFs from chromvar / chrombpnet predictions  ####
# TF = 'ELK3'
# #srt = srt[,srt$celltype2 == 'TAMs']

# pdf (file.path ('Plots',paste0(TF,'expression_dotplot.pdf')))
# DotPlot (srt, features = TF, group.by = 'celltype2')
# dev.off()

# module = 'cnmf.4'
# scrna_cnmf = srt@meta.data[, names(shared_cnmf_genes)]
# scrna_cnmf_sample = lapply (unique (srt$sampleID), function(x) scrna_cnmf[srt$sampleID == x,module,drop=F])
# scrna_cnmf_sample = lapply (scrna_cnmf_sample, function(x) head(rownames(x)[order(-x[,1])],20))
# srt$SPP1_high = ifelse (colnames(srt) %in% unlist (scrna_cnmf_sample), 'SPP1_high','SPP1_low')

# TF = c('ELK1','ELK3','KLF12','NFKB1','CTCF','CEBPA','NFYB','IRF1','NFIC','IRF8','ZNF76','ZEB1','MAF','SNAI1','ZN143')
# # pdf (file.path ('Plots',paste0(TF,'expression_dotplot.pdf')),width=12, height=8)
# # DotPlot (srt_tam, features = TF, group.by = 'SPP1_high')
# # VlnPlot (srt_tam, features = TF, group.by = 'SPP1_high',pt.size=0)
# # dev.off()	

# mod_sample_df = data.frame (
# 	module = srt$SPP1_high, 
# 	score = srt@meta.data[,module],
# 	sample = srt$sampleID)
# mod_sample_df = cbind(mod_sample_df)
# mod_sample_df = aggregate (
# 	t(srt@assays$RNA@data[rownames(srt) %in% TF,]),
# 	by = list(sampleID= mod_sample_df$sample,
# 		module = mod_sample_df$module), mean)

# mod_sample_df = gather (mod_sample_df, TF, expression, 3:ncol(mod_sample_df))
# bp = ggplot (mod_sample_df, aes (x = TF, y = expression, fill = module)) + 
# geom_boxplot() + gtheme

# pdf (file.path ('Plots','module_sample_boxplot.pdf'))
# bp
# dev.off()

# active_TFs = read.csv ('../scatac_ArchR/active_TFs.csv')

# reductionName = 'umap'
# fps = fp (srt, gene = c('TREM2','JUN','C1QA','C3','MAF','NFATC2','VCAN','APOE','SPP1','A2M'))
# fps = fp (srt, gene = active_TFs[[2]])
# pdf (file.path('Plots','markers_celltypes_umap.pdf'), width=52, height=50)
# wrap_plots (DimPlot (srt, group.by = 'sampleID', reduction=reductionName),DimPlot (srt, group.by = 'celltype2', reduction=reductionName))
# wrap_plots (fps)
# dev.off()

# pdf (file.path('Plots','modules_umap.pdf'), width=10, height=10)
# wrap_plots (fp (srt, gene = names(shared_cnmf_genes)))
# dev.off()


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
p_scatter <- ggplot(df_sub, aes(PC1, PC2, color = celltype)) +
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

p_density_y <- ggplot(df_sub, aes(PC2, fill = celltype)) +
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

pdf (file.path ('Plots','scrna_momac_regulon_pca.pdf'), width = 5,height=5)
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
