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

archp = loadArchRProject (file.path('..','scatac_ArchR'))

# Load scRNA ####
if (!file.exists ('srt.rds'))
	{
	#srt = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/reproduction2/scRNA/GSE190597_srt_tumor.rds')
	srt = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/srt.rds')
	
	# Import P14
	srt_p14 = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/MPM_naive_p14_analysis/_cellranger_raw_Filter_400_1000_25/no_harmony/srt.rds')
	srt_p14$celltype_simplified = srt_p14$celltype
	srt = merge (srt, srt_p14)
	srt$sampleID[srt$sampleID == 'MPM_naive_p14'] = 'P14'
	srt = srt[, srt$celltype_simplified %in% c('MonoMac','cDCs')]
	srt = srt[, srt$sampleID %in% sample_names]
	
	srt = NormalizeData (srt)
	sce = SingleCellExperiment (list(counts=srt@assays$RNA@counts, logcounts = srt@assays$RNA@data),
	rowData=rownames(srt)) 
	sce = modelGeneVar(sce)
	# remove batchy genes
	batchy_genes = c('RPL','RPS','MT-')
	sce = sce[!apply(sapply(batchy_genes, function(x) grepl (x, rownames(sce))),1,any),]
	nfeat = 5000
	vf = getTopHVGs (sce, n=nfeat)
	
	# remove ambient RNA genes
	vf = vf[!grepl ('SFTP',vf)]
	vf = vf[!grepl ('HB',vf)]
	VariableFeatures (srt) = vf
	srt <- ScaleData (srt)
	srt <- RunPCA (srt)
	
	batch = 'sampleID'
	reductionSave = 'harmony'
	srt = srt %>% 
	RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
	RunUMAP (reduction = reductionSave, dims = 1:15)

	srt = FindNeighbors (object = srt, reduction = reductionSave, dims = 1:15, k.param = 30,
                              verbose = TRUE)
	srt = FindClusters (srt, resolution = 1, verbose = T, n.start = 100)

	reductionName = 'umap'
	fps = fp (srt, gene = c('CD3D','KLRC1','FGFBP2','CD8A','CD4','FOXP3','CXCL13','CTLA4','PDCD1','CXCR5','MKI67'))
	dp = wrap_plots (DimPlot (srt, group.by = 'sampleID',reduction='umap'), DimPlot (srt, group.by = 'seurat_clusters', label=T,reduction='umap'),
		DimPlot (srt, group.by = 'celltype',reduction='umap'))
	pdf (file.path('Plots','markers_celltypes_umap.pdf'), width=12)
	wrap_plots (fps)
	dev.off()
	pdf (file.path('Plots','celltypes_umap.pdf'),width=15)
	wrap_plots (dp)
	dev.off()
	

metaGroupName = 'celltype2'
if (!file.exists ('DEG.rds'))
	{
	DefaultAssay(srt) = 'RNA'
	Idents(srt) = metaGroupName
	logfcThreshold=.5
	pvalAdjTrheshold = .05
	degClusters = FindAllMarkers (srt, max.cells.per.ident = 1000, min.pct = .1, logfc.threshold = logfcThreshold, verbose = T)
	saveRDS (degClusters, 'DEG.rds')
	write.csv (degClusters[degClusters$p_val_adj < pvalAdjTrheshold,], 'DEG.csv')
	} else {
	degClusters = readRDS ('DEG.rds')		
	}

message ('Generate Heatmap of top genes') # does not produce image when there are too many genes to plot
degClusters = degClusters[degClusters$p_val_adj < pvalAdjTrheshold, ]

top_genes = 5
	top_deg = degClusters[degClusters$avg_log2FC > 0,]
	top_deg = degClusters %>% group_by (cluster) %>% top_n (top_genes, avg_log2FC)
	#top_deg = top_deg %>% group_by (cluster) %>% top_n (top_genes, -p_val_adj)
	
# if (!is.na(is.integer(as.integer(srt@meta.data[,metaGroupName][1]))))
#  	{
#  	srt@meta.data[,metaGroupName] = factor (srt@meta.data[,metaGroupName], levels =unique(srt@meta.data[,metaGroupName])[order(as.numeric(as.character(unique(srt@meta.data[,metaGroupName]))))])
#  	} else {
# 	srt@meta.data[,metaGroupName] = factor (srt@meta.data[,metaGroupName], levels = unique (srt@meta.data[,metaGroupName]))
#  	}


#top_deg = top_deg[order (-top_deg$avg_log2FC),]
srt2 = ScaleData (srt, features = unique(top_deg$gene))

#srt2@meta.data[,metaGroupName] = factor(srt2@meta.data[,metaGroupName], levels = unique (srt2@meta.data[,metaGroupName]))
#srt2 = srt[, !is.na(srt@meta.data[,metaGroupName])]
#if (is.na(as.integer(levels(srt@meta.data[,metaGroupName])))) Idents(srt2) = factor (srt2@meta.data[,metaGroupName], levels =unique(srt2@meta.data[,metaGroupName])[order(as.numeric(as.character(unique(srt2@meta.data[,metaGroupName]))))])
#if (!is.na(is.integer(as.integer(srt@meta.data[,metaGroupName][1])))) Idents(srt2) = factor (srt2@meta.data[,metaGroupName], levels =unique(srt2@meta.data[,metaGroupName])[order(as.numeric(as.character(unique(srt2@meta.data[,metaGroupName]))))])
heat_p = DoHeatmap (srt2, features = unique(top_deg$gene))
png (file.path('Plots',paste0('top_',top_genes,'_heatmap.png')), height=2500, width=5300, res=400)
print (heat_p)
dev.off()
rm (srt2)

#### Generate Feature plots and dotplots of top genes ####
message ('Generate Feature plots of top genes')
top_deg = degClusters %>% arrange(-avg_log2FC) %>% group_by (cluster)  %>% slice(1:top_genes)
#top_deg = top_deg[order (-top_deg$avg_log2FC),]

feat_p = lapply (seq_along(top_deg$gene), function(x) FeaturePlot (srt, features = top_deg$gene[x], keep.scale = 'all', order=F,
  combine = FALSE, pt.size = .01, reduction = reductionName))
feat_p = unlist(feat_p, recursive=F)
for(i in 1:length(feat_p)) {
    feat_p[[i]] = feat_p[[i]] + theme_void() + NoLegend() + NoAxes() + ggtitle (paste0('CL',top_deg$cluster[i],':',top_deg$gene[i])) +
  scale_colour_gradientn (colours = viridis::turbo(100))
    }
  
png (file.path('Plots',paste0('top_',top_genes,'_fplots.png')), ((length(feat_p)*50) + 1500),((length(feat_p)*50) + 1500), res=300)
print (wrap_plots (feat_p))
dev.off()

dp = DotPlot(object = srt, features = rev(unique(top_deg$gene)), scale = T, group.by = metaGroupName) +
  	theme(axis.text.x = element_text(angle = 45, hjust=1), panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_line(colour = "gainsboro")) + 
    scale_color_gradientn(colours = rev(brewer.pal(11,"Spectral"))) +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)

png (file.path('Plots',paste0('top_',top_genes,'_dotplot.png')), ((length(unique(top_deg$gene))*90) + 500), ((length(unique (srt@meta.data[,metaGroupName]))*102) + 250), res=300)
dp
dev.off()

# Check ambient RNA ####
genes = c('HBB','SFTPC')
vp = VlnPlot (srt, features = genes, group.by ='sampleID')
pdf (file.path('Plots',paste0('check_ambientRNA_vlnplot.pdf')))
vp
dev.off()


saveRDS (srt, 'srt.rds')

} else {
srt = readRDS ('srt.rds')	
}


### Find cNMF shared between myeloid in each sample ####
#### Run cNMF ####
nfeat = 5000
force=F
k_list = c(5:10)
k_selections = c(5:10)
cores= 100

# Run cNMF only in samples matching scATAC-seq samples 
sams = unique (archp$Sample)
sams = sams[sams %in% unique(srt$sampleID)]
srt2 = srt[, srt$sampleID %in% sams]
#srt2 = srt2[, srt2$sampleID %in% names (table (srt2$sampleID)[table (srt2$sampleID) > 100])]
#srt2 = srt
sams = unique (srt2$sampleID)
for (sam in sams)
	{
	print(sam)	
	srt = srt2[,srt2$sampleID == sam]	
	cnmf_name = paste0('myeloid_',sam)
	cnmf_out = paste0('cNMF/cNMF_',cnmf_name,'_',paste0(k_list[1],'_',k_list[length(k_list)]),'_vf',nfeat)
	dir.create (file.path(cnmf_out,'Plots'), recursive=T)
	repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
	
	### RUN consensus NMF ####
	source (file.path ('..','..','git_repo','utils','cnmf_prepare_inputs.R'))
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
	source (file.path ('..','..','git_repo','utils','cnmf_format_spectra_files.R')) 
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
	ho
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
	# Remove cellcycle modules 
	#cellcycle_modules = c('cnmf.2','cnmf.7')
	#shared_cnmf_genes = shared_cnmf_genes[!names (shared_cnmf_genes) %in% cellcycle_modules]
new_cnmf_names = c(
	cnmf.2 = 'CC_G2M', 
	cnmf.1 = 'SPP1', 
	cnmf.3 = 'LILRA',
	cnmf.4 ='IM',
	cnmf.5 = 'TREM2',
	cnmf.6 = 'Mono',
	cnmf.7 = 'IFN_CXCLs',
	cnmf.8 = 'cDCs2',
	cnmf.9 = 'CC_S',
	cnmf.10 = 'cDCs'
	)
names (shared_cnmf_genes) = new_cnmf_names[names (shared_cnmf_genes)]
	saveRDS (shared_cnmf_genes, paste0 ('shared_cnmf_myeloid.rds'))	
	} else {	
	shared_cnmf_genes = readRDS ('shared_cnmf_myeloid.rds')
	}

shared_cnmf_genes = lapply (shared_cnmf_genes, function(x) head (x,50))

# Add modules to srt object and plot them on UMAPs ####
srt = ModScoreCor (
    seurat_obj = srt,
    geneset_list = shared_cnmf_genes,
    cor_threshold = NULL,
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'shared_cnmf', outdir = NULL)

remove_modules = c('cDCs2','CC_S','CC_G2M','LILRA') # remove monocyres cDC and CC modules. Consider re-inculding CC 
shared_cnmf_genes = shared_cnmf_genes[!names(shared_cnmf_genes) %in% remove_modules]
# Plot cnmfs 
reductionName = 'umap'
pdf (file.path ('Plots','cnmf_featplots.pdf'),10,10)
wrap_plots (fp (srt, names (shared_cnmf_genes)))
wrap_plots (fp (srt, c('TREM2','SPP1','C3','VCAN','NFKB1')))
dev.off()

### Annotate based on cnmfs ####
set.seed (123)
cnmf_scrna = srt@meta.data[, names(shared_cnmf_genes)]

km_cnmf = kmeans (scale(cnmf_scrna), centers=6) # double scale modules and cluster using k-means
#ha = HeatmapAnnotation (sample = srt$sampleID[, col=list(sample = palette_sample[sams]))
hm = Heatmap (t(scale(cnmf_scrna)), 
  col = palette_genescore_fun(cnmf_scrna), 
#  top_annotation = ha,
  clustering_distance_columns = 'pearson',
  clustering_distance_rows = 'pearson',
  show_column_dend = F,
  column_split = km_cnmf$cluster,
  #column_km=3,
  row_names_gp = gpar (fontsize = 8),
  column_names_gp = gpar (fontsize = 0),
  border=T)

pdf (file.path ('Plots','cnmf_scrna_scaled_heatmap.pdf'), height=1.5)
hm
dev.off()

### Annotate based on cnmfs ####
shared_cnmf_genes2 = shared_cnmf_genes

srt = ModScoreCor (
    seurat_obj = srt, 
    geneset_list = shared_cnmf_genes2, 
    cor_threshold = NULL, 
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'shared_cnmf2', outdir = NULL)

pdf (file.path ('Plots','cnmf_named_umap.pdf'))
DimPlot (srt, group.by = 'shared_cnmf2_r_max', reduction='umap')
dev.off()

sv()

### Plot meta cnmf coexpression ####
#mono_mac_cnmfs = c('TREM2','Monocytes','SPP1','IL1B','IFN','IM','C1Q')
ccomp_df = srt@meta.data[,names (shared_cnmf_genes)]
ccomp_df_cor = cor (ccomp_df, method = 'spearman')

pdf (file.path ('Plots','scRNA_meta_cnmf_coexpression_heatmap2.pdf'),4,width=5)
Heatmap (ccomp_df_cor, col = palette_expression_cor_fun (ccomp_df_cor), border=T)
dev.off()


library (hdWGCNA)
force = F
# table (srt$sampleID)
# srt$sampleID
#srt_tam = srt[,srt$celltype2 == 'TAMs']
if (!file.exists ('metacells.rds') | force)
	{
	srt <- SetupForWGCNA(
	  srt,
	  gene_select = "fraction", # the gene selection approach
	  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
	  wgcna_name = "metacells" # the name of the hdWGCNA experiment
	)
	# construct metacells  in each group
	srt <- MetacellsByGroups(
	  seurat_obj = srt,
	  group.by = c("sampleID"), # specify the columns in seurat_obj@meta.data to group by
	  reduction = 'umap', # select the dimensionality reduction to perform KNN on
	  k = 100, # nearest-neighbors parameter
	  max_shared = 50, # maximum number of shared cells between two metacells
	  ident.group = 'sampleID' # set the Idents of the metacell seurat object
	)
	
	# normalize metacell expression matrix:
	srt <- NormalizeMetacells (srt)
	metacells = GetMetacellObject (srt)
	saveRDS (metacells, 'metacells.rds')
	} else {
	metacells = readRDS ('metacells.rds')	
	}

# install.packages
# Correlate genes with cNMFs on metacells ####
vf = VariableFeatures (FindVariableFeatures (metacells, nfeat=20000))
metacells_assay = metacells@assays$RNA@layers$data
rownames (metacells_assay) = rownames (srt)
metacells_assay = metacells_assay[unique(c(vf)),]

metacells = ModScoreCor (
        seurat_obj = metacells, 
        geneset_list = shared_cnmf_genes, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'cNMF_', outdir = paste0(projdir,'Plots/'))

ccomp_df = metacells@meta.data[,c(names(shared_cnmf_genes))]

tc_cor = lapply (unique(metacells$sampleID), function(x)
	{
	metacells_assay_sample = t(metacells_assay[vf,metacells$sampleID == x])
	ccomp_df_sample = ccomp_df[metacells$sampleID == x,]
	res = cor (metacells_assay_sample, ccomp_df_sample, method = 'spearman')
	#rownames (res) = colnames (ccomp_df_sample)
	res
	})	
names (tc_cor) = unique(metacells$sampleID)
trem2_mod = tc_cor[[1]][,6]
trem2_mod = trem2_mod[order(-trem2_mod)]
genes = c('TREM2','SPP1','BACH2',
	'JUN','FOS','C1QA',
	'C1QB','STAT2','C3',
	'PRDM1','IRF1','CEBPZ','NFYC','SP2','IRF3',
	'PBX3','SP2','SMARCC1','JUNB','JDP2','NFEL2L2')
trem2_mod[genes]
#lapply (tc_cor, function(x) {x = x['cNMF19',]; head(x[order(-x)],10)})

# Check expression of ELK1 and other TFs from chromvar / chrombpnet predictions  ####
TF = 'ELK3'
#srt = srt[,srt$celltype2 == 'TAMs']

pdf (file.path ('Plots',paste0(TF,'expression_dotplot.pdf')))
DotPlot (srt, features = TF, group.by = 'celltype2')
dev.off()

module = 'cnmf.4'
scrna_cnmf = srt@meta.data[, names(shared_cnmf_genes)]
scrna_cnmf_sample = lapply (unique (srt$sampleID), function(x) scrna_cnmf[srt$sampleID == x,module,drop=F])
scrna_cnmf_sample = lapply (scrna_cnmf_sample, function(x) head(rownames(x)[order(-x[,1])],20))
srt$SPP1_high = ifelse (colnames(srt) %in% unlist (scrna_cnmf_sample), 'SPP1_high','SPP1_low')

TF = c('ELK1','ELK3','KLF12','NFKB1','CTCF','CEBPA','NFYB','IRF1','NFIC','IRF8','ZNF76','ZEB1','MAF','SNAI1','ZN143')
# pdf (file.path ('Plots',paste0(TF,'expression_dotplot.pdf')),width=12, height=8)
# DotPlot (srt_tam, features = TF, group.by = 'SPP1_high')
# VlnPlot (srt_tam, features = TF, group.by = 'SPP1_high',pt.size=0)
# dev.off()	

mod_sample_df = data.frame (
	module = srt$SPP1_high, 
	score = srt@meta.data[,module],
	sample = srt$sampleID)
mod_sample_df = cbind(mod_sample_df)
mod_sample_df = aggregate (
	t(srt@assays$RNA@data[rownames(srt) %in% TF,]),
	by = list(sampleID= mod_sample_df$sample,
		module = mod_sample_df$module), mean)

mod_sample_df = gather (mod_sample_df, TF, expression, 3:ncol(mod_sample_df))
bp = ggplot (mod_sample_df, aes (x = TF, y = expression, fill = module)) + 
geom_boxplot() + gtheme

pdf (file.path ('Plots','module_sample_boxplot.pdf'))
bp
dev.off()

active_TFs = read.csv ('../scatac_ArchR/active_TFs.csv')

reductionName = 'umap'
fps = fp (srt, gene = c('TREM2','JUN','C1QA','C3','MAF','NFATC2','VCAN','APOE','SPP1','A2M'))
fps = fp (srt, gene = active_TFs[[2]])
pdf (file.path('Plots','markers_celltypes_umap.pdf'), width=52, height=50)
wrap_plots (DimPlot (srt, group.by = 'sampleID', reduction=reductionName),DimPlot (srt, group.by = 'celltype2', reduction=reductionName))
wrap_plots (fps)
dev.off()

pdf (file.path('Plots','modules_umap.pdf'), width=10, height=10)
wrap_plots (fp (srt, gene = names(shared_cnmf_genes)))
dev.off()


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
reductionName = 'umap'
source (file.path(scrna_pipeline_dir, 'SCENIC_plots.R'))

auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

km = readRDS (file.path ('..','scatac_ArchR','TF_activity_modules.rds'))
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
ccomp$cnmf_celltypes = factor (ccomp$shared_cnmf_r_max, levels = rev(c('Mono','TREM2','SPP1','IFN_CXCLs','cDCs','IM')))
ccomp$module = srt$mod_2
ccomp = ccomp[!is.na(ccomp$cnmf_celltypes),]
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





