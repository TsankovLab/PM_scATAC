conda activate meso_scatac
use UGER # Add this before running R to be able to run cNMF scripts using UGER 

R
library (Seurat)
library (scran)
library (ggplot2)
library (RColorBrewer)
library(patchwork)
library (clusterProfiler)
library (hdWGCNA)
library(ComplexHeatmap)
library (fgsea)
library (circlize)
library (tidyverse)

set.seed(1234)

# Set project dir
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/Endothelial/scrna'
dir.create (file.path(projdir,'Plots'), recursive =T)
setwd (projdir)
source ('../../git_repo/utils/useful_functions.R')
source ('../../git_repo/utils/palettes.R')
source ('../../git_repo/utils/ggplot_aestetics.R')


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


# Load RNA
srt = readRDS (file.path('..','..','stroma','scrna','srt.rds'))
srt = srt[, srt$celltype == 'Endothelial']


	
	# Run Harmony
	
	srt = NormalizeData (srt)
	sce = SingleCellExperiment (list(counts=srt@assays$RNA@counts, logcounts = srt@assays$RNA@data),
	rowData=rownames(srt)) 
	sce = modelGeneVar(sce)
	# remove batchy genes
	batchy_genes = c('RPL','RPS','MT-')
	sce = sce[!apply(sapply(batchy_genes, function(x) grepl (x, rownames(sce))),1,any),]
	nfeat = 5000
	vf = getTopHVGs (sce, n=nfeat)
	VariableFeatures (srt) = vf
	srt <- ScaleData (srt)
	srt <- RunPCA (srt)
	batch = 'sampleID'
	reductionSave = 'harmony'
	srt = srt %>% 
	RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
	RunUMAP (reduction = reductionSave, dims = 1:15)

	#ref <- RunUMAP (ref, dims = 1:15)
	srt = FindNeighbors (object = srt, reduction = 'harmony', dims = 1:15, k.param = 30,
                              verbose = TRUE, force.recalc = T)
	srt = FindClusters (srt, resolution = .8, verbose = T, n.start = 100)

	pdf ('Plots/samples_umap.pdf', width=12)
	FeaturePlot (srt, feature = c('VWF','PECAM1','PLVAP','NOTCH4','COL4A1'), reduction = 'umap')
	DimPlot (srt, group.by = 'seurat_clusters', reduction = 'umap', label=T)
	DimPlot (srt, group.by = 'sampleID', reduction = 'umap', label=F)
	DimPlot (srt, group.by = 'celltype2', reduction = 'umap')
	dev.off()
	
	srt$celltype2='Vein'
	srt$celltype2[srt$seurat_clusters %in% c(2,3,7)] = 'COL4A1'
	srt$celltype2[srt$seurat_clusters %in% c(6,8)] = 'Artery'

	pdf ('Plots/samples_umap.pdf', width=12)
	FeaturePlot (srt, feature = c('VWF','PECAM1','PLVAP','NOTCH4','COL4A1'), reduction = 'umap')
	DimPlot (srt, group.by = 'seurat_clusters', reduction = 'umap', label=T)
	DimPlot (srt, group.by = 'sampleID', reduction = 'umap', label=F)
	DimPlot (srt, group.by = 'celltype2', reduction = 'umap')
	dev.off()

saveRDS (srt, 'srt.rds')	
