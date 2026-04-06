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
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/stroma/scrna'
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

# Load scRNA ####
if (!file.exists ('srt.rds'))
	{
	#srt = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/reproduction2/scRNA/GSE190597_srt_tumor.rds')
	srt = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/srt.rds')
	# Import P14
	srt_p14 = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_p14_analysis/_cellranger_raw_Filter_400_1000_25/no_harmony/srt.rds')
	srt_p14$celltype_simplified = srt_p14$celltype
	srt = merge (srt, srt_p14)
	srt$sampleID[srt$sampleID == 'MPM_naive_p14'] = 'P14'
	srt = srt[, srt$celltype_simplified %in% c('Endothelial','Fibroblasts','LEC','Mesothelium','SmoothMuscle')]
	srt = srt[, srt$sampleID %in% sample_names]
	
	# Import normal mesothelium RNA
	projdir_ref = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scRNA_meso/'
	ref = readRDS (paste0(projdir_ref, 'ref.rds'))	
	#subsample_ct = 1000
	#ref_bc = lapply (unique(ref$celltype), function(x) sample(colnames(ref)[ref$celltype == x], subsample_ct, replace=T))
	#names (ref_bc) = unique(ref$celltype)
	#cts = c('Mesothelium','Fibroblast')
	#ref = ref[,colnames (ref) %in% unlist(ref_bc[cts])]
	
	ref = NormalizeData (ref)
	sce = SingleCellExperiment (list(counts=ref@assays$RNA@counts, logcounts = ref@assays$RNA@data),
	rowData=rownames(ref)) 
	sce = modelGeneVar(sce)
	# remove batchy genes
	batchy_genes = c('RPL','RPS','MT-')
	sce = sce[!apply(sapply(batchy_genes, function(x) grepl (x, rownames(sce))),1,any),]
	nfeat = 2000
	vf = getTopHVGs (sce, n=nfeat)
	VariableFeatures (ref) = vf
	ref <- ScaleData (ref)
	ref <- RunPCA (ref)

	# Run Harmony
	batch = 'lungid'
	reductionSave = 'harmony'
	ref = ref %>% 
	RunHarmony (batch, plot_convergence = FALSE, reduction = 'pca', reduction.save= reductionSave) %>%
	RunUMAP (reduction = reductionSave, dims = 1:15)

	#ref <- RunUMAP (ref, dims = 1:15)
	ref = FindNeighbors (object = ref, reduction = 'harmony', dims = 1:15, k.param = 30,
                              verbose = TRUE, force.recalc = T)
	ref = FindClusters (ref, resolution = 10, verbose = T, n.start = 100)

	pdf ('Plots/ref_samples_umap.pdf', width=12)
	FeaturePlot (ref, feature = c('VWF','PECAM1','PLVAP','CLDN5'), reduction = 'umap')
	DimPlot (ref, group.by = 'seurat_clusters', reduction = 'umap', label=T)
	DimPlot (ref, group.by = 'celltype', reduction = 'umap')
	dev.off()
	ref$celltype[ref$seurat_clusters == 99] = 'Endothelial'
	ref = ref[,ref$celltype %in% c('Mesothelium','Fibroblast','SmoothMuscle','Endothelial')]

	#ref = ref[,ref$lungid %in% c('HU37','HU62')]
	ref$sampleID = ref$lungid
	srt = merge (srt, ref)
	
	srt[['SCT']] = NULL
	srt[['integrated']] = NULL	
	
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
	srt <- RunUMAP (srt, dims = 1:15)
	srt = FindNeighbors (object = srt, reduction = 'pca', dims = 1:15, k.param = 30,
                              verbose = TRUE, force.recalc = T)
	srt = FindClusters (srt, resolution = 1, verbose = T, n.start = 100)

	srt[['SCT']] = NULL
	srt[['integrated']] = NULL	
	srt$celltype [srt$celltype == 'PLVAP'] = 'Endothelial'
	srt$celltype [srt$celltype == 'Vein'] = 'Endothelial'
	srt$celltype [srt$celltype == 'Artery'] = 'Endothelial'
	srt$celltype [srt$celltype == 'Fibroblast'] = 'Fibroblasts'
	saveRDS (srt, 'srt.rds')
	
	pdf ('Plots/umap_samples.pdf', width=12)
	FeaturePlot (srt, feature = c('VWF','PECAM1','ACTA2','COL1A1'))
	wrap_plots (DimPlot (srt, group.by = 'sampleID'), DimPlot (srt, group.by = 'seurat_clusters', label=T),
		DimPlot (srt, group.by = 'celltype'))
	dev.off()
	} else {
	srt = readRDS ('srt.rds')	
	}
