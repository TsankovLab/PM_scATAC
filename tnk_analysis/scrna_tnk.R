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
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scrna'
dir.create (file.path(projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','PM_scATAC','utils','load_packages.R'))
source (file.path('..','..','PM_scATAC','utils','useful_functions.R'))
source (file.path('..','..','PM_scATAC','utils','ggplot_aestetics.R'))
source (file.path('..','..','PM_scATAC','utils','scATAC_functions.R'))
source (file.path('..','..','PM_scATAC','utils','palettes.R'))

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
	srt = srt[, srt$celltype_simplified %in% c('T_cells','NK')]
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
	srt = FindClusters (srt, resolution = .8, verbose = T, n.start = 100)

	reductionName = 'umap'
	fps = fp (srt, gene = c('CD3D','KLRC1','FGFBP2','CD8A','CD4','FOXP3','CXCL13','CTLA4','PDCD1','CXCR5','MKI67'))
	dp = wrap_plots (DimPlot (srt, group.by = 'sampleID'), DimPlot (srt, group.by = 'seurat_clusters', label=T),
		DimPlot (srt, group.by = 'celltype'))
	pdf (file.path('Plots','markers_celltypes_umap.pdf'), width=12)
	wrap_plots (fps)
	dev.off()
	pdf (file.path('Plots','celltypes_umap.pdf'),width=15)
	wrap_plots (dp)
	dev.off()
	

metaGroupName = 'seurat_clusters'
if (file.exists ('DEG.rds'))
	{
	DefaultAssay(srt) = 'RNA'
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
png (file.path('Plots',paste0('top_',top_genes,'_heatmap.png')), height=15500, width=3300, res=400)
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


# Annotate cells ####
srt$celltype2 = 0
srt$celltype2[srt$seurat_clusters == 11] = 'Tregs'
srt$celltype2[srt$seurat_clusters == 12] = 'Proliferating'
srt$celltype2[srt$seurat_clusters == 8] = 'NK_KLRC1'
srt$celltype2[srt$seurat_clusters == 8] = 'NK_KLRC1'
srt$celltype2[srt$seurat_clusters %in% c(2,9,7,6)] = 'CD8'
srt$celltype2[srt$seurat_clusters %in% c(0,3,5,1,10)] = 'CD4'
srt$celltype2[srt$seurat_clusters %in% c(4)] = 'NK_FGFBP2'

dp = DimPlot(srt, group.by = 'celltype2')
pdf (file.path('Plots',paste0('celltype_new_dimplot.pdf')))
dp
dev.off()


saveRDS (srt, 'srt.rds')

} else {
srt = readRDS ('srt.rds')	
}


pdf (file.path('Plots',paste0('Tregs_TFs.pdf')))
VlnPlot (srt, features = c('POU2F3','POU2F2'), group.by = 'celltype2')
dev.off()







