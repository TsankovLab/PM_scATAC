use UGER # Add this before running R to be able to run cNMF scripts using UGER 
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
library (infercnv)

set.seed(1234)

# Set project dir
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/main/scrna/'
dir.create (paste0(projdir,'/Plots/'), recursive =T)
setwd (projdir)
source ('../../PM_scATAC/useful_functions.R')
source ('../../PM_scATAC/palettes.R')
source ('../../PM_scATAC/ggplot_aestetics.R')
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
srt = readRDS ('srt.rds')	

### CNV analysis using inferCNV ####

library (gtools)
# Variables needed #
#org = 'mouse'
#projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/Romain_prj/' # project folder (use absolute path) where you have your seurat object saved as RDS
projdir_out = paste0('infercnv/') # output folder 
dir.create (projdir_out, showWarnings = FALSE)
# srt_obj  = 'seurat.obj.rds' # Specify the name of your seurat object saved as RDS 
#metaGroupName = 'seurat_clusters' # Specify the column in your meta data referring to the clustering (usually is 'seurat_clusters')
#metaGroupName2 = 'samples'
#cancer_clusters = c('0','4','6','7','8','11','14')
# cancer_clusters  = c("0", "1", "3", "9", "13", "18", "25", "26") # specify cancer clusters e.g. c('1','2','3')
# ref_clusters   = c("2", "4", "5", "6","7", "8", "10", "11","12", "14", "15", "16", "17","19", "20","21","22","23","24","27") # specify clusters to use as references e.g. c('5','6')
# subsample  = NULL
#chr_exclude = c('chrX','chrY', 'chrM')

#system (paste('mkdir -p',projdir_out))

# subset for malignant cells called from manual annotation and label trasnfer
#srt = srt [,srt$celltype  == 'Malignant' & srt$predicted.id %in% c('Sarcomatoid','Malignant')]
srt_cnv = srt[,srt$sampleID == 'P1']
srt_cnv$cnv_type = 'sample'

#load('/ahg/regevdata/projects/ICA_Lung/Maggie_fast/infercnv/10x_HU37_ResEpi.Rda')
#ResEpi<-UpdateSeuratObject(ResEpi)
#data.ref <- as.matrix(GetAssayData(object = ResEpi, slot = "counts")[,colnames(ResEpi)[seq(1,length(colnames(ResEpi)),20)]])

#ResEpi_sub <- ResEpi[,colnames(ResEpi)[seq(1,length(colnames(ResEpi)),20)]]
#data.ref <- as.matrix(GetAssayData(object = ResEpi, slot = "counts")[,colnames(ResEpi)[seq(1,length(colnames(ResEpi)),20)]])
ref = readRDS ('../../tumor_compartment/scrna/scRNA_meso.rds')
ref = ref[,ref$sampleID %in% c('HU37','HU62')]
ref$cnv_type = 'reference'
srt_cnv = merge (srt_cnv, ref)

#srt_cnv = merge (malig, ResEpi_sub)
srt_cnv = NormalizeData (object = srt_cnv, normalization.method = "LogNormalize")
srt_cnv = FindVariableFeatures (srt_cnv)
srt_cnv = ScaleData (srt_cnv)
srt_cnv = RunPCA (srt_cnv, npcs = 30, ndims.print = 1:5, nfeat.print = 5, verbose = FALSE)
srt_cnv = RunUMAP (srt_cnv, dims = 1:20)

# srt_cnv_obj  = 'seurat.obj.rds' # Specify the name of your seurat object saved as RDS 
  
#srt_cnv = readRDS (paste0(projdir,srt_cnv_obj)) # load seurat object 

# srt_cnv$malig_ref = ifelse (srt_cnv@meta.data[,metaGroupName] %in% cancer_clusters, 'malignant','reference')
message ('generate UMAP of malignant and reference clusters')
dim_p1 = DimPlot (srt_cnv, group.by = 'cnv_type')
png (paste0(projdir_out,'check_malig_clusters_ref_umap.png'), width=1500)
print (wrap_plots (dim_p1, ncol=3))
dev.off()

# load TxDbi object for the species and get genomic positions of genes
require (org.Hs.eg.db)
require (TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(txdb) = paste0('chr',1:22) # select chromosomes to include in the analysis
gene_regions = as.data.frame (genes (txdb))

# map entrez id to symbol ids
symbol = toTable (org.Hs.egSYMBOL)
gene_regions$symbol = symbol$symbol[match(gene_regions$gene_id, symbol$gene_id)]
gene_regions = gene_regions[gene_regions$symbol %in% rownames(srt_cnv),]

#Prepare files for inferCNV input
message ('Generate expression and annotation files')
#srt_cnv$orig.ident = srt_cnv@meta.data [,metaGroupName2]
malig_samples = unique(srt_cnv$cnv_type[srt_cnv$cnv_type != 'reference'])
        
exprMat = srt_cnv@assays$RNA@counts        
gene_regions2 = gene_regions[mixedorder(gene_regions$seqnames), ]
exprMat = exprMat[gene_regions2$symbol, ]

all (rownames(exprMat) == gene_regions2$symbol)
rownames (gene_regions2) = gene_regions2$symbol
colnames (gene_regions2) = NULL
message ('save gene_regions file')    
write.table (gene_regions2, paste0(projdir_out,'gene_regions.txt'), sep='\t')
#exprMat = exprMat[,rownames(annot_df)]
annot_df = srt_cnv@meta.data[,'cnv_type', drop=F]
annot_df[,1] = srt_cnv$cnv_type

colnames (annot_df) = NULL
message ('save annotation file')    
write.table (annot_df, paste0(projdir_out,'annot_df.txt'), sep='\t', col.names=NA)


library (biomaRt)

chr_exclude = c('chrX','chry')
infercnv_obj = CreateInfercnvObject (raw_counts_matrix=exprMat,
                                    annotations_file= paste0(projdir_out,'annot_df.txt'),
                                    delim="\t",
                                    gene_order_file=paste0(projdir_out,'gene_regions.txt'),
                                    ref_group_names=c('reference'),
                                    chr_exclude= chr_exclude
                                    )

infercnv_result = infercnv::run(infercnv_obj,
                           cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                           out_dir=projdir_out,  # dir is auto-created for storing outputs
                           cluster_by_groups=T,  # cluster
                           denoise=T,
                           HMM=F,
                           save_rds = F,
                           no_prelim_plot = T,
                           no_plot = T,
                           plot_probabilities = FALSE
                           )

#if (!is.null(custom_color_pal)) custom_color_pal = infercnv_pal
# plot_cnv (infercnv_result, out_dir = projdir_out_sam#, #title = paste0('inferCNV_'), 
  
#   #output_filename=paste0('inferCNV_',sam,'_subsample_',subsample)
#   )
saveRDS (infercnv_result, paste0(projdir_out, 'infercnv.results.obj.Rds'))
  
infercnv_result = readRDS (paste0(projdir_out, 'infercnv.results.obj.Rds'))
plot_cnv (infercnv_result, out_dir = projdir_out)#, #title = paste0('inferCNV_'), 

# Compute CNV load and identify sarcomatoid cluster 
cnv_mat = infercnv_result@expr.data
cnv_mat[cnv_mat > 1] = 2
cnv_load = apply (infercnv_result@expr.data,2, function(sum (abs(x)))

message (paste('read infercnv object of sample', x))    

lexpr.data_sample = log2(infercnv_result@expr.data) 
lexpr.data_sample = lexpr.data_sample[,-infercnv_result@reference_grouped_cell_indices[[1]]]
    #lexpr.data_sample = lexpr.data
cnvload = as.vector(colSums(abs(lexpr.data_sample)))
qnt = .7
cnvloadn = cnvload/quantile(cnvload, qnt)
cnvloadn = setNames (cnvloadn, colnames(lexpr.data_sample))
srt_cnv$cnv_load = 0
srt_cnv$cnv_load[match(names(cnvloadn),colnames(srt_cnv))] = cnvloadn
    
reductionName = 'umap'
srt_cnv <- FindNeighbors (object = srt_cnv)
srt_cnv <- FindClusters (object = srt_cnv, verbose = FALSE)
dim_p = DimPlot (object = srt_cnv, label = TRUE) + NoLegend()
pdf (paste0(projdir_out,'cnv_load.pdf'))
fp (srt_cnv, gene = 'cnv_load')
dim_p
dev.off()



cor.data = cor (lexpr.data_sample, rowMeans (lexpr.data_sample[, colnames(srt_cnv)[srt_cnv$seurat_clusters == '14']]))
srt_cnv$cnv_cor = setNames (cor.data[,1], rownames(cor.data))

pdf (paste0(projdir_out,'cnv_cor.pdf'))
    fp (srt_cnv, gene = 'cnv_cor')
    dev.off()


# subset for cells correlated to malignant cnv and show heatmap
cor_threshold = 0.5
high_cor_barcodes = list (
fibro = na.omit (colnames(srt_cnv)[srt_cnv$cnv_cor > cor_threshold & srt_cnv$seurat_clusters %in% c('13','5','2')]),
sarco = na.omit (colnames(srt_cnv)[srt_cnv$cnv_cor > cor_threshold & srt_cnv$seurat_clusters %in% c('14')]),
mesothelium = na.omit (colnames(srt_cnv)[srt_cnv$cnv_cor > cor_threshold & srt_cnv$seurat_clusters %in% c('12')]))

infercnv_result_subset = infercnv_result
infercnv_result_subset@expr.data = infercnv_result_subset@expr.data[,unlist(high_cor_barcodes)]

ha = HeatmapAnnotation (cnv = rep (names(high_cor_barcodes), sapply (high_cor_barcodes, function(x) length(x))), which='row')
pdf (paste0(projdir_out,'cnv_subset_heatmap.pdf'))
Heatmap (
	t(infercnv_result_subset@expr.data), 
	left_annotation = ha, 
	cluster_columns=F,
	row_names_gp = gpar(fontsize = 2),
	column_names_gp = gpar(fontsize = 0)
	)
dev.off()

# Map on umap the highly correlated cells to malignant cells
cnv_fibro_high = c(
'p786neg_GGTGTTACATACTTTC-1',
'p786neg_CTAGGTATCCATATGG-1',
'p786neg_TACTTCACAAGCTGTT-1'
)
srt_cnv$cnv_high_cor = 0
srt_cnv$cnv_high_cor[match (cnv_fibro_high, colnames(srt_cnv))] = 'high_cnv_cor_fibro'
pdf (paste0(projdir_out,'cnv_subset_umap.pdf'))
DimPlot (srt_cnv, group.by = 'cnv_high_cor', cells.highlight = cnv_fibro_high)
dev.off()



