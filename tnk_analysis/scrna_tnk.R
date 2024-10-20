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

set.seed(1234)

# Set project dir
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scrna'
dir.create (file.path(projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

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
	srt = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/srt.rds')
	
	# Import P14
	srt_p14 = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/MPM_naive_p14_analysis/_cellranger_raw_Filter_400_1000_25/no_harmony/srt.rds')
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
	srt = FindClusters (srt, resolution = 1, verbose = T, n.start = 100)

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
srt$celltype2[srt$RNA_snn_res.0.8 == 11] = 'Tregs'
srt$celltype2[srt$RNA_snn_res.0.8 == 12] = 'Proliferating'
srt$celltype2[srt$RNA_snn_res.0.8 == 8] = 'NK_KLRC1'
srt$celltype2[srt$RNA_snn_res.0.8 %in% c(2,9,7,6)] = 'CD8'
srt$celltype2[srt$RNA_snn_res.1 %in% c(6)] = 'CD8_exhausted'
srt$celltype2[srt$RNA_snn_res.0.8 %in% c(0,3,5,1,10)] = 'CD4'
srt$celltype2[srt$RNA_snn_res.0.8 %in% c(4)] = 'NK_FGFBP2'

dp = DimPlot(srt, group.by = 'celltype2')
pdf (file.path('Plots',paste0('celltype_new_dimplot.pdf')))
dp
DimPlot (srt, group.by = 'RNA_snn_res.1', label=T)
dev.off()


saveRDS (srt, 'srt.rds')

} else {
srt = readRDS ('srt.rds')	
}


pdf (file.path('Plots',paste0('Tregs_TFs.pdf')))
VlnPlot (srt, features = c('POU2F3','POU2F2'), group.by = 'celltype2')
dev.off()

pdf (file.path('Plots',paste0('Tregs_TFs.pdf')))
VlnPlot (srt, features = c('PRF1','GZMB','KLRC1'), group.by = 'celltype2')
dev.off()

pdf (file.path('Plots',paste0('NK_T_EXT.pdf')), width=12)
VlnPlot (srt, features = c('SPRY1','SFTPB','NR4A3','NR2F2'), group.by = 'celltype2', split.by='sampleID')
dev.off()


pdf (file.path('Plots',paste0('NK_T_EXT_TFs_featureplot.pdf')))
reductionName = 'umap'
fp (srt, gene = c('NR4A3','SPRY1','NR2F2'))
dev.off()


dp = DotPlot (srt,
  features = c('NR4A1','NR4A2','NR4A3','TOX','TOX2','IRF1'),
  #col = palette_gene_expression2,
  group.by = 'celltype2') + gtheme_italic
pdf (file.path('Plots','TNK_exhaustion_markers_dotplot.pdf'), height=2.6, width=4.6)
dp
dev.off()  

dp = DotPlot (srt,
  features = c('GZMB', 'PRF1', 'GZMH', 'GZMK', 'KLRC1','NCAM1'),
  #col = palette_gene_expression2,
  group.by = 'celltype2') + gtheme_italic
pdf (file.path('Plots','TNK_exhaustion_markers_dotplot.pdf'), height=2.6, width=5)
dp
dev.off()  

dp = DotPlot (srt,
  features = c('NR4A2','JUND','RUNX3','EOMES','JUN','JUNB','FOSB','REL','IRF1'),
  #col = palette_gene_expression2,
  group.by = 'celltype2') + gtheme_italic
pdf (file.path('Plots','TNK_exhaustion_markers_dotplot.pdf'), height=2.6, width=5)
dp
dev.off()  


# Run correlation to target TF to identify possible regulators ####
library (hdWGCNA)
selected_TF = 'NR4A2'

if (!file.exists ('metacells.rds'))
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
	  group.by = c("celltype2"), # specify the columns in seurat_obj@meta.data to group by
	  reduction = 'umap', # select the dimensionality reduction to perform KNN on
	  k = 50, # nearest-neighbors parameter
	  max_shared = 25, # maximum number of shared cells between two metacells
	  ident.group = 'celltype2' # set the Idents of the metacell seurat object
	)
	
	# normalize metacell expression matrix:
	srt <- NormalizeMetacells (srt)
	metacells = GetMetacellObject (srt)
	saveRDS (metacells, 'metacells.rds')
	} else {
	metacells = readRDS ('metacells.rds')	
	}

ct = c('CD8','CD8_exhausted','CD4','NK_FGFBP2','NK_KLRC1','Tregs')
selected_TF = 'NR4A2'
cor_ct = lapply (ct, function(x) 
	{
	metacells_assay = metacells@assays$RNA@layers$data
	rownames (metacells_assay) = rownames(srt)
	metacells_assay_sample = metacells_assay[,metacells$celltype2 == x]
	tc_cor = t(cor (metacells_assay_sample[selected_TF,], t(metacells_assay_sample), method='spearman'))
	tc_cor = setNames (tc_cor[,1], rownames (tc_cor))
	tc_cor[is.na(tc_cor)] = 0
	head (tc_cor,5000)
	#tc_cor
	tc_cor[order(-tc_cor)]
	})
names(cor_ct) = ct
head (cor_ct[['CD8_exhausted']]['CTLA4'],100)
head (cor_ct[['NK_KLRC1']]['RUNX3'],100)

cor_ct_df = do.call (cbind, cor_ct)
cor_ct_df = cor (cor_ct_df)
rownames(cor_ct_df) = ct
colnames(cor_ct_df) = ct
cor_ct_df

pdf (file.path ('Plots',paste0(selected_TF,'_cor.pdf')))

# Correlate to inhibitory molecules in CD8 and NK  ####
# Import TFs from ATAC
tfs = readRDS (file.path ('..','scatac_ArchR','Annotations','Motif-Matches-In-Peaks.rds'))
tfs = colnames(tfs)
tfs = gsub ('_.*','',tfs)
tfs = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", tfs)

selected_TF = 'PDCD1'
ct = 'CD8_exhausted'
cor_cd8 = lapply (ct, function(x) 
	{
	metacells_assay = metacells@assays$RNA@layers$data
	rownames (metacells_assay) = rownames(srt)
	metacells_assay_sample = metacells_assay[,metacells$celltype2 == x]
	tc_cor = t(cor (metacells_assay_sample[selected_TF,], t(metacells_assay_sample), method='spearman'))
	tc_cor = setNames (tc_cor[,1], rownames (tc_cor))
	tc_cor[is.na(tc_cor)] = 0
	head (tc_cor,5000)
	#tc_cor
	tc_cor = tc_cor[tfs]
	tc_cor[order(-tc_cor)]
	})
head (cor_cd8[[1]],100)
cor_cd8[[1]]['NR4A2']
cor_cd8[[1]] = na.omit (cor_cd8[[1]])

selected_TF = 'KLRC1'
ct = 'NK_KLRC1'
cor_nk = lapply (ct, function(x) 
	{
	metacells_assay = metacells@assays$RNA@layers$data
	rownames (metacells_assay) = rownames(srt)
	metacells_assay_sample = metacells_assay[,metacells$celltype2 == x]
	tc_cor = t(cor (metacells_assay_sample[selected_TF,], t(metacells_assay_sample), method='spearman'))
	tc_cor = setNames (tc_cor[,1], rownames (tc_cor))
	tc_cor[is.na(tc_cor)] = 0
	head (tc_cor,5000)
	#tc_cor
	tc_cor = tc_cor[tfs]
	tc_cor[order(-tc_cor)]
	})
cor_nk[[1]]
cor_nk[[1]]['NR4A2']
cor_nk[[1]] = na.omit (cor_nk[[1]])
intersect (head (names(cor_nk[[1]]),100), head (names(cor_cd8[[1]]),100))

# Export table
nk_cd8_ext_cor = data.frame (
	row.names = names(cor_cd8[[1]]),
	cd8 = cor_cd8[[1]],
	nk = cor_nk[[1]][names(cor_cd8[[1]])])
write.csv (nk_cd8_ext_cor, 'nk_cd8_ext_cor.csv')

nk_cd8_ext_cor$mean_cor = rowMeans (nk_cd8_ext_cor)
head (nk_cd8_ext_cor[order(-nk_cd8_ext_cor$mean_cor),],20)




#### Run cNMF ####
nfeat = 5000
force=F
k_list = c(5:30)
k_selections = c(5:30)
cores= 100
cnmf_name = 'TNK'
cnmf_out = paste0('cNMF/cNMF_',cnmf_name,'_',paste0(k_list[1],'_',k_list[length(k_list)]),'_vf',nfeat)
dir.create (file.path(cnmf_out,'Plots'), recursive=T)
repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'

### RUN consensus NMF ####
source (file.path ('..','..','git_repo','utils','cnmf_prepare_inputs.R')) 


### Import and format spectra files ####
k_selection = 19
source (file.path ('..','..','git_repo','utils','cnmf_format_spectra_files.R')) 
sapply (cnmf_spectra_unique, function(x) 'KLRC1' %in% x)
cnmf_spectra_unique[[6]]


### Run pathway enrichment on modules

#### Run GSEA enrichment on each cluster DE markers ####
genesets = lapply(cnmf_spectra_unique, function(x) head (x, 100))
enricher_universe = VariableFeatures(srt)
gmt_annotations = c(
#'c2.cp.kegg.v7.1.symbol.gmt',
#'c5.bp.v7.1.symbol.gmt',
'h.all.v7.4.symbols.gmt'
#'h.all.v7.1.symbol.gmt'
)
force=T
if (!file.exists (file.path('Pathway_Enrichment_clusters.rds')) | force)
	{
	# GSEA analysis on DEG per cluster
	EnrichRResAll = list()
	for (ann in gmt_annotations)
		{
		gmt.file = file.path ('..','..','git_repo','files', ann)
		pathways = read.gmt (gmt.file)
		#pathways = gmtPathways (gmt.file)
		message (paste('Compute enrichment per cluster using annotation:', ann))
		EnrichRResCluster = list()
		for (i in seq_along(genesets))
			{
			message (paste ('EnrichR running geneset',i))	
			geneset = genesets[[i]]
			egmt <- enricher(geneset, TERM2GENE=pathways, universe = enricher_universe)
			EnrichRResCluster[[i]] = egmt@result
			}
		EnrichRResAll[[ann]] = EnrichRResCluster
		}
	names(EnrichRResAll) = gmt_annotations
	saveRDS (EnrichRResAll, file.path('Pathway_Enrichment_clusters.rds'))	
	} else {
	EnrichRResAll = readRDS (file.path('Pathway_Enrichment_clusters.rds'))
	}

# Plot fgsea enrichments
pvalAdjTrheshold = 0.05
top_pathways = 50
EnrichRes_dp = lapply (EnrichRResAll, function(x) dotGSEA (
	enrichmentsTest_list = x, 
	type = 'enrich', 
	padj_threshold = pvalAdjTrheshold, 
	top_pathways= top_pathways,
	cluster_rows=T,
	cluster_cols=T))
#gmt_annotations = 'c5.bp.v7.1.symbol.gmt'
for (ann in gmt_annotations)
	{
	pdf (file.path('Plots', paste0('Pathway_Enrichment_',ann,'_dotplot.pdf')), width=8, height=6)
	print (EnrichRes_dp[[ann]])
	dev.off()
	}







# Generate heatmap of ext TF with associated spectra ####
# Import ext TFs
ext_tfs = read.csv (file.path('..','scatac_ArchR','top_TF_CD8_NK_dual_ext_TF_activity_RNA.csv'))

pdf (file.path('Plots','tf_ext_spectra_heatmap.pdf'), width=4)
Heatmap (cnmf_spectra[rownames(cnmf_spectra) %in% ext_tfs[,1],],
	col =palette_expression, 
	row_names_gp = gpar(fontsize = 5),
	border=T)#, row_names = gpar (fontsize=5))
dev.off()











srt_orig = srt
#### Run cNMF on NK KLRC1+ and CD8 exhausted ####
srt = srt_orig[, srt_orig$celltype2 == 'NK_KLRC1' & srt_orig$sampleID %in% 'P1']
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scrna'
setwd (projdir)
dir.create (file.path('NK_KLRC1','Plots'), recursive=TRUE)
setwd ('NK_KLRC1')
projdir = file.path (projdir,'NK_KLRC1')

#### Run cNMF ####
nfeat = 5000
force=F
k_list = c(5:30)
k_selections = c(5:30)
cores= 100
cnmf_name = 'NK_KLRC1'


repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'

### RUN consensus NMF ####
source (file.path ('..','..','..','git_repo','utils','cnmf_prepare_inputs.R')) 

### Import and format spectra files ####
k_selection = 5
source (file.path ('..','..','..','git_repo','utils','cnmf_format_spectra_files.R')) 



#### Run cNMF on CD8 exhausted ####
srt = srt2
srt2 = srt
srt = srt2[, srt$celltype2 == 'CD8_exhausted' & srt$sampleID %in% 'P1']
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scrna'
setwd (projdir)
projdir = file.path (projdir,'CD8_exhausted')
dir.create (file.path('CD8_exhausted','Plots'), recursive=TRUE)
setwd ('CD8_exhausted')

#### Run cNMF ####
nfeat = 5000
force=F
k_list = c(5:30)
k_selections = c(5:30)
cores= 100
cnmf_name = 'CD8_exhausted'
cnmf_out = paste0('cNMF/cNMF_',cnmf_name,'_',paste0(k_list[1],'_',k_list[length(k_list)]),'_vf',nfeat)
dir.create (file.path(cnmf_out,'Plots'), recursive=T)
repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'

### RUN consensus NMF ####
source (file.path ('..','..','..','git_repo','utils','cnmf_prepare_inputs.R')) 

### Import and format spectra files ####
k_selection = 5
source (file.path ('..','..','..','git_repo','utils','cnmf_format_spectra_files.R')) 




# Compare cNMFs ####
k_selection = 10
cnmf_nk = readRDS (paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scrna/NK_KLRC1/cnmf_genelist_',k_selection,'_nfeat_5000.rds'))
names (cnmf_nk) = paste0('nk_',names(cnmf_nk))
sapply (cnmf_nk, function(x) 'NR4A2' %in% x)
cnmf_cd8ext = readRDS (paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scrna/CD8_exhausted/cnmf_genelist_',k_selection,'_nfeat_5000.rds'))
names (cnmf_cd8ext) = paste0('cd8ext_',names(cnmf_cd8ext))
sapply (cnmf_cd8ext, function(x) 'NR4A2' %in% x)

cnmf_nk = lapply (cnmf_nk, function(x) head (x, 300))
cnmf_cd8ext = lapply (cnmf_cd8ext, function(x) head (x, 300))
cnmf_cd8ext[[6]][cnmf_cd8ext[[6]] %in% cnmf_nk[[4]]]
ov_mat = ovmat (c(cnmf_nk, cnmf_cd8ext), compare_lists = list(names(cnmf_nk),names(cnmf_cd8ext)))

projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scrna'
setwd (projdir)
pdf (file.path ('Plots','compare_cnmfs_NK_KLRC1_CD8exhausted.pdf'),5,5)
ov_mat
dev.off()


srt = srt2
dp = DotPlot (srt,
  features = c('NR4A2','SPRY1','CTNNB1','NFATC1','ID2','NFKB2'),
  #col = palette_gene_expression2,
  group.by = 'celltype2') + gtheme_italic
pdf (file.path('Plots','TNK_exhaustion_markers_dotplot.pdf'), height=2.6, width=4.6)
dp
dev.off()  

dp = DotPlot (srt,
  features = c(rownames(srt)[grep('RUNX', rownames(srt))]),
  #col = palette_gene_expression2,
  group.by = 'celltype2') + gtheme_italic
pdf (file.path('Plots','RUNX_dotplot.pdf'), height=2.6, width=7.6)
dp
dev.off()  


