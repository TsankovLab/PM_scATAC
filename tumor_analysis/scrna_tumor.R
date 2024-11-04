conda activate meso_scatac
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
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scrna/'
dir.create (paste0(projdir,'/Plots/'), recursive =T)
setwd (projdir)
source ('../../git_repo/utils/useful_functions.R')
source ('../../git_repo/utils/palettes.R')
source ('../../git_repo/utils/ggplot_aestetics.R')
source ('../../git_repo/utils/palettes.R')

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
	#srt = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/reproduction2/scRNA/srt_tumor.rds')
	srt = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/srt.rds')
	srt$celltype_simplified_mal = as.character (srt$celltype_simplified)
	srt$celltype_simplified_mal[srt$celltype_simplified == 'Malignant'] = paste0(srt$sampleID[srt$celltype_simplified == 'Malignant'], '_', srt$celltype_simplified[srt$celltype_simplified == 'Malignant'])
	# Import P14
	srt_p14 = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_p14_analysis/_cellranger_raw_Filter_400_1000_25/no_harmony/srt.rds')
	srt_p14$celltype_simplified = srt_p14$celltype
	srt = merge (srt, srt_p14)
	srt$sampleID[srt$sampleID == 'MPM_naive_p14'] = 'P14'
	srt = srt[, srt$celltype_simplified == 'Malignant']
	srt = srt[, srt$sampleID %in% sample_names]
	
	# Import normal mesothelium RNA
	projdir_ref = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scRNA_meso/'
	ref = readRDS (paste0(projdir_ref, 'ref.rds'))	
	#subsample_ct = 1000
	#ref_bc = lapply (unique(ref$celltype), function(x) sample(colnames(ref)[ref$celltype == x], subsample_ct, replace=T))
	#names (ref_bc) = unique(ref$celltype)
	#cts = c('Mesothelium','Fibroblast')
	#ref = ref[,colnames (ref) %in% unlist(ref_bc[cts])]
	ref = ref[,ref$celltype == 'Mesothelium']
	ref = ref[,ref$lungid %in% c('HU37','HU62')]
	ref$sampleID = ref$lungid
	ref$celltype_simplified = 'Normal_mesothelium'
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
	saveRDS (srt, 'srt.rds')
	
	pdf ('Plots/umap_samples.pdf', width=12)
	wrap_plots (DimPlot (srt, group.by = 'sampleID'), DimPlot (srt, group.by = 'seurat_clusters', label=T))
	dev.off()
	} else {
	srt = readRDS ('srt.rds')	
	}


#### Import bulk RNA subtype signatures to define malignant score ####
# Bueno ####
figures.dir = '/ahg/regevdata/projects/lungCancerBueno/Results/TCGAsubtypes/LUAD_LUSC_MESO/Bueno/v2/'
markers.res1.2 = get(load(paste0(figures.dir,'markers.res1.2.Rda'))) 
markers.res1.2.trimmed = markers.res1.2[markers.res1.2$avg_logFC > 0.2,]

top_bueno_genes = markers.res1.2.trimmed %>% group_by(cluster) %>% top_n(20, avg_logFC)

#top20genes = markers.res1.2.trimmed.2 %>% group_by(cluster) %>% top_n(10, avg_logFC)
top_bueno_genes = split (top_bueno_genes$gene, top_bueno_genes$cluster)

if (!all (names(top_bueno_genes) %in% colnames (srt@meta.data)))
	{
	srt = ModScoreCor (
	        seurat_obj = srt, 
	        geneset_list = top_bueno_genes, 
	        cor_threshold = NULL, 
	        pos_threshold = NULL, # threshold for fetal_pval2
	        listName = 'Bueno_', outdir = paste0(projdir,'Plots/'))

	}
  
#### Run cNMF ####
nfeat = 5000
force=F
k_list = c(5:30)
k_selections = c(5:30)
k_selection = 25
cores= 100

cnmf_name = 'scrna_tumor'
repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'

### RUN consensus NMF ####
source (file.path ('..','..','git_repo','utils','cnmf_prepare_inputs.R')) 

### Import and format spectra files ####
k_selection = 25
cnmf_name = 'scrna_tumor'
source (file.path ('..','..','git_repo','utils','cnmf_format_spectra_files.R')) 

top_nmf_genes = 50
cnmf_spectra_unique = lapply (cnmf_spectra_unique, function(x) head (x, top_nmf_genes))

# Add module score of cNMF modules ####
if (!all (names(cnmf_spectra_unique) %in% colnames (srt@meta.data)))
	{
	srt = ModScoreCor (
	        seurat_obj = srt, 
	        geneset_list = cnmf_spectra_unique, 
	        cor_threshold = NULL, 
	        pos_threshold = NULL, # threshold for fetal_pval2
	        listName = 'cNMF_', outdir = paste0(projdir,'Plots/'))

	}

srt$cNMF = 'cNMF'
ccomp_df = srt@meta.data[,c(names(cnmf_spectra_unique),'sampleID'), drop=FALSE]
      #ccomp_df = aggregate (ccomp_df, by=as.list(srt_wgcna@meta.data[,metaGroupNames,drop=F]), mean)    
bp1 = lapply (names(cnmf_spectra_unique), function(x) {
            ggplot (ccomp_df, aes_string (x= 'sampleID', y= x)) +
        #geom_violin (trim=TRUE, aes_string (fill = metaGroupNames[3])) +
        #geom_violin (aes_string(fill = metaGroupNames[3])) +
        geom_boxplot(width=0.5, color="black", alpha=0.2) +
        #geom_bar (stats='identity') +
        #geom_jitter (color="black", size=0.4, alpha=0.9) +
        theme_classic() + 
        #scale_fill_manual (values= module_pal) + 
        ggtitle (paste(x,'mod score')) + 
        theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + NoLegend()
      })  
png (file.path(cnmf_out,'Plots',paste0('cNMF_module_scores_boxplots_',k_selection,'nfeat_',nfeat,'.png')),5000,5000,res=300)
print (wrap_plots (bp1))
dev.off()


### FIGURE 3B - Plot boxplots of sarcomatoid cnmf ####
ccomp_df = srt@meta.data
ccomp_df = aggregate (ccomp_df[,c('cNMF20')], by=as.list(srt@meta.data[,'sampleID',drop=F]), 'mean')
rownames(ccomp_df) = ccomp_df[,1]
ccomp_df = ccomp_df[,-1, drop=F]
srt$sampleID = factor (srt$sampleID, levels = rownames(arrange (ccomp_df[,'x', drop=FALSE], -ccomp_df[,'x', drop=FALSE])))

ccomp_df = srt@meta.data
box = ggplot (ccomp_df, aes_string (x= 'sampleID', y= 'cNMF20')) +
  geom_violin (trim=TRUE, aes_string (fill = 'sampleID'),size=2,
    width=1, scale='width',
    linewidth = .2, alpha=0.7) +
  geom_boxplot (aes_string(fill = 'sampleID'),
    linewidth = .2,
    width=0.2,
    outlier.alpha = 0.2,
    outlier.size = .5,
     size=0.3, alpha=0.7
     ) +
  scale_fill_manual (values= palette_sample) +
  gtheme
  
png(file.path('Plots','sarc_signatures_order_boxplot.png'),1200,600, res=300) #width = 10, height = 11,
box
dev.off()



#### FIGURE S2A - Generate Neftel diagram using the four subtypes from bueno ####
# Generate GBM cell states as in paper Nefter et al
# On y axis = D = max(SCopc,SCnpc) - max(SCac,SCmes)
max_sarc = pmax(srt$Sarcomatoid, srt$`Biphasic-S`)
max_epit = pmax(srt$Epithelioid, srt$`Biphasic-E`)
srt$y_axis <- log2(abs(max_sarc - max_epit) + 1)
srt$y_axis[max_epit > max_sarc] <- -1 * srt$y_axis[max_epit > max_sarc]

srt$x_axis = 0
srt$x_axis[srt$y_axis > 0] = log2(abs(srt$Sarcomatoid - srt$`Biphasic-S`) + 1)[srt$y_axis > 0]
srt$x_axis[srt$y_axis > 0 & srt$Sarcomatoid > srt$`Biphasic-S`] <- -1 * srt$x_axis[srt$y_axis > 0 & srt$Sarcomatoid > srt$`Biphasic-S`]

srt$x_axis[srt$y_axis < 0] = log2(abs(srt$Epithelioid - srt$`Biphasic-E`) + 1)[srt$y_axis < 0]
srt$x_axis[srt$y_axis < 0 & srt$`Biphasic-E` > srt$Epithelioid] <- -1 * srt$x_axis[srt$y_axis < 0 & srt$`Biphasic-E` > srt$Epithelioid]

p2l = lapply (levels (srt$sampleID), function(x) 
  {    
  df = srt@meta.data[srt$sampleID == x,]
  tot_cells = nrow(df)    
  ggplot(df, aes(x=x_axis, y=y_axis, fill=y_axis, color=y_axis)) +
  geom_point(alpha=.8, shape=21, stroke=.05, linewidth=3) +
  xlim (c(-2,2)) + ylim (c(-2,2)) +
  #xlab ('Sarcomatoid - Biphasic-S / Epithelioid - Biphasic-E') +
  #ylab ('Sarcomatoid/Biphasic-S - Epithelioid/Biphasic-E') +#+ scale_fill_viridis() + scale_color_viridis()
  geom_vline(xintercept = 0,linetype = 'dashed', size=.1) +
  geom_hline(yintercept = 0,linetype = 'dashed', size=.1) +
  scale_fill_gradientn(colors=rev(palette_sample)) +
  scale_color_gradientn(colors=rev(palette_sample)) +
  #scale_fill_manual (values= sample_colors) +
  #scale_color_manual (values= sample_colors) +
  annotate("text", x = -1.1, y = 1.8, label = paste0("Sarco ",round(sum(df$x_axis < 0 & df$y_axis > 0) / tot_cells *100,1),'%') , size=2.5) +
  annotate("text", x = 1.2, y = 1.8, label = paste0("Bi-S ",round(sum(df$x_axis > 0 & df$y_axis > 0)/ tot_cells *100,1),'%'), size=2.5) +
  annotate("text", x = -1.2, y = -1.8, label = paste0("Bi-E ",round(sum(df$x_axis < 0 & df$y_axis < 0)/ tot_cells *100,1),'%'), size=2.5) +
  annotate("text", x = 1.2, y = -1.8, label = paste0("Epit ",round(sum(df$x_axis > 0 & df$y_axis < 0)/ tot_cells *100,1),'%'), size=2.5) + 
  xlab('') +
  ylab('') + 
  ggtitle (x) + 
  theme_void() + 
  NoLegend()
  })
png (file.path('Plots','neftel_diagram_on_malignant_cells_per_sample.png'),2200,2200, res=300)
wrap_plots (p2l, ncol=4)
dev.off()

reductionName = 'umap'
  umap_df = data.frame (srt[[reductionName]]@cell.embeddings, srt@meta.data[,c(names(cnmf_spectra_unique))])
  umap_p1 = lapply (names(cnmf_spectra_unique), function(x) ggplot(data = umap_df) + 
  geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], color = x), size = .1) + 
  scale_colour_gradientn (colours = rev(brewer.pal (n = 11, name = "RdBu")),limits=c(-max (abs(umap_df[,x])), max (abs(umap_df[,x])))) +
  ggtitle (x) + 
  #facet_wrap (as.formula(paste("~", metaGroupNames[3]))) + 
  theme_classic() +
  theme_void())
  
png (file.path(cnmf_out,'Plots',paste0('cNMF_module_scores_umaps_',k_selection,'nfeat_',nfeat,'.png')),10000,10000,res=300)
print (wrap_plots (umap_p1))
dev.off()

# Export sample order based on sarcomatoid score (cNMF20)
ccomp_df = srt@meta.data[,'cNMF20']
ccomp_df = aggregate (ccomp_df, by=as.list(srt@meta.data[,'sampleID',drop=F]), 'mean')
ccomp_df = ccomp_df[order (ccomp_df$x),]

write.csv (ccomp_df, 'cnmf20_sarcomatoid_sample_order.csv')

# # Run GSEA enrichment on each cluster DE markers
enricher_universe = vf
#do.fgsea = TRUE
gmt_annotations = c(
'h.all.v7.4.symbol.gmt',#,
'c5.bp.v7.1.symbol.gmt',
'c3.tft.v7.1.symbol.gmt'
)
gmt_annotation = gmt_annotations[3]
if (!file.exists (paste0('cNMF_normalized/',cnmf_out, '/EnrichR_cNMF_module_genes_k_',k_selection,'_top_nmf_genes_',top_nmf_genes,'_ann_',gmt_annotation,'.rds')) | force)
  {
    gmt.file = paste0 ('../../PM_scATAC/files/',gmt_annotation)
    pathways = read.gmt (gmt.file)
    EnrichRResCluster = list()
    for (i in names(cnmf_spectra_unique))
      {
      message (paste ('EnrichR running module',i)) 
      sig_genes = cnmf_spectra_unique[[i]]
      if (!any (sig_genes %in% pathways$gene)) next
      egmt = enricher (sig_genes, TERM2GENE=pathways, universe = enricher_universe)
      egmt@result$ID = stringr::str_trunc (egmt@result$ID, 50)
      EnrichRResCluster[[i]] = egmt@result
      }
    EnrichRResAll = EnrichRResCluster
    }
  saveRDS (EnrichRResAll, paste0(cnmf_out, '/EnrichR_cNMF_module_genes_k_',k_selection,'_top_nmf_genes_',top_nmf_genes,'_ann_',gmt_annotation,'.rds'))
  } else {
  EnrichRResAll = readRDS (paste0(cnmf_out, '/EnrichR_cNMF_module_genes_k_',k_selection,'_top_nmf_genes_',top_nmf_genes,'_ann_',gmt_annotation,'.rds'))
  }
  pvalAdjTrheshold = 0.05
  top_pathways = 10
  EnrichRes_dp = dotGSEA (enrichmentsTest_list = EnrichRResAll, type = 'enrich', padj_threshold = pvalAdjTrheshold, top_pathways= top_pathways)
  pdf (paste0(cnmf_out, '/Plots/EnrichR_',i,'_k_',k_selection,'_top_nmf_genes_',top_nmf_genes,'_ann_',gmt_annotation,'_dotplots.pdf'), width = 8 + length(length (cnmf_spectra_unique))/10, height = 15)
  print (EnrichRes_dp)
  dev.off()




# Make correlation network using TFs from scatac analysis
library (hdWGCNA)

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
	  group.by = c("sampleID"), # specify the columns in seurat_obj@meta.data to group by
	  reduction = 'umap', # select the dimensionality reduction to perform KNN on
	  k = 50, # nearest-neighbors parameter
	  max_shared = 25, # maximum number of shared cells between two metacells
	  ident.group = 'sampleID' # set the Idents of the metacell seurat object
	)
	
	# normalize metacell expression matrix:
	srt <- NormalizeMetacells (srt)
	metacells = GetMetacellObject (srt)
	saveRDS (metacells, 'metacells.rds')
	} else {
	metacells = readRDS ('metacells.rds')	
	}

# metacells = metacells[activeTFs[[2]], ]
# cor_TF_l = list()
# for (sam in unique(metacells$sampleID))
#   {
#   cor_TF_l[[sam]] = cor (t(as.matrix(metacells@assays$RNA@layers$data[,metacells$sampleID == sam])))
#   rownames (cor_TF_l[[sam]]) = activeTFs[[2]]
#   colnames (cor_TF_l[[sam]]) = activeTFs[[2]]
#   cor_TF_l[[sam]][is.na(cor_TF_l[[sam]])] = 0
#   cor_TF_l[[sam]] = Heatmap (cor_TF_l[[sam]], name = sam,
#     row_names_gp = gpar(fontsize = 5),
#     column_names_gp = gpar(fontsize = 5))
#   }

# pdf (paste0 ('Plots/selected_TF_exp_corr_heatmaps.pdf'), width = 8,height=9)
# cor_TF_l
# dev.off()

#
#metacells = GetMetacellObject (srt)

### Find enriched pathways in genes correlated with selected TFs ####
selected_TF = readRDS (file.path('..','scatac_ArchR','selected_TF.rds'))
vf = VariableFeatures (FindVariableFeatures (metacells, nfeat=10000))
metacells_assay = metacells@assays$RNA@layers$data
rownames (metacells_assay) = rownames(srt)
metacells_assay = metacells_assay[unique(c(vf, selected_TF)),]

enricher_universe = vf
#do.fgsea = TRUE
gmt_annotations = c(
'h.all.v7.4.symbols.gmt'#,#,
#'c5.bp.v7.1.symbol.gmt'
)

sams = unique(metacells$sampleID)[!unique(metacells$sampleID) %in% c('HU37','HU62')]
gmt_annotation = gmt_annotations[1]
force = T
top_genes= 100
if (!file.exists(paste0('EnrichR_activeTF_cor_top_genes_ann_',gmt_annotation,'_top_genes_',top_genes,'.rds')) | force)
	{
	TF_cor_sample = list()
	for (sam in sams)
	  {
	  metacells_assay_sample = as.matrix(metacells_assay[,metacells$sampleID == sam])
	  TF_cor_sample[[sam]] = lapply (selected_TF, function(x) 
	  	{
	  	tc_cor = t(cor (metacells_assay_sample[x,], t(metacells_assay_sample)))
	  	tc_cor = setNames (tc_cor[,1], rownames (tc_cor))
	  	tc_cor[is.na(tc_cor)] = 0
	  	tc_cor = tc_cor[order(-tc_cor)]
	  	gmt.file = file.path ('..','..','git_repo','files',gmt_annotation)
	  	pathways = read.gmt (gmt.file)
	  	#pathways = split (pathways$gene, pathways$term)
	    message (paste ('EnrichR running module',x)) 
	    egmt = enricher(head (names(tc_cor),top_genes), TERM2GENE=pathways, universe = enricher_universe)
	    egmt@result
	    })
	   #EnrichRResAll[[ann]] = EnrichRResCluster    
	    # fgseaRes = fgseaMultilevel (pathways, 
		# 			tc_cor#, 
		# 			#minSize=15, 
		# 			#maxSize=1500,
		# 			#BPPARAM = NULL
		# 			)
		# fgseaResCol = collapsePathways (fgseaRes, stats = tc_cor, pathway = pathways)
		# fgseaRes[fgseaRes$pathway %in% fgseaResCol$mainPathways]
		}
	saveRDS (TF_cor_sample, paste0('EnrichR_activeTF_cor_top_genes_ann_',gmt_annotation,'_top_genes_',top_genes,'.rds'))
	} else {
	TF_cor_sample = readRDS (paste0('EnrichR_activeTF_cor_top_genes_ann_',gmt_annotation,'_top_genes_',top_genes,'.rds'))
	}

gmt.file = file.path ('..','..','git_repo','files',gmt_annotation)
pathways = read.gmt (gmt.file)
pathway_names = as.character (unique(pathways$term))

TF_cor_sample_mat = lapply(TF_cor_sample, function(x)
	{
	pathway_mat = do.call (cbind, lapply (x, function(y) 
		{
		y = y[pathway_names, 'p.adjust', drop=F]
		y$p.adjust = y$p.adjust < 0.05
		}))
	})

TF_cor_sum = Reduce ('+', TF_cor_sample_mat)
rownames (TF_cor_sum) = pathway_names
colnames (TF_cor_sum) = selected_TF
TF_cor_sum[is.na(TF_cor_sum)] = 0
#TF_cor_sum[TF_cor_sum < 3] = 0
TF_cor_sum = TF_cor_sum[rowSums (TF_cor_sum) > 0, ]
TF_cor_sum = TF_cor_sum[,colSums (TF_cor_sum) > 0 ]

# Export enrichment table ####
saveRDS (TF_cor_sum, 'enrichment_pathways_TFs.rds')

# Add correlation to bulk RNA sarcomatoid score ####
tf_sarc_cor = read.csv ('../../bulkRNA_meso/activeTF_sarcomatoid_correlation.csv', row.names = 1)
rownames (TF_cor_sum) = sub ('HALLMARK_','',rownames (TF_cor_sum))
hm = Heatmap (
		t(TF_cor_sum),
		column_names_rot =45, 
		col =palette_enrichment, 
		cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.0f", t(TF_cor_sum)[i, j]), x, y, gp = gpar(fontsize = 10, col='white'))
},
		border=T)
hm2 = Heatmap (
	tf_sarc_cor[colnames(TF_cor_sum),],
	column_split = colnames (tf_sarc_cor),#,column_names_rot = 45, 
	cluster_columns = F,
	col =palette_module_correlation_fun, 
	border=T)
pdf (file.path('Plots','pathway_enrich_TF_cor.pdf'), width=14, height=14)
draw (hm)# + hm2,  padding = unit(c(2, 60, 2, 2), "mm"))
dev.off()

### Export genes correlated to TFs ####
TF_cor_sample2 = list()
for (tf in selected_TF)
	{
	TF_cor_sample = list()
	for (sam in samplesID)
	  {
	  metacells_assay_sample = metacells_assay[,metacells$sampleID == sam]	 
	  tc_cor = t(cor (metacells_assay_sample[tf,], t(metacells_assay_sample)))
	  TF_cor_sample[[sam]] = tc_cor
	  }
	 TF_cor_sample2[[tf]] = do.call (cbind, TF_cor_sample)
	 TF_cor_sample2[[tf]][is.na(TF_cor_sample2[[tf]])] = 0
	 colnames (TF_cor_sample2[[tf]]) = samplesID
	}
saveRDS (TF_cor_sample2, 'TF_cor_genes_per_sample.rds')

# pdf (paste0 ('Plots/selected_TF_exp_corr_heatmaps.pdf'), width = 8,height=9)
# cor_TF_l
# dev.off()

# Correlate TFs with cNMFs on metacells ####
vf = VariableFeatures (FindVariableFeatures (metacells, nfeat=10000))
metacells_assay = metacells@assays$RNA@layers$data
rownames (metacells_assay) = rownames (srt)
metacells_assay = metacells_assay[unique(c(vf, selected_TF)),]

metacells = ModScoreCor (
        seurat_obj = metacells, 
        geneset_list = cnmf_spectra_unique, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'cNMF_', outdir = paste0(projdir,'Plots/'))

ccomp_df = metacells@meta.data[,c(names(cnmf_spectra_unique))]

tc_cor = lapply (unique(metacells$sampleID), function(x)
	{
	metacells_assay_sample = metacells_assay[selected_TF,metacells$sampleID == x]
	ccomp_df_sample = ccomp_df[metacells$sampleID == x,]
	res = sapply (selected_TF, function (y) cor (metacells_assay_sample[y,], ccomp_df_sample, method = 'spearman'))
	rownames (res) = colnames (ccomp_df_sample)
	res
	})
names (tc_cor) = unique(metacells$sampleID)
#lapply (tc_cor, function(x) {x = x['cNMF19',]; head(x[order(-x)],10)})


cor_TF = tc_cor[[1]]
cor_TF[] = tapply(unlist(tc_cor), rep(seq(length(tc_cor[[1]])),length(tc_cor)), FUN=function(x) median(x,na.rm=T))

dim (cor_TF)

pdf ('Plots/cor_nmf_TF.pdf',width=15, height=3)
Heatmap (cor_TF, 
	column_names_gp = gpar(fontsize = 5),
	row_names_gp = gpar(fontsize = 5),
	clustering_distance_rows = 'pearson',
	clustering_distance_columns='pearson')
dev.off()

tc_cor_long = lapply (seq_along(tc_cor), function(x){tmp = gather (as.data.frame(tc_cor[[x]]['cNMF20',,drop=F]),TF, expression); tmp$sample = names(tc_cor)[x];tmp})
tc_cor_long = do.call (rbind, tc_cor_long)
tf_median = sapply (unique (tc_cor_long$TF), function(x) {median(tc_cor_long[tc_cor_long$TF == x,'expression'],na.rm=T)})
tf_median = tf_median[order (-tf_median)]
tc_cor_long$TF = factor (tc_cor_long$TF, levels = names (tf_median))
#filter_TF = names(tf_median[tf_median > quantile (tf_median)[4] | tf_median < quantile (tf_median)[2]])
tc_cor_long_box = ggplot (tc_cor_long, aes (x = TF, y = expression)) + 
geom_boxplot(outlier.size = 0.1) + 
theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) + 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf (paste0('Plots/distribution_TF_cor_cNMF.pdf'), width=10)
tc_cor_long_box
dev.off()

filter_TF = c('SOX9','HIC1','MESP1','SNAI2','TWIST1','OLIG3','PITX2','SMARCC2','MEF2A')
tc_cor_long = tc_cor_long[tc_cor_long$TF %in% filter_TF,]

tc_cor_long_box = ggplot (tc_cor_long, aes (x = TF, y = expression)) + 
geom_boxplot(outlier.size = 0.1) + 
geom_point() + 
theme_classic() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=8))
pdf (paste0('Plots/distribution_TF_cor_cNMF_filtered.pdf'), width=4, height=3)
tc_cor_long_box
dev.off()

# Order cells per samples along SOX9 expression and plot the rest of TF expression together
mMat = metacells@assays$RNA@layers$data[rownames(metacells) %in% selected_TF,]
rownames(mMat) = rownames(metacells)[rownames(metacells) %in% selected_TF]

traj_sample = list()
TF = 'SOX9'
for (sam in unique(metacells$sampleID))
    {
    mMat_ordered_sample = mMat[,metacells$sampleID == sam]
    mMat_ordered_sample = mMat_ordered_sample[, order(mMat_ordered_sample['SOX9',])]
    mMat_ordered_sample_sclaled = as.matrix(mMat_ordered_sample)
    mMat_ordered_sample_sclaled = t(scale(t(as.matrix(mMat_ordered_sample_sclaled))))
    mMat_ordered_sample_sclaled[is.na(mMat_ordered_sample_sclaled)] = 0
    ha = HeatmapAnnotation (tf = mMat_ordered_sample_sclaled[TF,])
    traj_sample[[sam]] = Heatmap (
    	top_annotation = ha,
      mMat_ordered_sample_sclaled[!rownames(mMat_ordered_sample_sclaled) %in% TF,], 
      col = viridis::plasma(100), 
      cluster_columns=F, name = sam,
      row_names_gp = gpar(fontsize = 4))
    }

pdf ('Plots/sarc_trajectory_per_sample.pdf', height=10)
traj_sample
dev.off()


# Run cNMF x TF correlation across all metacells ####
ccomp_df = metacells@meta.data[,c(names(cnmf_spectra_unique))]
metacells_assay = metacells@assays$RNA@layers$data
rownames (metacells_assay) = rownames(srt)
metacells_assay = metacells_assay[unique(c(vf, selected_TF)),]

res = sapply (selected_TF, function (y) cor (metacells_assay[y,], ccomp_df, method = 'spearman'))
rownames (res) = colnames (ccomp_df)

pdf ('Plots/cor_nmf_TF_across_all.pdf',width=15, height=3)
Heatmap (res, 
	column_names_gp = gpar(fontsize = 5),
	row_names_gp = gpar(fontsize = 5),
	clustering_distance_rows = 'pearson',
	clustering_distance_columns='pearson')
dev.off()


# Run cNMF x TF correlation across all metacells subsampled ####
ccomp_df = metacells@meta.data[,c(names(cnmf_spectra_unique))]
metacells_assay = metacells@assays$RNA@layers$data
rownames (metacells_assay) = rownames(srt)
metacells_assay = metacells_assay[unique(c(vf, selected_TF)),]

subsample = function(metagroupname, size) {
	tmp = split (metagroupname, metagroupname)
	tmp = lapply (tmp, function (s) 
		{
		if (!size > length (s))	
			{
			names(s) %in% names(sample (s, size, replace = F))
			} else {
			rep (TRUE, length (s))	
		}})
	do.call (c, tmp)
	}

subsampled = subsample (metacells$sampleID, 5)
table (metacells$sampleID, subsampled)

res = sapply (selected_TF, function (y) cor (metacells_assay[y,subsampled], ccomp_df[subsampled,], method = 'spearman'))
rownames (res) = colnames (ccomp_df)

# Export cNMF TF correlation matrix ####
saveRDS (res, 'cNMF_TF_correlation_subsampled.rds')
pdf ('Plots/cor_nmf_TF_across_all_subsampled.pdf',width=15, height=3)
Heatmap (res, 
	column_names_gp = gpar(fontsize = 5),
	row_names_gp = gpar(fontsize = 5),
	clustering_distance_rows = 'pearson',
	clustering_distance_columns='pearson')
dev.off()


# Order metacells by sarcomatoid score ####
mMat = metacells@assays$RNA@layers$data[rownames(metacells) %in% selected_TF,]
rownames(mMat) = rownames(metacells)[rownames(metacells) %in% selected_TF]
mMat = mMat[,order (metacells$cNMF20)]
mMat = mMat[colnames(res)[order(res['cNMF20',])],]


top_10_sarc_TF = head (colnames(res)[order(res['cNMF20',])],20)
top_10_sarc_TF = c(top_10_sarc_TF,head (colnames(res)[order(-res['cNMF20',])], 20))
mMat_scaled = t(scale(t(as.matrix(mMat))))[top_10_sarc_TF,]
mMat_scaled[mMat_scaled > abs(min(mMat_scaled))] = abs (min (mMat_scaled))
ha = HeatmapAnnotation (sample = metacells$sampleID[order (metacells$cNMF20)], col = list(sample = palette_sample))
sarc_traj = Heatmap (mMat_scaled, 
	top = ha,
  col = palette_expression, 
  cluster_columns=F, 
  cluster_rows=F,
  row_names_gp = gpar(fontsize = 14))

png ('Plots/sarc_trajectory_across_sample.png', height=600)
sarc_traj
dev.off()



# Make PCA of cNMF TF correlations
cor_TF_pca = prcomp (cor_TF)
# Example using the mtcars dataset


# Standardize the data: mean = 0, SD = 1

# Run PCA
pca_result <- prcomp(t(cor_TF), center = TRUE, scale. = TRUE)

# Summary of the PCA
summary(pca_result)

# Plot PCA
pdf ('Plots/cNMF_TF_correlation_pca.pdf')
fviz_pca_var(pca_result,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
fviz_pca_ind(pca_result,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
dev.off()

srt_TF = CreateSeuratObject (counts = cor_TF, data = cor_TF)
srt_TF[['RNA']]$scale.data = cor_TF
VariableFeatures(srt_TF) = rownames (cor_TF)
srt_TF = RunPCA (srt_TF)
srt_TF = RunUMAP (srt_TF, n.components = 3, n.neighbors = 20)

cor_TF_umap <- umap(cor_TF, config=umap.defaults)







# Correlate TFs against each other on metacells ####
selected_TF = readRDS (file.path('..','scatac_ArchR','selected_TF.rds'))
TF_dev_order = 
sams = c('P1','P11','P12','P13','P4','P5','P8')

metacells_assay = metacells@assays$RNA@layers$data
rownames (metacells_assay) = rownames (srt)
metacells_assay = metacells_assay[selected_TF,]

tc_cor = lapply (sams, function(x)
	{
	cor(as.matrix(t(metacells_assay[,metacells$sampleID == x])), method='spearman')
	})
names (tc_cor) = sams
#lapply (tc_cor, function(x) {x = x['cNMF19',]; head(x[order(-x)],10)})


cor_TF = tc_cor[[1]]
cor_TF[] = tapply(unlist(tc_cor), rep(seq(length(tc_cor[[1]])),length(tc_cor)), FUN=function(x) mean(x,na.rm=T))

saveRDS (cor_TF, 'selected_TF_metacells_scrna.rds')

dim (cor_TF)
diag(cor_TF) = 0
pdf ('Plots/cor_selected_TF_metacells_heatmap.pdf',width=5, height=5)
Heatmap (cor_TF, 
	column_names_gp = gpar(fontsize = 5),
	row_names_gp = gpar(fontsize = 5),
	cluster_rows=F,
	cluster_columns=F,
	clustering_distance_rows = 'pearson',
	clustering_distance_columns='pearson',
	col = palette_expression_cor_fun(cor_TF))
dev.off()

