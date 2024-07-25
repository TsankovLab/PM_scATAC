
####### ANALYSIS of P11 TUMOR #######
set.seed(1234)
source (file.path ('..','..','PM_scATAC','load_packages.R')

projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_scrna_P11'
dir.create (file.path(projdir,'Plots'), recursive =T)
setwd (projdir)
projdir_scatac = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR'

#devtools::install_github("immunogenomics/presto") #needed for DAA
source ('../../PM_scATAC/useful_functions.R')
source ('../../PM_scATAC/ggplot_aestetics.R')
source ('../../PM_scATAC/scATAC_functions.R')
source ('../../PM_scATAC/palettes.R')

set.seed(1234)
addArchRThreads (threads = 8) 
addArchRGenome ("Hg38")

# Load scATAC ####
if (!file.exists ('Save-ArchR-Project.rds')) 
  { source (file.path('..','..','PM_scATAC','scatac_tumor_create_ArchRobj.R'))
  } else {
 archp = loadArchRProject (projdir_scatac)   
  }

# Load scRNA ####
srt = readRDS (file.path('..','scrna','srt.rds'))
#sarc_order = read.csv ('../scrna/cnmf20_sarcomatoid_sample_order.csv', row.names=1)

### Test correlation of chr18 mega hubs regions from P11 with other genes ####

# Load metacells object ####
if (!file.exists (file.path('..','scrna','metacells.rds'))
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

# Load P11 megahubs regions ####
region = readRDS ('../scatac_ArchR/P11_chr18_region.rds')

# Make gene modules overlapping megahubs regions in P11 ####
all_genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes_in_region = all_genes$gene_id[subjectHits(findOverlaps (region, all_genes))]
genes_in_region = list(chr18_q23 = as.data.frame(org.Hs.egSYMBOL)[match (genes_in_region, as.data.frame(org.Hs.egSYMBOL)[,1]),'symbol'])

metacells = ModScoreCor (
    seurat_obj = metacells, 
    geneset_list = genes_in_region, 
    cor_threshold = NULL, 
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'mut', outdir = paste0(projdir,'Plots/'))

# check how it looks on UMAP ####
srt = ModScoreCor (
    seurat_obj = srt, 
    geneset_list = genes_in_region, 
    cor_threshold = NULL, 
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'mut', outdir = paste0(projdir,'Plots/'))

# Show mega hubs region and other P11 specific genes ####
reductionName = 'umap'
pdf ('Plots/chr18_q23_umap.pdf', width =10, height = 5)
wrap_plots (
	DimPlot (srt, group.by = 'sampleID'), 
	fp (srt, 'chr18_q23', reduction = reductionName)[[1]],
	fp (srt, 'MEF2A', reduction = reductionName)[[1]],
	fp (srt, 'HOXC6', reduction = reductionName)[[1]],
	fp (srt, 'HOXB13', reduction = reductionName)[[1]])
dev.off()


### Correlate megahubs region with all other genes in P11 metacells ####
ccomp_df = metacells@meta.data[,'chr18_q23', drop=F]
metacells_assay = metacells@assays$RNA@layers$data
rownames (metacells_assay) = rownames(srt)

res = sapply (unique(metacells$sampleID), function (y) cor (t(metacells_assay[,metacells$sampleID == y]), ccomp_df[metacells$sampleID == y,1, drop=F], method = 'pearson'))
rownames (res) = rownames (metacells_assay)
colnames (res) = unique (metacells$sampleID)

res = as.data.frame (res)
res = res[order (-res$P11),]
res = res[!rownames(res) %in% genes_in_region[[1]], ]
head (res)

y = 'P11'
pdf (file.path ('Plots','chr18_q23_CNDP2_scatter.pdf'))
plot (x = metacells_assay['CNDP2',metacells$sampleID == y], y = ccomp_df[metacells$sampleID == y,1])
dev.off()

### Restrict above analysis only on metacells high for megahubs in P11 ####
metacells_assay_P11s = metacells_assay[,metacells$sampleID == y & ccomp_df[,1] > 0.4]
ccomp_df_P11s = ccomp_df[metacells$sampleID == y & ccomp_df[,1] > 0.4,]

res_p11s = cor (t(metacells_assay_P11s), ccomp_df_P11s, method= 'pearson')
res_p11s = res_p11s[order (-res_p11s[,1]),]
res_p11s = res_p11s[!names(res_p11s) %in% genes_in_region[[1]]]

write.csv (res_p11s, 'correlated_genes_p11s_chr18_q23.csv')

pdf (file.path ('Plots','chr18_q23_ZKSCAN5_scatter.pdf'))
plot (x = metacells_assay_P11s['ZKSCAN5',], y = ccomp_df_P11s)
dev.off()

# Run DEG between high megahubs cells and low within P11 #### 
srt_p11 = srt[,srt$sampleID == 'P11']
pdf (file.path ('Plots','chr18_q23_in_P11_cells_hist.pdf'))
hist (srt_p11$chr18_q23)
dev.off()

srt_p11$megahubs = ifelse (srt_p11$chr18_q23 > 0.3, 'high','low')

deg_res = FindMarkers (srt_p11, ident.1 = 'high', ident.2 = 'low', group.by='megahubs')
head (deg_res,10)


logfcThreshold = 2
pvalAdjTrheshold = 0.05
deg_res$sig = ifelse (abs(deg_res$avg_log2FC) > logfcThreshold & deg_res$p_val_adj < pvalAdjTrheshold, 1,0)
deg_res$sig = deg_res$sig * sign (deg_res$avg_log2FC)
deg_res$alpha = ifelse (deg_res$sig != '0', .5,.2)
deg_res$sig[deg_res$sig == 1] = 'High'
deg_res$sig[deg_res$sig == -1] = 'Low'

# res_filtered = deg_res[abs(deg_res$avg_log2FC) > logfcThreshold & deg_res$p_val_adj < pvalAdjTrheshold,]
# res_filtered = head (rownames(res_filtered)[order (-abs(res_filtered$avg_log2FC))],20)
deg_res$labels = ''
deg_res[genes_in_region[[1]],'labels'] = genes_in_region[[1]]
deg_res = na.omit (deg_res)
vp = ggplot (deg_res, aes(x=avg_log2FC, y= -log10(p_val_adj))) +
    geom_point(shape=19, aes (color = sig, alpha=alpha),size=1) +
    geom_vline(xintercept = logfcThreshold, linetype="dotted", 
                color = "grey20", size=1) +
    geom_vline(xintercept = -logfcThreshold, linetype="dotted", 
                color = "grey20", size=1) +
    geom_hline(yintercept = -log10(pvalAdjTrheshold), linetype="dotted", 
                color = "grey20", size=1) + 
    geom_text_repel (size=2, data = deg_res, aes(label = labels),segment.size=.2) + 
    ggtitle ('Hubs differential accessibility') +
    #geom_label_repel (size=2,max.overlaps=10000, data = deg2_cl, aes(label = show_genes), color='red') + 
    scale_color_manual (values = c("0"='grey20',"Low"='green',"High"='red')) + 
    scale_fill_manual (values = c("0"='grey20',"Low"='green',"High"='red')) + 
    theme_light()

pdf (file.path ('Plots', 'P11_DEG_megahubs_high_low_volcano.pdf'),height=3,width=5)
vp
dev.off()

### Check raw counts of megahubs genes in high megahubs cells ####
p11_mh_high_mat = as.data.frame (AverageExpression (srt_p11, group.by = 'megahubs')[[1]])
p11_mh_high_mat = p11_mh_high_mat[,'high', drop=F]
p11_mh_high_mat = p11_mh_high_mat[order(-p11_mh_high_mat$high),, drop=F]

which (rownames(p11_mh_high_mat) %in% genes_in_region[[1]])
megahubs_gene_exp = data.frame (row.names = rownames(p11_mh_high_mat), expression = p11_mh_high_mat$high, megahubs = rownames(p11_mh_high_mat) %in% genes_in_region[[1]])
megahubs_gene_exp$rank = 1:nrow(megahubs_gene_exp)
megahubs_gene_exp$rank = factor (megahubs_gene_exp$rank, levels = megahubs_gene_exp$rank)
megahubs_gene_exp$label = ''
megahubs_gene_exp$label[megahubs_gene_exp$megahubs == TRUE] = rownames(megahubs_gene_exp)[megahubs_gene_exp$megahubs == TRUE]
megahubs_gene_exp$alpha = ifelse (megahubs_gene_exp$megahubs, 0.9,0.2)
megahubs_gene_exp = megahubs_gene_exp[megahubs_gene_exp$expression != 0, ]
megahubs_gene_exp = head (megahubs_gene_exp,100)
vp = ggplot (megahubs_gene_exp, aes(x=rank, y = expression, fill= megahubs)) +
    geom_bar(stat = 'identity', aes (fill = megahubs) ) +
    geom_text_repel (
    	size=2, segment.curvature = 0.1,
    	data = megahubs_gene_exp, 
    	aes(label = label),
    	segment.size=.2,
    	nudge_y = 50) + 
    theme_classic()
    #ggtitle ('Hubs differential accessibility') +
    #geom_label_repel (size=2,max.overlaps=10000, data = deg2_cl, aes(label = show_genes), color='red') + 
    #scale_color_manual (values = c("0"='grey20',"Low"='green',"High"='red')) + 
    #scale_fill_manual (values = c("0"='grey20',"Low"='green',"High"='red')) + 

pdf (file.path ('Plots', 'P11_high_cells_megahubs_gene_rank.pdf'),height=3,width=4)
vp
dev.off()





### Correlate HOXB13 with all other genes in P11 metacells ####
gene = 'HOXB13'
metacells_assay = metacells@assays$RNA@layers$data
rownames (metacells_assay) = rownames(srt)
metacells@meta.data[,gene] = metacells_assay[gene, ]

### Restrict above analysis only on metacells high for megahubs in P11 ####
y = 'P11'
summary (metacells@meta.data[,gene])
metacells_assay_P11s = metacells_assay[,metacells$sampleID == y & metacells@meta.data[,gene] > 0.03]
res_p11s = cor (t(metacells_assay_P11s), t(metacells_assay_P11s)[,gene], method= 'pearson')
res_p11s = res_p11s[order (-res_p11s[,1]),]
#res_p11s = res_p11s[!names(res_p11s) %in% genes_in_region[[1]]]

write.csv (res_p11s, 'correlated_genes_p11s_HOXB13.csv')

pdf (file.path ('Plots','chr18_q23_ZKSCAN5_scatter.pdf'))
plot (x = metacells_assay_P11s['ZKSCAN5',], y = ccomp_df_P11s)
dev.off()

# Check Proliferation index of HOX + cluster vs others ####
pdf (file.path ('Plots','cellcycle_fplot.pdf'))
FeaturePlot (srt, feature = 'cc')
dev.off()


### Plot heatmap of all HOX genes across tumors ####
gene = rownames(srt)[grep ('^HOX', rownames(srt))]

ha = HeatmapAnnotation (sample = srt$sampleID)
hm = Heatmap (as.matrix(srt@assays$RNA@data[gene,]), top_annotation=ha, 
	column_names_gp = gpar(fontsize = 0) )

pdf (file.path ('Plots','HOX_genes_cells_heatmap.pdf'), width=5, height=4)
hm
dev.off()



### Investigate P11 heterogeneity ####
# Identify differential TFs between HOX+ and HOX- clusters ####
library (presto)

#sample_names_rna = c('P1','P14','P13','P3','P12','P5','P11','P4','P8','P14','HU37','HU62')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat_P11 = mMat[,archp$Sample2 == 'P11']
P11_clusters = archp$Clusters[archp$Sample2 == 'P11']
P11_clusters = ifelse (P11_clusters == 'C14','P11_small','P11_large')
p11_dev_rna = wilcoxauc (mMat_P11 , y = P11_clusters)
p11_dev_rna = p11_dev_rna[p11_dev_rna$group == 'P11_large',]
head (dev[order(p11_dev_rna$logFC),],20)
p11_dev_rna[p11_dev_rna$feature == 'NFATC1',]
# Make volcano plot showing also scRNA RNA expression
pdf (paste0('Plots/scrna_clusters_umap.pdf'), 15,15)
wrap_plots (DimPlot (srt, group.by = 'seurat_clusters'), DimPlot (srt, group.by = 'sampleID'))
dev.off()

ps = log2(as.data.frame (AverageExpression (srt, features = p11_dev_rna$feature, group.by = 'seurat_clusters')[[1]]) +1)
#ps = ps[, colnames(ps) %in% sample_names_rna]
ps_P11_diff = ps[, 'g14', drop=F] - ps[, 'g9', drop=F]

p11_dev_rna$rna_diff = NA
p11_dev_rna$rna_diff = ps_P11_diff[p11_dev_rna$feature,1]
p11_dev_rna$rna_diff[is.na(p11_dev_rna$rna_diff)] = 0

logfcThreshold = .1
pvalAdjTrheshold = 0.05
p11_dev_rna$sig = ifelse (abs(p11_dev_rna$logFC) > logfcThreshold & p11_dev_rna$padj < pvalAdjTrheshold, 1,0)
p11_dev_rna$sig = p11_dev_rna$sig * sign (p11_dev_rna$logFC)
p11_dev_rna$sig[p11_dev_rna$sig == -1] = 'HOX+'
p11_dev_rna$sig[p11_dev_rna$sig == 1] = 'HOX-'

p11_dev_rna$rna_sign = ifelse (abs(p11_dev_rna$logFC) > logfcThreshold & p11_dev_rna$padj < pvalAdjTrheshold, 1,0)
p11_dev_rna$rna_sign = p11_dev_rna$rna_sign * sign (p11_dev_rna$rna_diff)
p11_dev_rna$rna_sign[p11_dev_rna$rna_sign == -1] = 'HOX+'
p11_dev_rna$rna_sign[p11_dev_rna$rna_sign == 1] = 'HOX-'

res_filtered = p11_dev_rna[abs(p11_dev_rna$logFC) > logfcThreshold & p11_dev_rna$padj < pvalAdjTrheshold,]
res_filtered = head (res_filtered$feature[order (-abs(res_filtered$logFC))],20)
p11_dev_rna$labels = ''
p11_dev_rna$labels[match (res_filtered, p11_dev_rna$feature)] = res_filtered
vp = ggplot (p11_dev_rna, aes(x=logFC, y= -log10(padj))) +
    geom_point(shape=21, aes (fill = sig, color = rna_sign, size = abs(rna_diff)), alpha=.5) +
    geom_vline(xintercept = logfcThreshold, linetype="dotted", 
                color = "grey20", size=1) +
    geom_vline(xintercept = -logfcThreshold, linetype="dotted", 
                color = "grey20", size=1) +
    geom_hline(yintercept = -log10(pvalAdjTrheshold), linetype="dotted", 
                color = "grey20", size=1) + 
    geom_text_repel (size=2, data = p11_dev_rna, aes(label = labels),segment.size=.2) + 
    ggtitle ('Hubs differential accessibility') +
    #geom_label_repel (size=2,max.overlaps=10000, data = deg2_cl, aes(label = show_genes), color='red') + 
    scale_color_manual (values = c("0"='grey77',"HOX-"='green',"HOX+"='red')) + 
    scale_fill_manual (values = c("0"='grey77',"HOX-"='green',"HOX+"='red')) + 
    theme_light()

pdf (file.path ('Plots', 'P11_TF_volcano.pdf'),height=3,width=5)
vp
dev.off()

# Check most expressed HOX genes from the TF co-deviation module analysis above ####
HOX_module = 3
hox_dev_rna = p11_dev_rna[p11_dev_rna$feature %in% names(km$cluster)[km$cluster == HOX_module],]
hox_dev_rna = hox_dev_rna[order (hox_dev_rna$rna_diff),]
head (hox_dev_rna[grep ('HOX', hox_dev_rna$feature),])

