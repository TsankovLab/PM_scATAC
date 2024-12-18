conda activate meso_scatac
R

####### ANALYSIS of P11 TUMOR #######
set.seed(1234)

projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scrna_P11'
projdir_scatac = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR'

dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

# Set # of threads and genome reference ####
addArchRThreads(threads = 1) 
addArchRGenome("hg38")

# Load scATAC ####
 archp = loadArchRProject (projdir_scatac)   
  

# Load scRNA ####
if (!file.exists ('srt.rds'))
{
srt = readRDS (file.path ('..','..','main','scrna','srt.rds'))
srt = srt[,srt$sampleID == 'P11']
srt = NormalizeData (object = srt, normalization.method = "LogNormalize", scale.factor = 10000)

nfeat = 3000
srt = FindVariableFeatures (srt, selection.method = "vst", nfeat = nfeat)
srt = ScaleData (srt)	
srt = RunPCA (srt, features = VariableFeatures (object = srt), npcs = ifelse(ncol(srt) <= 30,ncol(srt)-1,30), ndims.print = 1:5, nfeat.print = 5, verbose = FALSE)

reductionName = 'umap'
reductionSave = 'pca'
reductionGraphKnn = 'RNA_knn'
reductionGraphSnn = 'RNA_snn' 

srt = RunUMAP (object = srt, reduction = reductionSave, dims = 1:30)
srt = FindNeighbors (srt)
srt = FindClusters (srt, res=.6)
pdf (file.path ('Plots','celltypes_umap.pdf'))
DimPlot (srt, group.by = 'celltype_simplified', reduction = 'umap')
DimPlot (srt, group.by = 'seurat_clusters', reduction = 'umap', label=T)
fp (srt, gene = 'HOXB13',reduction= reductionName)
dev.off()
saveRDS (srt, 'srt.rds')
} else {
srt = readRDS ('srt.rds')    
}



#### Run cNMF ####
nfeat = 5000
force=F
k_list = c(5:30)
k_selections = c(5:30)
k_selection = 25
cores= 100

cnmf_name = 'scrna_p11'
repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'

### RUN consensus NMF ####
source (file.path ('..','..','git_repo','utils','cnmf_prepare_inputs.R')) 

### Import and format spectra files ####
k_selection = 25
cnmf_name = 'scrna_p11'
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
ccomp_df = srt@meta.data[,c(names(cnmf_spectra_unique),'sampleID','celltype_simplified'), drop=FALSE]
      #ccomp_df = aggregate (ccomp_df, by=as.list(srt_wgcna@meta.data[,metaGroupNames,drop=F]), mean)    
bp1 = lapply (names(cnmf_spectra_unique), function(x) {
            ggplot (ccomp_df, aes_string (x= 'celltype_simplified', y= x)) +
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

  
png (file.path('Plots',paste0('cNMF_module_scores_boxplots_',k_selection,'nfeat_',nfeat,'.png')),5000,5000,res=300)
print (wrap_plots (bp1))
dev.off()

metaGroupNames = c('celltype_simplified')
  umap_df = data.frame (srt[[reductionName]]@cell.embeddings, srt@meta.data[,c(names(cnmf_spectra_unique),metaGroupNames)])
  umap_p1 = lapply (names(cnmf_spectra_unique), function(x) ggplot(data = umap_df) + 
  geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], color = x), size = .1) + 
  scale_colour_gradientn (colours = rev(brewer.pal (n = 11, name = "RdBu")),limits=c(-max (abs(umap_df[,x])), max (abs(umap_df[,x])))) +
  ggtitle (x) + 
  #facet_wrap (as.formula(paste("~", metaGroupNames[3]))) + 
  theme_classic() +
  theme_void())
  
pdf (file.path('Plots',paste0('cNMF_module_scores_umap_',k_selection,'.pdf')),24,25)
print (wrap_plots (umap_p1))
dev.off()

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
    gmt.file = file.path ('..','..','git_repo','files',gmt_annotation)
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


## Rerun cNMF but only in HOXB13+ cluster ####

#### Run cNMF ####
nfeat = 5000
force=F
k_list = c(5:30)
k_selections = c(5:30)
k_selection = 10
cores= 100

cnmf_name = 'scrna_p11_HOX'
repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'

### RUN consensus NMF ####
srt2 = srt
srt = srt[,srt$seurat_clusters == 1]
source (file.path ('..','..','git_repo','utils','cnmf_prepare_inputs.R')) 

### Import and format spectra files ####
k_selection = 10
cnmf_name = 'scrna_p11_HOX'
source (file.path ('..','..','git_repo','utils','cnmf_format_spectra_files.R')) 

top_nmf_genes = 200
cnmf_spectra_unique = lapply (cnmf_spectra_unique, function(x) head (x, top_nmf_genes))
write.csv (cnmf_spectra_unique, file.path (cnmf_out,'cnmf_list.csv'))

# Add module score of cNMF modules ####
force = T
if (!all (names(cnmf_spectra_unique) %in% colnames (srt@meta.data)) | force)
    {
    srt = ModScoreCor (
            seurat_obj = srt, 
            geneset_list = cnmf_spectra_unique, 
            cor_threshold = NULL, 
            pos_threshold = NULL, # threshold for fetal_pval2
            listName = 'cNMF_', outdir = paste0(projdir,'Plots/'))

    }

srt$cNMF = 'cNMF'
nfeat = 3000
srt = FindVariableFeatures (srt, selection.method = "vst", nfeat = nfeat)
srt = ScaleData (srt)   
srt = RunPCA (srt, features = VariableFeatures (object = srt), npcs = ifelse(ncol(srt) <= 30,ncol(srt)-1,30), ndims.print = 1:5, nfeat.print = 5, verbose = FALSE)

reductionName = 'umap'
reductionSave = 'pca'
reductionGraphKnn = 'RNA_knn'
reductionGraphSnn = 'RNA_snn' 

srt = RunUMAP (object = srt, reduction = reductionSave, dims = 1:30)
srt = FindNeighbors (srt)
srt = FindClusters (srt, res=.6)

ccomp_df = srt@meta.data[,c(names(cnmf_spectra_unique),'sampleID','celltype_simplified'), drop=FALSE]
      #ccomp_df = aggregate (ccomp_df, by=as.list(srt_wgcna@meta.data[,metaGroupNames,drop=F]), mean)    
bp1 = lapply (names(cnmf_spectra_unique), function(x) {
            ggplot (ccomp_df, aes_string (x= 'celltype_simplified', y= x)) +
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

metaGroupNames = c('celltype_simplified')
  umap_df = data.frame (srt[[reductionName]]@cell.embeddings, srt@meta.data[,c(names(cnmf_spectra_unique),metaGroupNames)])
  umap_p1 = lapply (names(cnmf_spectra_unique), function(x) ggplot(data = umap_df) + 
  geom_point (mapping = aes_string (x = colnames(umap_df)[1], y= colnames(umap_df)[2], color = x), size = .1) + 
  scale_colour_gradientn (colours = rev(brewer.pal (n = 11, name = "RdBu")),limits=c(-max (abs(umap_df[,x])), max (abs(umap_df[,x])))) +
  ggtitle (x) + 
  #facet_wrap (as.formula(paste("~", metaGroupNames[3]))) + 
  theme_classic() +
  theme_void())
  
pdf (file.path(cnmf_out,'Plots',paste0('cNMF_module_scores_umap_',k_selection,'_HOX.pdf')),width=10,7)
print (wrap_plots (umap_p1))
dev.off()

enricher_universe = rownames(srt)
#do.fgsea = TRUE
gmt_annotations = c(
'h.all.v7.4.symbols.gmt',#,
'c5.bp.v7.1.symbol.gmt',
'c3.tft.v7.1.symbol.gmt'
)
gmt_annotation = gmt_annotations[2]
if (!file.exists (paste0('cNMF_normalized/',cnmf_out, '/EnrichR_cNMF_module_genes_k_',k_selection,'_top_nmf_genes_',top_nmf_genes,'_ann_',gmt_annotation,'.rds')) | force)
  {
    gmt.file = file.path ('..','..','git_repo','files',gmt_annotation)
    pathways = read.gmt (gmt.file)
    EnrichRResCluster = list()
    for (i in names(cnmf_spectra_unique))
      {
      message (paste ('EnrichR running module',i)) 
      sig_genes = cnmf_spectra_unique[[i]]
      if (!any (sig_genes %in% pathways$gene)) next
      egmt = enricher (sig_genes, TERM2GENE=pathways, universe = enricher_universe)
      egmt@result$ID = stringr::str_trunc (egmt@result$ID, 100)
      EnrichRResCluster[[i]] = egmt@result
      }
    EnrichRResAll = EnrichRResCluster
    saveRDS (EnrichRResAll, paste0(cnmf_out, '/EnrichR_cNMF_module_genes_k_',k_selection,'_top_nmf_genes_',top_nmf_genes,'_ann_',gmt_annotation,'.rds'))
    } else {
  EnrichRResAll = readRDS (paste0(cnmf_out, '/EnrichR_cNMF_module_genes_k_',k_selection,'_top_nmf_genes_',top_nmf_genes,'_ann_',gmt_annotation,'.rds'))
  }
  pvalAdjTrheshold = 0.05
  top_pathways = 10
  EnrichRes_dp = dotGSEA (enrichmentsTest_list = EnrichRResAll, type = 'enrich', padj_threshold = pvalAdjTrheshold, top_pathways= top_pathways)
  pdf (paste0(cnmf_out, '/Plots/EnrichR_',i,'_k_',k_selection,'_top_nmf_genes_',top_nmf_genes,'_ann_',gmt_annotation,'_dotplots.pdf'), width = 12 + length(length (cnmf_spectra_unique))/10, height = 15)
  print (EnrichRes_dp)
  dev.off()

# Restore original seurat object ####
srt = srt2


### Test correlation of chr18 mega hubs regions from P11 with other genes ####

# Load metacells object ####
metacells = readRDS (file.path('..','scrna','metacells.rds'))

# Load P11 megahubs regions ####
region = readRDS (file.path ('..','scatac_ArchR','P11_chr18_region.rds'))

# Make gene modules overlapping megahubs regions in P11 ####
all_genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes_in_region = all_genes$gene_id[subjectHits(findOverlaps (region, all_genes))]
genes_in_region = list(chr18_q23 = as.data.frame(org.Hs.egSYMBOL)[match (genes_in_region, as.data.frame(org.Hs.egSYMBOL)[,1]),'symbol'])

metacells = ModScoreCor (
    seurat_obj = metacells, 
    geneset_list = genes_in_region, 
    cor_threshold = NULL, 
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'mut', outdir = NULL)

# check how it looks on UMAP ####
srt = ModScoreCor (
    seurat_obj = srt, 
    geneset_list = genes_in_region, 
    cor_threshold = NULL, 
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'mut', outdir = NULL)

# Show mega hubs region and other P11 specific genes ####
reductionName = 'umap'
pdf (file.path('Plots','chr18_q23_umap.pdf'), width =10, height = 5)
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
rownames (metacells_assay) = rownames(metacells)

res = sapply (unique(metacells$sampleID3), function (y) cor (t(metacells_assay[,metacells$sampleID3 == y]), ccomp_df[metacells$sampleID3 == y,1, drop=F], method = 'pearson'))
rownames (res) = rownames (metacells_assay)
colnames (res) = unique (metacells$sampleID3)

res = as.data.frame (res)
res = res[order (-res$P11),]
res = res[!rownames(res) %in% genes_in_region[[1]], ]
head (res)

y = 'P11'
pdf (file.path ('Plots','chr18_q23_CNDP2_scatter.pdf'))
plot (x = metacells_assay['CNDP2',metacells$sampleID3 == y], y = ccomp_df[metacells$sampleID3 == y,1])
dev.off()

### Restrict above analysis only on metacells high for megahubs in P11 ####
metacells_assay_P11s = metacells_assay[,metacells$sampleID3 == y & ccomp_df[,1] > 0.4]
ccomp_df_P11s = ccomp_df[metacells$sampleID3 == y & ccomp_df[,1] > 0.4,]

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
rownames (metacells_assay) = rownames(metacells)
metacells@meta.data[,gene] = metacells_assay[gene, ]

### Restrict above analysis only on metacells high for HOXB13 in P11 ####
y = 'P11_HOX'
summary (metacells@meta.data[,gene])
metacells_assay_P11s = metacells_assay[,metacells$sampleID3 == y & metacells@meta.data[,gene] > 0.5]


pdf (file.path ('Plots','HOXB13_metacells_scatter.pdf'))
plot (x = metacells_assay_P11s['HOXB13',], y = metacells_assay_P11s['NFATC1',])
dev.off()

# Double check the metacells selection chosed ####
HOXB13_cells_selection = unlist(lapply (
	metacells$cells_merged[metacells$sampleID3 == y & 
	metacells@meta.data[,gene] > 0.5], function(x) unlist(strsplit(x, ','))))
srt$HOXB13_cells_selection = colnames(srt) %in% HOXB13_cells_selection

pdf (file.path ('Plots','HOXB13_cells_selection_umap.pdf'))
DimPlot (srt, group.by = 'HOXB13_cells_selection')
fp (srt, c('MBP','NFATC1'))
dev.off()




res_p11s = cor (t(metacells_assay_P11s), t(metacells_assay_P11s)[,gene], method= 'pearson')
res_p11s = res_p11s[order (-res_p11s[,1]),]
res_p11s = res_p11s[!is.na(res_p11s)]
#res_p11s = res_p11s[!names(res_p11s) %in% genes_in_region[[1]]]

write.csv (res_p11s, 'correlated_genes_p11s_HOXB13.csv')

# Run Enrichment on correlated genes ####
fetal_sigs = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/guccione_prj/final_fetal_sigs.csv')
fetal_sigs = as.list (fetal_sigs)
writeGMT (fetal_sigs, '/ahg/regevdata/projects/ICA_Lung/Bruno/DBs/GSEA_gs/human/fetal_sigs.gmt')

ranked_vector = res_p11s
source (file.path('..','..','git_repo','utils','fGSEA_enrichment.R'))

    message ('Print enrichment plots for each signficant pathway and cell type')
    pw = 'GO_POSITIVE_REGULATION_OF_MORPHOGENESIS_OF_AN_EPITHELIUM'
    pw = 'GO_EMBRYONIC_DIGIT_MORPHOGENESIS'
    ep = plotEnrichment(pathways[[pw]],
             ranked_vector) + 
                     labs(title='HOXB3 correlated genes')
                
    
    pdf (paste0('Plots/enrichment_plots.pdf'),5,3)
    print(ep)
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



# CRC fetal signature score across samples ####
fetal_sigs = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/guccione_prj/final_fetal_sigs.csv')
fetal_sigs = as.list (fetal_sigs)
srt = ModScoreCor (
            seurat_obj = srt, 
            geneset_list = fetal_sigs, 
            cor_threshold = NULL, 
            pos_threshold = NULL, # threshold for fetal_pval2
            listName = 'fetal_', outdir = paste0(projdir,'Plots/'))

ccomp_df = srt@meta.data[,c(names(fetal_sigs),'sampleID'), drop=FALSE]
      #ccomp_df = aggregate (ccomp_df, by=as.list(srt_wgcna@meta.data[,metaGroupNames,drop=F]), mean)    
bp1 = lapply (names(fetal_sigs), function(x) {
            ggplot (ccomp_df, aes_string (x= 'sampleID', y= x)) +
        #geom_violin (trim=TRUE, aes_string (fill = metaGroupNames[3])) +
        geom_violin (aes_string(fill = 'sampleID')) +
        geom_boxplot(width=0.5, color="black", alpha=0.2) +
        #geom_bar (stats='identity') +
        #geom_jitter (color="black", size=0.4, alpha=0.9) +
        theme_classic() + 
        scale_fill_manual (values= palette_sample) + 
        ggtitle (paste(x,'mod score')) + 
        theme (axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + NoLegend()
      })

  
png (file.path('Plots',paste0('fetal_score_sampleID_boxplot.png')),5000,5000,res=300)
print (wrap_plots (bp1))
dev.off()

srt = FindClusters (srt)

png (file.path('Plots',paste0('fetal_score_sampleID_featp.png')),height=2000,5000,res=300)
reductionName = 'umap'
wrap_plots (
    fp (srt, gene = 'fetal2', reduction = reductionName)[[1]],
    fp (srt, gene = 'fetal3', reduction = reductionName)[[1]],
    fp (srt, gene = 'fetal4', reduction = reductionName)[[1]],
    fp (srt, gene = 'TACSTD2', reduction = reductionName)[[1]],
    fp (srt, gene = 'HOXB13', reduction = reductionName)[[1]],
    DimPlot (srt, group.by = 'sampleID', cols = palette_sample,label=T),
    DimPlot (srt, group.by = 'seurat_clusters', cols = palette_sample,label=T))
dev.off()


# De novo marker discovery ####
#srt$PVLAP_deg = ifelse (srt$celltype == 'COL4A1','PLVAP','rest')
org='human'
enricher_universe = 'all'
logfcThreshold = .25
pvalAdjTrheshold = 0.01
metaGroupName = 'seurat_clusters'
#metaGroupName = c('celltype_simplified')
top_pathways = 10
top_genes = 5
force = F
source (file.path('..','..','git_repo','utils','DEG_standard.R'))








### Correlate HOX genes with CNV in P11 HOX+ cluster ####
p11cnv = readRDS (file.path('..','..','per_sample_QC_signac','CNV_analysis','CNV_LFC_GC_P11_ws_1e+07_ss_5e+06.rds'))
p11cnv_mat = assays(p11cnv)$counts
p11cnv_mat = scale (p11cnv_mat)

if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
# # Get deviation matrix ####
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name

# Get active TFs ####
selected_TF = read.csv ('Active_TFs.csv', row.names=1)
p11mMat = mMat[, archp$Clusters == 'C14']
#p11mMat = mMat[, archp$Sample2 == 'P11']
colnames(p11mMat) = sapply (colnames(p11mMat), function(x) unlist(strsplit(x, '#'))[2])
p11cnv_mat = p11cnv_mat[,colnames (p11mMat)]

#p11mMat = p11mMat[grepl ('^', rownames(p11mMat)), ]

# correlate
p11mMat = t(p11mMat)
p11cnv_mat = t(p11cnv_mat)

cnv_hox_cor = cor (as.matrix(p11mMat), as.matrix(p11cnv_mat), method='spearman')

hm = Heatmap  (cnv_hox_cor, cluster_columns = F, 
  column_names_gp = gpar(fontsize = 4),
  row_names_gp = gpar(fontsize = 4))
pdf (file.path ('Plots','cnv_HOX_P11_cor_heatmap.pdf'), width=20, height=30)
hm
dev.off()

min_cor = cnv_hox_cor['HOXB13',]
min_cor = min_cor[order(min_cor)]

pdf (file.path ('Plots','cnv_TF_cor.pdf'))
plot (p11mMat[,'HOXB13'], p11cnv_mat[,names(min_cor)[1]])
dev.off()
cor (p11mMat[,'HOXB13'], p11cnv_mat[,names(min_cor)[1]], method = 'spearman')

cnv_tf_cor = p11cnv_mat[,names(min_cor)[1]]
names(cnv_tf_cor) = paste0('P11#',names(cnv_tf_cor))
archp$p11_hox_cnv = cnv_tf_cor[match(rownames(archp@cellColData), names(cnv_tf_cor))]
archp$p11_hox_cnv[is.na(archp$p11_hox_cnv)] = 0
p <- plotEmbedding(
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = 'p11_hox_cnv', 
    embedding = "UMAP",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)

pdf (file.path('Plots','cnv_TF_cor_umap.pdf'))
p
dev.off()









