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
set.seed(1234)

# Set project dir
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scrna/'
dir.create (paste0(projdir,'/Plots/'), recursive =T)
setwd (projdir)
source ('../../PM_scATAC/useful_functions.R')

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
	srt = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/reproduction2/scRNA/srt_tumor.rds')
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
	
	srt[['SCT']] = NULL
	srt[['integrated']] = NULL	
	saveRDS (srt, 'scRNA_meso.rds')
	
	pdf ('Plots/umap_samples.pdf')
	DimPlot (srt, group.by = 'sampleID')
	dev.off()
	} else {
	srt = readRDS ('srt.rds')	
	}
  
#### Run cNMF ####
nfeat = 5000
force=F
k_list = c(5:30)
#k_list=10
k_selections = c(5:30)
#nfeat = 10000
cores= 100

### RUN consensus NMF ####
# conda create --yes --channel bioconda --channel conda-forge --channel defaults python=3.7 fastcluster matplotlib numpy palettable pandas scipy 'scikit-learn>=1.0' pyyaml 'scanpy>=1.8' -p /ahg/regevdata/projects/ICA_Lung/Bruno/conda/cnmf && conda clean --yes --all # Create environment, cnmf_env, containing required packages
# conda activate cnmf
# pip install cnmf
cnmf_name = 'scrna_tumor'
cnmf_out = paste0('cNMF/cNMF_',cnmf_name,'_',paste0(k_list[1],'_',k_list[length(k_list)]),'_vf',nfeat)
dir.create (paste0(cnmf_out,'/Plots/'), recursive=T)

# Extract and save count and data matrices from seurat object	
count_mat = t(srt@assays$RNA@counts[vf,])
norm_mat = t(srt@assays$RNA@data[vf,])
if (!file.exists (paste0('cNMF/counts_nmf_',nfeat,'.txt')) | force) write.table (count_mat, paste0('cNMF/counts_nmf_',nfeat,'.txt'), sep='\t', col.names = NA)
if (!file.exists (paste0('cNMF/norm_nmf_',nfeat,'.txt')) | force) write.table (norm_mat, paste0('cNMF/norm_nmf_',nfeat,'.txt'), sep='\t', col.names = NA)

## Format k_list variable to be passed correctly to bash script
message ('submit cNMF job')
if(length(k_list) > 1) 
	{
	k_list_formatted = paste (k_list, collapse=' ')	
	k_list_formatted = shQuote (k_list_formatted)
	}


# Run cNMF prepare script
system (paste0('chmod +x ../../PM_scATAC/cnmf_prepare_job.sh'), wait=FALSE) 
system (paste0('bash ','../../PM_scATAC/cnmf_prepare_job.sh ', projdir,' ', k_list_formatted,' ', nfeat, ' ', cnmf_out, ' ', cores), wait=TRUE)

# Submitting cNMF factorization job
system (paste0('chmod +x ../../PM_scATAC/cnmf_factorization_parallel.sh'), wait=FALSE)
system (paste0 ('qsub -t 1-',cores,' ','../../PM_scATAC/cnmf_factorization_parallel.sh ', projdir,' ', cnmf_out, ' ', cores))

# Run script to combine K iterations generated in previous script
system (paste0('chmod +x ../../PM_scATAC/cnmf_combine_job.sh'), wait=FALSE)
system (paste0('bash ','../../PM_scATAC/cnmf_combine_job.sh ', projdir, ' ',cnmf_out), wait=TRUE) # combine cnmf factors

k_selection = 25
# Run cNMF bash script
system (paste0('chmod +x ','../../PM_scATAC/cnmf_consensus_job.sh'), wait=FALSE)
system (paste0('qsub ','../../PM_scATAC/cnmf_consensus_job.sh ', projdir, ' ', cnmf_out,' ', k_selection), wait=FALSE)

# Read in NMF results ####
cnmf_spectra = read.table (paste0(cnmf_out,'/cnmf/cnmf.spectra.k_',k_selection,'.dt_0_3.consensus.txt'))

# Assign genes uniquely to cNMF modules based on spectra values
cnmf_spectra = t(cnmf_spectra)
max_spectra = apply (cnmf_spectra, 1, which.max)

top_nmf_genes = 50
cnmf_spectra_unique = lapply (1:ncol(cnmf_spectra), function(x) 
      {
      tmp = cnmf_spectra[names(max_spectra[max_spectra == x]),x,drop=F]
      tmp = tmp[order(-tmp[,1]),,drop=F]
      rownames (tmp) = gsub ('\\.','-', rownames (tmp))
      head(rownames(tmp),top_nmf_genes)
      })
names(cnmf_spectra_unique) = paste0('cNMF',seq_along(cnmf_spectra_unique))

saveRDS (cnmf_spectra_unique, paste0('cnmf_genelist_',k_selection,'_nfeat_',nfeat,'.rds'))
write.csv (patchvecs(cnmf_spectra_unique), paste0('cnmf_genelist_',k_selection,'_nfeat_',nfeat,'.csv'))

srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = cnmf_spectra_unique, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'cNMF_', outdir = paste0(projdir,'Plots/'))

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

  
png (paste0(cnmf_out,'/Plots/cNMF_module_scores_boxplots_',k_selection,'nfeat_',nfeat,'.png'),5000,5000,res=300)
print (wrap_plots (bp1))
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
  
png (paste0(cnmf_out,'/Plots/cNMF_module_scores_umaps_',k_selection,'nfeat_',nfeat,'.png'),10000,10000,res=300)
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
'c5.bp.v7.1.symbol.gmt'
)
gmt_annotation = gmt_annotations[2]
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
activeTFs = read.csv ('../scatac_ArchR/Active_TFs.csv')

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
  k = 30, # nearest-neighbors parameter
  max_shared = 20, # maximum number of shared cells between two metacells
  ident.group = 'sampleID' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
srt <- NormalizeMetacells(srt)
metacells = GetMetacellObject (srt)

metacells = metacells[activeTFs[[2]], ]
cor_TF_l = list()
for (sam in unique(metacells$sampleID))
  {
  cor_TF_l[[sam]] = cor (t(as.matrix(metacells@assays$RNA@layers$data[,metacells$sampleID == sam])))
  rownames (cor_TF_l[[sam]]) = activeTFs[[2]]
  colnames (cor_TF_l[[sam]]) = activeTFs[[2]]
  cor_TF_l[[sam]][is.na(cor_TF_l[[sam]])] = 0
  cor_TF_l[[sam]] = Heatmap (cor_TF_l[[sam]], name = sam,
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 5))
  }

pdf (paste0 ('Plots/selected_TF_exp_corr_heatmaps.pdf'), width = 8,height=9)
cor_TF_l
dev.off()

#
metacells = GetMetacellObject (srt)
vf = VariableFeatures (FindVariableFeatures (metacells, nfeat=10000))
metacells_assay = metacells@assays$RNA@layers$data
rownames (metacells_assay) = rownames(srt)
metacells_assay = metacells_assay[unique(c(vf, activeTFs[[2]])),]

enricher_universe = vf
#do.fgsea = TRUE
gmt_annotations = c(
'h.all.v7.4.symbols.gmt',#,
'c5.bp.v7.1.symbol.gmt'
)
gmt_annotation = gmt_annotations[1]

if (!file.exists(paste0('EnrichR_activeTF_cor_top_genes_,_ann_',gmt_annotation,'.rds')))
	{
	TF_cor_sample = list()
	for (sam in unique(metacells$sampleID))
	  {
	  metacells_assay_sample = metacells_assay[,metacells$sampleID == sam]
	  TF_cor_sample[[sam]] = lapply (activeTFs[[2]], function(x) 
	  	{
	  	tc_cor = t(cor (metacells_assay_sample[x,], t(metacells_assay_sample)))
	  	tc_cor = setNames (tc_cor[,1], rownames (tc_cor))
	  	tc_cor[is.na(tc_cor)] = 0
	  	tc_cor = tc_cor[order(-tc_cor)]
	  	gmt.file = paste0 ('../../PM_scATAC/files/',gmt_annotation)
	  	pathways = read.gmt (gmt.file)
	  	#pathways = split (pathways$gene, pathways$term)
	    message (paste ('EnrichR running module',x)) 
	    egmt <- enricher(head (names(tc_cor),100), TERM2GENE=pathways, universe = enricher_universe)
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
	saveRDS (TF_cor_sample, paste0('EnrichR_activeTF_cor_top_genes_,_ann_',gmt_annotation,'.rds'))
	} else {
	TF_cor_sample = readRDS (paste0('EnrichR_activeTF_cor_top_genes_,_ann_',gmt_annotation,'.rds'))
	}
pvalAdjTrheshold = 0.05
top_pathways = 10
TF_cor_sample = lapply (TF_cor_sample, function(x) {names(x) = activeTFs[[2]]; x})
EnrichRes_dp = lapply (TF_cor_sample, function(x) dotGSEA (enrichmentsTest_list = x, type = 'enrich', padj_threshold = pvalAdjTrheshold, top_pathways= top_pathways))
pdf (paste0('Plots/EnrichR_top_nmf_genes_ann_',gmt_annotation,'_dotplots.pdf'), width = 20, height = 15)
print (EnrichRes_dp)
dev.off()
  

pdf (paste0 ('Plots/selected_TF_exp_corr_heatmaps.pdf'), width = 8,height=9)
cor_TF_l
dev.off()


metacells = ModScoreCor (
        seurat_obj = metacells, 
        geneset_list = cnmf_spectra_unique, 
        cor_threshold = NULL, 
        pos_threshold = NULL, # threshold for fetal_pval2
        listName = 'cNMF_', outdir = paste0(projdir,'Plots/'))

ccomp_df = metacells@meta.data[,c(names(cnmf_spectra_unique))]

tc_cor = lapply (unique(metacells$sampleID), function(x)
	{
	metacells_assay_sample = metacells_assay[activeTFs[[2]],metacells$sampleID == x]
	ccomp_df_sample = ccomp_df[metacells$sampleID == x,]
	res = sapply (activeTFs[[2]], function (y) cor (metacells_assay_sample[y,], ccomp_df_sample, method = 'spearman'))
	rownames (res) = colnames (ccomp_df_sample)
	res
	})
names (tc_cor) = unique(metacells$sampleID)
#lapply (tc_cor, function(x) {x = x['cNMF19',]; head(x[order(-x)],10)})


cor_TF <- tc_cor[[1]]
cor_TF[] <- tapply(unlist(tc_cor), rep(seq(length(tc_cor[[1]])),length(tc_cor)), FUN=function(x) median(x,na.rm=T))

dim (cor_TF)

pdf ('Plots/cor_nmf_TF.pdf',width=10, height=3)
Heatmap (cor_TF, 
	column_names_gp = gpar(fontsize = 5),
	row_names_gp = gpar(fontsize = 5),
	clustering_distance_rows = 'pearson',
	clustering_distance_columns='pearson')
dev.off()



# Order cells per samples along SOX9 expression and plot the rest of TF expression together
mMat = metacells@assays$RNA@layers$data[rownames(metacells) %in% activeTFs[[2]],]
rownames(mMat) = rownames(metacells)[rownames(metacells) %in% activeTFs[[2]]]

traj_sample = list()
for (sam in unique(metacells$sampleID))
    {
    mMat_ordered_sample = mMat[,metacells$sampleID == sam]
    mMat_ordered_sample = mMat_ordered_sample[, order(mMat_ordered_sample['SOX9',])]
    mMat_ordered_sample_sclaled = as.matrix(mMat_ordered_sample)
    mMat_ordered_sample_sclaled[is.na(mMat_ordered_sample_sclaled)] = 0
    traj_sample[[sam]] = Heatmap (
      as.matrix(mMat_ordered_sample_sclaled), 
      col = viridis::plasma(100), 
      cluster_columns=F, name = sam,
      row_names_gp = gpar(fontsize = 4))
    }

pdf ('Plots/sarc_trajectory_per_sample.pdf', height=10)
traj_sample
dev.off()








