use UGER # Add this before running R to be able to run cNMF scripts using UGER 
library (Seurat)
library (scran)
library (ggplot2)
library (RColorBrewer)
library(patchwork)
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
	saveRDS (srt, 'scRNA_meso.rds')
	
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
	
	pdf ('Plots/umap_samples.pdf')
	DimPlot (srt, group.by = 'sampleID')
	dev.off()
	}
  
#### Run cNMF ####
force=F
k_list = c(5:30)
#k_list=10
k_selections = c(5:30)
#nfeat = 10000
cores= 100

cnmf_name='scran_hvg'

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

k = 25
# Run cNMF bash script
system (paste0('chmod +x ','../../PM_scATAC/cnmf_consensus_job.sh'), wait=FALSE)
system (paste0('qsub ','../../PM_scATAC/cnmf_consensus_job.sh ', projdir, ' ', cnmf_out,' ', k), wait=FALSE)

# Read in NMF results ####
cnmf_spectra = read.table (paste0(cnmf_out,'/cnmf/cnmf.spectra.k_',k,'.dt_0_3.consensus.txt'))

# Assign genes uniquely to cNMF modules based on spectra values
cnmf_spectra = t(cnmf_spectra)
max_spectra = apply (cnmf_spectra, 1, which.max)

top_nmf_genes = Inf
cnmf_spectra_unique = lapply (1:ncol(cnmf_spectra), function(x) 
      {
      tmp = cnmf_spectra[names(max_spectra[max_spectra == x]),x,drop=F]
      tmp = tmp[order(-tmp[,1]),,drop=F]
      rownames (tmp) = gsub ('\\.','-', rownames (tmp))
      head(rownames(tmp),top_nmf_genes)
      })
names(cnmf_spectra_unique) = paste0('cNMF',seq_along(cnmf_spectra_unique))

saveRDS (cnmf_spectra_unique, paste0('cnmf_genelist_',k,'_nfeat_',nfeat,'.rds'))
write.csv (patchvecs(cnmf_spectra_unique), paste0('cnmf_genelist_',k,'_nfeat_',nfeat,'.csv'))

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

  
png (paste0(cnmf_out,'/Plots/cNMF_module_scores_boxplots_',k,'nfeat_',nfeat,'.png'),5000,5000,res=300)
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
  
png (paste0(cnmf_out,'/Plots/cNMF_module_scores_umaps_',k,'nfeat_',nfeat,'.png'),10000,10000,res=300)
print (wrap_plots (umap_p1))
dev.off()


