use UGER # Add this before running R to be able to run cNMF scripts using UGER 
library (Seurat)
library (scran)
set.seed(1234)

# Set project dir
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scrna/'
dir.create (paste0(projdir,'/Plots/'), recursive =T)
setwd (projdir)

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
srt = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/reproduction2/scRNA/srt_tumor.rds')
srt$celltype_simplified_mal = as.character (srt$celltype_simplified)
srt$celltype_simplified_mal[srt$celltype_simplified == 'Malignant'] = paste0(srt$sampleID[srt$celltype_simplified == 'Malignant'], '_', srt$celltype_simplified[srt$celltype_simplified == 'Malignant'])
# Import P14
srt_p14 = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_p14_analysis/_cellranger_raw_Filter_400_1000_25/no_harmony/srt.rds')
srt_p14$celltype_simplified = srt_p14$celltype
srt = merge (srt, srt_p14[intersect(rownames(srt_p14), rownames(srt)),])
srt$sampleID[srt$sampleID == 'MPM_naive_p14'] = 'P14'
srt = srt[, srt$celltype_simplified == 'Malignant']
srt = srt[, srt$sampleID %in% sample_names]
srt = NormalizeData (srt)

sce = SingleCellExperiment (list(counts=srt@assays$RNA@counts, logcounts = srt@assays$RNA@data),
rowData=rownames(srt)) 
sce = modelGeneVar(sce)
# remove batchy genes
batchy_genes = c('RPL','RPS','MT-')
sce = sce[!apply(sapply(batchy_genes, function(x) grepl (x, rownames(sce))),1,any),]
nfeat = 3000
vf = getTopHVGs(sce, n=nfeat)
VariableFeatures (srt) = vf


  

#### Run cNMF ####
force=F
k_list = c(5:30)
#k_list=10
k_selections = c(5:30)
#nfeat = 10000
cores= 100
metacells = FALSE
min_cells = 20 # default 100
metacells_k = 30
max_shared=20
metacells_groups = 'sampleID'
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

system (paste0('chmod +x ../../PM_scATAC/cnmf_prepare_job.sh'), wait=FALSE) 
system (paste0('chmod +x ../../PM_scATAC/cnmf_factorization_parallel.sh'), wait=FALSE)

# Run cNMF prepare script
system (paste0('bash ','../../PM_scATAC/cnmf_prepare_job.sh ', projdir,' ', k_list_formatted,' ', nfeat, ' ', cnmf_out, ' ', cores), wait=TRUE)

# Submitting cNMF factorization job
system (paste0 ('qsub -t 1-',cores,' ','../../PM_scATAC/cnmf_factorization_parallel.sh ', projdir,' ', cnmf_out, ' ', cores))


# Run script to combine K iterations generated in previous script
message ('combine cnmf factors')
system (paste0('bash ','../../PM_scATAC/cnmf_combine_job.sh ', projdir, ' ',cnmf_out), wait=TRUE)

for (i in k_selections)
	{
	if (file.exists (paste0(cnmf_out,'/cnmf/cnmf.k_selection.png')) & !file.exists (paste0(projdir, cnmf_out,'/cnmf/cnmf.spectra.k_',i,'.dt_0_3.consensus.txt')) | force)
		{
		message ('submit cNMF consensus job')
		# Make sh file executable
		system (paste0('chmod +x ',scrna_pipeline_dir,'cnmf_consensus_job.sh'), wait=FALSE)
		
		# Run cNMF bash script
		#system (paste0('qsub -t 1-4 ',scrna_pipeline_dir,'cnmf_prepare_job.sh ', projdir,' ', k_list,' ', k_selection,' ', nfeat, ' ', cnmf_out), wait=FALSE)
		system (paste0('qsub ',scrna_pipeline_dir,'cnmf_consensus_job.sh ', projdir, ' ', cnmf_out,' ', i), wait=FALSE)
		}
	}

if (all (file.exists (paste0(cnmf_out,'/cnmf/cnmf.spectra.k_',k_selections,'.dt_0_3.consensus.txt'))))
	{
	message ('clean files')	
	system (paste0('rm -r ', projdir,'CNMF_factorization*'))
	system (paste0('rm -r ', projdir,'CNMF_consensus*'))
	system (paste0('rm -r ', projdir, cnmf_out,'/cnmf/cnmf.gene_spectra_tpm*'))
	system (paste0('rm -r ', projdir, cnmf_out,'/cnmf/cnmf.gene_spectra_score*'))
	system (paste0('rm -r ', projdir, cnmf_out,'/cnmf/cnmf.clustering*'))
	system (paste0('rm -r ', projdir, cnmf_out,'/cnmf/cnmf_tmp'))
	}
