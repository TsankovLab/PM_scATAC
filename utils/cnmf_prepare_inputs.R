require (scran)

if ('layers' %in% slotNames (srt@assays$RNA))
	{
	sce = SingleCellExperiment (list(counts=srt@assays$RNA@layers$counts, logcounts = srt@assays$RNA@layers$data),
	rowData=rownames(srt))
	rownames(sce) = rownames (srt) 
	} else {
	sce = SingleCellExperiment (list(counts=srt@assays$RNA@counts, logcounts = srt@assays$RNA@data),
	rowData=rownames(srt)) 
	}
	sce = modelGeneVar(sce)
	# remove batchy genes
	batchy_genes = c('RPL','RPS','MT-')
	sce = sce[!apply(sapply(batchy_genes, function(x) grepl (x, rownames(sce))),1,any),]
	vf = getTopHVGs (sce, n=nfeat)

cnmf_out = paste0('cNMF/cNMF_',cnmf_name,'_',paste0(k_list[1],'_',k_list[length(k_list)]),'_vf',nfeat)
dir.create (file.path(cnmf_out,'Plots'), recursive=T)

# Export count matrices ####
if (!file.exists (file.path('cNMF',paste0('counts_nmf_',nfeat,'_',cnmf_name,'.txt'))) | force) 
	{
	count_mat = t(srt@assays$RNA@counts[vf,])
	write.table (count_mat, (file.path('cNMF',paste0('counts_nmf_',nfeat,'_',cnmf_name,'.txt'))), sep='\t', col.names = NA)
	}
if (!file.exists (file.path('cNMF',paste0('norm_nmf_',nfeat,'_',cnmf_name,'.txt'))) | force) 
	{
	norm_mat = t(srt@assays$RNA@data[vf,])	
	write.table (norm_mat, (file.path('cNMF',paste0('norm_nmf_',nfeat,'_',cnmf_name,'.txt'))), sep='\t', col.names = NA)	
	}
	
### RUN consensus NMF ####
# conda create --yes --channel bioconda --channel conda-forge --channel defaults python=3.7 fastcluster matplotlib numpy palettable pandas scipy 'scikit-learn>=1.0' pyyaml 'scanpy>=1.8' -p /ahg/regevdata/projects/ICA_Lung/Bruno/conda/cnmf && conda clean --yes --all # Create environment, cnmf_env, containing required packages
# conda activate cnmf
# pip install cnmf
if (!all(file.exists(file.path(cnmf_out,'cnmf',paste0('cnmf.spectra.k_',k_selections,'.dt_0_3.consensus.txt')))))
	{
	message ('Run cNMF...')	
	# Extract and save count and data matrices from seurat object	
	## Format k_list variable to be passed correctly to bash script
	message ('submit cNMF job')
	if(length(k_list) > 1) 
		{
		k_list_formatted = paste (k_list, collapse=' ')	
		k_list_formatted = shQuote (k_list_formatted)
		}
	
	
	# Run cNMF prepare script
	system (paste0('chmod +x ',file.path(repodir,'utils','cnmf_master.sh')), wait=FALSE) 
	system (paste0('bash ',file.path(repodir,'utils','cnmf_master.sh '), projdir,' ',cnmf_out,' ',repodir,' ',nfeat,' ',k_list_formatted,' ', cores,' ', cnmf_name), wait=TRUE)
	} else {
	message ('cnmf spectra files aready generated!')	
	}	
