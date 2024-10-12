require (scran)

sce = SingleCellExperiment (list(counts=srt@assays$RNA@counts, logcounts = srt@assays$RNA@data),
	rowData=rownames(srt)) 
	sce = modelGeneVar(sce)
	# remove batchy genes
	batchy_genes = c('RPL','RPS','MT-')
	sce = sce[!apply(sapply(batchy_genes, function(x) grepl (x, rownames(sce))),1,any),]
	nfeat = 5000
	vf = getTopHVGs (sce, n=nfeat)

cnmf_name = 'CD8_exhausted'
cnmf_out = paste0('cNMF/cNMF_',cnmf_name,'_',paste0(k_list[1],'_',k_list[length(k_list)]),'_vf',nfeat)
dir.create (file.path(cnmf_out,'Plots'), recursive=T)

# Export count matrices ####
if (!file.exists (file.path('cNMF',paste0('counts_nmf_',nfeat,'.txt'))) | force) 
	{
	count_mat = t(srt@assays$RNA@counts[vf,])
	write.table (count_mat, (file.path('cNMF',paste0('counts_nmf_',nfeat,'.txt'))), sep='\t', col.names = NA)
	}
if (!file.exists (file.path('cNMF',paste0('norm_nmf_',nfeat,'.txt'))) | force) 
	{
	norm_mat = t(srt@assays$RNA@data[vf,])	
	write.table (norm_mat, (file.path('cNMF',paste0('norm_nmf_',nfeat,'.txt'))), sep='\t', col.names = NA)	
	}
	
### RUN consensus NMF ####
# conda create --yes --channel bioconda --channel conda-forge --channel defaults python=3.7 fastcluster matplotlib numpy palettable pandas scipy 'scikit-learn>=1.0' pyyaml 'scanpy>=1.8' -p /ahg/regevdata/projects/ICA_Lung/Bruno/conda/cnmf && conda clean --yes --all # Create environment, cnmf_env, containing required packages
# conda activate cnmf
# pip install cnmf
if (!all(file.exists(file.path(cnmf_out,'cnmf',paste0('cnmf.spectra.k_',k_selections,'.dt_0_3.consensus.txt')))))
	{
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
	system (paste0('bash ',file.path(repodir,'utils','cnmf_master.sh '), projdir,' ',cnmf_out,' ',repodir,' ',nfeat,' ',k_list_formatted,' ', cores), wait=TRUE)
	
	} else {
	cnmf_spectra = read.table (paste0(cnmf_out,'/cnmf/cnmf.spectra.k_',k_selection,'.dt_0_3.consensus.txt'))
	}

# Format NMF results ####
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

saveRDS (cnmf_spectra_unique, paste0('cnmf_genelist_',k_selection,'_nfeat_',nfeat,'.rds'))
write.csv (patchvecs(cnmf_spectra_unique), paste0('cnmf_genelist_',k_selection,'_nfeat_',nfeat,'.csv'))

top_nmf_genes = 200
cnmf_spectra_unique = lapply (1:ncol(cnmf_spectra), function(x) 
      {
      tmp = cnmf_spectra[names(max_spectra[max_spectra == x]),x,drop=F]
      tmp = tmp[order(-tmp[,1]),,drop=F]
      rownames (tmp) = gsub ('\\.','-', rownames (tmp))
      head(rownames(tmp),top_nmf_genes)
      })
names(cnmf_spectra_unique) = paste0('cNMF',seq_along(cnmf_spectra_unique))
saveRDS (cnmf_spectra_unique, paste0(cnmf_out,'/cnmf_genelist_',k_selection,'.rds'))
