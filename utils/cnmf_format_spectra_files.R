
# Format spectra NMF results ####
cnmf_spectra = read.table (paste0(cnmf_out,'/cnmf/cnmf.spectra.k_',k_selection,'.dt_0_3.consensus.txt'))

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

message ('cnmf results in cnmf_spectra_unique object')