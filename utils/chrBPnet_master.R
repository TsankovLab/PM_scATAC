# Variables
# chromBPdir
# repodir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/
# celltype
# fold_numbers

message ('Submit job for chromBPnet training model')	
# Extract and save count and data matrices from seurat object	
## Format k_list variable to be passed correctly to bash script

system (paste0('chmod +x ',file.path(repodir,'utils','chrBPnet_training.sh')), wait=FALSE) 

if (!all (file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),'evaluation','overall_report.html')))
	fold_numbers_remained = fold_numbers[!file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),'evaluation','overall_report.html'))]	
	{
	for (fold_number in fold_numbers_remained)
		{
		message (paste0('submit job for training fold ',fold_number))
		system (paste0('chmod +x ',file.path(repodir, 'utils','chrBPnet_training.sh')))
		system (paste0('bsub <',file.path(repodir,'utils','chrBPnet_training.sh'), chromBPdir,' ',grefdir,' ',repodir,' ',celltype,' ',fold_number), wait=FALSE)
		}
	} else {
	message ('all training model folds have been trained!')	
	}


