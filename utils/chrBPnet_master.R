# Variables
# chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet'
# repodir='/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
# grefdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet'
# celltype='CD8'
# fold_numbers = c(0,1,2,3,4)

message ('Submit job for chromBPnet training model')	

# Make script executable
system (paste0('chmod +x ',file.path(repodir,'utils','chrBPnet_training.sh')), wait=FALSE) 

# Check if models have been already generated if not submit jobs 
if (!all (file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),'evaluation','overall_report.html'))))
	{
	fold_numbers_remained = fold_numbers[!file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),'evaluation','overall_report.html'))]	
	for (fold_number in fold_numbers_remained)
		{
		message (paste0('submit job for training fold ',fold_number))
		#system (paste0('chmod +x ',file.path(repodir, 'utils','chrBPnet_training.sh')))
		#paste("bsub bash -c \"/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f3.sh, projdir)
		#system (paste("bsub bash -c \", file.path(repodir,utils','chrBPnet_training.sh '), chromBPdir,' ',grefdir,' ',repodir,' ',celltype,' ',fold_number,'"'), wait=FALSE)
		command <- paste("bsub -J", paste0(celltype,'_chrBP_model'), "-P acc_Tsankov_Normal_Lung -q gpu -n 8 -W 48:00 -gpu num=1 -R v100 -R rusage[mem=32000] -R span[hosts=1] -o",file.path(chromBPdir,paste0('chormBPtraining_',celltype,'.out')), "-eo" ,file.path(chromBPdir,paste0('chormBPtraining_',celltype,'.err')),file.path(repodir,'utils','chrBPnet_training.sh'))
		args <- paste(chromBPdir, grefdir, celltype, fold_number)
		system (paste(command, args))
		}
	} else {
	message ('all training model folds have been trained!')	
	}

# Check if all folds have been completed again #
data.frame (fold = fold_numbers, completed = file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),'evaluation','overall_report.html')))

message ('Submit job for chromBPnet contribution scores')	
system (paste0('chmod +x ',file.path(repodir,'utils','chrBPnet_contribution.sh')), wait=FALSE) 

if (!all (file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),'contribution_scores',paste0(celltype,'_contribution_scores.counts_scores.bw')))))
	{
	fold_numbers_remained = fold_numbers[!file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),'contribution_scores',paste0(celltype,'_contribution_scores.counts_scores.bw')))]	
	for (fold_number in fold_numbers_remained)
		{
		message (paste0('submit job for contribution fold ',fold_number))
		#system (paste0('chmod +x ',file.path(repodir, 'utils','chrBPnet_training.sh')))
		#paste("bsub bash -c \"/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f3.sh, projdir)
		#system (paste("bsub bash -c \", file.path(repodir,utils','chrBPnet_training.sh '), chromBPdir,' ',grefdir,' ',repodir,' ',celltype,' ',fold_number,'"'), wait=FALSE)
		command <- paste("bsub -J", paste0(celltype,'_chrBP_contribution'), "-P acc_Tsankov_Normal_Lung -q gpu -n 8 -W 48:00 -gpu num=1 -R v100 -R rusage[mem=32000] -R span[hosts=1] -o",file.path(chromBPdir,paste0('chormBPcontribution_',celltype,'.out')), "-eo" ,file.path(chromBPdir,paste0('chormBPcontribution_',celltype,'.err')),file.path(repodir,'utils','chrBPnet_contribution.sh'))
		args <- paste(chromBPdir, grefdir, celltype, fold_number)
		system (paste(command, args))
		}
	} else {
	message ('all contribution score folds have been computed!')	
	}



