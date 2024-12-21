# Variables
# chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet'
# chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/chromBPnet'
# repodir='/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
# grefdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet'
# celltype='SPP1'
# celltype='Mono'
# celltype='C1Q'
# celltype='P23__epithelioid'
# celltype = 'C2'
# celltype = 'C2_test'
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
		command <- paste("bsub -J", paste0(celltype,'_cBP'), "-P acc_Tsankov_Normal_Lung -q gpu -n 8 -W 48:00 -gpu num=2 -R v100 -R rusage[mem=32000] -R span[hosts=1] -o",file.path(chromBPdir,paste0('chromBPtraining_',celltype,'.out')), "-eo" ,file.path(chromBPdir,paste0('chromBPtraining_',celltype,'.err')),file.path(repodir,'utils','chrBPnet_training.sh'))
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

if (!all (file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),paste0(celltype,'_contribution_scores.counts_scores.bw')))))
	{
	fold_numbers_remained = fold_numbers[!file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),paste0(celltype,'_contribution_scores.counts_scores.bw')))]	
	for (fold_number in fold_numbers_remained)
		{
		message (paste0('submit job for contribution fold ',fold_number))
		#system (paste0('chmod +x ',file.path(repodir, 'utils','chrBPnet_training.sh')))
		#paste("bsub bash -c \"/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f3.sh, projdir)
		#system (paste("bsub bash -c \", file.path(repodir,utils','chrBPnet_training.sh '), chromBPdir,' ',grefdir,' ',repodir,' ',celltype,' ',fold_number,'"'), wait=FALSE)
		command <- paste("bsub -J", paste0(celltype,'_cBPcntr'), "-P acc_Tsankov_Normal_Lung -q gpu -n 8 -W 72:00 -gpu num=1 -R v100 -R rusage[mem=32000] -R span[hosts=1] -o",file.path(chromBPdir,paste0('chromBPcontribution_',celltype,'.out')), "-eo" ,file.path(chromBPdir,paste0('chromBPcontribution_',celltype,'.err')),file.path(repodir,'utils','chrBPnet_contribution.sh'))
		args <- paste(chromBPdir, grefdir, celltype, fold_number)
		system (paste(command, args))
		}
	} else {
	message ('all contribution score folds have been computed!')	
	}


system (paste0('chmod +x ',file.path(repodir,'utils','TFmodisco.sh')), wait=FALSE) 

if (!all (file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),'modisco','modisco_results.h5'))))
	{
	fold_numbers_remained = fold_numbers[!file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),'TF_modisco','_modisco_results.h5'))]
	# fold_numbers_remained = c(0:4)
	for (fold_number in fold_numbers_remained)
		{
		message (paste0('submit job for TFmodisco fold ',fold_number))
		#fold_number = 4
		#system (paste0('chmod +x ',file.path(repodir, 'utils','chrBPnet_training.sh')))
		#paste("bsub bash -c \"/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f3.sh, projdir)
		#system (paste("bsub bash -c \", file.path(repodir,utils','chrBPnet_training.sh '), chromBPdir,' ',grefdir,' ',repodir,' ',celltype,' ',fold_number,'"'), wait=FALSE)
		command <- paste ("bsub -J", paste0(celltype,'_TFmd'), "-P acc_Tsankov_Normal_Lung -q premium -n 8 -W 96:00 -R rusage[mem=64000] -R span[hosts=1] -o",file.path(chromBPdir,paste0('TFmodisco_',celltype,'.out')), "-eo" ,file.path(chromBPdir,paste0('TFmodisco_',celltype,'.err')),file.path(repodir,'utils','TFmodisco.sh'))
		args <- paste(chromBPdir,celltype, fold_number)
		system (paste(command, args))
		}
	}

system (paste0('chmod +x ',file.path(repodir,'utils','finemo_motif_calls.sh')), wait=FALSE) 

if (!all (file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),'modisco_profile','modisco_results.h5'))))
	{
	fold_numbers_remained = fold_numbers[!file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),'TF_modisco','_modisco_results.h5'))]	
	# fold_numbers_remained = c(0:4)
	for (fold_number in fold_numbers_remained)
		{
		message (paste0('submit job for finemo fold ',fold_number))
		#fold_number = 4
		#system (paste0('chmod +x ',file.path(repodir, 'utils','chrBPnet_training.sh')))
		#paste("bsub bash -c \"/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f3.sh, projdir)
		#system (paste("bsub bash -c \", file.path(repodir,utils','chrBPnet_training.sh '), chromBPdir,' ',grefdir,' ',repodir,' ',celltype,' ',fold_number,'"'), wait=FALSE)
		command <- paste("bsub -J", paste0(celltype,'_finemo'), "-P acc_Tsankov_Normal_Lung  -q gpu -n 1 -W 12:00 -gpu num=1 -R v100 -R rusage[mem=32000] -R span[hosts=1] -o",file.path(chromBPdir,paste0('finemo_',celltype,'.out')), "-eo" ,file.path(chromBPdir,paste0('finemo_',celltype,'.err')),file.path(repodir,'utils','finemo_motif_calls.sh'))
		args <- paste(chromBPdir,celltype, fold_number)
		system (paste(command, args))
		}
	}





message ('Submit job for chromBPnet training model and contribution scores')	

# Make script executable
system (paste0('chmod +x ',file.path(repodir,'utils','chrBPnet_training_new.sh')), wait=FALSE) 

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
		command <- paste("bsub -J", paste0(celltype,'_cBP'), "-P acc_Tsankov_Normal_Lung -q gpu -n 8 -W 96:00 -gpu num=2 -R h100nvl -R rusage[mem=32000] -R span[hosts=1] -o",file.path(chromBPdir,paste0('chromBPtraining_',celltype,'.out')), "-eo" ,file.path(chromBPdir,paste0('chromBPtraining_',celltype,'.err')),file.path(repodir,'utils','chrBPnet_training_new.sh'))
		args <- paste(chromBPdir, grefdir, celltype, fold_number)
		system (paste(command, args))
		}
	} else {
	message ('all training model folds have been trained!')	
	}



### Combine contribution files together 
library (rhdf5)
file_path = 'chromBPnet/P23__epithelioid_model/fold_0/P23__epithelioid_contribution_scores.counts_scores.h5'
h5ls(file_path)






command <- paste ("bsub -J", paste0(celltype,'_cBPm'), 
	"-P acc_Tsankov_Normal_Lung -q premium -n 8 -W 96:00 -R rusage[mem=32000] -R span[hosts=1] -o",
	file.path(chromBPdir,paste0('cBP_master_',celltype,'.out')), "-e" ,
	file.path(chromBPdir,paste0('cBP_master_',celltype,'.err')),
	file.path(repodir,'utils','chromBPnet_master.sh'))
args <- paste(chromBPdir, grefdir, repodir, celltype) 
system (paste0('chmod +x ',file.path(repodir,'utils','chromBPnet_master.sh')), wait=FALSE) 
system (paste0('chmod +x ',file.path(repodir,'utils','chrBPnet_training_new.sh')), wait=FALSE) 
system (paste(command, args))

