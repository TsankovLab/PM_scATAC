# Variables
# chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet'
# repodir='/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
# grefdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet'
# celltype='CD8'
# fold_numbers = c(0,1,2,3,4)

message ('Submit job for chromBPnet training model')	
# Extract and save count and data matrices from seurat object	
## Format k_list variable to be passed correctly to bash script

system (paste0('chmod +x ',file.path(repodir,'utils','chrBPnet_training.sh')), wait=FALSE) 

if (!all (file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),'evaluation','overall_report.html'))))
	fold_numbers_remained = fold_numbers[!file.exists(file.path (chromBPdir, paste0(celltype,'_model'), paste0('fold_',fold_numbers),'evaluation','overall_report.html'))]	
	{
	for (fold_number in fold_numbers_remained)
		{
		message (paste0('submit job for training fold ',fold_number))
		#system (paste0('chmod +x ',file.path(repodir, 'utils','chrBPnet_training.sh')))
		#paste("bsub bash -c \"/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/tnk_analysis/chrombnet_NKT_NK_KLRC1_f3.sh, projdir)
		#system (paste("bsub bash -c \", file.path(repodir,utils','chrBPnet_training.sh '), chromBPdir,' ',grefdir,' ',repodir,' ',celltype,' ',fold_number,'"'), wait=FALSE)
		command <- paste("bsub -J", paste0(celltype,'_chrBP_model'), "-P acc_Tsankov_Normal_Lung -q gpu -n 8 -W 48:00 -gpu num=1 -R v100 -R rusage[mem=32000] -R span[hosts=1] -o","/sc/arion/projects/Tsankov_Normal_Lung/Bruno/output_logs/job_output.out", "-eo" ,paste0(celltype,'_stderr'))
		args <- paste(chromBPdir, grefdir, celltype, fold_number)
		system (paste(command, args))
		)
		}
	} else {
	message ('all training model folds have been trained!')	
	}


