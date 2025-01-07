conda activate meso_scatac
R

set.seed(1234)

####### ANALYSIS of NKT compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

# Load functions for hub detection ####
source (file.path('..','..','git_repo','utils','knnGen.R'))
source (file.path('..','..','git_repo','utils','addCoax.R'))
source (file.path('..','..','git_repo','utils','Hubs_finder.R'))
source (file.path('..','..','git_repo','utils','hubs_track.R'))
#source (file.path('..','..','git_repo','utils','scATAC_functions.R'))

# Set # of threads and genome reference ####
addArchRThreads(threads = 1) 
addArchRGenome("hg38")

archp = loadArchRProject (projdir)

### Run chromBPnet bias model or provide folder of bias_model.h5 file ####
archp$TNK_cells = 'TNK_cells'

### Call peaks with MACS2 by metaGroupName ####
metaGroupName = 'TNK_cells'
source ('../../git_repo/chromBPnet_call_peaks.R')

# Note: you need to have a bias model for each of the fold used in the no bias model
chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet'
repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
grefdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet'
celltype = 'TNK_cells'

command <- paste ("bsub -J", paste0(celltype,'_cBP_bias'), 
	"-P acc_Tsankov_Normal_Lung -q premium -n 8 -W 96:00 -R rusage[mem=32000] -R span[hosts=1] -o",
	file.path(chromBPdir,paste0('cBP_master_bias_',celltype,'.out')), "-e" ,
	file.path(chromBPdir,paste0('cBP_master_bias_',celltype,'.err')),
	file.path(repodir,'utils','chromBPnet_bias_master.sh'))
args <- paste(chromBPdir, grefdir, repodir, celltype) 
system (paste0('chmod +x ',file.path(repodir,'utils','chromBPnet_bias_master.sh')), wait=FALSE)

system (paste(command, args))



# Run no bias chromBPnet model for each NKT cell subtype ####

### Call peaks with MACS2 by metaGroupName ####
metaGroupName = 'celltype2'
source ('../../git_repo/utils/chromBPnet_call_peaks.R')

chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet'
repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
grefdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet'
biasdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet/NKT_cells/bias_model'

celltypes = unique (as.character(archp@cellColData[, metaGroupName]))

for (celltype in celltypes)
	{
	command <- paste ("bsub -J", paste0(celltype,'_cBPnet'), 
		"-P acc_Tsankov_Normal_Lung -q premium -n 8 -W 96:00 -R rusage[mem=32000] -R span[hosts=1] -o",
		file.path(chromBPdir,paste0('cBP_master_',celltype,'.out')), "-e" ,
		file.path(chromBPdir,paste0('cBP_master_',celltype,'.err')),
		file.path(repodir,'utils','chromBPnet_master.sh'))
	args <- paste(chromBPdir, grefdir, repodir, celltype, biasdir) 
	system (paste0('chmod +x ',file.path(repodir,'utils','chromBPnet_master.sh')), wait=FALSE)
	
	system (paste(command, args))
	}







