conda activate meso_scatac
R

set.seed(1234)

####### ANALYSIS of tumor compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR'
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

### Call peaks with MACS2 by metaGroupName ####
archp$celltype_revised_sample = paste0(archp$celltype_lv1, '_', archp$Sample)

#metaGroupName = 'celltype_lv1'
metaGroupName = 'celltype_revised_sample'
#archp = archp[archp$celltype_lv1 %in% c('Fibroblasts','Malignant','Mesothelium')]
archp = archp[archp$celltype_revised_sample %in% c('Fibroblasts_P1')]
source ('../../git_repo/utils/chromBPnet_call_peaks.R')

# Run no bias chromBPnet model for each NKT cell subtype ####
chromBPdir = file.path (projdir,'chromBPnet')
dir.create (chromBPdir)
repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
grefdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet'
biasdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet/NKT_cells/bias_model'

celltypes = unique (as.character(archp@cellColData[, metaGroupName]))
#celltypes='B_cells'

celltypes = c('Fibroblasts','Malignant','Mesothelium')
celltypes = c('Fibroblasts_P1')
for (celltype in celltypes)
	{
	command <- paste ("bsub -J", paste0(celltype,'_CBPmaster'), 
		"-P acc_Tsankov_Normal_Lung -q premium -n 1 -W 72:00 -R rusage[mem=64000] -R span[hosts=1] -o",
		file.path(chromBPdir,paste0('cBP_master_',celltype,'.out')), "-e" ,
		file.path(chromBPdir,paste0('cBP_master_',celltype,'.err')),
		file.path(repodir,'utils','chromBPnet_master.sh'))
	args <- paste (chromBPdir, grefdir, repodir, celltype, biasdir) 
	system (paste0('chmod +x ',file.path(repodir,'utils','chromBPnet_master.sh')), wait=FALSE)
	
	system (paste(command, args))
	}
