source(local({a <- grep("^--file=", commandArgs(FALSE), value = TRUE); for (d in c(if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))), getwd())) { while (!file.exists(file.path(d, "config.R")) && d != dirname(d)) d <- dirname(d); if (file.exists(file.path(d, "config.R"))) break }; Sys.getenv("MESO_CONFIG", file.path(d, "config.R"))}))   # data / results / reference roots: see config.R at the repo root
conda activate meso_scatac
R

set.seed(1234)

####### ANALYSIS of tumor compartment #######
projdir = paste0(DATA_DIR, '/main/scatac_ArchR')
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','git_repo','utils','load_packages.R'))
source (file.path('..','git_repo','utils','useful_functions.R'))
source (file.path('..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','git_repo','utils','palettes.R'))

# Load functions for hub detection ####
source (file.path('..','git_repo','utils','knnGen.R'))
source (file.path('..','git_repo','utils','addCoax.R'))
source (file.path('..','git_repo','utils','Hubs_finder.R'))
source (file.path('..','git_repo','utils','hubs_track.R'))
#source (file.path('..','..','git_repo','utils','scATAC_functions.R'))

# Set # of threads and genome reference ####
addArchRThreads(threads = 1) 
addArchRGenome("hg38")

archp = loadArchRProject (projdir)

### Call peaks with MACS2 by metaGroupName ####
archp$celltype_revised_sample = paste0(archp$celltype_lv1, '_', archp$Sample)

#metaGroupName = 'celltype_lv1'
metaGroupName = 'celltype_lv1'
celltypes = unique (as.character(archp@cellColData[, metaGroupName]))
archp = archp[archp$celltype_lv1 %in% celltypes]
metaGroupName = 'celltype_lv1'
source ('../../git_repo/utils/chromBPnet_call_peaks.R')

# Run no bias chromBPnet model for each NKT cell subtype ####
chromBPdir = paste0(SCRATCH_DIR, '/chromBPnet')
#chromBPdir = 
dir.create (chromBPdir)
repodir = REPO_DIR
grefdir = paste0(REF_DIR, '/chromBPnet')
biasdir = paste0(DATA_DIR, '/NKT_cells/scatac_ArchR/chromBPnet/NKT_cells/bias_model')

#celltypes = 'B_cells'
#celltypes = c('pDCs')
celltypes = c('Mesothelium_P1','Fibroblasts_P1')
for (celltype in celltypes)
	{
	command <- paste ("bsub -J", paste0(celltype,'_CBPmaster'), 
		paste0("-P ", LSF_ACCOUNT, " -q premium -n 1 -W 72:00 -R rusage[mem=64000] -R span[hosts=1] -o"),
		file.path(chromBPdir,paste0('cBP_master_',celltype,'.out')), "-e" ,
		file.path(chromBPdir,paste0('cBP_master_',celltype,'.err')),
		file.path(repodir,'utils','chromBPnet_master.sh'))
	args <- paste (chromBPdir, grefdir, repodir, celltype, biasdir) 
	system (paste0('chmod +x ',file.path(repodir,'utils','chromBPnet_master.sh')), wait=FALSE)
	
	system (paste(command, args))
	}
