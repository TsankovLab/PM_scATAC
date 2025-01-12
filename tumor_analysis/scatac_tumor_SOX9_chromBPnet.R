conda activate meso_scatac
R

set.seed(1234)

####### ANALYSIS of tumor compartment #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR'
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


# Assign high vs low sc-S cells per sample ####
metaGroupName = 'sarc_score'
mods = as.data.frame (t(scale(t(archp@cellColData[,grep('cNMF',colnames(archp@cellColData))]))))
mods$Sample = archp$Sample2
mods = mods[order(-mods$cNMF20),]

sams = unique (archp$Sample2)
table (archp$Sample2)
sams = c('P10','P12','P23','P4','P5')
sarc_score = lapply (sams, function(sam) 
	{
	x = mods[mods$Sample == sam,'cNMF20', drop=F]
	barcodes_high = data.frame(
		barcode = rownames(x)[x[[1]] > quantile(x[[1]],0.8)],
		score='high')
	barcodes_low = data.frame(
		barcode= rownames(x)[x[[1]] < quantile(x[[1]],0.2)],
		score = 'low')
	rbind (barcodes_high, barcodes_low)
	})
sarc_score_df = do.call (rbind, sarc_score)
archp$sarc_score = 'mid'
archp$sarc_score[match(sarc_score_df$barcode, rownames(archp@cellColData))] = sarc_score_df$score

pdf (file.path ('Plots','sarc_score_groups.pdf'))
 umap_p1 = plotEmbedding (
 	ArchRProj = archp, 
 	colorBy = "cellColData",
   name = 'sarc_score', embedding = "UMAP")
umap_p1 
dev.off()



### Call peaks with MACS2 by metaGroupName ####
archp$sarc_score_sample = paste0(archp$sarc_score,'_', archp$Sample2)
archp2 = archp
archp = archp[archp$Sample2 %in% sams]
archp = archp[! grepl ('mid',archp$sarc_score_sample)]
metaGroupName = 'sarc_score_sample'
source ('../../git_repo/utils/chromBPnet_call_peaks.R')


# Run no bias chromBPnet model for each NKT cell subtype ####
metaGroupName = 'sarc_score_sample'

chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet'
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







