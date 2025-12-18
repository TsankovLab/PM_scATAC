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

### Call peaks with MACS2 by metaGroupName ####
archp2 = archp
sams = c('P23')
archp = archp[archp$Sample2 %in% sams]
archp$SOX9_cluster = ''
archp$SOX9_cluster[archp$Clusters2 == 'C4'] = 'SOX9_low_P23_v2'
archp$SOX9_cluster[archp$Clusters3 == 'C19'] = 'SOX9_high_P23_v2'

## Use SOX9 regulon tails ####
ccomp = archp@cellColData[,c('SCENIC_SOX9','Sample3')]
ccomp_P23 = ccomp[ccomp$Sample3 == 'P23',]
alpha <- 0.05  # total tail probability
lower <- quantile(ccomp_P23$SCENIC_SOX9, alpha / 2, na.rm = TRUE)
upper <- quantile(ccomp_P23$SCENIC_SOX9, 1 - alpha / 2, na.rm = TRUE)
ccomp_P23$breaks = 'not_selected'
ccomp_P23$breaks[ccomp_P23$SCENIC_SOX9 <= lower] = 'SOX9_regulon_low'
ccomp_P23$breaks[ccomp_P23$SCENIC_SOX9 >= upper] = 'SOX9_regulon_high'

archp = archp[match(rownames(ccomp_P23)[ccomp_P23$breaks %in% c('SOX9_regulon_low','SOX9_regulon_high')], rownames(archp@cellColData))]
archp$SOX9_regulon = ccomp_P23$breaks[match (rownames(archp@cellColData), rownames(ccomp_P23))]  
#archp = archp[! grepl ('mid',archp$sarc_score_sample)]
# archp2$sarc_score_sample2 = ''
# archp2$sarc_score_sample2[match(rownames(archp), rownames(archp2))] = archp$sarc_score_sample2

pdf (file.path ('Plots','SOX9_regulon_umap.pdf'))
plotEmbedding (ArchRProj = archp, labelMeans = F, 
  colorBy = "cellColData", name = 'SOX9_regulon',
  #name = "SOX9_cluster", 
  #pal = palette_sample,
   embedding = "UMAP")
dev.off()

# metaGroupName = 'SOX9_cluster'
# archp = archp[archp$SOX9_cluster %in% c('SOX9_low_P23_v2','SOX9_high_P23_v2')]
metaGroupName = 'SOX9_regulon'
source ('../../git_repo/utils/chromBPnet_call_peaks.R')


# Run no bias chromBPnet model for each cell subtype ####
#metaGroupName = 'sarc_score_sample2'

chromBPdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet'
repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
grefdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet'
biasdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet/NKT_cells/bias_model'

#celltypes = unique (as.character(archp@cellColData[, metaGroupName]))
#celltypes = c('SOX9_low_P23_v2','SOX9_high_P23_v2')
celltypes = c('SOX9_regulon_low','SOX9_regulon_high')
#celltypes = 'SOX9_high_P1'
for (celltype in celltypes)
	{
	command <- paste ("bsub -J", paste0(celltype,'_CBPmaster'), 
		"-P acc_Tsankov_Normal_Lung -q premium -n 8 -W 96:00 -R rusage[mem=64000] -R span[hosts=1] -o",
		file.path(chromBPdir,paste0('cBP_master_',celltype,'.out')), "-e" ,
		file.path(chromBPdir,paste0('cBP_master_',celltype,'.err')),
		file.path(repodir,'utils','chromBPnet_master.sh'))
	args <- paste(chromBPdir, grefdir, repodir, celltype, biasdir) 
	system (paste0('chmod +x ',file.path(repodir,'utils','chromBPnet_master.sh')), wait=FALSE)

	system (paste(command, args))
	}






# chromBPdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet
# repodir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo
# grefdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet
# biasdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet/NKT_cells/bias_model
# celltype=SOX9_low_P23
chromBPdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet
repodir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo
grefdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet
biasdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet/NKT_cells/bias_model
