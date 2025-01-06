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
archp$TNK_cells = 'TNK_cells'
### Preprocess files for launching chromBPnet bias model ####

# Export fragments by meta group
metaGroupName = 'TNK_cells'
if (!exists ('fragments')) fragments = unlist(getFragmentsFromProject (archp))

force = F
for (metagroup in unique(as.character(archp@cellColData[,metaGroupName])))
if (!file.exists(paste0('fragments_',metagroup,'.tsv')) | force)
    {  
    #fragments = ReadFragments(fragment_paths[sam], cutSite = FALSE)
    fragments_metagroup = as.data.frame (fragments[fragments$RG %in% rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) == metagroup]], row.names=NULL)
    fragments_metagroup$strand = '.'
    fragments_metagroup = fragments_metagroup[,c(1:3,5)]
    write.table (fragments_metagroup, file.path('chromBPnet',paste0('fragments_',metagroup,'.tsv')), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
    }

# Call peaks with MACS2 following ENCODE specifications 

# pval_thresh=0.01

# smooth_window=150 # default
# shiftsize=$(( -$smooth_window/2 ))
# fragment_file = 
# macs2 callpeak \
#     -t $fragment_file -f BED -n $prefix -g 1.3e+8 -p $pval_thresh \
#    --shift $shiftsize  --extsize $smooth_window --nomodel --keep-dup all --call-summits
# Set parameters
smooth_window <- 150
shiftsize <- -smooth_window / 2
prefix <- metagroup # Replace with your desired prefix
pval_thresh <- 0.01 # Replace with your p-value threshold

for (metagroup in unique(as.character(archp@cellColData[,metaGroupName])))
	{
	fragment_file <- file.path('chromBPnet',paste0('fragments_',metagroup,'.tsv')) # Replace with actual file path
	# Construct the command
	command <- paste(
	  "macs2 callpeak",
	  "-t", shQuote(fragment_file),
	  "-f BED",
	  "-n", shQuote(prefix),
	  "-g 1.3e+8",
	  "-p", pval_thresh,
	  "--shift", shiftsize,
	  "--extsize", smooth_window,
	  "--nomodel --keep-dup all --call-summits",
	  "--outdir", paste0('chromBPnet/','MACS2_',metagroup)
	)
	
	# Run the command
	system(command)
	}

# Read back MACS2 output and cap peaks to no more than 200K ####
NPEAKS=200000 # capping number of peaks called from MACS2
for (metagroup in unique(as.character(archp@cellColData[,metaGroupName])))
	{
	peaks = read.table (file.path ('chromBPnet',paste0('MACS2_',metagroup), paste0(metagroup,'_peaks.narrowPeak')), sep='\t')	
	peaks = head (peaks[order(-peaks[,7]),], NPEAKS)
	write.table (peaks, file.path('chromBPnet',paste0('MACS2_',metagroup), paste0(metagroup,'_peaks_capped.narrowPeak')), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
	}




