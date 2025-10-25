if (!exists ('fragments')) fragments = unlist(getFragmentsFromProject (archp))
dir.create ('chromBPnet')

force = F
for (metagroup in unique(as.character(archp@cellColData[,metaGroupName])))
if (!file.exists(paste0('fragments_',metagroup,'.tsv')) | force)
    { 
    message (paste('export fragment file for:', metagroup)) 
    #fragments = ReadFragments(fragment_paths[sam], cutSite = FALSE)
    fragments_metagroup = as.data.frame (fragments[fragments$RG %in% rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) == metagroup]], row.names=NULL)
    fragments_metagroup$strand = '.'
    fragments_metagroup = fragments_metagroup[,c(1:3,5)]
    write.table (fragments_metagroup, file.path(chromBPdir,paste0('fragments_',metagroup,'.tsv')), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
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
pval_thresh <- 0.05 # Replace with your p-value threshold

for (metagroup in unique(as.character(archp@cellColData[,metaGroupName])))
	{
	prefix <- metagroup # Replace with your desired prefix
	message (paste('call peaks with MACS2 for:', metagroup))	
	fragment_file <- file.path(chromBPdir,paste0('fragments_',metagroup,'.tsv')) # Replace with actual file path
	# Construct the command
	command <- paste(
	  "macs2 callpeak",
	  "-t", shQuote(fragment_file),
	  "-f BED",
	  "-n", shQuote(prefix),
	  "-g 2.7e+09",
	  "-p", pval_thresh,
	  "--shift", shiftsize,
	  "--extsize", smooth_window,
	  "--nomodel --keep-dup all --call-summits",
	  "--outdir", paste0(chromBPdir,'/','MACS2_',metagroup)
	)
	
	# Run the command
	system(command)
	}

# Read back MACS2 output and cap peaks to no more than 200K ####
NPEAKS=200000 # capping number of peaks called from MACS2
message ('Export peak files')	
for (metagroup in unique(as.character(archp@cellColData[,metaGroupName])))
	{
	peaks = read.table (file.path (chromBPdir,paste0('MACS2_',metagroup), paste0(metagroup,'_peaks.narrowPeak')), sep='\t')	
	peaks = head (peaks[order(-peaks[,7]),], NPEAKS)
	write.table (peaks, file.path(chromBPdir,paste0('MACS2_',metagroup), paste0(metagroup,'_peaks_capped.narrowPeak')), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
	}
message ('done!')	


