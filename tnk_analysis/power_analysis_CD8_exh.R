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
srt = readRDS (file.path ('..','scrna','srt.rds'))

### Call peaks with MACS2 by metaGroupName ####
#metaGroupName = 'TNK_cells'
metaGroupName = 'celltype2'

subsampling = c(20,50,100,150,200,250,300,350,400,450,500,550,600,1000,1500,2000)
#subsampling = c(100,1100,2100,3100,4100,5100,6100)
subsampled_barcodes_l= list()
for (sub_n in subsampling)
	{
	archp = loadArchRProject (projdir)	
	archp = archp[archp$celltype2 == 'CD8_exhausted']
	subsampled_barcodes = rownames(archp) [sample(1:length(rownames(archp)),sub_n)]
	archp = archp[rownames(archp) %in% subsampled_barcodes]	
	subsampled_barcodes_l[[sub_n]] = subsampled_barcodes
	if (!exists ('fragments')) fragments = unlist(getFragmentsFromProject (archp))
	dir.create ('power_analysis')	

	force = F
	for (metagroup in unique(as.character(archp@cellColData[,metaGroupName])))
	if (!file.exists(file.path('power_analysis',paste0('fragments_',metagroup,'_',sub_n,'.tsv'))) | force)
    	{ 
    	message (paste('export fragment file for:', metagroup)) 
    	#fragments = ReadFragments(fragment_paths[sam], cutSite = FALSE)
    	fragments_metagroup = as.data.frame (fragments[fragments$RG %in% rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) == metagroup]], row.names=NULL)
    	fragments_metagroup$strand = '.'
    	fragments_metagroup = fragments_metagroup[,c(1:3,5)]
    	write.table (fragments_metagroup, file.path('power_analysis',paste0('fragments_',metagroup,'_',sub_n,'.tsv')), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
    	}


smooth_window <- 150
shiftsize <- -smooth_window / 2
pval_thresh <- 0.05 # Replace with your p-value threshold

for (metagroup in unique(as.character(archp@cellColData[,metaGroupName])))
	{
	prefix <- metagroup # Replace with your desired prefix
	message (paste('call peaks with MACS2 for:', metagroup))	
	fragment_file <- file.path('power_analysis',paste0('fragments_',metagroup,'_',sub_n,'.tsv')) # Replace with actual file path
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
	  "--outdir", paste0('power_analysis/','MACS2_',metagroup,'_',sub_n)
	)
	
	# Run the command
	system(command)
	}
}
saveRDS (subsampled_barcodes_l, 'subsampled_barcodes.rds')

# Read back MACS2 output and cap peaks to no more than 200K ####
NPEAKS=200000 # capping number of peaks called from MACS2


subsampling = c(20,50,100,150,200,250,300,350,400,450,500,550,600,1000,1500,2000)
#subsampling = c(100,1100,2100,3100,4100,5100,6100)
res_l = list()
peaks_l = list()
metaGroupName = 'celltype2'

for (sub_n in subsampling)
	{
message ('read in Peaks and format peak matrix')
for (metagroup in 'CD8_exhausted')
	{
	peaks = read.table (file.path ('power_analysis',paste0('MACS2_',metagroup,'_',sub_n), paste0(metagroup,'_peaks.narrowPeak')), sep='\t')	
	peaks_l[[sub_n]] = length(peaks)
	colnames (peaks) = c('chr','start','end')
	peaks = makeGRangesFromDataFrame(peaks)

	# Define the desired width
	fixed_width <- 500

	# Create a list of GRanges with the desired width
	peaks <- resize(peaks, width = fixed_width, fix = "center") # Resize the ranges to center on the given width

	# Generate matrix of fragment counts of hubs x barcodes ####
	archp = loadArchRProject (projdir)
	if(!exists('pmat')) pmat = getMatrixFromProject (archp, 'PeakMatrix')

	archp = archp[rownames(archp) %in% c(subsampled_barcodes_l[[sub_n]], rownames(archp)[archp$celltype2 == 'CD8'])]
  	# fragments = unlist (getFragmentsFromProject (
    # 	ArchRProj = archp))
  hubsCell_mat = matrix (ncol = length(rownames(archp@cellColData)), nrow = length (peaks))
  colnames (hubsCell_mat) = rownames(archp@cellColData)
  #rownames (hubsCell_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (rownames(archp@cellColData)))
  
  
  pmat2 = pmat[,rownames(archp)]
  pmat2 = pmat2[unique(queryHits (findOverlaps (pmat2, peaks)))]
  pmat2 = assay(pmat2)
  # for (cell in rownames(archp@cellColData)) 
  #   {
  #   pb$tick()  
  #   fragments_in_cell = fragments[as.character(fragments$RG) %in% cell]  
  #   fragments_in_cell = countOverlaps (peaks, fragments_in_cell)
  #   hubsCell_mat[,cell] = fragments_in_cell
  #   }
  #all (colnames (pmat2) == rownames(archp@cellColData))  
  pmat2 = t(t(pmat2) * (10^6 / archp$ReadsInTSS)) # scale
  #saveRDS (hubsCell_mat, file.path (hubs_dir,paste0('hubs_cells_mat.rds')))
  # } else {
  # hubsCell_mat = readRDS (file.path (hubs_dir,paste0('hubs_cells_mat.rds')))  
  # }
#pmat2 = as.data.frame (pmat2)

# Compute differential hub accessibility DHA ####
library (presto)
metaGroupName = 'celltype2'

message ('compute DAP')
res = wilcoxauc (pmat2, as.character(archp@cellColData[,metaGroupName]))
res = res[res$group == 'CD8_exhausted' & res$logFC > 0 & res$pval < 0.05,]
res_l [[sub_n]] = res

	}	
}
saveRDS (res_l, 'power_analysis_subsampling.rds')
res_l = readRDS ('power_analysis_subsampling.rds')
res_df = res_l[subsampling]
res_df = data.frame (value = sapply (res_df, nrow), cells = paste0(subsampling,'_cells'))
res_df$cells = factor (res_df$cells , levels = res_df$cells)
res_df$type = 'DAP'

peaks_l = list()
for (sub_n in subsampling)
	{
	peaks = read.table (file.path ('power_analysis',paste0('MACS2_',metagroup,'_',sub_n), paste0(metagroup,'_peaks.narrowPeak')), sep='\t')	
	peaks_l[[sub_n]] = nrow(peaks)	
	}
res_df1 = rbind (res_df, data.frame (value = unlist(peaks_l[subsampling]),cells= res_df$cells,  type = 'n_peaks'))
res_df2 = cbind (res_df, data.frame (value = unlist(peaks_l[subsampling]),cells= res_df$cells,  type = 'n_peaks'))
colnames(res_df2)[1] = 'DAP'
colnames(res_df2)[4] = 'n_peaks'
res_df2 = res_df2[,-2]
res_df2$DAP_ratio = res_df2$DAP / res_df2$n_peaks

# First, compute a scaling factor to align the line with the barplot visually
# Create DAP_ratio in res_df2
res_df2$DAP_ratio <- res_df2$DAP / res_df2$n_peaks

# Compute scaling factor to make the ratio visible alongside DAP bars
scale_factor <- max(res_df1$value) / max(res_df2$DAP_ratio)

pdf(file.path('power_analysis', 'DAP_cells_barplot.pdf'), height = 3, width=5)

ggplot() + 
  # Stacked barplot for DAP
  geom_bar(data = res_df1, aes(x = cells, y = value, fill = type), 
           stat = 'identity') + 
  
  # Scaled line plot for DAP/n_peaks ratio
  geom_line(data = res_df2, aes(x = cells, y = DAP_ratio * scale_factor, group = 1), 
            color = "red", size = 0.5) +
  geom_point(data = res_df2, aes(x = cells, y = DAP_ratio * scale_factor), 
             color = "red", size = 1.5) +
  
  # Dual y-axis with inverse transformation
  scale_y_continuous(
    name = "DAP",
    sec.axis = sec_axis(~ . / scale_factor, name = "DAP / n_peaks")
  ) +
  
  gtheme +  # your custom theme
  labs(x = "Cells", title = "DAP vs detected peaks") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


