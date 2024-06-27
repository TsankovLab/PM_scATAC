
cd /ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/
wget -r http://catlas.org/catlas_downloads/humanbrain/cCREs/
wget -r http://catlas.org/catlas_downloads/humantissues/Peaks/

# Read ENCODE collection for peak annotation
H3K4me3_primary_filepath = '/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/ENCODE/H3K4me3_ChIP-seq/primary_cells/'
encode_H3K4me3_primary = readRDS (paste0(H3K4me3_primary_filepath, 'H3K4me3_primary_cells_bed_list.rds'))

# ### chromVAR analysis
archp = addBgdPeaks (archp, force= FALSE)
archp = addPeakAnnotations (ArchRProj = archp, 
     regions = encode_H3K4me3_primary, name = "ENCODE_H3K4me3")

archp = addDeviationsMatrix (
  ArchRProj = archp, 
  peakAnnotation = "ENCODE_H3K4me3",
  force = FALSE
)

if (!any (ls() == 'encSE')) encSE = ArchR::getMatrixFromProject (archp, useMatrix = 'ENCODE_H3K4me3Matrix')
  encSE = encSE[, archp$cellNames]
  encMat = as.matrix (assays (encSE)[[1]])

dev_sample = lapply (unique(archp$Sample2), function(x) {
	tmp = encMat[,archp$cellNames [archp$Sample2 == x]]
	tmp_max = apply (tmp, 2, function(y) which.max(y))
	data.frame (
		celltype = rownames (tmp)[tmp_max],
		barcode = colnames (tmp)
		)
})

prop_sample = lapply (dev_sample, function(x) as.data.frame (proportions (table (x$celltype))))

# ----- This section prepare a dataframe for labels ---- #
# Get the name and the y position of each label
# calculate the ANGLE of the labels
prop_sample = lapply (prop_sample, function(x)
	{
	number_of_bar <- nrow(x)
	x$id = 1:number_of_bar
	angle <- 90 - 360 * (x$id-0.5) / number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
	x$hjust<-ifelse( angle < -90, 1, 0)
	x$angle<-ifelse(angle < -90, angle+180, angle)
	x
	})
 
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
 
# flip angle BY to make them readable
# ----- ------------------------------------------- ---- #
 
p <- lapply (seq_along(prop_sample), function(x) ggplot(prop_sample[[x]], aes(x=as.factor(Var1), y=Freq)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", fill=alpha("blue", 0.3)) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-1,1) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(4,4), "cm")     # This remove unnecessary margin around plot
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar (start = 0) + 
  geom_text(data=prop_sample[[x]], aes(x=id, y=Freq+0.01, label=Var1, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= prop_sample[[x]]$angle, inherit.aes = FALSE ) +
  ggtitle (unique(archp$Sample2)[x]))


pdf (file.path ('Plots', 'encode_deviations_sample_barplot.pdf'))
p
dev.off()
