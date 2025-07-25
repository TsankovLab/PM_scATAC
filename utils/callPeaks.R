
archp = addGroupCoverages (
  ArchRProj = archp, 
  groupBy = metaGroupName,  
  force = TRUE,
  minCells= 20, # I think this should be set corresponding to the smallest cluster in the group or lower
  maxCells = 500,
  minReplicates = 2,
  sampleRatio = 0.8,
  useLabels = TRUE)

pdf()
archp = addReproduciblePeakSet (
    archp,
    groupBy= metaGroupName,
    peakMethod = 'Macs2',
    reproducibility = as.character(peak_reproducibility),
    maxPeaks = 500000, 
    minCells=20,
    force =TRUE) # I think this should be set corresponding to the smallest cluster in the group or lower
dev.off()
archp = addPeakMatrix (archp)
  
archp = saveArchRProject (archp, load=TRUE)
  
metaGroupNames = c('TSSEnrichment','nFrags','ReadsInTSS','FRIP')  
pdf()
  umap_p12 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
   name = x, embedding = "UMAP"))
dev.off()
pdf (file.path('Plots','qc_umap_after_filtering.pdf'), 15,15)
print (wrap_plots (umap_p12, ncol=5))
dev.off()
