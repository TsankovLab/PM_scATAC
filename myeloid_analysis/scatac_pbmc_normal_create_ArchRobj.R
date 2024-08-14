

  #setwd (projdir)
  fragment_paths = paste0 (file.path('/broad','hptmp','bgiotti','Tmp','Greenleaf','PBMC',
    list.files ('/broad/hptmp/bgiotti/Tmp/Greenleaf/PBMC', pattern = '.arrow')))
  sample_names = list.files ('/broad/hptmp/bgiotti/Tmp/Greenleaf/PBMC', pattern = '.arrow')
  sample_names = sapply (sample_names, function(x) unlist(strsplit (x, '\\.'))[[1]])
  meta = read.csv ('/broad/hptmp/bgiotti/Tmp/Greenleaf/PBMC/greenleaf_pbmc_bm.txt')
  meta$barcode_archr = paste0(meta$sample, '#', meta$barcode)
#ArrowFiles_dir = list.files ('/broad/hptmp/bgiotti/Tmp/Greenleaf/PBMC', pattern = '.arrow')
# if (!all (paste0(sample_names,'.arrow') %in% list.files(ArrowFiles_dir)))
#   {
#   ArrowFiles = createArrowFiles (inputFiles = fragment_paths,
#   sampleNames = sample_names,
#   minTSS = 4, #Dont set this too high because you can always increase later
#   minFrags = 1000,
#   maxFrags = Inf,
#   addTileMat = TRUE,
#   addGeneScoreMat = TRUE,
#   force = TRUE,
#   subThreading = T
#   )
#   } else {
  #ArrowFiles = paste0(ArrowFiles_dir, paste0(sample_names,'.arrow'))
  # }

archp = ArchRProject (
  ArrowFiles = fragment_paths, 
  outputDirectory = projdir,
  copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
  
  
archp$celltype = meta$cell_type[match(rownames(archp@cellColData), meta$barcode_archr)]

# Dimensionality reduction and clustering
varfeat = 25000
LSI_method=2
archp = addIterativeLSI (ArchRProj = archp, 
  useMatrix = "TileMatrix", name = "IterativeLSI",
  force=TRUE, LSIMethod=LSI_method,
  varFeatures = varfeat)
archp = addClusters (input = archp, resolution = .8,
  reducedDims = "IterativeLSI", maxClusters = 100,
  force = TRUE)
archp = addUMAP (ArchRProj = archp, 
  reducedDims = "IterativeLSI",
  force = TRUE)
archp = addTSNE (ArchRProj = archp, 
  reducedDims = "IterativeLSI",
  force = TRUE)

umap_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
 name = "Sample", embedding = "UMAP", pal = palette_sample)
umap_p2 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters",
   embedding = "UMAP")
umap_p3 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "celltype",
   embedding = "UMAP")
umap_p4 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "TSSEnrichment",
   embedding = "UMAP")


pdf (file.path('Plots','celltype_umap.pdf'))
print (umap_p1)
print (umap_p2)
print (umap_p3)
print (umap_p4)
dev.off()

