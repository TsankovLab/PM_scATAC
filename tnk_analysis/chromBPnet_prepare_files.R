# Export bed files of peaks per cell types and fragment files #### 
dir.create ('chromBPnet')

# Export fragments subset by meta group ####
fragments_l = list()
metaGroupName = 'celltype2'

if (!exists ('fragments')) fragments = unlist(getFragmentsFromProject (archp))
if (!file.exists (file.path('chromBPnet',paste0('fragments_NKT_cells.tsv')))) write.table (fragments, file.path('chromBPnet',paste0('fragments_NKT_cells.tsv')), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

for (metagroup in as.character(archp@cellColData[,metaGroupName]))
  {
  if (!file.exists(paste0('fragments_',metagroup,'.tsv')))
    {  
    #fragments = ReadFragments(fragment_paths[sam], cutSite = FALSE)
    fragments_metagroup = fragments[fragments$RG %in% rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) == metagroup]]
    write.table (fragments_metagroup, file.path('chromBPnet',paste0('fragments_',metagroup,'.tsv')), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
    }
  }

# Export peak sets by meta group ####
peaksets = lapply (unique(as.character(archp@cellColData[,metaGroupName])),
  function(x)
  {
   tmp = readRDS (file.path('PeakCalls',paste0(x, '-reproduciblePeaks.gr.rds')))
   tmp = extendGR(gr = tmp, upstream = 1500 - 250, downstream = 1500 - 250)

   tmp = data.frame (tmp)
   tmp = tmp[,1:10]
   tmp[[10]] = 1500
   tmp[,4:9] = '.'
   write.table (tmp, row.names =FALSE, col.names=FALSE, quote=FALSE, sep='\t',file.path('chromBPnet',paste0('peakset_',x,'.bed')))
  })

# also export merged peakset
ps = getPeakSet(archp)
ps = data.frame (extendGR(gr = ps, upstream = 1500 - 250, downstream = 1500 - 250))
ps = ps[head(order (-ps$score),250000),]
ps = ps[,1:10]
ps[[10]] = 1500
ps[,4:9] = '.'
write.table (ps, row.names =FALSE, col.names=FALSE, quote=FALSE, sep='\t',file.path('chromBPnet',paste0('peakset_all.bed')))
