
#### Import Hi-ChiP data from GSE168881 ####
library (HiCExperiment)
library (HiTC)

samples = c('SRR13961069','SRR13961070','SRR13961071','SRR13961072','SRR13961073','SRR13961074')

hic = list()
for (sam in samples)
  {
  matrix_files = paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/HiC_Texhaustion_GSE168881_hg38/hic_results/matrix/',sam,'/raw/10000/',sam,'_10000.matrix') # raw
  #matrix_files = paste0('/sc/arion/scratch/giottb01/HiCPro_matrices/hic_results/matrix/',sam,'/iced/10000/',sam,'_10000_iced.matrix') # iced normalized
  bed_files = paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/HiC_Texhaustion_GSE168881_hg38/hic_results/matrix/',sam,'/raw/10000/',sam,'_10000_abs.bed')
  bed = read.table (bed_files)
  
  hic[[sam]]<-importC(matrix_files,
               bed_files)
  }

E14exp = list()
E14norm.binned = list()

for (sam in samples)
  {
  # Focus on chr2:156312347−156492348
  
  #chr2 156292347−156642348
  E14subset = extractRegion (hic[[sam]]$chr2chr2, c(1,2),
  chr="chr2", from=156292347, to=156642348)
  ## Binning of 5C interaction map
  E14subset.binned <- binningC (E14subset, binsize=10000, method="median", step=3)
  pdf()
  E14exp[[sam]] <- getExpectedCounts (E14subset.binned, method="loess", stdev=TRUE, plot=TRUE)
  dev.off()
  E14norm.binned[[sam]] <- normPerExpected (E14subset.binned, method="loess", stdev=TRUE)
  }

projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/HiChIP'
dir.create (file.path(projdir,'Plots'))

#E14norm.binned = E14exp

pdf (file.path (projdir,'Plots','HiC_map_iced_normalized2.pdf'),height=3,width=4)
lapply (E14norm.binned, function(x) mapC(x,
#mapC(E14norm.binned[[samples[1]]] - E14norm.binned[[samples[6]]],
#tracks=list(RefSeqGene=gene, CTCF=ctcf),
maxrange=10))
dev.off()



samples2 = c('SRR13961072','SRR13961073')
# Focus on chr2:156312347−156492348

#chr2 156292347−156642348
E14subset = extractRegion(hic[[samples2[1]]]$chr2chr2, c(1,2),
chr="chr2", from=156292347, to=156642348)

E14subset2 = extractRegion(hic[[samples2[2]]]$chr2chr2, c(1,2),
chr="chr2", from=156292347, to=156642348)

## Binning of 5C interaction map
E14subset.binned <- binningC (E14subset, binsize=10000, method="median", step=3)
E14subset2.binned <- binningC (E14subset2, binsize=10000, method="median", step=3)

E14norm.binned <- normPerExpected (E14subset.binned, method="loess", stdev=TRUE)
E14norm2.binned <- normPerExpected (E14subset2.binned, method="loess", stdev=TRUE)

E14norm3.binned = E14norm.binned
E14norm.binned@intdata[is.na(E14norm.binned@intdata)] = 0
E14norm2.binned@intdata[is.na(E14norm2.binned@intdata)] = 0
E14norm3.binned@intdata = E14norm2.binned@intdata - E14norm.binned@intdata

pdf (file.path (projdir,'Plots','HiC_map_diff.pdf'),height=3,width=4)
mapC (E14norm3.binned, maxrange=10)
#mapC(E14norm.binned[[samples[1]]] - E14norm.binned[[samples[6]]],
#tracks=list(RefSeqGene=gene, CTCF=ctcf),

dev.off()


