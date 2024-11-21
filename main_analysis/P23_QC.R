library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library (EnsDb.Hsapiens.v86)

plan("multicore", workers = 16)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

####### START ANALYSIS #######
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/P23'
dir.create (file.path(projdir,'Plots'), recursive =T)
setwd (projdir)

sample_names = c(
'P23'
)

fragment_paths = c(
'/sc/arion/projects/Tsankov_Normal_Lung/data/meso_polyICLC/23Mesothelioma/cellranger_output/ALTS04_P22_0_v1/fragments.tsv.gz')
names (fragment_paths) = sample_names

metadata_paths = c(
'/sc/arion/projects/Tsankov_Normal_Lung/data/meso_polyICLC/23Mesothelioma/cellranger_output/ALTS04_P22_0_v1/singlecell.csv')
names (metadata_paths) = sample_names

pmat_paths = c(
'/sc/arion/projects/Tsankov_Normal_Lung/data/meso_polyICLC/23Mesothelioma/cellranger_output/ALTS04_P22_0_v1/raw_peak_bc_matrix.h5')
names (pmat_paths) = sample_names

# peak region paths
pbed_paths = c(
'/sc/arion/projects/Tsankov_Normal_Lung/data/meso_polyICLC/23Mesothelioma/cellranger_output/ALTS04_P22_0_v1/peaks.bed')
names (pbed_paths) = sample_names


sgn_l = list()
for (sam in sample_names)
  {
  counts <- Read10X_h5(filename = pmat_paths[sam])
  
  metadata <- read.csv(
    file = metadata_paths[sam],
    header = TRUE,
    row.names = 1
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = fragment_paths[sam],
    min.cells = 10,
    min.features = 200
  )
  
  sgn <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )  
  sgn_l[sam] = sgn  
  }

for (sam in sample_names)
  {
  # extract gene annotations from EnsDb
  annotations <- GetGRangesFromEnsDb (ensdb = EnsDb.Hsapiens.v86)
  
  # change to UCSC style since the data was mapped to hg19
  seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
  genome (annotations) <- "hg38"
  
  # add the gene information to the object
  Annotation (sgn_l[[sam]]) <- annotations
  
  # compute nucleosome signal score per cell
  sgn_l[[sam]] <- NucleosomeSignal (object = sgn_l[[sam]])
  
  # compute TSS enrichment score per cell
  sgn_l[[sam]] <- TSSEnrichment (object = sgn_l[[sam]], fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks
  sgn_l[[sam]]$pct_reads_in_peaks <- sgn_l[[sam]]$peak_region_fragments / sgn_l[[sam]]$passed_filters * 100
  sgn_l[[sam]]$blacklist_ratio <- sgn_l[[sam]]$blacklist_region_fragments / sgn_l[[sam]]$peak_region_fragments
  }

for (sam in sample_names)
  {
  sgn_l[[sam]] <- subset(
    x = sgn_l[[sam]],
    subset = nCount_peaks > 3000 &
      nCount_peaks < 30000 &
      pct_reads_in_peaks > 15 &
      blacklist_ratio < 0.05 &
      nucleosome_signal < 4 &
      TSS.enrichment > 3
  )
  sgn_l[[sam]]
  }

# Merge signac objects ####
# read in peak region sets
peak_regions = lapply (pbed_paths, function(x) makeGRangesFromDataFrame(read.table(
  file = x,
  col.names = c("chr", "start", "end"))
))
#combined.peaks = reduce(c(peak_regions[[1]],peak_regions[[2]],peak_regions[[3]],peak_regions[[4]]))
combined.peaks = peak_regions[[1]]
# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks


# Process merged sample
DefaultAssay (sgn) = 'peaks'
sgn = RunTFIDF(sgn)
sgn = FindTopFeatures(sgn, min.cutoff = 'q0')
sgn = RunSVD(sgn)
sgn <- RunUMAP(sgn, dims = 2:50, reduction = 'lsi')
sgn <- FindNeighbors (object = sgn, reduction = 'lsi', dims = 2:30)
sgn <- FindClusters (object = sgn, verbose = FALSE, algorithm = 3,res = 2)

pdf ('Plots/normal_pleura_dimplot.pdf')
DimPlot (sgn, group.by = 'dataset', pt.size = 0.1)
dev.off()


# Add GeneScore matrix
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  
# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome (annotations) <- "hg38"

# add the gene information to the object
Annotation (sgn) <- annotations
  
  
gene.activities <- GeneActivity (sgn)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
sgn[['RNA']] <- CreateAssayObject(counts = gene.activities)
sgn <- NormalizeData(
object = sgn,
assay = 'RNA',
normalization.method = 'LogNormalize',
scale.factor = median(sgn$nCount_RNA))


saveRDS ('sgn.rds')


DefaultAssay(sgn) <- 'RNA'
genes = c('KRT19', 'KRT5',
  'AXL', 'SOX9', 'SOX6','SNAI2', 
  'CD3D', 'LYZ','PECAM1','COL1A1',
  'VWF','KLRC1','GNLY','SFTPA1',
  'COL1A2','SOX5','EPCAM','CALB2',
  'ITLN1','TEAD1')
fp = FeaturePlot(
  object = sgn,
  features = genes,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
library (patchwork)
pdf (file.path ('Plots','markers_fplots.pdf'),16,16)
wrap_plots (fp)
dev.off()

### RNA integration
ref = readRDS (file.path ('..','main','scrna','srt.rds'))
# Subset distal and proximal lung reference
# ref = NormalizeData (object = ref, normalization.method = "LogNormalize")
# ref = FindVariableFeatures (ref)
# ref = ScaleData (ref)
# ref = RunPCA (ref, npcs = 30, ndims.print = 1:5, nfeat.print = 5, verbose = FALSE)
# ref = RunUMAP (ref, dims = 1:20)
# ref = FindNeighbors (object = ref, reduction = 'pca', dims = 1:15, k.param = 30,
#                             verbose = TRUE, force.recalc = T)

if (!file.exists ('labeltransfer_preditions.csv'))
	{
	DefaultAssay (sgn) = 'RNA'
	transfer.anchors <- FindTransferAnchors(
	reference = ref,
	query = sgn,
	reduction = 'cca'
	)
	predicted.labels <- TransferData(
	anchorset = transfer.anchors,
	refdata = ref$celltype_simplified,
	weight.reduction = sgn[['lsi']],
	dims = 2:30)
	write.csv (predicted.labels, 'labeltransfer_preditions.csv')
	sgn <- AddMetaData(object = sgn, metadata = predicted.labels)
	} else {
	predicted.labels = read.csv ('labeltransfer_preditions.csv')	
	sgn <- AddMetaData(object = sgn, metadata = predicted.labels)
	}
 

# label doublets
sgn$doublets = ifelse (sgn$seurat_clusters %in% c('14','23','16','26','30','25'), 'doublet','singlet')

dm = DimPlot(
object = sgn,
group.by = 'predicted.id',
label = TRUE,
repel = TRUE)

dm2 = DimPlot(
object = sgn,
group.by = 'seurat_clusters',
label = TRUE,
repel = TRUE)

fp = FeaturePlot(
  object = sgn,
  features = 'prediction.score.max',
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)

dm3 = DimPlot(
object = sgn,
group.by = 'doublets',
label = TRUE,
repel = TRUE)

pdf (paste0('Plots/label_transfer_labels_umap.pdf'),15,10)
wrap_plots (dm, dm2, fp,dm3, ncol=2)
dev.off()

sgn = sgn[,sgn$doublets != 'doublet']
write.csv (sgn$predicted.id, 'barcode_annotation.csv')
saveRDS (sgn, 'sgn.rds')

saveRDS (sgn,'signac_normal.rds')


### Read in ArchR to generate arrow ####
fragment_paths = '/sc/arion/projects/Tsankov_Normal_Lung/data/meso_polyICLC/23Mesothelioma/cellranger_output/ALTS04_P22_0_v1/fragments.tsv.gz'
sample_names = 'P23'
set.seed (1234)
addArchRThreads (threads = 8) 
addArchRGenome ("Hg38")

      #setwd (projdir)  
      ArrowFiles = createArrowFiles (inputFiles = fragment_paths,
      sampleNames = sample_names,
      minTSS = 4, #Dont set this too high because you can always increase later
      minFrags = 1000,
      maxFrags = Inf,
      addTileMat = TRUE,
      addGeneScoreMat = TRUE,
      force = TRUE,
      subThreading = T
      )

projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/P23'
  archp = ArchRProject (
    ArrowFiles = ArrowFiles, 
    outputDirectory = projdir,
    copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  )
saveArchRProject (archp)




