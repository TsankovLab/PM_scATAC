library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library (EnsDb.Hsapiens.v86)

plan("multicore", workers = 16)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

####### START ANALYSIS #######
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/'
dir.create (paste0(projdir,'Plots/'), recursive =T)
setwd (projdir)

sample_names = c(
'RPL_280_neg_1',
'RPL_280_neg_2',
'RPL_Epi_1',
'RPL_Epi_2'#,
)

fragment_paths = c(
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz',
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz',
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz',
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz'#,
#"/ahg/regevdata/projects/ICA_Lung/10x/200116/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01396_AlexTsankov/aDcd45n/outs/fragments.tsv.gz"
)
names (fragment_paths) = sample_names

metadata_paths = c(
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/singlecell.csv',
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/singlecell.csv',
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/singlecell.csv',
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/singlecell.csv'#,
#"/ahg/regevdata/projects/ICA_Lung/10x/200116/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01396_AlexTsankov/aDcd45n/outs/fragments.tsv.gz"
)
names (metadata_paths) = sample_names

pmat_paths = c(
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/raw_peak_bc_matrix.h5',
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/raw_peak_bc_matrix.h5',
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/raw_peak_bc_matrix.h5',
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/raw_peak_bc_matrix.h5'#,
#"/ahg/regevdata/projects/ICA_Lung/10x/200116/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01396_AlexTsankov/aDcd45n/outs/fragments.tsv.gz"
)
names (pmat_paths) = sample_names

# peak region paths
pbed_paths = c(
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/peaks.bed',
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/peaks.bed',
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/peaks.bed',
'/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/peaks.bed'#,
#"/ahg/regevdata/projects/ICA_Lung/10x/200116/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01396_AlexTsankov/aDcd45n/outs/fragments.tsv.gz"
)
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
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  
  # change to UCSC style since the data was mapped to hg19
  seqlevelsStyle(annotations) <- "UCSC"
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
combined.peaks = reduce(c(peak_regions[[1]],peak_regions[[2]],peak_regions[[3]],peak_regions[[4]]))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks


# # load metadata
md <- lapply(metadata_paths, function(x) read.table(
  file = x,
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]) # remove the first row

# # perform an initial filtering of low count cells
# md <- lapply (md, function(x) x[x$passed_filters > 500, ])

# create fragment objects
frags <- lapply(sample_names, function(x) CreateFragmentObject(
  path = fragment_paths[[x]],
  cells = colnames(sgn_l[[x]])
))
names (frags) = sample_names

pmat = lapply (sample_names, function(x) FeatureMatrix(
  fragments = frags[[x]],
  features = combined.peaks,
  cells = colnames(sgn_l[[x]])
))
names (pmat) = sample_names

sng_normal_assay <- lapply (sample_names, function(x) CreateChromatinAssay(pmat[[x]], fragments = frags[[x]]))
names (sng_normal_assay) = sample_names
sng_normal <- lapply (sample_names, function(x) CreateSeuratObject(sng_normal_assay[[x]], assay = "ATAC"))
names (sng_normal) = sample_names

sng_normal <- lapply (sample_names, function(x) {sng_normal[[x]]$dataset = x; sng_normal[[x]]})
names(sng_normal) = sample_names

# merge all datasets, adding a cell ID to make sure cell names are unique
sgn <- merge(
  x = sng_normal[[1]],
  y = sng_normal[2:length(sng_normal)],
  add.cell.ids = sample_names
)
sgn [["ATAC"]]

# Pass the metadata back to the merged object to apply cell filtering the same as for tumor samples
md = lapply (sample_names, function(x) {rownames(md[[x]]) = paste0(x, '_', rownames(md[[x]])); md[[x]]})
md = do.call (rbind, md)
sgn = AddMetaData (sgn, md)

# Process merged sample
sgn = RunTFIDF(sgn)
sgn = FindTopFeatures(sgn, min.cutoff = 'q0')
sgn = RunSVD(sgn)
sgn <- RunUMAP(sgn, dims = 2:50, reduction = 'lsi')

pdf ('Plots/normal_pleura_dimplot.pdf')
DimPlot (sgn, group.by = 'dataset', pt.size = 0.1)
dev.off()


# Add GeneScore matrix
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- "UCSC"
genome (annotations) <- "hg38"

# add the gene information to the object
Annotation (sgn) <- annotations
  
  
gene.activities <- GeneActivity(sgn)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
sgn[['RNA']] <- CreateAssayObject(counts = gene.activities)
sgn <- NormalizeData(
object = sgn,
assay = 'RNA',
normalization.method = 'LogNormalize',
scale.factor = median(sgn$nCount_RNA))

### RNA integration
ref_distal = load ('/ahg/regevdata/projects/ICA_Lung/Results/10x_All_191217/AllDno49_2/10x_integratedU.Rda')
ref_distal = gcdata.integrated
ref_distal = UpdateSeuratObject (ref_distal)

# ref_proximal = load ('/ahg/regevdata/projects/ICA_Lung/Results/10x_All_191217/AllP/10x_All_191217v3_integratedU.Rda')
# ref_proximal = gcdata.integrated
# ref_proximal = UpdateSeuratObject (ref_proximal)

ref = ref_distal

# Subset distal and proximal lung reference
ref = NormalizeData (object = ref, normalization.method = "LogNormalize")
ref = FindVariableFeatures (ref)
ref = ScaleData (ref)
ref = RunPCA (ref, npcs = 30, ndims.print = 1:5, nfeat.print = 5, verbose = FALSE)
ref = RunUMAP (ref, dims = 1:20)
ref = FindNeighbors (object = ref, reduction = 'pca', dims = 1:15, k.param = 30,
                            verbose = TRUE, force.recalc = T)

if (!file.exists ('scatac_normal_labeltransfer_preditions.csv'))
	{
	DefaultAssay (sgn) = 'RNA'
	transfer.anchors <- FindTransferAnchors(
	reference = ref,
	query = sgn,
	reduction = 'cca'
	)
	predicted.labels <- TransferData(
	anchorset = transfer.anchors,
	refdata = ref$celltype,
	weight.reduction = sgn[['lsi']],
	dims = 2:30)
	write.csv (predicted.labels, 'scatac_normal_labeltransfer_preditions.csv')
	sgn <- AddMetaData(object = sgn, metadata = predicted.labels)
	} else {
	predicted.labels = read.csv ('scatac_normal_labeltransfer_preditions.csv')	
	sgn <- AddMetaData(object = sgn, metadata = predicted.labels)
	}
 


dm = DimPlot(
object = sgn,
group.by = 'predicted.id',
label = TRUE,
repel = TRUE)

pdf (paste0('Plots/RNA_integration_normal_lung_umap.pdf'),10,10)
dm
dev.off()

saveRDS (sgn,'signac_normal.rds')


