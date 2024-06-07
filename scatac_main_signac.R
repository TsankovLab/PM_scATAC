library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library (EnsDb.Hsapiens.v86)

plan("multicore", workers = 16)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

####### START ANALYSIS #######
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_signac'
dir.create (file.path(projdir,'Plots'), recursive =T)
setwd (projdir)


metadata_paths = c(
# Tumor  
'/ahg/regevdata/projects/lungCancerBueno/10x/191121/scATAC_Pt_mesothelioma_CD45_neg_cellranger_atac_v1.2/138_ATACseq_CD45_neg_Lung_ATAC/outs/singlecell.csv',
'/ahg/regevdata/projects/lungCancerBueno/10x/200128/scATAC_Pt811_mesothelioma_CD45pos_neg_cellranger_atac_v1.2/161_ATACseq_Pt811_mesothelioma_CD45pos_CD45neg_ATAC/outs/singlecell.csv',
'/ahg/regevdata/projects/lungCancerBueno/10x/200331/10X_Single_Cell_ATAC/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01218_Robby_AlexTsankov/202_ATAC_826CD45pos_826CD45neg/outs/singlecell.csv',  
'/ahg/regevdata/projects/lungCancerBueno/10x/200721/scATAC_Pt846/10X_Single_Cell_RNA/TD01729_AlexTsankov/846-MesoPool-CD45-pos-CD45-neg-nuclei/outs/singlecell.csv',
'/ahg/regevdata/projects/lungCancerBueno/10x/200911/scATAC_Pt848/10X_Single_Cell_ATAC/TD01814_AlexTsankov/848/outs/singlecell.csv',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_10/230510/P10_scATAC/cellranger_output/ALTS03_Zhao6ATAC_0_v1/singlecell.csv',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_11/230714/ZHAO8mesotheliomaATAC/cellranger_output/ALTS04_Zhao8ATAC_0_v1/singlecell.csv',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_12/230718/ZHAO9mesotheliomaATAC/cellranger_output/ALTS04_Zhao9ATAC_0_v1/singlecell.csv',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_13/231018/ZHAO12mesotheliomaATAC/cellranger_output/ALTS04_Zhao12ATAC_0_v1/singlecell.csv',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_14/240109/ZHAO13mesotheliomaATAC/cellranger_output/ALTS04_Zhao13ATAC_0_v1/singlecell.csv'#,
# Normal
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz',
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz',
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz',
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz'#,
#"/ahg/regevdata/projects/ICA_Lung/10x/200116/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01396_AlexTsankov/aDcd45n/outs/fragments.tsv.gz"
)   
names (metadata_paths) = sample_names

fragment_paths = c(
'/ahg/regevdata/projects/lungCancerBueno/10x/191121/scATAC_Pt_mesothelioma_CD45_neg_cellranger_atac_v1.2/138_ATACseq_CD45_neg_Lung_ATAC/outs/fragments.tsv.gz',
'/ahg/regevdata/projects/lungCancerBueno/10x/200128/scATAC_Pt811_mesothelioma_CD45pos_neg_cellranger_atac_v1.2/161_ATACseq_Pt811_mesothelioma_CD45pos_CD45neg_ATAC/outs/fragments.tsv.gz',
'/ahg/regevdata/projects/lungCancerBueno/10x/200331/10X_Single_Cell_ATAC/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01218_Robby_AlexTsankov/202_ATAC_826CD45pos_826CD45neg/outs/fragments.tsv.gz',  
'/ahg/regevdata/projects/lungCancerBueno/10x/200721/scATAC_Pt846/10X_Single_Cell_RNA/TD01729_AlexTsankov/846-MesoPool-CD45-pos-CD45-neg-nuclei/outs/fragments.tsv.gz',
'/ahg/regevdata/projects/lungCancerBueno/10x/200911/scATAC_Pt848/10X_Single_Cell_ATAC/TD01814_AlexTsankov/848/outs/fragments.tsv.gz',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_10/230510/P10_scATAC/cellranger_output/ALTS03_Zhao6ATAC_0_v1/fragments.tsv.gz',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_11/230714/ZHAO8mesotheliomaATAC/cellranger_output/ALTS04_Zhao8ATAC_0_v1/fragments.tsv.gz',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_12/230718/ZHAO9mesotheliomaATAC/cellranger_output/ALTS04_Zhao9ATAC_0_v1/fragments.tsv.gz',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_13/231018/ZHAO12mesotheliomaATAC/cellranger_output/ALTS04_Zhao12ATAC_0_v1/fragments.tsv.gz',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_14/240109/ZHAO13mesotheliomaATAC/cellranger_output/ALTS04_Zhao13ATAC_0_v1/fragments.tsv.gz'#,
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz',
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz',
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz',
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz'#,
#"/ahg/regevdata/projects/ICA_Lung/10x/200116/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01396_AlexTsankov/aDcd45n/outs/fragments.tsv.gz"
)
names (fragment_paths) = sample_names

pbed_paths = c(
# Tumor  
'/ahg/regevdata/projects/lungCancerBueno/10x/191121/scATAC_Pt_mesothelioma_CD45_neg_cellranger_atac_v1.2/138_ATACseq_CD45_neg_Lung_ATAC/outs/peaks.bed',
'/ahg/regevdata/projects/lungCancerBueno/10x/200128/scATAC_Pt811_mesothelioma_CD45pos_neg_cellranger_atac_v1.2/161_ATACseq_Pt811_mesothelioma_CD45pos_CD45neg_ATAC/outs/peaks.bed',
'/ahg/regevdata/projects/lungCancerBueno/10x/200331/10X_Single_Cell_ATAC/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01218_Robby_AlexTsankov/202_ATAC_826CD45pos_826CD45neg/outs/peaks.bed',  
'/ahg/regevdata/projects/lungCancerBueno/10x/200721/scATAC_Pt846/10X_Single_Cell_RNA/TD01729_AlexTsankov/846-MesoPool-CD45-pos-CD45-neg-nuclei/outs/peaks.bed',
'/ahg/regevdata/projects/lungCancerBueno/10x/200911/scATAC_Pt848/10X_Single_Cell_ATAC/TD01814_AlexTsankov/848/outs/peaks.bed',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_10/230510/P10_scATAC/cellranger_output/ALTS03_Zhao6ATAC_0_v1/peaks.bed',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_11/230714/ZHAO8mesotheliomaATAC/cellranger_output/ALTS04_Zhao8ATAC_0_v1/peaks.bed',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_12/230718/ZHAO9mesotheliomaATAC/cellranger_output/ALTS04_Zhao9ATAC_0_v1/peaks.bed',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_13/231018/ZHAO12mesotheliomaATAC/cellranger_output/ALTS04_Zhao12ATAC_0_v1/peaks.bed',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_14/240109/ZHAO13mesotheliomaATAC/cellranger_output/ALTS04_Zhao13ATAC_0_v1/peaks.bed'#,
# Normal
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz',
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz',
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz',
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz'#,
#"/ahg/regevdata/projects/ICA_Lung/10x/200116/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01396_AlexTsankov/aDcd45n/outs/fragments.tsv.gz"
)
names (pbed_paths) = sample_names

# load processed meso sample signac objects for merging ####
sgn_l = readRDS ('../../per_sample_QC_signac/signac_list.rds')

# Merge signac objects ####

# read in peak region sets
file.exists(pbed_paths)
peak_regions = lapply (pbed_paths, function(x) makeGRangesFromDataFrame (read.table(
  file = x,
  col.names = c("chr", "start", "end"))
))
#combined.peaks = GenomicRanges::reduce (x = c(peak_regions[[1]],peak_regions[[2]],peak_regions[[3]],peak_regions[[4]]))
combined.peaks = GenomicRanges::reduce (Reduce('c', peak_regions))

# Filter out bad peaks based on length
peakwidths = width (combined.peaks)
combined.peaks = combined.peaks[peakwidths  < 10000 & peakwidths > 20]

# # load metadata
md <- lapply (metadata_paths, function(x) read.table(
  file = x,
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]) # remove the first row

# # perform an initial filtering of low count cells
# md <- lapply (md, function(x) x[x$passed_filters > 500, ])

# create fragment objects
frags <- lapply(sample_names, function(x) CreateFragmentObject (
  path = fragment_paths[[x]],
  cells = colnames(sgn_l[[x]])
))
names (frags) = sample_names

pmat = lapply (sample_names, function(x) FeatureMatrix (
  fragments = frags[[x]],
  features = combined.peaks,
  cells = colnames(sgn_l[[x]])
))
names (pmat) = sample_names

sng_assay <- lapply (sample_names, function(x) CreateChromatinAssay (pmat[[x]], fragments = frags[[x]]))
names (sng_assay) = sample_names
sng <- lapply (sample_names, function(x) CreateSeuratObject (sng_assay[[x]], assay = "ATAC"))
names (sng) = sample_names

sng <- lapply (sample_names, function(x) {sng[[x]]$dataset = x; sng[[x]]})
names(sng) = sample_names

# merge all datasets, adding a cell ID to make sure cell names are unique
sgn <- merge (
  x = sng[[1]],
  y = sng[2:length(sng)],
  add.cell.ids = sample_names
)
sgn[["ATAC"]]

# # Pass the metadata back to the merged object to apply cell filtering the same as for tumor samples
# md = lapply (sample_names, function(x) {rownames(md[[x]]) = paste0(x, '_', rownames(md[[x]])); md[[x]]})
# md = do.call (rbind, md)
# sgn = AddMetaData (sgn, md)
saveRDS (sgn, 'sgn.rds')


### Remove bad quality cells ####
write.csv (unlist (lapply (sgn_l, function(x) x$celltype)), 'celltype_annotation.csv')
sgn$celltype = unname (unlist (lapply (sgn_l, function(x) x$celltype)))
sgn = sgn [,sgn$celltype != 'bad_quality']

# Process merged sample ####
sgn = RunTFIDF (sgn)
sgn = FindTopFeatures (sgn, min.cutoff = 'q0')
sgn = RunSVD (sgn)
sgn = RunUMAP (sgn, dims = 2:50, reduction = 'lsi')
sgn = FindNeighbors(object = sgn, reduction = 'lsi', dims = 2:30)
sgn = FindClusters (object = sgn, verbose = FALSE, algorithm = 3)

pdf ('Plots/meso_dimplot_bad_quality_removed.pdf')
DimPlot (sgn, group.by = 'dataset', pt.size = 0.1)
DimPlot (sgn, group.by = 'seurat_clusters', pt.size = 0.1)
DimPlot (sgn, group.by = 'celltype', pt.size = 0.1)
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



