conda activate meso_scatac
R

set.seed(1234)

packages = c(
  'ArchR',
  'Signac',
  'Seurat',
  'biovizBase',
  'ggplot2',
  'patchwork',
  'scATACutils',
  'SummarizedExperiment',
  'epiAneufinder',
  'JASPAR2020',
  'TFBSTools',
  'TxDb.Hsapiens.UCSC.hg38.knownGene',
  'EnsDb.Hsapiens.v86',
  'gplots',
  'regioneR',
  'ComplexHeatmap',
  'BSgenome.Hsapiens.UCSC.hg38')
lapply(packages, require, character.only = TRUE)

####### START ANALYSIS #######
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/per_sample_QC_signac'
dir.create (file.path(projdir,'Plots'), recursive =T)
setwd (projdir)

# Load function (JGranja)
source (file.path('..','PM_scATAC','scATAC_functions.R'))


# combine fragments files from different channels of the normal RPL 
#system ('cat /ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz \
#/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz \
#/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz \
#/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz > RPL_normal_fragments.tsv.gz')

# Fix the fragment file of multiome sample by removing the header 
#system ('zcat atac_fragments.tsv.gz | grep -v ^\# | bgzip > atac_fragments_fixed.tzv.gz')
sample_names = c(
# Tumor  
'P1', # p786
'P4', # p811
'P8', # p826
'P3', # p846
'P5', #'p848'
'P10', # p10
'P11', # p11
'P12', # p12
'P13', # p13
'P14'#,# p14
# Normal
# 'RPL_280_neg_1',
# 'RPL_280_neg_2',
# 'RPL_Epi_1',
# 'RPL_Epi_2'#,
#'cf_distal'
)



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
pmat_paths = c(
# Tumor  
'/ahg/regevdata/projects/lungCancerBueno/10x/191121/scATAC_Pt_mesothelioma_CD45_neg_cellranger_atac_v1.2/138_ATACseq_CD45_neg_Lung_ATAC/outs/raw_peak_bc_matrix.h5',
'/ahg/regevdata/projects/lungCancerBueno/10x/200128/scATAC_Pt811_mesothelioma_CD45pos_neg_cellranger_atac_v1.2/161_ATACseq_Pt811_mesothelioma_CD45pos_CD45neg_ATAC/outs/raw_peak_bc_matrix.h5',
'/ahg/regevdata/projects/lungCancerBueno/10x/200331/10X_Single_Cell_ATAC/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01218_Robby_AlexTsankov/202_ATAC_826CD45pos_826CD45neg/outs/raw_peak_bc_matrix.h5',  
'/ahg/regevdata/projects/lungCancerBueno/10x/200721/scATAC_Pt846/10X_Single_Cell_RNA/TD01729_AlexTsankov/846-MesoPool-CD45-pos-CD45-neg-nuclei/outs/raw_peak_bc_matrix.h5',
'/ahg/regevdata/projects/lungCancerBueno/10x/200911/scATAC_Pt848/10X_Single_Cell_ATAC/TD01814_AlexTsankov/848/outs/raw_peak_bc_matrix.h5',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_10/230510/P10_scATAC/cellranger_output/ALTS03_Zhao6ATAC_0_v1/raw_peak_bc_matrix.h5',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_11/230714/ZHAO8mesotheliomaATAC/cellranger_output/ALTS04_Zhao8ATAC_0_v1/raw_peak_bc_matrix.h5',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_12/230718/ZHAO9mesotheliomaATAC/cellranger_output/ALTS04_Zhao9ATAC_0_v1/raw_peak_bc_matrix.h5',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_13/231018/ZHAO12mesotheliomaATAC/cellranger_output/ALTS04_Zhao12ATAC_0_v1/raw_peak_bc_matrix.h5',
'/ahg/regevdata/projects/lungCancerBueno/10x/MPM_polyICLC/patient_14/240109/ZHAO13mesotheliomaATAC/cellranger_output/ALTS04_Zhao13ATAC_0_v1/raw_peak_bc_matrix.h5'#,
# Normal
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_1/RPL_280_neg_1/outs/fragments.tsv.gz',
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_280_neg_2/RPL_280_neg_2/outs/fragments.tsv.gz',
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_1/RPL_Epi_1/outs/fragments.tsv.gz',
# '/ahg/regevdata/projects/ICA_Lung/10x_scatac/cellranger1.2_count_hg38/RPL_Epi_2/RPL_Epi_2/outs/fragments.tsv.gz'#,
#"/ahg/regevdata/projects/ICA_Lung/10x/200116/wangy33.u.hpc.mssm.edu/10X_Single_Cell_RNA/TD01396_AlexTsankov/aDcd45n/outs/fragments.tsv.gz"
)
names (pmat_paths) = sample_names

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
  seqlevelsStyle (annotations) <- "UCSC"
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

qc_ds = list()
for (sam in sample_names)
  {
  qc_ds[[sam]] = DensityScatter(sgn_l[[sam]], x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE) + ggtitle (sam)
  }
pdf (file.path('Plots','density_scatter.pdf'))
qc_ds
dev.off()

qc_tss = list()
for (sam in sample_names)
  {
  sgn_l[[sam]]$high.tss <- ifelse(sgn_l[[sam]]$TSS.enrichment > 3, 'High', 'Low')
  qc_tss[[sam]] = TSSPlot(sgn_l[[sam]], group.by = 'high.tss') + NoLegend() + ggtitle (sam)
  }
pdf (file.path('Plots','TSSPlot.pdf'))
qc_tss
dev.off()


qc_fr = list()  
for (sam in sample_names)
  {
  sgn_l[[sam]]$nucleosome_group <- ifelse(sgn_l[[sam]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  qc_fr[[sam]] = FragmentHistogram(object = sgn_l[[sam]], group.by = 'nucleosome_group') + ggtitle (sam)
  }
pdf (file.path('Plots','FragmentHistograms.pdf'))
qc_fr
dev.off()

qc_vln = list()
for (sam in sample_names)
  {
  qc_vln[[sam]] = VlnPlot(
    object = sgn_l[[sam]],
    features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
    pt.size = 0.1,
    ncol = 5
  ) + ggtitle (sam)
  }

pdf (file.path('Plots','QC_violin.pdf'), width=10, height=5)
qc_vln
dev.off()

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


for (sam in sample_names)
  {
  sgn_l[[sam]] <- RunTFIDF (sgn_l[[sam]])
  sgn_l[[sam]] <- FindTopFeatures (sgn_l[[sam]], min.cutoff = 'q0')
  sgn_l[[sam]] <- RunSVD (sgn_l[[sam]])
  }

dp = list()
for (sam in sample_names)
  {
  sgn_l[[sam]] <- RunUMAP (object = sgn_l[[sam]], reduction = 'lsi', dims = 2:30)
  sgn_l[[sam]] <- FindNeighbors (object = sgn_l[[sam]], reduction = 'lsi', dims = 2:30)
  sgn_l[[sam]] <- FindClusters (object = sgn_l[[sam]], verbose = FALSE, algorithm = 3)
  dp[[sam]] = DimPlot (object = sgn_l[[sam]], label = TRUE) + NoLegend()
  }

pdf (file.path('Plots','per_sample_umaps.pdf'), width=10, height=5)
dp
dev.off()

fp = list()

for (sam in sample_names)
  {
  gene.activities <- GeneActivity (sgn_l[[sam]])
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  sgn_l[[sam]][['RNA']] <- CreateAssayObject (counts = gene.activities)
  sgn_l[[sam]] <- NormalizeData (
  object = sgn_l[[sam]],
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(sgn_l[[sam]]$nCount_RNA))
  }

celltype_markers = read.csv (paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/highlevel_MPM_markers.csv'))
colnames (celltype_markers)[1] = 'gene'
for (sam in sample_names)
  {    
  DefaultAssay (sgn_l[[sam]]) <- 'RNA'
  fp[[sam]] = FeaturePlot (
    object = sgn_l[[sam]],
    features = celltype_markers[[1]],
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol = 3) 
  }

pdf (file.path('Plots','per_sample_fplots.pdf'), width=15, height=30)
fp
dev.off()


# Get fragments from each sample
fragments_l = list()
for (sam in sample_names)
  {
  if (!file.exists(paste0('fragments_',sam,'.tsv')))
    {  
    #fragments = ReadFragments(fragment_paths[sam], cutSite = FALSE)
    fragments = data.table::fread (cmd = paste0("zcat < ", fragment_paths[sam])) %>% 
    data.frame()
    colnames(fragments)[1:4] = c('seqnames','start','end','RG')
    fragments = fragments[fragments$RG %in% colnames (sgn_l[[sam]])]
    write.table (fragments, paste0('fragments_',sam,'.tsv'), sep='\t', row.names=FALSE)
    } else {
    fragments = read.table (paste0('fragments_',sam,'.tsv'), sep='\t', header=T)  
    fragments = GenomicRanges::makeGRangesFromDataFrame (fragments, keep.extra.columns = TRUE)
    fragments_l[[sam]] = fragments
    }
  }

#fragments_l = lapply (sample_names, function(x) fragments_l[[x]][fragments_l[[x]]$RG %in% colnames(sgn_l[[x]])])

# Get Granges of blacklist regions 
blacklist = toGRanges (paste0('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/blacklisted_regions/ENCODE_blacklist/',"hg38-blacklist.v2.bed")) # taken from https://github.com/Boyle-Lab/Blacklist/tree/master/lists
projdir_cnv = 'CNV_analysis'
dir.create (projdir_cnv)

# Get GRanges of bins excluding black list regions
ws = 10e6
ss = 5e6
ws = 10e5
ss = 10e5
if (!file.exists (file.path (projdir_cnv, paste0 ('windows_',ws,'_',ss,'.rds'))))
  {
  windows = makeWindows (genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist,
    windowSize = ws, slidingSize = ss)
  saveRDS (windows, file.path (projdir_cnv, paste0 ('windows_',ws,'_',ss,'.rds')))
  } else {
  windows = readRDS (file.path (projdir_cnv, paste0 ('windows_',ws,'_',ss,'.rds')))
  }

deleted_chr = c('chr4','chr22','chr13')

print_mat = F
force=F
# Loop per sample and run CNV analysis
for (sam in sample_names)
  {
  # Run CNV analysis
  #force=FALSE
  if (!file.exists (file.path(projdir_cnv, paste0('CNV_LFC_GC_',sam,'_ws_',ws,'_ss_',ss,'.rds'))) | force)
    {
    cnaObj = scCNA (windows, fragments_l[[sam]], neighbors = 100, LFC = 1.5, FDR = 0.1, force = FALSE, remove = c("chrX","chrM","chrY"))
    saveRDS(cnaObj, file.path(projdir_cnv, paste0('CNV_LFC_GC_',sam,'_ws_',ws,'_ss_',ss,'.rds')))
    } else {
    message ('cnaObj found! loading...')
    cnaObj = readRDS (file.path(projdir_cnv, paste0('CNV_LFC_GC_',sam,'_ws_',ws,'_ss_',ss,'.rds')))
    }
  # Get bins matching deleted chromosomes to assess ####
    deleted_chr_mean = sapply (deleted_chr, function(x) 
      colMeans (cnaObj@assays@data@listData[['z']][which(as.character(seqnames(cnaObj@rowRanges)) == x),]))
  
  # Attach CNV load and other CNV-related metrics to signac metadata per sample  
  meta.data = deleted_chr_mean
  rownames(meta.data) = colnames(cnaObj)
  # meta.data = data.frame (
  # row.names = colnames (cnaObj),
  # #cnvload_z = rowSums (abs (t(cnaObj@assays@data@listData[['z']]))),
  # #cnvload_lg = log2(rowSums (abs (t(cnaObj@assays@data@listData[['z']])))),
  # cnvload_counts = rowSums(abs(t(cnaObj@assays@data@listData[['counts']]))),
  # cCNV_score = log2(apply (t(cnaObj@assays@data@listData[['z']]),1, function (x) cCNV_score (x))+1),
  

  
  #rownames(meta.data_df2) = sapply (rownames(meta.data_df2), function(x) unlist(strsplit (x,'\\.'))[2])
  sgn_l[[sam]]@meta.data = sgn_l[[sam]]@meta.data[,!colnames(sgn_l[[sam]]@meta.data) %in% colnames(meta.data)]
  sgn_l[[sam]]@meta.data = cbind (sgn_l[[sam]]@meta.data, meta.data[match(colnames(sgn_l[[sam]]), rownames(meta.data)),])

  # Print CNV heatmap
  if (print_mat)
    {
    mat_type = 'z'
    cnaObj_mat = t(cnaObj@assays@data@listData[[mat_type]])
    colnames (cnaObj_mat) = paste0(seqnames(rowRanges(cnaObj)),':', ranges(rowRanges(cnaObj)))
    rownames (cnaObj_mat) = colnames (cnaObj)
    row_ann = as.character(sgn_l[[sam]]@meta.data[,'seurat_clusters']) [match (rownames(cnaObj_mat),colnames(sgn_l[[sam]]))]
    cnaObj_mat = cnaObj_mat[order(row_ann),]
    row_ann = row_ann[order(row_ann)]
    cnaObj_mat[is.na(cnaObj_mat)] = 0
    cnaObj_mat[is.infinite(cnaObj_mat)]= 0
    ha = HeatmapAnnotation (bar=row_ann,
      which='row')
    #col_fun = colorRamp2(c(min(cnaObj_z), 2,max(cnaObj_z)), c("blue", "white", "red"))  
    row_chr = gsub ('\\:.*','',colnames(cnaObj_mat))
    row_chr = factor (row_chr, levels = unique (row_chr))
    png (file.path('Plots', paste0('GL_method_',sam,'_',mat_type,'_gdn_heatmap.png')), width=1500,height=800)
    print (Heatmap (cnaObj_mat, cluster_rows=T,#col=col_fun,
     cluster_columns=F, left_annotation = ha, row_km = 2,
    row_names_gp = gpar(fontsize = 0.01),column_names_gp = gpar(fontsize = 0.1),
    column_split = row_chr,column_gap = unit(2, "mm"),
    cluster_column_slices = F,row_title_gp = gpar(fontsize = 1),column_title_rot=90))
    dev.off()
    }
  }

cnv_cols = c('cnvload_z','cnvload_lg','cnvload_counts','cCNV_score')
cnv_cols = c('chr22','chr13','chr4')
fp = list()
for (sam in sample_names)
  {
  fp[[sam]] = FeaturePlot (sgn_l[[sam]], feature = cnv_cols, combine=FALSE)
  for (ft in 1:length(fp[[sam]])) {fp[[sam]][[ft]] = fp[[sam]][[ft]] + scale_colour_gradientn (colours = viridis::turbo(100),na.value="white")}
  }
pdf (file.path(projdir,'Plots','CNV_score_per_sample_umap.pdf'),width=10,height=5)
lapply (sample_names, function(x) print (wrap_plots (fp[[x]])))
dev.off()



### RNA integration
srt = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/reproduction2/scRNA/srt_tumor.rds')
srt$celltype_simplified_mal = as.character (srt$celltype_simplified)
srt$celltype_simplified_mal[srt$celltype_simplified == 'Malignant'] = paste0(srt$sampleID[srt$celltype_simplified == 'Malignant'], '_', srt$celltype_simplified[srt$celltype_simplified == 'Malignant'])
# Import P14
srt_p14 = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/MPM_naive_p14_analysis/_cellranger_raw_Filter_400_1000_25/no_harmony/srt.rds')
srt_p14$celltype_simplified = srt_p14$celltype
srt = merge (srt, srt_p14)
srt$sampleID[srt$sampleID == 'MPM_naive_p14'] = 'P14'

predicted.labels = list()
for (sam in sample_names)
  {
  srt_sub = srt[, srt$sampleID == sam]
  srt_sub = NormalizeData (object = srt_sub, normalization.method = "LogNormalize", scale.factor = 10000)
  srt_sub = FindVariableFeatures (srt_sub)
  srt_sub = ScaleData (srt_sub, features = VariableFeatures (object=srt_sub))  
  srt_sub = RunPCA (srt_sub, features = VariableFeatures (object = srt_sub), npcs = 30, ndims.print = 1:5, nfeat.print = 5, verbose = FALSE)
  srt_sub = RunUMAP (object = srt_sub, reduction = 'pca', dims = 1:15)
  # Run denovo clustering on non-adjusted reductions
  srt_sub = FindNeighbors (object = srt_sub, reduction = 'pca', dims = 1:15, k.param = 30,
                              verbose = TRUE, force.recalc = T)
 
  DefaultAssay (sgn_l[[sam]]) = 'RNA'  
  transfer.anchors = FindTransferAnchors (
  reference = srt_sub,
  query = sgn_l[[sam]],
  reduction = 'cca'
  )
  predicted.labels = TransferData (
  anchorset = transfer.anchors,
  refdata = srt_sub$celltype_simplified,
  weight.reduction = sgn_l[[sam]][['lsi']],
  dims = 2:30)
  sgn_l[[sam]] = AddMetaData (object = sgn_l[[sam]], metadata = predicted.labels)
  }

dm = list()
for (sam in sample_names)
  {
  dm[[sam]] = DimPlot(
  object = sgn_l[[sam]],
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + ggtitle(sam)
  }

pdf (file.path('Plots','RNA_integration_x_sample_umap.pdf'),10,10)
dm
dev.off()

# Dont separate malignant annotation by sample
sgn_l = lapply (sgn_l, function(x) {x$predicted.id2 = x$predicted.id; x$predicted.id2[grep('Malignant',x$predicted.id2)] = 'Malignant'; x})

dm = list()
for (sam in sample_names)
  {
  dm[[sam]] = DimPlot(
  object = sgn_l[[sam]],
  group.by = 'predicted.id2',
  label = FALSE,
  repel = TRUE) + ggtitle(sam)
  }

pdf (file.path('Plots','RNA_integration_sample_umap2.pdf'),7,7)
dm
dev.off()


# Run EpiAneuFinder to detect malignant clusters
# Remember to set multiple cores for speeding up computations. Otherwise it take > 24h per sample to run!
ws=1e7
for (sam in sample_names)
  {
  dir.create (paste0(paste0('epiAneufinder_results_',sam,'_ws_',ws)))
  epiAneufinder (input=paste0('fragments_',sam,'.tsv'), #Enter path to your fragments.tsv file or the folder containing bam files
              outdir=paste0('epiAneufinder_results_',sam,'_ws_',ws), #Path to the directory where results should be written 
              blacklist='/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/blacklisted_regions/ENCODE_blacklist/hg38-blacklist.v2.bed', #Path to bed file that contains the blacklisted regions of your genome
              windowSize=ws, 
              genome="BSgenome.Hsapiens.UCSC.hg38", #Substitute with relevant BSgenome
              exclude=c('chrX','chrY','chrM'), 
              reuse.existing=FALSE,
              title_karyo=paste0("Karyogram of ",sam," data"), 
              ncores=12,
              minFrags=10000,
              minsizeCNV=0,
              k=4,
              plotKaryo=TRUE)
  }

# Add cnv load to each signac sample objects
for (sam in sample_names)
  {
  eaf_tab = read.table (paste0('epiAneufinder_results_',sam,'/epiAneufinder_results/results_table.tsv'))
  eaf_tab[eaf_tab==1] = -1
  eaf_tab[eaf_tab==0] = 1
  eaf_tab[eaf_tab==-1] = 0
  cnv_load = colSums (abs(eaf_tab[,4:ncol(eaf_tab)]))
  names (cnv_load) = gsub ('cell.','',names(cnv_load))
  names (cnv_load) = gsub ('\\.','-',names(cnv_load))
  sgn_l[[sam]]$cnv_load_eaf = unname (cnv_load[match(colnames(sgn_l[[sam]]), names(cnv_load))])
  sgn_l[[sam]]$cnv_load_eaf_log = log10(unname (cnv_load[match(colnames(sgn_l[[sam]]), names(cnv_load))]))  
  }

# Plot CNV load on UMAP per sample
fp= list()
for (sam in sample_names)
  {
  fp[[sam]] = FeaturePlot (sgn_l[[sam]], feature = c('cnv_load_eaf','cnv_load_eaf_log'), combine=FALSE)
  for (ft in 1:length(fp[[sam]])) {fp[[sam]][[ft]] = fp[[sam]][[ft]] + scale_colour_gradientn (colours = viridis::turbo(100),na.value="white")}
  }
pdf (file.path('Plots','epianeufinder_fplots.pdf'))
fp
dev.off()

# Find DA gene scores per cluster per sample
da_genescore = list()
for (sam in sample_names)
  {
  DefaultAssay(sgn_l[[sam]]) = 'RNA'  
  Idents (sgn_l[[sam]]) = 'seurat_clusters'
  da_genescore[[sam]] = FindAllMarkers(
    object = sgn_l[[sam]],
    test.use = 'wilcox',
    latent.vars = 'nCount_peaks')  
  }

fp = list()
for (sam in sample_names)
  {
  DefaultAssay(sgn_l[[sam]]) = 'RNA'    
  da_gs = da_genescore[[sam]]
  da_genes = do.call (rbind, lapply (split (da_gs, da_gs$cluster), function(x) head (x$gene,10)))
  fp[[sam]] = FeaturePlot (sgn_l[[sam]], feature = da_genes, combine=FALSE)
  for (ft in 1:length(fp[[sam]])) {fp[[sam]][[ft]] = fp[[sam]][[ft]] + scale_colour_gradientn (colours = viridis::plasma(100),na.value="white")}
  png (file.path ('Plots',paste0('DA_gene_score_',sam,'_fplots.png')), 13500, 13500)
  print (wrap_plots (fp[[sam]]))
  dev.off()
  }


# saveRDS (sgn_l, 'signac_list.rds')


# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
pfm_names = data.frame (
  ID = unlist(lapply(pfm, function(x) x@ID)),
  name = unlist(lapply(pfm, function(x) x@name)))

for (sam in sample_names[1:5])
  {
  DefaultAssay(sgn_l[[sam]]) = 'peaks' 

  # Remove peaks in scaffolds that generated error in P10 when computing chromVAR
  seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "UCSC"
  main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
  keep.peaks <- as.logical(seqnames(granges(sgn_l[[sam]])) %in% main.chroms)
  sgn_l[[sam]][["peaks"]] <- subset(sgn_l[[sam]][["peaks"]], features = rownames(sgn_l[[sam]][["peaks"]])[keep.peaks])

  # add motif information
  sgn_l[[sam]] <- AddMotifs(
    object = sgn_l[[sam]],
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm
  )
  sgn_l[[sam]] <- RunChromVAR(
    object = sgn_l[[sam]],
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  }

saveRDS (sgn_l, 'signac_list.rds')


da_chromvar = list()
for (sam in sample_names)
  {
  DefaultAssay(sgn_l[[sam]]) = 'chromvar'    
  Idents (sgn_l[[sam]]) = 'seurat_clusters'
  da_chromvar[[sam]] = FindAllMarkers(
    object = sgn_l[[sam]],
    test.use = 'wilcox'#,
    #latent.vars = 'nCount_peaks'
    )  
  }

pdf (file.path('Plots','DA_chromvar_heatmaps.pdf'),width=5, height=20)
for (sam in sample_names)
  {
    top_tf = unique(unlist(lapply(split (da_chromvar[[sam]], da_chromvar[[sam]]$cluster), function(x) head(x[,c('gene')],10))))
    da_chromvar_wide = do.call(rbind, lapply(split (da_chromvar[[sam]], da_chromvar[[sam]]$cluster), function(x) x[,c('avg_log2FC','cluster','gene')]))
    da_chromvar_wide = reshape(da_chromvar_wide, idvar = "gene", timevar = "cluster", direction = "wide")
    colnames(da_chromvar_wide) = gsub ('avg_log2FC.','',colnames(da_chromvar_wide))
    rownames(da_chromvar_wide) = da_chromvar_wide$gene
    da_chromvar_wide = da_chromvar_wide[top_tf,]
    da_chromvar_wide = da_chromvar_wide[,-1]
    da_chromvar_wide[is.na(da_chromvar_wide)] = 0
    rownames(da_chromvar_wide) = pfm_names$name[match(rownames(da_chromvar_wide), pfm_names$ID)]
    print (Heatmap (da_chromvar_wide, cluster_rows=F,cluster_columns=F, col = palette_chromvar_fun))
  }
dev.off()


# Further refine cell annotation per sample ####
sgn_l = lapply (sgn_l, function(x) {x$celltype = x$predicted.id; x})

#sample_annotation_P1 
sgn_l[['P1']]$celltype[sgn_l[['P1']]$seurat_clusters == c(5)] = 'Malignant'

#sample_annotation_P10 
sgn_l[['P10']]$celltype[sgn_l[['P10']]$seurat_clusters %in% c(1,5)] = 'Malignant'
sgn_l[['P10']]$celltype[sgn_l[['P10']]$seurat_clusters %in% c(9,14,8,17)] = 'bad_quality'


#sample_annotation_P11
sgn_l[['P11']]$celltype[sgn_l[['P11']]$seurat_clusters %in% c(9)] = 'bad_quality'


#sample_annotation_P12
sgn_l[['P12']]$celltype[sgn_l[['P12']]$seurat_clusters %in% c(7,11,8,12,10)] = 'bad_quality'


#sample_annotation_P14
sgn_l[['P14']]$celltype[!sgn_l[['P14']]$seurat_clusters %in% c(0,1,5,6,11) & sgn_l[['P14']]$celltype == 'T_cells' ] = 'bad_quality'
sgn_l[['P14']]$celltype[sgn_l[['P14']]$seurat_clusters == 12] = 'Malignant'

# Plot refined cell annotation per sample
dm = list()
for (sam in sample_names)
  {
  dm[[sam]] = DimPlot(
  object = sgn_l[[sam]],
  group.by = 'celltype',
  label = FALSE,
  repel = TRUE) + ggtitle(sam)
  }

pdf (file.path('Plots','celltype_sample_umap3.pdf'),7,7)
dm
dev.off()

# Export cell type annotation
write.csv (do.call (rbind, lapply (sgn_l, function(x) as.data.frame(x$celltype))), 'cell_annotation.csv')




