
conda activate meso_scatac
R

set.seed(1234)

# Load last istance
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/AS_human_lung_scatac/analysis/pbmc_myeloid/'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

# Load utils functions palettes and packages ####
source (file.path('..','..','PM_scATAC','utils','load_packages.R'))
source (file.path('..','..','PM_scATAC','utils','useful_functions.R'))
source (file.path('..','..','PM_scATAC','utils','ggplot_aestetics.R'))
source (file.path('..','..','PM_scATAC','utils','scATAC_functions.R'))
source (file.path('..','..','PM_scATAC','utils','palettes.R'))

# Load functions for hub detection ####
source (file.path('..','PM_scATAC','utils','knnGen.R'))
source (file.path('..','PM_scATAC','utils','addCoax.R'))
source (file.path('..','PM_scATAC','utils','Hubs_finder.R'))
source (file.path('..','PM_scATAC','utils','hubs_track.R'))

# Set # of threads and genome reference ####
addArchRThreads(threads = 8) 
addArchRGenome("hg38")


#### 
meta_pbmc = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/AS_human_lung_scatac/analysis/pbmc/meta_pbmc.csv')

fragmentsDir = '/broad/hptmp/bgiotti/AS_human_lung_scATAC/fragment_files/'
fragment_filenames = list.files (fragmentsDir, pattern = '*.tsv.gz$')
sampleID = sapply (fragment_filenames, function(x) unlist(strsplit(x, '_'))[1])

meta_pbmc = meta_pbmc[meta_pbmc$sample_ID %in% sampleID, ]
meta_pbmc = meta_pbmc[meta_pbmc$tissue == 'PBMC',]

fragment_filenames = fragment_filenames[sapply (meta_pbmc$sample_ID, function(x) grep (paste0('^',x), fragment_filenames))]
fragment_filenames = fragment_filenames[!grepl('.tbi',fragment_filenames)]
fragment_filepaths_pbmc = paste0(fragmentsDir, fragment_filenames)
sample_names_pbmc = paste0(meta_pbmc$patient_ID,'P')


meta = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/AS_human_lung_scatac/sample_annots.csv')
fr_atac = list.files ('/broad/hptmp/bgiotti/AS_human_lung_scATAC/fragment_files')
fr_atac = sapply (fr_atac, function(x) unlist(strsplit (x, '_'))[1])
table (fr_atac %in% meta$sample_ID)
meta = meta[meta$sample_ID %in% fr_atac, ]
meta = meta[meta$library_chemistry != 'V0_multiome',]
meta = meta[meta$tissue %in% c('Normal','Tumor'), ]
patient_atac = paste0(meta$patient_ID, ifelse(meta$tissue == 'Normal', 'N','T'))
fragments_atac = file.path('/broad/hptmp/bgiotti/AS_human_lung_scATAC/fragment_files',paste0(meta$sample_ID,'_fragments.tsv.gz'))
fragment_filepaths= c(
  '/broad/hptmp/bgiotti/AS_human_lung_scATAC/multiome/MIME30_Lu927N_0_v1/atac_fragments.tsv.gz',
  '/broad/hptmp/bgiotti/AS_human_lung_scATAC/multiome/MIME30_Lu941N_0_v1/atac_fragments.tsv.gz',
  '/broad/hptmp/bgiotti/AS_human_lung_scATAC/multiome/MIME30_Lu952N_0_v1/atac_fragments.tsv.gz',
  '/broad/hptmp/bgiotti/AS_human_lung_scATAC/multiome/MIME30_Lu954N_0_v1/atac_fragments.tsv.gz',
  '/broad/hptmp/bgiotti/AS_human_lung_scATAC/multiome/MIME30_Lu979N_0_v1/atac_fragments.tsv.gz',
  '/broad/hptmp/bgiotti/AS_human_lung_scATAC/multiome/MIME30_Lu927T_0_v1/atac_fragments.tsv.gz',
  '/broad/hptmp/bgiotti/AS_human_lung_scATAC/multiome/MIME30_Lu941T_0_v1/atac_fragments.tsv.gz',
  '/broad/hptmp/bgiotti/AS_human_lung_scATAC/multiome/MIME30_Lu952T_0_v1/atac_fragments.tsv.gz',
  '/broad/hptmp/bgiotti/AS_human_lung_scATAC/multiome/MIME30_Lu954T_0_v1/atac_fragments.tsv.gz',
  '/broad/hptmp/bgiotti/AS_human_lung_scATAC/multiome/MIME30_Lu979T_0_v1/atac_fragments.tsv.gz',
  fragments_atac,
  fragment_filepaths_pbmc
  )
sample_names = c(
  'Lu927N',
  'Lu941N',
  'Lu952N',
  'Lu954N',
  'Lu979N',
  'Lu927T',
  'Lu941T',
  'Lu952T',
  'Lu954T',
  'Lu979T',
  patient_atac,
  sample_names_pbmc
  )

write.csv (data.frame (sampleID = sample_names, fragment_path = fragment_filepaths), '/ahg/regevdata/projects/ICA_Lung/Bruno/AS_human_lung_scatac/metadata.csv')

sample_select = !grepl ('T$',sample_names)
sample_names = sample_names[sample_select]
fragment_filepaths = fragment_filepaths[sample_select]


####### START ANALYSIS - take only P and N #######
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/AS_human_lung_scatac/analysis/pbmc/'
system (paste('mkdir -p',projdir))
setwd (projdir)

ArrowFiles = createArrowFiles (inputFiles = fragment_filepaths,
  sampleNames = sample_names,
  filterTSS = 8, #Dont set this too high because you can always increase later
  filterFrags = 3000,
  minFrags = 1000,
  maxFrags = 100000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = FALSE
)

archp = ArchRProject (
  ArrowFiles = ArrowFiles, 
  outputDirectory = projdir,
  copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

archp = saveArchRProject (archp, load= T)

# Sample-level QC plots
p1 = plotGroups(
    ArchRProj = archp, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
)
p2 = plotGroups(
    ArchRProj = archp, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)

p3 = plotFragmentSizes (ArchRProj = archp)
p4 = plotTSSEnrichment (ArchRProj = archp)
sample_size = data.frame (table(archp$Sample))
sample_size = sample_size [match(sample_names,sample_size$Var1),]
sample_size$Var1 = factor (sample_size$Var1, levels=sample_size$Var1)
sample_size_p = ggplot (data=sample_size, 
  aes(x=Var1, y=Freq,
    fill = Var1)) +
  geom_bar(stat="identity") +
  scale_fill_manual (values=paletteDiscrete(sample_size$Var1)) +
  theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1))
plotPDF (p1,p2,p3,p4, sample_size_p, name = "QC-sample-FragSizes-TSSProfile.pdf", 
  ArchRProj = archp, addDOC = FALSE, width = 5, height = 5)

# Scatterplot per sample of Log10 unique Fragments x TSS Enrichment
tss_fr_p = list()
for (sam in unique (archp$Sample))
{
sam_tss_fr = getCellColData (archp, select = c("log10(nFrags)", "TSSEnrichment"))
sam_tss_fr = sam_tss_fr[archp$Sample == sam,]
tss_fr_p[[sam]] = ggPoint (
    x = sam_tss_fr[,1], 
    y = sam_tss_fr[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(sam_tss_fr[,1], probs = 0.99)),
    ylim = c(0, quantile(sam_tss_fr[,2], probs = 0.99))) + 
  geom_hline (yintercept = min(archp$TSSEnrichment), lty = "dashed") + 
  geom_vline (xintercept = log10(min(archp$nFrags)), lty = "dashed") +
  ggtitle (sam)
}
plotPDF (tss_fr_p, name = "QC-sample-FragSizes-TSSProfile_scatter.pdf", 
  ArchRProj = archp, addDOC = FALSE, width = 5, height = 5)

# Further filter 
#archp = archp[archp$TSSEnrichment > 8 & archp$nFrags > 3000]

# Dimensionality reduction and clustering
set.seed (1234)
archp = addIterativeLSI (ArchRProj = archp, 
  useMatrix = "TileMatrix", name = "IterativeLSI",
  force=TRUE)
archp = addClusters(input = archp, 
  reducedDims = "IterativeLSI",
  force = TRUE)
archp = addUMAP(ArchRProj = archp, 
  reducedDims = "IterativeLSI",
  force = TRUE)

umap_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
 name = "Sample", embedding = "UMAP")
umap_p2 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters",
   embedding = "UMAP")

plotPDF (umap_p1,umap_p2,
 name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = archp, addDOC = FALSE, width = 5, height = 5)

archp = saveArchRProject (ArchRProj = archp, load = TRUE, dropCells = TRUE)

### Harmony ###
archp = addHarmony (
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample", force=TRUE
)
archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony", name='UMAP_H',
    force = TRUE)
archp = addClusters (input = archp,
    reducedDims = "Harmony",
    name='Clusters_H',
    force = TRUE)

metaGroupName='Clusters_H'
metaGroupNames = c('TSSEnrichment','nFrags',metaGroupName,'Sample')  
umap_p1 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
name = x, embedding = "UMAP_H"))
  
pdf (paste0(projdir,'Plots/sample_clusters_harmony_umap.pdf'), 15,15)
wrap_plots (umap_p1, ncol=4)
dev.off()

metaGroupName='Clusters'
metaGroupNames = c('TSSEnrichment','nFrags','ReadsInTSS',metaGroupName,'Sample')  
umap_p2 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
name = x, embedding = "UMAP"))
  
pdf (paste0(projdir,'Plots/sample_clusters_umap.pdf'), 15,15)
wrap_plots (umap_p2, ncol=4)
dev.off()

# Get markers for gene score
immune_markers = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/gene_sets/scRNA_immune_markers_humanLUAD_Samarth_Assaf.csv')
immune_markers = immune_markers [immune_markers$group %in% c('Neutrophil','TRMac','IM','DC2','DC1','pDC','mregDC','CD14 mono','NK','Mast cell','Mgk','B/Plasma',' T cell','Treg','MoMac'),]

#markers = c('WT1','ITLN1','COL1A1','PECAM1','LYZ','CD3D','MSLN','KRT18','KRT5','VIM')
#gs = getMatrixFromProject (archp, useMatrix = 'GeneScoreMatrix')
if (!any (ls() == 'gsSE')) gsSE = getMatrixFromProject (archp, useMatrix = 'GeneScoreMatrix')
gsSE = gsSE[, archp$cellNames]

gs_mat = assays (gsSE)[[1]]
rownames (gs_mat) = gsSE@elementMetadata$name
gs_mat = gs_mat[,rownames(archp@cellColData)]
gs_mat = as.data.frame (t(gs_mat[rownames(gs_mat) %in% immune_markers$gene,]))
gs_mat$clusters = archp@cellColData[,metaGroupName]
gs_mat = aggregate (.~clusters, data= gs_mat, FUN = mean)
rownames(gs_mat) = gs_mat[,1]
gs_mat = gs_mat[,-1]
GeneScore_markers_hm = Heatmap (t(gs_mat), 
        #row_labels = rownames(mat_GeneScore),
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        cluster_rows = T,
        cluster_columns=T,#col = pals_heatmap[[1]],
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 6),
        rect_gp = gpar(col = "white", lwd = .5),
        border=TRUE, column_title = 'Immune markers GeneScore', column_title_gp = gpar(fontsize = 16)
        #right_annotation = motif_ha
        )

GeneScore_markers_scaled_hm = Heatmap (t(scale(gs_mat)), 
        #row_labels = rownames(mat_GeneScore),
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        cluster_rows = T,
        cluster_columns=T,#col = pals_heatmap[[1]],
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 6),
        rect_gp = gpar(col = "white", lwd = .5),
        border=TRUE, column_title = paste("scaled"), column_title_gp = gpar(fontsize = 16)
        #right_annotation = motif_ha
        )

metaGroupName='Clusters_H'
umap_p1 = lapply (metaGroupName, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
name = x, embedding = "UMAP_H"))

pdf (paste0('Plots/DAG_AS_markers_heatmap.pdf'), width=15, height=20)
GeneScore_markers_scaled_hm
dev.off()

projdir2= '/ahg/regevdata/projects/ICA_Lung/Bruno/AS_human_lung_scatac/analysis/NT_myeloid3/'
setwd (projdir2)
archp2 = loadArchRProject (projdir2)
myeloids = rownames(archp2)[as.logical(as.character(archp2@cellColData$status == 'normal'))]

archp$myeloid_normal = rownames(archp) %in% myeloids
archp$tissue = ifelse (grepl ('N$', archp$Sample), 'N','P')
archp$myeloid = archp@cellColData[,metaGroupName] %in% c('C24','C22','C23','C16','C25','C13','C14')
umap_p2 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "myeloid_normal",
   embedding = "UMAP_H", highlightCells = myeloids)
umap_p1 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "tissue",
   embedding = "UMAP_H")
umap_p3 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "myeloid",
   embedding = "UMAP_H")

plotPDF (umap_p2,umap_p1, umap_p3,
 name = "myeloid_normal_UMAP.pdf",
        ArchRProj = archp, addDOC = FALSE, width = 5, height = 5)



# Select myeloids cells ####
metaGroupName = 'Clusters_H'
archp = archp[archp$Clusters_H %in% archp@cellColData[,metaGroupName] %in% c('C24','C22','C23','C16','C25','C13','C14')],

# Dimensionality reduction and clustering
set.seed (1234)
archp = addIterativeLSI (ArchRProj = archp, 
  useMatrix = "TileMatrix", name = "IterativeLSI",
  force=TRUE)
archp = addClusters(input = archp, 
  reducedDims = "IterativeLSI",
  force = TRUE)
archp = addUMAP(ArchRProj = archp, 
  reducedDims = "IterativeLSI",
  force = TRUE)

umap_p1 = plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
 name = "Sample", embedding = "UMAP")
umap_p2 = plotEmbedding (ArchRProj = archp, 
  colorBy = "cellColData", name = "Clusters",
   embedding = "UMAP")

plotPDF (umap_p1,umap_p2,
 name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = archp, addDOC = FALSE, width = 5, height = 5)


# Remove further doublet clusters ####
archp = archp[!archp$Clusters_H %in% c('C2','C3','C4')]

set.seed (1234)
archp = addIterativeLSI (ArchRProj = archp, 
  useMatrix = "TileMatrix", name = "IterativeLSI",
  force=TRUE)

# Harmony ####
archp = addHarmony (
    ArchRProj = archp,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample", force=TRUE
)
archp = addUMAP (ArchRProj = archp, 
    reducedDims = "Harmony", name='UMAP_H',
    force = TRUE)
archp = addClusters (input = archp,
    reducedDims = "Harmony", resolution = 0.2,
    name='Clusters_H',
    force = TRUE)

metaGroupName='Clusters_H'
metaGroupNames = c('TSSEnrichment','nFrags',metaGroupName,'Sample')  
umap_p1 = lapply (metaGroupNames, function(x) plotEmbedding (ArchRProj = archp, colorBy = "cellColData",
name = x, embedding = "UMAP_H"))
  
pdf (paste0(projdir, 'Plots/sample_clusters_harmony_umap.pdf'), 15,15)
wrap_plots (umap_p1, ncol=4)
dev.off()


# Check for doublets ####
meso_markers = c('C1QA','APOE','IL1B','CD3D')
archp = addImputeWeights (archp)
p <- plotEmbedding(
    ArchRProj = archp,
    colorBy = "GeneScoreMatrix", 
    name = meso_markers, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp)
)

pdf (file.path(projdir, 'Plots','myeloid_markers_fplots.pdf'), width = 18, height = 15)
wrap_plots (p, ncol=3)
dev.off()





