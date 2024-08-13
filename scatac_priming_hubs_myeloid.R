conda activate meso_scatac
use UGER
R

set.seed(1234)

packages = c(
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
  'ArchR',
  'BSgenome.Hsapiens.UCSC.hg38',
  'tidyverse',
  'ggrepel',
  'RColorBrewer')
lapply(packages, require, character.only = TRUE)

####### ANALYSIS of Myeloid compartment #######
projdir = '/ahg/regevdata/projects/ICA_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/priming_hubs'
dir.create (file.path (projdir,'Plots'), recursive =T)
setwd (projdir)

#devtools::install_github("immunogenomics/presto") #needed for DAA
source (file.path('..','..','PM_scATAC','useful_functions.R'))
source (file.path('..','..','PM_scATAC','ggplot_aestetics.R'))
source (file.path('..','..','PM_scATAC','scATAC_functions.R'))
source (file.path('..','..','PM_scATAC','palettes.R'))

set.seed (1234)
addArchRThreads (threads = 8) 
addArchRGenome ("Hg38")

archp_tissue = loadArchRProject (file.path('..'))
archp_pbmc_t = loadArchRProject (file.path('..','..','..','..','AS_human_lung_scatac','analysis','pbmc_myeloid'))
archp_pbmc_h = loadArchRProject (file.path ('..','..','pbmc_normal','scatac_ArchR'))
archp_pbmc_h = archp_pbmc_h[!is.na(archp_pbmc_h$celltype)]


# Export bigiwg files ####
metaGroupName = 'celltype_simplified'
exp_bigwig = T
if (exp_bigwig)
  {
  getGroupBW(
    ArchRProj = archp_pbmc_t,
    groupBy = metaGroupName,
    normMethod = "ReadsInTSS",
    tileSize = 100,
    maxCells = 1000,
    ceiling = 4,
    verbose = TRUE,
    threads = getArchRThreads(),
    logFile = createLogFile("getGroupBW")
  )
  }
# # Check for doublets ####
meso_markers = c('C1QA','APOE','IL1B','CD3D')
archp_pbmc_t = addImputeWeights (archp_pbmc_t)
p <- plotEmbedding(
    ArchRProj = archp_pbmc_t,
    colorBy = "GeneScoreMatrix", 
    name = meso_markers, 
    size=1,
    embedding = "UMAP_H",
    pal = palette_expression,
    imputeWeights = getImputeWeights(archp_pbmc_t)
)

# # Run genescore DAG ####
# archp = archp_pbmc_t
# archp$Clusters_H_pbmc_T = archp_pbmc_t$Clusters_H
# metaGroupName = "Clusters_H_pbmc_T"
# force=FALSE
# if(!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force) source (file.path('..','..','PM_scATAC/DAG.R'))

pdf (file.path('Plots','pbmc_tumor_myeloid_markers_fplots.pdf'), width = 18, height = 15)
wrap_plots (p, ncol=3)
dev.off()



# umap_tissue = plotEmbedding (ArchRProj = archp_tissue, colorBy = "cellColData",
#  name = "celltype", embedding = "UMAP")
# umap_pbmc_h = plotEmbedding (ArchRProj = archp_pbmc_h, 
#   colorBy = "cellColData", name = "Sample2",
#    embedding = "UMAP")
# umap_pbmc_t = plotEmbedding (ArchRProj = archp_pbmc_t, 
#   colorBy = "cellColData", name = "Sample2",
#    embedding = "UMAP_H")
  
#   pdf (file.path('Plots','Clusters_umaps.pdf'),5,5)
#   print (umap_tissue)
#   print (umap_pbmc_h)
#   print (umap_pbmc_t)
#   dev.off()


# Import hubs from tissue myeloid analysis ####
metaGroupName = "Clusters_H"
cor_cutoff = 0.3
#max_dist = 12500
max_dist = 12500
min_peaks = 5
dgs = 0
hubs_dir = paste0 ('../hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks)
hubs_obj = readRDS (file.path(hubs_dir,'global_hubs_obj.rds'))  

# Generate matrix of fragment counts of hubs x metagroup per dataset ####
archp_obj = c('archp_pbmc_t','archp_pbmc_h','archp_tissue')
metaGroupNames = c('celltype_simplified','celltype','celltype')
#metaGroupName = 'Clusters_H'

hubs_mat_c = list()
if (!file.exists(file.path ('hubs_mat_combined.rds'))
  {
  for (i in seq(archp_obj))
    {
    archp = get (archp_obj[i])  
    metaGroupName = metaGroupNames[i]  
    fragments = unlist (getFragmentsFromProject (
      ArchRProj = archp))   
  hubsSample_mat = matrix (ncol = length(unique(archp@cellColData[,metaGroupName])), nrow = length(hubs_obj$hubsCollapsed))
  colnames (hubsSample_mat) = unique(archp@cellColData[,metaGroupName])
  rownames (hubsSample_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (unique(archp@cellColData[,metaGroupName])))
  for (sam in unique(archp@cellColData[,metaGroupName]))
    {
    pb$tick()  
    fragments_in_sample = fragments[fragments$RG %in% rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) == sam]]  
    fragments_in_sample_in_hubs = countOverlaps (hubs_obj$hubsCollapsed, fragments_in_sample)
    hubsSample_mat[,sam] = fragments_in_sample_in_hubs
    }
  frags_in_sample = sapply (unique(archp@cellColData[,metaGroupName]), function(x) sum (archp$nFrags[as.character(archp@cellColData[,metaGroupName]) == x]))
  hubsSample_mat = t(t(hubsSample_mat) * (10^6 / frags_in_sample)) # scale
  hubs_mat_c[[i]] = hubsSample_mat
  #saveRDS (hubsSample_mat, file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))
  }
  saveRDS (hubs_mat_c, 'hubs_mat_combined.rds') 
  } else {
  hubs_mat_c = readRDS ('hubs_mat_combined.rds')
  }

hubs_mat_df = as.data.frame (do.call (cbind, hubs_mat_c))
#colnames (hubs_mat_df)[1:23] =paste0('tumor_pbmc_',colnames (hubs_mat_df)[1:23]) 
hubs_mat_df = hubs_mat_df[,!colnames(hubs_mat_df) %in% 
c('Early.Eryth','CD8.CM','Unk',
  'GMP.Neut','GMP','CD4.M','B',
  'Pre.B','NK','Late.Eryth',
  'CLP.1','CMP.LMPP','CD8.EM',
  'CD4.N1','HSC','CD4.N2','Plasma',
  'CD8.N','CLP.2')]
ha = HeatmapAnnotation (size = anno_barplot(width (hubs_obj$hubsCollapsed), gp = gpar(color = "red"), height =  unit(8, "mm")))
hm = Heatmap (
  scale (t(hubs_mat_df)),
#  top_annotation = ha, 
  column_names_gp = gpar(fontsize = 3),
  show_column_dend = F,
  #row_dend_width = unit(5,'mm'),
  row_dend_side = 'left',
  col = rev(palette_hubs_accessibility),
  border=T,
  name = 'Hubs')
pdf (file.path ('Plots',paste0('hubs_combined_heatmap.pdf')), height=8.2, width=150)
hm
dev.off()




### Do it per sample in monocytes to run differential hub expression between healthy and tumor PBMC ####

# Subset first for monocytes only ####
#archp_tissue_mono = archp_tissue[archp_tissue$celltype == 'Monocytes']
archp_pbmc_t_mono = archp_pbmc_t[archp_pbmc_t$celltype_simplified %in% c('Mye_CD14 Mono','Mye_CD16 Mono')]
archp_pbmc_h_mono = archp_pbmc_h[archp_pbmc_h$celltype %in% c('CD14.Mono.1','CD14.Mono.2')]
#archp_tissue_mono$celltype_sample = paste0(archp_tissue_mono$Sample,'_', archp_tissue_mono$celltype)
archp_pbmc_t_mono$celltype_sample = paste0(archp_pbmc_t_mono$Sample,'_', archp_pbmc_t_mono$celltype_simplified)
archp_pbmc_h_mono$celltype_sample = paste0(archp_pbmc_h_mono$Sample,'_', archp_pbmc_h_mono$celltype)

# Filter out low cells samples ####
min_cells = 30
archp_pbmc_t_mono = archp_pbmc_t_mono[archp_pbmc_t_mono$celltype_sample %in% names(table (archp_pbmc_t_mono$celltype_sample)[table (archp_pbmc_t_mono$celltype_sample) > min_cells])]
archp_pbmc_h_mono = archp_pbmc_h_mono[archp_pbmc_h_mono$celltype_sample %in% names(table (archp_pbmc_h_mono$celltype_sample)[table (archp_pbmc_h_mono$celltype_sample) > min_cells])]
archp_obj = c('archp_pbmc_t_mono','archp_pbmc_h_mono')
metaGroupNames = c('celltype_sample','celltype_sample')

hubs_mat_c = list()
if (!file.exists(file.path ('hubs_mat_pbmc_sample_combined.rds'))
  {
  for (i in seq(archp_obj))
    {
    archp = get (archp_obj[i])  
    metaGroupName = metaGroupNames[i]  
    fragments = unlist (getFragmentsFromProject (
      ArchRProj = archp))   
  hubsSample_mat = matrix (ncol = length(unique(archp@cellColData[,metaGroupName])), nrow = length(hubs_obj$hubsCollapsed))
  colnames (hubsSample_mat) = unique(archp@cellColData[,metaGroupName])
  rownames (hubsSample_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (unique(archp@cellColData[,metaGroupName])))
  for (sam in unique(archp@cellColData[,metaGroupName]))
    {
    pb$tick()  
    fragments_in_sample = fragments[fragments$RG %in% rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) == sam]]  
    fragments_in_sample_in_hubs = countOverlaps (hubs_obj$hubsCollapsed, fragments_in_sample)
    hubsSample_mat[,sam] = fragments_in_sample_in_hubs
    }
  frags_in_sample = sapply (unique(archp@cellColData[,metaGroupName]), function(x) sum (archp$nFrags[as.character(archp@cellColData[,metaGroupName]) == x]))
  hubsSample_mat = t(t(hubsSample_mat) * (10^6 / frags_in_sample)) # scale
  hubs_mat_c[[i]] = hubsSample_mat
  #saveRDS (hubsSample_mat, file.path (hubs_dir,paste0('hubs_sample_',metaGroupName,'_mat.rds')))
  }
  saveRDS (hubs_mat_c, 'hubs_mat_pbmc_sample_combined.rds') 
  } else {
  hubs_mat_c = readRDS ('hubs_mat_pbmc_sample_combined.rds')
  }

hubs_mat_df = do.call (cbind, hubs_mat_c)
t_res = sapply (seq(nrow(hubs_mat_df)), function (x) 
  t.test (hubs_mat_c[[1]][x, grep ('CD14', colnames(hubs_mat_c[[1]]))], hubs_mat_c[[2]][x, grep ('CD14', colnames(hubs_mat_c[[2]]))])$p.value)
lfc_res = sapply (seq(nrow(hubs_mat_df)), function (x) 
  log2(mean(hubs_mat_c[[1]][x, grep ('CD14', colnames(hubs_mat_c[[1]]))])+1) - log2(mean(hubs_mat_c[[2]][x, grep ('CD14', colnames(hubs_mat_c[[2]]))])+1))

t_res_df = data.frame (row.names = hubs_obj$hubs_id, pval = t_res, padj = p.adjust (t_res), LFC= lfc_res)
t_res_df = t_res_df[order(t_res_df$padj),]
head (t_res_df)
hubs_mat_c_14 = hubs_mat_c[]
t.test 

logfcThreshold = 1
pvalAdjTrheshold = 0.05
t_res_df$sig = ifelse (abs(t_res_df$LFC) > logfcThreshold & t_res_df$padj < pvalAdjTrheshold, 1,0)
t_res_df$sig = t_res_df$sig * sign (t_res_df$LFC)
t_res_df$sig[t_res_df$sig == 1] = 'tumor PBMC'
t_res_df$sig[t_res_df$sig == -1] = 'healthy PBMC'

# t_res_df$rna_sign = ifelse (abs(p11_dev_rna$logFC) > logfcThreshold & p11_dev_rna$padj < pvalAdjTrheshold, 1,0)
# t_res_df$rna_sign = p11_dev_rna$rna_sign * sign (p11_dev_rna$rna_diff)
# t_res_df$rna_sign[p11_dev_rna$rna_sign == -1] = 'HOX+'
# t_res_df$rna_sign[p11_dev_rna$rna_sign == 1] = 'HOX-'

res_filtered = t_res_df[abs(t_res_df$LFC) > logfcThreshold & t_res_df$padj < pvalAdjTrheshold,]
res_filtered = head (rownames(res_filtered)[order (-abs(res_filtered$LFC))],20)
t_res_df$labels = ''
t_res_df$labels[match (res_filtered, rownames(t_res_df))] = res_filtered
vp = ggplot (t_res_df, aes(x=LFC, y= -log10(padj))) +
    geom_point(shape=21, aes (fill = sig, color = sig), alpha=.5, size = 1) +
    geom_vline(xintercept = logfcThreshold, linetype="dotted", 
                color = "grey20", size=.5) +
    geom_vline(xintercept = -logfcThreshold, linetype="dotted", 
                color = "grey20", size=.5) +
    geom_hline(yintercept = -log10(pvalAdjTrheshold), linetype="dotted", 
                color = "grey20", size=.5) + 
    geom_text_repel (size=2, data = t_res_df, aes(label = labels),segment.size=.2) + 
    ggtitle ('Hubs differential accessibility') +
    #geom_label_repel (size=2,max.overlaps=10000, data = deg2_cl, aes(label = show_genes), color='red') + 
    scale_color_manual (values = c("0"='grey77',"healthy PBMC"='green',"tumor PBMC"='red')) + 
    scale_fill_manual (values = c("0"='grey77',"healthy PBMC"='green',"tumor PBMC"='red')) + 
    gtheme_no_rot

pdf (file.path ('Plots', 'PBMC_hubs_volcano.pdf'),height=3,width=4)
vp
dev.off()









