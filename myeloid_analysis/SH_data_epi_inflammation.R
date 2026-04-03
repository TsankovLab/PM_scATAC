conda activate meso_scatac
R

# Load packages
packages = c('ArchR','rGREAT','igraph','karyoploteR','Seurat','ggpubr','org.Hs.eg.db',
  'TxDb.Hsapiens.UCSC.hg38.knownGene','EnsDb.Hsapiens.v86','gplots','ggrepel','patchwork','ComplexHeatmap',
  'RColorBrewer','ggrepel','gplots', 'rstatix','dplyr','tidyr')
lapply(packages, require, character.only = TRUE)

#devtools::install_github("immunogenomics/presto") #needed for DAA
#source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/useful_functions.R')
meso_repo = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
source (file.path (meso_repo, 'utils','useful_functions.R'))
source (file.path (meso_repo, 'utils','scATAC_functions.R'))

# Load functions for hub detection
# source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/cooltools/hubs_tools/knnGen.R')
# source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/cooltools/hubs_tools/addCoAx.R')
# source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/cooltools/hubs_tools/Hubs_finder.R')
# source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/cooltools/hubs_tools/hubs_track.R')
# source ('/ahg/regevdata/projects/ICA_Lung/Bruno/scripts/scrna_pipeline/useful_functions.R')
projdir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/NTP_multiome'
setwd (projdir)


set.seed(1234)
addArchRThreads(threads =1) 
addArchRGenome("hg38")

# Load utils functions palettes and packages ####
source (file.path('..','..','..','mesothelioma','scATAC_PM','git_repo','utils','load_packages.R'))
source (file.path('..','..','..','mesothelioma','scATAC_PM','git_repo','utils','useful_functions.R'))
source (file.path('..','..','..','mesothelioma','scATAC_PM','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','..','mesothelioma','scATAC_PM','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','..','mesothelioma','scATAC_PM','git_repo','utils','palettes.R'))


# load objects 
archp = loadArchRProject (projdir)
srt = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/Samarth_multiome_lung_TN_analysis/_cellranger_filtered_Filter_400_800_20/sampleID_harmony_nCount_RNA_regressed/myeloid_subset/sampleID_harmony_nCount_RNA_regressed/cleaned_subset/sampleID_harmony_nCount_RNA_regressed/srt.rds')


#write.table (rownames (archp@cellColData), 'multiome_tumor_cells_scatac.txt')
#write.table (rownames (srt@meta.data), 'multiome_tumor_cells_scrna.txt')
barcode_scatac = read.table ('multiome_tumor_cells_scatac.txt')[[1]]
barcode_scatac2 = gsub ('fragments.tsv.gz#','',barcode_scatac)
barcode_scrna = read.table ('multiome_tumor_cells_scrna.txt')[[1]]
archp = archp[barcode_scatac]
all (rownames(archp@cellColData) == barcode_scatac)
archp@cellColData = cbind (archp@cellColData, srt@meta.data[match(barcode_scatac2,colnames(srt)),])


# Subset objects
archp = archp[rownames(archp@cellColData) %in% barcode_scatac]
srt = srt [, colnames (srt) %in% barcode_scrna]

pdf (file.path ('Plots','celltypes_umap.pdf'))
DimPlot (srt, group.by = 'celltype_simplified')
 plotEmbedding (
    ArchRProj = archp, 
    colorBy = "cellColData", 
    name = 'celltype_simplified', 
    embedding = "UMAP_H",
    #pal = palette_expression,
    imputeWeights = getImputeWeights(archp))
dev.off()

shared_cnmf_genes = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/shared_cnmf_myeloid.rds')
remove_modules = c('Stress','CC','TAM_unknown','TAM_CXCLs')
shared_cnmf_genes = shared_cnmf_genes[!names (shared_cnmf_genes) %in% remove_modules]

#### Run SCENIC ####
scrna_pipeline_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scrna_pipeline'
force = TRUE
org = 'human'
motif_window = 'tss500bp'#'10kbp'
scenic_name = 'monomac_programs_tumor'
genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=5000))
source (file.path(scrna_pipeline_dir, 'SCENIC.R'))

# # Run SCENIC plots ####
# srt$myeloid = 'monomac_programs_tumor'
# motif_window = 'tss500bp'#'10kbp'
# genes.keep = VariableFeatures (FindVariableFeatures (srt, nfeat=5000))
# metaGroupNames = c('sampleID','shared_cnmf2_r_max','myeloid')
# reductionName = 'umap'
# source (file.path(scrna_pipeline_dir, 'SCENIC_plots.R'))



# ### Find cNMF shared between myeloid in each sample ####
# #### Run cNMF ####
# nfeat = 5000
# force=F
# k_list = c(10:20)
# k_selections = c(10:20)
# cores= 100

# cnmf_name = paste0('monomac_programs_tumor')
# #cnmf_out = paste0('cNMF/cNMF_',cnmf_name,'_',paste0(k_list[1],'_',k_list[length(k_list)]),'_vf',nfeat)
# #dir.create (file.path(cnmf_out,'Plots'), recursive=T)
# # repodir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo'
# scrna_pipeline_dir = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/scrna_pipeline'
  
# ### RUN consensus NMF ####
# org='human'
# force=F
# source (file.path (scrna_pipeline_dir,'cnmf_prepare_inputs.R'))     
# srt

# ### Import and format spectra files ####
# k_selection = 10
# cnmf_name = paste0('myeloid')  
# source (file.path (scrna_pipeline_dir,'cnmf_format_spectra_files.R')) 
# #cnmf_spectra_unique = lapply (cnmf_spectra_unique, function(x) head(x,50))


# new_cnmf_names = c(
#   cnmf.1 = 'TREM2', 
#   cnmf.3 = 'Monocytes',
#   cnmf.4 ='SPP1',
#   cnmf.5 = 'IL1B',
#   #cnmf.6 = 'cDCs',
#   cnmf.8 = 'IFN',
#   cnmf.9 = 'IM',
#   cnmf.10 = 'C1Q'
#   )
# names (shared_cnmf_genes) = new_cnmf_names[names (shared_cnmf_genes)]
# shared_cnmf_genes = shared_cnmf_genes[!is.na(names(shared_cnmf_genes))]
# pdf (file.path('Plots','overlap_cnmf_myeloids.pdf'))
# ovmat (c(cnmf_spectra_unique, shared_cnmf_genes), compare_lists = list(names (cnmf_spectra_unique), names(shared_cnmf_genes)), df=F)
# dev.off()

# cnmf_spectra_unique2 = list(
#   IL1B = cnmf_spectra_unique[['cNMF1']][cnmf_spectra_unique[['cNMF1']] %in% shared_cnmf_genes[['IL1B']]],
#   TREM2 = cnmf_spectra_unique[['cNMF8']][cnmf_spectra_unique[['cNMF8']] %in% shared_cnmf_genes[['TREM2']]],
#   Monocytes = cnmf_spectra_unique[['cNMF2']][cnmf_spectra_unique[['cNMF2']] %in% shared_cnmf_genes[['Monocytes']]],
#   SPP1 = cnmf_spectra_unique[['cNMF5']][cnmf_spectra_unique[['cNMF5']] %in% shared_cnmf_genes[['SPP1']]],
#   #C1Q = cnmf_spectra_unique[['cNMF3']][cnmf_spectra_unique[['cNMF3']] %in% shared_cnmf_genes[['C1Q']]],
#   IM = cnmf_spectra_unique[['cNMF3']][cnmf_spectra_unique[['cNMF3']] %in% shared_cnmf_genes[['IM']]]
#   )

# ### Add module scores of harmonized modules ####

# Add modules to srt object and plot them on UMAPs ####
# shared_cnmf_genes = cnmf_spectra_unique2
#shared_cnmf_atac = lapply (shared_cnmf_genes, function(x) x[x %in% getFeatures (archp)])
srt = ModScoreCor (
    seurat_obj = srt, 
    geneset_list = shared_cnmf_genes, 
    cor_threshold = NULL, 
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'shared_cnmf', outdir = NULL)

archp@cellColData = cbind (archp@cellColData, srt@meta.data[match(barcode_scatac2,colnames(srt)),])

#remove_modules = c('cnmf.3','cnmf.6','cnmf.7','cnmf.5') # remove monocyres cDC and CC modules. Consider re-inculding CC 
#shared_cnmf_MAC = shared_cnmf_atac[!names(shared_cnmf_atac) %in% remove_modules]

# force = T
# if (!all (names (shared_cnmf_genes) %in% colnames (archp@cellColData)) | force)
#   {
#   archp@cellColData = archp@cellColData[,!colnames(archp@cellColData) %in% names (shared_cnmf_genes)]
#   archp = addModuleScore (
#       ArchRProj = archp,
#       useMatrix = 'GeneScoreMatrix',
#       name = '',
#       features = shared_cnmf_genes,
#       nBin = 25,
#       nBgd = 100,
#       seed = 1,
#       threads = getArchRThreads(),
#       logFile = createLogFile("addModuleScore")
#     )
#   colnames (archp@cellColData) = gsub ('^\\.','',colnames(archp@cellColData))
#   }


### Plot regulon score of TFs found in km2 along with km2 average score from atac #####

# auc_mtx <- read.csv(file.path('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/NTP_multiome/SCENIC/vg_5000_mw_tss500bp/monomac_programs', 'auc_mtx.csv'), header=T)
# rownames (auc_mtx) = auc_mtx[,1]
# auc_mtx = auc_mtx[,-1]
# colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))
# #auc_mtx = auc_mtx[, colnames(auc_mtx) %in% names(km$cluster[km$cluster == 2])]
# auc_mtx = auc_mtx[, colnames(auc_mtx) %in% AP1]
# auc_mtx_avg = aggregate (auc_mtx, by=as.list(srt@meta.data[,'celltype2',drop=F]), mean)
# rownames (auc_mtx_avg) = auc_mtx_avg[,1]
# auc_mtx_avg = auc_mtx_avg[,-1]
# #auc_mtx_avg = auc_mtx_avg[!rownames(auc_mtx_avg) %in% c('IL1B','cDCs'),]
# auc_mtx_avg_scaled = as.data.frame (scale (auc_mtx_avg))
# auc_mtx_avg_scaled$celltype = rownames(auc_mtx_avg_scaled)
# auc_mtx_avg_scaled_l = gather (auc_mtx_avg_scaled, TF, score,1:(ncol(auc_mtx_avg_scaled)-1))
# #auc_mtx_avg_scaled_l$celltype = factor (auc_mtx_avg_scaled_l$celltype, levels = c('Monocytes','SPP1','TREM2','IFN','C1Q','IM'))


# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = assays (mSE)[[1]]
# mMat = mMat[,rownames(archp@cellColData)]
# rownames (mMat) = rowData (mSE)$name
# mMat = as.matrix(mMat)#[selected_TF,])
# #atac_mod_2 = as.data.frame(t(mMat[names(km$cluster[km$cluster == 2]),]))
# atac_mod_2 = as.data.frame(t(mMat[AP1,]))
# atac_mod_2 = aggregate (atac_mod_2, by=list(celltype = archp$celltype2), mean)
# #atac_mod_2 = atac_mod_2[atac_mod_2$celltype != 'cDCs',]
# atac_mod_2 = atac_mod_2 %>%
#     mutate_if(is.numeric, scale)
# atac_mod_2 = gather (atac_mod_2, TF, score,2:(ncol(atac_mod_2)))

# ### Plot  TF activity
# atac_mod_2_summary <- atac_mod_2 %>%
#   group_by(celltype) %>%
#   dplyr::  summarize(
#     mean_score = mean(score, na.rm = TRUE),
#     sd_score = sd(score, na.rm = TRUE),
#     n = n(),
#     se = sd_score / sqrt(n)  # standard error
#   )
# #atac_mod_2_summary$celltype = factor (atac_mod_2_summary$celltype, levels = c('Monocytes','SPP1','TREM2','C1Q','IFN','IM'))

# gp = ggplot(atac_mod_2_summary, aes(x = celltype, y = mean_score, group = 1)) +
#   geom_line(color = "darkred", size = 1.5) +
#   geom_ribbon(aes(ymin = mean_score - se, ymax = mean_score + se),
#               fill = "darkred", alpha = 0.2) +
#   theme_minimal() +
#   labs(
#     x = "Celltype",
#     y = "Mean Score",
#     title = "Mean Score per Celltype with Standard Error Shading"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ### Plot SCENIC regulon scores

# rna_mod_2_summary <- auc_mtx_avg_scaled_l %>%
#   group_by(celltype) %>%
#   dplyr::  summarize(
#     mean_score = mean(score, na.rm = TRUE),
#     sd_score = sd(score, na.rm = TRUE),
#     n = n(),
#     se = sd_score / sqrt(n)  # standard error
#   )
# #rna_mod_2_summary$celltype = factor (rna_mod_2_summary$celltype, levels = c('Monocytes','SPP1','TREM2','IFN','C1Q','IM'))


# gp = ggplot() +
#   # First dataset (normal scale, left axis)
#   geom_line(data = atac_mod_2_summary, aes(x = celltype, y = mean_score, group=1), color = "darkred", size = .5) +
#   geom_ribbon(data = atac_mod_2_summary, aes(x = celltype, ymin = mean_score - se, ymax = mean_score + se, group=1), 
#               fill = "darkred", alpha = 0.2) +

#   # Second dataset (scaled, right axis)
#   geom_line(data = rna_mod_2_summary, aes(x = celltype, y = mean_score, group=1), color = "navyblue", size = .5) +
#   geom_ribbon(data = rna_mod_2_summary, aes(x = celltype, ymin = mean_score - se, ymax = mean_score + se, group=1),
#               fill = "navyblue", alpha = 0.2) +

#   theme_minimal() +
#   labs(
#     x = "Celltype",
#     title = "Overlayed Mean Score Lines with Separate Y Axes"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


# pdf (file.path ('Plots','inflamed_module_atac_rna_lineplot.pdf'),width=5,height=5)
# gp
# dev.off()






### Plot correlation of regulon score of TFs found in km2 along with correlation of km2 average score from atac #####
km = readRDS ("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/TF_activity_modules.rds")

#AP1 = names (km$cluster[km$cluster == 1])

# AP1 = c('JUNB',
# 'FOSL2',
# 'JUN',
# 'SMARCC1',
# 'FOSL1',
# 'JUND',
# 'FOS',
# 'JDP2',
# 'BACH1',
# 'FOSB',
# 'BATF',
# 'NFE2',
# 'NFE2L2')



# #if (all (names(new_cnmf_names) %in% colnames(srt@meta.data))) colnames(srt@meta.data)[match(names(new_cnmf_names), colnames(srt@meta.data))] = new_cnmf_names
# scenic_name = 'monomac_programs'
# auc_mtx = read.csv(file.path(paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/NTP_multiome/SCENIC/vg_5000_mw_tss500bp/',scenic_name), 'auc_mtx.csv'), header=T)

# rownames (auc_mtx) = auc_mtx[,1]
# auc_mtx = auc_mtx[,-1]
# colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

# regulon_TFs_in_modules = list(
#   km1 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 1])],
#   km2 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 2])]
#   )
# # auc_mtx = auc_mtx[, colnames(auc_mtx) %in% names(km$cluster[km$cluster == 2])]
# auc_mtx = auc_mtx[, colnames(auc_mtx) %in% regulon_TFs_in_modules$km1]
# auc_mtx = auc_mtx[compartment_scrna == compartment_selection,]
# #colnames(srt@meta.data)[colnames(srt@meta.data) == 'Mono'] = 'Monocytes'
# #cnmf_mods = c('IL1B','Monocytes','IFN','C1Q','IM','SPP1','TREM2')
# rownames(auc_mtx) = gsub ('\\.','-',rownames(auc_mtx))
# all (rownames(auc_mtx) == colnames(srt))#2[,compartment == compartment_selection]))

# #cnmf_mods = colnames(srt@meta.data)[new_cnmf_names]
# auc_mtx_cor = as.data.frame (cor (auc_mtx, t(scale(t(srt@meta.data[,names(shared_cnmf_genes)]))), method='spearman'))
# auc_mtx_cor$TF = rownames(auc_mtx_cor)

# auc_mtx_avg_scaled_l = gather (auc_mtx_cor, celltype, score,1:(ncol(auc_mtx_cor)-1))
# #auc_mtx_avg_scaled_l$celltype = factor (auc_mtx_avg_scaled_l$celltype, levels = c('IL1B','Monocytes','SPP1','TREM2','IFN','C1Q','IM'))


# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = scale(assays (mSE)[[1]])
# rownames (mMat) = rowData (mSE)$name
# mMat = mMat[,rownames(archp@cellColData)]
# mMat = as.matrix(mMat)#[selected_TF,])
# #atac_mod_2 = as.data.frame(t(mMat[names(km$cluster[km$cluster == 2]),]))
# atac_mod_2 = as.data.frame(t(mMat[AP1,]))


# #if (all (names(new_cnmf_names) %in% colnames(archp@cellColData))) colnames(archp@cellColData)[match(names(new_cnmf_names), colnames(archp@cellColData))] = new_cnmf_names

# #compartment_selection = 'tumor'
# compartment <- regmatches(rownames(archp@cellColData), regexpr("[TNP]", rownames(archp@cellColData)))
# atac_mod_2 = atac_mod_2
# atac_mod_2 = as.data.frame (cor (atac_mod_2, t(scale(t(as.data.frame(archp@cellColData[,names(shared_cnmf_genes)])))), method='spearman'))
# atac_mod_2$TF = rownames(atac_mod_2)
# atac_mod_2 = gather (atac_mod_2, celltype, score,1:(ncol(atac_mod_2)-1))
# atac_mod_2$celltype = factor (atac_mod_2$celltype, levels = c('IL1B','Monocytes','SPP1','TREM2','IFN','C1Q','IM'))

# ### Plot  TF activity
# atac_mod_2_summary <- atac_mod_2 %>%
#   group_by(celltype) %>%
#   dplyr::  summarize(
#     mean_score = mean(score, na.rm = TRUE),
#     sd_score = sd(score, na.rm = TRUE),
#     n = n(),
#     se = sd_score / sqrt(n)  # standard error
#   )
# atac_mod_2_summary$celltype = factor (atac_mod_2_summary$celltype, levels = c('IL1B','Monocytes','SPP1','TREM2','C1Q','IFN','IM'))

# gp = ggplot(atac_mod_2_summary, aes(x = celltype, y = mean_score, group = 1)) +
#   geom_line(color = "darkred", size = 1.5) +
#   geom_ribbon(aes(ymin = mean_score - se, ymax = mean_score + se),
#               fill = "darkred", alpha = 0.2) +
#   theme_minimal() +
#   labs(
#     x = "Celltype",
#     y = "Mean Score",
#     title = "Mean Score per Celltype with Standard Error Shading"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ### Plot SCENIC regulon scores
# rna_mod_2_summary <- auc_mtx_avg_scaled_l %>%
#   group_by(celltype) %>%
#   dplyr::  summarize(
#     mean_score = mean(score, na.rm = TRUE),
#     sd_score = sd(score, na.rm = TRUE),
#     n = n(),
#     se = sd_score / sqrt(n)  # standard error
#   )
# rna_mod_2_summary$celltype = factor (rna_mod_2_summary$celltype, levels = c('IL1B','Monocytes','SPP1','TREM2','IFN','C1Q','IM'))


# gp = ggplot() +
#   # First dataset (normal scale, left axis)
#   geom_line(data = atac_mod_2_summary, aes(x = celltype, y = mean_score, group=1), color = "darkred", size = .5) +
#   geom_ribbon(data = atac_mod_2_summary, aes(x = celltype, ymin = mean_score - se, ymax = mean_score + se, group=1), 
#               fill = "darkred", alpha = 0.2) +

#   # Second dataset (scaled, right axis)
#   geom_line(data = rna_mod_2_summary, aes(x = celltype, y = mean_score, group=1), color = "navyblue", size = .5) +
#   geom_ribbon(data = rna_mod_2_summary, aes(x = celltype, ymin = mean_score - se, ymax = mean_score + se, group=1),
#               fill = "navyblue", alpha = 0.2) +

#   theme_minimal() +
#   ggtitle (compartment) +
#   labs(
#     x = "Celltype",
#     title = "Overlayed Mean Score Lines with Separate Y Axes"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey44", size = .5) 


# pdf (file.path ('Plots',paste0('inflamed_module_atac_rna_cor_',compartment_selection,'_lineplot3.pdf')),width=5,height=5)
# gp
# dev.off()


# # Show all TFs included in inflammation module ####
# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name
# mMat = as.matrix(mMat)#[selected_TF,])
# atac_mod_2 = as.data.frame(t(mMat[regulon_TFs_in_modules$km2,]))

# new_cnmf_names = c(
#   cnmf.1 = 'TREM2', 
#   cnmf.3 = 'Monocytes',
#   cnmf.4 ='SPP1',
#   cnmf.5 = 'IL1B',
#   #cnmf.6 = 'cDCs',
#   cnmf.8 = 'IFN',
#   cnmf.9 = 'IM',
#   cnmf.10 = 'C1Q'
#   )

# #if (all (names(new_cnmf_names) %in% colnames(archp@cellColData))) colnames(archp@cellColData)[match(names(new_cnmf_names), colnames(archp@cellColData))] = new_cnmf_names

# atac_mod_2 = as.data.frame (cor (atac_mod_2, t(scale(t(as.data.frame(archp@cellColData[,new_cnmf_names]))))))
# #atac_mod_2$TF = rownames(atac_mod_2)
# palette_deviation = paletteer::paletteer_c("ggthemes::Red-Black-White Diverging",100)
# hm = Heatmap (
#     atac_mod_2[,c('IL1B','Monocytes','SPP1','TREM2','IFN','C1Q','IM')],
# #    right_annotation = ha2,
#     column_names_rot =45, 
#     row_names_gp = gpar(fontsize = 5),
#     column_names_gp = gpar(fontsize = 6),
#     col = rev(as.character(palette_deviation)), 
#     cluster_rows=T,
#     cluster_columns = F,
#     border=T
# #rect_gp = gpar (col = "white", lwd = 1)
# )


# ### Check AP1 on the RNA side ####
# scenic_name = 'monomac_programs'
# auc_mtx = read.csv(file.path(paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/NTP_multiome/SCENIC/vg_5000_mw_tss500bp/',scenic_name), 'auc_mtx.csv'), header=T)
# rownames (auc_mtx) = auc_mtx[,1]
# auc_mtx = auc_mtx[,-1]
# colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

# auc_mtx = auc_mtx[, colnames(auc_mtx) %in% regulon_TFs_in_modules$km2]
# auc_mtx = auc_mtx[compartment_scrna == compartment_selection,]
# #colnames(srt@meta.data)[colnames(srt@meta.data) == 'Mono'] = 'Monocytes'
# #cnmf_mods = c('IL1B','Monocytes','IFN','C1Q','IM','SPP1','TREM2')
# rownames(auc_mtx) = gsub ('\\.','-',rownames(auc_mtx))
# all (rownames(auc_mtx) == colnames(srt))#2[,compartment == compartment_selection]))

# cnmf_mods = c('IL1B','Monocytes','IFN','C1Q','IM','SPP1','TREM2')
# auc_mtx_cor = as.data.frame (cor (auc_mtx, srt@meta.data[,cnmf_mods]))
# #auc_mtx_cor$TF = rownames(auc_mtx_cor)

# palette_expression_correlation = paletteer::paletteer_c("ggthemes::Green-Blue-White Diverging",100)
# hm2 = Heatmap (
#     auc_mtx_cor[,c('IL1B','Monocytes','SPP1','TREM2','IFN','C1Q','IM')],
#     #right_annotation = ha2,
#     column_names_rot =45, 
#     row_names_gp = gpar(fontsize = 5),
#     column_names_gp = gpar(fontsize = 6),
#     col = rev(as.character(palette_expression_correlation)), 
#     cluster_rows=T,
#     cluster_columns = F,
#     border=T
# #rect_gp = gpar (col = "white", lwd = 1)
# )


# pdf (file.path ('Plots','inflammation_module_atac_rna_TFs_cor_heatmap.pdf'), width = 3.6,height=3)
# hm + hm2
# dev.off()






### Run TF correlation to identify TF modules across TNK cells ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
#mMat = scale(as.matrix(mMat))#[selected_TF,])

# # Filter by RNA expression ####
metaGroupName = 'celltype_SM'
min_exp = .1
active_TFs = exp_genes (srt, rownames(mMat), min_exp = min_exp, metaGroupName)
mMat = t(scale (mMat[active_TFs, ]))
mMat_cor = cor (as.matrix(mMat), method = 'spearman')

set.seed(123)
centers=2
km = kmeans (mMat_cor, centers=centers)
if (!file.exists ('TF_activity_modules.rds')) saveRDS (km, 'TF_activity_modules.rds')
write.csv (patchvecs (split (names(km$cluster),km$cluster)), 'regMye_modules.csv')

genes_highlight = c(
'IKZF1
HIVEP3
HIVEP1
NFKB2
RELB
RELA
NFKB1
REL
NFE2L2
NFE2
BATF
JDP2
FOSB
SMARCC1
JUND
FOSL1
BACH1
JUNB
FOSL2
FOS
JUN')

genes_highlight2 = c(
'NFE2L2
NFE2
BATF
JDP2
FOSB
SMARCC1
JUND
FOSL1
BACH1
JUNB
FOSL2
FOS
JUN')

genes_highlight2 = c(
'FOSL1
FOSL2
BACH1
PPARG
NFKB2
KLF12
HIVEP3
SMAD1
NFKB1
REL
RUNX1
SNAI1
RUNX2
NFAT5')
# saveRDS (genes_highlight,'inflammation_TFs.rds')

#genes_highlight2 = unlist(strsplit(genes_highlight2,'\n'))
genes_highlight = unlist(strsplit(genes_highlight2,'\n'))

ha2 = rowAnnotation (foo = anno_mark(at = match(genes_highlight,colnames(mMat_cor)), 
    labels = genes_highlight, labels_gp = gpar(fontsize = 7)))

pdf (file.path ('Plots','TF_modules_heatmap3.pdf'), width = 4,height=3)
cor_mMat_hm = draw (Heatmap (mMat_cor,# row_km=15,
  right_annotation = ha2,
  #left_annotation = ha,
  #rect_gp = gpar(type = "none"),
  clustering_distance_rows='euclidean' ,
  clustering_distance_columns = 'euclidean', 
  col=palette_deviation_cor_fun, 
  row_split = km$cluster,
  column_split = km$cluster,
  #row_km=2, 
  #column_km=2,
#  right_annotation = ha,
  border=T,
#   ,
  row_names_gp = gpar(fontsize = 0), 
  column_names_gp = gpar(fontsize = 0)
# cell_fun = function(j, i, x, y, w, h, fill) {# THIS DOESNT WORK NEED TO USE LAYER_FUN
#         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#         }}
  ))
  # ,
  # cell_fun = function(j, i, x, y, w, h, fill) {
  #       if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
  #           grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#        }}))
dev.off()

pdf (file.path ('Plots','TF_modules_heatmap.pdf'), width = 4, height=3)
cor_mMat_hm
dev.off()

# ### Plot UMAP using chromvar TF ####
# mMat = readRDS ('mMat_scaled_active.rds')
# library (uwot)
# set.seed(42)  # for reproducibility

# samples_to_keep = c('Lu952T_fragments.tsv.gz','Lu979T_fragments.tsv.gz')
# mMat_sub = mMat[archp$Sample %in% samples_to_keep & !archp$celltype2 %in% 'Tumor_Mye_CD16 Mono',]
# umap_result <- umap(mMat_sub ,n_components=15, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
# umap_df = as.data.frame(umap_result)
# umap_df$celltypes = archp[archp$Sample %in% samples_to_keep & !archp$celltype2 %in% 'Tumor_Mye_CD16 Mono']$celltype2
# umap_df$FRIP = log2 (archp[archp$Sample %in% samples_to_keep & !archp$celltype2 %in% 'Tumor_Mye_CD16 Mono']$nFrags+1)
# umap_df$Sample = archp[archp$Sample %in% samples_to_keep & !archp$celltype2 %in% 'Tumor_Mye_CD16 Mono']$Sample
# umap_df$mod_2 = archp[archp$Sample %in% samples_to_keep & !archp$celltype2 %in% 'Tumor_Mye_CD16 Mono']$mod_2
# library(ggplot2)

# metaGroupNames1 = c('FRIP','mod_2')
# metaGroupNames2 = c('celltypes','Sample')

# sp = lapply (metaGroupNames1, function(x) 
#     ggplot (umap_df, aes_string ('V1', 'V2', color = x)) +
#     geom_point(size = 2) +
#     scale_color_viridis_b() +
#     labs(title = paste0("UMAP of ",x), x = "UMAP1", y = "UMAP2") +
#     theme_minimal())

# #metaGroupNames = c('celltypes','FRIP','Sample','mod_2')

# sp2 = lapply (metaGroupNames2, function(x) 
#     ggplot (umap_df, aes_string('V1', 'V2', color = x)) +
#     geom_point(size = 2) +
#     #scale_color_viridis_b() +
#     labs(title = paste0("UMAP of ",x), x = "UMAP1", y = "UMAP2") +
#     theme_minimal())

# pdf (file.path ('Plots',paste0('chromvar_mod2_umap.pdf')),width=15,height=15)
# wrap_plots (c(sp, sp2))
# dev.off()

# samples = unique(archp$Sample)
# for (sam in samples)
#   {
#   umap_result <- umap (mMat[as.character(archp@cellColData$Sample) == sam,] , n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
#   library(ggplot2)
  
#   umap_df = as.data.frame(umap_result)
#   umap_df$celltypes = archp$celltype2[as.character(archp@cellColData$Sample) == sam]
#   umap_df$FRIP = log2 (archp$nFrags+1)[as.character(archp@cellColData$Sample) == sam]
#   umap_df$Sample = archp$Sample[as.character(archp@cellColData$Sample) == sam]
#   umap_df$inflm = archp$mod_2[as.character(archp@cellColData$Sample) == sam]
#   sp = lapply (unique(archp$celltype2), function(x) 
#   ggplot (umap_df[umap_df$celltypes == x,], aes(V1, V2, color = celltypes)) +
#     geom_point(size = 2) +
#     labs(title = "UMAP of Iris Data", x = "UMAP1", y = "UMAP2") +
#     theme_minimal())
  
#   sp2 = ggplot (umap_df, aes(V1, V2, color = celltypes)) +
#     geom_point(size = 2) +
#     labs(title = "UMAP of Iris Data", x = "UMAP1", y = "UMAP2") +
#     theme_minimal()
#   sp3 = ggplot (umap_df, aes(V1, V2, color = FRIP)) +
#     geom_point(size = 2) +
#     labs(title = "UMAP of Iris Data", x = "UMAP1", y = "UMAP2") +
#     theme_minimal()
#   sp4 = ggplot (umap_df, aes(V1, V2, color = Sample)) +
#     geom_point(size = 2) +
#     labs(title = "UMAP of Iris Data", x = "UMAP1", y = "UMAP2") +
#     theme_minimal()
#   spq2 = ggplot (umap_df, aes(V1, V2, color = inflm)) +
#   geom_point(size = 2) +
#   scale_color_viridis_b() +
#   labs(title = "UMAP of Iris Data", x = "UMAP1", y = "UMAP2") +
#   theme_minimal()  
  
#   pdf (file.path ('Plots',paste0('chromvar_sample_',sam,'_umap.pdf')),width=15,height=15)
#   print(wrap_plots (sp2, sp4))
#   print(wrap_plots (sp))
#   print (spq2)
#   dev.off()
#   }


# shared_cnmf_genes = readRDS ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scrna/shared_cnmf_myeloid.rds')
# shared_cnmf_genes = lapply (shared_cnmf_genes, function(x) head (x,50))

# # Add modules to srt object and plot them on UMAPs ####
# remove_modules = c('cDCs2','CC_S','LILRA') # remove monocyres cDC and CC modules. Consider re-inculding CC 
# shared_cnmf_genes = shared_cnmf_genes[!names(shared_cnmf_genes) %in% remove_modules]
# srt = ModScoreCor (
#     seurat_obj = srt,
#     geneset_list = shared_cnmf_genes,
#     cor_threshold = NULL,
#     pos_threshold = NULL, # threshold for fetal_pval2
#     listName = 'shared_cnmf', outdir = NULL)
# scatac_barcodes = rownames (archp@cellColData)
# scatac_barcodes = gsub ('fragments.tsv.gz#','', scatac_barcodes)
# archp$shared_cnmf = srt$shared_cnmf_r_max[match(scatac_barcodes, rownames (srt@meta.data))]

# reductionName = 'orig.ident_harmony_umap'
# pdf (file.path ('Plots','cnmf_modules_umaps.pdf'))
# wrap_plots (fp (srt, names(shared_cnmf_genes), reduction=reductionName))
# dev.off()

# ### Check inflammation score across TAMs ####
# palette_myeloid
# mod_df = data.frame (
#   #celltype = archp$cnmf_celltypes,
#   #sample = archp$Sample,
#   #Infl_module = archp$mod_2,
#   Infl_module = archp$mod_2)
# mod_df = aggregate(archp$mod_2, by=list(
#   celltype = archp$celltype2,
#   sample = archp$Sample), FUN=mean)
# head (mod_df)
# df_order = mod_df %>% 
# group_by (celltype) %>% 
# summarize (avg_module = median(x)) %>% 
# arrange(avg_module)
# mod_df$celltype = factor (mod_df$celltype, levels = rev(df_order$celltype))

# bp = ggplot (mod_df, aes (x = celltype, y = x, fill=celltype)) +
# vlp + 
# bxpv + 
# scale_fill_manual (values = palette_myeloid) +
# #geom_point (position='identity', alpha=.3, color="grey44", size=1) +
# gtheme

# pdf (file.path ('Plots','celltype_infl_module_sample_boxplots.pdf'),2.3,width=3)
# bp
# dev.off()





# # Try with ridge plots ####
library (ggridges)
library (ggplot2)
library (viridis)
library (tidyr)
#library(hrbrthemes)


# cnmf_scatac = as.data.frame (scale(scale(t(archp@cellColData[,names(shared_cnmf)]))))

# if (!file.exists ('cnmf_scatac.rds')) saveRDS (cnmf_scatac, 'cnmf_scatac.rds')

# cap = 3
# cnmf_scatac_cap = cnmf_scatac
# cnmf_scatac_cap[cnmf_scatac_cap > cap] = cap
# cnmf_scatac_cap[cnmf_scatac_cap < -cap] = -cap
# set.seed (123)
# km_cnmf = kmeans (t(scale(t(cnmf_scatac_cap))), centers=6) # double scale modules and cluster using k-means
# ha = HeatmapAnnotation (sample = archp$Sample, col=list(sample = palette_sample))
# hm = Heatmap (t(cnmf_scatac_cap), 
#   col = palette_genescore_fun(cnmf_scatac_cap), 
#   top_annotation = ha,
#   clustering_distance_columns = 'pearson',
#   clustering_distance_rows = 'pearson',
#   show_column_dend = F,
#   column_split = km_cnmf$cluster,
#   #column_km=3,
#   row_names_gp = gpar (fontsize = 8),
#   column_names_gp = gpar (fontsize = 0),
#   border=T)

# pdf (file.path ('Plots','cnmf_scatac_scaled_only_MAC2_heatmap.pdf'), height=1.5)
# hm
# dev.off()

# # Show TAM modules in UMAPs ####
# archp = addImputeWeights (archp)
# pdf()
# p <- plotEmbedding (
#     ArchRProj = archp, 
#     colorBy = "cellColData", 
#     name = names (shared_cnmf), 
#     embedding = "UMAP_H",
#     pal = palette_expression,
#     imputeWeights = getImputeWeights(archp)
# )
# # p2 <- plotEmbedding (
# #     ArchRProj = archp, 
# #     colorBy = "GeneScoreMatrix", 
# #     name = c('C3','FOSB','CSF1R','A2M','NAIP','ADAP2','FOLR2','MALAT1' ),
# #     embedding = "UMAP_H",
# #     pal = palette_expression,
# #     imputeWeights = getImputeWeights(archp)
# # )
# dev.off()
# pdf (file.path ('Plots','shared_cnmf_TAMs_fplots.pdf'),14,14)
# wrap_plots (p, ncol=6)
# # p2
# dev.off()


# # Re-Annotate based on cnmf clustering ####
# all (names(km_cnmf$cluster) == rownames(archp@cellColData))
# archp$cnmf_cluster = paste0('cnmf_cluster_',km_cnmf$cluster)
# archp$celltype_lv2 = archp$cnmf_cluster
# archp$celltype_lv2[archp$celltype_lv2 == 'cnmf_cluster_1'] = 'Mono_CD16'
# archp$celltype_lv2[archp$celltype_lv2 == 'cnmf_cluster_2'] = 'TAM_MARCO'
# archp$celltype_lv2[archp$celltype_lv2 == 'cnmf_cluster_3'] = 'TAM_TREM2'
# archp$celltype_lv2[archp$celltype_lv2 == 'cnmf_cluster_4'] = 'Mono_CD14'
# archp$celltype_lv2[archp$celltype_lv2 == 'cnmf_cluster_5'] = 'cDCs'
# archp$celltype_lv2[archp$celltype_lv2 == 'cnmf_cluster_6'] = 'TAM_interstitial'
# #archp$celltype_lv2[archp$celltype_lv2 == 'cnmf_cluster_7'] = 'TAM_CXCLs'
# #archp$celltype_lv2[archp$celltype_lv2 == 'cnmf_cluster_6'] = 'TREM2'


# Plot
ccomp = as.data.frame (archp@cellColData)
#ccomp = ccomp[ccomp$cnmf_celltypes %in% c('cDCs'),]
ccomp$cnmf_celltypes = NA
ccomp$cnmf_celltypes[ccomp$celltype2 == 'Tumor_Mono CD14'] = 'Mono_CD14'
ccomp$cnmf_celltypes[ccomp$celltype2 == 'Tumor_Mono CD16'] = 'Mono_CD16'
ccomp$cnmf_celltypes[ccomp$celltype2 == 'Tumor_Mye_DC2'] = 'cDCs'
ccomp$cnmf_celltypes[ccomp$celltype2 == 'Tumor_Mye_Mac_TREM2'] = 'TREM2'
ccomp$cnmf_celltypes[ccomp$celltype2 == 'Tumor_Mye_Mac_SPP1'] = 'SPP1'
ccomp$cnmf_celltypes[ccomp$celltype2 == 'Tumor_Mye_Mac_IM'] = 'IM'
ccomp$cnmf_celltypes = factor (ccomp$cnmf_celltypes, levels = rev(c('Mono','TREM2','SPP1','IFN_CXCLs','cDCs','IM')))
ccomp$module = archp$mod_2
ccomp = ccomp[!is.na(ccomp$cnmf_celltypes),]
rp <- ggplot(ccomp, aes(x = module, y = cnmf_celltypes, fill = ..x..)) +
  geom_density_ridges_gradient(
  scale = 3,
  rel_min_height = 0.01,
  linewidth = 0.4,
  color='white',
  alpha = 0.3
) +

  scale_fill_viridis_c(option = "C") +  # Optional: nice color gradient
  theme_ridges() +                      # Optional: clean ridge plot theme
  theme(legend.position = "right")     # Adjust legend position
#   theme_classic() + facet_wrap (~sample, ncol=5)
pdf (file.path ('Plots','cnmf_inflammation_module_ridge_plots.pdf'), width = 5,height=3)
rp
dev.off()



# Add modules to srt object and plot them on UMAPs ####
srt = ModScoreCor (
    seurat_obj = srt,
    geneset_list = shared_cnmf_genes,
    cor_threshold = NULL,
    pos_threshold = NULL, # threshold for fetal_pval2
    listName = 'shared_cnmf', outdir = NULL)


# Ridge plots using SCENIC regulons from paired scRNA-seq
scenic_name = 'monomac_programs'
auc_mtx = read.csv(file.path(paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/NTP_multiome/SCENIC/vg_5000_mw_tss500bp/',scenic_name), 'auc_mtx.csv'), header=T)

rownames (auc_mtx) = auc_mtx[,1]
rownames(auc_mtx) = gsub ('\\.','-', rownames(auc_mtx))
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

km = readRDS (file.path ('TF_activity_modules.rds'))
regulon_TFs_in_modules = list(
  km1 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 1])],
  km2 = colnames(auc_mtx)[colnames(auc_mtx) %in% names(km$cluster[km$cluster == 2])]
  )
auc_mtx = auc_mtx[, colnames(auc_mtx) %in% regulon_TFs_in_modules$km2]



srt$mod_2 = unname(rowMeans (auc_mtx))[match(colnames(srt), rownames(auc_mtx))]

# # Try with ridge plots ####
library (ggridges)
library (ggplot2)
library (viridis)
library (tidyr)
#library(hrbrthemes)

# Plot
ccomp = as.data.frame (srt@meta.data)
ccomp$cnmf_celltypes = ccomp$shared_cnmf_r_max
#ccomp = ccomp[ccomp$cnmf_celltypes %in% c('cDCs'),]
#ccomp = as.data.frame (archp@cellColData)
#ccomp = ccomp[ccomp$cnmf_celltypes %in% c('cDCs'),]
# ccomp$cnmf_celltypes = NA
# ccomp$cnmf_celltypes[ccomp$celltype2 == 'Tumor_Mono CD14'] = 'Mono'
# ccomp$cnmf_celltypes[ccomp$celltype2 == 'Tumor_Mye_DC2'] = 'cDCs'
# ccomp$cnmf_celltypes[ccomp$celltype2 == 'Tumor_Mye_Mac_TREM2'] = 'TREM2'
# ccomp$cnmf_celltypes[ccomp$celltype2 == 'Tumor_Mye_Mac_SPP1'] = 'SPP1'
# ccomp$cnmf_celltypes[ccomp$celltype2 == 'Tumor_Mye_Mac_IM'] = 'IM'
# ccomp$cnmf_celltypes = factor (ccomp$cnmf_celltypes, levels = rev(c('Mono','TREM2','SPP1','IFN_CXCLs','cDCs','IM')))

cell_subsets_order = c("TAM_interstitial","TAM_MARCO","cDCs","TAM_TREM2","TAM_CXCLs","Mono_CD16","Mono_CD14")
ccomp$cnmf_celltypes = factor (ccomp$shared_cnmf_r_max, levels = cell_subsets_order)
ccomp$module = srt$mod_2

rp <- ggplot(ccomp, aes(x = module, y = cnmf_celltypes, fill = ..x..)) +
  geom_density_ridges_gradient(
  scale = 3,
  rel_min_height = 0.01,
  linewidth = 0.4,
  color='white',
  alpha = 0.3
) +

  scale_fill_gradientn (colors = palette_expression) +  # Optional: nice color gradient
  theme_ridges() +                      # Optional: clean ridge plot theme
  theme(legend.position = "right")     # Adjust legend position
#   theme_classic() + facet_wrap (~sample, ncol=5)
pdf (file.path ('Plots','cnmf_inflammation_module_ridge_scrna_plots.pdf'), width = 7,height=3)
rp
dev.off()






### Run TF correlation to identify TF modules across TNK cells ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
#mMat = scale(as.matrix(mMat))#[selected_TF,])

km = readRDS (file.path ('TF_activity_modules.rds'))
tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[rownames(mMat) %in% names(km$cluster[km$cluster == x]),]))
# tf_modules = c(tf_modules, list(AP1 = colMeans (t(scale(mMat))[rownames(t(mMat)) %in% c('JUN','FOSB','FOS','BACH1','SMARCC1','FOSL2','JUND','JDP2','BATF'),])))
#tf_module_infl = colMeans (t)


names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = do.call (cbind, tf_modules)
archp@cellColData = archp@cellColData[!colnames(archp@cellColData) %in% paste0('mod_',unique(km$cluster))]
archp@cellColData = cbind (archp@cellColData, tf_modules) 

ccomp = as.data.frame (archp@cellColData)
ccomp$module = ccomp$mod_2
rownames (ccomp) = gsub ('fragments.tsv.gz#', '', rownames (ccomp))
ccomp$cnmf_celltypes = srt$shared_cnmf_r_max[match (rownames (ccomp),rownames (srt@meta.data))]
ccomp$cnmf_celltypes = factor (ccomp$cnmf_celltypes, levels = cell_subsets_order)
rp <- ggplot(ccomp, aes(x = module, y = cnmf_celltypes, fill = ..x..)) +
  geom_density_ridges_gradient(
  scale = 3,
  rel_min_height = 0.01,
  linewidth = 0.4,
  color='white',
  alpha = 0.3
) +

  scale_fill_viridis_c(option = "C") +  # Optional: nice color gradient
  theme_ridges() +                      # Optional: clean ridge plot theme
  theme(legend.position = "right")     # Adjust legend position
#   theme_classic() + facet_wrap (~sample, ncol=5)
pdf (file.path ('Plots','cnmf_inflammation_module_ridge_plots.pdf'), width = 7,height=3)
rp
dev.off()






# ### Try using metacells for each cnmf ####
# sams = c('P1','P10','P11','P12','P13','P14','P23','P5') # Select samples that have at least 100 endothelial cells
# df = list()

# library (zoo)
# bin_width <- 30   # Number of observations per bin
# overlap <- 30
# cnmf_mods = c('Mono','cDCs','TREM2','SPP1','IFN_CXCLs','IM')
# #archp$Sample2 = 'sample'
# #sams = 'sample'

# for (sam in sams)
# {
#   cnmf_l = list()
#   for (cnmf_mod in cnmf_mods)
#   {
#   metacells_order = cnmf_scatac[,cnmf_mod][archp$Sample == sam]
#   metacells_order = order (-metacells_order)
#   mMat_sam = mMat[archp$Sample == sam, ]
#   mMat_sam = mMat_sam[metacells_order,]
#   cnmf_sam = cnmf_scatac[,cnmf_mod][archp$Sample == sam][metacells_order]
  
#  cnmf_l[[cnmf_mod]] <- cor (
#   rollapply (mMat_sam, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"),
#   rollapply (cnmf_sam, width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"), 
#   method = 'spearman')  
#   }
# df[[sam]] = as.data.frame (do.call (cbind, cnmf_l))
# colnames (df[[sam]]) = cnmf_mods
# #df[[sam]]$sample = sam
# }

# # Convert list to 3D array: rows x cols x n_matrices
# array_3d <- array(unlist(df), dim = c(nrow(df[[1]]), ncol(df[[1]]), length(df)))

# # Apply median across the third dimension (i.e., across matrices)
# median_matrix <- apply(array_3d, c(1, 2), median)
# rownames (median_matrix) = rownames (df[[1]])
# colnames (median_matrix) = colnames (df[[1]])
# top_5 = unlist(lapply (cnmf_mods, function(x) head (rownames(median_matrix[order(-median_matrix[,x]),]),5)))
# top_5 = unlist(lapply (cnmf_mods, function(x) head (rownames(df[[1]][order(-df[[1]][,x]),]),5)))

# DAM_hm = Heatmap (df[[1]][top_5,cnmf_mods], 
#           #row_labels = colnames (mMat_mg),
#           #column_title = paste('top',top_genes),
#           clustering_distance_columns = 'euclidean',
#           clustering_distance_rows = 'euclidean',
#           cluster_rows = F,
#           #col = pals_heatmap[[5]],
#           cluster_columns=F,#col = pals_heatmap[[1]],
#           row_names_gp = gpar(fontsize = 8),
#           column_names_gp = gpar(fontsize = 8),
#           column_names_rot = 45,
#           name = 'chromVAR',
#           #rect_gp = gpar(col = "white", lwd = .5),
#           border=TRUE,
#           col = rev(palette_deviation)#,
#           #width = unit(2, "cm")
#           #right_annotation = motif_ha
#           )

# pdf (file.path ('Plots','cnmf_clusters_DAM_metacells_heatmap.pdf'), width = 3,height=4)
# draw (DAM_hm)
# dev.off()


# # Try with PCA on mMat ####
library (ggridges)
library (ggplot2)
library (viridis)
library (tidyr)
#library(hrbrthemes)
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name


# # Filter by RNA expression ####
metaGroupName = 'celltype_simplified'
min_exp = .1
active_TFs = exp_genes (srt, rownames(mMat), min_exp = min_exp, metaGroupName)
mMat = t(scale (mMat[active_TFs, ]))
library(uwot)
library(ggplot2)
library(patchwork)


tail (rownames(archp@cellColData))
tail (colnames(srt))

#srt = srt[,gsub ('fragments.tsv.gz#','',rownames(archp@cellColData))]

p <- prcomp(mMat, center = TRUE, scale. = TRUE)
plot_df <- data.frame(p$x, celltype = archp$shared_cnmf_r_max)


# df_sub <- plot_df[plot_df$celltype %in% 
#                     c("Mye_Mac_TREM2","Mye_CD14 Mono","Mye_Mac_SPP1",
#                       "Mye_Mac_IM"), ]


df_sub <- plot_df[plot_df$celltype %in% 
                    c("Mono_CD14","TAM_interstitial","TAM_TREM2",
                      "TAM_MARCO","TAM_CXCLs"), ]



#df_sub = plot_df
#-------------------------
# 1. MAIN SCATTER
#-------------------------
p_scatter <- ggplot(df_sub, aes(PC1, PC2, color = celltype)) +
  geom_point(alpha = 0.4, size = 0.1) +
  geom_density_2d(size = 0.4) +
  scale_color_manual(values = palette_myeloid) +
  # scale_x_continuous(limits = x_range, expand = c(0,0)) +
  # scale_y_continuous(limits = y_range, expand = c(0,0)) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )


p_density_x <- ggplot(df_sub, aes(PC1, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) + gtheme_no_rot

p_density_y <- ggplot(df_sub, aes(PC2, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) + gtheme_no_rot

pdf (file.path ('Plots','scrna_momac_deviation_pca.pdf'), width = 5,height=6)
p_density_x + plot_spacer() + p_scatter + p_density_y + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4),
    guides = "collect"      # <- collect legends into one
  ) &
  theme(
    legend.position = "bottom",  # or "top"
    legend.box = "vertical"
  )

dev.off()




scenic_name = 'monomac_programs'
auc_mtx = read.csv(file.path(paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/NTP_multiome/SCENIC/vg_5000_mw_tss500bp/',scenic_name), 'auc_mtx.csv'), header=T)

rownames (auc_mtx) = auc_mtx[,1]
rownames(auc_mtx) = gsub ('\\.','-', rownames(auc_mtx))
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

auc_mtx = auc_mtx[rownames(auc_mtx) %in% colnames(srt),]
p <- prcomp(auc_mtx, center = TRUE, scale. = TRUE)
plot_df <- data.frame(p$x, celltype = srt$shared_cnmf_r_max[rownames(auc_mtx)])


# df_sub <- plot_df[plot_df$celltype %in% 
#                     c("Mye_Mac_TREM2","Mye_CD14 Mono","Mye_Mac_SPP1",
#                       "Mye_Mac_IM"), ]


df_sub <- plot_df[plot_df$celltype %in% 
                    c("Mono_CD14","TAM_interstitial","TAM_TREM2",
                      "TAM_MARCO","TAM_CXCLs"), ]



#df_sub = plot_df
#-------------------------
# 1. MAIN SCATTER
#-------------------------
p_scatter <- ggplot(df_sub, aes(PC1, PC2, color = celltype)) +
  geom_point(alpha = 0.4, size = 0.1) +
  geom_density_2d(size = 0.4) +
  scale_color_manual(values = palette_myeloid) +
  # scale_x_continuous(limits = x_range, expand = c(0,0)) +
  # scale_y_continuous(limits = y_range, expand = c(0,0)) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )


p_density_x <- ggplot(df_sub, aes(PC1, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) + gtheme_no_rot

p_density_y <- ggplot(df_sub, aes(PC2, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) + gtheme_no_rot

pdf (file.path ('Plots','scrna_momac_regulon_pca.pdf'), width = 5,height=6)
p_density_x + plot_spacer() + p_scatter + p_density_y + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4),
    guides = "collect"      # <- collect legends into one
  ) &
  theme(
    legend.position = "bottom",  # or "top"
    legend.box = "vertical"
  )

dev.off()






#-------------------------
# 1. MAIN SCATTER use only NFKB1 and JUNB
#-------------------------

plot_df <- data.frame(auc_mtx, celltype = srt$shared_cnmf_r_max[rownames(auc_mtx)])

df_sub <- plot_df[plot_df$celltype %in% 
                    c("Mono_CD14","TAM_interstitial","TAM_TREM2",
                      "TAM_MARCO","TAM_CXCLs"), ]


p_scatter <- ggplot(df_sub, aes(NFKB1, JUNB, color = celltype)) +
  geom_point(alpha = 0.4, size = 0.1) +
  geom_density_2d(size = 0.4) +
  scale_color_manual(values = palette_myeloid) +
  # scale_x_continuous(limits = x_range, expand = c(0,0)) +
  # scale_y_continuous(limits = y_range, expand = c(0,0)) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )


p_density_x <- ggplot(df_sub, aes(NFKB1, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) + gtheme_no_rot

p_density_y <- ggplot(df_sub, aes(JUNB, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) + gtheme_no_rot

pdf (file.path ('Plots','scrna_momac_regulon_NFKB1_JUNB_scatter.pdf'), width = 5,height=6)
p_density_x + plot_spacer() + p_scatter + p_density_y + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4),
    guides = "collect"      # <- collect legends into one
  ) &
  theme(
    legend.position = "bottom",  # or "top"
    legend.box = "vertical"
  )

dev.off()






#-------------------------
# 1. MAIN SCATTER use only NFKB1 and JUNB
#-------------------------

plot_df <- data.frame(auc_mtx, celltype = srt$celltype_simplified[rownames(auc_mtx)])

select_cells = c('Mye_CD14 Mono','Mye_CD14 Mono_activated','Mye_Mac_IM','Mye_Mac_HIF1A','Mye_Mac_SPP1','Mye_Mac_TRM','Mye_Mac_transition','Mye_Mac_TREM2')
df_sub <- plot_df[plot_df$celltype %in% 
                    select_cells, ]


p_scatter <- ggplot(df_sub, aes(NFKB1, JUNB, color = celltype)) +
  geom_point(alpha = 0.4, size = 0.1) +
  geom_density_2d(size = 0.4) +
  #scale_color_manual(values = palette_myeloid) +
  # scale_x_continuous(limits = x_range, expand = c(0,0)) +
  # scale_y_continuous(limits = y_range, expand = c(0,0)) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )


p_density_x <- ggplot(df_sub, aes(NFKB1, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  #scale_fill_manual(values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) + gtheme_no_rot

p_density_y <- ggplot(df_sub, aes(JUNB, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  #scale_fill_manual(values = palette_myeloid) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) + gtheme_no_rot

pdf (file.path ('Plots','scrna_momac_regulon_NFKB1_JUNB_samarth_anno_scatter.pdf'), width = 8,height=8)
p_density_x + plot_spacer() + p_scatter + p_density_y + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4),
    guides = "collect"      # <- collect legends into one
  ) &
  theme(
    legend.position = "bottom",  # or "top"
    legend.box = "vertical"
  )

dev.off()








# #-------------------------
# # 0. RUN UMAP
# #-------------------------
# set.seed(123)  # for reproducibility
# library(uwot)
# umap_res <- umap(
#   t(mMat),
#   n_neighbors = 30,
#   min_dist = 0.3,
#   metric = "cosine",
#   scale = TRUE
# )

# plot_df <- data.frame(
#   UMAP1 = umap_res[,1],
#   UMAP2 = umap_res[,2],
#   celltype = archp$shared_cnmf_r_max
# )

# df_sub <- plot_df[plot_df$celltype %in%
#                     c("Mono_CD14","TAM_interstitial","TAM_TREM2",
#                       "TAM_MARCO","TAM_CXCLs"), ]


# #-------------------------
# # 1. MAIN SCATTER
# #-------------------------
# p_scatter <- ggplot(df_sub, aes(UMAP1, UMAP2, color = celltype)) +
#   geom_point(alpha = 0.4, size = 0.1) +
#   geom_density_2d(size = 0.4) +
#   scale_color_manual(values = palette_myeloid) +
#   gtheme_no_rot +
#   theme(
#     legend.position = "right",
#     plot.margin = margin(0,0,0,0)
#   )


# #-------------------------
# # 2. TOP DENSITY (UMAP1)
# #-------------------------
# p_density_x <- ggplot(df_sub, aes(UMAP1, fill = celltype)) +
#   geom_density(alpha = 0.3, color = NA) +
#   scale_fill_manual(values = palette_myeloid) +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     axis.title.x = element_blank(),
#     axis.text.x  = element_blank(),
#     axis.ticks.x = element_blank()
#   ) +
#   gtheme_no_rot


# #-------------------------
# # 3. RIGHT DENSITY (UMAP2)
# #-------------------------
# p_density_y <- ggplot(df_sub, aes(UMAP2, fill = celltype)) +
#   geom_density(alpha = 0.3, color = NA) +
#   scale_fill_manual(values = palette_myeloid) +
#   coord_flip() +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     axis.title.y = element_blank(),
#     axis.text.y  = element_blank(),
#     axis.ticks.y = element_blank()
#   ) +
#   gtheme_no_rot


# #-------------------------
# # 4. SAVE PDF
# #-------------------------
# pdf(file.path("Plots","scatac_cnmf_inflammation_module_umap_ridge_plots.pdf"),
#     width = 5, height = 5)

# p_density_x + plot_spacer() + p_scatter + p_density_y +
#   plot_layout(
#     ncol = 2, nrow = 2,
#     widths  = c(4, 1),
#     heights = c(1, 4),
#     guides = "collect"
#   ) &
#   theme(
#     legend.position = "bottom",
#     legend.box = "horizontal"
#   )

# dev.off()



# umap_df <- data.frame(UMAP1 = umap_res[,1],
#                       UMAP2 = umap_res[,2])

# cors_UMAP1 <- apply(t(mMat), 2, function(x) cor(x, umap_df$UMAP1))
# cors_UMAP2 <- apply(t(mMat), 2, function(x) cor(x, umap_df$UMAP2))

# loadings_like <- data.frame(peak = colnames(mMat),
#                             cor_UMAP1 = cors_UMAP1,
#                             cor_UMAP2 = cors_UMAP2)



# #-------------------------
# # 0. RUN UMAP
# #-------------------------
# set.seed(123)  # for reproducibility
# library(uwot)
# umap_res <- umap(
#   auc_mtx,
#   n_neighbors = 30,
#   min_dist = 0.3,
#   metric = "cosine",
#   scale = TRUE
# )

# plot_df <- data.frame(
#   UMAP1 = umap_res[,1],
#   UMAP2 = umap_res[,2],
#   celltype = srt$shared_cnmf_r_max
# )

# df_sub <- plot_df[plot_df$celltype %in%
#                     c("Mono_CD14","TAM_interstitial","TAM_TREM2",
#                       "TAM_MARCO","TAM_CXCLs"), ]


# #-------------------------
# # 1. MAIN SCATTER
# #-------------------------
# p_scatter <- ggplot(df_sub, aes(UMAP1, UMAP2, color = celltype)) +
#   geom_point(alpha = 0.4, size = 0.1) +
#   geom_density_2d(size = 0.4) +
#   scale_color_manual(values = palette_myeloid) +
#   gtheme_no_rot +
#   theme(
#     legend.position = "right",
#     plot.margin = margin(0,0,0,0)
#   )


# #-------------------------
# # 2. TOP DENSITY (UMAP1)
# #-------------------------
# p_density_x <- ggplot(df_sub, aes(UMAP1, fill = celltype)) +
#   geom_density(alpha = 0.3, color = NA) +
#   scale_fill_manual(values = palette_myeloid) +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     axis.title.x = element_blank(),
#     axis.text.x  = element_blank(),
#     axis.ticks.x = element_blank()
#   ) +
#   gtheme_no_rot


# #-------------------------
# # 3. RIGHT DENSITY (UMAP2)
# #-------------------------
# p_density_y <- ggplot(df_sub, aes(UMAP2, fill = celltype)) +
#   geom_density(alpha = 0.3, color = NA) +
#   scale_fill_manual(values = palette_myeloid) +
#   coord_flip() +
#   theme_classic() +
#   theme(
#     legend.position = "none",
#     axis.title.y = element_blank(),
#     axis.text.y  = element_blank(),
#     axis.ticks.y = element_blank()
#   ) +
#   gtheme_no_rot


# #-------------------------
# # 4. SAVE PDF
# #-------------------------
# pdf(file.path("Plots","scrna_cnmf_inflammation_module_umap_ridge_plots.pdf"),
#     width = 5, height = 5)

# p_density_x + plot_spacer() + p_scatter + p_density_y +
#   plot_layout(
#     ncol = 2, nrow = 2,
#     widths  = c(4, 1),
#     heights = c(1, 4),
#     guides = "collect"
#   ) &
#   theme(
#     legend.position = "bottom",
#     legend.box = "horizontal"
#   )

# dev.off()




### Filter TFs by genescore correlation ####
seGroupMotif <- getGroupSE (ArchRProj = archp, useMatrix = "MotifMatrix", groupBy = "Clusters_H")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGSM_MM <- correlateMatrices(
    ArchRProj = archp,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
#corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.4 & corGSM_MM$padj < 0.05 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
active_TF = sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
corGSM_MM = corGSM_MM[corGSM_MM$GeneScoreMatrix_name %in% active_TFs,]
corGSM_MM = data.frame(corGSM_MM)

pdf ()
p = ggplot(corGSM_MM, aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )+
  # Add labels for TFRegulator == YES
  geom_text_repel(
    data = subset(corGSM_MM, TFRegulator == "YES"),
    aes(label = MotifMatrix_name),
    size = 3,
    max.overlaps = 20
  )
dev.off()

pdf (file.path ('Plots','positive_expression_TFs.pdf'), width = 4,height=4)
p
dev.off()






















# Compare MDM score between deviation and genescore ####
mdm_sig = c('FOSL1','FOSL2','BACH1','PPARG','NFKB2','KLF12','HIVEP3','SMAD1','NFKB1','REL','RUNX1','SNAI1','RUNX2','NFAT5')
mdm_sig2 = c('FOSL1','FOSL2','BACH1','NFKB2','NFKB1','REL')
myeloid_order = c('Mye_CD14 Mono_activated','Mye_CD14 Mono','Mye_Mac_TREM2','Mye_Mac_SPP1','Mye_Mac_IM','Mye_Mac_TRM')
myeloid_order = c('Mono_CD14','TAM_CXCLs','TAM_TREM2','TAM_MARCO','TAM_interstitial')

# Add metacolumns of average TF modules activity ####
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = scale(assays (mSE)[[1]])
rownames (mMat) = rowData (mSE)$name
#mMat = as.data.frame (t(scale(as.matrix(mMat)[mdm_sig,])))
mMat = as.data.frame (t(as.matrix(mMat)[mdm_sig2,]))
mMat = aggregate (mMat, by = list(celltype_lv3 = archp$shared_cnmf_r_max), mean)
rownames (mMat) = mMat[,1]
mMat = mMat[,-1]
mMat = as.data.frame (scale(mMat[myeloid_order,]))
pdf (file.path ('Plots','MDM_score_activity_heatmap.pdf'),2.5, 2.5)
Heatmap (t(mMat), 
  col = rev(palette_deviation), 
  cluster_columns=F,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  border=T)
dev.off()

#mMat = as.data.frame (scale (mMat))
mMat$celltype = rownames (mMat)
mMat = gather (mMat, TF, expression, 1:(ncol(mMat)-1))
mMat$type = 'deviation'

# Get genescore ####
if (!exists('gsSE')) gsSE = fetch_mat (archp, 'GeneScore')
gsMat = scale(assays (gsSE)[[1]])
rownames (gsMat) = rowData (gsSE)$name
#gsMat = scale(gsMat[rownames (gsMat) %in% mdm_sig,])
gsMat = gsMat[rownames (gsMat) %in% mdm_sig2,]
gsMat = t(gsMat)
gsMat = as.data.frame (gsMat)
gsMat = aggregate (gsMat, by = list(celltype_lv3 = archp$shared_cnmf_r_max), mean)
rownames (gsMat) = gsMat[,1]
gsMat = gsMat[,-1]
gsMat = as.data.frame (scale(gsMat[myeloid_order,]))
pdf (file.path ('Plots','MDM_score_genescore_heatmap.pdf'),2.5,2.5)
Heatmap (t(gsMat), 
  col = rev(palette_deviation), 
  cluster_columns=F,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  border=T)
dev.off()

#gsMat = as.data.frame (scale (gsMat))
gsMat$celltype = rownames (gsMat)
gsMat = gather (gsMat, TF, expression, 1:(ncol(gsMat)-1))
gsMat$type = 'genescore'


### Use scRNA ####
#gsMat = gather (gsMat, TF, expression, 2:ncol(gsMat))
ps = log2(as.data.frame (AverageExpression (srt, features = mdm_sig2, group.by = 'shared_cnmf_r_max')[[1]]) +1)
colnames (ps) = gsub ('-','_',colnames(ps))
ps = as.data.frame (scale(t(ps)))
ps$celltype = rownames(ps)
ps = gather (ps, TF, expression, 1:(ncol(ps)-1))
ps$type = 'expression'

combined = rbind (mMat, gsMat, ps)


df = combined
df$celltype <- factor(df$celltype, levels = unique(df$celltype))
#df$TF = paste0(df$TF, df$type)
df = df[df$celltype %in% myeloid_order,]
df$celltype = factor (df$celltype,myeloid_order)

p <- list()

for (g in mdm_sig2) {
  p[[g]] <- ggplot(
    df[df$TF == g, ],
    aes(
      x = celltype,
      y = expression,
      group = type,          # ← key change
      color = type
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    ggtitle(g) + gtheme + facet_wrap (~type, scales = 'free_y')
}

pdf (file.path ('Plots','deviation_vs_genescore_MDM_sig_linePlot.pdf'), width=14,height=7)
wrap_plots (p)
dev.off()


library(dplyr)
df = combined
#df = df[df$TF %in% c('FOSL1','FOSL2','BACH1','PPARG','NFKB2','NFKB1'),]
df = df[df$celltype %in% myeloid_order,]
df$celltype = factor (df$celltype,myeloid_order)

sum_df <- df %>%
  group_by(celltype, type) %>%
  summarise(
    median_expr = median(expression, na.rm = TRUE),
    mean_expr   = mean(expression, na.rm = TRUE),
    q25         = quantile(expression, 0.25, na.rm = TRUE),
    q75         = quantile(expression, 0.75, na.rm = TRUE),
    sd          = sd(expression, na.rm = TRUE),
    n           = n(),
    sem         = sd / sqrt(n),
    .groups = "drop"
  )
p =ggplot(
  sum_df,
  aes(
    x = celltype,
    y = mean_expr,
    group = type,
    color = type,
    fill  = type
  )
) +
  geom_ribbon(
    aes(ymin = q25, ymax = q75),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  labs(
    y = "Mean TF expression",
    x = "Cell type"
  ) + facet_wrap (~type, scales = 'free_y') + gtheme

pdf (file.path ('Plots','deviation_vs_genescore_MDM_sig_mean_linePlot.pdf'), width = 5,height=3)
p
dev.off()







### Make trajectory of metacells based on NFKB deviation and compare with genescore and expression ####

library(dplyr)
library(zoo)
library(ggplot2)

###---------------------------------------------------------
### 1. Load inputs and define samples
###---------------------------------------------------------

selected_TF <- c('FOSL1','FOSL2','BACH1','NFKB2','NFKB1','REL')
select_cells = c('Mye_CD14 Mono','Mye_CD14 Mono_activated','Mye_Mac_IM','Mye_Mac_HIF1A','Mye_Mac_SPP1','Mye_Mac_TRM','Mye_Mac_transition','Mye_Mac_TREM2')
#sarc_module <- "cNMF20"

archp_meta <- as.data.frame(archp@cellColData)
archp_meta = archp_meta[archp_meta$celltype_simplified %in% select_cells,]

# Filter samples: remove normal lung, low-cell samples, and outliers
sams <- unique(archp_meta$Sample)
#sams <- sams[!sams %in% c("normal1","normal2","normal3","normal4","P3","P13","P11_HOX")]


###---------------------------------------------------------
### 2. Extract and scale motif deviations (TF activity)
###---------------------------------------------------------

# Fetch motif matrix only if not already loaded
if (!exists ('mSE')) mSE <- fetch_mat(archp, "Motif")
mMat <- scale (assays(mSE)[[1]])
rownames(mMat) <- rowData(mSE)$name

# Keep only selected TFs and z-score normalize
#mMat <- scale(mMat[selected_TF,])
mMat <- mMat[selected_TF,rownames(archp_meta)]
mMat <- t(mMat)

# Split by sample
mMat <- lapply(sams, function(s) mMat[archp_meta$Sample == s,])
names(mMat) <- sams


###---------------------------------------------------------
### 3. Extract and scale GeneScore matrix
###---------------------------------------------------------

if (!exists ('gsSE')) gsSE <- fetch_mat(archp, "GeneScore")
gsMat <- scale(assays(gsSE)[[1]])
rownames(gsMat) <- rowData(gsSE)$name

# Keep only selected TFs and z-score normalize
gsMat <- gsMat[selected_TF,rownames(archp_meta)]
gsMat <- t(gsMat)

# Split by sample
gsMat <- lapply(sams, function(s) gsMat[archp_meta$Sample == s,])
names(gsMat) <- sams




###---------------------------------------------------------
### 3. Extract and scale Expression matrix
###---------------------------------------------------------

eMat <- srt@assays$RNA@data
#eMat <- scale(eSE)
archp_meta2 = archp_meta
rownames (archp_meta2) = gsub ('fragments.tsv.gz#','',rownames(archp_meta2))
# Keep only selected TFs and z-score normalize
eMat <- eMat[selected_TF,rownames(archp_meta2)]
eMat <- t(eMat)

# Split by sample
eMat <- lapply(sams, function(s) eMat[archp_meta2$Sample == s,])
names(eMat) <- sams

### Take also regulon scores 
scenic_name = 'monomac_programs'
auc_mtx = read.csv(file.path(paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/NTP_multiome/SCENIC/vg_5000_mw_tss500bp/',scenic_name), 'auc_mtx.csv'), header=T)

rownames (auc_mtx) = auc_mtx[,1]
rownames(auc_mtx) = gsub ('\\.','-', rownames(auc_mtx))
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

auc_mtx = auc_mtx[rownames(archp_meta2),colnames(auc_mtx) %in% selected_TF]
# Split by sample
auc_mtx <- lapply(sams, function(s) auc_mtx[archp_meta2$Sample == s,])
names(auc_mtx) <- sams

###---------------------------------------------------------
### 5. Order cells by sarcomatoid module + smooth (rolling mean)
###---------------------------------------------------------

bin_width <- 10   # Window size
overlap   <- 10   # Step size
sarc_module = 'NFKB1'
# Helper function: ordered → rolling average
rollavg_df <- function(mat, order_by, bin_width, overlap) {
  mat_ordered <- mat[order(-order_by),]
  as.data.frame(lapply(as.data.frame(mat_ordered), function(col)
    zoo::rollapply(
      col, 
      width = bin_width, 
      FUN = mean, 
      by = overlap, 
      partial = TRUE, 
      align = "left"
    )
  ))
}

# TF activity (motif)
mMat_ordered_avg <- lapply(sams, function(s)
  rollavg_df(mMat[[s]], mMat[[s]][,sarc_module], bin_width, overlap)
)
names (mMat_ordered_avg) = sams
# GeneScores
gsMat_ordered_avg <- lapply(sams, function(s) {
  df <- rollavg_df(gsMat[[s]], mMat[[s]][,sarc_module], bin_width, overlap)
#  df[is.na(df)] <- 0  # Handle NaNs (e.g. HIC2 issues)
  df
})
names (gsMat_ordered_avg) = sams

eMat_ordered_avg <- lapply(sams, function(s) {
  df <- rollavg_df(eMat[[s]], mMat[[s]][,sarc_module], bin_width, overlap)
#  df[is.na(df)] <- 0  # Handle NaNs (e.g. HIC2 issues)
  df
})
names (eMat_ordered_avg) = sams

regMat_ordered_avg <- lapply(sams, function(s) {
  df <- rollavg_df(auc_mtx[[s]], mMat[[s]][,sarc_module], bin_width, overlap)
#  df[is.na(df)] <- 0  # Handle NaNs (e.g. HIC2 issues)
  df
})
names (regMat_ordered_avg) = sams

library(ggplot2)

gene = 'NFKB1'
plot_nfkb1 <- function(mat_list, gene) {
  p <- list()

  for (nm in names(mat_list)) {
    df_long <- data.frame(
      x = seq_len(nrow(mat_list[[nm]])),
      NFKB1 = mat_list[[nm]][, gene]
    )

    p[[nm]] <- ggplot(df_long, aes(x = x, y = NFKB1)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2) +
      theme_bw() +
      labs(
        title = nm,
        x = "Row index",
        y = gene
      )
  }
  p
}

p_mMat <- plot_nfkb1(mMat_ordered_avg,gene)
p_gsMat <- plot_nfkb1(gsMat_ordered_avg,gene)
p_eMat <- plot_nfkb1(eMat_ordered_avg,gene)
p_regMat <- plot_nfkb1(regMat_ordered_avg,gene)

pdf (file.path('Plots',paste0(gene,'_pseudotime_lineplot.pdf')))
wrap_plots(p_mMat)
wrap_plots(p_gsMat)
wrap_plots(p_eMat)
wrap_plots(p_regMat)
dev.off()










### Take also regulon scores 
scenic_name = 'monomac_programs'
auc_mtx = read.csv(file.path(paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/NTP_multiome/SCENIC/vg_5000_mw_tss500bp/',scenic_name), 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
rownames(auc_mtx) = gsub ('\\.','-', rownames(auc_mtx))
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

auc_mtx = auc_mtx[rownames(archp_meta2),colnames(auc_mtx) %in% selected_TF]


# Fetch motif matrix only if not already loaded
if (!exists ('mSE')) mSE <- fetch_mat(archp, "Motif")
mMat <- scale (assays(mSE)[[1]])
rownames(mMat) <- rowData(mSE)$name

# Keep only selected TFs and z-score normalize
#mMat <- scale(mMat[selected_TF,])
mMat <- mMat[selected_TF,rownames(archp_meta)]
mMat <- t(mMat)

cor_mat = cor (mMat, auc_mtx)

pdf (file.path ('Plots','regulon_deviation_heatmap.pdf'))
Heatmap (cor_mat, col = palette_deviation_fun (cor_mat))
dev.off()

mMat2 = mMat
rownames(mMat2) = gsub ('fragments.tsv.gz#','',rownames(mMat2))
srt$NFKB1_dev = mMat2[,'NFKB1']
srt$NFKB1_regulon = auc_mtx[,'NFKB1'][match(colnames(srt), rownames(auc_mtx))]
srt$FOSL1_dev = mMat2[,'FOSL1']
srt$FOSL1_regulon = auc_mtx[,'FOSL1'][match(colnames(srt), rownames(auc_mtx))]
srt$FOSL2_dev = mMat2[,'FOSL2']
srt$FOSL2_regulon = auc_mtx[,'FOSL2'][match(colnames(srt), rownames(auc_mtx))]

pdf (file.path ('Plots','regulon_deviation_heatmap.pdf'),4,width=6)
reductionName = 'sampleID_harmony_umap'
wrap_plots (fp (srt, gene = c('NFKB1_regulon','NFKB1_dev','NFKB1')))
wrap_plots (fp (srt, gene = c('FOSL1_regulon','FOSL1_dev','FOSL1')))
wrap_plots (fp (srt, gene = c('FOSL2_regulon','FOSL2_dev','FOSL2','TREM2')))
DimPlot (srt, group.by = 'celltype_simplified')
DimPlot (srt, group.by = 'shared_cnmf_r_max')
dev.off()















 ### Plot meta cnmf coexpression ####
#mono_mac_cnmfs = c('TREM2','Monocytes','SPP1','IL1B','IFN','IM','C1Q')
### Take also regulon scores 
scenic_name = 'monomac_programs'
auc_mtx = read.csv(file.path(paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/NTP_multiome/SCENIC/vg_5000_mw_tss500bp/',scenic_name), 'auc_mtx.csv'), header=T)

rownames (auc_mtx) = auc_mtx[,1]
rownames(auc_mtx) = gsub ('\\.','-', rownames(auc_mtx))
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

library (hdWGCNA)
force = F
# table (srt$sampleID)
# srt$sampleID
#srt_tam = srt[,srt$celltype2 == 'TAMs']
select_cells = c('Mye_CD14 Mono','Mye_CD14 Mono_activated','Mye_Mac_IM','Mye_Mac_HIF1A','Mye_Mac_SPP1','Mye_Mac_TRM','Mye_Mac_transition','Mye_Mac_TREM2')
if (!file.exists ('metacells.rds') | force)
  {
  srt_sub = srt[,srt$celltype_simplified %in% select_cells]
  srt_sub = srt_sub[,srt_sub$sampleID %in% c('Lu952T')]

reductionName = 'umap'
reductionSave = 'pca'
reductionGraphKnn = 'RNA_knn'
reductionGraphSnn = 'RNA_snn' 
sigPCs=30
srt_sub  = RunUMAP (object = srt_sub , reduction = reductionSave, dims = 1:sigPCs)
  
pdf (file.path ('Plots','Lu952T_umap.pdf'))
DimPlot (srt_sub, reduction = 'umap', group.by = 'celltype_simplified')
dev.off()

  srt_sub
  srt_sub <- SetupForWGCNA(
    srt_sub,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "metacells" # the name of the hdWGCNA experiment
  )
  # construct metacells  in each group
  srt_sub <- MetacellsByGroups(
    seurat_obj = srt_sub,
    group.by = c("shared_cnmf_r_max"), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'umap', # select the dimensionality reduction to perform KNN on
    k = 10, # nearest-neighbors parameter
    max_shared = 7, # maximum number of shared cells between two metacells
    ident.group = 'shared_cnmf_r_max' # set the Idents of the metacell seurat object
  )
  
  # normalize metacell expression matrix:
  srt_sub <- NormalizeMetacells (srt_sub)
  metacells = GetMetacellObject (srt_sub)
  saveRDS (metacells, 'metacells.rds')
  } else {
  metacells = readRDS ('metacells.rds') 
  }

auc_mtx_meta = lapply (metacells$cells_merged, function(x) {
  tmp = unlist(strsplit(x, ','))
  colMeans (auc_mtx[tmp,])})
auc_mtx_meta = as.data.frame(do.call (cbind, auc_mtx_meta))


if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = scale(assays (mSE)[[1]])
rownames (mMat) = rowData (mSE)$name
mMat2 = t(mMat)
rownames(mMat2) = gsub ('fragments.tsv.gz#','',rownames(mMat2))
mMat_meta = lapply (metacells$cells_merged, function(x) {
  tmp = unlist(strsplit(x, ','))
  colMeans (mMat2[tmp,])})
mMat_meta = as.data.frame(do.call (cbind, mMat_meta))
head (mMat_meta)

all (colnames (auc_mtx_meta) == colnames(mMat_meta))

metacells_mat = metacells@assays$RNA@layers$data
rownames (metacells_mat) = rownames(metacells)
colnames (metacells_mat) = colnames (metacells)
auc_mtx_meta2 = auc_mtx_meta[rownames(auc_mtx_meta) %in% selected_TF,]
mMat_meta2 = mMat_meta[selected_TF, ]
metacells_mat = metacells_mat[c('FOSL1','FOSL2','NFKB1','NFKB2','REL','TREM2','SPP1','C1QA','C3','MAF','CCL2','EREG','VCAN','IL1B'),]
colnames(metacells_mat) == colnames (mMat_meta2)
hm = Heatmap (t(scale(t(auc_mtx_meta2))), name = 'regulon')
hm2 = Heatmap (t(scale(t(mMat_meta2))), name = 'motif')
hm3 = Heatmap ((t(scale(t(metacells_mat)))), name = 'expression')

pdf (file.path ('Plots','metacells_regulons_heatmap.pdf'), width=12,height=7)
hm %v% hm3 %v% hm2
dev.off()

metacells_mat2 = as.data.frame (metacells_mat)
mMat_meta2 = as.data.frame (mMat_meta2)
mMat_meta2$type = 'deviation'
mMat_meta2$TF = rownames (mMat_meta2)
metacells_mat2$type = 'expression'
metacells_mat2$TF = rownames (metacells_mat2)
auc_mtx_meta2 = as.data.frame (auc_mtx_meta2)
auc_mtx_meta2$type = 'regulon'
auc_mtx_meta2$TF = rownames (auc_mtx_meta2)
metacells_comb = rbind (auc_mtx_meta2, mMat_meta2, metacells_mat2)
#metacells_comb = gather (metacells_comb, metacell, score, 1:(ncol(metacells_comb)-2))

library(dplyr)
library(tidyr)

gene = 'FOSL2'
tf_long <- metacells_comb %>%
  #rownames_to_column("TF") %>%
  pivot_longer(
    cols = colnames (metacells_mat),
    names_to = "cell",
    values_to = "value"
  )
nfkb1_regulon <- tf_long %>%
  filter(TF == gene, type == "regulon") %>%
  select(cell, nfkb1_regulon = value)

nfkb1_deviation <- tf_long %>%
  filter(TF == gene, type == "deviation") %>%
  select(cell, nfkb1_deviation = value)

trem2_df <- as.data.frame (t(metacells_mat))
trem2_df$cell = rownames(trem2_df)


plot_df <- trem2_df %>%
  left_join(nfkb1_regulon,  by = "cell") %>%
  left_join(nfkb1_deviation, by = "cell")


grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "IL1B"
)
rp <- ggplot(
  plot_df,
  aes(
    x = nfkb1_deviation,
    y = nfkb1_regulon
  )
) +
  geom_point(aes(color = IL1B), alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste0(gene," deviation score"),
    y = paste0(gene," regulon"),
    title = paste0(gene," deviation vs regulon on IL1B (Mono)")
  ) + gtheme

grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "MAF"
)
dp <- ggplot(
  plot_df,
  aes(
    x = nfkb1_deviation,
    y = nfkb1_regulon
  )
) +
  geom_point(aes(color = MAF), alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste0(gene,"  deviation score"),
    y = paste0(gene," regulon"),
    title = paste0(gene," deviation vs regulon on MAF (IM)")
  ) + gtheme



grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "TREM2"
)

p_main <- ggplot(
  plot_df,
  aes(
    x = nfkb1_deviation,
    y = nfkb1_regulon,
    color = TREM2
  )
) +
  geom_point(alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste0(gene," deviation score"),
    y = paste0(gene," regulon")#,
    #title = "NFKB1 deviation vs regulon colored by TREM2"
  ) + gtheme_no_rot

pdf (file.path ('Plots',paste0(gene,'regulon_deviation_scatter.pdf')), width=7, height=3.5)
wrap_plots (rp, dp, p_main)
dev.off()





grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "IL1B"
)
rp <- ggplot(
  plot_df,
  aes_string(
    x = 'nfkb1_deviation',
    y = gene
  )
) +
  geom_point(aes(color = IL1B), alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste0(gene," deviation score"),
    y = paste0(gene," expression"),
    title = paste0(gene," deviation vs regulon on IL1B (Mono)")
  ) + gtheme

grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "MAF"
)
dp <- ggplot(
  plot_df,
  aes_string(
    x = 'nfkb1_deviation',
    y = gene
  )
) +
  geom_point(aes(color = MAF), alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste0(gene," deviation score"),
    y = paste0(gene," expression"),
    title = paste0(gene," deviation vs regulon on MAF (IM)")
  ) + gtheme



grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "TREM2"
)

p_main <- ggplot(
  plot_df,
  aes_string(
    x = 'nfkb1_deviation',
    y = gene,
    color = 'TREM2'
  )
) +
  geom_point(alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste0(gene," deviation score"),
    y = paste0(gene," expression")#,
    #title = "NFKB1 deviation vs regulon colored by TREM2"
  ) + gtheme_no_rot

pdf (file.path ('Plots',paste0(gene,'_expression_deviation_scatter.pdf')), width=7, height=3.5)
wrap_plots (rp, dp, p_main)
dev.off()







































library(dplyr)
library(tidyr)
for (m in mdm_sig2)
{

tf_long <- metacells_comb %>%
  #rownames_to_column("TF") %>%
  pivot_longer(
    cols = starts_with("Lu952T_"),
    names_to = "cell",
    values_to = "value"
  )
nfkb1_regulon <- tf_long %>%
  filter(TF == m, type == "regulon") %>%
  select(cell, nfkb1_regulon = value)

nfkb1_deviation <- tf_long %>%
  filter(TF == m, type == "deviation") %>%
  select(cell, nfkb1_deviation = value)

trem2_df <- as.data.frame (t(metacells_mat))
trem2_df$cell = rownames(trem2_df)


plot_df <- trem2_df %>%
  left_join(nfkb1_regulon,  by = "cell") %>%
  left_join(nfkb1_deviation, by = "cell")


grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "IL1B"
)
rp <- ggplot(
  plot_df,
  aes(
    x = nfkb1_deviation,
    y = nfkb1_regulon
  )
) +
  geom_point(aes(color = IL1B), alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste(m,"deviation score"),
    y = paste(m,"regulon"),
    title = paste(m,"deviation vs regulon on IL1B (Mono)")
  ) + gtheme

grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "MAF"
)
dp <- ggplot(
  plot_df,
  aes(
    x = nfkb1_deviation,
    y = nfkb1_regulon
  )
) +
  geom_point(aes(color = MAF), alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste(m,"deviation score"),
    y = paste(m,"regulon"),
    title = paste(m,"deviation vs regulon on MAF (IM)")
  ) + gtheme



grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "TREM2"
)

p_main <- ggplot(
  plot_df,
  aes(
    x = nfkb1_deviation,
    y = nfkb1_regulon,
    color = TREM2
  )
) +
  geom_point(alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste(m,"deviation score"),
    y = paste(m,"regulon")#,
    #title = "NFKB1 deviation vs regulon colored by TREM2"
  ) + gtheme_no_rot

pdf (file.path ('Plots',paste(m,'regulon_deviation_TREM2_scatter.pdf')), width=7, height=3.5)
print(wrap_plots (rp, dp, p_main))
dev.off()
}


library(dplyr)
library(tidyr)
for (m in mdm_sig2)
{

tf_long <- metacells_comb %>%
  #rownames_to_column("TF") %>%
  pivot_longer(
    cols = starts_with("Lu952T_"),
    names_to = "cell",
    values_to = "value"
  )
nfkb1_regulon <- tf_long %>%
  filter(TF == m, type == "expression") %>%
  select(cell, nfkb1_regulon = value)

nfkb1_deviation <- tf_long %>%
  filter(TF == m, type == "deviation") %>%
  select(cell, nfkb1_deviation = value)

trem2_df <- as.data.frame (t(metacells_mat))
trem2_df$cell = rownames(trem2_df)


plot_df <- trem2_df %>%
  left_join(nfkb1_regulon,  by = "cell") %>%
  left_join(nfkb1_deviation, by = "cell")


grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "IL1B"
)
rp <- ggplot(
  plot_df,
  aes(
    x = nfkb1_deviation,
    y = nfkb1_regulon
  )
) +
  geom_point(aes(color = IL1B), alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste(m,"deviation score"),
    y = paste(m,"expression"),
    title = paste(m,"deviation vs expression on IL1B (Mono)")
  ) + gtheme

grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "MAF"
)
dp <- ggplot(
  plot_df,
  aes(
    x = nfkb1_deviation,
    y = nfkb1_regulon
  )
) +
  geom_point(aes(color = MAF), alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste(m,"deviation score"),
    y = paste(m,"expression"),
    title = paste(m,"deviation vs expression on MAF (IM)")
  ) + gtheme



grad_scale <- scale_color_gradientn(
  colours = rev(palette_genescore),
  name = "TREM2"
)

p_main <- ggplot(
  plot_df,
  aes(
    x = nfkb1_deviation,
    y = nfkb1_regulon,
    color = TREM2
  )
) +
  geom_point(alpha = 0.6, size = 2) +
  grad_scale +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    x = paste(m,"deviation score"),
    y = paste(m,"expression")#,
    #title = "NFKB1 deviation vs regulon colored by TREM2"
  ) + gtheme_no_rot

pdf (file.path ('Plots',paste(m,'expression_deviation_TREM2_scatter.pdf')), width=7, height=3.5)
print(wrap_plots (rp, dp, p_main))
dev.off()
}




# # Try with PCA on mMat ####
library (ggridges)
library (ggplot2)
library (viridis)
library (tidyr)
#library(hrbrthemes)
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name


# # Filter by RNA expression ####
metaGroupName = 'celltype_simplified'
min_exp = .1
active_TFs = exp_genes (srt, rownames(mMat), min_exp = min_exp, metaGroupName)
mMat = t(scale (mMat[active_TFs, ]))
library(uwot)
library(ggplot2)
library(patchwork)


metacells_mat = metacells@assays$RNA@layers$data
rownames (metacells_mat) = rownames(metacells)
colnames (metacells_mat) = colnames (metacells)
metacells_mat = metacells_mat[active_TFs, ]
#srt = srt[,gsub ('fragments.tsv.gz#','',rownames(archp@cellColData))]
srt_sub = srt@assays$RNA@data[active_TFs,]
srt_sub = srt_sub[rowSums(srt_sub) > 0,]
#srt_sub = srt@assays$RNA@data[active_TFs,]
p <- prcomp(t(srt_sub), center = TRUE, scale. = TRUE)
plot_df <- data.frame(p$x, celltype = srt$shared_cnmf_r_max, sampleID=srt$sampleID)


# df_sub <- plot_df[plot_df$celltype %in% 
#                     c("Mye_Mac_TREM2","Mye_CD14 Mono","Mye_Mac_SPP1",
#                       "Mye_Mac_IM"), ]


df_sub <- plot_df[plot_df$celltype %in% 
                    c("Mono_CD14","TAM_interstitial","TAM_TREM2",
                      "TAM_MARCO","TAM_CXCLs"), ]



#df_sub = plot_df
#-------------------------
# 1. MAIN SCATTER
#-------------------------
p_scatter <- ggplot(df_sub, aes(PC1, PC2, color = celltype)) +
  geom_point(alpha = 0.4, size = 0.1) +
  geom_density_2d(size = 0.4) +
  scale_color_manual(values = palette_myeloid) +
  # scale_x_continuous(limits = x_range, expand = c(0,0)) +
  # scale_y_continuous(limits = y_range, expand = c(0,0)) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )


p_density_x <- ggplot(df_sub, aes(PC1, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
 scale_fill_manual(values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) + gtheme_no_rot

p_density_y <- ggplot(df_sub, aes(PC2, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) + gtheme_no_rot

pdf (file.path ('Plots','scrna_TF_expression_pca.pdf'), width = 5,height=6)
p_density_x + plot_spacer() + p_scatter + p_density_y + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4),
    guides = "collect"      # <- collect legends into one
  ) &
  theme(
    legend.position = "bottom",  # or "top"
    legend.box = "vertical"
  )

dev.off()

a = p$rotation[,1]
a[order(-a)]




all (rownames (mMat) == rownames (archp@cellColData))
p <- prcomp(mMat, center = TRUE, scale. = TRUE)
p <- prcomp(mMat, center = TRUE, scale. = TRUE)
plot_df <- data.frame(p$x, celltype = archp$shared_cnmf_r_max)


# df_sub <- plot_df[plot_df$celltype %in% 
#                     c("Mye_Mac_TREM2","Mye_CD14 Mono","Mye_Mac_SPP1",
#                       "Mye_Mac_IM"), ]


df_sub <- plot_df[plot_df$celltype %in% 
                    c("Mono_CD14","TAM_interstitial","TAM_TREM2",
                      "TAM_MARCO","TAM_CXCLs"), ]



#df_sub = plot_df
#-------------------------
# 1. MAIN SCATTER
#-------------------------
p_scatter <- ggplot(df_sub, aes(PC1, PC2, color = celltype)) +
  geom_point(alpha = 0.4, size = 0.1) +
  geom_density_2d(size = 0.4) +
  scale_color_manual(values = palette_myeloid) +
  # scale_x_continuous(limits = x_range, expand = c(0,0)) +
  # scale_y_continuous(limits = y_range, expand = c(0,0)) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )


p_density_x <- ggplot(df_sub, aes(PC1, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) + gtheme_no_rot

p_density_y <- ggplot(df_sub, aes(PC2, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) + gtheme_no_rot

pdf (file.path ('Plots','dev_momac_pca.pdf'), width = 5,height=6)
p_density_x + plot_spacer() + p_scatter + p_density_y + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4),
    guides = "collect"      # <- collect legends into one
  ) &
  theme(
    legend.position = "bottom",  # or "top"
    legend.box = "vertical"
  )

dev.off()

a = p$rotation[,1]
a[order(-a)]





# Check activity of NFKB1 regulon ####
srt = srt[, srt$sampleID %in% c('Lu952T','Lu979T')]
scenic_name = 'monomac_programs'
auc_mtx = read.csv(file.path(paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/AS_human_lung_scatac/analysis/NTP_multiome/SCENIC/vg_5000_mw_tss500bp/',scenic_name), 'auc_mtx.csv'), header=T)
rownames (auc_mtx) = auc_mtx[,1]
rownames(auc_mtx) = gsub ('\\.','-', rownames(auc_mtx))
auc_mtx = auc_mtx[,-1]
colnames (auc_mtx) = sub ('\\.\\.\\.','', colnames(auc_mtx))

#auc_mtx = auc_mtx['NFKB1',]

auc_mtx = aggregate (auc_mtx, by = list(
  celltype = srt$shared_cnmf_r_max[match(rownames(auc_mtx), colnames(srt))],
  sample = srt$sampleID[match(rownames(auc_mtx), colnames(srt))]), FUN = mean)
auc_mtx = gather (auc_mtx, TF, score, 3:(ncol(auc_mtx)))

gene = 'FOSL2'
cell_subsets_order = c("TAM_interstitial","TAM_MARCO","TAM_TREM2","TAM_CXCLs","Mono_CD14")
auc_mtx_sub = auc_mtx[auc_mtx$TF == gene,]
auc_mtx_sub = auc_mtx_sub[auc_mtx_sub$celltype %in% cell_subsets_order, ]
auc_mtx_sub$celltype = factor (auc_mtx_sub$celltype, levels = rev(cell_subsets_order))

bp = ggplot (auc_mtx_sub, 
  aes(x = celltype, y = score, fill = celltype)) + 
  #vlp +
  geom_boxplot (
    linewidth = .2,
    width=.5,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.9
     ) +
  scale_fill_manual (values = palette_myeloid) + 
  gtheme

pdf (file.path ('Plots',paste0(gene,'_regulon_boxplots_celltypes.pdf')), width = 4, height=4)
bp
dev.off()

# Check activity of NFKB1 expression ####
pb = as.data.frame (log2 (AverageExpression (srt, feature = gene, group.by = c('shared_cnmf_r_max','sampleID'))[[1]]+1))
pb$TF = gene
pb = gather (pb, celltype, score, 1:((ncol(pb)-1)))
pb$sample = sapply (pb$celltype, function(x) unlist(strsplit(x, '_'))[2])
pb$celltype = sapply (pb$celltype, function(x) unlist(strsplit(x, '_'))[1])
pb$celltype = gsub ('-','_',pb$celltype)
pb = pb[pb$celltype %in% cell_subsets_order,]
pb$celltype = factor (pb$celltype, levels = rev(cell_subsets_order))

bp = ggplot (pb, 
  aes(x = celltype, y = score, fill = celltype)) + 
  #vlp +
  geom_boxplot (
    linewidth = .2,
    width=.5,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.9
     ) +
  scale_fill_manual (values = palette_myeloid) + 
  gtheme

pdf (file.path ('Plots',paste0(gene, '_expression_boxplots_celltypes.pdf')), width = 4, height=4)
bp
dev.off()


# Check activity of NFKB1 deviation ####
cell_subsets_order = c("TAM_interstitial","TAM_MARCO","TAM_TREM2","TAM_CXCLs","Mono_CD14")
if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = mMat[rownames (mMat) %in% gene,,drop=F]
mMat = t(mMat)
mMat = as.data.frame (mMat)
mMat = aggregate (mMat, by = list(celltype = archp$shared_cnmf_r_max,
  sample = archp$Sample), mean)
mMat = mMat[mMat$celltype %in% cell_subsets_order,]
mMat$celltype = factor (mMat$celltype, levels = rev(cell_subsets_order))
bp = ggplot (mMat, 
  aes_string(x = 'celltype', y = gene, fill = 'celltype')) + 
  #vlp +
  geom_boxplot (
    linewidth = .2,
    width=.5,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.9
     ) +
  scale_fill_manual (values = palette_myeloid) + 
  gtheme

pdf (file.path ('Plots',paste0(gene,'_deviation_boxplots_celltypes.pdf')), width = 4, height=4)
bp
dev.off()


# Check activity of NFKB1 genescore ####
if (!exists('gsSE')) gsSE = fetch_mat (archp, 'GeneScore')
gsMat = scale(assays (gsSE)[[1]])
rownames (gsMat) = rowData (gsSE)$name
#gsMat = scale(gsMat[rownames (gsMat) %in% mdm_sig,])
gsMat = gsMat[rownames (gsMat) %in% gene,,drop=F]
gsMat = t(gsMat)
gsMat = as.data.frame (gsMat)
gsMat = aggregate (gsMat, by = list(celltype = archp$shared_cnmf_r_max,
  sample = archp$Sample), mean)
gsMat = gsMat[gsMat$celltype %in% cell_subsets_order,]
gsMat$celltype = factor (gsMat$celltype, levels = rev(cell_subsets_order))

bp = ggplot (gsMat, 
  aes_string(x = 'celltype', y = gene, fill = 'celltype')) + 
  #vlp +
  geom_boxplot (
    linewidth = .2,
    width=.5,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.9
     ) +
  scale_fill_manual (values = palette_myeloid) + 
  gtheme

pdf (file.path ('Plots',paste0(gene,'_genescore_boxplots_celltypes.pdf')), width = 4, height=4)
bp
dev.off()




