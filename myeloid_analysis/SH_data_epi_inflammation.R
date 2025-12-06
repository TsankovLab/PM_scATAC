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





# #### Identify TF regulators correlated to each scRNA cnmf #####
# if (!exists('mSE')) mSE = fetch_mat (archp, 'Motif')
# mMat = assays (mSE)[[1]]
# rownames (mMat) = rowData (mSE)$name

# # Filter by RNA expression ####
# metaGroupName = 'celltype_SM'
# active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)
# mMat = t(scale (mMat[active_TFs, ]))
# if (!file.exists ('mMat_scaled_active.rds')) saveRDS (mMat, 'mMat_scaled_active.rds')

# #mMat = scale(as.matrix(mMat))#[selected_TF,])
# mMat = readRDS ('mMat_scaled_active.rds')
# mMat_cor = cor (as.matrix(mMat), method = 'spearman')

# set.seed(1234)
# centers=2
# km = kmeans (mMat_cor, centers=centers)
# if (!file.exists ('TF_activity_modules.rds')) saveRDS (km, 'TF_activity_modules.rds')

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
# ha2 = rowAnnotation (foo = anno_mark(at = match(AP1,colnames(mMat_cor)), 
#     labels = AP1, labels_gp = gpar(fontsize = 7)))

# pdf (file.path ('Plots','TF_modules_heatmap3.pdf'), width = 4,height=3)
# cor_mMat_hm = draw (Heatmap (mMat_cor,# row_km=15,
#   right_annotation = ha2,
#   #left_annotation = ha,
#   #rect_gp = gpar(type = "none"),
#   clustering_distance_rows='euclidean' ,
#   clustering_distance_columns = 'euclidean', 
#   col=palette_deviation_cor_fun, 
#   row_split = km$cluster,
#   column_split = km$cluster,
#   #row_km=2, 
#   #column_km=2,
# #  right_annotation = ha,
#   border=T,
# #   ,
#   row_names_gp = gpar(fontsize = 0), 
#   column_names_gp = gpar(fontsize = 0)
# # cell_fun = function(j, i, x, y, w, h, fill) {# THIS DOESNT WORK NEED TO USE LAYER_FUN
# #         if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
# #             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
# #         }}
#   ))
#   # ,
#   # cell_fun = function(j, i, x, y, w, h, fill) {
#   #       if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#   #           grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
# #        }}))
# dev.off()

# pdf (file.path ('Plots','TF_modules_heatmap3.pdf'), width = 4.6, height=3)
# cor_mMat_hm
# dev.off()

# #km = readRDS ("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/TF_activity_modules.rds")

# #AP1 = names (km$cluster[km$cluster == 2])
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

# all (rownames(archp@cellColData) == rownames(mMat))
# archp$mod_2 = rowMeans (mMat[,AP1])

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

pdf (file.path ('Plots','scatac_cnmf_inflammation_module_ridge_plots.pdf'), width = 5,height=6)
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

auc_mtx = auc_mtx[colnames(srt),]
p <- prcomp(auc_mtx, center = TRUE, scale. = TRUE)
plot_df <- data.frame(p$x, celltype = srt$shared_cnmf_r_max)


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
# 0. RUN UMAP
#-------------------------
set.seed(123)  # for reproducibility
library(uwot)
umap_res <- umap(
  t(mMat),
  n_neighbors = 30,
  min_dist = 0.3,
  metric = "cosine",
  scale = TRUE
)

plot_df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  celltype = archp$shared_cnmf_r_max
)

df_sub <- plot_df[plot_df$celltype %in%
                    c("Mono_CD14","TAM_interstitial","TAM_TREM2",
                      "TAM_MARCO","TAM_CXCLs"), ]


#-------------------------
# 1. MAIN SCATTER
#-------------------------
p_scatter <- ggplot(df_sub, aes(UMAP1, UMAP2, color = celltype)) +
  geom_point(alpha = 0.4, size = 0.1) +
  geom_density_2d(size = 0.4) +
  scale_color_manual(values = palette_myeloid) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )


#-------------------------
# 2. TOP DENSITY (UMAP1)
#-------------------------
p_density_x <- ggplot(df_sub, aes(UMAP1, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  gtheme_no_rot


#-------------------------
# 3. RIGHT DENSITY (UMAP2)
#-------------------------
p_density_y <- ggplot(df_sub, aes(UMAP2, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  gtheme_no_rot


#-------------------------
# 4. SAVE PDF
#-------------------------
pdf(file.path("Plots","scatac_cnmf_inflammation_module_umap_ridge_plots.pdf"),
    width = 5, height = 5)

p_density_x + plot_spacer() + p_scatter + p_density_y +
  plot_layout(
    ncol = 2, nrow = 2,
    widths  = c(4, 1),
    heights = c(1, 4),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

dev.off()



umap_df <- data.frame(UMAP1 = umap_res[,1],
                      UMAP2 = umap_res[,2])

cors_UMAP1 <- apply(t(mMat), 2, function(x) cor(x, umap_df$UMAP1))
cors_UMAP2 <- apply(t(mMat), 2, function(x) cor(x, umap_df$UMAP2))

loadings_like <- data.frame(peak = colnames(mMat),
                            cor_UMAP1 = cors_UMAP1,
                            cor_UMAP2 = cors_UMAP2)



#-------------------------
# 0. RUN UMAP
#-------------------------
set.seed(123)  # for reproducibility
library(uwot)
umap_res <- umap(
  auc_mtx,
  n_neighbors = 30,
  min_dist = 0.3,
  metric = "cosine",
  scale = TRUE
)

plot_df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  celltype = srt$shared_cnmf_r_max
)

df_sub <- plot_df[plot_df$celltype %in%
                    c("Mono_CD14","TAM_interstitial","TAM_TREM2",
                      "TAM_MARCO","TAM_CXCLs"), ]


#-------------------------
# 1. MAIN SCATTER
#-------------------------
p_scatter <- ggplot(df_sub, aes(UMAP1, UMAP2, color = celltype)) +
  geom_point(alpha = 0.4, size = 0.1) +
  geom_density_2d(size = 0.4) +
  scale_color_manual(values = palette_myeloid) +
  gtheme_no_rot +
  theme(
    legend.position = "right",
    plot.margin = margin(0,0,0,0)
  )


#-------------------------
# 2. TOP DENSITY (UMAP1)
#-------------------------
p_density_x <- ggplot(df_sub, aes(UMAP1, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  gtheme_no_rot


#-------------------------
# 3. RIGHT DENSITY (UMAP2)
#-------------------------
p_density_y <- ggplot(df_sub, aes(UMAP2, fill = celltype)) +
  geom_density(alpha = 0.3, color = NA) +
  scale_fill_manual(values = palette_myeloid) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  gtheme_no_rot


#-------------------------
# 4. SAVE PDF
#-------------------------
pdf(file.path("Plots","scrna_cnmf_inflammation_module_umap_ridge_plots.pdf"),
    width = 5, height = 5)

p_density_x + plot_spacer() + p_scatter + p_density_y +
  plot_layout(
    ncol = 2, nrow = 2,
    widths  = c(4, 1),
    heights = c(1, 4),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

dev.off()




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
