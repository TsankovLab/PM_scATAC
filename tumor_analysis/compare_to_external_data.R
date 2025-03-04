
# ### check correlation of sarcomatoid score to external celltypes ####
# # Read ENCODE collection for peak annotation
# H3K4me3_primary_filepath = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/ENCODE/H3K4me3_ChIP-seq/primary_cells/'
# encode_H3K4me3_primary = readRDS (paste0(H3K4me3_primary_filepath, 'H3K4me3_primary_cells_bed_list.rds'))

# archp = addBgdPeaks (archp, force= TRUE)
# archp = addPeakAnnotations (ArchRProj = archp, force=T,
#      regions = encode_H3K4me3_primary, name = "ENCODE_H3K4me3")

# archp = addDeviationsMatrix (
#   ArchRProj = archp, 
#   peakAnnotation = "ENCODE_H3K4me3",
#   force = TRUE
# )
# # Compute correlation of sarcomatoid cNMF external celltypes ####
# sams = unique(archp$Sample3)
# sams = sams[!sams %in% c('normal_pleura','P3','P13','P11_HOX')] # Remove normal and low number samples and outliers
# cnmf_mat = as.matrix(archp@cellColData[,grep ('cNMF', colnames(archp@cellColData))])
# cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp$Sample3 == x, ])))
# names (cnmf_mat) = sams

# # # Get ENCODE matrix ####
# if (!exists('enSE')) enSE = fetch_mat (archp, 'ENCODE_H3K4me3')
# all (colnames(enSE) == rownames(archp))
# encode_Mat = assays (enSE)[[1]]
# rownames (encode_Mat) = rownames(rowData (enSE))

# # # Get external matrix ####
# # if (!exists('fSE')) fSE = fetch_mat (archp, 'scATAC_datasets')
# # all (colnames(fSE) == rownames(archp))
# # fetal_Mat = assays (fSE)[[1]]
# # rownames (fetal_Mat) = rownames(rowData (fSE))

# ext_mat = encode_Mat
# remove_immune = c(
#   rownames (ext_mat)[grep ('CD4',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('lymph',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('CD8',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('T-',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('\\.T\\.cell',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('killer',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('\\.B\\.cell',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('NK',rownames(ext_mat), ignore.case =F)],
#   rownames (ext_mat)[grep ('Treg',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('Naive.B',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('Naive.T',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('Mast',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('Memory.B',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('Mac',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('Memory.B',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('Plasma',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('Monocyte',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('Myelo',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('DC',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('CD5',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('myeloid',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('B cell',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('neutrophil',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('T cell',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('ILC',rownames(ext_mat), ignore.case =T)],
#   rownames (ext_mat)[grep ('mononuclear',rownames(ext_mat), ignore.case =F)],
#   rownames (ext_mat)[grep ('Eryth',rownames(ext_mat), ignore.case =F)])

# ext_mat = ext_mat[!rownames(ext_mat)%in% remove_immune,]



# ext_mat = lapply (sams, function(x) scale(ext_mat[,archp$Sample3 == x]))
# names (ext_mat) = sams

# cnmf_ext_l = list()
# for (sam in sams)
#   {
#   cnmf_ext_l[[sam]] = cor (t(cnmf_mat[[sam]]), t(ext_mat[[sam]]), method = 'spearman')
#   }

# cnmf_ext_array <- simplify2array(cnmf_ext_l)
# #any(lapply(corTF_array, function(x) any(is.na(x))))
# # Take element-wise median
# median_matrix <- apply(cnmf_ext_array, c(1, 2), median)

# cor_cnmf_ext = Heatmap (t(median_matrix[rownames(median_matrix) != 'sarcomatoid.cNMF9',]),
#   row_names_gp = gpar(fontsize = 8),
#   #rect_gp = gpar(type = "none"),
#   column_names_gp = gpar(fontsize = 7, fontface='italic'),
#   #col = palette_deviation_cor_fun
#   )

# pdf (paste0 ('Plots/selected_fetal_cnmf_corr_heatmaps.pdf'), width = 7,height=5)
# cor_cnmf_ext
# dev.off()

# # Also plot boxplots ordered by correlation to sarcomatoid score ####
# sarc_module = 'cNMF20'
# top_ext_ct = median_matrix[sarc_module,]
# top_ext_ct = top_ext_ct[order (-top_ext_ct)]
# top_ext_ct = names(c(head (top_ext_ct, 20)))#,tail (top_ext_ct, 20)))
# sarc_tf = lapply (sams, function(sam) cor (t(ext_mat[[sam]][top_ext_ct,]), t(cnmf_mat[[sam]][cnmf_module,,drop=F]), method = 'spearman'))
# sarc_tf_df = do.call (cbind,sarc_tf)

# sarc_tf_df = as.data.frame (sarc_tf_df)
# #sarc_tf_df2$TF = rownames(sarc_tf_df2)
# colnames(sarc_tf_df) = sams
# sarc_tf_df$ENCODE_celltype = factor (top_ext_ct, levels = rev(top_ext_ct))
# sarc_tf_df = gather (sarc_tf_df, scS, score, 1:(ncol(sarc_tf_df)-1))
# bp = ggplot (sarc_tf_df, aes (x = score, y = ENCODE_celltype)) + 
#   geom_boxplot (alpha=.8, outlier.shape = NA, fill = 'darkred', color='grey22') + 
#   geom_jitter(width = 0.2, alpha = 0.4, color = "grey22", size=.5) + 
#   gtheme_no_rot + 
#   geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1)

# pdf (paste0 ('Plots/sarcomatoid_score_ENCODE_boxplots2.pdf'), width = 6,height=5)
# bp
# dev.off()

# Add epithelioid gene signature 
epit = list(epit = 'PDZK1IP1
EFNA5
RAET1E
TGM1
CDON
ARHGAP44
EPHB6
MEIS2
NRG4
PARD6B
PLLP
RERG
ANXA9
SOX6
PROCR
C1S
BCO2
MSLN
SLC9A3R1
KLHL31
AGFG2
GALNT9
ITLN1
GAL3ST2
CFB
SGPP2
CFI
KLK11
NFIA
SELENBP1
CARNS1
FLRT3
COBL
GHR
TJP3
PRR15
TPD52L1
CLDN15
C1RL')
epit = strsplit(epit[[1]],'\\\n')

  archp = addModuleScore (
      ArchRProj = archp,
      useMatrix = 'GeneScoreMatrix',
      name = '',
      features = epit,
      nBin = 25,
      nBgd = 100,
      seed = 1,
      threads = getArchRThreads(),
      logFile = createLogFile("addModuleScore")
    )

# Define sarc - epithelioid module ####
archp$cNMF_sarc_epit = archp$cNMF20 - archp$cNMF6
archp$cNMF_sarc_epit = as.numeric(archp$cNMF_sarc_epit)

### Import H3K27Ac data from primary cells from ENCODE ####
encode_h3k27ac_path = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/ENCODE/H3K27Ac'
encode_h3k27ac_files = file.path(encode_h3k27ac_path, list.files(encode_h3k27ac_path,pattern = '.bed'))
if (!file.exists(file.path (encode_h3k27ac_path,'encode_h3k27ac.rds')))
  {
  encode_h3k27ac_gr = lapply (encode_h3k27ac_files, function(x) 
    {
    tmp = read.table (x)
    colnames(tmp) = c('chr','start','end')
    makeGRangesFromDataFrame (tmp)
    })
  names (encode_h3k27ac_gr) = sub('.bed','',list.files(encode_h3k27ac_path,pattern = '.bed'))
  metadata = read.table (file.path(encode_h3k27ac_path,'metadata.tsv'), sep='\t', header=T)
  metadata = metadata[metadata$File.assembly == 'GRCh38',]
  encode_h3k27ac_gr = encode_h3k27ac_gr[names (encode_h3k27ac_gr) %in% metadata$File.accession]
  names (encode_h3k27ac_gr) = paste0(metadata$Biosample.type,'|', metadata$Biosample.term.name)[match(names(encode_h3k27ac_gr), metadata$File.accession)]
  encode_h3k27ac_gr = encode_h3k27ac_gr[!duplicated (names(encode_h3k27ac_gr))]
  saveRDS (encode_h3k27ac_gr, file.path (encode_h3k27ac_path,'encode_h3k27ac.rds'))
  } else {
  encode_h3k27ac_gr = readRDS (file.path(encode_h3k27ac_path,'encode_h3k27ac.rds'))
  }

archp = addBgdPeaks (archp, force= FALSE)
archp = addPeakAnnotations (ArchRProj = archp, force=T,
     regions = encode_h3k27ac_gr, name = "ENCODE_H3K27Ac")

archp = addDeviationsMatrix (
  ArchRProj = archp, 
  peakAnnotation = "ENCODE_H3K27Ac",
  force = FALSE
  )

saveArchRProject (archp)


# Compute correlation of sarcomatoid cNMF external celltypes ####
sams = unique(archp$Sample3)
sams = sams[!sams %in% c('normal1','P3','P13','P11_HOX')] # Remove normal and low number samples and outliers
cnmf_mat = archp@cellColData[,grep ('cNMF', colnames(archp@cellColData))]
cnmf_mat$cNMF_sarc_epit = as.numeric (cnmf_mat$cNMF_sarc_epit)
cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp$Sample3 == x, ])))
names (cnmf_mat) = sams

# # Get ENCODE matrix ####
if (!exists('enSE')) enSE = fetch_mat (archp, 'ENCODE_H3K27Ac')
all (colnames(enSE) == rownames(archp))
encode_Mat = assays (enSE)[[1]]
rownames (encode_Mat) = rownames(rowData (enSE))

# # Get external matrix ####
# if (!exists('fSE')) fSE = fetch_mat (archp, 'scATAC_datasets')
# all (colnames(fSE) == rownames(archp))
# fetal_Mat = assays (fSE)[[1]]
# rownames (fetal_Mat) = rownames(rowData (fSE))

ext_mat = encode_Mat
remove_immune = c(
  rownames (ext_mat)[grep ('CD4',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('lymph',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('CD8',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('T-',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('\\.T\\.cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('killer',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('\\.B\\.cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('NK',rownames(ext_mat), ignore.case =F)],
  rownames (ext_mat)[grep ('Treg',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Naive.B',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Naive.T',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Mast',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Memory.B',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Mac',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Memory.B',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Plasma',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Monocyte',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Myelo',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('DC',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('CD5',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('myeloid',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('B cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('neutrophil',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('T cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('ILC',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('mononuclear',rownames(ext_mat), ignore.case =F)],
  rownames (ext_mat)[grep ('Eryth',rownames(ext_mat), ignore.case =F)])

ext_mat = ext_mat[!rownames(ext_mat)%in% remove_immune,]



ext_mat = lapply (sams, function(x) scale(ext_mat[,archp$Sample3 == x]))
names (ext_mat) = sams

cnmf_ext_l = list()
for (sam in sams)
  {
  cnmf_ext_l[[sam]] = cor (t(cnmf_mat[[sam]]), t(ext_mat[[sam]]), method = 'spearman')
  }

cnmf_ext_array <- simplify2array(cnmf_ext_l)
#any(lapply(corTF_array, function(x) any(is.na(x))))
# Take element-wise median
median_matrix <- apply(cnmf_ext_array, c(1, 2), median)
median_matrix = median_matrix[, sapply(colnames(median_matrix), function(x) nchar(x) < 100)]

cor_cnmf_ext = Heatmap (t(median_matrix[rownames(median_matrix) != 'sarcomatoid.cNMF9',]),
  row_names_gp = gpar(fontsize = 8),
  #rect_gp = gpar(type = "none"),
  column_names_gp = gpar(fontsize = 7, fontface='italic'),
  #col = palette_deviation_cor_fun
  )

pdf (paste0 ('Plots/selected_ENCODE_H3K27_cnmf_corr_heatmaps.pdf'), width = 7,height=25)
cor_cnmf_ext
dev.off()

# Also plot boxplots ordered by correlation to sarcomatoid score ####
sarc_module = 'cNMF20'
top_ext_ct = median_matrix[sarc_module,]
top_ext_ct = top_ext_ct[order (-top_ext_ct)]
top_ext_ct2 = sapply (names(top_ext_ct), function(x) unlist(strsplit (x, '\\|'))[1])
top_ext_ct = names (unlist(unname(lapply (split(top_ext_ct2, top_ext_ct2),function(x) c(head(x,5),tail(x,5))))))
#top_ext_ct = names(c(head (top_ext_ct, Inf)))#,tail (top_ext_ct, 20)))
sarc_tf = lapply (sams, function(sam) cor (t(ext_mat[[sam]][top_ext_ct,]), t(cnmf_mat[[sam]][sarc_module,,drop=F]), method = 'spearman'))
sarc_tf_df = do.call (cbind,sarc_tf)

sarc_tf_df = as.data.frame (sarc_tf_df)
#sarc_tf_df2$TF = rownames(sarc_tf_df2)
colnames(sarc_tf_df) = sams
sarc_tf_df$ENCODE_celltype = factor (top_ext_ct, levels = rev(top_ext_ct))
sarc_tf_df = gather (sarc_tf_df, scS, score, 1:(ncol(sarc_tf_df)-1))
sarc_tf_df$biotype = sapply (as.character(sarc_tf_df$ENCODE_celltype), function(x) unlist(strsplit (x, '\\|'))[1])
#sarc_tf_df = do.call (rbind, lapply(split(sarc_tf_df, sarc_tf_df$biotype), function(x) head (x,10)))
bp = lapply (unique (sarc_tf_df$biotype), function(x) ggplot (sarc_tf_df[sarc_tf_df$biotype == x,], aes (x = score, y = ENCODE_celltype)) + 
  geom_jitter(width = 0.2, alpha = 0.4, color = "grey22", size=.5) + 
  gtheme + 
  geom_boxplot (aes (fill = ENCODE_celltype),
    linewidth = .2,
    width=.6,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=1, alpha=0.4
     ) +
  scale_fill_manual (values = rev(paletteer::paletteer_c("grDevices::Purple-Blue",length(unique(sarc_tf_df[sarc_tf_df$biotype == x,]$ENCODE_celltype))))) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = .4) + NoLegend() + coord_flip())
  

pdf (paste0 ('Plots/sarcomatoid_score_ENCODE_H3K27_boxplots2.pdf'), width = 11,height=5)
wrap_plots (bp,ncol=4)
dev.off()



# Read in peak files from scATAC studies ####
projects = c('Tsankov_lung')
projects_peaks = lapply (seq_along(projects), function(x) {
  bed_files = list.files (file.path('..','..','tumor_compartment','all_tissues_ArchR',projects[x],'PeakCalls'), pattern = '.rds')
  grlist = lapply (seq_along(bed_files), 
    function(y) readRDS (file.path('..','..','tumor_compartment','all_tissues_ArchR',projects[x],'PeakCalls',bed_files[y])))
names (grlist) = paste0(projects[x], '_', sapply (bed_files, function(z) unlist(strsplit (z, '-'))[1]))
grlist
})
projects_peaks = unlist (projects_peaks, recursive=F)
# ### chromVAR analysis
archp = addBgdPeaks (archp, force= FALSE)
archp = addPeakAnnotations (ArchRProj = archp, 
     regions = projects_peaks, name = "Normal_Lung")

archp = addDeviationsMatrix (
  ArchRProj = archp, 
  peakAnnotation = "Normal_Lung",
  force = FALSE
)


# Compute correlation of sarcomatoid cNMF external celltypes ####
sams = unique(archp$Sample3)
sams = sams[!sams %in% c('normal1','P3','P13','P11_HOX')] # Remove normal and low number samples and outliers
cnmf_mat = archp@cellColData[,grep ('cNMF', colnames(archp@cellColData))]
cnmf_mat$cNMF_sarc_epit = as.numeric (cnmf_mat$cNMF_sarc_epit)
cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp$Sample3 == x, ])))
names (cnmf_mat) = sams

# # Get Normal Lung matrix ####
if (!exists('nlSE')) nlSE = fetch_mat (archp, 'Normal_Lung')
all (colnames(nlSE) == rownames(archp))
normal_Mat = assays (nlSE)[[1]]
rownames (normal_Mat) = rownames(rowData (nlSE))

# # Get external matrix ####
# if (!exists('fSE')) fSE = fetch_mat (archp, 'scATAC_datasets')
# all (colnames(fSE) == rownames(archp))
# fetal_Mat = assays (fSE)[[1]]
# rownames (fetal_Mat) = rownames(rowData (fSE))

ext_mat = normal_Mat
remove_immune = c(
  rownames (ext_mat)[grep ('CD4',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('lymph',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('CD8',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('T-',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('\\.T\\.cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('killer',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('\\.B\\.cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('NK',rownames(ext_mat), ignore.case =F)],
  rownames (ext_mat)[grep ('Treg',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Naive.B',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Naive.T',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Mast',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Memory.B',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Mac',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Memory.B',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Plasma',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Monocyte',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('Myelo',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('DC',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('CD5',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('myeloid',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('B cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('neutrophil',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('T cell',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('ILC',rownames(ext_mat), ignore.case =T)],
  rownames (ext_mat)[grep ('mononuclear',rownames(ext_mat), ignore.case =F)],
  rownames (ext_mat)[grep ('Eryth',rownames(ext_mat), ignore.case =F)])

ext_mat = ext_mat[!rownames(ext_mat)%in% remove_immune,]



ext_mat = lapply (sams, function(x) scale(ext_mat[,archp$Sample3 == x]))
names (ext_mat) = sams

cnmf_ext_l = list()
for (sam in sams)
  {
  cnmf_ext_l[[sam]] = cor (t(cnmf_mat[[sam]]), t(ext_mat[[sam]]), method = 'spearman')
  }

cnmf_ext_array <- simplify2array(cnmf_ext_l)
#any(lapply(corTF_array, function(x) any(is.na(x))))
# Take element-wise median
median_matrix <- apply(cnmf_ext_array, c(1, 2), median)
median_matrix = median_matrix[, sapply(colnames(median_matrix), function(x) nchar(x) < 100)]

cor_cnmf_ext = Heatmap (t(median_matrix[rownames(median_matrix) != 'sarcomatoid.cNMF9',]),
  row_names_gp = gpar(fontsize = 8),
  #rect_gp = gpar(type = "none"),
  column_names_gp = gpar(fontsize = 7, fontface='italic'),
  #col = palette_deviation_cor_fun
  )

pdf (paste0 ('Plots/selected_normal_lung_cnmf_corr_heatmaps.pdf'), width = 7,height=25)
cor_cnmf_ext
dev.off()

# Also plot boxplots ordered by correlation to sarcomatoid score ####
sarc_module = 'cNMF20'
sarc_module = 'cNMF_sarc_epit'
top_ext_ct = median_matrix[sarc_module,]
top_ext_ct = top_ext_ct[order (-top_ext_ct)]
#top_ext_ct2 = sapply (names(top_ext_ct), function(x) unlist(strsplit (x, '\\|'))[1])
#top_ext_ct = c(head(top_ext_ct,5),tail(top_ext_ct,5))
#top_ext_ct = names(c(head (top_ext_ct, Inf)))#,tail (top_ext_ct, 20)))
sarc_tf = lapply (sams, function(sam) cor (t(ext_mat[[sam]][names(top_ext_ct),]), t(cnmf_mat[[sam]][sarc_module,,drop=F]), method = 'spearman'))
sarc_tf_df = do.call (cbind,sarc_tf)

sarc_tf_df = as.data.frame (sarc_tf_df)
#sarc_tf_df2$TF = rownames(sarc_tf_df2)
colnames(sarc_tf_df) = sams
names (top_ext_ct) = sub ('Tsankov_lung_','',names(top_ext_ct))
sarc_tf_df$normal_celltype = factor (names(top_ext_ct), levels = rev(unique(rev(names(top_ext_ct)))))
sarc_tf_df = gather (sarc_tf_df, scS, score, 1:(ncol(sarc_tf_df)-1))
#sarc_tf_df$biotype = sapply (as.character(sarc_tf_df$ENCODE_celltype), function(x) unlist(strsplit (x, '\\|'))[1])
#sarc_tf_df = do.call (rbind, lapply(split(sarc_tf_df, sarc_tf_df$biotype), function(x) head (x,10)))
palette_expression_disc = paletteer::paletteer_c("grDevices::Purple-Blue",length(unique(sarc_tf_df$normal_celltype)))
head (sarc_tf_df)

bp = ggplot (sarc_tf_df, aes (x = score, y = normal_celltype)) + 
  #geom_boxplot (alpha=.8, outlier.shape = NA, fill = 'darkred', color='grey22') + 
  geom_jitter(width = 0.2, alpha = 0.4, color = "grey22", size=.5) + 
  geom_boxplot (aes (fill = normal_celltype),
    linewidth = .2,
    width=.6,
    outlier.alpha = 0.2,
    outlier.size = 0,
     size=1, alpha=0.6
     ) +
  scale_fill_manual(values = palette_expression_disc) +
  gtheme + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = .4)+
  coord_flip()
  

pdf (paste0 ('Plots/sarcomatoid_score_normal_lung_boxplots4.pdf'), width = 6,height=3)
bp
dev.off()








