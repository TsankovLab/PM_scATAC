
### check correlation of sarcomatoid score to external celltypes ####
# Read ENCODE collection for peak annotation
H3K4me3_primary_filepath = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/ENCODE/H3K4me3_ChIP-seq/primary_cells/'
encode_H3K4me3_primary = readRDS (paste0(H3K4me3_primary_filepath, 'H3K4me3_primary_cells_bed_list.rds'))

archp = addBgdPeaks (archp, force= TRUE)
archp = addPeakAnnotations (ArchRProj = archp, force=T,
     regions = encode_H3K4me3_primary, name = "ENCODE_H3K4me3")

archp = addDeviationsMatrix (
  ArchRProj = archp, 
  peakAnnotation = "ENCODE_H3K4me3",
  force = TRUE
)
# Compute correlation of sarcomatoid cNMF external celltypes ####
sams = unique(archp$Sample3)
sams = sams[!sams %in% c('normal_pleura','P3','P13','P11_HOX')] # Remove normal and low number samples and outliers
cnmf_mat = as.matrix(archp@cellColData[,grep ('cNMF', colnames(archp@cellColData))])
cnmf_mat = lapply (sams, function(x) scale(t(cnmf_mat[archp$Sample3 == x, ])))
names (cnmf_mat) = sams

# # Get ENCODE matrix ####
if (!exists('enSE')) enSE = fetch_mat (archp, 'ENCODE_H3K4me3')
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

cor_cnmf_ext = Heatmap (t(median_matrix[rownames(median_matrix) != 'sarcomatoid.cNMF9',]),
  row_names_gp = gpar(fontsize = 8),
  #rect_gp = gpar(type = "none"),
  column_names_gp = gpar(fontsize = 7, fontface='italic'),
  #col = palette_deviation_cor_fun
  )

pdf (paste0 ('Plots/selected_fetal_cnmf_corr_heatmaps.pdf'), width = 7,height=5)
cor_cnmf_ext
dev.off()

# Also plot boxplots ordered by correlation to sarcomatoid score ####
sarc_module = 'cNMF20'
top_ext_ct = median_matrix[sarc_module,]
top_ext_ct = top_ext_ct[order (-top_ext_ct)]
top_ext_ct = names(c(head (top_ext_ct, 20)))#,tail (top_ext_ct, 20)))
sarc_tf = lapply (sams, function(sam) cor (t(ext_mat[[sam]][top_ext_ct,]), t(cnmf_mat[[sam]][cnmf_module,,drop=F]), method = 'spearman'))
sarc_tf_df = do.call (cbind,sarc_tf)

sarc_tf_df = as.data.frame (sarc_tf_df)
#sarc_tf_df2$TF = rownames(sarc_tf_df2)
colnames(sarc_tf_df) = sams
sarc_tf_df$ENCODE_celltype = factor (top_ext_ct, levels = rev(top_ext_ct))
sarc_tf_df = gather (sarc_tf_df, scS, score, 1:(ncol(sarc_tf_df)-1))
bp = ggplot (sarc_tf_df, aes (x = score, y = ENCODE_celltype)) + 
  geom_boxplot (alpha=.8, outlier.shape = NA, fill = 'darkred', color='grey22') + 
  geom_jitter(width = 0.2, alpha = 0.4, color = "grey22", size=.5) + 
  gtheme_no_rot + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1)

pdf (paste0 ('Plots/sarcomatoid_score_ENCODE_boxplots2.pdf'), width = 6,height=5)
bp
dev.off()
