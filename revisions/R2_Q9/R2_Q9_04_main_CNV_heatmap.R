# R2_Q9 -- CNV heatmap of the malignant cells per sample (scatac_main_ArchR.R, Plots/GL_method_z_CNV_heatmap2.pdf)
# with the improved sample palette in the row annotation. The per-sample CNV objects saved by the original
# script are read back; nothing is recomputed.
source ('R2_Q9_common.R')
archp = load_project ('main/scatac_ArchR')
projdir_cnv = file.path (root, 'per_sample_QC_signac', 'CNV_analysis')
ws = 10e6
ss = 5e6

sample_names = unique (archp$Sample)
cnaObj_l = list()
for (sam in sample_names) cnaObj_l[[sam]] = readRDS (file.path(projdir_cnv, paste0('CNV_LFC_GC_',sam,'_ws_',ws,'_ss_',ss,'.rds')))

# Print CNV heatmap
malignant_cells = c(rownames(archp@cellColData)[archp$celltype_lv1 == 'Malignant'])

cnaMat_obj_l = list()
for (sam in sample_names)
    {
    mat_type = 'z'
    cnaObj_mat = t(cnaObj_l[[sam]]@assays@data@listData[[mat_type]])
    colnames (cnaObj_mat) = paste0(seqnames(rowRanges(cnaObj_l[[sam]])),':', ranges(rowRanges(cnaObj_l[[sam]])))
    bc = colnames(cnaObj_l[[sam]])
    if (any (grep ('\\#', colnames(cnaObj_l[[sam]])))) bc = sapply (colnames(cnaObj_l[[sam]]), function(x) unlist(strsplit(x, '\\#'))[2])
    rownames (cnaObj_mat) = paste0(sam,'#',bc)
    cnaObj_mat = cnaObj_mat[rownames(cnaObj_mat) %in% malignant_cells,]
    cnaObj_mat[is.na(cnaObj_mat)] = 0
    cnaObj_mat[is.infinite(cnaObj_mat)]= 0
    cnaMat_obj_l[[sam]] = cnaObj_mat
    }
cnaMat_obj_avg = lapply (cnaMat_obj_l, function(x) colMeans (x))
cnaMat_obj_avg = do.call (cbind, cnaMat_obj_avg)
cnaMat_obj_avg = cnaMat_obj_avg[,sample_names]
# row_chr is not defined in the published script: chromosome of each window, in genome order
row_chr = sapply (rownames (cnaMat_obj_avg), function(x) unlist (strsplit (x, ':'))[1])
row_chr = factor (row_chr, levels = unique (row_chr))

for (v in c('new', 'old'))
  {
  ha = rowAnnotation (sample = colnames(cnaMat_obj_avg), col = list(sample = list (new = palette_sample, old = palette_sample_old)[[v]][colnames(cnaMat_obj_avg)]))
  pdf (file.path (outdir, paste0('main_GL_method_',mat_type,'_CNV_heatmap2', ifelse (v == 'old', '_old_palette', ''), '.pdf')), width=6,height=3)
      print (Heatmap (t(cnaMat_obj_avg),
      cluster_rows=F,
      cluster_columns=F,
      left_annotation = ha,
      row_names_side = 'left',
      column_title_side = 'bottom',
      border=T,
      row_names_gp = gpar (fontsize = 8),
      column_names_gp = gpar(fontsize = 0),
      column_title_gp = gpar(fontsize = 7),
      column_split = row_chr,
      row_split = ifelse ('normal1' == colnames(cnaMat_obj_avg),'normal','tumor'),
      column_gap = unit(.5, "mm"),
      cluster_column_slices = F,
      row_title_gp = gpar(fontsize = 1),
      column_title_rot=45))
  dev.off()
  }
