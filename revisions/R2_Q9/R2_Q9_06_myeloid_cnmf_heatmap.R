# R2_Q9 -- heatmap of the TAM programme scores per myeloid cell with the sample annotation
# (scatac_myeloid_ArchR.R, Plots/cnmf_scatac_scaled_heatmap2.pdf), improved sample palette.
# The double-scaled module scores saved by the original script (cnmf_scatac.rds) are read back.
source ('R2_Q9_common.R')
proj = readRDS (file.path (root, 'myeloid_cells', 'scatac_ArchR', 'Save-ArchR-Project.rds'))
cnmf_scatac = readRDS (file.path (root, 'myeloid_cells', 'scatac_ArchR', 'cnmf_scatac.rds'))
stopifnot (all (rownames (cnmf_scatac) == rownames (proj@cellColData)))
sample_of_cell = proj$Sample

cap = 3
cnmf_scatac_cap = cnmf_scatac
cnmf_scatac_cap[cnmf_scatac_cap > cap] = cap
cnmf_scatac_cap[cnmf_scatac_cap < -cap] = -cap
set.seed (123)
km_cnmf = kmeans (cnmf_scatac_cap, centers=7) # double scale modules and cluster using k-means
for (v in c('new', 'old'))
  {
  ha = HeatmapAnnotation (sample = sample_of_cell, col=list(sample = list (new = palette_sample, old = palette_sample_old)[[v]][sort (unique (sample_of_cell))]))
  set.seed (123)
  hm = Heatmap (t(cnmf_scatac_cap),
    col = palette_genescore_fun(cnmf_scatac_cap),
    top_annotation = ha,
    clustering_distance_columns = 'pearson',
    clustering_distance_rows = 'pearson',
    show_column_dend = F,
    column_split = km_cnmf$cluster,
    row_names_gp = gpar (fontsize = 8),
    column_names_gp = gpar (fontsize = 0),
    border=T)
  pdf (file.path (outdir, paste0 ('myeloid_cnmf_scatac_scaled_heatmap2', ifelse (v == 'old', '_old_palette', ''), '.pdf')), height=1.5)
  draw (hm)
  dev.off()
  }
