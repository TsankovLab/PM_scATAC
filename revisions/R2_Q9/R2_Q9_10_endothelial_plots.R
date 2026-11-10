# R2_Q9 -- endothelial adult_vs_fetal_scatterplot_smoothed (scatac_endothelial_ArchR.R), improved palette:
#   fetal accessibility score vs fetal gene-score module, per sample (page 2 = old palette)
# Matrices are fetched from the endothelial project (read-only) and processed as in the original script.
source ('R2_Q9_common.R')
suppressPackageStartupMessages (library (zoo))
source (file.path (utils, 'scATAC_functions.R'))
archp = load_project ('Endothelial/scatac_ArchR')

fSE  = fetch_mat(archp, 'scATAC_datasets')
nSE2 = fetch_mat(archp, 'scATAC_normal2')
dense = function (se) as.matrix (assays (se)[[1]])

# Scatterplots of fetal accessibility score vs fetal gene-score module, rolling means ####
fMat = dense (fSE)
rownames (fMat) = gsub ('rawlins_fetal_lung_','', rownames (fMat))
nMat = dense (nSE2)
rownames (nMat) = gsub ('Tsankov_lung_','', rownames (nMat))
# Get endothelial subset module and scale them per cell to reduce seq-depth bias
modMat = as.data.frame (t(scale(t(as.matrix (archp@cellColData[, c('fetal.Artery','fetal.COL4A1','fetal.Vein')])))))
fnmMat = cbind (scale(t(fMat)), scale(t(nMat)))
fnmMat = fnmMat[,!duplicated(colnames(fnmMat))]

sams = c('P1','P10','P11','P14','P4') # Select samples that have at least 100 endothelial cells
frip = archp$FRIP
bin_width <- 30   # Number of observations per bin
overlap <- 30
mk_scatter = function (ext_dev, pal)
  {
  df = list()
  for (sam in sams)
    {
    fscore_order = modMat$fetal.COL4A1[archp$Sample == sam]
    fscore_order = order (fscore_order)
    df[[sam]] <- data.frame (
      fetal_ta = rollapply (fnmMat[archp$Sample == sam, ext_dev][fscore_order],width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"),
      fetal_gs = rollapply (modMat$fetal.COL4A1[archp$Sample == sam][fscore_order], width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"),
      FRIP = rollapply (frip[archp$Sample == sam][fscore_order], width = bin_width, FUN = mean, by = overlap, partial = TRUE, align = "left"),
      sample = sam)
    }
  df = do.call (rbind, df)
  ggplot(df, aes(x = fetal_ta, y = fetal_gs, color = sample)) +
    geom_point(
      alpha = 0.7,
      size = 2
    ) +
    geom_smooth(
      aes (group = sample),
      method = "lm",
      color = "dodgerblue",
      fill = "lightblue",
      se = FALSE,
      linetype = "solid",
      size = 1
    ) +
    geom_smooth(method=lm, fullrange=FALSE, se=FALSE)+
    facet_wrap (~sample, scales = 'free', ncol=length(sams), strip.position = "top") +
    gtheme_no_rot +
    scale_color_manual(values = pal)
  }
pdf (file.path (outdir, 'endothelial_adult_vs_fetal_scatterplot_smoothed.pdf'), width=9, height=4)
print (wrap_plots (list(mk_scatter ('Endothelial', palette_sample), mk_scatter ('Midcap', palette_sample)), nrow = 2))
print (wrap_plots (list(mk_scatter ('Endothelial', palette_sample_old), mk_scatter ('Midcap', palette_sample_old)), nrow = 2))
dev.off()
