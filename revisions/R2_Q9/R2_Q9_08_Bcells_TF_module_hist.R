# R2_Q9 -- distribution of the TF module activities per sample in B cells
# (scatac_Bcells_ArchR.R, Plots/TF_modules_dist_sample_hist.pdf), improved sample palette.
# The TF modules are rebuilt as in the original (k-means, 5 centres, on the Spearman correlation of the
# activity of TFs expressed in scRNA-seq), because the module scores are not stored in the project.
source ('R2_Q9_common.R')
suppressPackageStartupMessages (library (Seurat))
source (file.path (utils, 'useful_functions.R'))
source (file.path (utils, 'scATAC_functions.R'))
archp = load_project ('B_cells/scatac_ArchR')
srt = readRDS (file.path (root, 'B_cells', 'scrna', 'srt.rds'))

mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames (mMat) = rowData (mSE)$name
mMat = as.matrix(mMat)

# Filter by RNA expression ####
metaGroupName = 'seurat_clusters'
active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)
mMat = mMat[active_TFs, ]

mMat_cor = cor (as.matrix(t(scale(mMat))), method = 'spearman')

set.seed(1234)
centers=5
km = kmeans (mMat_cor, centers=centers)
print (table (km$cluster))

tf_modules = lapply (unique(km$cluster), function(x) colMeans (mMat[names(km$cluster[km$cluster == x]),]))
names (tf_modules) = paste0('mod_',unique(km$cluster))
tf_modules = do.call (cbind, tf_modules)
ccomp = cbind (as.data.frame (archp@cellColData), tf_modules)

mk = function (pal) wrap_plots (
  ggplot(ccomp, aes(x = mod_2, fill = Sample)) +
    geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
    scale_fill_manual(values= pal) + gtheme,
  ggplot(ccomp, aes(x = mod_1, fill = Sample)) +
    geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
    scale_fill_manual(values= pal) + gtheme)
pdf (file.path (outdir, 'Bcells_TF_modules_dist_sample_hist.pdf'),7,3)
print (mk (palette_sample))
print (mk (palette_sample_old))
dev.off()
