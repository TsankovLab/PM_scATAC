# R2_Q9 -- heatmap of the top differential TF activities per tumour cluster with the sample annotation
# (scatac_tumor_SOX9.R, Plots/DAM_clusters_Clusters_sample_heatmaps.pdf), improved sample palette.
# The differential-activity table saved by the original script (DAM_table_Clusters_sample.csv) is read back;
# the chromVAR deviations are fetched from the project, which is not modified.
source ('R2_Q9_common.R')
suppressPackageStartupMessages (library (stringr))
source (file.path (utils, 'scATAC_functions.R'))
archp = load_project ('tumor_compartment/scatac_ArchR')
archp_NN = archp[!archp$Sample3 %in% c('normal1')]

archp_NN$Clusters_sample = paste0(archp_NN$Clusters, '_',archp_NN$Sample2)
Clusters_sample = paste0(archp_NN$Clusters, '_',archp_NN$Sample2)
remove_low_clusters = !Clusters_sample %in% names(table (Clusters_sample)[table (Clusters_sample) < 10])
Clusters_sample = Clusters_sample[remove_low_clusters]
archp_NN = archp_NN[remove_low_clusters]
metaGroupName = "Clusters_sample"

DAM_df = read.csv (file.path (root, 'tumor_compartment', 'scatac_ArchR', paste0('DAM_table_',metaGroupName, '.csv')), row.names = 1)
top_genes = 3
DAM_df <- DAM_df %>%
group_by (comparison) %>%
slice_head(n = top_genes) %>%
ungroup() %>%
  mutate(
    second_term = str_extract(comparison, "(?<=_)P\\d+"),
    second_num  = as.numeric(str_extract(second_term, "\\d+"))
  ) %>%
  arrange(second_num) %>%
  select(-second_term, -second_num)

mSE = fetch_mat (archp, 'Motif')
mMat = assays (mSE)[[1]]
rownames(mMat) = rowData(mSE)$name
mMat = mMat[unique (DAM_df$gene), rownames(archp_NN@cellColData)]
mMat_mg = as.data.frame (as.matrix (t(mMat)))   # as.matrix: the deviation matrix comes back sparse in this session
mMat_mg$metaGroup = as.character (archp_NN@cellColData[,metaGroupName])
mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
rownames (mMat_mg) = mMat_mg[,1]
mMat_mg = mMat_mg[,-1]
order_p <- as.numeric(str_extract(rownames(mMat_mg), "(?<=_P)\\d+"))
mMat_mg <- mMat_mg[order(order_p), ]

ha = rowAnnotation (foo = anno_mark(at = match(head(DAM_df$gene[DAM_df$comparison=='C1_P11'],top_genes),colnames(mMat_mg)),
    labels = head(DAM_df$gene[DAM_df$comparison=='C1_P11'],top_genes), labels_gp = gpar(fontsize = 8, fontface = 'italic')))
split_raw <- sapply(rownames(mMat_mg), function(x) unlist(strsplit(x,"_"))[2])
p_nums <- as.numeric(str_extract(split_raw, "\\d+"))
unique_levels <- unique(split_raw[order(p_nums)])
split_factor <- factor(split_raw, levels = unique_levels)
print (table (split_raw))

for (v in c('new', 'old'))
  {
  ha2 = HeatmapAnnotation (sample = split_raw, col= list(sample = list (new = palette_sample, old = palette_sample_old)[[v]][unique_levels]))
  DAM_hm = Heatmap (t(scale(mMat_mg)),
          top_annotation = ha2,
          column_split = split_factor,
          right_annotation = ha,
          row_labels = colnames (mMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          cluster_columns=F,
          row_names_gp = gpar(fontsize = 0, fontface = 'italic'),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          name = 'TF activity',
          border=TRUE,
          col = palette_deviation_fun(scale(mMat_mg)))
  pdf (file.path (outdir, paste0('tumor_DAM_clusters_',metaGroupName,'_heatmaps', ifelse (v == 'old', '_old_palette', ''), '.pdf')), width = 6.5, height = 3.5)
  print(DAM_hm)
  dev.off()
  }
