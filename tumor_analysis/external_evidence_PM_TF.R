
### Look at overlap of TF with CNV from TCGA ####
## Generate circos plot for showing hubs on recurrent CNVs ####
# source script to load TCGA_CNV data
source (file.path('..','..','git_repo','tumor_analysis','compile_TCGA_CNV_data.R'))

# Find genomic loci of active TFs ####
all_genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
all_genes$gene = as.data.frame(org.Hs.egSYMBOL)[match (all_genes$gene_id, as.data.frame(org.Hs.egSYMBOL)[,1]),'symbol']
all_genes = as.data.frame (all_genes[all_genes$gene %in% selected_TF,])
all_genes = all_genes[, c('seqnames','start','end','gene')]


pdf (file.path ('Plots','CNV_TFs_circos.pdf'))
circos.initializeWithIdeogram (species = "hg38", chromosome.index= paste0('chr',1:22))
circos.genomicTrack (cnv_mat_avg,
    panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, bg.border = rep(0, 22),
            col = ifelse(value[[1]] > 0, "red", "blue"), border=NA, ...)
        #circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
})
# circos.genomicTrack (cnv_hubs_df_track,
#     panel.fun = function(region, value, ...) {
#         circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, 
#             col = ifelse(value[[1]] > 0, "red", "blue"), border=NA, ...)
#         #circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
# })
# circos.genomicHeatmap(cnv_hubs_df_track, col = col_fun, side = "inside", border = "white")
circos.genomicLabels(all_genes, labels.column = 4, side = "inside",
    col = 'grey22', line_col = 'grey22',padding = 0.6, cex=0.5)
dev.off()


# Also get CNV value for each TF ####
all_genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
all_genes$gene = as.data.frame(org.Hs.egSYMBOL)[match (all_genes$gene_id, as.data.frame(org.Hs.egSYMBOL)[,1]),'symbol']
all_genes = as.data.frame (all_genes[all_genes$gene %in% rownames(dev_res_df),])
all_genes = all_genes[, c('seqnames','start','end','gene')]

all_genes_gr = makeGRangesFromDataFrame (all_genes, keep.extra.columns = TRUE)
cnv_mat_avg_gr = makeGRangesFromDataFrame (cnv_mat_avg, keep.extra.columns = TRUE)
tf_idx = subjectHits (findOverlaps(all_genes_gr, cnv_mat_avg_gr))
tf_idx2 = queryHits (findOverlaps(all_genes_gr, cnv_mat_avg_gr))
cnv_mat_avg_gr = cnv_mat_avg_gr[tf_idx]
cnv_mat_avg_gr$gene = all_genes_gr$gene [tf_idx2]
cnv_mat_avg_df = as.data.frame (cnv_mat_avg_gr, row.names=NULL)
cnv_mat_avg_df = do.call (rbind, lapply (split (cnv_mat_avg_df, cnv_mat_avg_df$gene), function(x) data.frame (gene = x$gene[1], CNV_avg_log = mean (x$CNV_avg_log))))
write.csv (cnv_mat_avg_df, 'TF_CNV_TCGA.csv')

# Show differences in TF CNV score between up and down regulated in normal ####
metaGroupName = 'sampleID'
active_TFs = exp_genes (srt, rownames(mMat), min_exp = 0.1, metaGroupName)

cnv_mat_avg_df$selected_TFs = cnv_mat_avg_df$gene %in% selected_TF
#cnv_mat_avg_df = cnv_mat_avg_df[cnv_mat_avg_df$gene %in% active_TFs,]

bp = ggplot (cnv_mat_avg_df, aes (x = selected_TFs, y= CNV_avg_log)) + 
bxp + gtheme

stat.test = bp$data %>%
          wilcox_test (reformulate ('selected_TFs', 'CNV_avg_log')) %>%
          adjust_pvalue (method = "none") %>%
          add_significance ()
stat.test = stat.test %>% add_xy_position (x = 'selected_TFs', step.increase=0.5)

    bp = bp + stat_pvalue_manual (stat.test, 
      remove.bracket=FALSE,
      bracket.nudge.y = 0, 
      hide.ns = F,
      label = "p.adj") + 
      NoLegend()
    
pdf (file.path ('Plots','CNV_score_in_selected_TF.pdf'))
bp
dev.off()

### Check selected TF in TCGA ATAC ####
# Add deviations from TCGA ATAC bulk ####
tcga_dev = readRDS ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/TCGA_atac/tcga_meso_atac_deviations_clinical_info.rds')
tcga_mat = assays (tcga_dev)$z
tcga_mat = aggregate (t(tcga_mat), by = list(sampleID= colnames(tcga_mat)), mean)
rownames (tcga_mat) = tcga_mat[,1]
tcga_mat = tcga_mat[,-1]
tcga_mat = as.data.frame (t(tcga_mat))
rownames (tcga_mat) = unname (sapply (rownames(tcga_mat), function(x) unlist (strsplit (x, '_'))[3]))



# Plot heatmaps of candidate TFs ####
# Heatmap of TF deviations ####
column_split = ifelse (grepl ('normal_pleura', colnames(mMat_agg)), 'Normal','Tumor')
mMat_agg_tf = mMat_agg[selected_TF_ordered,]

TF_cluster_selected_hm = Heatmap (mMat_agg_tf,
        #right_annotation=tf_mark,
        column_split = column_split,
        cluster_rows = F, #km = 4, 
        name = 'z-score\ndeviations',
        column_gap = unit(.8, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T, 
        col = palette_deviation,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 5), 
        border=T,
        width = unit(2, "cm"))

#ps_order = ps[tf_name,sarc_order$sampleID]
# Add scRNA data ####
column_split_rna = ifelse (grepl ('HU', colnames(ps)), 'Normal','Tumor')
ps_tf = ps[selected_TF_ordered,]
TF_exp_selected_hm = Heatmap (ps_tf,
        #right_annotation=tf_mark,
        column_split = column_split_rna,
        cluster_rows = T, #km = 4, 
        name = 'expression',
        column_gap = unit(.5, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T, 
        col = palette_expression,
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 5), 
        border=T,
        width = unit(2, "cm"))

# Add regulons enriched for each TF ####
TF_regulons_hm = Heatmap (t(scale(t(nmf_TF_df[selected_TF_ordered,]))),
        #right_annotation=tf_mark,
#        column_split = column_split,
        cluster_rows = F, #km = 4, 
        name = 'expression',
        column_gap = unit(.2, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T,
        column_names_rot = 45, 
        col = viridis::inferno (100),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5), 
        border=F,
        width = unit(5, "cm"))

# Add correlation to bulk RNA sarcomatoid score ####
tf_sarc_cor = read.csv ('../../bulkRNA_meso/activeTF_sarcomatoid_correlation.csv', row.names = 1)
tf_sarc_cor_tf = tf_sarc_cor[selected_TF_ordered,]
TF_bulk_sarc_cor_hm = Heatmap (tf_sarc_cor_tf,
        #right_annotation=tf_mark,
        column_split = colnames (tf_sarc_cor),
        cluster_rows = F, #km = 4, 
        name = 'sarc_cor_bulk',
        column_gap = unit(.2, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T,
        column_names_rot = 45, 
        #col = palette_expression (100),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5), 
        border=F,
        width = unit(.6, "cm"))

# Add cNMF vs TF correlation from scRNA ####


cnmf_tf_cor_hm = Heatmap (t(cnmf_tf_cor)[selected_TF_ordered,],
        #right_annotation=tf_mark,
#        column_split = colnames (cnmf_tf_cor),
        cluster_rows = F, #km = 4, 
        name = 'sarc_cor_bulk',
        column_gap = unit(.2, "mm"),
        row_gap = unit(.2, "mm"),
        clustering_distance_rows = 'euclidean',
        clustering_distance_columns = 'euclidean',
        cluster_columns=T,
        column_names_rot = 45, 
        #col = palette_expression (100),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5), 
        border=F,
        width = unit(6, "cm"))


pdf (paste0 ('Plots/selected_TF_dev_exp_significant_bulkcor_heatmaps2.pdf'), width = 8,height=5)
draw (TF_cluster_selected_hm + TF_exp_selected_hm + cnmf_tf_cor_hm + TF_bulk_sarc_cor_hm)
dev.off()

# Import Viestra archetypes to add to the table below ####
motif_clusters = read.csv ('../../Viestra_motif_clustered.csv')
colnames (motif_clusters)[1] = 'Cluster'
head (motif_clusters)
motif_clusters$Cluster2 = sapply (motif_clusters$Motif, function(x) paste(unique(motif_clusters[motif_clusters$Motif == x,'Cluster']),collapse=','))

# Add correlation to bulk RNA sarcomatoid score ####
tf_sarc_cor = read.csv ('../../bulkRNA_meso/activeTF_sarcomatoid_correlation.csv', row.names = 1)
tf_sarc_cor_tf = tf_sarc_cor[selected_TF_ordered,]

# Export table ####
tf_cnv_tcga = read.csv ('TF_CNV_TCGA.csv')
ccle_exp = read.csv ('/ahg/regevdata/projects/ICA_Lung/Bruno/Public_data/CCLE/meso_scATAC_active_TFS.csv', row.names=1)
ccle_exp = ccle_exp[selected_TF,]
ccle_exp = ccle_exp[colnames(ccle_exp) != 'gene']
enr_path = t (readRDS ('../scrna/enrichment_pathways_TFs.rds'))
enr_path = enr_path[selected_TF,]
rna_res_sum_p = apply (do.call (cbind, lapply (rna_res, function(x) x[match(selected_TF, x$feature),c('padj'), drop=F])), 1, mean)
rna_res_sum_lfc = apply (do.call (cbind, lapply (rna_res, function(x) x[match(selected_TF, x$feature),c('logFC'), drop=F])),1, mean)
rna_res_sum_ae = apply (do.call (cbind, lapply (rna_res, function(x) x[match(selected_TF, x$feature),c('avgExpr'), drop=F])),1, mean)
rna_res_sum = data.frame (RNA_pvalue_mean = rna_res_sum_p, RNA_lfc_mean = rna_res_sum_lfc, RNA_avExpr_mean = rna_res_sum_ae)

dev_res_sum_p = apply (do.call (cbind, lapply (dev_res, function(x) x[match(selected_TF, x$feature),c('padj'), drop=F])), 1, mean)
dev_res_sum_lfc = apply (do.call (cbind, lapply (dev_res, function(x) x[match(selected_TF, x$feature),c('logFC'), drop=F])), 1, mean)
dev_res_sum_ae = apply (do.call (cbind, lapply (dev_res, function(x) x[match(selected_TF, x$feature),c('avgExpr'), drop=F])), 1, mean)
dev_res_sum = data.frame (DEV_pvalue_mean = dev_res_sum_p, DEV_lfc_mean = dev_res_sum_lfc, DEV_avExpr_mean = dev_res_sum_ae)



tf_table = cbind (
  dev_res_sum, 
  rna_res_sum, 
  tf_sarc_cor_tf[selected_TF,], 
  tf_cnv_tcga[match(selected_TF,tf_cnv_tcga$gene),'CNV_avg_log', drop=F], 
  motif_clusters[match (selected_TF, motif_clusters$Motif),'Cluster2',drop=F],
  ccle_exp,   
  enr_path)
rownames (tf_table) = selected_TF
write.csv (tf_table, 'candidate_TFs_table.csv')





