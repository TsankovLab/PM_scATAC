###############################################################
# Build data frame of deviation / gene score / expression
# differences between normal pleura and mesothelioma tumors
###############################################################

# ----------------------------
# Define sample groups
# ----------------------------
tumor_sams       <- c('P1','P10','P11','P12','P13','P14','P23','P3','P4','P5','P8')
normal_sams      <- c('normal1')

tumor_sams_rna   <- c('P1','P11','P12','P13','P14','P3','P4','P5','P8')
normal_sams_rna  <- 'normal_pleura'

###############################################################
# 1. Process chromVAR Deviation Matrix
###############################################################

# Fetch Motif matrix if missing
if (!exists('mSE')) mSE <- fetch_mat(archp, 'Motif')

archp_meta  <- as.data.frame(archp@cellColData)
metaGroupName <- 'Sample3'

# Extract deviation matrix and rename rows
mMat <- assays(mSE)[[1]]
rownames(mMat) <- rowData(mSE)$name

# Aggregate deviations per sample (average across cells)
mMat_agg <- as.data.frame(t(mMat))
mMat_agg$metaGroup <- as.character(archp_meta[, metaGroupName])
mMat_agg <- aggregate(. ~ metaGroup, mMat_agg, mean)

# Move sample IDs to rownames
rownames(mMat_agg) <- mMat_agg$metaGroup
mMat_agg <- t(mMat_agg[,-1])

# Keep only TFs present in the Seurat RNA object
mMat_agg <- mMat_agg[rownames(mMat_agg) %in% rownames(srt), ]


###############################################################
# 2. Process GeneScore Matrix
###############################################################

if (!exists('gsSE')) gsSE <- fetch_mat(archp, 'GeneScore')

gsMat <- assays(gsSE)[[1]]
rownames(gsMat) <- rowData(gsSE)$name

# Keep only TFs present in deviation matrix and scale
gsMat <- scale(gsMat[rownames(gsMat) %in% rownames(mMat_agg), ])

# Aggregate per sample
gsMat_agg <- as.data.frame(t(gsMat))
gsMat_agg$metaGroup <- as.character(archp_meta[, metaGroupName])
gsMat_agg <- aggregate(. ~ metaGroup, gsMat_agg, mean)

rownames(gsMat_agg) <- gsMat_agg$metaGroup
gsMat_agg <- t(gsMat_agg[,-1])


###############################################################
# 3. Process Seurat RNA Expression
###############################################################

DefaultAssay(srt) <- 'RNA'
metaGroupName <- 'sampleID3'
sample_names_rna <- c('P1','P3','P4','P5','P8','P11','P12','P13','P14','normal_pleura')

# Compute mean expression per sample
ps <- AverageExpression(srt, features = rownames(mMat_agg), group.by = metaGroupName)[[1]]
ps <- log2(ps + 1)

colnames(ps) <- gsub('-', '_', colnames(ps))
ps <- ps[, colnames(ps) %in% sample_names_rna]
ps <- ps[rownames(mMat_agg), ]


###############################################################
# 4. Combine Deviation + GS + RNA into summary DF
###############################################################

TF_diff_rna <- data.frame(
  tumor_dev  = apply(mMat_agg[, tumor_sams], 1, mean),
  normal_dev = mMat_agg[, normal_sams],
  tumor_gs   = apply(gsMat_agg[, tumor_sams], 1, mean)[rownames(mMat_agg)],
  normal_gs  = gsMat_agg[, normal_sams][rownames(mMat_agg)],
  tumor_rna  = apply(ps[, tumor_sams_rna], 1, mean),
  normal_rna = ps[, normal_sams_rna]
)

TF_diff_rna$dev_diff <- TF_diff_rna$tumor_dev - TF_diff_rna$normal_dev
TF_diff_rna$gs_diff  <- TF_diff_rna$tumor_gs  - TF_diff_rna$normal_gs
TF_diff_rna$rna_diff <- TF_diff_rna$tumor_rna - TF_diff_rna$normal_rna


###############################################################
# 5. Wilcoxon Test on chromVAR Deviations (presto)
###############################################################

library(presto)

# Recompute un-aggregated scaled deviations
mMat <- assays(mSE)[[1]]
rownames(mMat) <- rowData(mSE)$name
mMat <- scale(mMat)

# Ensure alignment
stopifnot(all(colnames(mMat) == rownames(archp@cellColData)))

# Build pairwise tumor-vs-normal comparison list
dev_comparisons <- setNames(
  lapply(tumor_sams, function(x) c(x, normal_sams)),
  tumor_sams
)

# Run Wilcoxon tests per tumor sample
dev_res <- lapply(dev_comparisons, function(x) {
  tmp <- wilcoxauc(mMat, y = archp$Sample3, groups_use = x)
  tmp <- tmp[tmp$group == x[1], ]
  tmp[tmp$logFC > 0, ]  # keep tumor-up TFs
})

# Collect union of all tumor-up TFs
pos_features <- unique(unlist(lapply(dev_res, function(x) x$feature)))

# Build p-value matrix for significant features
dev_res_df <- lapply(names(dev_comparisons), function(x) {
  tmp <- dev_res[[x]]
  rownames(tmp) <- tmp$feature
  tmp <- tmp[pos_features, "padj", drop = FALSE]
})
dev_res_df <- do.call(cbind, dev_res_df)

colnames(dev_res_df) <- names(dev_comparisons)
rownames(dev_res_df) <- pos_features


###############################################################
# 6. Feature Selection by Significance & Occurrence
###############################################################

pval_threshold       <- 0.05
occurrence_threshold <- 5

# Build logical significance matrix
dev_res_df_logic <- dev_res_df < pval_threshold
dev_res_df_logic[is.na(dev_res_df_logic)] <- FALSE

# Convert to numeric for rowSums
comb_res_df <- matrix(
  as.numeric(dev_res_df_logic),
  nrow = nrow(dev_res_df_logic),
  ncol = ncol(dev_res_df_logic)
)
rownames(comb_res_df) <- rownames(dev_res_df_logic)

# Keep TFs significant in at least N tumor samples
selected_TF <- rownames(comb_res_df)[rowSums(comb_res_df) >= occurrence_threshold]


###############################################################
# 7. Expression Filtering: keep TFs expressed in tumors
###############################################################

ps_tumor <- AverageExpression(
  srt[, srt$sampleID3 %in% tumor_sams_rna],
  features = selected_TF,
  group.by = metaGroupName
)[[1]]
ps_tumor <- log2(ps_tumor + 1)

min_exp <- 0.5
active_TFs <- rownames(ps_tumor)[apply(ps_tumor, 1, function(x) any(x > min_exp))]

# Keep only active in RNA
selected_TF <- intersect(selected_TF, active_TFs)

# Keep only tumor-up TFs
tf_tumor_pos <- rownames(TF_diff_rna)[TF_diff_rna$rna_diff > 0]
selected_TF <- intersect(selected_TF, tf_tumor_pos)


###############################################################
# 8. Order TFs by mean (dev_diff + rna_diff)
###############################################################

tf_order <- rownames(TF_diff_rna)[order(-rowMeans(TF_diff_rna[, c("dev_diff", "rna_diff")]))]
selected_TF_ordered <- tf_order[tf_order %in% selected_TF]

# Predefined guide genes
guides <- c('PITX1','TCF3','TEAD4','MEF2A','MEF2D','HMGA1','SOX9','TWIST1','BPTF')
TF_labels <- unique(c(head(selected_TF_ordered, 30), guides))


###############################################################
# 9. Annotate TFs for plotting
###############################################################

TF_diff_rna$label <- ''
TF_diff_rna$label[match(selected_TF, rownames(TF_diff_rna))] <- selected_TF

TF_diff_rna$label_top <- ''
TF_diff_rna$label_top[match(TF_labels, rownames(TF_diff_rna))] <- TF_labels

TF_diff_rna$color <- TF_diff_rna$label != ''
TF_diff_rna$alpha <- 0.6
TF_diff_rna$alpha[match(selected_TF, rownames(TF_diff_rna))] <- 0.8

TF_diff_rna$label_guides <- ifelse(rownames(TF_diff_rna) %in% guides, 'guides', 'rest')


# Direction categories (tumor-high vs normal-high)
TF_diff_rna$highest_rna  <- apply(TF_diff_rna[, c("tumor_rna","normal_rna")], 1, function(x) which.max(x))
TF_diff_rna$highest_atac <- apply(TF_diff_rna[, c("tumor_dev","normal_dev")], 1, function(x) which.max(x))

TF_diff_rna$highest_rna  <- ifelse(TF_diff_rna$highest_rna == 1, "tumor", "normal")
TF_diff_rna$highest_atac <- ifelse(TF_diff_rna$highest_atac == 1, "tumor", "normal")

TF_diff_rna$direction <- with(TF_diff_rna,
  ifelse(highest_rna == "tumor" & highest_atac == "tumor", "both_high",
  ifelse(highest_rna == "tumor" & highest_atac == "normal", "dev_low",
  ifelse(highest_rna == "normal" & highest_atac == "tumor", "rna_low",
         "both_low")))
)

# De-emphasize points not high in both modalities
TF_diff_rna$alpha[TF_diff_rna$direction != "both_high"] <- 0.6


###############################################################
# 10. Scatter Plot: Deviation diff vs RNA diff
###############################################################

diff_line <- 0

tf_diff_p <- ggplot(
  TF_diff_rna,
  aes(
    x = dev_diff,
    y = rna_diff,
    size = highest_rna,
    fill = direction,
    label = label_top
  )
) +
  geom_point(aes(color = color, alpha = alpha), shape = 21, stroke = 0.1) +
  scale_color_manual(values = c('FALSE' = "grey", 'TRUE' = "black")) +
  scale_fill_manual(values = c(
    both_low = "grey20",
    rna_low  = "grey20",
    dev_low  = "grey20",
    both_high = "darkred"
  )) +
  geom_vline(xintercept = diff_line, linetype = "dashed", color = "grey44", linewidth = .3) +
  geom_hline(yintercept = diff_line, linetype = "dashed", color = "grey44", linewidth = .3) +
  gtheme_no_rot +
  geom_text_repel(
    aes(color = label_guides, fill = label_guides),
    segment.size = .05,
    size = 1.8,
    max.overlaps = 100
  ) +
  xlab("activity difference") +
  ylab("RNA difference") +
  xlim(c(-0.2, 0.2)) +
  ylim(c(-0.6, 0.6))

pdf("Plots/Diff_normal_tumor_deviation_and_rna_scatterplot5.pdf", width = 5, height = 4)
print(tf_diff_p)
dev.off()

# Save selected TFs for later reuse
saveRDS(selected_TF, "selected_TF.rds")
selected_TF_df = cbind (dev_res_df[selected_TF, ],TF_diff_rna[selected_TF,])
write.csv (selected_TF_df, 'selected_TF_table.csv')
