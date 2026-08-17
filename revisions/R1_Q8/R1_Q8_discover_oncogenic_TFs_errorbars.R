###############################################################
# R1_Q8 : Reproduce Diff_normal_tumor_deviation_and_rna_scatterplot5.pdf
#          (from git_repo/tumor_analysis/discover_oncogenic_TFs.R)
#          and ADD error bars (+/- SEM across samples) for
#          TF activity (chromVAR deviation, x) and RNA (y).
#
# Faithful re-implementation of the original logic, but the per-sample
# aggregated matrices are retained so per-TF spread across samples can
# be turned into error bars.
###############################################################

set.seed(1234)

GITROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo"
BASE    <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment"
OUTDIR  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q8"
setwd(OUTDIR)
dir.create(file.path(OUTDIR, "Plots"), showWarnings = FALSE)

# ---- packages / helpers / aesthetics (same as original pipeline) ----
source(file.path(GITROOT, "utils", "load_packages.R"))
source(file.path(GITROOT, "utils", "scATAC_functions.R"))   # fetch_mat + hidden ArchR fns
source(file.path(GITROOT, "utils", "ggplot_aestetics.R"))   # gtheme_no_rot, etc.

addArchRThreads(threads = 1)
addArchRGenome("hg38")

# ----------------------------
# Define sample groups (verbatim from original)
# ----------------------------
tumor_sams       <- c('P1','P10','P11','P12','P13','P14','P23','P3','P4','P5','P8')
normal_sams      <- c('normal1')

tumor_sams_rna   <- c('P1','P11','P12','P13','P14','P3','P4','P5','P8')
normal_sams_rna  <- 'normal_pleura'

# ---- Load ArchR project & Seurat RNA object ----
archp <- loadArchRProject(file.path(BASE, "scatac_ArchR"), showLogo = FALSE)

srt <- readRDS(file.path(BASE, "scrna", "srt.rds"))
# The original merges HU37/HU62 into 'normal_pleura'
srt$sampleID3[srt$sampleID3 %in% c('HU37','HU62')] <- 'normal_pleura'

###############################################################
# 1. chromVAR Deviation Matrix -> per-sample aggregate
###############################################################
mSE <- fetch_mat(archp, 'Motif')

archp_meta    <- as.data.frame(archp@cellColData)
metaGroupName <- 'Sample3'

mMat <- assays(mSE)[[1]]
rownames(mMat) <- rowData(mSE)$name

# Aggregate deviations per sample (average across cells)  -> TF x sample
mMat_agg <- as.data.frame(t(mMat))
mMat_agg$metaGroup <- as.character(archp_meta[, metaGroupName])
mMat_agg <- aggregate(. ~ metaGroup, mMat_agg, mean)
rownames(mMat_agg) <- mMat_agg$metaGroup
mMat_agg <- t(mMat_agg[,-1])

# Keep only TFs present in the Seurat RNA object
mMat_agg <- mMat_agg[rownames(mMat_agg) %in% rownames(srt), ]

###############################################################
# 2. GeneScore Matrix -> per-sample aggregate
###############################################################
gsSE <- fetch_mat(archp, 'GeneScore')
gsMat <- assays(gsSE)[[1]]
rownames(gsMat) <- rowData(gsSE)$name
gsMat <- scale(gsMat[rownames(gsMat) %in% rownames(mMat_agg), ])

gsMat_agg <- as.data.frame(t(gsMat))
gsMat_agg$metaGroup <- as.character(archp_meta[, metaGroupName])
gsMat_agg <- aggregate(. ~ metaGroup, gsMat_agg, mean)
rownames(gsMat_agg) <- gsMat_agg$metaGroup
gsMat_agg <- t(gsMat_agg[,-1])

###############################################################
# 3. Seurat RNA Expression -> per-sample aggregate (log2 mean)
###############################################################
DefaultAssay(srt) <- 'RNA'
metaGroupName_rna <- 'sampleID3'
sample_names_rna  <- c('P1','P3','P4','P5','P8','P11','P12','P13','P14','normal_pleura')

ps <- AverageExpression(srt, features = rownames(mMat_agg), group.by = metaGroupName_rna)[[1]]
ps <- log2(ps + 1)
colnames(ps) <- gsub('-', '_', colnames(ps))
ps <- ps[, colnames(ps) %in% sample_names_rna]
ps <- ps[rownames(mMat_agg), ]

###############################################################
# 4. Combine Deviation + GS + RNA into summary DF  (verbatim)
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

# ---- NEW: per-TF spread across samples -> SEM for error bars ----
# Activity (deviation) spread is across the tumor samples; 'normal1' is a
# single sample so subtracting it shifts the mean but not the spread.
sem <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) return(NA_real_); sd(x) / sqrt(length(x)) }

dev_tumor_mat <- mMat_agg[, tumor_sams, drop = FALSE]
rna_tumor_mat <- ps[, tumor_sams_rna, drop = FALSE]

TF_diff_rna$dev_sem <- apply(dev_tumor_mat[rownames(TF_diff_rna), , drop = FALSE], 1, sem)  # x error
TF_diff_rna$rna_sem <- apply(rna_tumor_mat[rownames(TF_diff_rna), , drop = FALSE], 1, sem)  # y error
TF_diff_rna$dev_sd  <- apply(dev_tumor_mat[rownames(TF_diff_rna), , drop = FALSE], 1, function(x) sd(x[is.finite(x)]))
TF_diff_rna$rna_sd  <- apply(rna_tumor_mat[rownames(TF_diff_rna), , drop = FALSE], 1, function(x) sd(x[is.finite(x)]))

###############################################################
# 5. Wilcoxon Test on chromVAR Deviations (presto)  (verbatim)
###############################################################
library(presto)
mMat <- assays(mSE)[[1]]
rownames(mMat) <- rowData(mSE)$name
mMat <- scale(mMat)
stopifnot(all(colnames(mMat) == rownames(archp@cellColData)))

dev_comparisons <- setNames(lapply(tumor_sams, function(x) c(x, normal_sams)), tumor_sams)
dev_res <- lapply(dev_comparisons, function(x) {
  tmp <- wilcoxauc(mMat, y = archp$Sample3, groups_use = x)
  tmp <- tmp[tmp$group == x[1], ]
  tmp[tmp$logFC > 0, ]
})
pos_features <- unique(unlist(lapply(dev_res, function(x) x$feature)))
dev_res_df <- lapply(names(dev_comparisons), function(x) {
  tmp <- dev_res[[x]]; rownames(tmp) <- tmp$feature
  tmp[pos_features, "padj", drop = FALSE]
})
dev_res_df <- do.call(cbind, dev_res_df)
colnames(dev_res_df) <- names(dev_comparisons)
rownames(dev_res_df) <- pos_features

###############################################################
# 6. Feature Selection by Significance & Occurrence  (verbatim)
###############################################################
pval_threshold       <- 0.05
occurrence_threshold <- 5
dev_res_df_logic <- dev_res_df < pval_threshold
dev_res_df_logic[is.na(dev_res_df_logic)] <- FALSE
comb_res_df <- matrix(as.numeric(dev_res_df_logic), nrow = nrow(dev_res_df_logic), ncol = ncol(dev_res_df_logic))
rownames(comb_res_df) <- rownames(dev_res_df_logic)
selected_TF <- rownames(comb_res_df)[rowSums(comb_res_df) >= occurrence_threshold]

###############################################################
# 7. Expression Filtering: keep TFs expressed in tumors  (verbatim)
###############################################################
ps_tumor <- AverageExpression(
  srt[, srt$sampleID3 %in% tumor_sams_rna],
  features = selected_TF, group.by = metaGroupName_rna
)[[1]]
ps_tumor <- log2(ps_tumor + 1)
min_exp <- 0.5
active_TFs <- rownames(ps_tumor)[apply(ps_tumor, 1, function(x) any(x > min_exp))]
selected_TF <- intersect(selected_TF, active_TFs)
tf_tumor_pos <- rownames(TF_diff_rna)[TF_diff_rna$rna_diff > 0]
selected_TF  <- intersect(selected_TF, tf_tumor_pos)

###############################################################
# 8. Order TFs by mean (dev_diff + rna_diff)  (verbatim)
###############################################################
tf_order <- rownames(TF_diff_rna)[order(-rowMeans(TF_diff_rna[, c("dev_diff", "rna_diff")]))]
selected_TF_ordered <- tf_order[tf_order %in% selected_TF]
guides <- c('PITX1','TCF3','TEAD4','MEF2A','MEF2D','HMGA1','SOX9','TWIST1','BPTF')
TF_labels <- unique(c(head(selected_TF_ordered, 30), guides))

###############################################################
# 9. Annotate TFs for plotting  (verbatim)
###############################################################
TF_diff_rna$label <- ''
TF_diff_rna$label[match(selected_TF, rownames(TF_diff_rna))] <- selected_TF
TF_diff_rna$label_top <- ''
TF_diff_rna$label_top[match(TF_labels, rownames(TF_diff_rna))] <- TF_labels
TF_diff_rna$color <- TF_diff_rna$label != ''
TF_diff_rna$alpha <- 0.6
TF_diff_rna$alpha[match(selected_TF, rownames(TF_diff_rna))] <- 0.8
TF_diff_rna$label_guides <- ifelse(rownames(TF_diff_rna) %in% guides, 'guides', 'rest')

TF_diff_rna$highest_rna  <- apply(TF_diff_rna[, c("tumor_rna","normal_rna")], 1, function(x) which.max(x))
TF_diff_rna$highest_atac <- apply(TF_diff_rna[, c("tumor_dev","normal_dev")], 1, function(x) which.max(x))
TF_diff_rna$highest_rna  <- ifelse(TF_diff_rna$highest_rna == 1, "tumor", "normal")
TF_diff_rna$highest_atac <- ifelse(TF_diff_rna$highest_atac == 1, "tumor", "normal")
TF_diff_rna$direction <- with(TF_diff_rna,
  ifelse(highest_rna == "tumor" & highest_atac == "tumor", "both_high",
  ifelse(highest_rna == "tumor" & highest_atac == "normal", "dev_low",
  ifelse(highest_rna == "normal" & highest_atac == "tumor", "rna_low", "both_low"))))
TF_diff_rna$alpha[TF_diff_rna$direction != "both_high"] <- 0.6

###############################################################
# 10a. Scatter Plot: Deviation diff vs RNA diff  (EXACT reproduction)
###############################################################
diff_line <- 0
tf_diff_p <- ggplot(TF_diff_rna,
    aes(x = dev_diff, y = rna_diff, size = highest_rna, fill = direction, label = label_top)) +
  geom_point(aes(color = color, alpha = alpha), shape = 21, stroke = 0.1) +
  scale_color_manual(values = c('FALSE' = "grey", 'TRUE' = "black")) +
  scale_fill_manual(values = c(both_low = "grey20", rna_low = "grey20",
                               dev_low = "grey20", both_high = "darkred")) +
  geom_vline(xintercept = diff_line, linetype = "dashed", color = "grey44", linewidth = .3) +
  geom_hline(yintercept = diff_line, linetype = "dashed", color = "grey44", linewidth = .3) +
  gtheme_no_rot +
  geom_text_repel(aes(color = label_guides, fill = label_guides),
                  segment.size = .05, size = 1.8, max.overlaps = 100) +
  xlab("activity difference") + ylab("RNA difference") +
  xlim(c(-0.2, 0.2)) + ylim(c(-0.6, 0.6))

pdf(file.path("Plots", "Diff_normal_tumor_deviation_and_rna_scatterplot5.pdf"), width = 5, height = 4)
print(tf_diff_p)
dev.off()

###############################################################
# 10b. Same scatter WITH error bars (+/- SEM across samples)
#      Error bars drawn only for the selected (colored) TFs so the
#      plot stays legible; x = activity SEM, y = RNA SEM.
###############################################################
# Error bars only for the LABELLED genes (those carrying a text label).
eb_df <- TF_diff_rna[TF_diff_rna$label_top != '', ]

# Dot size encodes RNA expression level (mean log2 tumor expression),
# kept in a small range so points stay compact like the original figure.
tf_diff_eb <- ggplot(TF_diff_rna,
    aes(x = dev_diff, y = rna_diff, size = tumor_rna, fill = direction, label = label_top)) +
  # error bars first so points sit on top
  geom_errorbar(data = eb_df, inherit.aes = FALSE,
                aes(x = dev_diff, ymin = rna_diff - rna_sem, ymax = rna_diff + rna_sem),
                width = 0, linewidth = .2, color = "grey55", alpha = .8) +
  geom_errorbarh(data = eb_df, inherit.aes = FALSE,
                 aes(y = rna_diff, xmin = dev_diff - dev_sem, xmax = dev_diff + dev_sem),
                 height = 0, linewidth = .2, color = "grey55", alpha = .8) +
  geom_point(aes(color = color, alpha = alpha), shape = 21, stroke = 0.1) +
  scale_size_continuous(range = c(0.8, 3.5), name = "tumor RNA\n(log2 mean)") +
  scale_color_manual(values = c('FALSE' = "grey", 'TRUE' = "black")) +
  scale_fill_manual(values = c(both_low = "grey20", rna_low = "grey20",
                               dev_low = "grey20", both_high = "darkred")) +
  geom_vline(xintercept = diff_line, linetype = "dashed", color = "grey44", linewidth = .3) +
  geom_hline(yintercept = diff_line, linetype = "dashed", color = "grey44", linewidth = .3) +
  gtheme_no_rot +
  geom_text_repel(aes(color = label_guides, fill = label_guides),
                  segment.size = .05, size = 1.8, max.overlaps = 100) +
  xlab("activity difference") + ylab("RNA difference") +
  # coord_cartesian zooms WITHOUT dropping out-of-range error bars
  # (xlim()/ylim() would delete whole bars spilling past the limits).
  coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.6, 0.6))

pdf(file.path("Plots", "Diff_normal_tumor_deviation_and_rna_scatterplot5_errorbars.pdf"), width = 5, height = 4)
print(tf_diff_eb)
dev.off()

# Larger version for readability
pdf(file.path("Plots", "Diff_normal_tumor_deviation_and_rna_scatterplot5_errorbars_large.pdf"), width = 8, height = 6.5)
print(tf_diff_eb + theme(text = element_text(size = 12)) +
      geom_text_repel(aes(color = label_guides, fill = label_guides),
                      segment.size = .1, size = 3, max.overlaps = 100))
dev.off()

###############################################################
# Save tables
###############################################################
saveRDS(selected_TF, file.path(OUTDIR, "selected_TF.rds"))
write.csv(TF_diff_rna, file.path(OUTDIR, "TF_diff_rna_with_errorbars.csv"))
selected_TF_df <- cbind(dev_res_df[selected_TF, , drop = FALSE], TF_diff_rna[selected_TF, ])
write.csv(selected_TF_df, file.path(OUTDIR, "selected_TF_table_with_errorbars.csv"))

# Per-sample matrices used for the error bars (for transparency)
write.csv(as.data.frame(mMat_agg), file.path(OUTDIR, "deviation_per_sample.csv"))
write.csv(as.data.frame(ps),       file.path(OUTDIR, "rna_log2mean_per_sample.csv"))

cat("N selected TFs:", length(selected_TF), "\n")
cat("DONE\n")
