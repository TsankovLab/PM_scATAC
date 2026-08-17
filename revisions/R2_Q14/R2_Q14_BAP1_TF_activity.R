###############################################################
# R2_Q14 : TF activity (chromVAR) vs BAP1 status in the
#          malignant / tumor compartment.
#
# Reviewer 2, Q14: explore, in a systematic manner, whether TF
# activity differs between tumors WITH vs WITHOUT BAP1 loss, and
# whether TF activity increases/decreases in relation to low BAP1
# (gene score). Reported trend: RUNX1/RUNX2 higher in BAP1-retained.
#
# BAP1 status (per-sample genetic annotation):
#   LOST     : P3, P4, P5, P8, P10          (n = 5)
#   RETAINED : P1, P11, P12, P13, P14, P23  (n = 6)
#
# Data: ArchR object ../tumor_compartment/scatac_ArchR
#   - MotifMatrix  (chromVAR deviations + z)  -> TF activity
#   - GeneScoreMatrix (BAP1)                  -> BAP1 accessibility
# Sample label is in the `Sample` cellColData column.
#
# Two complementary, systematic read-outs:
#   (A) GROUP TEST  : per-sample pseudobulk TF activity, LOST vs
#                     RETAINED, sample as the unit of replication
#                     (Wilcoxon + t-test, BH-adjusted). This avoids
#                     the pseudoreplication that a per-cell test
#                     suffers here (P23 alone = 11,687 cells).
#   (B) CONTINUOUS  : per-sample correlation of each TF's activity
#                     with the sample's BAP1 gene score (Spearman).
#                     Directly answers "increase/decrease of TF
#                     activity in relation to low BAP1".
# A per-cell presto test is also run as a supplement (with caveat).
###############################################################

set.seed(1234)

GITROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo"
BASE    <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment"
OUTDIR  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q14"
setwd(OUTDIR)
dir.create(file.path(OUTDIR, "Plots"), showWarnings = FALSE)

suppressMessages({
  library(ArchR)
  library(SummarizedExperiment)
  library(ggplot2)
  library(ggrepel)
  library(presto)
})
source(file.path(GITROOT, "utils", "scATAC_functions.R"))   # fetch_mat + hidden ArchR fns
source(file.path(GITROOT, "utils", "ggplot_aestetics.R"))   # gtheme_no_rot

addArchRThreads(threads = 1)
addArchRGenome("hg38")

# ----------------------------------------------------------------
# Sample groups (BAP1 status from per-sample genetic annotation)
# ----------------------------------------------------------------
bap1_lost     <- c('P3', 'P4', 'P5', 'P8', 'P10')
bap1_retained <- c('P1', 'P11', 'P12', 'P13', 'P14', 'P23')
tumor_sams    <- c(bap1_lost, bap1_retained)

# ----------------------------------------------------------------
# Load project, restrict to the 11 annotated tumor samples
# ----------------------------------------------------------------
archp <- loadArchRProject(file.path(BASE, "scatac_ArchR"), showLogo = FALSE)
cc    <- as.data.frame(archp@cellColData)
cc$Sample <- as.character(cc$Sample)

keep_cells <- rownames(cc)[cc$Sample %in% tumor_sams]
archp <- archp[keep_cells, ]
cc    <- as.data.frame(archp@cellColData)
cc$Sample <- as.character(cc$Sample)
cc$BAP1_status <- ifelse(cc$Sample %in% bap1_lost, "BAP1_lost", "BAP1_retained")

cat("Cells per sample (tumor only):\n"); print(sort(table(cc$Sample)))
cat("\nCells per BAP1 status:\n");        print(table(cc$BAP1_status))

###############################################################
# 1. TF activity (chromVAR) : per-cell deviations, per-sample pseudobulk
#
# NOTE on the metric: this matches discover_oncogenic_TFs.R, where the
# plotted "TF activity" (dev_diff) is the RAW chromVAR *deviations*
# (assays(mSE)[[1]]) aggregated by mean per sample -- NOT the chromVAR z
# assay and NOT base-R scale()'d. Scaling in the original entered only the
# per-cell wilcoxauc significance step (mMat <- scale(mMat)); we mirror that
# exactly in the per-cell supplement (section 5) and leave the sample-level
# pseudobulk means unscaled, identical to the scatter axis.
###############################################################
mSE   <- fetch_mat(archp, 'Motif')                 # rowData$name cleaned
mDev  <- as.matrix(assays(mSE)[[1]])               # chromVAR deviations (assay 1), dense
rownames(mDev) <- rowData(mSE)$name
mDev  <- mDev[, archp$cellNames, drop = FALSE]
stopifnot(all(colnames(mDev) == rownames(cc)))

# drop TFs with any non-finite value across cells
finite_tf <- rownames(mDev)[apply(mDev, 1, function(x) all(is.finite(x)))]
mDev <- mDev[finite_tf, , drop = FALSE]
cat("\nN TFs (finite deviations):", nrow(mDev), "\n")

# per-sample pseudobulk = mean deviation across cells of that sample (TF x sample)
samp_of_cell <- cc$Sample
act_by_samp  <- t(apply(mDev, 1, function(v) tapply(v, samp_of_cell, mean)))
act_by_samp  <- act_by_samp[, tumor_sams, drop = FALSE]

###############################################################
# 2. BAP1 gene score : per-cell -> per-sample pseudobulk
###############################################################
gsSE  <- fetch_mat(archp, 'GeneScore')
gsMat <- as.matrix(assays(gsSE)[[1]])
rownames(gsMat) <- rowData(gsSE)$name
gsMat <- gsMat[, archp$cellNames, drop = FALSE]
stopifnot(all(colnames(gsMat) == rownames(cc)))

bap1_gs_cell <- as.numeric(gsMat["BAP1", ])
bap1_gs_samp <- tapply(bap1_gs_cell, samp_of_cell, mean)[tumor_sams]

cat("\nBAP1 gene score per sample (mean):\n")
print(round(sort(bap1_gs_samp), 3))

###############################################################
# 3. (A) GROUP TEST : LOST vs RETAINED at the SAMPLE level
###############################################################
grp <- ifelse(tumor_sams %in% bap1_lost, "lost", "retained")
names(grp) <- tumor_sams

group_stats <- data.frame(
  TF            = rownames(act_by_samp),
  mean_lost     = rowMeans(act_by_samp[, bap1_lost,     drop = FALSE]),
  mean_retained = rowMeans(act_by_samp[, bap1_retained, drop = FALSE]),
  stringsAsFactors = FALSE
)
# delta > 0  => higher activity in BAP1-LOST tumors
group_stats$delta_lost_minus_retained <- group_stats$mean_lost - group_stats$mean_retained

# per-TF two-group tests (sample as unit)
grp_vec <- factor(grp[colnames(act_by_samp)], levels = c("retained", "lost"))
p_wilcox <- apply(act_by_samp, 1, function(v)
  tryCatch(wilcox.test(v ~ grp_vec)$p.value, error = function(e) NA_real_))
p_ttest  <- apply(act_by_samp, 1, function(v)
  tryCatch(t.test(v ~ grp_vec)$p.value,      error = function(e) NA_real_))

group_stats$p_wilcox    <- p_wilcox[group_stats$TF]
group_stats$p_ttest     <- p_ttest[group_stats$TF]
group_stats$padj_wilcox <- p.adjust(group_stats$p_wilcox, method = "BH")
group_stats$padj_ttest  <- p.adjust(group_stats$p_ttest,  method = "BH")
group_stats <- group_stats[order(group_stats$p_ttest), ]

###############################################################
# 4. (B) CONTINUOUS : correlate TF activity with BAP1 gene score
#        across the 11 tumor samples (Spearman).
#        rho < 0  => TF activity goes UP as BAP1 goes DOWN.
###############################################################
cor_stats <- data.frame(TF = rownames(act_by_samp), stringsAsFactors = FALSE)
cs <- t(apply(act_by_samp, 1, function(v) {
  ct <- suppressWarnings(cor.test(v, bap1_gs_samp[colnames(act_by_samp)],
                                  method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value)
}))
cor_stats$rho_bap1  <- cs[cor_stats$TF, "rho"]
cor_stats$p_spearman <- cs[cor_stats$TF, "p"]
cor_stats$padj_spearman <- p.adjust(cor_stats$p_spearman, method = "BH")
# also Pearson for reference
csp <- t(apply(act_by_samp, 1, function(v) {
  ct <- suppressWarnings(cor.test(v, bap1_gs_samp[colnames(act_by_samp)],
                                  method = "pearson"))
  c(r = unname(ct$estimate), p = ct$p.value)
}))
cor_stats$pearson_r <- csp[cor_stats$TF, "r"]
cor_stats$pearson_p <- csp[cor_stats$TF, "p"]
cor_stats <- cor_stats[order(cor_stats$p_spearman), ]

###############################################################
# 5. Supplement : per-CELL presto Wilcoxon (LOST vs RETAINED)
#    Flagged as pseudoreplicated (cells not independent within a
#    tumor); provided because it mirrors the reviewer's observation.
###############################################################
# mirror discover_oncogenic_TFs.R exactly: base-R scale() on the deviations
# (column-wise, i.e. per cell) before the per-cell Wilcoxon.
mDev_scaled <- scale(mDev)
mDev_scaled <- mDev_scaled[apply(mDev_scaled, 1, function(x) all(is.finite(x))), , drop = FALSE]
cell_res <- wilcoxauc(mDev_scaled, y = cc$BAP1_status)
cell_res <- cell_res[cell_res$group == "BAP1_lost", ]   # logFC>0 => higher in lost
cell_res <- cell_res[order(cell_res$pval),
                     c("feature","auc","logFC","pval","padj")]
colnames(cell_res) <- c("TF","auc_cell","logFC_cell","pval_cell","padj_cell")

###############################################################
# 6. Assemble a master table + RUNX highlight
###############################################################
runx <- c("RUNX1","RUNX2","RUNX3")
master <- merge(group_stats, cor_stats, by = "TF")
master <- merge(master, cell_res, by = "TF", all.x = TRUE)
master$is_runx <- master$TF %in% runx
master <- master[order(master$p_ttest), ]

write.csv(master,       file.path(OUTDIR, "BAP1_TF_activity_master_table.csv"), row.names = FALSE)
write.csv(group_stats,  file.path(OUTDIR, "group_test_lost_vs_retained.csv"),   row.names = FALSE)
write.csv(cor_stats,    file.path(OUTDIR, "correlation_TFactivity_vs_BAP1genescore.csv"), row.names = FALSE)
write.csv(as.data.frame(act_by_samp), file.path(OUTDIR, "TFactivity_deviation_per_sample.csv"))
write.csv(data.frame(Sample = names(bap1_gs_samp),
                     BAP1_status = ifelse(names(bap1_gs_samp) %in% bap1_lost, "lost", "retained"),
                     BAP1_genescore = as.numeric(bap1_gs_samp)),
          file.path(OUTDIR, "BAP1_genescore_per_sample.csv"), row.names = FALSE)

cat("\n=== RUNX summary ===\n")
print(master[master$is_runx,
      c("TF","mean_lost","mean_retained","delta_lost_minus_retained",
        "p_ttest","p_wilcox","rho_bap1","p_spearman")])

###############################################################
# 7. PLOTS
###############################################################
lab_thr_p   <- 0.05      # label TFs below this (group test)
lab_top_n   <- 25

## 7a. Volcano-style: group delta vs -log10(p_ttest) --------------
gp <- group_stats
gp$neglogp <- -log10(gp$p_ttest)
gp$sig <- ifelse(gp$p_ttest < lab_thr_p,
                 ifelse(gp$delta_lost_minus_retained > 0, "up_in_lost", "up_in_retained"),
                 "ns")
gp$sig[gp$TF %in% runx] <- "RUNX"
gp$is_runx <- gp$TF %in% runx
# label: significant ones (top by p) + RUNX
lab_set <- unique(c(head(gp$TF[gp$p_ttest < lab_thr_p], lab_top_n), runx))
gp$lab  <- ifelse(gp$TF %in% lab_set, gp$TF, "")

p_volcano <- ggplot(gp, aes(x = delta_lost_minus_retained, y = neglogp)) +
  geom_hline(yintercept = -log10(lab_thr_p), linetype = "dashed", color = "grey60", linewidth = .3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = .3) +
  geom_point(aes(color = sig, size = is_runx), alpha = .8, stroke = 0) +
  scale_color_manual(values = c(up_in_lost = "#c0392b", up_in_retained = "#2471a3",
                                ns = "grey75", RUNX = "black"),
                     name = NULL) +
  scale_size_manual(values = c(`FALSE` = 1.1, `TRUE` = 2.6), guide = "none") +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.4, max.overlaps = 100,
                           segment.size = .2, min.segment.length = 0) +
  gtheme_no_rot +
  xlab("TF activity difference  (BAP1-lost  -  BAP1-retained)") +
  ylab("-log10  p (t-test, sample-level)") +
  ggtitle("TF activity: BAP1-lost vs BAP1-retained tumors")

pdf(file.path("Plots", "volcano_group_lost_vs_retained.pdf"), width = 7, height = 5.5)
print(p_volcano); dev.off()

## 7b. Correlation with BAP1 gene score: rho vs -log10(p) ----------
cp <- cor_stats
cp$neglogp <- -log10(cp$p_spearman)
cp$sig <- ifelse(cp$p_spearman < lab_thr_p,
                 ifelse(cp$rho_bap1 < 0, "up_as_BAP1_down", "down_as_BAP1_down"), "ns")
cp$sig[cp$TF %in% runx] <- "RUNX"
cp$is_runx <- cp$TF %in% runx
lab_set2 <- unique(c(head(cp$TF[cp$p_spearman < lab_thr_p], lab_top_n), runx))
cp$lab   <- ifelse(cp$TF %in% lab_set2, cp$TF, "")

p_cor <- ggplot(cp, aes(x = rho_bap1, y = neglogp)) +
  geom_hline(yintercept = -log10(lab_thr_p), linetype = "dashed", color = "grey60", linewidth = .3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = .3) +
  geom_point(aes(color = sig, size = is_runx), alpha = .8, stroke = 0) +
  scale_color_manual(values = c(up_as_BAP1_down = "#c0392b", down_as_BAP1_down = "#2471a3",
                                ns = "grey75", RUNX = "black"), name = NULL) +
  scale_size_manual(values = c(`FALSE` = 1.1, `TRUE` = 2.6), guide = "none") +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.4, max.overlaps = 100,
                           segment.size = .2, min.segment.length = 0) +
  gtheme_no_rot +
  xlab("Spearman rho: TF activity vs BAP1 gene score (across tumor samples)") +
  ylab("-log10  p (Spearman)") +
  ggtitle("TF activity in relation to BAP1 accessibility\n(rho < 0: activity up as BAP1 down)")

pdf(file.path("Plots", "correlation_TFactivity_vs_BAP1.pdf"), width = 7, height = 5.5)
print(p_cor); dev.off()

## 7c. BAP1 gene score by status (validation) ---------------------
bdf <- data.frame(Sample = names(bap1_gs_samp),
                  BAP1_genescore = as.numeric(bap1_gs_samp),
                  status = factor(ifelse(names(bap1_gs_samp) %in% bap1_lost,
                                         "BAP1-lost", "BAP1-retained"),
                                  levels = c("BAP1-retained","BAP1-lost")))
p_bap1 <- ggplot(bdf, aes(status, BAP1_genescore, fill = status)) +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = .5) +
  geom_jitter(width = .12, size = 2) +
  ggrepel::geom_text_repel(aes(label = Sample), size = 2.6, max.overlaps = 100) +
  scale_fill_manual(values = c(`BAP1-retained` = "#2471a3", `BAP1-lost` = "#c0392b"),
                    guide = "none") +
  gtheme_no_rot + xlab("") + ylab("BAP1 gene score (mean per sample)") +
  ggtitle("BAP1 accessibility validates genetic status")

pdf(file.path("Plots", "BAP1_genescore_by_status.pdf"), width = 4.2, height = 4)
print(p_bap1); dev.off()

## 7d. RUNX focused panels ---------------------------------------
# per-sample activity by status
runx_present <- intersect(runx, rownames(act_by_samp))
rdf <- do.call(rbind, lapply(runx_present, function(tf) {
  data.frame(TF = tf, Sample = colnames(act_by_samp),
             activity = as.numeric(act_by_samp[tf, ]),
             bap1_gs  = as.numeric(bap1_gs_samp[colnames(act_by_samp)]),
             status = ifelse(colnames(act_by_samp) %in% bap1_lost, "BAP1-lost", "BAP1-retained"),
             stringsAsFactors = FALSE)
}))
rdf$status <- factor(rdf$status, levels = c("BAP1-retained","BAP1-lost"))

p_runx_box <- ggplot(rdf, aes(status, activity, fill = status)) +
  geom_boxplot(width = .55, outlier.shape = NA, alpha = .5) +
  geom_jitter(width = .12, size = 1.8) +
  facet_wrap(~ TF, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(`BAP1-retained` = "#2471a3", `BAP1-lost` = "#c0392b"),
                    guide = "none") +
  gtheme_no_rot + xlab("") + ylab("TF activity (mean chromVAR deviation per sample)") +
  ggtitle("RUNX TF activity by BAP1 status") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf(file.path("Plots", "RUNX_activity_by_status.pdf"), width = 7, height = 3.6)
print(p_runx_box); dev.off()

p_runx_scatter <- ggplot(rdf, aes(bap1_gs, activity, color = status)) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              color = "grey55", linewidth = .4, linetype = "dashed") +
  geom_point(size = 2.2) +
  ggrepel::geom_text_repel(aes(label = Sample), size = 2.3, max.overlaps = 100, show.legend = FALSE) +
  facet_wrap(~ TF, nrow = 1, scales = "free") +
  scale_color_manual(values = c(`BAP1-retained` = "#2471a3", `BAP1-lost` = "#c0392b"), name = NULL) +
  gtheme_no_rot + xlab("BAP1 gene score (per sample)") +
  ylab("TF activity (chromVAR deviation)") +
  ggtitle("RUNX activity vs BAP1 accessibility across tumors")

pdf(file.path("Plots", "RUNX_activity_vs_BAP1genescore.pdf"), width = 8.5, height = 3.6)
print(p_runx_scatter); dev.off()

cat("\nDONE\n")
