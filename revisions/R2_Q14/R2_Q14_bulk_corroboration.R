###############################################################
# R2_Q14 (bulk corroboration) : RUNX vs BAP1 in bulk RNA-seq
#
# Corroborate the scATAC finding (RUNX1/2/3 chromVAR activity is higher
# in BAP1-retained tumors / drops with low BAP1 accessibility) in three
# independent bulk mesothelioma RNA-seq cohorts used elsewhere in the
# study: Bueno, TCGA, MESOMICS.
#
# Data objects (built by git_repo/bulk_analysis.R):
#   bulk_RNA_studies.rds          -> meso_bulk_l  (log2(expr+1), gene x sample)
#   bulk_RNA_studies_metadata.rds -> meso_bulk_meta_l (subtype, etc.)
#
# Bulk has expression, not chromVAR activity, so the analog read-outs are:
#   (A) CONTINUOUS  : Spearman correlation of RUNX1/2/3 (and cofactor CBFB)
#                     expression vs BAP1 expression across samples, per
#                     cohort. scATAC predicts POSITIVE correlation
#                     (RUNX down as BAP1 down).
#   (B) SUBTYPE-ADJUSTED : lm(RUNX ~ BAP1 + subtype) BAP1 term, since
#                     histology confounds both BAP1 status and RUNX.
#   (C) GROUP TEST  :
#         - MESOMICS: real BAP1 status from IHC.BAP1 (lost vs retained)
#                     -> RUNX expression, Wilcoxon. (Direction of the IHC
#                     label is verified against BAP1 mRNA in-script.)
#         - TCGA / Bueno: no BAP1 status column -> BAP1-mRNA tertiles
#                     (low = proxy for loss) vs RUNX expression, Wilcoxon.
###############################################################

set.seed(1234)

GITROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo"
BULKDIR <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/bulkRNA_meso"
OUTDIR  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q14"
setwd(OUTDIR)
dir.create(file.path(OUTDIR, "Plots"), showWarnings = FALSE)

suppressMessages({
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(reshape2)
})
source(file.path(GITROOT, "utils", "ggplot_aestetics.R"))   # gtheme, gtheme_no_rot

bl   <- readRDS(file.path(BULKDIR, "bulk_RNA_studies.rds"))
meta <- readRDS(file.path(BULKDIR, "bulk_RNA_studies_metadata.rds"))

studies    <- c("bueno", "tcga", "mesomics")
runx_genes <- c("RUNX1", "RUNX2", "RUNX3")
tf_genes   <- c(runx_genes, "CBFB")          # include RUNX cofactor CBFB

subtype_cols <- c(Epithelioid = "#2471a3", `Biphasic-E` = "#5dade2",
                  Biphasic = "#f5b041", `Biphasic-S` = "#e59866",
                  Sarcomatoid = "#c0392b")

# ----------------------------------------------------------------
# Assemble a tidy per-sample table (BAP1 + TF expr + subtype) per cohort
# ----------------------------------------------------------------
get_tab <- function(s) {
  m  <- bl[[s]]
  md <- meta[[s]]
  # align metadata to expression columns
  if (all(colnames(m) %in% rownames(md))) {
    md <- md[colnames(m), , drop = FALSE]
  } else if ("Sample" %in% colnames(md)) {
    md <- md[match(colnames(m), md$Sample), , drop = FALSE]
  }
  df <- data.frame(sample = colnames(m),
                   cohort = s,
                   subtype = as.character(md$subtype),
                   BAP1 = as.numeric(m["BAP1", ]),
                   stringsAsFactors = FALSE)
  for (g in tf_genes) df[[g]] <- if (g %in% rownames(m)) as.numeric(m[g, ]) else NA_real_
  # mesomics BAP1 IHC status (if present)
  df$IHC_BAP1 <- if ("IHC.BAP1" %in% colnames(md)) as.character(md$IHC.BAP1) else NA_character_
  df
}
tabs <- setNames(lapply(studies, get_tab), studies)

###############################################################
# (A) + (B) : correlation & subtype-adjusted association
###############################################################
cor_rows <- list()
for (s in studies) {
  df <- tabs[[s]]
  for (g in tf_genes) {
    x <- df$BAP1; y <- df[[g]]
    ok <- is.finite(x) & is.finite(y)
    xo <- x[ok]; yo <- y[ok]; sub <- df$subtype[ok]
    ct <- suppressWarnings(cor.test(xo, yo, method = "spearman"))
    # subtype-adjusted: lm(TF ~ BAP1 + subtype); report BAP1 term
    lm_p <- NA_real_; lm_beta <- NA_real_
    if (length(unique(sub[!is.na(sub)])) > 1) {
      fit <- try(lm(yo ~ xo + factor(sub)), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        cf <- summary(fit)$coefficients
        if ("xo" %in% rownames(cf)) { lm_beta <- cf["xo","Estimate"]; lm_p <- cf["xo","Pr(>|t|)"] }
      }
    }
    cor_rows[[paste(s,g)]] <- data.frame(
      cohort = s, gene = g, n = length(xo),
      spearman_rho = unname(ct$estimate), spearman_p = ct$p.value,
      lm_BAP1_beta = lm_beta, lm_BAP1_p_adjSubtype = lm_p)
  }
}
cor_tab <- do.call(rbind, cor_rows)
cor_tab$spearman_padj <- p.adjust(cor_tab$spearman_p, method = "BH")
rownames(cor_tab) <- NULL
write.csv(cor_tab, file.path(OUTDIR, "bulk_RUNX_BAP1_correlation.csv"), row.names = FALSE)

cat("=== RUNX/CBFB vs BAP1 expression correlation (per cohort) ===\n")
print(cor_tab[, c("cohort","gene","n","spearman_rho","spearman_p","lm_BAP1_p_adjSubtype")], digits = 3)

###############################################################
# (C) MESOMICS : BAP1 IHC status group test  (real BAP1 loss)
###############################################################
msm <- tabs[["mesomics"]]
# map IHC label -> lost / retained; drop mixed + NA
msm$BAP1_status <- NA_character_
msm$BAP1_status[msm$IHC_BAP1 == "YES"] <- "prelim_YES"
msm$BAP1_status[msm$IHC_BAP1 == "NO"]  <- "prelim_NO"
msm_g <- msm[!is.na(msm$BAP1_status), ]

# verify direction against BAP1 mRNA: the group with LOWER BAP1 mRNA is "lost"
bap1_by_lab <- tapply(msm_g$BAP1, msm_g$BAP1_status, mean, na.rm = TRUE)
lost_label  <- names(which.min(bap1_by_lab))     # lower BAP1 mRNA = lost
cat("\n=== MESOMICS IHC.BAP1 direction check (mean BAP1 mRNA) ===\n")
print(round(bap1_by_lab, 3))
cat("=> IHC label with lower BAP1 mRNA (=BAP1 lost):", lost_label,
    "( YES=", round(bap1_by_lab["prelim_YES"],3),
    "NO=", round(bap1_by_lab["prelim_NO"],3), ")\n")
msm_g$BAP1_status <- ifelse(msm_g$BAP1_status == lost_label, "BAP1-lost", "BAP1-retained")
msm_g$BAP1_status <- factor(msm_g$BAP1_status, levels = c("BAP1-retained","BAP1-lost"))
cat("MESOMICS group sizes:\n"); print(table(msm_g$BAP1_status))

ihc_rows <- list()
for (g in tf_genes) {
  wt <- suppressWarnings(wilcox.test(msm_g[[g]] ~ msm_g$BAP1_status))
  ihc_rows[[g]] <- data.frame(
    gene = g,
    mean_retained = mean(msm_g[[g]][msm_g$BAP1_status=="BAP1-retained"], na.rm=TRUE),
    mean_lost     = mean(msm_g[[g]][msm_g$BAP1_status=="BAP1-lost"],     na.rm=TRUE),
    wilcox_p = wt$p.value)
}
ihc_tab <- do.call(rbind, ihc_rows)
ihc_tab$delta_lost_minus_retained <- ihc_tab$mean_lost - ihc_tab$mean_retained
write.csv(ihc_tab, file.path(OUTDIR, "mesomics_IHC_BAP1_RUNX_grouptest.csv"), row.names = FALSE)
cat("\n=== MESOMICS IHC-BAP1 RUNX group test ===\n"); print(ihc_tab, digits = 3)

###############################################################
# (C') TCGA / Bueno : BAP1-mRNA tertile group test (proxy)
###############################################################
tert_rows <- list()
for (s in c("tcga","bueno")) {
  df <- tabs[[s]]
  q <- quantile(df$BAP1, c(1/3, 2/3), na.rm = TRUE)
  df$BAP1_tert <- ifelse(df$BAP1 <= q[1], "BAP1-low",
                  ifelse(df$BAP1 >= q[2], "BAP1-high", NA))
  d2 <- df[!is.na(df$BAP1_tert), ]
  d2$BAP1_tert <- factor(d2$BAP1_tert, levels = c("BAP1-high","BAP1-low"))
  tabs[[s]] <- df
  for (g in tf_genes) {
    wt <- suppressWarnings(wilcox.test(d2[[g]] ~ d2$BAP1_tert))
    tert_rows[[paste(s,g)]] <- data.frame(
      cohort = s, gene = g,
      mean_BAP1high = mean(d2[[g]][d2$BAP1_tert=="BAP1-high"], na.rm=TRUE),
      mean_BAP1low  = mean(d2[[g]][d2$BAP1_tert=="BAP1-low"],  na.rm=TRUE),
      wilcox_p = wt$p.value)
  }
}
tert_tab <- do.call(rbind, tert_rows)
tert_tab$delta_low_minus_high <- tert_tab$mean_BAP1low - tert_tab$mean_BAP1high
write.csv(tert_tab, file.path(OUTDIR, "tcga_bueno_BAP1tertile_RUNX_grouptest.csv"), row.names = FALSE)
cat("\n=== TCGA/Bueno BAP1-mRNA tertile RUNX group test ===\n"); print(tert_tab, digits = 3)

###############################################################
# PLOTS
###############################################################
cohort_titles <- c(bueno = "Bueno", tcga = "TCGA", mesomics = "MESOMICS")

## Fig 1 : RUNX (+CBFB) vs BAP1 expression scatter, facet gene x cohort ----
long <- do.call(rbind, lapply(studies, function(s) {
  df <- tabs[[s]]
  do.call(rbind, lapply(tf_genes, function(g) data.frame(
    cohort = cohort_titles[s], gene = g,
    BAP1 = df$BAP1, TF = df[[g]], subtype = df$subtype)))
}))
long$cohort <- factor(long$cohort, levels = cohort_titles[studies])
long$gene   <- factor(long$gene, levels = tf_genes)
long <- long[is.finite(long$BAP1) & is.finite(long$TF), ]

p_scatter <- ggplot(long, aes(x = BAP1, y = TF)) +
  geom_point(aes(color = subtype), size = .7, alpha = .75) +
  geom_smooth(method = "lm", se = FALSE, color = "grey25", linewidth = .4) +
  stat_cor(method = "spearman", size = 2.4, color = "black",
           label.x.npc = "left", label.y.npc = "top") +
  facet_grid(gene ~ cohort, scales = "free") +
  scale_color_manual(values = subtype_cols, na.value = "grey70", name = "subtype") +
  gtheme_no_rot +
  xlab("BAP1 expression (log2)") + ylab("TF expression (log2)") +
  ggtitle("RUNX / CBFB vs BAP1 expression across bulk mesothelioma cohorts")

pdf(file.path("Plots", "bulk_RUNX_vs_BAP1_scatter.pdf"), width = 8.5, height = 9)
print(p_scatter); dev.off()

## Fig 2 : MESOMICS IHC-BAP1 status boxplots (RUNX1/2/3 + BAP1 validation) ----
box_df <- do.call(rbind, lapply(c("BAP1", runx_genes), function(g)
  data.frame(gene = g, value = msm_g[[g]], status = msm_g$BAP1_status)))
box_df$gene <- factor(box_df$gene, levels = c("BAP1", runx_genes))

p_ihc <- ggplot(box_df, aes(status, value, fill = status)) +
  geom_boxplot(width = .55, outlier.shape = NA, alpha = .55) +
  geom_jitter(width = .12, size = .6, alpha = .5) +
  stat_compare_means(method = "wilcox.test", size = 2.4,
                     label = "p.format", label.x.npc = "center") +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(`BAP1-retained` = "#2471a3", `BAP1-lost` = "#c0392b"),
                    guide = "none") +
  gtheme_no_rot + xlab("") + ylab("expression (log2)") +
  ggtitle("MESOMICS: RUNX expression by BAP1 IHC status (BAP1 = validation)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf(file.path("Plots", "mesomics_IHC_BAP1_RUNX_boxplots.pdf"), width = 8, height = 3.4)
print(p_ihc); dev.off()

## Fig 3 : summary heatmap of Spearman rho (gene x cohort) ----
hm <- cor_tab
hm$sig <- ifelse(hm$spearman_p < 0.001, "***",
          ifelse(hm$spearman_p < 0.01, "**",
          ifelse(hm$spearman_p < 0.05, "*", "")))
hm$cohort <- factor(cohort_titles[hm$cohort], levels = cohort_titles[studies])
hm$gene   <- factor(hm$gene, levels = rev(tf_genes))

p_heat <- ggplot(hm, aes(cohort, gene, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = .6) +
  geom_text(aes(label = sprintf("%.2f%s", spearman_rho, sig)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-1,1), name = "Spearman rho\n(TF vs BAP1)") +
  gtheme_no_rot + xlab("") + ylab("") +
  ggtitle("RUNX/CBFB ~ BAP1 expression correlation\n(* p<0.05, ** p<0.01, *** p<0.001)")

pdf(file.path("Plots", "bulk_RUNX_BAP1_correlation_heatmap.pdf"), width = 4.6, height = 3.6)
print(p_heat); dev.off()

cat("\nDONE\n")
