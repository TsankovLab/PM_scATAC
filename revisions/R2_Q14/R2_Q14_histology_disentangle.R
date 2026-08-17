###############################################################
# R2_Q14 : Disentangle BAP1-status effect from HISTOLOGY
#
# Concern: in mesothelioma BAP1 loss is enriched in EPITHELIOID and rare
# in SARCOMATOID tumors, while RUNX1/RUNX2 are mesenchymal/EMT TFs higher
# in sarcomatoid. So "RUNX low when BAP1 lost" could be a pure histology
# artifact. This script (1) quantifies the BAP1<->histology association
# and (2) tests whether the RUNX~BAP1 relationship survives adjustment for
# histology, in both the bulk cohorts and the scATAC samples.
#
# Histology is represented two ways, consistently with the study:
#   - discrete subtype (Epithelioid / Biphasic / Sarcomatoid)
#   - continuous SARCOMATOID SCORE = mean expression of the cNMF20 program
#     (`sarc_score` already in the bulk metadata; per-sample scores for the
#      scATAC samples in cnmf20_sarcomatoid_sample_order.csv)
###############################################################

set.seed(1234)

GITROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo"
BULKDIR <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/bulkRNA_meso"
TCOMP   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment"
OUTDIR  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q14"
setwd(OUTDIR)
dir.create(file.path(OUTDIR, "Plots"), showWarnings = FALSE)

suppressMessages({ library(ggplot2); library(ggpubr); library(patchwork) })
source(file.path(GITROOT, "utils", "ggplot_aestetics.R"))

bl   <- readRDS(file.path(BULKDIR, "bulk_RNA_studies.rds"))
meta <- readRDS(file.path(BULKDIR, "bulk_RNA_studies_metadata.rds"))
studies    <- c("bueno", "tcga", "mesomics")
runx_genes <- c("RUNX1", "RUNX2", "RUNX3")
tf_genes   <- c(runx_genes, "CBFB")
subtype_cols <- c(Epithelioid = "#2471a3", `Biphasic-E` = "#5dade2",
                  Biphasic = "#f5b041", `Biphasic-S` = "#e59866",
                  Sarcomatoid = "#c0392b")

# partial Spearman (rho of x,y controlling z): Pearson corr of rank residuals
partial_spearman <- function(x, y, z) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  if (sum(ok) < 5) return(c(rho = NA, p = NA, n = sum(ok)))
  rx <- rank(x[ok]); ry <- rank(y[ok]); rz <- rank(z[ok])
  ex <- resid(lm(rx ~ rz)); ey <- resid(lm(ry ~ rz))
  ct <- suppressWarnings(cor.test(ex, ey, method = "pearson"))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

# ----------------------------------------------------------------
# Tidy per-sample tables incl subtype + sarc_score
# ----------------------------------------------------------------
get_tab <- function(s) {
  m  <- bl[[s]]; md <- meta[[s]]
  if (all(colnames(m) %in% rownames(md))) md <- md[colnames(m), , drop = FALSE] else
    if ("Sample" %in% colnames(md)) md <- md[match(colnames(m), md$Sample), , drop = FALSE]
  df <- data.frame(sample = colnames(m), cohort = s,
                   subtype = as.character(md$subtype),
                   BAP1 = as.numeric(m["BAP1", ]), stringsAsFactors = FALSE)
  for (g in tf_genes) df[[g]] <- if (g %in% rownames(m)) as.numeric(m[g, ]) else NA_real_
  df$IHC_BAP1  <- if ("IHC.BAP1" %in% colnames(md)) as.character(md$IHC.BAP1) else NA_character_
  # sarc_score: use metadata column if present, else compute from cNMF20 program
  if ("sarc_score" %in% colnames(md)) df$sarc_score <- as.numeric(md$sarc_score) else df$sarc_score <- NA_real_
  df
}
tabs <- setNames(lapply(studies, get_tab), studies)

# fallback: compute sarc_score from cNMF20 gene list if any cohort lacks it
if (any(sapply(tabs, function(d) all(is.na(d$sarc_score))))) {
  nmf <- read.csv(file.path(TCOMP, "scrna", "cnmf_genelist_25_nfeat_5000.csv"),
                  stringsAsFactors = FALSE)
  sarc_genes <- head(nmf[["cNMF20"]][!is.na(nmf[["cNMF20"]])], 20)
  for (s in studies) if (all(is.na(tabs[[s]]$sarc_score))) {
    m <- bl[[s]]; g <- intersect(sarc_genes, rownames(m))
    tabs[[s]]$sarc_score <- colMeans(m[g, , drop = FALSE], na.rm = TRUE)
  }
}

###############################################################
# PART 1 : Does BAP1 loss track with histology?  (bulk)
###############################################################
p1_rows <- list()
for (s in studies) {
  df <- tabs[[s]]
  # BAP1 mRNA vs subtype
  kw <- suppressWarnings(kruskal.test(BAP1 ~ factor(subtype), data = df))$p.value
  epi <- df$BAP1[df$subtype == "Epithelioid"]; sar <- df$BAP1[df$subtype == "Sarcomatoid"]
  ews_p <- if (length(epi) >= 3 & length(sar) >= 3)
    suppressWarnings(wilcox.test(epi, sar))$p.value else NA_real_
  # BAP1 mRNA vs continuous sarc score
  cs <- suppressWarnings(cor.test(df$BAP1, df$sarc_score, method = "spearman"))
  p1_rows[[s]] <- data.frame(cohort = s,
    BAP1_vs_subtype_KW_p = kw,
    mean_BAP1_epithelioid = mean(epi, na.rm = TRUE),
    mean_BAP1_sarcomatoid = mean(sar, na.rm = TRUE),
    epi_vs_sarc_wilcox_p = ews_p,
    BAP1_vs_sarcscore_rho = unname(cs$estimate),
    BAP1_vs_sarcscore_p = cs$p.value)
}
# MESOMICS: real IHC status vs subtype (Fisher) and vs sarc score
msm <- tabs[["mesomics"]]
msm$BAP1_stat <- ifelse(msm$IHC_BAP1 == "NO", "BAP1-lost",
                 ifelse(msm$IHC_BAP1 == "YES", "BAP1-retained", NA))  # NO=lost (verified by mRNA)
msm_g <- msm[!is.na(msm$BAP1_stat) & !is.na(msm$subtype), ]
ihc_subtype_tab <- table(msm_g$BAP1_stat, msm_g$subtype)
fisher_p_msm <- suppressWarnings(fisher.test(ihc_subtype_tab)$p.value)

p1 <- do.call(rbind, p1_rows)
write.csv(p1, file.path(OUTDIR, "histology_BAP1_association.csv"), row.names = FALSE)
cat("=== PART 1: BAP1 vs histology (bulk) ===\n"); print(p1, digits = 3)
cat("\nMESOMICS IHC-BAP1 status x subtype:\n"); print(ihc_subtype_tab)
cat("Fisher p (IHC-BAP1 vs subtype):", signif(fisher_p_msm, 3), "\n")

###############################################################
# PART 2 : Does RUNX~BAP1 survive histology adjustment?  (bulk)
###############################################################
p2_rows <- list()
for (s in studies) {
  df <- tabs[[s]]
  epi <- df[df$subtype == "Epithelioid", ]
  for (g in tf_genes) {
    ov  <- suppressWarnings(cor.test(df$BAP1, df[[g]], method = "spearman"))
    ep  <- if (nrow(epi) >= 8)
      suppressWarnings(cor.test(epi$BAP1, epi[[g]], method = "spearman")) else
      list(estimate = NA, p.value = NA)
    ps  <- partial_spearman(df[[g]], df$BAP1, df$sarc_score)  # RUNX~BAP1 | sarc
    # lm adjusting continuous sarc score, and interaction test
    fit  <- try(lm(df[[g]] ~ df$BAP1 + df$sarc_score), silent = TRUE)
    fiti <- try(lm(df[[g]] ~ df$BAP1 * df$sarc_score), silent = TRUE)
    lm_b <- lm_p <- int_p <- NA_real_
    if (!inherits(fit,"try-error")) { cf <- summary(fit)$coefficients
      if ("df$BAP1" %in% rownames(cf)) { lm_b <- cf["df$BAP1","Estimate"]; lm_p <- cf["df$BAP1","Pr(>|t|)"] } }
    if (!inherits(fiti,"try-error")) { cfi <- summary(fiti)$coefficients
      rn <- grep(":", rownames(cfi), value = TRUE); if (length(rn)) int_p <- cfi[rn[1],"Pr(>|t|)"] }
    p2_rows[[paste(s,g)]] <- data.frame(cohort = s, gene = g,
      rho_overall = unname(ov$estimate), p_overall = ov$p.value,
      n_epithelioid = nrow(epi), rho_epithelioid = unname(ep$estimate), p_epithelioid = ep$p.value,
      partial_rho_sarc = unname(ps["rho"]), partial_p_sarc = unname(ps["p"]),
      lm_BAP1_beta_adjSarc = lm_b, lm_BAP1_p_adjSarc = lm_p, interaction_p = int_p)
  }
}
p2 <- do.call(rbind, p2_rows); rownames(p2) <- NULL
write.csv(p2, file.path(OUTDIR, "RUNX_BAP1_histology_adjusted.csv"), row.names = FALSE)
cat("\n=== PART 2: RUNX~BAP1 histology-adjusted (bulk) ===\n")
print(p2[, c("cohort","gene","rho_overall","rho_epithelioid","partial_rho_sarc","lm_BAP1_p_adjSarc")], digits = 3)

###############################################################
# PART 3 : same disentangling in the scATAC samples
#   reads the saved per-sample activity + BAP1 gene score + sarc score
###############################################################
act <- read.csv(file.path(OUTDIR, "TFactivity_deviation_per_sample.csv"),
                row.names = 1, check.names = FALSE)
gs  <- read.csv(file.path(OUTDIR, "BAP1_genescore_per_sample.csv"), stringsAsFactors = FALSE)
sc  <- read.csv(file.path(TCOMP, "scrna", "cnmf20_sarcomatoid_sample_order.csv"),
                stringsAsFactors = FALSE)          # columns: sampleID, x (sarc score)
sc  <- setNames(sc$x, sc$sampleID)

sa  <- data.frame(sample = gs$Sample, BAP1_status = gs$BAP1_status,
                  BAP1_gs = gs$BAP1_genescore, stringsAsFactors = FALSE)
sa$sarc_score <- sc[sa$sample]
for (g in runx_genes) sa[[g]] <- if (g %in% rownames(act)) as.numeric(act[g, sa$sample]) else NA_real_
sa_c <- sa[is.finite(sa$sarc_score), ]           # samples with a sarc score (drops P10,P23)

cat("\n=== PART 3: scATAC samples (sarc score available) ===\n")
print(sa_c[, c("sample","BAP1_status","BAP1_gs","sarc_score", runx_genes)], digits = 3)

# does BAP1 status / gene score track sarc score in scATAC?
bap1_sarc <- suppressWarnings(cor.test(sa_c$BAP1_gs, sa_c$sarc_score, method = "spearman"))
sarc_by_status <- tapply(sa_c$sarc_score, sa_c$BAP1_status, mean)
cat("\nBAP1 gene score vs sarc score (scATAC): rho =", round(bap1_sarc$estimate,3),
    " p =", signif(bap1_sarc$p.value,3), "\n")
cat("mean sarc score by BAP1 status:\n"); print(round(sarc_by_status,3))

p3_rows <- list()
for (g in runx_genes) {
  ov <- suppressWarnings(cor.test(sa_c$BAP1_gs, sa_c[[g]], method = "spearman"))
  ps <- partial_spearman(sa_c[[g]], sa_c$BAP1_gs, sa_c$sarc_score)
  p3_rows[[g]] <- data.frame(gene = g, n = nrow(sa_c),
    rho_activity_vs_BAP1gs = unname(ov$estimate), p = ov$p.value,
    partial_rho_sarc = unname(ps["rho"]), partial_p_sarc = unname(ps["p"]))
}
p3 <- do.call(rbind, p3_rows); rownames(p3) <- NULL
write.csv(p3, file.path(OUTDIR, "scATAC_histology_disentangle.csv"), row.names = FALSE)
write.csv(sa, file.path(OUTDIR, "scATAC_per_sample_RUNX_BAP1_sarc.csv"), row.names = FALSE)
cat("\nscATAC RUNX activity ~ BAP1 gene score, raw vs partial(|sarc):\n"); print(p3, digits = 3)

###############################################################
# PLOTS
###############################################################
cohort_titles <- c(bueno = "Bueno", tcga = "TCGA", mesomics = "MESOMICS")

## P1a: BAP1 mRNA by subtype (does BAP1 track histology) ----
long_b <- do.call(rbind, lapply(studies, function(s) {
  df <- tabs[[s]]; data.frame(cohort = cohort_titles[s], subtype = df$subtype, BAP1 = df$BAP1) }))
long_b <- long_b[!is.na(long_b$subtype), ]
long_b$subtype <- factor(long_b$subtype,
  levels = c("Epithelioid","Biphasic-E","Biphasic","Biphasic-S","Sarcomatoid"))
p_bap1_sub <- ggplot(long_b, aes(subtype, BAP1, fill = subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = .6, linewidth = .2) +
  geom_jitter(width = .15, size = .4, alpha = .4) +
  facet_wrap(~ cohort, scales = "free", nrow = 1) +
  scale_fill_manual(values = subtype_cols, guide = "none") +
  gtheme_no_rot + xlab("") + ylab("BAP1 expression (log2)") +
  ggtitle("Does BAP1 expression track histology?") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
pdf(file.path("Plots", "BAP1_by_subtype_bulk.pdf"), width = 9, height = 3.4)
print(p_bap1_sub); dev.off()

## P2a: RUNX~BAP1 within EPITHELIOID only ----
long_e <- do.call(rbind, lapply(studies, function(s) {
  df <- tabs[[s]]; df <- df[df$subtype == "Epithelioid", ]
  do.call(rbind, lapply(runx_genes, function(g) data.frame(
    cohort = cohort_titles[s], gene = g, BAP1 = df$BAP1, TF = df[[g]]))) }))
long_e$cohort <- factor(long_e$cohort, levels = cohort_titles[studies])
long_e$gene   <- factor(long_e$gene, levels = runx_genes)
long_e <- long_e[is.finite(long_e$BAP1) & is.finite(long_e$TF), ]
p_epi <- ggplot(long_e, aes(BAP1, TF)) +
  geom_point(size = .7, alpha = .6, color = "#2471a3") +
  geom_smooth(method = "lm", se = FALSE, color = "grey25", linewidth = .4) +
  stat_cor(method = "spearman", size = 2.3) +
  facet_grid(gene ~ cohort, scales = "free") +
  gtheme_no_rot + xlab("BAP1 expression (log2)") + ylab("RUNX expression (log2)") +
  ggtitle("RUNX vs BAP1 within EPITHELIOID tumors only")
pdf(file.path("Plots", "RUNX_vs_BAP1_epithelioid_only.pdf"), width = 8.5, height = 7)
print(p_epi); dev.off()

## P2b: heatmap overall vs epithelioid vs partial(|sarc) rho ----
hm <- rbind(
  data.frame(gene = p2$gene, cohort = p2$cohort, kind = "overall",        rho = p2$rho_overall,       p = p2$p_overall),
  data.frame(gene = p2$gene, cohort = p2$cohort, kind = "epithelioid",    rho = p2$rho_epithelioid,   p = p2$p_epithelioid),
  data.frame(gene = p2$gene, cohort = p2$cohort, kind = "partial|sarc",   rho = p2$partial_rho_sarc,  p = p2$partial_p_sarc))
hm$sig <- ifelse(is.na(hm$p), "", ifelse(hm$p < 0.001, "***",
           ifelse(hm$p < 0.01, "**", ifelse(hm$p < 0.05, "*", ""))))
hm$cohort <- factor(cohort_titles[hm$cohort], levels = cohort_titles[studies])
hm$gene   <- factor(hm$gene, levels = rev(tf_genes))
hm$kind   <- factor(hm$kind, levels = c("overall","epithelioid","partial|sarc"))
p_hm <- ggplot(hm, aes(cohort, gene, fill = rho)) +
  geom_tile(color = "white", linewidth = .5) +
  geom_text(aes(label = ifelse(is.na(rho), "NA", sprintf("%.2f%s", rho, sig))), size = 2.6) +
  facet_wrap(~ kind, nrow = 1) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-1,1), name = "rho\n(RUNX~BAP1)") +
  gtheme_no_rot + xlab("") + ylab("") +
  ggtitle("RUNX~BAP1 correlation: overall vs epithelioid-only vs partial (adjusting sarc score)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
pdf(file.path("Plots", "RUNX_BAP1_histology_adjusted_heatmap.pdf"), width = 9, height = 3.4)
print(p_hm); dev.off()

## P3a: scATAC - BAP1 gene score vs sarc score, and sarc by BAP1 status ----
sa_c$BAP1_status <- factor(sa_c$BAP1_status, levels = c("retained","lost"))
sp1 <- ggplot(sa_c, aes(BAP1_gs, sarc_score)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey55", linewidth = .4, linetype = "dashed") +
  geom_point(aes(color = BAP1_status), size = 2.4) +
  ggrepel::geom_text_repel(aes(label = sample), size = 2.4, max.overlaps = 100) +
  stat_cor(method = "spearman", size = 2.6) +
  scale_color_manual(values = c(retained = "#2471a3", lost = "#c0392b"),
                     labels = c(retained = "BAP1-retained", lost = "BAP1-lost"), name = NULL) +
  gtheme_no_rot + xlab("BAP1 gene score (scATAC)") + ylab("sarcomatoid score (cNMF20)") +
  ggtitle("scATAC: BAP1 vs histology axis")
pdf(file.path("Plots", "scATAC_BAP1_vs_sarcscore.pdf"), width = 4.8, height = 3.8)
print(sp1); dev.off()

cat("\nDONE\n")
