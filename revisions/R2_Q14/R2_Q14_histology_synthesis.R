###############################################################
# R2_Q14 : Synthesis figure — BAP1 (expression / gene score) in relation to
#          RUNX activity is governed by HISTOLOGY and the SARCOMATOID SCORE.
#
# One coherent multi-panel figure that shows, in both modalities:
#   (i)   BAP1 tracks the histology axis (subtype / sarc score)
#   (ii)  RUNX tracks the same histology axis (opposite direction)
#   (iii) therefore the BAP1<->RUNX relationship is, to a large extent, a
#         sarcomatoid-score gradient (the confound made visible by colouring the
#         BAP1-vs-RUNX scatter by the continuous sarc score).
#
# Bulk (Bueno / TCGA / MESOMICS)  -> BAP1 & RUNX = mRNA expression
# scATAC (tumor_compartment)      -> BAP1 = gene score, RUNX = chromVAR activity
###############################################################

set.seed(1234)

GITROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo"
BULKDIR <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/bulkRNA_meso"
TCOMP   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment"
OUTDIR  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q14"
setwd(OUTDIR)
dir.create(file.path(OUTDIR, "Plots"), showWarnings = FALSE)

suppressMessages({ library(ggplot2); library(ggpubr); library(patchwork); library(ggrepel) })
source(file.path(GITROOT, "utils", "ggplot_aestetics.R"))

studies      <- c("bueno", "tcga", "mesomics")
cohort_title <- c(bueno = "Bueno", tcga = "TCGA", mesomics = "MESOMICS")
runx_genes   <- c("RUNX1", "RUNX2", "RUNX3")
tf_genes     <- c(runx_genes, "CBFB")
RUNX_MAIN    <- "RUNX2"          # representative RUNX for the bulk scatter panels
subtype_cols <- c(Epithelioid = "#2471a3", `Biphasic-E` = "#5dade2",
                  Biphasic = "#f5b041", `Biphasic-S` = "#e59866",
                  Sarcomatoid = "#c0392b")
status_cols  <- c(retained = "#2471a3", lost = "#c0392b")
SARC_GRAD    <- function(name) scale_color_viridis_c(option = "viridis", name = name)  # sequential = one hue

bl   <- readRDS(file.path(BULKDIR, "bulk_RNA_studies.rds"))
meta <- readRDS(file.path(BULKDIR, "bulk_RNA_studies_metadata.rds"))

# ---- tidy per-sample bulk tables (subtype + sarc_score), same as disentangle ----
get_tab <- function(s) {
  m  <- bl[[s]]; md <- meta[[s]]
  if (all(colnames(m) %in% rownames(md))) md <- md[colnames(m), , drop = FALSE] else
    if ("Sample" %in% colnames(md)) md <- md[match(colnames(m), md$Sample), , drop = FALSE]
  df <- data.frame(sample = colnames(m), cohort = s,
                   subtype = as.character(md$subtype),
                   BAP1 = as.numeric(m["BAP1", ]), stringsAsFactors = FALSE)
  for (g in tf_genes) df[[g]] <- if (g %in% rownames(m)) as.numeric(m[g, ]) else NA_real_
  df$sarc_score <- if ("sarc_score" %in% colnames(md)) as.numeric(md$sarc_score) else NA_real_
  df
}
tabs <- setNames(lapply(studies, get_tab), studies)
if (any(sapply(tabs, function(d) all(is.na(d$sarc_score))))) {
  nmf <- read.csv(file.path(TCOMP, "scrna", "cnmf_genelist_25_nfeat_5000.csv"),
                  stringsAsFactors = FALSE)
  sarc_genes <- head(nmf[["cNMF20"]][!is.na(nmf[["cNMF20"]])], 20)
  for (s in studies) if (all(is.na(tabs[[s]]$sarc_score))) {
    m <- bl[[s]]; g <- intersect(sarc_genes, rownames(m))
    tabs[[s]]$sarc_score <- colMeans(m[g, , drop = FALSE], na.rm = TRUE)
  }
}
bulk <- do.call(rbind, lapply(studies, function(s) {
  d <- tabs[[s]]; d$cohort <- cohort_title[s]; d }))
bulk$cohort  <- factor(bulk$cohort, levels = cohort_title[studies])
bulk$subtype <- factor(bulk$subtype,
  levels = c("Epithelioid","Biphasic-E","Biphasic","Biphasic-S","Sarcomatoid"))
bulk <- bulk[is.finite(bulk$sarc_score), ]

# ---- scATAC per-sample: BAP1 gene score, RUNX chromVAR activity, sarc score ----
act <- read.csv(file.path(OUTDIR, "TFactivity_deviation_per_sample.csv"),
                row.names = 1, check.names = FALSE)
gs  <- read.csv(file.path(OUTDIR, "BAP1_genescore_per_sample.csv"), stringsAsFactors = FALSE)
sc  <- read.csv(file.path(TCOMP, "scrna", "cnmf20_sarcomatoid_sample_order.csv"),
                stringsAsFactors = FALSE)
sc  <- setNames(sc$x, sc$sampleID)
sa  <- data.frame(sample = gs$Sample, BAP1_status = gs$BAP1_status,
                  BAP1_gs = gs$BAP1_genescore, stringsAsFactors = FALSE)
sa$sarc_score <- sc[sa$sample]
for (g in runx_genes) sa[[g]] <- if (g %in% rownames(act)) as.numeric(act[g, sa$sample]) else NA_real_
sa$RUNX_activity <- rowMeans(sa[, runx_genes], na.rm = TRUE)   # mean RUNX1/2/3 chromVAR deviation
sa <- sa[is.finite(sa$sarc_score), ]                          # 9 samples (P10,P23 lack a sarc score)
sa$BAP1_status <- factor(sa$BAP1_status, levels = c("retained","lost"))

# ================================================================
# PANELS
# ================================================================
th_lab <- theme(plot.title = element_text(size = 9, face = "bold"),
                plot.subtitle = element_text(size = 7.5),
                legend.key.size = unit(3.2, "mm"),
                legend.text = element_text(size = 7),
                legend.title = element_text(size = 7.5))

## ---- Row 1 : BULK (mRNA) --------------------------------------
# A. BAP1 expression rises along the histology axis (sarc score), by subtype
pA <- ggplot(bulk, aes(sarc_score, BAP1)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey55", linewidth = .35, linetype = "dashed") +
  geom_point(aes(color = subtype), size = .8, alpha = .75) +
  stat_cor(method = "spearman", size = 2.2, label.y.npc = "top") +
  facet_wrap(~ cohort, scales = "free", nrow = 1) +
  scale_color_manual(values = subtype_cols, name = "subtype", drop = TRUE, na.translate = FALSE) +
  gtheme_no_rot + th_lab +
  xlab("sarcomatoid score (cNMF20)") + ylab("BAP1 expression (log2)") +
  ggtitle("A  BAP1 tracks histology (bulk)",
          "BAP1 mRNA increases with sarcomatoid score / is lowest in epithelioid")

# B. RUNX expression rises along the SAME axis (opposite biology to BAP1)
pB <- ggplot(bulk, aes(sarc_score, .data[[RUNX_MAIN]])) +
  geom_smooth(method = "lm", se = FALSE, color = "grey55", linewidth = .35, linetype = "dashed") +
  geom_point(aes(color = subtype), size = .8, alpha = .75) +
  stat_cor(method = "spearman", size = 2.2, label.y.npc = "top") +
  facet_wrap(~ cohort, scales = "free", nrow = 1) +
  scale_color_manual(values = subtype_cols, name = "subtype", drop = TRUE, na.translate = FALSE) +
  gtheme_no_rot + th_lab +
  xlab("sarcomatoid score (cNMF20)") + ylab(paste0(RUNX_MAIN, " expression (log2)")) +
  ggtitle(paste0("B  ", RUNX_MAIN, " also tracks histology (bulk)"),
          "RUNX rises with sarcomatoid score — the same axis BAP1 follows")

# C. THE CONFOUND: BAP1 vs RUNX, coloured by continuous sarc score
pC <- ggplot(bulk, aes(BAP1, .data[[RUNX_MAIN]])) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30", linewidth = .35) +
  geom_point(aes(color = sarc_score), size = .9, alpha = .85) +
  stat_cor(method = "spearman", size = 2.2, label.y.npc = "top") +
  facet_wrap(~ cohort, scales = "free", nrow = 1) +
  SARC_GRAD("sarc\nscore") +
  gtheme_no_rot + th_lab +
  xlab("BAP1 expression (log2)") + ylab(paste0(RUNX_MAIN, " expression (log2)")) +
  ggtitle(paste0("C  BAP1 vs ", RUNX_MAIN, " is a sarcomatoid-score gradient (bulk)"),
          "the positive BAP1-RUNX correlation lines up with the sarc-score colour")

## ---- Row 2 : scATAC (gene score + chromVAR activity) ----------
lab_repel <- ggrepel::geom_text_repel(aes(label = sample), size = 2.3, max.overlaps = 100,
                                      seed = 1, min.segment.length = 0)
# D. BAP1 gene score tracks the histology axis
pD <- ggplot(sa, aes(sarc_score, BAP1_gs)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey55", linewidth = .35, linetype = "dashed") +
  geom_point(aes(color = BAP1_status), size = 2.6) + lab_repel +
  stat_cor(method = "spearman", size = 2.3) +
  scale_color_manual(values = status_cols, labels = c(retained = "BAP1-retained", lost = "BAP1-lost"),
                     name = NULL) +
  gtheme_no_rot + th_lab +
  xlab("sarcomatoid score (cNMF20)") + ylab("BAP1 gene score (scATAC)") +
  ggtitle("D  BAP1 gene score tracks histology (scATAC)",
          "BAP1-retained scATAC samples are the more sarcomatoid ones")

# E. RUNX chromVAR activity tracks the histology axis
pE <- ggplot(sa, aes(sarc_score, RUNX_activity)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey55", linewidth = .35, linetype = "dashed") +
  geom_point(aes(color = BAP1_status), size = 2.6) + lab_repel +
  stat_cor(method = "spearman", size = 2.3, label.x.npc = "right", label.y.npc = 0.12, hjust = 1) +
  scale_color_manual(values = status_cols, labels = c(retained = "BAP1-retained", lost = "BAP1-lost"),
                     name = NULL) +
  gtheme_no_rot + th_lab +
  xlab("sarcomatoid score (cNMF20)") + ylab("RUNX1/2/3 activity (chromVAR)") +
  ggtitle("E  RUNX activity tracks histology (scATAC)",
          "mean RUNX chromVAR deviation rises with sarcomatoid score")

# F. THE CONFOUND in scATAC: RUNX activity vs BAP1 gene score, coloured by sarc score
pF <- ggplot(sa, aes(BAP1_gs, RUNX_activity)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30", linewidth = .35) +
  geom_point(aes(color = sarc_score), size = 2.8) + lab_repel +
  stat_cor(method = "spearman", size = 2.3) +
  SARC_GRAD("sarc\nscore") +
  gtheme_no_rot + th_lab +
  xlab("BAP1 gene score (scATAC)") + ylab("RUNX1/2/3 activity (chromVAR)") +
  ggtitle("F  BAP1 vs RUNX activity follows the sarc score (scATAC)",
          "the BAP1-gs / RUNX-activity trend co-varies with the sarc-score colour")

# ================================================================
# ASSEMBLE
# ================================================================
fig <- (pA | pB | pC) / (pD | pE | pF) +
  plot_annotation(
    title = "R2_Q14 — BAP1 vs RUNX activity is governed by the epithelioid/sarcomatoid axis",
    subtitle = paste("Top: bulk RNA (Bueno/TCGA/MESOMICS), BAP1 & RUNX = expression.",
                     "Bottom: scATAC (n=9 tumors), BAP1 = gene score, RUNX = chromVAR activity."),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 8.5)))

pdf(file.path("Plots", "R2_Q14_histology_synthesis.pdf"), width = 15, height = 8.2)
print(fig); dev.off()

# ---- also drop the two "confound" panels as a compact stand-alone figure ----
fig2 <- (pC + ggtitle("Bulk RNA", paste0("BAP1 vs ", RUNX_MAIN, " expression, coloured by sarc score"))) /
        (pF + ggtitle("scATAC", "BAP1 gene score vs RUNX activity, coloured by sarc score"))
pdf(file.path("Plots", "R2_Q14_confound_BAP1_RUNX_by_sarcscore.pdf"), width = 11, height = 7.4)
print(fig2); dev.off()

# ---- numbers behind the figure ----
num <- list()
for (s in studies) {
  d <- tabs[[s]]
  num[[paste0(s,"_BAP1~sarc")]]      <- suppressWarnings(cor.test(d$BAP1, d$sarc_score, method="spearman"))
  num[[paste0(s,"_RUNX2~sarc")]]     <- suppressWarnings(cor.test(d[[RUNX_MAIN]], d$sarc_score, method="spearman"))
  num[[paste0(s,"_BAP1~RUNX2")]]     <- suppressWarnings(cor.test(d$BAP1, d[[RUNX_MAIN]], method="spearman"))
}
num[["scatac_BAP1gs~sarc"]]   <- suppressWarnings(cor.test(sa$BAP1_gs, sa$sarc_score, method="spearman"))
num[["scatac_RUNXact~sarc"]]  <- suppressWarnings(cor.test(sa$RUNX_activity, sa$sarc_score, method="spearman"))
num[["scatac_BAP1gs~RUNXact"]]<- suppressWarnings(cor.test(sa$BAP1_gs, sa$RUNX_activity, method="spearman"))
summ <- do.call(rbind, lapply(names(num), function(k)
  data.frame(comparison = k, rho = unname(num[[k]]$estimate), p = num[[k]]$p.value)))
write.csv(summ, file.path(OUTDIR, "histology_synthesis_correlations.csv"), row.names = FALSE)
cat("=== synthesis correlations ===\n"); print(summ, digits = 3)
cat("\nDONE\n")
