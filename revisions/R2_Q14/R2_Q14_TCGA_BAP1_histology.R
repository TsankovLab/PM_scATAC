###############################################################################
# R2_Q14 -- TCGA MESO: is BAP1 genetic inactivation enriched by histology, and does BAP1
#           expression (mRNA, protein) track histology?
#
# Histology: TCGA HISTOLOGICAL_DIAGNOSIS (57 epithelioid, 23 biphasic, 2 sarcomatoid, 5 NOS;
#   NOS excluded).  With only 2 sarcomatoid tumours the main contrast is epithelioid vs
#   non-epithelioid (biphasic + sarcomatoid); the continuous cNMF20 sarcomatoid score is used
#   as a second, sample-size-independent histology axis.
# BAP1 genetics, two definitions:
#   ours   tcga_genetic/TCGA_BAP1_genetic_status.csv (public PanCancer calls):
#          ALTERED 26 / HEMIZYGOUS 22 / VUS 5 / WILD-TYPE 34
#   paper  Hmeljak et al. 2018 tab 2A integrated call (74 tumours):
#          inactivated 36 / possibly 6 / no inactivation 32
# Expression: BAP1 mRNA (bulk_RNA_studies.rds, log2) and BAP1 protein (RPPA, 63 tumours).
# Question for expression: does any histology association survive once genetic status is
# in the model (i.e. is it just the genetic enrichment seen through expression)?
# Reference: MESOMICS clinical BAP1 IHC by histology.
#
# Output: TCGA_BAP1_histology_genetics.csv, TCGA_BAP1_histology_expression.csv,
#         Plots/TCGA_BAP1_histology.pdf
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R2_Q14")); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))

G    <- read.csv("tcga_genetic/TCGA_BAP1_genetic_status.csv", stringsAsFactors = FALSE)
meta <- readRDS(file.path(ROOT, "bulkRNA_meso", "bulk_RNA_studies_metadata.rds"))$tcga
PAP  <- read.csv("paper_BAP1_status_vs_ours.csv", stringsAsFactors = FALSE)
G$histology <- as.character(meta$subtype[match(G$sample, rownames(meta))])
G$sarc      <- as.numeric(meta$sarc_score[match(G$sample, rownames(meta))])
G$paper     <- PAP$BAP1_status[match(G$tcga_id, PAP$tcga_id)]
G$hist2     <- ifelse(G$histology == "Epithelioid", "Epithelioid",
                ifelse(G$histology %in% c("Biphasic", "Sarcomatoid"), "Non-epithelioid", NA))
H <- G[!is.na(G$hist2), ]
cat("tumours with histology:", nrow(H), "\n"); print(table(H$histology))
cat("\ngenetic class x histology:\n"); print(addmargins(table(H$class, factor(H$histology, c("Epithelioid","Biphasic","Sarcomatoid")))))
cat("\npaper status x histology:\n"); print(addmargins(table(H$paper, factor(H$histology, c("Epithelioid","Biphasic","Sarcomatoid")), useNA = "ifany")))

## ---- genetic enrichment ------------------------------------------------------------
fisher_row <- function(label, pos, neg) {
  d <- H[pos | neg, ]; alt <- pos[pos | neg]
  tb <- table(factor(alt, c(TRUE, FALSE)), factor(d$hist2, c("Epithelioid", "Non-epithelioid")))
  ft <- fisher.test(tb)
  w  <- suppressWarnings(wilcox.test(d$sarc[alt], d$sarc[!alt]))
  data.frame(definition = label, n_altered = sum(alt), n_reference = sum(!alt),
             pct_altered_epithelioid = 100 * tb["TRUE", "Epithelioid"] / sum(tb[, "Epithelioid"]),
             pct_altered_nonepithelioid = 100 * tb["TRUE", "Non-epithelioid"] / sum(tb[, "Non-epithelioid"]),
             odds_ratio_epithelioid = unname(ft$estimate), OR_low = ft$conf.int[1], OR_high = ft$conf.int[2],
             fisher_p = ft$p.value, sarc_score_altered = median(d$sarc[alt]), sarc_score_reference = median(d$sarc[!alt]),
             sarc_wilcox_p = w$p.value)
}
cl <- H$class; pp <- H$paper
GEN <- rbind(
  fisher_row("ours: ALTERED vs WILD-TYPE", cl == "ALTERED", cl == "WILD-TYPE"),
  fisher_row("ours: any lesion (ALTERED+HEMIZYGOUS+VUS) vs WILD-TYPE", cl != "WILD-TYPE", cl == "WILD-TYPE"),
  fisher_row("ours: ALTERED vs all other", cl == "ALTERED", cl != "ALTERED"),
  fisher_row("paper: inactivated vs no inactivation", pp %in% "inactivated", pp %in% "no_inactivation"),
  fisher_row("paper: inactivated+possibly vs no inactivation", pp %in% c("inactivated", "possibly_inactivated"), pp %in% "no_inactivation"),
  fisher_row("ours: truncating mutation vs WILD-TYPE", H$truncating %in% TRUE, cl == "WILD-TYPE"),
  fisher_row("ours: deep deletion vs WILD-TYPE", H$gistic %in% -2, cl == "WILD-TYPE"))
write.csv(GEN, "TCGA_BAP1_histology_genetics.csv", row.names = FALSE)
cat("\n=== BAP1 genetic inactivation vs histology (TCGA) ===\n")
print(transform(GEN, pct_altered_epithelioid = round(pct_altered_epithelioid), pct_altered_nonepithelioid = round(pct_altered_nonepithelioid),
                odds_ratio_epithelioid = round(odds_ratio_epithelioid, 2), OR_low = round(OR_low, 2), OR_high = round(OR_high, 2),
                fisher_p = signif(fisher_p, 2), sarc_score_altered = round(sarc_score_altered, 2),
                sarc_score_reference = round(sarc_score_reference, 2), sarc_wilcox_p = signif(sarc_wilcox_p, 2)), row.names = FALSE)

## ---- expression vs histology -----------------------------------------------------------
EXP <- list()
for (v in c("BAP1_mRNA", "BAP1_protein_rppa")) {
  d <- H[!is.na(H[[v]]), ]
  w  <- wilcox.test(d[[v]][d$hist2 == "Epithelioid"], d[[v]][d$hist2 == "Non-epithelioid"])
  cs <- suppressWarnings(cor.test(d[[v]], d$sarc, method = "spearman"))
  f0 <- summary(lm(d[[v]] ~ d$hist2))$coefficients
  cls <- factor(d$class, c("WILD-TYPE", "VUS", "HEMIZYGOUS", "ALTERED"))
  f1 <- summary(lm(d[[v]] ~ cls + d$hist2))$coefficients
  f2 <- summary(lm(d[[v]] ~ cls + d$sarc))$coefficients
  wt <- d[d$class == "WILD-TYPE", ]
  ww <- if (sum(wt$hist2 == "Non-epithelioid") >= 3) wilcox.test(wt[[v]][wt$hist2 == "Epithelioid"], wt[[v]][wt$hist2 == "Non-epithelioid"])$p.value else NA
  EXP[[v]] <- data.frame(readout = v, n = nrow(d),
    median_epithelioid = median(d[[v]][d$hist2 == "Epithelioid"]), median_nonepithelioid = median(d[[v]][d$hist2 == "Non-epithelioid"]),
    wilcox_p = w$p.value, rho_sarc = unname(cs$estimate), rho_sarc_p = cs$p.value,
    beta_nonepi_unadjusted = f0[2, 1], p_nonepi_unadjusted = f0[2, 4],
    beta_nonepi_adj_genetics = f1["d$hist2Non-epithelioid", 1], p_nonepi_adj_genetics = f1["d$hist2Non-epithelioid", 4],
    beta_sarc_adj_genetics = f2["d$sarc", 1], p_sarc_adj_genetics = f2["d$sarc", 4],
    n_wildtype = nrow(wt), wildtype_only_wilcox_p = ww)
}
EXP <- do.call(rbind, EXP); write.csv(EXP, "TCGA_BAP1_histology_expression.csv", row.names = FALSE)
cat("\n=== BAP1 expression vs histology (TCGA) ===\n")
print(t(transform(EXP, median_epithelioid = round(median_epithelioid, 3), median_nonepithelioid = round(median_nonepithelioid, 3),
                  wilcox_p = signif(wilcox_p, 2), rho_sarc = round(rho_sarc, 3), rho_sarc_p = signif(rho_sarc_p, 2),
                  beta_nonepi_unadjusted = round(beta_nonepi_unadjusted, 3), p_nonepi_unadjusted = signif(p_nonepi_unadjusted, 2),
                  beta_nonepi_adj_genetics = round(beta_nonepi_adj_genetics, 3), p_nonepi_adj_genetics = signif(p_nonepi_adj_genetics, 2),
                  beta_sarc_adj_genetics = round(beta_sarc_adj_genetics, 3), p_sarc_adj_genetics = signif(p_sarc_adj_genetics, 2),
                  wildtype_only_wilcox_p = signif(wildtype_only_wilcox_p, 2))))
cat("\nBAP1 mRNA median by histology x genetic class:\n")
print(round(tapply(H$BAP1_mRNA, list(H$class, H$histology), median), 2))
cat("n:\n"); print(table(H$class, H$histology))

## ---- reference: MESOMICS clinical IHC ----------------------------------------------------
mm <- readRDS(file.path(ROOT, "bulkRNA_meso", "bulk_RNA_studies_metadata.rds"))$mesomics
ih <- as.character(mm$IHC.BAP1); hs <- as.character(mm$subtype); ok <- ih %in% c("NO", "YES")
tb <- table(lost = ih[ok] == "NO", epithelioid = hs[ok] == "Epithelioid")
ftm <- fisher.test(tb[c("TRUE", "FALSE"), c("TRUE", "FALSE")])
cat(sprintf("\nreference MESOMICS IHC: BAP1 lost in %.0f%% epithelioid vs %.0f%% non-epithelioid | OR %.2f, Fisher p = %.2g\n",
            100 * tb["TRUE", "TRUE"] / sum(tb[, "TRUE"]), 100 * tb["TRUE", "FALSE"] / sum(tb[, "FALSE"]), ftm$estimate, ftm$p.value))
print(table(IHC_lost = ih[ok] == "NO", histology = hs[ok]))

## ---- figure --------------------------------------------------------------------------------
H$histology <- factor(H$histology, c("Epithelioid", "Biphasic", "Sarcomatoid"))
H$class <- factor(H$class, c("ALTERED", "HEMIZYGOUS", "VUS", "WILD-TYPE"))
cc <- c(ALTERED = "#b2182b", HEMIZYGOUS = "#ef8a62", VUS = "#fddbc7", `WILD-TYPE` = "#67a9cf")
p1 <- ggplot(H, aes(histology, fill = class)) + geom_bar(position = "fill", width = .7) +
  geom_text(stat = "count", aes(label = after_stat(count)), position = position_fill(vjust = .5), size = 2.6) +
  scale_fill_manual(values = cc, name = "BAP1 genetics") + scale_y_continuous(labels = scales::percent) +
  gtheme_no_rot + xlab(NULL) + ylab("fraction of tumours") +
  ggtitle(sprintf("ALTERED vs WT, epithelioid vs not: Fisher p = %.2g", GEN$fisher_p[1]))
p2 <- ggplot(H, aes(histology, BAP1_mRNA)) + geom_boxplot(outlier.shape = NA, width = .6, fill = "grey93") +
  geom_jitter(aes(colour = class), width = .15, size = 1.5) + scale_colour_manual(values = cc, guide = "none") +
  gtheme_no_rot + xlab(NULL) + ylab("BAP1 mRNA (log2)") + ggtitle(sprintf("mRNA: epi vs non-epi p = %.2g", EXP$wilcox_p[1]))
p3 <- ggplot(H[!is.na(H$BAP1_protein_rppa), ], aes(histology, BAP1_protein_rppa)) + geom_boxplot(outlier.shape = NA, width = .6, fill = "grey93") +
  geom_jitter(aes(colour = class), width = .15, size = 1.5) + scale_colour_manual(values = cc, guide = "none") +
  gtheme_no_rot + xlab(NULL) + ylab("BAP1 protein (RPPA)") + ggtitle(sprintf("protein: epi vs non-epi p = %.2g", EXP$wilcox_p[2]))
p4 <- ggplot(H, aes(sarc, BAP1_mRNA, colour = class)) + geom_point(size = 1.6) +
  scale_colour_manual(values = cc, guide = "none") + gtheme_no_rot +
  xlab("sarcomatoid score (cNMF20)") + ylab("BAP1 mRNA (log2)") + ggtitle(sprintf("rho = %.2f", EXP$rho_sarc[1]))
ggsave("Plots/TCGA_BAP1_histology.pdf", (p1 | p2) / (p3 | p4), width = 10, height = 8)
cat("\nDONE\n")
