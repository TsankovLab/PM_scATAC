###############################################################################
# R2_Q14 -- do the paper's IRF8 and YY1 TARGET genes (supplement tabs 2C, 2D) move
#           with BAP1 status in TCGA?
#
# Inferred TF activity is, at heart, the collective behaviour of a TF's targets.  A
# plain approximation -- mean z-score of the published target genes per tumour -- is
# tested here with both labels:
#   paper  : tab 2A BAP1_status (inactivated vs no_inactivation)
#   ours   : genetic class, ALTERED vs WILD-TYPE
# and with/without the sarcomatoid score.  This is NOT the Osmanbeyoglu model (no
# regression on the TF-target prior, no proteomics); it asks only whether the
# target-level signal is visible at all.
#
# IRF8 targets are expected to be dominated by immune genes, so the score is also
# related to PTPRC (CD45) to see how much of any effect is infiltration.
#
# Output: paper_target_scores.csv, Plots/paper_IRF8_YY1_target_scores.pdf
###############################################################################
suppressPackageStartupMessages({ library(readxl); library(ggplot2) })
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R2_Q14"))
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
X <- "21598290cd180804-sup-205173_2_supp_5073240_pg14sl.xlsx"

tg <- list(YY1  = as.data.frame(read_excel(X, sheet = "2C_YY1 target genes")),
           IRF8 = as.data.frame(read_excel(X, sheet = "2D_IRF8 target genes")))
tg <- lapply(tg, function(d) unique(as.character(d[[grep("gene", names(d), ignore.case = TRUE)[1]]])))

E    <- as.matrix(readRDS(file.path(ROOT, "bulkRNA_meso", "bulk_RNA_studies.rds"))$tcga)
meta <- readRDS(file.path(ROOT, "bulkRNA_meso", "bulk_RNA_studies_metadata.rds"))$tcga[colnames(E), ]
E    <- E[apply(E, 1, function(x) sd(x) > 0 && mean(x > 0) > 0.5), ]
Z    <- t(scale(t(E)))
for (k in names(tg)) cat(sprintf("%s targets: %d listed, %d measured in TCGA\n",
                                 k, length(tg[[k]]), sum(tg[[k]] %in% rownames(Z))))
cat(sprintf("overlap of IRF8 and YY1 target lists: %d genes\n", length(intersect(tg$IRF8, tg$YY1))))

S <- data.frame(sample = colnames(E), tcga_id = gsub("\\.", "-", colnames(E)),
                IRF8_target_score = colMeans(Z[intersect(tg$IRF8, rownames(Z)), ]),
                YY1_target_score  = colMeans(Z[intersect(tg$YY1,  rownames(Z)), ]),
                IRF8_mRNA = E["IRF8", ], YY1_mRNA = E["YY1", ], PTPRC = E["PTPRC", ],
                sarc = as.numeric(meta$sarc_score), stringsAsFactors = FALSE)

A <- as.data.frame(read_excel(X, sheet = "2A_BAP1 status"), check.names = FALSE)
S$paper_status <- A$BAP1_status[match(S$tcga_id, A[[1]])]
G <- read.csv("tcga_genetic/TCGA_BAP1_genetic_status.csv", stringsAsFactors = FALSE)
S$our_class <- G$class[match(S$sample, G$sample)]
cat("\npaper status categories:\n"); print(table(S$paper_status, useNA = "ifany"))
S$paper_inact <- ifelse(is.na(S$paper_status), NA, S$paper_status != "no_inactivation")
S$our_altered <- ifelse(S$our_class == "ALTERED", TRUE, ifelse(S$our_class == "WILD-TYPE", FALSE, NA))
write.csv(S, "paper_target_scores.csv", row.names = FALSE)

cat("\n=== target scores and mRNA by BAP1 label (positive diff = higher when BAP1 inactivated/altered) ===\n")
res <- do.call(rbind, lapply(c("paper_inact", "our_altered"), function(lab)
  do.call(rbind, lapply(c("IRF8_target_score","YY1_target_score","IRF8_mRNA","YY1_mRNA","PTPRC"), function(v){
    d <- S[!is.na(S[[lab]]), ]
    g <- d[[lab]]
    w  <- suppressWarnings(wilcox.test(d[[v]][g], d[[v]][!g]))
    f1 <- summary(lm(d[[v]] ~ g))$coefficients
    f2 <- summary(lm(d[[v]] ~ g + d$sarc))$coefficients
    data.frame(label = lab, n_inact = sum(g), n_other = sum(!g), variable = v,
               diff = mean(d[[v]][g]) - mean(d[[v]][!g]), p_wilcox = w$p.value,
               p_lm = f1["gTRUE", 4], p_lm_sarc_adj = f2["gTRUE", 4], stringsAsFactors = FALSE)
  }))))
print(transform(res, diff = round(diff, 3), p_wilcox = signif(p_wilcox, 3), p_lm = signif(p_lm, 3),
                p_lm_sarc_adj = signif(p_lm_sarc_adj, 3)), row.names = FALSE)

cat("\n=== what drives the IRF8 target score? (Spearman, all tumours) ===\n")
for (v in c("PTPRC", "IRF8_mRNA", "sarc", "YY1_target_score"))
  cat(sprintf("IRF8 target score vs %-17s rho = %.3f\n", v,
              cor(S$IRF8_target_score, S[[v]], method = "spearman", use = "complete.obs")))
cat(sprintf("YY1 target score  vs %-17s rho = %.3f\n", "PTPRC",
            cor(S$YY1_target_score, S$PTPRC, method = "spearman")))
fit <- summary(lm(IRF8_target_score ~ paper_inact + PTPRC, data = S))$coefficients
cat(sprintf("IRF8 target score ~ paper status + PTPRC: status p = %.3g (PTPRC p = %.3g)\n",
            fit["paper_inactTRUE", 4], fit["PTPRC", 4]))

L <- rbind(
  data.frame(label = "paper status", grp = ifelse(S$paper_inact, "inactivated", "no inactivation"),
             IRF8 = S$IRF8_target_score, YY1 = S$YY1_target_score),
  data.frame(label = "our genetic class", grp = ifelse(S$our_altered, "ALTERED", "WILD-TYPE"),
             IRF8 = S$IRF8_target_score, YY1 = S$YY1_target_score))
L <- L[!is.na(L$grp), ]
L <- rbind(data.frame(L[, 1:2], set = "IRF8 targets", score = L$IRF8),
           data.frame(L[, 1:2], set = "YY1 targets",  score = L$YY1))
p <- ggplot(L, aes(grp, score)) +
  geom_boxplot(outlier.shape = NA, width = .6, linewidth = .3, fill = "grey92") +
  geom_jitter(width = .15, size = .8, alpha = .7) +
  facet_grid(set ~ label, scales = "free") +
  gtheme_no_rot + xlab(NULL) + ylab("mean z-score of target genes") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  ggtitle("TCGA MESO: paper IRF8 / YY1 target-gene scores by BAP1 label")
pdf("Plots/paper_IRF8_YY1_target_scores.pdf", width = 5.6, height = 5.2)
print(p); dev.off()
cat("\nDONE\n")
