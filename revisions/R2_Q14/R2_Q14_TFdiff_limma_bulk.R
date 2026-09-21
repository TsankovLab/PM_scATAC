###############################################################################
# R2_Q14 -- the same histology-adjusted BAP1 comparison, in bulk RNA.
#
#     expression ~ BAP1_status + sarc_score      (limma, per cohort)
#
# BAP1 status is CLINICAL.  Only MESOMICS carries one: IHC.BAP1, routine BAP1
# immunohistochemistry (NO = nuclear staining lost, YES = retained).  Bueno and
# TCGA have no BAP1 annotation of any kind in their metadata -- only arm-level
# chr3p calls (TCGA STATUS_3P, Bueno FISH.chrom3), which are molecular inference
# of exactly the sort excluded here -- so they cannot enter this test and are
# reported as unavailable rather than filled in from expression.
#
# Direction is verified from the data rather than assumed: IHC "NO" tumours have
# lower BAP1 mRNA (5.34 vs 6.07), so NO = lost.
#
# The same confounding is present as in the scATAC cohort and in the same
# direction -- BAP1-lost tumours are LESS sarcomatoid (sarc score 4.74 vs 5.32) --
# which is why the adjustment is carried over.
#
# Bulk is the better-powered test of the same question: ~106 tumours instead of 9,
# so where scATAC can only rank, this can actually reject.  It measures a different
# thing (TF mRNA rather than motif accessibility), so agreement is corroboration and
# disagreement is not necessarily contradiction.
#
# Input : bulkRNA_meso/bulk_RNA_studies.rds, bulk_RNA_studies_metadata.rds
#         BAP1_TF_limma_sarcadjusted.csv (the scATAC result, for the comparison)
# Output: BAP1_bulk_limma_sarcadjusted.csv, BAP1_bulk_vs_scATAC.csv,
#         Plots/BAP1_bulk_limma_volcano.pdf, Plots/BAP1_bulk_vs_scATAC.pdf
###############################################################################
suppressPackageStartupMessages({ library(limma); library(ggplot2); library(ggrepel) })

ROOT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUTDIR <- file.path(ROOT, "git_repo_claude", "R2_Q14")
BULK   <- file.path(ROOT, "bulkRNA_meso")
setwd(OUTDIR); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))

RUNX <- c("RUNX1", "RUNX2", "RUNX3")
IRF  <- paste0("IRF", 1:9)
NTOP <- 20

bl   <- readRDS(file.path(BULK, "bulk_RNA_studies.rds"))
meta <- readRDS(file.path(BULK, "bulk_RNA_studies_metadata.rds"))

for (s in c("bueno", "tcga"))
  cat(sprintf("%s: no clinical BAP1 annotation (%s) -- excluded\n", s,
              if (any(grepl("BAP1", colnames(meta[[s]]), ignore.case = TRUE)))
                "BAP1 columns present but not a status call" else "no BAP1 column"))

## ---- MESOMICS ------------------------------------------------------------------
E  <- as.matrix(bl$mesomics)
md <- meta$mesomics[match(colnames(E), meta$mesomics$Sample), , drop = FALSE]
ih <- as.character(md$IHC.BAP1)

## only the two unambiguous IHC calls; "YES in MMS/NO in MME" is discordant between
## blocks of the same tumour and cannot be assigned.
ok <- ih %in% c("NO", "YES") & is.finite(md$sarc_score)
cat(sprintf("\nMESOMICS: %d tumours with an unambiguous IHC call and a sarc score",
            sum(ok)))
cat(sprintf(" (dropped %d: %s)\n", sum(!ok),
            paste(names(table(ih[!ok], useNA = "ifany")), collapse = ", ")))
E  <- E[, ok, drop = FALSE]; md <- md[ok, , drop = FALSE]
BAP1 <- factor(ifelse(ih[ok] == "NO", "lost", "retained"),
               levels = c("retained", "lost"))
SARC <- as.numeric(md$sarc_score)
print(table(BAP1))
cat(sprintf("BAP1 mRNA  : %.2f lost vs %.2f retained  (direction check)\n",
            mean(E["BAP1", BAP1 == "lost"]), mean(E["BAP1", BAP1 == "retained"])))
cat(sprintf("sarc score : %.2f lost vs %.2f retained  (the confound)\n",
            mean(SARC[BAP1 == "lost"]), mean(SARC[BAP1 == "retained"])))

## drop genes that are flat or largely undetected
keep <- apply(E, 1, function(x) all(is.finite(x)) && sd(x) > 0 && mean(x > 0) > 0.5)
E <- E[keep, , drop = FALSE]
cat("genes tested:", nrow(E), "\n")

## ---- fits --------------------------------------------------------------------------
## trend = TRUE: bulk log-expression has a mean-variance trend, so the empirical Bayes
## prior is fitted as a function of expression level rather than as one constant.
fit_adj <- eBayes(lmFit(E, model.matrix(~ BAP1 + SARC)), trend = TRUE, robust = TRUE)
fit_raw <- eBayes(lmFit(E, model.matrix(~ BAP1)),        trend = TRUE, robust = TRUE)

adj  <- topTable(fit_adj, coef = "BAP1lost", number = Inf, sort.by = "none")
raw  <- topTable(fit_raw, coef = "BAP1lost", number = Inf, sort.by = "none")
hist <- topTable(fit_adj, coef = "SARC",     number = Inf, sort.by = "none")
i    <- match(rownames(adj), rownames(raw)); j <- match(rownames(adj), rownames(hist))

RES <- data.frame(gene = rownames(adj),
  logFC_adj = adj$logFC, t_adj = adj$t, P_adj = adj$P.Value, FDR_adj = adj$adj.P.Val,
  logFC_unadj = raw$logFC[i], t_unadj = raw$t[i], P_unadj = raw$P.Value[i],
  FDR_unadj = raw$adj.P.Val[i],
  t_sarc = hist$t[j], P_sarc = hist$P.Value[j], FDR_sarc = hist$adj.P.Val[j],
  stringsAsFactors = FALSE)
RES <- RES[order(RES$P_adj), ]
write.csv(RES, "BAP1_bulk_limma_sarcadjusted.csv", row.names = FALSE)

cat("\n=== genes passing, MESOMICS ===\n")
for (a in c(0.01, 0.05, 0.10)) cat(sprintf(
  "FDR < %.2f : BAP1 adjusted %5d | BAP1 unadjusted %5d | histology %5d\n",
  a, sum(RES$FDR_adj < a), sum(RES$FDR_unadj < a), sum(RES$FDR_sarc < a)))

cat(sprintf("\n=== top %d genes, histology-adjusted ===\n", NTOP))
print(head(RES[, c("gene","logFC_adj","t_adj","P_adj","FDR_adj","t_unadj","t_sarc")], NTOP),
      row.names = FALSE, digits = 3)

cat("\n=== RUNX ===\n")
print(RES[RES$gene %in% RUNX, c("gene","logFC_adj","t_adj","P_adj","FDR_adj",
                                "logFC_unadj","t_unadj","P_unadj","t_sarc","P_sarc")],
      row.names = FALSE, digits = 3)
cat("\n=== IRF ===\n")
ir <- RES[RES$gene %in% IRF, c("gene","logFC_adj","t_adj","P_adj","FDR_adj","t_unadj","t_sarc")]
print(ir[order(ir$gene), ], row.names = FALSE, digits = 3)

sets <- list(IRF = which(rownames(E) %in% IRF), RUNX = which(rownames(E) %in% RUNX))
sets <- sets[lengths(sets) > 1]
cam  <- cameraPR(fit_adj$t[, "BAP1lost"], sets)
cat("\n=== cameraPR on the adjusted BAP1 t ===\n"); print(cam, digits = 3)

## ---- does bulk agree with scATAC? ---------------------------------------------------
## Motif name to gene symbol is one-to-one for most cisBP entries, so the shared names
## can be compared directly.  Motif accessibility and mRNA are different measurements;
## this asks whether they point the same way, not whether they agree in magnitude.
if (file.exists("BAP1_TF_limma_sarcadjusted.csv")) {
  at <- read.csv("BAP1_TF_limma_sarcadjusted.csv", stringsAsFactors = FALSE)
  CMP <- merge(RES[, c("gene","t_adj","P_adj","FDR_adj")],
               at[, c("TF","t_adj","P_adj")], by.x = "gene", by.y = "TF",
               suffixes = c("_bulk", "_atac"))
  ct <- suppressWarnings(cor.test(CMP$t_adj_bulk, CMP$t_adj_atac, method = "spearman"))
  write.csv(CMP[order(CMP$P_adj_bulk), ], "BAP1_bulk_vs_scATAC.csv", row.names = FALSE)
  cat(sprintf("\n=== bulk mRNA vs scATAC motif activity, %d shared TF names ===\n", nrow(CMP)))
  cat(sprintf("Spearman of the adjusted BAP1 t: rho = %.3f, p = %.3g\n",
              ct$estimate, ct$p.value))
  cat(sprintf("agree in sign: %d of %d (%.0f%%)\n",
              sum(sign(CMP$t_adj_bulk) == sign(CMP$t_adj_atac)), nrow(CMP),
              100 * mean(sign(CMP$t_adj_bulk) == sign(CMP$t_adj_atac))))
  CMP$lab <- ifelse(CMP$gene %in% RUNX |
                    rank(CMP$P_adj_bulk, ties.method = "min") <= 12, CMP$gene, NA)
  pc <- ggplot(CMP, aes(t_adj_atac, t_adj_bulk)) +
    geom_hline(yintercept = 0, linewidth = .2, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
    geom_point(colour = "grey70", size = 1, alpha = .8) +
    geom_point(data = subset(CMP, !is.na(lab)), colour = "grey20", size = 1.4) +
    geom_text_repel(aes(label = lab), size = 2.2, colour = "grey20",
                    max.overlaps = 40, seed = 1, min.segment.length = 0) +
    gtheme_no_rot +
    xlab("t, scATAC motif activity (histology-adjusted)") +
    ylab("t, MESOMICS bulk mRNA (histology-adjusted)") +
    ggtitle(sprintf("BAP1 effect, two modalities (%d shared TFs, rho = %.2f)",
                    nrow(CMP), ct$estimate))
  pdf(file.path("Plots", "BAP1_bulk_vs_scATAC.pdf"), width = 5.4, height = 4.4)
  print(pc); dev.off()
}

## ---- volcano -------------------------------------------------------------------------
RES$rank <- rank(RES$P_adj, ties.method = "min")
RES$lab  <- ifelse(RES$rank <= NTOP, RES$gene, NA)
lim <- max(abs(RES$logFC_adj))
pv <- ggplot(RES, aes(logFC_adj, -log10(P_adj))) +
  geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = .3, colour = "grey55") +
  geom_point(colour = "grey78", size = .8, alpha = .6) +
  geom_point(data = subset(RES, FDR_adj < 0.05), colour = "grey45", size = .9) +
  geom_point(data = subset(RES, rank <= NTOP), colour = "grey15", size = 1.4) +
  geom_text_repel(aes(label = lab), size = 2.2, colour = "grey15", max.overlaps = 40,
                  seed = 1, min.segment.length = 0, segment.size = .2) +
  scale_x_continuous(limits = c(-lim, lim)) +
  gtheme_no_rot +
  xlab("BAP1 effect, lost - retained (log2, histology-adjusted)") +
  ylab(expression(-log[10]~"moderated"~italic(P))) +
  ggtitle(sprintf("MESOMICS bulk RNA ~ BAP1(IHC) + sarcomatoid score (n = %d, %d genes)",
                  ncol(E), nrow(E)))
pdf(file.path("Plots", "BAP1_bulk_limma_volcano.pdf"), width = 6.4, height = 4.8)
print(pv); dev.off()

cat("\nDONE\n")
