###############################################################################
# R2_Q14 -- TCGA MESO, genetic BAP1 ALTERED vs WILD-TYPE, TF GENES ONLY.
#
#     TF mRNA ~ BAP1_status + sarc_score       (limma, trend + robust eBayes)
#
# Same labels and model as R2_Q14_TCGA_genetic_BAP1_limma.R, but the test is run on
# the TF genes alone rather than genome-wide.  That changes two things, both on
# purpose:
#   * the empirical Bayes prior is estimated from TF genes only, and
#   * the BH correction is over the TF panel (~680 genes), not ~17,500,
# so it is a focused test of the TF question rather than TFs read out of a
# genome-wide list.
#
# TF genes = the TFs behind the 869 chromVAR motifs of our scATAC cohort, i.e. the
# panel used in every scATAC-bulk comparison in R2_Q14, so results line up with them.
#
# chr3p TFs are flagged: 10 of 26 ALTERED tumours are BAP1 deep deletions, which also
# remove neighbouring 3p genes.
#
# Input : tcga_genetic/TCGA_BAP1_genetic_status.csv, bulkRNA_meso/*.rds,
#         BAP1_TF_limma_sarcadjusted.csv (the TF panel)
# Output: TCGA_genetic_BAP1_TFonly.csv, Plots/TCGA_genetic_BAP1_TFonly_volcano.pdf
###############################################################################
suppressPackageStartupMessages({
  library(limma); library(ggplot2); library(ggrepel); library(ArchR)
})
addArchRThreads(1); addArchRGenome("hg38")

ROOT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUTDIR <- file.path(ROOT, "git_repo_claude", "R2_Q14")
BULK   <- file.path(ROOT, "bulkRNA_meso")
setwd(OUTDIR); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
CEN3 <- 91.1e6; NTOP <- 20

## ---- data & labels -------------------------------------------------------------
G    <- read.csv("tcga_genetic/TCGA_BAP1_genetic_status.csv", stringsAsFactors = FALSE)
E    <- as.matrix(readRDS(file.path(BULK, "bulk_RNA_studies.rds"))$tcga)
meta <- readRDS(file.path(BULK, "bulk_RNA_studies_metadata.rds"))$tcga[colnames(E), , drop = FALSE]
cls  <- G$class[match(colnames(E), G$sample)]
sarc <- as.numeric(meta$sarc_score)
keep <- cls %in% c("ALTERED", "WILD-TYPE") & is.finite(sarc)
E    <- E[, keep, drop = FALSE]
BAP1 <- factor(ifelse(cls[keep] == "ALTERED", "lost", "retained"), levels = c("retained","lost"))
SARC <- sarc[keep]

## ---- TF genes only --------------------------------------------------------------
TFS <- setdiff(unique(read.csv("BAP1_TF_limma_sarcadjusted.csv", stringsAsFactors = FALSE)$TF), "BAP1")
E   <- E[rownames(E) %in% TFS, , drop = FALSE]
E   <- E[apply(E, 1, function(x) all(is.finite(x)) && sd(x) > 0 && mean(x > 0) > 0.5), , drop = FALSE]
cat("tumours:", ncol(E), "| ALTERED", sum(BAP1 == "lost"), "| WILD-TYPE", sum(BAP1 == "retained"), "\n")
cat("TF panel:", length(TFS), "names |", nrow(E), "TF genes measured and expressed in TCGA\n")

gn  <- getGeneAnnotation()$genes
loc <- data.frame(gene = mcols(gn)$symbol, chr = as.character(seqnames(gn)),
                  start = start(gn), stringsAsFactors = FALSE)
loc <- loc[!duplicated(loc$gene), ]
i3  <- match(rownames(E), loc$gene)
chr3p <- !is.na(i3) & loc$chr[i3] == "chr3" & loc$start[i3] < CEN3

## ---- fits ------------------------------------------------------------------------
fa <- eBayes(lmFit(E, model.matrix(~ BAP1 + SARC)), trend = TRUE, robust = TRUE)
fu <- eBayes(lmFit(E, model.matrix(~ BAP1)),        trend = TRUE, robust = TRUE)
a  <- topTable(fa, coef = "BAP1lost", number = Inf, sort.by = "none")
u  <- topTable(fu, coef = "BAP1lost", number = Inf, sort.by = "none")
h  <- topTable(fa, coef = "SARC",     number = Inf, sort.by = "none")
R  <- data.frame(TF = rownames(a), chr3p = chr3p,
                 logFC = a$logFC, t = a$t, P = a$P.Value, FDR = a$adj.P.Val,
                 t_unadjusted = u$t[match(rownames(a), rownames(u))],
                 FDR_unadjusted = u$adj.P.Val[match(rownames(a), rownames(u))],
                 t_sarc = h$t[match(rownames(a), rownames(h))],
                 FDR_sarc = h$adj.P.Val[match(rownames(a), rownames(h))],
                 stringsAsFactors = FALSE)
R <- R[order(R$P), ]
write.csv(R, "TCGA_genetic_BAP1_TFonly.csv", row.names = FALSE)

cat(sprintf("\nempirical Bayes prior df (TF genes only): %.1f\n", median(fa$df.prior)))
cat("\n=== TF genes passing (FDR over the TF panel) ===\n")
for (x in c(0.01, 0.05, 0.10)) cat(sprintf(
  "FDR < %.2f : adjusted %3d (chr3p %d; up %d, down %d) | unadjusted %3d | histology %3d\n", x,
  sum(R$FDR < x), sum(R$FDR < x & R$chr3p), sum(R$FDR < x & R$t > 0), sum(R$FDR < x & R$t < 0),
  sum(R$FDR_unadjusted < x), sum(R$FDR_sarc < x)))
hit <- R$FDR < 0.05
if (sum(hit) && sum(R$chr3p)) {
  ft <- fisher.test(table(hit, R$chr3p))
  cat(sprintf("chr3p among FDR<0.05 TFs: %d of %d (%.1f%%) vs %.1f%% of the panel | OR %.1f, p = %.3g\n",
              sum(hit & R$chr3p), sum(hit), 100 * mean(R$chr3p[hit]), 100 * mean(R$chr3p),
              ft$estimate, ft$p.value))
}
cat(sprintf("\n=== all TFs at FDR < 0.05 (%d) ===\n", sum(hit)))
print(R[hit, c("TF","chr3p","logFC","t","P","FDR","t_unadjusted","t_sarc")],
      row.names = FALSE, digits = 3)

## ---- volcano -------------------------------------------------------------------------
R$rank <- rank(R$P, ties.method = "min")
R$lab  <- ifelse(R$rank <= NTOP, R$TF, NA)
lim <- max(abs(R$logFC))
pv <- ggplot(R, aes(logFC, -log10(P))) +
  geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = .3, colour = "grey55") +
  geom_point(colour = "grey80", size = 1, alpha = .7) +
  geom_point(data = subset(R, FDR < 0.05), colour = "grey40", size = 1.2) +
  geom_point(data = subset(R, FDR < 0.05 & chr3p), colour = "#b2182b", size = 1.3) +
  geom_point(data = subset(R, rank <= NTOP), shape = 21, fill = NA, colour = "grey5", size = 1.9) +
  geom_text_repel(aes(label = lab), size = 2.3, colour = "grey10", max.overlaps = 40,
                  seed = 1, min.segment.length = 0, segment.size = .2) +
  scale_x_continuous(limits = c(-lim, lim)) +
  gtheme_no_rot +
  xlab("BAP1 altered - wild-type (log2, adjusted for sarcomatoid score)") +
  ylab(expression(-log[10]~"moderated"~italic(P))) +
  ggtitle(sprintf("TCGA MESO, TF genes only (n = %d): %d altered vs %d wild-type",
                  nrow(R), sum(BAP1 == "lost"), sum(BAP1 == "retained")),
          subtitle = "dark = FDR < 0.05 over the TF panel; red = on chr3p")
pdf(file.path("Plots", "TCGA_genetic_BAP1_TFonly_volcano.pdf"), width = 6.4, height = 5)
print(pv); dev.off()

cat("\nDONE\n")
