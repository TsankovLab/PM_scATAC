###############################################################################
# R2_Q14 -- TCGA MESO: BAP1 genetically ALTERED vs WILD-TYPE, adjusted for histology.
#
#     expression ~ BAP1_status + sarc_score        (limma, trend + robust eBayes)
#
# Same model as the scATAC limma (R2_Q14_TFdiff_limma.R): the cNMF20 sarcomatoid
# score is in the design, so the BAP1 coefficient is the lost-vs-retained difference
# with histology held constant.
#
# Label: genetic, from R2_Q14_TCGA_BAP1_genetic_status.R (TCGA PanCancer calls):
#   ALTERED   (26) truncating mutation, BAP1 fusion or deep deletion
#   WILD-TYPE (34) no BAP1 lesion detected
#   HEMIZYGOUS (22) and VUS (5) are excluded -- neither is clearly lost or retained.
#
# chr3p.  10 of the 26 ALTERED tumours carry a DEEP DELETION of BAP1, which removes
# neighbouring 3p21 genes too, so some differences will be positional (copy number)
# rather than BAP1-dependent regulation.  chr3p genes are flagged, their enrichment
# among the hits is tested, and TF congruence is reported with chr3p excluded.
#
# Input : tcga_genetic/TCGA_BAP1_genetic_status.csv, bulkRNA_meso/*.rds,
#         BAP1_TF_limma_sarcadjusted.csv (scATAC motif, 9 tumours),
#         BAP1_TFgenescore_limma.csv (scATAC gene score),
#         BAP1_TF_gsfiltered_results.csv (scATAC gene-score-filtered motif panel)
# Output: TCGA_genetic_BAP1_limma.csv, TCGA_genetic_BAP1_TFpanel.csv,
#         TCGA_genetic_vs_scATAC_congruence.csv,
#         Plots/TCGA_genetic_BAP1_volcano.pdf, Plots/TCGA_genetic_vs_scATAC.pdf
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
RUNX <- c("RUNX1","RUNX2","RUNX3"); IRF <- paste0("IRF", 1:9)
CEN3 <- 91.1e6; NTOP <- 20

## ---- data ----------------------------------------------------------------------
G    <- read.csv("tcga_genetic/TCGA_BAP1_genetic_status.csv", stringsAsFactors = FALSE)
E    <- as.matrix(readRDS(file.path(BULK, "bulk_RNA_studies.rds"))$tcga)
meta <- readRDS(file.path(BULK, "bulk_RNA_studies_metadata.rds"))$tcga
meta <- meta[colnames(E), , drop = FALSE]
cat("expression:", nrow(E), "genes x", ncol(E), "tumours; value range",
    paste(round(range(E, na.rm = TRUE), 2), collapse = " to "), "\n")

cls  <- G$class[match(colnames(E), G$sample)]
sarc <- as.numeric(meta$sarc_score)
keep <- cls %in% c("ALTERED", "WILD-TYPE") & is.finite(sarc)
E <- E[, keep, drop = FALSE]
BAP1 <- factor(ifelse(cls[keep] == "ALTERED", "lost", "retained"), levels = c("retained","lost"))
SARC <- sarc[keep]
cat("\ntumours:", ncol(E), "| ALTERED (lost)", sum(BAP1 == "lost"),
    "| WILD-TYPE (retained)", sum(BAP1 == "retained"), "\n")
cat(sprintf("confounding: sarc score %.2f lost vs %.2f retained (Wilcoxon p = %.3g)\n",
            mean(SARC[BAP1 == "lost"]), mean(SARC[BAP1 == "retained"]),
            suppressWarnings(wilcox.test(SARC ~ BAP1)$p.value)))
cat("histology:\n"); print(table(BAP1, subtype = meta$subtype[keep], useNA = "ifany"))
cat(sprintf("BAP1 mRNA: %.2f lost vs %.2f retained (label check)\n",
            mean(E["BAP1", BAP1 == "lost"]), mean(E["BAP1", BAP1 == "retained"])))

E <- E[apply(E, 1, function(x) all(is.finite(x)) && sd(x) > 0 && mean(x > 0) > 0.5), , drop = FALSE]
cat("genes tested:", nrow(E), "\n")

## ---- chr3p annotation ------------------------------------------------------------
gn  <- getGeneAnnotation()$genes
loc <- data.frame(gene = mcols(gn)$symbol, chr = as.character(seqnames(gn)),
                  start = start(gn), stringsAsFactors = FALSE)
loc <- loc[!duplicated(loc$gene), ]
i3  <- match(rownames(E), loc$gene)
chr3p <- !is.na(i3) & loc$chr[i3] == "chr3" & loc$start[i3] < CEN3

## ---- fits ------------------------------------------------------------------------
fit_adj <- eBayes(lmFit(E, model.matrix(~ BAP1 + SARC)), trend = TRUE, robust = TRUE)
fit_raw <- eBayes(lmFit(E, model.matrix(~ BAP1)),        trend = TRUE, robust = TRUE)
a <- topTable(fit_adj, coef = "BAP1lost", number = Inf, sort.by = "none")
u <- topTable(fit_raw, coef = "BAP1lost", number = Inf, sort.by = "none")
h <- topTable(fit_adj, coef = "SARC",     number = Inf, sort.by = "none")
R <- data.frame(gene = rownames(a), chr3p = chr3p,
                logFC_adj = a$logFC, t_adj = a$t, P_adj = a$P.Value, FDR_adj = a$adj.P.Val,
                t_unadj = u$t[match(rownames(a), rownames(u))],
                FDR_unadj = u$adj.P.Val[match(rownames(a), rownames(u))],
                t_sarc = h$t[match(rownames(a), rownames(h))],
                FDR_sarc = h$adj.P.Val[match(rownames(a), rownames(h))],
                stringsAsFactors = FALSE)
R <- R[order(R$P_adj), ]
write.csv(R, "TCGA_genetic_BAP1_limma.csv", row.names = FALSE)

cat("\n=== genome-wide, ALTERED vs WILD-TYPE ===\n")
for (aa in c(0.01, 0.05, 0.10)) cat(sprintf(
  "FDR < %.2f : adjusted %4d (chr3p %3d) | unadjusted %4d | histology %4d\n", aa,
  sum(R$FDR_adj < aa), sum(R$FDR_adj < aa & R$chr3p),
  sum(R$FDR_unadj < aa), sum(R$FDR_sarc < aa)))
hit <- R$FDR_adj < 0.05
if (sum(hit) > 0) {
  ft <- fisher.test(table(hit, R$chr3p))
  cat(sprintf("chr3p among adjusted FDR<0.05 hits: %d of %d (%.1f%%) vs %.1f%% of all genes | OR %.1f, p = %.3g\n",
              sum(hit & R$chr3p), sum(hit), 100 * mean(R$chr3p[hit]), 100 * mean(R$chr3p),
              ft$estimate, ft$p.value))
  cat(sprintf("of the chr3p hits, down in BAP1-altered: %d of %d\n",
              sum(hit & R$chr3p & R$t_adj < 0), sum(hit & R$chr3p)))
}
cat(sprintf("\n=== top %d, histology-adjusted ===\n", NTOP))
print(head(R[, c("gene","chr3p","logFC_adj","t_adj","P_adj","FDR_adj","t_unadj","t_sarc")], NTOP),
      row.names = FALSE, digits = 3)
cat(sprintf("\n=== top %d excluding chr3p ===\n", NTOP))
print(head(R[!R$chr3p, c("gene","logFC_adj","t_adj","P_adj","FDR_adj","t_sarc")], NTOP),
      row.names = FALSE, digits = 3)

## ---- TF panel (the 869 motif names from our scATAC cohort) --------------------------
at  <- read.csv("BAP1_TF_limma_sarcadjusted.csv", stringsAsFactors = FALSE)
TFS <- setdiff(unique(at$TF), "BAP1")
Et  <- E[rownames(E) %in% TFS, , drop = FALSE]
ft2 <- eBayes(lmFit(Et, model.matrix(~ BAP1 + SARC)), trend = TRUE, robust = TRUE)
tt  <- topTable(ft2, coef = "BAP1lost", number = Inf, sort.by = "none")
TP  <- data.frame(TF = rownames(tt), logFC = tt$logFC, t = tt$t, P = tt$P.Value,
                  FDR_panel = tt$adj.P.Val, chr3p = R$chr3p[match(rownames(tt), R$gene)],
                  stringsAsFactors = FALSE)
TP <- TP[order(TP$P), ]
write.csv(TP, "TCGA_genetic_BAP1_TFpanel.csv", row.names = FALSE)
cat(sprintf("\n=== TF panel: %d of our TFs, FDR over the panel ===\n", nrow(TP)))
cat(sprintf("FDR < 0.05: %d | < 0.10: %d | nominal P < 0.05: %d\n",
            sum(TP$FDR_panel < .05), sum(TP$FDR_panel < .10), sum(TP$P < .05)))
print(head(TP, 15), row.names = FALSE, digits = 3)
cat("\nRUNX / IRF:\n")
print(TP[TP$TF %in% c(RUNX, IRF), ], row.names = FALSE, digits = 3)
print(cameraPR(ft2$t[, "BAP1lost"], list(IRF = which(rownames(Et) %in% IRF),
                                         RUNX = which(rownames(Et) %in% RUNX))), digits = 3)

## ---- congruence with scATAC --------------------------------------------------------
gsl <- read.csv("BAP1_TFgenescore_limma.csv", stringsAsFactors = FALSE)
gsf <- read.csv("BAP1_TF_gsfiltered_results.csv", stringsAsFactors = FALSE)
tc  <- TP[!TP$chr3p %in% TRUE, c("TF","t","P","FDR_panel")]
sets <- list(
  "scATAC motif activity (all TFs)"          = setNames(at[, c("TF","t_adj")],      c("TF","t_atac")),
  "scATAC gene score"                         = setNames(gsl[, c("TF","t_gs")],      c("TF","t_atac")),
  "scATAC motif, gene-score-filtered panel"  = setNames(gsf[, c("TF","t_motif")],   c("TF","t_atac")))
C <- do.call(rbind, lapply(names(sets), function(k){
  m <- merge(tc, sets[[k]], by = "TF")
  ct <- suppressWarnings(cor.test(m$t, m$t_atac, method = "spearman"))
  data.frame(scATAC_readout = k, n_TF = nrow(m), rho = unname(ct$estimate), p = ct$p.value,
             pct_same_sign = 100 * mean(sign(m$t) == sign(m$t_atac)), stringsAsFactors = FALSE)
}))
write.csv(C, "TCGA_genetic_vs_scATAC_congruence.csv", row.names = FALSE)
cat("\n=== congruence of the BAP1 effect: TCGA genetic label vs scATAC (chr3p excluded) ===\n")
print(transform(C, rho = round(rho, 3), p = signif(p, 3), pct_same_sign = round(pct_same_sign)),
      row.names = FALSE)

top <- head(gsf[order(gsf$P_motif), c("TF","t_motif","P_motif")], NTOP)
top <- merge(top, TP[, c("TF","t","FDR_panel")], by = "TF", all.x = TRUE)
top$tcga <- ifelse(is.na(top$t), "not measured",
             ifelse(sign(top$t) == sign(top$t_motif) & top$FDR_panel < 0.05, "SIG same direction",
              ifelse(sign(top$t) != sign(top$t_motif) & top$FDR_panel < 0.05, "SIG opposite",
               ifelse(sign(top$t) == sign(top$t_motif), "same sign, ns", "opposite sign, ns"))))
top <- top[order(top$P_motif), ]
cat(sprintf("\n=== top %d scATAC gene-score-filtered TFs, in TCGA (genetic label) ===\n", NTOP))
print(top, row.names = FALSE, digits = 3)
cat("\n"); print(table(top$tcga))

## ---- volcano -------------------------------------------------------------------------
R$rank <- rank(R$P_adj, ties.method = "min"); R$lab <- ifelse(R$rank <= NTOP, R$gene, NA)
lim <- max(abs(R$logFC_adj))
pv <- ggplot(R, aes(logFC_adj, -log10(P_adj))) +
  geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = .3, colour = "grey55") +
  geom_point(colour = "grey80", size = .7, alpha = .6) +
  geom_point(data = subset(R, FDR_adj < 0.05), colour = "grey45", size = .9) +
  geom_point(data = subset(R, chr3p & FDR_adj < 0.05), colour = "#b2182b", size = 1) +
  geom_point(data = subset(R, rank <= NTOP), colour = "grey10", size = 1.4) +
  geom_text_repel(aes(label = lab), size = 2.2, colour = "grey10", max.overlaps = 40,
                  seed = 1, min.segment.length = 0, segment.size = .2) +
  scale_x_continuous(limits = c(-lim, lim)) +
  gtheme_no_rot +
  xlab("BAP1 altered - wild-type (log2, adjusted for sarcomatoid score)") +
  ylab(expression(-log[10]~"moderated"~italic(P))) +
  ggtitle(sprintf("TCGA MESO, genetic BAP1: %d altered vs %d wild-type (red = chr3p, FDR<0.05)",
                  sum(BAP1 == "lost"), sum(BAP1 == "retained")))
pdf(file.path("Plots", "TCGA_genetic_BAP1_volcano.pdf"), width = 6.6, height = 4.9)
print(pv); dev.off()

pl <- lapply(names(sets), function(k){
  m <- merge(tc, sets[[k]], by = "TF"); r <- cor(m$t, m$t_atac, method = "spearman")
  m$lab <- ifelse(m$TF %in% top$TF, m$TF, NA)
  ggplot(m, aes(t_atac, t)) +
    geom_hline(yintercept = 0, linewidth = .2, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
    geom_point(colour = "grey60", size = 1, alpha = .8) +
    geom_text_repel(aes(label = lab), size = 2, colour = "grey20", max.overlaps = 30,
                    seed = 1, min.segment.length = 0) +
    gtheme_no_rot + xlab(paste("t,", k)) + ylab("t, TCGA mRNA (genetic BAP1)") +
    ggtitle(sprintf("n = %d TFs, rho = %.2f", nrow(m), r))
})
pdf(file.path("Plots", "TCGA_genetic_vs_scATAC.pdf"), width = 5.2, height = 4.3)
for (p in pl) print(p); dev.off()

cat("\nDONE\n")
