###############################################################################
# R2_Q14 -- TCGA MESO, genetic BAP1, TF genes only, WITHOUT histology correction,
#           plus a targeted check of IRF8 and EGR2.
#
#     TF mRNA ~ BAP1_status                  (limma, trend + robust eBayes)
#
# IRF8 and EGR2 were reported as higher in BAP1-altered tumours in the TCGA MESO
# study.  If they do not reproduce here, the likeliest reasons are differences in
# how "BAP1 altered" was defined and in the gene universe tested, so both are varied:
#   contrast A  ALTERED (26) vs WILD-TYPE (34)                       <- used so far
#   contrast B  ANY lesion (ALTERED + HEMIZYGOUS + VUS, 53) vs WILD-TYPE (34)
#   contrast C  ALTERED (26) vs ALL OTHER tumours (61)
# each fitted (i) TF genes only and (ii) genome-wide, unadjusted and adjusted for the
# sarcomatoid score.
#
# IRF8 is a myeloid / lymphoid TF, so in bulk tumour RNA its level largely reflects
# immune infiltration rather than expression in malignant cells.  PTPRC (CD45) is
# printed alongside as the infiltration reference.
#
# Input : tcga_genetic/TCGA_BAP1_genetic_status.csv, bulkRNA_meso/*.rds,
#         BAP1_TF_limma_sarcadjusted.csv (TF panel)
# Output: TCGA_genetic_BAP1_TFonly_unadjusted.csv, TCGA_BAP1_IRF8_EGR2_sensitivity.csv,
#         Plots/TCGA_genetic_BAP1_TFonly_unadjusted_volcano.pdf,
#         Plots/TCGA_BAP1_IRF8_EGR2_by_class.pdf
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
CEN3 <- 91.1e6; NTOP <- 20; GOI <- c("IRF8", "EGR2")

G    <- read.csv("tcga_genetic/TCGA_BAP1_genetic_status.csv", stringsAsFactors = FALSE)
E0   <- as.matrix(readRDS(file.path(BULK, "bulk_RNA_studies.rds"))$tcga)
meta <- readRDS(file.path(BULK, "bulk_RNA_studies_metadata.rds"))$tcga[colnames(E0), , drop = FALSE]
cls  <- G$class[match(colnames(E0), G$sample)]
sarc <- as.numeric(meta$sarc_score)
TFS  <- setdiff(unique(read.csv("BAP1_TF_limma_sarcadjusted.csv", stringsAsFactors = FALSE)$TF), "BAP1")
expressed <- apply(E0, 1, function(x) all(is.finite(x)) && sd(x) > 0 && mean(x > 0) > 0.5)

cat("=== are IRF8 and EGR2 testable? ===\n")
for (g in c(GOI, "PTPRC")) cat(sprintf(
  "%-5s in TF panel: %-5s | in TCGA matrix: %-5s | expressed filter: %-5s | median log2 %.2f | %% >0: %.0f%%\n",
  g, g %in% TFS, g %in% rownames(E0), isTRUE(expressed[g]),
  if (g %in% rownames(E0)) median(E0[g, ]) else NA,
  if (g %in% rownames(E0)) 100 * mean(E0[g, ] > 0) else NA))

gn  <- getGeneAnnotation()$genes
loc <- data.frame(gene = mcols(gn)$symbol, chr = as.character(seqnames(gn)),
                  start = start(gn), stringsAsFactors = FALSE)
loc <- loc[!duplicated(loc$gene), ]
is3p <- function(x){ i <- match(x, loc$gene); !is.na(i) & loc$chr[i] == "chr3" & loc$start[i] < CEN3 }

## ---- generic fitter --------------------------------------------------------------
contrasts <- list(
  A = list(lost = "ALTERED", ret = "WILD-TYPE",
           desc = "ALTERED vs WILD-TYPE"),
  B = list(lost = c("ALTERED", "HEMIZYGOUS", "VUS"), ret = "WILD-TYPE",
           desc = "any BAP1 lesion vs WILD-TYPE"),
  C = list(lost = "ALTERED", ret = c("WILD-TYPE", "HEMIZYGOUS", "VUS"),
           desc = "ALTERED vs all other"))
fit1 <- function(k, universe, adjust){
  ct <- contrasts[[k]]
  keep <- cls %in% c(ct$lost, ct$ret) & is.finite(sarc)
  grp  <- factor(ifelse(cls[keep] %in% ct$lost, "lost", "retained"), levels = c("retained","lost"))
  S    <- sarc[keep]
  genes <- names(expressed)[expressed]
  if (universe == "TF") genes <- intersect(genes, TFS)
  E <- E0[genes, keep, drop = FALSE]
  des <- if (adjust) model.matrix(~ grp + S) else model.matrix(~ grp)
  f <- eBayes(lmFit(E, des), trend = TRUE, robust = TRUE)
  t <- topTable(f, coef = "grplost", number = Inf, sort.by = "none")
  t$gene <- rownames(t)
  attr(t, "n") <- c(lost = sum(grp == "lost"), retained = sum(grp == "retained"))
  t
}

## ---- primary: contrast A, TF genes only, UNADJUSTED ---------------------------------
tA <- fit1("A", "TF", FALSE)
R  <- data.frame(TF = tA$gene, chr3p = is3p(tA$gene), logFC = tA$logFC, t = tA$t,
                 P = tA$P.Value, FDR = tA$adj.P.Val, stringsAsFactors = FALSE)
adjA <- fit1("A", "TF", TRUE)
R$t_sarcadjusted   <- adjA$t[match(R$TF, adjA$gene)]
R$FDR_sarcadjusted <- adjA$adj.P.Val[match(R$TF, adjA$gene)]
R <- R[order(R$P), ]
write.csv(R, "TCGA_genetic_BAP1_TFonly_unadjusted.csv", row.names = FALSE)

cat(sprintf("\n=== TF genes only, NO histology correction: %d altered vs %d wild-type, %d TFs ===\n",
            attr(tA, "n")["lost"], attr(tA, "n")["retained"], nrow(R)))
for (x in c(0.01, 0.05, 0.10)) cat(sprintf(
  "FDR < %.2f : %3d TFs (up %d, down %d, chr3p %d) | with correction: %d\n", x,
  sum(R$FDR < x), sum(R$FDR < x & R$t > 0), sum(R$FDR < x & R$t < 0),
  sum(R$FDR < x & R$chr3p), sum(R$FDR_sarcadjusted < x)))
cat(sprintf("agreement of t with vs without correction: Pearson %.3f\n",
            cor(R$t, R$t_sarcadjusted)))
cat("\nTFs at FDR < 0.05 WITHOUT correction:\n")
print(R[R$FDR < 0.05, c("TF","chr3p","logFC","t","P","FDR","t_sarcadjusted","FDR_sarcadjusted")],
      row.names = FALSE, digits = 3)
only_unadj <- R$TF[R$FDR < 0.05 & R$FDR_sarcadjusted >= 0.05]
only_adj   <- R$TF[R$FDR >= 0.05 & R$FDR_sarcadjusted < 0.05]
cat("\nsignificant only WITHOUT correction:", if (length(only_unadj)) paste(only_unadj, collapse = ", ") else "none", "\n")
cat("significant only WITH correction   :", if (length(only_adj)) paste(only_adj, collapse = ", ") else "none", "\n")

## ---- IRF8 / EGR2 across definitions --------------------------------------------------
SENS <- do.call(rbind, lapply(names(contrasts), function(k)
  do.call(rbind, lapply(c("TF", "genome"), function(u)
    do.call(rbind, lapply(c(FALSE, TRUE), function(adj){
      t <- fit1(k, u, adj)
      do.call(rbind, lapply(c(GOI, "PTPRC"), function(g){
        r <- t[t$gene == g, ]
        if (!nrow(r)) return(NULL)
        data.frame(contrast = k, definition = contrasts[[k]]$desc,
                   n_lost = attr(t, "n")["lost"], n_ret = attr(t, "n")["retained"],
                   universe = u, sarc_adjusted = adj, gene = g,
                   logFC = r$logFC, t = r$t, P = r$P.Value, FDR = r$adj.P.Val,
                   stringsAsFactors = FALSE)
      })) }))))))
rownames(SENS) <- NULL
write.csv(SENS, "TCGA_BAP1_IRF8_EGR2_sensitivity.csv", row.names = FALSE)
cat("\n=== IRF8, EGR2 (and PTPRC as immune reference) under each definition ===\n")
cat("    logFC > 0 = higher in BAP1-altered\n")
print(transform(SENS[, c("contrast","n_lost","n_ret","universe","sarc_adjusted","gene","logFC","t","P","FDR")],
                logFC = round(logFC, 3), t = round(t, 2), P = signif(P, 3), FDR = signif(FDR, 3)),
      row.names = FALSE)

cat("\n=== means by genetic class (log2 expression) ===\n")
M <- data.frame(class = cls, IRF8 = E0["IRF8", ], EGR2 = if ("EGR2" %in% rownames(E0)) E0["EGR2", ] else NA,
                PTPRC = E0["PTPRC", ])
print(aggregate(cbind(IRF8, EGR2, PTPRC) ~ class, M, function(x) round(mean(x), 2)), row.names = FALSE)
cat(sprintf("\nIRF8 vs PTPRC across all 87 tumours: Spearman rho = %.3f\n",
            cor(M$IRF8, M$PTPRC, method = "spearman")))
if (!all(is.na(M$EGR2)))
  cat(sprintf("EGR2 vs PTPRC across all 87 tumours: Spearman rho = %.3f\n",
              cor(M$EGR2, M$PTPRC, method = "spearman")))

## ---- volcano (contrast A, TF-only, unadjusted) ---------------------------------------
R$rank <- rank(R$P, ties.method = "min")
R$lab  <- ifelse(R$rank <= NTOP | R$TF %in% GOI, R$TF, NA)
lim <- max(abs(R$logFC))
pv <- ggplot(R, aes(logFC, -log10(P))) +
  geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = .3, colour = "grey55") +
  geom_point(colour = "grey80", size = 1, alpha = .7) +
  geom_point(data = subset(R, FDR < 0.05), colour = "grey40", size = 1.2) +
  geom_point(data = subset(R, FDR < 0.05 & chr3p), colour = "#b2182b", size = 1.3) +
  geom_point(data = subset(R, rank <= NTOP), shape = 21, fill = NA, colour = "grey5", size = 1.9) +
  geom_point(data = subset(R, TF %in% GOI), shape = 23, fill = "#f5b041", colour = "grey10", size = 2.6) +
  geom_text_repel(aes(label = lab), size = 2.3, colour = "grey10", max.overlaps = 40,
                  seed = 1, min.segment.length = 0, segment.size = .2) +
  scale_x_continuous(limits = c(-lim, lim)) +
  gtheme_no_rot +
  xlab("BAP1 altered - wild-type (log2, NO histology correction)") +
  ylab(expression(-log[10]~"moderated"~italic(P))) +
  ggtitle(sprintf("TCGA MESO, TF genes only (n = %d): %d altered vs %d wild-type",
                  nrow(R), attr(tA, "n")["lost"], attr(tA, "n")["retained"]),
          subtitle = "dark = FDR < 0.05; red = chr3p; orange diamonds = IRF8, EGR2 (reported in the TCGA paper)")
pdf(file.path("Plots", "TCGA_genetic_BAP1_TFonly_unadjusted_volcano.pdf"), width = 6.6, height = 5.1)
print(pv); dev.off()

## ---- IRF8 / EGR2 / PTPRC by genetic class --------------------------------------------
L <- do.call(rbind, lapply(c(GOI, "PTPRC"), function(g) if (g %in% rownames(E0))
  data.frame(gene = g, class = cls, expr = E0[g, ], stringsAsFactors = FALSE)))
L$class <- factor(L$class, levels = c("WILD-TYPE", "VUS", "HEMIZYGOUS", "ALTERED"))
L$gene  <- factor(L$gene, levels = c(GOI, "PTPRC"))
pb <- ggplot(L, aes(class, expr)) +
  geom_boxplot(outlier.shape = NA, width = .6, linewidth = .3, fill = "grey92") +
  geom_jitter(width = .15, size = .8, alpha = .7) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  gtheme_no_rot + xlab(NULL) + ylab("log2 expression") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  ggtitle("TCGA MESO: IRF8, EGR2 and PTPRC by genetic BAP1 class")
pdf(file.path("Plots", "TCGA_BAP1_IRF8_EGR2_by_class.pdf"), width = 7.5, height = 3.4)
print(pb); dev.off()

cat("\nDONE\n")
