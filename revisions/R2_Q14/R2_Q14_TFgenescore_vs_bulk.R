###############################################################################
# R2_Q14 -- repeat the scATAC BAP1 test on GENE SCORE instead of motif activity,
#           and re-test congruence with bulk RNA.
#
# Motivation.  The motif-activity result showed no congruence with any bulk cohort
# (rho -0.07 to -0.01, sign agreement at chance) while the three bulk cohorts agreed
# with each other at rho 0.49-0.64.  One candidate explanation is the MODALITY: a
# chromVAR deviation is family-level binding across thousands of peaks, whereas bulk
# measures the mRNA of one gene.  ArchR gene score is the ATAC readout that is meant
# to approximate expression of a specific gene, so if modality is the obstacle,
# swapping motif activity for gene score should recover some congruence.  If it does
# not, modality is not the explanation and power or the BAP1 label is.
#
#     gene score(TF, tumour) ~ BAP1_status + sarc_score      (limma, 9 tumours)
#
# Same cohort rules as the motif analysis: clinical BAP1, malignant cells only,
# normal1 / P11_HOX dropped, and tumours under MIN_CELLS malignant cells excluded
# (P3 21 cells, P13 18).
#
# TWO CAVEATS THAT ARE NOT SYMMETRIC WITH BULK.
#   1. Gene score inherits copy number.  A gene on a gained arm scores higher for
#      reasons unrelated to regulation -- measured at 3-7x on-arm inflation in R2_Q3.
#   2. BAP1 loss co-occurs with chr3p deletion, so chr3p TFs are expected to show a
#      positional gene-score difference by BAP1 status in exactly the way they did
#      in the bulk mRNA proxy.  chr3p TFs are flagged and congruence is reported with
#      and without them.
#
# Input : tumor_compartment/scatac_ArchR, BAP1_genescore_per_sample.csv,
#         scATAC_sarcscore_per_sample.csv, BAP1_TF_limma_sarcadjusted.csv (motifs),
#         BAP1_TFcongruence_bulk_molecular.csv (bulk TF results)
# Output: BAP1_TFgenescore_limma.csv, BAP1_TFgenescore_congruence.csv,
#         Plots/BAP1_TFgenescore_congruence.pdf
###############################################################################
suppressPackageStartupMessages({
  library(ArchR); library(limma); library(SummarizedExperiment)
  library(ggplot2); library(ggrepel)
})
addArchRThreads(4); addArchRGenome("hg38")

ROOT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUTDIR <- file.path(ROOT, "git_repo_claude", "R2_Q14")
TCOMP  <- file.path(ROOT, "tumor_compartment")
setwd(OUTDIR); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))

RUNX <- c("RUNX1","RUNX2","RUNX3"); IRF <- paste0("IRF", 1:9)
MIN_CELLS <- 50
DROP <- c("normal1", "P11_HOX")
CEN3 <- 91.1e6
CACHE <- "TFgenescore_per_sample.csv"

at  <- read.csv("BAP1_TF_limma_sarcadjusted.csv", stringsAsFactors = FALSE)
TFS <- unique(at$TF)

## ---- per-sample TF gene score (cached) ------------------------------------------
if (file.exists(CACHE)) {
  GS <- read.csv(CACHE, row.names = 1, check.names = FALSE); cat("gene score from cache\n")
} else {
  proj <- loadArchRProject(file.path(TCOMP, "scatac_ArchR"), showLogo = FALSE)
  proj@projectMetadata$outputDirectory <- file.path(OUTDIR, "archr_scratch")
  dir.create(file.path(OUTDIR, "archr_scratch"), showWarnings = FALSE)
  cd  <- getCellColData(proj)
  grp <- if ("Sample3" %in% colnames(cd)) "Sample3" else "Sample"
  keep <- rownames(cd)[!cd[[grp]] %in% DROP & !grepl("^RPL", cd[[grp]])]
  proj <- proj[keep, ]
  se <- getGroupSE(proj, useMatrix = "GeneScoreMatrix", groupBy = grp,
                   divideN = TRUE, scaleTo = 10000)
  G <- log2(as.matrix(assay(se)) + 1)
  rownames(G) <- rowData(se)$name
  G <- G[!duplicated(rownames(G)), , drop = FALSE]
  GS <- G[rownames(G) %in% TFS, , drop = FALSE]
  write.csv(round(GS, 5), CACHE)
  rm(proj, se, G); invisible(gc(FALSE))
}
cat("TF gene scores:", nrow(GS), "TFs x", ncol(GS), "tumours\n")

## ---- design ----------------------------------------------------------------------
gs <- read.csv("BAP1_genescore_per_sample.csv", stringsAsFactors = FALSE)
sc <- read.csv("scATAC_sarcscore_per_sample.csv", stringsAsFactors = FALSE)
sams <- intersect(colnames(GS), intersect(gs$Sample, sc$sample))
ncell <- setNames(sc$n_cells, sc$sample)[sams]
small <- sams[ncell < MIN_CELLS]
cat("excluded for <", MIN_CELLS, "cells:",
    if (length(small)) paste(sprintf("%s (%d)", small, ncell[small]), collapse = ", ")
    else "none", "\n")
sams <- sams[ncell >= MIN_CELLS]
E    <- as.matrix(GS[, sams, drop = FALSE])
E    <- E[apply(E, 1, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE]
BAP1 <- factor(gs$BAP1_status[match(sams, gs$Sample)], levels = c("retained", "lost"))
SARC <- sc$sarc_score_atac[match(sams, sc$sample)]
cat("tumours:", length(sams), "|", sum(BAP1 == "lost"), "lost,",
    sum(BAP1 == "retained"), "retained | TFs:", nrow(E), "\n")

fit <- eBayes(lmFit(E, model.matrix(~ BAP1 + SARC)), robust = TRUE)
a <- topTable(fit, coef = "BAP1lost", number = Inf, sort.by = "none")
h <- topTable(fit, coef = "SARC",     number = Inf, sort.by = "none")
gn  <- getGeneAnnotation()$genes
loc <- data.frame(gene = mcols(gn)$symbol, chr = as.character(seqnames(gn)),
                  start = start(gn), stringsAsFactors = FALSE)
loc <- loc[!duplicated(loc$gene), ]
R <- data.frame(TF = rownames(a), logFC_gs = a$logFC, t_gs = a$t,
                P_gs = a$P.Value, FDR_gs = a$adj.P.Val,
                t_sarc_gs = h$t[match(rownames(a), rownames(h))],
                stringsAsFactors = FALSE)
R$chr3p <- with(loc[match(R$TF, loc$gene), ], !is.na(chr) & chr == "chr3" & start < CEN3)
R <- R[order(R$P_gs), ]
write.csv(R, "BAP1_TFgenescore_limma.csv", row.names = FALSE)

cat("\n=== TFs passing, gene score ===\n")
for (aa in c(0.05, 0.10, 0.20)) cat(sprintf("FDR < %.2f : %d TFs\n", aa, sum(R$FDR_gs < aa)))
cat("\n=== top 15 ===\n")
print(head(R[, c("TF","logFC_gs","t_gs","P_gs","FDR_gs","t_sarc_gs","chr3p")], 15),
      row.names = FALSE, digits = 3)
cat("\n=== RUNX / IRF ===\n")
print(R[R$TF %in% c(RUNX, IRF), c("TF","logFC_gs","t_gs","P_gs","FDR_gs","t_sarc_gs")],
      row.names = FALSE, digits = 3)

## ---- congruence -------------------------------------------------------------------
bulk <- read.csv("BAP1_TFcongruence_bulk_molecular.csv", stringsAsFactors = FALSE)
tb <- function(s) setNames(bulk[bulk$cohort == s, c("TF","t")], c("TF", paste0("t_", s)))
mot <- setNames(at[, c("TF","t_adj")], c("TF","t_motif"))
cmp <- function(a, b, la, lb, drop3p){
  m <- merge(a, b, by = "TF")
  if (drop3p) m <- m[!m$TF %in% R$TF[R$chr3p], ]
  ct <- suppressWarnings(cor.test(m[[2]], m[[3]], method = "spearman"))
  data.frame(comparison = paste(la, "vs", lb),
             chr3p = if (drop3p) "excluded" else "included", n_TF = nrow(m),
             rho = unname(ct$estimate), p = ct$p.value,
             pct_same_sign = 100 * mean(sign(m[[2]]) == sign(m[[3]])),
             stringsAsFactors = FALSE)
}
gsv <- setNames(R[, c("TF","t_gs")], c("TF","t_genescore"))
S <- do.call(rbind, lapply(c(FALSE, TRUE), function(d) rbind(
  cmp(gsv, tb("bueno"),    "scATAC gene score", "bueno",    d),
  cmp(gsv, tb("tcga"),     "scATAC gene score", "tcga",     d),
  cmp(gsv, tb("mesomics"), "scATAC gene score", "mesomics", d),
  cmp(gsv, mot,            "scATAC gene score", "scATAC motif", d),
  cmp(mot, tb("bueno"),    "scATAC motif",      "bueno",    d),
  cmp(mot, tb("tcga"),     "scATAC motif",      "tcga",     d),
  cmp(mot, tb("mesomics"), "scATAC motif",      "mesomics", d))))
write.csv(S, "BAP1_TFgenescore_congruence.csv", row.names = FALSE)
cat("\n=== congruence: does gene score agree with bulk better than motif activity? ===\n")
print(transform(S, rho = round(rho, 3), p = signif(p, 3),
                pct_same_sign = round(pct_same_sign)), row.names = FALSE)

## ---- figure -------------------------------------------------------------------------
D <- Reduce(function(x, y) merge(x, y, by = "TF"),
            list(gsv, mot, tb("bueno"), tb("tcga"), tb("mesomics")))
D$chr3p <- D$TF %in% R$TF[R$chr3p]
mk <- function(x, y, xl, yl){
  d <- D; r <- cor(d[[x]][!d$chr3p], d[[y]][!d$chr3p], method = "spearman")
  d$lab <- ifelse(d$TF %in% c(RUNX, IRF) | d$chr3p, d$TF, NA)
  ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_hline(yintercept = 0, linewidth = .2, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
    geom_point(colour = "grey78", size = 1, alpha = .8) +
    geom_point(data = subset(d, chr3p), colour = "#b2182b", size = 1.6) +
    geom_text_repel(aes(label = lab), size = 2, colour = "grey25",
                    max.overlaps = 30, seed = 1, min.segment.length = 0) +
    gtheme_no_rot + xlab(xl) + ylab(yl) +
    ggtitle(sprintf("rho = %.2f (chr3p TFs red, excluded)", r))
}
pdf(file.path("Plots", "BAP1_TFgenescore_congruence.pdf"), width = 5, height = 4.2)
print(mk("t_genescore", "t_bueno",    "t, scATAC TF gene score", "t, Bueno mRNA"))
print(mk("t_genescore", "t_tcga",     "t, scATAC TF gene score", "t, TCGA mRNA"))
print(mk("t_genescore", "t_mesomics", "t, scATAC TF gene score", "t, MESOMICS mRNA"))
print(mk("t_genescore", "t_motif",    "t, scATAC TF gene score", "t, scATAC motif activity"))
dev.off()

cat("\nDONE\n")
