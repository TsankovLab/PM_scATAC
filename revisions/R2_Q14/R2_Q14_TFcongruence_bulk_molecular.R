###############################################################################
# R2_Q14 -- Bueno and TCGA, molecular BAP1 proxy, restricted to the TFs assayed in
#           our scATAC cohort, to test congruence with the scATAC result.
#
#     expression ~ BAP1_group + sarc_score      (limma, per cohort)
#     BAP1_group = bottom vs top tertile of BAP1 mRNA (molecular inference)
#     features   = only the TF genes matching our 869 chromVAR motifs
#
# Restricting to our TFs before correcting changes the multiple-testing burden from
# ~18,000 genes to ~650, which is the correct correction for a focused question and
# is far more sensitive than reading our TFs out of a genome-wide FDR list.
#
# MANDATORY chr3p CONTROL.  Splitting on BAP1 mRNA selects tumours with chr3p
# deletion (BAP1 is at 3p21.1), so every gene on the arm drops with it: of the 304
# genes reproducing across all three cohorts, 82 (28%) are chr3p against a 2.9%
# background, OR 13.4, p = 2.7e-53, all in the same direction.  Any TF on chr3p is
# therefore expected to "differ by BAP1 status" for purely positional reasons.
# chr3p TFs are flagged and every congruence statistic is reported twice, with and
# without them.  The without-chr3p number is the one that means anything.
#
# Congruence is measured against the scATAC histology-adjusted motif result, which
# is a different measurement (motif accessibility, not mRNA) in different tumours,
# so agreement is corroboration and disagreement is weak evidence either way --
# especially since nothing in the scATAC fit passes FDR.
#
# Input : bulkRNA_meso/*, BAP1_TF_limma_sarcadjusted.csv (scATAC, 9 tumours)
# Output: BAP1_TFcongruence_bulk_molecular.csv, BAP1_TFcongruence_summary.csv,
#         Plots/BAP1_TFcongruence.pdf
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

RUNX <- c("RUNX1", "RUNX2", "RUNX3"); IRF <- paste0("IRF", 1:9)
CEN3 <- 91.1e6                      # hg38 chr3 centromere
STUDIES <- c("bueno", "tcga")       # mesomics carried along as the reference cohort

## ---- our TF panel + chr3p annotation ---------------------------------------------
at  <- read.csv("BAP1_TF_limma_sarcadjusted.csv", stringsAsFactors = FALSE)
TFS <- unique(at$TF)
gn  <- getGeneAnnotation()$genes
loc <- data.frame(gene = mcols(gn)$symbol, chr = as.character(seqnames(gn)),
                  start = start(gn), stringsAsFactors = FALSE)
loc <- loc[!duplicated(loc$gene), ]
loc$chr3p <- loc$chr == "chr3" & loc$start < CEN3
cat("TFs assayed in our scATAC cohort:", length(TFS), "\n")
cat("of these, on chr3p:", sum(loc$chr3p[match(TFS, loc$gene)], na.rm = TRUE), "->",
    paste(TFS[which(loc$chr3p[match(TFS, loc$gene)])], collapse = ", "), "\n")

bl   <- readRDS(file.path(BULK, "bulk_RNA_studies.rds"))
meta <- readRDS(file.path(BULK, "bulk_RNA_studies_metadata.rds"))

## ---- per cohort, TF-restricted -------------------------------------------------------
run <- function(s){
  E  <- as.matrix(bl[[s]]); md <- meta[[s]]
  md <- if (all(colnames(E) %in% rownames(md))) md[colnames(E), , drop = FALSE]
        else md[match(colnames(E), md$Sample), , drop = FALSE]
  bap1 <- as.numeric(E["BAP1", ]); sarc <- as.numeric(md$sarc_score)
  q   <- quantile(bap1, c(1/3, 2/3), na.rm = TRUE)
  grp <- ifelse(bap1 <= q[1], "lost", ifelse(bap1 >= q[2], "retained", NA))
  ok  <- !is.na(grp) & is.finite(sarc)
  G <- factor(grp[ok], levels = c("retained", "lost")); S <- sarc[ok]
  ## restrict to our TFs BEFORE fitting, so the FDR is over the TF panel
  keep <- rownames(E) %in% setdiff(TFS, "BAP1")
  Ei <- E[keep, ok, drop = FALSE]
  keep2 <- apply(Ei, 1, function(x) all(is.finite(x)) && sd(x) > 0 && mean(x > 0) > 0.5)
  Ei <- Ei[keep2, , drop = FALSE]
  fit <- eBayes(lmFit(Ei, model.matrix(~ G + S)), trend = TRUE, robust = TRUE)
  a <- topTable(fit, coef = "Glost", number = Inf, sort.by = "none")
  h <- topTable(fit, coef = "S",     number = Inf, sort.by = "none")
  cat(sprintf("\n--- %s: %d tumours (%d lost / %d retained), %d of our TFs measured ---\n",
              s, ncol(Ei), sum(G == "lost"), sum(G == "retained"), nrow(Ei)))
  cat(sprintf("    sarc score %.2f lost vs %.2f retained\n",
              mean(S[G == "lost"]), mean(S[G == "retained"])))
  R <- data.frame(cohort = s, TF = rownames(a), logFC = a$logFC, t = a$t,
                  P = a$P.Value, FDR = a$adj.P.Val,
                  t_sarc = h$t[match(rownames(a), rownames(h))],
                  stringsAsFactors = FALSE)
  R$chr3p <- loc$chr3p[match(R$TF, loc$gene)]
  R$chr3p[is.na(R$chr3p)] <- FALSE
  for (aa in c(0.05, 0.10)) cat(sprintf(
    "    FDR < %.2f over the TF panel: %d TFs (%d of them chr3p)\n", aa,
    sum(R$FDR < aa), sum(R$FDR < aa & R$chr3p)))
  R
}
res <- lapply(setNames(c(STUDIES, "mesomics"), c(STUDIES, "mesomics")), run)
ALL <- do.call(rbind, res); rownames(ALL) <- NULL
write.csv(ALL, "BAP1_TFcongruence_bulk_molecular.csv", row.names = FALSE)

cat("\n=== RUNX and IRF across cohorts (TF-panel FDR) ===\n")
print(ALL[ALL$TF %in% c(RUNX, IRF), c("cohort","TF","logFC","t","P","FDR","t_sarc")],
      row.names = FALSE, digits = 3)

## ---- congruence -----------------------------------------------------------------------
sc <- at[, c("TF", "t_adj", "P_adj")]; names(sc) <- c("TF", "t_scatac", "P_scatac")
cmp <- function(a, b, la, lb, drop3p){
  m <- merge(a, b, by = "TF")
  if (drop3p) m <- m[!m$TF %in% ALL$TF[ALL$chr3p], ]
  ct <- suppressWarnings(cor.test(m[[2]], m[[3]], method = "spearman"))
  data.frame(comparison = paste(la, "vs", lb),
             chr3p = if (drop3p) "excluded" else "included",
             n_TF = nrow(m), rho = unname(ct$estimate), p = ct$p.value,
             pct_same_sign = 100 * mean(sign(m[[2]]) == sign(m[[3]])),
             stringsAsFactors = FALSE)
}
tt <- function(s) setNames(res[[s]][, c("TF", "t")], c("TF", paste0("t_", s)))
S <- do.call(rbind, lapply(c(FALSE, TRUE), function(d) rbind(
  cmp(tt("bueno"),    tt("tcga"),     "bueno",    "tcga",     d),
  cmp(tt("bueno"),    tt("mesomics"), "bueno",    "mesomics", d),
  cmp(tt("tcga"),     tt("mesomics"), "tcga",     "mesomics", d),
  cmp(tt("bueno"),    sc[, 1:2],      "bueno",    "scATAC",   d),
  cmp(tt("tcga"),     sc[, 1:2],      "tcga",     "scATAC",   d),
  cmp(tt("mesomics"), sc[, 1:2],      "mesomics", "scATAC",   d))))
write.csv(S, "BAP1_TFcongruence_summary.csv", row.names = FALSE)
cat("\n=== congruence of the histology-adjusted BAP1 effect, TF panel only ===\n")
print(transform(S, rho = round(rho, 3), p = signif(p, 3),
                pct_same_sign = round(pct_same_sign)), row.names = FALSE)

## do bueno and tcga agree with each other more than either agrees with scATAC?
cat("\nreading: bulk-vs-bulk rho is the ceiling this design can reach;\n")
cat("bulk-vs-scATAC is the cross-modality congruence being tested against it.\n")

## ---- figure -------------------------------------------------------------------------
D <- merge(merge(tt("bueno"), tt("tcga"), by = "TF"), sc[, 1:2], by = "TF")
D$chr3p <- D$TF %in% ALL$TF[ALL$chr3p]
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
    ggtitle(sprintf("rho = %.2f (chr3p TFs in red, excluded)", r))
}
p1 <- mk("t_bueno", "t_tcga",   "t, Bueno (molecular BAP1)", "t, TCGA (molecular BAP1)")
p2 <- mk("t_bueno", "t_scatac", "t, Bueno (molecular BAP1)", "t, scATAC motif (clinical BAP1)")
p3 <- mk("t_tcga",  "t_scatac", "t, TCGA (molecular BAP1)",  "t, scATAC motif (clinical BAP1)")
pdf(file.path("Plots", "BAP1_TFcongruence.pdf"), width = 5, height = 4.2)
print(p1); print(p2); print(p3); dev.off()

cat("\nDONE\n")
