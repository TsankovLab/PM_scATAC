###############################################################################
# R2_Q14 -- filter TF motif activity by gene-score concordance, then re-test
#           congruence with bulk RNA.
#
# Idea being tested.  A chromVAR deviation says a motif family's sites are more
# accessible; it does not say the TF is present.  Gene score is the ATAC readout of
# the TF's own locus.  Keeping only motifs whose BAP1 effect points the SAME WAY in
# gene score should retain the TFs plausibly active and discard family-level motif
# signal with no TF behind it -- and gene score is the readout that did agree with
# bulk (rho 0.26-0.33 vs -0.01 for motifs), so the filtered motif set should inherit
# some of that agreement if the idea is right.
#
# Both readouts are refitted with BOTH covariates, as asked:
#     ~ BAP1_status + sarc_score + chr3p_accessibility          (9 tumours)
#
# TWO WARNINGS, both load-bearing.
#
# 1. ADJUSTING FOR chr3p IS PARTLY OVER-ADJUSTMENT.  chr3p deletion is the usual
#    MECHANISM of BAP1 loss (BAP1 is at 3p21.1), so chr3p dosage is not a nuisance
#    variable sitting beside BAP1 status -- it is on the causal path.  Conditioning
#    on it removes real BAP1 effect along with the positional artefact.  The design
#    is also collinear (correlation of BAP1 status with chr3p accessibility is
#    printed below) and leaves only 5 residual degrees of freedom.  The no-chr3p
#    model is therefore reported alongside, and neither is "the" answer.
#
# 2. THE FILTER CAN MANUFACTURE ITS OWN RESULT.  Selecting motifs that agree in sign
#    with gene score, and then asking whether those motifs agree with bulk, conditions
#    on a variable that itself correlates with bulk.  Some improvement is guaranteed
#    by construction.  A permutation null is therefore built: the gene-score t values
#    are shuffled across TFs, the same sign filter is applied, and the congruence is
#    recomputed.  The observed value only means something if it exceeds that null.
#
# Input : BAP1_TFcongruence_bulk_molecular.csv (bulk), TFgenescore_per_sample.csv,
#         TFactivity_deviation_per_sample.csv, BAP1_genescore_per_sample.csv,
#         scATAC_sarcscore_per_sample.csv, allsample_scATAC_bap1_chr3p.csv
# Output: BAP1_TF_activity_gsfiltered.csv, BAP1_TF_gsfilter_congruence.csv,
#         Plots/BAP1_TF_gsfilter_congruence.pdf
###############################################################################
suppressPackageStartupMessages({
  library(limma); library(ggplot2); library(ggrepel); library(ArchR)
})
addArchRThreads(1); addArchRGenome("hg38")

ROOT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUTDIR <- file.path(ROOT, "git_repo_claude", "R2_Q14")
setwd(OUTDIR); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))

RUNX <- c("RUNX1","RUNX2","RUNX3"); IRF <- paste0("IRF", 1:9)
MIN_CELLS <- 50; CEN3 <- 91.1e6; NPERM <- 2000

## ---- design ------------------------------------------------------------------------
gs <- read.csv("BAP1_genescore_per_sample.csv", stringsAsFactors = FALSE)
sc <- read.csv("scATAC_sarcscore_per_sample.csv", stringsAsFactors = FALSE)
c3 <- read.csv("allsample_scATAC_bap1_chr3p.csv", stringsAsFactors = FALSE)
chr3p <- tapply(c3$chr3p_access, c3$sample, mean)

ACT <- read.csv("TFactivity_deviation_per_sample.csv", row.names = 1, check.names = FALSE)
GSC <- read.csv("TFgenescore_per_sample.csv",          row.names = 1, check.names = FALSE)

sams <- Reduce(intersect, list(colnames(ACT), colnames(GSC), gs$Sample, sc$sample,
                               names(chr3p)))
ncell <- setNames(sc$n_cells, sc$sample)[sams]
sams  <- sams[ncell >= MIN_CELLS]
BAP1  <- factor(gs$BAP1_status[match(sams, gs$Sample)], levels = c("retained","lost"))
SARC  <- sc$sarc_score_atac[match(sams, sc$sample)]
CHR3P <- as.numeric(chr3p[sams])
cat("tumours:", length(sams), "->", paste(sams, collapse = ", "), "\n")
cat(sprintf("BAP1 lost %d / retained %d\n", sum(BAP1 == "lost"), sum(BAP1 == "retained")))
cat(sprintf("collinearity: cor(BAP1 lost, chr3p accessibility) = %.3f  <- see warning 1\n",
            cor(as.integer(BAP1 == "lost"), CHR3P)))
cat(sprintf("             cor(sarc score, chr3p accessibility) = %.3f\n",
            cor(SARC, CHR3P)))

fitit <- function(M, use_chr3p){
  E <- as.matrix(M[, sams, drop = FALSE])
  E <- E[apply(E, 1, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE]
  des <- if (use_chr3p) model.matrix(~ BAP1 + SARC + CHR3P) else model.matrix(~ BAP1 + SARC)
  f <- eBayes(lmFit(E, des), robust = TRUE)
  t <- topTable(f, coef = "BAP1lost", number = Inf, sort.by = "none")
  data.frame(TF = rownames(t), t = t$t, P = t$P.Value, FDR = t$adj.P.Val,
             stringsAsFactors = FALSE)
}

## ---- the two readouts, both models -------------------------------------------------
mot3 <- fitit(ACT, TRUE);  motN <- fitit(ACT, FALSE)
gsc3 <- fitit(GSC, TRUE);  gscN <- fitit(GSC, FALSE)
cat(sprintf("\nmotif  : %d TFs | gene score: %d TFs\n", nrow(mot3), nrow(gsc3)))

gn  <- getGeneAnnotation()$genes
loc <- data.frame(gene = mcols(gn)$symbol, chr = as.character(seqnames(gn)),
                  start = start(gn), stringsAsFactors = FALSE)
loc <- loc[!duplicated(loc$gene), ]
is3p <- function(x) { i <- match(x, loc$gene)
  !is.na(loc$chr[i]) & loc$chr[i] == "chr3" & loc$start[i] < CEN3 }

## detection: the TF's own locus must actually be accessible, else "same direction"
## is comparing two noise terms.  Floor = median gene score across the panel.
mean_gs <- rowMeans(GSC[, sams, drop = FALSE])
FLOOR   <- median(mean_gs, na.rm = TRUE)

build <- function(mot, gsc, tag){
  D <- merge(setNames(mot, c("TF","t_motif","P_motif","FDR_motif")),
             setNames(gsc, c("TF","t_gs","P_gs","FDR_gs")), by = "TF")
  D$mean_gs   <- mean_gs[match(D$TF, names(mean_gs))]
  D$chr3p     <- is3p(D$TF)
  D$detected  <- !is.na(D$mean_gs) & D$mean_gs >= FLOOR
  D$same_sign <- sign(D$t_motif) == sign(D$t_gs)
  D$keep      <- D$detected & D$same_sign & !D$chr3p
  D$model     <- tag
  D
}
D3 <- build(mot3, gsc3, "BAP1 + sarc + chr3p")
DN <- build(motN, gscN, "BAP1 + sarc")
write.csv(rbind(D3, DN), "BAP1_TF_activity_gsfiltered.csv", row.names = FALSE)

for (D in list(D3, DN)) {
  cat(sprintf("\n--- %s ---\n", D$model[1]))
  cat(sprintf("  TFs paired               : %d\n", nrow(D)))
  cat(sprintf("  detected (gene score >= %.2f): %d\n", FLOOR, sum(D$detected)))
  cat(sprintf("  same direction in both   : %d of %d detected (%.0f%%, 50%% expected by chance)\n",
              sum(D$detected & D$same_sign), sum(D$detected),
              100 * mean(D$same_sign[D$detected])))
  cat(sprintf("  retained after chr3p drop: %d\n", sum(D$keep)))
}

## ---- congruence with bulk -------------------------------------------------------------
bulk <- read.csv("BAP1_TFcongruence_bulk_molecular.csv", stringsAsFactors = FALSE)
tb   <- function(s) setNames(bulk[bulk$cohort == s, c("TF","t")], c("TF", "t_bulk"))
COH  <- c("bueno","tcga","mesomics")

rho1 <- function(D, s, subset_keep){
  d <- merge(if (subset_keep) D[D$keep, ] else D, tb(s), by = "TF")
  if (nrow(d) < 10) return(c(rho = NA, n = nrow(d)))
  c(rho = suppressWarnings(cor(d$t_motif, d$t_bulk, method = "spearman")), n = nrow(d))
}
S <- do.call(rbind, lapply(list(D3, DN), function(D)
  do.call(rbind, lapply(COH, function(s){
    a <- rho1(D, s, FALSE); b <- rho1(D, s, TRUE)
    ## null: shuffle the gene-score t across TFs, refilter, recompute
    set.seed(1)
    nul <- replicate(NPERM, {
      Dp <- D; Dp$t_gs <- sample(D$t_gs)
      Dp$keep <- Dp$detected & (sign(Dp$t_motif) == sign(Dp$t_gs)) & !Dp$chr3p
      rho1(Dp, s, TRUE)["rho"] })
    data.frame(model = D$model[1], cohort = s,
               rho_all = a["rho"], n_all = a["n"],
               rho_filtered = b["rho"], n_filtered = b["n"],
               null_mean = mean(nul, na.rm = TRUE),
               null_q95  = quantile(nul, .95, na.rm = TRUE, names = FALSE),
               p_vs_null = mean(nul >= b["rho"], na.rm = TRUE),
               stringsAsFactors = FALSE) }))))
rownames(S) <- NULL
write.csv(S, "BAP1_TF_gsfilter_congruence.csv", row.names = FALSE)
cat("\n=== does the gene-score filter make motif activity agree with bulk? ===\n")
cat("    rho_all = every TF; rho_filtered = only gene-score-concordant, non-chr3p TFs\n")
cat("    null = same filter with the gene-score t shuffled (selection effect alone)\n")
print(transform(S, rho_all = round(rho_all, 3), rho_filtered = round(rho_filtered, 3),
                null_mean = round(null_mean, 3), null_q95 = round(null_q95, 3),
                p_vs_null = signif(p_vs_null, 3)), row.names = FALSE)

## ---- what survives the filter -------------------------------------------------------
K <- D3[D3$keep, ]; K <- K[order(K$P_motif), ]
cat(sprintf("\n=== top 20 gene-score-concordant TFs (%s) ===\n", D3$model[1]))
print(head(K[, c("TF","t_motif","P_motif","FDR_motif","t_gs","P_gs","mean_gs")], 20),
      row.names = FALSE, digits = 3)
cat("\nRUNX / IRF retained by the filter:",
    paste(intersect(K$TF, c(RUNX, IRF)), collapse = ", "), "\n")

## ---- figure ---------------------------------------------------------------------------
pl <- lapply(COH, function(s){
  d <- merge(D3, tb(s), by = "TF")
  r_all <- cor(d$t_motif, d$t_bulk, method = "spearman")
  dk <- d[d$keep, ]; r_k <- cor(dk$t_motif, dk$t_bulk, method = "spearman")
  d$lab <- ifelse(d$TF %in% c(RUNX, IRF), d$TF, NA)
  ggplot(d, aes(t_motif, t_bulk)) +
    geom_hline(yintercept = 0, linewidth = .2, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
    geom_point(colour = "grey82", size = .9, alpha = .7) +
    geom_point(data = dk, colour = "#1b4f72", size = 1.3) +
    geom_text_repel(aes(label = lab), size = 2, colour = "grey25",
                    max.overlaps = 30, seed = 1, min.segment.length = 0) +
    gtheme_no_rot +
    xlab("t, scATAC motif activity (BAP1, adjusted)") +
    ylab(sprintf("t, %s mRNA (BAP1, adjusted)", s)) +
    ggtitle(sprintf("%s: all TFs rho = %.2f | gene-score-concordant (blue, n = %d) rho = %.2f",
                    s, r_all, nrow(dk), r_k))
})
pdf(file.path("Plots", "BAP1_TF_gsfilter_congruence.pdf"), width = 5.6, height = 4.4)
for (p in pl) print(p); dev.off()

cat("\nDONE\n")
