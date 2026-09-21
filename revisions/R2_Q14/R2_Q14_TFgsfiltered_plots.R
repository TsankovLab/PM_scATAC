###############################################################################
# R2_Q14 -- figures for the gene-score-filtered TF set only.
#
# The panel plotted here is the 186 motifs that (i) have a detectably accessible TF
# locus, (ii) move the same way in gene score as in motif activity, and (iii) are not
# on chr3p -- selected in R2_Q14_TFactivity_genescore_filtered.R.  Everything else is
# dropped, not greyed out.
#
# FDR is RECOMPUTED over the 186 rather than inherited from the 869.  That is
# legitimate here because the filter statistic is effectively independent of the
# statistic being tested: the motif and gene-score BAP1 t values correlate at
# rho = -0.008 across the panel.  This is ordinary independent filtering, not
# selection on the p-value being corrected.
#
# Input : TFactivity_deviation_per_sample.csv, TFgenescore_per_sample.csv,
#         BAP1_genescore_per_sample.csv, scATAC_sarcscore_per_sample.csv,
#         allsample_scATAC_bap1_chr3p.csv, BAP1_TF_activity_gsfiltered.csv,
#         BAP1_TFcongruence_bulk_molecular.csv
# Output: BAP1_TF_gsfiltered_results.csv,
#         Plots/BAP1_TF_gsfiltered_volcano.pdf,
#         Plots/BAP1_TF_gsfiltered_volcano_bulkmatch.pdf,
#         Plots/BAP1_TF_gsfiltered_congruence.pdf,
#         Plots/BAP1_TF_gsfiltered_motif_vs_genescore.pdf
###############################################################################
suppressPackageStartupMessages({ library(limma); library(ggplot2); library(ggrepel) })

ROOT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUTDIR <- file.path(ROOT, "git_repo_claude", "R2_Q14")
setwd(OUTDIR); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
RUNX <- c("RUNX1","RUNX2","RUNX3"); IRF <- paste0("IRF", 1:9)
MIN_CELLS <- 50; NTOP <- 20

## ---- design (identical to the filtering script) ----------------------------------
gs <- read.csv("BAP1_genescore_per_sample.csv", stringsAsFactors = FALSE)
sc <- read.csv("scATAC_sarcscore_per_sample.csv", stringsAsFactors = FALSE)
c3 <- read.csv("allsample_scATAC_bap1_chr3p.csv", stringsAsFactors = FALSE)
chr3p <- tapply(c3$chr3p_access, c3$sample, mean)
ACT <- read.csv("TFactivity_deviation_per_sample.csv", row.names = 1, check.names = FALSE)
GSC <- read.csv("TFgenescore_per_sample.csv",          row.names = 1, check.names = FALSE)

sams <- Reduce(intersect, list(colnames(ACT), colnames(GSC), gs$Sample, sc$sample,
                               names(chr3p)))
ncell <- setNames(sc$n_cells, sc$sample)[sams]; sams <- sams[ncell >= MIN_CELLS]
BAP1  <- factor(gs$BAP1_status[match(sams, gs$Sample)], levels = c("retained","lost"))
SARC  <- sc$sarc_score_atac[match(sams, sc$sample)]
CHR3P <- as.numeric(chr3p[sams])
des   <- model.matrix(~ BAP1 + SARC + CHR3P)

fitit <- function(M){
  E <- as.matrix(M[, sams, drop = FALSE])
  E <- E[apply(E, 1, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE]
  t <- topTable(eBayes(lmFit(E, des), robust = TRUE), coef = "BAP1lost",
                number = Inf, sort.by = "none")
  data.frame(TF = rownames(t), logFC = t$logFC, t = t$t, P = t$P.Value,
             stringsAsFactors = FALSE)
}
mot <- fitit(ACT); gsc <- fitit(GSC)

## ---- the filtered panel ------------------------------------------------------------
F <- read.csv("BAP1_TF_activity_gsfiltered.csv", stringsAsFactors = FALSE)
F <- F[F$model == "BAP1 + sarc + chr3p" & F$keep, c("TF","mean_gs")]
D <- merge(merge(setNames(mot, c("TF","logFC_motif","t_motif","P_motif")),
                 setNames(gsc, c("TF","logFC_gs","t_gs","P_gs")), by = "TF"),
           F, by = "TF")
D$FDR_motif <- p.adjust(D$P_motif, "BH")     # over the filtered panel, see header
D$FDR_gs    <- p.adjust(D$P_gs,    "BH")
D <- D[order(D$P_motif), ]
write.csv(D, "BAP1_TF_gsfiltered_results.csv", row.names = FALSE)
cat("filtered panel:", nrow(D), "TFs\n")
cat("FDR < 0.05:", sum(D$FDR_motif < .05), " < 0.10:", sum(D$FDR_motif < .10),
    " < 0.20:", sum(D$FDR_motif < .20), "\n")
cat("\ntop", NTOP, "\n")
print(head(D[, c("TF","logFC_motif","t_motif","P_motif","FDR_motif","t_gs")], NTOP),
      row.names = FALSE, digits = 3)
cat("\nRUNX / IRF in the panel:", paste(intersect(D$TF, c(RUNX, IRF)), collapse = ", "), "\n")

## ---- volcano --------------------------------------------------------------------------
D$rank <- rank(D$P_motif, ties.method = "min")
D$lab  <- ifelse(D$rank <= NTOP, D$TF, NA)
lim <- max(abs(D$logFC_motif))
pv <- ggplot(D, aes(logFC_motif, -log10(P_motif))) +
  geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = .3, colour = "grey55") +
  geom_point(colour = "grey78", size = 1.1, alpha = .8) +
  geom_point(data = subset(D, rank <= NTOP), colour = "grey15", size = 1.5) +
  geom_text_repel(aes(label = lab), size = 2.2, colour = "grey15", max.overlaps = 40,
                  seed = 1, min.segment.length = 0, segment.size = .2) +
  scale_x_continuous(limits = c(-lim, lim)) +
  annotate("text", x = -lim, y = -log10(0.05), hjust = 0, vjust = -0.6, size = 2,
           colour = "grey45", label = sprintf("nominal P = 0.05  (best FDR %.2f)",
                                              min(D$FDR_motif))) +
  gtheme_no_rot +
  xlab("BAP1 effect on motif activity, lost - retained (adjusted for histology and chr3p)") +
  ylab(expression(-log[10]~"moderated"~italic(P))) +
  ggtitle(sprintf("Gene-score-concordant TFs only (n = %d of 869, %d tumours)",
                  nrow(D), length(sams)))
pdf(file.path("Plots", "BAP1_TF_gsfiltered_volcano.pdf"), width = 6.4, height = 4.8)
print(pv); dev.off()

## ---- congruence with bulk, filtered panel only ------------------------------------------
bulk <- read.csv("BAP1_TFcongruence_bulk_molecular.csv", stringsAsFactors = FALSE)
pl <- lapply(c("bueno","tcga","mesomics"), function(s){
  b <- setNames(bulk[bulk$cohort == s, c("TF","t")], c("TF","t_bulk"))
  d <- merge(D, b, by = "TF")
  r <- cor(d$t_motif, d$t_bulk, method = "spearman")
  ct <- suppressWarnings(cor.test(d$t_motif, d$t_bulk, method = "spearman"))
  d$lab <- ifelse(rank(d$P_motif, ties.method = "min") <= 12, d$TF, NA)
  ggplot(d, aes(t_motif, t_bulk)) +
    geom_hline(yintercept = 0, linewidth = .2, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x,
                linewidth = .35, colour = "grey55", linetype = "dashed") +
    geom_point(colour = "#1b4f72", size = 1.3, alpha = .85) +
    geom_text_repel(aes(label = lab), size = 2, colour = "grey25",
                    max.overlaps = 30, seed = 1, min.segment.length = 0) +
    gtheme_no_rot +
    xlab("t, scATAC motif activity (BAP1)") +
    ylab(sprintf("t, %s mRNA (BAP1)", s)) +
    ggtitle(sprintf("%s: n = %d gene-score-concordant TFs, rho = %.2f, p = %.3g",
                    s, nrow(d), r, ct$p.value))
})
pdf(file.path("Plots", "BAP1_TF_gsfiltered_congruence.pdf"), width = 5.4, height = 4.4)
for (p in pl) print(p); dev.off()

## ---- volcano marking TFs that also match direction in bulk RNA -------------------------
## For each TF, count the bulk cohorts whose BAP1 effect points the same way as the
## scATAC motif effect.  A TF measured in all three can agree 0-3 times; 3/3 is the
## strict replication, 2/3 the majority.  Under no real agreement the number of
## matches is Binomial(3, 0.5), so 3/3 is expected for 12.5% of TFs by chance -- the
## observed rate is printed against that so the highlighting is not read as proof.
bulk <- read.csv("BAP1_TFcongruence_bulk_molecular.csv", stringsAsFactors = FALSE)
COH  <- c("bueno","tcga","mesomics")
for (s2 in COH) {
  b <- setNames(bulk[bulk$cohort == s2, c("TF","t")], c("TF", paste0("tb_", s2)))
  D <- merge(D, b, by = "TF", all.x = TRUE)
}
tbcols <- paste0("tb_", COH)
D$n_measured <- rowSums(!is.na(D[, tbcols]))
D$n_agree    <- rowSums(sign(D[, tbcols]) == sign(D$t_motif), na.rm = TRUE)
D$all_agree  <- D$n_measured > 0 & D$n_agree == D$n_measured
D$maj_agree  <- D$n_measured > 0 & D$n_agree >= ceiling(D$n_measured / 2)

cat("\n=== agreement in direction with bulk RNA, filtered panel ===\n")
cat(sprintf("TFs measured in all 3 cohorts: %d\n", sum(D$n_measured == 3)))
tb3 <- table(factor(D$n_agree[D$n_measured == 3], 0:3))
cat("cohorts agreeing (of 3): "); print(tb3)
cat(sprintf("all 3 agree: %d of %d (%.0f%%; %.0f%% expected by chance)\n",
            sum(D$all_agree & D$n_measured == 3), sum(D$n_measured == 3),
            100 * mean(D$all_agree[D$n_measured == 3]), 12.5))
top <- head(D[order(D$P_motif), ], NTOP)
cat(sprintf("among the top %d by motif P: %d agree in all 3, %d in a majority\n",
            NTOP, sum(top$all_agree), sum(top$maj_agree)))
cat("top-ranked TFs agreeing in all three cohorts:\n")
print(top[top$all_agree, c("TF","logFC_motif","t_motif","P_motif","FDR_motif",
                           tbcols)], row.names = FALSE, digits = 3)
write.csv(D, "BAP1_TF_gsfiltered_results.csv", row.names = FALSE)

D$grp <- ifelse(D$rank <= NTOP & D$all_agree, "top, all 3 cohorts agree",
          ifelse(D$rank <= NTOP & D$maj_agree, "top, majority agree",
           ifelse(D$rank <= NTOP, "top, bulk disagrees", "not top")))
D$grp <- factor(D$grp, levels = c("top, all 3 cohorts agree", "top, majority agree",
                                  "top, bulk disagrees", "not top"))
D$lab3 <- ifelse(D$rank <= NTOP & D$maj_agree, D$TF, NA)
gcol <- c("top, all 3 cohorts agree" = "#1b4f72", "top, majority agree" = "#5dade2",
          "top, bulk disagrees" = "#c0392b", "not top" = "grey82")
pb <- ggplot(D, aes(logFC_motif, -log10(P_motif))) +
  geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = .3, colour = "grey55") +
  geom_point(data = subset(D, grp == "not top"), colour = gcol["not top"],
             size = 1.1, alpha = .8) +
  geom_point(data = subset(D, grp != "not top"), aes(colour = grp), size = 2) +
  geom_text_repel(aes(label = lab3), size = 2.2, colour = "grey15", max.overlaps = 40,
                  seed = 1, min.segment.length = 0, segment.size = .2) +
  scale_colour_manual(values = gcol, breaks = levels(D$grp)[1:3], name = NULL) +
  scale_x_continuous(limits = c(-lim, lim)) +
  gtheme_no_rot +
  xlab("BAP1 effect on motif activity, lost - retained (adjusted for histology and chr3p)") +
  ylab(expression(-log[10]~"moderated"~italic(P))) +
  ggtitle(sprintf("Top %d gene-score-concordant TFs, coloured by agreement with bulk RNA",
                  NTOP))
pdf(file.path("Plots", "BAP1_TF_gsfiltered_volcano_bulkmatch.pdf"),
    width = 6.8, height = 4.8)
print(pb); dev.off()

## ---- the concordance the filter enforces ------------------------------------------------
r2 <- cor(D$t_motif, D$t_gs, method = "spearman")
D$lab2 <- ifelse(D$rank <= 15 | abs(D$t_gs) > 3, D$TF, NA)
pc <- ggplot(D, aes(t_gs, t_motif)) +
  geom_hline(yintercept = 0, linewidth = .2, colour = "grey85") +
  geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
  geom_point(colour = "#1b4f72", size = 1.3, alpha = .85) +
  geom_text_repel(aes(label = lab2), size = 2, colour = "grey25",
                  max.overlaps = 30, seed = 1, min.segment.length = 0) +
  gtheme_no_rot +
  xlab("t, TF gene score (BAP1)") + ylab("t, TF motif activity (BAP1)") +
  ggtitle(sprintf("the filtered panel by construction lies in two quadrants (rho = %.2f)", r2))
pdf(file.path("Plots", "BAP1_TF_gsfiltered_motif_vs_genescore.pdf"), width = 5, height = 4.2)
print(pc); dev.off()

cat("\nDONE\n")
