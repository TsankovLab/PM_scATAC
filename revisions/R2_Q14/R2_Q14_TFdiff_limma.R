###############################################################################
# R2_Q14 -- differential TF activity by BAP1 status, adjusted for histology,
#           using limma.
#
# Same question as R2_Q14_TFdiff_BAP1_sarcadjusted.R, standard machinery.
#
# The model, one line, fitted to all 869 motifs at once:
#
#     ~ BAP1_status + sarc_score        (11 tumours)
#
# limma fits that linear model per motif and reports the BAP1 coefficient with the
# sarcomatoid score held constant -- exactly the covariate adjustment asked for; a
# covariate in the design matrix is subtracted from the effect of interest.
#
# Why limma rather than a per-motif t test.  With 11 tumours and 2 predictors each
# motif has 8 residual degrees of freedom, so its variance estimate is unstable and
# a motif that happens to look consistent gets a huge t by accident.  limma's
# empirical Bayes step shrinks every motif's variance toward the trend across all
# 869, so each test effectively borrows degrees of freedom from the rest of the
# dataset.  That is the entire reason limma is used for small-n designs, and it is
# what makes the moderated t trustworthy where an ordinary t is not.  robust = TRUE
# stops a handful of extreme-variance motifs from dragging that shared prior.
#
# Multiple testing is Benjamini-Hochberg over the 869 motifs (limma's adj.P.Val).
#
# Three fits, so the effect of the adjustment is visible rather than asserted:
#   unadjusted   ~ BAP1_status
#   adjusted     ~ BAP1_status + sarc_score      <- the reported result
#   histology    the sarc_score coefficient of the adjusted fit
#
# BAP1 status is CLINICAL, taken from the per-sample genetic annotation of the
# cohort (hardcoded in R2_Q14_BAP1_TF_activity.R).  It is not inferred from these
# data, and it is deliberately not revised to match chr3p accessibility or the BAP1
# gene score -- those are measurements to be explained by the label, not inputs to
# it.  Molecular readouts disagree with the clinical call for some tumours (P3 reads
# as chr3p-retained, P11 as subclonal); that is a finding, not a labelling error.
#
# Tumours with fewer than MIN_CELLS malignant cells are excluded.  A tumour's TF
# activity here is a mean over its malignant cells, so a tumour resting on ~20 cells
# contributes an estimate an order of magnitude noisier than one resting on thousands
# -- yet in a per-tumour model every tumour counts once.  P3 (21 cells) and P13 (18)
# are the two affected; the next smallest is P14 at 131, so the threshold is not
# finely balanced.
#
# Input : TFactivity_deviation_per_sample.csv, BAP1_genescore_per_sample.csv,
#         scATAC_sarcscore_per_sample.csv (also supplies malignant cell counts)
# Output: BAP1_TF_limma_sarcadjusted.csv, Plots/BAP1_TF_limma_volcano.pdf,
#         Plots/BAP1_TF_limma_family.pdf
###############################################################################
suppressPackageStartupMessages({ library(limma); library(ggplot2); library(ggrepel) })

ROOT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUTDIR <- file.path(ROOT, "git_repo_claude", "R2_Q14")
setwd(OUTDIR); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))

RUNX <- c("RUNX1", "RUNX2", "RUNX3")
IRF  <- paste0("IRF", 1:9)
MIN_CELLS <- 50      # drops P3 (21) and P13 (18); next smallest is P14 at 131

## ---- data ---------------------------------------------------------------------
act <- read.csv("TFactivity_deviation_per_sample.csv", row.names = 1, check.names = FALSE)
gs  <- read.csv("BAP1_genescore_per_sample.csv", stringsAsFactors = FALSE)
sc  <- read.csv("scATAC_sarcscore_per_sample.csv", stringsAsFactors = FALSE)

sams <- intersect(colnames(act), intersect(gs$Sample, sc$sample))
ncell <- setNames(sc$n_cells, sc$sample)[sams]
small <- sams[ncell < MIN_CELLS]
cat("malignant cells per tumour:\n"); print(sort(ncell))
cat("excluded for < ", MIN_CELLS, " cells: ",
    if (length(small)) paste(sprintf("%s (%d)", small, ncell[small]), collapse = ", ")
    else "none", "\n", sep = "")
sams <- sams[ncell >= MIN_CELLS]
E    <- as.matrix(act[, sams, drop = FALSE])          # motifs x tumours, as limma wants
E    <- E[apply(E, 1, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE]
BAP1 <- factor(gs$BAP1_status[match(sams, gs$Sample)], levels = c("retained", "lost"))
SARC <- sc$sarc_score_atac[match(sams, sc$sample)]
stopifnot(!any(is.na(BAP1)), !any(is.na(SARC)))

cat("\nmotifs:", nrow(E), " tumours:", ncol(E), "\n")
cat("  BAP1-lost    :", paste(sams[BAP1 == "lost"], collapse = ", "), "\n")
cat("  BAP1-retained:", paste(sams[BAP1 == "retained"], collapse = ", "), "\n")
cat(sprintf("\nconfounding to be removed: mean sarc score %.3f (lost) vs %.3f (retained)\n",
            mean(SARC[BAP1 == "lost"]), mean(SARC[BAP1 == "retained"])))

## ---- the fits --------------------------------------------------------------------
des_adj <- model.matrix(~ BAP1 + SARC)      # coef 2 = BAP1lost, coef 3 = SARC
des_raw <- model.matrix(~ BAP1)
cat("\ndesign (adjusted):\n"); print(des_adj)

fit_adj <- eBayes(lmFit(E, des_adj), robust = TRUE)
fit_raw <- eBayes(lmFit(E, des_raw), robust = TRUE)
cat(sprintf("\nempirical Bayes: prior df = %.1f added to each motif's %d residual df\n",
            median(fit_adj$df.prior), fit_adj$df.residual[1]))

adj  <- topTable(fit_adj, coef = "BAP1lost", number = Inf, sort.by = "none")
raw  <- topTable(fit_raw, coef = "BAP1lost", number = Inf, sort.by = "none")
hist <- topTable(fit_adj, coef = "SARC",     number = Inf, sort.by = "none")

RES <- data.frame(
  TF = rownames(adj),
  diff_lost_minus_retained_adj = adj$logFC, t_adj = adj$t,
  P_adj = adj$P.Value, FDR_adj = adj$adj.P.Val,
  diff_unadjusted = raw$logFC[match(rownames(adj), rownames(raw))],
  t_unadjusted    = raw$t[match(rownames(adj), rownames(raw))],
  P_unadjusted    = raw$P.Value[match(rownames(adj), rownames(raw))],
  FDR_unadjusted  = raw$adj.P.Val[match(rownames(adj), rownames(raw))],
  beta_sarc = hist$logFC[match(rownames(adj), rownames(hist))],
  t_sarc    = hist$t[match(rownames(adj), rownames(hist))],
  P_sarc    = hist$P.Value[match(rownames(adj), rownames(hist))],
  FDR_sarc  = hist$adj.P.Val[match(rownames(adj), rownames(hist))],
  is_runx = rownames(adj) %in% RUNX, is_irf = rownames(adj) %in% IRF,
  stringsAsFactors = FALSE)
RES <- RES[order(RES$P_adj), ]
write.csv(RES, "BAP1_TF_limma_sarcadjusted.csv", row.names = FALSE)

## ---- what is significant -----------------------------------------------------------
cat("\n=== how many motifs pass ===\n")
for (a in c(0.05, 0.10, 0.20)) cat(sprintf(
  "FDR < %.2f : BAP1 adjusted %3d | BAP1 unadjusted %3d | histology (sarc) %3d\n",
  a, sum(RES$FDR_adj < a), sum(RES$FDR_unadjusted < a), sum(RES$FDR_sarc < a)))
cat(sprintf("nominal P < 0.05 : adjusted %d | unadjusted %d | histology %d\n",
            sum(RES$P_adj < .05), sum(RES$P_unadjusted < .05), sum(RES$P_sarc < .05)))

cat("\n=== top 20 by the histology-adjusted BAP1 effect ===\n")
print(head(RES[, c("TF","diff_lost_minus_retained_adj","t_adj","P_adj","FDR_adj",
                   "t_unadjusted","P_unadjusted","t_sarc")], 20),
      row.names = FALSE, digits = 3)

cat("\n=== RUNX ===\n")
print(RES[RES$is_runx, c("TF","diff_lost_minus_retained_adj","t_adj","P_adj","FDR_adj",
                         "t_unadjusted","P_unadjusted","t_sarc","P_sarc")],
      row.names = FALSE, digits = 3)
cat("\n=== IRF ===\n")
irf <- RES[RES$is_irf, c("TF","diff_lost_minus_retained_adj","t_adj","P_adj","FDR_adj",
                         "t_unadjusted","P_unadjusted","t_sarc")]
print(irf[order(irf$TF), ], row.names = FALSE, digits = 3)

## limma's own competitive gene-set test.  Unlike a Wilcoxon it corrects for the
## correlation between motifs, which is exactly the problem with treating 869 cisBP
## motifs as independent.  Both sets are passed in ONE call: cameraPR returns no FDR
## column when given a single set, and sprintf on that NULL prints nothing at all.
sets <- list(IRF = which(rownames(E) %in% IRF), RUNX = which(rownames(E) %in% RUNX))
cam  <- cameraPR(fit_adj$t[, "BAP1lost"], sets)
cat("\n=== cameraPR: is the family shifted relative to all other motifs? ===\n")
cat("    (on the histology-adjusted BAP1 moderated t)\n")
print(cam, digits = 3)

## ---- volcano -----------------------------------------------------------------------
## No family highlighting: colour encoded RUNX/IRF membership, which invited reading
## the colour as a result.  One colour for everything; the only annotation is the
## top NTOP motifs by moderated P, which is what the plot is for.
NTOP <- 20
RES$rank <- rank(RES$P_adj, ties.method = "min")
RES$lab  <- ifelse(RES$rank <= NTOP, RES$TF, NA)
lim <- max(abs(RES$diff_lost_minus_retained_adj))
sig <- min(RES$FDR_adj) < 0.05

pv <- ggplot(RES, aes(diff_lost_minus_retained_adj, -log10(P_adj))) +
  geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = .3, colour = "grey55") +
  geom_point(colour = "grey78", size = 1, alpha = .75) +
  geom_point(data = subset(RES, rank <= NTOP), colour = "grey25", size = 1.4) +
  geom_text_repel(aes(label = lab), size = 2.2, colour = "grey20", max.overlaps = 40,
                  seed = 1, segment.size = .2, min.segment.length = 0) +
  scale_x_continuous(limits = c(-lim, lim)) +
  annotate("text", x = -lim, y = -log10(0.05), hjust = 0, vjust = -0.6, size = 2,
           colour = "grey45",
           label = if (sig) "nominal P = 0.05" else
                   "nominal P = 0.05  (no motif reaches FDR 0.05)") +
  gtheme_no_rot +
  xlab("BAP1 effect, lost - retained (limma, histology-adjusted)") +
  ylab(expression(-log[10]~"moderated"~italic(P))) +
  ggtitle(sprintf("limma ~ BAP1 + sarcomatoid score (n = %d tumours, %d motifs)",
                  ncol(E), nrow(E)))
pdf(file.path("Plots", "BAP1_TF_limma_volcano.pdf"), width = 6.4, height = 4.8)
print(pv); dev.off()

cat(sprintf("\ntop %d motifs labelled on the volcano:\n", NTOP))
print(head(RES[, c("TF","diff_lost_minus_retained_adj","t_adj","P_adj","FDR_adj")], NTOP),
      row.names = FALSE, digits = 3)

## the family view lives in its own figure, so the volcano stays unannotated
dd <- rbind(data.frame(set = "all other motifs",
                       t = RES$t_adj[!RES$is_runx & !RES$is_irf]),
            data.frame(set = "IRF",  t = RES$t_adj[RES$is_irf]),
            data.frame(set = "RUNX", t = RES$t_adj[RES$is_runx]))
dd$set <- factor(dd$set, levels = c("all other motifs", "IRF", "RUNX"))
pf <- ggplot(dd, aes(set, t, fill = set)) +
  geom_hline(yintercept = 0, linewidth = .2, colour = "grey85") +
  geom_violin(scale = "width", linewidth = .25, colour = "grey40") +
  geom_jitter(width = .12, size = .8, alpha = .7) +
  scale_fill_manual(values = c(`all other motifs` = "grey85", IRF = "#1b7837",
                               RUNX = "#7d3c98"), guide = "none") +
  gtheme_no_rot + xlab(NULL) + ylab("moderated t (histology-adjusted)") +
  ggtitle("BAP1 effect by TF family")
pdf(file.path("Plots", "BAP1_TF_limma_family.pdf"), width = 4.6, height = 3.6)
print(pf); dev.off()

cat("\nDONE\n")
