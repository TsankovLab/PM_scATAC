###############################################################################
# R2_Q14 -- differential TF activity, BAP1-lost vs BAP1-retained,
#           with the sarcomatoid (cNMF20) score regressed out.
#
# Why the adjustment is not optional.  BAP1 status and histology are confounded in
# this cohort: the BAP1-retained tumours are the sarcomatoid ones.  Any TF that
# simply tracks the epithelioid-sarcomatoid axis will therefore look "BAP1
# associated" in an unadjusted test.  The model below asks the only question worth
# asking -- does the TF differ by BAP1 status ONCE the histology axis is accounted
# for -- and reports the unadjusted result alongside so the size of the confounding
# is visible per TF rather than asserted.
#
#     activity_TF(sample) ~ BAP1_status + sarc_score          (11 tumours)
#
# The unit is the TUMOUR, not the cell.  A per-cell test on a sample-level variable
# is pseudoreplication: with 20,635 malignant cells it would return vanishing
# p-values for an effect that rests on 11 independent observations.
#
# Significance is EXACT, not parametric.  With 11 tumours and 2 predictors there are
# 8 residual degrees of freedom, and a t distribution is not to be trusted there.
# There are only choose(11,5) = 462 ways to assign "lost" to 5 of 11 tumours, so
# every one is enumerated and the null distribution of the BAP1 t statistic is built
# exactly, holding the sarc score fixed.  The reported p is the fraction of those
# 462 assignments giving |t| at least as large as observed; the true labelling is
# one of them, so p is never below 1/462 = 0.00216.
#
# All 462 fits for all 869 motifs are solved in closed form (the design matrix is
# 11 x 3), so the exact test costs seconds.
#
# Input : TFactivity_deviation_per_sample.csv  (869 motifs x 11 tumours, mean
#           chromVAR deviation per tumour, from R2_Q14_BAP1_TF_activity.R)
#         BAP1_genescore_per_sample.csv
#         scATAC_sarcscore_per_sample.csv  (R2_Q14_sarcscore_scATAC.R, all 11)
# Output: BAP1_TF_diff_sarcadjusted.csv, BAP1_TF_diff_interferon_set.csv,
#         Plots/BAP1_TF_diff_sarcadjusted.pdf, Plots/BAP1_TF_family_distribution.pdf,
#         Plots/BAP1_TF_volcano_sarcadjusted.pdf
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })

ROOT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUTDIR <- file.path(ROOT, "git_repo_claude", "R2_Q14")
setwd(OUTDIR); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
RUNX <- c("RUNX1", "RUNX2", "RUNX3")
## IRF family (all 9 are in the cisBP set), plus the wider interferon/NF-kB axis
## they act in, kept as a separate set so the two questions stay separate.
IRF <- paste0("IRF", 1:9)
IFN <- c(IRF, paste0("STAT", c(1, 2, 3, 4, "5A", "5B", 6)),
         "NFKB1", "NFKB2", "REL", "RELA", "RELB", "PRDM1")

## ---- data --------------------------------------------------------------------
act <- read.csv("TFactivity_deviation_per_sample.csv", row.names = 1, check.names = FALSE)
gs  <- read.csv("BAP1_genescore_per_sample.csv", stringsAsFactors = FALSE)
sc  <- read.csv("scATAC_sarcscore_per_sample.csv", stringsAsFactors = FALSE)

sams <- intersect(colnames(act), intersect(gs$Sample, sc$sample))
Y    <- t(as.matrix(act[, sams, drop = FALSE]))          # tumours x TFs
bap1 <- as.integer(gs$BAP1_status[match(sams, gs$Sample)] == "lost")   # 1 = lost
sarc <- sc$sarc_score_atac[match(sams, sc$sample)]
stopifnot(!any(is.na(bap1)), !any(is.na(sarc)))
keep <- apply(Y, 2, function(x) all(is.finite(x)) && sd(x) > 0)
Y    <- Y[, keep, drop = FALSE]
n    <- length(sams); nl <- sum(bap1)

cat("tumours:", n, "|", nl, "BAP1-lost,", n - nl, "retained\n")
cat("  lost    :", paste(sams[bap1 == 1], collapse = ", "), "\n")
cat("  retained:", paste(sams[bap1 == 0], collapse = ", "), "\n")
cat("motifs tested:", ncol(Y), "\n")

## how badly are BAP1 status and histology confounded?
cf <- suppressWarnings(cor.test(bap1, sarc, method = "spearman"))
cat(sprintf("\nconfounding: mean sarc score %.3f in BAP1-lost vs %.3f in retained (rho = %.3f, p = %.3g)\n",
            mean(sarc[bap1 == 1]), mean(sarc[bap1 == 0]), cf$estimate, cf$p.value))

## ---- closed-form OLS for every motif at once -----------------------------------
## Returns the t statistic of the BAP1 term (column 2 of X) for all motifs.
tstat <- function(b, adjust = TRUE){
  X  <- if (adjust) cbind(1, b, sarc) else cbind(1, b)
  XX <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XX)) return(rep(NA_real_, ncol(Y)))
  B  <- XX %*% crossprod(X, Y)                    # coef x motifs
  R  <- Y - X %*% B                               # residuals
  s2 <- colSums(R^2) / (nrow(X) - ncol(X))
  list(beta = B[2, ], t = B[2, ] / sqrt(s2 * XX[2, 2]))
}
obs_adj <- tstat(bap1, TRUE)
obs_raw <- tstat(bap1, FALSE)

## the sarc term itself, so the histology-driven motifs can be named
Xs <- cbind(1, bap1, sarc); XXs <- solve(crossprod(Xs))
Bs <- XXs %*% crossprod(Xs, Y)
s2s <- colSums((Y - Xs %*% Bs)^2) / (n - 3)
t_sarc <- Bs[3, ] / sqrt(s2s * XXs[3, 3])
p_sarc <- 2 * pt(-abs(t_sarc), df = n - 3)

## ---- exact permutation null ------------------------------------------------------
## Every way of calling 5 of the 11 tumours "lost", the true one included.
combos <- utils::combn(n, nl)
cat("\nexact null: enumerating all", ncol(combos), "assignments of", nl, "lost of", n, "\n")
null_adj <- matrix(NA_real_, ncol(combos), ncol(Y))
null_raw <- matrix(NA_real_, ncol(combos), ncol(Y))
for (i in seq_len(ncol(combos))){
  b <- integer(n); b[combos[, i]] <- 1L
  null_adj[i, ] <- tstat(b, TRUE)$t
  null_raw[i, ] <- tstat(b, FALSE)$t
}
pperm <- function(nullmat, obs)
  colMeans(abs(nullmat) >= matrix(abs(obs), nrow(nullmat), length(obs), byrow = TRUE),
           na.rm = TRUE)
p_adj <- pperm(null_adj, obs_adj$t)
p_raw <- pperm(null_raw, obs_raw$t)

RES <- data.frame(
  TF = colnames(Y),
  delta_lost_minus_retained = colMeans(Y[bap1 == 1, , drop = FALSE]) -
                              colMeans(Y[bap1 == 0, , drop = FALSE]),
  beta_unadjusted = obs_raw$beta, t_unadjusted = obs_raw$t, p_unadjusted = p_raw,
  beta_sarcadjusted = obs_adj$beta, t_sarcadjusted = obs_adj$t, p_sarcadjusted = p_adj,
  beta_sarc = Bs[3, ], t_sarc = t_sarc, p_sarc = p_sarc,
  is_runx = colnames(Y) %in% RUNX, is_irf = colnames(Y) %in% IRF,
  is_ifn  = colnames(Y) %in% IFN, stringsAsFactors = FALSE)
RES$fdr_unadjusted   <- p.adjust(RES$p_unadjusted,   "BH")
RES$fdr_sarcadjusted <- p.adjust(RES$p_sarcadjusted, "BH")
RES$fdr_sarc         <- p.adjust(RES$p_sarc,         "BH")
RES <- RES[order(RES$p_sarcadjusted, -abs(RES$t_sarcadjusted)), ]
write.csv(RES, "BAP1_TF_diff_sarcadjusted.csv", row.names = FALSE)

## ---- what survives ------------------------------------------------------------------
cat(sprintf("\nsmallest attainable p with %d permutations: %.4f\n", ncol(combos), 1/ncol(combos)))
for (a in c(0.01, 0.05)) cat(sprintf("raw p < %.2f : unadjusted %3d motifs | sarc-adjusted %3d\n",
    a, sum(RES$p_unadjusted < a), sum(RES$p_sarcadjusted < a)))
cat(sprintf("BH FDR < 0.10: unadjusted %d | sarc-adjusted %d | (sarc term itself: %d)\n",
    sum(RES$fdr_unadjusted < .1), sum(RES$fdr_sarcadjusted < .1), sum(RES$fdr_sarc < .1)))

cat("\n=== top 20 motifs by the sarc-adjusted BAP1 term ===\n")
print(head(RES[, c("TF","delta_lost_minus_retained","t_unadjusted","p_unadjusted",
                   "t_sarcadjusted","p_sarcadjusted","fdr_sarcadjusted","t_sarc")], 20),
      row.names = FALSE, digits = 3)

cat("\n=== RUNX ===\n")
print(RES[RES$is_runx, c("TF","delta_lost_minus_retained","t_unadjusted","p_unadjusted",
                         "t_sarcadjusted","p_sarcadjusted","t_sarc","p_sarc")],
      row.names = FALSE, digits = 3)

## how much does adjusting change things?
cat(sprintf("\ncorrelation of unadjusted and adjusted t: Pearson %.3f\n",
            cor(RES$t_unadjusted, RES$t_sarcadjusted)))
lost <- RES$TF[RES$p_unadjusted < 0.05 & RES$p_sarcadjusted >= 0.05]
gain <- RES$TF[RES$p_unadjusted >= 0.05 & RES$p_sarcadjusted < 0.05]
cat("motifs significant unadjusted but NOT after adjustment (histology-driven):",
    length(lost), "\n  ", paste(head(lost, 25), collapse = ", "), "\n")
cat("motifs significant only AFTER adjustment (histology was masking them):",
    length(gain), "\n  ", paste(head(gain, 25), collapse = ", "), "\n")

## ---- IRF / interferon distribution --------------------------------------------
## Two questions, kept apart: (1) where do the individual IRFs sit, (2) is the family
## as a whole shifted relative to the other 860 motifs?  A family-level shift is
## tested with a two-sided Wilcoxon of the t statistics against all other motifs --
## the motifs are not independent, so this is descriptive, not a formal family test.
for (setn in list(c("IRF", "is_irf"), c("interferon/NF-kB axis", "is_ifn"))) {
  ix <- RES[[setn[2]]]
  cat(sprintf("\n=== %s: %d motifs of %d ===\n", setn[1], sum(ix), nrow(RES)))
  for (col in c("t_unadjusted", "t_sarcadjusted", "t_sarc")) {
    w <- suppressWarnings(wilcox.test(RES[[col]][ix], RES[[col]][!ix]))
    cat(sprintf("  %-16s median %6.2f (set) vs %6.2f (rest) | Wilcoxon p = %.3g\n",
                col, median(RES[[col]][ix]), median(RES[[col]][!ix]), w$p.value))
  }
  cat(sprintf("  nominal p < 0.05: unadjusted %d/%d, sarc-adjusted %d/%d (expected ~%.1f by chance)\n",
              sum(ix & RES$p_unadjusted < .05), sum(ix),
              sum(ix & RES$p_sarcadjusted < .05), sum(ix), 0.05 * sum(ix)))
}
cat("\n=== the 9 IRF motifs, individually ===\n")
irf_tab <- RES[RES$is_irf, c("TF","delta_lost_minus_retained","t_unadjusted","p_unadjusted",
                             "t_sarcadjusted","p_sarcadjusted","t_sarc","p_sarc")]
print(irf_tab[order(irf_tab$TF), ], row.names = FALSE, digits = 3)
write.csv(RES[RES$is_ifn, ], "BAP1_TF_diff_interferon_set.csv", row.names = FALSE)

## ---- figure ----------------------------------------------------------------------
## Points are NOT coloured by whether adjustment changed their significance -- that
## encoding invited reading the colour as a result.  Everything is grey; only the two
## families being asked about are marked.
RES$fam <- ifelse(RES$is_runx, "RUNX", ifelse(RES$is_irf, "IRF", "other"))
RES$lab <- ifelse(RES$fam != "other", RES$TF, NA)
famcol  <- c(RUNX = "#7d3c98", IRF = "#1b7837", other = "grey80")

base <- function(x, y, xlab, ylab, ttl)
  ggplot(RES, aes(.data[[x]], .data[[y]])) +
    geom_hline(yintercept = 0, linewidth = .2, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
    geom_point(data = subset(RES, fam == "other"), colour = famcol["other"],
               size = 1, alpha = .7) +
    geom_point(data = subset(RES, fam != "other"), aes(colour = fam), size = 1.7) +
    geom_text_repel(data = subset(RES, fam != "other"), aes(label = lab, colour = fam),
                    size = 2.2, max.overlaps = 40, seed = 1, show.legend = FALSE) +
    scale_colour_manual(values = famcol, breaks = c("RUNX","IRF"), name = NULL) +
    gtheme_no_rot + xlab(xlab) + ylab(ylab) + ggtitle(ttl)

p1 <- base("t_unadjusted", "t_sarcadjusted",
           "t, BAP1 lost vs retained (unadjusted)",
           "t, BAP1 lost vs retained (sarc-adjusted)",
           sprintf("TF activity by BAP1 status, before and after adjusting for histology (n = %d tumours)", n)) +
  geom_abline(slope = 1, linetype = 2, linewidth = .3, colour = "grey60")

## the second view: is a motif about BAP1, or about histology?
p2 <- base("t_sarc", "t_sarcadjusted",
           "t, sarcomatoid score (histology axis)",
           "t, BAP1 lost vs retained (sarc-adjusted)",
           "BAP1 effect vs histology effect, same model")

pdf(file.path("Plots", "BAP1_TF_diff_sarcadjusted.pdf"), width = 7, height = 5)
print(p1); print(p2); dev.off()

## ---- volcano, histology-adjusted only ------------------------------------------
## x is the adjusted BAP1 coefficient (difference in mean chromVAR deviation,
## lost - retained, holding the cNMF20 score fixed); y is the EXACT permutation p.
## That p cannot go below 1/462 = 0.00216, so the top of the plot is a hard ceiling
## and not a cluster of especially strong hits -- the ceiling is drawn explicitly.
## No motif reaches BH FDR < 0.10, so the only guide line is nominal p = 0.05 and the
## plot is a description of the ranking, not a set of discoveries.
PFLOOR <- 1 / ncol(combos)
V <- RES
V$logp <- -log10(pmax(V$p_sarcadjusted, PFLOOR))
V$top  <- rank(V$p_sarcadjusted, ties.method = "min") <= 15 &
          abs(V$beta_sarcadjusted) >= quantile(abs(V$beta_sarcadjusted), .90)
V$lab  <- ifelse(V$fam != "other" | V$top, V$TF, NA)
lim <- max(abs(V$beta_sarcadjusted))

p4 <- ggplot(V, aes(beta_sarcadjusted, logp)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = .3, colour = "grey55") +
  geom_hline(yintercept = -log10(PFLOOR), linetype = 3, linewidth = .3, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = .2, colour = "grey85") +
  geom_point(data = subset(V, fam == "other"), colour = "grey80", size = 1, alpha = .7) +
  geom_point(data = subset(V, fam != "other"), aes(colour = fam), size = 1.7) +
  geom_text_repel(aes(label = lab, colour = fam), size = 2.1, max.overlaps = 40,
                  seed = 1, show.legend = FALSE, segment.size = .2) +
  scale_colour_manual(values = famcol, breaks = c("RUNX", "IRF"), name = NULL) +
  scale_x_continuous(limits = c(-lim, lim)) +
  annotate("text", x = -lim, y = -log10(PFLOOR), hjust = 0, vjust = -0.5, size = 2,
           colour = "grey45", label = sprintf("permutation floor p = %.4f", PFLOOR)) +
  annotate("text", x = -lim, y = -log10(0.05), hjust = 0, vjust = -0.5, size = 2,
           colour = "grey45", label = "p = 0.05 (nominal; no motif reaches FDR 0.10)") +
  gtheme_no_rot +
  xlab("BAP1 effect, lost - retained (histology-adjusted)") +
  ylab(expression(-log[10]~"exact permutation p")) +
  ggtitle(sprintf("Histology-adjusted differential TF activity by BAP1 status (n = %d tumours)", n))
pdf(file.path("Plots", "BAP1_TF_volcano_sarcadjusted.pdf"), width = 6.4, height = 4.8)
print(p4); dev.off()
cat("\nvolcano: higher in BAP1-lost", sum(V$beta_sarcadjusted > 0 & V$p_sarcadjusted < .05),
    "| higher in retained", sum(V$beta_sarcadjusted < 0 & V$p_sarcadjusted < .05),
    "at nominal p < 0.05\n")

## distribution of the family t statistics against the rest
dd <- rbind(
  data.frame(set = "all other motifs", t = RES$t_sarcadjusted[RES$fam == "other"]),
  data.frame(set = "IRF",  t = RES$t_sarcadjusted[RES$is_irf]),
  data.frame(set = "RUNX", t = RES$t_sarcadjusted[RES$is_runx]))
dd$set <- factor(dd$set, levels = c("all other motifs", "IRF", "RUNX"))
p3 <- ggplot(dd, aes(set, t, fill = set)) +
  geom_hline(yintercept = 0, linewidth = .2, colour = "grey85") +
  geom_violin(scale = "width", linewidth = .25, colour = "grey40") +
  geom_jitter(width = .12, size = .8, alpha = .7) +
  scale_fill_manual(values = c(`all other motifs` = "grey85", IRF = "#1b7837",
                               RUNX = "#7d3c98"), guide = "none") +
  gtheme_no_rot + xlab(NULL) +
  ylab("t, BAP1 lost vs retained (sarc-adjusted)") +
  ggtitle("distribution of the BAP1 effect by TF family")
pdf(file.path("Plots", "BAP1_TF_family_distribution.pdf"), width = 4.6, height = 3.6)
print(p3); dev.off()

cat("\nDONE\n")
