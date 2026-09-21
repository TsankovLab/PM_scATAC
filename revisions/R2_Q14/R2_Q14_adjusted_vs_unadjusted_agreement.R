###############################################################################
# R2_Q14 -- does histology adjustment of the scATAC BAP1 effect improve agreement with
#           TCGA, or worsen it?
#
# scATAC (malignant cells, 9 tumours; MIN_CELLS 50 drops P3, P13), limma per TF:
#   unadjusted  ~ BAP1          adjusted  ~ BAP1 + sarcomatoid score (cNMF20)
#   readouts    motif activity (chromVAR deviation), TF gene score, and the
#               gene-score-filtered motif panel (motifs with a detectably accessible TF
#               locus whose motif and gene-score effects share a sign); for the panel the
#               comparison is restricted to TFs passing the filter under BOTH models so
#               the two correlations are computed on the same TFs.
# TCGA targets:
#   TCGA mRNA (genetic label)      ALTERED vs WILD-TYPE, unadjusted t
#   TCGA mRNA (genetic, adjusted)  same, ~ BAP1 + sarc score
#   paper tab 2B                   affinity-regression TF activity estimate (unadjusted)
#   AR reproduction                our reproduction of tab 2B (unadjusted t)
# Statistics per readout x target, on the TFs shared by both scATAC fits and the target:
#   rho_unadj, rho_adj           Spearman(scATAC t, TCGA statistic)
#   delta = rho_adj - rho_unadj  paired bootstrap over TFs (5,000) -> 95% CI, p
#                                Meng-Rosenthal-Rubin z-test for correlated correlations
#                                (the two scATAC t vectors are themselves correlated)
#   perm_p                       exact label-permutation p for each rho (126 assignments)
# Sensitivity: TFs on chr3p removed (positional artefact of BAP1 deletion in TCGA mRNA).
# Caveat: motifs are not independent observations, so the bootstrap and z-test p-values
# are optimistic; the permutation p is exact for the label but not for TF dependence.
#
# Output: adjusted_vs_unadjusted_agreement.csv, Plots/adjusted_vs_unadjusted_agreement.pdf
###############################################################################
suppressPackageStartupMessages({ library(limma); library(ggplot2) })
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R2_Q14")); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
set.seed(1); MIN_CELLS <- 50; NBOOT <- 5000

gs  <- read.csv("BAP1_genescore_per_sample.csv", stringsAsFactors = FALSE)
sc  <- read.csv("scATAC_sarcscore_per_sample.csv", stringsAsFactors = FALSE)
ACT <- as.matrix(read.csv("TFactivity_deviation_per_sample.csv", row.names = 1, check.names = FALSE))
GSC <- as.matrix(read.csv("TFgenescore_per_sample.csv", row.names = 1, check.names = FALSE))
sams <- Reduce(intersect, list(colnames(ACT), colnames(GSC), gs$Sample, sc$sample))
sams <- sams[setNames(sc$n_cells, sc$sample)[sams] >= MIN_CELLS]
LOST <- gs$BAP1_status[match(sams, gs$Sample)] == "lost"
SARC <- sc$sarc_score_atac[match(sams, sc$sample)]
prep <- function(M) { M <- M[, sams, drop = FALSE]; M[apply(M, 1, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE] }
ACT <- prep(ACT); GSC <- prep(GSC)
detected <- rownames(GSC)[rowMeans(GSC) >= median(rowMeans(GSC))]
cat("scATAC tumours:", length(sams), "| lost:", paste(sams[LOST], collapse = " "), "\n")
cat(sprintf("cor(sarc score, BAP1 lost) in these tumours: %.2f\n", cor(SARC, LOST)))

fit_t <- function(M, lost, adjust) {
  grp <- factor(ifelse(lost, "lost", "retained"), levels = c("retained", "lost"))
  des <- if (adjust) model.matrix(~ grp + SARC) else model.matrix(~ grp)
  tt <- topTable(eBayes(lmFit(M, des), robust = TRUE), coef = "grplost", number = Inf, sort.by = "none")
  setNames(tt$t, rownames(tt))
}
scatac <- function(lost) {
  m_u <- fit_t(ACT, lost, FALSE); m_a <- fit_t(ACT, lost, TRUE)
  g_u <- fit_t(GSC, lost, FALSE); g_a <- fit_t(GSC, lost, TRUE)
  both <- intersect(names(m_u), names(g_u))
  keep <- both[both %in% detected & sign(m_u[both]) == sign(g_u[both]) & sign(m_a[both]) == sign(g_a[both])]
  list(motif = list(u = m_u, a = m_a), genescore = list(u = g_u, a = g_a),
       gs_filtered_motif = list(u = m_u[keep], a = m_a[keep]))
}

## ---- TCGA targets ------------------------------------------------------------------
tm <- read.csv("TCGA_genetic_BAP1_TFonly_unadjusted.csv", stringsAsFactors = FALSE)
pb <- read.csv("affinity_regression/prepared/paper_2B.csv", stringsAsFactors = FALSE)
ar <- read.csv("affinity_regression/results/primary_cv/bap1_tf_ttest.csv", stringsAsFactors = FALSE)
TARGETS <- list(
  "TCGA mRNA, genetic label"             = setNames(tm$t, tm$TF),
  "TCGA mRNA, genetic label, sarc-adj"   = setNames(tm$t_sarcadjusted, tm$TF),
  "paper tab 2B (TF activity)"           = setNames(pb$estimate_tf, pb$TF),
  "affinity-regression reproduction"     = setNames(ar$t, ar$TF))
CHR3P <- tm$TF[tm$chr3p %in% TRUE]

mrr <- function(r1, r2, r12, n) {                       # Meng, Rosenthal & Rubin 1992
  z1 <- atanh(r1); z2 <- atanh(r2); rb2 <- (r1^2 + r2^2) / 2
  f <- min(1, (1 - r12) / (2 * (1 - rb2))); h <- (1 - f * rb2) / (1 - rb2)
  z <- (z2 - z1) * sqrt((n - 3) / (2 * (1 - r12) * h)); c(z = z, p = 2 * pnorm(-abs(z)))
}
rhos <- function(S, target, exclude = character(0)) {
  tf <- setdiff(Reduce(intersect, list(names(S$u), names(S$a), names(target))), exclude)
  tf <- tf[is.finite(target[tf])]
  c(n = length(tf), u = cor(S$u[tf], target[tf], method = "spearman"),
    a = cor(S$a[tf], target[tf], method = "spearman"))
}

OBS <- scatac(LOST)
combos <- combn(length(sams), sum(LOST))
cat("exact permutation null:", ncol(combos), "label assignments\n")
PERM <- lapply(seq_len(ncol(combos)), function(i) scatac(seq_along(sams) %in% combos[, i]))

RES <- list()
for (rd in names(OBS)) for (tg in names(TARGETS)) for (ex in c("all TFs", "chr3p TFs removed")) {
  exclude <- if (ex == "all TFs") character(0) else CHR3P
  S <- OBS[[rd]]; target <- TARGETS[[tg]]
  tf <- setdiff(Reduce(intersect, list(names(S$u), names(S$a), names(target))), exclude); tf <- tf[is.finite(target[tf])]
  u <- S$u[tf]; a <- S$a[tf]; y <- target[tf]
  ru <- cor(u, y, method = "spearman"); ra <- cor(a, y, method = "spearman"); rua <- cor(u, a, method = "spearman")
  bt <- replicate(NBOOT, { i <- sample.int(length(tf), replace = TRUE)
                           cor(a[i], y[i], method = "spearman") - cor(u[i], y[i], method = "spearman") })
  z <- mrr(ru, ra, rua, length(tf))
  pr <- t(sapply(PERM, function(P) rhos(P[[rd]], target, exclude)[c("u", "a")]))
  RES[[length(RES) + 1]] <- data.frame(scATAC_readout = rd, TCGA_target = tg, tf_set = ex, n_TF = length(tf),
    rho_unadjusted = ru, rho_adjusted = ra, delta_adj_minus_unadj = ra - ru,
    delta_CI_low = quantile(bt, .025), delta_CI_high = quantile(bt, .975),
    delta_boot_p = 2 * min(mean(bt <= 0), mean(bt >= 0)), cor_unadj_vs_adj_t = rua,
    MRR_z = z["z"], MRR_p = z["p"],
    perm_p_unadjusted = mean(pr[, "u"] >= ru), perm_p_adjusted = mean(pr[, "a"] >= ra), row.names = NULL)
}
RES <- do.call(rbind, RES)
write.csv(RES, "adjusted_vs_unadjusted_agreement.csv", row.names = FALSE)
fmt <- function(d) {                                     # round whichever columns are present
  dig <- c(rho_unadjusted = 3, rho_adjusted = 3, delta_adj_minus_unadj = 3, delta_CI_low = 3, delta_CI_high = 3,
           cor_unadj_vs_adj_t = 2, MRR_z = 2, perm_p_unadjusted = 3, perm_p_adjusted = 3)
  for (k in intersect(names(dig), names(d))) d[[k]] <- round(d[[k]], dig[[k]])
  for (k in intersect(c("delta_boot_p", "MRR_p"), names(d))) d[[k]] <- signif(d[[k]], 2)
  d
}
cat("\n=== scATAC unadjusted vs histology-adjusted: agreement with TCGA (all TFs) ===\n")
print(fmt(RES[RES$tf_set == "all TFs", setdiff(names(RES), "tf_set")]), row.names = FALSE)
cat("\n=== sensitivity: chr3p TFs removed ===\n")
print(fmt(RES[RES$tf_set != "all TFs", c("scATAC_readout", "TCGA_target", "n_TF", "rho_unadjusted", "rho_adjusted",
                                         "delta_adj_minus_unadj", "delta_CI_low", "delta_CI_high", "delta_boot_p",
                                         "perm_p_unadjusted", "perm_p_adjusted")]), row.names = FALSE)
cat(sprintf("\nadjusted better in %d of %d comparisons (all TFs); significant (bootstrap p<0.05): %d better, %d worse\n",
    sum(RES$delta_adj_minus_unadj[RES$tf_set == "all TFs"] > 0), sum(RES$tf_set == "all TFs"),
    sum(RES$tf_set == "all TFs" & RES$delta_boot_p < .05 & RES$delta_adj_minus_unadj > 0),
    sum(RES$tf_set == "all TFs" & RES$delta_boot_p < .05 & RES$delta_adj_minus_unadj < 0)))

D <- RES[RES$tf_set == "all TFs", ]
D$TCGA_target <- factor(D$TCGA_target, levels = rev(names(TARGETS)))
L <- rbind(data.frame(D[, c("scATAC_readout", "TCGA_target")], model = "~BAP1", rho = D$rho_unadjusted, p = D$perm_p_unadjusted),
           data.frame(D[, c("scATAC_readout", "TCGA_target")], model = "~BAP1 + sarc", rho = D$rho_adjusted, p = D$perm_p_adjusted))
p <- ggplot(L, aes(rho, TCGA_target)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = .3) +
  geom_line(aes(group = TCGA_target), colour = "grey70") +
  geom_point(aes(colour = model, shape = p < 0.05), size = 2.8) +
  scale_colour_manual(values = c(`~BAP1` = "#e08214", `~BAP1 + sarc` = "#2c7fb8"), name = "scATAC model") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "label-perm p < 0.05") +
  geom_text(data = D, aes(x = pmax(rho_unadjusted, rho_adjusted) + 0.03, y = TCGA_target,
                          label = sprintf("delta %+.2f [%+.2f, %+.2f]", delta_adj_minus_unadj, delta_CI_low, delta_CI_high)),
            hjust = 0, size = 2.3, colour = "grey30") +
  facet_wrap(~ scATAC_readout, ncol = 1) + scale_x_continuous(expand = expansion(mult = c(.05, .45))) +
  gtheme_no_rot + xlab("Spearman rho, scATAC BAP1 t vs TCGA statistic") + ylab(NULL) +
  ggtitle("Histology adjustment of the scATAC BAP1 effect: agreement with TCGA",
          subtitle = "delta = adjusted - unadjusted, paired bootstrap 95% CI over TFs") +
  theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
ggsave("Plots/adjusted_vs_unadjusted_agreement.pdf", p, width = 8.5, height = 8)
cat("\nDONE\n")
