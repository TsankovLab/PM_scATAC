###############################################################################
# R2_Q14 -- how well do our scATAC readouts agree with the TCGA MESO paper's BAP1
#           inferred-TF-activity list (Hmeljak et al. 2018, supp tab 2B)?
#
# Paper list: t-test of affinity-regression TF activity, BAP1 inactivated vs wild-type
#             (estimate > 0 = higher in inactivated; confirmed by our reproduction).
# Our readouts, scATAC malignant cells, 9 tumours (MIN_CELLS 50 drops P3, P13):
#   motif activity   per-tumour mean chromVAR deviation   (TFactivity_deviation_per_sample.csv)
#   gene score       per-tumour log2 TF gene score        (TFgenescore_per_sample.csv)
#   gene-score-filtered motif panel: motifs whose BAP1 effect has the same sign in gene
#                    score and whose TF locus is detectably accessible (>= panel median)
# each fitted with limma  ~ BAP1                (as the paper: no histology term)
#                    and  ~ BAP1 + sarc score   (our histology-adjusted analysis)
# BAP1 label: our genetic annotation (lost P4 P5 P8 P10 vs retained P1 P11 P12 P14 P23).
#
# Agreement = Spearman between our moderated t and the paper's estimate, sign agreement,
# and recovery of the paper's 28 FDR<0.01 TFs.
# Calibration: EXACT label-permutation null -- all choose(9,4) = 126 assignments of
# "lost" are refitted, so the question "is this agreement more than chance at n = 9?" has
# an exact answer.
# References: our affinity-regression reproduction of 2B on TCGA (results/primary_cv) and
# TCGA TF mRNA with the genetic BAP1 label.
# Sensitivity: the paper names six TFs via TRANSFAC alias quirks that do not correspond to
# the gene symbol (POSTN, ACTR1A, ADD1, SF1, NF1, TTF1); agreement is repeated without them.
#
# Output: paper2B_vs_scATAC_agreement.csv, paper2B_vs_scATAC_TFlevel.csv,
#         Plots/paper2B_vs_scATAC_agreement.pdf
###############################################################################
suppressPackageStartupMessages({ library(limma); library(ggplot2); library(ggrepel); library(patchwork) })
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R2_Q14")); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
MIN_CELLS <- 50
MISMAP <- c("POSTN", "ACTR1A", "ADD1", "SF1", "NF1", "TTF1")

paper <- read.csv("affinity_regression/prepared/paper_2B.csv", stringsAsFactors = FALSE)
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
gs_detected <- rownames(GSC)[rowMeans(GSC) >= median(rowMeans(GSC))]
cat("tumours:", length(sams), "| lost:", paste(sams[LOST], collapse = " "), "| retained:", paste(sams[!LOST], collapse = " "), "\n")
cat("motif TFs:", nrow(ACT), "| gene-score TFs:", nrow(GSC), "| paper TFs:", nrow(paper),
    "| paper in motif:", sum(paper$TF %in% rownames(ACT)), "| paper in gene score:", sum(paper$TF %in% rownames(GSC)), "\n")

fit_t <- function(M, lost, adjust) {
  grp <- factor(ifelse(lost, "lost", "retained"), levels = c("retained", "lost"))
  des <- if (adjust) model.matrix(~ grp + SARC) else model.matrix(~ grp)
  tt <- topTable(eBayes(lmFit(M, des), robust = TRUE), coef = "grplost", number = Inf, sort.by = "none")
  data.frame(TF = rownames(tt), t = tt$t, P = tt$P.Value, stringsAsFactors = FALSE)
}
readouts <- function(lost, adjust) {
  m <- fit_t(ACT, lost, adjust); g <- fit_t(GSC, lost, adjust)
  f <- merge(m, g[, c("TF", "t")], by = "TF", suffixes = c("", "_gs"))
  f <- f[f$TF %in% gs_detected & sign(f$t) == sign(f$t_gs), c("TF", "t", "P")]
  list(motif = m, genescore = g, gs_filtered_motif = f)
}
agree <- function(d, exclude = character(0)) {
  m <- merge(d, paper, by = "TF"); m <- m[!m$TF %in% exclude, ]
  if (nrow(m) < 5) return(NULL)
  s <- m$p_adj_tf < 0.01; same <- sign(m$t) == sign(m$estimate_tf)
  ct <- suppressWarnings(cor.test(m$t, m$estimate_tf, method = "spearman"))
  data.frame(n_shared = nrow(m), rho = unname(ct$estimate), rho_p = ct$p.value,
             rho_signed_logp = cor(sign(m$t) * -log10(m$P), sign(m$estimate_tf) * -log10(m$p_tf), method = "spearman"),
             pct_same_sign = 100 * mean(same), n_paper_sig = sum(s),
             pct_same_sign_paper_sig = 100 * mean(same[s]),
             n_paper_sig_same_sign_P05 = sum(s & same & m$P < 0.05))
}

## ---- observed ----------------------------------------------------------------------
OBS <- list(); TFL <- list()
for (adj in c(FALSE, TRUE)) {
  r <- readouts(LOST, adj); lab <- if (adj) "adjusted (~BAP1 + sarc)" else "unadjusted (~BAP1)"
  for (k in names(r)) {
    for (ex in list(character(0), MISMAP)) {
      a <- agree(r[[k]], ex)
      OBS[[length(OBS) + 1]] <- cbind(data.frame(readout = k, model = lab,
                                                 tf_set = if (length(ex)) "excluding alias-quirk TFs" else "all shared"), a)
    }
    TFL[[paste(k, adj)]] <- setNames(r[[k]][, c("TF", "t", "P")], c("TF", paste0(k, if (adj) "_adj" else "_unadj", c("_t", "_P"))))
  }
}
ref_ar   <- read.csv("affinity_regression/results/primary_cv/bap1_tf_ttest.csv", stringsAsFactors = FALSE)
ref_ar   <- data.frame(TF = ref_ar$TF, t = ref_ar$t, P = ref_ar$p)
ref_tcga <- read.csv("TCGA_genetic_BAP1_TFonly_unadjusted.csv", stringsAsFactors = FALSE)[, c("TF", "t", "P")]
for (x in list(list("affinity-regression reproduction (TCGA)", ref_ar), list("TCGA TF mRNA, genetic BAP1", ref_tcga)))
  for (ex in list(character(0), MISMAP))
    OBS[[length(OBS) + 1]] <- cbind(data.frame(readout = x[[1]], model = "reference",
                                               tf_set = if (length(ex)) "excluding alias-quirk TFs" else "all shared"), agree(x[[2]], ex))
OBS <- do.call(rbind, OBS)

## ---- exact label-permutation null (126 assignments) -----------------------------------
combos <- combn(length(sams), sum(LOST)); cat("\nexact null:", ncol(combos), "label assignments\n")
NULL_R <- list()
for (i in seq_len(ncol(combos))) {
  lost <- seq_along(sams) %in% combos[, i]
  for (adj in c(FALSE, TRUE)) {
    r <- readouts(lost, adj); lab <- if (adj) "adjusted (~BAP1 + sarc)" else "unadjusted (~BAP1)"
    for (k in names(r)) { a <- agree(r[[k]])
      if (!is.null(a)) NULL_R[[length(NULL_R) + 1]] <- data.frame(readout = k, model = lab, rho = a$rho,
                                                                  n_rec = a$n_paper_sig_same_sign_P05, true = all(lost == LOST)) }
  }
}
NULL_R <- do.call(rbind, NULL_R)
OBS$null_mean <- OBS$null_q95 <- OBS$perm_p_rho <- OBS$perm_p_recovered <- NA
for (i in which(OBS$tf_set == "all shared" & OBS$model != "reference")) {
  z <- NULL_R[NULL_R$readout == OBS$readout[i] & NULL_R$model == OBS$model[i], ]
  OBS$null_mean[i] <- mean(z$rho); OBS$null_q95[i] <- quantile(z$rho, .95)
  OBS$perm_p_rho[i] <- mean(z$rho >= OBS$rho[i]); OBS$perm_p_recovered[i] <- mean(z$n_rec >= OBS$n_paper_sig_same_sign_P05[i])
}
write.csv(OBS, "paper2B_vs_scATAC_agreement.csv", row.names = FALSE)
cat("\n=== agreement with paper tab 2B ===\n")
show <- OBS[, c("readout", "model", "tf_set", "n_shared", "rho", "rho_p", "pct_same_sign", "n_paper_sig",
                "pct_same_sign_paper_sig", "n_paper_sig_same_sign_P05", "null_q95", "perm_p_rho", "perm_p_recovered")]
print(transform(show, rho = round(rho, 3), rho_p = signif(rho_p, 2), pct_same_sign = round(pct_same_sign),
                pct_same_sign_paper_sig = round(pct_same_sign_paper_sig), null_q95 = round(null_q95, 3),
                perm_p_rho = round(perm_p_rho, 3), perm_p_recovered = round(perm_p_recovered, 3)), row.names = FALSE)

## ---- TF-level table for the paper's significant TFs ------------------------------------
W <- Reduce(function(a, b) merge(a, b, by = "TF", all.x = TRUE), c(list(paper), TFL))
W <- merge(W, setNames(ref_ar[, c("TF", "t")], c("TF", "AR_repro_t")), by = "TF", all.x = TRUE)
W <- W[order(W$p_tf), ]
write.csv(W, "paper2B_vs_scATAC_TFlevel.csv", row.names = FALSE)
cat("\n=== paper FDR<0.01 TFs: sign of effect in each readout (+ = higher in BAP1-lost/inactivated) ===\n")
S <- W[W$p_adj_tf < 0.01, ]
sg <- function(t, p) ifelse(is.na(t), "  .", paste0(ifelse(t > 0, "+", "-"), ifelse(p < 0.05, "*", " ")))
print(data.frame(TF = S$TF, paper = ifelse(S$estimate_tf > 0, "+", "-"),
                 motif = sg(S$motif_unadj_t, S$motif_unadj_P), motif_adj = sg(S$motif_adj_t, S$motif_adj_P),
                 genescore = sg(S$genescore_unadj_t, S$genescore_unadj_P), gs_adj = sg(S$genescore_adj_t, S$genescore_adj_P),
                 AR_repro = ifelse(is.na(S$AR_repro_t), ".", ifelse(S$AR_repro_t > 0, "+", "-"))), row.names = FALSE)
cat("(* = P < 0.05 in that readout; . = TF not measured)\n")

## ---- figure -----------------------------------------------------------------------------
P1 <- OBS[OBS$tf_set == "all shared", ]
P1$label <- paste0(P1$readout, ifelse(P1$model == "reference", "", paste0("\n", P1$model)))
P1$label <- factor(P1$label, levels = rev(P1$label))
p1 <- ggplot(P1, aes(rho, label)) +
  geom_vline(xintercept = 0, linewidth = .3, colour = "grey60") +
  geom_errorbarh(data = subset(P1, !is.na(null_q95)), aes(xmin = -null_q95, xmax = null_q95), height = .25, colour = "grey75") +
  geom_point(aes(colour = model == "reference"), size = 2.6) +
  scale_colour_manual(values = c(`FALSE` = "#1b4f72", `TRUE` = "grey35"), guide = "none") +
  gtheme_no_rot + xlab("Spearman rho with paper tab 2B estimates") + ylab(NULL) +
  ggtitle("Agreement with the TCGA MESO BAP1 inferred-TF-activity list",
          subtitle = "grey bars: central 95% of the exact label-permutation null (126 assignments)") +
  theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
sc_panel <- function(col, ttl) {
  d <- W[!is.na(W[[col]]), ]; r <- cor(d[[col]], d$estimate_tf, method = "spearman")
  d$lab <- ifelse(d$p_adj_tf < 0.001, d$TF, NA)
  ggplot(d, aes(estimate_tf, .data[[col]])) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = .2) + geom_vline(xintercept = 0, colour = "grey85", linewidth = .2) +
    geom_point(colour = "grey65", size = 1.1) + geom_point(data = subset(d, p_adj_tf < 0.01), colour = "grey20", size = 1.4) +
    geom_text_repel(aes(label = lab), size = 2, max.overlaps = 30, seed = 1, min.segment.length = 0) +
    gtheme_no_rot + xlab("paper estimate") + ylab("our moderated t") + ggtitle(sprintf("%s (rho = %.2f)", ttl, r)) +
    theme(plot.title = element_text(size = 9))
}
p2 <- sc_panel("motif_unadj_t", "scATAC motif activity") | sc_panel("genescore_unadj_t", "scATAC gene score")
ggsave("Plots/paper2B_vs_scATAC_agreement.pdf", p1 / p2 + plot_layout(heights = c(1.1, 1)), width = 9, height = 9)
cat("\nDONE\n")
