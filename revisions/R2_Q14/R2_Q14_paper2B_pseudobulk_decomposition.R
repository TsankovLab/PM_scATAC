###############################################################################
# R2_Q14 -- does a BULK-like scATAC profile agree with the paper's BAP1 TF list better
#           than tumour cells alone, and is any bulk BAP1 difference compositional?
#
# Bulk RNA / RPPA measure a cell mixture.  From the per-tumour x compartment means
# (R2_Q14_compartment_extract.R) the cell-weighted pseudobulk of each tumour is rebuilt:
#   all cells       every compartment, weighted by its cell count
#   non-malignant   every compartment except Malignant (>= MIN_CELLS non-malignant cells)
# For chromVAR deviations the cell-weighted mean of compartment means is exactly the mean
# over all cells; gene score is averaged on the linear scale then log2(x+1).
# Each is tested ~ BAP1 (as the paper) and scored against tab 2B with an exact
# label-permutation null, alongside malignant-only from the compartment analysis.
#
# Shift-share decomposition (motif deviations, 5 main compartments):
#   bulk difference (lost - retained) ~
#     composition  sum_c (mean weight_lost_c - mean weight_retained_c) * mean activity_c
#   + within       sum_c  mean weight_c * (mean activity_lost_c - mean activity_retained_c)
#   + remainder    (tumour-to-tumour covariance of weights and activities)
# Compartment activity means use tumours with >= MIN_DECOMP cells in that compartment.
#
# Output: pseudobulk_paper2B_agreement.csv, pseudobulk_decomposition_paperTFs.csv,
#         Plots/pseudobulk_paper2B_agreement.pdf
###############################################################################
suppressPackageStartupMessages({ library(limma); library(ggplot2) })
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R2_Q14")); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
MIN_CELLS <- 50; MIN_DECOMP <- 20
MAIN <- c("Malignant", "TNK", "Myeloid", "B_Plasma", "Stroma")

paper <- read.csv("affinity_regression/prepared/paper_2B.csv", stringsAsFactors = FALSE)
bap1  <- read.csv("BAP1_genescore_per_sample.csv", stringsAsFactors = FALSE)
status <- setNames(bap1$BAP1_status, bap1$Sample)
G   <- read.csv("compartment_groups.csv", stringsAsFactors = FALSE)
DEV <- as.matrix(read.csv("compartment_motif_dev.csv", row.names = 1, check.names = FALSE))
GSC <- as.matrix(read.csv("compartment_genescore_TF.csv", row.names = 1, check.names = FALSE))
G <- G[G$sample %in% names(status) & G$group %in% colnames(DEV), ]

pseudobulk <- function(M, linearise, comps, min_cells) {
  sams <- sort(unique(G$sample))
  tot  <- sapply(sams, function(s) sum(G$n_cells[G$sample == s & G$compartment %in% comps]))
  sams <- sams[tot >= min_cells]
  out <- sapply(sams, function(s) {
    g <- G[G$sample == s & G$compartment %in% comps, ]
    X <- M[, g$group, drop = FALSE]; if (linearise) X <- 2^X - 1
    v <- as.vector(X %*% g$n_cells) / sum(g$n_cells)
    if (linearise) log2(v + 1) else v })
  rownames(out) <- rownames(M); out[apply(out, 1, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE]
}
fit_t <- function(M, lost) {
  grp <- factor(ifelse(lost, "lost", "retained"), levels = c("retained", "lost"))
  tt <- topTable(eBayes(lmFit(M, model.matrix(~ grp)), robust = TRUE), coef = "grplost", number = Inf, sort.by = "none")
  data.frame(TF = rownames(tt), t = tt$t, P = tt$P.Value, stringsAsFactors = FALSE)
}
agree <- function(d) {
  m <- merge(d, paper, by = "TF"); s <- m$p_adj_tf < 0.01; same <- sign(m$t) == sign(m$estimate_tf)
  c(n_shared = nrow(m), rho = suppressWarnings(cor(m$t, m$estimate_tf, method = "spearman")),
    pct_same_sign_paper_sig = 100 * mean(same[s]), n_recovered = sum(s & same & m$P < 0.05))
}
test_one <- function(M, label, readout) {
  lost <- status[colnames(M)] == "lost"; a <- agree(fit_t(M, lost))
  combos <- combn(ncol(M), sum(lost))
  nul <- t(apply(combos, 2, function(ix) agree(fit_t(M, seq_len(ncol(M)) %in% ix))))
  data.frame(profile = label, readout = readout, n_tumours = ncol(M), n_lost = sum(lost),
             tumours = paste(colnames(M), collapse = " "), n_perm = ncol(combos),
             rho = a["rho"], null_q95 = quantile(nul[, "rho"], .95), perm_p_rho = mean(nul[, "rho"] >= a["rho"]),
             pct_same_sign_paper_sig = a["pct_same_sign_paper_sig"], n_recovered = a["n_recovered"],
             perm_p_recovered = mean(nul[, "n_recovered"] >= a["n_recovered"]), row.names = NULL)
}
allc <- unique(G$compartment); nonm <- setdiff(allc, "Malignant")
RES <- rbind(
  test_one(pseudobulk(DEV, FALSE, allc, MIN_CELLS), "all cells (pseudobulk)", "motif"),
  test_one(pseudobulk(GSC, TRUE,  allc, MIN_CELLS), "all cells (pseudobulk)", "genescore"),
  test_one(pseudobulk(DEV, FALSE, nonm, MIN_CELLS), "non-malignant (pseudobulk)", "motif"),
  test_one(pseudobulk(GSC, TRUE,  nonm, MIN_CELLS), "non-malignant (pseudobulk)", "genescore"))
prev <- read.csv("compartment_paper2B_agreement.csv", stringsAsFactors = FALSE)
prev <- prev[prev$compartment == "Malignant" & prev$model == "~BAP1", ]
RES <- rbind(RES, data.frame(profile = "malignant cells only", readout = prev$readout, n_tumours = prev$n_tumours,
                             n_lost = 4, tumours = NA, n_perm = prev$n_perm, rho = prev$rho, null_q95 = prev$null_q95,
                             perm_p_rho = prev$perm_p_rho, pct_same_sign_paper_sig = prev$pct_same_sign_paper_sig,
                             n_recovered = prev$n_recovered, perm_p_recovered = prev$perm_p_recovered))
write.csv(RES, "pseudobulk_paper2B_agreement.csv", row.names = FALSE)
cat("=== agreement with paper tab 2B: whole-tissue vs non-malignant vs tumour-cell profiles (~BAP1) ===\n")
print(transform(RES[, setdiff(names(RES), "tumours")], rho = round(rho, 3), null_q95 = round(null_q95, 3),
                perm_p_rho = round(perm_p_rho, 3), pct_same_sign_paper_sig = round(pct_same_sign_paper_sig),
                perm_p_recovered = round(perm_p_recovered, 3)), row.names = FALSE)

## ---- shift-share decomposition ---------------------------------------------------------
Gm <- G[G$compartment %in% MAIN, ]; sams <- sort(unique(Gm$sample)); lost_s <- status[sams] == "lost"
Wt <- sapply(MAIN, function(cp) sapply(sams, function(s) { x <- Gm$n_cells[Gm$sample == s & Gm$compartment == cp]; if (length(x)) x else 0 }))
Wt <- Wt / rowSums(Wt)
tfs <- intersect(paper$TF, rownames(DEV))
act_mean <- function(cp, tf, which) {
  g <- Gm[Gm$compartment == cp & Gm$n_cells >= MIN_DECOMP & Gm$sample %in% sams[which], ]
  if (!nrow(g)) return(NA_real_); mean(DEV[tf, g$group])
}
DEC <- do.call(rbind, lapply(tfs, function(tf) {
  mL <- sapply(MAIN, act_mean, tf = tf, which = lost_s); mR <- sapply(MAIN, act_mean, tf = tf, which = !lost_s)
  mA <- sapply(MAIN, act_mean, tf = tf, which = rep(TRUE, length(sams)))
  wL <- colMeans(Wt[lost_s, , drop = FALSE]); wR <- colMeans(Wt[!lost_s, , drop = FALSE])
  comp <- sum((wL - wR) * mA, na.rm = TRUE)
  within_c <- ((wL + wR) / 2) * (mL - mR)
  pb <- as.vector(sapply(sams, function(s) { g <- Gm[Gm$sample == s, ]; sum(DEV[tf, g$group] * g$n_cells) / sum(g$n_cells) }))
  data.frame(TF = tf, bulk_diff = mean(pb[lost_s]) - mean(pb[!lost_s]), composition = comp,
             within_total = sum(within_c, na.rm = TRUE), t(setNames(within_c, paste0("within_", MAIN))),
             top_within_compartment = MAIN[which.max(abs(within_c))], check.names = FALSE)
}))
DEC$remainder <- DEC$bulk_diff - DEC$composition - DEC$within_total
DEC$composition_share <- abs(DEC$composition) / (abs(DEC$composition) + abs(DEC$within_total))
DEC <- merge(DEC, paper[, c("TF", "estimate_tf", "p_adj_tf")], by = "TF"); DEC <- DEC[order(DEC$p_adj_tf), ]
write.csv(DEC, "pseudobulk_decomposition_paperTFs.csv", row.names = FALSE)
cat("\n=== cell composition by BAP1 status used in the decomposition (mean fraction of the 5 main compartments) ===\n")
print(round(rbind(lost = colMeans(Wt[lost_s, ]), retained = colMeans(Wt[!lost_s, ])), 3))
cat("\n=== paper FDR<0.01 TFs: whole-tissue BAP1 difference split into composition vs within-compartment ===\n")
D <- DEC[DEC$p_adj_tf < 0.01, ]
print(data.frame(TF = D$TF, paper = ifelse(D$estimate_tf > 0, "+", "-"), bulk_diff = signif(D$bulk_diff, 2),
                 composition = signif(D$composition, 2), within = signif(D$within_total, 2),
                 comp_share = round(D$composition_share, 2), largest_within = D$top_within_compartment,
                 malignant_within = signif(D$within_Malignant, 2)), row.names = FALSE)
cat(sprintf("\nmedian composition share, paper FDR<0.01 TFs: %.2f | all paper TFs: %.2f\n",
            median(D$composition_share, na.rm = TRUE), median(DEC$composition_share, na.rm = TRUE)))
cat("largest within-compartment contributor, paper FDR<0.01 TFs:\n"); print(table(D$top_within_compartment))
cat(sprintf("does the bulk difference match the paper direction? %d of %d paper FDR<0.01 TFs\n",
            sum(sign(D$bulk_diff) == sign(D$estimate_tf)), nrow(D)))
cat(sprintf("does the COMPOSITION term alone match the paper direction? %d of %d\n",
            sum(sign(D$composition) == sign(D$estimate_tf)), nrow(D)))

RES$profile <- factor(RES$profile, levels = rev(c("malignant cells only", "non-malignant (pseudobulk)", "all cells (pseudobulk)")))
p <- ggplot(RES, aes(rho, profile)) + geom_vline(xintercept = 0, colour = "grey60", linewidth = .3) +
  geom_errorbarh(aes(xmin = -null_q95, xmax = null_q95), height = .2, colour = "grey80") +
  geom_point(size = 2.6, colour = "#1b4f72") + facet_wrap(~ readout) + gtheme_no_rot +
  xlab("Spearman rho with paper tab 2B estimates") + ylab(NULL) +
  ggtitle("Does a bulk-like scATAC profile reproduce the BAP1 TF list better than tumour cells?",
          subtitle = "grey bars: central 95% of the exact label-permutation null") +
  theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
ggsave("Plots/pseudobulk_paper2B_agreement.pdf", p, width = 8, height = 3.4)
cat("\nDONE\n")
