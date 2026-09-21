###############################################################################
# R2_Q14 -- are the bulk BAP1-associated TF hits (TCGA MESO supp tab 2B) coming from
#           tumour cells or from the microenvironment?
#
# For every scATAC compartment separately (main project; R2_Q14_compartment_extract.R):
#   BAP1 lost vs retained, limma per TF on per-tumour means, ~ BAP1 (as the paper) and
#   ~ BAP1 + sarcomatoid score; readouts = chromVAR motif deviation and TF gene score.
#   Agreement with tab 2B: Spearman(our t, paper estimate), sign agreement and recovery of
#   the paper's FDR<0.01 TFs, each against an EXACT label-permutation null over all ways of
#   assigning "lost" to that compartment's qualifying tumours.
# Two further views that bulk inference cannot separate from within-cell change:
#   specificity  which compartment carries each paper TF's signal at baseline (mean motif
#                deviation / gene score per compartment, z-scored across compartments)
#   composition  do BAP1-lost tumours differ in cell-type make-up? (fraction of all cells,
#                and of non-malignant cells, per compartment; Wilcoxon)
# Tumours need >= MIN_CELLS cells in a compartment to enter that compartment's test.
#
# Output: compartment_paper2B_agreement.csv, compartment_paper2B_TFlevel.csv,
#         compartment_TF_specificity.csv, compartment_composition_by_BAP1.csv,
#         Plots/compartment_paper2B_agreement.pdf, Plots/compartment_paperTF_specificity.pdf
###############################################################################
suppressPackageStartupMessages({ library(limma); library(ggplot2); library(patchwork) })
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R2_Q14")); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
MIN_CELLS <- 50
COMPS <- c("Malignant", "TNK", "Myeloid", "B_Plasma", "Stroma")

paper <- read.csv("affinity_regression/prepared/paper_2B.csv", stringsAsFactors = FALSE)
bap1  <- read.csv("BAP1_genescore_per_sample.csv", stringsAsFactors = FALSE)
sarc  <- read.csv("scATAC_sarcscore_per_sample.csv", stringsAsFactors = FALSE)
G     <- read.csv("compartment_groups.csv", stringsAsFactors = FALSE)
DEV   <- as.matrix(read.csv("compartment_motif_dev.csv", row.names = 1, check.names = FALSE))
GSC   <- as.matrix(read.csv("compartment_genescore_TF.csv", row.names = 1, check.names = FALSE))
status <- setNames(bap1$BAP1_status, bap1$Sample)
sarcsc <- setNames(sarc$sarc_score_atac, sarc$sample)

fit_t <- function(M, lost, S, adjust) {
  grp <- factor(ifelse(lost, "lost", "retained"), levels = c("retained", "lost"))
  des <- if (adjust) model.matrix(~ grp + S) else model.matrix(~ grp)
  tt <- topTable(eBayes(lmFit(M, des), robust = TRUE), coef = "grplost", number = Inf, sort.by = "none")
  data.frame(TF = rownames(tt), t = tt$t, P = tt$P.Value, stringsAsFactors = FALSE)
}
agree <- function(d) {
  m <- merge(d, paper, by = "TF"); s <- m$p_adj_tf < 0.01; same <- sign(m$t) == sign(m$estimate_tf)
  c(n_shared = nrow(m), rho = suppressWarnings(cor(m$t, m$estimate_tf, method = "spearman")),
    pct_same_sign_paper_sig = 100 * mean(same[s]), n_recovered = sum(s & same & m$P < 0.05))
}

AGR <- list(); TFL <- list(); INFO <- list()
for (cp in COMPS) {
  g <- G[G$compartment == cp & G$n_cells >= MIN_CELLS & G$sample %in% names(status) & G$sample %in% names(sarcsc), ]
  sams <- g$sample; lost <- status[sams] == "lost"; S <- sarcsc[sams]
  INFO[[cp]] <- data.frame(compartment = cp, n_tumours = length(sams), n_lost = sum(lost), n_retained = sum(!lost),
                           lost = paste(sams[lost], collapse = " "), retained = paste(sams[!lost], collapse = " "),
                           cells = sum(g$n_cells))
  if (sum(lost) < 2 || sum(!lost) < 2) { cat(cp, ": too few tumours per group, skipped\n"); next }
  mats <- list(motif = DEV[, g$group, drop = FALSE], genescore = GSC[, g$group, drop = FALSE])
  mats <- lapply(mats, function(M) { colnames(M) <- sams; M[apply(M, 1, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE] })
  combos <- combn(length(sams), sum(lost))
  for (rd in names(mats)) for (adj in c(FALSE, TRUE)) {
    if (adj && length(sams) < 5) next
    obs <- fit_t(mats[[rd]], lost, S, adj); a <- agree(obs)
    nul <- t(apply(combos, 2, function(ix) agree(fit_t(mats[[rd]], seq_along(sams) %in% ix, S, adj))))
    AGR[[length(AGR) + 1]] <- data.frame(compartment = cp, readout = rd,
      model = if (adj) "~BAP1 + sarc" else "~BAP1", n_tumours = length(sams), n_perm = ncol(combos),
      n_shared = a["n_shared"], rho = a["rho"], perm_p_rho = mean(nul[, "rho"] >= a["rho"]),
      null_q95 = quantile(nul[, "rho"], .95), pct_same_sign_paper_sig = a["pct_same_sign_paper_sig"],
      n_recovered = a["n_recovered"], perm_p_recovered = mean(nul[, "n_recovered"] >= a["n_recovered"]))
    TFL[[paste(cp, rd, adj)]] <- setNames(obs, c("TF", paste(cp, rd, if (adj) "adj" else "unadj", c("t", "P"), sep = "_")))
  }
}
INFO <- do.call(rbind, INFO); AGR <- do.call(rbind, AGR); rownames(AGR) <- NULL
cat("\n=== tumours per compartment (>=", MIN_CELLS, "cells) ===\n"); print(INFO, row.names = FALSE)
write.csv(AGR, "compartment_paper2B_agreement.csv", row.names = FALSE)
cat("\n=== agreement with paper tab 2B, by compartment ===\n")
print(transform(AGR, rho = round(rho, 3), perm_p_rho = round(perm_p_rho, 3), null_q95 = round(null_q95, 3),
                pct_same_sign_paper_sig = round(pct_same_sign_paper_sig), perm_p_recovered = round(perm_p_recovered, 3)),
      row.names = FALSE)

W <- Reduce(function(a, b) merge(a, b, by = "TF", all.x = TRUE), c(list(paper), TFL)); W <- W[order(W$p_tf), ]
write.csv(W, "compartment_paper2B_TFlevel.csv", row.names = FALSE)
S <- W[W$p_adj_tf < 0.01, ]
sg <- function(t, p) ifelse(is.na(t), " .", paste0(ifelse(t > 0, "+", "-"), ifelse(p < 0.05, "*", " ")))
tab <- data.frame(TF = S$TF, paper = ifelse(S$estimate_tf > 0, "+", "-"))
for (cp in COMPS) for (rd in c("motif", "genescore")) {
  tc <- paste(cp, rd, "unadj_t", sep = "_"); pc <- paste(cp, rd, "unadj_P", sep = "_")
  if (tc %in% names(S)) tab[[paste(substr(cp, 1, 6), if (rd == "motif") "mot" else "gs", sep = ".")]] <- sg(S[[tc]], S[[pc]])
}
cat("\n=== paper FDR<0.01 TFs: direction in each compartment, unadjusted (+ = higher in BAP1-lost; * P<0.05) ===\n")
print(tab, row.names = FALSE)

## ---- specificity: which compartment carries each paper TF at baseline ----------------------
spec <- function(M, label) {
  cm <- sapply(COMPS, function(cp) { gg <- G$group[G$compartment == cp & G$n_cells >= MIN_CELLS]
                                     rowMeans(M[, intersect(gg, colnames(M)), drop = FALSE]) })
  z <- t(scale(t(cm))); tfs <- intersect(paper$TF, rownames(z))
  z <- z[tfs, , drop = FALSE]; z <- z[rowSums(is.finite(z)) == ncol(z), , drop = FALSE]   # drop TFs flat across compartments
  data.frame(TF = rownames(z), readout = label, z, top_compartment = COMPS[max.col(z, ties.method = "first")],
             check.names = FALSE, row.names = NULL)
}
SP <- rbind(spec(DEV, "motif"), spec(GSC, "genescore"))
SP <- merge(SP, paper[, c("TF", "estimate_tf", "p_adj_tf")], by = "TF"); SP <- SP[order(SP$p_adj_tf, SP$readout), ]
write.csv(SP, "compartment_TF_specificity.csv", row.names = FALSE)
cat("\n=== paper FDR<0.01 TFs: compartment with the highest baseline signal ===\n")
sp2 <- reshape(SP[SP$p_adj_tf < 0.01, c("TF", "readout", "top_compartment")], idvar = "TF", timevar = "readout", direction = "wide")
print(sp2, row.names = FALSE)
cat("\ntop compartment among paper FDR<0.01 TFs (motif):\n"); print(table(SP$top_compartment[SP$readout == "motif" & SP$p_adj_tf < 0.01]))
cat("top compartment among all paper TFs (motif):\n"); print(table(SP$top_compartment[SP$readout == "motif"]))
cat("top compartment, gene score, paper FDR<0.01 TFs:\n"); print(table(SP$top_compartment[SP$readout == "genescore" & SP$p_adj_tf < 0.01]))

## ---- composition by BAP1 status --------------------------------------------------------------
cnt <- xtabs(n_cells ~ sample + compartment, G[G$sample %in% names(status), ])
frac_all <- prop.table(cnt, 1)
nonmal <- cnt[, setdiff(colnames(cnt), "Malignant"), drop = FALSE]; frac_nm <- prop.table(nonmal, 1)
lost_s <- status[rownames(cnt)] == "lost"
COMP <- do.call(rbind, lapply(COMPS, function(cp) {
  fa <- frac_all[, cp]; fn <- if (cp %in% colnames(frac_nm)) frac_nm[, cp] else rep(NA, nrow(cnt))
  data.frame(compartment = cp, mean_frac_all_lost = mean(fa[lost_s]), mean_frac_all_retained = mean(fa[!lost_s]),
             p_frac_all = suppressWarnings(wilcox.test(fa[lost_s], fa[!lost_s])$p.value),
             mean_frac_nonmalignant_lost = mean(fn[lost_s]), mean_frac_nonmalignant_retained = mean(fn[!lost_s]),
             p_frac_nonmalignant = if (all(is.na(fn))) NA else suppressWarnings(wilcox.test(fn[lost_s], fn[!lost_s])$p.value))
}))
write.csv(COMP, "compartment_composition_by_BAP1.csv", row.names = FALSE)
cat("\n=== cell-type composition by BAP1 status (all 11 tumours) ===\n")
print(transform(COMP, mean_frac_all_lost = round(mean_frac_all_lost, 3), mean_frac_all_retained = round(mean_frac_all_retained, 3),
                p_frac_all = signif(p_frac_all, 2), mean_frac_nonmalignant_lost = round(mean_frac_nonmalignant_lost, 3),
                mean_frac_nonmalignant_retained = round(mean_frac_nonmalignant_retained, 3), p_frac_nonmalignant = signif(p_frac_nonmalignant, 2)),
      row.names = FALSE)

## ---- figures ------------------------------------------------------------------------------------
A <- AGR; A$label <- paste(A$readout, A$model)
p1 <- ggplot(A, aes(rho, factor(compartment, rev(COMPS)))) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = .3) +
  geom_errorbarh(aes(xmin = -null_q95, xmax = null_q95), height = .25, colour = "grey80") +
  geom_point(aes(shape = perm_p_rho < 0.05), size = 2.4, colour = "#1b4f72") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), name = "perm p < 0.05") +
  facet_grid(model ~ readout) + gtheme_no_rot + xlab("Spearman rho with paper tab 2B estimates") + ylab(NULL) +
  ggtitle("Agreement of each scATAC compartment's BAP1 effect with the bulk TF-activity list",
          subtitle = "grey bars: central 95% of the exact label-permutation null")
ggsave("Plots/compartment_paper2B_agreement.pdf", p1, width = 8, height = 5.5)
H <- SP[SP$p_adj_tf < 0.01, ]
HL <- do.call(rbind, lapply(COMPS, function(cp) data.frame(TF = H$TF, readout = H$readout, compartment = cp, z = H[[cp]])))
HL$TF <- factor(HL$TF, levels = rev(unique(H$TF[order(-H$estimate_tf)])))
p2 <- ggplot(HL, aes(factor(compartment, COMPS), TF, fill = z)) + geom_tile(colour = "white") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", name = "z across\ncompartments") +
  facet_wrap(~ readout) + gtheme_no_rot + xlab(NULL) + ylab("paper FDR < 0.01 TFs (ordered by paper estimate)") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) + ggtitle("Where is each bulk BAP1-associated TF active at baseline?")
ggsave("Plots/compartment_paperTF_specificity.pdf", p2, width = 7, height = 6.5)
cat("\nDONE\n")
