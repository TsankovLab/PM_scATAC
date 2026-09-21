###############################################################################
# R2_Q14 -- is each compartment's number of BAP1 TF hits more than chance?
#
# Leave-one-out showed no hit depends on a single tumour; it does not show that a
# compartment yields more hits than random labels would.  Here every possible assignment of
# the BAP1-lost label within a compartment (all combinations) is refitted with the SAME
# adjusted model and the SAME hit rule (motif P < 0.01, gene score detected and concordant).
#   expected_false_hits  mean hits over the permuted (non-true) assignments
#   empirical_FDR        expected_false_hits / observed hits
#   perm_p               fraction of assignments with >= observed hits (true one included)
# Per motif: perm_p_motif = fraction of assignments giving P <= observed in the same direction.
# Output: BAP1_TF_by_compartment_permFDR.csv; adds empirical_FDR / perm_p_motif to the hits table.
###############################################################################
suppressPackageStartupMessages(library(limma))
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R2_Q14"))
MIN_CELLS <- 50; P_HIT <- 0.01; COMPS <- c("Malignant", "TNK", "Myeloid", "B_Plasma", "Stroma")
bap1 <- read.csv("BAP1_genescore_per_sample.csv", stringsAsFactors = FALSE)
sarc <- read.csv("scATAC_sarcscore_per_sample.csv", stringsAsFactors = FALSE)
G    <- read.csv("compartment_groups.csv", stringsAsFactors = FALSE)
DEV  <- as.matrix(read.csv("compartment_motif_dev.csv", row.names = 1, check.names = FALSE))
GSC  <- as.matrix(read.csv("compartment_genescore_TF.csv", row.names = 1, check.names = FALSE))
ALL  <- read.csv("BAP1_TF_by_compartment_all.csv", stringsAsFactors = FALSE)
status <- setNames(bap1$BAP1_status, bap1$Sample); sarcsc <- setNames(sarc$sarc_score_atac, sarc$sample)
symbol <- function(x) sub("\\.[0-9]+$", "", x)

tfit <- function(M, lost, S) {
  grp <- factor(ifelse(lost, "lost", "retained"), levels = c("retained", "lost"))
  tt <- topTable(eBayes(lmFit(M, model.matrix(~ grp + S)), robust = TRUE), coef = "grplost", number = Inf, sort.by = "none")
  list(t = setNames(tt$t, rownames(tt)), P = setNames(tt$P.Value, rownames(tt)))
}
SUM <- list(); MOT <- list()
for (cp in COMPS) {
  g <- G[G$compartment == cp & G$n_cells >= MIN_CELLS & G$sample %in% names(status) & G$sample %in% names(sarcsc), ]
  sams <- g$sample; lost <- status[sams] == "lost"; S <- sarcsc[sams]
  md <- DEV[, g$group, drop = FALSE]; colnames(md) <- sams
  gs <- GSC[, g$group, drop = FALSE]; colnames(gs) <- sams
  md <- md[apply(md, 1, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE]
  gs <- gs[apply(gs, 1, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE]
  detected <- names(which(rowMeans(gs) >= median(rowMeans(gs))))
  sym <- symbol(rownames(md)); ok_gs <- sym %in% detected & sym %in% rownames(gs)
  hits_of <- function(l) {
    m <- tfit(md, l, S); gg <- tfit(gs, l, S)
    gt <- rep(NA_real_, nrow(md)); gt[ok_gs] <- gg$t[sym[ok_gs]]
    h <- m$P < P_HIT & ok_gs & !is.na(gt) & sign(gt) == sign(m$t)
    list(n = sum(h), t = m$t, P = m$P)
  }
  obs <- hits_of(lost)
  combos <- combn(length(sams), sum(lost)); true_ix <- which(apply(combos, 2, function(ix) all(seq_along(sams) %in% ix == lost)))
  perm <- lapply(seq_len(ncol(combos)), function(i) hits_of(seq_along(sams) %in% combos[, i]))
  nh <- sapply(perm, `[[`, "n"); nh_null <- nh[-true_ix]
  SUM[[cp]] <- data.frame(compartment = cp, tumours = length(sams), lost = sum(lost), assignments = ncol(combos),
    observed_hits = obs$n, expected_false_hits = mean(nh_null), null_95pct = unname(quantile(nh_null, .95)),
    null_max = max(nh_null), empirical_FDR = min(1, mean(nh_null) / max(obs$n, 1)), perm_p = mean(nh >= obs$n))
  hm <- ALL$motif[ALL$compartment == cp & ALL$hit]
  for (mt in hm) {
    s0 <- sign(obs$t[mt])
    MOT[[length(MOT) + 1]] <- data.frame(compartment = cp, motif = mt,
      perm_p_motif = mean(sapply(perm, function(p) sign(p$t[mt]) == s0 && p$P[mt] <= obs$P[mt])))
  }
}
SUM <- do.call(rbind, SUM); MOT <- do.call(rbind, MOT)
write.csv(SUM, "BAP1_TF_by_compartment_permFDR.csv", row.names = FALSE)
cat("=== are there more hits than random BAP1 labels give? ===\n")
print(transform(SUM, expected_false_hits = round(expected_false_hits, 1), empirical_FDR = round(empirical_FDR, 2),
                perm_p = round(perm_p, 3)), row.names = FALSE)
cat("\nminimum attainable per-motif permutation p per compartment:",
    paste(sprintf("%s 1/%d", SUM$compartment, SUM$assignments), collapse = " | "), "\n")

HT <- read.csv("BAP1_TF_by_compartment_hits.csv", stringsAsFactors = FALSE)
fdr <- setNames(SUM$empirical_FDR, SUM$compartment)
HT$best_compartment_empirical_FDR <- sapply(strsplit(HT$hit_compartments, ";"), function(hc) min(fdr[hc]))
HT$perm_p_motif_min <- sapply(seq_len(nrow(HT)), function(i) {
  x <- MOT$perm_p_motif[MOT$motif == HT$motif[i]]; if (length(x)) min(x) else NA })
HT$credible <- HT$best_compartment_empirical_FDR <= 0.25
write.csv(HT, "BAP1_TF_by_compartment_hits.csv", row.names = FALSE)
cat("\n=== hits in compartments with empirical FDR <= 0.25, by localisation ===\n")
print(table(HT$localisation, credible = HT$credible))
cr <- HT[HT$credible & HT$localisation %in% c("TME-only", "TME, tumour trend"), ]
cat("\n=== credible TME regulators ===\n")
print(cr[order(cr$hit_compartments, cr$direction), c("motif", "localisation", "direction", "hit_compartments",
          "best_compartment_empirical_FDR", "perm_p_motif_min", "TCGA_mark", "composition_share", "TME_share_of_within")],
      row.names = FALSE, digits = 2)
cat("\nDONE\n")
