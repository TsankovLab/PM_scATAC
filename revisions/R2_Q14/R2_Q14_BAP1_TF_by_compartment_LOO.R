###############################################################################
# R2_Q14 -- leave-one-tumour-out robustness of the per-compartment BAP1 TF hits.
#
# TME compartments rest on 2-4 BAP1-lost tumours, so a hit can be one tumour.  For every
# hit (compartment x motif, R2_Q14_BAP1_TF_by_compartment.R) the adjusted limma model is
# refitted dropping each tumour in turn (only when >= 2 lost and >= 2 retained remain).
#   loo_same_sign   fraction of refits with the same sign
#   loo_P05         fraction of refits with the same sign AND P < 0.05
#   robust          loo_P05 >= 0.8 and every possible refit was run
#   max_influence   the tumour whose removal weakens the effect most
# Output: BAP1_TF_by_compartment_LOO.csv; adds robust column to the hits table.
###############################################################################
suppressPackageStartupMessages(library(limma))
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R2_Q14"))
MIN_CELLS <- 50; COMPS <- c("Malignant", "TNK", "Myeloid", "B_Plasma", "Stroma")
bap1 <- read.csv("BAP1_genescore_per_sample.csv", stringsAsFactors = FALSE)
sarc <- read.csv("scATAC_sarcscore_per_sample.csv", stringsAsFactors = FALSE)
G    <- read.csv("compartment_groups.csv", stringsAsFactors = FALSE)
DEV  <- as.matrix(read.csv("compartment_motif_dev.csv", row.names = 1, check.names = FALSE))
ALL  <- read.csv("BAP1_TF_by_compartment_all.csv", stringsAsFactors = FALSE)
status <- setNames(bap1$BAP1_status, bap1$Sample); sarcsc <- setNames(sarc$sarc_score_atac, sarc$sample)
H <- ALL[ALL$hit, c("compartment", "motif", "motif_t_adj", "motif_P_adj")]

OUT <- list()
for (cp in COMPS) {
  h <- H[H$compartment == cp, ]; if (!nrow(h)) next
  g <- G[G$compartment == cp & G$n_cells >= MIN_CELLS & G$sample %in% names(status) & G$sample %in% names(sarcsc), ]
  sams <- g$sample; lost <- status[sams] == "lost"
  M <- DEV[, g$group, drop = FALSE]; colnames(M) <- sams
  M <- M[apply(M, 1, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE]
  fits <- list(); dropped <- character(0)
  for (i in seq_along(sams)) {
    keep <- -i; l <- lost[keep]
    if (sum(l) < 2 || sum(!l) < 2) next
    grp <- factor(ifelse(l, "lost", "retained"), levels = c("retained", "lost"))
    tt <- topTable(eBayes(lmFit(M[, keep, drop = FALSE], model.matrix(~ grp + sarcsc[sams[keep]])), robust = TRUE),
                   coef = "grplost", number = Inf, sort.by = "none")
    fits[[sams[i]]] <- tt; dropped <- c(dropped, sams[i])
  }
  for (j in seq_len(nrow(h))) {
    mt <- h$motif[j]; s0 <- sign(h$motif_t_adj[j])
    tv <- sapply(fits, function(tt) tt[mt, "t"]); pv <- sapply(fits, function(tt) tt[mt, "P.Value"])
    OUT[[length(OUT) + 1]] <- data.frame(compartment = cp, motif = mt, t_full = h$motif_t_adj[j], P_full = h$motif_P_adj[j],
      n_tumours = length(sams), n_lost = sum(lost), loo_fits = length(fits), loo_possible = length(sams),
      loo_same_sign = mean(sign(tv) == s0), loo_P05 = mean(sign(tv) == s0 & pv < 0.05),
      weakest_t = tv[which.min(s0 * tv)], max_influence = names(tv)[which.min(s0 * tv)],
      skipped = paste(setdiff(sams, dropped), collapse = " "))
  }
}
OUT <- do.call(rbind, OUT)
OUT$robust <- OUT$loo_P05 >= 0.8 & OUT$loo_fits == OUT$loo_possible
write.csv(OUT, "BAP1_TF_by_compartment_LOO.csv", row.names = FALSE)
cat("=== leave-one-tumour-out robustness of hits ===\n")
print(do.call(rbind, lapply(split(OUT, OUT$compartment), function(d) data.frame(compartment = d$compartment[1],
  hits = nrow(d), n_lost = d$n_lost[1], loo_fits = d$loo_fits[1], of_possible = d$loo_possible[1],
  robust = sum(d$robust), median_loo_P05 = median(d$loo_P05),
  top_influential_tumour = names(sort(table(d$max_influence), decreasing = TRUE))[1]))), row.names = FALSE)
cat("\nnote: refits dropping a lost tumour are impossible when only 2 lost tumours remain (B_Plasma), so no B_Plasma hit can be 'robust'\n")

HT <- read.csv("BAP1_TF_by_compartment_hits.csv", stringsAsFactors = FALSE)
HT$robust_compartments <- sapply(seq_len(nrow(HT)), function(i) {
  o <- OUT[OUT$motif == HT$motif[i] & OUT$robust, ]; paste(o$compartment, collapse = ";") })
HT$robust_any <- nzchar(HT$robust_compartments)
write.csv(HT, "BAP1_TF_by_compartment_hits.csv", row.names = FALSE)
cat("\n=== robust hits by localisation ===\n"); print(table(HT$localisation, robust = HT$robust_any))
show <- HT[HT$robust_any, c("motif", "localisation", "direction", "hit_compartments", "robust_compartments", "TCGA_mark", "composition_share", "TME_share_of_within")]
show$composition_share <- round(show$composition_share, 2); show$TME_share_of_within <- round(show$TME_share_of_within, 2)
cat("\n=== robust hits ===\n"); print(show[order(show$localisation), ], row.names = FALSE)
cat("\nDONE\n")
