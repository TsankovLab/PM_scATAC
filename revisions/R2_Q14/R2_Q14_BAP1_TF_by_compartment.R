###############################################################################
# R2_Q14 -- BAP1-associated TF regulators in OUR scATAC data, resolved by cellular
#           compartment, marked where TCGA (BAP1 ALTERED vs WILD-TYPE) also supports them,
#           and flagged where the effect sits in the tumour microenvironment (TME), which a
#           bulk measurement cannot resolve.
#
# Data: main scATAC ArchR project, per tumour x compartment means over cells
#       (R2_Q14_compartment_extract.R).  Compartments with enough tumours per BAP1 group:
#       Malignant, TNK (T + NK), Myeloid, B_Plasma, Stroma (fibroblast/endothelial/SMC).
#       A tumour enters a compartment if it has >= MIN_CELLS cells there.
# Model (per compartment, limma, one value per tumour):
#   PRIMARY     ~ BAP1 + sarcomatoid score   (histology is confounded with BAP1 in our cohort)
#   alongside   ~ BAP1
#   readouts    chromVAR motif deviation (TF activity) and TF gene score.
# Hit in a compartment (pre-specified):
#   motif activity P < P_HIT in the adjusted model, AND the TF's gene score is detectably
#   accessible in that compartment (mean >= compartment median) and moves in the SAME
#   direction -- i.e. the TF itself is present and its locus agrees with its motif activity.
#   "strong" additionally requires motif FDR < FDR_STRONG within the compartment.
# TCGA mark: TF mRNA, ALTERED vs WILD-TYPE (TCGA_genetic_BAP1_TFonly_unadjusted.csv),
#   "TCGA-supported" = FDR < 0.05 in the same direction as the scATAC hit.
# Localisation of each hit motif:
#   tumour-intrinsic  hit in Malignant only
#   shared            hit in Malignant and >= 1 TME compartment, same direction
#   TME-only          hit in >= 1 TME compartment; in Malignant P > P_ABSENT or opposite sign
#   TME, tumour trend hit in TME; Malignant same direction with P_HIT <= P <= P_ABSENT
#   discordant        hits in opposite directions across compartments
# Composition: shift-share decomposition of each hit motif's whole-tissue BAP1 difference
#   (cell-weighted pseudobulk of the 5 compartments, 11 tumours) into cell-composition and
#   within-compartment parts -- a bulk profile mixes both and cannot separate them.
#
# Output: BAP1_TF_by_compartment_all.csv (every motif x compartment),
#         BAP1_TF_by_compartment_hits.csv (one row per hit motif, localisation, TCGA mark),
#         Plots/BAP1_TF_by_compartment_dotplot.pdf, Plots/BAP1_TF_by_compartment_summary.pdf
###############################################################################
suppressPackageStartupMessages({ library(limma); library(ggplot2); library(patchwork) })
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R2_Q14")); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
MIN_CELLS <- 50; P_HIT <- 0.01; P_ABSENT <- 0.2; FDR_STRONG <- 0.1; N_SHOW <- 70
COMPS <- c("Malignant", "TNK", "Myeloid", "B_Plasma", "Stroma"); TME <- setdiff(COMPS, "Malignant")

bap1 <- read.csv("BAP1_genescore_per_sample.csv", stringsAsFactors = FALSE)
sarc <- read.csv("scATAC_sarcscore_per_sample.csv", stringsAsFactors = FALSE)
G    <- read.csv("compartment_groups.csv", stringsAsFactors = FALSE)
DEV  <- as.matrix(read.csv("compartment_motif_dev.csv", row.names = 1, check.names = FALSE))
GSC  <- as.matrix(read.csv("compartment_genescore_TF.csv", row.names = 1, check.names = FALSE))
TCGA <- read.csv("TCGA_genetic_BAP1_TFonly_unadjusted.csv", stringsAsFactors = FALSE)
status <- setNames(bap1$BAP1_status, bap1$Sample); sarcsc <- setNames(sarc$sarc_score_atac, sarc$sample)
symbol <- function(x) sub("\\.[0-9]+$", "", x)

fit <- function(M, lost, S, adjust) {
  grp <- factor(ifelse(lost, "lost", "retained"), levels = c("retained", "lost"))
  des <- if (adjust) model.matrix(~ grp + S) else model.matrix(~ grp)
  tt <- topTable(eBayes(lmFit(M, des), robust = TRUE), coef = "grplost", number = Inf, sort.by = "none")
  data.frame(t = tt$t, P = tt$P.Value, FDR = tt$adj.P.Val, logFC = tt$logFC, row.names = rownames(tt))
}

## ---- per-compartment tests ----------------------------------------------------------
ALL <- list(); INFO <- list()
for (cp in COMPS) {
  g <- G[G$compartment == cp & G$n_cells >= MIN_CELLS & G$sample %in% names(status) & G$sample %in% names(sarcsc), ]
  sams <- g$sample; lost <- status[sams] == "lost"; S <- sarcsc[sams]
  md <- DEV[, g$group, drop = FALSE]; colnames(md) <- sams
  gs <- GSC[, g$group, drop = FALSE]; colnames(gs) <- sams
  md <- md[apply(md, 1, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE]
  gs <- gs[apply(gs, 1, function(x) all(is.finite(x)) && sd(x) > 0), , drop = FALSE]
  ma <- fit(md, lost, S, TRUE); mu <- fit(md, lost, S, FALSE); ga <- fit(gs, lost, S, TRUE); gu <- fit(gs, lost, S, FALSE)
  gs_mean <- rowMeans(gs); detected <- names(gs_mean)[gs_mean >= median(gs_mean)]
  sym <- symbol(rownames(ma))
  d <- data.frame(compartment = cp, motif = rownames(ma), symbol = sym,
                  n_tumours = length(sams), n_lost = sum(lost),
                  motif_t_adj = ma$t, motif_P_adj = ma$P, motif_FDR_adj = ma$FDR, motif_logFC_adj = ma$logFC,
                  motif_t_unadj = mu[rownames(ma), "t"], motif_P_unadj = mu[rownames(ma), "P"],
                  gs_t_adj = ga[sym, "t"], gs_P_adj = ga[sym, "P"], gs_t_unadj = gu[sym, "t"],
                  gs_detected = sym %in% detected, stringsAsFactors = FALSE)
  d$gs_concordant <- !is.na(d$gs_t_adj) & sign(d$gs_t_adj) == sign(d$motif_t_adj)
  d$hit <- d$motif_P_adj < P_HIT & d$gs_detected & d$gs_concordant
  d$strong <- d$hit & d$motif_FDR_adj < FDR_STRONG
  d$unadj_same_direction <- sign(d$motif_t_unadj) == sign(d$motif_t_adj)
  ALL[[cp]] <- d
  INFO[[cp]] <- data.frame(compartment = cp, tumours = length(sams), lost = sum(lost), retained = sum(!lost),
                           residual_df_adjusted = length(sams) - 3, motifs = nrow(d),
                           motif_P01 = sum(d$motif_P_adj < P_HIT), motif_FDR10 = sum(d$motif_FDR_adj < FDR_STRONG),
                           hits = sum(d$hit), strong_hits = sum(d$strong))
}
ALL <- do.call(rbind, ALL); INFO <- do.call(rbind, INFO); rownames(ALL) <- NULL
tm <- TCGA[match(ALL$symbol, TCGA$TF), ]
ALL$TCGA_t <- tm$t; ALL$TCGA_FDR <- tm$FDR; ALL$TCGA_chr3p <- tm$chr3p
ALL$TCGA_mark <- ifelse(is.na(ALL$TCGA_t), "not measured in TCGA",
                  ifelse(ALL$TCGA_FDR < 0.05 & sign(ALL$TCGA_t) == sign(ALL$motif_t_adj), "TCGA-supported (FDR<0.05)",
                   ifelse(ALL$TCGA_FDR < 0.05, "TCGA opposite (FDR<0.05)", "TCGA ns")))
write.csv(ALL, "BAP1_TF_by_compartment_all.csv", row.names = FALSE)
cat("=== per-compartment tests (adjusted model) ===\n"); print(INFO, row.names = FALSE)
cat(sprintf("(expected motifs at P < %.2f by chance per compartment: ~%.0f)\n", P_HIT, P_HIT * INFO$motifs[1]))

## ---- localisation of hit motifs ---------------------------------------------------------
W <- reshape(ALL[, c("motif", "symbol", "compartment", "motif_t_adj", "motif_P_adj", "hit", "strong")],
             idvar = c("motif", "symbol"), timevar = "compartment", direction = "wide")
hitcols <- paste0("hit.", COMPS); strongcols <- paste0("strong.", COMPS)
for (k in c(hitcols, strongcols)) W[[k]] <- W[[k]] %in% TRUE          # NA (motif not tested) -> FALSE
for (k in c(paste0("motif_t_adj.", COMPS), paste0("motif_P_adj.", COMPS))) W[[k]] <- as.numeric(W[[k]])
W <- W[rowSums(as.matrix(W[, hitcols])) > 0, ]
## row-wise with typed columns (apply() would coerce the data frame to character)
loc1 <- function(i) {
  hc <- COMPS[unlist(W[i, hitcols])]; dirs <- sign(unlist(W[i, paste0("motif_t_adj.", hc)]))
  if (length(unique(dirs)) > 1) return(c("discordant", paste(hc, collapse = ";"), "mixed"))
  dir <- if (dirs[1] > 0) "up in BAP1-lost" else "down in BAP1-lost"
  mal_hit <- "Malignant" %in% hc; tme_hit <- any(TME %in% hc)
  pm <- W[i, "motif_P_adj.Malignant"]; tmal <- W[i, "motif_t_adj.Malignant"]
  lab <- if (mal_hit && !tme_hit) "tumour-intrinsic" else if (mal_hit && tme_hit) "shared" else
         if (is.na(pm) || pm > P_ABSENT || sign(tmal) != dirs[1]) "TME-only" else "TME, tumour trend"
  c(lab, paste(hc, collapse = ";"), dir)
}
L <- t(vapply(seq_len(nrow(W)), loc1, character(3)))
W$localisation <- L[, 1]; W$hit_compartments <- L[, 2]; W$direction <- L[, 3]
W$strong_any <- rowSums(as.matrix(W[, strongcols])) > 0
tm <- TCGA[match(W$symbol, TCGA$TF), ]
W$TCGA_t <- tm$t; W$TCGA_FDR <- tm$FDR; W$TCGA_chr3p <- tm$chr3p
W$TCGA_mark <- ifelse(is.na(W$TCGA_t), "not measured",
               ifelse(W$TCGA_FDR < 0.05 & W$direction != "mixed" & ((W$TCGA_t > 0) == (W$direction == "up in BAP1-lost")), "TCGA-supported",
                ifelse(W$TCGA_FDR < 0.05, "TCGA opposite", "TCGA ns")))

## baseline specificity: compartment with the highest mean activity across all tumours
cm <- sapply(COMPS, function(cp) rowMeans(DEV[, G$group[G$compartment == cp & G$n_cells >= MIN_CELLS], drop = FALSE]))
W$baseline_top_compartment <- COMPS[max.col(cm[W$motif, , drop = FALSE], ties.method = "first")]

## composition vs within-compartment decomposition (whole-tissue pseudobulk, 5 compartments)
Gm <- G[G$compartment %in% COMPS & G$sample %in% names(status), ]; ts <- sort(unique(Gm$sample)); ls <- status[ts] == "lost"
Wt <- sapply(COMPS, function(cp) sapply(ts, function(s) { x <- Gm$n_cells[Gm$sample == s & Gm$compartment == cp]; if (length(x)) x else 0 }))
Wt <- Wt / rowSums(Wt); wL <- colMeans(Wt[ls, ]); wR <- colMeans(Wt[!ls, ])
am <- function(cp, mt, which) { g <- Gm[Gm$compartment == cp & Gm$n_cells >= 20 & Gm$sample %in% ts[which], ]
  if (!nrow(g)) NA_real_ else mean(DEV[mt, g$group]) }
dec <- t(sapply(W$motif, function(mt) {
  mA <- sapply(COMPS, am, mt = mt, which = rep(TRUE, length(ts)))
  mL <- sapply(COMPS, am, mt = mt, which = ls); mR <- sapply(COMPS, am, mt = mt, which = !ls)
  comp <- sum((wL - wR) * mA, na.rm = TRUE); within <- ((wL + wR) / 2) * (mL - mR)
  pb <- sapply(ts, function(s) { g <- Gm[Gm$sample == s, ]; sum(DEV[mt, g$group] * g$n_cells) / sum(g$n_cells) })
  c(bulk_diff = mean(pb[ls]) - mean(pb[!ls]), composition = comp, within_malignant = within[["Malignant"]],
    within_TME = sum(within[TME], na.rm = TRUE)) }))
W <- cbind(W, dec)
W$composition_share <- abs(W$composition) / (abs(W$composition) + abs(W$within_malignant) + abs(W$within_TME))
W$TME_share_of_within <- abs(W$within_TME) / (abs(W$within_malignant) + abs(W$within_TME))
W$bulk_resolvable <- W$localisation %in% c("tumour-intrinsic")
W$category_order <- match(W$localisation, c("tumour-intrinsic", "shared", "TME, tumour trend", "TME-only", "discordant"))
W <- W[order(W$category_order, -apply(abs(W[, paste0("motif_t_adj.", COMPS)]), 1, max, na.rm = TRUE)), ]
out <- W[, c("motif", "symbol", "localisation", "direction", "hit_compartments", "strong_any", "baseline_top_compartment",
             paste0("motif_t_adj.", COMPS), paste0("motif_P_adj.", COMPS), "TCGA_mark", "TCGA_t", "TCGA_FDR", "TCGA_chr3p",
             "bulk_diff", "composition", "within_malignant", "within_TME", "composition_share", "TME_share_of_within")]
write.csv(out, "BAP1_TF_by_compartment_hits.csv", row.names = FALSE)

cat("\n=== hit motifs by localisation ===\n"); print(table(W$localisation))
cat("\n=== hits per compartment x TCGA mark ===\n")
H <- ALL[ALL$hit, ]; print(table(H$compartment, H$TCGA_mark)[COMPS, , drop = FALSE])
cat("\n=== localisation x TCGA mark (hit motifs) ===\n"); print(table(W$localisation, W$TCGA_mark))
show <- function(d, n = 40) print(head(transform(d[, c("motif", "localisation", "direction", "hit_compartments", "strong_any",
  "baseline_top_compartment", "TCGA_mark", "TCGA_t", "composition_share", "TME_share_of_within")],
  TCGA_t = round(TCGA_t, 2), composition_share = round(composition_share, 2), TME_share_of_within = round(TME_share_of_within, 2)), n), row.names = FALSE)
cat("\n=== TCGA-supported hits ===\n"); show(W[W$TCGA_mark == "TCGA-supported", ])
cat("\n=== TME hits (TME-only / TME with tumour trend): regulators bulk RNA cannot localise ===\n")
show(W[W$localisation %in% c("TME-only", "TME, tumour trend"), ], 60)
cat("\n=== tumour-intrinsic and shared hits ===\n"); show(W[W$localisation %in% c("tumour-intrinsic", "shared"), ], 40)
cat("\nstrong hits (FDR <", FDR_STRONG, "):", sum(W$strong_any), "\n")

## ---- figures ------------------------------------------------------------------------------
sel <- W
if (nrow(sel) > N_SHOW) {                     # keep all TCGA-supported + strongest others
  keep <- sel$TCGA_mark == "TCGA-supported" | sel$strong_any
  rest <- sel[!keep, ]; rest <- rest[order(apply(rest[, paste0("motif_P_adj.", COMPS)], 1, min, na.rm = TRUE)), ]
  sel <- rbind(sel[keep, ], head(rest, max(0, N_SHOW - sum(keep))))
}
sel$row <- factor(sel$motif, levels = rev(sel$motif[order(sel$category_order, -apply(abs(sel[, paste0("motif_t_adj.", COMPS)]), 1, max, na.rm = TRUE))]))
D <- do.call(rbind, lapply(COMPS, function(cp) data.frame(row = sel$row, localisation = sel$localisation, compartment = cp,
  t = sel[[paste0("motif_t_adj.", cp)]], P = sel[[paste0("motif_P_adj.", cp)]], hit = sel[[paste0("hit.", cp)]])))
T2 <- data.frame(row = sel$row, localisation = sel$localisation, compartment = "TCGA mRNA",
                 t = sel$TCGA_t, lab = ifelse(sel$TCGA_mark == "TCGA-supported", "*", ifelse(sel$TCGA_mark == "TCGA opposite", "x", "")))
D$compartment <- factor(D$compartment, levels = c(COMPS, "TCGA mRNA")); T2$compartment <- factor(T2$compartment, levels = c(COMPS, "TCGA mRNA"))
lim <- quantile(abs(D$t), .98, na.rm = TRUE)
p1 <- ggplot() +
  geom_tile(data = T2, aes(compartment, row, fill = pmax(pmin(t, lim), -lim)), colour = "white") +
  geom_text(data = T2, aes(compartment, row, label = lab), size = 3, vjust = .75) +
  geom_point(data = D, aes(compartment, row, colour = pmax(pmin(t, lim), -lim), size = -log10(P))) +
  geom_point(data = subset(D, hit), aes(compartment, row, size = -log10(P)), shape = 21, colour = "black", stroke = .6) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey92", high = "#b2182b", limits = c(-lim, lim), name = "t (BAP1-lost\nvs retained)") +
  scale_fill_gradient2(low = "#2166ac", mid = "grey97", high = "#b2182b", limits = c(-lim, lim), guide = "none", na.value = "white") +
  scale_size_continuous(range = c(.4, 4), name = expression(-log[10]~P)) +
  facet_grid(localisation ~ ., scales = "free_y", space = "free_y") +
  gtheme_no_rot + xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 6),
        strip.text.y = element_text(angle = 0, size = 7)) +
  ggtitle("BAP1-associated TF activity by scATAC compartment (histology-adjusted)",
          subtitle = "ringed = hit (motif P<0.01, gene score detected and concordant); TCGA mRNA tile: * FDR<0.05 same direction, x opposite")
ggsave("Plots/BAP1_TF_by_compartment_dotplot.pdf", p1, width = 8, height = max(6, 0.14 * nrow(sel) + 2.5), limitsize = FALSE)

S1 <- as.data.frame(table(compartment = H$compartment, TCGA = ifelse(H$TCGA_mark == "TCGA-supported (FDR<0.05)", "TCGA-supported", "not supported / not measured")))
S1$compartment <- factor(S1$compartment, COMPS)
p2 <- ggplot(S1, aes(compartment, Freq, fill = TCGA)) + geom_col(width = .7) +
  scale_fill_manual(values = c(`TCGA-supported` = "#1b4f72", `not supported / not measured` = "grey80"), name = NULL) +
  gtheme_no_rot + xlab(NULL) + ylab("hit motifs") + ggtitle("Hits per compartment") + theme(axis.text.x = element_text(angle = 35, hjust = 1))
S2 <- as.data.frame(table(localisation = factor(W$localisation, c("tumour-intrinsic", "shared", "TME, tumour trend", "TME-only", "discordant")),
                          TCGA = ifelse(W$TCGA_mark == "TCGA-supported", "TCGA-supported", "not supported / not measured")))
p3 <- ggplot(S2, aes(localisation, Freq, fill = TCGA)) + geom_col(width = .7) +
  scale_fill_manual(values = c(`TCGA-supported` = "#1b4f72", `not supported / not measured` = "grey80"), guide = "none") +
  gtheme_no_rot + xlab(NULL) + ylab("hit motifs") + ggtitle("Where the hit is localised") + theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave("Plots/BAP1_TF_by_compartment_summary.pdf", (p2 | p3) + plot_layout(guides = "collect"), width = 9, height = 4)
cat("\nDONE\n")
