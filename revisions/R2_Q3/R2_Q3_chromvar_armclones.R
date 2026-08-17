###############################################################
# R2_Q3 : chromVAR TF activity across subclones -- REDONE ON THE CURRENT (arm-level)
#         CLONE CALLS.
#
# The earlier run (R2_Q3_chromvar_clones.R) used the 1-Mb binarized-TileMatrix clone
# labels in scATAC_cnv/<S>_subclones.csv. Those calls have since been superseded: the
# current definition is GC-corrected arm-level CNV + cohort-adjusted multimodal arm
# filter (top-10 arms) + hierarchical/silhouette, and the partitions differ materially
# (e.g. P11 was 1074/906, is now 617/1363; P4 was 2869/211, is now 340/2740). So every
# TF conclusion has to be re-derived against the labels actually in use.
#
# scATAC only (chromVAR is an ATAC assay). Samples = those subclonal in scATAC under
# the current calls: P1, P4, P5, P8, P10, P11, P12. Where k>2 the two largest robust
# clones are compared, as before.
#
# Per sample:
#   eta^2 per motif   = between-clone variance / total  -> how much TF activity the
#                       clone split explains
#   ARI               = chromVAR-only k=2 clustering vs the CNV clones (does TF
#                       activity independently recover the same partition?)
#   Wilcoxon + BH     = which motifs differ, and by how much (difference of mean z)
# Across samples:
#   recurrence        = how many samples a motif is strongly differential in, with the
#                       AP-1/bZIP + SMARCC1 accessibility axis FLAGGED as technical
#                       (it is anti-correlated with nFrags and floors everywhere --
#                       see R2_Q3_tead_vs_nfrags.R); do not read it as clonal biology.
#
# Outputs: chromvar_armclones_mirror.csv, chromvar_armclones_diff_<S>.csv,
#          chromvar_armclones_recurrence.csv, Plots/R2_Q3_chromvar_armclones.pdf
###############################################################
suppressMessages({ library(ArchR); library(cluster); library(ggplot2); library(patchwork) })
addArchRThreads(1); set.seed(1)
SC  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUT <- file.path(SC, "git_repo_claude", "R2_Q3"); setwd(OUT); dir.create("Plots", showWarnings = FALSE)
SRC <- "arm_cnv_zscore"; KEEP <- 10
ZCACHE <- "chromvar_z_cache.rds"
COL_AT <- "#2a78d6"; COL_HI <- "#1baf7a"; COL_TECH <- "#b0413e"
TECH <- c("JUN","JUNB","JUND","FOS","FOSB","FOSL1","FOSL2","BACH1","BACH2","NFE2",
          "SMARCC1","JDP2","BATF","BATF3","NFE2L1","NFE2L2","MAF","MAFF","MAFK","MAFG")

## ---------------- STEP 1: rebuild the current clone calls -------------------
## Reproduced INLINE rather than read from arm_cnv_mmadjfilt/, because that directory can
## hold whichever K ran last (its path depends on MODE, not K). Same recipe as the tree
## scripts: cohort-adjusted multimodal arm score -> top 10 arms -> silhouette caller.
call_clones_arm <- function(A, kmax = 6){
  ok <- apply(A, 2, function(z) all(is.finite(z))); A <- A[, ok, drop = FALSE]
  none <- list(k = 1, sil = NA, cl = setNames(rep(1L, ncol(A)), colnames(A)))
  if (ncol(A) < 40 || nrow(A) < 2) return(none)
  Ac <- A - rowMeans(A); d <- dist(t(Ac))
  hc <- hclust(d, "ward.D2"); best <- none
  for (k in 2:min(kmax, ncol(A) - 1)){
    cl <- cutree(hc, k); if (length(unique(cl)) < k) next
    if (any(table(cl) < max(20, 0.05 * ncol(A)))) next
    s <- mean(silhouette(cl, d)[, 3])
    if (is.na(best$sil) || s > best$sil) best <- list(k = k, sil = s, cl = cl) }
  if (is.na(best$sil) || best$sil < 0.10)
    best <- list(k = 1, sil = best$sil, cl = setNames(rep(1L, ncol(A)), colnames(A)))
  best }

MM <- read.csv("arm_level_multimodal_arms.csv", stringsAsFactors = FALSE)
MM$pct <- ave(MM$bic_gain, MM$sample, MM$modality, FUN = function(z){
  z[!is.finite(z)] <- min(z[is.finite(z)], na.rm = TRUE); rank(z) / length(z) })
MM$adj <- MM$pct - ave(MM$pct, MM$arm, MM$modality,
  FUN = function(z) if (length(z) > 1) (sum(z) - z) / (length(z) - 1) else 0)

CL <- list()
for (f in list.files(SRC, "_scATAC_arm\\.rds$", full.names = TRUE)){
  o <- readRDS(f); A <- o$A
  t <- MM[MM$sample == o$sample & MM$modality == "scATAC", ]
  v <- setNames(t$adj[match(rownames(A), t$arm)], rownames(A)); v[!is.finite(v)] <- -Inf
  keep <- names(sort(v, decreasing = TRUE))[1:min(KEEP, nrow(A))]
  cc <- call_clones_arm(A[rownames(A) %in% keep, , drop = FALSE])
  if (cc$k > 1) CL[[o$sample]] <- cc$cl }
SAMPLES <- sort(names(CL))
cat("scATAC-subclonal samples under the current calls:", paste(SAMPLES, collapse = ", "), "\n")

## keep the two largest robust clones (>= max(20, 5%)), as in the earlier analysis
robust2 <- function(assign){
  tb <- sort(table(assign), decreasing = TRUE); thr <- max(20, 0.05 * length(assign))
  rob <- names(tb)[tb >= thr][1:2]
  out <- rep(NA_character_, length(assign)); names(out) <- names(assign)
  out[assign == rob[1]] <- "sc1"; out[assign == rob[2]] <- "sc2"
  out[!is.na(out)] }

## ---------------- STEP 2: chromVAR TF activity ------------------------------
## MotifMatrix assay "z" = per-cell deviation z-score per motif (869 cisBP motifs).
## Row names carry a numeric suffix (TFAP2B_1) that is stripped to the TF symbol.
## Cached to disk: pulling this out of the ArchR project is the slow step.
if (file.exists(ZCACHE)){
  Z <- readRDS(ZCACHE); cat("chromVAR z from cache\n")
} else {
  proj <- loadArchRProject(file.path(SC, "tumor_compartment", "scatac_ArchR"), showLogo = FALSE)
  se <- getMatrixFromProject(proj, useMatrix = "MotifMatrix", verbose = FALSE)
  Z <- as.matrix(assays(se)$z); rownames(Z) <- rowData(se)$name
  rownames(Z) <- make.unique(sub("_[0-9]+$", "", rownames(Z)))
  saveRDS(Z, ZCACHE) }
cat("chromVAR z:", nrow(Z), "motifs x", ncol(Z), "cells\n")

ARI <- function(a, b){ t <- table(a, b); n <- sum(t)
  ai <- rowSums(t); bj <- colSums(t)
  e <- sum(choose(ai, 2)) * sum(choose(bj, 2)) / choose(n, 2)
  m <- 0.5 * (sum(choose(ai, 2)) + sum(choose(bj, 2)))
  (sum(choose(t, 2)) - e) / (m - e) }

## ---------------- STEP 3: per sample, clone1 vs clone2 ----------------------
## Where a sample has k>2, the two LARGEST robust clones are compared, as before.
## Three complementary readouts per motif:
##   diff  = difference of mean chromVAR z between the clones (effect size, signed)
##   eta2  = fraction of that motif's variance explained by the clone split
##   fdr   = BH-adjusted Wilcoxon; note this is inflated by cell count, so it ranks
##           motifs within a sample but must NOT be compared across samples
## ARI compares an independent chromVAR-only k=2 clustering with the CNV clones, i.e.
## does TF activity recover the same partition without being told about it.
metr <- list(); diffs <- list()
for (s in SAMPLES){
  lab <- robust2(CL[[s]])
  cells <- intersect(names(lab), colnames(Z))
  if (length(cells) < 40){ cat("  skip", s, "- too few cells in the motif matrix\n"); next }
  lab <- lab[cells]; zt <- Z[, cells, drop = FALSE]
  i1 <- which(lab == "sc1"); i2 <- which(lab == "sc2"); n1 <- length(i1); n2 <- length(i2)
  m1 <- rowMeans(zt[, i1, drop = FALSE]); m2 <- rowMeans(zt[, i2, drop = FALSE]); gm <- rowMeans(zt)
  eta2 <- (n1 * (m1 - gm)^2 + n2 * (m2 - gm)^2) / rowSums((zt - gm)^2)
  pw  <- apply(zt, 1, function(v) tryCatch(wilcox.test(v[i1], v[i2])$p.value, error = function(e) NA))
  fdr <- p.adjust(pw, "BH")
  d <- data.frame(sample = s, TF = rownames(zt), mean_sc1 = m1, mean_sc2 = m2,
                  diff = m2 - m1, eta2 = eta2, p = pw, fdr = fdr, row.names = NULL)
  d$abs_pct <- rank(abs(d$diff)) / nrow(d)                     # within-sample percentile
  d$tech <- d$TF %in% TECH
  diffs[[s]] <- d
  write.csv(d[order(d$fdr), ], sprintf("chromvar_armclones_diff_%s.csv", s), row.names = FALSE)

  tv <- order(apply(zt, 1, var), decreasing = TRUE)[seq_len(min(100, nrow(zt)))]
  km <- kmeans(scale(t(zt[tv, , drop = FALSE])), centers = 2, nstart = 10)$cluster
  bio <- d[!d$tech, ]
  metr[[s]] <- data.frame(sample = s, n_sc1 = n1, n_sc2 = n2,
    n_sig_fdr05 = sum(fdr < 0.05, na.rm = TRUE), mean_eta2 = mean(eta2, na.rm = TRUE),
    max_eta2 = max(eta2, na.rm = TRUE), ARI_chromvar_vs_CNV = ARI(km, lab),
    top_TF_any = d$TF[which.max(abs(d$diff))],
    top_TF_nontech = bio$TF[which.max(abs(bio$diff))],
    TEAD_best_pct = max(d$abs_pct[grepl("^TEAD", d$TF)]),
    stringsAsFactors = FALSE) }

metr <- do.call(rbind, metr); rownames(metr) <- NULL
write.csv(metr, "chromvar_armclones_mirror.csv", row.names = FALSE)
cat("\n=== does TF activity mirror the clone split? ===\n")
print(transform(metr, mean_eta2 = round(mean_eta2, 4), max_eta2 = round(max_eta2, 3),
                ARI_chromvar_vs_CNV = round(ARI_chromvar_vs_CNV, 3),
                TEAD_best_pct = round(TEAD_best_pct, 3)), row.names = FALSE)

## ---------------- STEP 4: recurrence across samples -------------------------
## n_strong counts the samples where a motif lands in that sample's top 5% by |diff| --
## a rank, so it is immune to samples having very different dynamic range. The AP-1/bZIP
## + SMARCC1 motifs are flagged `tech`: they are the chromVAR accessibility axis
## (anti-correlated with nFrags), and they score high in every sample regardless.
D <- do.call(rbind, diffs)
rec <- do.call(rbind, lapply(split(D, D$TF), function(x)
  data.frame(TF = x$TF[1], n_samples = nrow(x),
             n_sig = sum(x$fdr < 0.05, na.rm = TRUE),
             n_strong = sum(x$abs_pct >= 0.95, na.rm = TRUE),
             mean_absdiff = mean(abs(x$diff)), med_eta2 = median(x$eta2, na.rm = TRUE),
             mean_pct = mean(x$abs_pct), same_sign = abs(sum(sign(x$diff))) == nrow(x),
             tech = x$tech[1], stringsAsFactors = FALSE)))
rec <- rec[order(-rec$n_strong, -rec$mean_absdiff), ]
write.csv(rec, "chromvar_armclones_recurrence.csv", row.names = FALSE)
cat("\n=== most recurrently clone-variable motifs (n_strong = samples in the top 5%) ===\n")
print(head(rec[, c("TF","n_samples","n_sig","n_strong","mean_absdiff","med_eta2","tech")], 25), row.names = FALSE)
cat("\n--- same, EXCLUDING the AP-1/bZIP technical axis ---\n")
print(head(rec[!rec$tech, c("TF","n_samples","n_sig","n_strong","mean_absdiff","med_eta2")], 20), row.names = FALSE)

## ---------------- STEP 5: figure --------------------------------------------
top <- head(rec$TF[!rec$tech], 12)
hm <- D[D$TF %in% top, ]
hm$TF <- factor(hm$TF, levels = rev(top))
p1 <- ggplot(metr, aes(reorder(sample, ARI_chromvar_vs_CNV), ARI_chromvar_vs_CNV)) +
  geom_col(fill = COL_AT, width = 0.65) + coord_flip() +
  labs(x = NULL, y = "ARI: chromVAR clustering vs CNV clones",
       title = "Does TF activity recover the clone split?") +
  theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank())
p2 <- ggplot(hm, aes(sample, TF, fill = diff)) + geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#08519c", mid = "white", high = "#a50f15", name = "diff\n(sc2-sc1)") +
  labs(x = NULL, y = NULL, title = "Recurrently clone-variable motifs (technical AP-1/bZIP axis excluded)") +
  theme_bw(base_size = 8) + theme(panel.grid = element_blank())
p3 <- ggplot(metr, aes(reorder(sample, TEAD_best_pct), TEAD_best_pct)) +
  geom_col(fill = COL_HI, width = 0.65) + coord_flip() + ylim(0, 1) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40") +
  labs(x = NULL, y = "best TEAD motif, |diff| percentile within sample",
       title = "Hippo-TEAD") +
  theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank())

ggsave("Plots/R2_Q3_chromvar_armclones.pdf",
       (p1 / p3 | p2) + plot_layout(widths = c(1, 1.25)) +
         plot_annotation(title = "chromVAR TF activity across subclones - current arm-level clone calls",
           theme = theme(plot.title = element_text(face = "bold", size = 10))),
       width = 9.5, height = 5.4, device = cairo_pdf)
cat("\nDONE -> Plots/R2_Q3_chromvar_armclones.pdf\n")
