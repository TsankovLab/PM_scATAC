# R1_Q3 – Visium spatial, revisit: is SNAI2 / RUNX2 expression really organised around
# SOX9+ malignant spots, or would any seed do?
#
# The earlier scripts (proximity_v2, radius_*, imputed_full_repeat) report that SNAI2 and
# RUNX2 decrease with distance from SOX9+ malignant spots. Three things were not controlled:
#   (a) SOX9+ is ~20 spots per section (SOX9 detected, SCT > 0), i.e. rare detection events
#       that are more likely in deeply sequenced, tumour-dense spots -- where every gene is
#       higher;
#   (b) spots are spatially autocorrelated, so Spearman / Wilcoxon p-values over spots are
#       not calibrated;
#   (c) the MAGIC-imputed version smooths expression over similar spots, which by itself
#       makes any two genes co-vary in space.
# Here the same statistic is recomputed and placed against three nulls, per section:
#   null 1 "other seed genes": seeds = positive spots of genes detected in a similar number
#           of malignant spots as SOX9 (0.7-1.4x). Asks: is SOX9 special as a seed?
#   null 2 "random seed spots": the same number of malignant spots drawn at random, matched
#           to the sequencing-depth distribution of the SOX9+ spots (depth quintiles).
#   null 3 (MAGIC only) seeds = top-N imputed spots of random expressed genes, N as for SOX9.
# Statistics: rho = Spearman (distance to nearest seed, target expression) over malignant
#   non-seed spots (negative = target higher near seeds; the statistic of proximity_v2),
#   and delta = mean target in ring1 (<= 1.5 spot pitch from a seed) - mean in far (> 5 pitch).
# Also: co-detection in the same spot (logistic regression adjusted for depth and
# malignant fraction).
#
# Outputs: R1_Q3/spatial_SOX9_RUNX2_SNAI2_revisit_table.csv
#          R1_Q3/spatial_SOX9_RUNX2_SNAI2_revisit_codetection.csv
#          R1_Q3/Plots/spatial_SOX9_RUNX2_SNAI2_revisit.pdf

set.seed(1234)
suppressMessages({ library(Seurat); library(Matrix); library(ggplot2) })
projdir <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
outdir  <- file.path(projdir, 'git_repo_claude', 'R1_Q3')
SPATIAL <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/meso_spatial/spa_all_magic.rds'
TARGETS <- c('SNAI2', 'RUNX2'); N_GENES <- 300; N_PERM <- 999; MIN_POS <- 10

spa <- readRDS(SPATIAL)
res <- list(); nulls <- list(); codet <- list(); info <- list()

seed_stats <- function(D, seeds, tv, pitch){            # D = distance matrix, seeds = indices, tv = target vector
  d <- apply(D[, seeds, drop = FALSE], 1, min); keep <- setdiff(seq_along(d), seeds)
  near <- keep[d[keep] <= 1.5 * pitch]; far <- keep[d[keep] > 5 * pitch]
  c(rho = suppressWarnings(cor(d[keep], tv[keep], method = 'spearman')),
    delta = if (length(near) >= 5 && length(far) >= 5) mean(tv[near]) - mean(tv[far]) else NA)
}

for (sl in names(spa)) {
  srt <- spa[[sl]]
  pred <- GetAssayData(srt, assay = 'predictions', layer = 'data')
  mal <- as.numeric(pred['Epithelioid-1', ]) + as.numeric(pred['Epithelioid-2', ]) + as.numeric(pred['Ephitelioid-p', ])
  nor <- as.numeric(pred['Fibroblasts', ]) + as.numeric(pred['Tcells', ]) + as.numeric(pred['Myeloid', ]) +
         as.numeric(pred['Endothelial', ]) + as.numeric(pred['Bcells', ])
  coords <- GetTissueCoordinates(srt, image = names(srt@images)[1])
  cells <- intersect(rownames(coords), colnames(srt)[mal > nor])
  purity <- setNames(mal / (mal + nor), colnames(srt))[cells]
  S <- GetAssayData(srt, assay = 'SCT', layer = 'data')[, cells]
  depth <- log10(colSums(GetAssayData(srt, assay = 'Spatial', layer = 'counts')[, cells]))
  D <- as.matrix(dist(coords[cells, 1:2])); Dn <- D; diag(Dn) <- Inf; pitch <- median(apply(Dn, 1, min)); rm(Dn)
  npos <- rowSums(S > 0); sox <- which(S['SOX9', ] > 0)
  info[[sl]] <- data.frame(slide = sl, n_malignant = length(cells), SOX9_pos = length(sox),
    SNAI2_pos = sum(S['SNAI2', ] > 0), RUNX2_pos = sum(S['RUNX2', ] > 0),
    depth_SOX9pos = round(median(depth[sox]), 3), depth_other = round(median(depth[-sox]), 3),
    depth_p = if (length(sox) >= 3) signif(wilcox.test(depth[sox], depth[-sox])$p.value, 2) else NA,
    purity_SOX9pos = round(median(purity[sox]), 3), purity_other = round(median(purity[-sox]), 3))
  message(sl, ': malignant ', length(cells), ' | SOX9+ ', length(sox))

  ## co-detection in the same spot, adjusted for depth and malignant fraction
  for (t in TARGETS) {
    df <- data.frame(y = as.numeric(S[t, ] > 0), sox9 = as.numeric(S['SOX9', ] > 0), depth = depth, purity = purity)
    cd <- data.frame(slide = sl, target = t, n_both = sum(df$y & df$sox9), frac_in_SOX9pos = round(mean(df$y[df$sox9 == 1]), 3),
                     frac_in_SOX9neg = round(mean(df$y[df$sox9 == 0]), 3), OR_adj = NA, p_adj = NA)
    if (length(sox) >= MIN_POS && cd$n_both >= 1) {
      fit <- summary(glm(y ~ sox9 + depth + purity, family = binomial, data = df))$coefficients
      cd$OR_adj <- round(exp(fit['sox9', 1]), 2); cd$p_adj <- signif(fit['sox9', 4], 2) }
    codet[[paste(sl, t)]] <- cd
  }
  if (length(sox) < MIN_POS) next

  ## null 1: genes detected in a similar number of malignant spots
  cand <- setdiff(names(npos)[npos >= 0.7 * length(sox) & npos <= 1.4 * length(sox)], c('SOX9', TARGETS))
  genes <- sample(cand, min(N_GENES, length(cand)))
  ## null 2: random seeds matched to the depth distribution of the SOX9+ spots
  qs <- cut(depth, quantile(depth, 0:5 / 5), include.lowest = TRUE, labels = FALSE); need <- table(factor(qs[sox], 1:5))
  rand_seeds <- function() unlist(lapply(1:5, function(q) { i <- which(qs == q); i[sample.int(length(i), need[q])] }))
  ## MAGIC: SOX9 seeds = top-N imputed spots (N as above), null 3 = top-N of random expressed genes
  M <- GetAssayData(srt, assay = 'MAGIC_to_impute', layer = 'data'); M <- M[, cells]
  topN <- function(g) order(M[g, ], decreasing = TRUE)[seq_along(sox)]
  mg <- sample(setdiff(intersect(rownames(M), names(npos)[npos >= 0.05 * length(cells)]), c('SOX9', TARGETS)), N_GENES)

  for (t in TARGETS) {
    tv <- as.numeric(S[t, ]); tm <- as.numeric(M[t, ])
    obs  <- seed_stats(D, sox, tv, pitch)
    n1   <- t(sapply(genes, function(g) seed_stats(D, which(S[g, ] > 0), tv, pitch)))
    n2   <- t(replicate(N_PERM, seed_stats(D, rand_seeds(), tv, pitch)))
    obsM <- seed_stats(D, topN('SOX9'), tm, pitch)
    n3   <- t(sapply(mg, function(g) seed_stats(D, topN(g), tm, pitch)))
    pl <- function(null, o) (sum(null <= o, na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1)   # one-sided: more negative rho than null
    pg <- function(null, o) (sum(null >= o, na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1)   # one-sided: larger delta than null
    res[[paste(sl, t)]] <- data.frame(slide = sl, target = t, n_SOX9pos = length(sox),
      rho = round(obs['rho'], 3), rho_null_genes = round(median(n1[, 'rho'], na.rm = TRUE), 3), p_rho_genes = signif(pl(n1[, 'rho'], obs['rho']), 2),
      rho_null_depth = round(median(n2[, 'rho']), 3), p_rho_depth = signif(pl(n2[, 'rho'], obs['rho']), 2),
      delta = round(obs['delta'], 3), delta_null_genes = round(median(n1[, 'delta'], na.rm = TRUE), 3), p_delta_genes = signif(pg(n1[, 'delta'], obs['delta']), 2),
      delta_null_depth = round(median(n2[, 'delta'], na.rm = TRUE), 3), p_delta_depth = signif(pg(n2[, 'delta'], obs['delta']), 2),
      magic_rho = round(obsM['rho'], 3), magic_rho_null = round(median(n3[, 'rho'], na.rm = TRUE), 3), p_magic_rho = signif(pl(n3[, 'rho'], obsM['rho']), 2),
      frac_null_genes_rho_negative = round(mean(n1[, 'rho'] < 0, na.rm = TRUE), 2), n_null_genes = length(genes))
    nulls[[paste(sl, t)]] <- rbind(
      data.frame(slide = sl, target = t, null = 'other seed genes\n(SCT)', rho = n1[, 'rho'], obs = obs['rho']),
      data.frame(slide = sl, target = t, null = 'random seed spots,\ndepth-matched (SCT)', rho = n2[, 'rho'], obs = obs['rho']),
      data.frame(slide = sl, target = t, null = 'other seed genes\n(MAGIC)', rho = n3[, 'rho'], obs = obsM['rho']))
  }
}

INFO <- do.call(rbind, info); RES <- do.call(rbind, res); COD <- do.call(rbind, codet); rownames(RES) <- rownames(COD) <- rownames(INFO) <- NULL
cat('\n=== detection and depth ===\n'); print(INFO, row.names = FALSE)
cat('\n=== co-detection in the same spot ===\n'); print(COD, row.names = FALSE)
cat('\n=== proximity against the nulls ===\n'); print(t(RES))
write.csv(merge(INFO, RES, by = 'slide', all = TRUE), file.path(outdir, 'spatial_SOX9_RUNX2_SNAI2_revisit_table.csv'), row.names = FALSE)
write.csv(COD, file.path(outdir, 'spatial_SOX9_RUNX2_SNAI2_revisit_codetection.csv'), row.names = FALSE)

NUL <- do.call(rbind, nulls); lv <- unique(NUL$null); NUL$null <- factor(NUL$null, lv)
OBS <- unique(NUL[, c('slide', 'target', 'null', 'obs')])
p <- ggplot(NUL, aes(null, rho)) + geom_hline(yintercept = 0, linetype = 2, linewidth = .3) +
  geom_boxplot(fill = 'grey85', colour = 'grey35', width = .55, linewidth = .3, outlier.size = .2) +
  geom_point(data = OBS, aes(y = obs), colour = '#d62728', size = 2) +
  facet_grid(target ~ slide) +
  labs(x = NULL, y = 'Spearman rho: distance to nearest seed vs expression\n(negative = higher near seeds)',
       title = 'SNAI2 / RUNX2 around SOX9+ malignant spots (red) against control seeds (boxes)') +
  theme_classic(base_size = 8) + theme(axis.text.x = element_text(angle = 35, hjust = 1), strip.background = element_blank(),
                                        strip.text = element_text(face = 'bold'), plot.title = element_text(size = 8.5, face = 'bold'))
ggsave(file.path(outdir, 'Plots', 'spatial_SOX9_RUNX2_SNAI2_revisit.pdf'), p, width = 2 + 1.6 * length(unique(NUL$slide)), height = 4.6)
ggsave(file.path(outdir, 'Plots', 'spatial_SOX9_RUNX2_SNAI2_revisit.png'), p, width = 2 + 1.6 * length(unique(NUL$slide)), height = 4.6, dpi = 150)
message('done')
