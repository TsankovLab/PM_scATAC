###############################################################################
# STEP 10 -- is the P4 chr8q clone grouped together AMONG MALIGNANT SPOTS? (Visium)
#
# Step 9 maps the chr8q-amplified clone and gives Moran's I of the continuous score.
# This step quantifies the clustering of the clone LABEL itself, restricted to the
# malignant spots, so the answer cannot come from tumour-versus-stroma layout:
#
#   1. same spot calls as step 9 (malignant from the deconvolution weights; chr8q score
#      referenced to the normal spots; k-means k = 2 on the score within malignant spots)
#   2. neighbours = malignant spots that are physically adjacent on the Visium hexagonal
#      grid (centre distance <= 1.3 x the spot pitch; up to six per spot)
#   3. statistics, each against 999 shuffles of the clone labels over the malignant
#      spots only (positions and tumour geometry fixed):
#        - same-clone neighbour fraction: of the malignant neighbours of a chr8q-amp spot,
#          the fraction that is also chr8q-amp (random expectation = fraction of amp spots)
#        - assortativity: (same-label edge fraction - expected) / (1 - expected);
#          0 = clones intermixed at random, 1 = fully segregated
#        - patches: number of connected chr8q-amp patches and size of the largest one
#   4. two controls for "any expression score is smooth in tissue":
#        - the same pipeline on every other chromosome arm (k-means k = 2 on that arm's
#          score), so chr8q can be ranked against the genome
#        - the chr8q score after regressing out the malignant fraction of the spot
#          (tumour purity), which a copy-number signal in mixed spots necessarily tracks
#        - the chr8q score RELATIVE TO THE REST OF THE GENOME of the same spot (arm score
#          minus the score of all genes on the other autosomal arms). The raw score of
#          step 9 turned out to be shared by every arm on the fresh sections (Spearman
#          ~0.9 with any other arm, same two-way split), i.e. it is a genome-wide signal;
#          only the relative score is specific to chr8q
#
# Output: P4_visium_chr8q_clustering.csv        (one row per slide)
#         P4_visium_arm_clustering.csv          (every arm, every slide)
#         Plots/R2_Q3_P4_visium_chr8q_clustering.pdf           compact: observed vs shuffled, per section
#         Plots/R2_Q3_P4_visium_chr8q_clustering_relative.pdf  same, genome-relative chr8q score
#         Plots/R2_Q3_P4_visium_chr8q_clustering_full.pdf      all panels incl. the controls
###############################################################################
suppressMessages({ library(Seurat); library(EnsDb.Hsapiens.v86); library(ggplot2)
                   library(patchwork); library(Matrix) })
source("00_common.R")
set.seed(1)
slides   <- c("A1", "B1", "C1", "D1")
PREP     <- c(A1 = "fresh", B1 = "fresh", C1 = "snap-frozen", D1 = "snap-frozen")
N_PERM   <- 999
MIN_GENES_ARM <- 100         # arms with fewer detected genes are not scored

## ---- STEP 1: genes per chromosome arm ---------------------------------------
edb <- EnsDb.Hsapiens.v86
ge <- genes(edb, columns = c("gene_name","seq_name","gene_seq_start","gene_seq_end"),
            filter = SeqNameFilter(as.character(1:22)), return.type = "data.frame")
ge <- ge[ge$gene_name != "", ]
ge$arm <- armof(paste0("chr", ge$seq_name), (ge$gene_seq_start + ge$gene_seq_end) / 2e6)
arm_genes <- lapply(split(ge$gene_name, ge$arm), unique)

## ---- STEP 2: helpers ---------------------------------------------------------
get_spot_coords <- function(srt){                       # as in step 9
  img <- Images(srt)[1]; io <- srt@images[[img]]; ih <- nrow(io@image)
  co <- GetTissueCoordinates(srt, image = img)
  if ("imagerow" %in% colnames(co)) data.frame(x = as.numeric(co$imagecol), y = ih - as.numeric(co$imagerow))
  else data.frame(x = as.numeric(co[, 2]), y = ih - as.numeric(co[, 1]))
}
## adjacency among the given spots: neighbours on the hexagonal grid
hex_adjacency <- function(xy, pitch){
  d <- as.matrix(dist(xy)); A <- Matrix((d <= 1.3 * pitch & d > 0) * 1, sparse = TRUE); A
}
## edge statistics of a 0/1 label on adjacency A (each undirected edge counted twice,
## which cancels in every ratio)
edge_stats <- function(A, z){
  aa <- sum(z * (A %*% z)); deg <- rowSums(A); a_any <- sum(z * deg)
  ll <- sum((1 - z) * (A %*% (1 - z))); tot <- sum(deg)
  pa <- a_any / tot                                      # degree-weighted amp proportion
  same <- (aa + ll) / tot; expd <- pa^2 + (1 - pa)^2
  c(same_amp = aa / a_any, same_low = ll / (tot - a_any), assort = (same - expd) / (1 - expd))
}
## connected patches of the labelled spots
patches <- function(A, z){
  idx <- which(z == 1); if (!length(idx)) return(c(n_patch = 0, largest = 0))
  B <- A[idx, idx, drop = FALSE]; comp <- integer(length(idx)); k <- 0
  for (i in seq_along(idx)) if (comp[i] == 0){ k <- k + 1; q <- i; comp[i] <- k
    while (length(q)){ nb <- which(B[q[1], ] != 0); nb <- nb[comp[nb] == 0]; comp[nb] <- k; q <- c(q[-1], nb) } }
  c(n_patch = k, largest = max(table(comp)))
}
## observed value, permutation null (labels shuffled over the malignant spots), z and p
perm_test <- function(A, z, with_patches = FALSE){
  obs <- edge_stats(A, z); if (with_patches) obs <- c(obs, patches(A, z))
  nul <- replicate(N_PERM, { zp <- sample(z); r <- edge_stats(A, zp); if (with_patches) r <- c(r, patches(A, zp)); r })
  list(obs = obs, null = nul, mean = rowMeans(nul), sd = apply(nul, 1, sd),
       p_hi = (rowSums(nul >= obs) + 1) / (N_PERM + 1), p_lo = (rowSums(nul <= obs) + 1) / (N_PERM + 1))
}
split2 <- function(s){ km <- kmeans(scale(s), 2, nstart = 25)$cluster     # 1 = higher-score cluster
  as.numeric(km == which.max(tapply(s, km, mean))) }

## ---- STEP 3: per slide ---------------------------------------------------------
spa <- readRDS(SPATIAL)[slides]
## MAL_CALL: which spots count as malignant
##   sum (default)  tumour fractions outweigh the normal ones, the rule of step 9
##   label          the spot's label-transfer annotation (highest prediction score) is Malignant, using the
##                  label transfer from the P4 scRNA with the study's cell types (spa_all_rna_int_1123.rds:
##                  Malignant, SMC, Mesothelium, Myeloid, T, NK, B, Fibroblast, VEC, LEC)
MAL_CALL <- Sys.getenv("MAL_CALL", "sum"); SUF <- if (MAL_CALL == "sum") "" else paste0("_", MAL_CALL)
if (MAL_CALL == "label") LT <- readRDS(file.path(dirname(SPATIAL), "spa_all_rna_int_1123.rds"))[slides]
tab <- list(); arms_tab <- list(); nulls <- list(); nulls_rel <- list()
for (sl in slides){
  obj <- spa[[sl]]
  pred <- GetAssayData(obj, assay = "predictions", slot = "data")
  mal <- as.numeric(pred["Epithelioid-1", ]) + as.numeric(pred["Epithelioid-2", ]) +
         as.numeric(pred["Ephitelioid-p", ])
  nor <- as.numeric(pred["Fibroblasts", ]) + as.numeric(pred["Tcells", ]) +
         as.numeric(pred["Myeloid", ]) + as.numeric(pred["Endothelial", ]) +
         as.numeric(pred["Bcells", ])
  is_mal <- mal > nor; n_sum <- sum(is_mal)
  if (MAL_CALL == "label") {
    pl <- GetAssayData(LT[[sl]], assay = "predictions", layer = "data"); pl <- pl[rownames(pl) != "max", colnames(obj)]
    lab <- rownames(pl)[apply(pl, 2, which.max)]
    cat(sl, "label-transfer annotation:\n"); print(table(lab))
    cat(sl, "Malignant label vs step-9 rule:\n"); print(table(label_malignant = lab == "Malignant", sum_rule = is_mal))
    mal <- as.numeric(pl["Malignant", ]); nor <- 1 - mal; is_mal <- lab == "Malignant"
  }
  mi <- which(is_mal); purity <- (mal / (mal + nor))[mi]
  cnt <- tryCatch(GetAssayData(obj, assay = "Spatial", layer = "counts"),
                  error = function(e) GetAssayData(obj, assay = "Spatial", slot = "counts"))
  ln <- log1p(sweep(cnt, 2, colSums(cnt), "/") * 1e4)
  xy <- get_spot_coords(obj)
  dall <- as.matrix(dist(xy)); diag(dall) <- Inf; pitch <- median(apply(dall, 1, min)); rm(dall)
  A <- hex_adjacency(xy[mi, ], pitch)
  arm_score <- function(g){ g <- intersect(g, rownames(ln)); s <- colMeans(as.matrix(ln[g, , drop = FALSE]))
    (s - mean(s[!is_mal]))[mi] }
  ## arm score minus the score of all autosomal genes outside the arm, same spot
  auto <- intersect(unique(unlist(arm_genes)), rownames(ln)); auto_sum <- colSums(ln[auto, , drop = FALSE])
  rel_score <- function(g){ g <- intersect(g, auto); sa <- colSums(ln[g, , drop = FALSE])
    r <- sa / length(g) - (auto_sum - sa) / (length(auto) - length(g)); (r - mean(r[!is_mal]))[mi] }
  depth <- log10(colSums(cnt))[mi]

  ## 3a. chr8q, exactly the step 9 call
  s8 <- arm_score(arm_genes[["8q"]]); z8 <- split2(s8)
  pt <- perm_test(A, z8, with_patches = TRUE)
  nulls[[sl]] <- data.frame(slide = sl, null = unname(pt$null["same_amp", ]))
  ## 3b. chr8q after removing the malignant fraction of the spot
  s8adj <- resid(lm(s8 ~ purity)); z8adj <- split2(s8adj); pta <- perm_test(A, z8adj)
  zsc <- function(p, k) unname((p$obs[k] - p$mean[k]) / p$sd[k])
  ## 3b'. chr8q relative to the rest of the genome
  s8rel <- rel_score(arm_genes[["8q"]]); z8rel <- split2(s8rel); ptr <- perm_test(A, z8rel, with_patches = TRUE)
  tab[[sl]] <- data.frame(slide = sl, prep = unname(PREP[sl]), n_malignant = length(mi),
    mean_malignant_neighbours = round(mean(rowSums(A)), 2), n_amp = sum(z8), frac_amp = round(mean(z8), 3),
    same_clone_nb_amp = round(pt$obs["same_amp"], 3), expected_amp = round(pt$mean["same_amp"], 3),
    fold_amp = round(pt$obs["same_amp"] / pt$mean["same_amp"], 2), z_amp = round(zsc(pt, "same_amp"), 1),
    same_clone_nb_low = round(pt$obs["same_low"], 3), expected_low = round(pt$mean["same_low"], 3),
    assortativity = round(pt$obs["assort"], 3), assort_null_sd = round(pt$sd["assort"], 3), p_perm = pt$p_hi["assort"], p_same_amp = pt$p_hi["same_amp"],
    n_patches = pt$obs["n_patch"], n_patches_null = round(pt$mean["n_patch"], 1),
    largest_patch = pt$obs["largest"], largest_patch_frac = round(pt$obs["largest"] / sum(z8), 2),
    largest_patch_null = round(pt$mean["largest"], 1), p_largest = pt$p_hi["largest"],
    rho_score_purity = round(cor(s8, purity, method = "spearman"), 3),
    assortativity_purity_adj = round(pta$obs["assort"], 3), z_purity_adj = round(zsc(pta, "assort"), 1),
    p_perm_purity_adj = pta$p_hi["assort"], agree_raw_vs_adj = round(mean(z8 == z8adj), 3),
    rho_score_depth = round(cor(s8, depth, method = "spearman"), 3),
    ## chr8q relative to the rest of the genome
    rel_frac_amp = round(mean(z8rel), 3), rel_same_clone_nb_amp = round(ptr$obs["same_amp"], 3),
    rel_expected_amp = round(ptr$mean["same_amp"], 3), rel_assortativity = round(ptr$obs["assort"], 3),
    rel_z = round(zsc(ptr, "assort"), 1), rel_p_perm = ptr$p_hi["assort"], rel_p_same_amp = ptr$p_hi["same_amp"],
    rel_largest_patch = ptr$obs["largest"], rel_largest_patch_null = round(ptr$mean["largest"], 1), rel_p_largest = ptr$p_hi["largest"],
    rel_rho_depth = round(cor(s8rel, depth, method = "spearman"), 3), rel_rho_purity = round(cor(s8rel, purity, method = "spearman"), 3),
    rel_agree_with_raw_split = round(max(mean(z8rel == z8), mean(z8rel != z8)), 3))
  nulls_rel[[sl]] <- data.frame(slide = sl, null = unname(ptr$null["same_amp", ]))

  ## 3c. every arm through the same pipeline (raw and purity-adjusted)
  for (a in names(arm_genes)){
    if (length(intersect(arm_genes[[a]], rownames(ln)[rowSums(cnt > 0) >= 10])) < MIN_GENES_ARM) next
    s <- arm_score(arm_genes[[a]]); z <- split2(s); za <- split2(resid(lm(s ~ purity)))
    sr <- rel_score(arm_genes[[a]]); zr <- split2(sr)
    arms_tab[[paste(sl, a)]] <- data.frame(slide = sl, arm = a, frac_high = round(mean(z), 3),
      assortativity = unname(edge_stats(A, z)["assort"]), assortativity_purity_adj = unname(edge_stats(A, za)["assort"]),
      assortativity_rel = unname(edge_stats(A, zr)["assort"]), frac_high_rel = round(mean(zr), 3),
      rho_score_purity = cor(s, purity, method = "spearman"),
      ## does this arm's split mark the same territory as the chr8q split?
      cor_with_8q_score = cor(s, s8, method = "spearman"), agree_with_8q_split = max(mean(z == z8), mean(z != z8)))
  }
  cat(sl, "done\n")
}
TAB <- do.call(rbind, tab); rownames(TAB) <- NULL
ARM <- do.call(rbind, arms_tab); rownames(ARM) <- NULL
rk <- function(v, arm) rank(-v)[arm == "8q"]
TAB$rank_8q_among_arms <- sapply(slides, function(sl) with(ARM[ARM$slide == sl, ], rk(assortativity, arm)))
TAB$rank_8q_purity_adj <- sapply(slides, function(sl) with(ARM[ARM$slide == sl, ], rk(assortativity_purity_adj, arm)))
TAB$rank_8q_rel <- sapply(slides, function(sl) with(ARM[ARM$slide == sl, ], rk(assortativity_rel, arm)))
TAB$n_arms <- sapply(slides, function(sl) sum(ARM$slide == sl))
write.csv(TAB, paste0("P4_visium_chr8q_clustering", SUF, ".csv"), row.names = FALSE)
write.csv(ARM, paste0("P4_visium_arm_clustering", SUF, ".csv"), row.names = FALSE)
cat("=== chr8q clone clustering among malignant spots ===\n"); print(t(TAB))

## ---- STEP 4: figure ----------------------------------------------------------------
hist_panel <- function(nl, obs, expd, title){
  NUL <- do.call(rbind, nl); NUL$slide <- factor(NUL$slide, slides)
  lab <- data.frame(slide = factor(TAB$slide, slides), obs = obs,
    txt = sprintf("observed %.2f\nexpected %.2f", obs, expd))
  ggplot(NUL, aes(null)) + geom_histogram(binwidth = .01, fill = "grey70", colour = NA) +
    geom_vline(data = lab, aes(xintercept = obs), colour = "#d62728", linewidth = .7) +
    geom_text(data = lab, aes(x = obs, y = Inf, label = txt), hjust = 1.08, vjust = 1.2, size = 2.6, colour = "#d62728") +
    facet_wrap(~ slide, nrow = 1) + coord_cartesian(xlim = c(0, 1)) +
    labs(x = "fraction of malignant neighbours of a chr8q-amp spot that are chr8q-amp",
         y = sprintf("label shuffles (n = %d)", N_PERM), title = title) +
    theme_classic(base_size = 9) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
}
pA  <- hist_panel(nulls, TAB$same_clone_nb_amp, TAB$expected_amp, "A  chr8q score as in step 9 (raw): chr8q-amp spots neighbour each other, among malignant spots")
pA2 <- hist_panel(nulls_rel, TAB$rel_same_clone_nb_amp, TAB$rel_expected_amp, "B  chr8q relative to the rest of the genome of the same spot")
kinds <- c(assortativity = "raw score", assortativity_purity_adj = "malignant fraction regressed out", assortativity_rel = "relative to rest of genome")
AL <- do.call(rbind, lapply(names(kinds), function(k) data.frame(ARM[, c("slide", "arm")], value = ARM[[k]], kind = kinds[[k]])))
AL$kind <- factor(AL$kind, kinds); AL$slide <- factor(AL$slide, slides)
AL$is8q <- AL$arm == "8q"
pB <- ggplot(AL, aes(slide, value)) +
  geom_jitter(data = AL[!AL$is8q, ], width = .15, height = 0, size = .9, colour = "grey60", alpha = .8) +
  geom_point(data = AL[AL$is8q, ], size = 2.6, colour = "#d62728") +
  geom_text(data = AL[AL$is8q, ], aes(label = "8q"), colour = "#d62728", size = 2.8, hjust = -0.6) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = .3) + facet_wrap(~ kind, nrow = 1) +
  labs(x = NULL, y = "assortativity of the two-way split\n(0 = intermixed, 1 = segregated)",
       title = "C  chr8q against every other chromosome arm (grey), same pipeline") +
  theme_classic(base_size = 9) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
PL <- rbind(data.frame(slide = TAB$slide, what = "observed", value = TAB$largest_patch),
            data.frame(slide = TAB$slide, what = "shuffled labels", value = TAB$largest_patch_null))
pC <- ggplot(PL, aes(slide, value, fill = what)) + geom_col(position = position_dodge(.7), width = .65) +
  scale_fill_manual(values = c(observed = "#d62728", `shuffled labels` = "grey70"), name = NULL) +
  labs(x = NULL, y = "spots in the largest\nconnected chr8q-amp patch", title = "D  largest patch (raw score)") +
  theme_classic(base_size = 9) + theme(legend.position = "bottom")
fig <- pA / pA2 / (pB + pC + plot_layout(widths = c(3, 1)))
ggsave(paste0("Plots/R2_Q3_P4_visium_chr8q_clustering", SUF, "_full.pdf"), fig, width = 10, height = 9)

## compact version of the test: per section, the shuffled-label distribution as a box and the observed value as a
## point. Main file = chr8q score as in step 9; _relative = chr8q relative to the rest of the genome of the spot.
compact <- function(nl, obs, p){
  NUL <- do.call(rbind, nl); NUL$slide <- factor(NUL$slide, slides)
  O <- data.frame(slide = factor(TAB$slide, slides), obs = obs)
  ggplot(NUL, aes(slide, null)) +
    geom_boxplot(fill = "grey80", width = .5, linewidth = .25, outlier.size = .1, colour = "grey35") +
    geom_point(data = O, aes(y = obs), colour = "#d62728", size = 1.4) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
    labs(x = NULL, y = "chr8q-amp fraction of\nneighbours of amp spots",
         subtitle = sprintf("red = observed, box = shuffled\np <= %s, all sections", signif(max(p), 1))) +
    theme_classic(base_size = 7) + theme(plot.subtitle = element_text(size = 5.5, colour = "grey25"), plot.margin = margin(2, 3, 2, 2))
}
ggsave(paste0("Plots/R2_Q3_P4_visium_chr8q_clustering", SUF, ".pdf"), compact(nulls, TAB$same_clone_nb_amp, TAB$p_same_amp), width = 1.7, height = 1.5)
ggsave(paste0("Plots/R2_Q3_P4_visium_chr8q_clustering", SUF, ".png"), compact(nulls, TAB$same_clone_nb_amp, TAB$p_same_amp), width = 1.7, height = 1.5, dpi = 300)
## the relative-score version is drawn tall and narrow
pr <- compact(nulls_rel, TAB$rel_same_clone_nb_amp, TAB$rel_p_same_amp) + scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
  labs(subtitle = sprintf("red = observed\nbox = shuffled labels\np <= %s, all sections", signif(max(TAB$rel_p_same_amp), 1)))
ggsave(paste0("Plots/R2_Q3_P4_visium_chr8q_clustering", SUF, "_relative.pdf"), pr, width = 1.6, height = 3.2)
ggsave(paste0("Plots/R2_Q3_P4_visium_chr8q_clustering", SUF, "_relative.png"), pr, width = 1.6, height = 3.2, dpi = 300)
cat("\nDONE\n")
