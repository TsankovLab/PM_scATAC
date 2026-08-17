###############################################################
# R2_Q3 : SUBCLONE TREE -- one leaf per subclone per sample per modality, showing
#         scATAC <-> scRNA agreement, with the CNV arm driving each clone labelled.
#
# Companion to R2_Q3_crossmodal_congruence.R (same clone definitions) and the
# lead-in to the P4 Visium figure.
#
# Leaves    : every clone of every sample in BOTH modalities under the recommended
#             config (arm_cnv_zscore + cohort-adjusted multimodal filter + top-10 arms
#             + hierarchical/silhouette). Clonal sample -> one pooled leaf.
#             Config is reproduced INLINE -- arm_cnv_mmadjfilt/ on disk is the K=5 run.
#             UNSUPERVISED calls only: P4's chr8q-targeted clones are a second partition
#             of the SAME cells (scRNA c1 vs 8qAMP profile r = 1.00), so including them
#             double-counted P4 and forced within-assay cherries.
# Profile   : per-leaf mean arm score, z-scored against that sample's own arm
#             distribution, so scATAC log2FC and scRNA expression-difference share a scale.
# Tree      : 1-Pearson over the 39 arms, ward.D2.
# Driver    : per leaf, the arm that most separates it from the OTHER leaves of the same
#             sample+modality (the between-clone driver). For pooled/clonal leaves, the
#             arm that most separates it from the tree mean.
# Agreement : (a) leaf colour/strip by modality; (b) per sample, whether each leaf's
#             nearest neighbour in the tree is the same tumour's leaf from the OTHER
#             modality -- the cross-modal pairing that failed at 1-Mb bin level;
#             (c) the between-clone arm-delta |r| per sample, printed alongside.
#
# Outputs: subclone_tree_leaves.csv, subclone_tree_pairing.csv,
#          Plots/R2_Q3_subclone_tree_crossmodal.pdf  (single compact ggtree page)
#
# SUBCLONAL, BOTH-ASSAY SAMPLES ONLY -- clonal samples (P14, P23, P3) carry no clone
# structure, and single-modality samples (P10 scATAC-only, P13 scRNA-only) cannot be
# compared across assays; all are dropped from this figure.
###############################################################
suppressMessages({ library(cluster); library(ape); library(ggtree); library(ggplot2) })
set.seed(1)
SC <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(SC, "git_repo_claude", "R2_Q3")); dir.create("Plots", showWarnings = FALSE)

SRC <- "arm_cnv_zscore"; KEEP <- 10; MINLEAF <- 20; K_CLADE <- 6
COL_AT <- "#2a78d6"; COL_RN <- "#eb6834"; COL_HI <- "#1baf7a"
modcol <- c(scATAC = COL_AT, scRNA = COL_RN)

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
arm_delta <- function(A, cl){ tb <- sort(table(cl), decreasing = TRUE); if (length(tb) < 2) return(NULL)
  rowMeans(A[, cl == names(tb)[1], drop = FALSE]) - rowMeans(A[, cl == names(tb)[2], drop = FALSE]) }

## ---------------- STEP 1: clone calls (recommended config) ------------------
## Cohort-adjusted multimodal arm score: an arm's within-sample percentile MINUS its
## leave-one-out mean percentile across the other tumours of the same modality. Technical
## arms score high in every tumour and are demoted; a real subclonal arm scores high in
## one and survives. Top 10 arms then feed the silhouette clone caller.
MM <- read.csv("arm_level_multimodal_arms.csv", stringsAsFactors = FALSE)
MM$pct <- ave(MM$bic_gain, MM$sample, MM$modality, FUN = function(z){
  z[!is.finite(z)] <- min(z[is.finite(z)], na.rm = TRUE); rank(z) / length(z) })
MM$adj <- MM$pct - ave(MM$pct, MM$arm, MM$modality,
  FUN = function(z) if (length(z) > 1) (sum(z) - z) / (length(z) - 1) else 0)

fs <- list.files(SRC, "_arm\\.rds$", full.names = TRUE)
objs <- lapply(fs, readRDS); names(objs) <- sub("_arm\\.rds$", "", basename(fs))
armlv <- rownames(objs[[1]]$A)

res <- list(); DELTA <- list()
for (nm in names(objs)){
  o <- objs[[nm]]; A <- o$A
  t <- MM[MM$sample == o$sample & MM$modality == o$modality, ]
  v <- setNames(t$adj[match(rownames(A), t$arm)], rownames(A)); v[!is.finite(v)] <- -Inf
  keep <- names(sort(v, decreasing = TRUE))[1:min(KEEP, nrow(A))]
  keep <- rownames(A)[rownames(A) %in% keep]
  cc <- call_clones_arm(A[keep, , drop = FALSE])
  if (cc$k > 1) DELTA[[nm]] <- arm_delta(A, cc$cl)
  res[[nm]] <- list(sample = o$sample, modality = o$modality, A = A, k = cc$k,
                    sil = ifelse(is.na(cc$sil), 0, cc$sil), cl = cc$cl) }

## ---------------- STEP 2: one leaf per clone --------------------------------
## Each leaf's profile is that clone's mean arm score, z-scored against its OWN sample's
## arm distribution -- that is what puts scATAC log2FC and scRNA expression-difference on
## a common scale so they can share one tree.
prof <- list(); meta <- list()
addleaf <- function(nm, cells, tag, cls){
  o <- res[[nm]]; A <- o$A
  mu <- mean(rowMeans(A)); sdv <- sd(rowMeans(A))          # sample's own arm distribution
  leaf <- paste(o$sample, o$modality, tag)
  prof[[leaf]] <<- (rowMeans(A[, cells, drop = FALSE]) - mu) / sdv
  meta[[leaf]] <<- data.frame(leaf = leaf, sample = o$sample, modality = o$modality,
                              tag = tag, n_cells = length(cells), class = cls,
                              stringsAsFactors = FALSE) }

## keep only tumours that are subclonal in BOTH assays -- a single-modality sample
## (P10 scATAC-only, P13 scRNA-only) cannot contribute a cross-modal comparison.
subcl <- sapply(res, function(r) r$k > 1)
BOTH <- intersect(sub("_scATAC$", "", names(which(subcl))[grepl("_scATAC$", names(which(subcl)))]),
                  sub("_scRNA$",  "", names(which(subcl))[grepl("_scRNA$",  names(which(subcl)))]))
cat("tumours kept (subclonal in both assays):", paste(sort(BOTH), collapse = ", "), "\n")

for (nm in names(res)){
  o <- res[[nm]]
  if (o$k <= 1 || !(o$sample %in% BOTH)) next                                  # clonal sample -> not shown
  g <- split(names(o$cl), paste0("c", o$cl))
  g <- g[sapply(g, length) >= MINLEAF]
  for (i in names(g)) addleaf(nm, g[[i]], i, "subclonal") }

meta <- do.call(rbind, meta); rownames(meta) <- meta$leaf
leaves <- meta$leaf
Z <- t(sapply(prof[leaves], function(v) v[armlv])); colnames(Z) <- armlv

## ---------------- STEP 3: driver arm per leaf -------------------------------
## The arm that most separates this clone from its SIBLING in the same sample+assay
## (not from the tree as a whole) -- i.e. what actually distinguishes the two clones.
zall <- colMeans(Z, na.rm = TRUE)
meta$driver <- NA_character_; meta$driver_dz <- NA_real_
for (l in leaves){
  sib <- leaves[meta$sample == meta[l, "sample"] & meta$modality == meta[l, "modality"] &
                meta$class == meta[l, "class"] & leaves != l]
  ref <- if (length(sib)) colMeans(Z[sib, , drop = FALSE], na.rm = TRUE) else zall
  dz  <- Z[l, ] - ref; b <- names(which.max(abs(dz)))
  meta[l, "driver"] <- sprintf("%s %s", b, ifelse(dz[b] < 0, "loss", "gain"))
  meta[l, "driver_dz"] <- round(dz[b], 2) }
write.csv(data.frame(meta, round(Z, 3), check.names = FALSE),
          "subclone_tree_leaves.csv", row.names = FALSE)

## ---------------- STEP 4: tree + cross-modal pairing ------------------------
D  <- as.dist(1 - cor(t(Z), use = "pairwise.complete.obs"))
hc <- hclust(D, "ward.D2")
Dm <- as.matrix(D)

## cross-modal pairing: is a leaf's nearest neighbour the same tumour, other modality?
pair <- do.call(rbind, lapply(leaves, function(l){
  d <- Dm[l, ]; d[l] <- Inf; nn <- names(which.min(d))
  data.frame(leaf = l, sample = meta[l, "sample"], modality = meta[l, "modality"],
             class = meta[l, "class"], nn = nn,
             nn_same_tumour = meta[nn, "sample"] == meta[l, "sample"],
             nn_other_modality = meta[nn, "modality"] != meta[l, "modality"],
             stringsAsFactors = FALSE) }))
pair$cross_modal_pair <- pair$nn_same_tumour & pair$nn_other_modality
write.csv(pair, "subclone_tree_pairing.csv", row.names = FALSE)
psamp <- tapply(pair$cross_modal_pair, pair$sample, mean)

## The quantitative cross-modal statistic: clone labels are arbitrary across assays, so
## correlate the between-clone DELTA (largest clone minus second largest) over all 39
## arms instead. It answers "do the same arms demarcate the split in both assays".
both <- intersect(sub("_scATAC$", "", grep("_scATAC$", names(DELTA), value = TRUE)),
                  sub("_scRNA$",  "", grep("_scRNA$",  names(DELTA), value = TRUE)))
dcor <- setNames(sapply(both, function(s)
  abs(cor(DELTA[[paste0(s, "_scATAC")]][armlv], DELTA[[paste0(s, "_scRNA")]][armlv]))), both)

cat("\n=== leaves ===\n"); print(meta[, c("sample","modality","tag","n_cells","class","driver","driver_dz")], row.names = FALSE)
cat("\n=== cross-modal nearest-neighbour pairing, fraction of leaves per sample ===\n"); print(round(psamp, 2))
cat("\n=== between-clone arm-delta |r| ===\n"); print(round(dcor, 2))

## ---------------- STEP 5: circular figure -----------------------------------
phy <- as.phylo(hc)
dat <- data.frame(
  label    = meta$leaf,
  sample   = meta$sample,
  modality = meta$modality,
  clone    = meta$tag,
  driver   = meta$driver,
  n_cells  = meta$n_cells,
  paired   = pair$cross_modal_pair[match(meta$leaf, pair$leaf)],
  stringsAsFactors = FALSE)
## assay is carried by tip colour, so it is left out of the label
dat$tip <- sprintf("%s %s  %s", dat$sample, dat$clone, dat$driver)

## per-sample background wedge, drawn only where that tumour's leaves form a clade
sampu <- sort(unique(dat$sample))
sampcol <- setNames(colorRampPalette(c("#cfe3f5","#f6dfd0","#d8ecdd","#efe0ef",
                                       "#fdf2cc","#dfe6ef","#f7dede","#e3e9d5"))(length(sampu)), sampu)
mono <- list()
for (sm in sampu){
  tips <- dat$label[dat$sample == sm]
  nd <- tryCatch(ape::getMRCA(phy, tips), error = function(e) NA)
  if (is.na(nd)) next
  if (setequal(ape::extract.clade(phy, nd)$tip.label, tips)) mono[[sm]] <- nd }
cat(sprintf("\nsamples forming a clean clade (both modalities together): %s\n",
            paste(names(mono), collapse = ", ")))

## mutual cross-modal clone pairs -- each leaf is the other's nearest neighbour, same
## tumour, opposite assay. Stronger than the tumour merely forming a clade.
mut <- pair[pair$cross_modal_pair &
            pair$leaf == pair$nn[match(pair$nn, pair$leaf)], c("leaf", "nn")]
mut <- mut[!duplicated(t(apply(mut, 1, sort))), , drop = FALSE]
cat(sprintf("mutual cross-modal clone pairs: %d  (%s)\n", nrow(mut),
            paste(sprintf("%s <-> %s", mut$leaf, mut$nn), collapse = "; ")))

p <- ggtree(phy, layout = "fan", open.angle = 12, size = 0.45, colour = "grey35") %<+% dat
hl <- do.call(rbind, lapply(names(mono), function(sm)
        data.frame(node = mono[[sm]], sample = sm, stringsAsFactors = FALSE)))
if (!is.null(hl))
  p <- p + geom_hilight(data = hl, mapping = aes(node = node, fill = sample),
                        alpha = 0.55, extend = 1.15, show.legend = FALSE) +
           scale_fill_manual(values = sampcol)
p <- p +
  geom_tippoint(aes(colour = modality, size = n_cells), stroke = 0) +
  geom_tiplab(aes(label = tip, colour = modality, angle = angle), size = 1.85,
              offset = 0.03, show.legend = FALSE) +
  scale_colour_manual(values = modcol, name = NULL) +
  scale_size_continuous(range = c(1, 3), trans = "sqrt", name = NULL,
                        breaks = c(100, 5000), labels = c("100 cells", "5000")) +
  labs(title = sprintf("Subclone tree - %d clones, %d subclonal tumours",
                       nrow(dat), length(sampu)),
       subtitle = paste0("arm-level CNV, 1 - Pearson over 39 chromosome arms  |  leaf = sample / clone / driver arm; colour = assay\n",
                         "shaded wedge = the tumour's clones form one clade")) +
  theme_tree() +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"),
        legend.position = "inside", legend.position.inside = c(0.5, 0.5),
        legend.key.size = unit(2.8, "mm"), legend.text = element_text(size = 5.6),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.spacing.y = unit(0, "mm"), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(2, 2, 2, 2))
p <- p + xlim(-0.25, max(p$data$x, na.rm = TRUE) * 1.72)

ggsave("Plots/R2_Q3_subclone_tree_crossmodal.pdf", p, width = 7.2, height = 7.0,
       device = cairo_pdf)
cat("\nDONE -> Plots/R2_Q3_subclone_tree_crossmodal.pdf\n")
