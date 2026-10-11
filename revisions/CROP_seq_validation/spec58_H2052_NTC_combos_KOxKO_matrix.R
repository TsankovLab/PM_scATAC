# Spec 58: full KO x KO Spearman correlation matrix (rows = CROP-seq screen KO, columns = bulk NCI-H2052 KO; NTC - KO profiles) for each of the 15 combinations of the
# 4 candidate H2052 NTC samples (H2052_NTC_R3 + QC-dropped H28_NTC_R1/2/3), in each feature space (20 Cm modules, top 1000/2000/3000 variable genes). "Best diagonal":
#   diag_minus_off = mean(diagonal) - mean(off-diagonal); top1 = number of screen KOs whose best-matching bulk KO is the same KO; exact permutation p (all 8! relabelings of
#   the bulk KOs) for mean(diagonal). Ranking of the 15 combinations by diag_minus_off (mean over spaces). NB: 15 combinations x 4 spaces are tried, so the best is optimistic.
# needs: git_repo/files/scRNA_PM_modules.rds, spec_common.R, repro_common.R
# outputs: spec58_KOxKO_matrices.rds, spec58_diagonal_metrics.csv, fig3Q_KOxKO_matrix_{Cm,top1000,top2000,top3000}.{png,pdf}, fig3R_diagonal_summary.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3"); SGN <- -1; TOPNS <- c(1000, 2000, 3000)
b <- load_bulk(); meta0 <- b$meta; lg0 <- log2(b$cp10k + 1); a <- load_full_audited_bulk()
add <- rownames(a$meta)[a$meta$Cell_Line == "NCI-H2052" & a$meta$Gene_Control == "NTC" & !rownames(a$meta) %in% rownames(meta0)]
g <- intersect(rownames(lg0), rownames(a$lg)); meta1 <- rbind(meta0[, c("Sample_ID", "Cell_Line", "Gene_Control")], a$meta[add, c("Sample_ID", "Cell_Line", "Gene_Control")])
lg1 <- cbind(lg0[g, ], a$lg[g, add])[, rownames(meta1)]; cnt1 <- cbind(b$counts[g, ], a$cs[g, add])[, rownames(meta1)]
NTCS <- rownames(meta1)[meta1$Cell_Line == "NCI-H2052" & meta1$Gene_Control == "NTC"]
cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds"); top50 <- lapply(cm_genes, head, 50); CM <- names(top50)
srt <- suppressWarnings(AddModuleScore(load_screen(), features = top50, name = "Cm_", seed = 1)); sgrp <- as.character(srt$merged_call)
scm <- as.matrix(srt@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(scm) <- CM; Y <- GetAssayData(srt, layer = "data")
sb <- CreateSeuratObject(counts = as.matrix(cnt1)); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg1)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
bcm <- as.matrix(sb@meta.data[rownames(meta1), paste0("Cm_", seq_along(CM))]); colnames(bcm) <- CM
KOS <- intersect(ORDER, unique(sgrp)); KOS <- KOS[KOS %in% meta1$Gene_Control[meta1$Cell_Line == "NCI-H2052"]]; K <- length(KOS)
lv <- sapply(LINES, function(cl) matrixStats::rowVars(as.matrix(lg1[, rownames(meta1)[meta1$Cell_Line == cl], drop = FALSE]))); vord <- intersect(rownames(lg1)[order(rowMeans(lv), decreasing = TRUE)], rownames(Y))
spaces <- c("20 Cm modules", sprintf("top %d variable genes", TOPNS)); sp_tag <- c("Cm", paste0("top", TOPNS)); names(sp_tag) <- spaces
feat <- function(sp) if (sp == "20 Cm modules") list(S = scm, B = bcm) else { n <- as.integer(sub("top (\\d+) .*", "\\1", sp)); list(S = t(as.matrix(Y[vord[1:n], ])), B = t(lg1[vord[1:n], ])) }
med <- function(M, ids) apply(M[ids, , drop = FALSE], 2, median)
subsets <- unlist(lapply(seq_along(NTCS), function(k) combn(NTCS, k, simplify = FALSE)), recursive = FALSE); lab <- sapply(subsets, function(s) paste(sub("_NTC_", "", s), collapse = "+")); nn <- sapply(subsets, length)
ord <- order(nn, lab); subsets <- subsets[ord]; lab <- lab[ord]; nn <- nn[ord]
# all 8! permutations of the bulk KO labels (exact null for mean diagonal)
perms <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(seq_len(K))) })

mats <- list(); met <- list()
for (sp in spaces) { F <- feat(sp); Sd <- t(sapply(KOS, function(k) SGN * (med(F$S, which(sgrp == k)) - med(F$S, which(sgrp == "NTC")))))
  Bko <- t(sapply(KOS, function(k) med(F$B, rownames(meta1)[meta1$Cell_Line == "NCI-H2052" & meta1$Gene_Control == k])))
  for (i in seq_along(subsets)) { Bd <- SGN * sweep(Bko, 2, med(F$B, subsets[[i]]), "-")
    C <- cor(t(Sd), t(Bd), method = "spearman"); dimnames(C) <- list(screen = KOS, bulk = KOS); mats[[sp]][[lab[i]]] <- C
    dm <- mean(diag(C)); om <- mean(C[row(C) != col(C)]); pd <- apply(perms, 1, function(p) mean(C[cbind(seq_len(K), p)]))
    met[[length(met) + 1]] <- data.frame(space = sp, subset = lab[i], n_ntc = nn[i], diag_mean = dm, off_mean = om, diag_minus_off = dm - om,
       top1_rows = sum(apply(C, 1, which.max) == seq_len(K)), top1_cols = sum(apply(C, 2, which.max) == seq_len(K)), mean_row_rank = mean(sapply(seq_len(K), function(r) rank(-C[r, ])[r])), perm_p = mean(pd >= dm - 1e-12)) } }
saveRDS(mats, file.path(OUT, "spec58_KOxKO_matrices.rds")); met <- do.call(rbind, met); write.csv(met, file.path(OUT, "spec58_diagonal_metrics.csv"), row.names = FALSE)
summ <- do.call(rbind, lapply(split(met, met$subset), function(d) data.frame(subset = d$subset[1], n_ntc = d$n_ntc[1], mean_diag_minus_off = mean(d$diag_minus_off), min_over_spaces = min(d$diag_minus_off),
     mean_top1 = mean(d$top1_rows), best_perm_p = min(d$perm_p), mean_row_rank = mean(d$mean_row_rank), Cm = d$diag_minus_off[d$space == spaces[1]], top1000 = d$diag_minus_off[d$space == spaces[2]], top2000 = d$diag_minus_off[d$space == spaces[3]], top3000 = d$diag_minus_off[d$space == spaces[4]])))
summ <- summ[order(-summ$mean_diag_minus_off), ]; options(width = 230)
cat("\n=== 15 H2052 NTC combinations ranked by diagonal quality (diag_minus_off = mean(diag) - mean(off-diag), Spearman NTC-KO; mean over 4 spaces) ===\n"); print(format(summ, digits = 2), row.names = FALSE)
cat("\n=== per space: best combination by diag_minus_off, its top1 (of 8) and exact permutation p; R3-only and all-4 for reference ===\n")
for (sp in spaces) { d <- met[met$space == sp, ]; bst <- d[which.max(d$diag_minus_off), ]; cat(sprintf("%-24s best %-26s diag-off %+.3f top1 %d/8 perm_p %.4f | R3 only %+.3f (top1 %d, p %.3f) | all4 %+.3f (top1 %d, p %.3f)\n", sp, bst$subset, bst$diag_minus_off, bst$top1_rows, bst$perm_p,
   d$diag_minus_off[d$subset == "H2052R3"], d$top1_rows[d$subset == "H2052R3"], d$perm_p[d$subset == "H2052R3"], d$diag_minus_off[nn == 4][1], d$top1_rows[nn == 4][1], d$perm_p[nn == 4][1])) }

## ---- figures: one per feature space, a 4x4 grid (15 combos + a diagonal-metric bar panel)
for (sp in spaces) save_fig(paste0("fig3Q_KOxKO_matrix_", sp_tag[[sp]]), 15, 15.5, function() {
  layout(matrix(1:16, nrow = 4, byrow = TRUE)); par(oma = c(0, 0, 3.4, 0), family = "sans"); d <- met[met$space == sp, ]; best <- d$subset[which.max(d$diag_minus_off)]
  for (s in lab) { C <- mats[[sp]][[s]]; m <- d[d$subset == s, ]; par(mar = c(3.6, 3.6, 3.4, 0.5))
    image(seq_len(K), seq_len(K), t(C)[, K:1], col = DIV, zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(K)) for (j in seq_len(K)) text(j, K - i + 1, sprintf("%.1f", C[i, j]), cex = 0.5, col = INK)
    rect(seq_len(K) - 0.5, K - seq_len(K) + 0.5, seq_len(K) + 0.5, K - seq_len(K) + 1.5, border = "black", lwd = 2)
    axis(1, at = seq_len(K), labels = KOS, las = 2, tick = FALSE, cex.axis = 0.62, line = -0.6); axis(2, at = K:1, labels = KOS, las = 1, tick = FALSE, cex.axis = 0.62, line = -0.6)
    box(col = if (s == best) "#d62728" else MUTED, lwd = if (s == best) 3 else 0.6)
    mtext(s, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.62, col = if (s == best) "#d62728" else INK)
    mtext(sprintf("diag-off %+.2f | top1 %d/8 | p=%.3f", m$diag_minus_off, m$top1_rows, m$perm_p), side = 3, adj = 0, line = 0.4, cex = 0.55, col = INK2) }
  par(mar = c(3.6, 1, 3.4, 1)); plot.new(); text(0.02, 0.9, "rows = CROP-seq screen KO\ncolumns = bulk H2052 KO\nSpearman, NTC-KO\nblack squares = diagonal\n(matching KOs)\nred frame = best diagonal\n(diag - off-diag)", adj = c(0, 1), cex = 0.85)
  xs <- seq(0.05, 0.95, length.out = 42); rect(xs[-42], 0.12, xs[-1], 0.2, col = DIV[round(seq(1, 101, length.out = 41))], border = NA); text(c(0.05, 0.5, 0.95), 0.05, c("-1", "0", "+1"), cex = 0.7); text(0.5, 0.26, "Spearman rho", cex = 0.7)
  mtext(sprintf("Screen vs bulk H2052, KO x KO Spearman correlation, every H2052 NTC combination: %s", sp), outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 1.0) })

## ---- summary figure: diag - off-diag per combination x space (heatmap) + top1 counts
save_fig("fig3R_diagonal_summary", 11, 6.5, function() {
  W <- sapply(spaces, function(sp) sapply(summ$subset, function(s) met$diag_minus_off[met$space == sp & met$subset == s])); rownames(W) <- summ$subset; n <- nrow(W); vm <- max(abs(W))
  T1 <- sapply(spaces, function(sp) sapply(summ$subset, function(s) met$top1_rows[met$space == sp & met$subset == s]))
  par(mar = c(4.5, 15, 3.6, 1), family = "sans"); image(seq_len(ncol(W)), seq_len(n), t(W)[, n:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(n)) for (j in seq_len(ncol(W))) text(j, n - i + 1, sprintf("%+.2f  (%d/8)", W[i, j], T1[i, j]), cex = 0.68, col = INK)
  axis(1, at = seq_len(ncol(W)), labels = FALSE); text(seq_len(ncol(W)), 0.2, colnames(W), srt = 20, adj = 1, xpd = NA, cex = 0.8)
  text(0.35, n:1, rownames(W), adj = 1, xpd = NA, cex = 0.7, col = ifelse(rownames(W) == "H2052R3", "#eda100", ifelse(summ$n_ntc == 4, "#4a3aa7", INK)), font = ifelse(rownames(W) == "H2052R3" | summ$n_ntc == 4, 2, 1)); box(col = MUTED, lwd = 0.6)
  mtext("Diagonal quality per H2052 NTC combination (rows ranked by mean over spaces)", side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext("cell = mean(diagonal) - mean(off-diagonal) of the KO x KO Spearman matrix (screen vs bulk H2052, NTC-KO), and (number of screen KOs whose best bulk match is the same KO)", side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2) })
