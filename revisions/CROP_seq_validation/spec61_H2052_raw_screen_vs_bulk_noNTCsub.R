# Spec 61: screen vs bulk NCI-H2052 correlation WITHOUT NTC subtraction: per-feature median expression (log-norm. in the screen, log2 CP10K in bulk; AddModuleScore for the
# 20 Cm modules) of each KO group is correlated directly (Spearman) across features. KO x KO matrix (rows = screen KO, columns = bulk H2052 KO, NTC row/column included as
# a baseline reference; NTC baseline = all 4 H2052 NTC samples, or R3 only for the alternative), in 4 feature spaces (20 Cm modules, top 1000/2000/3000 variable genes).
# Diagonal metrics on the 8 KOs (NTC excluded): mean diagonal - mean off-diagonal, top1, mean row rank, exact permutation p; plus the same after double-centering
# (row+column effects removed) and per-KO own-vs-null (screen KO vs its own bulk KO minus median over the other 7 bulk KOs).
# needs: git_repo/files/scRNA_PM_modules.rds, spec_common.R, repro_common.R
# outputs: spec61_raw_KOxKO_matrices.rds, spec61_raw_diagonal_metrics.csv, spec61_raw_per_KO_own_vs_null.csv, fig3W_raw_KOxKO_matrix_noNTCsub.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3"); TOPNS <- c(1000, 2000, 3000)
b <- load_bulk(); meta0 <- b$meta; lg0 <- log2(b$cp10k + 1); a <- load_full_audited_bulk()
add <- rownames(a$meta)[a$meta$Cell_Line == "NCI-H2052" & a$meta$Gene_Control == "NTC" & !rownames(a$meta) %in% rownames(meta0)]
g <- intersect(rownames(lg0), rownames(a$lg)); meta1 <- rbind(meta0[, c("Sample_ID", "Cell_Line", "Gene_Control")], a$meta[add, c("Sample_ID", "Cell_Line", "Gene_Control")])
lg1 <- cbind(lg0[g, ], a$lg[g, add])[, rownames(meta1)]; cnt1 <- cbind(b$counts[g, ], a$cs[g, add])[, rownames(meta1)]
NTC4 <- rownames(meta1)[meta1$Cell_Line == "NCI-H2052" & meta1$Gene_Control == "NTC"]
cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds"); top50 <- lapply(cm_genes, head, 50); CM <- names(top50)
srt <- suppressWarnings(AddModuleScore(load_screen(), features = top50, name = "Cm_", seed = 1)); sgrp <- as.character(srt$merged_call)
scm <- as.matrix(srt@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(scm) <- CM; Y <- GetAssayData(srt, layer = "data")
sb <- CreateSeuratObject(counts = as.matrix(cnt1)); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg1)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
bcm <- as.matrix(sb@meta.data[rownames(meta1), paste0("Cm_", seq_along(CM))]); colnames(bcm) <- CM
KOS <- intersect(ORDER, unique(sgrp)); KOS <- KOS[KOS %in% meta1$Gene_Control[meta1$Cell_Line == "NCI-H2052"]]; K <- length(KOS); GR <- c("NTC", KOS)
lv <- sapply(LINES, function(cl) matrixStats::rowVars(as.matrix(lg1[, rownames(meta1)[meta1$Cell_Line == cl], drop = FALSE]))); vord <- intersect(rownames(lg1)[order(rowMeans(lv), decreasing = TRUE)], rownames(Y))
spaces <- c("20 Cm modules", sprintf("top %d variable genes", TOPNS))
feat <- function(sp) if (sp == "20 Cm modules") list(S = scm, B = bcm) else { n <- as.integer(sub("top (\\d+) .*", "\\1", sp)); list(S = t(as.matrix(Y[vord[1:n], ])), B = t(lg1[vord[1:n], ])) }
med <- function(M, ids) apply(M[ids, , drop = FALSE], 2, median)
perms <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(seq_len(K))) })
dc <- function(C) C - rowMeans(C)[row(C)] - colMeans(C)[col(C)] + mean(C)
h_ids <- function(k) if (k == "NTC") NTC4 else rownames(meta1)[meta1$Cell_Line == "NCI-H2052" & meta1$Gene_Control == k]

mats <- list(); met <- list(); pko <- list()
for (sp in spaces) { F <- feat(sp)
  S <- t(sapply(GR, function(k) med(F$S, which(sgrp == k)))); B <- t(sapply(GR, function(k) med(F$B, h_ids(k))))
  C <- cor(t(S), t(B), method = "spearman"); dimnames(C) <- list(screen = GR, bulk = GR); mats[[sp]] <- C; Ck <- C[KOS, KOS]; Dk <- dc(Ck)
  pv <- function(M) { dm <- mean(diag(M)); pd <- apply(perms, 1, function(p) mean(M[cbind(seq_len(K), p)])); mean(pd >= dm - 1e-12) }
  met[[length(met) + 1]] <- data.frame(space = sp, NTC_vs_NTC_rho = C["NTC", "NTC"], diag_mean = mean(diag(Ck)), off_mean = mean(Ck[row(Ck) != col(Ck)]), diag_minus_off = mean(diag(Ck)) - mean(Ck[row(Ck) != col(Ck)]),
     top1_rows = sum(apply(Ck, 1, which.max) == seq_len(K)), mean_row_rank = mean(sapply(seq_len(K), function(r) rank(-Ck[r, ])[r])), perm_p = pv(Ck), dc_diag_resid = mean(diag(Dk)), dc_top1 = sum(apply(Dk, 1, which.max) == seq_len(K)), dc_perm_p = pv(Dk))
  pko[[length(pko) + 1]] <- data.frame(space = sp, KO = KOS, own_rho = diag(Ck), null_median = sapply(seq_len(K), function(i) median(Ck[i, -i])), null_max = sapply(seq_len(K), function(i) max(Ck[i, -i])), rank_of_own_in_row = sapply(seq_len(K), function(i) rank(-Ck[i, ])[i])) }
met <- do.call(rbind, met); pko <- do.call(rbind, pko); pko$sep <- pko$own_rho - pko$null_median
saveRDS(mats, file.path(OUT, "spec61_raw_KOxKO_matrices.rds")); write.csv(met, file.path(OUT, "spec61_raw_diagonal_metrics.csv"), row.names = FALSE); write.csv(pko, file.path(OUT, "spec61_raw_per_KO_own_vs_null.csv"), row.names = FALSE)
options(width = 220); cat("\n=== RAW (no NTC subtraction) screen vs bulk H2052, Spearman across features ===\n"); print(format(met, digits = 2), row.names = FALSE)
cat("\nper-KO: own rho | median rho vs the 7 other bulk KOs | separation | rank of own among 8 (1 = best)\n")
for (sp in spaces) { cat(" ", sp, "\n"); print(format(pko[pko$space == sp, c("KO", "own_rho", "null_median", "sep", "rank_of_own_in_row")], digits = 3), row.names = FALSE) }
cat("\nrange of all 8x8 raw correlations per space (min, median, max):\n"); for (sp in spaces) { x <- mats[[sp]][KOS, KOS]; cat(sprintf("  %-24s %.2f, %.2f, %.2f\n", sp, min(x), median(x), max(x))) }

save_fig("fig3W_raw_KOxKO_matrix_noNTCsub", 15, 8, function() {
  layout(matrix(1:8, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.4, 0), family = "sans")
  for (sp in spaces) { C <- mats[[sp]]; m <- met[met$space == sp, ]; n <- nrow(C); par(mar = c(3.6, 3.6, 3.4, 0.5)); lo <- min(C); hi <- max(C)
    image(seq_len(n), seq_len(n), t(C)[, n:1], col = colorRampPalette(c("#f0efec", "#eb6834"))(100), zlim = c(min(0, lo), hi), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(n)) for (j in seq_len(n)) text(j, n - i + 1, sprintf("%.2f", C[i, j]), cex = 0.5, col = INK)
    rect(seq_len(n) - 0.5, n - seq_len(n) + 0.5, seq_len(n) + 0.5, n - seq_len(n) + 1.5, border = "black", lwd = 2); abline(h = n - 0.5, v = 1.5, col = "white", lwd = 3)
    axis(1, at = seq_len(n), labels = GR, las = 2, tick = FALSE, cex.axis = 0.62, line = -0.6); axis(2, at = n:1, labels = GR, las = 1, tick = FALSE, cex.axis = 0.62, line = -0.6); box(col = MUTED, lwd = 0.6)
    mtext(paste(sp, "- raw rho"), side = 3, adj = 0, line = 1.5, font = 2, cex = 0.72); mtext(sprintf("diag-off %+.3f | top1 %d/8 | perm p=%.3f", m$diag_minus_off, m$top1_rows, m$perm_p), side = 3, adj = 0, line = 0.4, cex = 0.58, col = INK2) }
  for (sp in spaces) { D <- dc(mats[[sp]][KOS, KOS]); m <- met[met$space == sp, ]; par(mar = c(3.6, 3.6, 3.4, 0.5)); vm <- max(abs(D))
    image(seq_len(K), seq_len(K), t(D)[, K:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(K)) for (j in seq_len(K)) text(j, K - i + 1, sprintf("%.2f", D[i, j]), cex = 0.5, col = INK)
    rect(seq_len(K) - 0.5, K - seq_len(K) + 0.5, seq_len(K) + 0.5, K - seq_len(K) + 1.5, border = "black", lwd = 2)
    axis(1, at = seq_len(K), labels = KOS, las = 2, tick = FALSE, cex.axis = 0.62, line = -0.6); axis(2, at = K:1, labels = KOS, las = 1, tick = FALSE, cex.axis = 0.62, line = -0.6); box(col = MUTED, lwd = 0.6)
    mtext(paste(sp, "- row/col effects removed"), side = 3, adj = 0, line = 1.5, font = 2, cex = 0.72); mtext(sprintf("diag resid %+.3f | top1 %d/8 | perm p=%.3f", m$dc_diag_resid, m$dc_top1, m$dc_perm_p), side = 3, adj = 0, line = 0.4, cex = 0.58, col = INK2) }
  mtext("Screen vs bulk H2052, KO x KO Spearman correlation WITHOUT NTC subtraction (per-feature median expression / module score); rows screen, columns bulk; NTC = all 4 H2052 NTC samples", outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 0.9) })
