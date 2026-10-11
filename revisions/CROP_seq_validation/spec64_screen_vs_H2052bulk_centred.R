# Spec 64: screen vs bulk NCI-H2052 KO x KO correlation using GROUP-CENTRED profiles (fig13c/spec11/spec62 convention) instead of NTC subtraction: for each
# platform separately, build the group (NTC + 8 KOs) x feature matrix of median expression / module score, then subtract each feature's OWN mean across that
# platform's 9 group profiles (removes the platform's own baseline -- no NTC involved). Spearman-correlate the centred screen groups against the centred bulk
# H2052 groups, per feature space (20 Cm modules, top 1000/2000/3000 variable genes). Diagonal metrics on the 8 KOs (NTC excluded from the diagonal test, kept
# as a row/column for reference): mean diag - mean off-diag, top1, mean row rank, exact permutation p (all 8! relabelings); plus after double-centering the
# resulting correlation matrix (row+column effects removed) as in spec60.
# needs: git_repo/files/scRNA_PM_modules.rds, spec_common.R, repro_common.R
# outputs: spec64_centred_KOxKO_matrices.rds, spec64_centred_diagonal_metrics.csv, fig3Z_screen_vs_H2052bulk_centred_KOxKO.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3"); TOPNS <- c(1000, 2000, 3000)
b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1)
cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds"); top50 <- lapply(cm_genes, head, 50); CM <- names(top50)
srt <- suppressWarnings(AddModuleScore(load_screen(), features = top50, name = "Cm_", seed = 1)); sgrp <- as.character(srt$merged_call)
scm <- as.matrix(srt@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(scm) <- CM; Y <- GetAssayData(srt, layer = "data")
sb <- CreateSeuratObject(counts = as.matrix(b$counts)); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
bcm <- as.matrix(sb@meta.data[rownames(meta), paste0("Cm_", seq_along(CM))]); colnames(bcm) <- CM
KOS <- intersect(ORDER, unique(sgrp)); KOS <- KOS[KOS %in% meta$Gene_Control[meta$Cell_Line == "NCI-H2052"]]; K <- length(KOS); GR <- c("NTC", KOS)
lv <- sapply(LINES, function(cl) matrixStats::rowVars(as.matrix(lg[, rownames(meta)[meta$Cell_Line == cl], drop = FALSE]))); vord <- intersect(rownames(lg)[order(rowMeans(lv), decreasing = TRUE)], rownames(Y))
spaces <- c("20 Cm modules", sprintf("top %d variable genes", TOPNS)); sp_tag <- c("Cm", paste0("top", TOPNS)); names(sp_tag) <- spaces
feat <- function(sp) if (sp == "20 Cm modules") list(S = scm, B = bcm) else { n <- as.integer(sub("top (\\d+) .*", "\\1", sp)); list(S = t(as.matrix(Y[vord[1:n], ])), B = t(lg[vord[1:n], ])) }
med <- function(M, ids) apply(M[ids, , drop = FALSE], 2, median)
group_matrix <- function(M, grp, ids_fn) t(sapply(GR, function(k) med(M, ids_fn(k))))   # groups (rows) x features
centre_cols <- function(G) G - matrix(colMeans(G), nrow(G), ncol(G), byrow = TRUE)      # subtract each FEATURE's mean across this platform's own 9 groups
h_ids <- function(k) rownames(meta)[meta$Cell_Line == "NCI-H2052" & meta$Gene_Control == k]
perms <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(seq_len(K))) })
dc <- function(C) C - rowMeans(C)[row(C)] - colMeans(C)[col(C)] + mean(C)

mats <- list(); met <- list()
for (sp in spaces) { F <- feat(sp)
  Gs <- centre_cols(group_matrix(F$S, sgrp, function(k) which(sgrp == k))); Gb <- centre_cols(group_matrix(F$B, meta$Gene_Control, h_ids))
  C <- cor(t(Gs), t(Gb), method = "spearman"); dimnames(C) <- list(screen = GR, bulk = GR); mats[[sp]] <- C; Ck <- C[KOS, KOS]; Dk <- dc(Ck)
  pv <- function(M) { dm <- mean(diag(M)); pd <- apply(perms, 1, function(p) mean(M[cbind(seq_len(K), p)])); mean(pd >= dm - 1e-12) }
  met[[length(met) + 1]] <- data.frame(space = sp, NTC_vs_NTC_rho = C["NTC", "NTC"], diag_mean = mean(diag(Ck)), off_mean = mean(Ck[row(Ck) != col(Ck)]), diag_minus_off = mean(diag(Ck)) - mean(Ck[row(Ck) != col(Ck)]),
     top1_rows = sum(apply(Ck, 1, which.max) == seq_len(K)), mean_row_rank = mean(sapply(seq_len(K), function(r) rank(-Ck[r, ])[r])), perm_p = pv(Ck), dc_diag_resid = mean(diag(Dk)), dc_top1 = sum(apply(Dk, 1, which.max) == seq_len(K)), dc_perm_p = pv(Dk)) }
met <- do.call(rbind, met); saveRDS(mats, file.path(OUT, "spec64_centred_KOxKO_matrices.rds")); write.csv(met, file.path(OUT, "spec64_centred_diagonal_metrics.csv"), row.names = FALSE)
options(width = 220); cat("\n=== screen vs bulk H2052, GROUP-CENTRED (own-platform baseline removed, no NTC subtraction) Spearman across features ===\n"); print(format(met, digits = 3), row.names = FALSE)
cat("\nrange of raw diag/off (before double-centering) vs after removing row+col effects -- same numbers, dc columns confirm no extra info beyond diag_minus_off (as in spec58/60)\n")

for (sp in spaces) { C <- mats[[sp]][KOS, KOS]; cat(sprintf("\n-- %s: own KO diag vs median-of-other-7 (per KO) --\n", sp))
  df <- data.frame(KO = KOS, own = diag(C), null_median = sapply(seq_len(K), function(i) median(C[i, -i]))); df$sep <- df$own - df$null_median; print(format(df, digits = 3), row.names = FALSE) }

save_fig("fig3Z_screen_vs_H2052bulk_centred_KOxKO", 15, 8, function() {
  layout(matrix(1:8, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.4, 0), family = "sans")
  for (sp in spaces) { C <- mats[[sp]]; m <- met[met$space == sp, ]; n <- nrow(C); par(mar = c(3.6, 3.6, 3.4, 0.5)); vm <- max(abs(C))
    image(seq_len(n), seq_len(n), t(C)[, n:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(n)) for (j in seq_len(n)) text(j, n - i + 1, sprintf("%.2f", C[i, j]), cex = 0.5, col = INK)
    rect(seq_len(n) - 0.5, n - seq_len(n) + 0.5, seq_len(n) + 0.5, n - seq_len(n) + 1.5, border = "black", lwd = 2); abline(h = n - 0.5, v = 1.5, col = "white", lwd = 3)
    axis(1, at = seq_len(n), labels = GR, las = 2, tick = FALSE, cex.axis = 0.62, line = -0.6); axis(2, at = n:1, labels = GR, las = 1, tick = FALSE, cex.axis = 0.62, line = -0.6); box(col = MUTED, lwd = 0.6)
    mtext(paste(sp, "- centred rho"), side = 3, adj = 0, line = 1.5, font = 2, cex = 0.72); mtext(sprintf("diag-off %+.3f | top1 %d/8 | perm p=%.3f", m$diag_minus_off, m$top1_rows, m$perm_p), side = 3, adj = 0, line = 0.4, cex = 0.58, col = INK2) }
  for (sp in spaces) { D <- dc(mats[[sp]][KOS, KOS]); m <- met[met$space == sp, ]; par(mar = c(3.6, 3.6, 3.4, 0.5)); vm <- max(abs(D))
    image(seq_len(K), seq_len(K), t(D)[, K:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(K)) for (j in seq_len(K)) text(j, K - i + 1, sprintf("%.2f", D[i, j]), cex = 0.5, col = INK)
    rect(seq_len(K) - 0.5, K - seq_len(K) + 0.5, seq_len(K) + 0.5, K - seq_len(K) + 1.5, border = "black", lwd = 2)
    axis(1, at = seq_len(K), labels = KOS, las = 2, tick = FALSE, cex.axis = 0.62, line = -0.6); axis(2, at = K:1, labels = KOS, las = 1, tick = FALSE, cex.axis = 0.62, line = -0.6); box(col = MUTED, lwd = 0.6)
    mtext(paste(sp, "- row/col effects removed"), side = 3, adj = 0, line = 1.5, font = 2, cex = 0.72); mtext(sprintf("diag resid %+.3f | top1 %d/8 | perm p=%.3f", m$dc_diag_resid, m$dc_top1, m$dc_perm_p), side = 3, adj = 0, line = 0.4, cex = 0.58, col = INK2) }
  mtext("Screen vs bulk H2052, KO x KO Spearman correlation on GROUP-CENTRED profiles (each platform's own baseline removed, no NTC subtraction; fig13c/spec62 convention); rows screen, columns bulk", outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 0.85)
})
