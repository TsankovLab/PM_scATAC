# Spec 65: extends spec64 (group-centred profiles, fig13c/spec62 convention -- each platform's own baseline removed by subtracting each feature's mean across
# that platform's own 9 group profiles, NO NTC subtraction) to ALL 4 bulk lines, not just NCI-H2052. Screen vs bulk-line KO x KO Spearman correlation, per line,
# per feature space (20 Cm modules, top 1000/2000/3000 variable genes). Diagonal metrics (mean diag - mean off-diag, top1, mean row rank, exact permutation p)
# per line x space, plus a combined summary and per-KO x line separation table.
# needs: git_repo/files/scRNA_PM_modules.rds, spec_common.R, repro_common.R
# outputs: spec65_centred_KOxKO_matrices_allbulk.rds, spec65_centred_diagonal_metrics_allbulk.csv, spec65_centred_per_KO_separation_allbulk.csv,
#          fig4G_screen_vs_allbulk_centred_KOxKO_{line}.{png,pdf} (one file per line), fig4H_centred_diagonal_summary_allbulk.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3"); TOPNS <- c(1000, 2000, 3000)
b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1)
cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds"); top50 <- lapply(cm_genes, head, 50); CM <- names(top50)
srt <- suppressWarnings(AddModuleScore(load_screen(), features = top50, name = "Cm_", seed = 1)); sgrp <- as.character(srt$merged_call)
scm <- as.matrix(srt@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(scm) <- CM; Y <- GetAssayData(srt, layer = "data")
sb <- CreateSeuratObject(counts = as.matrix(b$counts)); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
bcm <- as.matrix(sb@meta.data[rownames(meta), paste0("Cm_", seq_along(CM))]); colnames(bcm) <- CM
KOS <- intersect(ORDER, unique(sgrp)); K <- length(KOS); GR <- c("NTC", KOS)
lv <- sapply(LINES, function(cl) matrixStats::rowVars(as.matrix(lg[, rownames(meta)[meta$Cell_Line == cl], drop = FALSE]))); vord <- intersect(rownames(lg)[order(rowMeans(lv), decreasing = TRUE)], rownames(Y))
spaces <- c("20 Cm modules", sprintf("top %d variable genes", TOPNS)); sp_tag <- c("Cm", paste0("top", TOPNS)); names(sp_tag) <- spaces
feat <- function(sp) if (sp == "20 Cm modules") list(S = scm, B = bcm) else { n <- as.integer(sub("top (\\d+) .*", "\\1", sp)); list(S = t(as.matrix(Y[vord[1:n], ])), B = t(lg[vord[1:n], ])) }
med <- function(M, ids) apply(M[ids, , drop = FALSE], 2, median)
group_matrix <- function(M, ids_fn) t(sapply(GR, function(k) med(M, ids_fn(k))))
centre_cols <- function(G) G - matrix(colMeans(G), nrow(G), ncol(G), byrow = TRUE)
perms <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(seq_len(K))) })

mats <- list(); met <- list(); pko <- list()
for (cl in LINES) { ids_fn <- function(k) rownames(meta)[meta$Cell_Line == cl & meta$Gene_Control == k]
  ko_cl <- KOS[sapply(KOS, function(k) length(ids_fn(k)) > 0)]; ntc_ok <- length(ids_fn("NTC")) > 0
  for (sp in spaces) { F <- feat(sp)
    Gs <- centre_cols(group_matrix(F$S, function(k) which(sgrp == k)))
    Gb <- centre_cols(group_matrix(F$B, ids_fn))
    C <- cor(t(Gs), t(Gb), method = "spearman"); dimnames(C) <- list(screen = GR, bulk = GR); mats[[cl]][[sp]] <- C
    Ck <- C[KOS, KOS]; dm <- mean(diag(Ck)); om <- mean(Ck[row(Ck) != col(Ck)])
    pd <- apply(perms, 1, function(p) mean(Ck[cbind(seq_len(K), p)])); pp <- mean(pd >= dm - 1e-12)
    met[[length(met) + 1]] <- data.frame(line = cl, space = sp, NTC_vs_NTC_rho = if (ntc_ok) C["NTC", "NTC"] else NA, diag_mean = dm, off_mean = om, diag_minus_off = dm - om,
       top1_rows = sum(apply(Ck, 1, which.max) == seq_len(K)), mean_row_rank = mean(sapply(seq_len(K), function(r) rank(-Ck[r, ])[r])), perm_p = pp)
    pko[[length(pko) + 1]] <- data.frame(line = cl, space = sp, KO = KOS, own = diag(Ck), null_median = sapply(seq_len(K), function(i) median(Ck[i, -i])),
                                          sep = diag(Ck) - sapply(seq_len(K), function(i) median(Ck[i, -i]))) } }
met <- do.call(rbind, met); pko <- do.call(rbind, pko)
saveRDS(mats, file.path(OUT, "spec65_centred_KOxKO_matrices_allbulk.rds")); write.csv(met, file.path(OUT, "spec65_centred_diagonal_metrics_allbulk.csv"), row.names = FALSE); write.csv(pko, file.path(OUT, "spec65_centred_per_KO_separation_allbulk.csv"), row.names = FALSE)
options(width = 220); cat("\n=== screen vs EACH bulk line, GROUP-CENTRED Spearman (own-platform baseline removed, no NTC subtraction) ===\n"); print(format(met, digits = 3), row.names = FALSE)
cat("\nmean diag-off per line (across 4 spaces):\n"); print(round(sort(tapply(met$diag_minus_off, met$line, mean), decreasing = TRUE), 3))
cat("\nmean separation per KO (across 4 lines x 4 spaces):\n"); print(round(sort(tapply(pko$sep, pko$KO, mean), decreasing = TRUE), 3))
cat("\nper-KO separation, mean over the 4 feature spaces, by line:\n"); print(round(tapply(pko$sep, list(pko$KO, pko$line), mean)[names(sort(tapply(pko$sep, pko$KO, mean), decreasing = TRUE)), ], 3))

for (cl in LINES) save_fig(paste0("fig4G_screen_vs_allbulk_centred_KOxKO_", sub("NCI-|-211H", "", cl)), 15, 4.2, function() {
  layout(matrix(1:4, nrow = 1)); par(oma = c(0, 0, 3.2, 0), family = "sans")
  for (sp in spaces) { C <- mats[[cl]][[sp]]; m <- met[met$line == cl & met$space == sp, ]; n <- nrow(C); par(mar = c(3.6, 3.6, 3.2, 0.5)); vm <- max(abs(C))
    image(seq_len(n), seq_len(n), t(C)[, n:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(n)) for (j in seq_len(n)) text(j, n - i + 1, sprintf("%.2f", C[i, j]), cex = 0.5, col = INK)
    rect(seq_len(n) - 0.5, n - seq_len(n) + 0.5, seq_len(n) + 0.5, n - seq_len(n) + 1.5, border = "black", lwd = 2); abline(h = n - 0.5, v = 1.5, col = "white", lwd = 3)
    axis(1, at = seq_len(n), labels = GR, las = 2, tick = FALSE, cex.axis = 0.62, line = -0.6); axis(2, at = n:1, labels = GR, las = 1, tick = FALSE, cex.axis = 0.62, line = -0.6); box(col = MUTED, lwd = 0.6)
    mtext(sp, side = 3, adj = 0, line = 1.4, font = 2, cex = 0.75); mtext(sprintf("diag-off %+.3f | top1 %d/8 | p=%.3f", m$diag_minus_off, m$top1_rows, m$perm_p), side = 3, adj = 0, line = 0.4, cex = 0.58, col = INK2) }
  mtext(sprintf("Screen vs bulk %s, KO x KO Spearman on group-centred profiles (own-platform baseline removed, no NTC subtraction); rows screen, columns bulk", cl), outer = TRUE, side = 3, adj = 0, line = 1.7, font = 2, cex = 0.9)
})

save_fig("fig4H_centred_diagonal_summary_allbulk", 10, 6.2, function() {
  W <- sapply(spaces, function(sp) sapply(LINES, function(cl) met$diag_minus_off[met$line == cl & met$space == sp])); rownames(W) <- sub("NCI-|-211H", "", LINES); n <- nrow(W); vm <- max(abs(W))
  T1 <- sapply(spaces, function(sp) sapply(LINES, function(cl) met$top1_rows[met$line == cl & met$space == sp]))
  par(mar = c(4.5, 6, 3.6, 1), family = "sans"); image(seq_len(ncol(W)), seq_len(n), t(W)[, n:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(n)) for (j in seq_len(ncol(W))) text(j, n - i + 1, sprintf("%+.2f  (%d/8)", W[i, j], T1[i, j]), cex = 0.75, col = INK)
  axis(1, at = seq_len(ncol(W)), labels = FALSE); text(seq_len(ncol(W)), 0.3, colnames(W), srt = 20, adj = 1, xpd = NA, cex = 0.85)
  axis(2, at = n:1, labels = rownames(W), las = 1, tick = FALSE, cex.axis = 0.9, font = 2); box(col = MUTED, lwd = 0.6)
  mtext("Screen vs bulk (group-centred, no NTC subtraction): diagonal quality per line x feature space", side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext("cell = mean(diagonal) - mean(off-diagonal) of the KO x KO Spearman matrix, and (number of screen KOs whose best bulk match is the same KO)", side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2)
})
