# Spec 66: bulk-vs-bulk (no screen) agreement between all 6 pairs of the 4 cell lines, using the SAME group-centred profiles as spec65 (each line's own
# baseline removed by subtracting each feature's mean across that line's own 9 group profiles, NO NTC subtraction). For every line pair, KO x KO Spearman
# correlation matrix (20 Cm modules, top 1000/2000/3000 variable genes); per TF, own-KO (same KO in both lines) vs cross-KO (different KOs) correlation,
# shown as boxplots -- this is the bulk-only ceiling to compare against the screen-vs-bulk results (spec64/65).
# needs: git_repo/files/scRNA_PM_modules.rds, repro_common.R
# outputs: spec66_bulk_pairwise_KOxKO.csv, spec66_bulk_pairwise_diagonal_metrics.csv,
#          fig4I_bulk_pairwise_own_vs_cross_boxplot.{png,pdf}, fig4J_bulk_pairwise_diagonal_summary.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3"); TOPNS <- c(1000, 2000, 3000)
b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1)
cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds"); top50 <- lapply(cm_genes, head, 50); CM <- names(top50)
sb <- CreateSeuratObject(counts = as.matrix(b$counts)); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
bcm <- as.matrix(sb@meta.data[rownames(meta), paste0("Cm_", seq_along(CM))]); colnames(bcm) <- CM
KOS <- ORDER; K <- length(KOS); GR <- c("NTC", KOS)
lv <- sapply(LINES, function(cl) matrixStats::rowVars(as.matrix(lg[, rownames(meta)[meta$Cell_Line == cl], drop = FALSE]))); vord <- rownames(lg)[order(rowMeans(lv), decreasing = TRUE)]
spaces <- c("20 Cm modules", sprintf("top %d variable genes", TOPNS))
feat <- function(sp) if (sp == "20 Cm modules") bcm else { n <- as.integer(sub("top (\\d+) .*", "\\1", sp)); t(lg[vord[1:n], ]) }   # samples x features, matching med()'s M[ids, ] convention
med <- function(M, ids) apply(M[ids, , drop = FALSE], 2, median)
group_matrix <- function(M, cl) t(sapply(GR, function(k) med(M, rownames(meta)[meta$Cell_Line == cl & meta$Gene_Control == k])))
centre_cols <- function(G) G - matrix(colMeans(G), nrow(G), ncol(G), byrow = TRUE)
pairs <- combn(LINES, 2, simplify = FALSE); plab <- sapply(pairs, function(p) paste(sub("NCI-|-211H", "", p), collapse = "-"))

mats <- list(); rows_own <- list(); rows_cross <- list(); met <- list()
for (sp in spaces) { F <- feat(sp); G <- lapply(setNames(LINES, LINES), function(cl) centre_cols(group_matrix(F, cl)))
  for (pi in seq_along(pairs)) { a <- pairs[[pi]][1]; bb <- pairs[[pi]][2]
    C <- cor(t(G[[a]]), t(G[[bb]]), method = "spearman"); dimnames(C) <- list(GR, GR); mats[[sp]][[plab[pi]]] <- C
    Ck <- C[KOS, KOS]; dm <- mean(diag(Ck)); om <- mean(Ck[row(Ck) != col(Ck)])
    met[[length(met) + 1]] <- data.frame(pair = plab[pi], line_a = a, line_b = bb, space = sp, diag_mean = dm, off_mean = om, diag_minus_off = dm - om, top1 = sum(apply(Ck, 1, which.max) == seq_len(K)))
    for (k in KOS) rows_own[[length(rows_own) + 1]] <- data.frame(space = sp, pair = plab[pi], KO = k, rho = Ck[k, k])
    for (k in KOS) for (k2 in setdiff(KOS, k)) rows_cross[[length(rows_cross) + 1]] <- data.frame(space = sp, pair = plab[pi], KO_a = k, KO_b = k2, rho = Ck[k, k2]) } }
own <- do.call(rbind, rows_own); cross <- do.call(rbind, rows_cross); met <- do.call(rbind, met)
write.csv(rbind(cbind(own, type = "own"), setNames(cbind(cross[, c("space", "pair", "KO_a", "rho")], type = "cross"), c("space", "pair", "KO", "rho", "type"))), file.path(OUT, "spec66_bulk_pairwise_KOxKO.csv"), row.names = FALSE)
write.csv(met, file.path(OUT, "spec66_bulk_pairwise_diagonal_metrics.csv"), row.names = FALSE)
options(width = 200); cat("\n=== bulk-vs-bulk (group-centred, no NTC), diagonal quality per line pair x space ===\n"); print(format(met[, c("pair", "space", "diag_mean", "off_mean", "diag_minus_off", "top1")], digits = 3), row.names = FALSE)
cat("\nmean diag-off per line pair (across 4 spaces):\n"); print(round(sort(tapply(met$diag_minus_off, met$pair, mean), decreasing = TRUE), 3))
cat("\nown-KO vs cross-KO median rho per TF (pooled over line pairs x spaces):\n")
print(do.call(rbind, lapply(KOS, function(k) data.frame(KO = k, own_median = median(own$rho[own$KO == k]), cross_median = median(cross$rho[cross$KO_a == k]), sep = median(own$rho[own$KO == k]) - median(cross$rho[cross$KO_a == k])))))

save_fig("fig4I_bulk_pairwise_own_vs_cross_boxplot", 13, 9.4, function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.4, 0), family = "sans")
  for (sp in spaces) { d_own <- own[own$space == sp, ]; d_cross <- cross[cross$space == sp, ]
    ko_ord <- names(sort(sapply(split(d_own$rho, d_own$KO), median), decreasing = TRUE)); n_ko <- length(ko_ord)
    par(mar = c(5, 4.4, 2.6, 0.8)); plot(NA, xlim = c(0.5, n_ko + 0.5), ylim = c(-1, 1), xaxt = "n", xlab = "", ylab = "Spearman rho (bulk line vs bulk line)")
    axis(1, at = seq_len(n_ko), labels = FALSE); text(seq_len(n_ko), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), ko_ord, srt = 45, adj = 1, xpd = NA, cex = 0.85)
    for (i in seq_len(n_ko)) boxplot(d_cross$rho[d_cross$KO_a == ko_ord[i]], at = i - 0.18, add = TRUE, boxwex = 0.3, col = "grey85", border = "grey40", outline = FALSE, axes = FALSE)
    for (i in seq_len(n_ko)) boxplot(d_own$rho[d_own$KO == ko_ord[i]], at = i + 0.18, add = TRUE, boxwex = 0.3, col = adjustcolor(KO_COL[ko_ord[i]], 0.5), border = INK2, outline = FALSE, axes = FALSE)
    set.seed(1); for (i in seq_len(n_ko)) { dd <- d_own[d_own$KO == ko_ord[i], ]; points(jitter(rep(i + 0.18, nrow(dd)), amount = 0.08), dd$rho, pch = 16, cex = 1, col = INK2) }
    abline(h = 0, col = INK, lty = 2); if (sp == spaces[1]) legend("bottomleft", legend = c("own KO (same KO, both lines)", "cross-KO (null)"), col = c("grey30", "grey60"), pch = 15, bty = "n", cex = 0.65)
    mtext(sp, side = 3, adj = 0, line = 0.3, font = 2, cex = 0.85) }
  mtext("Bulk-vs-bulk (all 6 line pairs), group-centred KO profiles, no NTC subtraction: own-KO vs cross-KO Spearman rho per TF", outer = TRUE, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext("coloured = own KO (dots = the 6 line pairs); grey = same screen KO's profile vs every other bulk KO in that line pair (null)", outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.58, col = INK2)
})

save_fig("fig4J_bulk_pairwise_diagonal_summary", 9.5, 5.6, function() {
  W <- sapply(spaces, function(sp) sapply(plab, function(p) met$diag_minus_off[met$pair == p & met$space == sp])); rownames(W) <- plab; n <- nrow(W); vm <- max(abs(W))
  T1 <- sapply(spaces, function(sp) sapply(plab, function(p) met$top1[met$pair == p & met$space == sp]))
  par(mar = c(4.5, 8, 3.6, 1), family = "sans"); image(seq_len(ncol(W)), seq_len(n), t(W)[, n:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(n)) for (j in seq_len(ncol(W))) text(j, n - i + 1, sprintf("%+.2f  (%d/8)", W[i, j], T1[i, j]), cex = 0.72, col = INK)
  axis(1, at = seq_len(ncol(W)), labels = FALSE); text(seq_len(ncol(W)), 0.3, sub("top ", "top", sub(" variable genes", "vg", spaces)), srt = 20, adj = 1, xpd = NA, cex = 0.8)
  axis(2, at = n:1, labels = rownames(W), las = 1, tick = FALSE, cex.axis = 0.85, font = 2); box(col = MUTED, lwd = 0.6)
  mtext("Bulk-vs-bulk diagonal quality, all 6 line pairs (group-centred, no NTC subtraction)", side = 3, adj = 0, line = 1.9, font = 2, cex = 0.9)
  mtext("cell = mean(diagonal) - mean(off-diagonal) of the KO x KO Spearman matrix; (n) = KOs whose best match in the other line is the same KO", side = 3, adj = 0, line = 0.7, cex = 0.55, col = INK2)
})
