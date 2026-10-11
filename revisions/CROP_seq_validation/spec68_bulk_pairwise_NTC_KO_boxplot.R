# Spec 68: bulk-vs-bulk (no screen) agreement between all 6 pairs of the 4 cell lines, using NTC-SUBTRACTED profiles (KO - NTC median, per line -- the
# normalization spec67 found performs best for cross-platform KO-identity recovery), Spearman correlation, across the pipeline's usual feature-space sweep
# (20 Cm modules, top 1000/2000/3000 variable genes). One boxplot per TF: coloured = intra-KO (same KO, both lines; n=6 line pairs), grey = cross-KO (that
# KO's profile vs every OTHER KO in the other line; null). One-sided Wilcoxon (intra > cross) per KO per space.
# needs: git_repo/files/scRNA_PM_modules.rds, repro_common.R
# outputs: spec68_bulk_pairwise_NTC_KOxKO.csv, spec68_bulk_pairwise_NTC_KO_stats.csv, fig5C_bulk_pairwise_NTC_intra_vs_cross_boxplot.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3"); TOPNS <- c(1000, 2000, 3000)
b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1)
cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds"); top50 <- lapply(cm_genes, head, 50); CM <- names(top50)
sb <- CreateSeuratObject(counts = as.matrix(b$counts)); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
bcm <- as.matrix(sb@meta.data[rownames(meta), paste0("Cm_", seq_along(CM))]); colnames(bcm) <- CM
KOS <- ORDER; K <- length(KOS)
lv <- sapply(LINES, function(cl) matrixStats::rowVars(as.matrix(lg[, rownames(meta)[meta$Cell_Line == cl], drop = FALSE]))); vord <- rownames(lg)[order(rowMeans(lv), decreasing = TRUE)]
spaces <- c("20 Cm modules", sprintf("top %d variable genes", TOPNS))
feat <- function(sp) if (sp == "20 Cm modules") bcm else { n <- as.integer(sub("top (\\d+) .*", "\\1", sp)); t(lg[vord[1:n], ]) }   # samples x features
med <- function(M, ids) apply(M[ids, , drop = FALSE], 2, median)
ntc_profile <- function(M, cl) { mm <- meta[meta$Cell_Line == cl, ]; ntc <- rownames(mm)[mm$Gene_Control == "NTC"]
  t(sapply(KOS, function(k) { ko <- rownames(mm)[mm$Gene_Control == k]; med(M, ko) - med(M, ntc) })) }   # KO - NTC, KOs x features
pairs <- combn(LINES, 2, simplify = FALSE); plab <- sapply(pairs, function(p) paste(sub("NCI-|-211H", "", p), collapse = "-"))

rows_own <- list(); rows_cross <- list()
for (sp in spaces) { F <- feat(sp); P <- lapply(setNames(LINES, LINES), ntc_profile, M = F)
  for (pi in seq_along(pairs)) { a <- pairs[[pi]][1]; bb <- pairs[[pi]][2]
    C <- cor(t(P[[a]]), t(P[[bb]]), method = "spearman"); dimnames(C) <- list(KOS, KOS)
    for (k in KOS) rows_own[[length(rows_own) + 1]] <- data.frame(space = sp, pair = plab[pi], KO = k, rho = C[k, k])
    for (k in KOS) for (k2 in setdiff(KOS, k)) rows_cross[[length(rows_cross) + 1]] <- data.frame(space = sp, pair = plab[pi], KO = k, other_KO = k2, rho = C[k, k2]) } }
own <- do.call(rbind, rows_own); cross <- do.call(rbind, rows_cross)
write.csv(rbind(cbind(own, type = "intra-KO"), setNames(cbind(cross[, c("space", "pair", "KO", "rho")], type = "cross-KO"), c("space", "pair", "KO", "rho", "type"))),
          file.path(OUT, "spec68_bulk_pairwise_NTC_KOxKO.csv"), row.names = FALSE)

stats <- do.call(rbind, lapply(spaces, function(sp) do.call(rbind, lapply(KOS, function(k) {
  o <- own$rho[own$space == sp & own$KO == k]; cr <- cross$rho[cross$space == sp & cross$KO == k]
  wt <- suppressWarnings(wilcox.test(o, cr, alternative = "greater"))
  data.frame(space = sp, KO = k, n_intra = length(o), n_cross = length(cr), median_intra = median(o), median_cross = median(cr), sep = median(o) - median(cr), wilcox_p = wt$p.value)
}))))
write.csv(stats, file.path(OUT, "spec68_bulk_pairwise_NTC_KO_stats.csv"), row.names = FALSE)
options(width = 200); cat("\n=== bulk-vs-bulk, NTC-subtracted, Spearman: intra-KO vs cross-KO per TF per space (one-sided Wilcoxon, intra > cross) ===\n")
for (sp in spaces) { cat("\n", sp, "\n"); print(format(stats[stats$space == sp, c("KO", "n_intra", "n_cross", "median_intra", "median_cross", "sep", "wilcox_p")], digits = 3), row.names = FALSE) }
cat("\npooled (across spaces) median separation and Wilcoxon p per KO:\n")
print(do.call(rbind, lapply(KOS, function(k) { o <- own$rho[own$KO == k]; cr <- cross$rho[cross$KO == k]; wt <- suppressWarnings(wilcox.test(o, cr, alternative = "greater"))
  data.frame(KO = k, median_intra = median(o), median_cross = median(cr), sep = median(o) - median(cr), wilcox_p = wt$p.value) })))
cat(sprintf("\noverall (pooled all KOs x spaces): intra-KO median %.3f | cross-KO median %.3f | Wilcoxon p=%.3g\n", median(own$rho), median(cross$rho), wilcox.test(own$rho, cross$rho, alternative = "greater")$p.value))

save_fig("fig5C_bulk_pairwise_NTC_intra_vs_cross_boxplot", 13, 9.4, function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.4, 0), family = "sans")
  stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
  for (sp in spaces) { d_own <- own[own$space == sp, ]; d_cross <- cross[cross$space == sp, ]; s <- stats[stats$space == sp, ]
    ko_ord <- s$KO[order(-s$sep)]; n_ko <- length(ko_ord)
    par(mar = c(5, 4.4, 2.6, 0.8)); plot(NA, xlim = c(0.5, n_ko + 0.5), ylim = c(-1, 1.15), xaxt = "n", xlab = "", ylab = "Spearman rho (bulk line vs bulk line, KO-NTC diff)")
    axis(1, at = seq_len(n_ko), labels = FALSE); text(seq_len(n_ko), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), ko_ord, srt = 45, adj = 1, xpd = NA, cex = 0.85)
    for (i in seq_len(n_ko)) boxplot(d_cross$rho[d_cross$KO == ko_ord[i]], at = i - 0.18, add = TRUE, boxwex = 0.3, col = "grey85", border = "grey40", outline = FALSE, axes = FALSE)
    for (i in seq_len(n_ko)) boxplot(d_own$rho[d_own$KO == ko_ord[i]], at = i + 0.18, add = TRUE, boxwex = 0.3, col = adjustcolor(KO_COL[ko_ord[i]], 0.5), border = INK2, outline = FALSE, axes = FALSE)
    set.seed(1); for (i in seq_len(n_ko)) { dd <- d_own[d_own$KO == ko_ord[i], ]; points(jitter(rep(i + 0.18, nrow(dd)), amount = 0.08), dd$rho, pch = 16, cex = 1, col = INK2) }
    abline(h = 0, col = INK, lty = 2)
    for (i in seq_len(n_ko)) { p <- s$wilcox_p[s$KO == ko_ord[i]]; text(i, 1.08, stars(p), cex = 0.8, font = ifelse(p < 0.05, 2, 1)) }
    if (sp == spaces[1]) legend("bottomleft", legend = c("intra-KO (same KO, both lines)", "cross-KO (null)"), col = c("grey30", "grey60"), pch = 15, bty = "n", cex = 0.65)
    mtext(sp, side = 3, adj = 0, line = 0.3, font = 2, cex = 0.85) }
  mtext("Bulk-vs-bulk (all 6 line pairs), NTC-subtracted (KO-NTC) profiles, Spearman: intra-KO vs cross-KO per TF", outer = TRUE, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext("coloured = intra-KO (dots = the 6 line pairs); grey = same KO's profile vs every other bulk KO in the other line (null); *** p<0.001 ** p<0.01 * p<0.05 ns not sig. (one-sided Wilcoxon, intra > cross)", outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.55, col = INK2)
})
