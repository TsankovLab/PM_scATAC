# Spec 17: PCA of the CROP-seq guides (pseudobulk). One profile per SINGLE-guide call = mean log-normalised expression of its cells (guides with >= 10 cells, so every KO of the screen is present;
# multi-guide calls are excluded; NTC_1 and NTC_3 kept separate), over the screen's top variable genes. Two versions: raw, and after regressing the cell-cycle index (cc = S.Score + G2M.Score) out per gene, pooled over all cells used.
# Genes are centred and scaled across guides before prcomp. PCs are related to the guide's proliferating fraction (Seurat Phase S/G2M) and mean cc.
# needs: screen_export_figs.R (variable gene list)
# outputs: spec17_guide_PCs.csv, spec17_PC_correlations.csv, fig21_guide_PCA.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
srt <- load_screen(); md <- srt@meta.data; Yall <- GetAssayData(srt, layer = "data")
vg <- readLines(file.path(OUT, "screen_top5000_variable_genes.txt")); vg <- vg[vg %in% rownames(Yall)]; cat("variable genes:", length(vg), "\n")
tb <- table(md$crispr_calls); guides <- sort(names(tb)[tb >= 10], method = "radix")
gene_of <- function(x) sapply(strsplit(x, "\\|"), function(p) { q <- p[!grepl("^NTC", p)]; if (!length(q)) "NTC" else sub("_.*", "", q[1]) }); guides <- guides[gene_of(guides) %in% c(ORDER, EXTRA, "NTC") & !grepl("\\|", guides)]      # single-guide calls only
cells <- which(md$crispr_calls %in% guides); lab <- md$crispr_calls[cells]; Y <- Yall[vg, cells]; ccv <- md$cc[cells]; ph <- md$Phase[cells] != "G1"
G <- sparseMatrix(i = seq_along(lab), j = match(lab, guides), x = 1, dims = c(length(lab), length(guides)), dimnames = list(NULL, guides)); cnt <- colSums(G)
M <- as.matrix(Y %*% G) %*% diag(1 / cnt); colnames(M) <- guides; mc <- as.numeric(crossprod(G, ccv)) / cnt
b <- (as.numeric(Y %*% (ccv - mean(ccv))) / (length(ccv) - 1)) / var(ccv); Mcc <- M - outer(b, mc - mean(ccv))
info <- data.frame(guide = guides, KO = gene_of(guides), type = ifelse(grepl("\\|", guides), "multi-guide", "single guide"), cells = as.integer(cnt), prolif_frac = as.numeric(crossprod(G, ph)) / cnt, mean_cc = mc); print(format(info, digits = 3), row.names = FALSE)
pca <- function(X) { keep <- apply(X, 1, sd) > 0; p <- prcomp(t(X[keep, ]), center = TRUE, scale. = TRUE); list(p = p, ve = p$sdev^2 / sum(p$sdev^2)) }
P <- list(raw = pca(M), cc = pca(Mcc))
out <- do.call(rbind, lapply(names(P), function(v) data.frame(version = v, info, P[[v]]$p$x[, 1:5]))); write.csv(out, file.path(OUT, "spec17_guide_PCs.csv"), row.names = FALSE)
pc_cor <- do.call(rbind, lapply(names(P), function(v) do.call(rbind, lapply(1:5, function(k) { x <- P[[v]]$p$x[, k]; data.frame(version = v, PC = k, var_explained = P[[v]]$ve[k], rho_prolif_frac = cor(x, info$prolif_frac, method = "spearman"), rho_mean_cc = cor(x, info$mean_cc, method = "spearman")) }))))
write.csv(pc_cor, file.path(OUT, "spec17_PC_correlations.csv"), row.names = FALSE); options(width = 200); cat("\n"); print(format(pc_cor, digits = 3), row.names = FALSE)
# is PC separation by KO or by guide? mean distance in PC1-3 space: same-KO guide pairs vs different-KO pairs
for (v in names(P)) { X <- P[[v]]$p$x[, 1:3]; d <- as.matrix(dist(X)); ut <- which(upper.tri(d), arr.ind = TRUE); same <- info$KO[ut[, 1]] == info$KO[ut[, 2]]
  cat(sprintf("%s: mean PC1-3 distance same-KO guide pairs = %.1f (n=%d) vs different-KO pairs = %.1f (n=%d)\n", v, mean(d[ut][same]), sum(same), mean(d[ut][!same]), sum(!same))) }

## ---------------------------------------------------------------- figure
kcol <- KO_COL
save_fig("fig21_guide_PCA", 13, 11, function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4.5, 0.5), mar = c(4.4, 4.6, 3.4, 1), family = "sans")
  for (v in names(P)) for (pr in list(c(1, 2), c(2, 3))) { p <- P[[v]]; X <- p$p$x[, pr]; rg <- apply(X, 2, range); pad <- 0.15 * apply(rg, 2, diff)
    plot(X[, 1], X[, 2], xlim = rg[, 1] + c(-1, 1) * pad[1], ylim = rg[, 2] + c(-1, 1) * pad[2], pch = ifelse(info$type == "single guide", 21, 24), bg = ifelse(info$type == "single guide", adjustcolor(kcol[info$KO], 0.9), adjustcolor(kcol[info$KO], 0.3)), col = ifelse(info$type == "single guide", "white", kcol[info$KO]), lwd = 1.3, cex = 0.9 + 2.4 * sqrt(info$cells / max(info$cells)), bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.85,
         xlab = sprintf("PC%d (%.0f%% of variance)", pr[1], 100 * p$ve[pr[1]]), ylab = sprintf("PC%d (%.0f%%)", pr[2], 100 * p$ve[pr[2]])); abline(h = 0, v = 0, col = GRID)
    text(X[, 1], X[, 2], gsub("_", "-", info$guide), pos = 3, cex = 0.6, offset = 0.6)
    rho <- sapply(pr, function(k) cor(p$p$x[, k], info$prolif_frac, method = "spearman"))
    mtext(sprintf("%s, PC%d vs PC%d", ifelse(v == "raw", "raw pseudobulk", "cell-cycle regressed"), pr[1], pr[2]), side = 3, adj = 0, line = 1.5, font = 2, cex = 0.85)
    mtext(sprintf("Spearman with the guide's proliferating fraction: PC%d %+.2f, PC%d %+.2f", pr[1], rho[1], pr[2], rho[2]), side = 3, adj = 0, line = 0.5, cex = 0.55, col = INK2)
    if (v == "raw" && pr[1] == 1) { sz <- c(20, 100, 500, 1000); legend("bottomleft", paste(sz, "cells"), pch = 21, pt.bg = "grey60", col = "white", pt.cex = 0.9 + 2.4 * sqrt(sz / max(info$cells)), bty = "n", cex = 0.7, title = "dot size = cells per guide", title.adj = 0, y.intersp = 1.5) }
    if (v == "raw" && pr[1] == 2) legend("bottomleft", names(kcol)[names(kcol) %in% info$KO], pt.bg = kcol[names(kcol) %in% info$KO], pch = 21, col = "white", bty = "n", cex = 0.7, pt.cex = 1.3, ncol = 2, title = "target gene", title.adj = 0) }
  mtext("PCA of the CROP-seq guides (pseudobulk)", outer = TRUE, side = 3, adj = 0, line = 2.6, font = 2, cex = 1.05)
  mtext(sprintf("%d single guides with >= 10 cells covering all %d KOs + NTC (multi-guide calls excluded); top %d variable genes, centred and scaled across guides; point size = cells (small guides are noisy); colour = target gene", nrow(info), length(unique(info$KO)) - 1, length(vg)), outer = TRUE, side = 3, adj = 0, line = 1.3, cex = 0.6, col = INK2)
})
