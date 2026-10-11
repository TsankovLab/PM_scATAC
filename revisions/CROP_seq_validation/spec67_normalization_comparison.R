# Spec 67: head-to-head comparison of the two normalizations used throughout this analysis for screen-vs-bulk KO-identity recovery:
#   (A) NTC-SUBTRACTED: KO-NTC diff profile (median KO minus median NTC, within each platform) -- spec44/45/47/58 convention.
#   (B) GROUP-CENTRED:   each feature's mean across that platform's own 9 group (NTC+8 KO) profiles subtracted, NO NTC involved -- spec62/64/65/66 convention.
# Both are computed for the SAME platform pairs (screen vs each of the 4 bulk lines), SAME feature spaces (20 Cm modules, top 1000/2000/3000 variable genes),
# and scored with the SAME diagonal metric (mean diagonal - mean off-diagonal of the 8x8 KO x KO Spearman matrix, top1, exact permutation p) so the two
# normalizations are directly comparable, not just eyeballed from separate scripts.
# needs: git_repo/files/scRNA_PM_modules.rds, spec_common.R, repro_common.R
# outputs: spec67_normalization_comparison.csv, fig5A_normalization_comparison_diagoff.{png,pdf}, fig5B_normalization_comparison_per_KO.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3"); TOPNS <- c(1000, 2000, 3000)
b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1)
cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds"); top50 <- lapply(cm_genes, head, 50); CM <- names(top50)
srt <- suppressWarnings(AddModuleScore(load_screen(), features = top50, name = "Cm_", seed = 1)); sgrp <- as.character(srt$merged_call)
scm <- as.matrix(srt@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(scm) <- CM; Y <- GetAssayData(srt, layer = "data")
sb <- CreateSeuratObject(counts = as.matrix(b$counts)); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
bcm <- as.matrix(sb@meta.data[rownames(meta), paste0("Cm_", seq_along(CM))]); colnames(bcm) <- CM
KOS <- ORDER; K <- length(KOS)
lv <- sapply(LINES, function(cl) matrixStats::rowVars(as.matrix(lg[, rownames(meta)[meta$Cell_Line == cl], drop = FALSE]))); vord <- intersect(rownames(lg)[order(rowMeans(lv), decreasing = TRUE)], rownames(Y))
spaces <- c("20 Cm modules", sprintf("top %d variable genes", TOPNS))
feat <- function(sp) if (sp == "20 Cm modules") list(S = scm, B = bcm) else { n <- as.integer(sub("top (\\d+) .*", "\\1", sp)); list(S = t(as.matrix(Y[vord[1:n], ])), B = t(lg[vord[1:n], ])) }
med <- function(M, ids) apply(M[ids, , drop = FALSE], 2, median)
perms <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(seq_len(K))) })
diagoff <- function(Ck) { dm <- mean(diag(Ck)); om <- mean(Ck[row(Ck) != col(Ck)]); pd <- apply(perms, 1, function(p) mean(Ck[cbind(seq_len(K), p)]))
  list(diag_mean = dm, off_mean = om, diag_minus_off = dm - om, top1 = sum(apply(Ck, 1, which.max) == seq_len(K)), perm_p = mean(pd >= dm - 1e-12)) }

## -- (A) NTC-subtracted profiles: KO - NTC (screen: all cells; bulk: per line), 8 x features matrix
ntc_profile <- function(M, grp, ntc_ids) { n <- med(M, ntc_ids); t(sapply(KOS, function(k) med(M, which(grp == k)) - n)) }   # screen uses logical grp
ntc_profile_bulk <- function(M, cl) { mm <- meta[meta$Cell_Line == cl, ]; ntc <- rownames(mm)[mm$Gene_Control == "NTC"]
  if (!length(ntc)) return(matrix(NA, K, ncol(M), dimnames = list(KOS, colnames(M))))
  t(sapply(KOS, function(k) { ko <- rownames(mm)[mm$Gene_Control == k]; if (!length(ko)) return(rep(NA, ncol(M))); med(M, ko) - med(M, ntc) })) }
## -- (B) group-centred profiles: subtract each feature's mean across the platform's own 9 group (NTC+8 KO) medians
group_matrix <- function(M, ids_fn) t(sapply(c("NTC", KOS), function(k) med(M, ids_fn(k))))
centre_cols <- function(G) G - matrix(colMeans(G), nrow(G), ncol(G), byrow = TRUE)

res <- list()
for (sp in spaces) { F <- feat(sp)
  S_ntc <- ntc_profile(F$S, sgrp, which(sgrp == "NTC")); S_grp <- centre_cols(group_matrix(F$S, function(k) which(sgrp == k)))[KOS, ]
  for (cl in LINES) {
    B_ntc <- ntc_profile_bulk(F$B, cl); B_grp <- centre_cols(group_matrix(F$B, function(k) rownames(meta)[meta$Cell_Line == cl & meta$Gene_Control == k]))[KOS, ]
    if (!all(is.na(B_ntc))) { Ck <- cor(t(S_ntc), t(B_ntc), method = "spearman"); m <- diagoff(Ck)
      res[[length(res) + 1]] <- data.frame(norm = "NTC-subtracted", line = cl, space = sp, m, own_per_KO = I(list(setNames(diag(Ck), KOS)))) }
    Ck2 <- cor(t(S_grp), t(B_grp), method = "spearman"); m2 <- diagoff(Ck2)
    res[[length(res) + 1]] <- data.frame(norm = "Group-centred", line = cl, space = sp, m2, own_per_KO = I(list(setNames(diag(Ck2), KOS)))) }
}
res <- do.call(rbind, res); pko <- do.call(rbind, lapply(seq_len(nrow(res)), function(i) data.frame(norm = res$norm[i], line = res$line[i], space = res$space[i], KO = KOS, own = res$own_per_KO[[i]])))
write.csv(res[, setdiff(names(res), "own_per_KO")], file.path(OUT, "spec67_normalization_comparison.csv"), row.names = FALSE)
options(width = 220); cat("\n=== normalization comparison: mean diag-off per (normalization, line), averaged over 4 spaces ===\n")
print(round(tapply(res$diag_minus_off, list(res$line, res$norm), mean), 3))
cat("\noverall mean diag-off per normalization (across all lines x spaces):\n"); print(round(tapply(res$diag_minus_off, res$norm, mean), 3))
cat("\noverall mean top1 (of 8) per normalization:\n"); print(round(tapply(res$top1, res$norm, mean), 2))
cat("\nfraction of (line x space) combos with diag-off > 0, per normalization:\n"); print(round(tapply(res$diag_minus_off > 0, res$norm, mean), 2))
cat("\nfraction with perm p < 0.10, per normalization:\n"); print(round(tapply(res$perm_p < 0.10, res$norm, mean), 2))
cat("\npaired comparison, same (line, space): which normalization wins?\n")
w <- merge(res[res$norm == "NTC-subtracted", c("line", "space", "diag_minus_off")], res[res$norm == "Group-centred", c("line", "space", "diag_minus_off")], by = c("line", "space"), suffixes = c("_NTC", "_group"))
w$winner <- ifelse(w$diag_minus_off_NTC > w$diag_minus_off_group, "NTC-subtracted", "Group-centred"); print(w); cat("\nwins:\n"); print(table(w$winner))
cat("\nmean diag-off per KO, per normalization (pooled over lines x spaces):\n"); print(round(tapply(pko$own, list(pko$KO, pko$norm), mean) - as.vector(tapply(res$off_mean, res$norm, mean))[rep(1:2, each = 8)], 3))

save_fig("fig5A_normalization_comparison_diagoff", 11, 6.4, function() {
  layout(matrix(1:4, nrow = 1)); par(oma = c(0, 0, 3.4, 0), family = "sans")
  for (sp in spaces) { d <- res[res$space == sp, ]; par(mar = c(7, 4.4, 2.6, 0.8))
    at <- c(seq_along(LINES) * 3 - 1, seq_along(LINES) * 3); vals <- c(d$diag_minus_off[d$norm == "NTC-subtracted"][match(LINES, d$line[d$norm == "NTC-subtracted"])], d$diag_minus_off[d$norm == "Group-centred"][match(LINES, d$line[d$norm == "Group-centred"])])
    cols <- rep(c("#4a3aa7", "#eda100"), each = length(LINES))
    bp <- barplot(vals[order(rep(seq_along(LINES), 2))], space = 0, xaxt = "n", ylim = range(c(vals, 0)) + c(-0.02, 0.02) * sign(range(c(vals,0))), border = NA, col = cols[order(rep(seq_along(LINES), 2))], ylab = "diagonal - off-diagonal")
    abline(h = 0, col = INK); axis(1, at = (bp[c(TRUE, FALSE)] + bp[c(FALSE, TRUE)]) / 2, labels = sub("NCI-|-211H", "", LINES), las = 2, tick = FALSE, cex.axis = 0.8)
    if (sp == spaces[1]) legend("bottomleft", legend = c("NTC-subtracted", "Group-centred"), fill = c("#4a3aa7", "#eda100"), border = NA, bty = "n", cex = 0.7)
    mtext(sp, side = 3, adj = 0, line = 0.3, font = 2, cex = 0.82) }
  mtext("Screen vs bulk: NTC-subtracted vs Group-centred normalization, diagonal quality per line", outer = TRUE, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.95)
  mtext("cell = mean(diagonal) - mean(off-diagonal) of the 8x8 screen-vs-bulk KO x KO Spearman matrix", outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.58, col = INK2)
})

save_fig("fig5B_normalization_comparison_per_KO", 9.5, 5.6, function() {
  ko_ord <- names(sort(tapply(pko$own[pko$norm == "NTC-subtracted"], pko$KO[pko$norm == "NTC-subtracted"], mean), decreasing = TRUE)); n_ko <- length(ko_ord); par(mar = c(5, 4.6, 3.6, 1), family = "sans")
  plot(NA, xlim = c(0.5, n_ko + 0.5), ylim = c(-1, 1), xaxt = "n", xlab = "", ylab = "Spearman rho (screen vs bulk, own KO)")
  axis(1, at = seq_len(n_ko), labels = FALSE); text(seq_len(n_ko), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), ko_ord, srt = 45, adj = 1, xpd = NA)
  abline(h = 0, col = INK, lty = 2)
  for (i in seq_len(n_ko)) boxplot(pko$own[pko$KO == ko_ord[i] & pko$norm == "NTC-subtracted"], at = i - 0.18, add = TRUE, boxwex = 0.3, col = adjustcolor("#4a3aa7", 0.5), border = INK2, outline = FALSE, axes = FALSE)
  for (i in seq_len(n_ko)) boxplot(pko$own[pko$KO == ko_ord[i] & pko$norm == "Group-centred"], at = i + 0.18, add = TRUE, boxwex = 0.3, col = adjustcolor("#eda100", 0.5), border = INK2, outline = FALSE, axes = FALSE)
  legend("bottomleft", legend = c("NTC-subtracted", "Group-centred"), fill = c(adjustcolor("#4a3aa7", 0.5), adjustcolor("#eda100", 0.5)), border = INK2, bty = "n", cex = 0.65)
  mtext("Per-TF own-KO Spearman rho, screen vs bulk, across 4 lines x 4 spaces", side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext("boxes pool the 4 lines x 4 feature spaces per TF; own-KO rho only, not corrected for the off-diagonal null", side = 3, adj = 0, line = 0.9, cex = 0.58, col = INK2)
})
