# Spec 11: replicate vs group (KO) structure in the RUN1 bulk samples.
#  (a) within each cell line: correlation between every pair of samples over the top variable genes; distribution for replicates of the SAME KO (incl. NTC)
#      vs pairs of DIFFERENT KOs. Pearson and Spearman, raw and "centred" (gene-wise mean of the line's samples removed, i.e. deviations from the line mean).
#  (b) across cell lines: the same comparison per KO -- the same KO in two different lines vs different KOs in two lines (within each line pair).
#  (c) variance partition per line: fraction of top-variable-gene variance explained by KO group (adjusted R^2).
# Genes: per line the top 3000 by variance among genes with mean CP10K >= 0.1 in that line; across lines the top 3000 by mean within-line variance among genes expressed in all 4 lines.
# AUC = P(random same-KO pair has a higher correlation than random different-KO pair); 0.5 = no group structure. Pairs share samples, so AUC is descriptive.
# outputs: spec_repl_within_line.csv, spec_repl_across_lines.csv, spec_variance_partition.csv, fig13_within_line_boxplots, fig14_across_lines_per_KO
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R"))
set.seed(5)
b <- load_bulk(); cp10k <- b$cp10k; meta <- b$meta; lg <- log2(cp10k + 1); NG <- 3000; MINEXP <- 0.1     # top-variable genes are chosen among genes with mean CP10K >= 0.1 (~350 reads at 35M depth)
auc <- function(w, o) { if (!length(w) || !length(o)) return(NA_real_); r <- rank(c(w, o)); (sum(r[seq_along(w)]) - length(w) * (length(w) + 1) / 2) / (length(w) * length(o)) }
line_genes <- function(cl) { s <- rownames(meta)[meta$Cell_Line == cl]; ok <- rowMeans(cp10k[, s]) >= MINEXP; v <- rowVars(lg[ok, s]); names(sort(setNames(v, rownames(lg)[ok]), decreasing = TRUE))[seq_len(NG)] }
centre_line <- function(y) y - rowMeans(y)

## ---------------------------------------------------------------- (a) within each line
within <- list(); var_part <- list(); pairs_store <- list()
for (cl in LINES) {
  s <- rownames(meta)[meta$Cell_Line == cl]; g <- line_genes(cl); y <- lg[g, s]; grp <- meta[s, "Gene_Control"]
  for (method in c("pearson", "spearman")) for (cen in c(FALSE, TRUE)) {
    C <- cor(if (cen) centre_line(y) else y, method = method); ut <- which(upper.tri(C), arr.ind = TRUE)
    d <- data.frame(line = cl, method = method, centred = cen, a = s[ut[, 1]], b = s[ut[, 2]], r = C[ut], same = grp[ut[, 1]] == grp[ut[, 2]])
    d$ko <- ifelse(d$same, as.character(grp[ut[, 1]]), NA); pairs_store[[length(pairs_store) + 1]] <- d
    w <- d$r[d$same]; o <- d$r[!d$same]
    within[[length(within) + 1]] <- data.frame(line = cl, method = method, centred = cen, n_samples = length(s), n_same_pairs = length(w), n_diff_pairs = length(o), median_same = median(w), median_diff = median(o), delta_median = median(w) - median(o), AUC = auc(w, o))
  }
  # variance partition: adjusted R^2 of KO group per gene
  n <- length(s); k <- length(unique(grp)); mu <- rowMeans(y); gsum <- t(rowsum(t(y), grp)); ng <- as.numeric(table(grp)[colnames(gsum)]); gm <- sweep(gsum, 2, ng, "/")
  ssb <- as.numeric(((gm - mu)^2) %*% ng); sst <- rowSums((y - mu)^2); ssw <- sst - ssb; r2 <- ssb / sst; r2adj <- 1 - (ssw / (n - k)) / (sst / (n - 1))
  var_part[[cl]] <- data.frame(line = cl, n_samples = n, n_groups = k, median_R2 = median(r2), median_adj_R2 = median(r2adj), frac_variance_between_KO = sum(ssb) / sum(sst))
}
within <- do.call(rbind, within); var_part <- do.call(rbind, var_part); pairs_all <- do.call(rbind, pairs_store)
write.csv(within, file.path(OUT, "spec_repl_within_line.csv"), row.names = FALSE); write.csv(var_part, file.path(OUT, "spec_variance_partition.csv"), row.names = FALSE)
options(width = 200); cat("\n=== (a) within cell line: replicates of the same KO vs different KOs ===\n"); print(format(within, digits = 3), row.names = FALSE)
cat("\n=== (c) variance partition (top-variable genes): share explained by KO group ===\n"); print(format(var_part, digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- (b) across lines, per KO
ok_all <- Reduce(`&`, lapply(LINES, function(cl) rowMeans(cp10k[, rownames(meta)[meta$Cell_Line == cl]]) >= MINEXP))
vl <- sapply(LINES, function(cl) rowVars(lg[, rownames(meta)[meta$Cell_Line == cl]])); rownames(vl) <- rownames(lg)
gx <- names(sort(setNames(rowMeans(vl)[ok_all], rownames(lg)[ok_all]), decreasing = TRUE))[seq_len(NG)]
across <- list(); pool <- list()
for (method in c("pearson", "spearman")) for (cen in c(FALSE, TRUE)) {
  Y <- lapply(LINES, function(cl) { s <- rownames(meta)[meta$Cell_Line == cl]; y <- lg[gx, s]; if (cen) centre_line(y) else y }); names(Y) <- LINES
  for (i in 1:(length(LINES) - 1)) for (j in (i + 1):length(LINES)) {
    A <- LINES[i]; B <- LINES[j]; C <- cor(Y[[A]], Y[[B]], method = method); ga <- meta[colnames(Y[[A]]), "Gene_Control"]; gb <- meta[colnames(Y[[B]]), "Gene_Control"]
    same <- outer(ga, gb, "=="); zc <- (C - mean(C)) / sd(C)
    pool[[length(pool) + 1]] <- data.frame(method = method, centred = cen, pair = paste(A, B, sep = " | "), r = as.vector(C), z = as.vector(zc), same = as.vector(same), ko = ifelse(as.vector(same), as.vector(outer(ga, gb, function(x, y) x)), NA))
    for (k in intersect(ga, gb)) { w <- as.vector(C[ga == k, gb == k]); o <- c(as.vector(C[ga == k, gb != k]), as.vector(C[ga != k, gb == k]))
      across[[length(across) + 1]] <- data.frame(method = method, centred = cen, pair = paste(A, B, sep = " | "), KO = k, n_same = length(w), n_diff = length(o), mean_r_same = mean(w), mean_r_diff = mean(o), AUC = auc(w, o)) }
  }
}
across <- do.call(rbind, across); pool <- do.call(rbind, pool); write.csv(across, file.path(OUT, "spec_repl_across_lines.csv"), row.names = FALSE)
ko_auc <- aggregate(AUC ~ method + centred + KO, across, mean)
cat("\n=== (b) across lines: mean AUC of same-KO vs different-KO pairs, per KO (mean over line pairs; 0.5 = none) ===\n")
for (m in c("pearson", "spearman")) for (cn in c(FALSE, TRUE)) { d <- ko_auc[ko_auc$method == m & ko_auc$centred == cn, ]; cat(sprintf("%s, %s: ", m, ifelse(cn, "centred", "raw"))); cat(paste(sprintf("%s %.2f", d$KO, d$AUC)[order(-d$AUC)], collapse = " | "), "\n") }
ov <- sapply(split(pool, list(pool$method, pool$centred)), function(d) c(AUC = auc(d$z[d$same], d$z[!d$same]), n_same = sum(d$same), n_diff = sum(!d$same))); cat("\npooled across line pairs (z within line pair): "); print(round(t(ov), 3))

## ---------------------------------------------------------------- boxplot helper: one pair (coloured = same KO, grey = different KO) per group
pair_boxes <- function(same, diff, labels, cols, main, sub, stats, ylab) {
  k <- length(labels); at_s <- seq_len(k) * 3 - 1; at_d <- seq_len(k) * 3
  yl <- range(unlist(c(same, diff))); yl[2] <- yl[2] + 0.14 * diff(yl)
  plot(NA, xlim = c(0.3, 3 * k + 0.7), ylim = yl, xaxt = "n", xlab = "", ylab = ylab, bty = "n", cex.axis = 0.8, cex.lab = 0.8, las = 1)
  abline(h = pretty(yl), col = GRID, lwd = 0.6)
  boxplot(diff, at = at_d, add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.85, col = adjustcolor(GREY, 0.9), border = INK2, medcol = INK, medlwd = 1.8, whisklty = 1, staplelty = 1)
  boxplot(same, at = at_s, add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.85, col = adjustcolor(cols, 0.85), border = INK2, medcol = INK, medlwd = 1.8, whisklty = 1, staplelty = 1)
  for (i in seq_len(k)) points(at_d[i] + runif(length(diff[[i]]), -0.25, 0.25), diff[[i]], pch = 16, col = adjustcolor(INK2, 0.3), cex = 0.35)
  for (i in seq_len(k)) points(at_s[i] + runif(length(same[[i]]), -0.22, 0.22), same[[i]], pch = 21, bg = cols[i], col = INK, cex = 0.65, lwd = 0.5)
  axis(1, at = (at_s + at_d) / 2, labels = labels, tick = FALSE, cex.axis = 0.8, font.axis = 2, line = -0.3)
  text((at_s + at_d) / 2, yl[2], stats, cex = 0.6, col = INK2)
  mtext(main, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.85); mtext(sub, side = 3, adj = 0, line = 0.4, cex = 0.6, col = INK2)
}

## ---------------------------------------------------------------- figure (a): within each line
save_fig("fig13_within_line_boxplots", 14, 8.4, function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4.5, 0), mar = c(3.2, 4.6, 3.6, 0.8), family = "sans")
  for (method in c("pearson", "spearman")) for (cen in c(FALSE, TRUE)) {
    d <- pairs_all[pairs_all$method == method & pairs_all$centred == cen, ]
    same <- lapply(LINES, function(cl) d$r[d$line == cl & d$same]); diffv <- lapply(LINES, function(cl) d$r[d$line == cl & !d$same])
    st <- sapply(LINES, function(cl) { w <- within[within$line == cl & within$method == method & within$centred == cen, ]; sprintf("AUC %.2f\nn=%d | %d", w$AUC, w$n_same_pairs, w$n_diff_pairs) })
    pair_boxes(same, diffv, sub("NCI-", "", LINES), LINE_COL[LINES], sprintf("%s, %s", method, ifelse(cen, "centred on the line mean", "raw")),
               if (cen) "deviations from each line's mean profile" else "dominated by baseline expression (correlations all high)", st, paste(method, "r between two samples"))
  }
  mtext(sprintf("Replicates of the same KO vs different KOs within each cell line (RUN1; top %d variable genes)", NG), outer = TRUE, side = 3, adj = 0, line = 2.6, font = 2, cex = 1)
  mtext("coloured box = pairs of replicates of the same KO (incl. NTC), grey box = pairs of different KOs; dots = individual sample pairs; AUC = P(same-KO pair correlates higher than different-KO pair); n = same | different pairs", outer = TRUE, side = 3, adj = 0, line = 1.3, cex = 0.6, col = INK2)
})

## ---------------------------------------------------------------- figure (b): across lines, one pair per line pair + per-KO AUC
save_fig("fig14_across_lines_per_KO", 15, 10.5, function() {
  layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 1.15)); par(oma = c(0, 0, 4.5, 0), family = "sans")
  prs <- unique(pool$pair); lab_pr <- gsub("NCI-|-211H", "", gsub(" \\| ", " / ", prs))
  for (method in c("pearson", "spearman")) for (cen in c(FALSE, TRUE)) {
    d <- pool[pool$method == method & pool$centred == cen, ]; same <- lapply(prs, function(p) d$r[d$pair == p & d$same]); diffv <- lapply(prs, function(p) d$r[d$pair == p & !d$same])
    st <- sapply(seq_along(prs), function(i) sprintf("AUC %.2f", auc(same[[i]], diffv[[i]])))
    par(mar = c(3.4, 4.6, 3.6, 0.8)); pair_boxes(same, diffv, lab_pr, rep("#4a3aa7", length(prs)), sprintf("%s, %s", method, ifelse(cen, "centred on the line mean", "raw")),
               sprintf("same KO in two lines vs different KOs; all line pairs pooled: AUC %.2f", auc(d$z[d$same], d$z[!d$same])), st, paste(method, "r across lines"))
  }
  d <- across[across$method == "pearson" & across$centred, ]; kos <- sort(unique(d$KO), method = "radix"); pr2 <- sort(unique(d$pair)); M <- sapply(pr2, function(p) sapply(kos, function(k) { v <- d$AUC[d$pair == p & d$KO == k]; if (length(v)) v else NA }))
  par(mar = c(3.5, 7, 9.5, 2)); image(seq_along(pr2), seq_along(kos), t(ifelse(is.na(M), 0.5, pmin(pmax(M, 0), 1)))[, length(kos):1], col = DIV, zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_along(kos)) for (j in seq_along(pr2)) text(j, length(kos) - i + 1, ifelse(is.na(M[i, j]), "n/a", sprintf("%.2f", M[i, j])), cex = 0.7)
  axis(2, at = length(kos):1, labels = kos, las = 1, tick = FALSE, cex.axis = 0.85); mtext(gsub("NCI-|-211H", "", gsub(" \\| ", " / ", pr2)), side = 3, at = seq_along(pr2), las = 2, cex = 0.65, line = 0.3, adj = 0)
  kp <- ko_auc[ko_auc$method == "pearson" & ko_auc$centred, ]; kp <- kp[order(-kp$AUC), ]
  mtext(sprintf("per-KO AUC (Pearson, centred); mean over line pairs: %s", paste(sprintf("%s %.2f", kp$KO, kp$AUC), collapse = " ")), side = 3, adj = 0, line = 6.4, cex = 0.6, col = INK2)
  mtext("Per KO: same KO across two cell lines vs different KOs (AUC 0.5 = no KO-specific similarity across lines)", side = 3, adj = 0, line = 7.6, font = 2, cex = 0.8)
  mtext(sprintf("Is the KO signature shared across cell lines? (top %d variable genes, RUN1)", NG), outer = TRUE, side = 3, adj = 0, line = 2.6, font = 2, cex = 1)
  mtext("purple box = the same KO in the two lines of a pair, grey box = different KOs; top row raw, second row centred on each line's mean", outer = TRUE, side = 3, adj = 0, line = 1.3, cex = 0.6, col = INK2)
})

## ---------------------------------------------------------------- compact single-panel version for presentation (Spearman, centred on the line mean)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
save_fig("fig13c_within_line_spearman_centred_compact", 4.9, 3.8, res = 300, fun = function() {
  par(mar = c(2.4, 3.5, 4.2, 0.5), mgp = c(2.2, 0.6, 0), family = "sans")
  d <- pairs_all[pairs_all$method == "spearman" & pairs_all$centred, ]; k <- length(LINES)
  same <- lapply(LINES, function(cl) d$r[d$line == cl & d$same]); diffv <- lapply(LINES, function(cl) d$r[d$line == cl & !d$same])
  pw <- sapply(seq_len(k), function(i) wilcox.test(same[[i]], diffv[[i]], alternative = "greater")$p.value)
  cat("\nfig13c: one-sided Wilcoxon p (same-KO r > different-KO r), per line:\n"); print(setNames(round(pw, 4), sub("NCI-", "", LINES)))
  at_s <- seq_len(k) * 3 - 1; at_d <- seq_len(k) * 3; yl <- range(unlist(c(same, diffv))) + c(-0.03, 0.28)
  plot(NA, xlim = c(0.3, 3 * k + 0.7), ylim = yl, xaxt = "n", xlab = "", ylab = "Spearman r between samples", bty = "n", las = 1, cex.axis = 0.7, cex.lab = 0.75)
  abline(h = 0, col = GRID, lwd = 0.8)
  boxplot(diffv, at = at_d, add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.85, col = adjustcolor(GREY, 0.9), border = INK2, medcol = INK, medlwd = 1.6, whisklty = 1, staplelty = 1, lwd = 0.8)
  boxplot(same, at = at_s, add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.85, col = adjustcolor(LINE_COL[LINES], 0.85), border = INK2, medcol = INK, medlwd = 1.6, whisklty = 1, staplelty = 1, lwd = 0.8)
  for (i in seq_len(k)) points(at_d[i] + runif(length(diffv[[i]]), -0.25, 0.25), diffv[[i]], pch = 16, col = adjustcolor(INK2, 0.3), cex = 0.3)
  for (i in seq_len(k)) points(at_s[i] + runif(length(same[[i]]), -0.22, 0.22), same[[i]], pch = 21, bg = LINE_COL[LINES][i], col = INK, cex = 0.5, lwd = 0.4)
  axis(1, at = (at_s + at_d) / 2, labels = sub("NCI-", "", LINES), tick = FALSE, cex.axis = 0.75, font.axis = 2, line = -0.6)
  a <- sapply(LINES, function(cl) within$AUC[within$line == cl & within$method == "spearman" & within$centred])
  off <- diff(yl) * 0.07
  text((at_s + at_d) / 2, yl[2], sprintf("AUC %.2f", a), cex = 0.65, col = INK2)
  text((at_s + at_d) / 2, yl[2] - off, stars(pw), cex = 0.75, col = INK, font = ifelse(stars(pw) == "ns", 1, 2))
  mtext("Replicates cluster by KO within each cell line", side = 3, adj = 0, line = 2.9, font = 2, cex = 0.8)
  mtext(sprintf("Spearman r over the top %d variable genes, centred on the line mean", NG), side = 3, adj = 0, line = 1.95, cex = 0.55, col = INK2)
  mtext("coloured = same-KO pairs, grey = different-KO pairs (dots = sample pairs)", side = 3, adj = 0, line = 1.15, cex = 0.55, col = INK2)
  mtext("*** p<0.001, ** p<0.01, * p<0.05, ns = not sig. (one-sided Wilcoxon, same-KO > diff-KO)", side = 3, adj = 0, line = 0.35, cex = 0.5, col = INK2)
})
