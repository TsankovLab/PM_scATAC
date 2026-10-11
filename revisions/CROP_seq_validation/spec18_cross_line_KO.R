# Spec 18: how well does each bulk KO reproduce ACROSS cell lines? Group-level KO profiles (DESeq2 Wald statistic vs the line's NTC; cc-adjusted (statcc) or not (stat0)), all genes expressed in all 4 lines,
# centred across the KOs of each line (removes the shared response). For every pair of lines (6 pairs), C = correlation matrix (KOs of line A x KOs of line B): the same-KO entry vs the entries of the
# other KOs (row / column z, rank, AUC). Significance per KO: permutation of the KO labels within every line (5000 permutations), statistic = mean same-KO r over the line pairs; BH across KOs.
# needs: spec1_prepare_stats.R
# outputs: spec18_cross_line_pairs.csv, spec18_cross_line_summary.csv, fig22_cross_line_KO.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
set.seed(18); NPERM <- 5000
bs <- readRDS(file.path(OUT, "spec_bulk_stats.rds")); ot <- read.csv(file.path(OUT, "spec_ontarget_knockdown.csv"))
genes <- Reduce(intersect, lapply(bs[LINES], `[[`, "genes")); cat("genes expressed in all 4 lines:", length(genes), "\n")
prof <- function(what) setNames(lapply(LINES, function(cl) { m <- bs[[cl]][[what]][genes, , drop = FALSE]; m[!is.finite(m)] <- 0; m - rowMeans(m) }), LINES)
prs <- combn(LINES, 2, simplify = FALSE); auc <- function(w, o) { if (!length(w) || !length(o)) return(NA_real_); r <- rank(c(w, o)); (sum(r[seq_along(w)]) - length(w) * (length(w) + 1) / 2) / (length(w) * length(o)) }

analyse <- function(what, method) {
  X <- prof(what); Cs <- lapply(prs, function(p) cor(X[[p[1]]], X[[p[2]]], method = method)); names(Cs) <- sapply(prs, paste, collapse = " | ")
  rows <- list()
  for (i in seq_along(prs)) { C <- Cs[[i]]; shared <- intersect(rownames(C), colnames(C))
    for (k in shared) { rr <- C[k, ]; cc <- C[, k]; zr <- (rr[k] - mean(rr)) / sd(rr); zc <- (cc[k] - mean(cc)) / sd(cc)
      rows[[length(rows) + 1]] <- data.frame(stat = what, method = method, pair = names(Cs)[i], KO = k, r_same = C[k, k], mean_r_other = mean(c(rr[-match(k, names(rr))], cc[-match(k, names(cc))])), z = mean(c(zr, zc)),
                                             rank_row = sum(rr > rr[k]) + 1, rank_col = sum(cc > cc[k]) + 1, n_kos = length(rr)) } }
  pr <- do.call(rbind, rows); kos <- sort(unique(pr$KO), method = "radix")
  # permutation of KO labels within every line
  obs <- sapply(kos, function(k) mean(pr$r_same[pr$KO == k])); nl <- sapply(LINES, function(l) ncol(X[[l]]))
  idx <- lapply(Cs, function(C) list(r = rownames(C), c = colnames(C)))
  perm <- replicate(NPERM, { pm <- lapply(LINES, function(l) { k <- colnames(X[[l]]); setNames(sample(seq_along(k)), k) })   # label -> profile index
    names(pm) <- LINES; sapply(kos, function(k) { v <- c(); for (i in seq_along(prs)) { a <- prs[[i]][1]; b <- prs[[i]][2]; if (k %in% names(pm[[a]]) && k %in% names(pm[[b]])) v <- c(v, Cs[[i]][pm[[a]][[k]], pm[[b]][[k]]]) }; mean(v) }) })
  p <- sapply(seq_along(kos), function(j) (1 + sum(perm[j, ] >= obs[j])) / (NPERM + 1))
  sm <- do.call(rbind, lapply(seq_along(kos), function(j) { d <- pr[pr$KO == kos[j], ]; ci <- Cs[d$pair]; same_v <- d$r_same
    other_v <- unlist(lapply(seq_len(nrow(d)), function(m) { C <- ci[[m]]; c(C[kos[j], colnames(C) != kos[j]], C[rownames(C) != kos[j], kos[j]]) }))
    data.frame(stat = what, method = method, KO = kos[j], n_line_pairs = nrow(d), mean_r_same = mean(same_v), mean_r_other = mean(other_v), mean_z = mean(d$z), frac_pairs_rank1 = mean(d$rank_row == 1 & d$rank_col == 1), AUC = auc(same_v, other_v), perm_p = p[j]) }))
  sm$perm_p_BH <- p.adjust(sm$perm_p, "BH"); list(pairs = pr, summary = sm, Cs = Cs, X = X) }

res <- list(); for (w in c("statcc", "stat0")) for (m in c("pearson", "spearman")) res[[paste(w, m)]] <- analyse(w, m)
pairs_all <- do.call(rbind, lapply(res, `[[`, "pairs")); sum_all <- do.call(rbind, lapply(res, `[[`, "summary")); rownames(pairs_all) <- rownames(sum_all) <- NULL
write.csv(pairs_all, file.path(OUT, "spec18_cross_line_pairs.csv"), row.names = FALSE); write.csv(sum_all, file.path(OUT, "spec18_cross_line_summary.csv"), row.names = FALSE); options(width = 200)
main <- res[["statcc pearson"]]; sm <- main$summary[order(-main$summary$mean_z), ]
cat("\n=== per KO, cc-adjusted statistic, Pearson (sorted by mean z) ===\n"); print(format(sm[, -(1:2)], digits = 3), row.names = FALSE)
cat("\n=== sensitivity: mean z and BH p per KO for each statistic / correlation ===\n"); 
z_tab <- sapply(names(res), function(n) res[[n]]$summary$mean_z[match(sm$KO, res[[n]]$summary$KO)]); rownames(z_tab) <- sm$KO; print(round(z_tab, 2))
cat("\nKOs with BH p < 0.05:", paste(sapply(names(res), function(n) sprintf("%s: %s", n, paste(res[[n]]$summary$KO[res[[n]]$summary$perm_p_BH < 0.05], collapse = ","))), collapse = " | "), "\n")
eff <- aggregate(cbind(lfc, padj) ~ KO, ot[ot$cell_line %in% LINES, ], function(x) mean(x)); cat("\nmean on-target log2FC across lines (for reference):\n"); print(round(setNames(eff$lfc, eff$KO), 2))

## ---------------------------------------------------------------- figure
Xc <- main$X; kos <- sm$KO; kcol <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7", "#008300", "#b3400f", "#8fa88a", "#5a7d9a", "#7a7a7a")[seq_along(kos)], kos)
units <- do.call(rbind, lapply(kos, function(k) do.call(rbind, lapply(LINES, function(l) if (k %in% colnames(Xc[[l]])) data.frame(KO = k, line = l, id = paste(k, l)) else NULL)))); n <- nrow(units)
big <- cor(do.call(cbind, lapply(seq_len(n), function(i) Xc[[units$line[i]]][, units$KO[i]]))); dimnames(big) <- list(units$id, units$id); big[outer(units$line, units$line, "==")] <- NA
save_fig("fig22_cross_line_KO", 15, 9.4, function() {
  layout(matrix(1:2, nrow = 1), widths = c(1, 1.35)); par(oma = c(0, 0, 4, 0.5), family = "sans")
  par(mar = c(1, 9, 8, 1)); vmax <- 0.5; image(seq_len(n), seq_len(n), heat_z(pmax(pmin(big, vmax), -vmax)), col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
  arr <- which(is.na(t(t(big)[, n:1])), arr.ind = TRUE); rect(arr[, 1] - 0.5, arr[, 2] - 0.5, arr[, 1] + 0.5, arr[, 2] + 0.5, col = "#f4f4f2", border = NA)
  lab <- sub("NCI-", "", sub("-211H", "", units$line)); axis(2, at = n:1, labels = paste(units$KO, lab), las = 1, tick = FALSE, cex.axis = 0.5, line = -0.5); axis(3, at = seq_len(n), labels = paste(units$KO, lab), las = 2, tick = FALSE, cex.axis = 0.5, line = -0.5)
  starts <- cumsum(c(1, table(factor(units$KO, kos))))[seq_along(kos)]; ends <- cumsum(table(factor(units$KO, kos)))
  for (i in seq_along(kos)) rect(starts[i] - 0.5, n - ends[i] + 0.5, ends[i] + 0.5, n - starts[i] + 1.5, border = INK, lwd = 1.3)
  mtext("A  Same KO in two lines (boxed blocks) vs other KOs", side = 3, adj = 0, line = 6.3, font = 2, cex = 0.8); mtext("cc-adjusted KO profiles, centred on each line; grey = same line (not comparable)", side = 3, adj = 0, line = 5.3, cex = 0.55, col = INK2)
  par(mar = c(4.2, 5, 8, 1)); k <- length(kos); plot(NA, xlim = c(0.3, 3 * k + 0.7), ylim = c(-0.6, 1.02), xaxt = "n", xlab = "", ylab = "correlation between two lines (Pearson)", bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8); abline(h = 0, col = GRID)
  at_s <- seq_len(k) * 3 - 1; at_d <- seq_len(k) * 3; pr <- main$pairs
  for (i in seq_len(k)) { s <- pr$r_same[pr$KO == kos[i]]; ci <- main$Cs[pr$pair[pr$KO == kos[i]]]; o <- unlist(lapply(ci, function(C) c(C[kos[i], colnames(C) != kos[i]], C[rownames(C) != kos[i], kos[i]])))
    boxplot(o, at = at_d[i], add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.85, col = adjustcolor(GREY, 0.9), border = INK2, medcol = INK, medlwd = 1.5, whisklty = 1, staplelty = 1); points(at_d[i] + runif(length(o), -0.25, 0.25), o, pch = 16, cex = 0.3, col = adjustcolor(INK2, 0.35))
    boxplot(s, at = at_s[i], add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.85, col = adjustcolor(kcol[kos[i]], 0.85), border = INK2, medcol = INK, medlwd = 1.5, whisklty = 1, staplelty = 1); points(at_s[i] + runif(length(s), -0.2, 0.2), s, pch = 21, bg = kcol[kos[i]], col = INK, cex = 0.8, lwd = 0.5)
    st <- sm[sm$KO == kos[i], ]; text((at_s[i] + at_d[i]) / 2, 1.0, sprintf("z %+.2f\n%s", st$mean_z, ifelse(st$perm_p_BH < 0.05, "*", ifelse(st$perm_p < 0.05, "·", ""))), cex = 0.6, col = INK2) }
  axis(1, at = (at_s + at_d) / 2, labels = kos, tick = FALSE, cex.axis = 0.75, font.axis = 2, las = 2, line = -0.3)
  mtext("B  Cross-line correlation of each KO", side = 3, adj = 0, line = 6.3, font = 2, cex = 0.8); mtext("coloured = same KO in two lines (one dot per line pair), grey = the KO vs other KOs; z = same-KO r standardised against the other KOs; * = BH p < 0.05, dot = permutation p < 0.05", side = 3, adj = 0, line = 5.3, cex = 0.5, col = INK2)
  mtext("How well does each bulk KO reproduce across cell lines?", outer = TRUE, side = 3, adj = 0, line = 2.4, font = 2, cex = 1)
  mtext(sprintf("bulk RNA-seq RUN1, %d genes expressed in all 4 lines, cc-adjusted DESeq2 statistics, 6 line pairs; KOs ordered by mean z; permutation test (%d shuffles of KO labels within lines) on the mean same-KO correlation", length(genes), NPERM), outer = TRUE, side = 3, adj = 0, line = 1.1, cex = 0.58, col = INK2)
})
