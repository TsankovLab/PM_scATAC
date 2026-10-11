# Spec 10: best-method matrix (cell-cycle adjusted statistics, all expressed genes, centred across KOs) restricted to the reproducible TFs of spec9,
#          with each bulk KO represented (a) by the average of its replicates and (b) ONLY by the replicate best correlated with the matching screen KO.
# Replicate-level statistic: per line, log2 CP10K ~ cell-cycle index + KO group (pooled per-gene residual variance); for every KO replicate
#   t = (cc-adjusted log2 CP10K of the replicate - mean of the NTC replicates) / sqrt(pooled variance * (1 + 1/n_NTC)), the per-replicate analogue of the DESeq2 Wald statistic.
# Selection-bias control: the best replicate is picked separately for every cell of the matrix (the advantage selection alone gives to the matched column).
# needs: spec1_prepare_stats.R, spec9_reproducible_TFs.R
# outputs: spec_best_replicate_selected.csv, spec_best_replicate_stats.csv, fig12_best_replicate_reproducible_TFs.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_variants.R"))
ctx <- load_spec_context(LINES); ss <- ctx$ss; bs <- ctx$bs
KEEP <- read.csv(file.path(OUT, "spec_TF_reproducibility.csv")); KEEP <- KEEP$TF[KEEP$keep]; cat("reproducible TFs:", paste(KEEP, collapse = ", "), "\n")
genes <- Reduce(intersect, c(list(ss$genes), lapply(bs[LINES], `[[`, "genes")))
S <- prep_dataset(ss$tcc[genes, c(ORDER, EXTRA)], TRUE, FALSE)                       # screen: cc-regressed Welch t, centred across its KOs (as F4)

## ---------------------------------------------------------------- replicate-level statistics per line
rep_prof <- list()
for (cl in LINES) {
  m <- ctx$meta[ctx$meta$Cell_Line == cl, ]; s <- rownames(m); y <- ctx$lg[genes, s, drop = FALSE]; cc <- ctx$ccb$cc[match(s, ctx$ccb$sample)]
  g <- factor(m$Gene_Control, levels = c("NTC", sort(setdiff(unique(m$Gene_Control), "NTC"), method = "radix"))); X <- model.matrix(~ cc + g); q <- qr(X)
  res <- qr.resid(q, t(y)); df <- length(s) - q$rank; sigma2 <- pmax(colSums(res^2) / df, 1e-4); bcc <- qr.coef(q, t(y))["cc", ]
  yadj <- y - outer(bcc, cc - mean(cc)); ntc_s <- s[m$Gene_Control == "NTC"]; ntc <- rowMeans(yadj[, ntc_s, drop = FALSE]); ko_s <- s[m$Gene_Control != "NTC"]
  Tm <- (yadj[, ko_s, drop = FALSE] - ntc) / sqrt(sigma2 * (1 + 1 / length(ntc_s))); Tm[!is.finite(Tm)] <- 0
  Tm <- Tm - rowMeans(Tm)                                                            # centred across all KO replicates of the line
  attr(Tm, "ko") <- setNames(as.character(m[ko_s, "Gene_Control"]), ko_s); rep_prof[[cl]] <- Tm
}
# sanity check: the mean of a KO's replicate statistics vs the DESeq2 Wald statistic of that KO (cc-adjusted)
chk <- unlist(lapply(LINES, function(cl) { Tm <- rep_prof[[cl]]; ko <- attr(Tm, "ko"); sapply(intersect(KEEP, unique(ko)), function(k) { v <- rowMeans(Tm[, names(ko)[ko == k], drop = FALSE]); w <- bs[[cl]]$statcc[genes, k]; w[!is.finite(w)] <- 0; cor(v, w - rowMeans(bs[[cl]]$statcc[genes, ], na.rm = TRUE)) }) }))
cat(sprintf("sanity: replicate-mean t vs DESeq2 Wald (centred), correlation across (line, TF): median %.2f, range %.2f-%.2f\n", median(chk), min(chk), max(chk)))

## ---------------------------------------------------------------- matrices
rows_lines <- c("NCI-H2052", "MSTO-211H", "NCI-H2452", "NCI-H28")
avg_rows <- list(); best_rows <- list(); ctl_rows <- list(); lab <- list(); sel <- list()
for (cl in rows_lines) { Tm <- rep_prof[[cl]]; ko <- attr(Tm, "ko")
  for (k in KEEP) { ss_k <- names(ko)[ko == k]; if (!length(ss_k)) next
    cr <- sapply(ss_k, function(sm) sapply(KEEP, function(j) cor(Tm[, sm], S[, j])))               # KEEP x replicates
    if (length(ss_k) == 1) cr <- matrix(cr, ncol = 1, dimnames = list(KEEP, ss_k))
    avgv <- rowMeans(Tm[, ss_k, drop = FALSE]); avg_rows[[length(avg_rows) + 1]] <- sapply(KEEP, function(j) cor(avgv, S[, j]))
    b <- ss_k[which.max(cr[k, ])]; best_rows[[length(best_rows) + 1]] <- cr[, b]; ctl_rows[[length(ctl_rows) + 1]] <- apply(cr, 1, max)
    lab[[length(lab) + 1]] <- c(cl, k); sel[[length(sel) + 1]] <- data.frame(line = cl, TF = k, n_replicates = length(ss_k), best_replicate = b, r_best = cr[k, b], r_other_replicate = if (length(ss_k) > 1) cr[k, setdiff(ss_k, b)] else NA) } }
lab <- setNames(as.data.frame(do.call(rbind, lab), stringsAsFactors = FALSE), c("line", "ko")); sel <- do.call(rbind, sel)
mk <- function(l) { M <- do.call(rbind, l); dimnames(M) <- list(paste(lab$line, lab$ko), KEEP); M }
CA <- mk(avg_rows); CB <- mk(best_rows); CC <- mk(ctl_rows)
write.csv(sel, file.path(OUT, "spec_best_replicate_selected.csv"), row.names = FALSE); write.csv(CB, file.path(OUT, "spec_best_replicate_matrix.csv"))

sets <- list(`3 lines` = lab$line %in% THREE, `4 lines` = rep(TRUE, nrow(lab)), `NCI-H2052 only` = lab$line == "NCI-H2052")
stt <- do.call(rbind, lapply(names(sets), function(sn) rbind(
  data.frame(variant = "replicates averaged", rows = sn, evaluate_specificity(CA[sets[[sn]], , drop = FALSE], lab$ko[sets[[sn]]], cols = KEEP)),
  data.frame(variant = "best replicate", rows = sn, evaluate_specificity(CB[sets[[sn]], , drop = FALSE], lab$ko[sets[[sn]]], cols = KEEP)),
  data.frame(variant = "control: best replicate per cell", rows = sn, evaluate_specificity(CC[sets[[sn]], , drop = FALSE], lab$ko[sets[[sn]]], cols = KEEP)))))
write.csv(stt, file.path(OUT, "spec_best_replicate_stats.csv"), row.names = FALSE); options(width = 200); print(format(stt, digits = 3), row.names = FALSE)
rank_k <- sapply(seq_len(nrow(lab)), function(i) sum(CB[i, ] > CB[i, lab$ko[i]]) + 1)

## ---------------------------------------------------------------- figure
n <- nrow(CB); k <- ncol(CB); vmax <- ceiling(quantile(abs(CB), 0.99, names = FALSE) * 20) / 20
rl <- sprintf("%s  (%s)", lab$ko, sub(".*_", "", sel$best_replicate))
save_fig("fig12_best_replicate_reproducible_TFs", 8.6, 9.8, function() {
  par(mar = c(5, 10, 13, 6.5), family = "sans"); image(seq_len(k), seq_len(n), heat_z(pmax(pmin(CB, vmax), -vmax)), col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
  axis(3, at = seq_len(k), labels = sprintf("%s (%d)", KEEP, ss$n[KEEP]), las = 2, tick = FALSE, cex.axis = 0.95, line = -0.4, font.axis = 2); axis(2, at = n:1, labels = rl, las = 1, tick = FALSE, cex.axis = 0.82, line = -0.4)
  for (i in seq_len(n)) for (j in seq_len(k)) text(j, n - i + 1, sub("^([+-])0\\.", "\\1.", sprintf("%+.2f", CB[i, j])), cex = 0.75, col = if (abs(CB[i, j]) < 0.6 * vmax) INK2 else "white")
  for (i in seq_len(n)) { j <- match(lab$ko[i], KEEP); rect(j - 0.5, n - i + 0.5, j + 0.5, n - i + 1.5, border = INK, lwd = 2) }
  edges <- which(lab$line[-1] != lab$line[-n]); abline(h = n - edges + 0.5, col = "white", lwd = 5); bounds <- c(0, edges, n); par(xpd = NA)
  for (b in seq_along(rows_lines)) { y1 <- n - bounds[b + 1] + 0.5; y2 <- n - bounds[b] + 0.5; rect(-1.25, y1, -1.0, y2, col = LINE_COL[[rows_lines[b]]], border = NA); text(-1.65, (y1 + y2) / 2, rows_lines[b], srt = 90, font = 2, cex = 0.9, col = LINE_COL[[rows_lines[b]]]) }
  text(k + 0.7, n:1, sprintf("rank %d/%d", rank_k, k), adj = 0, cex = 0.75, col = ifelse(rank_k <= 1, "#b3400f", MUTED), font = ifelse(rank_k <= 1, 2, 1))
  at <- -1.65; g <- function(v, rs) stt[stt$variant == v & stt$rows == rs, ]
  mtext(sprintf("Best method, reproducible TFs only, best replicate per KO (%s)", paste(KEEP, collapse = ", ")), side = 3, adj = 0, line = 11.4, font = 2, cex = 0.9, at = at)
  mtext("rows = bulk KO, ONLY the replicate best correlated with the matching screen KO (replicate in brackets); columns = CROP-seq KO (cells); cc-adjusted statistics, centred across KOs", side = 3, adj = 0, line = 10.2, cex = 0.55, col = INK2, at = at)
  mtext(sprintf("best replicate: 3 lines z = %.2f (p = %.3f), 4 lines z = %.2f (p = %.3f), NCI-H2052 z = %.2f (p = %.3f)", g("best replicate", "3 lines")$z_matched, g("best replicate", "3 lines")$perm_p, g("best replicate", "4 lines")$z_matched, g("best replicate", "4 lines")$perm_p, g("best replicate", "NCI-H2052 only")$z_matched, g("best replicate", "NCI-H2052 only")$perm_p), side = 3, adj = 0, line = 9.2, cex = 0.55, col = INK2, at = at)
  mtext(sprintf("replicates averaged: 3 lines z = %.2f (p = %.3f) | selection-bias control (best replicate picked per cell): z = %.2f (p = %.3f)", g("replicates averaged", "3 lines")$z_matched, g("replicates averaged", "3 lines")$perm_p, g("control: best replicate per cell", "3 lines")$z_matched, g("control: best replicate per cell", "3 lines")$perm_p), side = 3, adj = 0, line = 8.2, cex = 0.55, col = INK2, at = at)
  par(xpd = FALSE); draw_colorbar(c(0.3, 0.7), c(-vmax, vmax), "Pearson r")
})
