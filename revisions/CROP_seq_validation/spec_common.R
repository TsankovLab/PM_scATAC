# Helpers for the specificity analysis (spec1-spec3): screen loading, per-KO statistics, exact permutation test.
# Sourced after repro_common.R.
suppressPackageStartupMessages({ library(Matrix) })

SCREEN_RDS <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CRISPR_cropseq_analysis/_cellranger_raw_Filter_400_800_25/no_harmony/srt_filtered.rds"
SCREEN_LINE <- "NCI-H2052"     # the screen's cell line: CRISPR_cropseq.R reads DepMap ACH-000153 (= NCI-H2052) for its viability panel

# screen object filtered exactly as in CRISPR_cropseq.R: unique guide, guides with >9 cells, guides merged per gene, no combinations
load_screen <- function() {
  suppressPackageStartupMessages(library(Seurat))
  srt <- readRDS(SCREEN_RDS)
  srt <- srt[, srt$unique_guide %in% TRUE]
  tb <- table(srt$crispr_calls); srt <- srt[, srt$crispr_calls %in% names(tb)[tb > 9]]
  mg <- list(
    MEF2D = c('MEF2D-1','MEF2D-3','MEF2D−3|NTC−3','MEF2D−1|MEF2D−3','MEF2D−1|NTC−3','MEF2D−1|MEF2D−3|NTC−3','MEF2D−1|NTC−1'),
    MEF2A = c('MEF2A−3','MEF2A−1','MEF2A−3|NTC−3','MEF2A−1|NTC−3','MEF2A−3|NTC−1'),
    PITX1 = c('PITX1−2','PITX1−3','PITX1−1','PITX1−2|PITX1−3','PITX1−1|PITX1−2','PITX1−1|PITX1−3','NTC-1|PITX1-2','NTC-3|PITX1-2','NTC-3|PITX1-3','NTC-1|PITX1-3'),
    BPTF = c('BPTF−3','BPTF−2'), HMGA1 = 'HMGA1−3', TEAD2 = 'TEAD2-1', NTC = c('NTC-1','NTC-3'), TCF3 = c('TCF3-1','TCF3-3'),
    SOX9 = 'SOX9-3', TWIST1 = c('TWIST1-1','TWIST1-3'), TEAD4 = c('TEAD4-2','TEAD4-3'))
  mg <- lapply(mg, function(x) gsub("[−–-]", "_", x))
  srt$merged_call <- srt$crispr_calls
  for (i in names(mg)) srt$merged_call[srt$merged_call %in% mg[[i]]] <- i
  srt[, !grepl("\\|", srt$merged_call)]
}

# per-KO statistics vs NTC from a sparse log-normalised matrix Y (genes x cells), optionally after regressing a covariate (cc) out pooled over all cells
# (as Seurat::ScaleData(vars.to.regress) does). Uses group sums only, so the residual matrix is never built.
screen_ko_stats <- function(Y, group, cov, regress) {
  Y <- as(Y, "CsparseMatrix"); Y2 <- Y; Y2@x <- Y2@x^2
  cbar <- mean(cov); b <- if (regress) { (as.numeric(Y %*% (cov - cbar)) / (ncol(Y) - 1)) / var(cov) } else rep(0, nrow(Y))
  grp <- function(idx) {
    n <- length(idx); my <- rowSums(Y[, idx, drop = FALSE]) / n; c_ <- cov[idx]; mc <- mean(c_)
    vy <- (rowSums(Y2[, idx, drop = FALSE]) - n * my^2) / (n - 1); vc <- var(c_)
    cyc <- (as.numeric(Y[, idx, drop = FALSE] %*% c_) - n * my * mc) / (n - 1)
    list(n = n, mean = my - b * mc, var = pmax(vy + b^2 * vc - 2 * b * cyc, 0))
  }
  ntc <- grp(which(group == "NTC")); kos <- setdiff(sort(unique(group), method = "radix"), "NTC")
  diff <- sapply(setNames(kos, kos), function(k) grp(which(group == k))$mean - ntc$mean)
  tst  <- sapply(setNames(kos, kos), function(k) { g <- grp(which(group == k)); (g$mean - ntc$mean) / sqrt(g$var / g$n + ntc$var / ntc$n) })
  allg <- sort(unique(group), method = "radix")          # absolute group means (NTC included), residualised on cov when regress = TRUE
  abs_mean <- sapply(setNames(allg, allg), function(k) grp(which(group == k))$mean + b * cbar)
  list(diff = diff, t = tst, abs = abs_mean, n = c(NTC = ntc$n, sapply(setNames(kos, kos), function(k) sum(group == k))))
}

# all permutations of 1:k as rows
perm_matrix <- function(k) {
  if (k == 1) return(matrix(1L, 1, 1))
  p <- perm_matrix(k - 1)
  do.call(rbind, lapply(seq_len(k), function(i) cbind(i, matrix(setdiff(seq_len(k), i)[p], nrow = nrow(p)))))
}
PERM_CACHE <- list()
get_perms <- function(k) { key <- as.character(k); if (is.null(PERM_CACHE[[key]])) PERM_CACHE[[key]] <<- perm_matrix(k); PERM_CACHE[[key]] }

# Specificity of a rows x ORDER matrix C (rows = bulk profiles, columns = screen KOs in ORDER).
# Statistic T = mean over rows of the matched entry after standardising each row across the 8 screen KOs (z-score), so shared row offsets/scales drop out.
# p = exact permutation test over all 8! assignments of screen KOs to the columns.
evaluate_specificity <- function(C, row_ko, cols = ORDER) {
  PERMS <- get_perms(length(cols))
  C <- C[, cols, drop = FALSE]; j <- match(row_ko, cols); keep <- !is.na(j); C <- C[keep, , drop = FALSE]; j <- j[keep]; n <- nrow(C)
  Z <- t(scale(t(C))); Z[!is.finite(Z)] <- 0
  t_obs <- mean(Z[cbind(seq_len(n), j)])
  t_perm <- rowMeans(sapply(seq_len(n), function(i) Z[i, PERMS[, j[i]]]))
  ranks <- sapply(seq_len(n), function(i) sum(C[i, ] > C[i, j[i]]) + 1)
  data.frame(n_rows = n, matched_r = mean(C[cbind(seq_len(n), j)]), other_r = mean(C[-((j - 1) * n + seq_len(n))]),
             z_matched = t_obs, perm_p = mean(t_perm >= t_obs), top1 = sum(ranks == 1), top3 = sum(ranks <= 3), median_rank = median(ranks))
}

## ---- profile -> rows x ORDER matrices (used by spec3 and spec4) ----
prep_dataset <- function(m, centre, scale) { m[!is.finite(m)] <- 0; if (centre) m <- m - rowMeans(m); if (scale) { r <- rms_scale(m); m <- m / r; m[!is.finite(m)] <- 0 }; m }

# rows (line, KO in ORDER) x ORDER matrix of correlations between bulk KO profiles and screen KO profiles
cross_matrix <- function(S, Bl, centre, scale, method = "pearson", cols = ORDER) {
  D <- c(list(S), Bl); D <- lapply(D, prep_dataset, centre = centre, scale = scale); S <- D[[1]]; rows <- list(); lab <- list()
  for (i in seq_along(Bl)) { m <- D[[i + 1]]; for (ko in ORDER) if (ko %in% colnames(m)) { rows[[length(rows) + 1]] <- sapply(cols, function(cc) cor(m[, ko], S[, cc], method = method)); lab[[length(lab) + 1]] <- c(names(Bl)[i], ko) } }
  list(C = do.call(rbind, rows), lab = setNames(as.data.frame(do.call(rbind, lab), stringsAsFactors = FALSE), c("line", "ko")))
}
# signature approach: top-N up / down genes of each screen KO (screen t) scored in each bulk KO profile (bulk stat)
signature_matrix <- function(S, Bl, N, centre) {
  D <- lapply(c(list(S), Bl), prep_dataset, centre = centre, scale = FALSE); S <- D[[1]]; rows <- list(); lab <- list()
  sig <- lapply(ORDER, function(k) { o <- order(S[, k]); list(dn = o[seq_len(N)], up = rev(o)[seq_len(N)]) }); names(sig) <- ORDER
  for (i in seq_along(Bl)) { m <- D[[i + 1]]; for (ko in ORDER) if (ko %in% colnames(m)) { rows[[length(rows) + 1]] <- sapply(ORDER, function(cc) mean(m[sig[[cc]]$up, ko]) - mean(m[sig[[cc]]$dn, ko])); lab[[length(lab) + 1]] <- c(names(Bl)[i], ko) } }
  list(C = do.call(rbind, rows), lab = setNames(as.data.frame(do.call(rbind, lab), stringsAsFactors = FALSE), c("line", "ko")))
}



## ---- shared heatmap: bulk KO rows (grouped by line, diagonal ordered) x screen KO columns (best -> worst shared TFs, then screen-only KOs) ----
# C: rows x k ordered shared-TF columns (matched cell outlined); Cx: rows x m screen-only columns (or NULL); lab: data.frame(line, ko) in display order;
# n_cells: named vector of cell counts for all columns; row_labels: character(n); mark: TFs to flag with a dot; lines: line order
draw_ordered_matrix <- function(name, C, Cx, lab, n_cells, row_labels, main, sub1, sub2, sub3 = NULL, mark = NULL, lines = LINES) {
  A <- if (is.null(Cx)) C else cbind(C, Cx); n <- nrow(A); k <- ncol(C); kk <- ncol(A); vmax <- max(0.05, ceiling(quantile(abs(A), 0.99, names = FALSE) * 20) / 20)
  rk <- sapply(seq_len(n), function(i) sum(C[i, ] > C[i, match(lab$ko[i], colnames(C))]) + 1); lines_in <- lines[lines %in% lab$line]
  save_fig(name, 4.6 + 0.85 * kk, max(5, 2.6 + 0.24 * n), function() {
    par(mar = c(4.4, 11, 9.5, 6.2), family = "sans"); image(seq_len(kk), seq_len(n), heat_z(pmax(pmin(A, vmax), -vmax)), col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
    lb <- sprintf("%s (%d)%s", colnames(A), n_cells[colnames(A)], ifelse(colnames(A) %in% mark, "  ●", "")); sh <- seq_len(kk) <= k
    axis(3, at = which(sh), labels = lb[sh], las = 2, tick = FALSE, cex.axis = 0.85, line = -0.4, font.axis = 2, col.axis = INK)
    if (kk > k) axis(3, at = which(!sh), labels = lb[!sh], las = 2, tick = FALSE, cex.axis = 0.85, line = -0.4, font.axis = 1, col.axis = MUTED)
    axis(2, at = n:1, labels = row_labels, las = 1, tick = FALSE, cex.axis = 0.7, line = -0.4)
    for (i in seq_len(n)) for (j in seq_len(kk)) text(j, n - i + 1, sub("^([+-])0\\.", "\\1.", sprintf("%+.2f", A[i, j])), cex = 0.62, col = if (abs(A[i, j]) < 0.6 * vmax) INK2 else "white")
    for (i in seq_len(n)) { j <- match(lab$ko[i], colnames(C)); rect(j - 0.5, n - i + 0.5, j + 0.5, n - i + 1.5, border = INK, lwd = 1.8) }
    if (kk > k) abline(v = k + 0.5, col = "white", lwd = 5)
    ed <- which(lab$line[-1] != lab$line[-n]); abline(h = n - ed + 0.5, col = "white", lwd = 4); bd <- c(0, ed, n); par(xpd = NA); colw <- par("pin")[1] / kk; xin <- function(inch) 0.5 - inch / colw
    for (b in seq_along(lines_in)) { l <- lines_in[b]; y1 <- n - bd[b + 1] + 0.5; y2 <- n - bd[b] + 0.5; rect(xin(1.45), y1, xin(1.2), y2, col = LINE_COL[[l]], border = NA); text(xin(1.75), (y1 + y2) / 2, l, srt = 90, font = 2, cex = 0.7, col = LINE_COL[[l]]) }
    text(kk + 0.7, n:1, sprintf("rank %d/%d", rk, k), adj = 0, cex = 0.62, col = ifelse(rk == 1, "#b3400f", MUTED), font = ifelse(rk == 1, 2, 1))
    mtext(main, side = 3, adj = 0, line = 8.5, font = 2, cex = 0.8, at = xin(1.85)); mtext(sub1, side = 3, adj = 0, line = 7.5, cex = 0.55, col = INK2, at = xin(1.85)); mtext(sub2, side = 3, adj = 0, line = 6.5, cex = 0.55, col = INK2, at = xin(1.85))
    if (!is.null(sub3)) mtext(sub3, side = 3, adj = 0, line = 5.6, cex = 0.5, col = INK2, at = xin(1.85)); par(xpd = FALSE); draw_colorbar(c(0.3, 0.7), c(-vmax, vmax), "Pearson r") })
}
