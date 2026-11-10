# Helpers for the specificity analysis: screen loading, per-KO statistics. Sourced after repro_common.R. Used by spec1_prepare_stats.R.
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
