# Spec 23: QC of the FULL AUDITED bulk sample universe (load_full_audited_bulk(), repro_common.R) -- the 8 KOs shared with the CROP-seq screen (+
#   NTC), every status (not just RUN1's "OK"), including samples entirely missing from the sample sheet (H2452_TWIST1_R2, H28_HMGA1_R2,
#   H28_NTC_R1/2/3) and the registered NTC reassignment (H28_NTC_R1/2/3 -> NCI-H2052, H2052_NTC_R1/2 -> NCI-H28). Same 76-sample set as
#   spec25_bulk_PCA_prefilter.R's fig37/38/39. Highlights VISUAL_OUTLIERS (repro_common.R), the samples flagged in that same PCA.
#   excludes SNAI2, RUNX1, RUNX2 (present in bulk but not in the CROP-seq panel) and NR4A2/WT (bulk-only controls not used elsewhere in this pipeline)
# A: library size (M reads); B: genes detected (count > 0); C: read-assignment breakdown from STAR's ReadsPerGene.out.tab (reverse-stranded column, matching the counts matrix used everywhere else)
# needs: nothing beyond the raw counts / STAR ReadsPerGene.out.tab files
# outputs: spec23_QC_RUN1.csv, fig27_QC_RUN1_screenKOs.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))

KOSET <- c(ORDER, "NTC")
fab <- load_full_audited_bulk(); counts <- fab$cs; meta <- fab$meta
m <- meta[order(match(meta$Cell_Line, LINES), match(meta$Gene_Control, KOSET), meta$Sample_ID), ]
s <- rownames(m); n <- length(s)
cat(sprintf("QC set: %d samples (full audited universe, every status), KOs: %s\n", n, paste(KOSET, collapse = ", ")))
cat(sprintf("visually-flagged samples included here (%d): %s\n", sum(m$Sample_ID %in% VISUAL_OUTLIERS), paste(m$Sample_ID[m$Sample_ID %in% VISUAL_OUTLIERS], collapse = ", ")))

# read-assignment stats from STAR's ReadsPerGene.out.tab (col 4 = reverse-stranded)
map_one <- function(sid) {
  f <- file.path(PROC, "counts", paste0(sid, "_ReadsPerGene.out.tab"))
  d <- read.delim(f, header = FALSE, colClasses = c("character", "NULL", "NULL", "integer"))
  st <- setNames(d$V4[1:4], d$V1[1:4]); c(st, assigned = sum(d$V4[-(1:4)]))
}
MAP <- t(sapply(s, map_one)); total <- rowSums(MAP)
genes_detected <- colSums(counts[, s, drop = FALSE] > 0)
qc <- data.frame(sample = s, sample_id = m$Sample_ID, line = m$Cell_Line, KO = factor(m$Gene_Control, KOSET), total_reads = total, libsize_M = total / 1e6, genes_detected = genes_detected,
                  pct_assigned = 100 * MAP[, "assigned"] / total, pct_other = 100 * (MAP[, "N_ambiguous"] + MAP[, "N_noFeature"] + MAP[, "N_multimapping"]) / total,
                  pct_unmapped = 100 * MAP[, "N_unmapped"] / total, highlight = m$Sample_ID %in% VISUAL_OUTLIERS, row.names = NULL)
write.csv(qc, file.path(OUT, "spec23_QC_RUN1.csv"), row.names = FALSE)
options(width = 200); cat("\nper cell line summary:\n")
print(do.call(rbind, lapply(LINES, function(cl) { d <- qc[qc$line == cl, ]
  data.frame(line = cl, n = nrow(d), median_M_reads = round(median(d$libsize_M), 2), median_genes = median(d$genes_detected), median_pct_assigned = round(median(d$pct_assigned), 1)) })), row.names = FALSE)

## ---------------------------------------------------------------- figure (compact, log-scaled library size; by cell line or by KO)
cat_col <- c(assigned = "#1baf7a", other_unassigned = "#eda100", unmapped = "#7a7a7a")
line_col <- LINE_COL[as.character(qc$line)]
ko_col_a <- adjustcolor(KO_COL[as.character(qc$KO)], alpha.f = 0.45)      # same KO palette as the other figures, more transparent
edges <- which(m$Cell_Line[-1] != m$Cell_Line[-n]); bounds <- c(0, edges, n)
grp_key <- paste(m$Cell_Line, m$Gene_Control); ko_edges <- which(grp_key[-1] != grp_key[-n]); ko_bounds <- c(0, ko_edges, n)
ko_lab <- sub("^\\S+ ", "", grp_key[c(ko_bounds[-length(ko_bounds)] + 1)])

draw_qc <- function(bar_col, by_ko) function() {
  if (by_ko) layout(matrix(1:4, nrow = 4), heights = c(0.18, 1, 1, 1.35)) else layout(matrix(1:3, nrow = 3), heights = c(1, 1, 1.35))
  par(oma = c(0.1, 0, 2.1, 0.3), family = "sans")
  yr <- function() if (isTRUE(par("ylog"))) 10^par("usr")[3:4] else par("usr")[3:4]     # rect()/segments() need data-scale y even on a log axis
  sep <- function(bp) { y <- yr(); for (e in edges) segments((bp[e] + bp[e + 1]) / 2, y[1], (bp[e] + bp[e + 1]) / 2, y[2], col = "white", lwd = 1.8, xpd = FALSE) }

  hl_col <- ifelse(qc$highlight, "#d62728", INK)     # red border = visually-flagged outlier (VISUAL_OUTLIERS, repro_common.R)
  if (by_ko) { kk <- levels(qc$KO); par(mar = c(0, 0.4, 0, 0.4)); plot.new()
    legend("center", kk, fill = adjustcolor(KO_COL[kk], alpha.f = 0.45), border = NA, bty = "n", horiz = TRUE, cex = 0.42, x.intersp = 0.35, seg.len = 0.5) }

  par(mar = c(0.2, 3.6, 2.6, 0.4)); bp <- barplot(qc$libsize_M, col = bar_col, border = hl_col, lwd = ifelse(qc$highlight, 1.6, 0.4), space = 0, log = "y", ylab = "library size (M reads)", las = 1, cex.axis = 0.55, cex.lab = 0.6, xaxt = "n")
  abline(h = median(qc$libsize_M), col = INK2, lty = 2, lwd = 0.7); sep(bp)
  if (any(qc$highlight)) text(bp[qc$highlight], qc$libsize_M[qc$highlight], "*", col = "#d62728", font = 2, cex = 1.1, pos = 3, offset = 0.1, xpd = NA)
  for (i in seq_along(LINES)) mtext(LINES[i], side = 3, at = (bp[bounds[i] + 1] + bp[bounds[i + 1]]) / 2, line = 0.2, font = 2, cex = 0.5, col = LINE_COL[[LINES[i]]])
  mtext("A  Library size (log scale)", side = 3, adj = 0, line = 1.3, font = 2, cex = 0.62)

  par(mar = c(0.2, 3.6, 0.9, 0.4)); bp <- barplot(qc$genes_detected, col = bar_col, border = hl_col, lwd = ifelse(qc$highlight, 1.6, 0.4), space = 0, ylab = "genes detected", las = 1, cex.axis = 0.55, cex.lab = 0.6, xaxt = "n")
  abline(h = median(qc$genes_detected), col = INK2, lty = 2, lwd = 0.7); sep(bp)
  mtext("B  Gene complexity", side = 3, adj = 0, line = 0.15, font = 2, cex = 0.62)

  par(mar = c(4.4, 3.6, 0.9, 0.4)); M <- t(as.matrix(qc[, c("pct_assigned", "pct_other", "pct_unmapped")])); rownames(M) <- names(cat_col)
  bp <- barplot(M, col = cat_col, border = hl_col, lwd = ifelse(qc$highlight, 1.6, 0.4), space = 0, ylab = "% reads", las = 1, cex.axis = 0.55, cex.lab = 0.6, xaxt = "n")
  sep(bp)
  if (any(qc$highlight)) text(bp[qc$highlight], 100, "*", col = "#d62728", font = 2, cex = 1.1, pos = 3, offset = 0.1, xpd = NA)
  text(bp[ko_bounds[-length(ko_bounds)] + 1] + diff(bp[1:2]) * (diff(ko_bounds) - 1) / 2, par("usr")[3] - 1.5, ko_lab, srt = 90, adj = 1, cex = 0.32, xpd = NA, col = INK2)
  legend(par("usr")[2] * 0.005, par("usr")[3] - 11, c("assigned", "other unassigned", "unmapped"), fill = cat_col, border = NA, bty = "n", horiz = TRUE, cex = 0.4, xpd = NA, x.intersp = 0.4, seg.len = 0.6)
  mtext("C  Read assignment", side = 3, adj = 0, line = 0.15, font = 2, cex = 0.62)

  mtext(sprintf("Bulk RNA-seq QC, full audited universe (%d samples, CROP-seq-matched KOs only)", n), outer = TRUE, side = 3, adj = 0, line = 1.05, font = 2, cex = 0.6)
  mtext(sprintf("SNAI2, RUNX1/2, NR4A2, WT excluded; colour = %s; dashed = median; red border/* = visually-flagged outlier (VISUAL_OUTLIERS); white = cell line", if (by_ko) "KO (grey = NTC)" else "cell line"),
        outer = TRUE, side = 3, adj = 0, line = 0.15, cex = 0.38, col = INK2)
}
save_fig("fig27_QC_RUN1_screenKOs", 5.25, 5.6, draw_qc(line_col, FALSE))
save_fig("fig27_QC_RUN1_screenKOs_byKO", 5.25, 5.75, draw_qc(ko_col_a, TRUE))
