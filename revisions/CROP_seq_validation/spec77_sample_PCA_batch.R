# Spec 77: sample distribution for batch effects -- all 76 audited bulk samples (load_full_audited_bulk(): the 8 screen KOs + NTC, incl. undocumented/mixed samples).
# A: PCA of all samples (top 2000 variable genes, log2 CP10K), coloured by the GENOTYPE-corrected cell line, shape = library index group (i5 index; the sheet has no
#    run/date column, so the i5 index group is the only batch proxy), non-clean samples labelled.
# B-E: PCA within each line (top 2000 variable genes of that line), coloured by i5 index group, square = NTC, every sample labelled (KO_rep), red label = flagged
#    in the audit (swapped / undocumented / mixed / strict-QC outlier). Also reports how much of each within-line PC is explained by the i5 group and by % reads assigned.
# needs: repro_common.R, profiles.R, spec71_sample_audit.csv (optional, for the QC flags)
# outputs: spec77_sample_PCA_batch.csv, fig7F_sample_PCA_batch.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
a <- load_full_audited_bulk(); meta <- a$meta; lg <- a$lg
sheet <- read.delim(file.path(PROC, "sample_sheet/sample_metadata_RUN1.tsv"))
idx <- rbind(sheet[, c("fastq_sample_id", "i5_index_num")], data.frame(fastq_sample_id = H28_NTC_NATIVE$fastq_sample_id, i5_index_num = H28_NTC_NATIVE$i5_index_num))
meta$i5 <- idx$i5_index_num[match(rownames(meta), idx$fastq_sample_id)]; meta$i5g <- ifelse(is.na(meta$i5), "unknown", paste0("i5 #", meta$i5))
aud_f <- file.path(OUT, "spec71_sample_audit.csv")
if (file.exists(aud_f)) { aud <- read.csv(aud_f); meta$qc_outlier <- aud$qc_outlier[match(rownames(meta), aud$sample)]; meta$pct_assigned <- aud$pct_assigned[match(rownames(meta), aud$sample)] } else { meta$qc_outlier <- FALSE; meta$pct_assigned <- NA }
meta$flagged <- meta$flag != "clean" | (!is.na(meta$qc_outlier) & meta$qc_outlier)
I5_LEV <- c("i5 #2", "i5 #3", "i5 #4", "i5 #17", "i5 #18", "unknown")
I5_COL <- setNames(c("#2a78d6", "#1baf7a", "#eda100", "#b3400f", "#e87ba4", "grey55"), I5_LEV); I5_PCH <- setNames(c(16, 17, 15, 18, 8, 4), I5_LEV)
cat("i5 index group by genotype line:\n"); print(table(meta$Cell_Line, factor(meta$i5g, I5_LEV)))

pca <- function(M, n = 2000) { v <- matrixStats::rowVars(M); M <- M[order(-v)[seq_len(min(n, nrow(M)))], ]; p <- prcomp(t(M - rowMeans(M))); list(x = p$x[, 1:4], ve = 100 * p$sdev[1:4]^2 / sum(p$sdev^2)) }
keep <- rowMeans(a$cp10k) >= 0.1
G <- pca(lg[keep, ])
P <- lapply(setNames(LINES, LINES), function(cl) { s <- rownames(meta)[meta$Cell_Line == cl]; ok <- rowMeans(a$cp10k[, s]) >= 0.1; pca(lg[ok, s]) })

r2 <- function(y, x) { if (length(unique(x[!is.na(x)])) < 2) return(NA); summary(lm(y ~ x))$r.squared }
tab <- do.call(rbind, lapply(LINES, function(cl) { s <- rownames(P[[cl]]$x); m <- meta[s, ]
  do.call(rbind, lapply(1:3, function(k) data.frame(line = cl, PC = paste0("PC", k), var_pct = P[[cl]]$ve[k], R2_i5_group = r2(P[[cl]]$x[, k], factor(m$i5g)), R2_pct_assigned = r2(P[[cl]]$x[, k], m$pct_assigned),
                                                     R2_NTC_vs_KO = r2(P[[cl]]$x[, k], m$Gene_Control == "NTC")))) }))
write.csv(tab, file.path(OUT, "spec77_sample_PCA_batch.csv"), row.names = FALSE)
options(width = 200); cat("\nvariance of each within-line PC explained by i5 index group / % reads assigned / NTC vs KO:\n"); print(format(tab, digits = 2), row.names = FALSE)

PFX <- c("NCI-H2052" = "H2052", "MSTO-211H" = "MSTO", "NCI-H2452" = "H2452", "NCI-H28" = "H28")
short <- function(s) { own <- startsWith(s, paste0(PFX[meta[s, "Cell_Line"]], "_")); ifelse(own, sub("^[^_]+_", "", s), sub("^MIXED_", "MIX_", s)) }   # full name when the name prefix is not the genotype line
save_fig("fig7F_sample_PCA_batch", 15, 9.2, function() {
  layout(matrix(c(1, 1, 2, 3, 1, 1, 4, 5), nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.2, 0), family = "sans")
  ## A: all samples
  par(mar = c(4.2, 4.2, 2.6, 1)); x <- G$x; m <- meta[rownames(x), ]
  plot(x[, 1], x[, 2], type = "n", xlab = sprintf("PC1 (%.0f%%)", G$ve[1]), ylab = sprintf("PC2 (%.0f%%)", G$ve[2]))
  points(x[, 1], x[, 2], pch = I5_PCH[m$i5g], col = LINE_COL[m$Cell_Line], cex = 1.5, lwd = 1.4)
  lab <- m$flagged | m$Gene_Control == "NTC"; text(x[lab, 1], x[lab, 2], rownames(x)[lab], pos = 3, cex = 0.55, col = ifelse(m$flagged[lab], "#d62728", INK2), xpd = NA)
  legend("topright", legend = LINES, col = LINE_COL[LINES], pch = 16, bty = "n", cex = 0.75, title = "genotype line")
  legend("bottomright", legend = I5_LEV, pch = I5_PCH, col = "grey30", bty = "n", cex = 0.7, title = "library index group")
  mtext("A  All 76 samples (colour = line by genotype, shape = i5 index group)", side = 3, adj = 0, line = 0.6, font = 2, cex = 0.85)
  ## B-E: within line
  for (cl in LINES) { x <- P[[cl]]$x; m <- meta[rownames(x), ]; par(mar = c(3.8, 3.8, 2.4, 0.8))
    plot(x[, 1], x[, 2], type = "n", xlab = sprintf("PC1 (%.0f%%)", P[[cl]]$ve[1]), ylab = sprintf("PC2 (%.0f%%)", P[[cl]]$ve[2]), cex.lab = 0.8, cex.axis = 0.75)
    points(x[, 1], x[, 2], pch = ifelse(m$Gene_Control == "NTC", 15, 16), col = I5_COL[m$i5g], cex = 1.4)
    text(x[, 1], x[, 2], short(rownames(x)), pos = 3, cex = 0.5, col = ifelse(m$flagged, "#d62728", INK2), font = ifelse(m$Gene_Control == "NTC", 2, 1), xpd = NA)
    li <- unique(m$i5g); legend("bottomright", legend = li, col = I5_COL[li], pch = 16, bty = "n", cex = 0.65)
    t1 <- tab[tab$line == cl & tab$PC == "PC1", ]
    mtext(sprintf("%s  (PC1 R2: i5 %.2f, %% assigned %.2f, NTC %.2f)", cl, t1$R2_i5_group, t1$R2_pct_assigned, t1$R2_NTC_vs_KO), side = 3, adj = 0, line = 0.5, font = 2, cex = 0.66, col = LINE_COL[[cl]]) }
  mtext("Bulk sample distribution: cell line and library batch (i5 index group)", outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 1.0)
  mtext("PCA on log2 CP10K, top 2000 variable genes (genes with mean CP10K >= 0.1); B-E: within-line PCA, colour = i5 index group, square = NTC, red label = flagged in the audit (swapped / undocumented / mixed / strict-QC outlier)", outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.58, col = INK2)
})

## ---------------------------------------------------------------- fig7G: per-sample expression distributions (log2 CP10K of genes expressed overall: mean CP10K >= 0.1, i.e. 10 CPM), one box per sample,
## grouped by genotype line (NTC first, then KO order), coloured by line; strip below = i5 index group; red label = flagged in the audit; dashed = line median of sample medians.
ord <- unlist(lapply(LINES, function(cl) { s <- rownames(meta)[meta$Cell_Line == cl]; s[order(factor(meta[s, "Gene_Control"], levels = c("NTC", ORDER)), s)] }))
ex <- rowMeans(a$cp10k) >= 0.1; E <- lg[ex, ord]
st <- data.frame(sample = ord, line = meta[ord, "Cell_Line"], i5 = meta[ord, "i5g"], median = apply(E, 2, median), q25 = apply(E, 2, quantile, 0.25), q75 = apply(E, 2, quantile, 0.75), pct_zero = colMeans(E == 0))
st$median_dev_from_line <- st$median - ave(st$median, st$line, FUN = median)
write.csv(st, file.path(OUT, "spec77_sample_expression_distribution.csv"), row.names = FALSE)
cat(sprintf("\n%d genes with mean CP10K >= 0.1; samples whose median log2 CP10K deviates > 0.15 from their line:\n", sum(ex))); print(format(st[abs(st$median_dev_from_line) > 0.15, ], digits = 3), row.names = FALSE)
save_fig("fig7G_sample_expression_boxplots", 16, 5.8, function() {
  n <- length(ord); par(mar = c(9.5, 4.2, 3.6, 0.8), family = "sans")
  bx <- boxplot(E, at = seq_len(n), outline = FALSE, col = adjustcolor(LINE_COL[st$line], 0.55), border = INK2, xaxt = "n", ylab = "log2 CP10K (expressed genes)", lwd = 0.6, boxwex = 0.75, xlim = c(0.5, n + 0.5))
  edges <- which(st$line[-1] != st$line[-n]); abline(v = edges + 0.5, col = "grey70", lwd = 1.2)
  for (cl in LINES) { i <- which(st$line == cl); segments(min(i) - 0.4, median(st$median[i]), max(i) + 0.4, median(st$median[i]), lty = 2, col = LINE_COL[[cl]], lwd = 1.2)
    mtext(cl, side = 3, at = mean(i), line = 0.2, font = 2, cex = 0.8, col = LINE_COL[[cl]]) }
  flag <- meta[ord, "flagged"]; lab <- sub("^MIXED_", "MIX_", ord)
  text(seq_len(n), par("usr")[3] - 0.05 * diff(par("usr")[3:4]), lab, srt = 90, adj = 1, xpd = NA, cex = 0.52, col = ifelse(flag, "#d62728", INK2), font = ifelse(meta[ord, "Gene_Control"] == "NTC", 2, 1))
  y0 <- par("usr")[3]; h <- 0.035 * diff(par("usr")[3:4]); rect(seq_len(n) - 0.5, y0 - h, seq_len(n) + 0.5, y0, col = I5_COL[st$i5], border = NA, xpd = NA)
  li <- I5_LEV[I5_LEV %in% st$i5]; legend("topright", legend = li, fill = I5_COL[li], border = NA, bty = "n", cex = 0.65, title = "i5 index group (strip)", ncol = 3)
  mtext("Per-sample expression distribution, all 76 audited bulk samples, grouped by genotype line", side = 3, adj = 0, line = 2.2, font = 2, cex = 0.95)
  mtext(sprintf("%d genes with mean CP10K >= 0.1 (10 CPM); box = IQR, whiskers 1.5 IQR; dashed = line median; bold label = NTC, red = flagged in the audit (swapped / undocumented / mixed / strict-QC outlier)", sum(ex)), side = 3, adj = 0, line = 1.2, cex = 0.58, col = INK2)
})
