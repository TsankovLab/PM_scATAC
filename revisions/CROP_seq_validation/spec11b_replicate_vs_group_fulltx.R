## ──────────────────────────────────────────────────────────────────────
## PARAMETERS
## ──────────────────────────────────────────────────────────────────────
# Full-transcriptome version of fig13c (spec11_replicate_vs_group.R): replicate vs KO-group structure in the RUN1 bulk samples, within each
# cell line. spec11 uses the top 3000 variable genes per line; here EVERY gene expressed in the line (mean CP10K >= MINEXP) is used, with no
# variance-based selection. Spearman r between every pair of samples on line-centred log2 CP10K (gene-wise mean of the line's samples removed);
# same-KO pairs (incl. NTC) vs different-KO pairs. AUC = P(random same-KO pair has a higher r than random different-KO pair).
OUT_DIR  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CROPseq_validation_bulkRNA/processed/de/perturbation_panel/figs_repro_R"
COMMON   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/CROP_seq_validation/repro_common.R"
MINEXP   <- 0.1       # genes kept per line: mean CP10K >= 0.1 in that line (~350 reads at 35M depth); same filter as spec11, no top-N cut
SEED     <- 5
FIG_NAME <- "fig13c_within_line_spearman_centred_compact_fulltx"

## ──────────────────────────────────────────────────────────────────────
## LIBRARIES
## ──────────────────────────────────────────────────────────────────────
Sys.setenv(CROP_OUT = OUT_DIR)   # repro_common.R reads CROP_OUT into OUT (used by save_fig)
source(COMMON)

## ──────────────────────────────────────────────────────────────────────
## DATA
## ──────────────────────────────────────────────────────────────────────
set.seed(SEED)
b <- load_bulk(); cp10k <- b$cp10k; meta <- b$meta; lg <- log2(cp10k + 1)

## ──────────────────────────────────────────────────────────────────────
## FUNCTIONS
## ──────────────────────────────────────────────────────────────────────
auc <- function(w, o) { if (!length(w) || !length(o)) return(NA_real_); r <- rank(c(w, o)); (sum(r[seq_along(w)]) - length(w) * (length(w) + 1) / 2) / (length(w) * length(o)) }
line_genes <- function(cl) { s <- rownames(meta)[meta$Cell_Line == cl]; rownames(cp10k)[rowMeans(cp10k[, s]) >= MINEXP] }   # all expressed genes, no variance cut
centre_line <- function(y) y - rowMeans(y)
stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

## ──────────────────────────────────────────────────────────────────────
## Figure: fig13c_within_line_spearman_centred_compact_fulltx.pdf
## ──────────────────────────────────────────────────────────────────────
within <- list(); pairs_store <- list(); n_genes <- setNames(integer(length(LINES)), LINES)
for (cl in LINES) {
  s <- rownames(meta)[meta$Cell_Line == cl]; g <- line_genes(cl); n_genes[cl] <- length(g); y <- centre_line(lg[g, s]); grp <- meta[s, "Gene_Control"]
  C <- cor(y, method = "spearman"); ut <- which(upper.tri(C), arr.ind = TRUE)
  d <- data.frame(line = cl, method = "spearman", centred = TRUE, a = s[ut[, 1]], b = s[ut[, 2]], r = C[ut], same = grp[ut[, 1]] == grp[ut[, 2]])
  d$ko <- ifelse(d$same, as.character(grp[ut[, 1]]), NA); pairs_store[[cl]] <- d
  w <- d$r[d$same]; o <- d$r[!d$same]
  within[[cl]] <- data.frame(line = cl, method = "spearman", centred = TRUE, n_genes = length(g), n_samples = length(s), n_same_pairs = length(w), n_diff_pairs = length(o),
                             median_same = median(w), median_diff = median(o), delta_median = median(w) - median(o), AUC = auc(w, o),
                             p_wilcox_greater = wilcox.test(w, o, alternative = "greater")$p.value)
}
within <- do.call(rbind, within); pairs_all <- do.call(rbind, pairs_store)
write.csv(within, file.path(OUT_DIR, "spec_repl_within_line_fulltx.csv"), row.names = FALSE)
write.csv(pairs_all, file.path(OUT_DIR, "spec_repl_within_line_fulltx_pairs.csv"), row.names = FALSE)
options(width = 200); cat("\n=== within cell line (full transcriptome): replicates of the same KO vs different KOs ===\n"); print(format(within, digits = 3), row.names = FALSE)

save_fig(FIG_NAME, 4.9, 3.8, res = 300, fun = function() {
  par(mar = c(2.4, 3.5, 4.2, 0.5), mgp = c(2.2, 0.6, 0), family = "sans")
  k <- length(LINES)
  same <- lapply(LINES, function(cl) pairs_all$r[pairs_all$line == cl & pairs_all$same]); diffv <- lapply(LINES, function(cl) pairs_all$r[pairs_all$line == cl & !pairs_all$same])
  pw <- within$p_wilcox_greater[match(LINES, within$line)]
  at_s <- seq_len(k) * 3 - 1; at_d <- seq_len(k) * 3; yl <- range(unlist(c(same, diffv))) + c(-0.03, 0.28)
  plot(NA, xlim = c(0.3, 3 * k + 0.7), ylim = yl, xaxt = "n", xlab = "", ylab = "Spearman r between samples", bty = "n", las = 1, cex.axis = 0.7, cex.lab = 0.75)
  abline(h = 0, col = GRID, lwd = 0.8)
  boxplot(diffv, at = at_d, add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.85, col = adjustcolor(GREY, 0.9), border = INK2, medcol = INK, medlwd = 1.6, whisklty = 1, staplelty = 1, lwd = 0.8)
  boxplot(same, at = at_s, add = TRUE, axes = FALSE, outline = FALSE, boxwex = 0.85, col = adjustcolor(LINE_COL[LINES], 0.85), border = INK2, medcol = INK, medlwd = 1.6, whisklty = 1, staplelty = 1, lwd = 0.8)
  for (i in seq_len(k)) points(at_d[i] + runif(length(diffv[[i]]), -0.25, 0.25), diffv[[i]], pch = 16, col = adjustcolor(INK2, 0.3), cex = 0.3)
  for (i in seq_len(k)) points(at_s[i] + runif(length(same[[i]]), -0.22, 0.22), same[[i]], pch = 21, bg = LINE_COL[LINES][i], col = INK, cex = 0.5, lwd = 0.4)
  axis(1, at = (at_s + at_d) / 2, labels = sub("NCI-", "", LINES), tick = FALSE, cex.axis = 0.75, font.axis = 2, line = -0.6)
  a <- within$AUC[match(LINES, within$line)]
  off <- diff(yl) * 0.07
  text((at_s + at_d) / 2, yl[2], sprintf("AUC %.2f", a), cex = 0.65, col = INK2)
  text((at_s + at_d) / 2, yl[2] - off, stars(pw), cex = 0.75, col = INK, font = ifelse(stars(pw) == "ns", 1, 2))
  mtext("Replicates cluster by KO within each cell line", side = 3, adj = 0, line = 2.9, font = 2, cex = 0.8)
  mtext(sprintf("Spearman r over all expressed genes (mean CP10K ≥ %g per line; n = %s), centred on the line mean", MINEXP, paste(format(range(n_genes), big.mark = ","), collapse = "–")),
        side = 3, adj = 0, line = 1.95, cex = 0.52, col = INK2)
  mtext("coloured = same-KO pairs, grey = different-KO pairs (dots = sample pairs)", side = 3, adj = 0, line = 1.15, cex = 0.55, col = INK2)
  mtext("*** p<0.001, ** p<0.01, * p<0.05, ns = not sig. (one-sided Wilcoxon, same-KO > diff-KO)", side = 3, adj = 0, line = 0.35, cex = 0.5, col = INK2)
})
