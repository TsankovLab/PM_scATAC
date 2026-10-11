# Spec 43: for each of the 8 shared KOs, take every cNMF module (k=30, all 30) whose usage is LOWER in that KO than in NTC in the
# screen, then check whether that same "down vs NTC" direction holds in bulk RNA-seq (AddModuleScore), both pooled across lines and
# per line individually. This is the per-KO complement to spec39/40's per-module sweep: instead of asking "which KO is highest/
# lowest for module X", it asks "of the modules this KO pushes DOWN in the screen, how many also read down in bulk?"
# needs: CRISPR_cropseq_analysis/.../cNMF_crispr_20_30_vf5000/cnmf/cnmf.usages.k_30.dt_0_3.consensus.txt,
#        cnmf.gene_spectra_score.k_30.dt_0_3.txt, spec_common.R (load_screen()), repro_common.R (load_full_audited_bulk())
# outputs: spec43_screen_down_modules_bulk_agreement.csv, spec43_per_KO_summary.csv,
#          fig3A_screen_down_modules_bulk_agreement.{png,pdf}, fig3A_per_line_agreement.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat))
CNMF_D <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CRISPR_cropseq_analysis/_cellranger_raw_Filter_400_800_25/no_harmony/cNMF/cNMF_crispr_20_30_vf5000/cnmf"
K <- "30"; MODULES <- as.character(1:30)

## ---------------------------------------------------------------- screen: per-module top50 genes + per-KO usage (same as spec39)
gss <- read.delim(file.path(CNMF_D, sprintf("cnmf.gene_spectra_score.k_%s.dt_0_3.txt", K)), row.names = 1, check.names = FALSE)
usages <- read.delim(file.path(CNMF_D, sprintf("cnmf.usages.k_%s.dt_0_3.consensus.txt", K)), row.names = 1, check.names = FALSE)
srt <- load_screen(); grp_scr <- setNames(as.character(srt$merged_call), colnames(srt))
common <- intersect(rownames(usages), names(grp_scr)); grp_scr <- grp_scr[common]; usages <- usages[common, ]

top_genes_list <- lapply(MODULES, function(m) { x <- sort(setNames(as.numeric(gss[m, ]), colnames(gss)), decreasing = TRUE); names(x)[1:50] })
names(top_genes_list) <- paste0("m", MODULES)

GROUPS <- c("NTC", ORDER)
screen_usage <- sapply(MODULES, function(m) sapply(GROUPS, function(grp) mean(usages[grp_scr == grp, m])))
rownames(screen_usage) <- GROUPS

## ---------------------------------------------------------------- significance of "lower than NTC" per KO x module: one-sided Wilcoxon (KO < NTC), BH-adjusted across all 8x30 tests
ntc_idx <- which(grp_scr == "NTC")
screen_p <- sapply(MODULES, function(m) sapply(ORDER, function(k) {
  ko_idx <- which(grp_scr == k)
  wilcox.test(usages[ko_idx, m], usages[ntc_idx, m], alternative = "less", exact = FALSE)$p.value
}))
rownames(screen_p) <- ORDER
screen_padj <- matrix(p.adjust(as.vector(screen_p), method = "BH"), nrow = nrow(screen_p), dimnames = dimnames(screen_p))
cat(sprintf("screen-down significance test (one-sided Wilcoxon, KO < NTC, BH-adjusted over %d tests): %d/%d KO x module pairs significant at padj<0.05\n",
            length(screen_padj), sum(screen_padj < 0.05), length(screen_padj)))

## ---------------------------------------------------------------- bulk: score all 30 modules (cell-level, so we can split by line)
fab <- load_full_audited_bulk(); meta <- fab$meta; lg <- fab$lg
options(Seurat.object.assay.version = "v3")   # AddModuleScore errors ("subscript out of bounds") on the default Assay5 class here; classic Assay works
srt_bulk <- CreateSeuratObject(counts = as.matrix(fab$cs)); srt_bulk <- SetAssayData(srt_bulk, slot = "data", new.data = as.matrix(lg))
srt_bulk <- AddModuleScore(srt_bulk, features = top_genes_list, name = "m", seed = 1)
bulk_score_cells <- as.matrix(srt_bulk@meta.data[, paste0("m", seq_along(MODULES))]); colnames(bulk_score_cells) <- MODULES
bulk_score_cells <- bulk_score_cells[rownames(meta), ]

bulk_pooled <- sapply(MODULES, function(m) sapply(GROUPS, function(grp) mean(bulk_score_cells[meta$Gene_Control == grp, m])))
rownames(bulk_pooled) <- GROUPS
bulk_per_line <- lapply(LINES, function(cl) {
  mm <- meta[meta$Cell_Line == cl, ]
  sapply(MODULES, function(m) sapply(GROUPS, function(grp) { s <- rownames(mm)[mm$Gene_Control == grp]; if (length(s)) mean(bulk_score_cells[s, m]) else NA }))
})
names(bulk_per_line) <- LINES
for (cl in LINES) rownames(bulk_per_line[[cl]]) <- GROUPS

## ---------------------------------------------------------------- for each KO, which screen modules are SIGNIFICANTLY down vs NTC (padj<0.05)? check bulk correspondence
res <- do.call(rbind, lapply(ORDER, function(k) {
  down_mods <- MODULES[screen_padj[k, ] < 0.05]
  if (!length(down_mods)) return(NULL)
  do.call(rbind, lapply(down_mods, function(m) {
    line_down <- sapply(LINES, function(cl) bulk_per_line[[cl]][k, m] < bulk_per_line[[cl]]["NTC", m])
    data.frame(KO = k, module = m, screen_NTC = screen_usage["NTC", m], screen_KO = screen_usage[k, m], screen_padj = screen_padj[k, m],
               bulk_NTC_pooled = bulk_pooled["NTC", m], bulk_KO_pooled = bulk_pooled[k, m],
               bulk_down_pooled = bulk_pooled[k, m] < bulk_pooled["NTC", m],
               n_lines_down = sum(line_down, na.rm = TRUE), n_lines_tested = sum(!is.na(line_down)),
               t(setNames(line_down, paste0("down_", LINES))))
  }))
}))
write.csv(res, file.path(OUT, "spec43_screen_down_modules_bulk_agreement.csv"), row.names = FALSE)

summ <- do.call(rbind, lapply(ORDER, function(k) {
  d <- res[res$KO == k, ]
  data.frame(KO = k, n_screen_down_modules = nrow(d), n_bulk_down_pooled = sum(d$bulk_down_pooled),
             frac_pooled = mean(d$bulk_down_pooled), mean_frac_lines_down = mean(d$n_lines_down / d$n_lines_tested))
}))
write.csv(summ, file.path(OUT, "spec43_per_KO_summary.csv"), row.names = FALSE)

options(width = 200)
cat("=== per KO: of the cNMF modules that read LOWER than NTC in the SCREEN, how many are also lower in BULK? ===\n")
print(format(summ, digits = 3), row.names = FALSE)
cat(sprintf("\noverall: %d screen-down module-KO pairs tested, %d (%.0f%%) also down in pooled bulk [chance = 50%%]\n",
            nrow(res), sum(res$bulk_down_pooled), 100 * mean(res$bulk_down_pooled)))
bt <- binom.test(sum(res$bulk_down_pooled), nrow(res), p = 0.5)
cat(sprintf("binomial test vs 50%% chance: p = %.4g\n", bt$p.value))
cat("\n=== detail ===\n"); print(format(res, digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- fig: per-KO bar, fraction of screen-down modules also down in pooled bulk
save_fig("fig3A_screen_down_modules_bulk_agreement", 8.5, 6, function() {
  ord <- summ[order(-summ$frac_pooled), ]
  cols <- KO_COL[ord$KO]
  par(mar = c(6, 4.8, 3.6, 1), family = "sans")
  bp <- barplot(ord$frac_pooled, names.arg = ord$KO, col = adjustcolor(cols, 0.85), border = INK2, las = 2, cex.names = 0.85,
                ylab = "fraction also down in bulk (pooled)", ylim = c(0, 1))
  abline(h = 0.5, col = INK, lty = 2)
  text(bp, ord$frac_pooled, sprintf("%d/%d", ord$n_bulk_down_pooled, ord$n_screen_down_modules), pos = 3, cex = 0.7, xpd = NA)
  mtext("Significant screen-down cNMF modules (padj<0.05): also down in bulk?", side = 3, adj = 0, line = 2.2, font = 2, cex = 0.95)
  mtext("dashed line = 50% chance; numbers = modules agreeing / modules tested for that KO", side = 3, adj = 0, line = 0.9, cex = 0.6, col = INK2)
})

## ---------------------------------------------------------------- fig: per-line breakdown
save_fig("fig3A_per_line_agreement", 9, 6, function() {
  frac_line <- sapply(ORDER, function(k) { d <- res[res$KO == k, ]; sapply(LINES, function(cl) mean(d[[make.names(paste0("down_", cl))]], na.rm = TRUE)) })
  n_r <- length(LINES); n_c <- length(ORDER)
  par(mar = c(6.5, 8, 3.4, 1), family = "sans")
  image(seq_len(n_c), seq_len(n_r), t(frac_line)[, n_r:1], col = colorRampPalette(c("#2a78d6", "#f0efec", "#eb6834"))(101), zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(n_r)) for (j in seq_len(n_c)) text(j, n_r - i + 1, sprintf("%.2f", frac_line[i, j]), cex = 0.75, col = INK)
  axis(1, at = seq_len(n_c), labels = ORDER, las = 2, tick = FALSE, cex.axis = 0.85, col.axis = INK2)
  axis(2, at = seq_len(n_r), labels = rev(LINES), las = 2, tick = FALSE, cex.axis = 0.8, col.axis = INK2)
  box(col = MUTED, lwd = 0.6)
  mtext("Per-line: fraction of each KO's SIGNIFICANT screen-down modules also down in that line's bulk", side = 3, adj = 0, line = 2.0, font = 2, cex = 0.9)
  mtext("0.5 = chance; cell = fraction of that KO's screen-identified down-modules confirmed down in that line", side = 3, adj = 0, line = 0.8, cex = 0.58, col = INK2)
  draw_colorbar(c(0.06, 0.30), c(0, 1), "fraction confirmed down", cex = 0.55)
})

## ---------------------------------------------------------------- fig: per-KO binary heatmap, rows = Screen + 4 lines, columns = that KO's screen-down modules
save_fig("fig3B_screen_down_modules_binary_per_KO", 15, 9, function() {
  layout(matrix(1:8, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3, 0), family = "sans")
  SRC <- c("Screen", LINES)
  for (k in ORDER) {
    d <- res[res$KO == k, ]; d <- d[order(as.integer(d$module)), ]
    n_c <- nrow(d); n_r <- length(SRC)
    M <- matrix(NA, nrow = n_r, ncol = n_c, dimnames = list(SRC, d$module))
    M["Screen", ] <- 1
    for (cl in LINES) M[cl, ] <- as.integer(d[[make.names(paste0("down_", cl))]])
    par(mar = c(4.5, 6.5, 3, 1))
    Mcol <- ifelse(is.na(M), 0.5, M)
    image(seq_len(n_c), seq_len(n_r), t(Mcol)[, n_r:1], col = c("#eb6834", "grey85", "#2a78d6"), zlim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(n_r)) for (j in seq_len(n_c)) if (!is.na(M[i, j])) text(j, n_r - i + 1, ifelse(M[i, j] == 1, "down", "up"), cex = 0.55, col = "white")
    axis(1, at = seq_len(n_c), labels = colnames(M), las = 2, tick = FALSE, cex.axis = 0.65, col.axis = INK2)
    axis(2, at = seq_len(n_r), labels = rev(SRC), las = 2, tick = FALSE, cex.axis = 0.75, col.axis = INK2)
    box(col = MUTED, lwd = 0.6)
    mtext(sprintf("%s (n=%d sig. down in screen)", k, n_c), side = 3, adj = 0, line = 0.6, font = 2, cex = 0.85, col = KO_COL[k])
  }
  mtext("For each KO's SIGNIFICANTLY screen-down cNMF modules (padj<0.05): binary down/up call per line", outer = TRUE, side = 3, adj = 0, line = 1.6, font = 2, cex = 1.05)
  mtext("blue = down vs NTC (matches screen direction); orange = up vs NTC (does not match); columns = cNMF module number; Screen row is down by construction",
        outer = TRUE, side = 3, adj = 0, line = 0.4, cex = 0.6, col = INK2)
})
