# Spec 45: fully self-contained recompute of the 20 "Cm" modules' screen-vs-bulk agreement, using the ACTUAL source gene sets
# (git_repo/files/scRNA_PM_modules.rds -- cNMF gene programs from the primary PM tumor scRNA-seq atlas, Giotti et al 2024) and the
# exact method found in the original generating code (git_repo/tumor_analysis/CRISPR_cropseq.R + utils/useful_functions.R's
# ModScoreCor(), a thin AddModuleScore wrapper; top-50 genes per module; per-group MEDIAN, not mean; diff = NTC - KO median).
# This replaces spec44's dependence on the legacy screen_cm_diff_matrix.csv / cm_sample_scores_AddModuleScore_bulk_RUN1.csv with a
# from-scratch computation on BOTH sides (screen via load_screen(), bulk via load_bulk()), enabling: (1) validation against the
# legacy files, (2) a true RAW (non-NTC-differenced) screen-vs-bulk comparison (impossible before -- no raw screen file existed).
# needs: git_repo/files/scRNA_PM_modules.rds, spec_common.R (load_screen()), repro_common.R (load_bulk())
# outputs: spec45_20Cm_recompute_vs_legacy_validation.csv, spec45_20Cm_diff_agreement.csv, spec45_20Cm_raw_agreement.csv,
#          fig3D_20Cm_recompute_validation.{png,pdf}, fig3D_20Cm_diff_vs_raw_correlation.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat))
PP <- file.path(PROC, "de", "perturbation_panel")
# DIFF_DIR=NTC_KO flips the sign convention of the diff-based comparison (NTC - KO, the legacy convention) on BOTH sides; default KO_NTC.
DIFF_DIR <- Sys.getenv("DIFF_DIR", "KO_NTC"); SGN <- if (DIFF_DIR == "NTC_KO") -1 else 1; DSUF <- if (DIFF_DIR == "NTC_KO") "_NTCminusKO" else ""; DLAB <- if (DIFF_DIR == "NTC_KO") "NTC-KO" else "KO-NTC"

cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds")
top50 <- lapply(cm_genes, head, 50)   # exact original convention: head(x, 50)
CM <- names(top50)
cat(sprintf("loaded %d Cm modules (top-50 genes each) from scRNA_PM_modules.rds\n", length(CM)))

options(Seurat.object.assay.version = "v3")   # AddModuleScore errors ("subscript out of bounds") on the default Assay5 class here; classic Assay works

## ---------------------------------------------------------------- SCREEN: recompute raw per-cell Cm scores from scratch
srt <- load_screen()
srt <- AddModuleScore(srt, features = top50, name = "Cm_", seed = 1)
screen_cell_scores <- as.matrix(srt@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(screen_cell_scores) <- CM
grp_scr <- as.character(srt$merged_call)

screen_groups <- unique(grp_scr)
screen_median <- t(sapply(screen_groups, function(g) apply(screen_cell_scores[grp_scr == g, , drop = FALSE], 2, median)))
rownames(screen_median) <- screen_groups
screen_diff_recompute <- -sweep(screen_median, 2, screen_median["NTC", ], "-")   # NTC - group, the original CRISPR_cropseq.R sign convention

## ---------------------------------------------------------------- validate against the legacy screen_cm_diff_matrix.csv (after its own sign flip back to "NTC-KO" for this comparison)
legacy_screen_diff <- read.csv(file.path(PP, "screen_cm_diff_matrix.csv"), row.names = 1, check.names = FALSE)  # already stored as NTC - KO
common_ko <- intersect(rownames(legacy_screen_diff), ORDER)
valid <- do.call(rbind, lapply(common_ko, function(k) {
  a <- as.numeric(legacy_screen_diff[k, CM]); b <- screen_diff_recompute[k, CM]
  data.frame(KO = k, pearson_r = cor(a, b), rmse = sqrt(mean((a - b)^2)))
}))
write.csv(valid, file.path(OUT, "spec45_20Cm_recompute_vs_legacy_validation.csv"), row.names = FALSE)
cat("\n=== validation: our from-scratch recomputed screen diff (NTC-KO median) vs the legacy screen_cm_diff_matrix.csv ===\n")
print(format(valid, digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- BULK: recompute raw per-sample Cm scores from scratch (same top-50 gene sets)
b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1)
srt_bulk <- CreateSeuratObject(counts = as.matrix(b$counts)); srt_bulk <- SetAssayData(srt_bulk, slot = "data", new.data = as.matrix(lg))
srt_bulk <- AddModuleScore(srt_bulk, features = top50, name = "Cm_", seed = 1)
bulk_cell_scores <- as.matrix(srt_bulk@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(bulk_cell_scores) <- CM
bulk_cell_scores <- bulk_cell_scores[rownames(meta), ]

bulk_median <- t(sapply(c("NTC", ORDER), function(g) apply(bulk_cell_scores[rownames(meta)[meta$Gene_Control == g], , drop = FALSE], 2, median)))
rownames(bulk_median) <- c("NTC", ORDER)
bulk_mean <- t(sapply(c("NTC", ORDER), function(g) colMeans(bulk_cell_scores[rownames(meta)[meta$Gene_Control == g], , drop = FALSE])))
rownames(bulk_mean) <- c("NTC", ORDER)

## validate bulk recompute (median, KO-NTC) vs legacy cm_sample_scores_AddModuleScore_bulk_RUN1.csv (mean, computed by spec44)
legacy_bulk <- read.csv(file.path(PP, "cm_sample_scores_AddModuleScore_bulk_RUN1.csv"), check.names = FALSE); rownames(legacy_bulk) <- legacy_bulk$sample
common_s <- intersect(rownames(meta), rownames(legacy_bulk))
cat(sprintf("\nbulk recompute vs legacy per-SAMPLE score correlation (sanity check, all Cm pooled, n=%d samples): r=%.3f\n",
            length(common_s), cor(as.vector(as.matrix(bulk_cell_scores[common_s, CM])), as.vector(as.matrix(legacy_bulk[common_s, CM])))))

## ---------------------------------------------------------------- (1) DIFF-based agreement (KO - NTC, median, fully self-contained, no legacy files)
screen_diff_konly <- SGN * sweep(screen_median[common_ko, CM, drop = FALSE], 2, screen_median["NTC", CM], "-")   # KO - NTC (no extra negation -- that was only needed for the NTC-KO legacy-matching validation above)
bulk_diff <- SGN * sweep(bulk_median[common_ko, CM, drop = FALSE], 2, bulk_median["NTC", CM], "-")               # KO - NTC

diff_res <- do.call(rbind, lapply(common_ko, function(k) {
  s <- screen_diff_konly[k, ]; bl <- bulk_diff[k, ]
  pr <- cor.test(s, bl); sr <- suppressWarnings(cor.test(s, bl, method = "spearman", exact = FALSE))
  data.frame(KO = k, pearson_r = pr$estimate, pearson_p = pr$p.value, spearman_rho = sr$estimate, spearman_p = sr$p.value)
}))
diff_res <- diff_res[order(-diff_res$pearson_r), ]
write.csv(diff_res, file.path(OUT, paste0("spec45_20Cm_diff_agreement", DSUF, ".csv")), row.names = FALSE)
options(width = 200)
cat("\n=== fully self-contained DIFF-based (KO-NTC, median) screen-vs-bulk agreement, 20 Cm modules ===\n")
print(format(diff_res, digits = 3), row.names = FALSE)
cat(sprintf("mean r=%.3f | significant (p<0.05): %d/%d\n", mean(diff_res$pearson_r), sum(diff_res$pearson_p < 0.05), nrow(diff_res)))

## ---------------------------------------------------------------- (2) RAW (no NTC diff) agreement -- now genuinely possible with true screen data
raw_res <- do.call(rbind, lapply(c("NTC", common_ko), function(k) {
  s <- screen_median[k, CM]; bl <- bulk_median[k, CM]
  pr <- cor.test(s, bl)
  data.frame(KO = k, pearson_r = pr$estimate, pearson_p = pr$p.value)
}))
raw_res <- raw_res[order(-raw_res$pearson_r), ]
write.csv(raw_res, file.path(OUT, "spec45_20Cm_raw_agreement.csv"), row.names = FALSE)
cat("\n=== fully self-contained RAW (no NTC diff) screen-vs-bulk agreement, 20 Cm modules -- NTC included as baseline ===\n")
print(format(raw_res, digits = 3), row.names = FALSE)
raw_ko_only <- raw_res[raw_res$KO != "NTC", ]
cat(sprintf("KOs only: mean r=%.3f | significant: %d/%d | NTC baseline r=%.3f (p=%.2e)\n",
            mean(raw_ko_only$pearson_r), sum(raw_ko_only$pearson_p < 0.05), nrow(raw_ko_only),
            raw_res$pearson_r[raw_res$KO == "NTC"], raw_res$pearson_p[raw_res$KO == "NTC"]))

## ---------------------------------------------------------------- (2b) RAW (no NTC diff), per-line breakdown -- bulk median per line, no NTC subtraction at all
RAW_LINES <- setdiff(LINES, "NCI-H2052")   # excluded throughout: only 1 NTC sample in the curated sheet, no replication for the NTC baseline row
bulk_median_line_raw <- lapply(RAW_LINES, function(cl) {
  mm <- meta[meta$Cell_Line == cl, ]
  t(sapply(c("NTC", common_ko), function(k) { s <- rownames(mm)[mm$Gene_Control == k]
    if (length(s) < 2) return(setNames(rep(NA, length(CM)), CM))
    apply(bulk_cell_scores[s, CM, drop = FALSE], 2, median) }))
})
names(bulk_median_line_raw) <- RAW_LINES; for (cl in RAW_LINES) rownames(bulk_median_line_raw[[cl]]) <- c("NTC", common_ko)

r_per_line_raw <- do.call(rbind, lapply(RAW_LINES, function(cl) do.call(rbind, lapply(c("NTC", common_ko), function(k) {
  bl <- bulk_median_line_raw[[cl]][k, ]; if (all(is.na(bl))) return(NULL)
  s <- screen_median[k, CM]
  data.frame(line = cl, KO = k, pearson_r = cor(s, bl, use = "complete.obs"))
}))))
write.csv(r_per_line_raw, file.path(OUT, "spec45_20Cm_raw_agreement_per_line.csv"), row.names = FALSE)
cat("\n=== per-line, per-KO Pearson r (screen vs that line's own bulk, RAW -- no NTC diff, 20 Cm modules) ===\n")
print(format(r_per_line_raw, digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- cross-KO "null" distribution: for each screen TF X, correlate its raw profile
## against every OTHER KO's raw bulk profile (Y != X, all lines) -- is the "own KO" correlation actually higher than mismatched pairs, or is the
## whole thing just generic architecture that would correlate with ANY bulk KO regardless of identity?
cross_pairs <- do.call(rbind, lapply(RAW_LINES, function(cl) do.call(rbind, lapply(common_ko, function(screen_k) {
  other_ko <- setdiff(common_ko, screen_k)
  do.call(rbind, lapply(other_ko, function(bulk_k) {
    bl <- bulk_median_line_raw[[cl]][bulk_k, ]; if (all(is.na(bl))) return(NULL)
    data.frame(line = cl, screen_KO = screen_k, bulk_KO = bulk_k, pearson_r = cor(screen_median[screen_k, CM], bl, use = "complete.obs"))
  }))
}))))
write.csv(cross_pairs, file.path(OUT, "spec45_20Cm_raw_cross_KO_pairs.csv"), row.names = FALSE)
cat(sprintf("\ncross-KO (mismatched) pairs: n=%d, median r=%.3f (own-KO median for comparison: %.3f)\n",
            nrow(cross_pairs), median(cross_pairs$pearson_r), median(r_per_line_raw$pearson_r[r_per_line_raw$KO != "NTC"])))

save_fig("fig3D_20Cm_raw_correlation_per_TF_boxplot", 9, 5.5, function() {
  ko_only_raw <- r_per_line_raw[r_per_line_raw$KO != "NTC", ]
  ko_ord <- names(sort(sapply(split(ko_only_raw$pearson_r, ko_only_raw$KO), median), decreasing = TRUE))
  ko_only_raw$KO <- factor(ko_only_raw$KO, levels = ko_ord)
  n_ko <- length(ko_ord)
  par(mar = c(5, 4.8, 3.6, 1), family = "sans")
  plot(NA, xlim = c(0.5, n_ko + 0.5), ylim = c(-1, 1), xaxt = "n", xlab = "", ylab = "Pearson r (raw screen vs raw bulk, 20 Cm modules)")
  axis(1, at = seq_len(n_ko), labels = FALSE)
  text(seq_len(n_ko), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), ko_ord, srt = 45, adj = 1, xpd = NA)
  # grey "other KO pairs" box, screen_KO == this TF, against every OTHER bulk KO (the null: would ANY bulk KO correlate this well?)
  for (i in seq_len(n_ko)) {
    d_null <- cross_pairs$pearson_r[cross_pairs$screen_KO == ko_ord[i]]
    boxplot(d_null, at = i - 0.18, add = TRUE, boxwex = 0.3, col = "grey85", border = "grey40", outline = FALSE, axes = FALSE)
  }
  # coloured "own KO" box
  for (i in seq_len(n_ko)) {
    d_own <- ko_only_raw$pearson_r[ko_only_raw$KO == ko_ord[i]]
    boxplot(d_own, at = i + 0.18, add = TRUE, boxwex = 0.3, col = adjustcolor(KO_COL[ko_ord[i]], 0.5), border = INK2, outline = FALSE, axes = FALSE)
  }
  set.seed(1); for (i in seq_len(n_ko)) { d <- ko_only_raw[ko_only_raw$KO == ko_ord[i], ]
    points(jitter(rep(i + 0.18, nrow(d)), amount = 0.06), d$pearson_r, pch = 16, cex = 1, col = LINE_COL[d$line]) }
  abline(h = 0, col = INK, lty = 2)
  ntc_r <- r_per_line_raw$pearson_r[r_per_line_raw$KO == "NTC"]
  abline(h = median(ntc_r, na.rm = TRUE), col = "grey40", lty = 3)
  present_lines <- sort(unique(ko_only_raw$line))
  legend("bottomleft", legend = c(present_lines, "own KO", "other KOs (null)", "NTC median"),
         col = c(LINE_COL[present_lines], "grey30", "grey60", "grey40"),
         pch = c(rep(16, length(present_lines)), 15, 15, NA), lty = c(rep(NA, length(present_lines) + 2), 3), bty = "n", cex = 0.6)
  mtext("Raw screen-to-bulk correlation (20 Cm modules, no NTC differencing), per TF", side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext("coloured = screen TF vs its OWN bulk KO; grey = screen TF vs every OTHER bulk KO (null); dotted = NTC baseline; NCI-H2052 excluded", side = 3, adj = 0, line = 0.9, cex = 0.6, col = INK2)
})

## ---------------------------------------------------------------- (3) per-line breakdown: same diff-based (median, KO-NTC) correlation, computed within each line's own bulk samples
bulk_median_line <- lapply(LINES, function(cl) {
  mm <- meta[meta$Cell_Line == cl, ]; ntc_l <- rownames(mm)[mm$Gene_Control == "NTC"]
  t(sapply(common_ko, function(k) { ko_l <- rownames(mm)[mm$Gene_Control == k]
    if (!length(ko_l) || length(ntc_l) < 2) return(setNames(rep(NA, length(CM)), CM))   # NCI-H2052 excluded: only 1 NTC sample in the curated sheet
    SGN * (apply(bulk_cell_scores[ko_l, CM, drop = FALSE], 2, median) - apply(bulk_cell_scores[ntc_l, CM, drop = FALSE], 2, median)) }))
})
names(bulk_median_line) <- LINES; for (cl in LINES) rownames(bulk_median_line[[cl]]) <- common_ko

r_per_line <- do.call(rbind, lapply(LINES, function(cl) do.call(rbind, lapply(common_ko, function(k) {
  bl <- bulk_median_line[[cl]][k, ]; if (all(is.na(bl))) return(NULL)
  s <- screen_diff_konly[k, ]
  data.frame(line = cl, KO = k, pearson_r = cor(s, bl, use = "complete.obs"),
             spearman_rho = suppressWarnings(cor(s, bl, method = "spearman", use = "complete.obs")))
}))))
write.csv(r_per_line, file.path(OUT, paste0("spec45_20Cm_agreement_per_line", DSUF, ".csv")), row.names = FALSE)
cat("\n=== per-line, per-KO Pearson r (screen vs that line's own bulk, median-based KO-NTC diff, 20 Cm modules) ===\n")
print(format(r_per_line, digits = 3), row.names = FALSE)
cat("\nmedian Pearson r per line:\n"); print(round(sapply(split(r_per_line$pearson_r, r_per_line$line), median), 3))

## ---------------------------------------------------------------- cross-KO "null" distribution for the DIFF-based (KO-NTC) comparison: for each
## screen TF X, correlate its KO-NTC diff profile against every OTHER KO's bulk KO-NTC diff profile (Y != X, all lines)
cross_pairs_diff <- do.call(rbind, lapply(LINES, function(cl) do.call(rbind, lapply(common_ko, function(screen_k) {
  other_ko <- setdiff(common_ko, screen_k)
  do.call(rbind, lapply(other_ko, function(bulk_k) {
    bl <- bulk_median_line[[cl]][bulk_k, ]; if (all(is.na(bl))) return(NULL)
    data.frame(line = cl, screen_KO = screen_k, bulk_KO = bulk_k, pearson_r = cor(screen_diff_konly[screen_k, ], bl, use = "complete.obs"))
  }))
}))))
write.csv(cross_pairs_diff, file.path(OUT, paste0("spec45_20Cm_diff_cross_KO_pairs", DSUF, ".csv")), row.names = FALSE)
cat(sprintf("\ncross-KO (mismatched) pairs, diff-based: n=%d, median r=%.3f (own-KO median for comparison: %.3f)\n",
            nrow(cross_pairs_diff), median(cross_pairs_diff$pearson_r), median(r_per_line$pearson_r)))

save_fig(paste0("fig3D_20Cm_correlation_per_TF_boxplot", DSUF), 9, 5.5, function() {
  ko_ord <- names(sort(sapply(split(r_per_line$pearson_r, r_per_line$KO), median), decreasing = TRUE))
  r_per_line$KO <- factor(r_per_line$KO, levels = ko_ord)
  n_ko <- length(ko_ord)
  par(mar = c(5, 4.8, 3.6, 1), family = "sans")
  plot(NA, xlim = c(0.5, n_ko + 0.5), ylim = c(-1, 1), xaxt = "n", xlab = "", ylab = paste0("Pearson r (screen vs bulk, median ", DLAB, " diff, 20 Cm modules)"))
  axis(1, at = seq_len(n_ko), labels = FALSE)
  text(seq_len(n_ko), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), ko_ord, srt = 45, adj = 1, xpd = NA)
  for (i in seq_len(n_ko)) {
    d_null <- cross_pairs_diff$pearson_r[cross_pairs_diff$screen_KO == ko_ord[i]]
    boxplot(d_null, at = i - 0.18, add = TRUE, boxwex = 0.3, col = "grey85", border = "grey40", outline = FALSE, axes = FALSE)
  }
  for (i in seq_len(n_ko)) {
    d_own <- r_per_line$pearson_r[r_per_line$KO == ko_ord[i]]
    boxplot(d_own, at = i + 0.18, add = TRUE, boxwex = 0.3, col = adjustcolor(KO_COL[ko_ord[i]], 0.5), border = INK2, outline = FALSE, axes = FALSE)
  }
  set.seed(1); for (i in seq_len(n_ko)) { d <- r_per_line[r_per_line$KO == ko_ord[i], ]
    points(jitter(rep(i + 0.18, nrow(d)), amount = 0.06), d$pearson_r, pch = 16, cex = 1, col = LINE_COL[d$line]) }
  abline(h = 0, col = INK, lty = 2)
  present_lines <- sort(unique(r_per_line$line))
  legend("bottomleft", legend = c(present_lines, "own KO", "other KOs (null)"),
         col = c(LINE_COL[present_lines], "grey30", "grey60"), pch = c(rep(16, length(present_lines)), 15, 15), bty = "n", cex = 0.6)
  mtext(paste0("Screen-to-bulk correlation (20 Cm modules, ", DLAB, " diff), per TF"), side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext("coloured = screen TF vs its OWN bulk KO diff; grey = screen TF vs every OTHER bulk KO diff (null); NCI-H2052 excluded (only 1 NTC sample)", side = 3, adj = 0, line = 0.9, cex = 0.6, col = INK2)
})

## ---------------------------------------------------------------- Spearman version of the per-TF boxplot (own-KO per line vs cross-KO null), same data as above
EXCL <- setdiff(strsplit(Sys.getenv("EXCL_LINES", ""), ",")[[1]], ""); ESUF <- if (length(EXCL)) paste0("_no", paste(sub("NCI-|-211H", "", EXCL), collapse = "_")) else ""   # optional: drop lines (env EXCL_LINES=NCI-H28)
r_per_line <- r_per_line[!r_per_line$line %in% EXCL, ]
cross_pairs_sp <- do.call(rbind, lapply(setdiff(LINES, EXCL), function(cl) do.call(rbind, lapply(common_ko, function(screen_k) {
  do.call(rbind, lapply(setdiff(common_ko, screen_k), function(bulk_k) {
    bl <- bulk_median_line[[cl]][bulk_k, ]; if (all(is.na(bl))) return(NULL)
    data.frame(line = cl, screen_KO = screen_k, bulk_KO = bulk_k, spearman_rho = suppressWarnings(cor(screen_diff_konly[screen_k, ], bl, method = "spearman", use = "complete.obs")))
  }))
}))))
write.csv(cross_pairs_sp, file.path(OUT, paste0("spec45_20Cm_diff_cross_KO_pairs_spearman", DSUF, ESUF, ".csv")), row.names = FALSE)
cat(sprintf("\nSpearman: own-KO median rho=%.3f | cross-KO (null) median rho=%.3f | separation=%+.3f\n", median(r_per_line$spearman_rho), median(cross_pairs_sp$spearman_rho), median(r_per_line$spearman_rho) - median(cross_pairs_sp$spearman_rho)))
sp_tab <- do.call(rbind, lapply(split(r_per_line, r_per_line$KO), function(d) data.frame(KO = d$KO[1], own_median_rho = median(d$spearman_rho), null_median_rho = median(cross_pairs_sp$spearman_rho[cross_pairs_sp$screen_KO == d$KO[1]]))))
sp_tab$separation <- sp_tab$own_median_rho - sp_tab$null_median_rho; sp_tab <- sp_tab[order(-sp_tab$own_median_rho), ]; print(format(sp_tab, digits = 3), row.names = FALSE)
write.csv(sp_tab, file.path(OUT, paste0("spec45_20Cm_spearman_per_TF_summary", DSUF, ESUF, ".csv")), row.names = FALSE)
save_fig(paste0("fig3D_20Cm_correlation_per_TF_boxplot", DSUF, "_spearman", ESUF), 9, 5.5, function() {
  ko_ord <- names(sort(sapply(split(r_per_line$spearman_rho, r_per_line$KO), median), decreasing = TRUE)); n_ko <- length(ko_ord)
  par(mar = c(5, 4.8, 3.6, 1), family = "sans")
  plot(NA, xlim = c(0.5, n_ko + 0.5), ylim = c(-1, 1), xaxt = "n", xlab = "", ylab = paste0("Spearman rho (screen vs bulk, ", DLAB, " diff)"))
  axis(1, at = seq_len(n_ko), labels = FALSE); text(seq_len(n_ko), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), ko_ord, srt = 45, adj = 1, xpd = NA)
  for (i in seq_len(n_ko)) boxplot(cross_pairs_sp$spearman_rho[cross_pairs_sp$screen_KO == ko_ord[i]], at = i - 0.18, add = TRUE, boxwex = 0.3, col = "grey85", border = "grey40", outline = FALSE, axes = FALSE)
  for (i in seq_len(n_ko)) boxplot(r_per_line$spearman_rho[r_per_line$KO == ko_ord[i]], at = i + 0.18, add = TRUE, boxwex = 0.3, col = adjustcolor(KO_COL[ko_ord[i]], 0.5), border = INK2, outline = FALSE, axes = FALSE)
  set.seed(1); for (i in seq_len(n_ko)) { d <- r_per_line[r_per_line$KO == ko_ord[i], ]; points(jitter(rep(i + 0.18, nrow(d)), amount = 0.06), d$spearman_rho, pch = 16, cex = 1, col = LINE_COL[d$line]) }
  abline(h = 0, col = INK, lty = 2); present_lines <- sort(unique(r_per_line$line))
  legend("bottomleft", legend = c(present_lines, "own KO", "other KOs (null)"), col = c(LINE_COL[present_lines], "grey30", "grey60"), pch = c(rep(16, length(present_lines)), 15, 15), bty = "n", cex = 0.6)
  mtext(paste0("Screen-to-bulk Spearman correlation (20 Cm modules, ", DLAB, " diff), per TF"), side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext(paste0("coloured = screen TF vs its OWN bulk KO diff; grey = screen TF vs every OTHER bulk KO diff (null); NCI-H2052 excluded (only 1 NTC sample)", if (length(EXCL)) paste0("; ", paste(EXCL, collapse = ", "), " removed") else ""), side = 3, adj = 0, line = 0.9, cex = 0.6, col = INK2)
})

## ---------------------------------------------------------------- fig: validation scatter (one panel per KO: legacy vs recomputed screen diff)
save_fig("fig3D_20Cm_recompute_validation", 10, 5.5, function() {
  layout(matrix(1:8, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3, 0), family = "sans")
  for (k in common_ko) {
    a <- as.numeric(legacy_screen_diff[k, CM]); bb <- screen_diff_recompute[k, CM]
    par(mar = c(4, 4, 2.5, 1))
    plot(a, bb, pch = 16, col = adjustcolor(KO_COL[k], 0.7), xlab = "legacy NTC-KO diff", ylab = "recomputed NTC-KO diff", cex.lab = 0.8)
    abline(0, 1, col = INK, lty = 2); mtext(sprintf("%s (r=%.2f)", k, cor(a, bb)), side = 3, font = 2, cex = 0.85, col = KO_COL[k])
  }
  mtext("Validation: our from-scratch AddModuleScore recompute vs the legacy screen_cm_diff_matrix.csv", outer = TRUE, side = 3, adj = 0, line = 1.2, font = 2, cex = 1.0)
})

## ---------------------------------------------------------------- fig: diff-based vs raw agreement, side by side
save_fig("fig3D_20Cm_diff_vs_raw_correlation", 9, 5.5, function() {
  layout(matrix(1:2, nrow = 1)); par(oma = c(0, 0, 3, 0), family = "sans")
  par(mar = c(6, 4.8, 2.5, 1))
  ord1 <- diff_res[order(-diff_res$pearson_r), ]
  barplot(ord1$pearson_r, names.arg = ord1$KO, col = ifelse(ord1$pearson_p < 0.05, KO_COL[ord1$KO], adjustcolor(KO_COL[ord1$KO], 0.4)),
          border = INK2, las = 2, ylab = "Pearson r", ylim = c(-1, 1), main = "KO - NTC diff (self-contained)")
  abline(h = 0, col = INK)
  par(mar = c(6, 4.8, 2.5, 1))
  ord2 <- raw_res[order(-raw_res$pearson_r), ]
  cols2 <- ifelse(ord2$KO == "NTC", "grey40", ifelse(ord2$pearson_p < 0.05, KO_COL[ord2$KO], adjustcolor(KO_COL[ord2$KO], 0.4)))
  barplot(ord2$pearson_r, names.arg = ord2$KO, col = cols2, border = INK2, las = 2, ylab = "Pearson r", ylim = c(-1, 1), main = "raw, no NTC diff")
  abline(h = 0, col = INK); abline(h = raw_res$pearson_r[raw_res$KO == "NTC"], col = "grey40", lty = 2)
  mtext("20 Cm modules, fully self-contained recompute: diff-based vs raw agreement", outer = TRUE, side = 3, adj = 0, line = 1.2, font = 2, cex = 1.0)
})
