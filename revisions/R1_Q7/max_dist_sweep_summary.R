# max_dist_sweep_summary.R
#
# Summary table of the maximum distance benchmark: hubs at 5 / 10 / 12.5 / 15 / 20 kb (max_dist_sweep_SE.csv + myeloid_SE_hub_enrichment.csv
# for 12.5 kb), per method x setting (TSS filter x minimum peaks) x distance:
#   score_on = SE recovery, mean -log10 FDR of the matched cell type | SE pairs (capped at the 99th percentile of the 12.5 kb table)
#   spec     = SE specificity, score_on minus the mean of the same cell types' non-matching benchmark SE sets
# Output: plot_max_dist_sweep_table.csv (input of plot_benchmark_summary.R, panel B)

script_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q7"
rd <- function(f) { p <- file.path(script_dir, f); if (!file.exists(p)) return(NULL); x <- read.csv(p, stringsAsFactors = FALSE); x[!duplicated(x[, c("condition", "celltype", "se_db")]), c("condition", "celltype", "se_db", "neg_log10_fdr")] }
base <- rd("myeloid_SE_hub_enrichment.csv"); cap <- unname(quantile(base$neg_log10_fdr, 0.99, na.rm = TRUE))
base$condition <- ifelse(grepl("_tss\\d+$", base$condition), base$condition, paste0(base$condition, "_tss0")); base$condition <- paste0(base$condition, "_d12.5")
res <- rbind(base, rd("max_dist_sweep_SE.csv")); res$score <- pmin(res$neg_log10_fdr, cap)
res$dist <- as.numeric(sub(".*_d([0-9.]+)$", "\\1", res$condition)); res$tss <- as.integer(sub(".*_tss(\\d+)_d.*", "\\1", res$condition)); res$mp <- as.integer(sub(".*_mp(\\d+).*", "\\1", res$condition)); res$method <- sub("_mp\\d+.*", "", res$condition)
methods <- c(CoAccess_r0.2 = "r 0.2", CoAccess_r0.3 = "r 0.3", CoAccess_r0.4 = "r 0.4", CoAccess_r0.5 = "r 0.5", CoAccess_r0.6 = "r 0.6", Stitching = "Stitch-AP", Stitching_CT = "Stitch CT"); res <- res[res$method %in% names(methods), ]
pairs <- data.frame(celltype = c("Alveolar", "Alveolar", "Fibroblasts", "SmoothMuscle", "Endothelial", "Myeloid", "T_cells", "B_cells"),
  se_db = c("AT1", "AT2", "ENCODE_Lung_primary_cell_fibroblast-of-lung", "Smoothmuscle", "Endothelial", "ENCODE_Blood_primary_cell_CD14-positive-monocyte", "Tcell", "ENCODE_Blood_primary_cell_Bcell"))
se_bench <- c(unique(pairs$se_db), "ENCODE_Mammary_gland_primary_cell_mammary-epithelial-cell", "ENCODE_Skin_primary_cell_keratinocyte")
sub <- res[res$celltype %in% pairs$celltype & res$se_db %in% se_bench, ]; sub$matched <- paste(sub$celltype, sub$se_db) %in% paste(pairs$celltype, pairs$se_db); sub$is_off <- !sub$matched & !(sub$celltype == "Alveolar" & sub$se_db %in% c("AT1", "AT2"))
grid <- unique(sub[, c("condition", "method", "dist", "mp", "tss")]); full <- merge(merge(grid, pairs), sub[sub$matched, c("condition", "celltype", "se_db", "score")], all.x = TRUE); full$score[is.na(full$score)] <- 0   # matched pair without DA hubs = 0
on <- aggregate(score ~ condition + method + dist + mp + tss, full, mean); off <- aggregate(score ~ condition, sub[sub$is_off, ], mean)
S <- merge(on, off, by = "condition", suffixes = c("_on", "_off"), all.x = TRUE); S$score_off[is.na(S$score_off)] <- 0; S$spec <- S$score_on - S$score_off
S$method <- factor(methods[S$method], rev(methods)); S$distf <- factor(paste(S$dist, "kb"), paste(sort(unique(S$dist)), "kb"))
S$setting <- factor(sprintf("%s\n>= %d peaks", ifelse(S$tss == 0, "no TSS filter", "TSS >= 2 kb"), S$mp), levels = as.vector(outer(3:5, c("no TSS filter", "TSS >= 2 kb"), function(m, t) sprintf("%s\n>= %d peaks", t, m))))
write.csv(S[order(S$tss, S$mp, S$method, S$dist), ], file.path(script_dir, "plot_max_dist_sweep_table.csv"), row.names = FALSE)
print(xtabs(round(score_on, 1) ~ method + distf, S[S$mp == 5 & S$tss == 0, ])); print(xtabs(round(spec, 1) ~ method + distf, S[S$mp == 5 & S$tss == 0, ]))
