###############################################################################
# STEP 11 -- TF motifs that differ between subclones, by a per-tumour test on the cells.
#
# For each tumour with two clones and each motif:
#   Wilcoxon rank-sum test of the per-cell chromVAR z between the cells of clone c1 and clone c2
#   effect = |AUC - 0.5|, AUC = probability that a random c2 cell has higher activity than a
#            random c1 cell (0 = identical distributions, 0.5 = complete separation). Unlike the
#            P value it does not grow with the number of cells, so tumours of very different
#            size (P10 889 cells, P23 11,331) can be combined. Unsigned: c1 / c2 carry no order.
#   P, and Benjamini-Hochberg FDR over the motifs of that tumour
# Across tumours: motifs ranked by the number of tumours with FDR < 0.05, then by the median -log10 P.
#
# chromVAR z  per-cell deviation z-score for each cisBP motif, from the MotifMatrix of the tumour-compartment
#             ArchR project; pulling it out is the slow step, so it is cached to chromvar_z_cache.rds.
# Input : epi_clone_labels.csv (step 3); ArchR project (first run only); tumour-compartment scRNA (first run only)
# Output: chromvar_z_cache.rds, epi_chromvar_wilcoxon_per_tumour.csv, epi_chromvar_wilcoxon_summary_pvalue.csv,
#         Plots/epi_TF_variability_wilcoxon.pdf/.png
###############################################################################
suppressMessages({ library(data.table); library(presto); library(ggplot2) })
source("00_common.R")
N_TOP <- as.integer(Sys.getenv("N_TOP", 25))
ZCACHE <- "chromvar_z_cache.rds"
## Row names carry a numeric suffix (TFAP2B_1); strip it to the TF symbol.
if (file.exists(ZCACHE)) { Z <- readRDS(ZCACHE) } else {
  suppressMessages(library(ArchR)); addArchRThreads(1)
  proj <- loadArchRProject(ARCHR, showLogo = FALSE)
  se <- getMatrixFromProject(proj, useMatrix = "MotifMatrix", verbose = FALSE)
  Z <- as.matrix(SummarizedExperiment::assays(se)$z)
  rownames(Z) <- make.unique(sub("_[0-9]+$", "", SummarizedExperiment::rowData(se)$name))
  saveRDS(Z, ZCACHE)
}
LAB <- fread("epi_clone_labels.csv")

res <- list()
for (S in SAMPLES){
  l <- LAB[sample == S]; v <- setNames(l$clone, l$cell)
  cells <- intersect(names(v), colnames(Z)); v <- v[cells]
  if (length(unique(v)) < 2 || min(table(v)) < 20) next
  w <- wilcoxauc(Z[, cells, drop = FALSE], v); w <- w[w$group == "c2", ]        # c2 vs c1; two-sided P
  res[[S]] <- data.frame(sample = S, TF = w$feature, n_c1 = sum(v == "c1"), n_c2 = sum(v == "c2"),
                         auc = w$auc, effect = abs(w$auc - 0.5), p = w$pval, fdr = p.adjust(w$pval, "BH"))
  cat(sprintf("%-4s c1 %5d  c2 %5d | motifs FDR < 0.05: %3d | max |AUC-0.5| %.3f\n", S, sum(v == "c1"), sum(v == "c2"),
              sum(res[[S]]$fdr < 0.05), max(res[[S]]$effect)))
}
R <- do.call(rbind, res)

## Display: motifs ranked by REPRODUCIBILITY = number of tumours in which the motif differs between the
## clones (FDR < 0.05), ties broken by the median -log10 P. P = 0 (underflow) is set to the smallest non-zero P.
##
## Expression filter, as in tumor_analysis/discover_oncogenic_TFs.R (min_exp 0.5): in a tumour with matched scRNA-seq a
## motif is kept only if its TF is expressed, log2(mean normalised expression of the tumour's malignant cells + 1) >= 0.5.
## The clones are defined on scATAC cells and cannot be told apart in scRNA, so the filter is on the whole tumour
## ("expressed in at least one clone" cannot be lower than that). Tumours without malignant scRNA cells are not filtered.
R <- as.data.table(R); pmin_nz <- min(R$p[R$p > 0]); R[, mlp := -log10(pmax(p, pmin_nz))]
MIN_EXP <- 0.5; EXPR_CACHE <- "tf_rna_expression_per_tumour.csv"
if (!file.exists(EXPR_CACHE)) {
  suppressMessages(library(Seurat))
  srt <- readRDS(file.path(dirname(dirname(getwd())), "tumor_compartment", "scrna", "srt.rds")); srt <- srt[, srt$celltype == "Malignant"]
  g <- intersect(unique(R$TF), rownames(srt))
  ps <- log2(as.matrix(AverageExpression(srt, features = g, group.by = "sampleID")[[1]]) + 1)
  E <- data.frame(TF = rownames(ps), ps, check.names = FALSE); fwrite(E, EXPR_CACHE)
  fwrite(data.frame(sample = names(table(srt$sampleID)), n_malignant_rna = as.integer(table(srt$sampleID))), "tf_rna_expression_ncells.csv")
}
E <- fread(EXPR_CACHE); NC <- fread("tf_rna_expression_ncells.csv")
rna_ok <- intersect(unique(R$sample), intersect(colnames(E), NC$sample[NC$n_malignant_rna >= 50]))
cat("tumours with matched malignant scRNA (filter applied):", rna_ok, "| not filtered:", setdiff(unique(R$sample), rna_ok), "\n")
R[, expr := NA_real_]
for (S in rna_ok) R[sample == S, expr := E[[S]][match(TF, E$TF)]]
R[, keep := !(sample %in% rna_ok) | (!is.na(expr) & expr >= MIN_EXP)]
cat("motif x tumour tests kept:", sum(R$keep), "of", nrow(R), "\n"); print(R[, .(kept = sum(keep), dropped = sum(!keep)), by = sample])
fwrite(R, "epi_chromvar_wilcoxon_per_tumour.csv")
K <- R[keep == TRUE]
PS <- K[, .(n_tested = .N, n_sig = sum(fdr < 0.05), median_mlp = median(mlp), q25_mlp = quantile(mlp, .25), q75_mlp = quantile(mlp, .75),
            median_effect = median(effect), tumours_sig = paste(sample[fdr < 0.05], collapse = ";")), by = TF][order(-n_sig, -median_mlp)]
fwrite(PS, "epi_chromvar_wilcoxon_summary_pvalue.csv")
cat("\nmotifs by number of tumours significant (after the expression filter):\n"); print(table(n_sig = PS$n_sig, n_tested = PS$n_tested))
print(head(PS, N_TOP), digits = 3)
top <- head(PS, N_TOP); top$TF <- factor(top$TF, rev(top$TF)); nS <- length(unique(R$sample))
pts <- K[TF %in% top$TF]; pts$TF <- factor(pts$TF, levels(top$TF)); pts$sig <- ifelse(pts$fdr < 0.05, "FDR < 0.05", "not significant")
top$nlab <- factor(top$n_sig, levels = 0:nS)
## transposed layout: TFs along the x axis (most reproducible first), short panel
lv <- as.character(top$TF)[order(-top$n_sig, -top$median_mlp)]; pts$TF <- factor(pts$TF, lv); top$TF <- factor(top$TF, lv)
pA <- ggplot(pts, aes(TF, mlp)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3, fill = "grey90", colour = "grey30") +   # box over the tumours
  geom_point(aes(shape = sig), colour = "#c0392b", size = 0.9, position = position_jitter(width = 0.12, height = 0, seed = 1)) +
  scale_shape_manual(values = c(`FDR < 0.05` = 16, `not significant` = 1), name = NULL) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 2, 5, 10, 20, 50)) +
  labs(x = NULL, y = "-log10 P\n(c1 vs c2)", title = "TF motifs differing between subclones") +
  theme_classic(base_size = 7.5) + theme(plot.title = element_text(face = "bold", size = 8), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                                         legend.position = "top", legend.key.size = unit(2.5, "mm"), legend.margin = margin(0, 0, 0, 0), plot.margin = margin(2, 3, 0, 2))
pN <- ggplot(top, aes(TF, n_sig)) + geom_col(fill = "grey35", width = 0.7) +
  scale_y_continuous(limits = c(0, nS), breaks = c(0, nS), expand = c(0, 0)) + labs(x = NULL, y = "tumours\nFDR<0.05") +
  theme_classic(base_size = 7.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5), plot.margin = margin(1, 3, 2, 2))
p <- patchwork::wrap_plots(pA, pN, ncol = 1, heights = c(3, 1))
ggsave("Plots/epi_TF_variability_wilcoxon.pdf", p, width = 4.6, height = 2.3)
ggsave("Plots/epi_TF_variability_wilcoxon.png", p, width = 4.6, height = 2.3, dpi = 220)

cat("DONE\n")
