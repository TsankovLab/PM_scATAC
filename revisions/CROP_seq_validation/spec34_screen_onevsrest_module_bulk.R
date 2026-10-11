# Spec 34: for each KO, identify its top DE genes in the CROP-seq screen TWO WAYS -- (A) KO vs NTC, (B) KO vs the REST of the guides
# (all other non-NTC KO groups pooled) -- then take the top 50 genes MOST UPREGULATED in the KO from EACH comparison (2 gene lists
# per KO) and check whether those same genes are ALSO more highly expressed in that KO's own BULK RNA-seq samples: (1) an aggregate
# module score, one-sided tested (bulk KO > bulk NTC) per line, and (2) a simple per-gene concordance rate (fraction of the 50 genes
# individually higher in bulk KO than bulk NTC, pooled across lines).
# needs: spec1_prepare_stats.R (spec_screen_stats.rds, for the vs-NTC comparison), spec_common.R (load_screen(), for vs-rest)
# outputs: spec34_screen_top50_upregulated_genes.csv, spec34_bulk_module_scores.csv, spec34_KO_module_bulk_test.csv,
#          spec34_gene_level_concordance.csv, fig2Q_screen_top50_module_bulk_heatmap.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Matrix))

fab <- load_full_audited_bulk(); bulk_genes <- rownames(fab$lg); meta <- fab$meta; lg <- fab$lg

## ---------------------------------------------------------------- (A) vs NTC: reuse the already-computed screen stats (spec1_prepare_stats.R)
ss <- readRDS(file.path(OUT, "spec_screen_stats.rds"))   # t0 = Welch t vs NTC, no cc regression (matches the rest of this pipeline's default)

## ---------------------------------------------------------------- (B) vs rest of the guides: KO's cells vs all OTHER non-NTC KO groups pooled
srt <- load_screen(); Y <- GetAssayData(srt, layer = "data"); grp <- as.character(srt$merged_call)
Y <- as(Y, "CsparseMatrix"); Y2 <- Y; Y2@x <- Y2@x^2
grp_stats <- function(idx) { n <- length(idx); my <- as.numeric(rowSums(Y[, idx, drop = FALSE]) / n)
  vy <- as.numeric((rowSums(Y2[, idx, drop = FALSE]) - n * my^2) / (n - 1)); list(n = n, mean = my, var = pmax(vy, 1e-6)) }
t_vs_rest <- list()
for (k in ORDER) {
  ko_idx <- which(grp == k); rest_idx <- which(grp != "NTC" & grp != k)
  gko <- grp_stats(ko_idx); grest <- grp_stats(rest_idx)
  t_vs_rest[[k]] <- setNames((gko$mean - grest$mean) / sqrt(gko$var / gko$n + grest$var / grest$n), rownames(Y))
}

## ---------------------------------------------------------------- top 50 UPREGULATED genes per KO, per comparison (restricted to genes also measured in bulk)
top50 <- list()
for (k in ORDER) {
  if (k %in% colnames(ss$t0)) { x <- ss$t0[, k]; x <- x[names(x) %in% bulk_genes & is.finite(x)]
    top50[[paste(k, "vsNTC")]] <- data.frame(KO = k, comparison = "vsNTC", gene = names(sort(x, decreasing = TRUE))[1:50], t = sort(x, decreasing = TRUE)[1:50]) }
  x <- t_vs_rest[[k]]; x <- x[names(x) %in% bulk_genes & is.finite(x)]
  top50[[paste(k, "vsRest")]] <- data.frame(KO = k, comparison = "vsRest", gene = names(sort(x, decreasing = TRUE))[1:50], t = sort(x, decreasing = TRUE)[1:50])
}
top50 <- do.call(rbind, top50); rownames(top50) <- NULL
write.csv(top50, file.path(OUT, "spec34_screen_top50_upregulated_genes.csv"), row.names = FALSE)
cat(sprintf("top50-upregulated gene lists: %d (KO x comparison)\n", length(unique(paste(top50$KO, top50$comparison)))))

## ---------------------------------------------------------------- score each (KO, comparison) module in bulk (all 76 samples), z-scored
Z <- (lg - rowMeans(lg)) / matrixStats::rowSds(as.matrix(lg)); Z[!is.finite(Z)] <- 0
mod_keys <- unique(paste(top50$KO, top50$comparison))
score <- sapply(mod_keys, function(mk) { g <- top50$gene[paste(top50$KO, top50$comparison) == mk]; colMeans(Z[g, , drop = FALSE]) })
rownames(score) <- rownames(meta)
write.csv(cbind(sample = rownames(score), meta[, c("Sample_ID", "Cell_Line", "Gene_Control")], round(score, 4)), file.path(OUT, "spec34_bulk_module_scores.csv"), row.names = FALSE)

## ---------------------------------------------------------------- (1) aggregate module test: bulk KO > bulk NTC, one-sided, per line + pooled
test_res <- list()
for (k in ORDER) for (comp in c("vsNTC", "vsRest")) {
  mk <- paste(k, comp); if (!(mk %in% colnames(score))) next
  pooled_ko <- c(); pooled_ntc <- c()
  for (cl in LINES) {
    m <- meta[meta$Cell_Line == cl, ]; s <- rownames(m)
    ko_v <- score[s[m$Gene_Control == k], mk]; ntc_v <- score[s[m$Gene_Control == "NTC"], mk]
    pooled_ko <- c(pooled_ko, ko_v); pooled_ntc <- c(pooled_ntc, ntc_v)
    if (length(ko_v) < 1 || length(ntc_v) < 2) next
    p <- if (length(ko_v) >= 2) wilcox.test(ko_v, ntc_v, alternative = "greater", exact = FALSE)$p.value else NA
    test_res[[length(test_res) + 1]] <- data.frame(KO = k, comparison = comp, line = cl, n_ko = length(ko_v), n_ntc = length(ntc_v),
                                                    mean_ko = mean(ko_v), mean_ntc = mean(ntc_v), p_greater = p)
  }
  p_pooled <- wilcox.test(pooled_ko, pooled_ntc, alternative = "greater", exact = FALSE)$p.value
  test_res[[length(test_res) + 1]] <- data.frame(KO = k, comparison = comp, line = "POOLED", n_ko = length(pooled_ko), n_ntc = length(pooled_ntc),
                                                  mean_ko = mean(pooled_ko), mean_ntc = mean(pooled_ntc), p_greater = p_pooled)
}
test_res <- do.call(rbind, test_res)
write.csv(test_res, file.path(OUT, "spec34_KO_module_bulk_test.csv"), row.names = FALSE)
options(width = 200); cat("\n=== are the screen's top50-upregulated genes ALSO more expressed in that KO's own bulk samples? (one-sided: bulk KO > bulk NTC) ===\n")
print(format(test_res[order(test_res$comparison, test_res$KO, test_res$line != "POOLED"), ], digits = 3), row.names = FALSE)
pooled_only <- test_res[test_res$line == "POOLED", ]
cat(sprintf("\nsignificant (p<0.05), POOLED across lines: vsNTC = %d/8, vsRest = %d/8\n",
            sum(pooled_only$p_greater[pooled_only$comparison == "vsNTC"] < 0.05), sum(pooled_only$p_greater[pooled_only$comparison == "vsRest"] < 0.05)))

## ---------------------------------------------------------------- (2) simple per-gene concordance: of the 50 genes, what fraction are individually higher in bulk KO vs bulk NTC (pooled across lines)?
conc <- list()
for (k in ORDER) for (comp in c("vsNTC", "vsRest")) {
  genes_k <- top50$gene[top50$KO == k & top50$comparison == comp]; if (!length(genes_k)) next
  ko_s <- rownames(meta)[meta$Gene_Control == k]; ntc_s <- rownames(meta)[meta$Gene_Control == "NTC"]
  higher <- rowMeans(lg[genes_k, ko_s, drop = FALSE]) > rowMeans(lg[genes_k, ntc_s, drop = FALSE])
  conc[[length(conc) + 1]] <- data.frame(KO = k, comparison = comp, n_genes = length(genes_k), n_higher_in_bulk = sum(higher), frac_higher = mean(higher))
}
conc <- do.call(rbind, conc); write.csv(conc, file.path(OUT, "spec34_gene_level_concordance.csv"), row.names = FALSE)
cat("\n=== per-gene concordance: fraction of the 50 screen-upregulated genes ALSO higher in bulk KO vs bulk NTC (pooled across lines; 0.5 = chance) ===\n")
print(format(conc[order(-conc$frac_higher), ], digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- heatmap: 16 modules (8 KO x 2 comparisons) x 76 bulk samples
save_fig("fig2Q_screen_top50_module_bulk_heatmap", 17, 9.5, function() {
  ord <- rownames(meta)[order(match(meta$Cell_Line, LINES), match(meta$Gene_Control, c("NTC", ORDER, setdiff(unique(meta$Gene_Control), c("NTC", ORDER)))))]
  row_ord <- as.vector(t(outer(ORDER, c("vsNTC", "vsRest"), paste)))
  M <- t(score[ord, row_ord, drop = FALSE]); n_r <- nrow(M); n_c <- ncol(M)
  par(mar = c(6.5, 7, 3.4, 1), family = "sans")
  vmax <- max(abs(M)); image(seq_len(n_c), seq_len(n_r), heat_z(M), col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
  cl_of <- meta[ord, "Cell_Line"]; grp_of <- meta[ord, "Gene_Control"]
  edges <- which(cl_of[-1] != cl_of[-n_c]); for (e in edges) abline(v = e + 0.5, col = "white", lwd = 2)
  ko_edges <- which(paste(cl_of, grp_of)[-1] != paste(cl_of, grp_of)[-n_c]); for (e in ko_edges) abline(v = e + 0.5, col = adjustcolor("white", 0.6), lwd = 0.7)
  row_edges <- which(1:(n_r - 1) %% 2 == 0); for (e in row_edges) abline(h = n_r - e + 0.5, col = "white", lwd = 1.4)
  axis(1, at = seq_len(n_c), labels = meta[ord, "Sample_ID"], las = 2, tick = FALSE, cex.axis = 0.32, col.axis = INK2)
  for (ri in seq_len(n_r)) { rn <- rownames(M)[ri]; k <- sub(" .*", "", rn); text(0.2, n_r - ri + 1, rn, adj = 1, xpd = NA, cex = 0.62, font = 2, col = KO_COL[k]) }
  for (i in seq_along(LINES)) { at <- mean(which(cl_of == LINES[i])); mtext(LINES[i], side = 3, at = at, line = 0.4, font = 2, cex = 0.8, col = LINE_COL[[LINES[i]]]) }
  box(col = MUTED, lwd = 0.6)
  mtext("Screen's top50 upregulated genes per KO (vs NTC and vs rest-of-guides), scored in bulk RNA-seq -- all 76 samples", side = 3, adj = 0, line = 2.4, font = 2, cex = 0.95)
  mtext("rows = one module per KO per comparison (paired rows); columns = bulk samples grouped by line then KO", side = 3, adj = 0, line = 1.2, cex = 0.55, col = INK2)
  draw_colorbar(c(0.03, 0.20), c(-vmax, vmax), "module activity z-score", cex = 0.55)
})
