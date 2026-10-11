# Spec 41: which of the 30 cNMF modules (k=30, CROP-seq screen's own cNMF sweep) are enriched for EMT genes? Tests each module's
# top-N defining genes (by gene_spectra_score) for overrepresentation of the MSigDB Hallmark EPITHELIAL_MESENCHYMAL_TRANSITION gene
# set (hypergeometric test), background = the ~5000 overdispersed genes cNMF was actually run on (cnmf.overdispersed_genes.txt).
# needs: CRISPR_cropseq_analysis/.../cNMF_crispr_20_30_vf5000/cnmf/cnmf.gene_spectra_score.k_30.dt_0_3.txt,
#        cnmf.overdispersed_genes.txt, msigdbr (Hallmark gene sets)
# outputs: spec41_cNMF_EMT_enrichment.csv, fig2Y_cNMF_EMT_enrichment.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages({ library(msigdbr) })
CNMF_D <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CRISPR_cropseq_analysis/_cellranger_raw_Filter_400_800_25/no_harmony/cNMF/cNMF_crispr_20_30_vf5000/cnmf"
K <- "30"; MODULES <- as.character(1:30); TOPN <- as.integer(Sys.getenv("EMT_TOPN", "100"))

gss <- read.delim(file.path(CNMF_D, sprintf("cnmf.gene_spectra_score.k_%s.dt_0_3.txt", K)), row.names = 1, check.names = FALSE)
universe <- readLines(file.path(CNMF_D, "cnmf.overdispersed_genes.txt"))
universe <- intersect(universe, colnames(gss))
cat(sprintf("universe (cNMF overdispersed genes, present in gene_spectra_score): %d genes\n", length(universe)))

emt <- msigdbr(species = "Homo sapiens", category = "H")
emt <- emt[emt$gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ]
emt_genes <- intersect(unique(emt$gene_symbol), universe)
cat(sprintf("Hallmark EMT gene set: %d genes total, %d present in the cNMF universe\n", length(unique(emt$gene_symbol)), length(emt_genes)))

top_genes_list <- lapply(MODULES, function(m) {
  x <- sort(setNames(as.numeric(gss[m, universe]), universe), decreasing = TRUE)
  names(x)[1:TOPN]
})
names(top_genes_list) <- MODULES

res <- do.call(rbind, lapply(MODULES, function(m) {
  g <- top_genes_list[[m]]; hit_genes <- intersect(g, emt_genes); k <- length(hit_genes)
  N <- length(universe); K_ <- length(emt_genes); n <- length(g)
  p <- phyper(k - 1, K_, N - K_, n, lower.tail = FALSE)
  expected <- n * K_ / N
  data.frame(module = m, n_top = n, n_EMT_genes_in_set = K_, overlap = k, expected = round(expected, 2),
             fold_enrichment = round(k / expected, 2), p_value = p, hit_genes = paste(hit_genes, collapse = ";"))
}))
res$padj <- p.adjust(res$p_value, method = "BH")
res <- res[order(res$p_value), ]
write.csv(res, file.path(OUT, sprintf("spec41_cNMF_EMT_enrichment_top%d.csv", TOPN)), row.names = FALSE)

options(width = 200)
cat(sprintf("\n=== EMT (Hallmark) enrichment across all 30 cNMF modules (k=%s, top%d genes per module) ===\n", K, TOPN))
print(format(res[, c("module", "overlap", "expected", "fold_enrichment", "p_value", "padj")], digits = 3), row.names = FALSE)
sig <- res[res$padj < 0.05, ]
cat(sprintf("\n%d/%d modules significantly EMT-enriched (BH padj < 0.05):\n", nrow(sig), nrow(res)))
if (nrow(sig)) for (i in seq_len(nrow(sig))) cat(sprintf("  module %s: %d/%d genes overlap EMT (%.1fx enrichment, padj=%.2e) -- %s\n",
    sig$module[i], sig$overlap[i], sig$n_top[i], sig$fold_enrichment[i], sig$padj[i], sig$hit_genes[i]))

## ---------------------------------------------------------------- figure: -log10(padj) per module, EMT hits highlighted
save_fig(sprintf("fig2Y_top%d_cNMF_EMT_enrichment", TOPN), 11, 6, function() {
  ord <- res[order(as.numeric(res$module)), ]
  cols <- ifelse(ord$padj < 0.05, "#eb6834", adjustcolor(INK2, 0.4))
  par(mar = c(5, 4.5, 3.6, 1), family = "sans")
  bp <- barplot(-log10(ord$padj), names.arg = paste0("m", ord$module), col = cols, border = INK2, las = 2, cex.names = 0.7,
                ylab = "-log10(BH-adjusted p), Hallmark EMT enrichment")
  abline(h = -log10(0.05), col = INK, lty = 2)
  text(bp, -log10(ord$padj), ifelse(ord$padj < 0.05, ord$overlap, ""), pos = 3, cex = 0.65, xpd = NA)
  mtext(sprintf("Which cNMF modules (k=%s) are enriched for Hallmark EMT genes?", K), side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext(sprintf("hypergeometric test, top%d genes/module vs %d-gene cNMF universe (%d Hallmark EMT genes in universe); dashed line = padj 0.05; numbers = overlap count",
                TOPN, length(universe), length(emt_genes)), side = 3, adj = 0, line = 0.9, cex = 0.58, col = INK2)
})
