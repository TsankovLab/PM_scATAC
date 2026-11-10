# 01_nr4a2_deseq2.R: NR4A2 enhancer CRISPRa (guides NR4A2_1, NR4A2_2) vs SCRAMBLE (1-3) in Jurkat cells, from the STAR gene counts (raw/04-05).
# Reconstruction of the analysis behind processed/de/NR4A2_1_2_vs_SCRAMBLE (2026-09-15; the original script was not kept). Verified to reproduce the saved
# dds.rds, DESeq2_results_full.csv and the GSEA input exactly (R 4.4.1, DESeq2 1.46.0).
#   counts   gene_counts_matrix_stranded_reverse.tsv (STAR column 4, reverse-stranded library), the 5 samples; genes with >= 10 reads summed over them
#   design   ~ cc_score + condition (SCRAMBLE reference); DESeq() defaults; Wald test of condition_NR4A2_vs_SCRAMBLE
#   results  results(alpha = 0.05) (independent filtering targeted at padj 0.05); log2FC and lfcSE from lfcShrink(type = "normal", contrast = NR4A2 vs SCRAMBLE),
#            stat / pvalue / padj from the unshrunk Wald test; symbol from gene_id_to_symbol.tsv; rows ordered by padj (NA last, ties in
#            DESeq2 row order)
#   cc_score per-sample cell-cycle covariate. The values are those stored in the original dds.rds: their recipe was not recorded and is not one of the
#            standard ones (Seurat S + G2M, mean z of the Seurat S / G2M genes, PC1 of those genes: best match r = 0.98). They are fixed here so the
#            published results reproduce exactly.
#   the ranked gene list for GSEA is built from DESeq2_results_full.csv by 02_nr4a2_gsea.py
# outputs (DE_DIR): dds.rds, DESeq2_results_full.csv, DESeq2_results_padj0.05.csv, normalized_counts.csv
# Usage: [DE_DIR=...] [PROC=...] Rscript 01_nr4a2_deseq2.R      (module load R/4.4.1)
suppressPackageStartupMessages(library(DESeq2))
PROC <- Sys.getenv("PROC", "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CROPseq_validation_bulkRNA/processed")
DE_DIR <- Sys.getenv("DE_DIR", file.path(PROC, "de", "NR4A2_1_2_vs_SCRAMBLE")); dir.create(file.path(DE_DIR, "gsea"), recursive = TRUE, showWarnings = FALSE)
S5 <- c("NR4A2_1", "NR4A2_2", "SCRAMBLE_1", "SCRAMBLE_2", "SCRAMBLE_3")
CC <- c(NR4A2_1 = -0.081381560609297979, NR4A2_2 = 0.36272718958156569, SCRAMBLE_1 = 0.73812092269449414, SCRAMBLE_2 = -0.64887578974909987,
        SCRAMBLE_3 = -0.37059076191765838)   # full double precision, as stored in the original dds.rds
cnt <- read.delim(file.path(PROC, "counts", "gene_counts_matrix_stranded_reverse.tsv"), row.names = 1, check.names = FALSE)
sym <- read.delim(file.path(PROC, "counts", "gene_id_to_symbol.tsv"), header = FALSE, col.names = c("gene_id", "symbol"))
x <- as.matrix(cnt[, S5]); x <- x[rowSums(x) >= 10, ]; cat("genes with >= 10 reads:", nrow(x), "\n")
cd <- data.frame(condition = factor(c("NR4A2", "NR4A2", "SCRAMBLE", "SCRAMBLE", "SCRAMBLE"), levels = c("SCRAMBLE", "NR4A2")), cc_score = CC[S5], row.names = S5)
dds <- DESeq(DESeqDataSetFromMatrix(x, cd, ~ cc_score + condition), quiet = TRUE); saveRDS(dds, file.path(DE_DIR, "dds.rds"))
res <- results(dds, name = "condition_NR4A2_vs_SCRAMBLE", alpha = 0.05)
shr <- lfcShrink(dds, contrast = c("condition", "NR4A2", "SCRAMBLE"), res = res, type = "normal", quiet = TRUE)
out <- data.frame(gene_id = rownames(res), symbol = sym$symbol[match(rownames(res), sym$gene_id)], baseMean = res$baseMean, log2FoldChange = shr$log2FoldChange,
                  lfcSE = shr$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj); out <- out[order(out$padj), ]
write.csv(out, file.path(DE_DIR, "DESeq2_results_full.csv"), row.names = FALSE)
write.csv(out[!is.na(out$padj) & out$padj < 0.05, ], file.path(DE_DIR, "DESeq2_results_padj0.05.csv"), row.names = FALSE)
write.csv(counts(dds, normalized = TRUE), file.path(DE_DIR, "normalized_counts.csv"))
cat(sprintf("padj < 0.05: %d genes (%d up, %d down)\n", sum(out$padj < 0.05, na.rm = TRUE), sum(out$padj < 0.05 & out$log2FoldChange > 0, na.rm = TRUE),
            sum(out$padj < 0.05 & out$log2FoldChange < 0, na.rm = TRUE)))
