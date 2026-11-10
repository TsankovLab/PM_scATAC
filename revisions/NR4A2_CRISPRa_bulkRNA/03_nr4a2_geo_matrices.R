# 03_nr4a2_geo_matrices.R: processed data files for GEO, NR4A2 enhancer CRISPRa (NR4A2_1, NR4A2_2) vs SCRAMBLE (1-3), Jurkat, the 5 samples of the analysis.
#   raw counts         STAR GeneCounts, reverse-stranded column (gene_counts_matrix_stranded_reverse.tsv), all GENCODE v46 genes
#   normalised counts  raw counts / DESeq2 size factors (median-of-ratios, estimated on the genes with >= 10 reads over the 5 samples, as dds.rds of 01;
#                      the same normalisation as the differential expression), all genes, not log-transformed
#   both: tab-separated, gzipped, columns gene_id (Ensembl, with version), gene_symbol (GENCODE v46 gene_name), then the 5 samples
# needs: 01_nr4a2_deseq2.R output (dds.rds); outputs (DE_DIR/GEO): NR4A2_CRISPRa_Jurkat_raw_counts.tsv.gz, NR4A2_CRISPRa_Jurkat_DESeq2_normalized_counts.tsv.gz
# Usage: [DE_DIR=...] [PROC=...] Rscript 03_nr4a2_geo_matrices.R      (module load R/4.4.1)
suppressPackageStartupMessages(library(DESeq2))
PROC <- Sys.getenv("PROC", "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CROPseq_validation_bulkRNA/processed")
DE_DIR <- Sys.getenv("DE_DIR", file.path(PROC, "de", "NR4A2_1_2_vs_SCRAMBLE")); GD <- file.path(DE_DIR, "GEO"); dir.create(GD, showWarnings = FALSE)
dds <- readRDS(file.path(DE_DIR, "dds.rds")); S5 <- colnames(dds); sf <- sizeFactors(dds)
cnt <- as.matrix(read.delim(file.path(PROC, "counts", "gene_counts_matrix_stranded_reverse.tsv"), row.names = 1, check.names = FALSE)[, S5])
sym <- read.delim(file.path(PROC, "counts", "gene_id_to_symbol.tsv"), header = FALSE, col.names = c("gene_id", "symbol"))
stopifnot(identical(cnt[rownames(dds), ], counts(dds)))          # same counts as the DE object
nrm <- sweep(cnt, 2, sf[S5], "/"); stopifnot(isTRUE(all.equal(nrm[rownames(dds), ], counts(dds, normalized = TRUE))))
wr <- function(M, f, digits) { d <- data.frame(gene_id = rownames(M), gene_symbol = sym$symbol[match(rownames(M), sym$gene_id)], if (is.null(digits)) M else round(M, digits), check.names = FALSE)
  con <- gzfile(file.path(GD, f), "w"); write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con) }
wr(cnt, "NR4A2_CRISPRa_Jurkat_raw_counts.tsv.gz", NULL); wr(nrm, "NR4A2_CRISPRa_Jurkat_DESeq2_normalized_counts.tsv.gz", 4)
cat("genes:", nrow(cnt), "| samples:", paste(S5, collapse = ", "), "\nsize factors:", paste(sprintf("%s %.4f", S5, sf[S5]), collapse = ", "), "\nwritten:", GD, "\n")
