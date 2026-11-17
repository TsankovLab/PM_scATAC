# 00_nr4a2_count_matrix.R: gene x sample count matrix of the NR4A2 enhancer CRISPRa experiment (Jurkat; NR4A2_1, NR4A2_2 vs SCRAMBLE_1-3) from the STAR
# GeneCounts files of raw/03_star_align.lsf: column 4 of <sample>_ReadsPerGene.out.tab (reads on the reverse strand of the gene: reverse-stranded library),
# STAR summary rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous) dropped; genes in GENCODE v46 order.
# output: DE_DIR/NR4A2_CRISPRa_counts.tsv (gene_id, then the 5 samples), the input of 01_nr4a2_deseq2.R and 03_nr4a2_geo_matrices.R
# Usage: [COUNTS_DIR=...] [DE_DIR=...] Rscript 00_nr4a2_count_matrix.R
source(local({a <- grep("^--file=", commandArgs(FALSE), value = TRUE); for (d in c(if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))), getwd())) { while (!file.exists(file.path(d, "config.R")) && d != dirname(d)) d <- dirname(d); if (file.exists(file.path(d, "config.R"))) break }; Sys.getenv("MESO_CONFIG", file.path(d, "config.R"))}))   # data / results / reference roots: see config.R at the repo root
PROC <- Sys.getenv("PROC", paste0(DATA_DIR, "/CROPseq_validation_bulkRNA/processed"))
COUNTS_DIR <- Sys.getenv("COUNTS_DIR", file.path(PROC, "counts")); DE_DIR <- Sys.getenv("DE_DIR", file.path(PROC, "de", "NR4A2_1_2_vs_SCRAMBLE")); dir.create(DE_DIR, recursive = TRUE, showWarnings = FALSE)
S5 <- c("NR4A2_1", "NR4A2_2", "SCRAMBLE_1", "SCRAMBLE_2", "SCRAMBLE_3")
rd <- function(s) { x <- read.delim(file.path(COUNTS_DIR, paste0(s, "_ReadsPerGene.out.tab")), header = FALSE, stringsAsFactors = FALSE); x <- x[!startsWith(x$V1, "N_"), ]; setNames(x$V4, x$V1) }
L <- lapply(setNames(S5, S5), rd); stopifnot(all(sapply(L, function(v) identical(names(v), names(L[[1]])))))
M <- do.call(cbind, L); write.table(data.frame(gene_id = rownames(M), M, check.names = FALSE), file.path(DE_DIR, "NR4A2_CRISPRa_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("genes:", nrow(M), "| reads assigned to genes:", paste(sprintf("%s %.1f M", S5, colSums(M) / 1e6), collapse = ", "), "\n")
