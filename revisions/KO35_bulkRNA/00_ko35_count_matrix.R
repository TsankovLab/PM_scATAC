# 00_ko35_count_matrix.R: processed count matrices of the bulk RNA-seq KO cohort: 35 transcription-factor knockout samples (6 KOs: TEAD4, TWIST1, PITX1, SOX9, TEAD2, BPTF)
# + 8 non-targeting controls (NTC), mesothelioma lines MSTO-211H, NCI-H2452, NCI-H28; 43 samples, from the STAR GeneCounts files of raw/03_star_align.lsf.
#   sample       sample name (FASTQ / STAR files; <line>_<KO>_<replicate>) = column name of the matrices
#   raw counts   column 4 of <sample>_ReadsPerGene.out.tab (reverse-stranded library), STAR summary rows (N_*) dropped, all GENCODE v46 genes
#   normalised   TMM log2 CPM (edgeR: calcNormFactors on all 43 samples, cpm(log = TRUE, prior.count = 2)), all genes
#   annotation   sample, cell line, knocked-out gene (NTC = non-targeting control), replicate
# outputs (OUT_DIR): KO35_raw_counts.tsv.gz, KO35_TMM_log2CPM.tsv.gz (gene_id, gene_symbol, then the samples), KO35_samples.tsv, KO35_samples.csv
# Usage: [COUNTS_DIR=...] [OUT_DIR=...] Rscript 00_ko35_count_matrix.R      (R with edgeR)
source(local({a <- grep("^--file=", commandArgs(FALSE), value = TRUE); for (d in c(if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))), getwd())) { while (!file.exists(file.path(d, "config.R")) && d != dirname(d)) d <- dirname(d); if (file.exists(file.path(d, "config.R"))) break }; Sys.getenv("MESO_CONFIG", file.path(d, "config.R"))}))   # data / results / reference roots: see config.R at the repo root
suppressPackageStartupMessages(library(edgeR))
PROC <- Sys.getenv("PROC", paste0(DATA_DIR, "/CROPseq_validation_bulkRNA/04_processed_counts"))
COUNTS_DIR <- Sys.getenv("COUNTS_DIR", file.path(PROC, "star_counts")); OUT_DIR <- Sys.getenv("OUT_DIR", file.path(PROC, "KO35")); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
LINES <- c(MSTO = "MSTO-211H", H2452 = "NCI-H2452", H28 = "NCI-H28")
KO <- c("MSTO_TEAD4_R1", "MSTO_TEAD4_R2", "MSTO_TWIST1_R2", "MSTO_TWIST1_R3", "MSTO_PITX1_R1", "MSTO_PITX1_R3", "MSTO_SOX9_R2", "MSTO_SOX9_R3", "MSTO_TEAD2_R2", "MSTO_TEAD2_R3", "MSTO_BPTF_R1", "MSTO_BPTF_R3",
        "H2452_TEAD4_R1", "H2452_TEAD4_R3", "H2452_TWIST1_R1", "H2452_PITX1_R1", "H2452_PITX1_R3", "H2452_SOX9_R2", "H2452_SOX9_R3", "H2452_TEAD2_R1", "H2452_TEAD2_R2", "H2452_BPTF_R1", "H2452_BPTF_R2",
        "H28_TEAD4_R1", "H28_TEAD4_R3", "H28_TWIST1_R1", "H28_TWIST1_R3", "H28_PITX1_R1", "H28_PITX1_R3", "H28_SOX9_R1", "H28_SOX9_R2", "H28_TEAD2_R1", "H28_TEAD2_R2", "H28_BPTF_R1", "H28_BPTF_R3")
NTC <- c("MSTO_NTC_R1", "MSTO_NTC_R2", "MSTO_NTC_R3", "H2452_NTC_R1", "H2452_NTC_R2", "H2452_NTC_R3", "H28_NTC_R1", "H28_NTC_R2")
S <- c(KO, NTC)
rd <- function(s) { x <- read.delim(file.path(COUNTS_DIR, paste0(s, "_ReadsPerGene.out.tab")), header = FALSE, stringsAsFactors = FALSE); x <- x[!startsWith(x$V1, "N_"), ]; setNames(x$V4, x$V1) }
L <- lapply(setNames(S, S), rd); stopifnot(all(sapply(L, function(v) identical(names(v), names(L[[1]]))))); M <- do.call(cbind, L)
sym <- read.delim(file.path(COUNTS_DIR, "gene_id_to_symbol.tsv"), header = FALSE, col.names = c("gene_id", "symbol"))
Y <- cpm(calcNormFactors(DGEList(M)), log = TRUE, prior.count = 2)
wr <- function(X, f, digits) { d <- data.frame(gene_id = rownames(X), gene_symbol = sym$symbol[match(rownames(X), sym$gene_id)], if (is.null(digits)) X else round(X, digits), check.names = FALSE)
  con <- gzfile(file.path(OUT_DIR, f), "w"); write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con) }
wr(M, "KO35_raw_counts.tsv.gz", NULL); wr(Y, "KO35_TMM_log2CPM.tsv.gz", 4)
p <- strsplit(S, "_")
A <- data.frame(sample = S, cell_line = unname(LINES[sub("_.*", "", S)]), knockout = sapply(p, `[`, 2), replicate = sapply(p, `[`, 3))
write.table(A, file.path(OUT_DIR, "KO35_samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); write.csv(A, file.path(OUT_DIR, "KO35_samples.csv"), row.names = FALSE, quote = FALSE)
cat("genes:", nrow(M), "| samples:", ncol(M), "\n"); print(table(A$knockout, A$cell_line))
