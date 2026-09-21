suppressMessages({ library(Seurat) })
set.seed(1234)
PROC <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CROPseq_validation_bulkRNA/processed"
OUT <- Sys.getenv("CROP_OUT", unset = "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CROPseq_validation_bulkRNA/processed/de/perturbation_panel/figs_repro")

counts <- read.delim(file.path(PROC, "counts/gene_counts_matrix_stranded_reverse.tsv"), row.names = 1, check.names = FALSE)
sym <- read.delim(file.path(PROC, "counts/gene_id_to_symbol.tsv"), header = FALSE, col.names = c("gene_id", "symbol"))
sym <- sym[sym$gene_id %in% rownames(counts) & !is.na(sym$symbol) & sym$symbol != "", ]
cs <- rowsum(as.matrix(counts[sym$gene_id, ]), group = sym$symbol)   # sum counts per symbol
cat("symbol-level matrix:", dim(cs), "\n")

meta <- read.delim(file.path(PROC, "sample_sheet/sample_metadata_RUN1.tsv"))
meta <- meta[meta$status == "OK" & !meta$Cell_Line %in% c("Jurkat", "Lily"), ]
rownames(meta) <- meta$fastq_sample_id

g <- Seurat::cc.genes
cat("S genes present:", sum(g$s.genes %in% rownames(cs)), "/", length(g$s.genes),
    " G2M present:", sum(g$g2m.genes %in% rownames(cs)), "/", length(g$g2m.genes), "\n")

out <- list()
for (cl in c("MSTO-211H", "NCI-H2052", "NCI-H2452", "NCI-H28")) {
  samples <- meta$fastq_sample_id[meta$Cell_Line == cl]
  obj <- CreateSeuratObject(counts = cs[, samples], min.cells = 0, min.features = 0)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- CellCycleScoring(obj, s.features = g$s.genes, g2m.features = g$g2m.genes, set.ident = FALSE)
  d <- obj@meta.data[, c("S.Score", "G2M.Score", "Phase")]
  d$cc <- d$S.Score + d$G2M.Score
  d$sample <- rownames(d); d$Cell_Line <- cl; d$Gene_Control <- meta[rownames(d), "Gene_Control"]
  out[[cl]] <- d
}
res <- do.call(rbind, out); rownames(res) <- NULL
write.csv(res, file.path(OUT, "bulk_cellcycle_per_sample_RUN1.csv"), row.names = FALSE)
print(aggregate(cc ~ Cell_Line + Gene_Control, res, mean), digits = 3)
