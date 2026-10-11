# Spec 71b (step 1 of 3): regions and BAM list for the SNP-based cell-line identity check of the bulk RNA-seq samples.
# Regions = exons (GENCODE v46, autosomes) of the NGENES protein-coding genes with the highest MINIMUM-over-lines mean CPM, so that every line has
# coverage at the same positions (genotype concordance is only comparable on sites covered in both samples). Samples = every BAM whose name carries
# one of the 4 mesothelioma line prefixes (incl. SNAI2/RUNX1/RUNX2 KOs and the MIXED_H28_* samples: more reference samples per line, same cost).
# env: GENO_DIR (work folder, default scratch work/genotype), BAM_DIR, NGENES (default 3000)
# needs: counts/gene_counts_matrix_stranded_reverse.tsv, gencode.v46.annotation.gtf
# outputs (in GENO_DIR): regions.bed, bam_list.txt, sample_names.txt
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R"))
suppressPackageStartupMessages(library(data.table))
GENO_DIR <- Sys.getenv("GENO_DIR", unset = "/sc/arion/scratch/giottb01/CROPseq_validation_bulkRNA/work/genotype"); dir.create(GENO_DIR, showWarnings = FALSE, recursive = TRUE)
BAM_DIR  <- Sys.getenv("BAM_DIR", unset = "/sc/arion/scratch/giottb01/CROPseq_validation_bulkRNA/work/bam")
GTF <- "/sc/arion/projects/Tsankov_Normal_Lung/genome_assemblies/human/gencode.v46.annotation.gtf"; NGENES <- as.integer(Sys.getenv("NGENES", unset = "3000"))

ids <- list.files(BAM_DIR); ids <- ids[grepl("^(H2052|MSTO|H2452|H28)_|^MIXED_H28_", ids)]
bam <- file.path(BAM_DIR, ids, paste0(ids, "_Aligned.sortedByCoord.out.bam")); ok <- file.exists(bam); ids <- ids[ok]; bam <- bam[ok]
writeLines(bam, file.path(GENO_DIR, "bam_list.txt")); writeLines(ids, file.path(GENO_DIR, "sample_names.txt"))
message(length(ids), " BAMs listed")

counts <- as.matrix(read.delim(file.path(PROC, "counts/gene_counts_matrix_stranded_reverse.tsv"), row.names = 1, check.names = FALSE))
counts <- counts[, intersect(ids, colnames(counts))]; cpm <- sweep(counts, 2, colSums(counts), "/") * 1e6
pfx <- sub("^MIXED_", "", colnames(cpm)); pfx <- sub("_.*", "", pfx)                       # as-named line prefix; only used to balance coverage
minline <- apply(sapply(unique(pfx), function(p) rowMeans(cpm[, pfx == p, drop = FALSE])), 1, min)

g <- fread(cmd = sprintf("awk -F'\\t' '$3==\"exon\" && $1 ~ /^chr[0-9]+$/' %s | grep 'gene_type \"protein_coding\"'", GTF), header = FALSE, sep = "\t", select = c(1, 4, 5, 9))
setnames(g, c("chr", "start", "end", "attr")); g[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", attr)]
top <- names(sort(minline[names(minline) %in% g$gene_id], decreasing = TRUE))[seq_len(NGENES)]
e <- unique(g[gene_id %in% top, .(chr, start = start - 1L, end)]); setorder(e, chr, start, end)
e[, grp := cumsum(c(TRUE, start[-1] > cummax(end)[-.N])), by = chr]; bed <- e[, .(start = min(start), end = max(end)), by = .(chr, grp)][, grp := NULL]   # merge overlapping exons
bed <- bed[order(as.integer(sub("chr", "", chr)), start)]
fwrite(bed, file.path(GENO_DIR, "regions.bed"), sep = "\t", col.names = FALSE)
message(sprintf("regions: %d merged exons of %d genes, %.2f Mb (min-over-lines mean CPM of the last gene: %.1f)", nrow(bed), length(top), sum(bed$end - bed$start) / 1e6, min(minline[top])))
