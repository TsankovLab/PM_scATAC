###############################################################################
# STEP 5 -- do the two annotations recover the same tissue composition?
#
# The third, weakest-assumption test.  It needs no integration and no per-cell
# correspondence: for each of the 10 patients profiled by both assays, take the
# fraction of cells each annotation assigns to each cell type and compare them.
#
# This is the question a reviewer actually cares about for downstream biology --
# if the chromatin annotation says a tumour is 60% malignant and 20% T cells and
# the transcriptomic annotation of the same tumour says the same, the annotations
# are interchangeable at the level any composition analysis uses them.
#
# It is deliberately a weak test of the annotation and a strong test of the assay:
# disagreement here can come from the annotation OR from the dissociation and
# nucleus-extraction biases that differ between scATAC and scRNA, and the two
# cannot be separated with these data.  Reported as such.
#
# Input : atac_cells.csv, rna_cells.csv
# Output: composition_by_patient.csv, composition_concordance.csv
###############################################################################
source("00_common.R")

A <- read.csv("atac_cells.csv", stringsAsFactors = FALSE)
R <- read.csv("rna_cells.csv",  stringsAsFactors = FALSE)
A <- A[A$sample %in% SAMPLES_SHARED, ]; R <- R[R$sample %in% SAMPLES_SHARED, ]

frac <- function(d, lab){
  tb <- table(d$sample, d[[lab]])
  as.data.frame.table(prop.table(tb, 1), responseName = "frac",
                      stringsAsFactors = FALSE)
}
fa <- frac(A, "atac_label"); names(fa)[1:2] <- c("sample", "celltype")
fr <- frac(R, "rna_label");  names(fr)[1:2] <- c("sample", "celltype")
CP <- merge(fa, fr, by = c("sample", "celltype"), all = TRUE,
            suffixes = c("_atac", "_rna"))
CP[is.na(CP)] <- 0
CP$n_atac <- as.integer(round(CP$frac_atac * table(A$sample)[CP$sample]))
CP$n_rna  <- as.integer(round(CP$frac_rna  * table(R$sample)[CP$sample]))
CP <- CP[order(CP$sample, -CP$frac_atac), ]
write.csv(CP, "composition_by_patient.csv", row.names = FALSE)

cat("=== composition, 10 shared patients x", length(unique(CP$celltype)), "cell types ===\n")
cat("pairs:", nrow(CP), "\n\n")

## overall agreement across all patient x cell type pairs
p_all <- cor(CP$frac_atac, CP$frac_rna, method = "pearson")
s_all <- cor(CP$frac_atac, CP$frac_rna, method = "spearman")
cat(sprintf("all pairs           : Pearson r = %.3f | Spearman rho = %.3f | n = %d\n",
            p_all, s_all, nrow(CP)))

## per patient: does the ranking of cell types agree within a tumour?
per_pt <- do.call(rbind, lapply(split(CP, CP$sample), function(x)
  data.frame(sample = x$sample[1], n_types = nrow(x),
             pearson  = cor(x$frac_atac, x$frac_rna),
             spearman = cor(x$frac_atac, x$frac_rna, method = "spearman"),
             max_abs_diff = max(abs(x$frac_atac - x$frac_rna)),
             worst_type = x$celltype[which.max(abs(x$frac_atac - x$frac_rna))],
             stringsAsFactors = FALSE)))
cat("\n=== per patient ===\n")
print(transform(per_pt, pearson = round(pearson, 3), spearman = round(spearman, 3),
                max_abs_diff = round(max_abs_diff, 3)), row.names = FALSE)

## per cell type: is a type systematically over- or under-represented in one assay?
per_ct <- do.call(rbind, lapply(split(CP, CP$celltype), function(x)
  data.frame(celltype = x$celltype[1], n_patients = nrow(x),
             mean_atac = mean(x$frac_atac), mean_rna = mean(x$frac_rna),
             mean_diff = mean(x$frac_atac - x$frac_rna),
             pearson = if (sd(x$frac_atac) > 0 && sd(x$frac_rna) > 0)
                         cor(x$frac_atac, x$frac_rna) else NA,
             p_paired = tryCatch(wilcox.test(x$frac_atac, x$frac_rna,
                                             paired = TRUE)$p.value,
                                 error = function(e) NA),
             stringsAsFactors = FALSE)))
per_ct <- per_ct[order(-abs(per_ct$mean_diff)), ]
per_ct$fdr <- p.adjust(per_ct$p_paired, "BH")
cat("\n=== per cell type (mean fraction over the 10 patients) ===\n")
print(transform(per_ct, mean_atac = round(mean_atac, 4), mean_rna = round(mean_rna, 4),
                mean_diff = round(mean_diff, 4), pearson = round(pearson, 3),
                p_paired = signif(p_paired, 2), fdr = signif(fdr, 2)), row.names = FALSE)

write.csv(rbind(
  data.frame(level = "all", id = "all pairs", n = nrow(CP),
             pearson = p_all, spearman = s_all, stringsAsFactors = FALSE),
  data.frame(level = "patient", id = per_pt$sample, n = per_pt$n_types,
             pearson = per_pt$pearson, spearman = per_pt$spearman),
  data.frame(level = "celltype", id = per_ct$celltype, n = per_ct$n_patients,
             pearson = per_ct$pearson, spearman = NA)),
  "composition_concordance.csv", row.names = FALSE)

cat("\nDONE\n")
