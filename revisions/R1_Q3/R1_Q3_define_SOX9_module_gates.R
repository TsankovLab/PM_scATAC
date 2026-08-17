#!/usr/bin/env Rscript
# R1_Q3_define_SOX9_module_gates.R
#
# Define conservative SOX9+ / SOX9- cell sets from the SOX9-program module score
# (top-20 P4 cluster-32 genes) on the P23 ATAC tumor cells, for rebuilding
# purity-gated ChromBPNet pseudobulks. Gating is on the module score across ALL
# tumor cells (not the C9/C4 cluster labels), so high-program cells misassigned
# to C4 are included and low-program cells in C9 are excluded.
#
# Gates:
#   SOX9+ (conservative) : module >= POS  (default 0.20)
#   SOX9- (clean)        : module <= NEG  (default -0.10)
#   middle band excluded to maximize contrast.

suppressMessages(library(data.table)); setDTthreads(1)
OUT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
POS <- 0.20; NEG <- -0.10

m <- fread(file.path(OUT, "SOX9_program_module_P23_percell.csv"))
p <- fread(file.path(OUT, "SOX9_C9_purity_percell.csv"))
d <- merge(m, p[, .(cell, gene_SOX9, motif_SOX9_z)], by = "cell")

d[, gate := fifelse(module >= POS, "SOX9pos",
             fifelse(module <= NEG, "SOX9neg", "excluded"))]

summ <- d[, .(n = .N,
              pct_fromC9   = round(100*mean(SOX9_cluster=="SOX9_high_P23"),1),
              median_module= round(median(module),3),
              median_SOX9gene = round(median(gene_SOX9),3),
              pct_SOX9gene_pos= round(100*mean(gene_SOX9>0),1),
              median_motifz   = round(median(motif_SOX9_z,na.rm=TRUE),3)),
          by = gate][order(-median_module)]
cat(sprintf("Gates: SOX9pos module>=%.2f , SOX9neg module<=%.2f\n\n", POS, NEG))
print(summ)

# write barcode lists (ArchR cellNames = Sample#barcode) + labelled table
for (g in c("SOX9pos","SOX9neg")) {
  bc <- d[gate==g, cell]
  writeLines(bc, file.path(OUT, sprintf("SOX9_%s_module_barcodes.txt", g)))
  # bare barcode (strip Sample# prefix) for fragment-file matching
  writeLines(sub("^.*#", "", bc), file.path(OUT, sprintf("SOX9_%s_module_barcodes_bare.txt", g)))
}
fwrite(d[, .(cell, SOX9_cluster, module, gene_SOX9, motif_SOX9_z, gate)],
       file.path(OUT, "SOX9_module_gated_cells.csv"))

cat("\nExample cellNames:\n"); print(head(d$cell, 3))
cat(sprintf("\nWrote barcode lists + SOX9_module_gated_cells.csv to %s\n", OUT))
cat(sprintf("  SOX9pos: %d cells   SOX9neg: %d cells   excluded: %d\n",
            sum(d$gate=="SOX9pos"), sum(d$gate=="SOX9neg"), sum(d$gate=="excluded")))
