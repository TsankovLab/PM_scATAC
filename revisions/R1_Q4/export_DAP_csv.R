# export_DAP_csv.R  (R1_Q4)
# Reads the already-computed DAP RDS files and exports significant-peak CSV tables.
# Run on login node: Rscript export_DAP_csv.R

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
})

out_dir        <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/meso_vs_pbmc"
marker_targets <- c("CD14_Mono", "CD16_Mono", "DC")

message("Exporting DA peak tables (FDR < 0.05, Log2FC > 1)...")
for (tgt in marker_targets) {
  rds_out <- file.path(out_dir, paste0("DAP_", tgt, "_vs_rest.rds"))
  csv_out <- file.path(out_dir, paste0("DAP_", tgt, "_sig.csv"))
  if (!file.exists(rds_out)) { message("  Missing: ", basename(rds_out)); next }

  mkr <- readRDS(rds_out)
  col <- colnames(mkr)[1]   # e.g. "CD14_Mono"
  rd  <- as.data.frame(rowData(mkr))
  res <- data.frame(
    Log2FC = as.numeric(assays(mkr)[["Log2FC"]][[col]]),
    FDR    = as.numeric(assays(mkr)[["FDR"]][[col]]),
    region = paste0(rd$seqnames, ":", rd$start, "-", rd$end),
    stringsAsFactors = FALSE
  )
  sig <- res[!is.na(res$FDR) & res$FDR < 0.05 & res$Log2FC > 1, ]
  sig <- sig[order(-sig$Log2FC), ]
  message(sprintf("  %s: %d significant peaks", tgt, nrow(sig)))
  if (nrow(sig) > 0)
    write.csv(sig, csv_out, row.names = FALSE, quote = FALSE)
}
message("Done.")
