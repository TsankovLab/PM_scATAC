###############################################################################
# STEP 1 -- per-cell fragments for every scATAC tumour, straight out of the
#           tumour-compartment ArchR project.
#
# epiAneufinder reads a headerless BED-like .tsv (chr, start, end, barcode). The only
# thing that matters here is WHICH cells go in: the ArchR project at
# tumor_compartment/scatac_ArchR holds the malignant cells only, so taking every cell
# of a sample from that project is exactly the malignant population -- no separate
# barcode list, no chance of the filter and the object disagreeing.
#
# The project is loaded once for all samples; loading is the slow part. Existing
# outputs are skipped, so the script is safe to re-run.
#
# Runtime ~10 min, ~11 GB of output.
# Output : frags/<S>_fragments_archr.tsv, frags/extraction_summary.csv
###############################################################################
suppressMessages({ library(ArchR); library(GenomicRanges); library(data.table) })
source("00_common.R")
addArchRThreads(4)
dir.create("frags", showWarnings = FALSE)
AUT <- paste0("chr", 1:22)

proj <- loadArchRProject(ARCHR, showLogo = FALSE)
cd <- getCellColData(proj)
cat("project cells:", nrow(cd), "\n"); print(table(cd$Sample))

## tumour samples only -- the project also holds RPL_* normal reference cells
ALL <- grep("^P[0-9]+$", sort(unique(as.character(cd$Sample))), value = TRUE)
cat("tumour samples:", paste(ALL, collapse = ", "), "\n")

summ <- list()
for (s in ALL) {
  f <- file.path("frags", sprintf("%s_fragments_archr.tsv", s))
  if (file.exists(f)) { cat(sprintf("%-5s exists, skipping\n", s)); next }
  cells <- rownames(cd)[cd$Sample == s]
  gr <- unlist(getFragmentsFromProject(proj, cellNames = cells), use.names = FALSE)
  gr <- gr[as.character(seqnames(gr)) %in% AUT]
  ## ArchR tags each fragment with "<sample>#<barcode>"; epiAneufinder wants the
  ## bare barcode in column 4
  bc <- sub("^.*#", "", as.character(mcols(gr)$RG))
  dt <- data.table(seqnames = as.character(seqnames(gr)), start = start(gr),
                   end = end(gr), barcode = bc)
  fwrite(dt, f, sep = "\t", col.names = FALSE)
  pc <- dt[, .N, by = barcode]$N
  summ[[s]] <- data.table(sample = s, cells = length(cells), frags = nrow(dt),
                          median_per_cell = as.integer(median(pc)),
                          cells_ge_minfrags = sum(pc >= MINFRAGS))
  cat(sprintf("%-5s %5d cells | %10d frags | median/cell %6d | >=%d frags: %5d\n",
              s, length(cells), nrow(dt), as.integer(median(pc)), MINFRAGS, sum(pc >= MINFRAGS)))
  rm(gr, dt); invisible(gc())
}
if (length(summ)) fwrite(rbindlist(summ), "frags/extraction_summary.csv")
cat("\nDONE\n")
