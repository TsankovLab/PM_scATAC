###############################################################################
# STEP 2 -- epiAneufinder on one sample.  Usage:  Rscript 02_run_epianeufinder.R P4
#
# All parameters come from 00_common.R so every sample is called identically and the
# 9 result tables share the same 298-window genome grid (required by step 3, which
# intersects the bins across samples).
#
# Runtime 10-60 min per sample depending on depth (P23 is the slow one).
# Output : out_5Mb/<S>/epiAneufinder_results/{results_table.tsv, cnv_calls.rds, Karyogram.png, ...}
###############################################################################
suppressMessages(library(epiAneufinder))
source("00_common.R")
S <- commandArgs(TRUE)[1]
stopifnot(!is.na(S), nzchar(S))
dir.create("out_5Mb", showWarnings = FALSE)
OUT <- file.path("out_5Mb", S)
cat("sample", S, "->", OUT, "\n")

epiAneufinder(
  input       = sprintf("frags/%s_fragments_archr.tsv", S),
  outdir      = OUT,
  blacklist   = BLACKLIST,
  windowSize  = WINDOW,
  genome      = GENOME,
  exclude     = EXCLUDE,
  minFrags    = MINFRAGS,
  ncores      = NCORES,
  k           = KVAL,
  gc_correction = TRUE,
  plotKaryo   = TRUE,
  title_karyo = sprintf("%s malignant - epiAneufinder %g Mb", S, WINDOW/1e6))

cat("\nDONE", S, "\n")
