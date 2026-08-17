#!/usr/bin/env Rscript
# R1_Q3_mp3_SOX9_overlap.R
#
# Does the HDMA merged motif Average_65__merged_pattern_3 (mp3; putative
# C/EBP / bZIP palindrome, NOT TOMTOM-called SOX) co-localize with SOX9
# chromVAR motif hits from the ArchR tumor_compartment (SOX9 P23) project?
#
# SOX9 "hits" = genomic positions of the SOX9 motif within ArchR peaks
#   (ArchR Annotations/Motif-Positions-In-Peaks.rds, i.e. the matchMotifs
#    output that chromVAR deviations are computed from).
# mp3 instances = FiNeMo hit calls of Average_65__merged_pattern_3 inside the
#   ChromBPNet SOX9-high and SOX9-low peaksets.
#
# Tests (per ChromBPNet cell type and pooled):
#   (1) instance-level: fraction of mp3 instances within +/-FLANK bp of a SOX9 motif
#   (2) peak-level: are mp3+ peaks enriched for SOX9 motifs vs all model peaks (Fisher)
# mp0 (Average_65__merged_pattern_0, strong C/EBP) and the full merged-motif
# hit set are included as comparison anchors.

suppressMessages({
  library(GenomicRanges)
  library(data.table)
})

P     <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
MESO  <- "/sc/arion/scratch/giottb01/meso_SOX9_P23_chromBPnet"
OUT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
FLANK <- 50L   # bp window around a SOX9 motif centre for instance-level overlap

CTS <- c(high = "SOX9_high_P23", low = "SOX9_low_P23")

## ---- 1. SOX9 chromVAR motif positions from ArchR -------------------------
message("Loading ArchR Motif-Positions-In-Peaks.rds ...")
pos <- readRDS(file.path(P, "Annotations", "Motif-Positions-In-Peaks.rds"))
sox9_names <- grep("SOX9", names(pos), value = TRUE, ignore.case = TRUE)
message("SOX9 motif entries in annotation: ", paste(sox9_names, collapse = ", "))
stopifnot(length(sox9_names) >= 1)
sox9 <- unlist(GRangesList(pos[sox9_names]), use.names = FALSE)
sox9 <- reduce(resize(sox9, width = width(sox9) + 2L * FLANK, fix = "center"))
seqlevelsStyle(sox9) <- "UCSC"
message("SOX9 motif positions (merged, +/-", FLANK, "bp): ", length(sox9))

## ---- helpers --------------------------------------------------------------
read_peaks <- function(ct) {
  f <- file.path(MESO, ct, paste0(ct, "_peakset_all_no_blacklist.bed"))
  d <- fread(f, header = FALSE, select = 1:3,
             col.names = c("chr", "start", "end"))
  gr <- GRanges(d$chr, IRanges(d$start + 1L, d$end))  # BED 0-based -> 1-based
  seqlevelsStyle(gr) <- "UCSC"
  gr
}

read_hits <- function(ct, motif_regex) {
  f <- file.path(MESO, ct, "no_bias_model", "finemo_out_merged_counts", "hits.tsv")
  d <- fread(f, header = TRUE)
  d <- d[grepl(motif_regex, motif_name)]
  gr <- GRanges(d$chr, IRanges(d$start + 1L, d$end))
  seqlevelsStyle(gr) <- "UCSC"
  gr
}

# one analysis row: instance overlap + peak-level Fisher vs all peaks
analyse <- function(label, ct, hits_gr, peaks_gr, sox9_gr) {
  n_hits <- length(hits_gr)
  # instance-level
  inst_ov <- sum(overlapsAny(hits_gr, sox9_gr))
  # peak-level: which peaks are SOX9+, which host this motif
  peaks_sox9 <- overlapsAny(peaks_gr, sox9_gr)
  peaks_mot  <- overlapsAny(peaks_gr, hits_gr)
  a <- sum(peaks_mot  & peaks_sox9)   # motif+ & SOX9+
  b <- sum(peaks_mot  & !peaks_sox9)  # motif+ & SOX9-
  c <- sum(!peaks_mot & peaks_sox9)   # motif- & SOX9+
  d <- sum(!peaks_mot & !peaks_sox9)  # motif- & SOX9-
  ft <- fisher.test(matrix(c(a, b, c, d), 2), alternative = "two.sided")
  data.table(
    motif            = label,
    celltype         = ct,
    n_instances      = n_hits,
    inst_overlap_SOX9= inst_ov,
    inst_pct_SOX9    = round(100 * inst_ov / max(n_hits, 1), 2),
    n_motif_peaks    = a + b,
    motifpeaks_SOX9pos = a,
    motifpeaks_pct_SOX9 = round(100 * a / max(a + b, 1), 2),
    baseline_pct_SOX9 = round(100 * sum(peaks_sox9) / length(peaks_gr), 2),
    odds_ratio       = round(unname(ft$estimate), 3),
    fisher_p         = signif(ft$p.value, 3)
  )
}

## ---- 2. run per cell type -------------------------------------------------
MOTIFS <- list(
  "Average_65__mp3 (C/EBP-palindrome / putative)" = "Average_65__merged_pattern_3$|Average_65__merged_pattern_3\\b",
  "Average_65__mp0 (C/EBP, control)"              = "Average_65__merged_pattern_0\\b",
  "ALL merged_counts motifs (baseline)"           = "."
)

res <- list()
for (k in names(CTS)) {
  ct <- CTS[[k]]
  message("\n=== ", ct, " ===")
  peaks <- read_peaks(ct)
  for (lab in names(MOTIFS)) {
    rgx  <- MOTIFS[[lab]]
    hits <- if (rgx == ".") {
      f <- file.path(MESO, ct, "no_bias_model", "finemo_out_merged_counts", "hits.tsv")
      d <- fread(f, header = TRUE)
      gr <- GRanges(d$chr, IRanges(d$start + 1L, d$end)); seqlevelsStyle(gr) <- "UCSC"; gr
    } else read_hits(ct, rgx)
    res[[paste(ct, lab)]] <- analyse(lab, ct, hits, peaks, sox9)
    message(sprintf("  %-46s n=%d  inst%%SOX9=%.1f  peak%%SOX9=%.1f (base %.1f) OR=%.2f p=%.1e",
                    lab, tail(res, 1)[[1]]$n_instances,
                    tail(res, 1)[[1]]$inst_pct_SOX9,
                    tail(res, 1)[[1]]$motifpeaks_pct_SOX9,
                    tail(res, 1)[[1]]$baseline_pct_SOX9,
                    tail(res, 1)[[1]]$odds_ratio,
                    tail(res, 1)[[1]]$fisher_p))
  }
}

out <- rbindlist(res)
fwrite(out, file.path(OUT, "mp3_SOX9_overlap_results.csv"))
message("\nWrote ", file.path(OUT, "mp3_SOX9_overlap_results.csv"))
print(out)
