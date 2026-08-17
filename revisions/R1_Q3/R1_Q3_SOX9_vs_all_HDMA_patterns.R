#!/usr/bin/env Rscript
# R1_Q3_SOX9_vs_all_HDMA_patterns.R
#
# Systematic test: which HDMA-merged ChromBPNet patterns (all merged_counts +
# merged_profile modisco motifs, called by FiNeMo) significantly overlap the
# SOX9 chromVAR motif (SOX9_756) hit regions from the ArchR tumor_compartment
# (SOX9 P23) project?
#
# For each pattern, per ChromBPNet cell type (SOX9_high_P23 / SOX9_low_P23):
#   * instance-level : fraction of FiNeMo instances within +/-FLANK bp of a SOX9 motif
#   * peak-level     : Fisher's exact enrichment of (pattern+ peak x SOX9+ peak)
#                      across all model peaks; OR + BH-adjusted p
# SOX-annotated patterns act as an internal positive control.

suppressMessages({ library(GenomicRanges); library(data.table) })

P     <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
MESO  <- "/sc/arion/scratch/giottb01/meso_SOX9_P23_chromBPnet"
OUT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
FLANK <- 50L
CTS   <- c("SOX9_high_P23", "SOX9_low_P23")

## ---- SOX9 chromVAR motif positions ---------------------------------------
message("Loading ArchR Motif-Positions-In-Peaks.rds ...")
pos  <- readRDS(file.path(P, "Annotations", "Motif-Positions-In-Peaks.rds"))
snm  <- grep("SOX9", names(pos), value = TRUE, ignore.case = TRUE)
message("SOX9 motif entries: ", paste(snm, collapse = ", "))
sox9 <- unlist(GRangesList(pos[snm]), use.names = FALSE)
sox9 <- reduce(resize(sox9, width(sox9) + 2L * FLANK, fix = "center"))
seqlevelsStyle(sox9) <- "UCSC"
message("SOX9 regions (+/-", FLANK, "bp, reduced): ", length(sox9))

## ---- TOMTOM annotation (TF label per pattern) ----------------------------
ann <- rbindlist(lapply(c("counts", "profile"), function(ty) {
  d <- fread(file.path(MESO, paste0("modisco_merged_", ty), "compiled_tomtom",
                       "modisco_compiled.tsv"))
  data.table(type = ty, pattern = d$pattern, TF = d$match0, qval = d$qval0)
}))
strip <- function(x) sub("^(pos|neg)_patterns\\.", "", x)

## ---- per cell type --------------------------------------------------------
all_res <- list()
for (ct in CTS) {
  message("\n=== ", ct, " ===")
  bed <- fread(file.path(MESO, ct, paste0(ct, "_peakset_all_no_blacklist.bed")),
               header = FALSE, select = 1:4,
               col.names = c("chr", "start", "end", "peak_name"))
  peaks <- GRanges(bed$chr, IRanges(bed$start + 1L, bed$end), peak_name = bed$peak_name)
  seqlevelsStyle(peaks) <- "UCSC"
  peak_sox9 <- data.table(peak_name = peaks$peak_name,
                          sox9 = overlapsAny(peaks, sox9))
  n_peaks  <- nrow(peak_sox9)
  n_sox9pk <- sum(peak_sox9$sox9)
  message(sprintf("  peaks=%d  SOX9+ peaks=%d (%.2f%%)",
                  n_peaks, n_sox9pk, 100 * n_sox9pk / n_peaks))

  hits <- rbindlist(lapply(c("counts", "profile"), function(ty) {
    h <- fread(file.path(MESO, ct, "no_bias_model", paste0("finemo_out_merged_", ty),
                         "hits.tsv"),
               select = c("chr", "start", "end", "motif_name", "peak_name"))
    h[, type := ty]; h
  }))
  # instance-level SOX9 overlap for every hit (one overlap call)
  hg <- GRanges(hits$chr, IRanges(hits$start + 1L, hits$end)); seqlevelsStyle(hg) <- "UCSC"
  hits[, in_sox9 := overlapsAny(hg, sox9)]
  hits[, sox9_peak := peak_sox9$sox9[match(peak_name, peak_sox9$peak_name)]]

  # per-pattern aggregation
  agg <- hits[, .(
      n_instances   = .N,
      inst_in_sox9  = sum(in_sox9),
      n_motif_peaks = uniqueN(peak_name),
      motif_sox9_pk = uniqueN(peak_name[sox9_peak])
    ), by = .(type, motif_name)]

  # Fisher per pattern (pattern+ x SOX9+ over all peaks)
  fish <- function(a, m) {           # a = motif&SOX9 peaks, m = motif peaks
    b <- m - a; c <- n_sox9pk - a; d <- n_peaks - m - c
    ft <- fisher.test(matrix(c(a, b, c, d), 2))
    c(OR = unname(ft$estimate), p = ft$p.value)
  }
  ff <- t(mapply(fish, agg$motif_sox9_pk, agg$n_motif_peaks))
  agg[, `:=`(
      inst_pct_sox9   = round(100 * inst_in_sox9 / n_instances, 3),
      peak_pct_sox9   = round(100 * motif_sox9_pk / n_motif_peaks, 2),
      baseline_pct    = round(100 * n_sox9pk / n_peaks, 2),
      odds_ratio      = round(ff[, "OR"], 3),
      fisher_p        = ff[, "p"]
  )]
  agg[, padj := p.adjust(fisher_p, "BH")]
  agg[, celltype := ct]
  agg[, pattern := strip(motif_name)]
  agg <- merge(agg, ann, by = c("type", "pattern"), all.x = TRUE, sort = FALSE)
  all_res[[ct]] <- agg
}

res <- rbindlist(all_res, fill = TRUE)
setorder(res, celltype, -odds_ratio)
fwrite(res, file.path(OUT, "SOX9_vs_all_HDMA_patterns.csv"))
message("\nWrote ", file.path(OUT, "SOX9_vs_all_HDMA_patterns.csv"))

## ---- concise console summary ---------------------------------------------
sig <- res[padj < 0.05 & odds_ratio > 1]
message("\nPatterns significantly ENRICHED for SOX9 overlap (padj<0.05, OR>1): ",
        nrow(sig), " of ", nrow(res), " pattern x celltype tests")

show <- function(ct) {
  message("\n--- ", ct, ": top 15 by odds ratio ---")
  x <- res[celltype == ct][order(-odds_ratio)][1:15,
        .(type, pattern, TF, n_instances, inst_pct_sox9, peak_pct_sox9,
          baseline_pct, odds_ratio, padj = signif(padj, 2))]
  print(x)
}
invisible(lapply(CTS, show))

message("\n--- SOX/SRY-annotated patterns (positive control) across both cell types ---")
print(res[grepl("SOX|SRY", TF), .(celltype, type, pattern, TF, inst_pct_sox9,
        peak_pct_sox9, baseline_pct, odds_ratio, padj = signif(padj, 2))][order(-odds_ratio)])
