#!/usr/bin/env Rscript
# R1_Q3_chromVAR_motif_vs_all_HDMA.R   <motif_query>  <tag>
#
# Generalized version of the SOX9 analysis: which HDMA-merged ChromBPNet
# patterns (all merged_counts + merged_profile, called by FiNeMo) significantly
# overlap the chromVAR motif hit regions of a chosen TF (e.g. RUNX2) from the
# ArchR tumor_compartment (SOX9 P23) project.
#
#   motif_query : regex matched against names() of the ArchR Motif positions
#                 (e.g. "RUNX2", "SOX9")
#   tag         : short label used for the output CSV
#
# Per ChromBPNet cell type (SOX9_high_P23 / SOX9_low_P23):
#   instance-level : fraction of FiNeMo instances within +/-FLANK bp of a TF motif
#   peak-level     : Fisher enrichment of (pattern+ peak x TF+ peak), BH-adjusted

suppressMessages({ library(GenomicRanges); library(data.table) })

args  <- commandArgs(trailingOnly = TRUE)
QUERY <- if (length(args) >= 1) args[[1]] else "SOX9"
TAG   <- if (length(args) >= 2) args[[2]] else QUERY

P     <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
MESO  <- "/sc/arion/scratch/giottb01/meso_SOX9_P23_chromBPnet"
OUT   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
FLANK <- 50L
CTS   <- c("SOX9_high_P23", "SOX9_low_P23")

## ---- TF chromVAR motif positions -----------------------------------------
message("Loading ArchR Motif-Positions-In-Peaks.rds ...")
pos <- readRDS(file.path(P, "Annotations", "Motif-Positions-In-Peaks.rds"))
tfn <- grep(QUERY, names(pos), value = TRUE, ignore.case = TRUE)
message("Motif entries matching '", QUERY, "': ", paste(tfn, collapse = ", "))
if (length(tfn) == 0) stop("No motif in annotation matches query '", QUERY, "'")
tf <- unlist(GRangesList(pos[tfn]), use.names = FALSE)
tf <- reduce(resize(tf, width(tf) + 2L * FLANK, fix = "center"))
seqlevelsStyle(tf) <- "UCSC"
message(TAG, " regions (+/-", FLANK, "bp, reduced): ", length(tf))

## ---- TOMTOM annotation ----------------------------------------------------
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
  peak_tf <- data.table(peak_name = peaks$peak_name, tf = overlapsAny(peaks, tf))
  n_peaks <- nrow(peak_tf); n_tfpk <- sum(peak_tf$tf)
  message(sprintf("  peaks=%d  %s+ peaks=%d (%.2f%%)", n_peaks, TAG, n_tfpk,
                  100 * n_tfpk / n_peaks))

  hits <- rbindlist(lapply(c("counts", "profile"), function(ty) {
    h <- fread(file.path(MESO, ct, "no_bias_model", paste0("finemo_out_merged_", ty),
                         "hits.tsv"),
               select = c("chr", "start", "end", "motif_name", "peak_name"))
    h[, type := ty]; h
  }))
  hg <- GRanges(hits$chr, IRanges(hits$start + 1L, hits$end)); seqlevelsStyle(hg) <- "UCSC"
  hits[, in_tf := overlapsAny(hg, tf)]
  hits[, tf_peak := peak_tf$tf[match(peak_name, peak_tf$peak_name)]]

  agg <- hits[, .(
      n_instances   = .N,
      inst_in_tf    = sum(in_tf),
      n_motif_peaks = uniqueN(peak_name),
      motif_tf_pk   = uniqueN(peak_name[tf_peak])
    ), by = .(type, motif_name)]

  fish <- function(a, m) {
    b <- m - a; c <- n_tfpk - a; d <- n_peaks - m - c
    ft <- fisher.test(matrix(c(a, b, c, d), 2)); c(OR = unname(ft$estimate), p = ft$p.value)
  }
  ff <- t(mapply(fish, agg$motif_tf_pk, agg$n_motif_peaks))
  agg[, `:=`(
      inst_pct_tf  = round(100 * inst_in_tf / n_instances, 3),
      peak_pct_tf  = round(100 * motif_tf_pk / n_motif_peaks, 2),
      baseline_pct = round(100 * n_tfpk / n_peaks, 2),
      odds_ratio   = round(ff[, "OR"], 3),
      fisher_p     = ff[, "p"])]
  agg[, padj := p.adjust(fisher_p, "BH")]
  agg[, celltype := ct]
  agg[, pattern := strip(motif_name)]
  agg <- merge(agg, ann, by = c("type", "pattern"), all.x = TRUE, sort = FALSE)
  all_res[[ct]] <- agg
}

res <- rbindlist(all_res, fill = TRUE)
setorder(res, celltype, -odds_ratio)
outf <- file.path(OUT, paste0(TAG, "_vs_all_HDMA_patterns.csv"))
fwrite(res, outf)
message("\nWrote ", outf)

sig <- res[padj < 0.05 & odds_ratio > 1]
message("\n", TAG, ": patterns significantly ENRICHED (padj<0.05, OR>1): ",
        nrow(sig), " of ", nrow(res), " pattern x celltype tests")
message(TAG, ": instance-level overlap distribution per celltype:")
print(res[, .(median = median(inst_pct_tf), mean = round(mean(inst_pct_tf), 2),
              max = max(inst_pct_tf)), by = celltype])
