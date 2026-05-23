# run_stitching_ct_enrichment.R
#
# Adds a "Stitching_CT" condition to myeloid_SE_hub_enrichment.csv:
#   For each cell type, stitch that cell type's own PeakCalls peakset
#   (gap <= 12,500 bp, >= 5 peaks) and run hypergeometric SE enrichment
#   using the full CT peakset as background.
#
# Unlike the existing "Stitching" condition (union peakset + Wilcoxon DA),
# hubs here are cell-type-specific by construction — no DA step needed.
#
# Appends results to myeloid_SE_hub_enrichment.csv and regenerates plots.
# Runs on login node (~2 min, no fragment loading required).

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
})
addArchRGenome("hg38")
addArchRThreads(threads = 1)

script_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q7"
main_dir    <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR"
peakcalls   <- file.path(main_dir, "PeakCalls")
se_rds      <- file.path(main_dir, "SE_regions_SE_database.rds")
se_path2    <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/SE_db2_from_Wooseung"
chain_file  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/liftover/hg38ToHg19.over.chain"
csv_out     <- file.path(script_dir, "myeloid_SE_hub_enrichment.csv")

MAX_DIST          <- 12500L
MIN_PEAKS_VALUES  <- c(3L, 4L, 5L)
MIN_DIST_TSS      <- 0L
REMOVE_CHR        <- c("chrX", "chrY")
FDR_CUT           <- 0.05

# ---- load project, SE databases, chain --------------------------------------
message("Loading ArchR project...")
proj         <- loadArchRProject(main_dir, showLogo = FALSE)
ct_levels    <- sort(unique(as.character(proj@cellColData[, "celltype_lv1"])))
all_peaks    <- getPeakSet(proj)
peaks_distal <- all_peaks[
  !as.character(seqnames(all_peaks)) %in% REMOVE_CHR &
  all_peaks$distToGeneStart >= MIN_DIST_TSS
]
message(sprintf("  %d distal peaks (TSS >= %d bp, no chrX)",
                length(peaks_distal), MIN_DIST_TSS))

message("Loading SE databases (hg19)...")
SE_db1 <- readRDS(se_rds)
SE_db2 <- lapply(list.files(se_path2, pattern = "\\.bed$", full.names = TRUE), function(f) {
  y <- read.table(f, header = TRUE, sep = "\t")[, 1:3]
  colnames(y) <- c("chr", "start", "end")
  makeGRangesFromDataFrame(y)
})
names(SE_db2) <- sub("\\.bed$", "", basename(list.files(se_path2, pattern = "\\.bed$")))
SE_all <- c(SE_db1, SE_db2)
message(sprintf("  %d SE sets", length(SE_all)))

ch <- import.chain(chain_file)

.liftover <- function(gr) {
  seqlevelsStyle(gr) <- "UCSC"
  unlist(liftOver(gr, ch))
}

.phyper_enrich <- function(hub_h19, bg_h19, se_gr) {
  q <- sum(overlapsAny(hub_h19, se_gr, ignore.strand = TRUE))
  m <- sum(overlapsAny(bg_h19,  se_gr, ignore.strand = TRUE))
  n <- length(bg_h19) - m
  k <- length(hub_h19)
  list(q = q, m = m, n = n, k = k,
       frac_hub = q / max(k, 1),
       frac_bg  = m / max(length(bg_h19), 1),
       pval     = phyper(q = q, m = m, n = n, k = k, lower.tail = FALSE))
}

# ---- per-celltype stitching enrichment (sweep over MIN_PEAKS_VALUES) ---------
all_rows <- list()

for (mp in MIN_PEAKS_VALUES) {
  cond_label <- sprintf("Stitching_CT_mp%d", mp)
  message(sprintf("\nRunning %s enrichment per cell type...", cond_label))

  for (ct in ct_levels) {
    pk_file <- file.path(peakcalls, paste0(ct, "-reproduciblePeaks.gr.rds"))
    if (!file.exists(pk_file)) {
      message(sprintf("  %s: PeakCalls file not found, skipping", ct))
      next
    }

    peaks_ct <- readRDS(pk_file)
    if ("distToGeneStart" %in% colnames(mcols(peaks_ct))) {
      peaks_ct <- peaks_ct[peaks_ct$distToGeneStart >= MIN_DIST_TSS]
    } else {
      peaks_ct <- peaks_ct[overlapsAny(peaks_ct, peaks_distal, ignore.strand = TRUE)]
    }
    peaks_ct <- peaks_ct[!as.character(seqnames(peaks_ct)) %in% REMOVE_CHR]

    stitched <- reduce(peaks_ct, min.gapwidth = MAX_DIST, ignore.strand = TRUE)
    hits     <- findOverlaps(peaks_ct, stitched)
    npeak    <- tabulate(subjectHits(hits), nbins = length(stitched))
    stitched <- stitched[npeak >= mp]

    if (length(stitched) == 0) {
      message(sprintf("  %s / %s: 0 stitching hubs, skipping", cond_label, ct))
      next
    }

    hubs_h19 <- .liftover(stitched)
    bg_h19   <- .liftover(peaks_ct)

    message(sprintf("  %s / %s: %d peaks -> %d hubs -> %d after liftover",
                    cond_label, ct, length(peaks_ct), length(stitched), length(hubs_h19)))

    if (length(hubs_h19) == 0) next

    rows <- lapply(names(SE_all), function(se_name) {
      r <- .phyper_enrich(hubs_h19, bg_h19, SE_all[[se_name]])
      data.frame(condition   = cond_label,
                 celltype    = ct,
                 se_db       = se_name,
                 n_da_hubs   = r$k,
                 q_overlap   = r$q,
                 frac_hub_SE = r$frac_hub,
                 frac_bg_SE  = r$frac_bg,
                 pval        = r$pval,
                 stringsAsFactors = FALSE)
    })
    all_rows <- c(all_rows, rows)
  }
}

new_results <- do.call(rbind, all_rows)
new_results <- new_results %>%
  group_by(condition, celltype) %>%
  mutate(fdr           = p.adjust(pval, "BH"),
         neg_log10_fdr = -log10(pmax(fdr, 1e-300)),
         enriched      = fdr < FDR_CUT) %>%
  ungroup()

# ---- append to existing CSV (idempotent) ------------------------------------
existing  <- read.csv(csv_out, stringsAsFactors = FALSE)
ct_labels <- unique(new_results$condition)
existing  <- existing[!existing$condition %in% ct_labels, ]
combined  <- bind_rows(existing, new_results)
write.csv(combined, csv_out, row.names = FALSE, quote = FALSE)
message(sprintf("\nAppended %d rows (%s) to %s",
                nrow(new_results), paste(ct_labels, collapse=", "), basename(csv_out)))

# summary
sig <- new_results %>% filter(enriched) %>% count(condition, celltype, name = "n_sig")
message("\nSignificant SE sets per condition x cell type:")
print(as.data.frame(sig[order(sig$condition, -sig$n_sig), ]), row.names = FALSE)

message("\nDone. Re-run plot_SE_hub_comparison.R to update figures.")
