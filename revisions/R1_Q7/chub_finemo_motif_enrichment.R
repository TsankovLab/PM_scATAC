# chub_finemo_motif_enrichment.R
#
# 1. Load (or build) co-accessibility hubs at r=0.3, min_peaks=3
#    Hub GRanges cached as hub_gr_CoAccess_r0.3_mp3.rds
#    Hub x cell CPM matrix cached as hub_cell_mat_CoAccess_r0.3_mp3.rds
# 2. Wilcoxon DA (presto) per celltype_lv2 -> DA hub GRanges per CT
# 3. Hypergeometric enrichment: FiNeMo motifs over-represented in DA hub peaks
#
# Output: chub_finemo_motif_enrichment.csv
# Submit: bsub < submit_chub_finemo_motif_enrichment.sh

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(Matrix)
  library(igraph)
  library(presto)
  library(dplyr)
  library(data.table)
})
addArchRGenome("hg38")
addArchRThreads(threads = 1)

source("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils/Hubs_finder.R")

# ---- paths -------------------------------------------------------------------
script_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q7"
main_dir    <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR"
finemo_base <- "/sc/arion/scratch/giottb01/Xmen/meso"
mat_cache   <- file.path(script_dir, "hub_cell_matrices_v2")

# ---- parameters --------------------------------------------------------------
COR_CUTOFF   <- 0.3
MIN_PEAKS    <- 3L
MAX_DIST     <- 12500L
MIN_DIST_TSS <- 0L
REMOVE_CHR   <- c("chrX", "chrY")
META_GROUP   <- "celltype_lv2"
FDR_CUT      <- 0.05
LFC_CUT      <- 1
FINEMO_TYPE  <- "finemo_out_counts"

COND_LBL   <- sprintf("CoAccess_r%.1f_mp%d", COR_CUTOFF, MIN_PEAKS)
GR_PATH    <- file.path(mat_cache, paste0("hub_gr_",       COND_LBL, ".rds"))
MAT_PATH   <- file.path(mat_cache, paste0("hub_cell_mat_", COND_LBL, ".rds"))

# ---- load ArchR project ------------------------------------------------------
message("Loading ArchR project...")
proj      <- loadArchRProject(main_dir, showLogo = FALSE)
cells     <- rownames(proj@cellColData)
celltypes <- as.character(proj@cellColData[, META_GROUP])
ct_levels <- unique(celltypes[!is.na(celltypes)])
message(sprintf("  %d cells, %d cell types (lv2)", nCells(proj), length(ct_levels)))

# ---- load or build hub GRanges -----------------------------------------------
if (file.exists(GR_PATH)) {
  message("Loading cached hub GRanges: ", basename(GR_PATH))
  hub_gr <- readRDS(GR_PATH)
} else {
  message("Hub GRanges cache missing — running hubs_finder at r=", COR_CUTOFF, "...")
  all_peaks    <- getPeakSet(proj)
  peaks_distal <- all_peaks[
    !as.character(seqnames(all_peaks)) %in% REMOVE_CHR &
    all_peaks$distToGeneStart >= MIN_DIST_TSS
  ]
  cA <- as.data.frame(
    getCoAccessibility(proj, corCutOff = COR_CUTOFF, resolution = 1,
                       returnLoops = FALSE),
    row.names = NULL)
  message(sprintf("  %d co-accessibility pairs at r >= %.1f", nrow(cA), COR_CUTOFF))
  hub_info <- hubs_finder(proj,
                          coAccessibility = cA,
                          cor_cutoff      = COR_CUTOFF,
                          max_dist        = MAX_DIST,
                          min_dist_tss    = MIN_DIST_TSS,
                          min_peaks       = MIN_PEAKS,
                          remove_chr      = REMOVE_CHR)
  hub_gr <- hub_info$hubsCollapsed
  saveRDS(hub_gr, GR_PATH)
  message(sprintf("  %d hubs detected; saved to %s", length(hub_gr), basename(GR_PATH)))
}
hub_ids <- names(hub_gr)
message(sprintf("  %d hubs loaded", length(hub_gr)))

# ---- load or build hub x cell matrix -----------------------------------------
if (file.exists(MAT_PATH)) {
  message("Loading cached hub x cell matrix: ", basename(MAT_PATH))
  hub_mat <- readRDS(MAT_PATH)
} else {
  message("Matrix cache missing — loading fragments to build matrix...")
  fragments <- unlist(getFragmentsFromProject(ArchRProj = proj))
  frag_rg   <- as.character(fragments$RG)
  cell_idx  <- setNames(seq_along(cells), cells)
  n_frags   <- proj$nFrags

  ovl    <- findOverlaps(hub_gr, fragments, ignore.strand = TRUE)
  q_hits <- queryHits(ovl)
  s_hits <- subjectHits(ovl)
  j      <- cell_idx[frag_rg[s_hits]]
  keep   <- !is.na(j)
  hub_mat <- sparseMatrix(
    i = q_hits[keep], j = j[keep], x = 1,
    dims = c(length(hub_gr), length(cells)),
    dimnames = list(hub_ids, cells)
  )
  hub_mat <- t(t(hub_mat) * (1e6 / n_frags))
  saveRDS(hub_mat, MAT_PATH)
  message(sprintf("  Saved to %s", basename(MAT_PATH)))
  rm(fragments, frag_rg, cell_idx, ovl, q_hits, s_hits, j, keep)
  gc()
}

# ---- DA hubs per cell type ---------------------------------------------------
message("\nComputing DA hubs per celltype_lv2...")
da_res <- wilcoxauc(log2(hub_mat + 1), celltypes)
da_res <- da_res[da_res$padj < FDR_CUT & da_res$logFC > LFC_CUT, ]
da_cts <- unique(da_res$group)
message(sprintf("  %d DA hub-CT pairs across %d cell types",
                nrow(da_res), length(da_cts)))

# ---- load FiNeMo data for all cell types -------------------------------------
message("\nPre-loading FiNeMo data...")

.load_finemo <- function(ct) {
  bed_file <- file.path(finemo_base, ct, paste0(ct, "_peakset_all_no_blacklist.bed"))
  occ_file <- file.path(finemo_base, ct, "no_bias_model", FINEMO_TYPE,
                        "report", "motif_occurrences.tsv")
  qc_file  <- file.path(finemo_base, ct, "no_bias_model", FINEMO_TYPE, "peaks_qc.tsv")
  gb_file  <- file.path(finemo_base, ct, "no_bias_model",
                        paste0(ct, "_finemo_counts_to_genome_browser.tsv"))

  if (!all(file.exists(bed_file, occ_file, qc_file, gb_file))) {
    message("  Skipping ", ct, " — missing FiNeMo files")
    return(NULL)
  }

  bed <- fread(bed_file, header = FALSE, select = 1:4)
  setnames(bed, c("chr", "start", "end", "peak_name"))
  peak_gr <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE,
                                      seqnames.field = "chr",
                                      start.field    = "start",
                                      end.field      = "end")
  peak_gr$peak_name <- bed$peak_name

  occ       <- fread(occ_file, header = TRUE)
  setnames(occ, 1, "peak_id")
  qc        <- fread(qc_file, header = TRUE, select = c("peak_id", "peak_name"))
  occ_named <- merge(occ, qc, by = "peak_id", sort = FALSE)
  motif_cols <- setdiff(names(occ), "peak_id")
  mat        <- as.matrix(occ_named[, ..motif_cols]) > 0
  rownames(mat) <- occ_named$peak_name

  gb <- fread(gb_file, header = FALSE, select = c(4, 6))
  setnames(gb, c("tf_name", "motif_name"))
  tf_map <- unique(gb[, .(motif_name, tf_name)])

  message("  Loaded: ", ct, " (", nrow(mat), " peaks, ", ncol(mat), " motifs)")
  list(peak_gr = peak_gr, motif_mat = mat, tf_map = tf_map)
}

# Discover CTs with FiNeMo output and intersect with ArchR lv2 labels
finemo_dirs <- basename(list.dirs(finemo_base, recursive = FALSE))
finemo_dirs <- finemo_dirs[!grepl("^MACS2_|^fragments_", finemo_dirs)]
overlap_cts <- intersect(ct_levels, finemo_dirs)
message(sprintf("  %d cell types overlap between ArchR lv2 and FiNeMo: %s",
                length(overlap_cts), paste(sort(overlap_cts), collapse = ", ")))

finemo_data <- setNames(lapply(overlap_cts, .load_finemo), overlap_cts)
finemo_data <- Filter(Negate(is.null), finemo_data)
overlap_cts <- names(finemo_data)
message(sprintf("  %d cell types successfully loaded", length(overlap_cts)))

# ---- hypergeometric enrichment -----------------------------------------------
message("\nRunning FiNeMo motif enrichment per cell type...")

.phyper_motif <- function(in_hub, has_hit) {
  k <- sum(in_hub)
  m <- sum(has_hit)
  n <- length(has_hit) - m
  q <- sum(in_hub & has_hit)
  list(q = q, m = m, n = n, k = k,
       frac_hub = if (k > 0) q / k else 0,
       frac_bg  = m / length(has_hit),
       pval     = phyper(q = q, m = m, n = n, k = k, lower.tail = FALSE))
}

all_rows <- list()

for (ct in overlap_cts) {
  finemo <- finemo_data[[ct]]

  # DA hub IDs for this cell type
  da_ids <- da_res$feature[da_res$group == ct]
  if (length(da_ids) == 0) {
    message(sprintf("  %s: 0 DA hubs, skipping", ct))
    next
  }

  hub_idx   <- match(da_ids, hub_ids)
  hub_idx   <- hub_idx[!is.na(hub_idx)]
  da_hub_gr <- hub_gr[hub_idx]

  # FiNeMo peaks overlapping DA hub regions
  in_hub <- overlapsAny(finemo$peak_gr, da_hub_gr, ignore.strand = TRUE)
  n_hub  <- sum(in_hub)
  n_bg   <- length(in_hub)
  message(sprintf("  %s: %d DA hubs -> %d / %d FiNeMo peaks in hubs",
                  ct, length(da_hub_gr), n_hub, n_bg))
  if (n_hub == 0) next

  peak_names <- finemo$peak_gr$peak_name
  mat_sub    <- finemo$motif_mat[match(peak_names, rownames(finemo$motif_mat)), ,
                                 drop = FALSE]

  # per-CT tf_name lookup
  tf_map_ct <- setNames(finemo$tf_map$tf_name, finemo$tf_map$motif_name)

  rows <- lapply(colnames(mat_sub), function(m) {
    res <- .phyper_motif(in_hub, mat_sub[, m])
    data.frame(celltype    = ct,
               motif_name  = m,
               tf_name     = unname(tf_map_ct[m]),
               n_hub_peaks = res$k,
               n_bg_peaks  = n_bg,
               q_overlap   = res$q,
               m_hits      = res$m,
               frac_hub    = res$frac_hub,
               frac_bg     = res$frac_bg,
               pval        = res$pval,
               stringsAsFactors = FALSE)
  })
  all_rows <- c(all_rows, rows)
}

# ---- FDR correction ---------------------------------------------------------
message("\nFinalizing results...")

results <- bind_rows(all_rows)

# tf_name is already per-CT (added inside the loop above)
results <- results %>%
  group_by(celltype) %>%
  mutate(fdr           = p.adjust(pval, "BH"),
         neg_log10_fdr = -log10(pmax(fdr, 1e-300)),
         enriched      = fdr < FDR_CUT) %>%
  ungroup() %>%
  select(celltype, motif_name, tf_name, n_hub_peaks, n_bg_peaks,
         q_overlap, m_hits, frac_hub, frac_bg, pval, fdr, neg_log10_fdr, enriched)

out_csv <- file.path(script_dir, "chub_finemo_motif_enrichment.csv")
write.csv(results, out_csv, row.names = FALSE, quote = FALSE)
message("Wrote: ", out_csv)

sig_summary <- results %>%
  filter(enriched) %>%
  count(celltype, name = "n_sig_motifs")
message("\nSignificant motifs (FDR < 0.05) per cell type:")
print(as.data.frame(sig_summary[order(sig_summary$celltype), ]), row.names = FALSE)

message("\nDone.")
