# myeloid_SE_hub_enrichment.R
#
# Computes Super Enhancer enrichment in hub regions across:
#   1. Co-accessibility hubs (r = 0.2–0.6, min_peaks = 3/4/5, tss_dist = 0/2000)
#   2. Stitching hubs — consensus peak set reduced at 12.5 kb gap
#   3. Stitching_CT hubs — per-cell-type PeakCalls reduced at 12.5 kb gap
#
# Condition labels: {type}_mp{N}_tss{D}  e.g. CoAccess_r0.3_mp3_tss0
# For tss=0: falls back to legacy hub_cell_mat cache files (no _tss suffix).
#
# Parallel design (mclapply, fork-safe on Linux):
#   Phase 1a — hubs_finder: parallel over all conditions (N_CORES)
#              memory-light; fragments not touched here
#   Phase 1b — matrix building: limited parallelism (N_CORES_MAT=2)
#              memory-heavy (findOverlaps on fragments); only uncached conditions
#   Phase 2  — DA wilcoxauc + SE enrichment: parallel over conditions (N_CORES)
#              each worker reads its own matrix from cache; not held in parent
#   Phase 3  — Stitching_CT: parallel over (mp, tss, ct) triplets (N_CORES)
#
# Output: myeloid_SE_hub_enrichment.csv
# Submit: bsub < submit_myeloid_SE_hub_enrichment.sh

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(Matrix)
  library(rtracklayer)
  library(igraph)
  library(presto)
  library(dplyr)
  library(parallel)
})
addArchRGenome("hg38")
addArchRThreads(threads = 1)

source("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils/Hubs_finder.R")

# ---- paths ------------------------------------------------------------------
script_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q7"
main_dir    <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR"
peakcalls   <- file.path(main_dir, "PeakCalls")
se_rds      <- file.path(main_dir, "SE_regions_SE_database.rds")
se_path2    <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/SE_db2_from_Wooseung"
chain_file  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/liftover/hg38ToHg19.over.chain"
mat_cache   <- file.path(script_dir, "hub_cell_matrices_v2")
dir.create(mat_cache, showWarnings = FALSE)

# ---- parameters -------------------------------------------------------------
COR_CUTOFFS          <- c(0.2, 0.3, 0.4, 0.5, 0.6)
MAX_DIST             <- 12500L
MIN_PEAKS_VALUES     <- c(3L, 4L, 5L)
MIN_DIST_TSS_VALUES  <- c(0L, 2000L)
REMOVE_CHR           <- c("chrX", "chrY")
META_GROUP           <- "celltype_lv1"
FDR_CUT              <- 0.05
LFC_CUT              <- 1

N_CORES     <- min(as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "4")), 8L)
N_CORES_MAT <- 2L  # matrix building hits fragments; keep to 2 to avoid OOM
message(sprintf("N_CORES=%d  N_CORES_MAT=%d", N_CORES, N_CORES_MAT))

# ---- load main ArchR project ------------------------------------------------
message("Loading main ArchR project...")
proj <- loadArchRProject(main_dir, showLogo = FALSE)
cells     <- rownames(proj@cellColData)
celltypes <- as.character(proj@cellColData[, META_GROUP])
ct_levels <- unique(celltypes[!is.na(celltypes)])
message(sprintf("  %d cells, %d peaks, %d cell types",
                nCells(proj), length(getPeakSet(proj)), length(ct_levels)))

all_peaks <- getPeakSet(proj)

# ---- load SE databases (hg19) -----------------------------------------------
message("Loading SE databases...")
SE_db1 <- readRDS(se_rds)
SE_db2 <- lapply(list.files(se_path2, pattern = "\\.bed$", full.names = TRUE), function(f) {
  y <- read.table(f, header = TRUE, sep = "\t")[, 1:3]
  colnames(y) <- c("chr", "start", "end")
  makeGRangesFromDataFrame(y)
})
names(SE_db2) <- sub("\\.bed$", "", basename(list.files(se_path2, pattern = "\\.bed$")))
SE_all <- c(SE_db1, SE_db2)
message(sprintf("  %d SE sets (%d db1 + %d db2)", length(SE_all), length(SE_db1), length(SE_db2)))
ch <- import.chain(chain_file)

# ---- load all fragments once ------------------------------------------------
message("Loading all fragments (shared across all conditions)...")
fragments <- unlist(getFragmentsFromProject(ArchRProj = proj))
message(sprintf("  %s fragments loaded", format(length(fragments), big.mark = ",")))
frag_rg  <- as.character(fragments$RG)
cell_idx <- setNames(seq_along(cells), cells)
n_frags  <- proj$nFrags

# ---- helpers ----------------------------------------------------------------

# Safely combine a list that may contain NULLs or mclapply try-errors.
.safe_bind <- function(lst) {
  dfs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, lst)
  if (length(dfs) == 0) return(NULL)
  bind_rows(dfs)
}

# Builds hub x cell CPM matrix; writes to cache_path.
# Returns matrix (used by caller or discarded if only caching is needed).
.hub_cell_mat <- function(hub_gr, hub_ids, cache_path) {
  if (file.exists(cache_path)) {
    message("    cache hit: ", basename(cache_path))
    return(invisible(cache_path))  # signal cached; caller loads if needed
  }
  message(sprintf("    Computing %d hubs x %d cells matrix...",
                  length(hub_gr), length(cells)))
  ovl    <- findOverlaps(hub_gr, fragments, ignore.strand = TRUE)
  q_hits <- queryHits(ovl)
  s_hits <- subjectHits(ovl)
  j      <- cell_idx[frag_rg[s_hits]]
  keep   <- !is.na(j)
  mat    <- sparseMatrix(
    i = q_hits[keep], j = j[keep], x = 1,
    dims = c(length(hub_gr), length(cells)),
    dimnames = list(hub_ids, cells)
  )
  mat <- t(t(mat) * (1e6 / n_frags))
  saveRDS(mat, cache_path)
  message(sprintf("    Saved: %s", basename(cache_path)))
  invisible(cache_path)
}

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
       frac_hub = q / k, frac_bg = m / length(bg_h19),
       pval = phyper(q = q, m = m, n = n, k = k, lower.tail = FALSE))
}

.filter_bg <- function(gr, min_tss = 0L) {
  if ("distToGeneStart" %in% colnames(mcols(gr))) {
    gr <- gr[gr$distToGeneStart >= min_tss]
  } else {
    local_distal <- all_peaks[
      !as.character(seqnames(all_peaks)) %in% REMOVE_CHR &
      all_peaks$distToGeneStart >= min_tss
    ]
    gr <- gr[overlapsAny(gr, local_distal, ignore.strand = TRUE)]
  }
  gr[!as.character(seqnames(gr)) %in% REMOVE_CHR]
}

.stitching_hubs <- function(min_peaks, min_tss = 0L) {
  peaks <- all_peaks[
    !as.character(seqnames(all_peaks)) %in% REMOVE_CHR &
    all_peaks$distToGeneStart >= min_tss
  ]
  stitched <- reduce(peaks, min.gapwidth = MAX_DIST, ignore.strand = TRUE)
  hits     <- findOverlaps(peaks, stitched)
  npeak    <- tabulate(subjectHits(hits), nbins = length(stitched))
  stitched <- stitched[npeak >= min_peaks]
  hub_ids  <- paste0("ST_HUB", seq_along(stitched))
  message(sprintf("    -> %d stitching hubs (min_peaks=%d, tss>=%d bp)",
                  length(stitched), min_peaks, min_tss))
  list(hub_gr = stitched, hub_ids = hub_ids)
}

# For tss=0, fall back to legacy cache name (without _tss0) to skip recompute.
.resolve_cache_path <- function(lbl, tss_dist) {
  new_path <- file.path(mat_cache, paste0("hub_cell_mat_", lbl, ".rds"))
  if (!file.exists(new_path) && tss_dist == 0L) {
    old_lbl  <- sub("_tss0$", "", lbl)
    old_path <- file.path(mat_cache, paste0("hub_cell_mat_", old_lbl, ".rds"))
    if (file.exists(old_path)) {
      message("    [tss=0] reusing legacy cache: ", basename(old_path))
      return(old_path)
    }
  }
  new_path
}

# ---- pre-compute co-accessibility once --------------------------------------
message("\nPre-computing co-accessibility at r >= ", min(COR_CUTOFFS), "...")
cA_precomputed <- as.data.frame(
  getCoAccessibility(proj, corCutOff = min(COR_CUTOFFS), resolution = 1,
                     returnLoops = FALSE),
  row.names = NULL)
message(sprintf("  %d pairs fetched", nrow(cA_precomputed)))

# ---- build parameter grids --------------------------------------------------
coacc_grid <- do.call(rbind, lapply(MIN_DIST_TSS_VALUES, function(tss) {
  do.call(rbind, lapply(MIN_PEAKS_VALUES, function(mp) {
    data.frame(tss_dist = tss, mp = mp, r = COR_CUTOFFS, stringsAsFactors = FALSE)
  }))
}))

stitch_grid <- do.call(rbind, lapply(MIN_DIST_TSS_VALUES, function(tss) {
  data.frame(tss_dist = tss, mp = MIN_PEAKS_VALUES, stringsAsFactors = FALSE)
}))

# ---- Phase 1a: parallel hubs_finder (memory-light; no matrix ops) -----------
message(sprintf("\nPhase 1a: hubs_finder across %d CoAccess + %d Stitching conditions (%d cores)...",
                nrow(coacc_grid), nrow(stitch_grid), N_CORES))

coacc_hub_list <- mclapply(seq_len(nrow(coacc_grid)), function(i) {
  r <- coacc_grid$r[i]; mp <- coacc_grid$mp[i]; tss_dist <- coacc_grid$tss_dist[i]
  lbl      <- sprintf("CoAccess_r%.1f_mp%d_tss%d", r, mp, tss_dist)
  hub_info <- hubs_finder(proj,
                           coAccessibility = cA_precomputed,
                           cor_cutoff      = r,
                           max_dist        = MAX_DIST,
                           min_dist_tss    = tss_dist,
                           min_peaks       = mp,
                           remove_chr      = REMOVE_CHR)
  if (is.null(hub_info)) return(NULL)
  hub_gr  <- hub_info$hubsCollapsed
  hub_ids <- names(hub_gr)
  list(lbl        = lbl,
       tss_dist   = tss_dist,
       hub_gr     = hub_gr,
       hub_ids    = hub_ids,
       cache_path = .resolve_cache_path(lbl, tss_dist))
}, mc.cores = N_CORES, mc.preschedule = FALSE)

stitch_hub_list <- mclapply(seq_len(nrow(stitch_grid)), function(i) {
  mp <- stitch_grid$mp[i]; tss_dist <- stitch_grid$tss_dist[i]
  lbl     <- sprintf("Stitching_mp%d_tss%d", mp, tss_dist)
  st_info <- .stitching_hubs(min_peaks = mp, min_tss = tss_dist)
  list(lbl        = lbl,
       tss_dist   = tss_dist,
       hub_gr     = st_info$hub_gr,
       hub_ids    = st_info$hub_ids,
       cache_path = .resolve_cache_path(lbl, tss_dist))
}, mc.cores = N_CORES, mc.preschedule = FALSE)

all_hub_list <- Filter(Negate(is.null),
                       c(coacc_hub_list, stitch_hub_list))
message(sprintf("  -> %d hub sets built", length(all_hub_list)))

# ---- Phase 1b: matrix building — limited cores to avoid OOM ----------------
# Each findOverlaps(hub_gr, fragments) uses ~10-15 GB; base memory ~70 GB.
# N_CORES_MAT=2 keeps peak usage well within 128 GB.
uncached <- Filter(function(x) !file.exists(x$cache_path), all_hub_list)
message(sprintf("\nPhase 1b: %d cached, %d to compute (N_CORES_MAT=%d)...",
                length(all_hub_list) - length(uncached),
                length(uncached), N_CORES_MAT))

if (length(uncached) > 0) {
  mclapply(uncached, function(entry) {
    .hub_cell_mat(entry$hub_gr, entry$hub_ids, entry$cache_path)
  }, mc.cores = N_CORES_MAT, mc.preschedule = FALSE)
}

# Build conditions index: hub_gr + hub_ids + cache_path (matrix NOT held in parent)
conditions <- setNames(
  lapply(all_hub_list, function(x) x[c("hub_gr", "hub_ids", "cache_path")]),
  sapply(all_hub_list, `[[`, "lbl")
)
message(sprintf("  -> %d conditions ready", length(conditions)))

# ---- Phase 2: DA + SE enrichment (parallel over conditions) -----------------
# Each worker loads its own matrix from cache → parent holds no matrices.
# wilcoxauc called once per condition (all CTs in one shot, not per-CT).
message(sprintf("\nPhase 2: DA + SE enrichment across %d conditions (%d cores)...",
                length(conditions), N_CORES))

valid_cts <- ct_levels[file.exists(
  file.path(peakcalls, paste0(ct_levels, "-reproduciblePeaks.gr.rds")))]

enrich_results_list <- mclapply(names(conditions), function(cond) {
  entry    <- conditions[[cond]]
  tss_cond <- as.integer(sub(".*_tss(\\d+)$", "\\1", cond))

  info_mat <- readRDS(entry$cache_path)

  da_all <- wilcoxauc(log2(info_mat + 1), celltypes)
  da_all <- da_all[da_all$padj < FDR_CUT & da_all$logFC > LFC_CUT, ]
  rm(info_mat); gc()

  cond_rows <- list()
  for (ct in valid_cts) {
    bg_h19 <- .liftover(.filter_bg(
      readRDS(file.path(peakcalls, paste0(ct, "-reproduciblePeaks.gr.rds"))),
      min_tss = tss_cond))

    da_ids  <- da_all$feature[da_all$group == ct]
    if (length(da_ids) == 0) next

    hub_idx  <- match(da_ids, entry$hub_ids)
    hub_idx  <- hub_idx[!is.na(hub_idx)]
    hubs_gr  <- entry$hub_gr[hub_idx]
    hubs_h19 <- .liftover(hubs_gr)
    if (length(hubs_h19) == 0) next

    message(sprintf("    %s / %s: %d DA hubs -> %d after liftover",
                    cond, ct, length(hubs_gr), length(hubs_h19)))

    rows <- lapply(names(SE_all), function(se_name) {
      r <- .phyper_enrich(hubs_h19, bg_h19, SE_all[[se_name]])
      data.frame(condition = cond, celltype = ct, se_db = se_name,
                 n_da_hubs = r$k, q_overlap = r$q,
                 frac_hub_SE = r$frac_hub, frac_bg_SE = r$frac_bg,
                 pval = r$pval, stringsAsFactors = FALSE)
    })
    cond_rows <- c(cond_rows, rows)
  }
  cond_rows
}, mc.cores = N_CORES, mc.preschedule = FALSE)

all_results <- unlist(enrich_results_list, recursive = FALSE)
da_results  <- .safe_bind(all_results)
if (!is.null(da_results)) {
  da_results <- da_results %>%
    mutate(n_da_hubs = as.integer(n_da_hubs),
           q_overlap = as.integer(q_overlap)) %>%
    group_by(condition, celltype) %>%
    mutate(fdr           = p.adjust(pval, "BH"),
           neg_log10_fdr = -log10(pmax(fdr, 1e-300)),
           enriched      = fdr < FDR_CUT) %>%
    ungroup()
}

# ---- Phase 3: Stitching_CT (parallel over mp × tss × ct) -------------------
message(sprintf("\nPhase 3: Stitching_CT enrichment (%d cores)...", N_CORES))

ct_grid <- do.call(rbind, lapply(MIN_DIST_TSS_VALUES, function(tss) {
  do.call(rbind, lapply(MIN_PEAKS_VALUES, function(mp) {
    data.frame(tss_dist = tss, mp = mp, ct = valid_cts, stringsAsFactors = FALSE)
  }))
}))

ct_results_list <- mclapply(seq_len(nrow(ct_grid)), function(i) {
  mp       <- ct_grid$mp[i]
  tss_dist <- ct_grid$tss_dist[i]
  ct       <- ct_grid$ct[i]
  cond_label <- sprintf("Stitching_CT_mp%d_tss%d", mp, tss_dist)

  peaks_ct <- .filter_bg(
    readRDS(file.path(peakcalls, paste0(ct, "-reproduciblePeaks.gr.rds"))),
    min_tss = tss_dist)
  stitched <- reduce(peaks_ct, min.gapwidth = MAX_DIST, ignore.strand = TRUE)
  hits     <- findOverlaps(peaks_ct, stitched)
  npeak    <- tabulate(subjectHits(hits), nbins = length(stitched))
  stitched <- stitched[npeak >= mp]
  if (length(stitched) == 0) return(NULL)

  hubs_h19 <- .liftover(stitched)
  bg_h19   <- .liftover(peaks_ct)
  if (length(hubs_h19) == 0) return(NULL)

  message(sprintf("    %s / %s: %d peaks -> %d hubs -> %d after liftover",
                  cond_label, ct, length(peaks_ct), length(stitched), length(hubs_h19)))

  lapply(names(SE_all), function(se_name) {
    r <- .phyper_enrich(hubs_h19, bg_h19, SE_all[[se_name]])
    data.frame(condition = cond_label, celltype = ct, se_db = se_name,
               n_da_hubs = r$k, q_overlap = r$q,
               frac_hub_SE = r$frac_hub, frac_bg_SE = r$frac_bg,
               pval = r$pval, stringsAsFactors = FALSE)
  })
}, mc.cores = N_CORES, mc.preschedule = FALSE)

ct_rows    <- unlist(ct_results_list, recursive = FALSE)
ct_results <- .safe_bind(ct_rows)
if (!is.null(ct_results)) {
  ct_results <- ct_results %>%
    mutate(n_da_hubs = as.integer(n_da_hubs),
           q_overlap = as.integer(q_overlap)) %>%
    group_by(condition, celltype) %>%
    mutate(fdr           = p.adjust(pval, "BH"),
           neg_log10_fdr = -log10(pmax(fdr, 1e-300)),
           enriched      = fdr < FDR_CUT) %>%
    ungroup()
}

# ---- write combined results --------------------------------------------------
results <- bind_rows(da_results, ct_results)
write.csv(results, file.path(script_dir, "myeloid_SE_hub_enrichment.csv"),
          row.names = FALSE, quote = FALSE)

sig_summary <- results %>%
  filter(enriched) %>%
  count(condition, celltype, name = "n_sig_SE")
message("\nSignificant SE sets (FDR < 0.05) per condition x celltype:")
print(as.data.frame(sig_summary[order(sig_summary$celltype, sig_summary$condition), ]),
      row.names = FALSE)

message("\nDone. Run plot_SE_hub_comparison.R to generate the figure.")
