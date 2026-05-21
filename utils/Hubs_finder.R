#' hubs_finder
#'
#' Detects co-accessibility hubs from an ArchR project.
#'
#' Logic follows the simplified sensitivity-analysis version (R1_Q7):
#'   1. Filter peaks to distal (>= min_dist_tss) non-blacklisted chromosomes
#'   2. Get co-accessibility pairs (or accept pre-computed to avoid re-fetching)
#'   3. Retain pairs on the same chromosome within max_dist of each other
#'   4. Build undirected graph; find connected components (igraph, C-level)
#'   5. Keep components with >= min_peaks peaks; collapse to hull GRanges
#'
#' Key differences from Hubs_finder_old.R:
#'   - No FDR or VarQuantile filters (those are strongly threshold-dependent)
#'   - Explicit max_dist spatial cap (comparable to stitching gap)
#'   - Global peak set rather than per-group PeakCalls
#'   - Accepts pre-computed coAccessibility to avoid re-reading Arrow files
#'     when sweeping cor_cutoff (major speedup for sensitivity analyses)
#'   - Vectorised hub collapsing (no per-hub lapply / peak-index lookup)
#'
#' @param ArchRProj      ArchRProject object
#' @param coAccessibility Pre-computed result of getCoAccessibility()
#'                        (data.frame or GInteractions). If NULL, computed
#'                        internally. Pass the result at the lowest desired
#'                        cor_cutoff and re-use for all higher thresholds.
#' @param cor_cutoff     Minimum co-accessibility correlation (default 0.4)
#' @param max_dist       Maximum distance (bp) between peak midpoints (default 12500)
#' @param min_dist_tss   Minimum distance from TSS to keep a peak (default 2000)
#' @param min_peaks      Minimum peaks per hub (default 5)
#' @param remove_chr     Chromosomes to exclude (default "chrX")
#' @param cores          Cores for parallel steps (default 1)
#'
#' @return A list with:
#'   $hubsCollapsed  GRanges of hub regions (min start – max end per component)
#'   $peakLinks      GRanges of all filtered peak-pair links (for arc plots)
#'   $n_hubs         Integer count of hubs
#'   $params         List of the parameters used
#'
#' @export
hubs_finder <- function(
  ArchRProj,
  coAccessibility = NULL,
  cor_cutoff      = 0.4,
  max_dist        = 12500L,
  min_dist_tss    = 2000L,
  min_peaks       = 5L,
  remove_chr      = "chrX",
  cores           = 1L
) {
  message(sprintf(
    "--* Hubs Finder *-- cor>=%.2f | max_dist=%d bp | min_dist_tss=%d bp | min_peaks=%d",
    cor_cutoff, max_dist, min_dist_tss, min_peaks
  ))

  # ---- 1. Peak set: distal + allowed chromosomes --------------------------------
  peaks     <- getPeakSet(ArchRProj)
  peaks$idx <- seq_along(peaks)
  keep_mask <- !as.character(seqnames(peaks)) %in% remove_chr &
               peaks$distToGeneStart >= min_dist_tss
  keep_idx  <- which(keep_mask)
  message(sprintf("  Peaks: %d / %d pass TSS / chr filters",
                  length(keep_idx), length(peaks)))

  # ---- 2. Co-accessibility pairs -----------------------------------------------
  # Accept pre-computed object to avoid re-reading Arrow files.
  # If sweeping multiple cor_cutoff values, call getCoAccessibility() once at
  # the lowest threshold, store the result, and pass it here for all others.
  if (is.null(coAccessibility)) {
    message("  Fetching co-accessibility from ArchR project...")
    cA <- as.data.frame(
      getCoAccessibility(ArchRProj,
                         corCutOff   = cor_cutoff,
                         resolution  = 1,
                         returnLoops = FALSE),
      row.names = NULL
    )
  } else {
    cA <- as.data.frame(coAccessibility, row.names = NULL)
    cA <- cA[cA$correlation >= cor_cutoff, , drop = FALSE]
  }
  message(sprintf("  Co-accessibility pairs at cor >= %.2f: %d", cor_cutoff, nrow(cA)))

  if (nrow(cA) == 0) {
    message("  No pairs above threshold — returning NULL")
    return(NULL)
  }

  # ---- 3. Spatial filters: distal peaks + same chr + within max_dist -----------
  cA <- cA[cA$queryHits   %in% keep_idx &
           cA$subjectHits %in% keep_idx, , drop = FALSE]

  peak_chr <- as.character(seqnames(peaks))
  peak_mid <- (start(peaks) + end(peaks)) / 2L

  same_chr <- peak_chr[cA$queryHits] == peak_chr[cA$subjectHits]
  in_range <- abs(peak_mid[cA$queryHits] - peak_mid[cA$subjectHits]) <= max_dist
  cA <- cA[same_chr & in_range, , drop = FALSE]
  message(sprintf("  Pairs after spatial filter (same chr, <= %d bp): %d",
                  max_dist, nrow(cA)))

  if (nrow(cA) == 0) {
    message("  No pairs passed spatial filter — returning NULL")
    return(NULL)
  }

  # ---- 4. Connected components via igraph (C implementation, fast) --------------
  all_nodes <- sort(unique(c(cA$queryHits, cA$subjectHits)))
  g <- igraph::graph_from_data_frame(
    d        = data.frame(from = cA$queryHits, to = cA$subjectHits),
    directed = FALSE,
    vertices = data.frame(name = all_nodes)
  )
  comp      <- igraph::components(g)
  comp_size <- tabulate(comp$membership)
  large_ids <- which(comp_size >= min_peaks)
  message(sprintf("  Components: %d total, %d with >= %d peaks",
                  length(comp_size), length(large_ids), min_peaks))

  if (length(large_ids) == 0) {
    message("  No hubs pass min_peaks filter — returning NULL")
    return(NULL)
  }

  # ---- 5. Vectorised hub collapsing --------------------------------------------
  # membership is named by node name (character peak index), valued by component ID
  node_ids  <- as.integer(names(comp$membership))
  comp_ids  <- comp$membership

  # Keep only large-component nodes and remap to 1..n_hubs
  keep_nodes  <- comp_ids %in% large_ids
  node_ids    <- node_ids[keep_nodes]
  hub_id_raw  <- comp_ids[keep_nodes]
  hub_id      <- match(hub_id_raw, large_ids)   # 1..n_hubs

  # Extract peak metadata for all hub peaks at once (one data.frame access)
  pk_df        <- as.data.frame(peaks[node_ids], row.names = NULL)
  pk_df$hub_id <- hub_id

  # Collapse: one row per hub via split + vectorised rbind
  hub_list <- lapply(split(pk_df, pk_df$hub_id), function(x) {
    data.frame(
      seqnames = x$seqnames[1],
      start    = min(x$start),
      end      = max(x$end),
      gene     = paste(unique(x$nearestGene), collapse = "-"),
      n_peaks  = nrow(x),
      stringsAsFactors = FALSE
    )
  })
  hub_df        <- do.call(rbind, hub_list)
  hub_gr        <- makeGRangesFromDataFrame(hub_df, keep.extra.columns = TRUE)
  names(hub_gr) <- paste0("HUB", seq_along(hub_gr))
  message(sprintf("  -> %d hubs detected", length(hub_gr)))

  # ---- 6. Peak-pair links GRanges (for arc / genome-browser plots) -------------
  s1 <- as.data.frame(resize(peaks[cA$queryHits],  1L, "center"),
                      row.names = NULL)[, c("seqnames", "start")]
  s2 <- as.data.frame(resize(peaks[cA$subjectHits], 1L, "center"),
                      row.names = NULL)[, "start", drop = FALSE]
  colnames(s2) <- "end"
  lnk <- data.frame(s1, end = s2$end, value = cA$correlation,
                    stringsAsFactors = FALSE)
  # Ensure start <= end
  lnk$end   <- pmax(lnk$start, lnk$end)
  lnk$start <- pmin(lnk$start, lnk$end)
  links_gr  <- makeGRangesFromDataFrame(lnk, keep.extra.columns = TRUE)

  # ---- Return ------------------------------------------------------------------
  list(
    hubsCollapsed = hub_gr,
    peakLinks     = links_gr,
    n_hubs        = length(hub_gr),
    params        = list(
      cor_cutoff   = cor_cutoff,
      max_dist     = max_dist,
      min_dist_tss = min_dist_tss,
      min_peaks    = min_peaks,
      remove_chr   = remove_chr
    )
  )
}
