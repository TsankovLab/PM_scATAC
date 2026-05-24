# myeloid_SE_hub_enrichment.R
#
# Computes Super Enhancer enrichment in hub regions across three hub types:
#   1. Co-accessibility hubs (r = 0.2–0.6, min_peaks = 3/4/5)
#      DA hubs per cell type via Wilcoxon (presto); hypergeometric SE test
#   2. Stitching hubs — consensus peak set reduced at 12.5 kb gap
#      Same DA + hypergeometric pipeline
#   3. Stitching_CT hubs — per-cell-type PeakCalls reduced at 12.5 kb gap
#      No DA needed; hypergeometric SE test directly
#
# Output: myeloid_SE_hub_enrichment.csv
# Plot:   run plot_SE_hub_comparison.R after this script completes
#
# Submit: bsub < submit_myeloid_SE_hub_enrichment.sh

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(Matrix)
  library(rtracklayer)
  library(igraph)
  library(presto)
  library(dplyr)
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
dir.create(mat_cache, showWarnings=FALSE)

# ---- parameters -------------------------------------------------------------
COR_CUTOFFS       <- c(0.2, 0.3, 0.4, 0.5, 0.6)
MAX_DIST          <- 12500L
MIN_PEAKS_VALUES  <- c(3L, 4L, 5L)
MIN_DIST_TSS      <- 0L
REMOVE_CHR   <- c("chrX", "chrY")
META_GROUP   <- "celltype_lv1"
FDR_CUT      <- 0.05
LFC_CUT      <- 1


# ---- load main ArchR project ------------------------------------------------
message("Loading main ArchR project...")
proj <- loadArchRProject(main_dir, showLogo=FALSE)
cells     <- rownames(proj@cellColData)
celltypes <- as.character(proj@cellColData[, META_GROUP])
ct_levels <- unique(celltypes[!is.na(celltypes)])
message(sprintf("  %d cells, %d peaks, %d cell types",
                nCells(proj), length(getPeakSet(proj)), length(ct_levels)))

# TSS-filtered distal peak set (used for all hub detection and background filtering)
all_peaks    <- getPeakSet(proj)
peaks_distal <- all_peaks[
  !as.character(seqnames(all_peaks)) %in% REMOVE_CHR &
  all_peaks$distToGeneStart >= MIN_DIST_TSS
]
message(sprintf("  %d distal peaks after TSS filter (>= %d bp, no chrX)",
                length(peaks_distal), MIN_DIST_TSS))


# ---- load SE databases (hg19) -----------------------------------------------
message("Loading SE databases...")
SE_db1 <- readRDS(se_rds)
SE_db2 <- lapply(list.files(se_path2, pattern="\\.bed$", full.names=TRUE), function(f) {
  y <- read.table(f, header=TRUE, sep="\t")[, 1:3]
  colnames(y) <- c("chr","start","end")
  makeGRangesFromDataFrame(y)
})
names(SE_db2) <- sub("\\.bed$", "", basename(list.files(se_path2, pattern="\\.bed$")))
SE_all <- c(SE_db1, SE_db2)
message(sprintf("  %d SE sets (%d db1 + %d db2)",
                length(SE_all), length(SE_db1), length(SE_db2)))
ch <- import.chain(chain_file)


# ---- load all fragments once ------------------------------------------------
message("Loading all fragments (done once, shared across conditions)...")
fragments <- unlist(getFragmentsFromProject(ArchRProj=proj))
message(sprintf("  %s fragments loaded", format(length(fragments), big.mark=",")))

# fast cell-index lookup; convert RG Rle to plain character vector once
frag_rg  <- as.character(fragments$RG)
cell_idx <- setNames(seq_along(cells), cells)
n_frags  <- proj$nFrags   # total fragments per cell for CPM scaling


# ---- helpers ----------------------------------------------------------------

# Build hub x cell CPM matrix (sparse) from hub GRanges + all fragments
.hub_cell_mat <- function(hub_gr, hub_ids, cache_path) {
  if (file.exists(cache_path)) {
    message("    Loading cached matrix...")
    return(readRDS(cache_path))
  }
  message(sprintf("    Computing hub x cell matrix (%d hubs x %d cells)...",
                  length(hub_gr), length(cells)))
  ovl    <- findOverlaps(hub_gr, fragments, ignore.strand=TRUE)
  q_hits <- queryHits(ovl)
  s_hits <- subjectHits(ovl)
  j      <- cell_idx[frag_rg[s_hits]]
  keep   <- !is.na(j)
  mat    <- sparseMatrix(
    i=q_hits[keep], j=j[keep], x=1,
    dims=c(length(hub_gr), length(cells)),
    dimnames=list(hub_ids, cells)
  )
  # CPM scaling
  mat <- t(t(mat) * (1e6 / n_frags))
  saveRDS(mat, cache_path)
  message(sprintf("    Saved to %s", basename(cache_path)))
  mat
}

# Wilcoxon DA: return hub IDs enriched in `celltype`
.da_hub_ids <- function(mat, celltype) {
  res <- wilcoxauc(log2(mat + 1), celltypes)
  res <- res[res$group == celltype & res$padj < FDR_CUT & res$logFC > LFC_CUT, ]
  res$feature
}

# LiftOver GRanges to hg19
.liftover <- function(gr) {
  seqlevelsStyle(gr) <- "UCSC"
  unlist(liftOver(gr, ch))
}

# Hypergeometric test against one SE set
.phyper_enrich <- function(hub_h19, bg_h19, se_gr) {
  q <- sum(overlapsAny(hub_h19, se_gr, ignore.strand=TRUE))
  m <- sum(overlapsAny(bg_h19,  se_gr, ignore.strand=TRUE))
  n <- length(bg_h19) - m
  k <- length(hub_h19)
  list(q=q, m=m, n=n, k=k,
       frac_hub=q/k, frac_bg=m/length(bg_h19),
       pval=phyper(q=q, m=m, n=n, k=k, lower.tail=FALSE))
}

# Filter background peak GRanges to distal peaks (>= MIN_DIST_TSS from TSS)
.filter_bg <- function(gr) {
  if ("distToGeneStart" %in% colnames(mcols(gr))) {
    gr <- gr[gr$distToGeneStart >= MIN_DIST_TSS]
  } else {
    gr <- gr[overlapsAny(gr, peaks_distal, ignore.strand = TRUE)]
  }
  gr[!as.character(seqnames(gr)) %in% REMOVE_CHR]
}

# hubs_finder() is sourced from git_repo/utils/Hubs_finder.R

# Derive stitching hub GRanges (uses TSS-filtered distal peaks)
.stitching_hubs <- function(min_peaks) {
  peaks    <- peaks_distal   # already TSS-filtered + chrX-removed
  stitched <- reduce(peaks, min.gapwidth=MAX_DIST, ignore.strand=TRUE)
  hits     <- findOverlaps(peaks, stitched)
  npeak    <- tabulate(subjectHits(hits), nbins=length(stitched))
  stitched <- stitched[npeak >= min_peaks]
  hub_ids  <- paste0("ST_HUB", seq_along(stitched))
  message(sprintf("    -> %d stitching hubs", length(stitched)))
  list(hub_gr=stitched, hub_ids=hub_ids)
}


# ---- build conditions -------------------------------------------------------
message("\nBuilding hub sets and hub x cell matrices...")

conditions <- list()

# Pre-compute co-accessibility once at the lowest threshold; reuse for all cutoffs
message("  Pre-computing co-accessibility at r >= ", min(COR_CUTOFFS), "...")
cA_precomputed <- as.data.frame(
  getCoAccessibility(proj, corCutOff=min(COR_CUTOFFS), resolution=1, returnLoops=FALSE),
  row.names=NULL)
message(sprintf("  %d pairs fetched", nrow(cA_precomputed)))

for (mp in MIN_PEAKS_VALUES) {
  message(sprintf("\n  --- min_peaks = %d ---", mp))

  for (r in COR_CUTOFFS) {
    lbl <- sprintf("CoAccess_r%.1f_mp%d", r, mp)
    message(sprintf("  %s...", lbl))
    hub_info   <- hubs_finder(proj,
                               coAccessibility = cA_precomputed,
                               cor_cutoff      = r,
                               max_dist        = MAX_DIST,
                               min_dist_tss    = MIN_DIST_TSS,
                               min_peaks       = mp,
                               remove_chr      = REMOVE_CHR)
    if (is.null(hub_info)) next
    hub_gr  <- hub_info$hubsCollapsed
    hub_ids <- names(hub_gr)
    cache_path <- file.path(mat_cache, paste0("hub_cell_mat_", lbl, ".rds"))
    mat        <- .hub_cell_mat(hub_gr, hub_ids, cache_path)
    conditions[[lbl]] <- list(hub_gr=hub_gr, hub_ids=hub_ids, mat=mat)
  }

  # Stitching with this min_peaks
  lbl <- sprintf("Stitching_mp%d", mp)
  message(sprintf("  %s...", lbl))
  st_info    <- .stitching_hubs(min_peaks = mp)
  cache_path <- file.path(mat_cache, paste0("hub_cell_mat_", lbl, ".rds"))
  mat_st     <- .hub_cell_mat(st_info$hub_gr, st_info$hub_ids, cache_path)
  conditions[[lbl]] <- list(hub_gr=st_info$hub_gr, hub_ids=st_info$hub_ids, mat=mat_st)
}


# ---- per-celltype SE enrichment loop ----------------------------------------
message("\nRunning DA + SE enrichment per condition and cell type...")

all_results <- list()

for (cond in names(conditions)) {
  info <- conditions[[cond]]
  message(sprintf("  Condition: %s (%d hubs)", cond, length(info$hub_ids)))

  for (ct in ct_levels) {
    pk_file <- file.path(peakcalls, paste0(ct, "-reproduciblePeaks.gr.rds"))
    if (!file.exists(pk_file)) next
    bg_h19 <- .liftover(.filter_bg(readRDS(pk_file)))

    # DA hubs for this cell type
    da_ids <- .da_hub_ids(info$mat, ct)
    if (length(da_ids) == 0) {
      message(sprintf("    %s / %s: 0 DA hubs", cond, ct))
      next
    }
    hub_idx  <- match(da_ids, info$hub_ids)
    hub_idx  <- hub_idx[!is.na(hub_idx)]
    hubs_gr  <- info$hub_gr[hub_idx]
    hubs_h19 <- .liftover(hubs_gr)
    message(sprintf("    %s / %s: %d DA hubs -> %d after liftover",
                    cond, ct, length(hubs_gr), length(hubs_h19)))
    if (length(hubs_h19) == 0) next

    rows <- lapply(names(SE_all), function(se_name) {
      r <- .phyper_enrich(hubs_h19, bg_h19, SE_all[[se_name]])
      data.frame(condition=cond, celltype=ct, se_db=se_name,
                 n_da_hubs=r$k, q_overlap=r$q,
                 frac_hub_SE=r$frac_hub, frac_bg_SE=r$frac_bg,
                 pval=r$pval, stringsAsFactors=FALSE)
    })
    all_results <- c(all_results, rows)
  }
}

da_results <- do.call(rbind, all_results)
rownames(da_results) <- NULL
da_results <- da_results %>%
  group_by(condition, celltype) %>%
  mutate(fdr           = p.adjust(pval, "BH"),
         neg_log10_fdr = -log10(pmax(fdr, 1e-300)),
         enriched      = fdr < FDR_CUT) %>%
  ungroup()


# ---- Stitching_CT: per-cell-type PeakCalls, no DA step ----------------------
message("\nRunning Stitching_CT enrichment per cell type...")

ct_rows <- list()
for (mp in MIN_PEAKS_VALUES) {
  cond_label <- sprintf("Stitching_CT_mp%d", mp)
  message(sprintf("  %s...", cond_label))

  for (ct in ct_levels) {
    pk_file <- file.path(peakcalls, paste0(ct, "-reproduciblePeaks.gr.rds"))
    if (!file.exists(pk_file)) next

    peaks_ct <- .filter_bg(readRDS(pk_file))
    stitched <- reduce(peaks_ct, min.gapwidth = MAX_DIST, ignore.strand = TRUE)
    hits     <- findOverlaps(peaks_ct, stitched)
    npeak    <- tabulate(subjectHits(hits), nbins = length(stitched))
    stitched <- stitched[npeak >= mp]
    if (length(stitched) == 0) next

    hubs_h19 <- .liftover(stitched)
    bg_h19   <- .liftover(peaks_ct)
    if (length(hubs_h19) == 0) next

    message(sprintf("    %s / %s: %d peaks -> %d hubs -> %d after liftover",
                    cond_label, ct, length(peaks_ct), length(stitched), length(hubs_h19)))

    rows <- lapply(names(SE_all), function(se_name) {
      r <- .phyper_enrich(hubs_h19, bg_h19, SE_all[[se_name]])
      data.frame(condition=cond_label, celltype=ct, se_db=se_name,
                 n_da_hubs=r$k, q_overlap=r$q,
                 frac_hub_SE=r$frac_hub, frac_bg_SE=r$frac_bg,
                 pval=r$pval, stringsAsFactors=FALSE)
    })
    ct_rows <- c(ct_rows, rows)
  }
}

ct_results <- do.call(rbind, ct_rows)
ct_results <- ct_results %>%
  group_by(condition, celltype) %>%
  mutate(fdr           = p.adjust(pval, "BH"),
         neg_log10_fdr = -log10(pmax(fdr, 1e-300)),
         enriched      = fdr < FDR_CUT) %>%
  ungroup()


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
