# max_dist_sweep_SE.R
#
# Super-enhancer (SE) enrichment of hubs when the maximum distance parameter is changed, same pipeline as myeloid_SE_hub_enrichment.R
# (12.5 kb there):
#   CoAccess_r{R}_mp{N}_tss{D}_d{K}  co-accessibility hubs (hubs_finder, peak pairs at most K kb apart, midpoint to midpoint)
#   Stitching_mp{N}_tss{D}_d{K}      consensus peaks merged when the gap is <= K kb, DA hubs per cell type
#   Stitching_CT_mp{N}_tss{D}_d{K}   per-cell-type reproducible peaks merged the same way (all hubs of the cell type)
# K = MAX_DISTS (5, 10, 15, 20 kb; 12.5 kb is in myeloid_SE_hub_enrichment.csv). Hub x cell CPM matrices are cached in hub_cell_matrices_v2.
# Output: max_dist_sweep_SE.csv        Submit: bsub < submit_max_dist_sweep_SE.sh
# env: MAX_DISTS (5000,10000,15000,20000), N_CORES_MAT (3)

suppressPackageStartupMessages({ library(ArchR); library(GenomicRanges); library(Matrix); library(rtracklayer); library(igraph); library(presto); library(dplyr); library(parallel) })
addArchRGenome("hg38"); addArchRThreads(threads = 1)
source("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils/Hubs_finder.R")
script_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q7"
main_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR"; peakcalls <- file.path(main_dir, "PeakCalls")
se_rds     <- file.path(main_dir, "SE_regions_SE_database.rds"); se_path2 <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/SE_db2_from_Wooseung"
chain_file <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/liftover/hg38ToHg19.over.chain"; mat_cache <- file.path(script_dir, "hub_cell_matrices_v2")
MAX_DISTS <- as.integer(strsplit(Sys.getenv("MAX_DISTS", "5000,10000,15000,20000"), ",")[[1]]); COR_CUTOFFS <- c(0.2, 0.3, 0.4, 0.5, 0.6); MIN_PEAKS_VALUES <- c(3L, 4L, 5L); MIN_DIST_TSS_VALUES <- c(0L, 2000L)
REMOVE_CHR <- c("chrX", "chrY"); META_GROUP <- "celltype_lv1"; FDR_CUT <- 0.05; LFC_CUT <- 1
N_CORES <- min(as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "4")), 8L); N_CORES_MAT <- as.integer(Sys.getenv("N_CORES_MAT", "3"))

proj <- loadArchRProject(main_dir, showLogo = FALSE); cells <- rownames(proj@cellColData); celltypes <- as.character(proj@cellColData[, META_GROUP]); ct_levels <- unique(celltypes[!is.na(celltypes)]); all_peaks <- getPeakSet(proj)
SE_db1 <- readRDS(se_rds); f2 <- list.files(se_path2, pattern = "\\.bed$", full.names = TRUE)
SE_db2 <- setNames(lapply(f2, function(f) { y <- read.table(f, header = TRUE, sep = "\t")[, 1:3]; colnames(y) <- c("chr", "start", "end"); makeGRangesFromDataFrame(y) }), sub("\\.bed$", "", basename(f2))); SE_all <- c(SE_db1, SE_db2); ch <- import.chain(chain_file)
.liftover <- function(gr) { seqlevelsStyle(gr) <- "UCSC"; unlist(liftOver(gr, ch)) }
.phyper_enrich <- function(hub_h19, bg_h19, se_gr) { q <- sum(overlapsAny(hub_h19, se_gr, ignore.strand = TRUE)); m <- sum(overlapsAny(bg_h19, se_gr, ignore.strand = TRUE)); n <- length(bg_h19) - m; k <- length(hub_h19)
  list(q = q, k = k, frac_hub = q / k, frac_bg = m / length(bg_h19), pval = phyper(q = q, m = m, n = n, k = k, lower.tail = FALSE)) }
.filter_bg <- function(gr, min_tss = 0L) { gr <- gr[gr$distToGeneStart >= min_tss]; gr[!as.character(seqnames(gr)) %in% REMOVE_CHR] }
.stitch <- function(peaks, mp, d) { st <- reduce(peaks, min.gapwidth = d, ignore.strand = TRUE); np <- tabulate(subjectHits(findOverlaps(peaks, st)), nbins = length(st)); st[np >= mp] }
.safe_bind <- function(lst) { dfs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, lst); if (!length(dfs)) NULL else bind_rows(dfs) }
.row <- function(lbl, ct, h, bg, n_tot) .safe_bind(lapply(names(SE_all), function(s) { r <- .phyper_enrich(h, bg, SE_all[[s]]); data.frame(condition = lbl, celltype = ct, se_db = s, n_da_hubs = r$k, q_overlap = r$q, frac_hub_SE = r$frac_hub, frac_bg_SE = r$frac_bg, pval = r$pval, n_hubs_total = n_tot) }))
valid_cts <- ct_levels[file.exists(file.path(peakcalls, paste0(ct_levels, "-reproduciblePeaks.gr.rds")))]
kb <- function(d) sub("\\.0$", "", sprintf("%.1f", d / 1000))

## ---- hub sets ----
cA <- as.data.frame(getCoAccessibility(proj, corCutOff = min(COR_CUTOFFS), resolution = 1, returnLoops = FALSE), row.names = NULL)
cg <- expand.grid(r = COR_CUTOFFS, mp = MIN_PEAKS_VALUES, tss = MIN_DIST_TSS_VALUES, d = MAX_DISTS)
co <- mclapply(seq_len(nrow(cg)), function(i) { x <- cg[i, ]; lbl <- sprintf("CoAccess_r%.1f_mp%d_tss%d_d%s", x$r, x$mp, x$tss, kb(x$d))
  h <- hubs_finder(proj, coAccessibility = cA, cor_cutoff = x$r, max_dist = x$d, min_dist_tss = x$tss, min_peaks = x$mp, remove_chr = REMOVE_CHR); if (is.null(h)) return(NULL)
  list(lbl = lbl, tss = x$tss, hub_gr = h$hubsCollapsed, hub_ids = names(h$hubsCollapsed)) }, mc.cores = N_CORES, mc.preschedule = FALSE)
sg <- expand.grid(mp = MIN_PEAKS_VALUES, tss = MIN_DIST_TSS_VALUES, d = MAX_DISTS)
stl <- lapply(seq_len(nrow(sg)), function(i) { x <- sg[i, ]; lbl <- sprintf("Stitching_mp%d_tss%d_d%s", x$mp, x$tss, kb(x$d)); pk <- all_peaks[!as.character(seqnames(all_peaks)) %in% REMOVE_CHR & all_peaks$distToGeneStart >= x$tss]; st <- .stitch(pk, x$mp, x$d)
  list(lbl = lbl, tss = x$tss, hub_gr = st, hub_ids = paste0("ST_HUB", seq_along(st))) })
hubs <- Filter(Negate(is.null), c(co, stl)); for (i in seq_along(hubs)) hubs[[i]]$cache <- file.path(mat_cache, paste0("hub_cell_mat_", hubs[[i]]$lbl, ".rds"))
message(sprintf("%d hub sets: %s", length(hubs), paste(sprintf("%s=%d", sapply(hubs, `[[`, "lbl"), sapply(hubs, function(h) length(h$hub_gr))), collapse = ", ")))
todo <- Filter(function(x) !file.exists(x$cache), hubs)
if (length(todo)) { message("Loading fragments..."); fragments <- unlist(getFragmentsFromProject(ArchRProj = proj)); frag_rg <- as.character(fragments$RG); cell_idx <- setNames(seq_along(cells), cells); n_frags <- proj$nFrags
  mclapply(todo, function(e) { ovl <- findOverlaps(e$hub_gr, fragments, ignore.strand = TRUE); j <- cell_idx[frag_rg[subjectHits(ovl)]]; k <- !is.na(j)
    mat <- sparseMatrix(i = queryHits(ovl)[k], j = j[k], x = 1, dims = c(length(e$hub_gr), length(cells)), dimnames = list(e$hub_ids, cells)); saveRDS(t(t(mat) * (1e6 / n_frags)), e$cache); message("    saved ", basename(e$cache)); NULL }, mc.cores = N_CORES_MAT, mc.preschedule = FALSE)
  rm(fragments, frag_rg); gc() }
## ---- DA hubs + SE enrichment ----
res1 <- mclapply(hubs, function(e) { m <- readRDS(e$cache); da <- wilcoxauc(log2(m + 1), celltypes); da <- da[da$padj < FDR_CUT & da$logFC > LFC_CUT, ]; rm(m); out <- list()
  for (ct in valid_cts) { ids <- da$feature[da$group == ct]; if (!length(ids)) next; h <- .liftover(e$hub_gr[match(ids, e$hub_ids)]); if (!length(h)) next
    out[[ct]] <- .row(e$lbl, ct, h, .liftover(.filter_bg(readRDS(file.path(peakcalls, paste0(ct, "-reproduciblePeaks.gr.rds"))), e$tss)), length(e$hub_gr)) }
  .safe_bind(out) }, mc.cores = N_CORES, mc.preschedule = FALSE)
tg <- expand.grid(mp = MIN_PEAKS_VALUES, tss = MIN_DIST_TSS_VALUES, d = MAX_DISTS, ct = valid_cts, stringsAsFactors = FALSE)
res2 <- mclapply(seq_len(nrow(tg)), function(i) { x <- tg[i, ]; lbl <- sprintf("Stitching_CT_mp%d_tss%d_d%s", x$mp, x$tss, kb(x$d)); pk <- .filter_bg(readRDS(file.path(peakcalls, paste0(x$ct, "-reproduciblePeaks.gr.rds"))), x$tss); st <- .stitch(pk, x$mp, x$d); if (!length(st)) return(NULL)
  .row(lbl, x$ct, .liftover(st), .liftover(pk), length(st)) }, mc.cores = N_CORES, mc.preschedule = FALSE)
res <- bind_rows(.safe_bind(res1), .safe_bind(res2)) %>% group_by(condition, celltype) %>% mutate(fdr = p.adjust(pval, "BH"), neg_log10_fdr = -log10(pmax(fdr, 1e-300)), enriched = fdr < FDR_CUT) %>% ungroup()
write.csv(res, file.path(script_dir, "max_dist_sweep_SE.csv"), row.names = FALSE, quote = FALSE); message("Done: max_dist_sweep_SE.csv")
