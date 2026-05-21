# myeloid_SE_hub_enrichment.R
#
# Compares Super Enhancer enrichment in cHub regions across:
#   - Co-accessibility hubs at r = 0.2, 0.3, 0.4, 0.5, 0.6
#   - Proximity stitching hubs (gap <= 12,500bp, >= 5 peaks)
#
# Exactly replicates scatac_main_hub_analysis.R logic for each condition:
#   1. Build hub regions (hubsCollapsed GRanges)
#   2. Count fragments per cell overlapping each hub → hub x cell matrix
#      (using sparse Matrix for efficiency; same result as cell-by-cell loop)
#   3. Scale by total fragments per cell (CPM)
#   4. Wilcoxon DA test per cell type (presto)  →  DA hubs (padj<0.05, logFC>1)
#   5. Hypergeometric test: DA hubs vs PeakCalls/{celltype} background per SE set
#
# Hub x cell matrices are cached in git_repo_claude/hub_cell_matrices/ so they
# are not recomputed on reruns.
#
# Outputs (git_repo_claude/):
#   myeloid_SE_hub_enrichment.csv
#   plot_myeloid_SE_hub_enrichment.pdf   — heatmap per condition (SE x celltype)
#   plot_myeloid_SE_hub_summary.pdf      — bar + scatter summary
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
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})
addArchRGenome("hg38")
addArchRThreads(threads = 1)

source("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils/Hubs_finder.R")

# ---- paths ------------------------------------------------------------------
script_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude"
main_dir    <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR"
peakcalls   <- file.path(main_dir, "PeakCalls")
se_rds      <- file.path(main_dir, "SE_regions_SE_database.rds")
se_path2    <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/SE_db2_from_Wooseung"
chain_file  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/liftover/hg38ToHg19.over.chain"
mat_cache   <- file.path(script_dir, "hub_cell_matrices_v2")
dir.create(mat_cache, showWarnings=FALSE)

# ---- parameters -------------------------------------------------------------
COR_CUTOFFS  <- c(0.2, 0.3, 0.4, 0.5, 0.6)
MAX_DIST     <- 12500L
MIN_PEAKS    <- 5L
MIN_DIST_TSS <- 2000L
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
.stitching_hubs <- function() {
  peaks    <- peaks_distal   # already TSS-filtered + chrX-removed
  stitched <- reduce(peaks, min.gapwidth=MAX_DIST, ignore.strand=TRUE)
  hits     <- findOverlaps(peaks, stitched)
  npeak    <- tabulate(subjectHits(hits), nbins=length(stitched))
  stitched <- stitched[npeak >= MIN_PEAKS]
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

for (r in COR_CUTOFFS) {
  lbl <- sprintf("CoAccess_r%.1f", r)
  message(sprintf("  %s...", lbl))
  hub_info   <- hubs_finder(proj,
                             coAccessibility = cA_precomputed,
                             cor_cutoff      = r,
                             max_dist        = MAX_DIST,
                             min_dist_tss    = MIN_DIST_TSS,
                             min_peaks       = MIN_PEAKS,
                             remove_chr      = REMOVE_CHR)
  if (is.null(hub_info)) next
  hub_gr  <- hub_info$hubsCollapsed
  hub_ids <- names(hub_gr)
  cache_path <- file.path(mat_cache, paste0("hub_cell_mat_", lbl, ".rds"))
  mat        <- .hub_cell_mat(hub_gr, hub_ids, cache_path)
  conditions[[lbl]] <- list(hub_gr=hub_gr, hub_ids=hub_ids, mat=mat)
}

# Stitching: derive hubs, compute/load matrix
message("  Stitching...")
st_info    <- .stitching_hubs()
cache_path <- file.path(mat_cache, "hub_cell_mat_Stitching.rds")
mat_st     <- .hub_cell_mat(st_info$hub_gr, st_info$hub_ids, cache_path)
conditions[["Stitching"]] <- list(hub_gr=st_info$hub_gr, hub_ids=st_info$hub_ids, mat=mat_st)


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

results <- do.call(rbind, all_results)
rownames(results) <- NULL

results <- results %>%
  group_by(condition, celltype) %>%
  mutate(fdr           = p.adjust(pval, "BH"),
         neg_log10_fdr = -log10(pmax(fdr, 1e-300)),
         enriched      = fdr < FDR_CUT) %>%
  ungroup()

write.csv(results, file.path(script_dir, "myeloid_SE_hub_enrichment.csv"),
          row.names=FALSE, quote=FALSE)

sig_summary <- results %>%
  filter(enriched) %>%
  count(condition, celltype, name="n_sig_SE")
message("\nSignificant SE sets (FDR < 0.05) per condition x celltype:")
print(as.data.frame(sig_summary[order(sig_summary$celltype, sig_summary$condition),]),
      row.names=FALSE)


# ---- Plot 1: heatmap per condition (SE x celltype) --------------------------
message("\nGenerating plots...")

cond_order <- c(sprintf("CoAccess_r%.1f", COR_CUTOFFS), "Stitching")
cond_order <- cond_order[cond_order %in% unique(results$condition)]
ct_order   <- c("Malignant","Mesothelium","Alveolar","Fibroblasts","SmoothMuscle",
                "Endothelial","Myeloid","T_cells","NK","B_cells","Plasma","pDCs")
ct_order   <- ct_order[ct_order %in% unique(results$celltype)]
fdr_thresh <- -log10(FDR_CUT)

pdf(file.path(script_dir, "plot_myeloid_SE_hub_enrichment.pdf"), width=7, height=5)
for (cond in cond_order) {
  sub_df <- results[results$condition == cond, ]
  mat <- tapply(sub_df$neg_log10_fdr,
                list(sub_df$se_db, sub_df$celltype), mean)
  mat[is.na(mat)] <- 0
  mat <- mat[, intersect(ct_order, colnames(mat)), drop=FALSE]
  mat[mat < fdr_thresh] <- 0
  col_fun <- colorRamp2(c(0, fdr_thresh, max(mat, fdr_thresh+0.01)),
                        c("white","#fee8c8","#e34a33"))
  draw(Heatmap(mat,
    name            = "-log10(FDR)",
    col             = col_fun,
    cluster_columns = FALSE,
    cluster_rows    = TRUE,
    show_row_dend   = FALSE,
    row_names_gp    = gpar(fontsize=6),
    column_names_gp = gpar(fontsize=9),
    column_names_rot= 45,
    border          = TRUE,
    column_title    = cond,
    column_title_gp = gpar(fontsize=10, fontface="bold")
  ))
}
dev.off()


# ---- Plot 2 & 3: summary (bar + scatter) ------------------------------------
cond_colours <- setNames(
  c(colorRampPalette(c("#1a9850","#d73027"))(length(COR_CUTOFFS)), "#4575b4"),
  cond_order
)

bar_df <- results %>%
  count(condition, celltype, wt=as.integer(enriched)) %>%
  rename(n_sig=n) %>%
  mutate(condition=factor(condition, levels=cond_order),
         celltype =factor(celltype,  levels=ct_order))

hub_counts <- data.frame(
  condition = names(conditions),
  n_hubs    = sapply(names(conditions), function(x) length(conditions[[x]]$hub_ids)),
  stringsAsFactors=FALSE
)

line_df <- sig_summary %>%
  left_join(hub_counts, by="condition") %>%
  mutate(condition=factor(condition, levels=cond_order),
         celltype =factor(celltype,  levels=ct_order))

p_bar <- ggplot(bar_df, aes(x=condition, y=n_sig, fill=condition)) +
  geom_col(width=0.7) +
  geom_text(aes(label=n_sig), vjust=-0.3, size=2.5) +
  facet_wrap(~celltype, nrow=2) +
  scale_fill_manual(values=cond_colours, guide="none") +
  labs(title="Significant SE sets per hub method (FDR < 0.05)",
       x=NULL, y="# significant SE sets") +
  theme_bw(base_size=9) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
        strip.text=element_text(size=8))

p_line <- ggplot(line_df, aes(x=n_hubs, y=n_sig_SE, colour=condition)) +
  geom_point(size=2.5) +
  ggrepel::geom_text_repel(aes(label=condition), size=2.3,
                             colour="black", max.overlaps=20) +
  scale_colour_manual(values=cond_colours, guide="none") +
  facet_wrap(~celltype, nrow=2) +
  labs(title="Hub count vs SE enrichment per cell type",
       x="Number of hubs", y="# significant SE sets (FDR < 0.05)") +
  theme_bw(base_size=9) +
  theme(strip.text=element_text(size=8))

pdf(file.path(script_dir, "plot_myeloid_SE_hub_summary.pdf"), width=10, height=6)
print(p_bar)
print(p_line)
dev.off()

message("Done.")
