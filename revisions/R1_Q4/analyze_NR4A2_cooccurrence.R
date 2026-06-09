# analyze_NR4A2_cooccurrence.R  (R1_Q4)
#
# Tests whether NR4A2 motif sites in CD8_exhausted and NK_KLRC1 peaks
# co-occur with AP-1 motifs more than expected by chance — supporting
# the hypothesis that AP-1 co-binding dominates ChromBPNet contribution
# scores at those sites and suppresses the NR4A2 MoDISco signal.
#
# Also exports NR4A2 motif positions (per cell type) as BED files for
# the Python contribution score pileup script.
#
# Outputs → R1_Q4/output/NR4A2_analysis/
#   NR4A2_AP1_cooccurrence_barplot.pdf
#   NR4A2_top_comotifs_heatmap.pdf
#   NR4A2_cooccurrence_stats.csv
#   NR4A2_sites_<celltype>.bed   (for Python pileup)

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(paletteer)
})

addArchRGenome("hg38")
addArchRThreads(4)

proj_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR"
out_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/NR4A2_analysis"

CT_FOCUS <- c("CD8_exhausted", "NK_KLRC1")
CT_ALL   <- c("CD4", "CD8", "CD8_exhausted", "NK_FGFBP2", "NK_KLRC1", "Tregs")

# ── Load project ───────────────────────────────────────────────────────────────
message("Loading TNK ArchR project...")
proj <- loadArchRProject(proj_dir, showLogo = FALSE)
message("  Cells: ", nCells(proj), "  colData: ", paste(colnames(proj@cellColData), collapse=", "))

# ── Load motif match matrix (peaks × motifs, binary) ──────────────────────────
message("Loading motif match matrix...")
motif_se  <- readRDS(file.path(proj_dir, "Annotations", "Motif-Matches-In-Peaks.rds"))
motif_mat <- as.matrix(assay(motif_se))   # peaks × motifs (0/1)
peak_ct   <- rownames(motif_se)            # cell type label per peak

message("  Peaks: ", nrow(motif_mat), "  Motifs: ", ncol(motif_mat))
message("  Peak cell types: ", paste(sort(unique(peak_ct)), collapse=", "))

# ── Motif column indices ───────────────────────────────────────────────────────
nr4a2_col <- "NR4A2_686"
ap1_cols  <- grep("^(FOSL1|FOSL2|FOS_|JUNB|JUND|JUN_|BATF_|BATF3)", colnames(motif_mat),
                  value = TRUE, perl = TRUE)
message("  NR4A2 col: ", nr4a2_col)
message("  AP-1 cols (", length(ap1_cols), "): ", paste(ap1_cols, collapse=", "))

# Composite AP-1 flag: peak has ANY AP-1 motif
has_ap1  <- rowSums(motif_mat[, ap1_cols, drop = FALSE]) > 0
has_nr4a2 <- motif_mat[, nr4a2_col] == 1

# ── Co-occurrence statistics per cell type ────────────────────────────────────
message("Computing co-occurrence statistics...")
stats_list <- lapply(CT_ALL, function(ct) {
  idx       <- which(peak_ct == ct)
  n_peaks   <- length(idx)
  n_nr4a2   <- sum(has_nr4a2[idx])
  n_ap1     <- sum(has_ap1[idx])
  n_both    <- sum(has_nr4a2[idx] & has_ap1[idx])
  pct_ap1_given_nr4a2 <- if (n_nr4a2 > 0) n_both / n_nr4a2 * 100 else NA
  pct_ap1_bg          <- n_ap1 / n_peaks * 100

  # Fisher's exact test: NR4A2 × AP-1
  ct_tbl <- matrix(c(
    n_both,
    n_nr4a2 - n_both,
    n_ap1   - n_both,
    n_peaks - n_nr4a2 - n_ap1 + n_both
  ), nrow = 2, dimnames = list(NR4A2 = c("yes","no"), AP1 = c("yes","no")))
  ftest  <- fisher.test(ct_tbl, alternative = "greater")

  data.frame(celltype = ct, n_peaks, n_nr4a2, n_ap1, n_both,
             pct_ap1_given_nr4a2, pct_ap1_bg,
             OR = ftest$estimate, pval = ftest$p.value,
             stringsAsFactors = FALSE)
})
stats_df <- do.call(rbind, stats_list)
stats_df$padj  <- p.adjust(stats_df$pval, method = "BH")
stats_df$log2OR <- log2(stats_df$OR)

message("\nCo-occurrence stats:")
print(stats_df[, c("celltype","n_nr4a2","pct_ap1_given_nr4a2","pct_ap1_bg","OR","pval")])
write.csv(stats_df, file.path(out_dir, "NR4A2_cooccurrence_stats.csv"), row.names = FALSE)

# ── Barplot: % AP-1+ among NR4A2+ peaks by cell type ─────────────────────────
ct_order <- stats_df$celltype[order(stats_df$pct_ap1_given_nr4a2, decreasing = TRUE)]
plot_df  <- stats_df |>
  select(celltype, `NR4A2+ with AP-1` = pct_ap1_given_nr4a2,
         `Background AP-1 rate` = pct_ap1_bg) |>
  pivot_longer(-celltype, names_to = "metric", values_to = "pct")
plot_df$celltype <- factor(plot_df$celltype, levels = ct_order)

pal <- c("NR4A2+ with AP-1" = "#D73027", "Background AP-1 rate" = "#4292C6")
p_bar <- ggplot(plot_df, aes(x = celltype, y = pct, fill = metric)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = pal, name = NULL) +
  labs(x = NULL, y = "% peaks with AP-1 motif",
       title = "AP-1 co-occurrence at NR4A2+ peaks per cell type") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
        legend.position = "top", plot.title = element_text(size = 10))

ggsave(file.path(out_dir, "NR4A2_AP1_cooccurrence_barplot.pdf"), p_bar,
       width = 6, height = 4)
message("Barplot saved.")

# ── Top co-occurring motifs at NR4A2+ peaks per cell type ────────────────────
message("Computing top co-occurring motifs at NR4A2+ peaks...")
top_motifs <- lapply(CT_FOCUS, function(ct) {
  idx_nr4a2 <- which(peak_ct == ct & has_nr4a2)
  if (length(idx_nr4a2) < 5) return(NULL)
  submat  <- motif_mat[idx_nr4a2, , drop = FALSE]
  freqs   <- colMeans(submat)
  freqs   <- sort(freqs, decreasing = TRUE)
  # Strip numeric suffix
  clean   <- gsub("_\\d+$", "", names(freqs))
  names(freqs) <- clean
  head(freqs, 40)
})
names(top_motifs) <- CT_FOCUS

# Heatmap of top 30 unique motifs across focus cell types
valid_top <- top_motifs[!sapply(top_motifs, is.null)]
all_top   <- unique(unlist(lapply(valid_top, names)))[1:min(30, length(unique(unlist(lapply(valid_top, names)))))]
hm_mat <- do.call(rbind, lapply(CT_FOCUS, function(ct) {
  v   <- top_motifs[[ct]]
  out <- setNames(rep(0, length(all_top)), all_top)
  if (!is.null(v)) out[intersect(all_top, names(v))] <- v[intersect(all_top, names(v))]
  out
}))
rownames(hm_mat) <- CT_FOCUS

col_fun <- colorRamp2(c(0, 0.3, 1), c("white", "#FDAE61", "#D73027"))
ht <- Heatmap(hm_mat,
  name = "Motif freq.\nat NR4A2+ peaks",
  col  = col_fun,
  cluster_rows = TRUE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_title = "Top co-occurring motifs at NR4A2+ peaks",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  width = unit(3, "cm"), height = unit(length(all_top) * 0.45, "cm")
)

pdf(file.path(out_dir, "NR4A2_top_comotifs_heatmap.pdf"),
    width = 5, height = max(4, length(all_top) * 0.4 + 2))
draw(ht)
dev.off()
message("Co-motif heatmap saved.")

# ── Export NR4A2 motif positions as BED for pileup ────────────────────────────
message("Exporting NR4A2 motif positions...")
motif_pos <- readRDS(file.path(proj_dir, "Annotations", "Motif-Positions-In-Peaks.rds"))
peaks_gr  <- getPeakSet(proj)
peaks_gr$celltype <- peak_ct   # attach cell type label

for (ct in CT_FOCUS) {
  ct_peaks   <- peaks_gr[peaks_gr$celltype == ct]
  nr4a2_gr   <- motif_pos[["NR4A2_686"]]
  # Keep NR4A2 sites overlapping this cell type's peaks
  ol         <- findOverlaps(nr4a2_gr, ct_peaks)
  nr4a2_ct   <- nr4a2_gr[queryHits(ol)]
  nr4a2_ct   <- unique(nr4a2_ct)
  # Write BED: chr, start-1 (0-based), end, name, score, strand
  bed_df <- data.frame(
    chr    = as.character(seqnames(nr4a2_ct)),
    start  = start(nr4a2_ct) - 1L,
    end    = end(nr4a2_ct),
    name   = paste0("NR4A2_", seq_along(nr4a2_ct)),
    score  = round(nr4a2_ct$score, 3),
    strand = as.character(strand(nr4a2_ct))
  )
  bed_path <- file.path(out_dir, paste0("NR4A2_sites_", ct, ".bed"))
  write.table(bed_df, bed_path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  message("  ", ct, ": ", nrow(bed_df), " NR4A2 sites → ", bed_path)
}

# Also export FOSL2 sites (AP-1 positive control for pileup)
for (ct in CT_FOCUS) {
  ct_peaks   <- peaks_gr[peaks_gr$celltype == ct]
  fosl2_col  <- grep("^FOSL2_", names(motif_pos), value = TRUE)[1]
  fosl2_gr   <- motif_pos[[fosl2_col]]
  ol         <- findOverlaps(fosl2_gr, ct_peaks)
  fosl2_ct   <- unique(fosl2_gr[queryHits(ol)])
  bed_df <- data.frame(
    chr = as.character(seqnames(fosl2_ct)),
    start = start(fosl2_ct) - 1L,
    end   = end(fosl2_ct),
    name  = paste0("FOSL2_", seq_along(fosl2_ct)),
    score = round(fosl2_ct$score, 3),
    strand = as.character(strand(fosl2_ct))
  )
  bed_path <- file.path(out_dir, paste0("FOSL2_sites_", ct, ".bed"))
  write.table(bed_df, bed_path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  message("  ", ct, ": ", nrow(bed_df), " FOSL2 sites → ", bed_path)
}

message("Done.")
