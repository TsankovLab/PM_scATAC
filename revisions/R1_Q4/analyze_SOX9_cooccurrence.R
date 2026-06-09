# analyze_SOX9_cooccurrence.R  (R1_Q4)
#
# Tests whether SOX9 motif sites in SOX9_regulon_high P23 tumor peaks
# co-occur with TEAD and AP-1 motifs — and whether this co-binding could
# explain why ChromBPNet/MoDISco does not recover a SOX9 motif pattern
# despite strong SCENIC regulon activity.
#
# Peak set: MACS2 peaks from the SOX9_regulon_high ChromBPNet run
# (the exact peaks used to train the model).
# Motif annotations: from tumor ArchR Motif-Matches-In-Peaks.rds
# (overlapping peaks are matched by genomic coordinates).
#
# Additional hypotheses tested:
#   1. SOX9 sites co-occur with TEAD (YAP/TAZ pathway in mesothelioma)
#   2. SOX9 sites co-occur with AP-1
#   3. SOX9 motif frequency at SOX9_high peaks vs genome-wide background
#
# Outputs → R1_Q4/output/SOX9_analysis/
#   SOX9_comotif_barplot.pdf
#   SOX9_top_comotifs_heatmap.pdf
#   SOX9_cooccurrence_stats.csv
#   SOX9_sites_SOX9_regulon_high.bed    (for pileup)
#   TEAD_sites_SOX9_regulon_high.bed    (positive control for pileup)

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(paletteer)
  library(SummarizedExperiment)
})

addArchRGenome("hg38")
addArchRThreads(4)

proj_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR"
macs2_dir <- "/sc/arion/scratch/giottb01/SOX9_chromBPnet/MACS2_SOX9_regulon_high"
out_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/output/SOX9_analysis"

# ── Load tumor ArchR (for motif annotations and peak set) ─────────────────────
message("Loading tumor ArchR project...")
proj <- loadArchRProject(proj_dir, showLogo = FALSE)
message("  Cells: ", nCells(proj))

# ── Load MACS2 SOX9_high peak calls (the actual ChromBPNet training peaks) ────
message("Loading MACS2 SOX9_regulon_high peaks...")
narrow_path <- file.path(macs2_dir, "SOX9_regulon_high_peaks.narrowPeak")
narrow <- read.table(narrow_path, sep = "\t",
                     col.names = c("chr","start","end","name","score",
                                   "strand","fc","neglog10p","neglog10q","summit"))
sox9_high_peaks <- GRanges(seqnames = narrow$chr,
                            ranges   = IRanges(narrow$start + 1, narrow$end))
message("  SOX9_high peaks: ", length(sox9_high_peaks))

# ── Load tumor motif match matrix (all tumor peaks × 869 motifs) ──────────────
message("Loading tumor motif match matrix...")
motif_se  <- readRDS(file.path(proj_dir, "Annotations", "Motif-Matches-In-Peaks.rds"))
motif_mat <- as.matrix(assay(motif_se))
peaks_all <- rowRanges(motif_se)
message("  All tumor peaks: ", nrow(motif_mat))

# ── Overlap: find tumor peaks that overlap SOX9_high MACS2 peaks ──────────────
ol       <- findOverlaps(peaks_all, sox9_high_peaks)
sox9_idx <- unique(queryHits(ol))
bg_idx   <- setdiff(seq_len(nrow(motif_mat)), sox9_idx)
message("  Tumor peaks overlapping SOX9_high: ", length(sox9_idx))
message("  Background peaks: ", length(bg_idx))

# ── Motif columns of interest ─────────────────────────────────────────────────
sox9_col  <- "SOX9_756"
tead_cols <- grep("^TEAD", colnames(motif_mat), value = TRUE)
ap1_cols  <- grep("^(FOSL1|FOSL2|FOS_|JUNB|JUND|JUN_|BATF_|BATF3)",
                  colnames(motif_mat), value = TRUE, perl = TRUE)
message("  TEAD cols: ", paste(tead_cols, collapse=", "))
message("  AP-1 cols: ", paste(ap1_cols, collapse=", "))

has_sox9 <- motif_mat[, sox9_col] == 1
has_tead <- rowSums(motif_mat[, tead_cols, drop=FALSE]) > 0
has_ap1  <- rowSums(motif_mat[, ap1_cols,  drop=FALSE]) > 0

# ── Co-occurrence statistics: SOX9_high vs background peaks ───────────────────
cooccur_stats <- function(idx, label_fg, label_bg) {
  fg <- list(n = length(idx),
             sox9 = sum(has_sox9[idx]),
             tead = sum(has_tead[idx]),
             ap1  = sum(has_ap1[idx]),
             sox9_tead = sum(has_sox9[idx] & has_tead[idx]),
             sox9_ap1  = sum(has_sox9[idx] & has_ap1[idx]))
  fg$pct_sox9 <- fg$sox9 / fg$n * 100
  fg$pct_tead_given_sox9 <- if (fg$sox9 > 0) fg$sox9_tead / fg$sox9 * 100 else NA
  fg$pct_ap1_given_sox9  <- if (fg$sox9 > 0) fg$sox9_ap1  / fg$sox9 * 100 else NA
  fg$pct_tead_bg <- sum(has_tead[idx]) / fg$n * 100
  as.data.frame(c(group = label_fg, fg))
}

stats_high <- cooccur_stats(sox9_idx, "SOX9_regulon_high", "rest")
stats_bg   <- cooccur_stats(bg_idx,   "background", "rest")
stats_df   <- rbind(stats_high, stats_bg)
message("\nCo-occurrence summary:")
print(stats_df[, c("group","n","sox9","pct_sox9","sox9_tead","pct_tead_given_sox9","sox9_ap1","pct_ap1_given_sox9")])

# Fisher's test: SOX9 × TEAD at SOX9_high peaks
make_ftab <- function(flag1, flag2, idx) {
  matrix(c(sum( flag1[idx] &  flag2[idx]),
            sum(!flag1[idx] &  flag2[idx]),
            sum( flag1[idx] & !flag2[idx]),
            sum(!flag1[idx] & !flag2[idx])), nrow=2)
}
ft_tead <- fisher.test(make_ftab(has_sox9, has_tead, sox9_idx), alternative="greater")
ft_ap1  <- fisher.test(make_ftab(has_sox9, has_ap1,  sox9_idx), alternative="greater")
message("Fisher SOX9×TEAD: OR=", round(ft_tead$estimate,2), " p=", signif(ft_tead$p.value,3))
message("Fisher SOX9×AP1 : OR=", round(ft_ap1$estimate,2),  " p=", signif(ft_ap1$p.value,3))

write.csv(
  data.frame(test    = c("SOX9×TEAD", "SOX9×AP1"),
             OR      = c(ft_tead$estimate, ft_ap1$estimate),
             pval    = c(ft_tead$p.value,  ft_ap1$p.value),
             CI_low  = c(ft_tead$conf.int[1], ft_ap1$conf.int[1])),
  file.path(out_dir, "SOX9_cooccurrence_stats.csv"), row.names=FALSE
)

# ── Barplot: co-occurrence rates ───────────────────────────────────────────────
bar_df <- data.frame(
  group   = rep(c("SOX9_regulon_high", "Background"), 3),
  motif   = c("TEAD (at SOX9+ peaks)", "TEAD (at SOX9+ peaks)",
              "AP-1 (at SOX9+ peaks)", "AP-1 (at SOX9+ peaks)",
              "SOX9 motif freq.", "SOX9 motif freq."),
  pct     = c(as.numeric(stats_high$pct_tead_given_sox9),
              as.numeric(stats_bg$pct_tead_given_sox9),
              as.numeric(stats_high$pct_ap1_given_sox9),
              as.numeric(stats_bg$pct_ap1_given_sox9),
              as.numeric(stats_high$pct_sox9),
              as.numeric(stats_bg$pct_sox9))
)

p_bar <- ggplot(bar_df, aes(x = motif, y = pct, fill = group)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
  scale_fill_manual(values = c("SOX9_regulon_high"="#D73027","Background"="#4292C6"),
                    name = NULL) +
  labs(x = NULL, y = "% peaks", title = "Motif co-occurrence at SOX9+ peaks") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 9),
        legend.position = "top", plot.title = element_text(size = 10))

ggsave(file.path(out_dir, "SOX9_comotif_barplot.pdf"), p_bar, width = 5, height = 4)

# ── Top co-occurring motifs at SOX9+ SOX9_high peaks ─────────────────────────
sox9_pos_idx <- sox9_idx[has_sox9[sox9_idx]]
message("SOX9+ peaks in SOX9_high set: ", length(sox9_pos_idx))

if (length(sox9_pos_idx) >= 5) {
  freqs    <- colMeans(motif_mat[sox9_pos_idx, , drop=FALSE])
  freqs_bg <- colMeans(motif_mat[bg_idx[has_sox9[bg_idx]], , drop=FALSE])
  top40    <- names(sort(freqs, decreasing=TRUE))[1:min(40, length(freqs))]
  hm_mat   <- rbind("SOX9_high"  = freqs[top40],
                     "Background" = freqs_bg[top40])
  clean_names <- gsub("_\\d+$", "", top40)
  colnames(hm_mat) <- clean_names

  col_fun <- colorRamp2(c(0, 0.3, 1), c("white", "#FDAE61", "#D73027"))
  ht <- Heatmap(hm_mat,
    name = "Motif freq.\nat SOX9+ peaks",
    col  = col_fun,
    cluster_rows    = FALSE,
    cluster_columns = TRUE,
    row_names_gp    = gpar(fontsize = 9, fontface = "bold"),
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 45,
    column_title    = "Top 40 motifs at SOX9+ peaks (SOX9_high vs background)",
    column_title_gp = gpar(fontsize = 9, fontface = "bold"),
    width  = unit(14, "cm"), height = unit(1.8, "cm")
  )
  pdf(file.path(out_dir, "SOX9_top_comotifs_heatmap.pdf"), width = 14, height = 4)
  draw(ht)
  dev.off()
  message("Heatmap saved.")
}

# ── Export SOX9 and TEAD motif positions within SOX9_high peaks as BED ────────
message("Exporting motif positions...")
motif_pos <- readRDS(file.path(proj_dir, "Annotations", "Motif-Positions-In-Peaks.rds"))

export_bed <- function(motif_name, out_name, peaks_gr) {
  gr  <- motif_pos[[motif_name]]
  if (is.null(gr)) { message("  [WARN] motif not found: ", motif_name); return() }
  ol  <- findOverlaps(gr, peaks_gr)
  sub <- unique(gr[queryHits(ol)])
  bed <- data.frame(
    chr    = as.character(seqnames(sub)),
    start  = start(sub) - 1L,
    end    = end(sub),
    name   = paste0(gsub("_\\d+$","",motif_name), "_", seq_along(sub)),
    score  = round(sub$score, 3),
    strand = as.character(strand(sub))
  )
  path <- file.path(out_dir, paste0(out_name, ".bed"))
  write.table(bed, path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  message("  ", out_name, ": ", nrow(bed), " sites → ", path)
}

export_bed("SOX9_756", "SOX9_sites_SOX9_regulon_high", sox9_high_peaks)

# Export the TEAD motif with most hits as positive control
tead_counts <- sapply(tead_cols, function(col) sum(motif_mat[sox9_idx, col]))
best_tead   <- names(which.max(tead_counts))
export_bed(best_tead, "TEAD_sites_SOX9_regulon_high", sox9_high_peaks)

message("Done.")
