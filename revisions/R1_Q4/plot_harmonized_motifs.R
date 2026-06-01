# plot_harmonized_motifs.R
#
# Generates a two-page PDF summarising the harmonised ChromBPNet motifs
# for momac, resident, and PBMC_CD14_Mono:
#
#   Page 1 – Panel A (seqlet composition bar chart) + Panel B (peak abundance heatmap)
#   Page 2 – Panel C (CWM logos for the top N merged motifs)
#
# Inputs (all under modisco_merged_profile/):
#   compiled/modisco_compiled.tsv   — motif x cell-type seqlet table
#   compiled/modisco_compiled.h5    — CWMs (contrib_scores 30×4 per pattern)
#   compiled_tomtom/tomtom/*.tsv    — TOMTOM best-match per merged pattern
#   {ct}/no_bias_model/finemo_out_profile/hits.tsv  — per-cell-type FI-NeMo hits
#
# Output: plot_harmonized_motifs.pdf  (in the same directory as this script)

suppressPackageStartupMessages({
  library(rhdf5)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggseqlogo)
  library(patchwork)
  library(scales)
})
setDTthreads(1)

# ---- paths ------------------------------------------------------------------
CHROMBPDIR <- "/sc/arion/scratch/giottb01/scatac_meso/revisions/R1_Q4"
MERGED_DIR <- file.path(CHROMBPDIR, "modisco_merged_profile")
COMPILED_H5  <- file.path(MERGED_DIR, "compiled/modisco_compiled.h5")
COMPILED_TSV <- file.path(MERGED_DIR, "compiled/modisco_compiled.tsv")
TOMTOM_DIR   <- file.path(MERGED_DIR, "compiled_tomtom/tomtom")
SCRIPT_DIR   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4"
CELLTYPES    <- c("momac", "resident", "PBMC_CD14_Mono")
MODEL_HEAD   <- "profile"
TOP_N        <- 40   # motifs shown in bar chart + abundance heatmap
TOP_LOGO     <- 30   # motifs shown in logo panel

CT_COLORS  <- c(momac = "#E64B35", resident = "#4DBBD5", PBMC_CD14_Mono = "#00A087")
CT_LABELS  <- c(momac = "Momac", resident = "Resident", PBMC_CD14_Mono = "PBMC Mono")

# ---- 1. Compiled modisco table ----------------------------------------------
message("Loading compiled modisco table...")
compiled <- fread(COMPILED_TSV, sep="\t") %>% as_tibble()
colnames(compiled) <- c("merged_pattern", "component_celltype", "pattern_class", "pattern", "n_seqlets")

motif_total <- compiled %>%
  group_by(merged_pattern, pattern_class) %>%
  summarise(total_seqlets = sum(n_seqlets), .groups="drop") %>%
  arrange(desc(total_seqlets))

motif_ct <- compiled %>%
  group_by(merged_pattern, component_celltype) %>%
  summarise(n_seqlets = sum(n_seqlets), .groups="drop")

# ---- 2. TOMTOM best match per merged pattern --------------------------------
message("Loading TOMTOM matches...")
tomtom_files <- list.files(TOMTOM_DIR, pattern="\\.tomtom\\.tsv$", full.names=TRUE)

strip_tf <- function(id) {
  # "JUN_MA0488.1" -> "JUN",  "ATF2_HUMAN.H11MO.0.B" -> "ATF2"
  id <- sub("_HUMAN\\..*$", "", id, ignore.case=TRUE)
  id <- sub("_MOUSE\\..*$", "", id, ignore.case=TRUE)
  id <- sub("_[A-Z]{2}[0-9]+\\..*$", "", id)
  id
}

tomtom_best <- lapply(tomtom_files, function(f) {
  pat <- sub("\\.tomtom\\.tsv$", "", basename(f))
  d   <- tryCatch(fread(f, sep="\t") %>% as_tibble(), error=function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || !("q-value" %in% colnames(d)))
    return(tibble(merged_pattern=pat, tf_name=NA_character_, tomtom_qval=NA_real_))
  d <- d %>% rename(qval=`q-value`) %>% filter(!is.na(qval), !is.na(Target_ID))
  if (nrow(d) == 0)
    return(tibble(merged_pattern=pat, tf_name=NA_character_, tomtom_qval=NA_real_))
  best <- d %>% slice_min(qval, n=1, with_ties=FALSE)
  tibble(merged_pattern=pat, tf_name=strip_tf(best$Target_ID), tomtom_qval=best$qval)
})
tomtom_df <- bind_rows(tomtom_best)

# ---- 3. Motif metadata and short labels -------------------------------------
motif_meta <- motif_total %>%
  left_join(tomtom_df, by="merged_pattern") %>%
  mutate(
    cluster_id = sub("__merged_pattern_", "_p", sub("Average_", "Avg", merged_pattern)),
    label = case_when(
      !is.na(tf_name) & tomtom_qval < 0.05 ~ paste0(tf_name, " [", cluster_id, "]"),
      TRUE ~ cluster_id
    )
  )

top_motifs  <- motif_meta %>% head(TOP_N)
top_pat     <- top_motifs$merged_pattern
label_map   <- setNames(top_motifs$label, top_pat)

# ---- 4. FI-NeMo peak abundance per cell type --------------------------------
message("Loading FI-NeMo hits (this may take a moment)...")

peak_abundance <- lapply(CELLTYPES, function(ct) {
  hits_f <- file.path(CHROMBPDIR, ct, "no_bias_model",
                      paste0("finemo_out_", MODEL_HEAD), "hits.tsv")
  pk_f   <- file.path(CHROMBPDIR, ct, paste0(ct, "_peakset_all_no_blacklist.bed"))
  if (!file.exists(hits_f)) { message("  Missing hits: ", hits_f); return(NULL) }

  n_peaks <- fread(pk_f, header=FALSE, select=1L) %>% nrow()
  message(sprintf("  %s: %d peaks", ct, n_peaks))

  hits <- fread(hits_f, sep="\t") %>% as_tibble()
  colnames(hits)[1:13] <- c("chr","start","end","start_unt","end_unt","motif_name",
                             "hit_coef","hit_coef_global","hit_corr","hit_importance",
                             "strand","peak_name","peak_id")

  # parse pattern class and pattern number from motif_name
  hits <- hits %>%
    mutate(
      pattern_class = sub("\\..*$", "", motif_name),
      pattern       = sub("^[^.]+\\.", "", motif_name)
    )

  # join to compiled table to get merged_pattern
  ct_map <- compiled %>%
    filter(component_celltype == ct) %>%
    select(merged_pattern, pattern_class, pattern)

  hits <- hits %>% left_join(ct_map, by=c("pattern_class","pattern"))

  hits %>%
    filter(!is.na(merged_pattern)) %>%
    group_by(merged_pattern) %>%
    summarise(
      n_peaks_hit   = n_distinct(peak_id),
      mean_importance = mean(hit_importance),
      .groups = "drop"
    ) %>%
    mutate(frac_peaks = n_peaks_hit / n_peaks, celltype = ct)
}) %>% bind_rows()

# ---- 5. CWMs from h5 --------------------------------------------------------
message("Extracting CWMs from h5...")

get_cwm <- function(mp, pcls) {
  key <- paste0(pcls, "/", mp, "/contrib_scores")
  tryCatch({
    cwm <- h5read(COMPILED_H5, key)  # rhdf5 returns (4, 30): rows=ACGT, cols=positions
    rownames(cwm) <- c("A","C","G","T")
    cwm
  }, error=function(e) NULL)
}

logo_pats <- motif_meta %>% head(TOP_LOGO)

cwm_list <- lapply(seq_len(nrow(logo_pats)), function(i) {
  get_cwm(logo_pats$merged_pattern[i], logo_pats$pattern_class[i])
})
names(cwm_list) <- logo_pats$merged_pattern
cwm_list <- cwm_list[!sapply(cwm_list, is.null)]

# ---- 6. Panel A: seqlet composition bar chart --------------------------------
message("Building Panel A: seqlet composition...")

comp_df <- motif_ct %>%
  filter(merged_pattern %in% top_pat) %>%
  mutate(
    merged_pattern    = factor(merged_pattern, levels=rev(top_pat)),
    component_celltype = factor(component_celltype, levels=CELLTYPES)
  )

pA <- ggplot(comp_df, aes(x=n_seqlets, y=merged_pattern, fill=component_celltype)) +
  geom_col(orientation="y", width=0.8) +
  scale_fill_manual(values=CT_COLORS, labels=CT_LABELS, name=NULL) +
  scale_x_continuous(labels=comma_format(), expand=expansion(mult=c(0,0.05))) +
  scale_y_discrete(labels=label_map) +
  labs(x="Number of seqlets", y=NULL, title="A  Seqlet composition") +
  theme_minimal(base_size=9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y  = element_text(size=7),
    legend.position = "bottom",
    legend.text = element_text(size=8),
    plot.title = element_text(face="bold", size=10)
  )

# ---- 7. Panel B: peak abundance heatmap -------------------------------------
message("Building Panel B: peak abundance...")

abund_wide <- peak_abundance %>%
  filter(merged_pattern %in% top_pat) %>%
  select(merged_pattern, celltype, frac_peaks) %>%
  complete(merged_pattern = top_pat, celltype = CELLTYPES, fill=list(frac_peaks=0)) %>%
  mutate(
    merged_pattern = factor(merged_pattern, levels=rev(top_pat)),
    celltype       = factor(celltype, levels=CELLTYPES)
  )

pB <- ggplot(abund_wide, aes(x=celltype, y=merged_pattern, fill=frac_peaks)) +
  geom_tile(color="white", linewidth=0.4) +
  scale_fill_gradient(
    low="grey97", high="#C0392B",
    name="Fraction\npeaks",
    labels=percent_format(accuracy=1),
    limits=c(0, NA)
  ) +
  scale_y_discrete(labels=label_map) +
  scale_x_discrete(labels=CT_LABELS) +
  labs(x=NULL, y=NULL, title="B  Peak-level abundance") +
  theme_minimal(base_size=9) +
  theme(
    axis.text.x  = element_text(angle=35, hjust=1, size=8),
    axis.text.y  = element_text(size=7),
    panel.grid   = element_blank(),
    legend.position = "right",
    plot.title   = element_text(face="bold", size=10)
  )

# ---- 8. Panel C: CWM logos --------------------------------------------------
message("Building Panel C: CWM logos...")

avail_cwm  <- names(cwm_list)
logo_order <- logo_pats$merged_pattern[logo_pats$merged_pattern %in% avail_cwm]

logo_label_map <- setNames(logo_pats$label, logo_pats$merged_pattern)
class_map      <- setNames(logo_pats$pattern_class, logo_pats$merged_pattern)

logo_plots <- lapply(logo_order, function(mp) {
  cwm   <- cwm_list[[mp]]
  lbl   <- logo_label_map[mp]
  cls   <- class_map[mp]
  pcol  <- if (grepl("^pos", cls)) "#2166AC" else "#B2182B"
  ggseqlogo(cwm, method="custom", seq_type="dna") +
    ggtitle(lbl) +
    theme(
      plot.title      = element_text(size=6, hjust=0.5, color=pcol, face="bold"),
      axis.text       = element_blank(),
      axis.title      = element_blank(),
      axis.ticks      = element_blank(),
      plot.margin     = margin(2, 2, 2, 2, "mm"),
      panel.background = element_blank()
    )
})

message(sprintf("  %d logos generated out of %d requested", length(logo_plots), length(logo_order)))
if (length(logo_plots) == 0) stop("No CWMs extracted — check h5 key format")

ncols <- 5
pC <- wrap_plots(logo_plots, ncol=ncols) +
  plot_annotation(
    title   = "C  CWM logos (top merged motifs by seqlet count)",
    caption = "Blue title = pos_patterns; Red title = neg_patterns",
    theme   = theme(
      plot.title   = element_text(face="bold", size=11),
      plot.caption = element_text(size=8, color="grey50")
    )
  )

# ---- 9. Save PDF ------------------------------------------------------------
outpdf <- file.path(SCRIPT_DIR, "plot_harmonized_motifs.pdf")
message("Saving ", outpdf)

pdf(outpdf, width=18, height=22)
print(
  (pA | pB) +
    plot_layout(widths=c(2.5, 1)) +
    plot_annotation(
      title   = "Harmonised ChromBPNet motifs: momac | resident | PBMC_CD14_Mono",
      subtitle = sprintf("Top %d merged motifs (of %d total) by seqlet count - profile model head",
                         TOP_N, nrow(motif_total)),
      theme   = theme(
        plot.title    = element_text(face="bold", size=12),
        plot.subtitle = element_text(size=9, color="grey40")
      )
    )
)
print(pC)
dev.off()
message("Done.")
