# R1_Q3 – F23: GSEA enrichment of per-sample co-expression with SOX9/RUNX2/SNAI2
#
# Input: F17_coexpr_all_genes_table.csv (Spearman r per gene x sample x anchor)
# For each anchor x sample: rank all genes by Spearman r → run fgsea
# Gene sets: Hallmark + GO BP
#
# Outputs:
#   F23a_gsea_per_sample_heatmap_{anchor}_{gs}.pdf  – NES heatmap samples x pathways
#   F23b_gsea_consistent_pathways.pdf               – pathways enriched in ≥50% samples
#   F23_gsea_per_sample_results.csv                 – full table

set.seed(1234)

projdir   <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
files_dir <- file.path(projdir, 'git_repo', 'files')
outdir    <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))
library(fgsea); library(patchwork); library(tidyr); library(dplyr)

anchor_genes   <- c('SOX9', 'RUNX2', 'SNAI2')
padj_threshold <- 0.1   # relaxed for per-sample (smaller n)
top_pw         <- 20    # top pathways to show in heatmap

##############################################################################
# Load gene sets
##############################################################################
load_gmt <- function(path) {
  lines <- readLines(path)
  setNames(lapply(lines, function(l) {
    parts <- strsplit(l, '\t')[[1]]
    parts[3:length(parts)]
  }), sapply(lines, function(l) strsplit(l, '\t')[[1]][1]))
}

message('Loading gene sets...')
hallmark <- load_gmt(file.path(files_dir, 'h.all.v7.4.symbols.gmt'))
gobp     <- load_gmt(file.path(files_dir, 'c5.bp.v7.1.symbol.gmt'))
names(hallmark) <- gsub('HALLMARK_', '', names(hallmark))
names(hallmark) <- gsub('_', ' ', names(hallmark))
names(gobp)     <- gsub('GOBP_', '', names(gobp))
names(gobp)     <- gsub('_', ' ', names(gobp))
gs_list <- list(Hallmark = hallmark, GOBP = gobp)

##############################################################################
# Load F17 per-sample co-expression
##############################################################################
message('Loading F17 per-sample correlations...')
f17 <- read.csv(file.path(outdir, 'F17_coexpr_all_genes_table.csv'),
                stringsAsFactors = FALSE)

# Sample order (low → high scS-score)
sarc_order    <- read.csv(file.path(projdir, 'tumor_compartment', 'scrna',
                                     'cnmf20_sarcomatoid_sample_order.csv'), row.names = 1)
tumor_sams    <- sarc_order$sampleID[!sarc_order$sampleID %in% c('HU37','HU62','normal_pleura')]
sample_levels <- tumor_sams[order(sarc_order$x[sarc_order$sampleID %in% tumor_sams])]
# keep only samples present in F17
sample_levels <- sample_levels[sample_levels %in% unique(f17$sample)]

##############################################################################
# Run fgsea per anchor x sample x gene set
##############################################################################
message('Running fgsea per sample...')
all_gsea <- list()

for (anc in anchor_genes) {
  for (sam in sample_levels) {
    sub <- f17[f17$anchor == anc & f17$sample == sam & !is.na(f17$rho), ]
    if (nrow(sub) < 50) {
      message('  Skipping ', anc, ' / ', sam, ' (n=', nrow(sub), ')')
      next
    }
    # Ranked vector: all genes sorted by Spearman r
    ranked <- setNames(sub$rho, sub$gene)
    ranked <- sort(ranked[!duplicated(names(ranked))], decreasing = TRUE)
    message('  ', anc, ' / ', sam, ' : n=', length(ranked))

    for (gs_name in names(gs_list)) {
      res <- tryCatch(
        fgseaMultilevel(gs_list[[gs_name]], ranked,
                        minSize = 10, maxSize = 500,
                        nPermSimple = 5000, eps = 0),
        error = function(e) { message('  fgsea error: ', e$message); NULL }
      )
      if (is.null(res)) next
      all_gsea[[paste(anc, sam, gs_name, sep = '|')]] <- data.frame(
        anchor   = anc,
        sample   = sam,
        gs       = gs_name,
        pathway  = res$pathway,
        NES      = res$NES,
        pval     = res$pval,
        padj     = res$padj,
        stringsAsFactors = FALSE
      )
    }
  }
}

gsea_df <- do.call(rbind, Filter(Negate(is.null), all_gsea))
write.csv(gsea_df, file.path(outdir, 'F23_gsea_per_sample_results.csv'),
          row.names = FALSE)
message('Total GSEA results: ', nrow(gsea_df))

##############################################################################
# F23a – NES heatmap: samples x pathways per anchor x gene set
##############################################################################
pal_div <- rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100))

for (gs_name in names(gs_list)) {
  hm_plots <- lapply(anchor_genes, function(anc) {
    sub <- gsea_df[gsea_df$anchor == anc & gsea_df$gs == gs_name, ]
    if (nrow(sub) == 0) return(NULL)

    # Select top pathways: those significant (padj < threshold) in ≥1 sample,
    # ranked by mean |NES| across samples
    sig_pws <- sub %>%
      group_by(pathway) %>%
      summarise(
        n_sig    = sum(padj < padj_threshold, na.rm = TRUE),
        mean_abs = mean(abs(NES), na.rm = TRUE),
        mean_NES = mean(NES, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      filter(n_sig >= 1) %>%
      slice_max(order_by = mean_abs, n = top_pw)

    if (nrow(sig_pws) == 0) return(NULL)

    hm_df <- sub[sub$pathway %in% sig_pws$pathway, ]
    hm_df$sample  <- factor(hm_df$sample, levels = sample_levels)
    # Order pathways by mean NES
    pw_order <- sig_pws$pathway[order(sig_pws$mean_NES)]
    hm_df$pathway <- factor(hm_df$pathway, levels = pw_order)

    # Mark significance
    hm_df$sig_label <- ifelse(hm_df$padj < padj_threshold,
                              ifelse(hm_df$padj < 0.01, '**', '*'), '')

    ggplot(hm_df, aes(x = sample, y = pathway, fill = NES)) +
      geom_tile(color = 'white', linewidth = 0.3) +
      geom_text(aes(label = sig_label), size = 2.5, vjust = 0.75) +
      scale_fill_gradientn(colours = pal_div, limits = c(-3, 3),
        oob = scales::squish, na.value = 'grey85', name = 'NES') +
      labs(x = 'Sample (low→high scS-score)', y = NULL,
           title = paste0(anc, ' – ', gs_name)) +
      theme_minimal(base_size = 8) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 35, hjust = 1),
            axis.text.y = element_text(size = 7))
  })
  hm_plots <- Filter(Negate(is.null), hm_plots)
  if (length(hm_plots) == 0) next

  max_pws <- max(sapply(hm_plots, function(p) length(unique(p$data$pathway))))
  pdf(file.path(outdir, 'Plots',
      paste0('F23a_gsea_per_sample_heatmap_', gs_name, '.pdf')),
      width = 4 + length(sample_levels) * 0.8,
      height = 3 + max_pws * 0.35 * length(anchor_genes))
  print(wrap_plots(hm_plots, ncol = 1))
  dev.off()
  message('F23a written for ', gs_name)
}

##############################################################################
# F23b – Consistent pathways: enriched in ≥50% of samples per anchor
##############################################################################
n_min <- max(2, floor(length(sample_levels) / 2))

consistent_plots <- lapply(anchor_genes, function(anc) {
  lapply(names(gs_list), function(gs_name) {
    sub <- gsea_df[gsea_df$anchor == anc & gsea_df$gs == gs_name, ]
    if (nrow(sub) == 0) return(NULL)

    consist <- sub %>%
      group_by(pathway) %>%
      summarise(
        n_sig       = sum(padj < padj_threshold, na.rm = TRUE),
        mean_NES    = mean(NES, na.rm = TRUE),
        direction   = ifelse(mean_NES > 0, 'Positive', 'Negative'),
        .groups = 'drop'
      ) %>%
      filter(n_sig >= n_min) %>%
      slice_max(order_by = abs(mean_NES), n = 20)

    if (nrow(consist) == 0) return(NULL)
    consist$pathway <- factor(consist$pathway,
      levels = consist$pathway[order(consist$mean_NES)])

    ggplot(consist, aes(x = mean_NES, y = pathway, fill = direction)) +
      geom_bar(stat = 'identity', alpha = 0.85) +
      geom_vline(xintercept = 0, linewidth = 0.4, color = 'grey50') +
      scale_fill_manual(values = c(Positive = '#B2182B', Negative = '#2166AC')) +
      labs(x = 'Mean NES', y = NULL,
           title = paste0(anc, ' – ', gs_name),
           subtitle = paste0('Sig. in ≥', n_min, '/', length(sample_levels),
                             ' samples (padj<', padj_threshold, ')')) +
      theme_minimal(base_size = 8) +
      theme(panel.grid.major.y = element_blank(),
            axis.text.y = element_text(size = 7),
            legend.position = 'none')
  })
})

consist_flat <- Filter(Negate(is.null), do.call(c, consistent_plots))

if (length(consist_flat) > 0) {
  pdf(file.path(outdir, 'Plots', 'F23b_gsea_consistent_pathways.pdf'),
      width = 5 * length(names(gs_list)),
      height = 5 * length(anchor_genes))
  print(wrap_plots(consist_flat, ncol = length(names(gs_list))))
  dev.off()
  message('F23b written')
}

message('F23 done.')
