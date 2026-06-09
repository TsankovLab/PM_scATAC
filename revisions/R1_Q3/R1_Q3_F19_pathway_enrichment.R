# R1_Q3 – F19: Pathway enrichment on top co-expressed genes with SOX9/RUNX2/SNAI2
#
# Uses:
#   - Metacell co-expression (F17_coexpr_all_genes_table.csv) — scRNA-seq ranked list
#   - Bulk RNA co-expression (F18_metacell_bulk_cor_table.csv) — bulk ranked list
#
# Methods:
#   1. GSEA (fgsea) — full ranked gene list ordered by Spearman r
#   2. Fisher's exact test — top 200 positively correlated genes vs gene sets
#
# Gene sets: Hallmark (h.all.v7.4) + GO Biological Process (c5.bp.v7.1)
#
# Outputs:
#   F19a_gsea_scrna_barplot.pdf         – GSEA NES per anchor (scRNA)
#   F19b_gsea_bulk_barplot.pdf          – GSEA NES per anchor (bulk)
#   F19c_fisher_scrna_dotplot.pdf       – Fisher enrichment (scRNA top 200)
#   F19d_fisher_bulk_dotplot.pdf        – Fisher enrichment (bulk top 200)
#   F19e_scrna_vs_bulk_NES_scatter.pdf  – NES comparison scRNA vs bulk
#   F19_gsea_results.csv                – full GSEA table
#   F19_fisher_results.csv              – full Fisher table

set.seed(1234)

projdir  <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
files_dir <- file.path(projdir, 'git_repo', 'files')
outdir   <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))
library(fgsea); library(patchwork); library(tidyr); library(dplyr)

anchor_genes   <- c('SOX9', 'RUNX2', 'SNAI2')
top_fisher     <- 200   # top N genes for Fisher test
padj_threshold <- 0.05

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

# Clean names
names(hallmark) <- gsub('HALLMARK_', '', names(hallmark))
names(hallmark) <- gsub('_', ' ', names(hallmark))
names(gobp)     <- gsub('GOBP_', '', names(gobp))
names(gobp)     <- gsub('_', ' ', names(gobp))

gs_list <- list(Hallmark = hallmark, GOBP = gobp)
message('Hallmark: ', length(hallmark), ' | GOBP: ', length(gobp))

##############################################################################
# Load co-expression ranked lists
##############################################################################
message('Loading metacell co-expression (F17)...')
f17 <- read.csv(file.path(outdir, 'F17_coexpr_all_genes_table.csv'),
                stringsAsFactors = FALSE)

# Median rho across samples per anchor x gene (scRNA ranked list)
scrna_ranked <- f17 %>%
  group_by(anchor, gene) %>%
  summarise(rho = median(rho, na.rm = TRUE), .groups = 'drop')

message('Loading bulk co-expression (F18)...')
f18 <- read.csv(file.path(outdir, 'F18_metacell_bulk_cor_table.csv'),
                stringsAsFactors = FALSE)

# Median bulk rho across cohorts per anchor x gene
bulk_ranked <- f18 %>%
  select(anchor, gene, bulk_rho_median) %>%
  rename(rho = bulk_rho_median) %>%
  distinct()

##############################################################################
# Helper: create named ranked vector for fgsea
##############################################################################
make_ranked <- function(df, anc) {
  sub  <- df[df$anchor == anc & !is.na(df$rho), ]
  sub  <- sub[!duplicated(sub$gene), ]
  vec  <- setNames(sub$rho, sub$gene)
  sort(vec, decreasing = TRUE)
}

##############################################################################
# Helper: run fgsea on one ranked vector + one gene set collection
##############################################################################
run_gsea <- function(ranked_vec, pathways, label) {
  message('  GSEA: ', label, ' (n=', length(ranked_vec), ')')
  res <- tryCatch(
    fgseaMultilevel(pathways, ranked_vec, minSize = 10, maxSize = 500,
                    nPermSimple = 10000, eps = 0),
    error = function(e) { message('  fgsea error: ', e$message); NULL }
  )
  if (is.null(res)) return(NULL)
  res$label <- label
  as.data.frame(res[, c('pathway','pval','padj','NES','size','label')])
}

##############################################################################
# Helper: Fisher's exact test
##############################################################################
run_fisher <- function(top_genes, background, pathways, label) {
  message('  Fisher: ', label, ' (top=', length(top_genes), ')')
  bg_n     <- length(background)
  top_n    <- length(top_genes)
  do.call(rbind, lapply(names(pathways), function(pw) {
    pw_genes  <- intersect(pathways[[pw]], background)
    pw_n      <- length(pw_genes)
    if (pw_n < 5) return(NULL)
    overlap   <- length(intersect(top_genes, pw_genes))
    # Contingency: overlap | top_not_pw | pw_not_top | neither
    mat <- matrix(c(overlap, top_n - overlap,
                    pw_n - overlap, bg_n - top_n - pw_n + overlap),
                  nrow = 2)
    ft  <- fisher.test(mat, alternative = 'greater')
    data.frame(pathway = pw, overlap = overlap, pw_size = pw_n,
               pval = ft$p.value, OR = ft$estimate, label = label,
               stringsAsFactors = FALSE)
  }))
}

##############################################################################
# Run enrichments
##############################################################################
gsea_results   <- list()
fisher_results <- list()

for (anc in anchor_genes) {
  message('\n=== ', anc, ' ===')

  scrna_vec <- make_ranked(scrna_ranked, anc)
  bulk_vec  <- make_ranked(bulk_ranked, anc)
  background <- union(names(scrna_vec), names(bulk_vec))

  top_scrna <- names(head(scrna_vec[scrna_vec > 0], top_fisher))
  top_bulk  <- names(head(bulk_vec[bulk_vec > 0], top_fisher))

  for (gs_name in names(gs_list)) {
    gset <- gs_list[[gs_name]]

    # GSEA
    gsea_results[[paste(anc, gs_name, 'scRNA')]] <-
      cbind(anchor = anc, gs = gs_name, modality = 'scRNA',
            run_gsea(scrna_vec, gset, paste(anc, gs_name, 'scRNA')))
    gsea_results[[paste(anc, gs_name, 'bulk')]] <-
      cbind(anchor = anc, gs = gs_name, modality = 'bulk',
            run_gsea(bulk_vec,  gset, paste(anc, gs_name, 'bulk')))

    # Fisher
    fisher_results[[paste(anc, gs_name, 'scRNA')]] <-
      cbind(anchor = anc, gs = gs_name, modality = 'scRNA',
            run_fisher(top_scrna, background, gset,
                       paste(anc, gs_name, 'scRNA')))
    fisher_results[[paste(anc, gs_name, 'bulk')]] <-
      cbind(anchor = anc, gs = gs_name, modality = 'bulk',
            run_fisher(top_bulk,  background, gset,
                       paste(anc, gs_name, 'bulk')))
  }
}

gsea_df   <- do.call(rbind, Filter(Negate(is.null), gsea_results))
fisher_df <- do.call(rbind, Filter(Negate(is.null), fisher_results))
gsea_df$padj   <- as.numeric(gsea_df$padj)
gsea_df$NES    <- as.numeric(gsea_df$NES)
fisher_df$pval <- as.numeric(fisher_df$pval)
fisher_df$padj <- p.adjust(fisher_df$pval, method = 'BH')

write.csv(gsea_df,   file.path(outdir, 'F19_gsea_results.csv'),   row.names = FALSE)
write.csv(fisher_df, file.path(outdir, 'F19_fisher_results.csv'), row.names = FALSE)

##############################################################################
# Plotting helpers
##############################################################################
pal_div <- rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100))
modality_colors <- c(scRNA = '#E66101', bulk = '#5E3C99')

plot_gsea_barplot <- function(gsea_sub, title, top_n = 15) {
  sig <- gsea_sub[!is.na(gsea_sub$padj) & gsea_sub$padj < padj_threshold, ]
  if (nrow(sig) == 0) {
    sig <- gsea_sub[!is.na(gsea_sub$pval), ]
    sig <- sig[order(sig$pval), ][seq_len(min(top_n, nrow(sig))), ]
    sig$padj_label <- 'ns'
  } else {
    sig <- rbind(
      head(sig[order(-sig$NES), ], top_n),
      head(sig[order(sig$NES),  ], top_n)
    )
    sig <- sig[!duplicated(sig$pathway), ]
    sig$padj_label <- ifelse(sig$padj < 0.001, '***',
                      ifelse(sig$padj < 0.01,  '**', '*'))
  }
  sig$pathway <- factor(sig$pathway, levels = sig$pathway[order(sig$NES)])
  ggplot(sig, aes(x = NES, y = pathway, fill = NES)) +
    geom_bar(stat = 'identity', alpha = 0.85) +
    geom_text(aes(label = padj_label, x = ifelse(NES > 0, NES + 0.05, NES - 0.05)),
              hjust = ifelse(sig$NES > 0, 0, 1), size = 2.5) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = 'grey50') +
    scale_fill_gradientn(colours = pal_div, limits = c(-3, 3),
                         oob = scales::squish, name = 'NES') +
    labs(x = 'NES', y = NULL, title = title) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 7))
}

plot_fisher_dotplot <- function(fisher_sub, title, top_n = 15) {
  sig <- fisher_sub[!is.na(fisher_sub$padj) & fisher_sub$padj < padj_threshold, ]
  if (nrow(sig) == 0)
    sig <- fisher_sub[order(fisher_sub$pval), ][seq_len(min(top_n, nrow(sig))), ]
  sig <- head(sig[order(sig$pval), ], top_n)
  if (nrow(sig) == 0) return(ggplot() + labs(title = paste0(title, ' (no sig.)')) + gtheme)
  sig$pathway <- factor(sig$pathway, levels = sig$pathway[order(sig$OR)])
  sig$log10p  <- -log10(sig$pval)
  ggplot(sig, aes(x = log10p, y = pathway, size = overlap, color = log10p)) +
    geom_point(alpha = 0.8) +
    scale_color_gradientn(colours = c('grey80','#B2182B'), name = '-log10(p)') +
    scale_size_continuous(range = c(2, 7), name = 'Overlap') +
    geom_vline(xintercept = -log10(0.05), linetype = 'dashed',
               color = 'grey50', linewidth = 0.4) +
    labs(x = '-log10(p-value)', y = NULL, title = title) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 7))
}

##############################################################################
# F19a/b – GSEA barplots per anchor x modality
##############################################################################
for (gs_name in names(gs_list)) {
  for (mod in c('scRNA', 'bulk')) {
    sub    <- gsea_df[gsea_df$gs == gs_name & gsea_df$modality == mod, ]
    plots  <- lapply(anchor_genes, function(anc)
      plot_gsea_barplot(sub[sub$anchor == anc, ],
                        paste0(anc, ' – ', mod, ' (', gs_name, ')')))
    fname  <- paste0('F19_gsea_', gs_name, '_', mod, '_barplot.pdf')
    pdf(file.path(outdir, 'Plots', fname), width = 5 * length(anchor_genes), height = 7)
    print(wrap_plots(plots, ncol = length(anchor_genes)))
    dev.off()
    message('Wrote ', fname)
  }
}

##############################################################################
# F19c/d – Fisher dotplots per anchor x modality
##############################################################################
for (gs_name in names(gs_list)) {
  for (mod in c('scRNA', 'bulk')) {
    sub   <- fisher_df[fisher_df$gs == gs_name & fisher_df$modality == mod, ]
    plots <- lapply(anchor_genes, function(anc)
      plot_fisher_dotplot(sub[sub$anchor == anc, ],
                          paste0(anc, ' – ', mod, ' (', gs_name, ')')))
    fname <- paste0('F19_fisher_', gs_name, '_', mod, '_dotplot.pdf')
    pdf(file.path(outdir, 'Plots', fname), width = 5 * length(anchor_genes), height = 6)
    print(wrap_plots(plots, ncol = length(anchor_genes)))
    dev.off()
    message('Wrote ', fname)
  }
}

##############################################################################
# F19e – NES comparison: scRNA vs bulk (per pathway per anchor)
##############################################################################
for (gs_name in names(gs_list)) {
  gsea_wide <- gsea_df[gsea_df$gs == gs_name, ] %>%
    select(anchor, modality, pathway, NES, padj) %>%
    pivot_wider(names_from = modality, values_from = c(NES, padj)) %>%
    filter(!is.na(NES_scRNA) & !is.na(NES_bulk))

  sig_either <- gsea_wide[
    (!is.na(gsea_wide$padj_scRNA) & gsea_wide$padj_scRNA < padj_threshold) |
    (!is.na(gsea_wide$padj_bulk)  & gsea_wide$padj_bulk  < padj_threshold), ]

  nes_plots <- lapply(anchor_genes, function(anc) {
    sub <- gsea_wide[gsea_wide$anchor == anc, ]
    if (nrow(sub) < 3) return(NULL)
    sig_sub <- sig_either[sig_either$anchor == anc, ]
    r_val   <- cor(sub$NES_scRNA, sub$NES_bulk, use = 'complete.obs', method = 'pearson')
    p_sub   <- ggplot(sub, aes(x = NES_scRNA, y = NES_bulk)) +
      geom_point(alpha = 0.25, size = 0.8, color = 'grey50') +
      geom_point(data = sig_sub, aes(color = NES_scRNA > 0),
                 size = 1.5, alpha = 0.8) +
      ggrepel::geom_text_repel(data = head(sig_sub[order(-abs(sig_sub$NES_scRNA +
                                                               sig_sub$NES_bulk)), ], 15),
        aes(label = pathway), size = 2, max.overlaps = 20,
        segment.size = 0.2) +
      geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.5) +
      geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey70', linewidth = 0.3) +
      geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey70', linewidth = 0.3) +
      scale_color_manual(values = c(`TRUE` = '#B2182B', `FALSE` = '#2166AC')) +
      annotate('text', x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3,
               label = sprintf('r=%.3f', r_val)) +
      labs(x = 'NES (scRNA metacells)', y = 'NES (bulk RNA)',
           title = paste0(anc, ' – ', gs_name)) +
      gtheme + NoLegend()
    p_sub
  })
  nes_plots <- Filter(Negate(is.null), nes_plots)

  pdf(file.path(outdir, 'Plots',
      paste0('F19e_NES_scRNA_vs_bulk_', gs_name, '.pdf')),
      width = 5 * length(nes_plots), height = 5)
  print(wrap_plots(nes_plots, ncol = length(nes_plots)))
  dev.off()
  message('Wrote F19e NES scatter for ', gs_name)
}

message('\nF19 done.')
