# R1_Q3 – F11: Correlation of SOX9 / RUNX2 / SNAI2 with immune infiltration
#
# Strategy:
#   Mesomics  – use pre-computed quanTIseq cell fractions (from SamplesOverview)
#   Bueno/TCGA – derive immune scores from canonical marker gene signatures
#
# Outputs:
#   F11a_immune_cor_heatmap.pdf     – Spearman r heatmap: anchor x immune cell type
#   F11b_immune_cor_scatter_msm.pdf – scatter: anchor vs quanTIseq fractions (Mesomics)
#   F11_immune_cor_table.csv        – full correlation table

set.seed(1234)

projdir  <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
bulk_dir <- file.path(projdir, 'bulkRNA_meso')
outdir   <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))

library(tidyr)
library(dplyr)
library(patchwork)
library(ggpubr)

anchor_genes <- c('SOX9', 'RUNX2', 'SNAI2')
studies      <- c('bueno', 'tcga', 'mesomics')
study_labels <- c(bueno = 'Bueno et al.', tcga = 'TCGA', mesomics = 'Mesomics')

meso_bulk_l    <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies.rds'))
meso_bulk_meta <- readRDS(file.path(bulk_dir, 'bulk_RNA_studies_metadata.rds'))

# ---- Canonical immune gene signatures ------------------------------------
immune_sigs <- list(
  CD8_T_cells    = c('CD8A', 'CD8B', 'GZMB', 'PRF1', 'GZMA'),
  CD4_T_cells    = c('CD4', 'IL7R', 'TCF7'),
  Tregs           = c('FOXP3', 'IL2RA', 'CTLA4', 'IKZF2'),
  B_cells         = c('CD19', 'MS4A1', 'CD79A', 'CD79B'),
  NK_cells        = c('KLRD1', 'NCR1', 'NCAM1', 'GNLY'),
  M1_macrophages  = c('CD68', 'CXCL10', 'TNF', 'IL12B'),
  M2_macrophages  = c('CD163', 'MRC1', 'IL10', 'TGFB1'),
  Monocytes       = c('CD14', 'FCGR3A', 'CSF1R'),
  Dendritic_cells = c('ITGAX', 'CLEC9A', 'FLT3', 'THBD'),
  Total_immune    = c('PTPRC')   # CD45 — pan-leukocyte
)

##############################################################################
# 1. Compute marker-based immune scores for Bueno and TCGA
##############################################################################
compute_immune_scores <- function(bulk_mat) {
  scores <- lapply(names(immune_sigs), function(ct) {
    genes <- immune_sigs[[ct]][immune_sigs[[ct]] %in% rownames(bulk_mat)]
    if (length(genes) == 0) return(rep(NA, ncol(bulk_mat)))
    colMeans(bulk_mat[genes, , drop = FALSE], na.rm = TRUE)
  })
  scores_df <- as.data.frame(do.call(cbind, scores))
  colnames(scores_df) <- names(immune_sigs)
  scores_df
}

##############################################################################
# 2. Mesomics: use quanTIseq fractions + marker-based scores
##############################################################################
msm_meta <- meso_bulk_meta[['mesomics']]

quantiseq_cols <- c(
  'B.cells.quanTIseq', 'Macrophages.M1.quanTIseq', 'Macrophages.M2.quanTIseq',
  'Monocytes.quanTIseq', 'Neutrophils.quanTIseq', 'NK.cells.quanTIseq',
  'T.cells.CD4.quanTIseq', 'T.cells.CD8.quanTIseq', 'Tregs.quanTIseq',
  'Dendritic.cells.quanTIseq', 'Total.Immune.cells.quanTIseq'
)
# handle space in column name
colnames(msm_meta) <- gsub(' ', '.', colnames(msm_meta))
quantiseq_cols_present <- quantiseq_cols[quantiseq_cols %in% colnames(msm_meta)]

msm_immune <- msm_meta[, quantiseq_cols_present, drop = FALSE]
colnames(msm_immune) <- gsub('\\.quanTIseq', '', colnames(msm_immune))
colnames(msm_immune) <- gsub('\\.', '_', colnames(msm_immune))
msm_immune <- as.data.frame(lapply(msm_immune, as.numeric))

##############################################################################
# 3. Compute Spearman correlations per study
##############################################################################
cor_all <- list()

for (study in studies) {
  bulk <- meso_bulk_l[[study]]

  if (study == 'mesomics') {
    immune_df <- msm_immune
  } else {
    immune_df <- compute_immune_scores(bulk)
  }

  for (anchor in anchor_genes) {
    if (!anchor %in% rownames(bulk)) next
    anchor_expr <- as.numeric(bulk[anchor, ])

    cor_vec <- sapply(colnames(immune_df), function(ct) {
      vals <- as.numeric(immune_df[, ct])
      if (sum(!is.na(vals)) < 10) return(NA)
      cor(anchor_expr, vals, method = 'spearman', use = 'complete.obs')
    })

    cor_all[[paste(anchor, study, sep = '|')]] <- data.frame(
      anchor     = anchor,
      cell_type  = names(cor_vec),
      rho        = cor_vec,
      cohort     = study_labels[study],
      data_source = ifelse(study == 'mesomics', 'quanTIseq', 'Marker genes'),
      stringsAsFactors = FALSE
    )
  }
}

cor_df <- do.call(rbind, cor_all)
rownames(cor_df) <- NULL
cor_df <- cor_df[!is.na(cor_df$rho), ]

write.csv(cor_df, file.path(outdir, 'F11_immune_cor_table.csv'), row.names = FALSE)

##############################################################################
# F11a – Heatmap: Spearman r (anchor x immune cell type) per cohort
##############################################################################
cor_df$anchor  <- factor(cor_df$anchor,  levels = anchor_genes)
cor_df$cohort  <- factor(cor_df$cohort,  levels = unname(study_labels))

# Order cell types by mean absolute rho
ct_order <- cor_df %>%
  group_by(cell_type) %>%
  summarise(mean_abs = mean(abs(rho), na.rm = TRUE)) %>%
  arrange(mean_abs) %>%
  pull(cell_type)
cor_df$cell_type <- factor(cor_df$cell_type, levels = ct_order)

f11a <- ggplot(cor_df, aes(x = anchor, y = cell_type, fill = rho)) +
  geom_tile(color = 'white', linewidth = 0.5) +
  geom_text(aes(label = ifelse(abs(rho) > 0.2, sprintf('%.2f', rho), '')),
    size = 2.5, color = 'black') +
  scale_fill_gradientn(
    colours = rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100)),
    limits  = c(-1, 1), name = 'Spearman r') +
  facet_wrap(~ cohort, ncol = length(studies)) +
  labs(x = NULL, y = NULL,
       title = 'Correlation of SOX9 / RUNX2 / SNAI2 with immune infiltration') +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = 'bold'),
        panel.grid = element_blank(),
        axis.text.x = element_text(face = 'italic'))

pdf(file.path(outdir, 'Plots', 'F11a_immune_cor_heatmap.pdf'), width = 11, height = 5)
print(f11a)
dev.off()

##############################################################################
# F11b – Scatter: anchor vs quanTIseq fractions in Mesomics
##############################################################################
bulk_msm  <- meso_bulk_l[['mesomics']]
msm_df    <- cbind(t(bulk_msm[anchor_genes[anchor_genes %in% rownames(bulk_msm)], ]),
                   msm_immune)
msm_df    <- as.data.frame(msm_df)

# Pick the most immunologically informative cell types
key_cts <- intersect(c('T_cells_CD8', 'T_cells_CD4', 'Macrophages_M2',
                        'NK_cells', 'B_cells', 'Tregs', 'Total_Immune_cells'),
                     colnames(msm_df))

scatter_plots <- list()
for (anchor in anchor_genes) {
  if (!anchor %in% colnames(msm_df)) next
  for (ct in key_cts) {
    scatter_plots[[paste(anchor, ct)]] <- ggplot(msm_df,
        aes_string(x = anchor, y = ct)) +
      geom_point(shape = 21, fill = 'steelblue', size = 1.5, alpha = 0.7, stroke = 0.3) +
      geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.5) +
      stat_cor(method = 'spearman', size = 2.5, label.x.npc = 'left', label.y.npc = 'top') +
      labs(x = paste0(anchor, ' (log2)'),
           y = gsub('_', ' ', ct),
           title = paste0(anchor, ' vs ', gsub('_', ' ', ct))) +
      gtheme
  }
}

n_ct <- length(key_cts)
pdf(file.path(outdir, 'Plots', 'F11b_immune_cor_scatter_mesomics.pdf'),
    width = 3.5 * n_ct, height = 3.5 * length(anchor_genes))
print(wrap_plots(scatter_plots, ncol = n_ct, byrow = TRUE))
dev.off()

##############################################################################
# F11c/d – Scatter: anchor vs marker-based immune scores (Bueno + TCGA)
##############################################################################
make_scatter_cohort <- function(study, study_label) {
  bulk        <- meso_bulk_l[[study]]
  immune_df   <- compute_immune_scores(bulk)
  # Keep only cell types with non-trivial range (max > 0)
  immune_df   <- immune_df[, apply(immune_df, 2, function(x) max(x, na.rm=TRUE) > 0),
                            drop = FALSE]
  key_cts_mk  <- colnames(immune_df)

  df <- cbind(
    t(bulk[anchor_genes[anchor_genes %in% rownames(bulk)], , drop = FALSE]),
    immune_df
  )
  df <- as.data.frame(df)

  plots <- list()
  for (anchor in anchor_genes) {
    if (!anchor %in% colnames(df)) next
    for (ct in key_cts_mk) {
      plots[[paste(anchor, ct)]] <- ggplot(df, aes_string(x = anchor, y = ct)) +
        geom_point(shape = 21, fill = 'tomato3', size = 1.2, alpha = 0.6, stroke = 0.3) +
        geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 0.5) +
        stat_cor(method = 'spearman', size = 2.5,
                 label.x.npc = 'left', label.y.npc = 'top') +
        labs(x = paste0(anchor, ' (log2)'),
             y = paste0(gsub('_', ' ', ct), '\n(marker genes)'),
             title = paste0(anchor, ' vs ', gsub('_', ' ', ct))) +
        gtheme
    }
  }
  list(plots = plots, n_ct = length(key_cts_mk))
}

bueno_sc  <- make_scatter_cohort('bueno', 'Bueno et al.')
tcga_sc   <- make_scatter_cohort('tcga',  'TCGA')

pdf(file.path(outdir, 'Plots', 'F11c_immune_cor_scatter_bueno.pdf'),
    width = 3.5 * bueno_sc$n_ct, height = 3.5 * length(anchor_genes))
print(wrap_plots(bueno_sc$plots, ncol = bueno_sc$n_ct, byrow = TRUE))
dev.off()

pdf(file.path(outdir, 'Plots', 'F11d_immune_cor_scatter_tcga.pdf'),
    width = 3.5 * tcga_sc$n_ct, height = 3.5 * length(anchor_genes))
print(wrap_plots(tcga_sc$plots, ncol = tcga_sc$n_ct, byrow = TRUE))
dev.off()

message('F11 done.')
