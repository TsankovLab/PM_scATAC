# Compute Pearson and Spearman correlation of every SCENIC regulon vs scS-score
# for All cells pooled + per sample

set.seed(1234)

projdir    <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM'
scenic_dir <- file.path(projdir, 'tumor_compartment', 'scrna', 'SCENIC',
                         'vg_5000_mw_tss500bp', 'tumor_programs')
scrna_dir  <- file.path(projdir, 'tumor_compartment', 'scrna')
outdir     <- file.path(projdir, 'git_repo_claude', 'R1_Q3')

source(file.path(projdir, 'git_repo', 'utils', 'load_packages.R'))
source(file.path(projdir, 'git_repo', 'utils', 'ggplot_aestetics.R'))
source(file.path(projdir, 'git_repo', 'utils', 'palettes.R'))
library(patchwork); library(tidyr); library(dplyr)

tumor_prefixes <- c('p11','p12','p13','p811','p826','p848')

# ---- Load AUC matrix -------------------------------------------------------
message('Loading AUC matrix...')
auc <- read.csv(file.path(scenic_dir, 'auc_mtx.csv'), row.names = 1,
                check.names = FALSE)
colnames(auc) <- gsub('\\(\\+\\)|\\(\\-\\)', '', colnames(auc))
auc$sample <- sapply(rownames(auc), function(x) strsplit(x, '_')[[1]][1])
auc <- auc[auc$sample %in% tumor_prefixes, ]
regulon_cols <- colnames(auc)[colnames(auc) != 'sample']

# ---- Compute scS-score per cell --------------------------------------------
message('Loading cNMF20 genes and expression...')
cnmf    <- readRDS(file.path(scrna_dir, 'cnmf_genelist_25_nfeat_5000.rds'))
sarc_genes <- head(cnmf[['cNMF20']], 20)

expr <- read.table(file.path(scenic_dir, 'expr_vg_5000.txt'),
                   header = TRUE, row.names = 1, sep = '\t', check.names = FALSE)
expr <- expr[rownames(expr) %in% rownames(auc), , drop = FALSE]
genes_present <- sarc_genes[sarc_genes %in% colnames(expr)]
scS  <- rowMeans(expr[, genes_present, drop = FALSE], na.rm = TRUE)

common_cells <- intersect(rownames(auc), names(scS))
df <- data.frame(
  scS    = scS[common_cells],
  sample = auc[common_cells, 'sample'],
  auc[common_cells, regulon_cols, drop = FALSE],
  check.names = FALSE
)
df$sample <- factor(df$sample, levels = tumor_prefixes)

# ---- Compute both Pearson and Spearman per regulon per sample --------------
compute_cors <- function(data_sub) {
  sapply(regulon_cols, function(reg) {
    x <- data_sub$scS
    y <- data_sub[[reg]]
    if (var(y, na.rm = TRUE) == 0 | var(x, na.rm = TRUE) == 0) return(c(NA, NA))
    c(pearson  = cor(x, y, method = 'pearson',  use = 'complete.obs'),
      spearman = cor(x, y, method = 'spearman', use = 'complete.obs'))
  })
}

samples_all <- c('All', tumor_prefixes)

cor_list <- lapply(samples_all, function(sam) {
  sub <- if (sam == 'All') df else df[df$sample == sam, ]
  n   <- nrow(sub)
  mat <- compute_cors(sub)   # 2 x n_regulons matrix
  data.frame(
    sample   = sam,
    n_cells  = n,
    regulon  = regulon_cols,
    pearson  = mat['pearson', ],
    spearman = mat['spearman', ],
    row.names = NULL
  )
})
cor_df <- do.call(rbind, cor_list)

# Save full table
write.csv(cor_df, file.path(outdir, 'F14_regulon_scS_pearson_spearman.csv'),
          row.names = FALSE)

# ---- Print SOX9 results clearly --------------------------------------------
sox9 <- cor_df[cor_df$regulon == 'SOX9', ]
cat('\nSOX9 regulon vs scS-score:\n')
cat(sprintf('%-8s %8s %12s %12s\n', 'Sample', 'N cells', 'Pearson r', 'Spearman r'))
cat(strrep('-', 46), '\n')
for (i in seq_len(nrow(sox9))) {
  cat(sprintf('%-8s %8d %12.3f %12.3f\n',
    sox9$sample[i], sox9$n_cells[i], sox9$pearson[i], sox9$spearman[i]))
}

# ---- Print all-sample top regulons by Spearman ----------------------------
all_cors <- cor_df[cor_df$sample == 'All', ]
all_cors <- all_cors[order(-all_cors$spearman), ]
cat('\nTop 15 regulons (Spearman, all cells):\n')
cat(sprintf('%-12s %12s %12s\n', 'Regulon', 'Pearson r', 'Spearman r'))
cat(strrep('-', 38), '\n')
for (i in seq_len(min(15, nrow(all_cors)))) {
  cat(sprintf('%-12s %12.3f %12.3f\n',
    all_cors$regulon[i], all_cors$pearson[i], all_cors$spearman[i]))
}

# ---- Figure: heatmap of Spearman r (all regulons x samples) ---------------
pal_div <- rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 100))

for (method in c('pearson', 'spearman')) {
  cor_wide <- cor_df %>%
    select(sample, regulon, all_of(method)) %>%
    pivot_wider(names_from = sample, values_from = all_of(method))

  # Order regulons by All correlation
  reg_order <- cor_wide$regulon[order(-cor_wide$All, na.last = TRUE)]
  cor_long <- pivot_longer(cor_wide, -regulon, names_to = 'sample', values_to = 'r')
  cor_long$regulon <- factor(cor_long$regulon, levels = reg_order)
  cor_long$sample  <- factor(cor_long$sample, levels = samples_all)

  hm <- ggplot(cor_long, aes(x = sample, y = regulon, fill = r)) +
    geom_tile(color = 'white', linewidth = 0.3) +
    geom_text(data = cor_long[!is.na(cor_long$r) & abs(cor_long$r) > 0.15, ],
      aes(label = sprintf('%.2f', r)), size = 2.2) +
    scale_fill_gradientn(colours = pal_div, limits = c(-1, 1),
      na.value = 'grey85', name = paste0(tools::toTitleCase(method), ' r')) +
    labs(x = NULL, y = NULL,
         title = paste0('SCENIC regulon vs scS-score (', method, ')')) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(size = 7, face = 'italic'),
          axis.text.x = element_text(angle = 30, hjust = 1))

  pdf(file.path(outdir, 'Plots',
      paste0('F14d_regulons_scS_', method, '_heatmap.pdf')), width = 6, height = 9)
  print(hm)
  dev.off()
}

message('Done.')
