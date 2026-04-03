# Differential Peaks in SOX9 high vs SOX9 low cells in P23 ####
archp = addClusters (input = archp, resolution = 1.5,
  reducedDims = "IterativeLSI", name = 'Clusters2',
  maxClusters = 100,
  force = TRUE)

archp = addClusters (input = archp, resolution = 20,
  reducedDims = "IterativeLSI", name = 'Clusters3',
  maxClusters = 100,
  force = TRUE)

force=T
metaGroupName = 'Clusters2'
if (!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force)
  {
  DAG_list = getMarkerFeatures (
    ArchRProj = archp, 
    testMethod = "wilcoxon",
          useGroups = "C9",
          bgdGroups = "C4",
    k=100,
    binarize = FALSE,
    useMatrix = "GeneScoreMatrix",
    groupBy = metaGroupName
  #  useSeqnames="z"
  )
  saveRDS (DAG_list, paste0('DAG_',metaGroupName,'.rds'))
  } else {
  DAG_list = readRDS (paste0('DAG_',metaGroupName,'.rds'))
  }

pdf(file.path('Plots','P23_Clusters2_MA_plot.pdf'), width=5,height=5)
dma <- markerPlot (seMarker = DAG_list, name = 'C9', cutOff = "FDR <= 0.01", plotAs = "MA")
dev.off()


# check SOX9 deviation ####
archp_P23 = archp[archp$Sample %in% c('P23')]
archp_P23 = archp_P23[archp_P23@embeddings$UMAP[[1]][[2]] < 0]
tf_markers = c('SOX9','SNAI2','RUNX2','SOX6','JUND','FOS')
markerMotifs = getFeatures (archp_P23, select = paste(tf_markers, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs = grep ("z:", markerMotifs, value = TRUE)
archp_P23 = addImputeWeights (archp_P23)

pdf()
TF_p = plotEmbedding(
    ArchRProj = archp_P23, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    pal = rev(palette_deviation),
    imputeWeights = getImputeWeights(archp_P23)
)

TF_p2 = plotEmbedding(
    ArchRProj = archp_P23, 
    colorBy = "GeneScoreMatrix", 
    name = c('SOX9','SNAI2','RUNX2','RUNX1','NFIC','SOX6'), 
    embedding = "UMAP",
    pal = rev(palette_genescore),
    imputeWeights = getImputeWeights(archp_P23)
)

TF_p3 = plotEmbedding(
    ArchRProj = archp_P23, 
    colorBy = "GeneScoreMatrix", 
    name = c('SOX9','SNAI2','RUNX2','RUNX1','NFIC','SOX6'), 
    embedding = "UMAP",
    pal = rev(palette_genescore),
    imputeWeights = NULL
)

umap_p1 = plotEmbedding (ArchRProj = archp_P23, labelMeans = T, 
  colorBy = "cellColData", name = "Clusters2", 
  #pal = palette_sample,
   embedding = "UMAP")
umap_p2 = plotEmbedding (ArchRProj = archp_P23, labelMeans = T, 
  colorBy = "cellColData", name = "Clusters3", 
  #pal = palette_sample,
   embedding = "UMAP")

umap_p3 = plotEmbedding (ArchRProj = archp_P23, labelMeans = F, 
  colorBy = "cellColData", name = "cNMF20", 
  pal = rev(palette_genescore),
  imputeWeights = getImputeWeights(archp_P23),
  #pal = palette_sample,
   embedding = "UMAP")
dev.off()

pdf (file.path('Plots','SOX9_dev_P23_fplot.pdf'), width = 14, height = 12)
wrap_plots (TF_p)
wrap_plots (TF_p2)
wrap_plots (TF_p3)
wrap_plots (umap_p1, umap_p2, umap_p3)
dev.off()



# Run GSEA enrichment analysis #### 
#!!! To run this analysis load only ArcHR and clusterprofiler packages !!!!
library (fgsea)    
options(warn = 0)
ps = getPeakSet (archp)
metaGroupName = 'Clusters2'
gmt_annotations = c(
'h.all.v7.4.symbol.gmt',#,
'c5.bp.v7.1.symbol.gmt',
'c3.tft.v7.1.symbol.gmt'
)

gmt.file = paste0 ('../../git_repo/files/c5.bp.v7.1.symbol.gmt')
gmt.file = paste0 ('../../git_repo/files/h.all.v7.4.symbols.gmt')

pathways = clusterProfiler::read.gmt (gmt.file)

pathways = split(pathways$gene, pathways$term)
#pathways = gmtPathways (gmt.file)
DAG_list = readRDS (paste0('DAG_',metaGroupName,'.rds'))
DAG_res = do.call (cbind, (assays(DAG_list)))
colnames (DAG_res) = names(assays(DAG_list))
rownames(DAG_res) = rowData(DAG_list)$name

ranked_genes = setNames (-log10(DAG_res$FDR) * sign(DAG_res$Log2FC), rownames(DAG_res))
ranked_genes = na.omit (ranked_genes)
#peak_genes = peak_genes[!duplicated(names(peak_genes))]
#names (peak_genes) = 
ranked_genes = ranked_genes[!duplicated(names(ranked_genes))]
ranked_genes = ranked_genes[!is.na(names(ranked_genes))]
library (BiocParallel)
BiocParallel::register(BiocParallel::SerialParam())

### fgsea throws a BiocParallel error when I load all packages including clusterProfiler...try avoiding loading packages except ArchR and fgsea
#peak_genes2 = setNames(order(peak_genes), names(peak_genes))
fgseaRes = fgsea (pathways, 
          ranked_genes,#, 
          minSize=15, 
        #  scoreType='pos',
          maxSize=1500,
          nproc=1,
          nPermSimple=100000,
          BPPARAM = NULL
          )
pvalAdjTrheshold = 0.05
top_pathways=5
#fgseaRes$padj = fgseaRes$pval
fgseaResAll_dp = dotGSEA (
  list(fgseaRes), 
  padj_threshold = pvalAdjTrheshold, 
  type = 'fgsea',
  top_pathways = top_pathways,
  cluster_rows=F,
  cluster_cols=F)

pdf (file.path ('Plots','fgsea_dotplot.pdf'), width=7, height=3)
fgseaResAll_dp
dev.off()

### Compare DEG SOX9 high vs SOX9 low to list of genes from Yang et al paper ####
rankby = 'pval_signedLFC'

sox9_gs = read.csv ('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/gene_sets/Sox9_targets_mouse.csv')
D0 = head(sox9_gs$symbol[order(-sox9_gs$D0)],500)
W1 = head(sox9_gs$symbol[order(-sox9_gs$W1)],500)
W2 = head(sox9_gs$symbol[order(-sox9_gs$W2)],500)
W6 = head(sox9_gs$symbol[order(-sox9_gs$W6)],500)
W12 = head(sox9_gs$symbol[order(-sox9_gs$W12)],500)
#sox9_gs_class = split (sox9_gs$symbol, sox9_gs$class)
sox9_gs_class = list (D0 = D0, W1 = W1, W2 = W2, W6 = W6, W12=W12)
pathways = sox9_gs_class
pathways = lapply (pathways, function(x) mouseMan (x))
pathways = lapply (pathways, function(x) x$HGNC.symbol)

fgseaRes_fuchs = fgsea (pathways, 
          ranked_genes,#, 
          minSize=15, 
        #  scoreType='pos',
          maxSize=5000,
          nproc=1,
          nPermSimple=100000,
          BPPARAM = NULL
          )

pvalAdjTrheshold = 0.05
top_pathways=5
fgseaRes_fuchs$padj = fgseaRes_fuchs$padj
fgseaResAll_dp = dotGSEA (
  list(fgseaRes_fuchs), 
  padj_threshold = pvalAdjTrheshold, 
  type = 'fgsea',
  top_pathways = top_pathways,
  cluster_rows=F,
  cluster_cols=F)

pdf (file.path ('Plots','fgsea_Yang_dotplot.pdf'), width=4, height=3)
fgseaResAll_dp
dev.off()


#### Make boxplots of gene modules
gs_feat = getFeatures (archp, 'GeneScoreMatrix')
pathways = lapply (pathways, function(x) x[x %in% gs_feat])
archp = addModuleScore (
  ArchRProj = archp,
  useMatrix = 'GeneScoreMatrix',
  name = "Sox9_mouse",
  features = pathways,
  nBin = 25,
  nBgd = 100,
  seed = 1,
  threads = getArchRThreads(),
  logFile = createLogFile("addModuleScore")
)

ccomp = as.data.frame (archp@cellColData)[, c(paste0('Sox9_mouse.',names (sox9_gs_class)),'Clusters2')]
ccomp = ccomp[ccomp$Clusters2 %in% c('C4','C9'),]
library (tidyr)
ccomp = gather (ccomp, module, score, 1:(ncol(ccomp)-1))
ccomp$module = factor (ccomp$module, levels = paste0('Sox9_mouse.',c('D0','W1','W2','W6','W12')))
ccomp$Clusters2 = factor (ccomp$Clusters2, levels = c('C4','C9'))
gp = ggplot (ccomp, aes (x = Clusters2, y = score), alpha=.5) + 
  geom_violin (aes (fill = Clusters2), 
    trim=TRUE,size=2,
    width=1,
    scale='width',
    linewidth = .2, alpha=0.7) + 
  geom_boxplot (aes (fill = Clusters2),
    linewidth = .2,
    width=.4,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.7) + 
  gtheme + 
  scale_fill_manual (values = c(C4 = 'blue',C9='purple')) + 
  facet_wrap (~module, ncol = 5, scale='free')

stat.test = gp$data %>% group_by (module) %>%
wilcox_test (reformulate ('Clusters2', 'score')) %>%
adjust_pvalue (method = "fdr") %>%
add_significance ()
stat.test = stat.test %>% add_xy_position (x = 'Clusters2', step.increase=.4)
gp = gp + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
bracket.nudge.y = 0, hide.ns = T,
label = "p.adj.signif") + NoLegend()

pdf (file.path ('Plots','Yang_SOX9_module_scores_boxplot.pdf'), width = 5, height=2.4)
gp
dev.off()






# Import chromBPnet finemo motifs ####
#archp_P1 = archp[archp$Clusters %in% c('C2') & archp$Sample == 'P1']
library ('universalmotif')

ps = getPeakSet (archp)

chromBPdir = 'chromBPnet'

chrombpnet_counts = list()
celltypes = c('SOX9_low_P23','SOX9_high_P23')#,'SOX9_high_P1')
#celltypes = c('SOX9_high_P1')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_counts[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_counts_to_genome_browser.tsv')))
  gr = makeGRangesFromDataFrame (chrombpnet_counts[[celltype]], keep.extra.columns=T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  chrombpnet_counts[[celltype]]$peak_type = ps$peakType[findOverlaps(gr, ps, select='first')]
  }


chrombpnet_profile = list()
#celltypes = c('Mesothelium','Malignant')
for (celltype in celltypes)
  {
  message (paste0('reading finemo output for ', celltype))  
  chrombpnet_profile[[celltype]] = read.table (file.path (chromBPdir, celltype,'no_bias_model',paste0(celltype, '_finemo_profile_to_genome_browser.tsv')))
  gr = makeGRangesFromDataFrame (chrombpnet_profile[[celltype]], keep.extra.columns=T, seqnames.field = 'V1', start.field = 'V2', end.field = 'V3')
  chrombpnet_profile[[celltype]]$peak_type = ps$peakType[findOverlaps(gr, ps, select='first')]
  }

chrombpnet_profile = lapply (chrombpnet_profile, function(x) x[x$V5 != 'NaN_NaN_NaN',])


top_n <- 5
n <- length(chrombpnet_counts )

bp_list <- lapply(seq_len(n), function(i) {
  tbl <- table(chrombpnet_counts [[i]]$V4)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_counts[[i]]$V5[chrombpnet_counts [[i]]$V4 == tf][1]
  })
  
  data.frame(
    Freq = proportions(top_tbl),
    TF   = tf_names,
    direction = directions,
    type = rep(celltypes[[i]], length(top_tbl))
  )
})

bp_df <- do.call(rbind, bp_list)

# Make neg values negative
bp_df <- bp_df %>%
  mutate(Freq = ifelse(direction == "neg", -Freq.Freq, Freq.Freq))

# Create custom ordering per type
bp_df <- bp_df %>%
  group_by(type, direction) %>%
  mutate(
    TF_order = ifelse(direction == "pos",
                      rank(-Freq, ties.method = "first"),  # descending
                      rank(Freq, ties.method = "first"))   # ascending for neg (opposite)
  ) %>%
  ungroup()

# Build a combined factor: ensures pos stack from bottom up, neg from top down
bp_df <- bp_df %>%
  arrange(type, direction, TF_order)

bp_df$TF_id <- paste(bp_df$TF, bp_df$type, sep = "_")
bp_df$TF_id <- factor(bp_df$TF_id, levels = unique(bp_df$TF_id))
bp_df$type = factor (bp_df$type, levels = celltypes)
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF)) ) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle(paste0("Top ",top_n," TFs (counts head)")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

  

pdf (file.path ('Plots', 'TF_abundance_counts_barplot.pdf'),6,width=6.5)
bp
dev.off()


top_n <- 5
n <- length(chrombpnet_profile)

bp_list <- lapply(seq_len(n), function(i) {
  tbl <- table(chrombpnet_profile[[i]]$V4)
  tbl_sorted <- sort(tbl, decreasing = TRUE)
  top_tbl <- head(tbl_sorted, top_n)
  
  tf_names <- names(top_tbl)
  directions <- sapply(tf_names, function(tf) {
    chrombpnet_profile[[i]]$V5[chrombpnet_profile[[i]]$V4 == tf][1]
  })
  
  data.frame(
    Freq = proportions(top_tbl),
    TF   = tf_names,
    direction = directions,
    type = rep(celltypes[[i]], length(top_tbl))
  )
})

bp_df <- do.call(rbind, bp_list)

# Make neg values negative
bp_df <- bp_df %>%
  mutate(Freq = ifelse(direction == "neg", -Freq.Freq, Freq.Freq))

# Create custom ordering per type
bp_df <- bp_df %>%
  group_by(type, direction) %>%
  mutate(
    TF_order = ifelse(direction == "pos",
                      rank(-Freq, ties.method = "first"),  # descending
                      rank(Freq, ties.method = "first"))   # ascending for neg (opposite)
  ) %>%
  ungroup()

# Build a combined factor: ensures pos stack from bottom up, neg from top down
bp_df <- bp_df %>%
  arrange(type, direction, TF_order)

bp_df$TF_id <- paste(bp_df$TF, bp_df$type, sep = "_")
bp_df$TF_id <- factor(bp_df$TF_id, levels = unique(bp_df$TF_id))
bp_df$type = factor (bp_df$type, levels = celltypes)
# Plot stacked bars
bp <- ggplot(bp_df, aes(x = type, y = Freq, fill = TF_id)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = paletteer_d("palettesForR::LaTeX", length(bp_df$TF))) +
  theme_minimal(base_size = 14) +
  ylab("Proportion of counts") +
  xlab("Cell type") +
  ggtitle(paste0("Top ",top_n," TFs (profile head)")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1)

  

pdf (file.path ('Plots', 'TF_abundance_profile_barplot.pdf'),6,width=6.5)
bp
dev.off()



