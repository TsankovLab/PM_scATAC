
#### Identify bins of divergent ATAC coverage and gene expression ####
# Get GRanges of bins excluding black list regions
ws = 10e5
ss = 10e5
ws = 10e4
ss = 10e4

# blacklist = toGRanges (paste0('/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Public_data/blacklisted_regions/ENCODE_blacklist/',"hg38-blacklist.v2.bed")) # taken from https://github.com/Boyle-Lab/Blacklist/tree/master/lists
# if (!file.exists (file.path (paste0 ('windows_',ws,'_',ss,'.rds'))))
#   {
#   windows = makeWindows (genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist,
#     windowSize = ws, slidingSize = ss)
#   saveRDS (windows, file.path (paste0 ('windows_',ws,'_',ss,'.rds')))
#   } else {
#   windows = readRDS (file.path (paste0 ('windows_',ws,'_',ss,'.rds')))
#   }
genome = BSgenome.Hsapiens.UCSC.hg38  
chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
chromSizes = chromSizes[seqnames(chromSizes) %in% paste0('chr',1:22),]
chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")  
windows <- unlist(slidingWindows(x = chromSizes, width = ws, step = ss))

if (!exists('fragments')) fragments = getFragmentsFromProject (archp)

celltypes = c('CD8_exhausted','CD4','CD8','Tregs','NK_KLRC1','NK_FGFBP2')
metaGroupName = 'celltype2'
metaGroupName2= 'P10'
metaGroupName2= unique (srt$sampleID)

# Bin rna counts ####
hg38_genes = genes (TxDb.Hsapiens.UCSC.hg38.knownGene)
symbols = toTable (org.Hs.egSYMBOL)
hg38_genes$symbol = symbols$symbol[match (hg38_genes$gene_id, symbols$gene_id)]
genes_in_bins = findOverlaps (windows,hg38_genes)

srt_sub = srt[,srt$sampleID %in% metaGroupName2]
pseudobulk_counts = lapply (celltypes, function(x) 
  {
  count_mat = srt_sub@assays$RNA@counts[, srt_sub@meta.data[,metaGroupName] %in% x]
  count_mat[count_mat == 0] = NA
  rowMeans(count_mat, na.rm=TRUE)
  })

pseudobulk_counts = do.call (cbind,pseudobulk_counts)
pseudobulk_counts[is.na(pseudobulk_counts)] = NA
colnames (pseudobulk_counts) = paste0(celltypes,'_scrna')
pseudobulk_counts = as.data.frame (pseudobulk_counts)
pseudobulk_counts = pseudobulk_counts[hg38_genes$symbol,]
pseudobulk_counts_mean = lapply(split(subjectHits(genes_in_bins), queryHits(genes_in_bins)), function(x) colSums(na.omit(pseudobulk_counts[x,, drop=F])))
pseudobulk_counts_mean = do.call (rbind, pseudobulk_counts_mean)

elementMetadata (windows) = cbind (elementMetadata (windows), pseudobulk_counts_mean[match(seq(windows),unique(queryHits(genes_in_bins))),])

# Bin atac fragments ####
fragments_sub = unlist(fragments[metaGroupName2])
fragments_sub = lapply (celltypes, function(x) countOverlaps (windows, fragments_sub[fragments_sub$RG %in% rownames(archp@cellColData)[as.character(archp@cellColData[,metaGroupName]) %in% x]]))
scatac_bins = do.call (cbind,fragments_sub)
colnames (scatac_bins) = paste0(celltypes,'_scatac')
elementMetadata (windows) = cbind (elementMetadata(windows), scatac_bins)

## Generate heatmap of atac vs rna bins ####
scatac_bins_cor = cor(as.matrix(elementMetadata (windows)[,grep('scatac',colnames(elementMetadata(windows)))]))
scrna_bins_cor = cor(as.matrix(elementMetadata (windows)[,grep('scrna',colnames(elementMetadata(windows)))]), use='complete.obs')

hm_dev = draw (Heatmap (scatac_bins_cor, 
  rect_gp = gpar(type = "none"),
  clustering_distance_rows='pearson' ,
  clustering_distance_columns = 'pearson', 
  col=rev(palette_deviation_correlation), 
  border=F,
  ,
  cell_fun = function(j, i, x, y, w, h, fill) {
        if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }}))

pdf (file.path('Plots','celltypes_similarity_bins_deviations_heatmap.pdf'),width = 4,height=3.2)
print (hm_dev)
dev.off()

hm = draw (Heatmap (scrna_bins_cor, 
  rect_gp = gpar(type = "none"),
  cluster_columns = F,
  cluster_rows = F,
  #clustering_distance_rows='pearson' ,
  #clustering_distance_columns = 'pearson', 
  col=palette_expression_correlation, 
  border=F,
  ,
  cell_fun = function(j, i, x, y, w, h, fill) {
        if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
        }}))

pdf (file.path('Plots','celltypes_similarity_bins_expression_heatmap.pdf'),width = 4,height=3.2)
print (hm)
dev.off()


# Fit regression line and compute residuals to compare scRNA vs scATAC coverage of NK KLRC1 and CD8 exhausted ####
cor_df = as.data.frame (windows[,c('CD8_exhausted_scrna','NK_KLRC1_scrna','CD8_exhausted_scatac','NK_KLRC1_scatac')], row.names=NULL)
cor_df$region = as.character(windows)
cor_df$region_short = paste0(seqnames(windows),':',seq(windows))

cor_df[is.na(cor_df)] = 0
model_scrna <- lm(CD8_exhausted_scrna ~ NK_KLRC1_scrna, data=cor_df)
standard_res_scrna <- rstandard (model_scrna)
model_scatac <- lm (CD8_exhausted_scatac ~ NK_KLRC1_scatac, data=cor_df)
standard_res_scatac <- rstandard (model_scatac)

cor_df$residuals_scrna = standard_res_scrna
cor_df$residuals_scatac = standard_res_scatac

cor_df2 = as.data.frame (windows[,c('CD8_exhausted_scrna','NK_KLRC1_scrna','CD8_exhausted_scatac','NK_KLRC1_scatac')], row.names=NULL)
cor_df$residuals_scrna[is.na(cor_df2$CD8_exhausted_scrna)] = NA

top_dif_residuals = head (cor_df$region_short[order(-(abs(cor_df$residuals_scrna) - abs(cor_df$residuals_scatac)))],10)
cor_df$top_regions = ''
cor_df$top_regions[match(top_dif_residuals,cor_df$region_short)] = top_dif_residuals

library(ggpointdensity)
library(viridis)
library (ggnewscale)

cor_df$rna_counts = rowMeans (as.data.frame(elementMetadata(windows)[,grep('scrna',colnames(elementMetadata(windows)))]))
cor_df$scatac_counts = rowMeans (as.data.frame(elementMetadata(windows)[,grep('scatac',colnames(elementMetadata(windows)))]))
cor_df$region = factor (cor_df$region, levels = cor_df$region)

scale_to_01 <- function(x) {
  na_ind = is.na(x)
  x[na_ind] = 0
  scaled = (x - min(x)) / (max(x) - min(x))
  scaled[na_ind] = NA
  return (scaled)
}

cor_df$rna_counts_scaled = scale_to_01(cor_df$rna_counts)
cor_df$scatac_counts_scaled = scale_to_01(cor_df$scatac_counts)

corp1 = ggplot (cor_df, 
  aes (x = rna_counts, y=abs(residuals_scrna))) + 
  geom_point(size=0.5, alpha=0.3) + gtheme_no_rot

corp2 = ggplot (cor_df, 
  aes (x = residuals_scatac, 
    y = residuals_scrna,
    size= rna_counts)) + 
#geom_point (aes (size = rna_counts), size=.2, alpha = 0.3, color='red') +
geom_pointdensity(aes(size=rna_counts), alpha=.5) +
  scale_color_viridis() + gtheme_no_rot


corp = ggplot () + 
geom_point (data = cor_df, aes (x = region, 
  y = residuals_scatac, size=scatac_counts_scaled), 
  alpha = 0.5, color='grey') +
geom_point (data = cor_df, aes(x = region, 
  y = residuals_scrna, size = rna_counts_scaled), 
alpha = 0.5, color='darkgreen') + 
#scale_colour_gradientn (colours = rev(viridis::turbo(100))) +
#new_scale_colour() +
geom_text_repel (data = cor_df, 
  aes (x = region, 
  y = residuals_scrna,label=top_regions), size=2, max.overlaps=10000) + 
scale_size(range = c(0, 1)) +
#facet_wrap(~seqnames, scales="free_x", ncol=23) +
#scale_colour_gradientn (colours = rev(viridis::turbo(100))) +
#geom_text(position = position_dodge(width = 1), aes(x=seqnames, y=-30)) +
theme(
    axis.text.x = element_blank(), # Remove x-axis labels
    axis.title.x = element_blank(), # Optional: Remove x-axis title
    axis.ticks.x = element_blank(),
    panel.spacing = unit(1, "lines") # Adjust space between groups
  )
pdf (file.path ('Plots',paste0('correlation_',paste(celltype_pair,collapse='_'),'_scatac_scrna2.pdf')))
corp1
corp2
dev.off()

pdf (file.path ('Plots',paste0('correlation_',paste(celltype_pair,collapse='_'),'_scatac_scrna.pdf')), width=8)
corp
dev.off()
