# Load utils functions palettes and packages ####
source (file.path('..','..','git_repo','utils','load_packages.R'))
source (file.path('..','..','git_repo','utils','useful_functions.R'))
source (file.path('..','..','git_repo','utils','ggplot_aestetics.R'))
source (file.path('..','..','git_repo','utils','scATAC_functions.R'))
source (file.path('..','..','git_repo','utils','palettes.R'))

set.seed(1234)


srt = readRDS ('GSM9326198_CROPseq_srt.rds')


### INFERCNV ####
library (infercnv)
ref = readRDS ('seurat_distal_normal_lung.rds')
normal_mesothelium_scrna_barcodes = read.csv (file.path ('..','..','git_repo','files','normal_lung_mesothelium_scrna_barcodes.csv'))
ref = ref[, colnames(ref) %in% normal_mesothelium_scrna_barcodes$x]  # Subset for only normal mesothelial cells

projdir_out = 'infercnv_CROPSEQ' # output folder 
dir.create (projdir_out, showWarnings = FALSE)

srt_cnv = srt
srt_cnv$cnv_type = srt_cnv$sampleID
ref$cnv_type = 'reference'
srt_cnv = merge (srt_cnv, ref)

# load TxDbi object for the species and get genomic positions of genes
require (org.Hs.eg.db)
require (TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(txdb) = paste0('chr',1:22) # select chromosomes to include in the analysis
gene_regions = as.data.frame (genes (txdb))

# map entrez id to symbol ids
symbol = toTable (org.Hs.egSYMBOL)
gene_regions$symbol = symbol$symbol[match(gene_regions$gene_id, symbol$gene_id)]
gene_regions = gene_regions[gene_regions$symbol %in% rownames(srt_cnv),]

message ('Generate expression and annotation files')
#srt_cnv$orig.ident = srt_cnv@meta.data [,metaGroupName2]
malig_samples = unique(srt_cnv$cnv_type[srt_cnv$cnv_type != 'reference'])
        
exprMat = srt_cnv@assays$RNA@counts        
gene_regions2 = gene_regions[mixedorder(gene_regions$seqnames), ]
exprMat = exprMat[gene_regions2$symbol, ]

all (rownames(exprMat) == gene_regions2$symbol)
rownames (gene_regions2) = gene_regions2$symbol
colnames (gene_regions2) = NULL
message ('save gene_regions file')    
write.table (gene_regions2, file.path(projdir_out,'gene_regions.txt'), sep='\t')
#exprMat = exprMat[,rownames(annot_df)]
annot_df = srt_cnv@meta.data[,'cnv_type', drop=F]
annot_df[,1] = srt_cnv$cnv_type

colnames (annot_df) = NULL
message ('save annotation file')    
write.table (annot_df, file.path(projdir_out,'annot_df.txt'), sep='\t', col.names=NA)

chr_exclude = c('chrX','chry')
infercnv_obj = CreateInfercnvObject (raw_counts_matrix=exprMat,
                                    annotations_file= file.path(projdir_out,'annot_df.txt'),
                                    delim="\t",
                                    gene_order_file=file.path(projdir_out,'gene_regions.txt'),
                                    ref_group_names=c('reference'),
                                    chr_exclude= chr_exclude
                                    )

infercnv_result = infercnv::run(infercnv_obj,
                           cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                           out_dir=projdir_out,  # dir is auto-created for storing outputs
                           cluster_by_groups=T,  # cluster
                           denoise=T,
                           HMM=F,
                           save_rds = F,
                           no_prelim_plot = T,
                           no_plot = T,
                           plot_probabilities = FALSE
                           )


# Identify cells with one gRNA called or in combination with NTC ####
x <- srt$crispr_calls

keep <- sapply(strsplit(x, "\\|"), function(parts) {
  # Rule 1: single gene
  if (length(parts) == 1) return(TRUE)
  
  # Check if NTC is present
  ntc_present <- any(grepl("^NTC", parts))
  
  if (ntc_present) {
    # Remove NTC from list
    other_genes <- parts[!grepl("^NTC", parts)]
    # After removing NTC, get unique base gene names
    gene_names <- sub("_\\d+$", "", other_genes)
    # Keep only if there is exactly ONE unique non-NTC gene
    return(length(unique(gene_names)) == 1)
  }
  
  # Rule 3: all same gene name ignoring _number
  gene_names <- sub("_\\d+$", "", parts)
  if (length(unique(gene_names)) == 1) return(TRUE)
  
  return(FALSE)
})
srt$unique_guide = keep
ccomp = as.data.frame (table (srt$unique_guide))
bp = ggplot (ccomp, aes (x = Var1, y = Freq, fill = Var1)) + geom_bar (stat = 'identity') + 
scale_fill_manual (values = c(`TRUE` = 'brown',`FALSE` = 'grey')) + gtheme
bp$data$Var1 = factor (bp$data$Var1, levels = c('TRUE','FALSE'))
pdf (file.path ('Plots','fraction_unique_guides_barplot.pdf'), width=3, height=3)
bp
dev.off()


### Filter non-unique guides ####
srt <- srt[,keep]


# # Remove cells with no guides ####
res2 = as.data.frame (table (srt$crispr_calls))
res2 = res2[order(-res2$Freq),]
res2$Var1 = factor (res2$Var1, levels = unique (res2$Var1))
bp = ggplot (res2, aes (x = Var1, y = Freq)) + geom_bar (stat = 'identity') + gtheme

pdf (file.path ('Plots','gRNA_abundances_not_filtered_barplot.pdf'), width=93, height=20)
bp
dev.off()


res2 = as.data.frame (table (srt$crispr_calls))
res2 = res2[order(-res2$Freq),]
res2$Var1 = factor (res2$Var1, levels = unique (res2$Var1))
bp = ggplot (res2, aes (x = Var1, y = Freq)) + geom_bar (stat = 'identity') + gtheme

pdf (file.path ('Plots','gRNA_abundances_barplot.pdf'), width=13, height=5)
bp
dev.off()

### Further remove guides with less than n cells ####
# Show genes left in after also filtering by lowest amount of cells 
cells_threshold = 9
gene_left = as.character(res2$Var1[res2$Freq > cells_threshold])
srt = srt[,srt$crispr_calls %in% gene_left]

pdf (file.path ('Plots','gRNA_abundances_barplot.pdf'), width=10, height=5)
bp + geom_hline(yintercept = cells_threshold)
dev.off()

### Check consistency of guides ####
var_feat_avg = log2(AverageExpression (srt, group.by = 'crispr_calls', features = VariableFeatures (srt))[[1]]+1)
#var_feat_avg = t(scale(t(var_feat_avg)))
var_feat_avg_cor = cor(as.matrix(var_feat_avg), method='pearson')
colnames (var_feat_avg_cor) = gsub ('-','_', colnames(var_feat_avg_cor))
#ha = HeatmapAnnotation (cells = anno_barplot(as.numeric(table (srt$crispr_calls)[colnames(var_feat_avg_cor)])))
hm = Heatmap (var_feat_avg_cor[rownames(var_feat_avg_cor) != 'BPTF-2|BPTF-3',colnames(var_feat_avg_cor) != 'BPTF-2|BPTF-3'], 
 # top_annotation = ha, 
  clustering_distance_columns = 'pearson',
    clustering_distance_rows = 'pearson',
    row_names_gp = gpar (fontsize = 5), 
    column_names_gp = gpar (fontsize = 5),
    border=T
    )

pdf (file.path ('Plots','guides_cor_heatmap.pdf'),height = 4,width = 4.5)
hm
dev.off()

# Merge guides ####
merged_guides = list (
  MEF2D = c('MEF2D-1','MEF2D-3','MEF2D−3|NTC−3','MEF2D−1|MEF2D−3','MEF2D−1|NTC−3','MEF2D−1|MEF2D−3|NTC−3','MEF2D−1|NTC−1'),
  MEF2A = c('MEF2A−3','MEF2A−1','MEF2A−3|NTC−3','MEF2A−1|NTC−3','MEF2A−3|NTC−1'),
  PITX1 = c('PITX1−2','PITX1−3','PITX1−1','PITX1−2|PITX1−3','PITX1−1|PITX1−2','PITX1−1|PITX1−3','NTC-1|PITX1-2','NTC-3|PITX1-2','NTC-3|PITX1-3','NTC-1|PITX1-3'),
  BPTF = c('BPTF−3','BPTF−2'),
  HMGA1 = 'HMGA1−3',
  TEAD2 = 'TEAD2-1',
  NTC = c('NTC-1','NTC-3'),
  TCF3 = c('TCF3-1','TCF3-3'),
  SOX9 = 'SOX9-3',
  TWIST1 = c('TWIST1-1','TWIST1-3'),
  TEAD4 = c('TEAD4-2','TEAD4-3'))
merged_guides = lapply (merged_guides, function(x) gsub("[−–-]", "_", x))
srt$merged_call = srt$crispr_calls

for (i in names(merged_guides)) srt$merged_call[srt$merged_call %in% merged_guides[[i]]] = i
table (srt$merged_call)  

srt = srt[, !grepl ('\\|', srt$merged_call)]


### QC per guide ####
ccomp_df = as.data.frame (table(srt$merged_call))
cc_p2 = ggplot (ccomp_df, aes (x= Var1, y= Freq)) +
        geom_bar (position="stack", stat="identity") +
        theme (axis.text.x = element_text (angle = 90, vjust = 0.5, hjust=1)) + 
        ggtitle (paste('Tot cells',ncol(srt))) + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
srt$nFeature_RNAL = log10 (srt$nFeature_RNA)
srt$nCount_RNAL = log10 (srt$nCount_RNA)
vln_p = VlnPlot (srt, features = c("nFeature_RNAL", "nCount_RNAL", "percent.mt"), combine=T,group.by = 'merged_call',pt.size = 0, ncol = 3)
png (file.path("Plots","QC_nFeat_nCount_m.percent_vlnPlot.png"), 3000, 1000, res=300)
print (cc_p2 | vln_p[[1]] | vln_p[[2]] | vln_p[[3]]) + plot_layout (widths=c(1,2,2,2))
dev.off()




### Show changes in variable genes ####
var_feat_avg = log2(AverageExpression (srt, group.by = 'merged_call', features = VariableFeatures(FindVariableFeatures(srt, nfeat = 5000)))[[1]]+1)
#var_feat_avg = t(scale(t(var_feat_avg)))
#var_feat_avg_cor = cor(as.matrix(var_feat_avg), method='pearson')
colnames (var_feat_avg) = gsub ('-','_', colnames(var_feat_avg))
var_scaled = scale(t(as.matrix(var_feat_avg)))
var_scaled[is.na(var_scaled)] = 0
var_scaled = var_scaled[,apply(var_scaled,2,var) > 0]
#var_scaled = var_scaled[!rownames(var_scaled) %in% c('HMGA1','TEAD2'),]
ncells = table (srt$merged_call)
#ncells = ncells[!names(ncells) %in% c('HMGA1','TEAD2')]

ha = HeatmapAnnotation (cells = anno_barplot(as.numeric(ncells[rownames(var_scaled)])), which='row')
hm = Heatmap (var_scaled, right_annotation = ha, 
  row_names_side = 'left',
  row_split = ifelse (rownames(var_scaled) == 'NTC','control','KO'),
  clustering_distance_columns = 'pearson',
  , clustering_distance_rows = 'pearson', 
  column_names_gp=gpar (fontsize=0),
  row_names_gp=gpar (fontsize=8),
  col = rev(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 100)),
  border=T)
pdf (file.path ('Plots','most_variable_genes_heatmap.pdf'),width=4,height=2.5)
hm
dev.off()


#### Show expression of guides in relative genes ####
guides = unique(srt$merged_call)
guides = guides[!guides %in% 'NTC']
srt$celltype = 'crispr'
gdot_p2 = geneDot (
  seurat_obj = srt, 
  gene = unname(unlist(guides)),
  y = 'celltype',
  x = 'merged_call',
  scale.data = T,
  include_NA = FALSE,
  min_expression = 0,
  returnDF = F,
  plotcol = as.character(rev(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 100))))
gdot_p2$data = gdot_p2$data %>%
  filter(x_axis == "NTC" | gene == x_axis)
gdot_p2$data$x_axis = ifelse(gdot_p2$data$x_axis == 'NTC','control', 'guide')

dp = DotPlot (srt, group.by = 'merged_call', features = c('SOX9'), scale =F)
dp_data1 = dp$data[dp$data$id %in% c('NTC', 'SOX9'),]
dp_data2 = gdot_p2$data[gdot_p2$data$y_axis %in% c('NTC', 'SOX9'),]
pdf (file.path('Plots','guides_expression_dotplot.pdf'), height=3.5, width=3)
#gdot_p2
gdot_p2
dev.off()


### Boxplot cc genes ####
ccomp = srt@meta.data[,c('cc','merged_call','S.Score','G2M.Score'), drop=F]
bxpv = geom_boxplot (
    linewidth = .3,
    width=.1,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.8, alpha=.6
     ) 

vlp = geom_violin (trim=TRUE,size=2,
    width=1,
    scale='width',
    linewidth = .3, alpha=.6)

ccomp2 = aggregate (ccomp[,'cc', drop=F], by = list(merged_call = srt$merged_call), mean)
ccomp$merged_call = factor (ccomp$merged_call, levels = ccomp2$merged_call[order(-ccomp2$cc)])
  

gp1 = ggplot (ccomp, aes (x = merged_call, y = cc)) + geom_violin (
  aes (fill = merged_call),
  trim=TRUE,size=2,
    width=1,
    scale='width',
    linewidth = .3, alpha=.6) + 
geom_boxplot (
    linewidth = .3,
    width=.1,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.8, alpha=.6
     )  + gtheme +
scale_fill_manual (values = c(NTC = 'grey', TWIST1='brown',
  TCF3 = 'brown', PITX1 = 'brown', TEAD4 = 'brown', TEAD2 = 'brown',
  MEF2A = 'brown', MEF2D = 'brown', SOX9 = 'brown', BPTF = 'brown', HMGA1 = 'brown'))

stat.test <- ccomp %>%
  wilcox_test(reformulate('merged_call', 'cc')) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test = stat.test[stat.test$group1 == 'NTC',]  
stat.test <- stat.test %>% add_xy_position (x = "merged_call", step.increase=0.3)
gp1 = gp1 + stat_pvalue_manual (stat.test, remove.bracket=FALSE,
  bracket.nudge.y = -1, hide.ns = TRUE,
   label = "p.adj.signif") 
# gp2 = ggplot (ccomp, aes (x = merged_call, y = G2M.Score)) + geom_boxplot()
# gp3 = ggplot (ccomp, aes (x = merged_call, y = S.Score)) + geom_boxplot()
pdf (file.path ('Plots','cellcycle_boxplots.pdf'), height=3, width = 4)
gp1
dev.off()


#### Check viability CRISPR assay from DepMap database ####
# https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2026Q1&filename=CRISPRGeneEffect.csv
ccle_cri = read.csv ('CRISPRGeneEffect.csv')
rownames (ccle_cri) = ccle_cri[,1]
ccle_cri = as.data.frame (t(ccle_cri[,-1]))
ccle_cri_sub = ccle_cri[,'ACH-000153', drop=F]
row_names = sapply (rownames (ccle_cri_sub), function(x) unlist(strsplit(x, '\\...'))[1])
ccle_cri_sub = ccle_cri_sub[!duplicated(row_names),, drop=F]
rownames (ccle_cri_sub) = row_names [!duplicated(row_names)]
#ccle_cri_sub = as.data.frame (t(ccle_cri_sub))

ccle_res = ccle_cri_sub[unique(srt$merged_call), , drop=F]
ccle_res = na.omit (ccle_res)
colnames (ccle_res)[1] = 'ACH000153'
ccle_res$TF = rownames(ccle_res)
bp = ggplot (ccle_res, aes (x = TF, y = ACH000153)) + geom_bar (stat= 'identity', fill = 'brown') + gtheme
bp$data$TF = factor (bp$data$TF, levels = levels (gp1$data$merged_call)[levels (gp1$data$merged_call) != 'NTC']) 
pdf (file.path ('Plots','CCLE_CRISPR_screens_barplot.pdf'),4,height=3)
bp
dev.off()

### Plot Cms scores from PM scRNA-seq Giotti et al 2024. ####
cnmf_spectra_unique = readRDS (file.path ('..','..','git_repo','files','scRNA_PM_modules.rds'))

# compute module scores ####
if (!all (names(cnmf_spectra_unique) %in% colnames(srt@meta.data))) {
srt = ModScoreCor (
        seurat_obj = srt, 
        geneset_list = lapply(cnmf_spectra_unique, function(x) head(x,50)),
        cor_threshold = NULL, 
        pos_threshold = NULL, # 
        listName = 'Cm', outdir = paste0(projdir,'Plots/'))
}

df = as.data.frame (srt@meta.data) 
df$merged_call[is.na (df$merged_call)] = 'NT'
df = df[!is.na(df$merged_call),]
#top_kds = names(head(table (df$crispr_calls)[order(-table(df$crispr_calls))],50))
#df = df[df$crispr_calls %in% top_kds,]

# Extract Cm columns
cm_cols <- grep("^Cm\\d+", colnames(df), value = TRUE)

# Get mean expression for NTC_1
ntc_df1 <- df[df$merged_call == "NTC", ]
#ntc_df3 <- df[df$merged_call == "NTC_3", ]
ntc_means1 <- colMeans(ntc_df1[, cm_cols], na.rm = TRUE)
ntc_means1 = apply (ntc_df1[, cm_cols], 2, median)
#ntc_means3 <- colMeans(ntc_df3[, cm_cols], na.rm = TRUE)
#ntc_means = rowMeans (cbind(ntc_means1, ntc_means3))
ntc_means = ntc_means1

# Get all unique crispr_calls (including NTC_1)
all_calls <- unique(df$merged_call)
# Initialize matrix
diff_matrix <- matrix(nrow = length(all_calls), ncol = length(cm_cols))
rownames(diff_matrix) <- all_calls
colnames(diff_matrix) <- cm_cols

# Compute difference: NTC_1 mean - crispr_call mean
for (call in all_calls) {
  group_df <- df[df$merged_call == call, ]
  #group_means <- colMeans(group_df[, cm_cols], na.rm = TRUE)
  group_means = apply (group_df[, cm_cols], 2, median)
  diff_matrix[call, ] <- ntc_means - group_means
}

# Initialize p-value matrix
pval_matrix <- matrix(NA, nrow = length(all_calls), ncol = length(cm_cols))
rownames(pval_matrix) <- all_calls
colnames(pval_matrix) <- cm_cols

# Loop over calls and Cm columns
for (call in all_calls) {
  
  group_df <- df[df$merged_call == call, ]
  
  for (cm in cm_cols) {
    control_vals <- df[df$merged_call %in% c("NTC"), cm]
    group_vals   <- group_df[[cm]]
    
    # Check if both groups have >1 unique values
    if (length(unique(control_vals)) > 1 && length(unique(group_vals)) > 1) {
      pval_matrix[call, cm] <- wilcox.test(control_vals, group_vals, exact = FALSE)$p.value
    } else {
      pval_matrix[call, cm] <- NA
    }
  }
}
#diff_matrix[diff_matrix > 0.2] = 0.2
#diff_matrix[diff_matrix < -0.2] = -0.2

diff_matrix = diff_matrix[rownames (diff_matrix) != 'NTC',]
pval_matrix = pval_matrix[rownames (pval_matrix) != 'NTC',]
diff_matrix[diff_matrix > min(c(abs(min(diff_matrix)), abs(max (diff_matrix))))] = min(c(abs(min(diff_matrix)), abs(max (diff_matrix))))
diff_matrix[diff_matrix < -min(c(abs(min(diff_matrix)), abs(max (diff_matrix))))] = -min(c(abs(min(diff_matrix)), abs(max (diff_matrix))))
ha = HeatmapAnnotation (cells = anno_barplot(as.vector(unname(table (srt$merged_call)[rownames(diff_matrix)]))), which= 'row')
x <- seq(-1, 1, length.out = 100)
# Make palette: blue (-1) -> white (0) -> red (1)
pal <- colorRampPalette(rev(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 100)))
# Generate 100 interpolated colors
cols <- pal(100)
pdf (file.path ('Plots','knock_downs_heatmap2.pdf'), height=2.7, width=5)
Heatmap (diff_matrix, row_names_gp=gpar(fontsize = 7), right_annotation=ha,
  column_names_rot=45,
  border=T, row_names_side = 'left',
  column_names_gp = gpar (fontsize=7),
  #col = cols,
  col = rev(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 100)),
  #col = col_fun,
  cell_fun = function (j, i, x, y, width, height, fill) 
            {
           if (pval_matrix[i, j] < 0.001)
              {
               grid.text("***", x, y, just='center', vjust=.8,
                gp = gpar(fontsize = 5, col='black'))
              } else {
              if(pval_matrix[i, j] < 0.01)
                  {
                  grid.text("**", x, y, just='center', vjust=.8,
                  gp = gpar(fontsize = 5, col='black'))   
                  } else {
                  if(pval_matrix[i, j] < 0.05)
                    {
                    grid.text("*", x, y, just='center', vjust=.8,
                    gp = gpar(fontsize = 5, col='black'))         
                    }}}
      })
dev.off()
