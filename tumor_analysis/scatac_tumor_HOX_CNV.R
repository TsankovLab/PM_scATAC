

### Correlate HOX deviation with CNVs ####
p11cnv = readRDS (file.path('..','..','per_sample_QC_signac','CNV_analysis','CNV_LFC_GC_P11_ws_1e+07_ss_5e+06.rds'))
p11cnv_mat = assays(p11cnv)$counts
p11cnv_mat = scale (p11cnv_mat)

#if (!exists('gsSE')) gsSE = fetch_mat(archp, 'GeneScore')
if (!exists('mSE')) mSE = fetch_mat(archp, 'Motif')
mMat = scale(assays (mSE)[[1]])
rownames (mMat) = rowData (mSE)$name


# # Get deviation matrix ####
#gsMat = assays (gsSE)[[1]]
#rownames (gsMat) = rowData (gsSE)$name

# Get active TFs ####
selected_TF = read.csv ('Active_TFs.csv', row.names=1)
#p11mMat = gsMat[, archp$Clusters == 'C1']
p11mMat = mMat[, archp$Clusters == 'C1']
#p11mMat = mMat[, archp$Sample2 == 'P11']
colnames(p11mMat) = sapply (colnames(p11mMat), function(x) unlist(strsplit(x, '#'))[2])
p11cnv_mat = p11cnv_mat[,colnames (p11mMat)]

p11mMat = p11mMat[grepl ('^HOX', rownames(p11mMat)) | grepl ('^TEAD', rownames(p11mMat)), ]

# correlate
p11mMat = t(p11mMat)
p11cnv_mat = t(p11cnv_mat)

cnv_hox_cor = cor (as.matrix(p11mMat), as.matrix(p11cnv_mat), method='spearman')

hm = Heatmap  (cnv_hox_cor, cluster_columns = F, 
  column_names_gp = gpar(fontsize = 0), column_split = factor (sapply (colnames(cnv_hox_cor), function(x) unlist(strsplit(x,':'))[1]),levels=unique(sapply (colnames(cnv_hox_cor), function(x) unlist(strsplit(x,':'))[1]))),
  row_names_gp = gpar(fontsize = 4),
  column_title_rot = 45)
pdf (file.path ('Plots','cnv_HOX_P11_cor_gs_heatmap.pdf'), width=10, height=4)
hm
dev.off()

min_cor = cnv_hox_cor['HOXB13',]
min_cor = min_cor[order(min_cor)]

pdf (file.path ('Plots','cnv_TF_cor_gs.pdf'))
plot (p11mMat[,'HOXB13'], p11cnv_mat[,names(min_cor)[1]])
dev.off()
cor (p11mMat[,'HOXB13'], p11cnv_mat[,names(min_cor)[1]], method = 'spearman')





### Look at overlap of HOXB13 with specific CNVs from TCGA ####
## Generate circos plot for showing hubs on recurrent CNVs ####
# source script to load TCGA_CNV data
source (file.path('..','..','git_repo','tumor_analysis','compile_TCGA_CNV_data.R'))
meso_CNV_gr_hg38_cnv = readRDS ('TCGA_CNV_hg38.rds')    

# Compute bins ####
ws = 1e7
ss = 5e6
blacklist = toGRanges (file.path ('..','..','git_repo','files',"hg38-blacklist.v2.bed")) # taken from https://github.com/Boyle-Lab/Blacklist/tree/master/lists

if (!file.exists (paste0('windows_',ws,'_',ss,'.rds')))
  {
  windows = makeWindows (genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist,
    windowSize = ws, slidingSize = ss)
  saveRDS (windows, paste0('windows_',ws,'_',ss,'.rds'))
  } else {
  windows = readRDS (paste0('windows_',ws,'_',ss,'.rds'))  
  }
windows = windows[!grepl ('chrX',seqnames(windows))]
windows = windows[!grepl ('chrY',seqnames(windows))]
pb = progress::progress_bar$new(total = length(meso_CNV_gr_hg38_cnv))
cnv_mat = matrix (nrow = length(windows), ncol = length(meso_CNV_gr_hg38_cnv))
colnames (cnv_mat) = names (meso_CNV_gr_hg38_cnv)
rownames (cnv_mat) = as.character(windows)
for (sam in names(meso_CNV_gr_hg38_cnv))
  {
  pb$tick()  
  ov = as.data.frame (findOverlaps (windows, meso_CNV_gr_hg38_cnv[[sam]]))
  ov_sum1 = split (ov, ov$queryHits)
  ov_sum2 = unlist(lapply (ov_sum1, function(x) mean (meso_CNV_gr_hg38_cnv[[sam]]$Segment_Mean[x$subjectHits])))
  cnv_mat[as.numeric(names(ov_sum1)),sam] = ov_sum2
  }
  


# load bulkRNA cohorts ####
gene = 'HOXB13'
study = 'tcga'
meso_bulk_l = readRDS (file.path ('..','..','bulkRNA_meso','bulk_RNA_studies.rds'))
meso_bulk_meta_l = readRDS (file.path ('..','..','bulkRNA_meso','bulk_RNA_studies_metadata.rds'))

bulk_gene = meso_bulk_l[[study]][gene,]
names (bulk_gene) = gsub ('\\.','-', names (bulk_gene))
colnames(cnv_mat) = sapply(colnames(cnv_mat), function(x) unlist(strsplit(x, '01A-'))[1])
colnames(cnv_mat) = paste0(colnames(cnv_mat), '01')
bulk_gene = bulk_gene[colnames(cnv_mat)]
cnv_mat[is.na(cnv_mat)] = 0
cnv_mat = cnv_mat[!grepl ('chrX',rownames(cnv_mat)),]
cnv_mat = cnv_mat[!grepl ('chrY',rownames(cnv_mat)),]
ha = rowAnnotation (gene = bulk_gene)
pdf (file.path('Plots',paste0(gene,'_CNV_TCGA_data_heatmap.pdf')),width=12, height=4)
Heatmap (t(cnv_mat), clustering_distance_rows = 'pearson',
  column_split =seqnames(windows),
  column_title_rot=45,
  right_annotation = ha, 
  cluster_columns=F,
  column_names_gp=gpar (fontsize=0), 
  row_names_gp=gpar (fontsize=0))
dev.off()

# correlate gene with CNV ####
gene_cnv_pvalue = p.adjust(sapply (1:nrow(cnv_mat), function(x) cor.test (cnv_mat[x,], bulk_gene, method = 'spearman')$p.value))
gene_cnv_cor = cor (t(cnv_mat), bulk_gene, method='spearman')
res = data.frame (cor = gene_cnv_cor, pvalue = gene_cnv_pvalue)
head (res[order(res[,2]), ,drop=F],20)
table (res$pvalue)

