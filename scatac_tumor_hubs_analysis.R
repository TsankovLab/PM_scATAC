
# Load functions for hub detection
source ('../../PM_scATAC/knnGen.R')
source ('../../PM_scATAC/addCoax.R')
source ('../../PM_scATAC/Hubs_finder.R')
source ('../../PM_scATAC/hubs_track.R')

### Run Hub detection and correlate with sarcomatoid TFs
k= 50
metaGroupName = 'Sample2'
KNNs = readRDS (paste0 ('KNN_',metaGroupName,'_k_',k,'.rds'))


###--- Hubs analysis --###
metaGroupName = "Sample2"
cor_cutoff = 0.2
max_dist = 12500
min_peaks = 5
dgs = 2000
projdir_hubs = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks,'/')

system (paste('mkdir -p',projdir_hubs))
system (paste('mkdir -p',paste0(projdir_hubs,'Plots')))

run_coax = FALSE
if (run_coax)
  {
  archp = addCoAx (
    archp, 
    KNNs,
    maxDist = max_dist)
  }

#archp = saveArchRProject (archp, load=T)
force=F
if (!file.exists (paste0(projdir_hubs,'global_hubs_obj.rds')) | force)
  {
  hubs_obj = hubs_finder (
    ArchRProj = archp, 
    group_by = NULL,
    cor_cutoff = cor_cutoff,
    #select_group = metaGroup_df$metaGroup,
    cor_FDR = 1, 
    cor_var = 0,
    min_peaks = min_peaks,
    macs_score = 1,
    dgs = dgs,
    cores=1
    ) 
  saveRDS (hubs_obj, paste0(projdir_hubs,'global_hubs_obj.rds'))
  write.table (as.data.frame (hubs_obj$hubsCollapsed), paste0(projdir_hubs, 'hub_regions.bed'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  } else {
  hubs_obj = readRDS (paste0(projdir_hubs,'global_hubs_obj.rds'))  
  }



# Export bed file of hub regions

# Plot KNN groups on UMAP
KNNs_df = lapply (seq_along(KNNs), function(x) data.frame (
  cell = KNNs[[x]], 
  group=paste0('KNN_',x), 
  group2 = unlist(strsplit (names(KNNs)[x],'KNN'))[1]))

KNNs_df = do.call (rbind, KNNs_df)
archp$knn_groups = KNNs_df$group[match (archp$cellNames, KNNs_df$cell)]
umap_knn = plotEmbedding (ArchRProj = archp, embedding = 'UMAP',
  colorBy = "cellColData", name = 'knn_groups',plotAs ='hex',
    baseSize=0, labelMeans=FALSE) + NoLegend() 
pdf (paste0(projdir_hubs, 'Plots/knn_test.pdf'), height=20, width=20)
umap_knn
dev.off()

fragments = unlist(getFragmentsFromProject (
  ArchRProj = archp))  
  
###--- Generate matrces of hubs x cells / knn from collapsed non-redundant hubs to identify differential hubs  ---###
metaGroup = as.character (archp@cellColData[,metaGroupName])
metaGroup_df = data.frame (barcode = rownames(archp@cellColData), metaGroup = metaGroup)

metaGroup3 = names (sapply (seq_along(KNNs), function(x)
  table(metaGroup_df$metaGroup[match (KNNs[[x]],metaGroup_df$barcode)])))

# Generate matrix of fragment counts of peaks x cluster
peaksKnn_mat = matrix (ncol = length(KNNs), nrow = length(hubs_obj$peaksMerged))
pb =progress::progress_bar$new(total = length (KNNs)) 
for (i in 1:length(KNNs)) {
    pb$tick()  
    peaks_fr = fragments[fragments$RG %in% KNNs[[i]]]  
    fragments_hubs = countOverlaps (hubs_obj$peaksMerged, peaks_fr)
    peaksKnn_mat[,i] = fragments_hubs
}
peaksKnn_mat = as.data.frame (peaksKnn_mat)
peaksKnn_mat$hubs_id = paste0(rep (hubs_obj$hubs_id,sapply(hubs_obj$hubsMerged,nrow)))

# Aggregate peaks x knn matrix to hub level taking mean
#gm_mean = function (x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) } # geometric mean, doesnt work as good as normal mean i think
hubsKnn_mat = aggregate (. ~ hubs_id, data = peaksKnn_mat, FUN = mean)
rownames (hubsKnn_mat) = apply (data.frame(as.data.frame (hubs_obj$hubsCollapsed), hubs_id = hubs_obj$hubs_id)[,c(7,1,2,3,6)], 1,
 function(x) paste(x, collapse = '_')) [match (hubsKnn_mat$hubs_id, hubs_obj$hubs_id)]
rownames (hubsKnn_mat) = gsub (' ','', rownames(hubsKnn_mat))
hubsKnn_mat = hubsKnn_mat[,-1]

# Normalise by seq depth
cellpool_nFrags = sapply(KNNs,function(v) sum(archp$ReadsInTSS[archp$cellNames %in% v]))
hubsKnn_mat = t(t(hubsKnn_mat) * (10^6 /cellpool_nFrags)) # scale
colnames (hubsKnn_mat) = metaGroup3
saveRDS (hubsKnn_mat, paste0 (projdir_hubs, 'hubs_knn_matrix_global.rds'))
hubsKnn_mat = readRDS (paste0 (projdir_hubs, 'hubs_knn_matrix_global.rds'))

# hubsClust_mat = as.data.frame (t (hubsKnn_mat))
# hubsClust_mat$metaGroup = metaGroup3
# hubsClust_mat = aggregate (.~ metaGroup, data = hubsClust_mat, FUN = mean)
# rownames (hubsClust_mat) = hubsClust_mat[,1]
# hubsClust_mat = hubsClust_mat[,-1]
# saveRDS (hubsClust_mat, paste0(projdir_hubs, 'hubsClust_global_mat.rds'))
# hubsClust_mat = readRDS (paste0(projdir_hubs, 'hubsClust_global_mat.rds'))


# Correlate with SOX9 genescore
if (!any (ls() == 'gsSE')) gsSE = ArchR::getMatrixFromProject (archp, useMatrix = 'GeneScoreMatrix')
gsSE = gsSE[, archp$cellNames]
gsMat = assays (gsSE)[[1]]
rownames (gsMat) = rowData (gsSE)$name

head (KNNs_df)
gsMat = gsMat[c(tf_name_selected, 'AXL','VIM','CALB2','ITLN1'),]
gsMat = as.data.frame (t(gsMat))
gsMat = gsMat[unlist(KNNs),]
gsMat$cellGroup = KNNs_df$group
gsMat = aggregate (.~ cellGroup, data = gsMat, FUN = mean)
rownames (gsMat) = gsMat[,1]
gsMat = gsMat[,-1]

hubsKnn_mat = t(hubsKnn_mat)
hubsKnn_mat = as.data.frame (hubsKnn_mat)
rownames (hubsKnn_mat) = rownames (gsMat)
TF_hub_cor = lapply (unique(KNNs_df$group2), function(x) cor (gsMat[unique(KNNs_df$group[KNNs_df$group2 == x]),'SOX9'], hubsKnn_mat[unique(KNNs_df$group[KNNs_df$group2 == x]),], method = 'spearman'))
names (TF_hub_cor) = unique(KNNs_df$group2)

TF_hub_cor2 = as.data.frame (t(TF_hub_cor[['P1']]))
TF_hub_cor2 = TF_hub_cor2[order (-TF_hub_cor2$V1),, drop=F]
head (TF_hub_cor2, 10)
rownames (TF_hub_cor2) = sapply (rownames(TF_hub_cor2), function(x) unlist(strsplit(x, '_'))[1])


KNNs_df2 = KNNs_df[!duplicated(KNNs_df$group),]
table (KNNs_df2$group2)
a = lapply (unique(KNNs_df$group2), function(x) cor (gsMat[unique(KNNs_df$group[KNNs_df$group2 == x]),'SNAI2'], gsMat[unique(KNNs_df$group[KNNs_df$group2 == x]),'SOX9'], method='spearman'))
names(a) = unique(KNNs_df$group2)






nfeat=5000
  k=25
  cnmf_list = readRDS (paste0('../scrna/cnmf_genelist_',k,'_nfeat_',nfeat,'.rds'))
  cnmf_list = lapply (cnmf_list, function(x) head (x,50))
gsMat = assays (gsSE)[[1]]
rownames (gsMat) = rowData (gsSE)$name
gsMat = gsMat[rownames (gsMat) %in% c(cnmf_list[[20]],'SOX9'),]
gsMat = as.data.frame (t(gsMat))

gsMat = gsMat[unlist(KNNs),]
gsMat$cellGroup = KNNs_df$group
gsMat = aggregate (.~ cellGroup, data = gsMat, FUN = mean)
rownames (gsMat) = gsMat[,1]
gsMat = gsMat[,-1]
gsMat_p1 = gsMat[unique(KNNs_df$group[KNNs_df$group2 == 'P1']),]

cor (gsMat_p1[,'TWIST1'],gsMat_p1[,'SNAI2'], method='spearman')

pdf ('Plots/cor_genescore_genes.pdf')
plot (gsMat[,'SOX9'], gsMat[,'VIM'])
dev.off 




# Check enrichment of TF in hubs using fgsea
TF = 'MESP1'
hubs_id = sapply (rownames(TF_hub_cor2), function(x) unlist(strsplit (x, '_'))[1])
tf_matches = getMatches (archp)
TF = colnames(tf_matches)[grep (TF, colnames(tf_matches))]


hub_peaks = rep (hubs_obj$hubs_id, sapply (hubs_obj$hubsClusters[[1]], nrow))
hubs_obj$peaksMerged$hub_id = hub_peaks
hubs_peaks_idx = subjectHits (findOverlaps (hubs_obj$peaksMerged,rowRanges(tf_matches)))
tf_matches = tf_matches[hubs_peaks_idx, colnames(tf_matches) == TF]
tf_matches_hit = which (assays (tf_matches)[[1]][,1])
hub_peaks_cor = setNames (TF_hub_cor2[hubs_obj$peaksMerged$hub_id,]$V1, paste0('h',1:nrow(tf_matches))) 
#tf_matches_peaks = unique(queryHits (findOverlaps (hubs_obj$hubsCollapsed[match(TF_hubs_id, hubs_obj$hubs_id)], rowRanges(tf_matches)[tf_matches_hit])))
pathways = list(SOX9 = paste0('h',tf_matches_hit))

library (fgsea)
hub_peaks_cor = na.omit (hub_peaks_cor)
fgseaRes = fgseaMultilevel (pathways, 
          hub_peaks_cor#, 
          #minSize=15, 
          #maxSize=1500,
          #BPPARAM = NULL
          )


tf_match = getMatches (archp)
colnames(tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)
top_cor_hubs = head (sapply (rownames(TF_hub_cor2), function(x) unlist(strsplit (x, '_'))[1]),100)
top_cor_hubs_peaks = lapply(top_cor_hubs, function(x) getPeakSet(archp) [queryHits(findOverlaps(getPeakSet(archp), hubs_obj$hubsCollapsed[which(hubs_obj$hubs_id %in% x)]))])
top_cor_hubs_TF = lapply (top_cor_hubs_peaks, function(x) hyperMotif (
  selected_peaks = x,
  motifmatch = tf_match))
  
top_cor_hubs_TF = lapply (top_cor_hubs_TF, function(x) x[rownames(top_cor_hubs_TF[[1]]), ])
  top_cor_hubs_TF_df = do.call (cbind, top_cor_hubs_TF)
  top_cor_hubs_TF_df = top_cor_hubs_TF_df[, grep ('padj', colnames(top_cor_hubs_TF_df))]
  top_cor_hubs_TF_df[top_cor_hubs_TF_df > 0.05] = 1
  top_cor_hubs_TF_df = -log10(top_cor_hubs_TF_df)
  top_cor_hubs_TF_df = top_cor_hubs_TF_df[rowSums (top_cor_hubs_TF_df) != 0, ]
  top_cor_hubs_TF_df[sapply(top_cor_hubs_TF_df, is.infinite)] <- 300
  
  TF_ht = Heatmap (top_cor_hubs_TF_df, row_names_gp = gpar (fontsize=3), column_names_gp = gpar (fontsize=5))
  
  pdf (paste0('Plots/TF_top_cor_hubs_heatmap.pdf'),width = 15,height=15)
  print (TF_ht)
  dev.off()





# Generate matrix of fragment counts of hubs x barcodes ####
if (!file.exists('hubs_cells_mat.rds'))
  {
  hubsCell_mat = matrix (ncol = length(rownames(archp@cellColData)), nrow = length(hubs_obj$hubsCollapsed))
  colnames (hubsCell_mat) = rownames(archp@cellColData)
  rownames (hubsCell_mat) = hubs_obj$hubs_id
  pb =progress::progress_bar$new(total = length (rownames(archp@cellColData)))
  for (cell in rownames(archp@cellColData)) 
    {
    pb$tick()  
    fragments_in_cell = fragments[fragments$RG %in% cell]  
    fragments_in_cell_in_hubs = countOverlaps (hubs_obj$hubsCollapsed, fragments_in_cell)
    hubsCell_mat[,cell] = fragments_in_cell_in_hubs
    }
  hubsCell_mat = t(t(hubsCell_mat) * (10^6 /archp$nFrags)) # scale
  saveRDS (hubsCell_mat, 'hubs_cells_mat.rds')  
  } else {
  hubsCell_mat = readRDS ('hubs_cells_mat.rds')  
  }

hubsCell_mat = as.data.frame (hubsCell_mat)

# Compute differential hub accessibility DHA ####
library (presto)
hubsCell_mat_comp = hubsCell_mat[,archp$Sample2 %in% c('P1','normal_pleura')]
comparison = archp$Sample2[archp$Sample2 %in% c('P1','normal_pleura')]
hubsCell_mat_comp2 = hubsCell_mat[,archp$Sample2 %in% c('P1','P8','P5','P4')]
comparison2 = archp$Sample2[archp$Sample2 %in% c('P1','P8','P5','P4')] 
comparison2 = ifelse(comparison2 == 'P1','P1','epit')  
res = wilcoxauc (log2(hubsCell_mat_comp+1), comparison)
res2 = wilcoxauc (log2(hubsCell_mat_comp2+1), comparison2)

# comparison2 = archp$Sample2 == 'P1'
# res_P1_rest = wilcoxauc(log2(hubsCell_mat[,comparison2 %in% c('P1','normal_pleura')]+1), comparison2[comparison2 %in% c('P1','normal_pleura')])

res_p1 = res[res$group == 'P1',]
res_p1$region = as.character(hubs_obj$hubsCollapsed)
#res_p1 = res_p1[res_p1$padj < 0.05 & res_p1$logFC > 0 & res_p1$avgExpr > 3,]

# Intersect with hubs higher in P1 vs epithelioid samples ####
res2_p1 = res2[res2$group == 'P1',]
res2_p1$region = as.character(hubs_obj$hubsCollapsed)
res2_p1 = res2_p1[res2_p1$padj < 0.05 & res2_p1$logFC > 0 & res2_p1$avgExpr > 3,]
res_p1 = res_p1[res_p1$feature %in% res2_p1$feature,]


# intersect with correlated hubs to SOX9 ####
res_p1$SOX9_cor = TF_hub_cor2[res_p1$feature,]$V1
res_p1 = res_p1[order (-res_p1$SOX9_cor),]
res_p1 = res_p1[res_p1$SOX9_cor > 0.1, ]

tf_match = getMatches (archp)
colnames(tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
tf_match = tf_match[queryHits(findOverlaps(rowRanges(tf_match),hubs_obj$hubsCollapsed)),]
top_cor_hubs = res_p1$feature
top_cor_hubs_peaks = lapply(top_cor_hubs, function(x) getPeakSet(archp) [queryHits(findOverlaps(getPeakSet(archp), hubs_obj$hubsCollapsed[which(hubs_obj$hubs_id %in% x)]))])
top_cor_hubs_TF = lapply (top_cor_hubs_peaks, function(x) hyperMotif (
  selected_peaks = x,
  motifmatch = tf_match))
names (top_cor_hubs_TF) = top_cor_hubs

top_cor_hubs_TF = lapply (top_cor_hubs_TF, function(x) x[rownames(top_cor_hubs_TF[[1]]), ])
top_cor_hubs_TF_df = do.call (cbind, top_cor_hubs_TF)
top_cor_hubs_TF_df = top_cor_hubs_TF_df[, grep ('padj', colnames(top_cor_hubs_TF_df))]
top_cor_hubs_TF_df[top_cor_hubs_TF_df > 0.05] = 1
top_cor_hubs_TF_df = -log10(top_cor_hubs_TF_df)
top_cor_hubs_TF_df = top_cor_hubs_TF_df[rowSums (top_cor_hubs_TF_df) != 0, ]
top_cor_hubs_TF_df[sapply(top_cor_hubs_TF_df, is.infinite)] <- 300

TF_ht = Heatmap (top_cor_hubs_TF_df, row_names_gp = gpar (fontsize=3), column_names_gp = gpar (fontsize=5))

pdf (paste0('Plots/TF_top_DHA_P1_hubs_heatmap.pdf'),width = 40,height=25)
print (TF_ht)
dev.off()

top_cor_hubs_TF_df['SOX9',]
res_p1[res_p1$feature == 'HUB22404',]
# all peaks of top hubs
top_cor_hubs = res_p1$feature[1:1000]
top_cor_hubs_peaks = getPeakSet(archp)[queryHits(findOverlaps(getPeakSet(archp), hubs_obj$hubsCollapsed[which(hubs_obj$hubs_id %in% top_cor_hubs)]))]
top_cor_hubs_TF = hyperMotif (
  selected_peaks = top_cor_hubs_peaks,
  motifmatch = tf_match)
head (top_cor_hubs_TF[order (top_cor_hubs_TF$padj),],100)






library (presto)
pSE = getMatrixFromProject (archp, useMatrix = 'PeakMatrix')
pMat = assays (pSE)[[1]]
hubs_peaks_idx1 = subjectHits (findOverlaps (hubs_obj$peaksMerged,rowRanges(pSE)))
pMat = pMat[hubs_peaks_idx1,]

pMat_comp = pMat[,archp$Sample2 %in% c('P1','normal_pleura')]
comparison = archp$Sample2[archp$Sample2 %in% c('P1','normal_pleura')]
pMat_comp2 = pMat[,archp$Sample2 %in% c('P1','P8','P5','P4')]
comparison2 = archp$Sample2[archp$Sample2 %in% c('P1','P8','P5','P4')] 
comparison2 = ifelse(comparison2 == 'P1','P1','epit')  
res = wilcoxauc (pMat_comp, comparison)
res2 = wilcoxauc (log2(pMat_comp2), comparison2)

# comparison2 = archp$Sample2 == 'P1'
# res_P1_rest = wilcoxauc(log2(hubsCell_mat[,comparison2 %in% c('P1','normal_pleura')]+1), comparison2[comparison2 %in% c('P1','normal_pleura')])

res_p1 = res[res$group == 'P1',]
res_p1 = res_p1[res_p1$avgExpr != 0,]
#res_p1$region = as.character(hubs_obj$hubsCollapsed)

# Try using fgsea ####
TF = 'MESP1'
tf_matches = getMatches (archp)
TF = colnames(tf_matches)[grep (TF, colnames(tf_matches))]

hubs_peaks_idx2 = subjectHits (findOverlaps (hubs_obj$peaksMerged, rowRanges(tf_matches)))
tf_matches = tf_matches[hubs_peaks_idx2, colnames(tf_matches) == TF]
tf_matches_hit = which (assays (tf_matches)[[1]][,1])
hub_peaks_rank = setNames (-log10(res_p1$pval) * res_p1$logFC, paste0('h',1:nrow(res_p1)))
#tf_matches_peaks = unique(queryHits (findOverlaps (hubs_obj$hubsCollapsed[match(TF_hubs_id, hubs_obj$hubs_id)], rowRanges(tf_matches)[tf_matches_hit])))
pathways = list(SOX9 = paste0('h',tf_matches_hit))

library (fgsea)
#hub_peaks_rank = na.omit (hub_peaks_rank)
fgseaRes = fgseaMultilevel (pathways, 
          hub_peaks_rank#, 
          #minSize=15, 
          #maxSize=1500,
          #BPPARAM = NULL
          )


tf_match = getMatches (archp)
colnames(tf_match) = sapply (colnames (tf_match), function(x) unlist(strsplit(x,'_'))[1])
bg_peakSet = rowRanges (tf_match)
top_cor_hubs = head (sapply (rownames(TF_hub_cor2), function(x) unlist(strsplit (x, '_'))[1]),100)
top_cor_hubs_peaks = lapply(top_cor_hubs, function(x) getPeakSet(archp) [queryHits(findOverlaps(getPeakSet(archp), hubs_obj$hubsCollapsed[which(hubs_obj$hubs_id %in% x)]))])
top_cor_hubs_TF = lapply (top_cor_hubs_peaks, function(x) hyperMotif (
  selected_peaks = x,
  motifmatch = tf_match))
  
top_cor_hubs_TF = lapply (top_cor_hubs_TF, function(x) x[rownames(top_cor_hubs_TF[[1]]), ])
  top_cor_hubs_TF_df = do.call (cbind, top_cor_hubs_TF)
  top_cor_hubs_TF_df = top_cor_hubs_TF_df[, grep ('padj', colnames(top_cor_hubs_TF_df))]
  top_cor_hubs_TF_df[top_cor_hubs_TF_df > 0.05] = 1
  top_cor_hubs_TF_df = -log10(top_cor_hubs_TF_df)
  top_cor_hubs_TF_df = top_cor_hubs_TF_df[rowSums (top_cor_hubs_TF_df) != 0, ]
  top_cor_hubs_TF_df[sapply(top_cor_hubs_TF_df, is.infinite)] <- 300
  
  TF_ht = Heatmap (top_cor_hubs_TF_df, row_names_gp = gpar (fontsize=3), column_names_gp = gpar (fontsize=5))
  
  pdf (paste0('Plots/TF_top_cor_hubs_heatmap.pdf'),width = 15,height=15)
  print (TF_ht)
  dev.off()






# # Try using clusters as KNN to improve accuracy of correlations ####
# # DOESNT WORK! ####
# fragments = unlist(getFragmentsFromProject (
#   ArchRProj = archp))  
  
# ###--- Generate matrces of hubs x cells / knn from collapsed non-redundant hubs to identify differential hubs  ---###
# metaGroupName = 'Clusters'
# metaGroup = as.character (archp@cellColData[,metaGroupName])
# metaGroup_df = data.frame (barcode = rownames(archp@cellColData), metaGroup = metaGroup)

# addArchRThreads (threads = 1) 
# peaks_sb = getGroupSE (archp, useMatrix = 'PeakMatrix', groupBy = 'Clusters')
# peaks_sb_mat = assays (peaks_sb)[[1]]
# gene_sb = getGroupSE (archp, useMatrix = 'GeneScoreMatrix', groupBy = 'Clusters')
# gene_sb_mat = assays (gene_sb)[[1]]
# rownames (gene_sb_mat) = rowData(gene_sb)$name

# cor (gene_sb_mat['SOX9',],gene_sb_mat['VIM',], method='spearman')

# # Correlate with SOX9 genescore
# if (!any (ls() == 'gsSE')) gsSE = ArchR::getMatrixFromProject (archp, useMatrix = 'GeneScoreMatrix')
# gsSE = gsSE[, archp$cellNames]
# gsMat = assays (gsSE)[[1]]
# rownames (gsMat) = rowData (gsSE)$name

# head (KNNs_df)
# gsMat = gsMat[c(tf_name_selected, 'AXL','VIM','CALB2','ITLN1'),]
# gsMat = as.data.frame (t(gsMat))
# gsMat = gsMat[unlist(KNNs),]
# gsMat$cellGroup = KNNs_df$group
# gsMat = aggregate (.~ cellGroup, data = gsMat, FUN = mean)
# rownames (gsMat) = gsMat[,1]
# gsMat = gsMat[,-1]

# hubsKnn_mat = t(hubsKnn_mat)
# hubsKnn_mat = as.data.frame (hubsKnn_mat)
# rownames (hubsKnn_mat) = rownames (gsMat)
# TF_hub_cor = lapply (unique(KNNs_df$group2), function(x) cor (gsMat[unique(KNNs_df$group[KNNs_df$group2 == x]),'SOX9'], hubsKnn_mat[unique(KNNs_df$group[KNNs_df$group2 == x]),], method = 'spearman'))
# names (TF_hub_cor) = unique(KNNs_df$group2)

# TF_hub_cor2 = as.data.frame (t(TF_hub_cor[['P1']]))
# TF_hub_cor2 = TF_hub_cor2[order (-TF_hub_cor2$V1),, drop=F]
# head (TF_hub_cor2, 10)
# rownames (TF_hub_cor2) = sapply (rownames(TF_hub_cor2), function(x) unlist(strsplit(x, '_'))[1])


# KNNs_df2 = KNNs_df[!duplicated(KNNs_df$group),]
# table (KNNs_df2$group2)
# a = lapply (unique(KNNs_df$group2), function(x) cor (gsMat[unique(KNNs_df$group[KNNs_df$group2 == x]),'MESP1'], gsMat[unique(KNNs_df$group[KNNs_df$group2 == x]),'SOX9'], method='spearman'))
# names(a) = unique(KNNs_df$group2)

# pdf ('Plots/cor_genescore_genes.pdf')
# plot (gsMat[,'SOX9'], gsMat[,'VIM'])
# dev.off 

# # Check enrichment of TF in hubs using fgsea
# TF = 'MESP1'
# hubs_id = sapply (rownames(TF_hub_cor2), function(x) unlist(strsplit (x, '_'))[1])
# tf_matches = getMatches (archp)
# TF = colnames(tf_matches)[grep (TF, colnames(tf_matches))]


# hub_peaks = rep (hubs_obj$hubs_id, sapply (hubs_obj$hubsClusters[[1]], nrow))
# hubs_obj$peaksMerged$hub_id = hub_peaks
# hubs_peaks_idx = subjectHits (findOverlaps (hubs_obj$peaksMerged,rowRanges(tf_matches)))
# tf_matches = tf_matches[hubs_peaks_idx, colnames(tf_matches) == TF]
# tf_matches_hit = which (assays (tf_matches)[[1]][,1])
# hub_peaks_cor = setNames (TF_hub_cor2[hubs_obj$peaksMerged$hub_id,]$V1, paste0('h',1:nrow(tf_matches))) 
# #tf_matches_peaks = unique(queryHits (findOverlaps (hubs_obj$hubsCollapsed[match(TF_hubs_id, hubs_obj$hubs_id)], rowRanges(tf_matches)[tf_matches_hit])))
# pathways = list(SOX9 = paste0('h',tf_matches_hit))

# library (fgsea)
# hub_peaks_cor = na.omit (hub_peaks_cor)
# fgseaRes = fgseaMultilevel (pathways, 
#           hub_peaks_cor#, 
#           #minSize=15, 
          #maxSize=1500,
          #BPPARAM = NULL
          )



# Compute differential hub accessibility DHA ####
library (presto)
hubsCell_mat_comp = hubsCell_mat[,archp$Sample2 %in% c('P1','normal_pleura')]
comparison = archp$Sample2[archp$Sample2 %in% c('P1','normal_pleura')]
hubsCell_mat_comp2 = hubsCell_mat[,archp$Sample2 %in% c('P1','P8','P5','P4')]
comparison2 = archp$Sample2[archp$Sample2 %in% c('P1','P8','P5','P4')] 
comparison2 = ifelse(comparison2 == 'P1','P1','epit')  
res = wilcoxauc (log2(hubsCell_mat_comp+1), comparison)
res2 = wilcoxauc (log2(hubsCell_mat_comp2+1), comparison2)

# comparison2 = archp$Sample2 == 'P1'
# res_P1_rest = wilcoxauc(log2(hubsCell_mat[,comparison2 %in% c('P1','normal_pleura')]+1), comparison2[comparison2 %in% c('P1','normal_pleura')])

res_p1 = res[res$group == 'P1',]
res_p1$region = as.character(hubs_obj$hubsCollapsed)
#res_p1 = res_p1[res_p1$padj < 0.05 & res_p1$logFC > 0 & res_p1$avgExpr > 3,]

# # Intersect with hubs higher in P1 vs epithelioid samples ####
# res2_p1 = res2[res2$group == 'P1',]
# #res2_p1$region = as.character(hubs_obj$hubsCollapsed)
# #res2_p1 = res2_p1[res2_p1$padj < 0.05 & res2_p1$logFC > 0 & res2_p1$avgExpr > 3,]
# res_p1 = res_p1[res_p1$feature %in% res2_p1$feature,]


# # intersect with correlated hubs to SOX9 ####
# res_p1$SOX9_cor = TF_hub_cor2[res_p1$feature,]$V1
# res_p1 = res_p1[order (-res_p1$SOX9_cor),]
# res_p1 = res_p1[res_p1$SOX9_cor > 0.1, ]

TF = 'SOX9'
tf_matches = getMatches (archp)
TF = colnames(tf_matches)[grep (TF, colnames(tf_matches))]

hub_peaks = rep (hubs_obj$hubs_id, sapply (hubs_obj$hubsClusters[[1]], nrow))
hubs_obj$peaksMerged$hub_id = hub_peaks

tf_matches = tf_matches[subjectHits(findOverlaps(hubs_obj$peaksMerged, rowRanges(tf_matches))),]
tf_matches = tf_matches[,colnames(tf_matches) == TF]
tf_matches_mat = as.data.frame(assays(tf_matches)[[1]])
tf_matches_mat$hubs_id = hubs_obj$peaksMerged$hub_id
tf_matches_mat = aggregate (. ~ hubs_id, data = tf_matches_mat, FUN = sum)
tf_matches_mat = tf_matches_mat[match(res_p1$feature, tf_matches_mat$hubs_id),]
tf_matches_mat[tf_matches_mat >0] = 1
hub_rank = setNames (res_p1$logFC, res_p1$feature)
pathways = list(res_p1$feature[tf_matches_mat[,2] > 0])
names (pathways) = TF
library (fgsea)

table (hub_rank[pathways[[1]]]>0)
table (res_p1$logFC>0)
#hub_peaks_rank = na.omit (hub_peaks_rank)
fgseaRes = fgseaMultilevel (pathways, 
          hub_rank#, 
          #minSize=15, 
          #maxSize=1500,
          #BPPARAM = NULL
          )

pdf (paste0('Plots/hubs_enrichmentplot_',TF,'.pdf'))
plotEnrichment(pathways[[TF]],
               hub_rank) + labs(title=TF)
dev.off()
