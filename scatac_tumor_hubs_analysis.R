
### Run Hub detection and correlate with sarcomatoid TFs
KNN = readRDS (paste0 ('KNN_',metaGroupName,'.rds'))


###--- Hubs analysis --###
metaGroupName = "Clusters"
cor_cutoff = 0.3
max_dist = 12500
min_peaks = 5
dgs = 0
projdir_hubs = paste0 ('hubs_obj_cor_',cor_cutoff,'_md_',max_dist,'_dgs_',dgs,'_min_peaks_',min_peaks,'/')
system (paste('mkdir -p',projdir_hubs))
system (paste('mkdir -p',paste0(projdir_hubs,'Plots')))


archp = addCoAx (
  archp, 
  KNN,
  maxDist = max_dist)

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
    cores=4
    ) 
  saveRDS (hubs_obj, paste0(projdir_hubs,'global_hubs_obj.rds'))
  } else {
  hubs_obj = readRDS (paste0(projdir_hubs,'global_hubs_obj.rds'))  
  }



# Export bed file of hub regions
write.table (as.data.frame (hubs_obj$hubsCollapsed), paste0(projdir_hubs, 'hub_regions.bed'), sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)

# Plot KNN groups on UMAP
KNNs_df = lapply (seq_along(KNNs), function(x) data.frame (cell = KNNs[[x]], group=paste0('KNN_',x)))
KNNs_df = do.call (rbind, KNNs_df)
archp$knn_groups = KNNs_df$group[match (archp$cellNames, KNNs_df$cell)]
umap_knn = plotEmbedding (ArchRProj = archp, embedding = 'UMAP',
  colorBy = "cellColData", name = 'knn_groups',plotAs ='hex',
    baseSize=0, labelMeans=FALSE) + NoLegend() 
pdf (paste0(projdir_hubs, 'Plots/knn_test.pdf'), height=20, width=20)
umap_knn
dev.off()

###--- Generate matrces of hubs x cells / knn from collapsed non-redundant hubs to identify differential hubs  ---###
metaGroup3 = names (sapply (seq_along(KNNs), function(x) 
  table(metaGroup_df$metaGroup2[match (KNNs[[x]],metaGroup_df$barcode)])))

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

hubsClust_mat = as.data.frame (t (hubsKnn_mat))
hubsClust_mat$metaGroup = metaGroup3
hubsClust_mat = aggregate (.~ metaGroup, data = hubsClust_mat, FUN = mean)
rownames (hubsClust_mat) = hubsClust_mat[,1]
hubsClust_mat = hubsClust_mat[,-1]
saveRDS (hubsClust_mat, paste0(projdir_hubs, 'hubsClust_global_mat.rds'))
hubsClust_mat = readRDS (paste0(projdir_hubs, 'hubsClust_global_mat.rds'))



