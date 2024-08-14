### Gene score based analysis ####
# Find DAG ####
force = FALSE
if (!file.exists (paste0('DAG_',metaGroupName,'.rds')) | force)
  {
  DAG_list = getMarkerFeatures (
    ArchRProj = archp, 
    testMethod = "wilcoxon",
          #useGroups = "ClusterA",
          #bgdGroups = "Clusters1B",
    binarize = FALSE,
    useMatrix = "GeneScoreMatrix",
    groupBy = metaGroupName
  #  useSeqnames="z"
  )
  listnames = colnames (DAG_list)
  DAG_list = lapply (1:ncol (DAG_list), function(x) 
    {
    df = DAG_list[,x]  
    df = do.call (cbind, (assays(df)))
    colnames(df) = names (assays(DAG_list))
    df$gene = rowData (DAG_list)$name
    df
    })
  names (DAG_list) = listnames
  saveRDS (DAG_list, paste0 ('DAG_',metaGroupName,'.rds'))    
  } else {
  DAG_list = readRDS (paste0('DAG_',metaGroupName,'.rds'))
  }

FDR_threshold = 1e-8
lfc_threshold = 1
top_genes = 20
DAG_top_list = DAG_list[sapply (DAG_list, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$Log2FC) > lfc_threshold,]) > 0)]
DAG_top_list = lapply (seq_along(DAG_top_list), function(x) {
  res = DAG_top_list[[x]]
  res = na.omit (res)
  res = res[res$FDR < FDR_threshold,]
  res = res[order (res$FDR), ]
  res = res[abs(res$Log2FC) > lfc_threshold,]
  res$comparison = names(DAG_top_list)[x]
  if (nrow(res) < top_genes) 
    {
    res
    } else {
    head (res,top_genes)
    }
  })
DAG_df = Reduce (rbind ,DAG_top_list)

if (!any (ls() == 'gsSE')) gsSE = ArchR::getMatrixFromProject (archp, useMatrix = 'GeneScoreMatrix')
gsSE = gsSE[, archp$cellNames]
gsMat = assays (gsSE)[[1]]
rownames (gsMat) = rowData (gsSE)$name
gsMat_mg = gsMat[rownames (gsMat) %in% DAG_df$gene, ]
gsMat_mg = as.data.frame (t(gsMat_mg))
gsMat_mg$metaGroup = as.character(archp@cellColData[,metaGroupName])
gsMat_mg = aggregate (.~ metaGroup, gsMat_mg, mean)
rownames (gsMat_mg) = gsMat_mg[,1]
gsMat_mg = gsMat_mg[,-1]
gsMat_mg = gsMat_mg[names(table (archp@cellColData[,metaGroupName])[table (archp@cellColData[,metaGroupName]) > 50]),]
DAG_hm = Heatmap (t(scale(gsMat_mg)), 
        row_labels = colnames (gsMat_mg),
        column_title = paste('top',top_genes),
        clustering_distance_columns = 'euclidean',
        clustering_distance_rows = 'euclidean',
        cluster_rows = T,
        #col = pals_heatmap[[5]],
        cluster_columns=T,#col = pals_heatmap[[1]],
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 4),
        rect_gp = gpar(col = "white", lwd = .5),
        border=TRUE
        #right_annotation = motif_ha
        )
         
  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (paste0('Plots/DAG_clusters_',metaGroupName,'_heatmaps.pdf'), width = 8, height = 50)
print (DAG_hm)
dev.off()
  