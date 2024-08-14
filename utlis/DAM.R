
### ChromVAR based analysis ####
  # Find DAM ####

  if (!file.exists (paste0('DAM_',metaGroupName,'.rds')) | force)
    {
    DAM_list = getMarkerFeatures (
      ArchRProj = archp, 
      testMethod = "wilcoxon",
            #useGroups = "ClusterA",
            #bgdGroups = "Clusters1B",
      binarize = FALSE,
      useMatrix = "MotifMatrix",
      groupBy = metaGroupName
    #  useSeqnames="z"
    )

    listnames = colnames (DAM_list)
    DAM_list = lapply (1:ncol (DAM_list), function(x) 
      {
      df = DAM_list[,x]  
      df = do.call (cbind, (assays(df)))
      colnames(df) = names (assays(DAM_list))
      df$gene = rowData (DAM_list)$name
      df
      })
    names (DAM_list) = listnames
    saveRDS (DAM_list, paste0 ('DAM_',metaGroupName,'.rds'))    
    } else {
    DAM_list = readRDS (paste0('DAM_',metaGroupName,'.rds'))
    }
  
  active_genes = corGSM_MM$MotifMatrix_name[corGSM_MM$cor > 0.1]
  DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% active_genes,])
  names (DAM_list2) = names (DAM_list)
  FDR_threshold = 1e-3
  meandiff_threshold = 0
  top_genes = 3
  DAM_top_list = DAM_list2[sapply (DAM_list2, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$MeanDiff) > meandiff_threshold,]) > 0)]
  DAM_top_list = lapply (seq_along(DAM_top_list), function(x) {
    res = DAM_top_list[[x]]
    #res = na.omit (res)
    res = res[res$FDR < FDR_threshold,]
    res = res[order (res$FDR), ]
    res = res[res$MeanDiff > meandiff_threshold,]
    res$comparison = names(DAM_top_list)[x]
    if (nrow(res) < top_genes) 
      {
      res
      } else {
      head (res,top_genes)
      }
    })
  DAM_df = Reduce (rbind ,DAM_top_list)
  
  
  DAM_df$gene = gsub ('_.*','',DAM_df$gene)
  DAM_df$gene = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", DAM_df$gene)
  mMat = assays (mSE)[[1]]
  rownames (mMat) = rowData (mSE)$name
  mMat_mg = mMat[DAM_df$gene, ]
  mMat_mg = as.data.frame (t(mMat_mg))
  mMat_mg$metaGroup = as.character (archp@cellColData[,metaGroupName])
  mMat_mg = aggregate (.~ metaGroup, mMat_mg, mean)
  rownames (mMat_mg) = mMat_mg[,1]
  mMat_mg = mMat_mg[,-1]
  mMat_mg = mMat_mg[names (DAM_list),]
  #mMat_mg = mMat_mg[names(table (archp@cellColData[,metaGroupName])[table (archp@cellColData[,metaGroupName]) > 50]),]
  DAM_hm = Heatmap (t(scale(mMat_mg)), 
          row_labels = colnames (mMat_mg),
          column_title = paste('top',top_genes),
          clustering_distance_columns = 'euclidean',
          clustering_distance_rows = 'euclidean',
          cluster_rows = F,
          #col = pals_heatmap[[5]],
          cluster_columns=F,#col = pals_heatmap[[1]],
          row_names_gp = gpar(fontsize = 8, fontface = 'italic'),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          name = 'chromVAR',
          #rect_gp = gpar(col = "white", lwd = .5),
          border=TRUE,
          col = rev(palette_deviation)

          #right_annotation = motif_ha
          )

  #DAG_grob = grid.grabExpr(draw(DAG_hm, column_title = 'DAG GeneScore2', column_title_gp = gpar(fontsize = 16)))
pdf (file.path ('Plots',paste0('DAM_clusters_',metaGroupName,'_heatmaps.pdf')), width = 3, height = 5)
print(DAM_hm)
dev.off()

