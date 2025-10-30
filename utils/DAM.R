
# Find DAM ####
### ChromVAR based analysis ####

DAM = function (
  ArchRProj = NULL,
  metaGroupName = NULL,
  FDR_threshold = 1e-2,
  meandiff_threshold = 0,
  top_genes=5,
  filter_by_scRNA=TRUE, # use 'srt' object
  min_exp=.1,
  force = FALSE) {

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
  
  # Clean TF names
  DAM_list = lapply (DAM_list, function(x)
       {
       x$gene = gsub ('_.*','',x$gene)
       x$gene = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", x$gene)
       x
       })
  

  if (filter_by_scRNA)
  {
  #Get active genes from RNA
  ps = log2(as.data.frame (AverageExpression (srt, 
  features = sapply (unique(unlist(lapply(DAM_list, function(x) x$gene))), function(x) unlist(strsplit (x, '_'))[1]), 
  group.by = metaGroupName)[[1]]) +1)
  ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
  active_TFs = rownames(ps)

  #active_genes = corGSM_MM$MotifMatrix_name[corGSM_MM$cor > 0.1]
  DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% active_TFs,])    
  } else {
  DAM_list2 = DAM_list  
  }
  names (DAM_list2) = names (DAM_list)
  DAM_top_list = DAM_list2[sapply (DAM_list2, function(x) nrow (x[x$FDR < FDR_threshold & abs(x$MeanDiff) > meandiff_threshold,]) > 0)]
  DAM_top_list = lapply (seq_along(DAM_top_list), function(x) {
    res = DAM_top_list[[x]]
    res = na.omit (res)
    res = res[res$FDR < FDR_threshold,]
    res = res[order (res$FDR), ]
#    res = res[order (-res$MeanDiff), ]
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
return (DAM_df)
}  

