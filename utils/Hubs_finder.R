#' hubs_finder
#'
#' @param ArchRProj 
#' @param projdir 
#' @param group_by Find Hubs within group type
#' @param min_cells 
#' @param max_gap Max link width (bp)
#' @param cor_cutoff Correlation cut-off
#' @param cor_FDR 
#' @param cor_var 
#' @param max_hub_width Max hub width
#' @param min_peaks Min peaks per hub
#' @param macs_score MACS2 score cut-off
#' @param dgs Distance from start site
#'
#' @return
#' @export
#'
#' @examples
hubs_finder = function (
  ArchRProj = NULL,
  projdir = getOutputDirectory (ArchRProj),
  group_by = NULL,
  select_group = NULL,
  min_cells = 50,
  cor_cutoff = 0.4,
  cor_FDR = 1e-10,
  cor_var = 0.3,
  min_peaks = 3,
  macs_score = 1,
  dgs = 2000,
  cores = 1,
  remove_chr = c('chrX')
  )
  {
  require (parallel)  
  cat (paste(
    '--* Hubs Finder *--','\n',
    'Find Hubs in',group_by,'\n',
    'Correlation cut-off:', cor_cutoff,'\n',
    'Min peaks per hub:',min_peaks,'\n',
    'MACS2 score cut-off:',macs_score,'\n',
    'Distance from start site:',dgs,'\n',
    'Number of cores:',cores,'\n')
    )
  
  if (!is.null(group_by))
    {
    cell_meta = as.data.frame (getCellColData (ArchRProj, group_by))[,1]
    groups = names (table(cell_meta) [table (cell_meta) > min_cells])
    groupsfiles = unlist(sapply (groups, function(x) list.files (file.path(projdir,'PeakCalls'))[grepl(paste0('^',x,'-'),list.files (file.path(projdir,'PeakCalls')))]))
    groups = groups[groups %in% names(groupsfiles)]
    message (paste('run hub finder for groups:', groups))
    if (!is.null(select_group)) groups = groups[groups %in% select_group]
    main_peakSet = getPeakSet (ArchRProj)
    main_peakSet$idx = 1:length (main_peakSet) 
    hubsClustersL = list()
    hubs_linksL = list()
  #fragments = unlist(getFragmentsFromProject (
  # ArchRProj = ArchRProj))  
  
    for (cluster in groups) 
      {
      message(paste('Running cell group:',cluster))  
      
      # get peak calls and fragments per cluster
      cluster_peakSet = readRDS (paste0(projdir,'/PeakCalls/',cluster,'-reproduciblePeaks.gr.rds')) 
      
      # Filter peak set 
      cluster_peakSet = cluster_peakSet [cluster_peakSet$distToGeneStart > dgs]
      cluster_peakSet = cluster_peakSet [cluster_peakSet$score > macs_score,]
      cluster_to_main_idx = subjectHits(findOverlaps(cluster_peakSet, main_peakSet))
      if (length(cluster_to_main_idx) < 1) 
        {
        message ('no peaks passed filters!')
        next
        }
      
      cA1 = getCoAccessibility (
        ArchRProj = ArchRProj,
        corCutOff = cor_cutoff,
        resolution = 1,
        returnLoops = FALSE)
        
      # Filter peak-pairs correlation
      cA1 <- cA1[cA1$FDR < cor_FDR,]
      if (nrow(cA1) < 1) 
        {
        message ('no coaccessed peaks passed FDR filter!')
        next
        }
      cA1 <- cA1[rowMins(cbind(cA1$VarQuantile1,cA1$VarQuantile2)) > cor_var,]
      if (nrow(cA1) < 1) 
        {
        message ('no coaccessed peaks passed cor variance filter!')
        next
        }
      cA1 = as.data.frame (cA1, row.names=NULL)
      cA1 = cA1[cA1$queryHits %in% cluster_to_main_idx,] # subset for cluster-specific peaks 
      cA1 = cA1[cA1$subjectHits %in% cluster_to_main_idx,] # subset for cluster-specific peaks 
    
      # Convert edge list to graph and select components (hubs)
      edge_list = cA1[,c('queryHits','subjectHits','correlation')]
      edge_list = as.data.frame (edge_list)
      colnames (edge_list)[3] = 'weight'
      hubs_grp = igraph::graph.data.frame (edge_list, directed = FALSE)
      components_hubs = igraph::components (hubs_grp)
      clusters_size = table (components_hubs$membership)
      clusters_size = clusters_size[order(clusters_size, decreasing=T)]
      
      # Filter by minimun number of peaks in hubs
      clusters_size = clusters_size[clusters_size >= min_peaks]
      cluster_ids = names (clusters_size)
      hubs = mclapply (cluster_ids,function(x){
        names(components_hubs$membership)[components_hubs$membership == x]
      }, mc.cores = cores)
      if (length(hubs) == 0) 
        {
        message('No hubs passing min_peaks filter!')
        next  
        }
      # Replace peak index with peak regions per hub
      hubs_exp = mclapply (hubs, function(x) {
      main_peakSet [as.numeric(x)] 
      }, mc.cores=cores) # error when more than 1 core is set
      names (hubs_exp) = 1:length(hubs_exp)
      hubsClustersL[[cluster]] = mclapply (hubs_exp, function(x) {
        as.data.frame(x, row.names=NULL)}, mc.cores=cores)
      # Make list of hub links
      hubs_eg = edge_list[edge_list$queryHits %in% unlist(hubs) | edge_list$subjectHits %in% unlist(hubs), ] 
      hubs_summits_1 = as.data.frame(resize (main_peakSet [hubs_eg$queryHits,], 1, "center"), row.names=NULL)[,c(1,2)]
      hubs_summits_2 = as.data.frame(resize (main_peakSet [hubs_eg$subjectHits,], 1, "center"),, row.names=NULL)[,c(1,2)]
      colnames(hubs_summits_2)[2] = 'end'
      hubs_links = data.frame (hubs_summits_1, end = hubs_summits_2[,2], value = hubs_eg$weight)
      hubs_links = transform (hubs_links, end = ifelse(end < start, start, end), start = ifelse(start > end, end, start))
      hubs_linksL[[cluster]] = makeGRangesFromDataFrame (hubs_links, keep.extra.columns=T)
      message (paste (length(hubsClustersL[[cluster]]), 'Hubs detected in',cluster))
      }
    cat (paste(
    'Hubs detected for cluster','\n'))
    print (data.frame (
      cluster = names(hubsClustersL), 
      hubs = sapply (hubsClustersL, length)))

    message('Format hubs list for downstream analyses')
    message ('1 - Unlist hubs from clusters')
    hubsMerged = do.call (c, hubsClustersL)
  
    # Find duplicated hubs
    message ('2 - Remove duplicated hubs')
    hub_dup = sapply (hubsMerged, function(x) paste(paste0(x$seqnames, x$start, x$end),collapse=""))
    hubsMerged = hubsMerged[!duplicated(hub_dup)]
    message (paste(sum (duplicated(hub_dup)), 'duplicated hubs removed'))

    if (!is.null(remove_chr)) 
      {
      hubsMerged = hubsMerged[!sapply (hubsMerged, function(x) any(x$seqnames %in% remove_chr))]
      hubs_linksL = lapply (hubs_linksL, function(x) x[!seqnames(x) %in% remove_chr])
      }

    # Collapse hubs
    message ('3 - collapse Hubs by genomc ranges')
    hubs_cluster_collapsed = mclapply (hubsMerged, function(x){
      data.frame (seqnames = x$seqnames[1], start=x$start[1], end=x$end[nrow(x)],
      gene = paste(unique(x$nearestGene), collapse= '-'))
      },mc.cores = cores)
    hubs_cluster_collapsed = do.call (rbind, hubs_cluster_collapsed)
    hubs_cluster_collapsed = makeGRangesFromDataFrame (hubs_cluster_collapsed, keep.extra.columns=T)
    
    message ('4 - merge all peaks from hubs')
    peaksMerged = do.call (rbind, hubsMerged)
    rownames (peaksMerged) = NULL
    peaksMerged = makeGRangesFromDataFrame (peaksMerged, keep.extra.columns=T)
    hubs_obj = list(
      hubsClusters = hubsClustersL,
      peaksMerged = peaksMerged,
      hubsMerged = hubsMerged, 
      hubsCollapsed = hubs_cluster_collapsed, 
      peakLinks = hubs_linksL,
      hubs_id = paste0('HUB',1:length(hubs_cluster_collapsed))
      )        
    return (hubs_obj)

    } else {
    hubsClustersL = list()
    hubs_linksL = list()  
    main_peakSet = getPeakSet (ArchRProj)  

    cA1 = getCoAccessibility (
        ArchRProj = ArchRProj,
        corCutOff = cor_cutoff,
        resolution = 1,
        returnLoops = FALSE)
        
      # Filter peak-pairs correlation
      message ('Filter based on FDR and variance')
      cA1 <- cA1[cA1$FDR < cor_FDR,]
      if (nrow(cA1) < 1) 
        {
        message ('no coaccessed peaks passed FDR filter!')
        next
        }
      cA1 <- cA1[rowMins(cbind(cA1$VarQuantile1,cA1$VarQuantile2)) > cor_var,]
      if (nrow(cA1) < 1) 
        {
        message ('no coaccessed peaks passed cor variance filter!')
        next
        }
      cA1 = as.data.frame (cA1, row.names=NULL)  
      # Convert edge list to graph and select components (hubs)
      message ('Convert edge list to graph')
      edge_list = cA1[,c('queryHits','subjectHits','correlation')]
      edge_list = as.data.frame (edge_list)
      colnames (edge_list)[3] = 'weight'
      hubs_grp = igraph::graph.data.frame (edge_list, directed = FALSE)
      components_hubs = igraph::components (hubs_grp)
      clusters_size = table (components_hubs$membership)
      clusters_size = clusters_size[order(clusters_size, decreasing=T)]
      
      # Filter by minimun number of peaks in hubs
      message ('filter by minimun number of peaks in hubs')
      clusters_size = clusters_size[clusters_size >= min_peaks]
      cluster_ids = names (clusters_size)
      message ('Get peak indexes from each igraph component')
      hubs = mclapply (cluster_ids,function(x){
        names(components_hubs$membership)[components_hubs$membership == x]
      }, mc.cores = cores)
      if (length(hubs) == 0) 
        {
        message('No hubs passing min_peaks filter!')
        next  
        }
      # Replace peak index with peak regions per hub
      message ('Replace peak indexes with peak regions per hub')
      hubs_exp = mclapply (hubs, function(x) {
      main_peakSet [as.numeric(x)] 
      }, mc.cores=cores) # error when more than 1 core is set
      names (hubs_exp) = 1:length(hubs_exp)
      hubsClustersL[['main']] = mclapply (hubs_exp, function(x) {
        as.data.frame(x, row.names=NULL)}, mc.cores=cores)
      # Make list of hub links
      message ('Generate hub links Granges')
      hubs_eg = edge_list[edge_list$queryHits %in% unlist(hubs) | edge_list$subjectHits %in% unlist(hubs), ] 
      hubs_summits_1 = as.data.frame(resize (main_peakSet [hubs_eg$queryHits,], 1, "center"), row.names=NULL)[,c(1,2)]
      hubs_summits_2 = as.data.frame(resize (main_peakSet [hubs_eg$subjectHits,], 1, "center"),, row.names=NULL)[,c(1,2)]
      colnames(hubs_summits_2)[2] = 'end'
      hubs_links = data.frame (hubs_summits_1, end = hubs_summits_2[,2], value = hubs_eg$weight)
      hubs_links = transform (hubs_links, end = ifelse(end < start, start, end), start = ifelse(start > end, end, start))
      hubs_linksL[['main']] = makeGRangesFromDataFrame (hubs_links, keep.extra.columns=T)
      message (paste (length(hubsClustersL[['main']]), 'Hubs detected'))
      }
    
    message('Format hubs list')
    
    hubsMerged = hubsClustersL[['main']]
    if (!is.null(remove_chr)) 
      {
      hubsMerged = hubsMerged[!sapply (hubsMerged, function(x) any(x$seqnames %in% remove_chr))]
      hubs_linksL = lapply (hubs_linksL, function(x) x[!seqnames(x) %in% remove_chr])
      }
    
    # Collapse hubs
    message ('1 - collapse Hubs by genomc ranges')
    hubs_cluster_collapsed = mclapply (hubsMerged, function(x){
      data.frame (seqnames = x$seqnames[1], start=x$start[1], end=x$end[nrow(x)],
      gene = paste(unique(x$nearestGene), collapse= '-'))
      },mc.cores = cores)
    hubs_cluster_collapsed = do.call (rbind, hubs_cluster_collapsed)
    hubs_cluster_collapsed = makeGRangesFromDataFrame (hubs_cluster_collapsed, keep.extra.columns=T)
    
    message ('2 - merge all peaks from hubs')
    peaksMerged = do.call (rbind, hubsMerged)
    rownames (peaksMerged) = NULL
    peaksMerged = makeGRangesFromDataFrame (peaksMerged, keep.extra.columns=T)
    hubs_obj = list(
      hubsClusters = hubsClustersL,
      peaksMerged = peaksMerged,
      hubsMerged = hubsMerged, 
      hubsCollapsed = hubs_cluster_collapsed, 
      peakLinks = hubs_linksL,
      hubs_id = paste0('HUB',1:length(hubs_cluster_collapsed))
      )
    message('hubs object created!')    
    return (hubs_obj)
    }

  #return (hubs_cluster)


