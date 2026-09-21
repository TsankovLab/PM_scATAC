#### Load hidden functions from ArchR ####
hidden_func_names = ls(getNamespace("ArchR"), all.names = TRUE)[grep ('^\\.', ls(getNamespace("ArchR"), all.names = TRUE))]
for (func_name in hidden_func_names) {
  assign(func_name, getFromNamespace(func_name, "ArchR"), envir = .GlobalEnv)
}

# Retreive matrix from ArchR object, fix barcode order and names
fetch_mat = function(archp = NULL, mat = 'Motif')
  {
  mat = ArchR::getMatrixFromProject (archp, useMatrix = paste0(mat,'Matrix'))
  mat = mat[, archp$cellNames]
  rowData(mat)$name = gsub ('_.*','',rowData(mat)$name)
  rowData(mat)$name = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", rowData(mat)$name)
  return (mat)
  }



pbknn = function(SE = NULL, KNN)
  {
  Mat = assay (SE)
  rownames (Mat) = rowData (SE)$name
  message ('generate pseudobulk KNN')
  pseudobulk_knn = do.call (rbind, lapply (KNN, function(x) rowMeans (Mat[,x])))
  return (pseudobulk_knn)
  }
  
hubMatsGen = function(hubs_obj, KNN, projdir_hubs, fragments, force = FALSE)
{
###--- Generate matrces of hubs x cells / knn from collapsed non-redundant hubs to identify differential hubs  ---###  
  # Generate matrix of fragment counts of peaks x cluster
    message ('Generate matrix of fragment counts of peaks x KNN')
    peaksKnn_mat = matrix (ncol = length(KNN), nrow = length(hubs_obj$peaksMerged))
    pb =progress::progress_bar$new(total = length (KNN)) 
    for (i in 1:length(KNN)) 
      {
      pb$tick()  
      peaks_fr = fragments[fragments$RG %in% KNN[[i]]]  
      fragments_hubs = countOverlaps (hubs_obj$peaksMerged, peaks_fr)
      peaksKnn_mat[,i] = fragments_hubs
      }
    peaksKnn_mat = as.data.frame (peaksKnn_mat)
    peaksKnn_mat$hubs_id = paste0(rep (hubs_obj$hubs_id,sapply(hubs_obj$hubsMerged,nrow)))
    
    # Aggregate peaks x knn matrix to hub level taking mean
    #gm_mean = function (x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) } # geometric mean, doesnt work as good as normal mean i think
    message ('Aggregate peaks x knn matrix to hub level taking mean')
    hubsKnn_mat = aggregate (. ~ hubs_id, data = peaksKnn_mat, FUN = mean)
    #rownames (hubsKnn_mat) = apply (data.frame(as.data.frame (hubs_obj$hubsCollapsed), hubs_id = hubs_obj$hubs_id)[,c(7,1,2,3,6)], 1,
    # function(x) paste(x, collapse = '_')) [match (hubsKnn_mat$hubs_id, hubs_obj$hubs_id)]
    #rownames (hubsKnn_mat) = gsub (' ','', rownames(hubsKnn_mat))
    rownames (hubsKnn_mat) = hubsKnn_mat$hubs_id
    hubsKnn_mat = hubsKnn_mat[,-1]
    colnames (hubsKnn_mat) = sub ('KNN.*','',names(KNN)) # add column names
    
    # Normalise by seq depth
    message ('Normalise by ReadsInTSS')
    cellpool_nFrags = sapply(KNN,function(v) sum(archp$nFrags[archp$cellNames %in% v]))
    hubsKnn_mat = t(t(hubsKnn_mat) * (10^6 /cellpool_nFrags)) # scale

    message ('Aggregate Hubs x knn to clusters')  
    hubsClust_mat = as.data.frame (t (hubsKnn_mat))
    hubsClust_mat$metaGroup = colnames (hubsKnn_mat)
    hubsClust_mat = aggregate (.~ metaGroup, data = hubsClust_mat, FUN = mean)
    rownames (hubsClust_mat) = hubsClust_mat[,1]
    hubsClust_mat = hubsClust_mat[,-1]
return (list (hubsKnnMat = hubsKnn_mat, hubsClusterMat = hubsClust_mat))
}

findDiffHubs = function (hubsKnn_mat, hubs_obj, tt_group1, tt_group2, projdir_hubs)
  {
  # Run pairwise Differential Accessible Hubs (DAH)
  #padj_threshold = 0.01
  DAH_l = list()
  message ('Run hub differential accessibility')  
  for (test in seq_along (tt_group1))
    {
    test_comparison = paste0 (paste(tt_group1[[test]],collapse='_'), '--', paste(tt_group2[[test]], collapse='_'))
    message (paste ('Running ttest for', test_comparison))
    # Run ttest on mesothelium vs fibroblasts
    t_test = NULL
    #hubsKnn_mat2 = hubsKnn_mat[!grepl ('chrX', rownames(hubsKnn_mat)),]
    hubsKnn_mat2 = hubsKnn_mat
    for (x in 1:nrow(hubsKnn_mat2))
    {
    if (var (unlist(hubsKnn_mat2[x,colnames(hubsKnn_mat2) %in% tt_group1[[test]]])) == 0 & 
      var (unlist(hubsKnn_mat2[x,colnames(hubsKnn_mat2) %in% tt_group2[[test]]])) == 0)
      {
      t_test[x] = NA
      } else {
      t_test[x] = t.test (unlist(hubsKnn_mat2[x,colnames(hubsKnn_mat2) %in% tt_group1[[test]]]),
      unlist(hubsKnn_mat2[x,colnames(hubsKnn_mat2) %in% tt_group2[[test]]]))$p.value
      }
    }
    LogFC = sapply (1:nrow(hubsKnn_mat2), function (x)  mean(unlist(log2(hubsKnn_mat2[x,colnames(hubsKnn_mat2) %in% tt_group1[[test]]]+1)))
    - mean(unlist(log2(hubsKnn_mat2[x,colnames(hubsKnn_mat2) %in% tt_group2[[test]]]+1))))
    maxMean = sapply (1:nrow(hubsKnn_mat2), function (x)  max(mean(unlist(hubsKnn_mat2[x,colnames(hubsKnn_mat2) %in% tt_group1[[test]]])),
    mean(unlist(hubsKnn_mat2[x,colnames(hubsKnn_mat2) %in% tt_group2[[test]]]))))
    padj = p.adjust(t_test, method = 'fdr', n = length(t_test))
     
    DAH_df = data.frame (
      hub_id= rownames (hubsKnn_mat2), 
      region= as.character(granges(hubs_obj$hubsCollapsed))[match (rownames(hubsKnn_mat2), hubs_obj$hubs_id)],
      padj = padj, 
      log2FC= LogFC,
      maxMean = maxMean, 
      gene = hubs_obj$hubsCollapsed$gene[match (rownames(hubsKnn_mat2), hubs_obj$hubs_id)],
      comparison = test_comparison
      )
    DAH_df = DAH_df[!is.na(DAH_df$padj),]
    DAH_df = DAH_df[order (DAH_df$log2FC, decreasing=T),]
    DAH_l[[test_comparison]] = DAH_df
    #message ('print csv of signficant hubs')  
    #write.csv (DAH_df, paste0(projdir_hubs, 'ttest_',tt_group1[test],'_vs_',tt_group2[test],'_global.csv'))
    }
   message ('done!')   
   return (DAH_l) 
   } 

peakToTF = function (regions, pSE, mSE, KNN, cores)
  {  
  pMat = assay (pSE) 
  pMat = pMat [,unlist(KNN)]
  mMat = assay (mSE)
  mMat = mMat[, unlist(KNN)]
  rownames (mMat) = sapply (rownames(mMat), function(x) unlist (strsplit (x, '_'))[1])
  message (paste('Check order of barcode match between pMat and mMat',all(colnames(pMat) == colnames (mMat))))
  message ('Factorize pMat by regions')
  pMatL = mclapply (seq_along(regions), function(x) pMat[queryHits (findOverlaps (rowRanges (pSE), regions[x])),,drop=F], mc.cores=cores)
  pMat_df = do.call (rbind, mclapply (pMatL, function (x)
    {
    x = as.data.frame (t(x))
    x$cellGroup = rep (names(KNN), sapply(KNN, length))
    x = aggregate (.~cellGroup, x, mean)
    rownames(x) = x[,1]
    x = x[,-1, drop=F]
    x = colMeans(t(x))
    }, mc.cores=cores))
  mMat = as.data.frame (t(mMat[,unlist(KNN)]))
  mMat$cellGroup = rep (names(KNN), sapply(KNN, length))
  mMat = aggregate (.~cellGroup, mMat, mean)
  rownames(mMat) = mMat[,1]
  mMat = mMat[,-1]
  message ('Compute correlation between TF deviations and aggregated peaks')  
  return (cor (mMat, t(pMat_df)))
  }


# Compute correlation of deviations vs genescore. Need an archr object KNN groups and a metaGroupName 
activeTF = function (
  archp = NULL, 
  knns = NULL, 
  metaGroupName = NULL, 
  motifSE = NULL, 
  #useMotifmatrix = 'Motif', 
  geneScoreSE = NULL,
  subsample = NULL,
  iter = 1000
  )
    {
    #if (!any (ls() == 'gsSE')) gsSE = getMatrixFromProject (archp, useMatrix = 'GeneScoreMatrix')
    #gsSE = gsSE[, archp$cellNames]
    #if (!any (ls() == 'mSE')) mSE = getMatrixFromProject (archp, useMatrix = 'MotifMatrix')
    #mSE = mSE[, archp$cellNames]
    #rowData(motifSE)$name = gsub ('_.*','',rowData(motifSE)$name) 
    #rowData(motifSE)$name = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", rowData(motifSE)$name)
    if (!is.null(subsample))
      {
      subsample_groups = lapply (seq(iter), function(z) 
        {
        do.call (c, unlist(lapply (unique (archp@cellColData[,metaGroupName]), function(x) 
          {
          if (length (names(knns)[grep(paste0('^',x), names(knns))]) >= subsample)
            subsampled_knn = sample(names(knns)[grep(paste0('^',x), names(knns))], subsample)    
            else subsampled_knn = names(knns)[grep(paste0('^',x), names(knns))]             
          knns[subsampled_knn]  
          })))
        })      
      }

    gsMat = assay (geneScoreSE)
    gsMat = gsMat[,rownames(archp@cellColData)]
    rownames (gsMat) = rowData (geneScoreSE)$name
    mMat = assay (motifSE)
    mMat = mMat[,rownames(archp@cellColData)]
    rownames (mMat) = rowData (motifSE)$name
    gene_int = intersect(rownames(gsMat) , rownames(mMat))
    gsMat = gsMat[gene_int, ]
    mMat = mMat [gene_int, ]

    message ('Compute correlation deviations vs genescore')
    mAvknn_res = lapply (subsample_groups, function(y) 
      {
      mAvknn = lapply (y, function(x) rowMeans (mMat[,x])) 
      do.call (cbind, mAvknn)
      })
    gsAvknn_res = lapply (subsample_groups, function(y) 
      {
      gsAvknn = lapply (y, function(x) rowMeans (gsMat[,x]))
      do.call (cbind, gsAvknn)
      })
    tf_gs_cor = as.data.frame (do.call (cbind, lapply (seq_along(mAvknn_res), function(x) sapply (gene_int, function(y) cor (mAvknn_res[[x]][y,], gsAvknn_res[[x]][y,])))))
    tf_gs_cor$gene = rownames(tf_gs_cor)
    tf_gs_cor = gather (tf_gs_cor, iteration, correlation, 1:(ncol(tf_gs_cor) - 1))  
    tf_gs_cor$correlation[is.na(tf_gs_cor$correlation)] = 0

    tf_gs_cor_random = as.data.frame (do.call (cbind, lapply (seq_along(mAvknn_res), function(x) sapply (gene_int, function(y) cor (mAvknn_res[[x]][y,sample(colnames(mAvknn_res[[x]]),ncol(mAvknn_res[[x]]))], gsAvknn_res[[x]][y,sample(colnames(gsAvknn_res[[x]]),ncol(gsAvknn_res[[x]]))])))))
    tf_gs_cor_random$gene = rownames(tf_gs_cor_random)
    tf_gs_cor_random = gather (tf_gs_cor_random, iteration, correlation, 1:(ncol(tf_gs_cor_random) - 1))  
    tf_gs_cor_random$correlation[is.na(tf_gs_cor_random$correlation)] = 0
      
    mAvknn_mean = as.data.frame (do.call (cbind, lapply (gene_int, function(x) rowMeans(sapply (seq_along(mAvknn_res), function(y) mAvknn_res[[y]][x,])))))
    mAvknn_min = as.data.frame (do.call (cbind, lapply (gene_int, function(x) rowMin(sapply (seq_along(mAvknn_res), function(y) mAvknn_res[[y]][x,])))))
    mAvknn_max = as.data.frame (do.call (cbind, lapply (gene_int, function(x) rowMax(sapply (seq_along(mAvknn_res), function(y) mAvknn_res[[y]][x,])))))
    gsAvknn_mean = as.data.frame (do.call (cbind, lapply (gene_int, function(x) rowMeans(sapply (seq_along(gsAvknn_res), function(y) gsAvknn_res[[y]][x,])))))
    gsAvknn_min = as.data.frame (do.call (cbind, lapply (gene_int, function(x) rowMin(sapply (seq_along(gsAvknn_res), function(y) gsAvknn_res[[y]][x,])))))
    gsAvknn_max = as.data.frame (do.call (cbind, lapply (gene_int, function(x) rowMax(sapply (seq_along(gsAvknn_res), function(y) gsAvknn_res[[y]][x,])))))
    colnames (mAvknn_mean) = gene_int
    colnames (mAvknn_min) = gene_int
    colnames (mAvknn_max) = gene_int
    colnames (gsAvknn_mean) = gene_int
    colnames (gsAvknn_min) = gene_int
    colnames (gsAvknn_max) = gene_int

    # Compute max TF Motif Delta (from ArchR)
    # seGroupMotif = getGroupSE(ArchRProj = archp, useMatrix = useMotifmatrix, groupBy = metaGroupName) 
    # seZ = seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
    # rowData(seZ)$maxDelta = lapply(seq_len(ncol(seZ)), function(x)
    #   {
    #   rowMaxs(assay(seZ) - assay(seZ)[,x])
    #   }) %>% Reduce("cbind", .) %>% rowMaxs
    # rowData(seZ)$name = gsub ('_.*','',rowData(seZ)$name)  
    # rowData(seZ)$name = gsub ("(NKX\\d)(\\d{1})$","\\1-\\2", rowData(seZ)$name)
    # tf_gs_cor$delta = rowData(seZ)[match(rownames(tf_gs_cor), rowData(seZ)$name), "maxDelta"]
    return (list (
      cor_df = tf_gs_cor, 
      motifAvKnn_mean = mAvknn_mean, 
      motifAvKnn_min = mAvknn_min,
      motifAvKnn_max = mAvknn_max,
      gsAvKnn_mean = gsAvknn_mean,
      gsAvKnn_min = gsAvknn_min,
      gsAvKnn_max = gsAvknn_max,
      permuted_cor_df = tf_gs_cor_random
      ))
    }


sigActiveTF = function (motifSE, geneScoreSE, KNN, nperm = 10000)
  {
  message('Harmonize motifSE and genescoreSE')
  gsMat = assay (geneScoreSE)
  rownames (gsMat) = rowData (geneScoreSE)$name
  mMat = assay (motifSE)
  rownames (mMat) = rowData (motifSE)$name
  gene_int = intersect(rownames(gsMat) , rownames(mMat))
  gsMat = gsMat[gene_int, ]
  mMat = mMat [gene_int, ]

  message ('Average deviations and genescores per KNN and correlate')
  mAvknn_res = do.call (rbind, lapply (KNN, function(x) rowMeans (mMat[,x])))
  gsAvknn_res = do.call (rbind, lapply (KNN, function(x) rowMeans (gsMat[,x])))
  tf_gs_cor = sapply (gene_int, function(x) cor (mAvknn_res[,x], gsAvknn_res[,x]))
  tf_gs_cor[is.na (tf_gs_cor)] = 0
  
  message('Permutation test')
  tf_gs_cor_perm=list()
  library (progress)
  pb = progress::progress_bar$new(total = nperm)
  for(i in 1:nperm)
    {
    pb$tick()   
    tf_gs_cor_perm[[i]] = sapply (gene_int, function(x) 
      cor (mAvknn_res[sample(seq(nrow(mAvknn_res)),nrow(mAvknn_res), replace=F),x],
          gsAvknn_res[sample(seq(nrow(gsAvknn_res)),nrow(gsAvknn_res), replace=F),x]))
    } 
  tf_gs_cor_perm_df = do.call (cbind, tf_gs_cor_perm)  
    #tf_gs_cor_perm = gather (tf_gs_cor_random, iteration, correlation, 1:(ncol(tf_gs_cor_random) - 1))  
    
  #  tf_gs_cor_perm$correlation[is.na(tf_gs_cor_random$correlation)] = 0
  result=rowSums(abs(tf_gs_cor_perm_df) >= abs(tf_gs_cor))/(nperm)
  return(list(pvalue = result[order(result)], correlations = tf_gs_cor, permutations = tf_gs_cor_perm_df, mknn = mAvknn_res, gsknn = gsAvknn_res))
  }  
  
  # test1 <- permutation.test(treatment, outcome, 10000)
  # hist(test1[[2]], breaks=50, col='grey', main="Permutation Distribution", las=1, xlab='')
  # abline(v=original, lwd=3, col="red")  



# Survival analysis function that accept a list of genes
SurV = function (surv_genes, counts_mat, Time, Status, clinical, thr='Median')
  {
  require(survival)
  require(survminer)
  k = 1
  passed_genes = NULL
  fit = list()
  pval_survival = NULL
  direction = NULL
  counts_mat = counts_mat[ rownames(counts_mat) %in% surv_genes,]
  for (i in 1:nrow(counts_mat))
    {
      clinical$gene_exp = unname(unlist(counts_mat[i,, drop=T]))
    
      lowExpr = summary(clinical$gene_exp)[thr]
      clinical$gene_bin = ifelse(clinical$gene_exp > lowExpr, 'High','Low')
      clinical$Time = Time      
     fit[[k]] = survfit (Surv (Time, Status) ~ gene_bin, clinical)
     if (dim(fit[[k]]) < 2) {
      print (paste0('Not enough obs in either group for gene:', rownames(counts_mat)[i]))
     next
     }
     diff = survdiff(formula = Surv(Time, Status) ~ gene_bin, data = clinical)
     pval_survival[k] = pchisq(diff$chisq, length(diff$n)-1, lower.tail = FALSE)
     direct_tmp = surv_median(fit[[k]])[which.min(surv_median(fit[[k]])[,2]),'strata']
     if (length(direct_tmp) == 0) {
      direction[k] = '-'
      } else {
     direction[k] = direct_tmp
      }
     passed_genes[k] = rownames(counts_mat)[i]
     k = k + 1
     }
  return (data.frame (gene = passed_genes,
  pval = pval_survival,
  direction = direction, stringsAsFactors=FALSE))
  }


# Function returns the Jaccard index and Jaccard distance
jaccard <- function(df, margin) {
  if (margin == 1 | margin == 2) {
    M_00 <- apply(df, margin, sum) == 0
    M_11 <- apply(df, margin, sum) == 2
    if (margin == 1) {
      df <- df[!M_00, ]
      JSim <- sum(M_11) / nrow(df)
    } else {
      df <- df[, !M_00]
      JSim <- sum(M_11) / length(df)
    }
    JDist <- 1 - JSim
    return(c(JSim = JSim, JDist = JDist))
  } else break
}

# Palette centered on 0 for correlation heatmaps
require (circlize)
pal_corr1 = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
pal_corr2 = colorRamp2(c(-0.5, 0, 0.5), c("green", "white", "red"))

# Co-accessibility networks 

# take a GRange object and collapse to one element with ranges spanning all 
# the previous ones
max_range = function (granges_obj,upstream = 20000, downstream=20000, grange=TRUE)
        {
         if (grange == TRUE) {
            df = as.data.frame (granges_obj, row.names=NULL)
            } else {
            df = granges_obj
            }
         min_start = min (df$start)
         min_start = min_start - upstream
         max_end = max (df$end)
         max_end = max_end + downstream
         seqnames = unique(df$seqnames)
         cluster = df$GroupReplicate 
        return (GRanges(seqnames = seqnames,
                  ranges = IRanges(min_start,max_end),
                   length=max_end - min_start
                  ))
        }

# extend ranges in GRanges object for plotTracks in ArchR
ext_range = function (granges_obj,upstream = 20000, downstream=20000)
        {
         df = as.data.frame (granges_obj, row.names=NULL)
         df$start = df$start - upstream
         df$end = df$end + downstream
        return (GRanges(seqnames = df$seqnames,
                  ranges = IRanges(df$start,df$end)
                  ))
        } 

#### FUNCTIONS FOR CHORMOSOME INTERACTION MAPS ###
# Normalise correlation matrix by subtraction of background matrix      
rotate = function(x) t(apply(x, 2, rev))
diag_mean = function (mat)
  {
  #mat = rotate(mat)  
  k = 1
  j = ncol(mat) 
  diag_mean = NULL
  while (k  < (nrow(mat)))
    {
    diag_mean[k] = mean (diag(mat[k:nrow(mat),1:j]))    
    k = k + 1
    j = j - 1
    }
  diag_mean[k] = mat[nrow(mat),j]   
  return (diag_mean)
  }

# get bg matrix
get_bg_mat = function (mat, diag_mean)
  {
  k = 1
  #mat = rotate (mat)
  j = ncol (mat)
  j2 = nrow (mat) 
  bg_mat = matrix (nrow = nrow(mat), ncol = ncol(mat))
  while (k  < (nrow(mat)))
    {
    diag(bg_mat[k:nrow(mat),1:j]) = (diag(mat[k:nrow(mat),1:j])) - diag_mean[k]
    diag(bg_mat[1:j2,k:ncol(mat)]) = (diag(mat[k:nrow(mat),1:j])) - diag_mean[k]        
    k = k + 1
    j = j - 1
    j2 = j2 - 1
    }
  bg_mat[nrow(mat),1] = mat[nrow(mat),1] - diag_mean[k]   
  bg_mat[1,ncol(mat)] = mat[1,ncol(mat)] - diag_mean[k]
  #return (rotate(rotate(rotate(bg_mat))))
  return (bg_mat)
    }


### TADs ANALYSIS ###
# Run SpectralTADs on different matrices 
coordCols = function(x, resolution) {
colnames(x) = seq(0, (nrow(x)-1) * resolution, resolution)
return(x)
}


### FUNCTION FOR CNV ANALYSIS FROM GRANJA ####

#Estimating Copy Number Variation in scATAC-seq
#05/02/19
#Cite Satpathy*, Granja*, et al. 
#Massively parallel single-cell chromatin landscapes of human immune 
#cell development and intratumoral T cell exhaustion (2019)
#Created by Jeffrey Granja

"%ni%" <- Negate("%in%")

countInsertions <- function (query, fragments, by = "RG"){
  #Count By Fragments Insertions
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
  )
  by = 'RG'
  overlapDF <- DataFrame(findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Calculate Overlap Stats
  inPeaks <- table(overlapDF$name)
  total <- table(mcols(inserts)[, by])
  total <- total[names(inPeaks)]
  frip <- inPeaks / total
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1], 
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)), 
    dims = c(length(query), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  total <- total[colnames(sparseM)]
  frip <- frip[colnames(sparseM)]
  out <- list(counts = sparseM, frip = frip, total = total)
  return(out)
}

makeWindows <- function(genome, blacklist, windowSize = 10e6, slidingSize = 2e6,
  gene_level_annotation = TxDb.Hsapiens.UCSC.hg38.knownGene){
  chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
  chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
  windows <- slidingWindows(x = chromSizes, width = windowSize, step = slidingSize) %>% unlist %>% .[which(width(.)==windowSize),]
  mcols(windows)$wSeq <- as.character(seqnames(windows))
    mcols(windows)$wStart <- start(windows)
    mcols(windows)$wEnd <- end(windows)
  message("Subtracting Blacklist...")
  windowsBL <- lapply(seq_along(windows), function(x){
      if(x %% 100 == 0){
        message(sprintf("%s of %s", x, length(windows)))
      }
      gr <- GenomicRanges::setdiff(windows[x,], blacklist)
      mcols(gr) <- mcols(windows[x,])
      return(gr)
    })
  names(windowsBL) <- paste0("w",seq_along(windowsBL))
  windowsBL <- unlist(GRangesList(windowsBL), use.names = TRUE)
  mcols(windowsBL)$name <- names(windowsBL)
  message("Adding Nucleotide Information...")
  windowSplit <- split (windowsBL, as.character(seqnames(windowsBL)))
  windowNuc <- lapply(seq_along(windowSplit), function(x){
    message(sprintf("%s of %s", x, length(windowSplit)))
      chrSeq <- Biostrings::getSeq(genome,chromSizes[which(as.character(seqnames(chromSizes))==names(windowSplit)[x])])
      grx <- windowSplit[[x]]
      aFreq <- alphabetFrequency(Biostrings::Views(chrSeq[[1]], ranges(grx)))
      mcols(grx)$GC <- rowSums(aFreq[, c("G","C")]) / rowSums(aFreq)
      mcols(grx)$AT <- rowSums(aFreq[, c("A","T")]) / rowSums(aFreq)
      return(grx)
    }) %>% GRangesList %>% unlist %>% sortSeqlevels %>% sort
  windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
  # get gene density
  gene_density = genes (gene_level_annotation)
  mcols(windowNuc)$gene_density = countOverlaps (windowNuc, gene_density)
  windowNuc
}

# scCNA <- function (windows, fragments, neighbors = 100, force = FALSE, remove = c("chrM","chrX","chrY"),
#   by_col = 'RG',prc_length=10, gene_density = FALSE, bgr_k = 'GC', LFC=1.5, FDR=0.1){
  
#   #Keep only regions in filtered chromosomes
#   windows <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
#   fragments <- GenomeInfoDb::keepStandardChromosomes(fragments, pruning.mode = "coarse")
#   windows <- windows[seqnames(windows) %ni% remove]
#   fragments <- fragments[seqnames(fragments) %ni% remove]

#   #Count Insertions in windows
#   message("Getting Counts...")
#   counts <- countInsertions (windows, fragments, by = by_col)[[1]]
#   message("Summarizing...")
#   windowSummary <- GRangesList()
#   countSummary <- matrix(nrow=length(unique(windows$name)), ncol = ncol(counts))
#   gene_density=NULL
#   for(x in seq_along(unique(mcols(windows)$name)))
#     {
#     if(x %% 100 == 0)
#       {
#       message(sprintf("%s of %s", x, length(unique(mcols(windows)$name))))
#       }
#     idx <- which(mcols(windows)$name == unique(mcols(windows)$name)[x])
#     wx <- windows[idx,]
#     wo <- GRanges(mcols(wx)$wSeq , ranges = IRanges(mcols(wx)$wStart, mcols(wx)$wEnd))[1,]
#     mcols(wo)$name <- mcols(wx)$name[1]
#     mcols(wo)$effectiveLength <- sum(width(wx))
#     mcols(wo)$percentEffectiveLength <- 100*sum(width(wx))/width(wo)
#     mcols(wo)$GC <- sum(mcols(wx)$GC * width(wx))/width(wo)
#     mcols(wo)$AT <- sum(mcols(wx)$AT * width(wx))/width(wo)
#     mcols(wo)$N <- sum(mcols(wx)$N * width(wx))/width(wo)
#     #mcols(wo)$gene_density = sum(as.numeric(mcols(wx)$gene_density) * (width(wx))/width(wo))
#     #gene_density = append (gene_density, sum(as.numeric(mcols(wx)$gene_density) * (width(wx))/width(wo)))
#     countSummary[x,] <- Matrix::colSums(counts[x,,drop=FALSE]) * (wo$percentEffectiveLength / 100)
#     windowSummary[[x]] <- wo
#     }
#   windowSummary = unlist(windowSummary)
#   which(is.na(countSummary), arr.ind=TRUE)
  
#   # Extract residuals
#   #model <- lm(countSummary ~ gene_density)
#   #countSummary <- resid(model)
  
#   #Keep only regions with less than 0.1% N
#   keep <- which (windowSummary$N < 0.1)
#   windowSummary <- windowSummary[keep,]
#   countSummary <- countSummary[keep,]
  
#   #Keep only regions with effective length at least 50% of the window length
#   low_width = windowSummary$percentEffectiveLength > prc_length
#   keep = which (low_width)
#   windowSummary = windowSummary[keep,]
#   countSummary = countSummary[keep,]
#   message (paste('Regions with lower than',prc_length,'% of window length removed:',table(low_width)[1]))
#   #countSummary = countSummary / windowSummary$gene_density # normalising for gene density does not work!!
  
#   # Determine the nearest neighbors by GC content
#   message("Computing Background...")
#   bdgMean = matrix (nrow=nrow(countSummary), ncol= ncol(countSummary))
#   bdgSd = matrix (nrow=nrow(countSummary), ncol= ncol(countSummary))
#   log2FC = matrix (nrow=nrow(countSummary), ncol= ncol(countSummary))
#   z = matrix (nrow=nrow(countSummary), ncol= ncol(countSummary))
#   pval = matrix (nrow=nrow(countSummary), ncol= ncol(countSummary))

#   for(x in seq_len(nrow(countSummary)))
#     {
#     if(x %% 100 == 0)
#       {
#       message(sprintf("%s of %s", x, nrow(countSummary)))
#       }
#   #   #Get Nearest Indices
#     knn = abs(windowSummary$GC[x] - windowSummary$GC)
#     idxNN = head(order(knn), neighbors + 1)
#     idxNN = idxNN[idxNN %ni% x]
#     # idxNN = intersect (idxNN, idxNNgd)
#     #Background
#     if(any(colMeans(countSummary[idxNN, ]) == 0))
#       {
#       if(force)
#         {
#         message("Warning! Background Mean = 0 Try a higher neighbor count or remove cells with 0 in colMins")
#         } else {
#         stop ("Background Mean = 0!")
#         }
#       }
#   bdgMean[x, ] =  apply (countSummary[idxNN,],2,function(x) weighted.mean (x, w= max(knn[idxNN])-knn[idxNN])) # added weighted mean
#     bdgSd[x, ] = matrixStats::colSds(countSummary[idxNN, ])
#     log2FC[x, ] = log2 ((countSummary[x, ] + 1e-5) / (bdgMean[x, ] + 1e-5)) # it was 1e-5 for both
#   z[x, ] = (countSummary[x,] - bdgMean[x, ]) / bdgSd[x, ]
#     pval[x, ] = 2*pnorm(-abs(z[x, ]))
#     }
#   padj = apply (pval, 2, function(x) p.adjust(x, method = "fdr"))
#   CNA = matrix (0, nrow=nrow(countSummary), ncol=ncol(countSummary))
#   CNA [which (log2FC >= LFC & padj <= FDR)] <- 1

#   se <- SummarizedExperiment::SummarizedExperiment(
#     assays = S4Vectors::SimpleList(
#         counts = countSummary,
#         z = z
#         ),
#     rowRanges = windowSummary
#   )
#   colnames(se) <- colnames(counts)
#   return(se)
# }

# Original functions

#"%ni%" <- Negate("%in%")
#
#countInsertions <- function (query, fragments, by = "RG"){
#  #Count By Fragments Insertions
#  inserts <- c(
#    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
#    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
#  )
#  by = 'RG'
#  overlapDF <- DataFrame(findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
#  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
#  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
#  #Calculate Overlap Stats
#  inPeaks <- table(overlapDF$name)
#  total <- table(mcols(inserts)[, by])
#  total <- total[names(inPeaks)]
#  frip <- inPeaks / total
#  #Summarize
#  sparseM <- Matrix::sparseMatrix(
#    i = overlapTDF[, 1], 
#    j = overlapTDF[, 4],
#    x = rep(1, nrow(overlapTDF)), 
#    dims = c(length(query), length(unique(overlapDF$name))))
#  colnames(sparseM) <- unique(overlapDF$name)
#  total <- total[colnames(sparseM)]
#  frip <- frip[colnames(sparseM)]
#  out <- list(counts = sparseM, frip = frip, total = total)
#  return(out)
#}
#
#makeWindows <- function(genome, blacklist, windowSize = 10e6, slidingSize = 2e6){
#  chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
#  chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
#  windows <- slidingWindows(x = chromSizes, width = windowSize, step = slidingSize) %>% unlist %>% .[which(width(.)==windowSize),]
#  mcols(windows)$wSeq <- as.character(seqnames(windows))
#    mcols(windows)$wStart <- start(windows)
#    mcols(windows)$wEnd <- end(windows)
#  message("Subtracting Blacklist...")
#  windowsBL <- lapply(seq_along(windows), function(x){
#      if(x %% 100 == 0){
#        message(sprintf("%s of %s", x, length(windows)))
#      }
#      gr <- GenomicRanges::setdiff(windows[x,], blacklist)
#      mcols(gr) <- mcols(windows[x,])
#      return(gr)
#    })
#  names(windowsBL) <- paste0("w",seq_along(windowsBL))
#  windowsBL <- unlist(GRangesList(windowsBL), use.names = TRUE)
#  mcols(windowsBL)$name <- names(windowsBL)
#  message("Adding Nucleotide Information...")
#  windowSplit <- split(windowsBL, as.character(seqnames(windowsBL)))
#  windowNuc <- lapply(seq_along(windowSplit), function(x){
#    message(sprintf("%s of %s", x, length(windowSplit)))
#      chrSeq <- Biostrings::getSeq(genome,chromSizes[which(as.character(seqnames(chromSizes))==names(windowSplit)[x])])
#      grx <- windowSplit[[x]]
#      aFreq <- alphabetFrequency(Biostrings::Views(chrSeq[[1]], ranges(grx)))
#      mcols(grx)$GC <- rowSums(aFreq[, c("G","C")]) / rowSums(aFreq)
#      mcols(grx)$AT <- rowSums(aFreq[, c("A","T")]) / rowSums(aFreq)
#      return(grx)
#    }) %>% GRangesList %>% unlist %>% sortSeqlevels %>% sort
#  windowNuc$N <- 1 - (windowNuc$GC + windowNuc$AT)
#  windowNuc
#}
#
scCNA = function(windows, fragments, neighbors = 100, LFC = 1.5, FDR = 0.1, force = FALSE, remove = c("chrM","chrX","chrY"),
 by_col = 'RG'){
 
 #Keep only regions in filtered chromosomes
 windows <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
 fragments <- GenomeInfoDb::keepStandardChromosomes(fragments, pruning.mode = "coarse")
 windows <- windows[as.character(seqnames(windows)) %ni% remove]
 fragments <- fragments[as.character(seqnames(fragments)) %ni% remove]

 #Count Insertions in windows
 message("Getting Counts...")
 counts <- countInsertions(windows, fragments, by = by_col)[[1]]
 message("Summarizing...")
 windowSummary <- GRangesList()
 countSummary <- matrix(nrow=length(unique(windows$name)), ncol = ncol(counts))
 for(x in seq_along(unique(mcols(windows)$name))){
   if(x %% 100 == 0){
     message(sprintf("%s of %s", x, length(unique(mcols(windows)$name))))
   }
   idx <- which(mcols(windows)$name == unique(mcols(windows)$name)[x])
   wx <- windows[idx,]
   wo <- GRanges(mcols(wx)$wSeq , ranges = IRanges(mcols(wx)$wStart, mcols(wx)$wEnd))[1,]
   mcols(wo)$name <- mcols(wx)$name[1]
   mcols(wo)$effectiveLength <- sum(width(wx))
   mcols(wo)$percentEffectiveLength <- 100*sum(width(wx))/width(wo)
   mcols(wo)$GC <- sum(mcols(wx)$GC * width(wx))/width(wo)
   mcols(wo)$AT <- sum(mcols(wx)$AT * width(wx))/width(wo)
   mcols(wo)$N <- sum(mcols(wx)$N * width(wx))/width(wo)
   countSummary[x,] <- Matrix::colSums(counts[idx,,drop=FALSE])
   windowSummary[[x]] <- wo
 }
 windowSummary <- unlist(windowSummary)
 
 #Keep only regions with less than 0.1% N
 keep <- which(windowSummary$N < 0.001) 
 windowSummary <- windowSummary[keep,]
 countSummary <- countSummary[keep,]
 
 #Now determine the nearest neighbors by GC content
 message("Computing Background...")
 bdgMean <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
 bdgSd <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
 log2FC <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
 z <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
 pval <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))

 for(x in seq_len(nrow(countSummary))){
   if(x %% 100 == 0){
     message(sprintf("%s of %s", x, nrow(countSummary)))
   }
   #Get Nearest Indices
   idxNN <- head(order(abs(windowSummary$GC[x] - windowSummary$GC)), neighbors + 1)
   idxNN <- idxNN[idxNN %ni% x]
   #Background
   if(any(colMeans(countSummary[idxNN, ])==0)){
     if(force){
       message("Warning! Background Mean = 0 Try a higher neighbor count or remove cells with 0 in colMins")
     }else{
       stop("Background Mean = 0!")
     }
   }
   bdgMean[x, ] <- colMeans(countSummary[idxNN, ])
   bdgSd[x, ] <- matrixStats::colSds(countSummary[idxNN, ])
   log2FC[x, ] <- log2((countSummary[x, ]+1e-5) / (bdgMean[x, ]+1e-5))
   z[x, ] <- (countSummary[x,] - bdgMean[x, ]) / bdgSd[x, ]
   pval[x, ] <- 2*pnorm(-abs(z[x, ]))
 }
 padj <- apply(pval, 2, function(x) p.adjust(x, method = "fdr"))
 CNA <- matrix(0, nrow=nrow(countSummary), ncol=ncol(countSummary))
 CNA[which(log2FC >= LFC & padj <= FDR)] <- 1

 se <- SummarizedExperiment(
   assays = SimpleList(
       CNA = CNA,
       counts = countSummary,
       log2FC = log2FC,
       padj = padj,
       pval = pval,
       z = z,
       bdgMean = bdgMean,
       bdgSd = bdgSd
     ),
   rowRanges = windowSummary
 )
 colnames(se) <- colnames(counts)
 rownames(se) = as.character (windowSummary)
 return(se)
}

################ CNV ANALYSIS END ###################

#--- Original code CNV inference ---#
#Estimating Copy Number Variation in scATAC-seq
#05/02/19
#Cite Satpathy*, Granja*, et al. 
#Massively parallel single-cell chromatin landscapes of human immune 
#cell development and intratumoral T cell exhaustion (2019)
#Created by Jeffrey Granja

#scCNA <- function(windows, fragments, neighbors = 100, LFC = 1.5, FDR = 0.1, force = FALSE, remove = c("chrM","chrX","chrY"),
#  by_col = 'RG2'){
#  
#  #Keep only regions in filtered chromosomes
#  windows   <- GenomeInfoDb::keepStandardChromosomes(windows, pruning.mode = "coarse")
#  fragments <- GenomeInfoDb::keepStandardChromosomes(fragments, pruning.mode = "coarse")
#  windows <- windows[seqnames(windows) %ni% remove]
#  fragments <- fragments[seqnames(fragments) %ni% remove]
#
#  #Count Insertions in windows
#  message("Getting Counts...")
#  counts <- countInsertions(windows, fragments, by = by_col)[[1]]
#  message("Summarizing...")
#  windowSummary <- GenomicRangesList()
#  countSummary <- matrix(nrow=length(unique(windows$name)), ncol = ncol(counts))
#  for(x in seq_along(unique(mcols(windows)$name))){
#    if(x %% 100 == 0){
#      message(sprintf("%s of %s", x, length(unique(mcols(windows)$name))))
#    }
#    idx <- which(mcols(windows)$name == unique(mcols(windows)$name)[x])
#    wx <- windows[idx,]
#    wo <- GRanges(mcols(wx)$wSeq , ranges = IRanges(mcols(wx)$wStart, mcols(wx)$wEnd))[1,]
#    mcols(wo)$name <- mcols(wx)$name[1]
#    mcols(wo)$effectiveLength <- sum(width(wx))
#    mcols(wo)$percentEffectiveLength <- 100*sum(width(wx))/width(wo)
#    mcols(wo)$GC <- sum(mcols(wx)$GC * width(wx))/width(wo)
#    mcols(wo)$AT <- sum(mcols(wx)$AT * width(wx))/width(wo)
#    mcols(wo)$N <- sum(mcols(wx)$N * width(wx))/width(wo)
#    countSummary[x,] <- Matrix::colSums(counts[idx,,drop=FALSE])
#    windowSummary[[x]] <- wo
#  }
#  windowSummary <- unlist(windowSummary)
#  
#  #Keep only regions with less than 0.1% N
#  keep <- which(windowSummary$N < 0.001) 
#  windowSummary <- windowSummary[keep,]
#  countSummary <- countSummary[keep,]
#  
#  #Now determine the nearest neighbors by GC content
#  message("Computing Background...")
#  bdgMean <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
#  bdgSd <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
#  log2FC <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
#  z <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
#  pval <- matrix(nrow=nrow(countSummary), ncol=ncol(countSummary))
#
#  for(x in seq_len(nrow(countSummary))){
#    if(x %% 100 == 0){
#      message(sprintf("%s of %s", x, nrow(countSummary)))
#    }
#    #Get Nearest Indices
#    idxNN <- head(order(abs(windowSummary$GC[x] - windowSummary$GC)), neighbors + 1)
#    idxNN <- idxNN[idxNN %ni% x]
#    #Background
#    if(any(colMeans(countSummary[idxNN, ])==0)){
#      if(force){
#        message("Warning! Background Mean = 0 Try a higher neighbor count or remove cells with 0 in colMins")
#      }else{
#        stop("Background Mean = 0!")
#      }
#    }
#    bdgMean[x, ] <- colMeans(countSummary[idxNN, ])
#    bdgSd[x, ] <- matrixStats::colSds(countSummary[idxNN, ])
#    log2FC[x, ] <- log2((countSummary[x, ]+1e-5) / (bdgMean[x, ]+1e-5))
#    z[x, ] <- (countSummary[x,] - bdgMean[x, ]) / bdgSd[x, ]
#    pval[x, ] <- 2*pnorm(-abs(z[x, ]))
#  }
#  padj <- apply(pval, 2, function(x) p.adjust(x, method = "fdr"))
#  CNA <- matrix(0, nrow=nrow(countSummary), ncol=ncol(countSummary))
#  CNA[which(log2FC >= LFC & padj <= FDR)] <- 1
#
#  se <- SummarizedExperiment(
#    assays = SimpleList(
#        CNA = CNA,
#        counts = countSummary,
#        log2FC = log2FC,
#        padj = padj,
#        pval = pval,
#        z = z,
#        bdgMean = bdgMean,
#        bdgSd = bdgSd
#      ),
#    rowRanges = windowSummary
#  )
#  colnames(se) <- colnames(counts)
#
#  return(se)
#}




#' Generate bin by cell matrix for scATAC-seq data
# From: https://github.com/seasoncloud/Basic_CNV_SNV/blob/main/scripts/Gen_bin_cell_atac.R
#'
#' @param bin_bed A matrix of the BED format for fixed bins across the genome (not the bed file for each peak region). 
#' The first three columns are 'chr', "start site", "end site". Each row is a bin region.
#' There is a helper function "Generate_bed.R" to help generate this. 
#' @param barcodes A matrix/ data.frame with barcodes for each cell in the first column.
#' @param path_to_fragments The path to the "fragments.tsv.gz" file with the same format as that from the Cell Ranger software.
#' @param out_path The path to save the output matrix.
#' 
#' @import rtracklayer
#' @return A bin by cell matrix as the raw_counts input for Alleloscope.
#'
#' @export
## generate bed files for the bins
Gen_bin_cell_atac=function(bin_bed=NULL, barcodes=NULL, path_to_fragments="./fragments.tsv.gz", out_path="./" ){
  # generate bin by cell matrix from fragment file
  barcodes=barcodes[,1]
  
  cat("Read fragment file...\n")
  fragments <- import.bed(path_to_fragments, extraCols = c( "type"="character", "type"="integer"))
  colnames(mcols(fragments)) <- c("barcode", "dup_counts")
  
  cat(paste0("Total ", length(fragments)," fragments.\n"))
  #length(unique(fragments$barcode))
  fragments_incell <- fragments[fragments$barcode %in% barcodes]
  
  cat(paste0("Total ", length(fragments_incell)," fragments in cell.\n")) # 236637826
  
  ## check bin bed chr column
  if(grepl("chr",bin_bed[1,1])){
    bin_bed[,1]=gsub("chr","",bin_bed[,1])
  }else{
    bin_bed=bin_bed
  }
  
  chr200k=bin_bed
  chr200k=chr200k[which(as.numeric(chr200k[,1]) %in% 1:22),]
  chr200k=chr200k[order(as.numeric(chr200k[,1]), as.numeric(chr200k[,2])),]
  bins=paste0('chr',chr200k[,1],'-',chr200k[,2],"-", chr200k[,3])
  query=GRanges(paste0('chr',chr200k[,1]), IRanges(chr200k[,2]+1,chr200k[,3]))
  
  ov=findOverlaps(query, fragments_incell )
  ov=as.matrix(ov)
  tmp=fragments_incell$barcode[ov[,2]]
  ov=cbind(ov,match(tmp, barcodes))
  
  cat("Generate bin-by-cell matrix...\n")
  mm=table(ov[,1],ov[,3])
  colnames(mm)=barcodes[as.numeric(colnames(mm))]
  rownames(mm)=bins[as.numeric(rownames(mm))]
  
  message("The bin-by-cell matrix has beed successfully generated!")
  saveRDS(mm,paste0(out_path,"/bin_cell_atac_fragments.rds"))
  return(mm)
}


#' Help generate customized bin_bed file of the BED forma as the input for the Gen_bin_cell_atac.R function
# From: https://github.com/seasoncloud/Basic_CNV_SNV/blob/main/scripts/Generate_bed.R 
#' @param chr_size A 22x2 matrix. The first column is 'chr' and the second column is the size (bps), Each row is a chromosome (from 1 to 22).
#' The matrix for hg19 and GRCh38 can be direcly downloaded from ./data-raw/ directory.
#' @param bin_res Numeric. The fixed bin size for generating the bed file across the chromosomes.
#' @param out_path The path to save the output bed file.
#' 
#' @return A bin_bed file across the chromosomes. 
#'
#' @export
## generate bed files for the bins
Generate_bed=function(chr_size=NULL, bin_res=200000, out_path='./'){
chr_bin=NULL
for(ii in 1:22){
  sub=chr_size[which(chr_size[,1]==paste0('chr',ii)),, drop=F]
  sub=as.numeric(c(0,sub[,2]) )
  minn=min(sub)
  maxx=max(sub)
  rr=maxx-minn
  nn=ceiling(rr/bin_res)
  mm=cbind(rep(ii, nn), (0:(nn-1))*bin_res, c((1:(nn-1))*bin_res, maxx))
  chr_bin=rbind(chr_bin, mm)
  print(ii)
}

return(chr_bin)
write.table(chr_bin,paste0(out_path, "./chrbin", bin_res, ".bed"), sep='\t', col.names=F, row.names=F, quote=F)

}



### ArchR function to identify differential markers fixed ###
# From https://github.com/GreenleafLab/ArchR/issues/562
.MarkersSC_fixed <- function(
  ArchRProj = NULL,
  groupBy = "Clusters",
  useGroups = NULL,
  bgdGroups = NULL,
  normBy = NULL,
  maxCells = 500,
  scaleTo = 10^4,
  bufferRatio = 0.8,
  bias = NULL,
  k = 100,
  threads = 1,
  binarize = FALSE,
  useSeqnames = NULL,
  testMethod = "wilcoxon",
  useMatrix = "GeneScoreMatrix",
  markerParams = list(),
  verbose = TRUE,
  logFile = NULL
){
  
  tstart <- Sys.time()
  
  #####################################################
  # Feature Info
  #####################################################
  ArrowFiles <- getArrowFiles(ArchRProj)
  featureDF <- .getFeatureDF(head(ArrowFiles, 2), useMatrix)
  matrixClass <- as.character(h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix, "/Info/Class")))
  
  .logThis(range(as.vector(table(paste0(featureDF$seqnames)))), "FeaturesPerSeqnames", logFile = logFile)
  
  isDeviations <- FALSE
  if(all(unique(paste0(featureDF$seqnames)) %in% c("z", "dev"))){
    isDeviations <- TRUE
  }
  
  .logThis(featureDF, "FeatureDF", logFile=logFile)
  .logMessage(paste0("MatrixClass = ", matrixClass), logFile=logFile)
  
  seqnames <- unique(as.vector(featureDF$seqnames))
  useSeqnames <- useSeqnames[useSeqnames %in% seqnames]
  if(length(useSeqnames)==0){
    useSeqnames <- NULL
  }
  
  if(!is.null(useSeqnames)){
    if(matrixClass == "Sparse.Assays.Matrix"){
      if(length(useSeqnames) == 1){
        featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
      }else{
        .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires 1 seqname!\n",
                    "Continuing with first seqname '", seqnames[1], "'!\n",
                    "If confused, try getSeqnames(ArchRProj, '", useMatrix,"'') to list out available seqnames for input!", verbose = verbose, logFile = logFile)
        useSeqnames <- seqnames[1]
        featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
      }
    }else{
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }
  }else{
    if(matrixClass == "Sparse.Assays.Matrix"){
      .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires 1 seqname!\n",
                  "Continuing with first seqname '", seqnames[1], "'!\n",
                  "If confused, try getSeqnames(ArchRProj, '", useMatrix,"'') to list out available seqnames for input!", verbose = verbose, logFile = logFile)
      useSeqnames <- seqnames[1]
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }
  }
  if(!(nrow(featureDF) > 1)){
    .logStop("Less than 1 feature is remaining in featureDF please check input!", logFile = logFile)
  }
  
  #####################################################
  # Match Bias Groups
  #####################################################
  .logDiffTime("Matching Known Biases", t1 = tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
  groups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  colDat <- getCellColData(ArchRProj)
  
  matchObj <- .matchBiasCellGroups(
    input = colDat, 
    groups = groups,
    useGroups = useGroups,
    bgdGroups = bgdGroups,
    bias = bias,
    k = k,
    n = maxCells,
    bufferRatio = bufferRatio,
    logFile = logFile
  )
  
  #####################################################
  # Pairwise Test Per Seqnames
  #####################################################
  mColSums <- tryCatch({
    suppressMessages(tmpColSum <- .getColSums(ArrowFiles, seqnames = featureDF$seqnames@values, useMatrix = useMatrix, threads = threads))
    tmpColSum[ArchRProj$cellNames]
    }, error = function(x){
    rep(1, nCells(ArchRProj))
  })
  if(all(mColSums==1) & is.null(normBy)){
    normBy <- "none"
  }
  
  if(is.null(normBy)){
    if(tolower(testMethod) == "binomial"){
      normFactors <- NULL
    }else{
      if(tolower(useMatrix) %in% c("tilematrix", "peakmatrix")){
        normBy <- "ReadsInTSS"
        normFactors <- getCellColData(ArchRProj, normBy, drop=FALSE)
        normFactors[,1] <- median(normFactors[,1]) / normFactors[,1]
      }else{
        normFactors <- scaleTo / mColSums
        normFactors <- DataFrame(normFactors)
      }
    }
  }else{
    if(tolower(normBy) == "none"){
      normFactors <- NULL
    }else{
      normFactors <- scaleTo / mColSums
      normFactors <- DataFrame(normFactors)
    }
  }
  
  if(!is.null(normFactors)){
    normFactors[,1] <- normFactors[,1] * (scaleTo / median(normFactors[names(mColSums), 1] * mColSums))
  }
  
  diffList <- .safelapply(seq_along(matchObj[[1]]), function(x){
    .logDiffTime(sprintf("Computing Pairwise Tests (%s of %s)", x, length(matchObj[[1]])), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    .testMarkerSC(
      ArrowFiles = ArrowFiles,
      matchObj = matchObj, 
      group = names(matchObj[[1]])[x], 
      testMethod = testMethod, 
      threads = 1, 
      useMatrix = useMatrix,
      featureDF = featureDF,
      normFactors = normFactors,
      binarize = binarize,
      logFile = logFile
    )
  }, threads = threads)
  
  .logDiffTime("Completed Pairwise Tests", tstart, addHeader = TRUE, verbose = verbose, logFile = logFile)
  
  #####################################################
  # Summarize Output
  #####################################################
  if(tolower(testMethod) == "wilcoxon"){
    pse <- SummarizedExperiment::SummarizedExperiment(
      assays = 
        SimpleList(
          Log2FC = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$log2FC)) %>% Reduce("cbind",.),
          Mean = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1)) %>% Reduce("cbind",.),
          FDR = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$fdr)) %>% Reduce("cbind",.),
          Pval = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$pval)) %>% Reduce("cbind",.),
          MeanDiff = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1 - diffList[[x]]$mean2)) %>% Reduce("cbind",.),
          AUC = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$auc)) %>% Reduce("cbind",.),
          MeanBGD = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean2)) %>% Reduce("cbind",.)
        ),
      rowData = featureDF
    )
  }else if(tolower(testMethod) == "ttest"){
    pse <- SummarizedExperiment::SummarizedExperiment(
      assays = 
        SimpleList(
          Log2FC = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$log2FC)) %>% Reduce("cbind",.),
          Mean = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1)) %>% Reduce("cbind",.),
          Variance = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$var1)) %>% Reduce("cbind",.),
          FDR = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$fdr)) %>% Reduce("cbind",.),
          Pval = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$pval)) %>% Reduce("cbind",.),
          MeanDiff = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1 - diffList[[x]]$mean2)) %>% Reduce("cbind",.),
          MeanBGD = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean2)) %>% Reduce("cbind",.),
          VarianceBGD = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$var2)) %>% Reduce("cbind",.)
        ),
      rowData = featureDF
    )
  }else if(tolower(testMethod) == "binomial"){
    pse <- SummarizedExperiment::SummarizedExperiment(
      assays = 
        SimpleList(
          Log2FC = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$log2FC)) %>% Reduce("cbind",.),
          Mean = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1)) %>% Reduce("cbind",.),
          FDR = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$fdr)) %>% Reduce("cbind",.),
          Pval = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$pval)) %>% Reduce("cbind",.),
          MeanDiff = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean1 - diffList[[x]]$mean2)) %>% Reduce("cbind",.),
          MeanBGD = lapply(seq_along(diffList), function(x) data.frame(x = diffList[[x]]$mean2)) %>% Reduce("cbind",.)
        ),
      rowData = featureDF
    )
  }else{
    stop("Error Unrecognized Method!")
  }
  colnames(pse) <- names(matchObj[[1]])
  
  metadata(pse)$MatchInfo <- matchObj
  
  if(isDeviations){
    assays(pse)[["Log2FC"]] <- NULL #This measure does not make sense with deviations matrices better to just remove
  }
  
  
  return(pse)
  
}

environment(.MarkersSC_fixed) <- asNamespace('ArchR')
assignInNamespace(".MarkersSC", .MarkersSC_fixed, ns='ArchR')





### FIX STUPID BUG WHEN DOUBLET FUNCTION DOESNT COMPLETE AND MESSES UP COLDATA PER SAMPLE ####

getMatrixFromProject <- function(
  ArchRProj = NULL,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = useSeqnames, name = "useSeqnames", valid = c("character","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getMatrixFromProject Input-Parameters", logFile = logFile)

  ArrowFiles <- getArrowFiles(ArchRProj)

  cellNames <- ArchRProj$cellNames

  avMat <- getAvailableMatrices(ArchRProj)
  if(useMatrix %ni% avMat){
    stop("useMatrix is not in Available Matrices see getAvailableMatrices")
  }

  seL <- .safelapply(seq_along(ArrowFiles), function(x){

    .logDiffTime(paste0("Reading ", useMatrix," : ", names(ArrowFiles)[x], "(",x," of ",length(ArrowFiles),")"), 
      t1 = tstart, verbose = FALSE, logFile = logFile)

    allCells <- .availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]

    if(length(allCells) != 0){

      o <- getMatrixFromArrow(
        ArrowFile = ArrowFiles[x],
        useMatrix = useMatrix,
        useSeqnames = useSeqnames,
        cellNames = allCells, 
        ArchRProj = ArchRProj,
        verbose = FALSE,
        binarize = binarize,
        logFile = logFile
      )

      .logDiffTime(paste0("Completed ", useMatrix," : ", names(ArrowFiles)[x], "(",x," of ",length(ArrowFiles),")"), 
        t1 = tstart, verbose = FALSE, logFile = logFile)

      o

    }else{

      NULL
      
    }

  }, threads = threads) 

  #ColData
  .logDiffTime("Organizing colData", t1 = tstart, verbose = verbose, logFile = logFile)
  cD <- lapply(seq_along(seL), function(x){
    colData(seL[[x]]) = colData(seL[[x]])[!colnames(colData(seL[[x]])) %in% c("DoubletEnrichment", "DoubletScore")]  ### THIS IS THE LINE I ADDED 
    colData(seL[[x]])
  }) %>% Reduce("rbind", .)
  
  #RowData
  .logDiffTime("Organizing rowData", t1 = tstart, verbose = verbose, logFile = logFile)
  rD1 <- rowData(seL[[1]])
  rD <- lapply(seq_along(seL), function(x){
    identical(rowData(seL[[x]]), rD1)
  }) %>% unlist %>% all
  if(!rD){
    stop("Error with rowData being equal for every sample!")
  }

  #RowRanges
  .logDiffTime("Organizing rowRanges", t1 = tstart, verbose = verbose, logFile = logFile)
  rR1 <- rowRanges(seL[[1]])
  rR <- lapply(seq_along(seL), function(x){
    identical(rowRanges(seL[[x]]), rR1)
  }) %>% unlist %>% all
  if(!rR){
    stop("Error with rowRanges being equal for every sample!")
  }

  #Assays
  nAssays <- names(assays(seL[[1]]))
  asy <- lapply(seq_along(nAssays), function(i){
    .logDiffTime(sprintf("Organizing Assays (%s of %s)", i, length(nAssays)), t1 = tstart, verbose = verbose, logFile = logFile)
    m <- lapply(seq_along(seL), function(j){
      assays(seL[[j]])[[nAssays[i]]]
    }) %>% Reduce("cbind", .)
    m
  }) %>% SimpleList()
  names(asy) <- nAssays
  
  .logDiffTime("Constructing SummarizedExperiment", t1 = tstart, verbose = verbose, logFile = logFile)
  if(!is.null(rR1)){
    se <- SummarizedExperiment(assays = asy, colData = cD, rowRanges = rR1)
    se <- sort(se)
  }else{
    se <- SummarizedExperiment(assays = asy, colData = cD, rowData = rD1)
  }
  rm(seL)
  gc()

  .logDiffTime("Finished Matrix Creation", t1 = tstart, verbose = verbose, logFile = logFile)

  se
  
}



# Calculate hypergeometric test of peaks with motifs
hyperMotif = function (selected_peaks, motifmatch)
  {
  bg_peakSet = rowRanges (motifmatch)
  bg_peak_motif_mat = assays (motifmatch)[[1]]
  selected_peaks_idx = queryHits (findOverlaps (bg_peakSet, selected_peaks))
  all_peaks_w_motif = colSums (bg_peak_motif_mat)
  selected_peaks_w_motif = colSums (bg_peak_motif_mat[selected_peaks_idx,])
  
  # Run hypergeometric test per celltype using the cluster-specific peaks for that celltypes as background
  hyper_res = list()
    
  # Run hypergeometric test per TF
  message (paste('Run hyper tests for all TFs'))
  hyper_tf = NULL
  pb = progress::progress_bar$new(total =  ncol(bg_peak_motif_mat))
  for (i in 1:ncol(bg_peak_motif_mat))
    {
    pb$tick()
    q = selected_peaks_w_motif[i]
    m = all_peaks_w_motif[i]
    n = nrow(bg_peak_motif_mat) - m
    k = length (selected_peaks)
    hyper_tf[i] = phyper (q, m, n, k, lower.tail = FALSE, log.p = FALSE)
    }
  hyper_res = data.frame (
      pval = hyper_tf, 
      padj = p.adjust (hyper_tf, method = 'fdr'),
      row.names = colnames (bg_peak_motif_mat))
  hyper_res = hyper_res[order (hyper_res$padj),]
  return (hyper_res)
  }


cCNV_score = function (x)
  {
   #x_pos = mean (unlist(lapply(split(x, cumsum(x < 0)), length)))
   x_pos = split(x, cumsum(x < 0))
   #x_pos = x_pos[sapply (x_pos, function(x) length(x) != 1)]
   x_pos = unlist(lapply (x_pos, length))
   x_pos = x_pos[order (-x_pos)]
   x_pos = mean(head (x_pos, round (length(x_pos) / 10)))
   #x_neg = mean (unlist(lapply(split(x, cumsum(x > 0)), length)))
   x_neg = split(x, cumsum(x > 0))
   #x_neg = x_neg[sapply (x_neg, function(x) length(x) != 1)]
   x_neg = unlist(lapply (x_neg, length))
   x_neg = x_neg[order (-x_neg)]
   x_neg = mean(head (x_neg, round(length(x_neg) / 10)))
   #x_neg = summary (unlist(lapply (x_neg, length)))[5]
  return (x_pos + abs(x_neg))
  }


subsetArchRProject_light = function(
  ArchRProject = NULL, 
  cells = NULL,
  projdir_new = NULL)
{
source_dir = getOutputDirectory(ArchRProject)
projdir_temp = file.path(source_dir, 'tempdir')
dir.create (file.path(projdir_new,'Plots'), recursive=T)
  copy_project_files <- function(source_dir, projdir_temp) {
    # Ensure both source_dir and projdir_temp are provided
    if (missing(source_dir) || missing(projdir_temp)) {
      stop("Both 'source_dir' and 'projdir_temp' must be specified.")
    }
    
    # Create the projdir_temp directory if it doesn't exist
    if (!dir.exists(projdir_temp)) {
      dir.create(projdir_temp, recursive = TRUE)
      message(paste("Created target directory:", projdir_temp))
    }
    
    # Construct paths for ArrowFiles and Save-ArchR-Project.rds
    arrow_files_path <- file.path(source_dir, "ArrowFiles")
    rds_file_path <- file.path(source_dir, "Save-ArchR-Project.rds")
    
    # Check if ArrowFiles directory exists
    if (!dir.exists(arrow_files_path)) {
      stop("ArrowFiles directory does not exist in the source directory.")
    }
    
    # Check if Save-ArchR-Project.rds file exists
    if (!file.exists(rds_file_path)) {
      stop("Save-ArchR-Project.rds file does not exist in the source directory.")
    }
    
    # Copy ArrowFiles directory
    arrow_copy_command <- sprintf("cp -r '%s' '%s/'", arrow_files_path, projdir_temp)
    system(arrow_copy_command)
    message(paste("Copied ArrowFiles to:", projdir_temp))
    
    # Copy Save-ArchR-Project.rds file
    rds_copy_command <- sprintf("cp '%s' '%s/'", rds_file_path, projdir_temp)
    system(rds_copy_command)
    message(paste("Copied Save-ArchR-Project.rds to:", projdir_temp))
  }
message ('copy arrow files and archr object to temp directory...')  
copy_project_files (source_dir, projdir_temp)
ArchRProject = loadArchRProject (projdir_temp, force=T)
saveArchRProject (ArchRProject[rownames(ArchRProject@cellColData) %in% cells], outputDirectory = projdir_new, dropCells=T)

unlink(file.path(projdir_temp, "*"), recursive = TRUE, force = TRUE)
  message(paste("All contents in", projdir_temp, "have been deleted."))
}





# Find DAM ####
### ChromVAR based analysis ####

DAM = function (
  ArchRProj = NULL,
  metaGroupName = NULL,
  FDR_threshold = 1e-2,
  meandiff_threshold = 0,
  top_genes=5,
  filter_by_scRNA=TRUE, # use 'srt' object
  seurat_obj = NULL,
  min_exp=.1,
  force = FALSE) {

  if (!file.exists (paste0('DAM_',metaGroupName,'.rds')) | force)
    {
    DAM_list = getMarkerFeatures (
      ArchRProj = ArchRProj, 
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
  
  DAM_list = lapply(DAM_list, function(res) {
  res = res[res$FDR < FDR_threshold, ]
  res = res[order(res$FDR), ]
  res = res[res$MeanDiff > meandiff_threshold, , drop = FALSE]
  if (nrow(res) == 0) return(NULL)
  res
  })

if (filter_by_scRNA)
    {
    message ('filtering using scRNA')
    ps = log2(as.data.frame (AverageExpression (seurat_obj, 
    features = sapply (unique(unlist(lapply(DAM_list, function(x) x$gene))), function(x) unlist(strsplit (x, '_'))[1]), 
    group.by = metaGroupName)[[1]]) +1)
    colnames(ps) = gsub ('-','_',colnames(ps))

    if (all(unique (ArchRProj@cellColData[,metaGroupName]) %in% names (table (seurat_obj@meta.data[,metaGroupName]))))
    {
    message ('filtering celltype-wise')
    active_TFs = unique(unlist(lapply (names (DAM_list), function(x) DAM_list[[x]]$gene[DAM_list[[x]]$gene %in% rownames(ps)[ps[,x] > min_exp]])))
    } else {
    message ('filtering across all cells')
    ps = ps[apply(ps, 1, function(x) any (x > min_exp)),]
    active_TFs = rownames(ps)  
    }

  #active_genes = corGSM_MM$MotifMatrix_name[corGSM_MM$cor > 0.1]
  DAM_list2 = lapply (DAM_list, function(x) x[x$gene %in% active_TFs,])    
  } else {
  DAM_list2 = DAM_list  
  }
  DAM_top_list = lapply (names(DAM_list2), function(x)
    {
     res = DAM_list2[[x]]
     res$comparison = x
    head (res, top_genes)
    })
    DAM_top_list = DAM_top_list[unlist(sapply (  DAM_top_list, function(x) !is.null (dim(x))))]
  DAM_df = Reduce (rbind ,DAM_top_list)
return (DAM_df)
}


# Plot ArchR gene score feature plots with minimal styling (gene name only).
# ArchR's plotEmbedding sets the title to a verbose string like
# "UMAP_H of Harmony colored by\nGeneScoreMatrix : CD3D"; this function
# extracts just the gene name and strips all other text and the legend.
#
# Usage:
#   archp <- addImputeWeights(archp)
#   plots <- plotFeaturePlots(archp, genes = c("CD3D","MS4A1"), embedding = "UMAP_H")
#   pdf("out.pdf", width=25, height=25); patchwork::wrap_plots(plots); dev.off()
plotFeaturePlots = function(
  ArchRProj,
  genes,
  embedding     = "UMAP_H",
  pal           = viridis::plasma(100),
  imputeWeights = getImputeWeights(ArchRProj)
) {
  genes <- genes[genes %in% getFeatures(ArchRProj, "GeneScoreMatrix")]

  clean_theme <- ggplot2::theme(
    axis.title      = ggplot2::element_blank(),
    axis.text       = ggplot2::element_blank(),
    axis.ticks      = ggplot2::element_blank(),
    axis.line       = ggplot2::element_blank(),
    plot.subtitle   = ggplot2::element_blank(),
    plot.caption    = ggplot2::element_blank(),
    legend.position = "none",
    panel.grid      = ggplot2::element_blank(),
    panel.border    = ggplot2::element_blank()
  )

  pdf(NULL)
  plots <- plotEmbedding(
    ArchRProj     = ArchRProj,
    colorBy       = "GeneScoreMatrix",
    name          = genes,
    embedding     = embedding,
    pal           = pal,
    imputeWeights = imputeWeights
  )
  dev.off()

  lapply(plots, function(p) {
    # ArchR title: "UMAP_H of Harmony colored by\nGeneScoreMatrix : GENE"
    p$labels$title    <- sub(".*: ", "", p$labels$title)
    p$labels$subtitle <- NULL
    p$labels$x        <- NULL
    p$labels$y        <- NULL
    p$labels$fill     <- NULL
    p$labels$colour   <- NULL
    p$labels$color    <- NULL
    p$labels$z        <- NULL
    p + clean_theme + ggplot2::guides(colour = "none", fill = "none", color = "none")
  })
}



###############################################################################
# plotTracksSideBySide
#
# Browser tracks for several genomic regions drawn as panels in ONE row, so that loci can be
# compared across cell types at a glance. Wraps ArchR::plotBrowserTrack and then rebuilds the
# layout, because stacking ArchR panels as-is spends most of the page on labels that repeat
# once per region (y-axis title, facet strip, coordinate axis).
#
# What it changes relative to a plain plotBrowserTrack:
#   - cell-type names drawn INSIDE each coverage facet (top left), not as a strip per panel,
#     and by default only on the left-most panel of each row
#   - gene bodies coloured by strand with the symbol centred above each gene
#   - coverage profiles optionally outlined
#   - panel frames and coordinate axes optionally removed; the region is named in the header
#   - whitespace between the stacked tracks collapsed
#   - loop-track colour key kept on the right-most panel only
#
# All of the above is done by temporarily rewriting ArchR internals (.bulkTracks, .geneTracks,
# .loopTracks, theme_ArchR). They are ALWAYS restored on exit, including on error, so nothing
# leaks into the calling session.
#
# ARGUMENTS
#  ArchRProj    ArchR project.
#  regions      GRanges to plot, or a character vector of hub ids (then `hubs` is required).
#               An `mcols(regions)$label` column, if present, is used for the panel header.
#  hubs         Optional list with `hubs_id`, `hubsMerged`, `hubsCollapsed`, `peakLinks`
#               (i.e. a global_hubs_obj.rds). Used to resolve hub ids to regions and, unless
#               overridden, to supply the peak feature track and the pairwise link track.
#  hub_bed      Optional data.frame/path of hub_regions.bed (chr,start,end,id,strand,genes),
#               used to resolve hub ids when `hubs` has no hubsCollapsed names.
#  groupBy      cellColData column defining the tracks (default 'celltype_lv1').
#  groups       character vector giving which groups to show AND their top-to-bottom order.
#  pal          named colour vector for the groups.
#  features     GRanges/GRangesList for the feature track; defaults to the hubs' own peaks.
#  loops        GRanges/GRangesList for the loop track; defaults to the hubs' pairwise links.
#  file         output PDF path. If NULL the assembled gtable(s) are returned invisibly.
#  perRow       wrap panels over several rows (NULL = a single row).
#  header       function(i, regions) returning the panel header, or NULL for id + coordinates.
#
# VALUE  the assembled gtable (one row) or a list of gtables (wrapped), invisibly.
#
# EXAMPLE
#   hubs = readRDS('hubs_obj_cor_0.3.../global_hubs_obj.rds')
#   bed  = read.table('hubs_obj_cor_0.3.../hub_regions.bed', sep='\t')
#   plotTracksSideBySide(archp, c('HUB9','HUB32','HUB85'), hubs = hubs, hub_bed = bed,
#                        groups = c('Malignant','Myeloid','T_cells','NK'),
#                        pal = palette_celltype_lv1, file = 'tracks.pdf')
###############################################################################
plotTracksSideBySide = function(
  ArchRProj,
  regions,
  hubs            = NULL,
  hub_bed         = NULL,
  groupBy         = 'celltype_lv1',
  groups          = NULL,
  pal             = NULL,
  features        = NULL,
  loops           = NULL,
  file            = NULL,
  perRow          = NULL,
  padding         = 0.05,
  plotSummary     = c('bulkTrack', 'featureTrack', 'loopTrack', 'geneTrack'),
  sizes           = c(10, 0.45, 1, 1.9),
  tileSize        = 200,
  minCells        = 25,
  panelWidth      = 2.1,
  panelHeight     = 8.1,
  labelWidth      = 2.6,
  genePlusColor   = 'brown',
  geneMinusColor  = 'black',
  geneLabelSize   = 3,
  geneLabelNudge  = 0.18,
  groupLabelSize  = 3.4,
  labelFirstOnly  = TRUE,
  loopPal         = NULL,
  outlineProfiles = TRUE,
  outlineWidth    = 0.05,
  panelBorders    = FALSE,
  showXaxis       = FALSE,
  baseSize        = 11,
  facetSize       = 11,
  headerSize      = 9.2,
  header          = NULL,
  caption         = NULL,
  verbose         = TRUE
){
  for (pkg in c('ArchR','GenomicRanges','ggplot2','ggrepel','grid','gtable'))
    if (!requireNamespace(pkg, quietly = TRUE)) stop('plotTracksSideBySide needs package: ', pkg)
  requireNamespace('ggrepel', quietly = TRUE)   # drawn lazily by the gene track

  ## ---- 1. resolve regions ---------------------------------------------------------------
  if (is.character(regions)) {
    if (is.null(hubs) && is.null(hub_bed)) stop('hub ids given but neither `hubs` nor `hub_bed` supplied')
    if (!is.null(hub_bed)) {
      bed = if (is.character(hub_bed)) utils::read.table(hub_bed, sep = '\t', stringsAsFactors = FALSE) else hub_bed
      colnames(bed)[1:6] = c('chr','start','end','id','strand','genes')
      miss = setdiff(regions, bed$id)
      if (length(miss)) warning('not found in hub_bed, skipped: ', paste(miss, collapse = ', '))
      regions = regions[regions %in% bed$id]
      b  = bed[match(regions, bed$id), ]
      gr = GenomicRanges::GRanges(b$chr, IRanges::IRanges(b$start, b$end))
      GenomicRanges::mcols(gr)$label = b$id
      GenomicRanges::mcols(gr)$genes = b$genes
    } else {
      idx = match(regions, hubs$hubs_id)
      if (anyNA(idx)) { warning('hub ids not in `hubs`, skipped: ',
                                paste(regions[is.na(idx)], collapse = ', ')); idx = idx[!is.na(idx)] }
      gr = hubs$hubsCollapsed[idx]
      GenomicRanges::mcols(gr) = NULL
      GenomicRanges::mcols(gr)$label = hubs$hubs_id[idx]
    }
    hub_ids = GenomicRanges::mcols(gr)$label
  } else {
    gr = regions
    if (is.null(GenomicRanges::mcols(gr)$label))
      GenomicRanges::mcols(gr)$label = sprintf('%s:%s-%s', as.character(GenomicRanges::seqnames(gr)),
                                               GenomicRanges::start(gr), GenomicRanges::end(gr))
    hub_ids = NULL
  }
  if (!length(gr)) stop('no regions left to plot')
  w  = GenomicRanges::width(gr)
  gr = GenomicRanges::resize(gr, w + 2 * round(padding * w), fix = 'center')
  n  = length(gr)

  ## ---- 2. feature and loop tracks from the hub object, unless supplied -------------------
  if (is.null(features) && !is.null(hubs) && !is.null(hub_ids)) {
    k = match(hub_ids, hubs$hubs_id)
    pk = do.call(c, lapply(k, function(j) {
      d = hubs$hubsMerged[[j]]
      GenomicRanges::GRanges(as.character(d$seqnames), IRanges::IRanges(d$start, d$end))
    }))
    features = GenomicRanges::GRangesList(`cHub peaks` = pk)
  }
  if (is.null(loops) && !is.null(hubs) && !is.null(hubs$peakLinks)) {
    lk = hubs$peakLinks[[1]]
    loops = GenomicRanges::GRangesList(`cHub links` = IRanges::subsetByOverlaps(lk, gr, type = 'within'))
  }
  if (!'looptrack'    %in% tolower(plotSummary)) loops = NULL
  if (!'featuretrack' %in% tolower(plotSummary)) features = NULL
  if (length(sizes) != length(plotSummary)) stop('`sizes` must have one value per `plotSummary` entry')

  ## ---- 3. groups and palette --------------------------------------------------------------
  cd = ArchR::getCellColData(ArchRProj)
  if (!groupBy %in% colnames(cd)) stop('groupBy not in cellColData: ', groupBy)
  if (is.null(groups)) groups = sort(unique(as.character(cd[[groupBy]])))
  groups = groups[groups %in% unique(as.character(cd[[groupBy]]))]
  if (!length(groups)) stop('none of `groups` present in ', groupBy)
  if (!is.null(pal)) {
    if (!all(groups %in% names(pal))) stop('`pal` is missing: ', paste(setdiff(groups, names(pal)), collapse = ', '))
    pal = pal[groups]
  }
  ArchRProj = ArchRProj[as.character(cd[[groupBy]]) %in% groups, ]

  ## ---- 4. temporarily rewrite the ArchR track builders ------------------------------------
  ## deparse() wraps at width.cutoff, so sources are rejoined with newlines and every pattern
  ## tolerates a line break inside the expression.
  ns = asNamespace('ArchR')
  orig = list(bulk  = ArchR:::.bulkTracks, gene = ArchR:::.geneTracks,
              loop  = ArchR:::.loopTracks, theme = ArchR::theme_ArchR)
  on.exit({
    assignInNamespace('.bulkTracks', orig$bulk,  ns = 'ArchR')
    assignInNamespace('.geneTracks', orig$gene,  ns = 'ArchR')
    assignInNamespace('.loopTracks', orig$loop,  ns = 'ArchR')
    assignInNamespace('theme_ArchR', orig$theme, ns = 'ArchR')
  }, add = TRUE)

  rewrite = function(fun, subs) {
    txt = paste(deparse(fun), collapse = '\n')
    for (u in subs) {
      if (!grepl(u[1], txt)) stop('plotTracksSideBySide: ArchR source pattern not found -> ', u[1],
                                  '\n(ArchR version may have changed; the styling patch needs updating)')
      txt = sub(u[1], u[2], txt)
    }
    f = eval(parse(text = txt)); environment(f) = ns; f
  }

  ## gene track: strand colours, symbol centred above its own gene body
  gt = rewrite(orig$gene, list(
    c('x\\s*=\\s*start,\\s*y\\s*=\\s*cluster,\\s*label\\s*=\\s*symbol',
      'x = (start + end)/2, y = cluster, label = symbol'),
    c('x\\s*=\\s*end,\\s*y\\s*=\\s*cluster,\\s*label\\s*=\\s*symbol',
      'x = (start + end)/2, y = cluster, label = symbol'),
    c('nudge_x\\s*=\\s*-0\\.01\\s*\\*\\s*\\(end\\(region\\)\\s*-\\s*start\\(region\\)\\)', 'nudge_x = 0'),
    c('nudge_x\\s*=\\s*\\+0\\.01\\s*\\*\\s*\\(end\\(region\\)\\s*-\\s*start\\(region\\)\\)', 'nudge_x = 0'),
    c('nudge_y\\s*=\\s*-0\\.25', paste0('nudge_y = ', geneLabelNudge)),
    c('nudge_y\\s*=\\s*0\\.25',  paste0('nudge_y = ', geneLabelNudge)),
    c('ggrepel::geom_label_repel', 'ggrepel::geom_text_repel')))
  formals(gt)$colorPlus = genePlusColor
  formals(gt)$colorMinus = geneMinusColor
  formals(gt)$labelSize = geneLabelSize
  assignInNamespace('.geneTracks', gt, ns = 'ArchR')

  ## loop track: plotBrowserTrack never passes `pal` down, so set the formal default
  if (is.null(loopPal))
    loopPal = grDevices::colorRampPalette(c('#F7F5FA', '#B9A7D4', '#6A51A3', '#3F007D'))(100)
  lt = orig$loop; formals(lt)$pal = loopPal
  assignInNamespace('.loopTracks', lt, ns = 'ArchR')

  ## panel frames come from theme_ArchR, not from plotBrowserTrack(borderWidth = )
  if (!panelBorders)
    assignInNamespace('theme_ArchR', function(...)
      orig$theme(...) + ggplot2::theme(panel.border = ggplot2::element_blank()), ns = 'ArchR')

  ## bulk track: optional outline, and the group name inside each facet
  bulk_variant = function(labelled) {
    subs = list()
    if (outlineProfiles)
      subs = c(subs, list(c('geom_area\\(stat\\s*=\\s*"identity"\\)',
        sprintf('geom_area(stat = "identity", colour = "black", linewidth = %s)', outlineWidth))))
    if (labelled)
      subs = c(subs, list(c('facet_wrap\\(facets\\s*=\\s*~group,\\s*strip\\.position\\s*=\\s*"right",\\s*ncol\\s*=\\s*1\\)',
        paste0('facet_wrap(facets = ~group, strip.position = "right", ncol = 1) + ',
               'geom_text(data = data.frame(group = unique(df$group)), ',
               'mapping = aes(x = -Inf, y = Inf, label = group), hjust = -0.07, vjust = 1.45, size = ',
               groupLabelSize, ', inherit.aes = FALSE, colour = "black")'))))
    if (!length(subs)) orig$bulk else rewrite(orig$bulk, subs)
  }

  ## ---- 5. render ----------------------------------------------------------------------------
  render = function(labelled) {
    assignInNamespace('.bulkTracks', bulk_variant(labelled), ns = 'ArchR')
    p = ArchR::plotBrowserTrack(ArchRProj, region = gr, groupBy = groupBy, useGroups = groups,
                                features = features, loops = loops, plotSummary = plotSummary,
                                sizes = sizes, tileSize = tileSize, minCells = minCells,
                                baseSize = baseSize, facetbaseSize = facetSize, pal = pal, title = '')
    ## a gtable IS a list, so is.list() cannot distinguish 'one region' from 'many':
    ## with a single region plotBrowserTrack returns the gtable itself
    if (inherits(p, 'gtable')) p = list(p)
    p
  }
  if (verbose) message('rendering ', n, ' regions x ', length(groups), ' groups',
                       if (labelFirstOnly) ' (x2 label variants)' else '')
  plain = render(FALSE)
  lab   = if (labelFirstOnly) render(TRUE) else plain

  ## ---- 6. trim and assemble ------------------------------------------------------------------
  ## column 5 = y-axis title, 6 = feature-row label, 8 = facet strip in ArchR's panel gtable
  COL = c(ylab = 5, axis_l = 6, strip = 8)
  blank_cols = function(g, cols) {
    cols = cols[!is.na(cols)]
    if (!length(cols)) return(g)
    for (i in which(g$layout$l %in% cols)) g$grobs[[i]] = grid::nullGrob()
    g$widths[cols] = grid::unit(0, 'cm'); g
  }
  blank_named = function(g, rx, drop = c('width','height')) {
    i = grep(rx, g$layout$name)
    if (length(i)) {
      for (k in i) g$grobs[[k]] = grid::nullGrob()
      if ('width'  %in% drop) g$widths[unique(g$layout$l[i])]  = grid::unit(0, 'cm')
      if ('height' %in% drop) g$heights[unique(g$layout$t[i])] = grid::unit(0, 'cm')
    }
    g
  }
  compact_rows = function(g, keep = grid::unit(0.02, 'cm')) {
    last = suppressWarnings(max(g$layout$b[grepl('^panel-1-', g$layout$name)]))
    if (!is.finite(last)) return(g)
    for (r in seq_len(nrow(g))) {
      if (r <= last) next
      occ = unique(g$layout$name[g$layout$t <= r & g$layout$b >= r])
      if (!length(setdiff(occ, 'background'))) g$heights[r] = keep
    }
    g
  }
  default_header = function(i) {
    lab_i = GenomicRanges::mcols(gr)$label[i]
    sprintf('%s\n%s:%s-%s', sub('^HUB', 'cHub', lab_i), as.character(GenomicRanges::seqnames(gr))[i],
            format(GenomicRanges::start(gr)[i], big.mark = ','),
            format(GenomicRanges::end(gr)[i], big.mark = ','))
  }
  retitle = function(g, txt) {
    i = which(g$layout$name == 'title')
    if (length(i)) {
      g$grobs[[i[1]]] = grid::textGrob(txt, gp = grid::gpar(fontsize = headerSize, lineheight = 1.2))
      g$heights[g$layout$t[i[1]]] = grid::unit(2.2, 'lines')
    }
    g
  }
  prep = function(g, k, m, idx) {
    g = retitle(g, if (is.null(header)) default_header(idx) else header(idx, gr))
    g = blank_named(g, '^guide-box-bottom$', drop = 'height')       # stray strand key
    if (k < m) g = blank_named(g, '^guide-box-right$', drop = 'width')  # loop key: last panel only
    if (!showXaxis) g = blank_named(g, '^(axis-b|xlab-b)', drop = 'height')
    g = compact_rows(g)
    g = blank_cols(g, c(COL[['ylab']], COL[['strip']], if (k > 1) COL[['axis_l']]))
    g$widths[c(1, ncol(g))] = grid::unit(0.13, 'cm')               # gutter between panels
    g
  }
  row_of = function(idx) {
    m = length(idx)
    tr = lapply(seq_along(idx), function(k)
      prep(if (k == 1) lab[[idx[k]]] else plain[[idx[k]]], k, m, idx[k]))
    ## cbind() on a single gtable falls through to the matrix method and returns a matrix,
    ## so the one-region case is returned directly
    if (length(tr) == 1) tr[[1]] else do.call(cbind, c(tr, list(size = 'max')))
  }

  chunks = if (is.null(perRow)) list(seq_len(n)) else split(seq_len(n), ceiling(seq_len(n) / perRow))
  rows   = lapply(chunks, row_of)

  ## ---- 7. draw ---------------------------------------------------------------------------------
  if (!is.null(file)) {
    wdt = panelWidth * max(lengths(chunks)) + labelWidth
    grDevices::pdf(file, width = wdt, height = panelHeight * length(rows) + 0.3)
    grid::grid.newpage()
    if (length(rows) == 1) grid::grid.draw(rows[[1]]) else {
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(length(rows), 1)))
      for (k in seq_along(rows)) {
        grid::pushViewport(grid::viewport(layout.pos.row = k, layout.pos.col = 1))
        grid::grid.draw(rows[[k]]); grid::popViewport()
      }
      grid::popViewport()
    }
    if (!is.null(caption))
      grid::grid.text(caption, x = grid::unit(6, 'mm'), y = grid::unit(3, 'mm'),
                      just = c('left','bottom'), gp = grid::gpar(fontsize = 8.5, col = 'grey25'))
    grDevices::dev.off()
    if (verbose) message('wrote ', file, '  (', round(wdt, 1), ' x ',
                         round(panelHeight * length(rows) + 0.3, 1), ' in)')
  }
  invisible(if (length(rows) == 1) rows[[1]] else rows)
}
