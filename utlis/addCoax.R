#Run This to accesss hidden functions
fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
for (i in seq_along(fn)) {
    tryCatch({
        eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
    }, error = function(x) {
    })
}

addCoAx <- function (
  ArchRProj = NULL,
  KNN = NULL,
  scaleDims = NULL, # this variable never appear in the original code 
  maxDist = 100000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  return_df = FALSE, # return GRanges of peak peak correlation instead of ArchR project obj
  seed = 1, 
  threads = getArchRThreads(),
  force = FALSE

  ){
  coax_parameters = paste0('KNN_number_',length(KNN),'_KNN_size_',length(KNN[[1]]),'_max_dist_',maxDist)
  cat (paste(
    '--* Coaccessibility *--','\n',
    'KNN number:',length(KNN),'\n',
    'KNN size:', length(KNN[[1]]),'\n',
    'Max distance:', max_dist,'\n'
    ))
  if (is.null(ArchRProj@projectMetadata[['coax_parameters']])) ArchRProj@projectMetadata[['coax_parameters']] = 'not_computed'
  if (ArchRProj@projectMetadata[['coax_parameters']] != coax_parameters | force)
    {
    message ('Running coaccessibility')  
    #Get Peak Set
    peakSet = getPeakSet(ArchRProj)

    #Check Chromosomes
    chri <- gtools::mixedsort(.availableChr(getArrowFiles(ArchRProj), subGroup = "PeakMatrix"))
    chrj <- gtools::mixedsort(unique(paste0(seqnames(getPeakSet(ArchRProj)))))
    stopifnot(identical(chri,chrj))

    #Create Ranges
    peakSummits <- resize (peakSet, 1, "center") # get peak summits from peakset
    peakWindows <- resize (peakSummits, maxDist * 2, "center") # get peak windows equal to user-specified max distance of links

    #Create Pairwise Things to Test
    o <- DataFrame (findOverlaps (peakSummits, peakWindows, ignore.strand = TRUE))
    o <- o[o[,1] != o[,2],] # filter overlap between peakSummit and peakWindow of the same peak
    o$seqnames <- seqnames(peakSet)[o[,1]]
    o$idx1 <- peakSet$idx[o[,1]]
    o$idx2 <- peakSet$idx[o[,2]]
    o$correlation <- -999.999
    o$Variability1 <- 0.000
    o$Variability2 <- 0.000

    #Peak Matrix ColSums
    cS <- .getColSums(getArrowFiles(ArchRProj), chri, verbose = FALSE, useMatrix = "PeakMatrix") # this get sum of peaks per cell
    gS <- sapply(seq_along(KNN), function(x) sum(cS[KNN[[x]]], na.rm=TRUE)) # sum of fragments of each knn

    # compute correlation between each peak pair per chromosome
    for(x in seq_along(chri)){
    message (paste('Computing co-accessibility in',chri[x]))  
      #Features
      featureDF <- mcols(peakSet)[BiocGenerics::which(seqnames(peakSet) == chri[x]),]
      featureDF$seqnames <- chri[x]

      #Group Matrix - this is to generate a knn x peak matrix
      groupMat <- .getGroupMatrix(
        ArrowFiles = getArrowFiles (ArchRProj), 
        featureDF = featureDF, 
        groupList = KNN, 
        useMatrix = "PeakMatrix",
        verbose = FALSE
      )    

      #Scale
      groupMat <- t(t(groupMat) / gS) * scaleTo # seq-depth normalisation

      if(log2Norm){
        groupMat <- log2(groupMat + 1)
      }

      #Correlations
      idx <- BiocGenerics::which(o$seqnames == chri[x])
      corVals <- rowCorCpp (idxX = o[idx,]$idx1, idxY = o[idx,]$idx2, X = as.matrix(groupMat), Y = as.matrix(groupMat))

      rowVars <- as.numeric(matrixStats::rowVars(groupMat))

      o[idx,]$correlation <- as.numeric(corVals)
      o[idx,]$Variability1 <- rowVars[o[idx,]$idx1]
      o[idx,]$Variability2 <- rowVars[o[idx,]$idx2]

    }
    
    o$idx1 <- NULL
    o$idx2 <- NULL
    o <- o[!is.na(o$correlation),]

    o$TStat <- (o$correlation / sqrt((max(1-o$correlation^2, 0.00000000000000001))/(length(KNN)-2))) #T-statistic P-value
    o$Pval <- 2*pt(-abs(o$TStat), length(KNN) - 2)
    o$FDR <- p.adjust(o$Pval, method = "fdr")
    o$VarQuantile1 <- .getQuantiles(o$Variability1)
    o$VarQuantile2 <- .getQuantiles(o$Variability2)

    mcols(peakSet) <- NULL
    o@metadata$peakSet <- peakSet

    metadata (ArchRProj@peakSet)$CoAccessibility <- o
    ArchRProj@projectMetadata[['coax_parameters']] = coax_parameters

    if (return_df == TRUE){
    return (o)  
    } else {
    return (ArchRProj)
    }
    } else {
    stop ('coaccessibility with same parameters already computed!')  
    }
}

