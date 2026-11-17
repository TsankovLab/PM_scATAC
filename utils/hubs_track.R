#' Plot an ArchR Region Track
#' 
#' This function will plot the coverage at an input region in the style of a browser track. It allows for normalization of the signal
#' which enables direct comparison across samples.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param region A `GRanges` region that indicates the region to be plotted. If more than one region exists in the `GRanges` object,
#' all will be plotted. If no region is supplied, then the `geneSymbol` argument can be used to center the plot window at the
#' transcription start site of the supplied gene.
#' @param groupBy A string that indicates how cells should be grouped. This string corresponds to one of the standard or
#' user-supplied `cellColData` metadata columns (for example, "Clusters"). Cells with the same value annotated in this metadata
#' column will be grouped together and the average signal will be plotted.
#' @param useGroups A character vector that is used to select a subset of groups by name from the designated `groupBy` column in
#' `cellColData`. This limits the groups to be plotted.
#' @param plotSummary A character vector containing the features to be potted. Possible values include "bulkTrack" (the ATAC-seq signal),
#' "scTrack" (scATAC-seq signal), "featureTrack" (i.e. the peak regions), "geneTrack" (line diagrams of genes with introns and exons shown. 
#' Blue-colored genes are on the minus strand and red-colored genes are on the plus strand), and "loopTrack" (links between a peak and a gene).
#' @param sizes A numeric vector containing up to 3 values that indicate the sizes of the individual components passed to `plotSummary`.
#' The order must be the same as `plotSummary`.
#' @param features A `GRanges` object containing the "features" to be plotted via the "featureTrack". This should be thought of as a
#' bed track. i.e. the set of peaks obtained using `getPeakSet(ArchRProj))`. 
#' @param loops A `GRanges` object containing the "loops" to be plotted via the "loopTrack".
#' This `GRanges` object start represents the center position of one loop anchor and the end represents the center position of another loop anchor. 
#' A "loopTrack" draws an arc between two genomic regions that show some type of interaction. This type of track can be used 
#' to display chromosome conformation capture data or co-accessibility links obtained using `getCoAccessibility()`. 
#' @param geneSymbol If `region` is not supplied, plotting can be centered at the transcription start site corresponding to the gene symbol(s) passed here.
#' @param useMatrix If supplied geneSymbol, one can plot the corresponding GeneScores/GeneExpression within this matrix. I.E. "GeneScoreMatrix"
#' @param log2Norm If supplied geneSymbol, Log2 normalize the corresponding GeneScores/GeneExpression matrix before plotting.
#' @param upstream The number of basepairs upstream of the transcription start site of `geneSymbol` to extend the plotting window.
#' If `region` is supplied, this argument is ignored.
#' @param downstream The number of basepairs downstream of the transcription start site of `geneSymbol` to extend the plotting window.
#' If `region` is supplied, this argument is ignored.
#' @param tileSize The numeric width of the tile/bin in basepairs for plotting ATAC-seq signal tracks. All insertions in a single bin will be summed.
#' @param minCells The minimum number of cells contained within a cell group to allow for this cell group to be plotted. This argument can be
#' used to exclude pseudo-bulk replicates generated from low numbers of cells.
#' @param normMethod The name of the column in `cellColData` by which normalization should be performed. The recommended and default value
#' is "ReadsInTSS" which simultaneously normalizes tracks based on sequencing depth and sample data quality.
#' @param threads The number of threads to use for parallel execution.
#' @param ylim The numeric quantile y-axis limit to be used for for "bulkTrack" plotting. If not provided, the y-axis limit will be c(0, 0.999).
#' @param pal A custom palette (see `paletteDiscrete` or `ArchRPalettes`) used to override coloring for groups.
#' @param baseSize The numeric font size to be used in the plot. This applies to all plot labels.
#' @param scTileSize The width of the tiles in scTracks. Larger numbers may make cells overlap more. Default is 0.5 for about 100 cells.
#' @param scCellsMax The maximum number of cells for scTracks.
#' @param borderWidth The numeric line width to be used for plot borders.
#' @param tickWidth The numeric line width to be used for axis tick marks.
#' @param facetbaseSize The numeric font size to be used in the facets (gray boxes used to provide track labels) of the plot.
#' @param geneAnnotation The `geneAnnotation` object to be used for plotting the "geneTrack" object. See `createGeneAnnotation()` for more info.
#' @param title The title to add at the top of the plot next to the plot's genomic coordinates.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
plotBrowserTrack2 <- function(
  ArchRProj = NULL, 
  region = NULL, 
  groupBy = "Clusters",
  useGroups = NULL, 
  sample_levels = NULL,
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack", "hubTrack",'hubregiontrack'),
  sizes = c(10, 1.5, 3, 4,3, 3),
  features = getPeakSet(ArchRProj),
  loops = getCoAccessibility(ArchRProj),
  loop_size = 0.5,
  hubs = NULL,
  hubs_regions = NULL,
  geneSymbol = NULL,
  useMatrix = NULL,
  log2Norm = TRUE,
  upstream = 50000,
  downstream = 50000,
  genelabelsize = 2,
  tileSize = 250, 
  minCells = 25,
  normMethod = "ReadsInTSS",
  threads = getArchRThreads(), 
  ylim = NULL,
  pal = NULL,
  baseSize = 7,
  scTileSize = 0.5,
  scCellsMax = 100,
  borderWidth = 0.3,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  verbose = TRUE,
  logFile = createLogFile("plotBrowserTrack2")
  ){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  .validInput(input = region, name = "region", valid = c("granges","null"))
  .validInput(input = groupBy, name = "groupBy", valid = "character")
  .validInput(input = useGroups, name = "useGroups", valid = c("character", "null"))
  .validInput(input = plotSummary, name = "plotSummary", valid = "character")
  .validInput(input = sizes, name = "sizes", valid = "numeric")
  .validInput(input = features, name = "features", valid = c("granges", "grangeslist", "null"))
  .validInput(input = loops, name = "loops", valid = c("granges", "grangeslist", "null"))
  .validInput(input = geneSymbol, name = "geneSymbol", valid = c("character", "null"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character", "null"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = upstream, name = "upstream", valid = c("integer"))
  .validInput(input = downstream, name = "downstream", valid = c("integer"))
  .validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  .validInput(input = minCells, name = "minCells", valid = c("integer"))
  .validInput(input = normMethod, name = "normMethod", valid = c("character"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = ylim, name = "ylim", valid = c("numeric", "null"))
  .validInput(input = pal, name = "pal", valid = c("palette", "null"))
  .validInput(input = baseSize, name = "baseSize", valid = "numeric")
  .validInput(input = scTileSize, name = "scTileSize", valid = "numeric")
  .validInput(input = scCellsMax, name = "scCellsMax", valid = "integer")
  .validInput(input = borderWidth, name = "borderWidth", valid = "numeric")
  .validInput(input = tickWidth, name = "tickWidth", valid = "numeric")
  .validInput(input = facetbaseSize, name = "facetbaseSize", valid = "numeric")
  geneAnnotation <- .validGeneAnnotation(geneAnnotation)
  .validInput(input = title, name = "title", valid = "character")

  tstart <- Sys.time()
  .startLogging(logFile=logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "plotBrowserTrack2 Input-Parameters", logFile = logFile)

  ##########################################################
  # Get Region Where Plot Will Occur (GenomicRanges)
  ##########################################################
  .logDiffTime("Validating Region", t1=tstart, verbose=verbose, logFile=logFile)
  if(is.null(region)){
    if(!is.null(geneSymbol)){
      region <- geneAnnotation$genes
      region <- region[which(tolower(mcols(region)$symbol) %in% tolower(geneSymbol))]
      region <- region[order(match(tolower(mcols(region)$symbol), tolower(geneSymbol)))]
      print(region)
      region <- resize(region, 1, "start")
      strand(region) <- "*"
      region <- extendGR(region, upstream = upstream, downstream = downstream)
    }
  }
  region <- .validGRanges(region)
  .logThis(region, "region", logFile = logFile)

  if(is.null(geneSymbol)){
    useMatrix <- NULL
  }

  if(!is.null(useMatrix)){
    featureMat <- .getMatrixValues(
      ArchRProj = ArchRProj,
      matrixName = useMatrix,
      name = mcols(region)$symbol
    )
    if(log2Norm){
      featureMat <- log2(featureMat + 1) 
    }
    featureMat <- data.frame(t(featureMat))
    featureMat$Group <- getCellColData(ArchRProj, groupBy, drop = FALSE)[rownames(featureMat), 1]
  }

  ggList <- lapply(seq_along(region), function(x){

    plotList <- list()

    ##########################################################
    # Bulk Tracks
    ##########################################################
    if("bulktrack" %in% tolower(plotSummary)){
      .logDiffTime(sprintf("Adding Bulk Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$bulktrack <- .bulkTracks(
        ArchRProj = ArchRProj, 
        region = region[x], 
        tileSize = tileSize, 
        groupBy = groupBy,
        sample_levels = sample_levels,
        threads = threads, 
        minCells = minCells,
        pal = pal,
        ylim = ylim,
        baseSize = baseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth,
        facetbaseSize = facetbaseSize,
        normMethod = normMethod,
        geneAnnotation = geneAnnotation,
        title = title,
        useGroups = useGroups,
        tstart = tstart,
        logFile = logFile) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    }
    
    ##########################################################
    # Bulk Tracks
    ##########################################################
    if("sctrack" %in% tolower(plotSummary)){
      .logDiffTime(sprintf("Adding SC Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$sctrack <- .scTracks(
        ArchRProj = ArchRProj, 
        region = region[x], 
        tileSize = tileSize, 
        groupBy = groupBy,
        threads = threads, 
        minCells = 5,
        maxCells = scCellsMax,
        pal = pal,
        baseSize = baseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth,
        scTileSize = scTileSize,
        facetbaseSize = facetbaseSize,
        geneAnnotation = geneAnnotation,
        title = title,
        useGroups = useGroups,
        tstart = tstart,
        logFile = logFile) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    }

    ##########################################################
    # Feature Tracks
    ##########################################################
    if("featuretrack" %in% tolower(plotSummary)){
      if(!is.null(features)){
        .logDiffTime(sprintf("Adding Feature Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
        plotList$featuretrack <- .featureTracks(
            features = features, 
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            title = "Peaks",
            logFile = logFile) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
      }
    }

    ##########################################################
    # hubregion Tracks
    ##########################################################
    if("hubregiontrack" %in% tolower(plotSummary)){
      if(!is.null(hubs_regions)){
        .logDiffTime(sprintf("Adding hubregion Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
        plotList$hubregiontrack <- .hubregionTracks(
            features = hubs_regions, 
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            title = "hub region",
            logFile = logFile) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
      }
    }

    ##########################################################
    # Loops Tracks
    ##########################################################
    if("looptrack" %in% tolower(plotSummary)){
      if(!is.null(loops)){
        .logDiffTime(sprintf("Adding Loop Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
        plotList$looptrack <- .loopTracks(
            loops = loops, 
            loop_size = loop_size,
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            hideY = TRUE,
            title = "Loops",
            logFile = logFile) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
      }
    }
     ##########################################################
    # Hubs Tracks
    ##########################################################
    if("hubtrack" %in% tolower(plotSummary)){
      if(!is.null(hubs)){
        .logDiffTime(sprintf("Adding Hubs Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
        plotList$hubTrack <- hubsTracks(
            hubs = hubs, 
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            hideY = TRUE,
            title = "Hubs",
            logFile = logFile) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
      }
    }

    ##########################################################
    # Gene Tracks
    ##########################################################
    if("genetrack" %in% tolower(plotSummary)){
      .logDiffTime(sprintf("Adding Gene Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$genetrack <- .geneTracks(
        geneAnnotation = geneAnnotation,
        labelSize=genelabelsize, 
        region = region[x], 
        facetbaseSize = facetbaseSize,
        title = "Genes",
        logFile = logFile) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    }

    ##########################################################
    # Time to plot
    ##########################################################
    plotSummary <- tolower(plotSummary)
    names(sizes) <- plotSummary
    sizes <- sizes[order(plotSummary)]
    plotSummary <- plotSummary[order(plotSummary)]

    # nullSummary <- unlist(lapply(seq_along(plotSummary), function(x) is.null(eval(parse(text=paste0("plotList$", plotSummary[x]))))))
    # if(any(nullSummary)){
    #   sizes <- sizes[-which(nullSummary)]
    # }
    sizes <- sizes[tolower(names(plotList))]

    if(!is.null(useMatrix)){

      suppressWarnings(.combinedFeaturePlot(
        plotList = plotList,
        log2Norm = log2Norm,
        featureMat = featureMat,
        feature = region[x]$symbol[[1]],
        useMatrix = useMatrix,
        pal = pal,
        sizes = sizes,
        baseSize = baseSize,
        facetbaseSize = facetbaseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth
      ))

    }else{

      .logThis(names(plotList), sprintf("(%s of %s) names(plotList)",x,length(region)), logFile=logFile)
      .logThis(sizes, sprintf("(%s of %s) sizes",x,length(region)), logFile=logFile)
      #.logThis(nullSummary, sprintf("(%s of %s) nullSummary",x,length(region)), logFile=logFile)
      .logDiffTime("Plotting", t1=tstart, verbose=verbose, logFile=logFile)
      
      tryCatch({
        suppressWarnings(ggAlignPlots(plotList = plotList, sizes=sizes, draw = FALSE))
      }, error = function(e){
        .logMessage("Error with plotting, diagnosing each element", verbose = TRUE, logFile = logFile)
        for(i in seq_along(plotList)){
          tryCatch({
            print(plotList[[i]])
          }, error = function(f){
            #.logError(f, fn = names(plotList)[i], info = "", errorList = NULL, logFile = logFile)
          })
        }
        .logError(e, fn = "ggAlignPlots", info = "", errorList = NULL, logFile = logFile)
      })

    }

  })

  if(!is.null(mcols(region)$symbol)){
    names(ggList) <- mcols(region)$symbol
  }else{
    if(length(ggList) == 1){
      ggList <- ggList[[1]]
    }
  }

  .endLogging(logFile=logFile)

  ggList

}
    
#######################################################
# Bulk Aggregated ATAC Track Methods
#######################################################
.bulkTracks <- function(
  ArchRProj = NULL, 
  region = NULL, 
  tileSize = 100, 
  minCells = 25,
  groupBy = "Clusters",
  sample_levels = NULL,
  useGroups = NULL,
  normMethod = "ReadsInTSS",
  threads = 1, 
  ylim = NULL,
  baseSize = 7,
  borderWidth = 0.3,
  tickWidth = 0,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  pal = NULL,
  tstart = NULL,
  verbose = FALSE,
  logFile = NULL
  ){

  .requirePackage("ggplot2", source = "cran")

  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  df <- .groupRegionSumArrows(
    ArchRProj = ArchRProj, 
    groupBy = groupBy, 
    normMethod = normMethod,
    useGroups = useGroups,
    minCells = minCells,
    region = region, 
    tileSize = tileSize, 
    threads = threads,
    verbose = verbose,
    logFile = logFile
  )
  .logThis(split(df, df[,3]), ".bulkTracks df", logFile = logFile)

  ######################################################
  # Plot Track
  ######################################################
  if(!is.null(ylim)){
    df$y[df$y < ylim[1]] <- ylim[1]
    df$y[df$y > ylim[2]] <- ylim[2]
  }else{
    ylim <- c(0,quantile(df$y, probs=c(0.999)))
    df$y[df$y < ylim[1]] <- ylim[1]
    df$y[df$y > ylim[2]] <- ylim[2]
  }
  uniqueGroups <- gtools::mixedsort(unique(paste0(df$group)))
  if(!is.null(useGroups)){
    uniqueGroups <- unique(useGroups)
  }
  df$group <- factor(df$group, levels = uniqueGroups)
  if (!is.null (sample_levels))  df$group = factor(df$group, levels = sample_levels[sample_levels %in% levels(df$group)])
  title <- paste0(as.character(seqnames(region)),":", start(region)-1, "-", end(region), " ", title)

  allGroups <- gtools::mixedsort(unique(getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = TRUE)))

  if(is.null(pal)){
    pal <- suppressWarnings(paletteDiscrete(values = allGroups))
  }
  
  #Plot Track
  p <- ggplot(df, aes_string("x","y", fill = "group")) + 
    geom_area(stat = "identity", color = "black", size = 0.2) + 
    facet_wrap(facets = ~group, strip.position = 'right', ncol = 1) +
    ylab(sprintf("Coverage\n(Norm. ATAC Signal Range (%s-%s) by %s)", round(min(ylim),2), round(max(ylim),2), normMethod)) +
    #scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
    scale_y_continuous(limits = ylim, expand = c(0,0)) +
    theme_ArchR(baseSize = baseSize,
                baseRectSize = borderWidth,
                baseLineSize = tickWidth,
                legendPosition = "right",
                axisTickCm = 0) +
    theme(panel.spacing= unit(0, "lines"),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text = element_text(
            size = facetbaseSize, 
            color = "black", 
            margin = margin(0,0,0,0, "cm")),
            strip.text.y = element_text(angle = 0),
          strip.background = element_rect(color="black")) +
    guides(fill = FALSE, colour = FALSE) + ggtitle(title)

  p

}

##############################################################################
# Create Average Tracks from Arrows
##############################################################################
.groupRegionSumArrows <- function(
  ArchRProj = NULL,
  useGroups = NULL,
  groupBy = NULL,
  region = NULL,
  tileSize = NULL,
  normMethod = NULL,
  verbose = FALSE,
  minCells = 25,
  maxCells = Inf, # Changed from 500 
  threads = NULL,
  logFile = NULL
  ){
  
  #Group Info
  cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(!is.null(minCells)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% names(table(cellGroups)[table(cellGroups) >= minCells]),,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  if(!is.null(useGroups)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% useGroups,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  tabGroups <- table(cellGroups)
  
  if(any(tabGroups > maxCells)){
    cellGroups2 <- getCellColData(ArchRProj, groupBy, drop = FALSE)
    splitGroups <- split(rownames(cellGroups2), cellGroups2[,1])
    useCells <- lapply(seq_along(splitGroups), function(x){
      if(length(splitGroups[[x]]) > maxCells){
        sample(splitGroups[[x]], maxCells)
      }else{
        splitGroups[[x]]
      }
    }) %>% unlist
    ArchRProj@cellColData <- ArchRProj@cellColData[useCells,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
    tabGroups <- table(cellGroups)
  }

  cellsBySample <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj, "Sample", drop = TRUE))
  groupsBySample <- split(cellGroups, getCellColData(ArchRProj, "Sample", drop = TRUE))
  uniqueGroups <- gtools::mixedsort(unique(cellGroups))
  
  #Tile Region
  regionTiles <- seq(trunc(start(region) / tileSize), trunc(end(region) / tileSize) + 1) * tileSize
  ArrowFiles <- getArrowFiles(ArchRProj)
  ArrowFiles <- ArrowFiles[names(cellsBySample)]

  groupMat <- .safelapply(seq_along(ArrowFiles), function(i){
    .logMessage(sprintf("Getting Region From Arrow Files %s of %s", i, length(ArrowFiles)), logFile = logFile)
    tryCatch({
      .regionSumArrows(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
    }, error = function(e){
      errorList <- list(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
      .logError(e, fn = ".groupRegionSumArrows", info = .sampleName(ArrowFiles[i]), errorList = errorList, logFile = logFile)
    })
  }, threads = threads) %>% Reduce("+" , .)

  #Plot DF
  df <- data.frame(which(groupMat > 0, arr.ind=TRUE))
  df$y <- groupMat[cbind(df[,1], df[,2])]

  #Minus 1 Tile Size
  dfm1 <- df
  dfm1$row <- dfm1$row - 1
  dfm1$y <- 0

  #Plus 1 Size
  dfp1 <- df
  dfp1$row <- dfp1$row + 1
  dfp1$y <- 0

  #Create plot DF
  df <- rbind(df, dfm1, dfp1)
  df <- df[!duplicated(df[,1:2]),]
  df <- df[df$row > 0,]
  df$x <- regionTiles[df$row]
  df$group <- uniqueGroups[df$col]

  #Add In Ends
  dfs <- data.frame(
    col = seq_along(uniqueGroups), 
    row = 1, 
    y = 0,
    x = start(region),
    group = uniqueGroups
  )

  dfe <- data.frame(
    col = seq_along(uniqueGroups),
    row = length(regionTiles),
    y = 0,
    x = end(region),
    group = uniqueGroups
  )
  
  #Final output
  plotDF <- rbind(df,dfs,dfe)
  plotDF <- df[order(df$group,df$x),]
  plotDF <- df[,c("x", "y", "group")]
  
  #Normalization 
  g <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(tolower(normMethod) == "readsintss"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "readsinpromoter"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "nfrags"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "ncells"){
      groupNormFactors <- table(g)
  }else if(tolower(normMethod) == "none"){
      groupNormFactors <- rep(10^4, length(g))
      names(groupNormFactors) <- g
  }else{
    stop("Norm Method Not Recognized : ", normMethod)
  }

  #Scale with Norm Factors
  scaleFactors <- 10^4 / groupNormFactors
  matchGroup <- match(paste0(plotDF$group), names(scaleFactors))
  plotDF$y <- plotDF$y * as.vector(scaleFactors[matchGroup])

  return(plotDF)

}

.regionSumArrows <- function(
  ArrowFile = NULL,
  region = NULL,
  regionTiles = NULL,
  tileSize = NULL,
  cellNames = NULL,
  cellGroups = NULL,
  uniqueGroups = NULL,
  logFile = NULL
  ){
  
  cellFragsRegion <- .getFragsFromArrow(
      ArrowFile = ArrowFile, 
      chr = paste0(seqnames(region)), 
      cellNames = cellNames, 
      out = "GRanges"
    ) %>% subsetByOverlaps(., region, ignore.strand = TRUE)
  
  #Starts
  ts <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ids <- which(ts > 0)
  
  #Ends
  te <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ide <- which(te > 0)
  
  #Match
  matchID <- S4Vectors::match(mcols(cellFragsRegion)$RG, cellNames)
  
  #Sparse Matrix
  mat <- Matrix::sparseMatrix(
    i = c(ts[ids], te[ide]),
    j = c(matchID[ids], matchID[ide]),
    x = rep(1,  length(ids) + length(ide)),
    dims = c(length(regionTiles), length(cellNames))
  )
  colnames(mat) <- cellNames
  
  mat@x[mat@x > 1] <- 1

  #Create Group Matrix
  groupMat <- matrix(0, nrow = length(regionTiles), ncol = length(uniqueGroups))
  colnames(groupMat) <- uniqueGroups
  uniqueGroups <- uniqueGroups[uniqueGroups %in% unique(cellGroups)]
  for(i in seq_along(uniqueGroups)){
    groupMat[,uniqueGroups[i]] <- Matrix::rowSums(mat[,which(cellGroups == uniqueGroups[i]),drop=FALSE])
  }

  return(groupMat)

}

#######################################################
# Gene Tracks
#######################################################
.geneTracks <- function(
  geneAnnotation = NULL, 
  region = NULL, 
  baseSize = 9, 
  borderWidth = 0.3, 
  title = "Genes",
  geneWidth = 1, 
  exonWidth = 1.5, 
  labelSize = 2,
  facetbaseSize,
  colorMinus = "grey22",
  colorPlus = "brown",
  logFile = NULL
  ){

  .requirePackage("ggplot2", source = "cran")
  .requirePackage("ggrepel", source = "cran")

  #only take first region
  region <- .validGRanges(region)
  region <- .subsetSeqnamesGR(region[1], as.character(seqnames(region[1])))

  genes <- sort(sortSeqlevels(geneAnnotation$genes), ignore.strand = TRUE)
  exons <- sort(sortSeqlevels(geneAnnotation$exons), ignore.strand = TRUE)
  genesO <- data.frame(subsetByOverlaps(genes, region, ignore.strand = TRUE))

  if(nrow(genesO) > 0){

    #Identify Info for Exons and Genes
    exonsO <- data.frame(subsetByOverlaps(exons, region, ignore.strand = TRUE))
    exonsO <- exonsO[which(exonsO$symbol %in% genesO$symbol),]
    genesO$facet = title
    genesO$start <- matrixStats::rowMaxs(cbind(genesO$start, start(region)))
    genesO$end <- matrixStats::rowMins(cbind(genesO$end, end(region)))

    #Collapse Iteratively
    #backwards iteration so that the last value chosen is the lowest cluster possible to fit in.
    genesO$cluster <- 0
    for(i in seq_len(nrow(genesO))){
      if(i==1){
        genesO$cluster[i] <- 1
      }else{
        for(j in seq_len(max(genesO$cluster))){
          jEnd <- rev(genesO$end)[match(rev(seq_len(max(genesO$cluster)))[j], rev(genesO$cluster))]
          if(genesO$start[i] > jEnd + median(genesO$width)){
            genesO$cluster[i] <- rev(genesO$cluster)[match(rev(seq_len(max(genesO$cluster)))[j],rev(genesO$cluster))]
          }
        }
        if(genesO$cluster[i]==0){
          genesO$cluster[i] <- genesO$cluster[i-1] + 1
        }
      }
    }
    exonsO$cluster <- genesO$cluster[match(exonsO$symbol, genesO$symbol)]
    pal <- c("-"=colorMinus,"+"=colorPlus,"*"=colorPlus)
    
    p <- ggplot(data = genesO, aes(color = strand, fill = strand)) +
      facet_grid(facet~.) +
      #################################################
      #Limits
      #################################################
      ylim(c(0.5, max(genesO$cluster) + 0.5)) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) + 
      #################################################
      #Segment for Not Minus Stranded
      #################################################
      geom_segment(data = genesO[which(as.character(genesO$strand)!="-"),], 
        aes(x = start, xend = end, y = cluster, yend = cluster, color = strand),size=geneWidth) +
      #################################################
      #Segment for Minus Stranded
      #################################################
      geom_segment(data = genesO[which(as.character(genesO$strand)=="-"),], 
        aes(x = end, xend = start, y = cluster, yend = cluster, color = strand),size=geneWidth) +
      #################################################
      #Segement for Exons
      #################################################
      geom_segment(data = exonsO, aes(x = start, xend = end, y = cluster, 
        yend = cluster, color = strand),size=exonWidth) +
      #################################################
      #Colors
      #################################################
      scale_color_manual(values = pal, guide = FALSE) + 
      scale_fill_manual(values = pal) +
      #################################################
      #Theme
      #################################################
      #theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      theme_void() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
      theme(legend.text = element_text(size = baseSize), strip.text.y = element_text(size = facetbaseSize, angle = 0)) +
      guides(fill = guide_legend(override.aes = list(colour = NA, shape = "c", size=3)), color = FALSE) + 
      theme(legend.position="bottom") +
      theme(legend.title=element_text(size=5), legend.text=element_text(size=7),
        legend.key.size = unit(0.75,"line"), legend.background = element_rect(color =NA), strip.background = element_blank())

    #Add Labels if There are Genes with this orientation!
    if(length(which(genesO$strand!="-")) > 0){
      p <- p + ggrepel::geom_text_repel(data=genesO[which(genesO$strand!="-"),], 
        aes(x = start, y = cluster, label = symbol, color = strand),fill = 'white', 
          segment.color = "grey", nudge_x = -0.01*(end(region) - start(region)), nudge_y = -0.25, 
          size = labelSize, direction = "x")
    }

    #Add Labels if There are Genes with this orientation!
    if(length(which(genesO$strand=="-")) > 0){
      p <- p + ggrepel::geom_text_repel(data=genesO[which(genesO$strand=="-"),], 
        aes(x = end, y = cluster, label = symbol, color = strand),,fill = 'white', 
          segment.color = "grey", nudge_x = +0.01*(end(region) - start(region)), nudge_y = 0.25, 
          size = labelSize, direction = "x")
    }

    p <- p + theme(legend.justification = c(0, 1), 
      legend.background = element_rect(colour = NA, fill = NA), legend.position="none")

  }else{

    #create empty plot
    df <- data.frame(facet = "GeneTrack", start = 0, end = 0, strand = "*", symbol = "none")
    pal <- c("*"=colorPlus)
    p <- ggplot(data = df, aes(start, end, fill = strand)) + geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_color_manual(values = pal) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  if(!is.ggplot(p)){
    .logError("geneTrack is not a ggplot!", fn = ".geneTracks", info = "", errorList = NULL, logFile = logFile)
  }

  return(p)

}

#######################################################
# Feature Tracks
#######################################################
.featureTracks <- function(
  features = NULL, 
  region = NULL, 
  title = "FeatureTrack", 
  pal = NULL,
  baseSize = 9, 
  facetbaseSize = NULL,
  featureWidth = 2, 
  borderWidth = 0.3, 
  hideX = FALSE, 
  hideY = FALSE,
  logFile = NULL
  ){

  .requirePackage("ggplot2", source = "cran")

  #only take first region
  region <- .validGRanges(region)
  region <- .subsetSeqnamesGR(region[1], as.character(seqnames(region[1])))

  if(!is.null(features)){

    if(!.isGRList(features)){
      features <- .validGRanges(features)
      featureList <- SimpleList(FeatureTrack = features)
      hideY <- TRUE
    }else{
      featureList <- features
      hideY <- FALSE
    }
    featureList <- featureList[rev(seq_along(featureList))]

    featureO <- lapply(seq_along(featureList), function(x){
      featurex <- featureList[[x]]
      namex <- names(featureList)[x]
      mcols(featurex) <- NULL
      sub <- subsetByOverlaps(featurex, region, ignore.strand = TRUE)
      if(length(sub) > 0){
        data.frame(sub, name = namex)
      }else{
        empty <- GRanges(as.character(seqnames(region[1])), ranges = IRanges(0,0))
        data.frame(empty, name = namex)
      }

    })

    featureO <- Reduce("rbind", featureO)
    
    .logThis(featureO, "featureO", logFile = logFile)

    featureO$facet <- title

    if(is.null(pal)){
      pal <- paletteDiscrete(set = "stallion", values = rev(unique(paste0(featureO$name))))
    }
    
    featureO$name <- factor(paste0(featureO$name), levels=names(featureList))

    p <- ggplot(data = featureO, aes(color = name)) +
      facet_grid(facet~.) +
      geom_segment(data = featureO, aes(x = start, xend = end, y = name, yend = name, color = name), size=featureWidth) +
      ylab("") + xlab("") + 
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      scale_color_manual(values = pal) +
      theme_void() + 
      theme(legend.text = element_text(size = baseSize)) + 
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = 0) +
      guides(color = FALSE, fill = FALSE) + theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank())

  }else{

    #create empty plot
    df <- data.frame(facet = "FeatureTrack", start = 0, end = 0, strand = "*", symbol = "none")
    p <- ggplot(data = df, aes(start, end)) + 
      geom_point() +
      facet_grid(facet~.) +
      #theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      theme_void() +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  if(hideX){
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  if(hideY){
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  if(!is.ggplot(p)){
    .logError("featureTrack is not a ggplot!", fn = ".featureTracks", info = "", errorList = NULL, logFile = logFile)
  }

  return(p)

}


#######################################################
# hub region Tracks
#######################################################
.hubregionTracks <- function(
  features = NULL, 
  region = NULL, 
  title = "hubregiontrack", 
  pal = NULL,
  baseSize = 9, 
  facetbaseSize = NULL,
  featureWidth = 2, 
  borderWidth = 0.3, 
  hideX = FALSE, 
  hideY = FALSE,
  logFile = NULL
  ){

  .requirePackage("ggplot2", source = "cran")

  #only take first region
  region <- .validGRanges(region)
  region <- .subsetSeqnamesGR(region[1], as.character(seqnames(region[1])))

  if(!is.null(features)){

    if(!.isGRList(features)){
      features <- .validGRanges(features)
      featureList <- SimpleList(FeatureTrack = features)
      hideY <- TRUE
    }else{
      featureList <- features
      hideY <- FALSE
    }
    featureList <- featureList[rev(seq_along(featureList))]

    featureO <- lapply(seq_along(featureList), function(x){
      featurex <- featureList[[x]]
      namex <- names(featureList)[x]
      mcols(featurex) <- NULL
      sub <- subsetByOverlaps(featurex, region, ignore.strand = TRUE)
      if(length(sub) > 0){
        data.frame(sub, name = namex)
      }else{
        empty <- GRanges(as.character(seqnames(region[1])), ranges = IRanges(0,0))
        data.frame(empty, name = namex)
      }

    })

    featureO <- Reduce("rbind", featureO)
    
    .logThis(featureO, "featureO", logFile = logFile)

    featureO$facet <- title

    if(is.null(pal)){
      pal <- paletteDiscrete(set = "stallion", values = rev(unique(paste0(featureO$name))))
    }
    
    featureO$name <- factor(paste0(featureO$name), levels=names(featureList))

    p <- ggplot(data = featureO, aes(color = name)) +
      facet_grid(facet~.) +
      geom_segment(data = featureO, aes(x = start, xend = end, y = name, yend = name, color = name), size=featureWidth) +
      ylab("") + xlab("") + 
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      scale_color_manual(values = pal) +
      theme_void() + 
      theme(legend.text = element_text(size = baseSize)) + 
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = 0) +
      guides(color = FALSE, fill = FALSE) + theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank())

  }else{

    #create empty plot
    df <- data.frame(facet = "hubregionTrack", start = 0, end = 0, strand = "*", symbol = "none")
    p <- ggplot(data = df, aes(start, end)) + 
      geom_point() +
      facet_grid(facet~.) +
      #theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      theme_void() +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  if(hideX){
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  if(hideY){
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  if(!is.ggplot(p)){
    .logError("hubregionTrack is not a ggplot!", fn = "hubregionTrack", info = "", errorList = NULL, logFile = logFile)
  }

  return(p)

}

#######################################################
# Loop Tracks
#######################################################
.loopTracks <- function(
  loops = NULL, 
  region = NULL, 
  title = "LoopTrack", 
  pal = NULL,
  baseSize = 9, 
  loop_size = NULL,
  facetbaseSize = 9,
  featureWidth = .5, 
  borderWidth = 0.4, 
  hideX = FALSE, 
  hideY = FALSE,
  logFile = NULL
  ){

  getArchDF <- function(lp, r = 100){
    angles <- seq(pi, 2*pi,length.out=100)
    rx <- (end(lp)-start(lp))/2
    rscale <- r * (rx/max(rx))
    cx <- start(lp) + rx
    if(is.null(mcols(lp)$value)){
      mcols(lp)$value <- 1
    }
    df <- lapply(seq_along(cx), function(z){
      xz <- rx[z]*cos(angles)+cx[z]
      dfz <- DataFrame(x=xz, y=rscale[z]*sin(angles), id=Rle(paste0("l",z)), value = mcols(lp)$value[z])
    }) %>% Reduce("rbind",.)
    return(df)
  }

  if(!is.null(loops)){

    if(is(loops, "GRanges")){
      loops <- SimpleList(Loops = loops)
    }else if(.isGRList(loops)){
    }else{
      stop("Loops is not a GRanges or a list of GRanges! Please supply valid input!")
    }

    valueMin <- min(unlist(lapply(loops, function(x) min(x$value))))
    valueMax <- max(unlist(lapply(loops, function(x) max(x$value))))

    loopO <- lapply(seq_along(loops), function(x){
       subLoops <- subsetByOverlaps(loops[[x]], region, ignore.strand = TRUE, type = "within") 
       if(length(subLoops)>0){
         dfx <- getArchDF(subLoops)
         dfx$name <- Rle(paste0(names(loops)[x]))
         dfx
       }else{
         NULL
       }
    }) %>% Reduce("rbind",.)
    .logThis(loopO, "loopO", logFile = logFile)

    testDim <- tryCatch({
      if(is.null(loopO)){
        FALSE
      }
      if(nrow(loopO) > 0){
        TRUE
      }else{
        FALSE
      }
    }, error = function(x){
      FALSE
    })

    if(testDim){

      loopO$facet <- title
      if(is.null(pal)){
        pal <- colorRampPalette(c("white",'purple',"black"))(100)
      }

      p <- ggplot(data = data.frame(loopO), aes(x = x, y = y, group = id, color = value)) + 
        geom_line(size=loop_size) +
        facet_grid(name ~ .) +
        ylab("") + 
        coord_cartesian(ylim = c(-100,0)) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        scale_color_gradientn(colors = pal, limits = c(valueMin, valueMax)) +
        theme(legend.text = element_text(size = baseSize)) +
        theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth, legendPosition = "right") +
        theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank(),
          legend.box.background = element_rect(color = NA)) +
        guides(color= guide_colorbar(barwidth = 0.75, barheight = 3))

    }else{

      #create empty plot
      df <- data.frame(facet = "LoopTrack", start = 0, end = 0, strand = "*", symbol = "none")
      p <- ggplot(data = df, aes(start, end)) + 
        geom_point() +
        facet_grid(facet~.) +
        theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
        theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

    }

  }else{

    #create empty plot
    df <- data.frame(facet = "LoopTrack", start = 0, end = 0, strand = "*", symbol = "none")
    p <- ggplot(data = df, aes(start, end)) + 
      geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  if(hideX){
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  if(hideY){
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  if(!is.ggplot(p)){
    .logError("loopTracks is not a ggplot!", fn = ".loopTracks", info = "", errorList = NULL, logFile = logFile)
  }

  return(p)

}

.subsetSeqnamesGR <- function(gr = NULL, names = NULL){
  .validInput(input = gr, name = "gr", valid = c("GRanges"))
  .validInput(input = names, name = "names", valid = c("character"))
  gr <- gr[which(as.character(seqnames(gr)) %in% names),]
  seqlevels(gr) <- as.character(unique(seqnames(gr)))
  return(gr)
}

hubsTracks <- function(
  hubs = NULL, 
  region = NULL, 
  title = "HubTrack", 
  pal = NULL,
  baseSize = 5, 
  facetbaseSize = 9,
  featureWidth = 2, 
  borderWidth = 0.1, 
  hideX = FALSE, 
  hideY = FALSE,
  logFile = NULL
  ){

  getArchDF <- function(lp, r = 100){
    angles <- seq(pi, 2*pi,length.out=100)
    rx <- (end(lp)-start(lp))/2
    rscale <- r * (rx/max(rx))
    cx <- start(lp) + rx
    if(is.null(mcols(lp)$value)){
      mcols(lp)$value <- 1
    }
    df <- lapply(seq_along(cx), function(z){
      xz <- rx[z]*cos(angles)+cx[z]
      dfz <- DataFrame(x=xz, y=rscale[z]*sin(angles), id=Rle(paste0("l",z)), value = mcols(lp)$value[z])
    }) %>% Reduce("rbind",.)
    return(df)
  }

  hubs <- SimpleList(hubs = hubs)
 
    valueMin <- min(unlist(lapply(hubs, function(x) min(x$value))))
    valueMax <- max(unlist(lapply(hubs, function(x) max(x$value))))

    hubsO <- lapply(seq_along(hubs), function(x){
       subHubs <- subsetByOverlaps(hubs[[x]], region, ignore.strand = TRUE, type = "within") 
       if(length(subHubs)>0){
         dfx <- getArchDF(subHubs)
         dfx$name <- Rle(paste0(names(hubs)[x]))
         dfx
       }else{
         NULL
       }
    }) %>% Reduce("rbind",.)
    .logThis(hubsO, "hubsO", logFile = logFile)

    testDim <- tryCatch({
      if(is.null(hubsO)){
        FALSE
      }
      if(nrow(hubsO) > 0){
        TRUE
      }else{
        FALSE
      }
    }, error = function(x){
      FALSE
    })

    if(testDim){

      hubsO$facet <- title
      if(is.null(pal)){
        pal <- colorRampPalette(c("white",'purple',"black"))(100)
      }

      p <- ggplot(data = data.frame(hubsO), aes(x = x, y = y, group = id, color = value)) + 
        geom_line(size=.5) +
        facet_grid(name ~ .) +
        ylab("") + 
        coord_cartesian(ylim = c(-100,0)) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        scale_color_gradientn(colors = pal, limits = c(valueMin, valueMax)) +
        guides(color= guide_colorbar(barwidth = 0.75, barheight = 3)) +
        #theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth, legendPosition = "right") +
        theme_void() +
        theme(
        legend.text = element_text(size = baseSize),
        #legendPosition = "right",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
      }
    return (p)
    }

#######################################################
# scATAC Track Methods
#######################################################

.scTracks <- function(
  ArchRProj = NULL,
  region = NULL,
  tileSize = 100,
  minCells = 5,
  maxCells = 100,
  groupBy = "Clusters",
  useGroups = NULL,
  threads = 1,
  baseSize = 7,
  scTileSize = 0.5,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  pal = NULL,
  tstart = NULL,
  verbose = FALSE,
  logFile = NULL
  ){

  .requirePackage("ggplot2", source = "cran")

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #Group Info
  cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(!is.null(minCells)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% names(table(cellGroups)[table(cellGroups) >= minCells]),,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  if(!is.null(useGroups)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% useGroups,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  tabGroups <- table(cellGroups)
  
  if(any(tabGroups > maxCells)){
    cellGroups2 <- getCellColData(ArchRProj, groupBy, drop = FALSE)
    splitGroups <- split(rownames(cellGroups2), cellGroups2[,1])
    useCells <- lapply(seq_along(splitGroups), function(x){
      if(length(splitGroups[[x]]) > maxCells){
        sample(splitGroups[[x]], maxCells)
      }else{
        splitGroups[[x]]
      }
    }) %>% unlist
    ArchRProj@cellColData <- ArchRProj@cellColData[useCells,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
    tabGroups <- table(cellGroups)
  }

  cellsBySample <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj, "Sample", drop = TRUE))
  groupsBySample <- split(cellGroups, getCellColData(ArchRProj, "Sample", drop = TRUE))
  uniqueGroups <- gtools::mixedsort(unique(cellGroups))
  
  #Tile Region
  regionTiles <- seq(trunc(start(region) / tileSize), trunc(end(region) / tileSize) + 1) * tileSize
  ArrowFiles <- getArrowFiles(ArchRProj)
  ArrowFiles <- ArrowFiles[names(cellsBySample)]

  groupMat <- .safelapply(seq_along(ArrowFiles), function(i){
    .logMessage(sprintf("Getting Region From Arrow Files %s of %s", i, length(ArrowFiles)), logFile = logFile)
    tryCatch({
      .regionSCArrows(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
    }, error = function(e){
      errorList <- list(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
      .logError(e, fn = ".groupRegionSCArrows", info = .sampleName(ArrowFiles[i]), errorList = errorList, logFile = logFile)
    })
  }, threads = threads) %>% Reduce("cbind" , .)

  groupDF <- data.frame(Matrix::summary(groupMat))
  groupDF$group <- getCellColData(ArchRProj, groupBy, drop = FALSE)[colnames(groupMat)[groupDF$j], 1]
  groupDF <- lapply(split(groupDF, groupDF$group), function(z){
    nz <- tabGroups[z$group[1]]
    nc <- length(unique(z$j))
    idx <- sort(sample(seq_len(nz), nc))
    idx[1] <- 1
    idx[length(idx)] <- nz
    z$y <- idx[match(z$j, unique(z$j))]
    z
  }) %>% Reduce("rbind", .)
  groupDF$bp <- regionTiles[groupDF$i]
  
  if(is.null(pal)){
    pal <- suppressWarnings(paletteDiscrete(values = names(tabGroups)))
  }

  nn <- paste0(names(tabGroups), ":", tabGroups)
  names(nn) <- names(tabGroups)
  groupDF$group2 <- nn[groupDF$group]
  names(pal) <- nn[names(pal)]

  title <- paste0(as.character(seqnames(region)),":", start(region)-1, "-", end(region), " ", title)
  
  #Re-Order
  groupDF$group2 <- factor(
    paste0(groupDF$group2), 
    levels = gtools::mixedsort(unique(paste0(groupDF$group2)))
  )

  p <- ggplot(groupDF, aes(x=bp, y=y, width = tileSize, fill = group2, color = group2)) + 
      geom_tile(size = scTileSize) + 
      facet_grid(group2 ~ ., scales="free_y") + 
      theme_ArchR() + 
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      ylab("Binarized SC Coverage") + 
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme_ArchR(baseSize = baseSize,
                  baseRectSize = borderWidth,
                  baseLineSize = tickWidth,
                  legendPosition = "right",
                  axisTickCm = 0.1) +
      theme(panel.spacing= unit(0, "lines"),
            axis.title.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.text = element_text(
              size = facetbaseSize, 
              color = "black", 
              margin = margin(0,0,0,0, "cm")),
              strip.text.y = element_text(angle = 0),
            strip.background = element_rect(color="black")) +
      guides(fill = FALSE, colour = FALSE) + ggtitle(title)

    p

}

.regionSCArrows <- function(
  ArrowFile = NULL,
  region = NULL,
  regionTiles = NULL,
  tileSize = NULL,
  cellNames = NULL,
  cellGroups = NULL,
  uniqueGroups = NULL,
  logFile = NULL
  ){
  
  cellFragsRegion <- .getFragsFromArrow(
      ArrowFile = ArrowFile, 
      chr = paste0(seqnames(region)), 
      cellNames = cellNames, 
      out = "GRanges"
    ) %>% subsetByOverlaps(., region, ignore.strand = TRUE)
  
  #Starts
  ts <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ids <- which(ts > 0)
  
  #Ends
  te <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ide <- which(te > 0)
  
  #Match
  matchID <- S4Vectors::match(mcols(cellFragsRegion)$RG, cellNames)
  
  #Sparse Matrix
  mat <- Matrix::sparseMatrix(
    i = c(ts[ids], te[ide]),
    j = as.vector(c(matchID[ids], matchID[ide])),
    x = rep(1,  length(ids) + length(ide)),
    dims = c(length(regionTiles), length(cellNames))
  )
  colnames(mat) <- cellNames
  
  mat@x[mat@x > 1] <- 1

  return(mat)

}



####################################
# Combined Feature Plot
####################################

.combinedFeaturePlot <- function(
  plotList = NULL,
  useMatrix = NULL,
  featureMat = NULL,
  log2Norm = FALSE,
  feature = NULL,
  pal = NULL,
  sizes = NULL,
  baseSize = NULL,
  facetbaseSize = NULL,
  borderWidth = NULL,
  tickWidth = NULL
  ){

  .requirePackage("patchwork", installInfo = "devtools::install_github('thomasp85/patchwork')")

  if(is.null(pal)){
    pal <- paletteDiscrete(values=featureMat$Group, set = "stallion")
  }

  if(log2Norm){
    title <- paste0("Log2 ", useMatrix, " : ", feature)
  }else{
    title <- paste0("Raw ", useMatrix, " : ", feature) 
  }

  featurePlot <- ggGroup(
      x = featureMat$Group,
      y = featureMat[,feature],
      groupOrder = gtools::mixedsort(paste0(unique(featureMat$Group))),
      pal = pal
    ) + 
    facet_wrap(x~., ncol=1,scales="free_y",strip.position="right") +
    guides(fill = FALSE, colour = FALSE) +
    theme_ArchR(baseSize = baseSize,
              baseRectSize = borderWidth,
              baseLineSize = tickWidth,
              legendPosition = "right",
              axisTickCm = 0.1) +
    theme(panel.spacing= unit(0, "lines"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_text(
          size = facetbaseSize, 
          color = "black", 
          margin = margin(0,0,0,0, "cm")),
          strip.text.y = element_text(angle = 0),
        strip.background = element_rect(color="black")) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    ggtitle(title)

  if(any(tolower(names(plotList)) %in% "bulktrack")){

    idx <- which(tolower(names(plotList)) == "bulktrack")
    
    p <- plotList[[idx]] + featurePlot
    
    plotList[idx] <- NULL
    
    for(i in seq_along(plotList)){
      p <- p + plotList[[i]] 
    }
    
    p <- p + plot_layout(
      ncol = 3,
      widths = c(3, 1, 0.2), 
      heights = sizes
    )

  }else{


    idx <- which(tolower(names(plotList)) == "sctrack")
    
    p <- plotList[[idx]] + featurePlot
    
    plotList[idx] <- NULL
    
    for(i in seq_along(plotList)){
      p <- p + plotList[[i]] 
    }
    
    p <- p + plot_layout(
      ncol = 3,
      widths = c(3, 1, 0.2), 
      heights = sizes
    )

  }

  p

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
