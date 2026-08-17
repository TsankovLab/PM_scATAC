suppressPackageStartupMessages({library(ArchR); library(chromVAR); library(Matrix); library(SummarizedExperiment); library(data.table); library(BSgenome.Hsapiens.UCSC.hg38)})
addArchRThreads(4)
PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
WD   <- "/sc/arion/scratch/giottb01/meso_hdma_chromvar"
proj <- loadArchRProject(PROJ, showLogo=FALSE)
ps <- getPeakSet(proj); seqlevelsStyle(ps) <- "UCSC"; N <- length(ps)

## peaks x patterns matches from FI-NeMo counts-head hits
h <- fread(file.path(WD,"hdma_counts_hits.tsv"), header=FALSE, col.names=c("chr","start","end","pat"))
gr <- GRanges(h$chr, IRanges(h$start+1, h$end))
ov <- findOverlaps(gr, ps); pats <- sort(unique(h$pat))
patf <- factor(h$pat[queryHits(ov)], levels=pats)
matches <- sparseMatrix(i=subjectHits(ov), j=as.integer(patf), x=1L, dims=c(N,length(pats)), dimnames=list(NULL,pats))>0
cat(sprintf("patterns=%d peaks=%d\n", length(pats), N))

## counts in peakSet order
pm <- getMatrixFromProject(proj, useMatrix="PeakMatrix"); pmGR <- rowRanges(pm)
ov2 <- findOverlaps(ps, pmGR, type="equal")
counts <- Matrix(0, N, ncol(pm), sparse=TRUE); counts[queryHits(ov2),] <- assay(pm)[subjectHits(ov2),]
colnames(counts) <- colnames(pm)

## drop empty peaks (chromVAR requirement); keep matches aligned
keep <- Matrix::rowSums(counts) > 0
cat(sprintf("peaks with >=1 fragment: %d of %d\n", sum(keep), N))
se <- SummarizedExperiment(assays=list(counts=counts[keep,]), rowRanges=ps[keep])
matches <- matches[keep,,drop=FALSE]
npk <- Matrix::colSums(matches); matches <- matches[, npk>=10, drop=FALSE]
cat(sprintf("patterns kept (>=10 peaks): %d\n", ncol(matches)))

## chromVAR own GC background
se <- addGCBias(se, genome=BSgenome.Hsapiens.UCSC.hg38)
bg <- getBackgroundPeaks(se, niterations=50)
dev <- computeDeviations(object=se, annotations=matches, background_peaks=bg)
z <- deviationScores(dev)
saveRDS(dev, file.path(WD,"hdma_chromvar_dev.rds"))
write.csv(t(z), file.path(OUT,"HDMA_pattern_chromvar_z_percell.csv"))

## summarize by SOX9 status
g <- read.csv(file.path(OUT,"SOX9_module_gated_cells.csv"), stringsAsFactors=FALSE)
clu <- g$SOX9_cluster[match(colnames(z), g$cell)]; hi<-clu=="SOX9_high_P23"; lo<-clu=="SOX9_low_P23"
summ <- data.frame(pattern=rownames(z), mean_z_high=round(rowMeans(z[,hi,drop=FALSE],na.rm=TRUE),3),
                   mean_z_low=round(rowMeans(z[,lo,drop=FALSE],na.rm=TRUE),3))
summ$delta <- round(summ$mean_z_high-summ$mean_z_low,3)
summ$wilcox_p <- apply(z,1,function(r) tryCatch(wilcox.test(r[hi],r[lo])$p.value, error=function(e) NA))
summ <- summ[order(-summ$delta),]
write.csv(summ, file.path(OUT,"HDMA_pattern_chromvar_by_SOX9.csv"), row.names=FALSE)
cat("\n=== most SOX9-HIGH-selective patterns ===\n"); print(head(summ,15), row.names=FALSE)
cat("\n=== most SOX9-LOW-selective ===\n"); print(tail(summ,10), row.names=FALSE)
cat("\n=== Average_50 mp2 & Average_65 mp3 ===\n")
print(summ[grepl("Average_50__merged_pattern_2$|Average_65__merged_pattern_3$", summ$pattern),], row.names=FALSE)
cat("\nDONE\n")
