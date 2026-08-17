#!/usr/bin/env Rscript
# R1_Q3_RUNXdimer_target_genes.R
# What genes do the RUNX-RUNX dimer (Average_50 mp2) peaks map to?
# Annotate the mp2 peaks (in scatac_ArchR_SOX9_P23) to (a) nearest gene (ArchR
# peakSet annotation) and (b) Peak2Gene-linked target genes (corr>0.45), output
# the gene lists, and attempt a GO enrichment of the linked targets.

suppressPackageStartupMessages({library(ArchR)})
addArchRThreads(4)
PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
proj <- loadArchRProject(PROJ, showLogo=FALSE)
ps <- getPeakSet(proj)

gr_bed <- function(f){ d<-read.table(f,sep="\t",col.names=c("chr","s","e")); GRanges(d$chr,IRanges(d$s+1,d$e)) }
mp2gr <- gr_bed(file.path(OUT,"Average50_mp2_hits.bed"))
peaks_mp2 <- sort(unique(queryHits(findOverlaps(ps, mp2gr))))
mp2 <- ps[peaks_mp2]
cat(sprintf("mp2 (RUNX-RUNX): %d hits -> %d unique ArchR peaks\n", length(mp2gr), length(mp2)))
cat("peakSet annotation columns: ", paste(colnames(mcols(mp2)),collapse=", "), "\n")

## (a) nearest-gene annotation
ng_col  <- intersect(c("nearestGene"), colnames(mcols(mp2)))
dist_col<- intersect(c("distToTSS","distToGeneStart"), colnames(mcols(mp2)))
pt_col  <- intersect(c("peakType"), colnames(mcols(mp2)))
ann <- data.frame(
  peak = paste0(seqnames(mp2),":",start(mp2),"-",end(mp2)),
  nearestGene = if(length(ng_col)) mcols(mp2)[[ng_col[1]]] else NA,
  distToTSS   = if(length(dist_col)) mcols(mp2)[[dist_col[1]]] else NA,
  peakType    = if(length(pt_col)) mcols(mp2)[[pt_col[1]]] else NA,
  stringsAsFactors=FALSE)
write.csv(ann, file.path(OUT,"RUNXdimer_mp2_peaks_annotated.csv"), row.names=FALSE)

cat("\n=== peakType distribution of mp2 peaks ===\n"); print(table(ann$peakType, useNA="ifany"))
near_genes <- sort(table(ann$nearestGene[!is.na(ann$nearestGene)]), decreasing=TRUE)
cat(sprintf("\nunique nearest genes: %d\n", length(near_genes)))
cat("=== top 30 nearest genes (by #mp2 peaks) ===\n"); print(head(near_genes,30))
writeLines(names(near_genes), file.path(OUT,"RUNXdimer_mp2_nearest_genes.txt"))

## (b) Peak2Gene-linked target genes
p2g_genes <- character(0)
tryCatch({
  p2g <- getPeak2GeneLinks(proj, corCutOff=0.45, returnLoops=FALSE)
  gs <- metadata(p2g)$geneSet
  hit <- p2g$idxATAC %in% peaks_mp2 & p2g$Correlation>0.45
  p2g_genes <- sort(unique(as.character(gs$name[p2g$idxRNA[hit]])))
  cat(sprintf("\nPeak2Gene-linked target genes (corr>0.45) for mp2 peaks: %d\n", length(p2g_genes)))
  print(head(p2g_genes, 50))
  writeLines(p2g_genes, file.path(OUT,"RUNXdimer_mp2_P2G_target_genes.txt"))
}, error=function(e) message("P2G step failed: ", conditionMessage(e)))

## (c) GO enrichment of the gene set (nearest genes; fall back gracefully)
geneset <- unique(c(names(near_genes), p2g_genes))
cat(sprintf("\nCombined unique gene set for enrichment: %d genes\n", length(geneset)))
tryCatch({
  suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db)})
  eg <- bitr(geneset, "SYMBOL","ENTREZID", org.Hs.eg.db)$ENTREZID
  ego <- enrichGO(eg, org.Hs.eg.db, ont="BP", pvalueCutoff=0.1, qvalueCutoff=0.2, readable=TRUE)
  if(!is.null(ego) && nrow(as.data.frame(ego))>0){
    df <- as.data.frame(ego)[,c("ID","Description","GeneRatio","p.adjust","Count")]
    write.csv(df, file.path(OUT,"RUNXdimer_mp2_GO_BP.csv"), row.names=FALSE)
    cat("\n=== top 20 GO:BP terms (RUNX-dimer target genes) ===\n"); print(head(df,20), row.names=FALSE)
  } else cat("\nGO enrichment: no significant terms\n")
}, error=function(e) message("GO enrichment skipped: ", conditionMessage(e)))
cat("\nDONE — wrote annotated peaks + nearest/P2G gene lists (+GO if available)\n")
