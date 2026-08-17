#!/usr/bin/env Rscript
# R1_Q3_RUNXdimer_P2G02.R
# Relaxed co-accessibility: map RUNX-RUNX (mp2) peaks to Peak2Gene targets at
# corr>0.2 (vs 0.45), then GO:BP + overlap with sarcomatoid cNMF20. More targets
# -> more power for the overlap/enrichment tests (at the cost of weaker links).

suppressPackageStartupMessages({library(ArchR)})
addArchRThreads(4)
PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
SCRN <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scrna"
CUT  <- 0.2

proj <- loadArchRProject(PROJ, showLogo=FALSE)
ps <- getPeakSet(proj)
gr_bed <- function(f){ d<-read.table(f,sep="\t",col.names=c("chr","s","e")); GRanges(d$chr,IRanges(d$s+1,d$e)) }
peaks_mp2 <- sort(unique(queryHits(findOverlaps(ps, gr_bed(file.path(OUT,"Average50_mp2_hits.bed"))))))

p2g <- getPeak2GeneLinks(proj, corCutOff=CUT, returnLoops=FALSE)
gs  <- metadata(p2g)$geneSet
universe <- unique(as.character(gs$name))
hit <- p2g$idxATAC %in% peaks_mp2 & p2g$Correlation>CUT
targets <- sort(unique(as.character(gs$name[p2g$idxRNA[hit]])))
cat(sprintf("P2G@%.2f: RUNX-dimer target genes = %d (universe %d)\n", CUT, length(targets), length(universe)))
writeLines(targets, file.path(OUT,"RUNXdimer_mp2_P2G02_target_genes.txt"))

## overlap with sarcomatoid cNMF20
spectra <- readRDS(file.path(SCRN,"cnmf_genelist_25_nfeat_5000.rds"))
sets <- list(scS_top20=head(spectra[["cNMF20"]],20), sarc_top100=head(spectra[["cNMF20"]],100),
             sarc_top200=head(spectra[["cNMF20"]],200))
cat("\n=== overlap RUNX-dimer P2G@0.2 targets vs sarcomatoid cNMF20 ===\n")
for(nm in names(sets)){
  s<-intersect(sets[[nm]],universe); hits<-intersect(targets,s)
  N<-length(universe);K<-length(s);n<-length(intersect(targets,universe));k<-length(hits)
  p<-phyper(k-1,K,N-K,n,lower.tail=FALSE); fold<-(k/n)/(K/N)
  cat(sprintf("%-12s: %d/%d in set | exp %.2f | fold=%.2f | p=%.3g | genes: %s\n",
              nm,k,n,n*K/N,fold,p,paste(sort(hits),collapse=",")))
}

## GO:BP (clusterProfiler, universe = P2G genes)
tryCatch({
  suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db)})
  eg<-bitr(targets,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID
  bg<-bitr(universe,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID
  ego<-enrichGO(eg,org.Hs.eg.db,ont="BP",universe=bg,pvalueCutoff=0.1,qvalueCutoff=0.2,readable=TRUE)
  if(!is.null(ego)&&nrow(as.data.frame(ego))>0){
    df<-as.data.frame(ego)[,c("ID","Description","GeneRatio","p.adjust","Count")]
    write.csv(df,file.path(OUT,"RUNXdimer_P2G02_GO_BP.csv"),row.names=FALSE)
    cat("\n=== GO:BP RUNX-dimer P2G@0.2 targets, top 20 ===\n"); print(head(df,20),row.names=FALSE)
  } else cat("\nGO: no significant terms\n")
}, error=function(e) message("GO skipped: ",conditionMessage(e)))
cat("\nDONE\n")
