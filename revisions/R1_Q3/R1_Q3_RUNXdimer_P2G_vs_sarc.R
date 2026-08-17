#!/usr/bin/env Rscript
# R1_Q3_RUNXdimer_P2G_vs_sarc.R
# Using the Peak2Gene-linked target genes of the RUNX-RUNX (mp2) peaks (no nearest-
# gene artifact): (1) GO:BP enrichment, (2) overlap/enrichment with the sarcomatoid
# scS program (cNMF20). Background universe = all P2G-eligible genes.

suppressPackageStartupMessages({library(ArchR)})
addArchRThreads(4)
PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
SCRN <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scrna"

p2g_genes <- readLines(file.path(OUT,"RUNXdimer_mp2_P2G_target_genes.txt"))
p2g_genes <- unique(p2g_genes[nzchar(p2g_genes)])
cat(sprintf("RUNX-dimer P2G target genes: %d\n", length(p2g_genes)))

## universe = all P2G-eligible genes (geneSet of the Peak2Gene links)
proj <- loadArchRProject(PROJ, showLogo=FALSE)
p2g  <- getPeak2GeneLinks(proj, corCutOff=0.45, returnLoops=FALSE)
universe <- unique(as.character(metadata(p2g)$geneSet$name))
cat(sprintf("P2G gene universe: %d genes\n", length(universe)))

## sarcomatoid scS program (cNMF20)
spectra <- readRDS(file.path(SCRN,"cnmf_genelist_25_nfeat_5000.rds"))
cat("cNMF programs available: ", paste(names(spectra),collapse=", "), "\n")
sarc_full <- spectra[["cNMF20"]]
sets <- list(scS_top20  = head(sarc_full,20),
             sarc_top100= head(sarc_full,100),
             sarc_top200= head(sarc_full,min(200,length(sarc_full))))

cat("\n=== overlap of RUNX-dimer P2G targets with sarcomatoid cNMF20 ===\n")
ovrows <- list()
for(nm in names(sets)){
  s <- intersect(sets[[nm]], universe)            # restrict program to testable universe
  hits <- intersect(p2g_genes, s)
  N <- length(universe); K <- length(s); n <- length(intersect(p2g_genes,universe)); k <- length(hits)
  p <- phyper(k-1, K, N-K, n, lower.tail=FALSE)
  fold <- (k/n)/(K/N)
  cat(sprintf("%-12s: %d/%d P2G targets in set | expected %.2f | fold=%.2f | hyperg p=%.3g\n",
              nm, k, n, n*K/N, fold, p))
  cat("   overlap genes: ", paste(sort(hits),collapse=", "), "\n")
  ovrows[[nm]] <- data.frame(set=nm, set_size=K, p2g_in_universe=n, overlap=k,
                             expected=round(n*K/N,2), fold=round(fold,2), hyperg_p=signif(p,3),
                             overlap_genes=paste(sort(hits),collapse=";"))
}
write.csv(do.call(rbind,ovrows), file.path(OUT,"RUNXdimer_P2G_vs_sarc_overlap.csv"), row.names=FALSE)

## GO:BP on the 152 P2G targets (universe-restricted)
tryCatch({
  suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db)})
  eg  <- bitr(p2g_genes, "SYMBOL","ENTREZID", org.Hs.eg.db)$ENTREZID
  bg  <- bitr(universe,  "SYMBOL","ENTREZID", org.Hs.eg.db)$ENTREZID
  ego <- enrichGO(eg, org.Hs.eg.db, ont="BP", universe=bg, pvalueCutoff=0.1, qvalueCutoff=0.2, readable=TRUE)
  if(!is.null(ego) && nrow(as.data.frame(ego))>0){
    df <- as.data.frame(ego)[,c("ID","Description","GeneRatio","p.adjust","Count")]
    write.csv(df, file.path(OUT,"RUNXdimer_P2G_GO_BP.csv"), row.names=FALSE)
    cat("\n=== GO:BP on RUNX-dimer P2G targets (universe = P2G genes), top 20 ===\n")
    print(head(df,20), row.names=FALSE)
  } else cat("\nGO: no significant terms\n")
}, error=function(e) message("GO skipped: ", conditionMessage(e)))
cat("\nDONE\n")
