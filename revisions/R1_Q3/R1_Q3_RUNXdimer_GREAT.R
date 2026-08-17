#!/usr/bin/env Rscript
# R1_Q3_RUNXdimer_GREAT.R
# GREAT-style enrichment of the RUNX-RUNX (Average_50 mp2) peaks, using the
# regulatory-domain binomial model (rGREAT local) with ALL ATAC peaks as the
# background (ATAC-appropriate; corrects the gene-desert / large-family bias that
# inflated mesenchyme/EMT + olfactory/Hb in the naive nearest-gene test).
# Tests GO:BP and the sarcomatoid cNMF20 program (+ SOX9 module as a control set).

suppressPackageStartupMessages({library(ArchR); library(rGREAT); library(org.Hs.eg.db)})
addArchRThreads(4)
PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
SCRN <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scrna"
TSS  <- "TxDb.Hsapiens.UCSC.hg38.knownGene"

proj <- loadArchRProject(PROJ, showLogo=FALSE)
ps <- getPeakSet(proj); ps <- keepStandardChromosomes(ps, pruning.mode="coarse")
gr_bed <- function(f){ d<-read.table(f,sep="\t",col.names=c("chr","s","e")); GRanges(d$chr,IRanges(d$s+1,d$e)) }
mp2hit <- gr_bed(file.path(OUT,"Average50_mp2_hits.bed"))
mp2 <- ps[unique(queryHits(findOverlaps(ps, mp2hit)))]
cat(sprintf("foreground mp2 peaks: %d ; background all peaks: %d\n", length(mp2), length(ps)))

## ---- GO:BP via GREAT (binomial), background = all peaks ----
res_go <- great(mp2, "GO:BP", TSS, background=ps, cores=4, verbose=FALSE)
tb <- getEnrichmentTable(res_go)
tb <- tb[order(tb$p_value), ]
write.csv(head(tb, 200), file.path(OUT,"RUNXdimer_GREAT_GO_BP.csv"), row.names=FALSE)
cat("\n=== GREAT GO:BP top 20 (binomial over regulatory domains, all-peaks bg) ===\n")
cols <- intersect(c("id","description","genome_fraction","observed_region_hits","fold_enrichment","p_value","p_adjust"), colnames(tb))
print(head(tb[,cols], 20), row.names=FALSE)
cat("\nEMT/mesenchyme terms in table:\n")
print(tb[grepl("mesenchym|epithelial to mesenchymal|mesenchyme", tb$description, ignore.case=TRUE), cols], row.names=FALSE)

## ---- sarcomatoid cNMF20 (+ SOX9 module control) as custom gene sets ----
spectra <- readRDS(file.path(SCRN,"cnmf_genelist_25_nfeat_5000.rds"))
to_entrez <- function(sym){ na.omit(mapIds(org.Hs.eg.db, unique(sym), "ENTREZID","SYMBOL")) }
modgenes <- read.csv(file.path(OUT,"cluster32_top20_genes.csv"))$gene
gs <- list(
  sarc_cNMF20_top50  = to_entrez(head(spectra[["cNMF20"]],50)),
  sarc_cNMF20_top100 = to_entrez(head(spectra[["cNMF20"]],100)),
  sarc_cNMF20_top200 = to_entrez(head(spectra[["cNMF20"]],200)),
  SOX9_module_top20  = to_entrez(modgenes))
res_cust <- great(mp2, gs, TSS, background=ps, cores=4, verbose=FALSE, min_gene_set_size=5)
tc <- getEnrichmentTable(res_cust)
write.csv(tc, file.path(OUT,"RUNXdimer_GREAT_sarc_custom.csv"), row.names=FALSE)
cat("\n=== GREAT enrichment of RUNX-dimer peaks for sarcomatoid / SOX9-module gene sets ===\n")
print(tc[, intersect(c("id","genome_fraction","observed_region_hits","fold_enrichment","p_value","p_adjust"), colnames(tc))], row.names=FALSE)
cat("\nDONE\n")
