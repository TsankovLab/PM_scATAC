#!/usr/bin/env Rscript
# R1_Q3_RUNXRUNX_at_SOX9module.R
# Are the RUNX-RUNX (Average_50 mp2) peaks specifically the cis-elements of the
# SOX9 program? Test whether mp2 peaks concentrate at the SOX9 + 20-module-gene
# cis-regulome (peaks +/-50kb of TSS + Peak2Gene links), using single-RUNX (mp0)
# and all peaks as baselines. mp2-vs-mp0 controls for "RUNX-ness" -> isolates the
# clustered/dimeric effect. Also per-module-gene mp2 vs mp0 peak counts.

suppressPackageStartupMessages({library(ArchR); library(ggplot2)})
addArchRThreads(4)
PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
modgenes <- read.csv(file.path(OUT,"cluster32_top20_genes.csv"))$gene
proj <- loadArchRProject(PROJ, showLogo=FALSE)
ps <- getPeakSet(proj); N <- length(ps)
genes <- getGenes(proj)

gr_bed <- function(f){ d<-read.table(f,sep="\t",col.names=c("chr","s","e")); GRanges(d$chr,IRanges(d$s+1,d$e)) }
peaks_mp0 <- sort(unique(queryHits(findOverlaps(ps, gr_bed(file.path(OUT,"Average50_mp0_hits.bed"))))))
peaks_mp2 <- sort(unique(queryHits(findOverlaps(ps, gr_bed(file.path(OUT,"Average50_mp2_hits.bed"))))))

## cis-regulome of a gene set: peaks +/-50kb of TSS  +  P2G-linked peaks
p2g <- tryCatch(getPeak2GeneLinks(proj, corCutOff=0.45, returnLoops=FALSE), error=function(e) NULL)
cis_of <- function(syms, flank=50000){
  g <- genes[!is.na(genes$symbol) & genes$symbol %in% syms]
  prox <- unique(queryHits(findOverlaps(ps, resize(resize(g,1,fix="start"), 2*flank+1, fix="center"))))
  pl <- integer(0)
  if(!is.null(p2g)){ gs<-metadata(p2g)$geneSet; pl<-unique(p2g$idxATAC[p2g$idxRNA %in% which(gs$name %in% syms) & p2g$Correlation>0.45]) }
  sort(unique(c(prox,pl)))
}
mod_cis  <- cis_of(modgenes)
sox9_cis <- cis_of("SOX9")

enr <- function(setidx, cisidx, label, setlabel){
  inset <- seq_len(N) %in% setidx; incis <- seq_len(N) %in% cisidx
  a<-sum(inset&incis); b<-sum(inset&!incis); c<-sum(!inset&incis); d<-sum(!inset&!incis)
  ft<-fisher.test(matrix(c(a,b,c,d),2), alternative="greater")
  data.frame(peak_set=setlabel, cis=label, n_peaks=length(setidx), n_in_cis=a,
             pct_in_cis=round(100*a/length(setidx),2), pct_bg=round(100*sum(incis)/N,2),
             odds_ratio=round(ft$estimate,2), p=signif(ft$p.value,3))
}
allpk <- seq_len(N)
res <- rbind(
  enr(allpk,     mod_cis,"SOX9_module_cis","all_peaks"),
  enr(peaks_mp0, mod_cis,"SOX9_module_cis","mp0_RUNX"),
  enr(peaks_mp2, mod_cis,"SOX9_module_cis","mp2_RUNX_RUNX"),
  enr(allpk,     sox9_cis,"SOX9_locus_cis","all_peaks"),
  enr(peaks_mp0, sox9_cis,"SOX9_locus_cis","mp0_RUNX"),
  enr(peaks_mp2, sox9_cis,"SOX9_locus_cis","mp2_RUNX_RUNX"))
cat("=== enrichment of peak sets in SOX9 cis-regulomes ===\n"); print(res, row.names=FALSE)
write.csv(res, file.path(OUT,"RUNXRUNX_SOX9module_enrichment.csv"), row.names=FALSE)

## direct mp2 vs mp0-only: is RUNX-RUNX more module-cis than single RUNX? (controls RUNX-ness)
mp0only <- setdiff(peaks_mp0, peaks_mp2)
a<-sum(peaks_mp2 %in% mod_cis); b<-length(peaks_mp2)-a
c<-sum(mp0only %in% mod_cis);   d<-length(mp0only)-c
ft<-fisher.test(matrix(c(a,b,c,d),2))
cat(sprintf("\nmp2(RUNX-RUNX) vs mp0-only(single RUNX) in module-cis: %.1f%% vs %.1f%%  OR=%.2f p=%.3g\n",
            100*a/length(peaks_mp2), 100*c/length(mp0only), ft$estimate, ft$p.value))

## per module gene: how many mp2 / mp0 peaks in its cis-regulome
perg <- do.call(rbind, lapply(c("SOX9",modgenes), function(g){
  ci <- cis_of(g)
  data.frame(gene=g, n_cis_peaks=length(ci),
             mp2_peaks=sum(peaks_mp2 %in% ci), mp0_peaks=sum(peaks_mp0 %in% ci))
}))
perg <- unique(perg); perg <- perg[order(-perg$mp2_peaks),]
cat("\n=== per gene: mp2 (RUNX-RUNX) and mp0 (RUNX) peaks in its cis-regulome ===\n"); print(perg, row.names=FALSE)
write.csv(perg, file.path(OUT,"RUNXRUNX_per_module_gene.csv"), row.names=FALSE)

## plot
res$peak_set <- factor(res$peak_set, levels=c("all_peaks","mp0_RUNX","mp2_RUNX_RUNX"))
p <- ggplot(res, aes(peak_set, pct_in_cis, fill=peak_set)) + geom_col() +
  geom_text(aes(label=sprintf("%.1f%%\nOR=%.1f",pct_in_cis,odds_ratio)), vjust=-0.1, size=2.8) +
  facet_wrap(~cis, scales="free_y") +
  scale_fill_manual(values=c("grey70","#4C78A8","#E45756"), guide="none") +
  labs(title="RUNX-RUNX (mp2) vs single-RUNX (mp0) peaks at the SOX9 cis-regulome",
       subtitle="% of each peak set falling in the cis-elements (+/-50kb TSS + P2G); OR vs all peaks",
       x=NULL, y="% peaks in cis-regulome") + theme_bw(base_size=11)
ggsave(file.path(OUT,"RUNXRUNX_SOX9module_enrichment.pdf"), p, width=8, height=4.5)
cat("\nWrote enrichment csv + per-gene csv + plot\n")
