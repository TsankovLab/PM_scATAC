#!/usr/bin/env Rscript
# R1_Q3_Average65_mp3_core.R
# Same battery as Average_50 mp2, for Average_65 . mp3 (merged_pattern_3):
#  (1) per-cell accessibility activity of mp3 peaks, by SOX9 status
#  (2) enrichment of mp3 peaks at the SOX9 module + SOX9-locus cis-regulome
#  (3) gene mapping: nearest gene + Peak2Gene targets (corr>0.45 and >0.2)
#  (4) overlap of P2G targets with sarcomatoid cNMF20; GO:BP on P2G@0.2 targets

suppressPackageStartupMessages({library(ArchR); library(Matrix)})
addArchRThreads(4)
PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
SCRN <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scrna"
BED  <- file.path(OUT,"Average65_mp3_hits.bed")
proj <- loadArchRProject(PROJ, showLogo=FALSE)
ps <- getPeakSet(proj); N <- length(ps); genes <- getGenes(proj)
gr_bed <- function(f){ d<-read.table(f,sep="\t",col.names=c("chr","s","e")); GRanges(d$chr,IRanges(d$s+1,d$e)) }
peaks <- sort(unique(queryHits(findOverlaps(ps, gr_bed(BED)))))
mp <- ps[peaks]
cat(sprintf("Average_65 mp3: %d hits -> %d ArchR peaks\n", length(gr_bed(BED)), length(peaks)))

## (1) per-cell activity by SOX9 status
pm <- getMatrixFromProject(proj, useMatrix="PeakMatrix")
pmGR <- rowRanges(pm); ov <- findOverlaps(pmGR, ps, type="equal")
ps2row <- integer(N); ps2row[subjectHits(ov)] <- queryHits(ov)
rows <- ps2row[peaks]; rows <- rows[rows>0]
M <- assay(pm); depth <- Matrix::colSums(M)
act <- Matrix::colSums(M[rows,,drop=FALSE])/length(rows)/depth*1e4
g <- read.csv(file.path(OUT,"SOX9_module_gated_cells.csv"), stringsAsFactors=FALSE)
clust <- g$SOX9_cluster[match(colnames(M), g$cell)]; gate <- g$gate[match(colnames(M), g$cell)]
df <- data.frame(cell=colnames(M), clust=clust, gate=gate, act_mp3_meanpp=act)
write.csv(df, file.path(OUT,"Average65_mp3_percell_activity.csv"), row.names=FALSE)
cat("\n=== per-cell mp3 activity (mean-per-peak x1e4) by SOX9 cluster ===\n")
print(aggregate(act_mp3_meanpp~clust, df, function(x) round(mean(x),4)))
cat("by gate:\n"); print(aggregate(act_mp3_meanpp~gate, df, function(x) round(mean(x),4)))

## (2) cis-regulome enrichment
p2g45 <- tryCatch(getPeak2GeneLinks(proj, corCutOff=0.45, returnLoops=FALSE), error=function(e) NULL)
cis_of <- function(syms, flank=50000, p2g=p2g45, cut=0.45){
  gg<-genes[!is.na(genes$symbol)&genes$symbol%in%syms]
  prox<-unique(queryHits(findOverlaps(ps, resize(resize(gg,1,fix="start"),2*flank+1,fix="center"))))
  pl<-integer(0); if(!is.null(p2g)){gs<-metadata(p2g)$geneSet; pl<-unique(p2g$idxATAC[p2g$idxRNA%in%which(gs$name%in%syms)&p2g$Correlation>cut])}
  sort(unique(c(prox,pl)))
}
modg <- read.csv(file.path(OUT,"cluster32_top20_genes.csv"))$gene
enr <- function(cis,label){ inset<-seq_len(N)%in%peaks; incis<-seq_len(N)%in%cis
  a<-sum(inset&incis);b<-sum(inset&!incis);c<-sum(!inset&incis);d<-sum(!inset&!incis)
  ft<-fisher.test(matrix(c(a,b,c,d),2),alternative="greater")
  data.frame(cis=label,n_peaks=length(peaks),n_in_cis=a,pct=round(100*a/length(peaks),2),pct_bg=round(100*sum(incis)/N,2),OR=round(ft$estimate,2),p=signif(ft$p.value,3)) }
cisr <- rbind(enr(cis_of(modg),"SOX9_module_cis"), enr(cis_of("SOX9"),"SOX9_locus_cis"))
cat("\n=== mp3 enrichment at SOX9 cis-regulome ===\n"); print(cisr,row.names=FALSE)
write.csv(cisr, file.path(OUT,"Average65_mp3_cis_enrichment.csv"), row.names=FALSE)

## (3) gene mapping
ann <- data.frame(peak=paste0(seqnames(mp),":",start(mp),"-",end(mp)),
                  nearestGene=mp$nearestGene, distToTSS=mp$distToTSS, peakType=mp$peakType)
write.csv(ann, file.path(OUT,"Average65_mp3_peaks_annotated.csv"), row.names=FALSE)
cat("\npeakType:\n"); print(table(ann$peakType))
targ <- function(cut){ p2g<-getPeak2GeneLinks(proj,corCutOff=cut,returnLoops=FALSE); gs<-metadata(p2g)$geneSet
  list(t=sort(unique(as.character(gs$name[p2g$idxRNA[p2g$idxATAC%in%peaks & p2g$Correlation>cut]]))), u=unique(as.character(gs$name))) }
t45<-targ(0.45); t02<-targ(0.2)
cat(sprintf("\nP2G targets: corr>0.45 = %d ; corr>0.2 = %d\n", length(t45$t), length(t02$t)))
writeLines(t45$t, file.path(OUT,"Average65_mp3_P2G045_genes.txt"))
writeLines(t02$t, file.path(OUT,"Average65_mp3_P2G02_genes.txt"))

## (4) sarc overlap (P2G@0.2) + GO
spectra <- readRDS(file.path(SCRN,"cnmf_genelist_25_nfeat_5000.rds")); uni<-t02$u
sets <- list(scS_top20=head(spectra[["cNMF20"]],20),sarc_top100=head(spectra[["cNMF20"]],100),sarc_top200=head(spectra[["cNMF20"]],200))
cat("\n=== Average_65 mp3 P2G@0.2 targets vs sarcomatoid cNMF20 ===\n")
for(nm in names(sets)){ s<-intersect(sets[[nm]],uni);h<-intersect(t02$t,s)
  K<-length(s);n<-length(intersect(t02$t,uni));k<-length(h);p<-phyper(k-1,K,length(uni)-K,n,lower.tail=FALSE)
  cat(sprintf("%-12s: %d/%d | exp %.2f | fold=%.2f | p=%.3g | %s\n",nm,k,n,n*K/length(uni),(k/n)/(K/length(uni)),p,paste(sort(h),collapse=","))) }
tryCatch({ suppressPackageStartupMessages({library(clusterProfiler);library(org.Hs.eg.db)})
  eg<-bitr(t02$t,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID; bg<-bitr(uni,"SYMBOL","ENTREZID",org.Hs.eg.db)$ENTREZID
  ego<-enrichGO(eg,org.Hs.eg.db,ont="BP",universe=bg,pvalueCutoff=0.1,qvalueCutoff=0.2,readable=TRUE)
  if(!is.null(ego)&&nrow(as.data.frame(ego))>0){df<-as.data.frame(ego)[,c("ID","Description","GeneRatio","p.adjust","Count")]
    write.csv(df,file.path(OUT,"Average65_mp3_P2G02_GO_BP.csv"),row.names=FALSE)
    cat("\n=== GO:BP Average_65 mp3 P2G@0.2 targets, top 20 ===\n"); print(head(df,20),row.names=FALSE)}
}, error=function(e) message("GO skipped: ",conditionMessage(e)))
cat("\nDONE\n")
