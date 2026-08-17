#!/usr/bin/env Rscript
# R1_Q3_RUNX2_SNAI2_regulate_SOX9.R
# Test the hypothesis that the bound drivers RUNX2 / SNAI2 sit upstream of SOX9
# and its module (with SOX9 a downstream readout). Three chromatin-level tests:
#  (1) Motif enrichment in the cis-regulatory elements (peaks within +/-50kb of TSS,
#      plus Peak2Gene-linked peaks) of SOX9 + the 20-gene module: are RUNX2/SNAI2
#      (and AP-1) motifs over-represented there vs all peaks? (Fisher)
#  (2) Per-cell: does RUNX2/SNAI2 chromVAR motif activity predict the SOX9 module
#      score (directional, consistent with upstream regulation)?
#  (3) Occupancy: footprint RUNX2 + SNAI2 in the gated cells (are they actually bound,
#      unlike SOX9?) -> handled by extending the footprint set (separate, optional).

suppressPackageStartupMessages({library(ArchR); library(ggplot2)})
addArchRThreads(4)
PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_module_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
MODF <- file.path(OUT, "cluster32_top20_genes.csv")
proj <- loadArchRProject(PROJ, showLogo = FALSE)

modgenes <- read.csv(MODF, stringsAsFactors=FALSE)$gene

## ---- peak set + motif matches ----
ps <- getPeakSet(proj)
ps$peakID <- paste0(seqnames(ps), ":", start(ps), "-", end(ps))
mt <- readRDS(file.path(PROJ, "Annotations/Motif-Matches-In-Peaks.rds"))
M  <- assay(mt, "matches")                       # peaks x motifs logical
mnames <- colnames(M)
pick <- function(p) grep(p, mnames, ignore.case=TRUE, value=TRUE)[1]
panel <- na.omit(setNames(c(pick("^SOX9_"),pick("^RUNX2_"),pick("^SNAI2_"),pick("^SNAI1_"),
                            pick("^FOSL2_"),pick("^JUN_"),pick("^TEAD1_"),pick("^NFIC_"),pick("^CTCF_")),
                          c("SOX9","RUNX2","SNAI2","SNAI1","FOSL2","JUN","TEAD1","NFIC","CTCF")))
cat("motif panel:\n"); print(panel)

## ---- target cis-elements: peaks within +/-50kb of module-gene TSS ----
genes <- getGenes(proj)
g <- genes[!is.na(genes$symbol) & genes$symbol %in% modgenes]
tss <- resize(g, 1, fix="start")
prox <- unique(queryHits(findOverlaps(ps, resize(tss, 100001, fix="center"))))
# add Peak2Gene-linked peaks for these genes if available
p2g_idx <- tryCatch({
  p2g <- getPeak2GeneLinks(proj, corCutOff=0.45, returnLoops=FALSE)
  gs <- metadata(p2g)$geneSet
  hit_genes <- which(gs$name %in% modgenes)
  unique(p2g$idxATAC[p2g$idxRNA %in% hit_genes & p2g$Correlation>0.45])
}, error=function(e){ message("P2G unavailable: ", conditionMessage(e)); integer(0) })
target <- sort(unique(c(prox, p2g_idx)))
cat(sprintf("module cis-peaks: %d proximal + %d P2G = %d unique  (of %d total peaks)\n",
            length(prox), length(p2g_idx), length(target), nrow(M)))

## ---- Fisher enrichment of each motif in target peaks vs all peaks ----
N <- nrow(M); inT <- seq_len(N) %in% target
res <- lapply(names(panel), function(tf){
  hasm <- as.logical(M[, panel[tf]])
  a<-sum(hasm & inT); b<-sum(hasm & !inT); c<-sum(!hasm & inT); d<-sum(!hasm & !inT)
  ft<-fisher.test(matrix(c(a,b,c,d),2), alternative="greater")
  data.frame(TF=tf, motif=panel[tf], n_target_with_motif=a, pct_target=round(100*a/length(target),1),
             pct_bg=round(100*sum(hasm)/N,1), odds_ratio=round(ft$estimate,2), p=signif(ft$p.value,3))
})
enr <- do.call(rbind, res); rownames(enr)<-NULL
enr <- enr[order(-enr$odds_ratio),]
cat("\n=== motif enrichment in SOX9-module cis-elements ===\n"); print(enr)
write.csv(enr, file.path(OUT,"SOX9module_cis_motif_enrichment.csv"), row.names=FALSE)

## ---- per-cell: does RUNX2/SNAI2 motif activity predict the SOX9 module? ----
mm <- getMatrixFromProject(proj, useMatrix="MotifMatrix"); Z <- assays(mm)[["z"]]
mod <- read.csv(file.path(OUT,"SOX9_program_module_P23_percell.csv"), stringsAsFactors=FALSE)
mcol <- grep("module|score", names(mod), ignore.case=TRUE, value=TRUE)[1]
ccol <- names(mod)[1]
modv <- setNames(mod[[mcol]], mod[[ccol]])
cells <- intersect(colnames(Z), names(modv)); modv <- modv[cells]
zrow <- function(tf){ m<-grep(paste0("^",tf,"_[0-9]"),rownames(Z),value=TRUE)[1]; if(is.na(m)) NULL else as.numeric(Z[m,cells]) }
cordf <- do.call(rbind, lapply(c("SOX9","RUNX2","SNAI2","SNAI1","FOSL2","JUN"), function(tf){
  z<-zrow(tf); if(is.null(z)) return(NULL)
  data.frame(TF=tf, pearson=round(cor(z,modv,use="complete.obs"),3),
             spearman=round(cor(z,modv,method="spearman",use="complete.obs"),3))
}))
cat("\n=== per-cell: TF motif deviation vs SOX9 module score ===\n"); print(cordf)
write.csv(cordf, file.path(OUT,"TFdeviation_vs_SOX9module_percell.csv"), row.names=FALSE)

## ---- plot: enrichment bar ----
enr$TF <- factor(enr$TF, levels=enr$TF)
p <- ggplot(enr, aes(TF, odds_ratio, fill=-log10(p+1e-300))) + geom_col() +
  geom_hline(yintercept=1, linetype=2) +
  geom_text(aes(label=sprintf("%.1f%%\np=%.0e", pct_target, p)), vjust=-0.2, size=2.6) +
  scale_fill_viridis_c(name="-log10 p") +
  labs(title="TF motif enrichment in SOX9-module cis-regulatory elements",
       subtitle=sprintf("peaks +/-50kb of %d module genes + P2G links; vs all peaks (Fisher)", length(unique(g$symbol))),
       y="odds ratio (enrichment)", x=NULL) + theme_bw(base_size=11)
ggsave(file.path(OUT,"SOX9module_cis_motif_enrichment.pdf"), p, width=7, height=4.5)
cat("\nWrote enrichment + per-cell correlation + plot\n")
