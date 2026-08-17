#!/usr/bin/env Rscript
# R1_Q3_allTF_cis_enrichment_GCmatched.R
# GC/accessibility-MATCHED motif enrichment in SOX9 cis-elements, fixing the
# GC/promoter artifact of the all-peaks background. Background = the chromVAR
# Background-Peaks (each peak's GC+accessibility-matched peers), so promoter/GC
# composition is controlled. Hypergeometric/Fisher per motif vs matched bg, BH.
# Target sets: SOX9 locus only, and the full SOX9 20-gene module.

suppressPackageStartupMessages({library(ArchR); library(Matrix); library(ggplot2)})
addArchRThreads(4)
PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_module_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
MODF <- file.path(OUT, "cluster32_top20_genes.csv")
proj <- loadArchRProject(PROJ, showLogo = FALSE)
modgenes <- read.csv(MODF, stringsAsFactors=FALSE)$gene

ps <- getPeakSet(proj)
mt <- readRDS(file.path(PROJ, "Annotations/Motif-Matches-In-Peaks.rds"))
M  <- assay(mt, "matches")                       # peaks x motifs (sparse logical)
N  <- nrow(M)
genes <- getGenes(proj)
message("Loading Background-Peaks (GC+accessibility matched) ...")
bgd_se <- readRDS(file.path(PROJ, "Background-Peaks.rds"))  # SummarizedExperiment
bgd <- as.matrix(SummarizedExperiment::assay(bgd_se, "bgdPeaks"))  # nPeaks x k matched indices
stopifnot(nrow(bgd) == N)
cat("Background-Peaks dim:", paste(dim(bgd),collapse=" x "), "\n")

cis_peaks <- function(syms, flank=100000) {
  g <- genes[!is.na(genes$symbol) & genes$symbol %in% syms]
  tss <- resize(g, 1, fix="start")
  prox <- unique(queryHits(findOverlaps(ps, resize(tss, 2*flank+1, fix="center"))))
  p2g <- tryCatch({
    pl <- getPeak2GeneLinks(proj, corCutOff=0.45, returnLoops=FALSE)
    gs <- metadata(pl)$geneSet
    unique(pl$idxATAC[pl$idxRNA %in% which(gs$name %in% syms) & pl$Correlation>0.45])
  }, error=function(e) integer(0))
  sort(unique(c(prox, p2g)))
}

enrich_matched <- function(target, label) {
  target <- target[target>=1 & target<=N]
  # GC+accessibility-matched background = matched peers of the target peaks
  bg <- unique(as.vector(bgd[target, ]))
  bg <- setdiff(bg, target)                       # don't let target be its own bg
  nT <- length(target); nB <- length(bg)
  aT <- Matrix::colSums(M[target, , drop=FALSE])  # motif in target
  aB <- Matrix::colSums(M[bg, , drop=FALSE])       # motif in matched bg
  OR <- (aT/(nT-aT)) / (aB/(nB-aB))
  # Fisher one-sided (greater) per motif via hypergeometric on the pooled target+bg
  white <- aT + aB; draws <- nT
  p <- phyper(aT-1, white, (nT+nB)-white, draws, lower.tail=FALSE)
  q <- p.adjust(p, "BH")
  df <- data.frame(motif=colnames(M), n_target=as.integer(aT),
                   pct_target=round(100*aT/nT,2), pct_matchedBg=round(100*aB/nB,2),
                   odds_ratio=round(OR,2), p=signif(p,3), qval=signif(q,3), target=label)
  df <- df[order(-df$odds_ratio), ]
  write.csv(df, file.path(OUT, sprintf("allTF_cis_GCmatched_%s.csv", label)), row.names=FALSE)
  cat(sprintf("\n=== %s : %d target / %d matched-bg peaks | top 20 (a>=5) ===\n", label, nT, nB))
  print(head(df[df$n_target>=5, c("motif","n_target","pct_target","pct_matchedBg","odds_ratio","qval")], 20))
  df
}

mod <- enrich_matched(cis_peaks(modgenes), "SOX9module")
loc <- enrich_matched(cis_peaks("SOX9"),   "SOX9locus")

## candidate ranks
cand <- c("SOX9_756","RUNX2_732","SNAI2_161","SNAI1_199","FOSL2_105","JUN_143","TEAD1_796","SMAD5_866","CTCF_177")
rankof <- function(df){ d<-df[df$n_target>=5,]; d$rank<-seq_len(nrow(d)); d[d$motif %in% cand, c("rank","motif","odds_ratio","qval","pct_target","pct_matchedBg")] }
cat("\n=== candidate ranks (GC-matched, SOX9module) ===\n"); print(rankof(mod))
cat("\n=== candidate ranks (GC-matched, SOX9locus) ===\n");  print(rankof(loc))

## plot top 15 module (GC-matched)
top <- head(mod[mod$n_target>=5 & mod$qval<0.1, ], 15)
if (nrow(top)>=1){
  top$motif <- factor(top$motif, levels=rev(top$motif))
  p <- ggplot(top, aes(odds_ratio, motif, fill=-log10(qval+1e-300))) + geom_col() +
    geom_vline(xintercept=1, linetype=2) + scale_fill_viridis_c(name="-log10 q") +
    labs(title="GC-matched motif enrichment: SOX9-module cis-elements",
         subtitle="background = chromVAR GC+accessibility-matched peaks (q<0.1)",
         x="odds ratio", y=NULL) + theme_bw(base_size=10)
  ggsave(file.path(OUT,"allTF_cis_GCmatched_SOX9module_top.pdf"), p, width=7, height=5)
}
cat("\nWrote allTF_cis_GCmatched_{SOX9module,SOX9locus}.csv (+plot)\n")
