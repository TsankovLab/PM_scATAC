#!/usr/bin/env Rscript
# R1_Q3_allTF_cis_enrichment_SOX9.R
# Unbiased motif enrichment across ALL motifs (~870) in the cis-regulatory
# elements of (A) the SOX9 gene locus only, and (B) the full SOX9 20-gene module.
# Hypergeometric (one-sided Fisher) per motif, BH-adjusted. Ranks which TFs'
# motifs populate SOX9's regulatory DNA, with no pre-selected panel.

suppressPackageStartupMessages({library(ArchR); library(Matrix); library(ggplot2)})
addArchRThreads(4)
PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_module_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
MODF <- file.path(OUT, "cluster32_top20_genes.csv")
proj <- loadArchRProject(PROJ, showLogo = FALSE)
modgenes <- read.csv(MODF, stringsAsFactors=FALSE)$gene

ps <- getPeakSet(proj)
mt <- readRDS(file.path(PROJ, "Annotations/Motif-Matches-In-Peaks.rds"))
M  <- assay(mt, "matches")              # peaks x motifs (sparse logical)
N  <- nrow(M); m1 <- Matrix::colSums(M)  # motif present, all peaks
genes <- getGenes(proj)

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

enrich_all <- function(target, label) {
  nT <- length(target); inT <- seq_len(N) %in% target
  a  <- Matrix::colSums(M[inT, , drop=FALSE])        # motif present & in target
  OR <- (a * (N - m1 - (nT - a))) / ((m1 - a) * (nT - a))
  p  <- phyper(a - 1, m1, N - m1, nT, lower.tail=FALSE)
  q  <- p.adjust(p, "BH")
  df <- data.frame(motif=names(a), n_target=as.integer(a),
                   pct_target=round(100*a/nT,2), pct_bg=round(100*m1/N,2),
                   odds_ratio=round(OR,2), p=signif(p,3), qval=signif(q,3))
  df <- df[order(-df$odds_ratio), ]
  df$target <- label
  write.csv(df, file.path(OUT, sprintf("allTF_cis_enrichment_%s.csv", label)), row.names=FALSE)
  cat(sprintf("\n=== %s : %d cis-peaks | top 20 enriched motifs (a>=5) ===\n", label, nT))
  print(head(df[df$n_target>=5, c("motif","n_target","pct_target","pct_bg","odds_ratio","qval")], 20))
  # where does SOX9 rank?
  rk <- which(df$motif[df$n_target>=5] == grep("^SOX9_", df$motif, value=TRUE)[1])
  cat(sprintf("SOX9 motif rank among enriched (a>=5): %s of %d\n",
              ifelse(length(rk), rk, "NA"), sum(df$n_target>=5)))
  df
}

sox9_only <- enrich_all(cis_peaks("SOX9"),     "SOX9locus")
module    <- enrich_all(cis_peaks(modgenes),   "SOX9module")

## ---- plot: top 15 enriched at the SOX9 locus ----
top <- head(sox9_only[sox9_only$n_target>=5, ], 15)
top$motif <- factor(top$motif, levels=rev(top$motif))
p <- ggplot(top, aes(odds_ratio, motif, fill=-log10(qval+1e-300))) + geom_col() +
  geom_vline(xintercept=1, linetype=2) +
  scale_fill_viridis_c(name="-log10 q") +
  labs(title="Top motifs enriched at the SOX9-locus cis-elements (all-TF scan)",
       subtitle=sprintf("peaks +/-100kb of SOX9 TSS + P2G; vs all peaks; %d cis-peaks", length(cis_peaks("SOX9"))),
       x="odds ratio (enrichment)", y=NULL) + theme_bw(base_size=10)
ggsave(file.path(OUT,"allTF_cis_enrichment_SOX9locus_top.pdf"), p, width=7, height=5)
cat("\nWrote allTF_cis_enrichment_{SOX9locus,SOX9module}.csv + plot\n")
