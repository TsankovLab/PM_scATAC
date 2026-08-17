###############################################################
# R2_Q14 : WITHIN-SAMPLE TF activity vs BAP1 (removes the histology confound).
#
# Rationale: bulk showed the BAP1<->RUNX link is largely histology (sarc score).
# To test BAP1 *independent of histology*, compare cells WITHIN a single tumor:
# BAP1-high vs BAP1-low cells. The reviewer suggested P8; inferCNV shows P8 is
# CLONALLY BAP1-deleted (no retained subclone), whereas P11 and P5 are genuinely
# subclonal for chr3p/BAP1. So we run the analysis systematically across every
# tumor that has a BAP1-accessible subpopulation, feature the subclonal ones, and
# report P8 with its (honest) clonality caveat.
#
# TF activity = chromVAR deviations (scaled), differential BAP1-high vs BAP1-low
# via presto::wilcoxauc across ALL motifs (not just RUNX).
###############################################################
suppressMessages({ library(ArchR); library(Matrix); library(presto)
                   library(ggplot2); library(ggrepel) })
addArchRThreads(threads = 1); addArchRGenome("hg38")

BASE   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment"
GITROOT<- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo"
OUTDIR <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q14"
setwd(OUTDIR); dir.create("Plots", showWarnings = FALSE)
CHR3_CENT <- 90e6
runx <- c("RUNX1","RUNX2","RUNX3","CBFB")

archp <- loadArchRProject(file.path(BASE,"scatac_ArchR"), showLogo=FALSE)
cd <- getCellColData(archp)

## ---- gene score: BAP1 + chr3p accessibility per cell ----
gsm <- getMatrixFromProject(archp, useMatrix="GeneScoreMatrix")
rd  <- rowData(gsm); gs <- assay(gsm); rownames(gs) <- rd$name
bap1_gs_all <- gs["BAP1", ]
chr3p_genes <- intersect(rd$name[as.character(rd$seqnames)=="chr3" & rd$start<CHR3_CENT], rownames(gs))
chr3p_acc_all <- colMeans(gs[chr3p_genes,,drop=FALSE])

## ---- chromVAR TF deviations (scaled, matching original wilcoxauc metric) ----
mSE <- getMatrixFromProject(archp, useMatrix="MotifMatrix")
mDev <- as.matrix(assays(mSE)[[1]])                 # deviations
rownames(mDev) <- sub("_.*","", rownames(mDev))     # clean TF names (drop _idx)
# align cells
common <- intersect(colnames(mDev), names(bap1_gs_all))
mDev <- mDev[, common]; bap1_gs <- bap1_gs_all[common]; chr3p_acc <- chr3p_acc_all[common]
samp <- cd[common,"Sample"]
mDev_sc <- t(scale(t(mDev)))                         # scale each TF across cells
mDev_sc[!is.finite(mDev_sc)] <- 0

## ---- per-sample BAP1-high vs BAP1-low differential TF activity ----
tum <- c("P1","P3","P4","P5","P8","P10","P11","P12","P13","P14","P23")
run_sample <- function(s) {
  idx <- which(samp == s); if (length(idx) < 60) return(NULL)
  b <- bap1_gs[idx]
  # split: BAP1-accessible (gs>0) vs BAP1-silent (gs==0); if too few silent, use median
  hi <- b > 0; lo <- b == 0
  if (sum(hi) < 25 || sum(lo) < 25) {               # fall back to top vs bottom tertile
    q <- quantile(b, c(1/3, 2/3)); hi <- b >= q[2]; lo <- b <= q[1]
    method <- "tertile"
  } else method <- "accessible_vs_silent"
  if (sum(hi) < 25 || sum(lo) < 25) return(NULL)
  y <- rep(NA, length(idx)); y[hi] <- "BAP1_high"; y[lo] <- "BAP1_low"
  keep <- !is.na(y)
  res <- wilcoxauc(mDev_sc[, idx[keep], drop=FALSE], y[keep])
  res <- res[res$group == "BAP1_high", ]
  data.frame(sample=s, method=method, n_high=sum(hi), n_low=sum(lo),
             TF=res$feature, auc=res$auc, logFC=res$logFC,
             pval=res$pval, padj=res$padj,
             chr3p_acc_hi=mean(chr3p_acc[idx][hi]), chr3p_acc_lo=mean(chr3p_acc[idx][lo]))
}
all_res <- do.call(rbind, lapply(tum, run_sample))
write.csv(all_res, "within_sample_TFactivity_BAP1_high_vs_low.csv", row.names=FALSE)
tested <- unique(all_res$sample)
cat("Samples tested (n_high/n_low):\n")
print(unique(all_res[,c("sample","method","n_high","n_low")]), row.names=FALSE)

## ---- summarize: which TFs differ, and is it consistent across samples? ----
# per-sample RUNX result
cat("\n=== RUNX within-sample (auc>0.5 = higher in BAP1-high) ===\n")
print(all_res[all_res$TF %in% runx, c("sample","TF","auc","logFC","pval","padj")],
      row.names=FALSE, digits=3)

# meta: mean signed effect (auc-0.5) across samples per TF, and how many samples sig
sig <- all_res$padj < 0.05
agg <- aggregate(cbind(eff = all_res$auc - 0.5) ~ TF, data=all_res, FUN=mean)
nsig <- aggregate(sig ~ TF, data=all_res, FUN=sum); colnames(nsig)[2] <- "n_sig_samples"
ndir <- aggregate((all_res$auc>0.5) ~ all_res$TF, FUN=function(x) max(sum(x), sum(!x)))
colnames(ndir) <- c("TF","n_consistent_dir")
meta <- Reduce(function(a,b) merge(a,b,by="TF"), list(agg, nsig, ndir))
meta$n_samples <- length(tested)
meta <- meta[order(-abs(meta$eff)), ]
write.csv(meta, "within_sample_TF_meta.csv", row.names=FALSE)
cat("\n=== top consistent TFs across within-sample tests ===\n")
print(head(meta[meta$n_sig_samples>=2, ], 25), row.names=FALSE, digits=3)

## ---- PLOTS ----
# (1) per-sample chr3p CNV clonality context (from probe csv)
cnv <- read.csv("persample_chr3p_bimodality.csv")
cnv$clonality <- ifelse(cnv$frac_retained_chr3p>0.25 & cnv$frac_lost_chr3p>0.25, "subclonal",
                 ifelse(cnv$frac_lost_bap1>0.4, "clonal-lost", "clonal-retained"))
scr <- read.csv("allsample_scRNA_chr3p_cnv.csv")
scr <- scr[scr$sample %in% tum, ]
scr$sample <- factor(scr$sample, levels=cnv$sample[order(cnv$chr3p_med)])
pV <- ggplot(scr, aes(sample, chr3p_norm, fill=sample)) +
  geom_violin(scale="width", linewidth=.2) +
  geom_hline(yintercept=1, linetype="dashed", color="grey40") +
  geom_hline(yintercept=0.92, linetype="dotted", color="red") +
  guides(fill="none") + theme_bw(base_size=10) +
  ylab("chr3p inferCNV (norm to genome median)") + xlab("") +
  ggtitle("P8 is clonally chr3p/BAP1-lost; P11/P5 are subclonal")
ggsave("Plots/within_P8_chr3p_clonality_violin.pdf", pV, width=7, height=3.6)

# (2) volcano for the subclonal samples (P11, P5) + P8
for (s in intersect(c("P11","P5","P8"), tested)) {
  d <- all_res[all_res$sample==s, ]; d$sig <- d$padj<0.05
  d$lab <- ifelse(d$TF %in% runx | (d$sig & abs(d$auc-0.5)>0.12), d$TF, NA)
  pv <- ggplot(d, aes(auc-0.5, -log10(pval))) +
    geom_point(aes(color=sig), size=.8) +
    scale_color_manual(values=c(`FALSE`="grey70",`TRUE`="#c0392b"), guide="none") +
    geom_vline(xintercept=0, color="grey50") +
    ggrepel::geom_text_repel(aes(label=lab), size=2.4, max.overlaps=30) +
    theme_bw(base_size=10) +
    xlab("AUC-0.5  (>0 = higher in BAP1-high cells)") + ylab("-log10 p") +
    ggtitle(paste0("Within ", s, " : TF activity BAP1-high vs BAP1-low"),
            paste0(unique(d$method), "  n_high=", d$n_high[1], " n_low=", d$n_low[1]))
  ggsave(paste0("Plots/within_", s, "_TF_volcano.pdf"), pv, width=5, height=4.5)
}

# (3) RUNX across samples
rd2 <- all_res[all_res$TF %in% runx, ]
rd2$sample <- factor(rd2$sample, levels=cnv$sample[order(cnv$chr3p_med)])
pR <- ggplot(rd2, aes(sample, auc-0.5, fill=TF)) +
  geom_col(position="dodge") + geom_hline(yintercept=0, color="grey40") +
  theme_bw(base_size=10) + ylab("AUC-0.5 (BAP1-high vs low)") + xlab("") +
  ggtitle("RUNX within-sample effect by tumor (no consistent BAP1-high/low direction)")
ggsave("Plots/within_sample_RUNX_by_tumor.pdf", pR, width=7, height=3.6)

cat("\nDONE\n")
