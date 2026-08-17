###############################################################
# R2_Q14 : headline figure for the within-sample / genetic BAP1 analysis
# (plot-only, rebuilt from saved CSVs — no ArchR/Seurat reload).
###############################################################
suppressMessages({ library(ggplot2); library(ggrepel); library(patchwork) })
OUTDIR <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q14"
setwd(OUTDIR)
runx <- c("RUNX1","RUNX2","RUNX3","CBFB")
th <- theme_bw(base_size=9) + theme(plot.title=element_text(size=9,face="bold"),
       plot.subtitle=element_text(size=7.5), legend.key.size=unit(3,"mm"),
       legend.text=element_text(size=7), legend.title=element_text(size=7.5))

## data
res  <- read.csv("within_sample_TFactivity_BAP1_high_vs_low.csv", stringsAsFactors=FALSE)
cnvs <- read.csv("persample_chr3p_bimodality.csv", stringsAsFactors=FALSE)
scr  <- read.csv("allsample_scRNA_chr3p_cnv.csv", stringsAsFactors=FALSE)
clon <- read.csv("scRNA_persample_clonality.csv", stringsAsFactors=FALSE)
cm   <- read.csv("crossmodal_scRNA_vs_scATAC_BAP1.csv", stringsAsFactors=FALSE)

## (A) chr3p clonality violin, ordered by median, honest title
ord <- cnvs$sample[order(cnvs$chr3p_med)]
scr2 <- scr[scr$sample %in% ord, ]; scr2$sample <- factor(scr2$sample, levels=ord)
pA <- ggplot(scr2, aes(sample, chr3p_norm, fill=sample)) +
  geom_violin(scale="width", linewidth=.2) +
  geom_hline(yintercept=1, linetype="dashed", color="grey40") +
  geom_hline(yintercept=0.92, linetype="dotted", color="red") +
  guides(fill="none") + th + ylab("chr3p inferCNV\n(norm. to genome median)") + xlab("") +
  ggtitle("A  Clonal vs subclonal BAP1 loss (scRNA inferCNV)",
          "P8/P5 clonally lost (unimodal <1); P11 subclonal (bimodal); rest retained")

## (B) cross-modal concordance
cm$class <- factor(cm$class)
rho <- suppressWarnings(cor.test(cm$frac_retained, cm$frac_bap1_accessible, method="spearman"))
pB <- ggplot(cm, aes(frac_retained, frac_bap1_accessible)) +
  geom_smooth(method="lm", se=FALSE, color="grey55", linetype="dashed", linewidth=.4) +
  geom_point(aes(color=class), size=2.6) + geom_text_repel(aes(label=sample), size=2.4) +
  scale_color_manual(values=c(`clonal-lost`="#c0392b",`clonal-retained`="#2471a3",
                     `mixed/intermediate`="#16a085",`SUBCLONAL`="#8e44ad"), name=NULL) +
  th + xlab("scRNA: fraction chr3p/BAP1-retained") +
  ylab("scATAC: fraction BAP1-accessible") +
  ggtitle("B  Cross-modal concordance of BAP1 status",
          sprintf("Spearman rho=%.2f, p=%.2f (n=9)", rho$estimate, rho$p.value))

## (C) within-P11 volcano (RUNX1 down in BAP1-high subclone)
d <- res[res$sample=="P11", ]; d$sig <- d$padj<0.05
d$lab <- ifelse(d$TF %in% runx | (d$sig & abs(d$auc-0.5)>0.09), d$TF, NA)
pC <- ggplot(d, aes(auc-0.5, -log10(pval))) +
  geom_point(aes(color=sig), size=.7) +
  scale_color_manual(values=c(`FALSE`="grey75",`TRUE`="#c0392b"), guide="none") +
  geom_vline(xintercept=0, color="grey50") +
  geom_text_repel(aes(label=lab), size=2.2, max.overlaps=25) +
  th + xlab("AUC-0.5  (>0 = higher in BAP1-high cells)") + ylab("-log10 p") +
  ggtitle("C  Within P11 (subclonal): TF activity BAP1-high vs -low",
          "hundreds of TFs shift; RUNX1 is LOWER in the BAP1-retained subclone")

## (D) RUNX within-sample effect by tumor (scATAC-tested samples only, no NA)
rd <- res[res$TF %in% runx, ]
sc_ord <- intersect(cnvs$sample[order(cnvs$chr3p_med)], unique(rd$sample))
sc_ord <- c(sc_ord, setdiff(unique(rd$sample), sc_ord))   # append P10,P23 (no scRNA)
rd$sample <- factor(rd$sample, levels=sc_ord); rd$TF <- factor(rd$TF, levels=runx)
pD <- ggplot(rd, aes(sample, auc-0.5, fill=TF)) +
  geom_col(position="dodge") + geom_hline(yintercept=0, color="grey40") +
  scale_fill_brewer(palette="Set2") +
  th + ylab("AUC-0.5 (BAP1-high vs low)") + xlab("") +
  ggtitle("D  RUNX within-sample effect by tumor",
          "no consistent BAP1-high/low direction across tumors")

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(title="R2_Q14 — Within-tumor BAP1-high vs -low TF activity (histology held constant)",
    theme=theme(plot.title=element_text(size=11,face="bold")))
ggsave("Plots/R2_Q14_within_sample_headline.pdf", fig, width=12, height=8)
# also refresh the two standalone fixed plots
ggsave("Plots/within_chr3p_clonality_violin.pdf", pA, width=7, height=3.4)
ggsave("Plots/within_sample_RUNX_by_tumor.pdf", pD, width=7, height=3.4)
cat("DONE\n")
