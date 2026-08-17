###############################################################
# R2_Q14 : inferCNV subclones (scRNA) <-> scATAC concordance.
#   (1) In scRNA, define per-cell BAP1 status from chr3p inferCNV
#       (retained vs lost) and VALIDATE it with BAP1 mRNA.
#   (2) Show which tumors are subclonal (mix of retained+lost) vs clonal.
#   (3) MATCH to scATAC: per-sample fraction BAP1-retained (scRNA CNV) vs
#       fraction BAP1-accessible (scATAC gene score).
###############################################################
suppressMessages({ library(Seurat); library(ggplot2); library(ggrepel) })

BASE   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment"
OUTDIR <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q14"
setwd(OUTDIR); dir.create("Plots", showWarnings = FALSE)

cnv <- read.csv("allsample_scRNA_chr3p_cnv.csv", stringsAsFactors=FALSE)  # cell,sample,chr3p_norm,bap1_norm
# per-cell scRNA CNV-based BAP1 status
cnv$cnv_status <- ifelse(cnv$bap1_norm > 0.98, "retained",
                  ifelse(cnv$bap1_norm < 0.90, "lost", "intermediate"))

## ---- BAP1 mRNA from the Seurat object ----
srt <- readRDS(file.path(BASE, "scrna", "srt.rds"))
DefaultAssay(srt) <- "RNA"
bap1_mrna <- tryCatch(as.numeric(GetAssayData(srt, assay="RNA", layer="data")["BAP1", ]),
              error=function(e) as.numeric(GetAssayData(srt, assay="RNA", slot="data")["BAP1", ]))
names(bap1_mrna) <- colnames(srt)
cnv$BAP1_mrna <- bap1_mrna[cnv$cell]
cat("cells with BAP1 mRNA matched:", sum(!is.na(cnv$BAP1_mrna)), "/", nrow(cnv), "\n")

## ---- (1) validate: BAP1 mRNA higher in CNV-retained cells (per subclonal sample) ----
val <- do.call(rbind, lapply(sort(unique(cnv$sample)), function(s){
  d <- cnv[cnv$sample==s & cnv$cnv_status %in% c("retained","lost") & is.finite(cnv$BAP1_mrna), ]
  if (sum(d$cnv_status=="retained")<15 || sum(d$cnv_status=="lost")<15)
    return(data.frame(sample=s, n_ret=sum(d$cnv_status=="retained"), n_lost=sum(d$cnv_status=="lost"),
                      mrna_ret=NA, mrna_lost=NA, wilcox_p=NA))
  wp <- suppressWarnings(wilcox.test(BAP1_mrna ~ cnv_status, data=d)$p.value)
  data.frame(sample=s, n_ret=sum(d$cnv_status=="retained"), n_lost=sum(d$cnv_status=="lost"),
    mrna_ret=mean(d$BAP1_mrna[d$cnv_status=="retained"]),
    mrna_lost=mean(d$BAP1_mrna[d$cnv_status=="lost"]), wilcox_p=wp)
}))
write.csv(val, "scRNA_CNV_BAP1mrna_validation.csv", row.names=FALSE)
cat("\n=== BAP1 mRNA in CNV-retained vs CNV-lost cells (per sample) ===\n")
print(val, row.names=FALSE, digits=3)

## ---- (2) per-sample clonality summary ----
clon <- do.call(rbind, lapply(sort(unique(cnv$sample)), function(s){
  d <- cnv[cnv$sample==s, ]
  data.frame(sample=s, n=nrow(d),
    frac_retained=mean(d$bap1_norm>0.98), frac_lost=mean(d$bap1_norm<0.90),
    frac_interm=mean(d$bap1_norm>=0.90 & d$bap1_norm<=0.98))
}))
clon$class <- ifelse(clon$frac_retained>0.25 & clon$frac_lost>0.25, "SUBCLONAL",
              ifelse(clon$frac_lost>0.4, "clonal-lost",
              ifelse(clon$frac_retained>0.6, "clonal-retained", "mixed/intermediate")))
write.csv(clon, "scRNA_persample_clonality.csv", row.names=FALSE)
cat("\n=== per-sample BAP1/chr3p clonality (scRNA) ===\n"); print(clon, row.names=FALSE, digits=3)

## ---- (3) cross-modal: scRNA frac-retained vs scATAC frac BAP1-accessible ----
atac <- read.csv("allsample_scATAC_bap1_chr3p.csv", stringsAsFactors=FALSE)
atac_frac <- do.call(rbind, lapply(sort(unique(atac$sample)), function(s){
  d <- atac[atac$sample==s, ]
  data.frame(sample=s, n_atac=nrow(d), frac_bap1_accessible=mean(d$BAP1_gs>0))
}))
cm <- merge(clon[,c("sample","frac_retained","class")], atac_frac, by="sample")
write.csv(cm, "crossmodal_scRNA_vs_scATAC_BAP1.csv", row.names=FALSE)
cat("\n=== cross-modal (scRNA CNV-retained frac vs scATAC BAP1-accessible frac) ===\n")
print(cm, row.names=FALSE, digits=3)
cmr <- suppressWarnings(cor.test(cm$frac_retained, cm$frac_bap1_accessible, method="spearman"))
cat("Spearman rho:", round(cmr$estimate,3), " p:", signif(cmr$p.value,3), "\n")

## ---- PLOTS ----
# validation: BAP1 mRNA by CNV status for the subclonal / mixed samples
subs <- clon$sample[clon$class %in% c("SUBCLONAL")]
pd <- cnv[cnv$sample %in% union(subs, c("P8","P5","P11")) &
          cnv$cnv_status %in% c("retained","lost") & is.finite(cnv$BAP1_mrna), ]
pd$cnv_status <- factor(pd$cnv_status, levels=c("lost","retained"))
pV <- ggplot(pd, aes(cnv_status, BAP1_mrna, fill=cnv_status)) +
  geom_violin(scale="width", linewidth=.2) + geom_boxplot(width=.15, outlier.shape=NA, alpha=.5) +
  facet_wrap(~sample, scales="free_y", nrow=1) +
  scale_fill_manual(values=c(lost="#c0392b", retained="#2471a3"), guide="none") +
  theme_bw(base_size=10) + xlab("chr3p/BAP1 inferCNV status") + ylab("BAP1 mRNA (lognorm)") +
  ggtitle("Validation: CNV-retained cells express more BAP1 mRNA")
ggsave("Plots/scRNA_CNV_BAP1mrna_validation.pdf", pV, width=8, height=3.4)

# cross-modal scatter
cm$class <- factor(cm$class)
pC <- ggplot(cm, aes(frac_retained, frac_bap1_accessible)) +
  geom_smooth(method="lm", se=FALSE, color="grey55", linetype="dashed", linewidth=.4) +
  geom_point(aes(color=class), size=3) + ggrepel::geom_text_repel(aes(label=sample), size=2.6) +
  theme_bw(base_size=10) +
  xlab("scRNA: fraction chr3p/BAP1-retained (inferCNV)") +
  ylab("scATAC: fraction BAP1-accessible (gene score>0)") +
  ggtitle("Cross-modal concordance of subclonal BAP1 status")
ggsave("Plots/crossmodal_scRNA_vs_scATAC.pdf", pC, width=5, height=4)

cat("\nDONE\n")
