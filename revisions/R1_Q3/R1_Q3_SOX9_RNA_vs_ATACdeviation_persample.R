#!/usr/bin/env Rscript
# R1_Q3_SOX9_RNA_vs_ATACdeviation_persample.R
# Rigorous cross-patient test of "does SOX9 expression track SOX9 motif activity?"
# Patient-matched, pseudobulk-per-sample (ATAC and RNA are different cells from
# the same patients; P23 has no scRNA so it drops out):
#   x = per-sample mean SOX9 scRNA expression (srt.rds, RNA data layer)
#   y = per-sample mean SOX9 chromVAR deviation z (full tumor ATAC scatac_ArchR)
# Correlate across matched samples; compare to RUNX2/AP-1 (real drivers); and
# PARTIAL-correlate controlling for the epithelioid->sarcomatoid scS axis to test
# whether any coupling is independent of the program both co-vary with.

suppressPackageStartupMessages({library(ArchR); library(Seurat); library(ggplot2); library(ggrepel)})
addArchRThreads(4)
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUT  <- file.path(ROOT, "git_repo_claude/R1_Q3")
ATAC <- file.path(ROOT, "tumor_compartment/scatac_ArchR")
SCRN <- file.path(ROOT, "tumor_compartment/scrna/srt.rds")
ORD  <- file.path(ROOT, "tumor_compartment/scrna/cnmf20_sarcomatoid_sample_order.csv")
GENES <- c("SOX9","RUNX2","FOSL2","JUN")

## ---- scS program axis + tumor sample list ----
sarc <- read.csv(ORD, stringsAsFactors = FALSE)
tumor_sams <- sarc$sampleID[!sarc$sampleID %in% c("HU37","HU62","normal_pleura")]
scS_axis <- setNames(sarc$x, sarc$sampleID)

## ---- ATAC: per-sample mean chromVAR deviation z ----
message("Loading ATAC project ...")
proj <- loadArchRProject(ATAC, showLogo = FALSE)
mm <- getMatrixFromProject(proj, useMatrix = "MotifMatrix")
Z  <- assays(mm)[["z"]]
samp <- setNames(proj$Sample, proj$cellNames)[colnames(Z)]
motif_of <- function(sym){ h<-grep(paste0("^",sym,"_[0-9]"), rownames(Z)); if(length(h)) rownames(Z)[h[1]] else NA }
atac_tab <- list()
for (g in GENES) {
  m <- motif_of(g); if (is.na(m)) next
  z <- as.numeric(Z[m, ])
  agg <- tapply(z, samp, mean, na.rm=TRUE)
  atac_tab[[g]] <- data.frame(sampleID=names(agg), gene=g, atac_dev=as.numeric(agg))
}
atac_tab <- do.call(rbind, atac_tab)
atac_n <- tapply(samp, samp, length)
message("ATAC samples: ", paste(names(atac_n), collapse=", "))

## ---- RNA: per-sample mean SOX9 (+others) expression ----
message("Loading scRNA srt.rds ...")
srt <- readRDS(SCRN)
srt <- srt[, srt$sampleID %in% tumor_sams]
expr <- GetAssayData(srt, assay="RNA", layer="data")
sid  <- as.character(srt$sampleID)
rna_tab <- list()
for (g in GENES) {
  if (!g %in% rownames(expr)) next
  e <- as.numeric(expr[g, ])
  agg <- tapply(e, sid, mean, na.rm=TRUE)
  rna_tab[[g]] <- data.frame(sampleID=names(agg), gene=g, rna_expr=as.numeric(agg))
}
rna_tab <- do.call(rbind, rna_tab)
rna_n <- tapply(sid, sid, length)
message("RNA samples: ", paste(names(rna_n), collapse=", "))

## ---- merge + correlate per gene ----
res <- list(); merged_all <- list()
for (g in GENES) {
  a <- atac_tab[atac_tab$gene==g, ]; r <- rna_tab[rna_tab$gene==g, ]
  m <- merge(a, r, by="sampleID")
  m <- m[m$sampleID %in% tumor_sams, ]           # matched tumor samples only
  m$n_rna <- rna_n[m$sampleID]; m$n_atac <- atac_n[m$sampleID]
  m$scS <- scS_axis[m$sampleID]; m$gene <- g
  if (nrow(m) < 4) next
  rp <- cor(m$rna_expr, m$atac_dev, method="pearson")
  rs <- cor(m$rna_expr, m$atac_dev, method="spearman")
  # partial correlation controlling for scS program axis
  pc <- with(m, {
    rxy<-cor(rna_expr,atac_dev); rxz<-cor(rna_expr,scS); ryz<-cor(atac_dev,scS)
    (rxy-rxz*ryz)/sqrt((1-rxz^2)*(1-ryz^2)) })
  res[[g]] <- data.frame(gene=g, n_samples=nrow(m), pearson=round(rp,3),
                         spearman=round(rs,3), partial_given_scS=round(pc,3))
  merged_all[[g]] <- m[, c("gene","sampleID","rna_expr","atac_dev","scS","n_rna","n_atac")]
  if (g=="SOX9") write.csv(m, file.path(OUT,"SOX9_RNA_vs_ATACdev_persample_table.csv"), row.names=FALSE)
}
summ <- do.call(rbind, res)
merged_all <- do.call(rbind, merged_all)
write.csv(merged_all, file.path(OUT,"RNA_vs_ATACdev_persample_allgenes_table.csv"), row.names=FALSE)

## ---- multi-panel: SOX9 (passenger) vs RUNX2/AP-1 (drivers) ----
ord <- c("SOX9","RUNX2","JUN","FOSL2"); ord <- ord[ord %in% summ$gene]
merged_all$gene <- factor(merged_all$gene, levels=ord)
labr <- setNames(sprintf("%s  raw r=%.2f | partial|scS=%.2f",
                         summ$gene, summ$pearson, summ$partial_given_scS), summ$gene)
merged_all$panel <- labr[as.character(merged_all$gene)]
merged_all$panel <- factor(merged_all$panel, levels=labr[ord])
pmw <- ggplot(merged_all, aes(rna_expr, atac_dev)) +
  geom_smooth(method="lm", se=TRUE, color="grey40", linewidth=.5) +
  geom_point(aes(size=n_atac, color=scS)) +
  ggrepel::geom_text_repel(aes(label=sampleID), size=2.4, max.overlaps=20) +
  scale_color_gradient2(low="#2166AC", mid="#FFFFBF", high="#D6604D",
                        midpoint=median(merged_all$scS)) +
  facet_wrap(~panel, scales="free", ncol=2) +
  labs(title="TF scRNA expression vs ATAC chromVAR deviation, per patient (n=9)",
       subtitle="Passenger (SOX9: corr vanishes after removing scS program axis) vs drivers (RUNX2/AP-1: corr survives)",
       x="mean scRNA expression", y="mean ATAC deviation z",
       size="n ATAC cells", color="scS axis") +
  theme_bw(base_size=11) + theme(strip.text=element_text(face="bold", size=9))
ggsave(file.path(OUT,"TF_RNA_vs_ATACdeviation_multipanel.pdf"), pmw, width=9, height=8)
cat("\n=== per-sample RNA expression vs ATAC chromVAR deviation (matched patients) ===\n")
print(summ)
write.csv(summ, file.path(OUT, "RNA_vs_ATACdeviation_persample_summary.csv"), row.names=FALSE)

## ---- scatter for SOX9 ----
m <- read.csv(file.path(OUT,"SOX9_RNA_vs_ATACdev_persample_table.csv"), stringsAsFactors=FALSE)
p <- ggplot(m, aes(rna_expr, atac_dev)) +
  geom_smooth(method="lm", se=TRUE, color="grey40", linewidth=.6) +
  geom_point(aes(size=n_rna, color=scS)) +
  ggrepel::geom_text_repel(aes(label=sampleID), size=3) +
  scale_color_gradient2(low="#2166AC", mid="#FFFFBF", high="#D6604D", midpoint=median(m$scS)) +
  labs(title="SOX9: scRNA expression vs ATAC chromVAR deviation (per patient)",
       subtitle=sprintf("Pearson r=%.2f  Spearman=%.2f  partial|scS=%.2f  (n=%d samples)",
                        cor(m$rna_expr,m$atac_dev), cor(m$rna_expr,m$atac_dev,method="spearman"),
                        summ$partial_given_scS[summ$gene=="SOX9"], nrow(m)),
       x="mean SOX9 scRNA expression", y="mean SOX9 ATAC deviation z",
       size="n RNA cells", color="scS axis") +
  theme_bw(base_size=12)
ggsave(file.path(OUT,"SOX9_RNA_vs_ATACdeviation_persample.pdf"), p, width=6, height=5)
cat("\nWrote summary CSV + per-sample table + scatter PDF\n")
