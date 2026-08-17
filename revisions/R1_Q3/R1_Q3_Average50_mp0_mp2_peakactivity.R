#!/usr/bin/env Rscript
# R1_Q3_Average50_mp0_mp2_peakactivity.R
# Per-cell accessibility activity in the P23 ArchR object of:
#   peaks containing FI-NeMo hits for Average_50 . mp2  vs
#   peaks containing FI-NeMo hits for Average_50 . mp0
# (HDMA merged-modisco sub-patterns; hits pooled across SOX9_high_P23 + SOX9_low_P23
# FI-NeMo merged_counts). ArchR peaks overlapping >=1 hit define each set; per-cell
# activity = depth-normalised accessibility (fraction-in-set and mean-per-peak).

suppressPackageStartupMessages({library(ArchR); library(Matrix); library(ggplot2)})
addArchRThreads(4)
PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
MP0  <- file.path(OUT, "Average50_mp0_hits.bed")   # pre-extracted chr start end
MP2  <- file.path(OUT, "Average50_mp2_hits.bed")

proj <- loadArchRProject(PROJ, showLogo = FALSE)
ps <- getPeakSet(proj)

gr_from_bed <- function(f){
  d <- read.table(f, sep="\t", header=FALSE, col.names=c("chr","start","end"))
  GRanges(d$chr, IRanges(d$start+1, d$end))           # bed 0-based -> 1-based
}
mp0 <- gr_from_bed(MP0); mp2 <- gr_from_bed(MP2)
peaks_mp0 <- sort(unique(queryHits(findOverlaps(ps, mp0))))
peaks_mp2 <- sort(unique(queryHits(findOverlaps(ps, mp2))))
cat(sprintf("Average_50 mp0: %d hits -> %d ArchR peaks\n", length(mp0), length(peaks_mp0)))
cat(sprintf("Average_50 mp2: %d hits -> %d ArchR peaks\n", length(mp2), length(peaks_mp2)))
cat(sprintf("overlap of the two peak sets: %d peaks\n", length(intersect(peaks_mp0,peaks_mp2))))

## ---- per-cell activity from PeakMatrix ----
pm_se <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
# rows of PeakMatrix follow the peakSet order within each seqnames group; align by GRanges
pmGR <- rowRanges(pm_se)
ov <- findOverlaps(pmGR, ps, type="equal")
# map peakSet index -> PeakMatrix row
ps2row <- integer(length(ps)); ps2row[subjectHits(ov)] <- queryHits(ov)
rows_mp0 <- ps2row[peaks_mp0]; rows_mp0 <- rows_mp0[rows_mp0>0]
rows_mp2 <- ps2row[peaks_mp2]; rows_mp2 <- rows_mp2[rows_mp2>0]
M <- assay(pm_se)                                   # peaks x cells (counts)
depth <- Matrix::colSums(M)                          # total reads-in-peaks per cell
norm <- function(rows){
  sub <- M[rows,,drop=FALSE]
  list(frac = Matrix::colSums(sub)/depth,            # fraction of in-peak signal in this set
       meanpp = Matrix::colSums(sub)/length(rows)/depth*1e4)  # mean per-peak, depth-norm (x1e4)
}
a0 <- norm(rows_mp0); a2 <- norm(rows_mp2)

cc <- getCellColData(proj)
grpcol <- grep("SOX9", colnames(cc), value=TRUE, ignore.case=TRUE)
grp <- if (length(grpcol)) as.character(cc[[grpcol[1]]]) else rep("P23", nCells(proj))
df <- data.frame(cell=rownames(cc), group=grp, depth=depth,
                 act_mp0_frac=a0$frac, act_mp2_frac=a2$frac,
                 act_mp0_meanpp=a0$meanpp, act_mp2_meanpp=a2$meanpp)
write.csv(df, file.path(OUT,"Average50_mp0_mp2_percell_activity.csv"), row.names=FALSE)

cat("\n=== per-cell activity summary (mean across cells) ===\n")
cat(sprintf("n cells = %d ; mp0 peaks = %d ; mp2 peaks = %d\n", nrow(df), length(rows_mp0), length(rows_mp2)))
cat(sprintf("fraction-in-set : mp0 = %.4f  mp2 = %.5f\n", mean(df$act_mp0_frac), mean(df$act_mp2_frac)))
cat(sprintf("mean-per-peak   : mp0 = %.4f  mp2 = %.4f  (depth-norm x1e4)\n", mean(df$act_mp0_meanpp), mean(df$act_mp2_meanpp)))
ct <- wilcox.test(df$act_mp2_meanpp, df$act_mp0_meanpp, paired=TRUE)
cat(sprintf("paired Wilcoxon mp2 vs mp0 (mean-per-peak): p=%.3g\n", ct$p.value))
if (length(grpcol)){
  cat("\nby group (mean-per-peak):\n")
  print(aggregate(cbind(act_mp0_meanpp,act_mp2_meanpp)~group, df, mean))
}

## ---- plots ----
library(reshape2)
pl <- melt(df[,c("group","act_mp0_meanpp","act_mp2_meanpp")], id.vars="group")
pl$variable <- factor(pl$variable, labels=c("Average_50 . mp0","Average_50 . mp2"))
p <- ggplot(pl, aes(variable, value, fill=variable)) +
  geom_violin(scale="width", alpha=.8, linewidth=.3) + geom_boxplot(width=.12, outlier.shape=NA, fill="white", linewidth=.3) +
  facet_wrap(~group, scales="free_y") +
  scale_fill_manual(values=c("#4C78A8","#E45756"), guide="none") +
  labs(title="Per-cell accessibility activity: Average_50 mp0 vs mp2 peak sets (P23)",
       subtitle=sprintf("mean per-peak depth-normalised accessibility; mp0=%d peaks, mp2=%d peaks",
                        length(rows_mp0), length(rows_mp2)),
       x=NULL, y="activity (mean per-peak, x1e4)") +
  theme_bw(base_size=11)
ggsave(file.path(OUT,"Average50_mp0_mp2_percell_activity.pdf"), p, width=8, height=4.5)
cat("\nWrote per-cell activity csv + plot\n")
