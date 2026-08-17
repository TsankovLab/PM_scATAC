#!/usr/bin/env Rscript
# R1_Q3_SOX9_chromVAR_by_gate.R
# Verify the chromVAR SOX9 motif deviation is higher in SOX9-high (module) cells
# than SOX9-low, on the SAME gated project used for ChromBPNet + footprinting.
# Reports per-gate mean/median deviation z-score + Wilcoxon, and a violin plot.

suppressPackageStartupMessages({library(ArchR); library(ggplot2)})
addArchRThreads(4)
PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_module_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
proj <- loadArchRProject(PROJ, showLogo = FALSE)

mm <- getMatrixFromProject(proj, useMatrix = "MotifMatrix")
# chromVAR stores two seqnames: 'deviations' and 'z' (z-scores). Use z.
rn <- rowData(mm)
sox_rows <- grep("SOX9", rownames(mm), ignore.case = TRUE, value = TRUE)
cat("SOX9 motif rows:\n"); print(sox_rows)

zmat <- assays(mm)[["z"]]
if (is.null(zmat)) zmat <- assays(mm)[["deviations"]]
gate <- proj$SOX9_module_gate
names(gate) <- proj$cellNames
# align cells
gate <- gate[colnames(zmat)]

res <- list()
for (sr in sox_rows) {
  z <- as.numeric(zmat[sr, ])
  hi <- z[gate == "SOX9pos"]; lo <- z[gate == "SOX9neg"]
  wt <- wilcox.test(hi, lo)
  cat(sprintf("\n== %s ==\n", sr))
  cat(sprintf("  SOX9-high (n=%d): mean z=%.3f  median z=%.3f\n", length(hi), mean(hi), median(hi)))
  cat(sprintf("  SOX9-low  (n=%d): mean z=%.3f  median z=%.3f\n", length(lo), mean(lo), median(lo)))
  cat(sprintf("  delta mean z (high-low) = %.3f ; Wilcoxon p = %.3g\n", mean(hi)-mean(lo), wt$p.value))
  res[[sr]] <- data.frame(motif=sr, n_high=length(hi), n_low=length(lo),
                          mean_high=mean(hi), mean_low=mean(lo),
                          median_high=median(hi), median_low=median(lo),
                          delta_mean=mean(hi)-mean(lo), wilcox_p=wt$p.value)
}
resdf <- do.call(rbind, res)
write.csv(resdf, file.path(OUT, "SOX9_chromVAR_deviation_by_gate.csv"), row.names = FALSE)

# ---- context: rank SOX9's high-vs-low shift against ALL motifs ----
# (is SOX9 a standout, or does the whole landscape shift with the gate?)
gpos <- gate == "SOX9pos"; gneg <- gate == "SOX9neg"
delta_all <- rowMeans(zmat[, gpos, drop=FALSE]) - rowMeans(zmat[, gneg, drop=FALSE])
delta_all <- sort(delta_all, decreasing = TRUE)
nmot <- length(delta_all)
sox_rank <- which(names(delta_all) == sox_rows[1])
cat(sprintf("\n--- SOX9 shift in context of all %d motifs ---\n", nmot))
cat(sprintf("SOX9_756 delta-mean-z = %.3f  -> rank %d / %d (top %.0f%%)\n",
            delta_all[sox_rows[1]], sox_rank, nmot, 100*sox_rank/nmot))
cat("Fraction of motifs with higher deviation in SOX9-high (delta>0): ",
    sprintf("%.0f%%\n", 100*mean(delta_all > 0)))
cat("\nTop 15 motifs more accessible in SOX9-high (co-shifting TFs):\n")
print(round(head(delta_all, 15), 3))
write.csv(data.frame(motif=names(delta_all), delta_mean_z=as.numeric(delta_all),
                     rank=seq_len(nmot)),
          file.path(OUT, "chromVAR_all_motifs_high_vs_low_delta.csv"), row.names = FALSE)

# violin for the primary SOX9 motif
sr <- sox_rows[1]
df <- data.frame(z = as.numeric(zmat[sr, ]), gate = gate)
df <- df[df$gate %in% c("SOX9pos","SOX9neg"), ]
df$gate <- factor(df$gate, levels = c("SOX9neg","SOX9pos"))
p <- ggplot(df, aes(gate, z, fill = gate)) +
  geom_violin(scale="width", trim=TRUE, alpha=.7) +
  geom_boxplot(width=.15, outlier.size=.3, alpha=.9) +
  scale_fill_manual(values=c(SOX9neg="#377EB8", SOX9pos="#E41A1C")) +
  labs(title=paste0("chromVAR deviation z: ", sr),
       subtitle="SOX9 module-gated P23 malignant cells", y="chromVAR deviation z", x=NULL) +
  theme_bw(base_size=12) + theme(legend.position="none")
ggsave(file.path(OUT, "SOX9_chromVAR_deviation_by_gate.pdf"), p, width=4.5, height=4.5)
cat("\nWrote SOX9_chromVAR_deviation_by_gate.{csv,pdf}\n")
