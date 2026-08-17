#!/usr/bin/env Rscript
# Distribution of SOX9+ vs SOX9- cells from the module-score gate.
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
setDTthreads(1)
OUT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
POS <- 0.20; NEG <- -0.10
d <- fread(file.path(OUT, "SOX9_module_gated_cells.csv"))
d[, gate := factor(gate, levels=c("SOX9neg","excluded","SOX9pos"))]
cols <- c(SOX9neg="#377EB8", excluded="grey75", SOX9pos="#E41A1C")

# 1. histogram of module score, colored by gate, with threshold lines
p1 <- ggplot(d, aes(module, fill=gate)) +
  geom_histogram(bins=60, colour="white", linewidth=0.05) +
  geom_vline(xintercept=c(NEG,POS), linetype="dashed", linewidth=0.4) +
  scale_fill_manual(values=cols,
     labels=c(sprintf("SOX9- (<=%.2f): %d",NEG,sum(d$gate=="SOX9neg")),
              sprintf("excluded: %d",sum(d$gate=="excluded")),
              sprintf("SOX9+ (>=%.2f): %d",POS,sum(d$gate=="SOX9pos")))) +
  labs(title="SOX9 program module score - P23 tumor cells (n=2416)",
       x="module score", y="cells", fill="gate") +
  theme_classic(base_size=11) + theme(legend.position=c(0.78,0.8))

# 2. density split by ORIGINAL cluster (shows overlap the gate corrects)
p2 <- ggplot(d, aes(module, fill=SOX9_cluster, colour=SOX9_cluster)) +
  geom_density(alpha=0.35, linewidth=0.6) +
  geom_vline(xintercept=c(NEG,POS), linetype="dashed", linewidth=0.4) +
  scale_fill_manual(values=c(SOX9_high_P23="#E41A1C",SOX9_low_P23="#377EB8")) +
  scale_colour_manual(values=c(SOX9_high_P23="#E41A1C",SOX9_low_P23="#377EB8")) +
  labs(title="Module by original cluster (C9 high vs C4 low)",
       x="module score", y="density", fill="cluster", colour="cluster") +
  theme_classic(base_size=11) + theme(legend.position=c(0.78,0.8))

# 3. gate composition: how many of each gate come from C9 vs C4
comp <- d[, .N, by=.(gate, SOX9_cluster)]
p3 <- ggplot(comp, aes(gate, N, fill=SOX9_cluster)) +
  geom_col(position="stack", width=0.7) +
  scale_fill_manual(values=c(SOX9_high_P23="#E41A1C",SOX9_low_P23="#377EB8")) +
  labs(title="Gate composition by original cluster", x=NULL, y="cells", fill="cluster") +
  theme_classic(base_size=11)

PDF <- file.path(OUT,"Plots","SOX9_module_gate_distribution.pdf")
pdf(PDF, width=8, height=9)
print(p1 / p2 / p3 + plot_layout(heights=c(1,1,0.9)))
dev.off()
cat("Wrote", PDF, "\n")
print(d[, .N, by=gate][order(gate)])
