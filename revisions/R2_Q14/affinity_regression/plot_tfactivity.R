## Figures for one affinity-regression result directory (results/<tag>):
##   compare_estimates.pdf   our BAP1 effect on inferred TF activity vs tab 2B
##   volcano_ours_vs_paper.pdf  Fig 2E-style volcano, ours and the paper side by side
##   boxplots_YY1_IRF8.pdf   Fig 2F/G equivalents
args <- commandArgs(trailingOnly = TRUE); dir <- args[1]
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(patchwork) })
source("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils/ggplot_aestetics.R")
R <- read.csv(file.path(dir, "bap1_tf_ttest.csv"), stringsAsFactors = FALSE)
R <- R[!is.na(R$paper_estimate), ]
rho <- cor(R$estimate, R$paper_estimate, method = "spearman")
R$lab <- ifelse(R$paper_padj < 0.01 & R$padj < 0.01 | R$TF %in% c("IRF8", "EGR2", "YY1"), R$TF, NA)
p1 <- ggplot(R, aes(paper_estimate, estimate)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = .2) + geom_vline(xintercept = 0, colour = "grey85", linewidth = .2) +
  geom_point(colour = "grey70", size = 1.2) +
  geom_point(data = subset(R, paper_padj < 0.01), colour = "grey25", size = 1.6) +
  geom_text_repel(aes(label = lab), size = 2.3, max.overlaps = 40, seed = 1, min.segment.length = 0) +
  gtheme_no_rot + xlab("paper (tab 2B): BAP1 effect on inferred TF activity") +
  ylab("reproduction: BAP1 effect on inferred TF activity") +
  ggtitle(sprintf("BAP1 effect on inferred TF activity: %d TFs, Spearman rho = %.2f", nrow(R), rho),
          subtitle = "dark points = paper FDR < 0.01") +
  theme(plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
ggsave(file.path(dir, "compare_estimates.pdf"), p1, width = 6, height = 5)

V <- rbind(data.frame(source = "reproduction", TF = R$TF, x = R$estimate, fdr = R$padj),
           data.frame(source = "paper (tab 2B)", TF = R$TF, x = R$paper_estimate, fdr = R$paper_padj))
V$sig <- V$fdr < 0.01
focus <- c("IRF8", "EGR2", "YY1")
both  <- R$TF[R$padj < 0.01 & R$paper_padj < 0.01]          # significant in BOTH analyses
V$lab <- ifelse((V$sig & V$TF %in% both) | V$TF %in% focus, V$TF, NA)
p2 <- ggplot(V, aes(x, -log10(fdr))) +
  geom_hline(yintercept = 2, linetype = 2, linewidth = .3, colour = "grey55") +
  geom_point(aes(colour = sig), size = 1.2) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#e67e22"), guide = "none") +
  geom_text_repel(aes(label = lab), size = 2.1, max.overlaps = 40, seed = 1, min.segment.length = 0) +
  facet_wrap(~ source, scales = "free") + gtheme_no_rot +
  xlab("mean inferred TF activity difference, BAP1 inactivated - wild-type") + ylab(expression(-log[10]~FDR)) +
  ggtitle(sprintf("Fig 2E: orange = FDR < 0.01; labelled = significant in both (%d) and IRF8/EGR2/YY1", length(both)))
ggsave(file.path(dir, "volcano_ours_vs_paper.pdf"), p2, width = 9, height = 4.4)

A <- read.csv(file.path(dir, "tf_activity.csv"), check.names = FALSE, stringsAsFactors = FALSE)
S <- read.csv(file.path(dir, "samples.csv"), stringsAsFactors = FALSE)
L <- do.call(rbind, lapply(c("YY1", "IRF8"), function(g) data.frame(TF = g, sample = names(A)[-1],
      activity = as.numeric(A[A$TF == g, -1]))))
L$status <- S$bap1_status[match(L$sample, S$sample)]
L <- L[L$status %in% c("inactivated", "no_inactivation"), ]
L$status <- factor(ifelse(L$status == "inactivated", "BAP1 inactivated", "BAP1 wild-type"), c("BAP1 wild-type", "BAP1 inactivated"))
pv <- sapply(c("YY1", "IRF8"), function(g) R$p[R$TF == g])
L$TF <- factor(sprintf("%s (t-test p = %.2g)", L$TF, pv[L$TF]), levels = sprintf("%s (t-test p = %.2g)", c("YY1","IRF8"), pv))
p3 <- ggplot(L, aes(status, activity)) + geom_boxplot(outlier.shape = NA, width = .55, fill = "grey92", linewidth = .3) +
  geom_jitter(width = .15, size = .9, alpha = .7) + facet_wrap(~ TF, scales = "free_y") +
  gtheme_no_rot + xlab(NULL) + ylab("inferred TF activity") + ggtitle("Fig 2F/G equivalents")
ggsave(file.path(dir, "boxplots_YY1_IRF8.pdf"), p3, width = 5.6, height = 3.6)
cat("figures written to", dir, "\n")
