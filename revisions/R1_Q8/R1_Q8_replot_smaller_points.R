###############################################################
# R1_Q8 : Re-plot the deviation-vs-RNA scatter with error bars,
#         using SMALLER points (matching the original figure size).
#         Reads the saved TF_diff_rna table -> no heavy recompute.
###############################################################

GITROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo"
OUTDIR  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q8"
setwd(OUTDIR)

library(ggplot2)
library(ggrepel)
source(file.path(GITROOT, "utils", "ggplot_aestetics.R"))   # gtheme_no_rot

TF_diff_rna <- read.csv(file.path(OUTDIR, "TF_diff_rna_with_errorbars.csv"), row.names = 1,
                        check.names = FALSE, stringsAsFactors = FALSE)
# restore empty-string labels (read.csv turns "" into NA)
TF_diff_rna$label[is.na(TF_diff_rna$label)]         <- ''
TF_diff_rna$label_top[is.na(TF_diff_rna$label_top)] <- ''

diff_line <- 0
# Error bars only for the LABELLED genes (those carrying a text label).
eb_df <- TF_diff_rna[TF_diff_rna$label_top != '', ]

# Dot size encodes RNA expression level (mean log2 tumor expression),
# kept in a small range so points stay compact like the original figure.
size_range <- c(0.8, 3.5)

base_plot <- function(text_size = 1.8, seg = .05) {
  ggplot(TF_diff_rna,
      aes(x = dev_diff, y = rna_diff, size = tumor_rna, fill = direction, label = label_top)) +
    geom_errorbar(data = eb_df, inherit.aes = FALSE,
                  aes(x = dev_diff, ymin = rna_diff - rna_sem, ymax = rna_diff + rna_sem),
                  width = 0, linewidth = .2, color = "grey55", alpha = .8) +
    geom_errorbarh(data = eb_df, inherit.aes = FALSE,
                   aes(y = rna_diff, xmin = dev_diff - dev_sem, xmax = dev_diff + dev_sem),
                   height = 0, linewidth = .2, color = "grey55", alpha = .8) +
    geom_point(aes(color = color, alpha = alpha), shape = 21, stroke = 0.1) +
    scale_size_continuous(range = size_range, name = "tumor RNA\n(log2 mean)") +
    scale_color_manual(values = c('FALSE' = "grey", 'TRUE' = "black")) +
    scale_fill_manual(values = c(both_low = "grey20", rna_low = "grey20",
                                 dev_low = "grey20", both_high = "darkred")) +
    geom_vline(xintercept = diff_line, linetype = "dashed", color = "grey44", linewidth = .3) +
    geom_hline(yintercept = diff_line, linetype = "dashed", color = "grey44", linewidth = .3) +
    gtheme_no_rot +
    geom_text_repel(aes(color = label_guides, fill = label_guides),
                    segment.size = seg, size = text_size, max.overlaps = 100) +
    xlab("activity difference") + ylab("RNA difference") +
    # coord_cartesian zooms WITHOUT dropping out-of-range error bars
    # (xlim()/ylim() would delete whole bars that spill past the limits,
    #  which is why e.g. CUX1's activity bar was missing before).
    coord_cartesian(xlim = c(-0.2, 0.2), ylim = c(-0.6, 0.6))
}

# 5x4 version (original dimensions)
pdf(file.path("Plots", "Diff_normal_tumor_deviation_and_rna_scatterplot5_errorbars.pdf"), width = 5, height = 4)
print(base_plot(text_size = 1.8, seg = .05))
dev.off()

# larger readable version, still with small points
pdf(file.path("Plots", "Diff_normal_tumor_deviation_and_rna_scatterplot5_errorbars_large.pdf"), width = 8, height = 6.5)
print(base_plot(text_size = 3, seg = .1) + theme(text = element_text(size = 12)))
dev.off()

cat("DONE\n")
