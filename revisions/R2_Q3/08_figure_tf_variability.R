###############################################################################
# STEP 8 -- TF variability figure.
#
#  A  per tumour: the ARI between an independent chromVAR clustering and the CNV
#     clones. Near zero everywhere = the clones are not distinct regulatory states.
#  B  the most variable motifs, median |between-clone difference| +/- IQR over the
#     tumours. Immune/interferon and the AP-1/bZIP technical axis are colour-coded.
#
# Input : epi_chromvar_recurrence.csv, epi_chromvar_mirror.csv (step 7)
# Output: Plots/epi_TF_variability.pdf
###############################################################################
suppressMessages({ library(ggplot2); library(patchwork) })
source("00_common.R")
COL_TECH <- "#b0413e"; COL_IMM <- "#7b3fa0"; COL_OTH <- "grey55"

rec  <- read.csv("epi_chromvar_recurrence.csv", stringsAsFactors = FALSE)
metr <- read.csv("epi_chromvar_mirror.csv",     stringsAsFactors = FALSE)
rec$class <- ifelse(rec$tech, "AP-1/bZIP (technical)",
             ifelse(rec$immune, "immune / interferon", "other TF"))
pal <- c(`AP-1/bZIP (technical)` = COL_TECH, `immune / interferon` = COL_IMM,
         `other TF` = COL_OTH)

top <- head(rec[order(-rec$median_absdiff), ], 22)
top$TF <- factor(top$TF, levels = rev(top$TF))

pA <- ggplot(metr, aes(reorder(sample, ARI_chromvar_vs_CNV), ARI_chromvar_vs_CNV)) +
  geom_col(aes(fill = sample == "P4"), width = 0.68, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = COL_HI, `FALSE` = "grey65")) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) + coord_flip() +
  labs(x = NULL, y = "ARI: chromVAR clustering vs epiAneufinder clones",
       title = "A  Does TF activity recover the clone split?",
       subtitle = "0 = no correspondence") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"))

pB <- ggplot(top, aes(median_absdiff, TF, colour = class)) +
  geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 0.45) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = "median |between-clone difference| in chromVAR z   (bars = IQR over tumours)",
       y = NULL,
       title = "B  Most variable TF motifs between epiAneufinder subclones",
       subtitle = "median across samples, so a TF must move in most tumours, not one") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 6.5),
        axis.text.y = element_text(size = 6.6),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"))

ggsave("Plots/epi_TF_variability.pdf", (pA | pB) + plot_layout(widths = c(1, 1.35)),
       width = 10, height = 5, device = cairo_pdf)
cat("DONE -> Plots/epi_TF_variability.pdf\n")
