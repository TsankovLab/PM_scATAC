###############################################################################
# STEP 5 -- one circular tree over all 18 subclones.
#
# Leaf   = one subclone.  Distance = 1 - Pearson between clone CNV profiles over the
#          shared 5 Mb bins; ward.D2.  Clones of the same tumour normally sit together,
#          so a leaf that does NOT is a clone whose CNV profile is genuinely distinct.
# Label  = sample / clone / the arm change that separates it from its sibling.
# Ring   = P4's chr8q-amplified clone, the one the Visium data localise in step 9.
#
# Input : epi_clone_profiles.rds (step 3)
# Output: Plots/epiAneufinder_circular_tree.pdf
###############################################################################
suppressMessages({ library(ape); library(ggtree); library(ggplot2) })
source("00_common.R")
P <- readRDS("epi_clone_profiles.rds")
Z <- P$Z; meta <- P$meta

D <- as.dist(1 - cor(t(Z))); phy <- as.phylo(hclust(D, "ward.D2"))
dat <- data.frame(label = meta$leaf, sample = meta$sample, n_cells = meta$n_cells,
                  tip = sprintf("%s  %s", meta$leaf, meta$driver), stringsAsFactors = FALSE)

p <- ggtree(phy, layout = "fan", open.angle = 14, size = 0.5, colour = "grey35") %<+% dat +
  geom_tippoint(aes(colour = sample, size = n_cells), stroke = 0) +
  geom_tiplab(aes(label = tip, colour = sample, angle = angle), size = 2.1,
              offset = 0.02, show.legend = FALSE) +
  scale_colour_manual(values = SAMPCOL, name = NULL) +
  scale_size_continuous(range = c(1.2, 4), trans = "sqrt", name = NULL,
                        breaks = c(100, 1000, 10000)) +
  labs(title = sprintf("epiAneufinder subclones, all scATAC tumours (%g Mb windows)", WINDOW/1e6),
       subtitle = paste0("leaf = sample / clone / dominant arm change | 1 - Pearson over ",
                         ncol(Z), " shared ", WINDOW/1e6, " Mb bins | green ring = P4 chr8q-amplified clone")) +
  theme_tree() +
  theme(plot.title = element_text(face = "bold", size = 10.5),
        plot.subtitle = element_text(size = 6.8, colour = "grey35"),
        legend.position = "inside", legend.position.inside = c(0.5, 0.5),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 6),
        legend.background = element_rect(fill = NA, colour = NA),
        plot.margin = margin(2, 2, 2, 2))

hl <- p$data[p$data$label == P$chr8q_leaf & p$data$isTip, ]
p <- p + geom_point(data = hl, aes(x = x, y = y), shape = 21, size = 5.2,
                    colour = COL_HI, fill = NA, stroke = 1.4) +
         xlim(-0.15, max(p$data$x, na.rm = TRUE) * 1.75)

ggsave("Plots/epiAneufinder_circular_tree.pdf", p, width = 7.6, height = 7.2, device = cairo_pdf)
cat("DONE -> Plots/epiAneufinder_circular_tree.pdf\n")
