###############################################################################
# STEP 9 -- the joint embedding figure.
#
#   A  joint UMAP coloured by ASSAY.  The honest reading of this panel is limited:
#      imputation pulls the ATAC cells toward the RNA reference by construction, so
#      overlap is expected and is not by itself evidence.  It is here to show that
#      no cell type forms an assay-only island.
#   B  the same embedding coloured by CELL TYPE -- chromatin label for the ATAC
#      cells, expression label for the RNA cells, assigned independently.  If the
#      two annotations name the same biology, the colours must agree across the two
#      assays without ever having been matched.
#   C  the same, split by assay, so the two annotations can be compared side by side
#      on identical coordinates.
#   D  modality mixing per cell type: fraction of each cell's 30 nearest joint-space
#      neighbours drawn from the other assay, against what the cohort composition
#      alone would give.
#   E  centroid geometry: distance between the two assays' centroids for a cell type,
#      over the median distance to other cell types.  Below 1 = organised by cell
#      type, not by assay.
#   F  why D and E disagree: the spread of the imputed ATAC cells around their own
#      centroid, relative to the real RNA cells.  TransferData imputes an
#      anchor-weighted average, so the ATAC cloud is roughly half as wide as the RNA
#      one and each ATAC cell's nearest neighbours are other ATAC cells even though
#      the cloud sits in the right place.  D is a density artefact; E is the result.
#
# Input : joint_embedding.csv, joint_mixing_by_celltype.csv,
#         joint_centroid_distance.csv, joint_dispersion.csv
# Output: Plots/R2Q5_joint_embedding.pdf and the panels separately
###############################################################################
source("00_common.R")
suppressMessages({ library(ggplot2); library(patchwork) })

th <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        legend.key.size = unit(3, "mm"))
thu <- th + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                  axis.title = element_blank(), plot.title = element_text(size = 8))

E  <- read.csv("joint_embedding.csv", stringsAsFactors = FALSE)
MX <- read.csv("joint_mixing_by_celltype.csv", stringsAsFactors = FALSE)
CD <- read.csv("joint_centroid_distance.csv", stringsAsFactors = FALSE)
SP <- read.csv("joint_dispersion.csv", stringsAsFactors = FALSE)
lv <- CELLTYPES[CELLTYPES %in% E$celltype]
set.seed(1); E <- E[sample(nrow(E)), ]

MODCOL <- c(scATAC = "#3b6ea5", scRNA = "#e08214")

pA <- ggplot(E, aes(UMAP1, UMAP2, colour = modality)) +
  geom_point(size = .05, stroke = 0, alpha = .5) +
  scale_colour_manual(values = MODCOL, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  labs(title = sprintf("joint CCA embedding: %s scATAC + %s scRNA cells",
                       format(sum(E$modality == "scATAC"), big.mark = ","),
                       format(sum(E$modality == "scRNA"), big.mark = ","))) + thu

pB <- ggplot(E, aes(UMAP1, UMAP2, colour = factor(celltype, lv))) +
  geom_point(size = .05, stroke = 0, alpha = .5) +
  scale_colour_manual(values = CTCOL, drop = FALSE, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1), ncol = 1)) +
  labs(title = "cell type — annotated independently in each assay") + thu

pC <- ggplot(E, aes(UMAP1, UMAP2, colour = factor(celltype, lv))) +
  geom_point(size = .05, stroke = 0, alpha = .5) +
  facet_wrap(~ modality) +
  scale_colour_manual(values = CTCOL, drop = FALSE, guide = "none") +
  labs(title = "the two annotations on identical coordinates") + thu

MX$celltype <- factor(MX$celltype, rev(lv))
pD <- ggplot(MX, aes(mixing_ratio, celltype, fill = modality)) +
  geom_col(position = position_dodge(width = .75), width = .7) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = .3, colour = "grey30") +
  scale_fill_manual(values = MODCOL, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, .05))) +
  labs(x = "neighbours from the other assay / expected from composition", y = NULL,
       title = "modality mixing — scATAC bars are ~0 (see F), not missing") + th +
  theme(plot.title = element_text(size = 8))

CD <- CD[order(CD$ratio), ]
pE <- ggplot(CD, aes(ratio, factor(celltype, CD$celltype), fill = celltype)) +
  geom_col(width = .7) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = .3, colour = "grey30") +
  scale_fill_manual(values = CTCOL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, .05))) +
  labs(x = "assay-to-assay distance / median distance to other cell types", y = NULL,
       title = "same cell type across assays vs across cell types") + th +
  theme(plot.title = element_text(size = 8))

SP <- SP[order(SP$shrinkage), ]
pF <- ggplot(SP, aes(shrinkage, factor(celltype, SP$celltype), fill = celltype)) +
  geom_col(width = .7) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = .3, colour = "grey30") +
  scale_fill_manual(values = CTCOL, guide = "none") +
  scale_x_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, .02))) +
  labs(x = "spread of imputed scATAC cells / spread of real scRNA cells", y = NULL,
       title = sprintf("imputation shrinks the ATAC cloud (median %.2f)",
                       median(SP$shrinkage))) + th +
  theme(plot.title = element_text(size = 8))

fig <- (pA | pB) / pC / (pD | pE) / (pF | plot_spacer()) +
  plot_annotation(tag_levels = "A",
    title = "R2_Q5  Joint CCA embedding of the scATAC and scRNA cells",
    subtitle = "CCA anchors -> RNA expression imputed for every ATAC cell -> shared PCA/UMAP. Cell type labels were assigned independently in each assay and never matched.",
    theme = theme(plot.title = element_text(size = 10, face = "bold"),
                  plot.subtitle = element_text(size = 7.5, colour = "grey30")))
ggsave("Plots/R2Q5_joint_embedding.pdf", fig, width = 11, height = 16,
       device = cairo_pdf, limitsize = FALSE)
for (n in c("pA","pB","pC","pD","pE","pF"))
  ggsave(sprintf("Plots/R2Q5_joint_panel_%s.pdf", sub("^p", "", n)), get(n),
         width = 5.5, height = 4.5, device = cairo_pdf)
cat("wrote Plots/R2Q5_joint_embedding.pdf and 6 panels\nDONE\n")
