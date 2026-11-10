###############################################################################
# STEP 9 -- the joint embedding figures.
#
#   A  joint UMAP coloured by ASSAY.  The honest reading of this panel is limited:
#      imputation pulls the ATAC cells toward the RNA reference by construction, so
#      overlap is expected and is not by itself evidence.  It is here to show that
#      no cell type forms an assay-only island.
#   B  the same embedding coloured by CELL TYPE -- chromatin label for the ATAC
#      cells, expression label for the RNA cells, assigned independently.  If the
#      two annotations name the same biology, the colours must agree across the two
#      assays without ever having been matched.
#   composition  for each scATAC cell type, the expression label of the nearest
#      scRNA cell in the joint UMAP (1-NN), as the fraction of its cells.
#
# Input : joint_embedding.csv (step 8)
# Output: Plots/R2Q5_joint_panel_A.pdf, Plots/R2Q5_joint_panel_B.pdf,
#         Plots/scatac_celltype_comp_in_scrna.pdf
###############################################################################
suppressMessages({ library(ggplot2); library(RColorBrewer); library(FNN); library(dplyr); library(scales) })
source("00_common.R")
source(file.path(SC, "git_repo", "utils", "palettes.R"))   # palette_celltype_lv1

th <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        legend.key.size = unit(3, "mm"))
thu <- th + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                  axis.title = element_blank(), plot.title = element_text(size = 8))

E  <- read.csv("joint_embedding.csv", stringsAsFactors = FALSE)
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

for (n in c("pA", "pB"))
  ggsave(sprintf("Plots/R2Q5_joint_panel_%s.pdf", sub("^p", "", n)), get(n),
         width = 5.5, height = 4.5, device = cairo_pdf)

## ---- scATAC cell types by the label of the nearest scRNA cell ----------------
# Split modalities
rna  <- E %>% filter(modality == "scRNA")
atac <- E %>% filter(modality == "scATAC")

# Find nearest scRNA cell for each scATAC cell based on UMAP1/UMAP2
nn <- get.knnx(
  data  = as.matrix(rna[, c("UMAP1", "UMAP2")]),
  query = as.matrix(atac[, c("UMAP1", "UMAP2")]),
  k = 1
)

# Assign nearest scRNA cell type
atac$closest_scRNA_celltype <- rna$celltype[nn$nn.index[, 1]]

# Composition
plot_df <- atac %>%
  count(celltype, closest_scRNA_celltype) %>%
  group_by(celltype) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Order bars by purity
purity_order <- plot_df %>%
  group_by(celltype) %>%
  summarise(purity = max(prop), .groups = "drop") %>%
  arrange(desc(purity)) %>%
  pull(celltype)

plot_df <- plot_df %>%
  mutate(celltype = factor(celltype, levels = purity_order)) %>%
  
  # Within EACH bar: largest -> smallest
  group_by(celltype) %>%
  arrange(desc(prop), .by_group = TRUE) %>%
  mutate(
    ymax = cumsum(prop),
    ymin = lag(ymax, default = 0)
  ) %>%
  ungroup() %>%
  
  # Numeric x position for geom_rect
  mutate(x = as.numeric(celltype))

# Plot
sp = ggplot(plot_df) +
  geom_rect(
    aes(
      xmin = x - 0.4,
      xmax = x + 0.4,
      ymin = ymin,
      ymax = ymax,
      fill = closest_scRNA_celltype
    )
  ) +
  scale_fill_manual(values = palette_celltype_lv1) +
  scale_x_continuous(
    breaks = seq_along(purity_order),
    labels = purity_order
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "scATAC cell type",
    y = "Proportion",
    fill = "Closest scRNA cell type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

pdf (file.path ('Plots','scatac_celltype_comp_in_scrna.pdf'), height=3)
print (sp)
dev.off()
cat("DONE\n")
