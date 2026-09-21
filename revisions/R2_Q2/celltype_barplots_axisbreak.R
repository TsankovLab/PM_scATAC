# celltype_barplots_axisbreak.R  (R2_Q2)
#
# Same two-panel cell-type composition figure as celltype_barplots.R, but the count
# panel now has a BROKEN COUNT AXIS.
#
# Why: P23 contributes 18,228 of the 49,849 cells, roughly 3.6x the next largest
# sample (P11, 5,052) and 96x the smallest (P3, 189).  On a single linear axis the
# other ten samples are squeezed into the left quarter of the panel and their
# composition cannot be read.  The break keeps the axis linear and honest within
# each segment -- unlike a log or sqrt axis, a stacked bar stays proportional, so
# segment lengths remain directly comparable.
#
# The break is drawn by plotting the same bars twice with different coord_cartesian
# x-limits and gluing the two panels with patchwork: a wide low-range panel
# (0 - BREAK_LO) and a narrow high-range panel (BREAK_HI - max) holding only the end
# of the P23 bar.  Nothing is dropped or rescaled; the middle of the axis is simply
# not drawn, marked with the usual pair of diagonal glyphs.
#
# Cell counts are cached to celltype_counts.csv on the first run so the figure can
# be re-styled without reloading the ArchR project.
#
# Output: celltype_barplots_axisbreak.pdf  (width FIG_W in, height FIG_W/2)
# Run via submit_celltype_barplots.sh (edit the script name) or ./run_break.sh

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(paletteer)
  library(circlize)
  library(RColorBrewer)
})

script_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q2"
utils_dir  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils"
main_dir   <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR"
setwd(script_dir)

source(file.path(utils_dir, "ggplot_aestetics.R"))
source(file.path(utils_dir, "useful_functions.R"))
source(file.path(utils_dir, "palettes.R"))

## ---- figure size -------------------------------------------------------------
## gtheme's point sizes are absolute, so shrinking the canvas alone overflows the
## legend and clips the n= labels.  Everything type-related is therefore derived
## from FIG_W: set the width and the rest follows, so the figure can be re-scaled
## from one place.  Reference is the 12 in canvas the original figure used.
FIG_W <- 4.2
FIG_H <- FIG_W / 2
SC    <- FIG_W / 12                     # scale relative to the original canvas
compact <- theme(
  axis.text          = element_text(size = 13 * SC),
  axis.title         = element_text(size = 17 * SC),
  legend.title       = element_text(size = 16 * SC),
  legend.text        = element_text(size = 13.5 * SC),
  legend.key.size    = unit(6.2 * SC, "mm"),
  legend.margin      = margin(0, 0, 0, 0),
  legend.box.spacing = unit(4 * SC, "mm"),
  axis.ticks         = element_line(linewidth = 0.8 * SC),
  panel.border       = element_blank())
LAB_SIZE <- 4.2 * SC                    # the n= annotations

CACHE <- "celltype_counts.csv"
sample_order <- rev(c("P1","P3","P4","P5","P8","P10","P11","P12","P13","P14","P23"))

## ---- counts (cached) ---------------------------------------------------------
if (file.exists(CACHE)) {
  counts_df <- read.csv(CACHE, stringsAsFactors = FALSE)
  message("counts from cache")
} else {
  suppressPackageStartupMessages(library(ArchR))
  addArchRGenome("hg38"); addArchRThreads(threads = 2)
  message("Loading main ArchR project...")
  archp <- loadArchRProject(main_dir, showLogo = FALSE)
  archp_meta <- as.data.frame(archp@cellColData)
  archp_meta$celltype_lv1 <- as.character(archp_meta$celltype_lv1)
  archp_meta$celltype_lv1[archp_meta$celltype_lv1 == "SmoothMuscle"] <- "Smooth Muscle"
  archp_meta$celltype_lv1 <- factor(archp_meta$celltype_lv1,
                                    levels = names(palette_celltype_lv1))
  archp_meta$Sample2 <- factor(archp_meta$Sample, levels = sample_order)
  counts_df <- cellComp(seurat_obj = archp_meta,
                        metaGroups = c("Sample2", "celltype_lv1"),
                        plot_as = "bar", pal = palette_celltype_lv1,
                        prop = FALSE, returnDF = TRUE)
  write.csv(counts_df, CACHE, row.names = FALSE)
}
counts_df$Sample2      <- factor(counts_df$Sample2, levels = sample_order)
counts_df$celltype_lv1 <- factor(counts_df$celltype_lv1, levels = names(palette_celltype_lv1))

totals_df <- aggregate(Freq ~ Sample2, data = counts_df, FUN = sum)
colnames(totals_df)[2] <- "total"
totals_df$Sample2 <- factor(totals_df$Sample2, levels = sample_order)

## ---- where to break ----------------------------------------------------------
## Low segment ends just above the largest sample that is NOT P23, so every other
## bar is drawn complete; the high segment is a narrow window around P23's total.
big      <- max(totals_df$total)                                  # P23
next_big <- max(totals_df$total[totals_df$total < big])           # P11
BREAK_LO <- ceiling(next_big * 1.20 / 500) * 500                  # 6100
BREAK_HI <- floor(big * 0.94 / 500) * 500                         # 17000
XMAX     <- ceiling(big * 1.07 / 500) * 500                       # 19500
W        <- c(3.6, 1)                                             # panel widths
cat(sprintf("break: low 0-%d, high %d-%d (largest %d = P23, next %d)\n",
            BREAK_LO, BREAK_HI, XMAX, big, next_big))

## ---- panel 1: proportions (unchanged) ----------------------------------------
prop_df <- counts_df
prop_df$prop <- ave(prop_df$Freq, prop_df$Sample2, FUN = function(x) x / sum(x))
bp_prop <- ggplot(prop_df, aes(x = prop, y = Sample2, fill = celltype_lv1)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = palette_celltype_lv1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  gtheme + compact +
  labs(x = "Proportion", y = NULL, fill = "Cell type") +
  theme(legend.position = "none")

## ---- panel 2: counts, broken axis --------------------------------------------
## The same layer stack twice, differing only in coord_cartesian(xlim).  Limits are
## set on the COORD, not the scale, so bars that leave the window are clipped rather
## than dropped and the P23 bar reads as continuing through the break.
bars <- function(xlim, labs_df, axis_y, brks, xlab) {
  p <- ggplot(counts_df, aes(x = Freq, y = Sample2, fill = celltype_lv1)) +
    geom_col(position = "stack", width = 0.8) +
    scale_fill_manual(values = palette_celltype_lv1) +
    scale_x_continuous(expand = c(0, 0), breaks = brks) +
    ## clip MUST stay on: with clip = "off" the P23 bar (18,228) is drawn past the
    ## end of the low panel and overflows across the rest of the figure.
    coord_cartesian(xlim = xlim, clip = "on") +
    gtheme + compact + labs(x = xlab, y = NULL, fill = "Cell type")
  if (nrow(labs_df))
    p <- p + geom_text(data = labs_df, aes(x = total, y = Sample2,
                                           label = paste0("n=", total)),
                       hjust = -0.15, size = LAB_SIZE, inherit.aes = FALSE)
  if (!axis_y)
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
}

## the n= label goes in whichever segment its bar ends in
lab_lo <- totals_df[totals_df$total <= BREAK_LO, ]
lab_hi <- totals_df[totals_df$total >  BREAK_LO, ]

## neither count panel repeats the sample names -- the proportion panel carries them
p_lo <- bars(c(0, BREAK_LO * 1.02), lab_lo, axis_y = FALSE,
             brks = seq(0, BREAK_LO, by = 2000), xlab = "Cell count") +
  theme(plot.margin = margin(2, 1, 2, 3) * SC * 3, legend.position = "none")
p_hi <- bars(c(BREAK_HI, XMAX * 1.14), lab_hi, axis_y = FALSE,
             brks = c(17500, 19000), xlab = NULL) +
  theme(plot.margin = margin(2, 2, 2, 1) * SC * 3, legend.position = "right",
        axis.title.x = element_blank())

## break glyphs: a pair of short diagonals at each inner edge.  They sit just INSIDE
## the panel (the discrete axis expands to 0.4, so 0.45-0.75 is the lowest strip of
## drawable space) because clipping has to stay on -- see above.
gl <- function(p, at, dx, side) {
  off <- c(-dx, dx)
  p + annotate("segment",
               x = at + off - dx * 0.9, xend = at + off + dx * 0.9,
               y = 0.42, yend = 0.95, linewidth = 0.9 * SC, colour = "grey20")
}
p_lo <- gl(p_lo, BREAK_LO * 0.995, dx = BREAK_LO * 0.016, side = "lo")
p_hi <- gl(p_hi, BREAK_HI * 1.004, dx = (XMAX - BREAK_HI) * 0.016, side = "hi")

## ---- combine and save ---------------------------------------------------------
## One patchwork, three panels, guides collected once.  The two count panels carry
## the widths that make the break: the high segment gets a quarter of the space of
## the low one, which is what un-squeezes the ten small samples.
combined <- bp_prop + p_lo + p_hi +
  plot_layout(widths = c(4.4, 3.4, 1.2), guides = "collect") &
  theme(legend.box.margin = margin(0, 0, 0, 4))

pdf(file.path(script_dir, "celltype_barplots_axisbreak.pdf"),
    height = FIG_H, width = FIG_W)
print(combined)
dev.off()

message("Done. Output: celltype_barplots_axisbreak.pdf")
