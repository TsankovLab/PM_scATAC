###############################################################################
# R2_Q14 -- volcano of BAP1-associated TF activity by compartment, linked to the malignant
#           result, with TCGA inferred-TF-activity agreement marked.
#
# One panel per compartment, shared axes:
#   Malignant   malignant volcano; hits coloured and labelled
#   Myeloid / Stroma / T-NK / B-plasma
#               faint malignant volcano as background; for each TF that is a hit in that
#               compartment an arrow runs from the TF's malignant point to its point in the
#               compartment (motif P < 0.01, gene score detected and concordant).
#               Dashed arrows: compartment hit count not above chance (permutation FDR > 0.25).
# x = BAP1 effect on motif activity (lost - retained, histology-adjusted); y = -log10 P.
# TCGA mark (black ring): the TF is in the TCGA MESO inferred-TF-activity table (Hmeljak et
#   al. 2018 tab 2B), moves in the SAME direction there (estimate > 0 = higher in BAP1-
#   inactivated) and has FDR < 0.05 there.  Only 141 TFs are in that table.
# Colours: project palette_celltype_lv1 (merged compartments take their largest member).
# Output: Plots/BAP1_TF_compartment_volcano.pdf, BAP1_TF_compartment_volcano_marks.csv
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(RColorBrewer); library(paletteer); library(circlize) })
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R2_Q14")); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
source(file.path(ROOT, "git_repo", "utils", "palettes.R"))
COMPS <- c("Malignant", "Myeloid", "Stroma", "TNK", "B_Plasma")
COL <- c(Malignant = palette_celltype_lv1[["Malignant"]], Myeloid = palette_celltype_lv1[["Myeloid"]],
         Stroma = palette_celltype_lv1[["Fibroblasts"]], TNK = palette_celltype_lv1[["T_cells"]],
         B_Plasma = palette_celltype_lv1[["B_cells"]])
NAME <- c(Malignant = "Malignant", Myeloid = "Myeloid", Stroma = "Stroma", TNK = "T/NK", B_Plasma = "B/plasma")
TCGA_FDR <- 0.05

A   <- read.csv("BAP1_TF_by_compartment_all.csv", stringsAsFactors = FALSE)
FDR <- read.csv("BAP1_TF_by_compartment_permFDR.csv", stringsAsFactors = FALSE)
P2B <- read.csv("affinity_regression/prepared/paper_2B.csv", stringsAsFactors = FALSE)
fdr <- setNames(FDR$empirical_FDR, FDR$compartment)
A$y <- -log10(A$motif_P_adj); A$x <- A$motif_logFC_adj
A$paper_est <- P2B$estimate_tf[match(A$symbol, P2B$TF)]; A$paper_fdr <- P2B$p_adj_tf[match(A$symbol, P2B$TF)]
A$tcga_same <- !is.na(A$paper_est) & A$paper_fdr < TCGA_FDR & sign(A$paper_est) == sign(A$x)

M <- A[A$compartment == "Malignant", ]
panel_lab <- sapply(COMPS, function(cp) {
  n  <- if (cp == "Malignant") sum(M$hit) else sum(A$hit & A$compartment == cp)
  nm <- sum(A$hit & A$compartment == cp & A$tcga_same)
  sprintf("%s: %d hits%s, %d TCGA-concordant", NAME[[cp]], n,
          if (fdr[[cp]] > 0.25) " (not above chance)" else "", nm) })

BG  <- do.call(rbind, lapply(COMPS, function(cp) data.frame(panel = cp, x = M$x, y = M$y)))
PTS <- A[A$hit, c("compartment", "motif", "x", "y", "tcga_same", "paper_est", "paper_fdr")]
PTS$panel <- PTS$compartment
SEG <- merge(PTS[PTS$compartment != "Malignant", ], M[, c("motif", "x", "y")], by = "motif", suffixes = c("", "_mal"))
SEG$dashed <- fdr[SEG$compartment] > 0.25
lvl <- setNames(panel_lab, COMPS)
for (nm in c("BG", "PTS", "SEG")) { d <- get(nm); d$panel <- factor(lvl[d$panel], levels = lvl); assign(nm, d) }

write.csv(A[A$hit, c("compartment", "motif", "symbol", "x", "motif_P_adj", "paper_est", "paper_fdr", "tcga_same")],
          "BAP1_TF_compartment_volcano_marks.csv", row.names = FALSE)
cat("TCGA-concordant hits (same direction, tab 2B FDR <", TCGA_FDR, "):\n")
print(A[A$hit & A$tcga_same, c("compartment", "motif", "x", "motif_P_adj", "paper_est", "paper_fdr")], row.names = FALSE, digits = 3)
cat("\nhits present in tab 2B but NOT concordant:\n")
print(A[A$hit & !is.na(A$paper_est) & !A$tcga_same, c("compartment", "motif", "x", "paper_est", "paper_fdr")], row.names = FALSE, digits = 3)

lim <- max(abs(PTS$x), abs(M$x)) * 1.05
p <- ggplot() +
  geom_hline(yintercept = 2, linetype = 3, linewidth = .3, colour = "grey60") +
  geom_vline(xintercept = 0, linewidth = .25, colour = "grey85") +
  geom_point(data = BG, aes(x, y), colour = "grey80", size = .6) +
  geom_segment(data = SEG, aes(x = x_mal, y = y_mal, xend = x, yend = y, colour = compartment, linetype = dashed),
               linewidth = .4, alpha = .75, arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  geom_point(data = SEG, aes(x_mal, y_mal), colour = COL[["Malignant"]], size = 1.1, alpha = .8) +
  geom_point(data = PTS, aes(x, y, colour = compartment), size = 2.4) +
  geom_point(data = subset(PTS, tcga_same), aes(x, y), shape = 21, colour = "black", fill = NA, size = 4, stroke = .8) +
  geom_text_repel(data = PTS, aes(x, y, label = motif, colour = compartment,
                                  fontface = ifelse(tcga_same, "bold", "plain")),
                  size = 2.7, max.overlaps = Inf, box.padding = .35, point.padding = .25, min.segment.length = 0,
                  segment.size = .2, segment.colour = "grey40", seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
  facet_wrap(~ panel, ncol = 3) +
  scale_x_continuous(limits = c(-lim, lim)) +
  gtheme_no_rot +
  theme(strip.text = element_text(size = 9, face = "bold"), strip.background = element_blank(),
        panel.border = element_rect(colour = "grey75", fill = NA, linewidth = .3),
        plot.title = element_text(size = 12), plot.subtitle = element_text(size = 9)) +
  xlab("BAP1 effect on motif activity, lost - retained (histology-adjusted)") + ylab(expression(-log[10]~italic(P))) +
  labs(title = "BAP1-associated TF activity by compartment",
       subtitle = paste0("Grey = all malignant motifs. Arrows run from a TF's malignant result to its hit in that compartment (motif P < 0.01, gene score concordant).\n",
                         "Black ring + bold label = same direction and FDR < ", TCGA_FDR,
                         " in the TCGA MESO inferred TF-activity table (tab 2B). Dashed = compartment hit count not above chance. Dotted line P = 0.01."))
ggsave("Plots/BAP1_TF_compartment_volcano.pdf", p, width = 15, height = 10)
cat("\nDONE\n")
