###############################################################################
# R1_Q8 -- make the metacell-parameter effect visible
#
# The 12 sweep plots are all genuinely different (12 distinct files) but look alike, because
# each shows 20 TFs x 2 readouts on a fixed -1..1 axis and the parameter shifts medians by
# only ~0.05-0.15. This condenses the same numbers into the comparison that is actually
# readable: what each parameter does to the score, and by how much.
#   A  median |score| vs width, one line per step   -> width matters, step does not
#   B  per-TF activity median across widths          -> which TFs move, and in which direction
#   C  the same TFs as boxplots at width 20 / 40 / 50, for a handful of TFs rather than all 20
# Reads the CSVs written by the two sweep scripts; nothing is recomputed.
#
# Output: Plots/binstep_comparison.pdf (+ .png)
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr) })
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R1_Q8"))
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))

G <- read.csv("sarc_TF_binstep_summary.csv")
G$step_lab <- ifelse(G$step == G$width, "no overlap (step = width)",
              ifelse(G$step == 1, "sliding (step = 1)", "50% overlap"))
pA <- ggplot(G, aes(width, median_abs, colour = step_lab, group = step_lab)) +
  geom_line(linewidth = .5) + geom_point(size = 2) +
  scale_colour_manual(values = c("no overlap (step = width)" = "#b2182b",
                                 "50% overlap" = "#ef8a62", "sliding (step = 1)" = "#2166ac"), name = NULL) +
  theme_classic(base_size = 9) +
  labs(x = "metacell width (cells)", y = "median |score| over all TFs",
       title = "A  Width drives the spread; step barely does",
       subtitle = "each point is one of the 12 sweep settings")

M <- read.csv("sarc_TF_binsweep_medians.csv")
A <- M[M$type == "activity", ]
ord <- A %>% filter(bin == 40) %>% arrange(desc(median)) %>% pull(TF)
A$TF <- factor(A$TF, levels = ord)
mv <- A %>% group_by(TF) %>% summarise(rng = diff(range(median))) %>% arrange(desc(rng))
A$hl <- ifelse(A$TF %in% head(mv$TF, 4), as.character(A$TF), "other")
pB <- ggplot(A, aes(factor(bin), median, group = TF)) +
  geom_line(aes(colour = hl, linewidth = hl != "other"), alpha = .8) +
  geom_point(aes(colour = hl), size = .9) +
  scale_colour_manual(values = c(setNames(c("#b2182b","#ef8a62","#2166ac","#66a61e"), head(mv$TF,4)),
                                 other = "grey80"), name = NULL) +
  scale_linewidth_manual(values = c(`FALSE` = .3, `TRUE` = .8), guide = "none") +
  theme_classic(base_size = 9) +
  labs(x = "metacell width (cells)", y = "median activity score",
       title = "B  Per-TF medians across widths",
       subtitle = sprintf("coloured = the 4 TFs that move most (range up to %.2f)", max(mv$rng)))

sel <- head(mv$TF, 6)
B <- M[M$TF %in% sel, ]
B$TF <- factor(B$TF, levels = sel)
pC <- ggplot(B, aes(factor(bin), median, fill = type)) +
  geom_col(position = position_dodge(.8), width = .7) +
  geom_hline(yintercept = 0, linewidth = .3) +
  facet_wrap(~ TF, nrow = 1) +
  scale_fill_manual(values = c(activity = "#AE123AFF", genescore = "#001260FF"), name = NULL) +
  theme_classic(base_size = 8) +
  labs(x = "metacell width (cells)", y = "median score",
       title = "C  The TFs most sensitive to the setting")
ggsave("Plots/binstep_comparison.pdf", pA / pB / pC, width = 9, height = 9)
ggsave("Plots/binstep_comparison.png", pA / pB / pC, width = 9, height = 9, dpi = 170)

cat("=== how much does each TF's activity median move across widths 20-50? ===\n")
print(head(as.data.frame(mv), 10), row.names = FALSE, digits = 3)
cat("\nmedian |score| by width (averaged over steps):\n")
print(aggregate(median_abs ~ width, G, mean), row.names = FALSE, digits = 3)
cat("median |score| by step type (averaged over widths):\n")
print(aggregate(median_abs ~ step_lab, G, mean), row.names = FALSE, digits = 3)
cat("\nDONE\n")
