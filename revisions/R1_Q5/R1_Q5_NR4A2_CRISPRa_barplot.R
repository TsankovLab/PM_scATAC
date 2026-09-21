###############################################################################
# R1_Q5 -- Barplot of the NR4A2 CRISPRa qPCR result, both normalisers.
#
# Plots the analysis set used in the source spreadsheet's final columns:
#   targeting  NR4A2 2, NR4A2 4, NR4A2 5
#   scramble   Scramble 1, Scramble 3, Scramble 6/24/26
# Values are recomputed from the raw Ct in the Excel file (not the pasted derived columns) and
# checked against the spreadsheet's own dCt / ddCt / 2^-ddCt, so the figure is reproducible
# from source.
#
#   dCt   = mean Ct(NR4A2) - mean Ct(reference)     per biological replicate
#   ddCt  = dCt - mean dCt(scramble)
#   fold  = 2^(-ddCt)
#   bar   = arithmetic mean of fold changes, matching the spreadsheet's "Average fold change"
#           (the geometric mean is also printed: it is the statistically preferable summary
#           for ratios, and it is what the dCt-scale test actually compares)
#   error = SEM of the fold changes across biological replicates; points = replicates
#   p     = two-sample equal-variance t-test on dCt, which reproduces the spreadsheet's
#           p = 0.003311 (beta-actin) and 0.003808 (18S) exactly; Welch is printed alongside
#
# NOTE the source file's 18S block mis-states ddCt for NR4A2 2 (-5.0021 instead of -3.8892),
# inflating its fold change to 32.05 and the reported 18S average to 22.6; recomputed here as
# 14.82 and 16.8. The p-values are unaffected.
#
# Output: Plots/NR4A2_CRISPRa_barplot.pdf, NR4A2_CRISPRa_barplot_values.csv
###############################################################################
suppressPackageStartupMessages({
  library(readxl); library(ggplot2); library(dplyr); library(RColorBrewer)
  library(paletteer); library(circlize)
})
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
setwd(file.path(ROOT, "git_repo_claude", "R1_Q5")); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
XLSX <- file.path(ROOT, "NR4A2_CRISPRa_Lillian", "7_15_26_qPCR_JK_CRISPRa.xlsx")
KEEP_T <- c("NR4A2 2", "NR4A2 4", "NR4A2 5")
KEEP_S <- c("Scramble 1", "Scramble 3", "Scramble 6/24")
GRP_COL <- c(Scramble = "grey55", `NR4A2 enhancer` = "#B2182B")

raw <- suppressMessages(read_excel(XLSX, sheet = 1, col_names = FALSE, .name_repair = "minimal"))
num <- function(r, cols) as.numeric(unlist(raw[r, cols]))
rows <- list(`NR4A2 1` = 2, `NR4A2 2` = 3, `NR4A2 3` = 4, `NR4A2 4` = 5, `NR4A2 5` = 6,
             `Scramble 1` = 9, `Scramble 2` = 10, `Scramble 3` = 11, `Scramble 6/24` = 12)
CT <- do.call(rbind, Map(function(id, r) data.frame(replicate = id,
  group = if (grepl("^NR4A2", id)) "NR4A2 enhancer" else "Scramble",
  NR4A2 = mean(num(r, 2:4)), `beta-actin` = mean(num(r, 5:7)), `18S` = mean(num(r, 8:10)),
  check.names = FALSE, stringsAsFactors = FALSE), names(rows), rows))
CT <- CT[CT$replicate %in% c(KEEP_T, KEEP_S), ]
CT$group <- factor(CT$group, levels = c("Scramble", "NR4A2 enhancer"))

D <- do.call(rbind, lapply(c("beta-actin", "18S"), function(n) {
  d <- data.frame(CT[, c("replicate", "group")], normaliser = n, dCt = CT$NR4A2 - CT[[n]])
  d$ddCt <- d$dCt - mean(d$dCt[d$group == "Scramble"]); d$fold <- 2^(-d$ddCt); d }))
D$normaliser <- factor(D$normaliser, levels = c("beta-actin", "18S"))

## Check against the spreadsheet's own derived columns. One cell disagrees: in the 18S block
## NR4A2 2 is given ddCt = -5.0021 (fold 32.05), but its own dCt (9.6967) minus the block's
## stated average scramble dCt (13.58589) is -3.8892 (fold 14.82) -- every other cell in that
## block is self-consistent, so this is an isolated formula error. It inflates the reported
## average 18S fold change from 16.8 to 22.6. It does NOT affect the p-values, which are
## computed on dCt and reproduce exactly. The recomputed value is used here.
ref <- data.frame(
  replicate  = rep(c(KEEP_T, KEEP_S), 2),
  normaliser = rep(c("beta-actin", "18S"), each = 6),
  fold_sheet = c(7.817881829, 19.24543592, 6.345693226, 1.286493747, 0.824797069, 0.942421557,
                 32.04686022, 28.62536576, 7.103625648, 1.268293716, 1.133320222, 0.695708839),
  stringsAsFactors = FALSE)
m <- merge(D, ref, by = c("replicate", "normaliser"))
m$diff <- m$fold - m$fold_sheet
bad <- m[abs(m$diff) > 1e-3, ]
cat("replicates matching the spreadsheet:", sum(abs(m$diff) <= 1e-3), "of", nrow(m), "\n")
if (nrow(bad)) {
  cat("MISMATCH with the source file (recomputed value used):\n")
  print(bad[, c("replicate", "normaliser", "dCt", "ddCt", "fold", "fold_sheet")], row.names = FALSE, digits = 6)
}
stopifnot(nrow(bad) <= 1)          # a single known bad cell; more would mean a parsing problem

SUM <- D %>% group_by(normaliser, group) %>%
  summarise(n = n(), mean_fold = mean(fold), sem = sd(fold) / sqrt(n()),
            geom_mean_fold = 2^(-mean(ddCt)), .groups = "drop")
ST <- D %>% group_by(normaliser) %>% summarise(
  p_student = t.test(dCt[group == "NR4A2 enhancer"], dCt[group == "Scramble"], var.equal = TRUE)$p.value,
  p_welch   = t.test(dCt[group == "NR4A2 enhancer"], dCt[group == "Scramble"])$p.value, .groups = "drop")
write.csv(merge(SUM, ST, by = "normaliser"), "NR4A2_CRISPRa_barplot_values.csv", row.names = FALSE)
cat("=== bar values ===\n"); print(as.data.frame(SUM %>% mutate(across(where(is.numeric), ~round(., 3)))), row.names = FALSE)
cat("\n=== tests on dCt ===\n"); print(as.data.frame(ST %>% mutate(across(where(is.numeric), ~signif(., 4)))), row.names = FALSE)

## ---- figure ---------------------------------------------------------------------------
BR <- merge(SUM, ST, by = "normaliser")
BR$lab <- sprintf("p = %s", signif(BR$p_student, 3))
brk <- BR %>% group_by(normaliser) %>%
  summarise(y = max(mean_fold + sem, max(D$fold[D$normaliser == normaliser[1]])) * 1.12,
            lab = lab[1], .groups = "drop")
p <- ggplot(SUM, aes(group, mean_fold, fill = group)) +
  geom_col(width = .6, colour = "grey25", linewidth = .3) +
  geom_errorbar(aes(ymin = pmax(mean_fold - sem, 0), ymax = mean_fold + sem), width = .18, linewidth = .4) +
  geom_point(data = D, aes(group, fold), inherit.aes = FALSE, size = 1.7,
             position = position_jitter(width = .1, seed = 1), colour = "grey20") +
  geom_segment(data = brk, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE, linewidth = .3) +
  geom_text(data = brk, aes(x = 1.5, y = y, label = lab), inherit.aes = FALSE,
            vjust = -0.45, size = 2.7) +
  scale_fill_manual(values = GRP_COL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, .16))) +
  facet_wrap(~ normaliser, nrow = 1, scales = "free_y") +
  gtheme_no_rot + labs(x = NULL, y = "NR4A2 fold change vs scramble",
    title = "CRISPRa of the NR4A2 distal enhancer in Jurkat cells",
    subtitle = "mean +/- SEM of 3 biological replicates per group; points = replicates; equal-variance t-test on dCt") +
  theme(plot.title = element_text(size = 10.5), plot.subtitle = element_text(size = 7.5))
ggsave("Plots/NR4A2_CRISPRa_barplot.pdf", p, width = 5.6, height = 3.6)
cat("\nDONE\n")
