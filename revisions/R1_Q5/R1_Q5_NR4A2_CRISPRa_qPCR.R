###############################################################################
# R1_Q5 -- Functional validation of the NR4A2 distal enhancer by CRISPRa in Jurkat cells.
#
# Reviewer 1 Q5 asks for functional validation of the NR4A2 distal enhancer. A dCas9-VP64
# (CRISPRa) assay targeting that element was run in Jurkat (JK) cells and NR4A2 mRNA measured
# by RT-qPCR against two normalisers (beta-actin, 18S).
#
# Source (read-only, never modified): NR4A2_CRISPRa_Lillian/7_15_26_qPCR_JK_CRISPRa.xlsx
#   rows  2-6   NR4A2-targeting gRNA, biological replicates 1-5   (7/15/26)
#   rows  9-11  scramble gRNA, biological replicates 1-3          (7/15/26)
#   row  12     scramble gRNA from a SEPARATE experiment          (6/24/26)
#   each row = technical triplicate for NR4A2 | beta-actin | 18S
#
# This script recomputes everything from the raw Ct values rather than reusing the
# spreadsheet's derived columns, because the sheet contains three different analyses of the
# same data (all replicates; a subset with NR4A2 1, NR4A2 3 and Scramble 2 dropped; that
# subset plus the 6/24 scramble) whose p-values range from 0.49 to 0.003. All three are
# reproduced here as a sensitivity analysis so the choice is explicit rather than implicit.
#
# Method (standard delta-delta-Ct, tested on the log2 scale):
#   dCt      = mean Ct(NR4A2) - mean Ct(normaliser)        per biological replicate
#              (technical triplicates averaged first; dCt is already log2 units)
#   ddCt     = dCt - mean dCt(scramble)
#   fold     = 2^(-ddCt)
#   test     Welch t-test on dCt (NOT on fold change: fold changes are log-normal, so a
#            t-test on them is mis-specified and inflates the influence of high values)
#   effect   fold change = 2^(mean dCt scramble - mean dCt targeting), 95% CI from the t
#            interval for that difference, back-transformed
#   Normalisers: beta-actin, 18S, and their mean Ct (= geometric mean of the two, the
#   standard multi-reference normalisation).
#
# Outputs: NR4A2_CRISPRa_replicate_level.csv   per-replicate Ct, dCt, fold change
#          NR4A2_CRISPRa_summary_stats.csv     every analysis variant x normaliser
#          NR4A2_CRISPRa_technical_QC.csv      technical triplicate spread
#          Plots/NR4A2_CRISPRa_qPCR.pdf        main figure
#          Plots/NR4A2_CRISPRa_sensitivity.pdf analysis-choice sensitivity
###############################################################################
suppressPackageStartupMessages({
  library(readxl); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(RColorBrewer); library(paletteer); library(circlize)
})
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUT  <- file.path(ROOT, "git_repo_claude", "R1_Q5")
setwd(OUT); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
XLSX <- file.path(ROOT, "NR4A2_CRISPRa_Lillian", "7_15_26_qPCR_JK_CRISPRa.xlsx")
GRP_COL <- c(`NR4A2 enhancer` = "#B2182B", Scramble = "grey45")

## ---- 1. raw Ct -----------------------------------------------------------------------
raw <- suppressMessages(read_excel(XLSX, sheet = 1, col_names = FALSE, .name_repair = "minimal"))
num <- function(r, cols) as.numeric(unlist(raw[r, cols]))
blocks <- list(
  list(rows = 2:6,  grp = "NR4A2 enhancer", batch = "7/15/26", id = paste0("NR4A2 ", 1:5)),
  list(rows = 9:11, grp = "Scramble",       batch = "7/15/26", id = paste0("Scramble ", 1:3)),
  list(rows = 12,   grp = "Scramble",       batch = "6/24/26", id = "Scramble 6/24"))
CT <- do.call(rbind, lapply(blocks, function(b) do.call(rbind, Map(function(r, id)
  data.frame(replicate = id, group = b$grp, batch = b$batch,
             NR4A2 = mean(num(r, 2:4)), bactin = mean(num(r, 5:7)), rRNA18S = mean(num(r, 8:10)),
             sd_NR4A2 = sd(num(r, 2:4)), sd_bactin = sd(num(r, 5:7)), sd_18S = sd(num(r, 8:10)),
             stringsAsFactors = FALSE), b$rows, b$id))))
CT$group <- factor(CT$group, levels = c("Scramble", "NR4A2 enhancer"))

## verify against the spreadsheet's own averaged-Ct block (rows 16-20, 23-25, 27)
chk <- c(num(16:20, 2), num(23:25, 2), num(27, 2))
stopifnot(max(abs(chk - CT$NR4A2)) < 1e-3)
cat("technical-triplicate means reproduce the spreadsheet's averaged Ct (max diff",
    signif(max(abs(chk - CT$NR4A2)), 2), ")\n\n")
write.csv(CT[, c("replicate","group","batch","sd_NR4A2","sd_bactin","sd_18S")],
          "NR4A2_CRISPRa_technical_QC.csv", row.names = FALSE)
cat("=== technical triplicate SD (Ct) ===\n")
print(CT %>% transmute(replicate, group, batch, NR4A2 = round(sd_NR4A2, 3),
                       bactin = round(sd_bactin, 3), `18S` = round(sd_18S, 3)), row.names = FALSE)

## ---- 2. dCt per normaliser ------------------------------------------------------------
CT$mean_ref <- (CT$bactin + CT$rRNA18S) / 2          # geometric mean of the two references
NORM <- c(`beta-actin` = "bactin", `18S` = "rRNA18S", `geometric mean of both` = "mean_ref")
D <- do.call(rbind, lapply(names(NORM), function(n)
  data.frame(CT[, c("replicate","group","batch")], normaliser = n,
             Ct_NR4A2 = CT$NR4A2, Ct_ref = CT[[NORM[[n]]]], dCt = CT$NR4A2 - CT[[NORM[[n]]]])))

## ---- 3. analysis variants (exactly the three in the spreadsheet) -----------------------
DROPPED <- c("NR4A2 1", "NR4A2 3", "Scramble 2")
## QC_FAIL: the only replicate whose technical triplicates disagree badly (SD 0.70 Ct for
## NR4A2, ~3x every other well); excluding it is defensible from the assay itself, unlike the
## other two replicates the source file drops.
QC_FAIL <- "NR4A2 1"
VARIANTS <- list(
  `all replicates`                 = function(d) d[d$batch == "7/15/26", ],
  `all replicates + 6/24 scramble` = function(d) d,
  `drop NR4A2 1 (technical QC)`    = function(d) d[!d$replicate %in% QC_FAIL, ],
  `spreadsheet subset`             = function(d) d[d$batch == "7/15/26" & !d$replicate %in% DROPPED, ],
  `spreadsheet subset + 6/24`      = function(d) d[!d$replicate %in% DROPPED, ])
PRIMARY <- "all replicates + 6/24 scramble"

analyse <- function(d, variant, n) {
  s <- d$dCt[d$group == "Scramble"]; t <- d$dCt[d$group == "NR4A2 enhancer"]
  tt <- t.test(t, s)                                   # Welch, on the log2 (dCt) scale
  wt <- suppressWarnings(wilcox.test(t, s))            # rank alternative (min attainable p is
                                                       # 2/choose(n1+n2,n1), so it is weak here)
  ddct <- mean(t) - mean(s)                            # negative = induction
  ci <- -rev(tt$conf.int)                              # CI for (scramble - targeting)
  data.frame(variant, normaliser = n, n_targeting = length(t), n_scramble = length(s),
             mean_dCt_scramble = mean(s), mean_dCt_targeting = mean(t), ddCt = ddct,
             fold_change = 2^(-ddct), CI_low = 2^(ci[1]), CI_high = 2^(ci[2]),
             p = tt$p.value, p_wilcoxon = wt$p.value, row.names = NULL)
}
RES <- do.call(rbind, lapply(names(VARIANTS), function(v) do.call(rbind, lapply(names(NORM), function(n) {
  d <- VARIANTS[[v]](D[D$normaliser == n, ]); analyse(d, v, n) }))))
write.csv(RES, "NR4A2_CRISPRa_summary_stats.csv", row.names = FALSE)
cat("\n=== NR4A2 CRISPRa vs scramble: every analysis variant ===\n")
print(RES %>% transmute(variant, normaliser, n = sprintf("%d vs %d", n_targeting, n_scramble),
        `fold change` = sprintf("%.1f (%.1f-%.1f)", fold_change, CI_low, CI_high),
        ddCt = round(ddCt, 2), p = signif(p, 3), p_wilcox = signif(p_wilcoxon, 3)), row.names = FALSE)

## ---- 4. per-replicate fold change, relative to the scramble mean of its variant --------
REP <- do.call(rbind, lapply(names(NORM), function(n) {
  d <- VARIANTS[[PRIMARY]](D[D$normaliser == n, ])
  d$fold <- 2^(-(d$dCt - mean(d$dCt[d$group == "Scramble"]))); d }))
REP$excluded <- REP$replicate %in% DROPPED
write.csv(REP, "NR4A2_CRISPRa_replicate_level.csv", row.names = FALSE)
cat("\n=== per-replicate fold change (primary variant:", PRIMARY, ") ===\n")
print(REP %>% filter(normaliser == "geometric mean of both") %>%
        transmute(replicate, group, batch, dCt = round(dCt, 2), fold = round(fold, 2),
                  `dropped in sheet` = ifelse(excluded, "yes", "")), row.names = FALSE)

## ---- 5. figures ------------------------------------------------------------------------
lab <- RES %>% filter(variant == PRIMARY) %>%
  mutate(txt = sprintf("%.1f-fold\np = %s", fold_change, format.pval(p, digits = 2, eps = 1e-4)))
pA <- ggplot(REP, aes(group, dCt, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = .55, fill = NA, colour = "grey70", linewidth = .3) +
  geom_point(aes(shape = batch), size = 2.2, position = position_jitter(width = .12, seed = 1)) +
  scale_colour_manual(values = GRP_COL, guide = "none") +
  scale_shape_manual(values = c(`7/15/26` = 16, `6/24/26` = 17), name = "experiment") +
  scale_y_reverse() + facet_wrap(~ normaliser, nrow = 1) + gtheme_no_rot +
  labs(x = NULL, y = expression(Delta*"Ct (NR4A2 - reference), reversed")) +
  theme(legend.position = "top")
pB <- ggplot(REP, aes(group, fold, colour = group)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = .3, colour = "grey55") +
  stat_summary(fun = mean, geom = "crossbar", width = .45, linewidth = .3, colour = "grey30") +
  geom_point(aes(shape = batch), size = 2.2, position = position_jitter(width = .12, seed = 1)) +
  geom_text(data = lab, aes(x = 1.5, y = Inf, label = txt), inherit.aes = FALSE,
            vjust = 1.2, size = 2.6, colour = "grey20", lineheight = .95) +
  scale_colour_manual(values = GRP_COL, guide = "none") +
  scale_shape_manual(values = c(`7/15/26` = 16, `6/24/26` = 17), guide = "none") +
  scale_y_continuous(trans = "log2", breaks = c(.25, .5, 1, 2, 4, 8, 16, 32),
                     expand = expansion(mult = c(.05, .25))) +
  facet_wrap(~ normaliser, nrow = 1) + gtheme_no_rot +
  labs(x = NULL, y = "NR4A2 fold change vs scramble")
ggsave("Plots/NR4A2_CRISPRa_qPCR.pdf", (pA / pB) +
  plot_annotation(title = "CRISPRa of the NR4A2 distal enhancer in Jurkat cells",
                  subtitle = sprintf("RT-qPCR, %s; Welch t-test on dCt; bars = mean", PRIMARY),
                  theme = theme(plot.title = element_text(size = 11), plot.subtitle = element_text(size = 8))),
  width = 8, height = 6.5)

S <- RES; S$variant <- factor(S$variant, levels = rev(names(VARIANTS)))
pS <- ggplot(S, aes(fold_change, variant, colour = p < 0.05)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = .3, colour = "grey55") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = .18, linewidth = .4) +
  geom_point(size = 2.2) +
  geom_text(aes(label = sprintf("p = %s", signif(p, 2))), hjust = -0.25, vjust = -0.9, size = 2.2, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey55"), name = "p < 0.05") +
  scale_x_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 16, 32), expand = expansion(mult = c(.08, .3))) +
  facet_wrap(~ normaliser, nrow = 1) + gtheme_no_rot +
  labs(x = "NR4A2 fold change vs scramble (95% CI)", y = NULL,
       title = "Sensitivity of the CRISPRa result to replicate inclusion",
       subtitle = paste("'spreadsheet subset' drops", paste(DROPPED, collapse = ", "), "as in the source file"),
       caption = "Welch t-test on dCt") +
  theme(legend.position = "top", plot.title = element_text(size = 10), plot.subtitle = element_text(size = 7.5))
ggsave("Plots/NR4A2_CRISPRa_sensitivity.pdf", pS, width = 9, height = 4)
cat("\nDONE\n")
