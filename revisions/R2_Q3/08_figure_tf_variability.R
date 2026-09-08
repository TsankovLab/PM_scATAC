###############################################################################
# STEP 8 -- TF variability figures.  Everything here is UNSIGNED: which clone is
# labelled c1 is arbitrary, so only the magnitude of the between-clone difference,
# |Δz|, is ever plotted.
#
# Figure 1, Plots/epi_TF_variability.pdf
#  A  per tumour: the ARI between an independent chromVAR clustering and the CNV
#     clones. Near zero everywhere = the clones are not distinct regulatory states.
#  B  the most variable motifs, median |Δz| +/- IQR over the tumours that split.
#  C  the same without any selection: all 869 motifs ranked within each tumour.
#     Panel B is a ranking and a ranking always has a top; panel C shows the curve
#     the top is being taken off, and where the technical and immune motifs sit on it.
#
# Figure 2, Plots/epi_TF_genescore.pdf -- the same contrast on a second ATAC readout.
#  D  per tumour: |Δz| (motif activity between the clones) against |Δ gene score| (the
#     TF's own locus between the same two clones). Both are unsigned between-clone
#     variabilities; the question is whether a motif that moves belongs to a TF whose
#     locus also moves.
#     Genes sitting on the ARM THAT DEFINES THE SPLIT are drawn ringed: a clone with an
#     extra copy of that arm has more fragments across all of it, so those genes gain
#     score for reasons that are not regulatory. That is visible as a shifted subset.
#     CAVEAT: gene score comes from the SAME fragments as the chromVAR deviations, so
#     this is an internal consistency check, not independent evidence.
#  E  the same collapsed to one point per motif (median over tumours), square, immune
#     drawn on top, and every TF in the TOP-RIGHT QUADRANT named -- the motifs that move
#     between the clones AND whose own locus moves. This is the headline panel and it is
#     written on its own to Plots/epi_TF_genescore.pdf; D and F go to the supporting file.
#  F  median |Δz| of the immune motifs against everything else, with a Wilcoxon test.
#     Computed on ALL 869 motifs over ALL the tumours that split (not the expression
#     subset of D/E), so it matches the Fisher test in step 7. Two different questions
#     give two different answers and both belong here: the immune DISTRIBUTION is
#     shifted up, while immune motifs are NOT over-represented among the largest movers.
#
# Input : epi_chromvar_recurrence.csv, epi_chromvar_mirror.csv,
#         epi_chromvar_diff_<S>.csv, archr_tf_genescore.csv (step 7)
# Output: Plots/epi_TF_variability.pdf, Plots/epi_TF_genescore.pdf (panel E alone),
#         Plots/epi_TF_genescore_supporting.pdf (panels D and F)
###############################################################################
suppressMessages({ library(ggplot2); library(patchwork); library(ggrepel) })
source("00_common.R")
COL_TECH <- "#b0413e"; COL_IMM <- "#7b3fa0"; COL_OTH <- "grey55"
QUAD_Q <- 0.85     # the top-right quadrant of panel E: this quantile on BOTH axes
pal <- c(`AP-1/bZIP (technical)` = COL_TECH, `immune / interferon` = COL_IMM,
         `other TF` = COL_OTH)

rec  <- read.csv("epi_chromvar_recurrence.csv", stringsAsFactors = FALSE)
metr <- read.csv("epi_chromvar_mirror.csv",     stringsAsFactors = FALSE)
EX   <- read.csv("archr_tf_genescore_clonediff.csv", stringsAsFactors = FALSE)
rec$class <- ifelse(rec$tech, "AP-1/bZIP (technical)",
             ifelse(rec$immune, "immune / interferon", "other TF"))

## per-sample motif tables, in the order the tumours are listed in SAMPLES
SPL <- metr$sample[order(match(metr$sample, SAMPLES))]
D <- do.call(rbind, lapply(SPL, function(S)
       read.csv(sprintf("epi_chromvar_diff_%s.csv", S), stringsAsFactors = FALSE)))
D$class <- rec$class[match(D$TF, rec$TF)]

theme_p <- theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6.5),
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text = element_text(size = 6.4),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"))

## ---- A ----------------------------------------------------------------------
pA <- ggplot(metr, aes(reorder(sample, ARI_chromvar_vs_CNV), ARI_chromvar_vs_CNV)) +
  geom_col(aes(fill = sample == "P4"), width = 0.68, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = COL_HI, `FALSE` = "grey65")) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) + coord_flip() +
  labs(x = NULL, y = "ARI: chromVAR clustering vs epiAneufinder clones",
       title = "A  Does TF activity recover the clone split?",
       subtitle = "0 = no correspondence") +
  theme_p + theme(legend.position = "none")

## ---- B ----------------------------------------------------------------------
top <- head(rec[order(-rec$median_absdiff), ], 22)
top$TF <- factor(top$TF, levels = rev(top$TF))
pB <- ggplot(top, aes(median_absdiff, TF, colour = class)) +
  geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 0.45) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = "median |Δz| between clones   (bars = IQR over tumours)", y = NULL,
       title = "B  Most variable TF motifs between epiAneufinder subclones",
       subtitle = "median across samples, so a TF must move in most tumours, not one") +
  theme_p + theme(axis.text.y = element_text(size = 6.6))

## ---- C: the whole ranking, per tumour ---------------------------------------
D <- do.call(rbind, lapply(split(D, D$sample), function(x){
       x <- x[order(-x$absdiff), ]; x$rank <- seq_len(nrow(x)); x }))
lab <- setNames(sprintf("%s   max |Δz| %.2f | mean η² %.4f", SPL,
                        sapply(SPL, function(S) max(D$absdiff[D$sample == S])),
                        metr$mean_eta2[match(SPL, metr$sample)]), SPL)
D$facet <- factor(lab[D$sample], levels = lab[SPL])
D <- D[order(D$class == "other TF", decreasing = TRUE), ]        # grey drawn first
markC <- do.call(rbind, lapply(split(D, D$sample), function(x) head(x[order(x$rank), ], 5)))

pC <- ggplot(D, aes(rank, absdiff)) +
  geom_point(aes(colour = class), size = 0.5, alpha = 0.7, stroke = 0) +
  geom_point(data = markC, colour = "black", shape = 21, size = 1.3, stroke = 0.4, fill = NA) +
  geom_text_repel(data = markC, aes(label = TF, colour = class), size = 1.9,
                  min.segment.length = 0, segment.size = 0.2, segment.colour = "grey55",
                  max.overlaps = 20, nudge_x = 90, show.legend = FALSE) +
  facet_wrap(~ facet, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = pal, name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.06)) +
  guides(colour = "none") +                   # same colour key as panel B, drawn once
  labs(x = "motifs, ranked within the tumour (1 = most variable, 869 = least)",
       y = "|Δz| between the two clones",
       title = "C  All 869 motifs, ranked within each tumour",
       subtitle = paste0("no direction: the clone labels carry no order. Circled and ",
                         "labelled = the 5 most variable in that tumour.\n",
                         "The curve falls away immediately in every tumour but P11")) +
  theme_p

fig1 <- ((pA | pB) + plot_layout(widths = c(1, 1.35))) / pC + plot_layout(heights = c(1, 0.88))
ggsave("Plots/epi_TF_variability.pdf", fig1, width = 10.5, height = 8.6, device = cairo_pdf)
cat("DONE -> Plots/epi_TF_variability.pdf\n")

###############################################################################
# Figure 2 -- variability against the gene score of the same TF in the same tumour.
###############################################################################
have <- intersect(SPL, unique(EX$sample))          # gene score covers every tumour
cat("tumours in the gene-score comparison:", paste(have, collapse = ", "), "\n")
E <- merge(D[D$sample %in% have, c("sample","TF","absdiff","eta2","class")],
           EX[, c("sample","gene","absdiff_gs","mean_score","arm","on_split_arm")],
           by.x = c("sample","TF"), by.y = c("sample","gene"))
cat(sprintf("motifs with a matched gene: %d of %d\n",
            length(unique(E$TF)), length(unique(D$TF))))

## do the two between-clone variabilities agree, and how much of the gene-score one is
## just the copy-number difference that defines the clones?
qs <- do.call(rbind, lapply(split(E, E$sample), function(x){
  t20 <- head(x[order(-x$absdiff), ], 20)
  data.frame(sample = x$sample[1],
             rho = cor(x$absdiff_gs, x$absdiff, method = "spearman"),
             med_gs_on  = median(x$absdiff_gs[x$on_split_arm]),
             med_gs_off = median(x$absdiff_gs[!x$on_split_arm]),
             top20_gs = median(t20$absdiff_gs)) }))
cat("\n=== |dz| vs |d gene score| ===\n")
print(transform(qs, rho = round(rho, 3), med_gs_on = round(med_gs_on, 3),
                med_gs_off = round(med_gs_off, 3), top20_gs = round(top20_gs, 3)),
      row.names = FALSE)

E <- E[order(E$class == "other TF", decreasing = TRUE), ]
elab <- setNames(sprintf("%s   rho %.2f   %s %.3f / off %.3f", have,
                         qs$rho[match(have, qs$sample)],
                         EX$arm[match(paste(have, TRUE), paste(EX$sample, EX$on_split_arm))],
                         qs$med_gs_on[match(have, qs$sample)],
                         qs$med_gs_off[match(have, qs$sample)]), have)
E$facet <- factor(elab[E$sample], levels = elab[have])
onarm <- E[E$on_split_arm, ]
markD <- do.call(rbind, lapply(split(E, E$sample), function(x) head(x[order(-x$absdiff), ], 6)))

pD <- ggplot(E, aes(absdiff_gs, absdiff)) +
  geom_point(aes(colour = class), size = 0.55, alpha = 0.7, stroke = 0) +
  geom_point(data = onarm, aes(fill = class), shape = 21, size = 1.5, stroke = 0.35,
             colour = "grey15", alpha = 0.95) +
  geom_point(data = markD, colour = "black", shape = 21, size = 1.4, stroke = 0.4, fill = NA) +
  geom_text_repel(data = markD, aes(label = TF, colour = class), size = 2,
                  min.segment.length = 0, segment.size = 0.2, segment.colour = "grey55",
                  max.overlaps = 25, show.legend = FALSE) +
  facet_wrap(~ facet, nrow = 1, scales = "free") +
  scale_colour_manual(values = pal, name = NULL) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  labs(x = "|Δ gene score| of the TF's own gene between the same two clones   log2(score + 1)",
       y = "|Δz| between the two clones",
       title = "D  Does the TF's own locus move between the clones as well?",
       subtitle = paste0("both axes are unsigned between-clone differences. Black rings = the ",
                        "6 most variable motifs of that tumour;\ngrey-outlined points = genes ",
                        "on the arm that DEFINES the split, which gain score from the extra copy ",
                        "alone\n(their median |Δ gene score| is given in each strip against the ",
                        "off-arm median). Same fragments as the chromVAR z: a consistency check.")) +
  theme_p

## ---- E / F: one point per motif, tumours collapsed to the median -------------
## medians over the tumours, and the on-arm genes dropped from the gene-score median so
## the copy-number effect does not leak into it
qq <- function(v, p) if (length(v)) quantile(v, p, names = FALSE) else NA_real_
MED <- do.call(rbind, lapply(split(E, E$TF), function(x){
  off <- x[!x$on_split_arm, ]
  data.frame(TF = x$TF[1], class = x$class[1], n = nrow(x), n_off = nrow(off),
             madiff = median(x$absdiff),
             mad_lo = qq(x$absdiff, .25),      mad_hi = qq(x$absdiff, .75),
             expr   = qq(off$absdiff_gs, .50),
             expr_lo = qq(off$absdiff_gs, .25), expr_hi = qq(off$absdiff_gs, .75),
             stringsAsFactors = FALSE) }))
MED <- MED[!is.na(MED$expr), ]
cat(sprintf("median |dz| vs median |d gene score| (on-arm genes excluded): Spearman %.3f\n",
            cor(MED$madiff, MED$expr, method = "spearman")))
MED <- MED[order(MED$class == "other TF", decreasing = TRUE), ]
imm <- MED[MED$class == "immune / interferon", ]
## ---- the top-right quadrant: motif moves AND locus moves --------------------
QX <- quantile(MED$expr, QUAD_Q); QY <- quantile(MED$madiff, QUAD_Q)
quad <- MED[MED$expr >= QX & MED$madiff >= QY, ]
quad <- quad[order(-quad$madiff), ]
cat(sprintf("\ntop-right quadrant (>= %.0fth pct on both: |d gene score| %.3f, |dz| %.3f): %d motifs\n",
            100 * QUAD_Q, QX, QY, nrow(quad)))
print(quad[, c("TF","class","expr","madiff")], row.names = FALSE, digits = 3)
write.csv(quad, "epi_TF_genescore_quadrant.csv", row.names = FALSE)

## one colour for everything: the point of this panel is the quadrant, not the classes.
## The AP-1 / immune split is still colour-coded in panels B, D and F.
COL_PT <- "grey62"; COL_QUAD <- "grey38"; COL_LAB <- "grey35"
pE <- ggplot(MED, aes(expr, madiff)) +
  annotate("rect", xmin = QX, xmax = Inf, ymin = QY, ymax = Inf,
           fill = "grey80", alpha = 0.28) +
  geom_vline(xintercept = QX, linetype = 2, colour = "grey45", linewidth = 0.3) +
  geom_hline(yintercept = QY, linetype = 2, colour = "grey45", linewidth = 0.3) +
  geom_point(colour = COL_PT, size = 0.8, alpha = 0.55, stroke = 0) +
  ## IQR over the tumours, drawn for the named motifs only -- the background cloud would
  ## be unreadable with 772 pairs of bars
  geom_segment(data = quad, aes(x = expr_lo, xend = expr_hi, y = madiff, yend = madiff),
               colour = "grey72", linewidth = 0.25) +
  geom_segment(data = quad, aes(x = expr, xend = expr, y = mad_lo, yend = mad_hi),
               colour = "grey72", linewidth = 0.25) +
  geom_point(data = quad, shape = 21, size = 2, stroke = 0.35,
             fill = COL_QUAD, colour = "grey15") +
  geom_text_repel(data = quad, aes(label = TF), colour = COL_LAB, size = 3.2,
                  min.segment.length = 0, segment.size = 0.2,
                  segment.colour = "grey55", max.overlaps = Inf, box.padding = 0.36,
                  point.padding = 0.2, force = 3, force_pull = 0.6, seed = 1) +
  scale_x_continuous(expand = expansion(mult = 0.16)) +
  scale_y_continuous(expand = expansion(mult = 0.16)) +
  labs(x = "median |\u0394 gene score| of the TF's own gene between the clones\n(genes on the arm that defines the split are excluded)",
       y = "median |\u0394z| of the TF motif between the clones",
       title = "TF variability between subclones: motif activity against the TF's own locus",
       subtitle = sprintf(paste0("one point per motif, median over the %d tumours that split ",
                          "(n = %d motifs with a matched gene). Spearman %.2f.\nShaded = the ",
                          "top-right quadrant, above the %.0fth percentile on BOTH axes ",
                          "(%d motifs, all named):\nthe motifs that move between the clones and ",
                          "whose TF locus moves as well.\nBars = IQR over the tumours: wide ",
                          "because the tumours differ several-fold in how far apart their ",
                          "clones sit,\nwhich is why the summary is a median."),
                          nrow(metr), nrow(MED), cor(MED$madiff, MED$expr, method = "spearman"),
                          100 * QUAD_Q, nrow(quad))) +
  theme_p + theme(aspect.ratio = 1, legend.position = "none")

ggsave("Plots/epi_TF_genescore.pdf", pE, width = 5.6, height = 6.5, device = cairo_pdf)
cat("DONE -> Plots/epi_TF_genescore.pdf\n")

## F uses every motif and every tumour that split, so it is comparable with step 7
wt <- wilcox.test(median_absdiff ~ I(class == "immune / interferon"),
                  data = rec[rec$class != "AP-1/bZIP (technical)", ])
cat(sprintf("immune vs other over all %d motifs and %d tumours: %.3f vs %.3f, Wilcoxon p = %.3g\n",
            nrow(rec), nrow(metr), median(rec$median_absdiff[rec$class == "immune / interferon"]),
            median(rec$median_absdiff[rec$class == "other TF"]), wt$p.value))

pF <- ggplot(rec, aes(class, median_absdiff, colour = class)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.16, size = 0.6, alpha = 0.45, stroke = 0) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_discrete(labels = function(x) sub(" / ", "/\n", sub(" \\(", "\n(", x))) +
  labs(x = NULL, y = "median |Δz| between clones",
       title = "F  Immune motifs shift a little; none of them lead",
       subtitle = sprintf(paste0("all %d motifs, all %d tumours that split. Immune %.3f vs\n",
                                 "%.3f for the other non-technical motifs (Wilcoxon\n",
                                 "p = %.4f), but immune motifs are NOT over-represented\n",
                                 "among the largest movers (top-decile Fisher OR 1.26, p = 0.59)"),
                          nrow(rec), nrow(metr),
                          median(rec$median_absdiff[rec$class == "immune / interferon"]),
                          median(rec$median_absdiff[rec$class == "other TF"]),
                          wt$p.value)) +
  theme_p + theme(legend.position = "none", axis.text.x = element_text(size = 6.5))

## supporting: the per-tumour version (D) and the immune comparison (F)
sup <- pD / pF + plot_layout(heights = c(1, 1))
ggsave("Plots/epi_TF_genescore_supporting.pdf", sup, width = 11.5, height = 8.6,
       device = cairo_pdf)
cat("DONE -> Plots/epi_TF_genescore_supporting.pdf\n")
