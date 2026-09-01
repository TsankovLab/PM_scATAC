###############################################################################
# STEP 6 -- is P4's split really the chr8q amplification?
#
# The single-cell half of the P4 argument (step 9 is the spatial half). Everything
# here comes out of the epiAneufinder calls alone -- no other CNV pipeline is involved.
#
#  1. per cell, the fraction of chr8q bins called GAIN, and the same for chr8p as an
#     internal control: a whole-chromosome or global artefact would move both arms
#  2. the unsupervised clones from step 3 vs those chr8q fractions
#  3. a TARGETED split (k-means, k = 2, on the chr8q fraction alone) as the reference
#     partition, and the ARI of the unsupervised clones against it -- i.e. did the
#     genome-wide clustering find the chr8q clone without being pointed at chr8?
#
# Input : out_5Mb/P4/..., epi_clone_labels.csv (step 3)
# Output: P4_chr8q_epiAneufinder.csv, Plots/P4_chr8q_epiAneufinder.pdf
###############################################################################
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("00_common.R")
S <- "P4"

e <- read_epi(S); M <- e$M; d <- e$d
is8q <- d$seq == "chr8" & d$start >  CEN8
is8p <- d$seq == "chr8" & d$start <= CEN8
cat(sprintf("%s: %d bins x %d cells | chr8q bins %d, chr8p bins %d\n",
            S, nrow(M), ncol(M), sum(is8q), sum(is8p)))
cat("call distribution (0=loss, 1=normal, 2=gain):\n"); print(round(prop.table(table(M)), 3))

## ---- 1. per-cell arm-level gain fractions -----------------------------------
f8q <- colMeans(M[is8q, , drop = FALSE] == 2)
f8p <- colMeans(M[is8p, , drop = FALSE] == 2)
fgw <- colMeans(M == 2)
cat(sprintf("\nchr8q gain fraction: mean %.3f | chr8p %.3f | genome-wide %.3f\n",
            mean(f8q), mean(f8p), mean(fgw)))

## ---- 2. the unsupervised clones ---------------------------------------------
LAB <- fread("epi_clone_labels.csv")[sample == S]
cl  <- setNames(LAB$clone, sub("^P4#", "", LAB$cell))
cl  <- cl[intersect(names(cl), colnames(M))]
cat("\nunsupervised clones:\n"); print(table(cl))
cat("\nchr8q gain fraction by unsupervised clone:\n")
print(round(sapply(split(f8q[names(cl)], cl), function(z)
        c(mean = mean(z), median = median(z))), 3))
cat("chr8p (control arm):\n")
print(round(sapply(split(f8p[names(cl)], cl), mean), 3))
w8q <- wilcox.test(f8q[names(cl)] ~ cl); w8p <- wilcox.test(f8p[names(cl)] ~ cl)
cat(sprintf("Wilcoxon  chr8q p = %.3g | chr8p p = %.3g\n", w8q$p.value, w8p$p.value))

## ---- 3. targeted chr8q split, and does the unsupervised one match it? -------
km  <- kmeans(scale(f8q), 2, nstart = 25)$cluster
hi  <- which.max(tapply(f8q, km, mean))
ref <- setNames(ifelse(km == hi, "chr8q-amp", "chr8q-low"), names(f8q))
tb  <- table(unsupervised = cl, targeted = ref[names(cl)])
cat("\n=== unsupervised clones vs the targeted chr8q split ===\n"); print(tb)
amp_clone <- names(which.max(tapply(f8q[names(cl)], cl, mean)))
prec <- tb[amp_clone, "chr8q-amp"] / sum(tb[amp_clone, ])
rec  <- tb[amp_clone, "chr8q-amp"] / sum(tb[, "chr8q-amp"])
ari  <- ARI(cl, ref[names(cl)])
cat(sprintf("ARI = %.3f | clone %s: precision %.3f, recall %.3f for chr8q-amp\n",
            ari, amp_clone, prec, rec))

write.csv(data.frame(sample = S, n_cells = ncol(M),
  n_amp_clone = sum(cl == amp_clone), frac_amp_clone = mean(cl == amp_clone),
  chr8q_amp = round(mean(f8q[names(cl)][cl == amp_clone]), 3),
  chr8q_other = round(mean(f8q[names(cl)][cl != amp_clone]), 3),
  chr8p_amp = round(mean(f8p[names(cl)][cl == amp_clone]), 3),
  chr8p_other = round(mean(f8p[names(cl)][cl != amp_clone]), 3),
  p_chr8q = signif(w8q$p.value, 3), p_chr8p = signif(w8p$p.value, 3),
  ARI_vs_targeted = round(ari, 3), precision = round(prec, 3), recall = round(rec, 3)),
  "P4_chr8q_epiAneufinder.csv", row.names = FALSE)

## ---- figure -----------------------------------------------------------------
df <- data.frame(cell = names(cl), clone = unname(cl),
                 chr8q = f8q[names(cl)], chr8p = f8p[names(cl)])
pal <- setNames(ifelse(sort(unique(cl)) == amp_clone, COL_HI, "grey60"), sort(unique(cl)))

long <- rbind(data.frame(clone = df$clone, arm = "chr8q (amplified)", v = df$chr8q),
              data.frame(clone = df$clone, arm = "chr8p (control)",   v = df$chr8p))
long$arm <- factor(long$arm, levels = c("chr8q (amplified)", "chr8p (control)"))
pA <- ggplot(long, aes(clone, v, fill = clone)) +
  geom_violin(scale = "width", colour = NA, alpha = 0.85) +
  geom_boxplot(width = 0.16, outlier.size = 0.25, fill = "white", linewidth = 0.3) +
  facet_wrap(~ arm) + scale_fill_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "fraction of bins called GAIN, per cell",
       title = "A  The P4 clone split is a chr8q event",
       subtitle = sprintf(paste0("chr8q %.2f vs %.2f (p = %.1g)\n",
                                 "chr8p control shifts the same way but %.0fx less (%.3f vs %.3f, p = %.1g)"),
                          mean(df$chr8q[df$clone != amp_clone]), mean(df$chr8q[df$clone == amp_clone]),
                          w8q$p.value,
                          diff(range(tapply(df$chr8q, df$clone, mean))) /
                            diff(range(tapply(df$chr8p, df$clone, mean))),
                          mean(df$chr8p[df$clone != amp_clone]), mean(df$chr8p[df$clone == amp_clone]),
                          w8p$p.value)) +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"),
        strip.background = element_rect(fill = "grey93", colour = NA),
        panel.grid.minor = element_blank())

## chr8 profile of each clone along the chromosome
i8 <- which(d$seq == "chr8")
pr <- do.call(rbind, lapply(sort(unique(cl)), function(g){
  cells <- names(cl)[cl == g]
  data.frame(clone = g, mb = d$start[i8] / 1e6,
             v = rowMeans(M[i8, cells, drop = FALSE] == 2) -
                 rowMeans(M[i8, cells, drop = FALSE] == 0)) }))
pB <- ggplot(pr, aes(mb, v, colour = clone)) +
  geom_vline(xintercept = CEN8 / 1e6, linetype = 2, colour = "grey45", linewidth = 0.3) +
  geom_line(linewidth = 0.7) + geom_point(size = 1) +
  scale_colour_manual(values = pal, name = NULL) +
  annotate("text", x = CEN8/1e6, y = Inf, label = " centromere", hjust = 0, vjust = 1.4,
           size = 2.1, colour = "grey45") +
  labs(x = "chr8 position (Mb)", y = "gain - loss",
       title = "B  chr8 CNV profile of the two P4 clones",
       subtitle = "the difference is concentrated distal to the centromere\n") +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 6.5, colour = "grey35"),
        legend.position = "inside", legend.position.inside = c(0.08, 0.85),
        legend.key.size = unit(3, "mm"),
        legend.background = element_rect(fill = NA, colour = NA),
        panel.grid.minor = element_blank())

ggsave("Plots/P4_chr8q_epiAneufinder.pdf", (pA | pB) + plot_layout(widths = c(1, 1.25)),
       width = 9.5, height = 3.6, device = cairo_pdf)
cat("\nDONE -> Plots/P4_chr8q_epiAneufinder.pdf\n")
