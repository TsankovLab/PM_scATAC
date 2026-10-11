###############################################################################
# R1_Q8 -- sarcomatoid-score TF boxplots WITHOUT the per-cell scaling of TF activity
#
# Identical to R1_Q8_sarc_TF_boxplots_paired.R (bins of 40 cells, step 40, same 9 samples,
# same TF panel, paired lines, n = label) except for one step.
#
# The original does  mMat <- scale(mMat[selected_TF, ])  on a TF x cell matrix, so scale()
# acts on COLUMNS: every cell is centred and divided by the sd of the selected TFs WITHIN
# that cell. That is a cell-wise composition normalisation, not a per-TF z-score, and it is
# not monotone in any single TF -- so removing it genuinely changes each sample's Spearman
# correlation rather than just rescaling it. Here the motif deviations are used as they come
# out of chromVAR. Gene scores keep the original scaling (SCALE_GENESCORE), so the comparison
# is against an unchanged reference; set it FALSE for the symmetric version.
#
# Output: sarc_TF_correlations_per_sample_noscale.csv,
#         Plots/sarcomatoid_score_TF_activity_boxplots_noscale.pdf (+ .png)
#         and a side-by-side scaled-vs-unscaled panel.
###############################################################################
suppressPackageStartupMessages({ library(ArchR); library(dplyr); library(zoo); library(ggplot2)
                                 library(RColorBrewer); library(paletteer); library(circlize) })
addArchRThreads(8); addArchRGenome("hg38")
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
TUM  <- file.path(ROOT, "tumor_compartment", "scatac_ArchR")
OUT  <- file.path(ROOT, "git_repo_claude", "R1_Q8")
setwd(OUT); dir.create("Plots", showWarnings = FALSE)
source(file.path(ROOT, "git_repo", "utils", "scATAC_functions.R"))
source(file.path(ROOT, "git_repo", "utils", "ggplot_aestetics.R"))
source(file.path(ROOT, "git_repo", "utils", "palettes.R"))

SARC <- "cNMF20"; BIN <- 40; STEP <- 40
SCALE_ACTIVITY <- FALSE   # was TRUE in the original: scale() there z-scores each CELL across
SCALE_GENESCORE <- TRUE   # the selected TFs, so dropping it changes the values, not just units
DROP <- c("normal1", "normal2", "normal3", "normal4", "P3", "P13", "P11_HOX")

archp <- loadArchRProject(TUM, showLogo = FALSE)
meta  <- as.data.frame(archp@cellColData)
sams  <- setdiff(unique(meta$Sample3), DROP)
selected_TF <- readRDS(file.path(TUM, "selected_TF.rds"))
top_sarc_TF <- readRDS(file.path(TUM, "top_sarc_TF.rds"))
cat("samples used:", length(sams), "->", paste(sams, collapse = ", "), "\n")
cat("TFs plotted:", length(top_sarc_TF), "\n")

prep <- function(kind, do_scale) {
  se <- fetch_mat(archp, kind)
  M  <- assays(se)[[1]]; rownames(M) <- rowData(se)$name
  M  <- M[selected_TF, ]
  ## scale() on this TF x cell matrix normalises each CELL across the selected TFs
  ## scale() also densifies; without it the matrix stays sparse, so densify explicitly
  M  <- if (do_scale) t(scale(M)) else as.matrix(t(M))
  lapply(setNames(sams, sams), function(s) M[meta$Sample3 == s, , drop = FALSE])
}
mMat  <- prep("Motif",     SCALE_ACTIVITY)
gsMat <- prep("GeneScore", SCALE_GENESCORE)
cn <- as.data.frame(t(scale(t(meta[, grep("cNMF", colnames(meta))]))))
cn <- lapply(setNames(sams, sams), function(s) cn[meta$Sample3 == s, , drop = FALSE])

rollavg <- function(mat, ord) {
  mat <- mat[order(ord), , drop = FALSE]
  as.data.frame(lapply(as.data.frame(mat), function(col)
    zoo::rollapply(col, width = BIN, FUN = mean, by = STEP, partial = TRUE, align = "left")))
}
cor_df <- function(matlist, label) {
  do.call(rbind, lapply(sams, function(s) {
    ord <- cn[[s]][, SARC]
    a <- rollavg(matlist[[s]], ord); b <- rollavg(cn[[s]], ord)
    if (label == "genescore") a[is.na(a)] <- 0
    d <- as.data.frame(cor(a, b[, SARC], method = "spearman"))
    names(d) <- "score"; d$TF <- rownames(d); d$sample <- s; d$type <- label; d
  }))
}
D <- rbind(cor_df(mMat, "activity"), cor_df(gsMat, "genescore"))
D <- D[D$TF %in% top_sarc_TF, ]
D$TF <- factor(D$TF, levels = top_sarc_TF)
D$type <- factor(D$type, levels = c("activity", "genescore"))
write.csv(D, "sarc_TF_correlations_per_sample_noscale.csv", row.names = FALSE)
cat("\nper-TF sample counts (should equal", length(sams), "):",
    paste(range(table(D$TF, D$type)), collapse = "-"), "\n")

## explicit dodged x, so the paired lines land on the dots
DODGE <- 0.9; W <- 0.8
D$x <- as.numeric(D$TF) + ifelse(D$type == "activity", -1, 1) * DODGE / 4
pal <- paletteer::paletteer_c("grDevices::Purple-Blue", length(top_sarc_TF))
TYPE_COL <- c(activity = "#AE123AFF", genescore = "#001260FF")
n_lab <- sprintf("n = %d samples", length(sams))

p <- ggplot(D, aes(x = x, y = score)) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed", linewidth = .35) +
  geom_boxplot(aes(group = interaction(TF, type), fill = TF, colour = type),
               width = W / 2, linewidth = .4, alpha = .6, outlier.shape = NA) +
  geom_line(aes(group = interaction(TF, sample)), colour = "grey45", linewidth = .18, alpha = .55) +
  geom_point(aes(colour = type), size = 1.1, alpha = .8) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_colour_manual(values = TYPE_COL, name = NULL) +
  scale_x_continuous(breaks = seq_along(top_sarc_TF), labels = top_sarc_TF,
                     expand = expansion(add = .6)) +
  annotate("text", x = length(top_sarc_TF) + .4, y = max(D$score, na.rm = TRUE),
           label = n_lab, hjust = 1, vjust = 1, size = 2.6, colour = "grey25") +
  gtheme + labs(x = "TF", y = "score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5))

ggsave("Plots/sarcomatoid_score_TF_activity_boxplots_noscale.pdf", p, width = 11, height = 3.6)
ggsave("Plots/sarcomatoid_score_TF_activity_boxplots_noscale.png", p, width = 11, height = 3.6, dpi = 200)

###############################################################################
## side-by-side against the scaled version, so the effect of dropping scale() is visible
###############################################################################
ref_f <- "sarc_TF_correlations_per_sample.csv"
if (file.exists(ref_f)) {
  R <- read.csv(ref_f); R$scaling <- "scaled (original)"
  N <- D[, c("score", "TF", "sample", "type")]; N$scaling <- "unscaled activity"
  CMP <- rbind(R[, names(N)], N)
  CMP$TF <- factor(CMP$TF, levels = top_sarc_TF)
  CMP$scaling <- factor(CMP$scaling, levels = c("scaled (original)", "unscaled activity"))

  act <- CMP[CMP$type == "activity", ]
  pc <- ggplot(act, aes(TF, score, colour = scaling)) +
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed", linewidth = .35) +
    geom_boxplot(aes(fill = scaling), position = position_dodge(.75), width = .6,
                 linewidth = .35, alpha = .35, outlier.shape = NA) +
    geom_point(position = position_dodge(.75), size = .9, alpha = .8) +
    scale_colour_manual(values = c("scaled (original)" = "grey40",
                                   "unscaled activity" = "#AE123AFF"), name = NULL) +
    scale_fill_manual(values = c("scaled (original)" = "grey70",
                                 "unscaled activity" = "#AE123AFF"), name = NULL) +
    gtheme + labs(x = "TF", y = "activity score", title = "TF activity: effect of dropping the per-cell scaling") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5))
  ggsave("Plots/sarc_TF_scaled_vs_unscaled.pdf", pc, width = 11, height = 3.6)
  ggsave("Plots/sarc_TF_scaled_vs_unscaled.png", pc, width = 11, height = 3.6, dpi = 200)

  md <- merge(aggregate(score ~ TF, act[act$scaling == "scaled (original)", ], median),
              aggregate(score ~ TF, act[act$scaling == "unscaled activity", ], median),
              by = "TF", suffixes = c(".scaled", ".unscaled"))
  md$delta <- md$score.unscaled - md$score.scaled
  cat("\n=== per-TF median activity score, scaled vs unscaled ===\n")
  print(md[order(-abs(md$delta)), ], row.names = FALSE, digits = 3)
  cat("\nmedian |delta|:", round(median(abs(md$delta)), 3),
      " | sign flips:", sum(sign(md$score.scaled) != sign(md$score.unscaled)), "of", nrow(md), "\n")
  ## do the small samples still tie exactly?
  ties <- sum(sapply(split(act[act$scaling == "unscaled activity", ], act$TF[act$scaling == "unscaled activity"]),
                     function(x) any(duplicated(x$score))))
  cat("panels with an exact tie (overplotted dots), unscaled:", ties, "of", length(top_sarc_TF), "\n")
}
cat("\nDONE\n")
