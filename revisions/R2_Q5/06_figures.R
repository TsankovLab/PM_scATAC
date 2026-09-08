###############################################################################
# STEP 6 -- the figure.
#
#   A,B  the scATAC UMAP, coloured by the chromatin annotation (A) and by the
#        transcriptomic annotation transferred onto the same cells (B).  Same
#        embedding, same cells, two independent annotations -- the visual form of
#        the question.
#   C    confusion matrix, row-normalised (chromatin label -> transferred RNA label).
#   D    recall per chromatin-defined cell type, with the cell count.
#   E    agreement as a function of the transfer's own confidence.
#   F    RNA marker sets scored in ATAC chromatin -- the integration-free test.
#   G    per-patient composition, chromatin vs transcriptomic annotation.
#   H    the one substantial disagreement (malignant -> fibroblast) against the
#        independently computed sarcomatoid score of each tumour (step 7).
#
# Input : atac_joined.csv, concordance_*.csv, confusion_matrix.csv,
#         marker_enrichment.csv, composition_by_patient.csv
# Output: Plots/R2Q5_annotation_concordance.pdf and the panels separately
###############################################################################
source("00_common.R")
suppressMessages({ library(ggplot2); library(patchwork); library(ggrepel) })

th <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        legend.key.size = unit(3, "mm"))

X  <- read.csv("atac_joined.csv", stringsAsFactors = FALSE)
S  <- read.csv("concordance_summary.csv", stringsAsFactors = FALSE)
CT <- read.csv("concordance_by_celltype.csv", stringsAsFactors = FALSE)
SB <- read.csv("concordance_by_score.csv", stringsAsFactors = FALSE)
CF <- as.matrix(read.csv("confusion_matrix.csv", row.names = 1, check.names = FALSE))
ME <- as.matrix(read.csv("marker_enrichment.csv", row.names = 1, check.names = FALSE))
CP <- read.csv("composition_by_patient.csv", stringsAsFactors = FALSE)
MD <- read.csv("malignant_disagreement_by_sample.csv", stringsAsFactors = FALSE)

lv <- CELLTYPES[CELLTYPES %in% union(X$atac_label, X$predictedGroup)]
ov <- S$agreement[S$subset == "all ATAC cells"]
kp <- S$kappa[S$subset == "all ATAC cells"]

## ---- A. the two annotations on the same embedding ------------------------------
umap <- function(col, title){
  d <- X[sample(nrow(X)), ]
  ggplot(d, aes(UMAP1, UMAP2, colour = factor(.data[[col]], lv))) +
    geom_point(size = .05, stroke = 0, alpha = .6) +
    scale_colour_manual(values = CTCOL, drop = FALSE, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1), ncol = 1)) +
    labs(title = title) + th +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), plot.title = element_text(size = 8))
}
pA1 <- umap("atac_label", "chromatin annotation (scATAC clusters)")
pA2 <- umap("predictedGroup", "transcriptomic annotation (transferred)") +
       theme(legend.position = "none")

## ---- B. confusion matrix ---------------------------------------------------------
P  <- 100 * prop.table(CF, 1)
dB <- as.data.frame.table(P, responseName = "pct", stringsAsFactors = FALSE)
names(dB)[1:2] <- c("atac", "rna")
pB <- ggplot(dB, aes(factor(rna, lv), factor(atac, rev(lv)), fill = pct)) +
  geom_tile(colour = "white", linewidth = .3) +
  geom_text(data = subset(dB, pct >= 5), aes(label = round(pct)), size = 2) +
  scale_fill_gradient(low = "white", high = "#3b6ea5", limits = c(0, 100),
                      name = "% of row") +
  labs(x = "transcriptomic label (transferred)", y = "chromatin label",
       title = sprintf("agreement %.1f%%, kappa %.2f", 100 * ov, kp)) + th +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 8))

## ---- C. recall per cell type ------------------------------------------------------
CT$lab <- sprintf("%s (%s)", CT$celltype, format(CT$n_atac, big.mark = ","))
CT <- CT[order(CT$recall), ]
pC <- ggplot(CT, aes(recall, factor(lab, CT$lab), fill = celltype)) +
  geom_col(width = .7) +
  geom_vline(xintercept = ov, linetype = 2, linewidth = .3, colour = "grey30") +
  scale_fill_manual(values = CTCOL, guide = "none") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, .02))) +
  labs(x = "fraction of the chromatin-defined cells given the same RNA label",
       y = NULL, title = "recall per cell type (dashed = cohort agreement)") + th +
  theme(plot.title = element_text(size = 8))

## ---- D. agreement vs transfer confidence -------------------------------------------
SB$score_bin <- factor(SB$score_bin, SB$score_bin)
pD <- ggplot(SB, aes(score_bin, agreement)) +
  geom_col(aes(y = 1), fill = "grey92", width = .75) +
  geom_col(fill = "#3b6ea5", width = .75) +
  geom_text(aes(label = sprintf("%.0f%%\nof cells", pct_of_cells)), y = .04,
            size = 1.9, colour = "grey25", vjust = 0) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, .02))) +
  labs(x = "transfer confidence (predictedScore)", y = "agreement",
       title = "disagreement is concentrated at low transfer confidence") + th +
  theme(plot.title = element_text(size = 8))

## ---- E. RNA marker sets scored in chromatin -----------------------------------------
dE <- as.data.frame.table(ME, responseName = "z", stringsAsFactors = FALSE)
names(dE)[1:2] <- c("atac", "markers")
le <- intersect(lv, unique(dE$atac))
lim <- max(abs(dE$z))
pE <- ggplot(dE, aes(factor(markers, le), factor(atac, rev(le)), fill = z)) +
  geom_tile(colour = "white", linewidth = .3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       limits = c(-lim, lim), name = "mean z\ngene score") +
  labs(x = "RNA marker set (top 100 genes)", y = "chromatin label",
       title = "no integration: RNA markers scored in ATAC chromatin") + th +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 8))

## ---- F. composition -------------------------------------------------------------------
r <- cor(CP$frac_atac, CP$frac_rna)
pF <- ggplot(CP, aes(frac_rna, frac_atac, colour = factor(celltype, CELLTYPES))) +
  geom_abline(slope = 1, linetype = 2, linewidth = .3, colour = "grey50") +
  geom_point(size = 1.1, alpha = .85) +
  scale_colour_manual(values = CTCOL, drop = FALSE, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 1.6), ncol = 1)) +
  labs(x = "fraction of cells, transcriptomic annotation",
       y = "fraction of cells, chromatin annotation",
       title = sprintf("composition of the 10 shared patients (r = %.3f)", r)) + th +
  theme(plot.title = element_text(size = 8))

## ---- G. the malignant disagreement vs sarcomatoid histology --------------------
md  <- MD[!is.na(MD$sarc_score), ]
rho <- cor(md$sarc_score, md$pct_to_fibroblast, method = "spearman")
pG <- ggplot(md, aes(sarc_score, pct_to_fibroblast)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .3, colour = "grey60",
              formula = y ~ x) +
  geom_point(aes(size = n_malignant), colour = "plum4", alpha = .85) +
  geom_text_repel(aes(label = sample), size = 2.6, min.segment.length = 0, seed = 1) +
  scale_size_continuous(range = c(1, 5), name = "malignant\ncells") +
  labs(x = "sarcomatoid score (scATAC, R2_Q14)",
       y = "% of malignant cells\ntransferred as Fibroblasts",
       title = sprintf("the malignant disagreement tracks sarcomatoid histology (rho = %.2f)", rho)) +
  th + theme(plot.title = element_text(size = 8))

## ---- assemble -----------------------------------------------------------------------
fig <- ((pA1 | pA2) / (pB | pC) / (pD | pE) / (pF | pG)) +
  plot_annotation(tag_levels = "A",
    title = "R2_Q5  Chromatin-based and transcriptome-based cell type annotation of the same tumours",
    subtitle = sprintf("%s scATAC cells, 11 patients | %s scRNA cells, 10 patients | unconstrained CCA label transfer",
                       format(nrow(X), big.mark = ","),
                       format(nrow(read.csv("rna_cells.csv")), big.mark = ",")),
    theme = theme(plot.title = element_text(size = 10, face = "bold"),
                  plot.subtitle = element_text(size = 8, colour = "grey30")))
ggsave("Plots/R2Q5_annotation_concordance.pdf", fig, width = 11, height = 15,
       device = cairo_pdf, limitsize = FALSE)

for (n in c("pA1","pA2","pB","pC","pD","pE","pF","pG"))
  ggsave(sprintf("Plots/R2Q5_panel_%s.pdf", sub("^p", "", n)), get(n),
         width = 5, height = 4, device = cairo_pdf)

cat("wrote Plots/R2Q5_annotation_concordance.pdf and 8 panels\n")
cat("DONE\n")
