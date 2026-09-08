###############################################################################
# STEP 7 -- what the one substantial disagreement actually is.
#
# Step 3 shows the agreement is near-perfect for every immune and stromal type and
# that essentially all of the loss is one cell type: 23.5% of the cells the
# chromatin annotation calls Malignant are transferred as Fibroblasts, and all 8
# LSI clusters whose majority labels disagree are Malignant clusters.
#
# That is not a random failure.  Sarcomatoid mesothelioma is spindle-cell and
# mesenchymal by definition; its cells are transcriptionally fibroblast-like, and a
# transfer trained on an RNA reference in which "Fibroblasts" are normal stromal
# fibroblasts has no other label to give them.  If that explanation is right, the
# per-tumour rate of Malignant -> Fibroblasts must track the tumour's sarcomatoid
# score, which was computed independently in R2_Q14 from the scATAC data.
#
# This is a falsifiable prediction about which tumours disagree, not a rescue of
# the agreement number: the disagreeing cells are still counted as disagreements
# everywhere else in this folder.
#
# Input : atac_joined.csv, R2_Q14/scATAC_per_sample_RUNX_BAP1_sarc.csv
# Output: malignant_disagreement_by_sample.csv, Plots/R2Q5_panel_G.pdf
###############################################################################
source("00_common.R")
suppressMessages({ library(ggplot2); library(ggrepel) })

SARC <- file.path(SC, "git_repo_claude", "R2_Q14",
                  "scATAC_per_sample_RUNX_BAP1_sarc.csv")
MINMAL <- 50    # tumours with fewer malignant cells are shown but not correlated

X <- read.csv("atac_joined.csv", stringsAsFactors = FALSE)
M <- X[X$atac_label == "Malignant", ]

D <- do.call(rbind, lapply(split(M, M$sample), function(d)
  data.frame(sample = d$sample[1], n_malignant = nrow(d),
             pct_kept_malignant = 100 * mean(d$predictedGroup == "Malignant"),
             pct_to_fibroblast  = 100 * mean(d$predictedGroup == "Fibroblasts"),
             pct_to_lymphoid    = 100 * mean(d$predictedGroup %in%
                                             c("T_cells","NK","B_cells","Plasma")),
             median_score = median(d$predictedScore), stringsAsFactors = FALSE)))

sc <- read.csv(SARC, stringsAsFactors = FALSE)
D  <- merge(D, sc[, c("sample", "sarc_score", "BAP1_status")], by = "sample", all.x = TRUE)
D  <- D[order(-D$pct_to_fibroblast), ]
write.csv(D, "malignant_disagreement_by_sample.csv", row.names = FALSE)
cat("=== malignant cells, per tumour: where does the transfer send them? ===\n")
print(transform(D, pct_kept_malignant = round(pct_kept_malignant, 1),
                pct_to_fibroblast = round(pct_to_fibroblast, 1),
                pct_to_lymphoid = round(pct_to_lymphoid, 1),
                median_score = round(median_score, 3),
                sarc_score = round(sarc_score, 3)), row.names = FALSE)

d <- D[!is.na(D$sarc_score), ]
dm <- d[d$n_malignant >= MINMAL, ]
cat(sprintf("\nsarcomatoid score vs %% transferred to Fibroblasts:\n"))
cat(sprintf("  all %d tumours with a score : Spearman rho = %.3f, p = %.3g\n",
            nrow(d), cor(d$sarc_score, d$pct_to_fibroblast, method = "spearman"),
            cor.test(d$sarc_score, d$pct_to_fibroblast, method = "spearman",
                     exact = FALSE)$p.value))
cat(sprintf("  the %d with >= %d malignant cells : Spearman rho = %.3f, p = %.3g\n",
            nrow(dm), MINMAL, cor(dm$sarc_score, dm$pct_to_fibroblast, method = "spearman"),
            cor.test(dm$sarc_score, dm$pct_to_fibroblast, method = "spearman",
                     exact = FALSE)$p.value))
cat("  (P10 and P23 have no sarcomatoid score in R2_Q14 and are excluded from the test)\n")

p <- ggplot(d, aes(sarc_score, pct_to_fibroblast)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = .3, colour = "grey60") +
  geom_point(aes(size = n_malignant), colour = "plum4", alpha = .85) +
  geom_text_repel(aes(label = sample), size = 2.6, min.segment.length = 0, seed = 1) +
  scale_size_continuous(range = c(1, 5), name = "malignant\ncells") +
  labs(x = "sarcomatoid score (scATAC, R2_Q14)",
       y = "% of chromatin-defined malignant cells\ntransferred as Fibroblasts",
       title = sprintf("the malignant disagreement tracks sarcomatoid histology (rho = %.2f)",
                       cor(d$sarc_score, d$pct_to_fibroblast, method = "spearman"))) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), plot.title = element_text(size = 8),
        legend.key.size = unit(3, "mm"))
ggsave("Plots/R2Q5_panel_G.pdf", p, width = 5, height = 4, device = cairo_pdf)

cat("\nDONE\n")
