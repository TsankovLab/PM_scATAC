# Spec 26: does our bulk RNA-seq proliferation index (cc = S.Score + G2M.Score, raw KO score -- NOT subtracted vs NTC) agree with DepMap's actual CRISPR viability
# screen for the same gene in the same cell line? DepMap CRISPRGeneEffect.csv (Chronos gene-effect score; negative = knockout suppresses
# growth/fitness, positive = knockout promotes growth) downloaded from the public DepMap release (Figshare, DepMap_24Q4_Public,
# article 27993248, file Model.csv + CRISPRGeneEffect.csv) -- extracted here to depmap_crispr_chronos_meso4lines.csv (checked into this
# repo, 4 rows x 11 genes) via the 4 ACH IDs matched by name in Model.csv: MSTO-211H=ACH-000335, NCI-H2052=ACH-000153 (matches the value
# already hardcoded in spec_common.R's SCREEN_LINE comment -- good cross-check), NCI-H2452=ACH-000092, NCI-H28=ACH-000648.
# All 11 bulk KOs compared (the 8 CROP-seq-shared KOs + SNAI2/RUNX1/RUNX2, all present in DepMap).
# needs: fig2_proliferation.R (fig2_bulk_cc_deltas.csv), depmap_crispr_chronos_meso4lines.csv (this dir)
# outputs: spec26_depmap_merged.csv, spec26_depmap_correlation_summary.csv, fig2G_bulk_proliferation_vs_depmap_viability.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R"))

LINES4 <- c("MSTO-211H", "NCI-H2052", "NCI-H2452", "NCI-H28")
bd <- read.csv(file.path(OUT, "fig2_bulk_cc_deltas.csv"))
dm <- read.csv(file.path(here, "depmap_crispr_chronos_meso4lines.csv"))
dm_long <- do.call(rbind, lapply(setdiff(colnames(dm), "cell_line"), function(g) data.frame(cell_line = dm$cell_line, KO = g, chronos = dm[[g]])))

m <- merge(bd, dm_long, by = c("cell_line", "KO")); m$cell_line <- factor(m$cell_line, LINES4); m <- m[order(m$cell_line, m$KO), ]
write.csv(m, file.path(OUT, "spec26_depmap_merged.csv"), row.names = FALSE)
cat(sprintf("merged: %d cell-line x KO pairs (%d KOs x %d lines)\n", nrow(m), length(unique(m$KO)), length(unique(m$cell_line))))

summ <- list()
for (cl in LINES4) {
  d <- m[m$cell_line == cl, ]; ct <- cor.test(d$mean_cc, d$chronos, method = "spearman", exact = FALSE)
  summ[[cl]] <- data.frame(cell_line = cl, n_KO = nrow(d), spearman_rho = unname(ct$estimate), p = ct$p.value, sign_agreement = mean(sign(d$mean_cc) == sign(d$chronos)))
}
ct <- cor.test(m$mean_cc, m$chronos, method = "spearman", exact = FALSE)
summ[["ALL (pooled)"]] <- data.frame(cell_line = "ALL (pooled)", n_KO = nrow(m), spearman_rho = unname(ct$estimate), p = ct$p.value, sign_agreement = mean(sign(m$mean_cc) == sign(m$chronos)))
summ <- do.call(rbind, summ); rownames(summ) <- NULL
write.csv(summ, file.path(OUT, "spec26_depmap_correlation_summary.csv"), row.names = FALSE)
cat("\nbulk proliferation index (raw KO mean_cc, NOT subtracted vs NTC) vs DepMap Chronos gene-effect score (negative = essential/growth-suppressing knockout):\n")
print(format(summ, digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- figure: scatter per line + pooled, same layout as fig2B
save_fig("fig2G_bulk_proliferation_vs_depmap_viability", 21, 4.6, function() {
  par(mfrow = c(1, 5), mar = c(4, 4.2, 4.2, 0.8), oma = c(0, 0, 2.2, 0), family = "sans")
  for (cl in LINES4) {
    d <- m[m$cell_line == cl, ]; s <- summ[summ$cell_line == cl, ]
    plot(d$chronos, d$mean_cc, pch = 21, bg = LINE_COL[[cl]], col = "white", cex = 1.8, bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8,
         xlab = "DepMap Chronos score (KO viability effect)", ylab = "bulk KO cell-cycle index (raw)", xlim = range(d$chronos) + c(-0.05, 0.12))
    abline(v = 0, col = GRID); text(d$chronos, d$mean_cc, d$KO, pos = 4, cex = 0.7, offset = 0.4)
    mtext(cl, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.8, col = LINE_COL[[cl]])
    mtext(sprintf("Spearman ρ=%+.2f (p=%.2f) · sign agreement %.0f%%", s$spearman_rho, s$p, 100 * s$sign_agreement), side = 3, adj = 0, line = 0.7, cex = 0.7, font = 2, col = LINE_COL[[cl]])
  }
  s <- summ[summ$cell_line == "ALL (pooled)", ]
  plot(m$chronos, m$mean_cc, pch = 21, bg = LINE_COL[as.character(m$cell_line)], col = "white", cex = 1.4, bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8,
       xlab = "DepMap Chronos score", ylab = ""); abline(v = 0, col = GRID)
  legend("topleft", LINES4, pt.bg = LINE_COL[LINES4], pch = 21, col = "white", bty = "n", cex = 0.7)
  mtext("all lines pooled", side = 3, adj = 0, line = 1.9, font = 2, cex = 0.8)
  mtext(sprintf("ρ=%+.2f (p=%.2f)", s$spearman_rho, s$p), side = 3, adj = 0, line = 0.7, cex = 0.7, font = 2)
  mtext("Does our bulk RNA-seq proliferation index agree with DepMap's actual CRISPR viability screen?", outer = TRUE, side = 3, adj = 0, font = 2, cex = 1)
  mtext("x: DepMap Chronos gene-effect score (negative = knockout suppresses growth); y: our bulk KO cell-cycle index, raw (NOT subtracted vs NTC), replicate-averaged; all 11 bulk KOs (incl. SNAI2/RUNX1/RUNX2)",
        outer = TRUE, side = 3, adj = 0, line = -1.1, cex = 0.55, col = INK2)
})
