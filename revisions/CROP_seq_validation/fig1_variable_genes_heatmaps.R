# Fig 1: KO pattern of the CROP-seq screen's top-5000 variable genes, CROP-seq vs bulk RNA-seq (RUN1)
# needs: screen_export_figs.R
# outputs: fig1_KO_pattern_top5000_variable_genes (stacked, all 4 lines)
#          fig1_side_by_side_CROPseq_vs_3lines_RUN1 (genes as rows, 3 lines, all KOs)
#          fig1_side_by_side_CROPseq_vs_3lines_RUN1_matchedKOs (bulk restricted to screen KOs, screen column order)
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "repro_common.R"))

scr   <- as.matrix(read.csv(file.path(OUT, "screen_var5000_log2avg_by_KO.csv"), row.names = 1, check.names = FALSE))
ncell <- read.csv(file.path(OUT, "screen_ncells_by_KO.csv")); ncell <- setNames(ncell$Freq, ncell$Var1)
b <- load_bulk(); cp10k <- b$cp10k; meta <- b$meta
cat("screen genes:", nrow(scr), "| in bulk:", sum(rownames(scr) %in% rownames(cp10k)), "\n")

zs <- zscore_by_gene(scr); zs <- zs[rowSums(abs(zs)) > 0, , drop = FALSE]
genes_all <- rownames(zs)[corr_order(zs)]                        # gene order = screen clustering
genes_bulk <- genes_all[genes_all %in% rownames(cp10k)]

# NTC first, then KOs clustered on their gene profile (rows of t(z))
order_ko <- function(z) { ko <- setdiff(colnames(z), "NTC"); c("NTC", ko[corr_order(t(z[, ko, drop = FALSE]))]) }
write.csv(data.frame(gene = genes_all), file.path(OUT, "fig1_gene_order_screen.csv"), row.names = FALSE)

## ------------------------------------------------------------ stacked version (rows = KO, columns = genes)
stacked <- function() {
  panels <- list(list(title = "CROP-seq screen", col = LINE_COL[["CROP-seq"]], z = zs[genes_all, ], n = ncell, lab = "cells"))
  for (cl in LINES) {
    g <- reindex_rows(group_log2avg(cp10k, meta, cl), genes_all)          # genes absent from bulk -> NA
    n_abs <- sum(is.na(g[, 1]))
    z <- (g - rowMeans(g, na.rm = TRUE)) / rowSds(g, na.rm = TRUE); z[!is.finite(z)] <- 0
    panels[[length(panels) + 1]] <- list(title = sprintf("%s bulk (RUN1)  [%d genes absent]", cl, n_abs), col = LINE_COL[[cl]], z = z,
                                         n = table(meta$Gene_Control[meta$Cell_Line == cl]), lab = "reps")
  }
  for (i in seq_along(panels)) { o <- order_ko(panels[[i]]$z); panels[[i]]$z <- panels[[i]]$z[, o, drop = FALSE]; panels[[i]]$n <- panels[[i]]$n[o] }
  save_fig("fig1_KO_pattern_top5000_variable_genes", 14, 0.24 * sum(sapply(panels, function(p) ncol(p$z))) + 4, function() {
    np <- length(panels); nk <- sapply(panels, function(p) ncol(p$z))
    layout(matrix(1:(2 * np), ncol = 2, byrow = TRUE), widths = c(40, 1.6), heights = nk + 2.2)
    par(oma = c(4, 0, 5, 1), family = "sans")
    for (p in panels) {
      k <- ncol(p$z)
      par(mar = c(0.3, 5.5, 1.6, 0.3))
      image(seq_len(nrow(p$z)), seq_len(k), pmax(pmin(p$z, 2.5), -2.5)[, k:1, drop = FALSE], col = DIV, zlim = c(-2.5, 2.5), axes = FALSE, xlab = "", ylab = "")   # z[gene, KO]; KO axis flipped so the first KO is on top
      axis(2, at = seq_len(k), labels = rev(colnames(p$z)), las = 1, tick = FALSE, cex.axis = 0.85, line = -0.6)
      abline(h = k - 0.5, col = "white", lwd = 3); box(col = MUTED, lwd = 0.6)
      mtext(p$title, side = 3, adj = 0, line = 0.2, font = 2, cex = 0.85, col = p$col)
      par(mar = c(0.3, 0.3, 1.6, 1))
      barplot(rev(as.numeric(p$n)), horiz = TRUE, col = p$col, border = NA, axes = FALSE, space = 0.3)
      axis(1, cex.axis = 0.6, lwd = 0.5, mgp = c(1, 0.3, 0)); mtext(p$lab, side = 1, line = 1.2, cex = 0.6, col = INK2)
    }
    mtext("KO pattern of the top-5000 variable genes: CROP-seq screen vs bulk RNA-seq", outer = TRUE, side = 3, adj = 0, line = 3, font = 2, cex = 1.15)
    mtext("columns = the screen's 5000 variable genes, same order in every panel (clustered on the screen); rows = KO groups (NTC on top, others clustered per panel); values = log2(mean CP10K+1), z-scored per gene across KO groups",
          outer = TRUE, side = 3, adj = 0, line = 1.6, cex = 0.6, col = INK2)
    draw_colorbar(c(0.10, 0.35), c(-2.5, 2.5), "gene z-score across KO groups")
  })
}

## ------------------------------------------------------------ side-by-side versions (rows = genes, columns = KO)
side_by_side <- function(name, matched, title, sub_lines) {
  scr_z <- zs[genes_bulk, , drop = FALSE]
  scr_z <- scr_z[, order_ko(scr_z), drop = FALSE]
  panels <- list(list(title = "CROP-seq screen", col = LINE_COL[["CROP-seq"]], z = scr_z, n = ncell))
  for (cl in THREE) {
    g <- group_log2avg(cp10k, meta, cl)[genes_bulk, , drop = FALSE]
    if (matched) g <- g[, intersect(colnames(scr_z), colnames(g)), drop = FALSE]
    z <- zscore_by_gene(g)
    z <- z[, if (matched) colnames(g) else order_ko(z), drop = FALSE]
    panels[[length(panels) + 1]] <- list(title = paste0(cl, " bulk (RUN1)"), col = LINE_COL[[cl]], z = z, n = table(meta$Gene_Control[meta$Cell_Line == cl]))
  }
  w <- sapply(panels, function(p) ncol(p$z))
  save_fig(name, 13.5, 11, function() {
    layout(matrix(c(1, 2, 3, 4, 5), nrow = 1), widths = c(w[1], 2.2, w[2:4]))
    par(oma = c(4.5, 0, 6, 0.5), family = "sans")
    for (i in seq_along(panels)) {
      p <- panels[[i]]
      par(mar = c(0.5, if (i == 1) 3 else 0.3, 7.5, 0.3))
      image(seq_len(ncol(p$z)), seq_len(nrow(p$z)), heat_z(pmax(pmin(p$z, 2.5), -2.5)),
            col = DIV, zlim = c(-2.5, 2.5), axes = FALSE, xlab = "", ylab = "")
      mtext(sprintf("%s (%d)", colnames(p$z), as.integer(p$n[colnames(p$z)])), side = 3, at = seq_len(ncol(p$z)), las = 2, cex = 0.62, line = 0.2, adj = 0)
      abline(v = 1.5, col = "white", lwd = 3); box(col = MUTED, lwd = 0.6)
      mtext(p$title, side = 3, adj = 0, line = 6.3, font = 2, cex = 0.85, col = p$col)
      if (i == 1) mtext(sprintf("%d top variable genes (screen order)", nrow(p$z)), side = 2, line = 1, cex = 0.7, col = INK2)
      if (i == 1) { plot.new() ; }                # gap panel
    }
    mtext(title, outer = TRUE, side = 3, adj = 0, line = 4.4, font = 2, cex = 1.05)
    for (k in seq_along(sub_lines)) mtext(sub_lines[k], outer = TRUE, side = 3, adj = 0, line = 3.3 - (k - 1) * 0.9, cex = 0.55, col = INK2)
    draw_colorbar(c(0.05, 0.25), c(-2.5, 2.5), "gene z-score across KO groups")
  })
}

stacked()
side_by_side("fig1_side_by_side_CROPseq_vs_3lines_RUN1", FALSE,
             "KO pattern of the top variable genes: CROP-seq (left) vs bulk RNA-seq, three cell lines, RUN1 (right)",
             c("rows = the screen's 5000 variable genes present in bulk, same order in every panel (clustered on the CROP-seq screen); columns = KO groups, NTC first, others clustered per panel;",
               "values = log2(mean CP10K + 1) per KO group, z-scored per gene across KO groups; number in brackets = cells (screen) or replicates (bulk). NCI-H28 excluded; outlier samples removed (RUN1)."))
side_by_side("fig1_side_by_side_CROPseq_vs_3lines_RUN1_matchedKOs", TRUE,
             "KO pattern of the top variable genes: CROP-seq (left) vs bulk RNA-seq, three cell lines, RUN1 (right)",
             c("rows = the screen's 5000 variable genes present in bulk, same order in every panel (clustered on the CROP-seq screen)",
               "columns = KO groups present in the screen, in the CROP-seq panel's order (bulk has no MEF2A/MEF2D); values = log2(mean CP10K + 1), z-scored per gene across the KO groups shown",
               "number in brackets = cells (screen) or replicates (bulk). NCI-H28 excluded; outlier samples removed (RUN1)."))
