# KO - NTC profiles for the correlation figures (fig3, fig4, fig5). Sourced after repro_common.R.
# needs: screen_export_figs.R (screen variable genes) and the earlier module-score outputs
#        screen_cm_diff_matrix.csv (screen, NTC - KO) and cm_sample_scores_AddModuleScore_bulk_RUN1.csv (bulk, per sample)

ORDER <- c("TCF3", "TEAD4", "TWIST1", "PITX1", "SOX9", "TEAD2", "BPTF", "HMGA1")   # KOs present in screen AND bulk, CROP-seq (fig1) order
EXTRA <- c("MEF2A", "MEF2D")                                                        # screen-only KOs

# returns, for the two feature sets ("mod" = 20 Cm modules, "gen" = screen's top variable genes present in bulk):
#   $avg[[dataset]] features x KO  (replicate-averaged KO - NTC;  dataset = "CROP-seq" or a cell line)
#   $rep[[cell line]] named list sample -> feature vector (single-replicate KO - NTC)
#   $meta  bulk sample metadata
load_profiles <- function(lines) {
  b <- load_bulk(); cp10k <- b$cp10k; meta <- b$meta
  scr_cm <- -as.matrix(read.csv(file.path(PP, "screen_cm_diff_matrix.csv"), row.names = 1, check.names = FALSE))   # stored NTC - KO -> KO - NTC
  cm_cols <- colnames(scr_cm)
  sc_df <- read.csv(file.path(PP, "cm_sample_scores_AddModuleScore_bulk_RUN1.csv"), check.names = FALSE)
  sc <- as.matrix(sc_df[, cm_cols]); rownames(sc) <- sc_df$sample
  scr_g <- as.matrix(read.csv(file.path(OUT, "screen_var5000_log2avg_by_KO.csv"), row.names = 1, check.names = FALSE))
  genes <- rownames(scr_g)[rownames(scr_g) %in% rownames(cp10k)]
  lg <- log2(cp10k[genes, , drop = FALSE] + 1)

  out <- list(mod = list(avg = list("CROP-seq" = t(scr_cm)), rep = list()),
              gen = list(avg = list("CROP-seq" = scr_g[genes, colnames(scr_g) != "NTC", drop = FALSE] - scr_g[genes, "NTC"]), rep = list()),
              meta = meta, genes = genes)
  for (cl in lines) {
    m <- meta[meta$Cell_Line == cl, ]
    ntc_s <- rownames(m)[m$Gene_Control == "NTC"]
    ko_s <- rownames(m)[m$Gene_Control != "NTC"]
    kos <- sort(setdiff(unique(m$Gene_Control), "NTC"), method = "radix")
    # modules
    ntc_m <- apply(sc[ntc_s, , drop = FALSE], 2, median)
    out$mod$avg[[cl]] <- sapply(setNames(kos, kos), function(g) colMeans(sc[rownames(m)[m$Gene_Control == g], , drop = FALSE]) - ntc_m)
    out$mod$rep[[cl]] <- setNames(lapply(ko_s, function(s) sc[s, ] - ntc_m), ko_s)
    # variable genes
    g_avg <- group_log2avg(cp10k, meta, cl)[genes, , drop = FALSE]
    out$gen$avg[[cl]] <- g_avg[, kos, drop = FALSE] - g_avg[, "NTC"]
    ntc_g <- log2(rowMeans(cp10k[genes, ntc_s, drop = FALSE]) + 1)
    out$gen$rep[[cl]] <- setNames(lapply(ko_s, function(s) lg[, s] - ntc_g), ko_s)
  }
  out
}

# cross-correlation of bulk KO profiles (rows: cell line x KO in ORDER) with CROP-seq KO profiles (columns: ORDER + EXTRA)
cross_corr <- function(N, lines, cols = c(ORDER, EXTRA)) {
  S <- N[["CROP-seq"]][, cols, drop = FALSE]; rows <- list(); lab <- list()
  for (cl in lines) for (ko in ORDER) if (ko %in% colnames(N[[cl]])) {
    rows[[length(rows) + 1]] <- sapply(cols, function(cc) cor(N[[cl]][, ko], S[, cc])); lab[[length(lab) + 1]] <- c(cl, ko)
  }
  lab <- as.data.frame(do.call(rbind, lab), stringsAsFactors = FALSE); names(lab) <- c("line", "ko")
  list(C = do.call(rbind, rows), lab = lab)
}

# two-panel heatmap: bulk KO (rows, grouped by cell line) x CROP-seq KO (columns); matched pair outlined
draw_cross_heatmaps <- function(name, title, sub_lines, res, panel_titles, row_labels, lines, vmax) {
  save_fig(name, 17.5, 9.2, function() {
    layout(matrix(1:2, nrow = 1)); par(oma = c(4.5, 0, 6, 0.5), family = "sans")
    for (p in 1:2) {
      C <- res[[p]]$C; lab <- res[[p]]$lab; n <- nrow(C); k <- ncol(C)
      par(mar = c(0.5, 16, 6.5, 0.5), xpd = FALSE)
      image(seq_len(k), seq_len(n), heat_z(pmax(pmin(C, vmax), -vmax)), col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
      axis(3, at = seq_len(k), labels = colnames(C), las = 2, tick = FALSE, cex.axis = 0.8, line = -0.4)
      axis(2, at = n:1, labels = row_labels[[p]], las = 1, tick = FALSE, cex.axis = 0.68, line = -0.4)
      for (i in seq_len(n)) for (j in seq_len(k)) {
        v <- C[i, j]; txt <- sub("^([+-])0\\.", "\\1.", sprintf("%+.2f", v))
        text(j, n - i + 1, txt, cex = 0.5, col = if (abs(v) < 0.6 * vmax) INK2 else "white")
      }
      for (i in seq_len(n)) { j <- match(lab$ko[i], colnames(C)); rect(j - 0.5, n - i + 0.5, j + 0.5, n - i + 1.5, border = INK, lwd = 1.6) }
      abline(v = length(ORDER) + 0.5, col = "white", lwd = 4)
      edges <- which(lab$line[-1] != lab$line[-n]); abline(h = n - edges + 0.5, col = "white", lwd = 4)
      bounds <- c(0, edges, n)
      par(xpd = NA)
      for (b in seq_along(lines)) {
        y1 <- n - bounds[b + 1] + 0.5; y2 <- n - bounds[b] + 0.5
        rect(-2.05, y1, -1.75, y2, col = LINE_COL[[lines[b]]], border = NA)
        text(-2.4, (y1 + y2) / 2, lines[b], srt = 90, adj = c(0.5, 0), font = 2, cex = 0.85, col = LINE_COL[[lines[b]]])
      }
      mtext(panel_titles[[p]][1], side = 3, adj = 0, line = 4.4, cex = 0.8); mtext(panel_titles[[p]][2], side = 3, adj = 0, line = 3.3, cex = 0.8)
    }
    mtext(title, outer = TRUE, side = 3, adj = 0, line = 4.3, font = 2, cex = 1.15)
    for (i in seq_along(sub_lines)) mtext(sub_lines[i], outer = TRUE, side = 3, adj = 0, line = 3 - (i - 1) * 0.95, cex = 0.62, col = INK2)
    draw_colorbar(c(0.36, 0.66), c(-vmax, vmax), "Pearson r between bulk KO profile and CROP-seq KO profile")
  })
}

cross_stats <- function(C, lab) {       # matched-vs-other summary of one cross-correlation matrix
  k <- length(ORDER); match_j <- match(lab$ko, ORDER)
  mt <- C[cbind(seq_len(nrow(C)), match_j)]
  msk <- matrix(TRUE, nrow(C), ncol(C)); msk[cbind(seq_len(nrow(C)), match_j)] <- FALSE; msk[, (k + 1):ncol(C)] <- FALSE
  ranks <- sapply(seq_len(nrow(C)), function(i) sum(C[i, seq_len(k)] > C[i, match_j[i]]) + 1)
  list(matched = mean(mt), other = mean(C[msk]), n_first = sum(ranks == 1), n = nrow(C), median_rank = median(ranks), ranks = ranks)
}
