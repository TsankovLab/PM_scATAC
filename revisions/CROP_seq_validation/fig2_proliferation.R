# Fig 2: proliferation index (cc = S.Score + G2M.Score) by KO, CROP-seq vs bulk RNA-seq (RUN1)
# needs: screen_export_figs.R, bulk_cellcycle_RUN1.R
# outputs: fig2_screen_cc_stats.csv, fig2_bulk_cc_deltas.csv, fig2B_summary.csv
#          fig2A_cell_cycle_index_by_KO            (screen violins + bulk bars/dots, 4 lines)
#          fig2B_cell_cycle_bulk_vs_screen         (bulk delta vs screen delta, shared KOs)
#          fig2C_proliferation_boxplots_CROPseq_vs_bulk_RUN1  (CROP-seq style boxplots, 3 lines, all bulk KOs)
#          fig2D_proliferation_side_by_side_CROPseq_vs_3lines_RUN1 (fig1 layout: shared KOs only, screen order)
#          fig2E_proliferation_combined_zscored_by_KO (all 4 bulk lines + screen combined into one box per KO, Δ vs NTC, per-source z-scored)
#          fig2F_proliferation_boxplots_bulk_all_lines_pooled (bulk only, one point per line (replicates averaged), raw index, no NTC subtraction, no CROP-seq, CROP-seq screen order)
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "repro_common.R"))
set.seed(1)

scr  <- read.csv(file.path(OUT, "screen_cellcycle_per_cell.csv"), row.names = 1)
bulk <- read.csv(file.path(OUT, "bulk_cellcycle_per_sample_RUN1.csv"))
bulk_full <- bulk     # unrestricted (all 4 lines, all bulk KOs) -- kept aside since `bulk` gets restricted to screen KOs further down for fig2D
LINES4 <- c("MSTO-211H", "NCI-H2052", "NCI-H2452", "NCI-H28")

## ---------------------------------------------------------------- screen: Wilcoxon vs NTC, BH across KOs
ntc <- scr$cc[scr$merged_call == "NTC"]
kos <- setdiff(sort(unique(scr$merged_call), method = "radix"), "NTC")
st <- do.call(rbind, lapply(kos, function(g) {
  x <- scr$cc[scr$merged_call == g]
  data.frame(KO = g, n = length(x), mean_cc = mean(x), delta_vs_NTC = mean(x) - mean(ntc), p = wilcox.test(x, ntc, exact = FALSE)$p.value)
}))
st$padj  <- p.adjust(st$p, "BH")
st$stars <- ifelse(st$padj < 0.001, "***", ifelse(st$padj < 0.01, "**", ifelse(st$padj < 0.05, "*", "")))
write.csv(st, file.path(OUT, "fig2_screen_cc_stats.csv"), row.names = FALSE)
print(format(st, digits = 4), row.names = FALSE)

## ---------------------------------------------------------------- bulk deltas vs NTC
bd <- do.call(rbind, lapply(LINES4, function(cl) {
  d <- bulk[bulk$Cell_Line == cl, ]; nm <- mean(d$cc[d$Gene_Control == "NTC"])
  do.call(rbind, lapply(sort(setdiff(unique(d$Gene_Control), "NTC"), method = "radix"), function(g) {
    v <- d$cc[d$Gene_Control == g]
    data.frame(cell_line = cl, KO = g, n = length(v), mean_cc = mean(v), delta_vs_NTC = mean(v) - nm)
  }))
}))
write.csv(bd, file.path(OUT, "fig2_bulk_cc_deltas.csv"), row.names = FALSE)

## ---------------------------------------------------------------- panel helpers
frame_panel <- function(k, yl, labels) {
  plot(NA, xlim = c(0.4, k + 0.6), ylim = yl, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "", xaxs = "i")
  abline(h = pretty(yl), col = GRID, lwd = 0.7)
  axis(2, las = 1, cex.axis = 0.7, col = GRID, col.axis = INK2, lwd = 0.7)
  axis(1, at = seq_len(k), labels = labels, las = 2, tick = FALSE, cex.axis = 0.62, col.axis = INK2, line = -0.4)
  abline(h = par("usr")[3], col = GRID)
}
ylab_cc <- function() mtext("cell-cycle index (S.Score + G2M.Score)", side = 2, line = 2.6, cex = 0.65)

scr_violin_panel <- function(order, title) {
  data <- lapply(order, function(g) scr$cc[scr$merged_call == g]); k <- length(order)
  ytop <- max(sapply(data, function(d) quantile(d, 0.997)))
  frame_panel(k, c(min(sapply(data, min)), ytop + 0.16), sprintf("%s (%d)", order, sapply(data, length)))
  for (i in seq_len(k)) violin_at(data[[i]], i, if (order[i] == "NTC") GREY else BROWN, 0.47)
  boxplot(data, at = seq_len(k), add = TRUE, boxwex = 0.14, outline = FALSE, axes = FALSE,
          col = adjustcolor("white", 0.7), border = INK2, medcol = INK, medlwd = 1.4, lwd = 0.8, whisklty = 1, staplelty = 1)
  for (i in seq_len(k)) if (order[i] != "NTC") text(i, ytop + 0.06, st$stars[st$KO == order[i]], cex = 0.9)
  abline(h = mean(ntc), col = MUTED, lty = 2, lwd = 0.9)
  ylab_cc(); mtext(title, side = 3, adj = 0, line = 0.6, font = 2, cex = 0.8, col = INK)
}
bulk_box_panel <- function(cl, order, title, col, ylab) {
  d <- bulk[bulk$Cell_Line == cl, ]; data <- lapply(order, function(g) d$cc[d$Gene_Control == g]); k <- length(order)
  r <- range(unlist(data))
  frame_panel(k, r + c(-1, 1) * 0.06 * diff(r), sprintf("%s (%d)", order, sapply(data, length)))
  boxplot(data, at = seq_len(k), add = TRUE, boxwex = 0.55, outline = FALSE, axes = FALSE,
          col = adjustcolor(ifelse(order == "NTC", GREY, BROWN), 0.6), border = INK2, medcol = INK, medlwd = 1.4, lwd = 0.8, whisklty = 1, staplelty = 1)
  for (i in seq_len(k)) points(i + runif(length(data[[i]]), -0.13, 0.13), data[[i]], pch = 21, bg = INK, col = "white", cex = 0.85)
  abline(h = mean(d$cc[d$Gene_Control == "NTC"]), col = MUTED, lty = 2, lwd = 0.9)
  if (ylab) ylab_cc()
  mtext(title, side = 3, adj = 0, line = 0.6, font = 2, cex = 0.8, col = col)
}
bulk_bar_panel <- function(cl, order) {
  d <- bulk[bulk$Cell_Line == cl, ]; data <- lapply(order, function(g) d$cc[d$Gene_Control == g]); k <- length(order)
  yl <- range(c(0, unlist(data))); yl <- yl + c(-1, 1) * 0.05 * diff(yl)
  frame_panel(k, yl, sprintf("%s (%d)", order, sapply(data, length)))
  for (i in seq_len(k)) {
    rect(i - 0.35, 0, i + 0.35, mean(data[[i]]), col = adjustcolor(if (order[i] == "NTC") GREY else LINE_COL[[cl]], 0.35), border = NA)
    points(i + seq(-0.12, 0.12, length.out = length(data[[i]])) * (length(data[[i]]) > 1), data[[i]], pch = 21, bg = INK, col = "white", cex = 0.85)
  }
  abline(h = mean(d$cc[d$Gene_Control == "NTC"]), col = MUTED, lty = 2, lwd = 0.9)
  mtext(paste0(cl, " bulk (RUN1)"), side = 3, adj = 0, line = 0.6, font = 2, cex = 0.8, col = LINE_COL[[cl]])
}
order_mean <- function(v) names(sort(v, decreasing = TRUE))                    # groups sorted by mean, high -> low
scr_means  <- tapply(scr$cc, scr$merged_call, mean)
bulk_order <- function(cl) { d <- bulk[bulk$Cell_Line == cl, ]; order_mean(tapply(d$cc, d$Gene_Control, mean)) }

## ---------------------------------------------------------------- Fig 2A
save_fig("fig2A_cell_cycle_index_by_KO", 21, 5, function() {
  so <- order_mean(scr_means); ords <- lapply(LINES4, bulk_order)
  layout(matrix(1:5, nrow = 1), widths = c(length(so) * 1.3, sapply(ords, length)))
  par(oma = c(0, 0, 3, 0), mar = c(8.5, 4.2, 2, 0.8), family = "sans")
  scr_violin_panel(so, "CROP-seq screen (single cells)")
  for (i in 1:4) bulk_bar_panel(LINES4[i], ords[[i]])
  mtext("Cell-cycle index by KO relative to NTC (sorted by mean; dashed = NTC mean)", outer = TRUE, side = 3, adj = 0, font = 2, cex = 1)
})

## ---------------------------------------------------------------- Fig 2B
m <- merge(bd, data.frame(KO = st$KO, screen_delta = st$delta_vs_NTC), by = "KO")
m$cell_line <- factor(m$cell_line, LINES4); m <- m[order(m$cell_line, m$KO), ]
summ <- list()
save_fig("fig2B_cell_cycle_bulk_vs_screen", 21, 4.6, function() {
  par(mfrow = c(1, 5), mar = c(4, 4.2, 4.2, 0.8), oma = c(0, 0, 2.2, 0), family = "sans")
  for (cl in LINES4) {
    d <- m[m$cell_line == cl, ]; ct <- cor.test(d$screen_delta, d$delta_vs_NTC, method = "spearman", exact = FALSE)
    ag <- mean(sign(d$screen_delta) == sign(d$delta_vs_NTC))
    summ[[cl]] <<- data.frame(cell_line = cl, n_KO = nrow(d), spearman_rho = unname(ct$estimate), p = ct$p.value, sign_agreement = ag)
    plot(d$screen_delta, d$delta_vs_NTC, pch = 21, bg = LINE_COL[[cl]], col = "white", cex = 1.8, bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8,
         xlab = "screen Δ cell-cycle index (KO − NTC)", ylab = "bulk Δ cell-cycle index (KO − NTC)", xlim = range(d$screen_delta) + c(-0.02, 0.07))
    abline(h = 0, v = 0, col = GRID); text(d$screen_delta, d$delta_vs_NTC, d$KO, pos = 4, cex = 0.7, offset = 0.4)
    mtext(cl, side = 3, adj = 0, line = 1.9, font = 2, cex = 0.8, col = LINE_COL[[cl]])
    mtext(sprintf("Spearman ρ=%+.2f (p=%.2f) · sign agreement %.0f%%", ct$estimate, ct$p.value, 100 * ag), side = 3, adj = 0, line = 0.7, cex = 0.7, font = 2, col = LINE_COL[[cl]])
  }
  ct <- cor.test(m$screen_delta, m$delta_vs_NTC, method = "spearman", exact = FALSE)
  summ[["ALL (pooled)"]] <<- data.frame(cell_line = "ALL (pooled)", n_KO = nrow(m), spearman_rho = unname(ct$estimate), p = ct$p.value,
                                        sign_agreement = mean(sign(m$screen_delta) == sign(m$delta_vs_NTC)))
  plot(m$screen_delta, m$delta_vs_NTC, pch = 21, bg = LINE_COL[as.character(m$cell_line)], col = "white", cex = 1.4, bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.8,
       xlab = "screen Δ cell-cycle index (KO − NTC)", ylab = ""); abline(h = 0, v = 0, col = GRID)
  legend("topleft", LINES4, pt.bg = LINE_COL[LINES4], pch = 21, col = "white", bty = "n", cex = 0.7)
  mtext("all lines pooled", side = 3, adj = 0, line = 1.9, font = 2, cex = 0.8)
  mtext(sprintf("ρ=%+.2f (p=%.2f)", ct$estimate, ct$p.value), side = 3, adj = 0, line = 0.7, cex = 0.7, font = 2)
  mtext("Does bulk reproduce the screen's proliferation shifts? (shared KOs only)", outer = TRUE, side = 3, adj = 0, font = 2, cex = 1)
})
summ <- do.call(rbind, summ); rownames(summ) <- NULL
write.csv(summ, file.path(OUT, "fig2B_summary.csv"), row.names = FALSE); print(format(summ, digits = 3), row.names = FALSE)

## ---------------------------------------------------------------- Fig 2C (all bulk KOs) and 2D (fig1 layout, screen KOs only)
prolif_boxes <- function(name, w, h, side_by_side) {
  save_fig(name, w, h, function() {
    so <- order_mean(scr_means); so <- c("NTC", setdiff(so, "NTC"))     # CROP-seq screen order (NTC first, then by screen mean) -- shared reference for every panel
    bo <- lapply(THREE, function(cl) {
      avail <- unique(bulk$Gene_Control[bulk$Cell_Line == cl]); shared <- intersect(so, avail); extra <- setdiff(avail, so)   # extra = bulk KOs not in the screen (RUNX1/RUNX2/SNAI2): keep them, append after the screen-ordered ones
      if (!length(extra)) return(shared)
      d <- bulk[bulk$Cell_Line == cl & bulk$Gene_Control %in% extra, ]; c(shared, order_mean(tapply(d$cc, d$Gene_Control, mean)))
    })
    wd <- c(length(so), if (side_by_side) 2.2, sapply(bo, length))
    layout(matrix(seq_along(wd), nrow = 1), widths = wd)
    par(oma = c(0, 0, 5, 0), mar = c(8.5, if (side_by_side) 3.8 else 4.2, 2.2, 0.6), family = "sans")
    scr_violin_panel(so, if (side_by_side) "CROP-seq screen" else "CROP-seq screen (cells; Wilcoxon vs NTC, FDR)")
    if (side_by_side) plot.new()
    for (i in seq_along(THREE))
      bulk_box_panel(THREE[i], bo[[i]], if (side_by_side) paste0(THREE[i], " bulk (RUN1)") else paste0(THREE[i], " bulk (RUN1; each dot = one replicate)"),
                     if (side_by_side) LINE_COL[[THREE[i]]] else INK, ylab = !side_by_side)
    if (side_by_side) {
      mtext("Proliferation index by KO: CROP-seq (left) vs bulk RNA-seq, three cell lines, RUN1 (right)", outer = TRUE, side = 3, adj = 0, line = 3, font = 2, cex = 1)
      mtext("cell-cycle index = S.Score + G2M.Score (same as CROP-seq); bulk restricted to KOs present in the screen; grey = NTC, brown = KO; NTC first, then in the CROP-seq order (sorted by screen mean);",
            outer = TRUE, side = 3, adj = 0, line = 1.7, cex = 0.55, col = INK2)
      mtext("bulk: each dot = one replicate, dashed = NTC mean, y-axis standardised per line; screen: violin + box, stars = Wilcoxon vs NTC (FDR); bracket = cells (screen) or replicates (bulk).",
            outer = TRUE, side = 3, adj = 0, line = 0.9, cex = 0.55, col = INK2)
    } else {
      mtext("Proliferation index by KO: CROP-seq vs bulk RNA-seq (S + G2M cell-cycle marker score, sorted by mean; grey = NTC, brown = KO)", outer = TRUE, side = 3, adj = 0, line = 2.5, font = 2, cex = 1)
      mtext("bulk KOs shared with the screen follow the CROP-seq screen order (NTC first, then its own mean, high to low); KOs not in the screen (RUNX1/RUNX2/SNAI2) appended after, sorted by their own bulk mean",
            outer = TRUE, side = 3, adj = 0, line = 1.3, cex = 0.55, col = INK2)
    }
  })
}
prolif_boxes("fig2C_proliferation_boxplots_CROPseq_vs_bulk_RUN1", 21, 5.2, FALSE)
bulk_screen_kos <- bulk[bulk$Gene_Control %in% c(setdiff(unique(scr$merged_call), "NTC"), "NTC"), ]     # 2D: only KOs present in the screen
bulk <- bulk_screen_kos
prolif_boxes("fig2D_proliferation_side_by_side_CROPseq_vs_3lines_RUN1", 19, 5.6, TRUE)

## ---------------------------------------------------------------- Fig 2E: proliferation combined into ONE boxplot per KO -- replicate-averaged per line,
## z-scored within each source's own KO distribution (so lines/screen are comparable despite different baseline dynamic ranges), all 4 bulk lines + screen
bd$z <- ave(bd$delta_vs_NTC, bd$cell_line, FUN = function(x) (x - mean(x)) / sd(x))
st$z <- (st$delta_vs_NTC - mean(st$delta_vs_NTC)) / sd(st$delta_vs_NTC)
comb <- rbind(data.frame(KO = bd$KO, source = bd$cell_line, z = bd$z, is_screen = FALSE),
              data.frame(KO = st$KO, source = "CROP-seq screen", z = st$z, is_screen = TRUE))
so2 <- setdiff(order_mean(scr_means), "NTC")
extra2 <- setdiff(unique(bd$KO), so2)
extra2_ord <- if (length(extra2)) order_mean(tapply(bd$delta_vs_NTC[bd$KO %in% extra2], bd$KO[bd$KO %in% extra2], mean)) else character(0)
ko_order2 <- c(so2, extra2_ord)

save_fig("fig2E_proliferation_combined_zscored_by_KO", 9.5, 5.2, function() {
  k <- length(ko_order2)
  par(mar = c(5.8, 4.4, 4.4, 0.6), family = "sans")
  dat <- lapply(ko_order2, function(g) comb$z[comb$KO == g])
  yl <- range(unlist(dat)); yl <- yl + c(-1, 1) * 0.10 * diff(yl)
  plot(NA, xlim = c(0.4, k + 0.6), ylim = yl, xaxt = "n", xlab = "", ylab = "z-scored proliferation index (Δ vs NTC, within source)", bty = "n", las = 1, cex.axis = 0.75, cex.lab = 0.7)
  abline(h = 0, col = GRID, lwd = 0.8)
  boxplot(dat, at = seq_len(k), add = TRUE, boxwex = 0.55, outline = FALSE, axes = FALSE, col = adjustcolor(BROWN, 0.3), border = INK2, medcol = INK, medlwd = 1.4, whisklty = 1, staplelty = 1, lwd = 0.8)
  for (i in seq_len(k)) {
    d <- comb[comb$KO == ko_order2[i], ]; scr_i <- d$is_screen
    points(i + runif(sum(!scr_i), -0.16, 0.16), d$z[!scr_i], pch = 21, bg = LINE_COL[d$source[!scr_i]], col = "white", cex = 1.5, lwd = 0.6)
    if (any(scr_i)) points(i, d$z[scr_i], pch = 23, bg = LINE_COL[["CROP-seq"]], col = INK, cex = 1.8, lwd = 1.3)
  }
  axis(1, at = seq_len(k), labels = ko_order2, las = 2, tick = FALSE, cex.axis = 0.72, font = 2, line = -0.4)
  legend("topright", c(LINES4, "CROP-seq screen"), pch = c(rep(21, 4), 23), pt.bg = LINE_COL[c(LINES4, "CROP-seq")], col = c(rep("white", 4), INK), pt.cex = 1.1, bty = "n", cex = 0.6)
  mtext("Proliferation index by KO: all bulk lines + CROP-seq screen, combined into one box per KO", side = 3, adj = 0, line = 3.0, font = 2, cex = 0.9)
  mtext("each round point = one cell line's replicate-averaged Δ(cell-cycle index vs NTC), z-scored within that line's own KO distribution; diamond = CROP-seq screen, z-scored the same way",
        side = 3, adj = 0, line = 1.8, cex = 0.55, col = INK2)
  mtext("ordered by the CROP-seq screen's own mean (high to low); KOs not in the screen (RUNX1/RUNX2/SNAI2) appended after, sorted by their own bulk mean; box = distribution across the up-to-5 sources",
        side = 3, adj = 0, line = 0.9, cex = 0.55, col = INK2)
})

## ---------------------------------------------------------------- Fig 2F: proliferation index, bulk only, one boxplot per KO -- one point per cell line (replicates averaged),
## raw index (NOT subtracted vs NTC), no CROP-seq shown, screen-shared KOs only (RUNX1/RUNX2/SNAI2 excluded), CROP-seq screen order.
## Statistics: paired t-test per KO, KO vs NTC mean_cc, paired by cell line (n = up to 4 lines); BH-adjusted across the 8 KOs.
bd_full <- do.call(rbind, lapply(LINES4, function(cl) {
  d <- bulk_full[bulk_full$Cell_Line == cl, ]
  do.call(rbind, lapply(sort(unique(d$Gene_Control), method = "radix"), function(g) data.frame(cell_line = cl, KO = g, mean_cc = mean(d$cc[d$Gene_Control == g]))))
}))
so2b <- intersect(so2, unique(bulk_full$Gene_Control))    # so2 (CROP-seq screen order, no NTC) restricted to KOs actually present in bulk -- drops MEF2A/MEF2D (screen-only, not in bulk); RUNX1/RUNX2/SNAI2 excluded separately (not in so2)
so3 <- c("NTC", so2b)
data3 <- lapply(so3, function(g) bd_full$mean_cc[bd_full$KO == g]); line3 <- lapply(so3, function(g) bd_full$cell_line[bd_full$KO == g])
k3 <- length(so3)

ntc_v <- setNames(bd_full$mean_cc[bd_full$KO == "NTC"], bd_full$cell_line[bd_full$KO == "NTC"])
ko_stats3 <- do.call(rbind, lapply(so2b, function(g) {
  ko_v <- setNames(bd_full$mean_cc[bd_full$KO == g], bd_full$cell_line[bd_full$KO == g]); cl <- intersect(names(ntc_v), names(ko_v))
  tt <- t.test(ko_v[cl], ntc_v[cl], paired = TRUE)
  data.frame(KO = g, n_lines = length(cl), mean_delta = mean(ko_v[cl] - ntc_v[cl]), t = unname(tt$statistic), p = tt$p.value)
}))
ko_stats3$padj <- p.adjust(ko_stats3$p, "BH")
ko_stats3$stars <- ifelse(ko_stats3$padj < 0.001, "***", ifelse(ko_stats3$padj < 0.01, "**", ifelse(ko_stats3$padj < 0.05, "*", "")))
write.csv(ko_stats3, file.path(OUT, "fig2F_bulk_pooled_stats.csv"), row.names = FALSE)
cat("\nFig 2F statistics: paired t-test, KO mean_cc vs NTC mean_cc, paired by cell line (BH-adjusted across 8 KOs):\n")
print(format(ko_stats3, digits = 3), row.names = FALSE)

save_fig("fig2F_proliferation_boxplots_bulk_all_lines_pooled", 2.4, 2.3, function() {
  par(mar = c(3.0, 2.2, 2.5, 0.15), family = "sans")
  r <- range(unlist(data3)); yl <- r + c(-1, 1) * 0.10 * diff(r); yl[2] <- yl[2] + 0.30 * diff(yl)
  frame_panel(k3, yl, so3)
  boxplot(data3, at = seq_len(k3), add = TRUE, boxwex = 0.5, outline = FALSE, axes = FALSE,
          col = adjustcolor(ifelse(so3 == "NTC", GREY, BROWN), 0.5), border = INK2, medcol = INK, medlwd = 0.8, lwd = 0.5, whisklty = 1, staplelty = 1)
  for (i in seq_len(k3)) points(i + runif(length(data3[[i]]), -0.16, 0.16), data3[[i]], pch = 21, bg = LINE_COL[as.character(line3[[i]])], col = "white", cex = 0.6, lwd = 0.25)
  for (i in seq_len(k3)) if (so3[i] != "NTC") text(i, yl[2] - 0.26 * diff(yl), ko_stats3$stars[ko_stats3$KO == so3[i]], cex = 0.45, font = 2)
  abline(h = median(bd_full$mean_cc[bd_full$KO == "NTC"]), col = MUTED, lty = 2, lwd = 0.5)
  mtext("cc index", side = 2, line = 1.3, cex = 0.32)
  legend("top", LINES4, pt.bg = LINE_COL[LINES4], pch = 21, col = "white", pt.cex = 0.42, bty = "n", cex = 0.24, ncol = 2, x.intersp = 0.4, y.intersp = 0.9, inset = c(0, -0.02))
  mtext("Proliferation by KO (bulk, avg/line)", side = 3, adj = 0, line = 1.35, font = 2, cex = 0.36)
})
