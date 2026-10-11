# Spec 106: the scatter of the 160 KO x Cm heatmap cells of fig9O (spec93), CROP-seq vs bulk NCI-H2052, as a stand-alone figure coloured by KO.
# Same values as fig9C / fig9O (spec81: median NTC - median KO module score, AddModuleScore top-50 genes):
#   A  the 160 cells as in the heatmaps (NTC - KO)      B  the 20 module averages over the 8 KOs (shared response)      C  KO-specific part: each cell minus its module average (the 160-cell scatter of fig9O)
# needs: spec81 (spec81_Cm_KO_NTCminusKO_screen_and_bulk.csv)
# outputs: spec106_Cm_heatmap_cells.csv, fig21_Cm_heatmap_cells_scatter.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
CM <- paste0("Cm", 1:20); x <- read.csv(file.path(OUT, "spec81_Cm_KO_NTCminusKO_screen_and_bulk.csv")); W <- function(s) { d <- x[x$source == s, ]; M <- tapply(d$NTC_minus_KO, list(d$KO, d$Cm), mean); M[ORDER, CM] }
S <- W("screen"); B <- W("NCI-H2052"); Ss <- sweep(S, 2, colMeans(S)); Bs <- sweep(B, 2, colMeans(B)); fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
D <- data.frame(KO = rep(ORDER, 20), Cm = rep(CM, each = 8), screen = as.vector(S), bulk = as.vector(B), screen_specific = as.vector(Ss), bulk_specific = as.vector(Bs)); write.csv(D, file.path(OUT, "spec106_Cm_heatmap_cells.csv"), row.names = FALSE)
ct <- function(a, b) cor.test(a, b, method = "spearman", exact = FALSE); r1 <- ct(D$screen, D$bulk); r2 <- ct(colMeans(S), colMeans(B)); r3 <- ct(D$screen_specific, D$bulk_specific)
cat(sprintf("full 160 cells: rho %.2f p %.3g | 20 module averages: rho %.2f p %.3g | KO-specific 160 cells: rho %.2f p %.3g\n", r1$estimate, r1$p.value, r2$estimate, r2$p.value, r3$estimate, r3$p.value))
save_fig("fig21_Cm_heatmap_cells_scatter", 15, 5.8, function() {
  layout(matrix(1:3, nrow = 1)); par(oma = c(0, 0, 3, 0.4), family = "sans", mar = c(4.6, 4.8, 3.4, 1)); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.85); mtext(b, side = 3, adj = 0, line = 0.4, cex = 0.62, col = INK2) }
  sc <- function(a, b, xl, yl, n_lab = 12) { plot(a, b, pch = 21, bg = adjustcolor(KO_COL[D$KO], 0.85), col = INK2, cex = 1.4, bty = "n", las = 1, xlab = xl, ylab = yl, cex.axis = 0.85); abline(h = 0, v = 0, lty = 2, col = GRID); abline(lm(b ~ a), col = INK2, lty = 3)
    far <- order(-abs(scale(a)) - abs(scale(b)))[seq_len(n_lab)]; text(a[far], b[far], D$Cm[far], pos = 3, cex = 0.6, col = KO_COL[D$KO[far]], xpd = NA) }
  sc(D$screen, D$bulk, "CROP-seq: Cm score, NTC - KO", "bulk H2052: Cm score, NTC - KO"); legend("topleft", ORDER, pch = 21, pt.bg = KO_COL[ORDER], col = INK2, bty = "n", cex = 0.7, pt.cex = 1.2)
  ttl("A  the 160 heatmap cells (8 KOs x 20 Cm)", sprintf("Spearman %.2f, p = %s", r1$estimate, fp(r1$p.value)))
  a <- colMeans(S); b <- colMeans(B); plot(a, b, pch = 21, bg = "grey45", col = INK2, cex = 1.6, bty = "n", las = 1, xlab = "CROP-seq: module average over the 8 KOs (NTC - KO)", ylab = "bulk H2052: module average over the 8 KOs (NTC - KO)", cex.axis = 0.85); abline(h = 0, v = 0, lty = 2, col = GRID); abline(lm(b ~ a), col = INK2, lty = 3); text(a, b, CM, pos = 3, cex = 0.62, col = INK2, xpd = NA)
  ttl("B  shared response: 20 module averages", sprintf("Spearman %.2f, p = %s", r2$estimate, fp(r2$p.value)))
  sc(D$screen_specific, D$bulk_specific, "CROP-seq: KO-specific module effect", "bulk H2052: KO-specific module effect"); ttl("C  KO-specific part: cell minus its module average (fig9O)", sprintf("Spearman %.2f, p = %s", r3$estimate, fp(r3$p.value)))
  mtext("NCI-H2052: the KO x Cm heatmap cells, CROP-seq vs bulk RNA-seq", outer = TRUE, side = 3, adj = 0, line = 1.5, font = 2, cex = 1)
  mtext("median NTC - median KO module score (AddModuleScore, top-50 genes), as in fig9C / fig9O; > 0 = module lower in the KO; colour = KO", outer = TRUE, side = 3, adj = 0, line = 0.3, cex = 0.65, col = INK2)
})
