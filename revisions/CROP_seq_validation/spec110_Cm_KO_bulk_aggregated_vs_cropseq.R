# Spec 110: bulk KO effect with replicates averaged and the 4 cell lines aggregated, against the same quantity in the CROP-seq screen.
# For every Cm module x KO (20 x 8 = 160 pairs) the effect is defined identically in both data sets:
#       delta = mean(the KO) - mean(the rest of that data set: NTC + the other 7 KOs)
# Bulk:  scores z-scored per module WITHIN each line; per line the replicates of each group are averaged first (so a line with 3 NTC
#        replicates does not outweigh one with 1), delta is computed per line, then the 4 line deltas are averaged with EQUAL weight.
#        Averaging replicates before the line delta is the point of this script -- spec109 pooled samples, which let replicate count
#        drive the result; here each line contributes exactly one number per pair.
# CROP:  cell-level scores z-scored per module across the cells of the 9 groups, group means, same delta.
# Because each line gets one delta, the sign consistency across lines (how many of the 4 agree with the CROP-seq sign) is also reported --
# that, not the pooled magnitude, is what separates a reproducible pair from one carried by a single line.
# needs: repro_common.R, profiles.R, spec72_profiles_S70.rds, spec75_screen_scores_by_cell.rds
# outputs: spec110_Cm_KO_bulk_vs_cropseq.csv, fig22_Cm_KO_bulk_aggregated_vs_cropseq.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
CM <- paste0("Cm", 1:20); GR <- c("NTC", KOS)
HILITE <- list(c("Cm20", "TEAD4"), c("Cm8", "TWIST1"), c("Cm9", "TWIST1"), c("Cm10", "TCF3"), c("Cm10", "TEAD2"), c("Cm1", "TEAD2"))

pr <- readRDS(file.path(OUT, "spec72_profiles_S70.rds")); fab <- load_full_audited_bulk(); meta <- fab$meta
ids <- pr$samples; line <- meta[ids, "Cell_Line"]; grp <- meta[ids, "Gene_Control"]
S <- pr$cm_scores$tmm[CM, ids, drop = FALSE]
Z <- S; for (cl in LINES) { j <- which(line == cl); Z[, j] <- t(scale(t(S[, j, drop = FALSE]))) }; Z[!is.finite(Z)] <- NA

## bulk: per line -> replicate means per group -> delta(KO vs rest) -> average the 4 lines
dl <- array(NA_real_, c(length(CM), length(KOS), length(LINES)), dimnames = list(CM, KOS, LINES))
for (cl in LINES) { j <- line == cl
  gm <- sapply(GR, function(g) rowMeans(Z[, j & grp == g, drop = FALSE], na.rm = TRUE))   # replicate mean per group
  for (k in KOS) dl[, k, cl] <- gm[, k] - rowMeans(gm[, setdiff(GR, k), drop = FALSE], na.rm = TRUE) }
B <- apply(dl, c(1, 2), mean, na.rm = TRUE)                                               # equal weight per line

## CROP-seq: same delta from cell-level scores
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); gs <- as.character(sc$group); keep <- gs %in% GR
zS <- scale(sc$score[keep, CM]); gk <- gs[keep]
mS <- sapply(GR, function(g) colMeans(zS[gk == g, , drop = FALSE]))
C <- sapply(KOS, function(k) mS[, k] - rowMeans(mS[, setdiff(GR, k), drop = FALSE]))
dimnames(C) <- list(CM, KOS)

nsign <- sapply(KOS, function(k) sapply(CM, function(m) sum(sign(dl[m, k, ]) == sign(C[m, k]), na.rm = TRUE)))
D <- data.frame(Cm = rep(CM, length(KOS)), KO = rep(KOS, each = length(CM)),
                cropseq = as.vector(C), bulk_agg = as.vector(B), n_lines_agree = as.vector(nsign),
                as.data.frame(matrix(aperm(dl, c(1, 2, 3)), ncol = length(LINES), dimnames = list(NULL, sub("NCI-|-211H", "", LINES)))))
write.csv(D, file.path(OUT, "spec110_Cm_KO_bulk_vs_cropseq.csv"), row.names = FALSE)
ct <- function(x, y) { r <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)); c(rho = unname(r$estimate), p = r$p.value) }
all_r <- ct(D$cropseq, D$bulk_agg); perKO <- t(sapply(KOS, function(k) ct(C[, k], B[, k])))
options(width = 200)
cat(sprintf("\n=== bulk (replicates averaged, 4 lines aggregated) vs CROP-seq, %d Cm x KO pairs ===\nSpearman rho = %.3f, p = %.3g\n", nrow(D), all_r[1], all_r[2]))
cat("\nper KO (over the 20 modules):\n"); print(round(perKO, 3))
cat("\nhow often do all 4 bulk lines agree in sign with CROP-seq?\n"); print(table(D$n_lines_agree))
cat("\n=== the 6 requested pairs ===\n")
h <- do.call(rbind, lapply(HILITE, function(x) D[D$Cm == x[1] & D$KO == x[2], ]))
print(format(h, digits = 3), row.names = FALSE)

fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
hid <- paste(D$Cm, D$KO) %in% sapply(HILITE, paste, collapse = " ")
save_fig("fig22_Cm_KO_bulk_aggregated_vs_cropseq", 13.5, 7.2, function() {
  layout(cbind(matrix(1, 2, 2), matrix(2:9, nrow = 2, byrow = TRUE)), widths = c(1, 1, rep(0.8, 4)))
  par(oma = c(0, 0, 3.4, 0.4), family = "sans", mar = c(4.4, 4.6, 3.0, 1))
  plot(D$cropseq, D$bulk_agg, pch = 21, bg = adjustcolor(KO_COL[D$KO], 0.85), col = INK2, cex = 1.3, bty = "n", las = 1,
       xlab = "CROP-seq: KO minus the rest (z)", ylab = "bulk: KO minus the rest, replicates averaged, 4 lines aggregated (z)", cex.axis = 0.85)
  abline(h = 0, v = 0, lty = 2, col = GRID); abline(lm(bulk_agg ~ cropseq, D), col = INK2, lty = 3)
  points(D$cropseq[hid], D$bulk_agg[hid], pch = 21, bg = adjustcolor(KO_COL[D$KO[hid]], 0.95), col = "black", cex = 2.0, lwd = 1.4)
  text(D$cropseq[hid], D$bulk_agg[hid], paste(D$Cm[hid], D$KO[hid]), pos = 3, cex = 0.62, font = 2, col = INK, xpd = NA)
  legend("topleft", KOS, pch = 21, pt.bg = KO_COL[KOS], col = INK2, bty = "n", cex = 0.72, pt.cex = 1.1)
  mtext("all 20 Cm x 8 KOs", side = 3, adj = 0, line = 1.5, font = 2, cex = 0.85)
  mtext(sprintf("Spearman %.2f, p = %s (%d pairs); ringed = the 6 requested pairs", all_r[1], fp(all_r[2]), nrow(D)), side = 3, adj = 0, line = 0.4, cex = 0.62, col = INK2)
  for (k in KOS[order(-perKO[, 1])]) { par(mar = c(3.6, 3.4, 3.0, 0.6))
    plot(C[, k], B[, k], pch = 21, bg = adjustcolor(KO_COL[[k]], 0.85), col = INK2, cex = 1.1, bty = "n", las = 1, xlab = "", ylab = "", cex.axis = 0.7)
    abline(h = 0, v = 0, lty = 2, col = GRID); abline(lm(B[, k] ~ C[, k]), col = INK2, lty = 3)
    text(C[, k], B[, k], sub("Cm", "", CM), pos = 3, cex = 0.5, col = INK2, xpd = NA)
    mtext(k, side = 3, adj = 0, line = 1.3, font = 2, cex = 0.78, col = KO_COL[[k]])
    mtext(sprintf("rho %.2f, p %s", perKO[k, 1], fp(perKO[k, 2])), side = 3, adj = 0, line = 0.4, cex = 0.58, col = INK2, font = ifelse(perKO[k, 2] < 0.05, 2, 1))
    mtext("CROP-seq", side = 1, line = 2.2, cex = 0.52); mtext("bulk (4 lines)", side = 2, line = 2.2, cex = 0.52) }
  mtext("KO effect on each Cm module: bulk with replicates averaged and cell lines aggregated, vs the CROP-seq screen", outer = TRUE, side = 3, adj = 0, line = 1.7, font = 2, cex = 0.95)
  mtext("effect = mean(KO) - mean(NTC + the other 7 KOs), computed the same way in both; bulk: z per module within line, replicate means, then the 4 line deltas averaged with equal weight", outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.6, col = INK2)
})
