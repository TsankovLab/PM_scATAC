# Spec 116: the three summary figures of the reproducibility investigation (spec114 reliability, spec115 pipeline grid, spec73-79 benchmarks).
# fig28  ACROSS BULK LINES -- same KO vs other KOs, in the configuration that separates the most TFs.
#        Gene level, not module level: spec73 (12,672 configurations) and spec115 (1,280) agree that the cross-line KO signal lives in the genes that respond to the KOs
#        (the 500 genes with the largest KO signal in the two lines NOT being compared, so the selection never sees the pair it is scored on), not in the 20 patient-derived Cm.
#        Configuration = the top of spec73 by marginal effects: S70, TMM logCPM, limma moderated t of (mean of the other KOs) - KO, centred over the KOs, Spearman.
#        A 32 x 32 KO-by-line correlation map ordered by KO | B per TF same-KO vs other-KO with KO-label permutation p | C robustness of each TF over the 256 gene-level
#        configurations of spec73 | D the same TFs scored on the 20 Cm (best cross-line Cm configuration of spec115) for contrast.
# fig29  CROP-SEQ vs BULK on the Cm, lower-in-KO focus (sign convention here: KO - reference, NEGATIVE = lower in the KO).
#        A screen vs 4-line bulk mean, KO - NTC | B the agreement split into the response SHARED by all KOs (20 module means) and the KO-SPECIFIC remainder, per line |
#        C the pairs the screen puts lower in the KO under BOTH references, ranked by how often the bulk lines agree across ALL processing choices of spec115
#        (640 S70 configurations) -- a pair that only reproduces under one pipeline is not counted as reproducible.
# fig30  WHERE THE REPRODUCIBILITY IS LOST -- the ladder from within-screen to across data sets, the batch nesting, the knockdown and the variance budget.
# needs: spec72 (S70), spec73 grids, spec74 tables, spec79 i5 blocks, spec114, spec115
# outputs: spec116_across_lines_per_KO.csv, spec116_lower_in_KO_robust.csv, spec116_ladder.csv, spec116_i5_block_similarity.csv, fig28/fig29/fig30 {png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R")); options(width = 230)
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000")); permsL <- make_perms(N_PERM); GR <- c("NTC", KOS); CM <- paste0("Cm", 1:20); short <- sub("NCI-|-211H", "", LINES)
fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2)); sp <- function(a, b) suppressWarnings(cor(as.vector(a), as.vector(b), method = "spearman", use = "pairwise.complete.obs"))
pr <- readRDS(file.path(OUT, "spec72_profiles_S70.rds")); fab <- load_full_audited_bulk(); meta <- fab$meta; cs <- fab$cs; q <- load_audit(); ids <- pr$samples; line <- meta[ids, "Cell_Line"]; grp <- meta[ids, "Gene_Control"]

## ================================================================ fig28: across lines, gene level
A <- array(NA_real_, c(length(pr$genes), 8, 4), dimnames = list(pr$genes, KOS, LINES)); sig <- matrix(NA_real_, length(pr$genes), 4, dimnames = list(pr$genes, LINES)); GM <- list()
for (cl in LINES) { s <- ids[line == cl]; y <- normalise(cs[pr$genes, s], "tmm")$y; e <- ko_effects(y, meta[s, "Gene_Control"]); A[, , cl] <- e$t_rest; sig[, cl] <- matrixStats::rowMaxs(abs(e$t), na.rm = TRUE)
  GM[[cl]] <- sapply(KOS, function(k) rowMeans(y[, meta[s, "Gene_Control"] == k, drop = FALSE])) }
CL <- cl_from_array(A, "signal", sig, n_sig = 500, centre = TRUE, method = "spearman"); ag <- agreement(CL, permsL); pk <- ag$per_KO[order(ag$per_KO$perm_p, -ag$per_KO$mean_intra), ]
k73 <- read.csv(file.path(OUT, "spec73_grid_per_KO.csv")); k73 <- k73[k73$set == "S70" & k73$space %in% c("signal500", "signal1000") & k73$prep == "centred", ]
pk$robust_frac <- tapply(k73$perm_p < 0.05, k73$KO, mean)[pk$KO]; pk$robust_median_z <- tapply(k73$perm_z, k73$KO, median)[pk$KO]
R115 <- readRDS(file.path(OUT, "spec115_profiles.rds")); G115 <- read.csv(file.path(OUT, "spec115_grid.csv")); bestCm <- G115[G115$set == "S70", ]; bestCm <- bestCm[which.max(bestCm$lines_z), ]
keyCm <- paste(bestCm$set, bestCm$norm, bestCm$score, bestCm$cov, bestCm$zline, bestCm$ref, bestCm$stat); agCm <- agreement(cl_from_array(R115$profiles[[keyCm]], centre = TRUE, method = "spearman"), permsL)
pk$Cm_perm_z <- agCm$per_KO$perm_z[match(pk$KO, agCm$per_KO$KO)]; pk$Cm_perm_p <- agCm$per_KO$perm_p[match(pk$KO, agCm$per_KO$KO)]
write.csv(pk, file.path(OUT, "spec116_across_lines_per_KO.csv"), row.names = FALSE)
cat("=== fig28: across the 4 bulk lines, gene level (S70 | tmm | t_rest | top-500 KO-signal genes | centred | Spearman) ===\n")
cat(sprintf("global: mean same-KO %.3f, other-KO %.3f, permutation p = %s (z %.2f); TFs with p < 0.05: %d of 8\n", ag$global[["diag_mean"]], ag$global[["off_mean"]], fp(ag$global[["perm_p"]]), ag$global[["perm_z"]], sum(pk$perm_p < 0.05)))
print(format(pk[, c("KO", "mean_intra", "median_cross", "n_pos_pairs", "perm_p", "perm_fdr", "perm_z", "robust_frac", "Cm_perm_z", "Cm_perm_p")], digits = 3), row.names = FALSE)
cat("best Cm-level configuration for cross-line agreement:", keyCm, sprintf("-> global z %.2f, p %s\n", agCm$global[["perm_z"]], fp(agCm$global[["perm_p"]])))

BIG <- matrix(NA_real_, 32, 32); lab <- expand.grid(line = LINES, KO = pk$KO, stringsAsFactors = FALSE); ix <- function(k, l) match(paste(k, l), paste(lab$KO, lab$line))
for (i in seq_along(PAIRS)) { p <- PAIRS[[i]]; C <- CL[[i]]; for (a in KOS) for (b in KOS) { BIG[ix(a, p[1]), ix(b, p[2])] <- C[a, b]; BIG[ix(b, p[2]), ix(a, p[1])] <- C[a, b] } }
save_fig("fig28_across_lines_best_configuration", 15.5, 8.6, function() {
  layout(matrix(c(1, 2, 1, 3, 1, 4), nrow = 2), widths = c(1.25, 1, 0.8), heights = c(1.5, 1)); par(family = "sans", oma = c(0, 0, 3.4, 0.4))
  par(mar = c(1.5, 7.5, 4.2, 1.5)); vm <- 0.35; n <- 32
  image(1:n, 1:n, t(pmax(pmin(BIG, vm), -vm))[, n:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = ""); na <- which(is.na(BIG), arr.ind = TRUE); rect(na[, 2] - 0.5, n - na[, 1] + 0.5, na[, 2] + 0.5, n - na[, 1] + 1.5, col = "#e9e8e4", border = NA)
  for (b in seq(4.5, 28.5, 4)) abline(h = n - b + 1, v = b, col = "white", lwd = 2); for (b in 0:7) rect(b * 4 + 0.5, n - b * 4 - 3.5, b * 4 + 4.5, n - b * 4 + 0.5, border = "black", lwd = 1.6)
  text(seq(2.5, 30.5, 4), n + 1.4, pk$KO, xpd = NA, font = 2, cex = 0.8, col = KO_COL[pk$KO]); text(-0.6, n - seq(2.5, 30.5, 4) + 1, pk$KO, xpd = NA, font = 2, cex = 0.8, col = KO_COL[pk$KO], adj = 1)
  text(1:n, n + 0.75, substr(sub("NCI-|-211H", "", lab$line), 1, 5), srt = 90, adj = 0, xpd = NA, cex = 0.36, col = INK2)
  mtext("A  KO-by-line Spearman of the KO profiles (black boxes = same KO, different lines; grey = same line, not compared)", side = 3, adj = 0, line = 2.9, font = 2, cex = 0.72)
  par(mar = c(4.2, 6.0, 3.4, 6.4)); o <- pk$KO; n8 <- 8; xr <- range(c(unlist(ag$cross), ag$intra), na.rm = TRUE)
  plot(NA, xlim = xr, ylim = c(n8 + 0.6, 0.4), yaxt = "n", bty = "n", las = 1, ylab = "", xlab = "Spearman between two bulk cell lines", cex.axis = 0.8); abline(v = 0, lty = 2, col = INK)
  for (i in seq_len(n8)) { j <- match(o[i], ag$per_KO$KO); v <- ag$intra[o[i], ]
    boxplot(ag$cross[[j]], at = i + 0.25, add = TRUE, horizontal = TRUE, boxwex = 0.5, col = "#d5d4cf", border = INK2, outline = FALSE, axes = FALSE)
    boxplot(v, at = i - 0.25, add = TRUE, horizontal = TRUE, boxwex = 0.5, col = adjustcolor(KO_COL[[o[i]]], 0.75), border = INK2, outline = FALSE, axes = FALSE); points(v, i - 0.25 + seq(-0.15, 0.15, length.out = length(v)), pch = 16, cex = 0.7, col = INK2)
    text(xr[2] + 0.03 * diff(xr), i, sprintf("p = %s%s", fp(pk$perm_p[i]), ifelse(pk$perm_fdr[i] < 0.01, " **", ifelse(pk$perm_fdr[i] < 0.05, " *", ""))), adj = 0, xpd = NA, cex = 0.74, font = ifelse(pk$perm_p[i] < 0.05, 2, 1)) }
  text(xr[1] - 0.04 * diff(xr), seq_len(n8), o, adj = 1, xpd = NA, font = 2, cex = 0.85, col = KO_COL[o])
  mtext(sprintf("B  same KO (colour) vs other KOs (grey): %d of 8 TFs p < 0.05; all p = %s", sum(pk$perm_p < 0.05), fp(ag$global[["perm_p"]])), side = 3, adj = 0, line = 1.6, font = 2, cex = 0.72, at = xr[1] - 0.12 * diff(xr))
  mtext(sprintf("KO-label permutation within lines (%s); * / ** BH FDR < 0.05 / 0.01", format(N_PERM, big.mark = ",")), side = 3, adj = 0, line = 0.5, cex = 0.55, col = INK2, at = xr[1] - 0.12 * diff(xr))
  par(mar = c(4.2, 5.4, 3.4, 1)); bp <- barplot(rev(pk$robust_frac), horiz = TRUE, col = adjustcolor(KO_COL[rev(o)], 0.85), border = NA, names.arg = rev(o), las = 1, cex.names = 0.75, xlim = c(0, 1), xlab = "fraction of 256 gene-level configurations with p < 0.05")
  abline(v = 0.05, lty = 2, col = "red"); mtext("C  robust to the analysis choices?", side = 3, adj = 0, line = 1.6, font = 2, cex = 0.72); mtext("spec73: normalisation x covariates x statistic x n genes x method; red = chance", side = 3, adj = 0, line = 0.5, cex = 0.52, col = INK2)
  m <- rbind(genes = pk$perm_z, Cm = pk$Cm_perm_z); bp <- barplot(m[, 8:1], beside = TRUE, horiz = TRUE, col = c("#4a3aa7", "#d5d4cf"), border = NA, names.arg = rev(o), las = 1, cex.names = 0.75, xlab = "permutation z of the same-KO correlation")
  abline(v = c(0, 1.645), lty = c(1, 2), col = c(INK, "red")); legend("bottomright", c("KO-signal genes", "20 Cm modules (best Cm configuration)"), fill = c("#4a3aa7", "#d5d4cf"), border = NA, bty = "n", cex = 0.62)
  mtext("D  genes vs Cm modules as the feature space", side = 3, adj = 0, line = 1.6, font = 2, cex = 0.72); mtext("red = one-sided p 0.05", side = 3, adj = 0, line = 0.5, cex = 0.52, col = INK2)
  mtext("Reproducibility of each TF KO across the 4 bulk cell lines, in the configuration that separates the most TFs", outer = TRUE, side = 3, adj = 0.01, line = 1.7, font = 2, cex = 1)
  mtext("S70, TMM logCPM, limma moderated t of (mean of the other KOs) - KO, KO-centred; 500 genes with the largest KO signal in the two lines not being compared; Spearman", outer = TRUE, side = 3, adj = 0.01, line = 0.5, cex = 0.62, col = INK2)
})

## ================================================================ fig29: screen vs bulk on the Cm, lower in the KO
kN <- "S70 tmm mz100 none FALSE ntc diff"; kR <- "S70 tmm mz100 none FALSE rest diff"; AN <- R115$profiles[[kN]]; AR <- R115$profiles[[kR]]; SN <- R115$screen[["FALSE ntc"]]; SR <- R115$screen[["FALSE rest"]]
dec <- do.call(rbind, lapply(c(LINES, "mean of 4 lines"), function(cl) { B <- if (cl %in% LINES) AN[, , cl] else apply(AN, c(1, 2), mean, na.rm = TRUE); bs <- rowMeans(B, na.rm = TRUE); ss <- rowMeans(SN)
  data.frame(line = cl, total = sp(SN, B), shared = sp(ss, bs), specific = sp(SN - ss, B - bs)) }))
cat("\n=== fig29: screen vs bulk, KO - NTC, split into the response shared by all KOs and the KO-specific remainder (mz100, tmm, no covariate) ===\n"); print(format(dec, digits = 3), row.names = FALSE)
key <- do.call(rbind, strsplit(names(R115$profiles), " ")); colnames(key) <- c("set", "norm", "score", "cov", "zline", "ref", "stat")
rob <- function(rf, f) { i <- which(key[, "set"] == "S70" & key[, "ref"] == rf); Reduce(`+`, lapply(R115$profiles[i], function(X) f(apply(X < 0, c(1, 2), sum, na.rm = TRUE)) * 1)) / length(i) }
LW <- data.frame(Cm = rep(CM, 8), KO = rep(KOS, each = 20), screen_ntc = as.vector(SN), screen_rest = as.vector(SR), all4_ntc = as.vector(rob("ntc", function(x) x == 4)), all4_rest = as.vector(rob("rest", function(x) x == 4)),
                 ge3_ntc = as.vector(rob("ntc", function(x) x >= 3)), ge3_rest = as.vector(rob("rest", function(x) x >= 3)))
for (l in 1:4) { LW[[paste0("ntc_", short[l])]] <- as.vector(AN[, , l]); LW[[paste0("rest_", short[l])]] <- as.vector(AR[, , l]) }
LW$robust <- pmin(LW$ge3_ntc, LW$ge3_rest); LW$robust_all4 <- pmin(LW$all4_ntc, LW$all4_rest)
LD <- LW[LW$screen_ntc < 0 & LW$screen_rest < 0, ]; LD <- LD[order(-LD$robust, -LD$robust_all4), ]; write.csv(LD, file.path(OUT, "spec116_lower_in_KO_robust.csv"), row.names = FALSE)
cat(sprintf("\n=== pairs lower in the KO in the screen under both references: %d. Ranked by the fraction of the 320 + 320 S70 configurations in which >= 3 of 4 bulk lines are also lower ===\n", nrow(LD)))
print(format(head(LD[, c("Cm", "KO", "screen_ntc", "screen_rest", "ge3_ntc", "ge3_rest", "all4_ntc", "all4_rest", paste0("rest_", short))], 14), digits = 2), row.names = FALSE)
cat("chance: >= 3 of 4 lines lower = 5/16 = 0.31; all 4 = 1/16 = 0.06\n")
save_fig("fig29_cropseq_vs_bulk_lower_in_KO", 15.5, 6.6, function() {
  layout(matrix(1:3, nrow = 1), widths = c(1, 0.8, 1.75)); par(family = "sans", oma = c(0, 0.6, 3.4, 0.4), mar = c(4.6, 4.8, 3.4, 1))
  Bm <- apply(AN, c(1, 2), mean, na.rm = TRUE); plot(SN, Bm, pch = 21, bg = adjustcolor(KO_COL[rep(KOS, each = 20)], 0.85), col = INK2, cex = 1.2, bty = "n", las = 1, xlab = "CROP-seq: KO - NTC (z)", ylab = "bulk: KO - NTC, mean of 4 lines (z)", cex.axis = 0.8)
  abline(h = 0, v = 0, lty = 2, col = GRID); abline(lm(as.vector(Bm) ~ as.vector(SN)), lty = 3, col = INK2); rect(par("usr")[1], par("usr")[3], 0, 0, col = adjustcolor("#2a78d6", 0.06), border = NA)
  text(par("usr")[1], par("usr")[3], " lower in the KO in both", adj = c(0, -0.6), cex = 0.62, col = "#2a78d6", font = 2); legend("topleft", KOS, pch = 21, pt.bg = KO_COL[KOS], col = INK2, bty = "n", cex = 0.62, pt.cex = 1)
  mtext("A  all 160 KO x Cm effects, NTC-referenced", side = 3, adj = 0, line = 1.5, font = 2, cex = 0.76); mtext(sprintf("Spearman %.2f", dec$total[5]), side = 3, adj = 0, line = 0.4, cex = 0.6, col = INK2)
  par(mar = c(6.4, 4.6, 3.4, 1)); m <- t(as.matrix(dec[, c("total", "shared", "specific")])); colnames(m) <- c(short, "mean")
  bp <- barplot(m, beside = TRUE, col = c("#898781", "#1baf7a", "#eb6834"), border = NA, las = 2, ylab = "Spearman, screen vs bulk", ylim = range(c(m, -0.3, 0.7)), cex.names = 0.8); abline(h = 0, col = INK, lwd = 0.6)
  legend("topright", c("all 160 effects", "shared by all KOs (20 module means)", "KO-specific remainder"), fill = c("#898781", "#1baf7a", "#eb6834"), border = NA, bty = "n", cex = 0.62)
  mtext("B  agreement = shared response", side = 3, adj = 0, line = 1.5, font = 2, cex = 0.76); mtext("the KO-specific part does not transfer", side = 3, adj = 0, line = 0.4, cex = 0.6, col = INK2)
  n <- min(14, nrow(LD)); d <- LD[seq_len(n), ]; cols <- c("screen_rest", paste0("rest_", short)); M <- as.matrix(d[, cols]); par(mar = c(4.6, 8.5, 4.6, 17)); vm <- 2
  image(seq_along(cols), seq_len(n), t(pmax(pmin(M, vm), -vm))[, n:1], col = DIV, zlim = c(-vm, vm), axes = FALSE, xlab = "", ylab = ""); abline(v = 1.5, col = "white", lwd = 4)
  for (i in seq_len(n)) for (j in seq_along(cols)) text(j, n - i + 1, sprintf("%+.2f", M[i, j]), cex = 0.6, col = INK)
  text(seq_along(cols), n + 0.8, c("CROP-seq", short), xpd = NA, font = 2, cex = 0.72); text(0.4, n:1, paste(d$Cm, d$KO), adj = 1, xpd = NA, font = 2, cex = 0.74, col = KO_COL[d$KO])
  xb <- length(cols) + 0.7; w <- 1.6; text(xb + w / 2, n + 0.8, ">= 3 of 4 lower", xpd = NA, font = 2, cex = 0.62)
  for (i in seq_len(n)) { y <- n - i + 1; rect(xb, y - 0.32, xb + w, y + 0.32, col = "#efefef", border = NA, xpd = NA); rect(xb, y, xb + w * d$ge3_rest[i], y + 0.32, col = "#eb6834", border = NA, xpd = NA); rect(xb, y - 0.32, xb + w * d$ge3_ntc[i], y, col = "#898781", border = NA, xpd = NA)
    text(xb + w + 0.1, y, sprintf("%.0f / %.0f%%", 100 * d$ge3_rest[i], 100 * d$ge3_ntc[i]), adj = 0, cex = 0.55, xpd = NA, col = INK2) }
  segments(xb + w * 5 / 16, 0.5, xb + w * 5 / 16, n + 0.5, lty = 2, col = "red", xpd = NA); legend(xb, 0.1, c("KO - rest", "KO - NTC"), fill = c("#eb6834", "#898781"), border = NA, bty = "n", cex = 0.55, xpd = NA, horiz = TRUE)
  mtext("C  lower in the KO in the screen: which hold in bulk whatever the processing?", side = 3, adj = 0, line = 3.2, font = 2, cex = 0.76, at = 0.5)
  mtext("cells: KO - rest (z; blue = lower in the KO). bars: % of 320 configurations per reference with >= 3 of 4 lines lower; red = chance (31%)", side = 3, adj = 0, line = 2.2, cex = 0.54, col = INK2, at = 0.5)
  mtext("CROP-seq vs bulk on the Cm modules: what reproduces, and the candidates that go down on knockout", outer = TRUE, side = 3, adj = 0.01, line = 1.7, font = 2, cex = 1)
  mtext("module score = mean within-line gene z over the top-100 module genes (the scoring with the best screen agreement in spec115), TMM, no covariate; NEGATIVE = lower in the KO", outer = TRUE, side = 3, adj = 0.01, line = 0.5, cex = 0.62, col = INK2)
})

## ================================================================ fig30: where the reproducibility is lost
scr <- read.csv(file.path(OUT, "spec114_screen_reliability.csv")); bl <- read.csv(file.path(OUT, "spec114_bulk_reliability.csv")); cei <- read.csv(file.path(OUT, "spec114_ceiling.csv"))
S <- pr$cm_scores$tmm[CM, ids]; Z <- S; for (cl in LINES) { j <- which(line == cl); Z[, j] <- t(scale(t(S[, j]))) }
EB <- lapply(setNames(LINES, LINES), function(cl) { M <- sapply(GR, function(g) rowMeans(Z[, line == cl & grp == g, drop = FALSE])); sapply(KOS, function(k) M[, k] - rowMeans(M[, setdiff(GR, k)])) })
xl <- sapply(PAIRS, function(p) sp(EB[[p[1]]], EB[[p[2]]])); obs <- cei$observed[grepl("specific", cei$kind)]; names(obs) <- cei$line[grepl("specific", cei$kind)]
lad <- data.frame(step = c("screen: half of the cells vs the other half", "bulk: replicate vs replicate, same line", "bulk: one line vs another line", "screen vs bulk, SAME line (H2052)", "screen vs bulk, other lines"),
  value = c(scr$half_rho[scr$what == "specific"], mean(bl$half_specific), mean(xl), obs[["NCI-H2052"]], mean(obs[LINES[-1]])),
  lo = c(NA, min(bl$half_specific), min(xl), NA, min(obs[LINES[-1]])), hi = c(NA, max(bl$half_specific), max(xl), NA, max(obs[LINES[-1]])))
write.csv(lad, file.path(OUT, "spec116_ladder.csv"), row.names = FALSE); cat("\n=== fig30: the ladder -- Spearman of the KO-specific KO x Cm effects (160 values) ===\n"); print(format(lad, digits = 3), row.names = FALSE)
i5 <- read.csv(file.path(OUT, "spec79_i5_blocks.csv"), stringsAsFactors = FALSE); bs <- list()
for (cl in LINES) { g <- GM[[cl]]; v <- order(matrixStats::rowVars(g), decreasing = TRUE)[1:2000]; C <- cor(g[v, ] - rowMeans(g[v, ]), method = "spearman"); b <- setNames(i5$i5[i5$line == cl], i5$group[i5$line == cl])[KOS]
  for (a in 1:7) for (c2 in (a + 1):8) if (!grepl("/", b[a]) && !grepl("/", b[c2])) bs[[length(bs) + 1]] <- data.frame(line = cl, KO1 = KOS[a], KO2 = KOS[c2], same_block = b[a] == b[c2], rho = C[a, c2]) }
bs <- do.call(rbind, bs); write.csv(bs, file.path(OUT, "spec116_i5_block_similarity.csv"), row.names = FALSE)
set.seed(1); bobs <- mean(bs$rho[bs$same_block]) - mean(bs$rho[!bs$same_block]); bnul <- replicate(5000, { d <- 0; x <- bs; for (cl in LINES) { i <- which(x$line == cl); x$same_block[i] <- sample(x$same_block[i]) }; mean(x$rho[x$same_block]) - mean(x$rho[!x$same_block]) }); bp_ <- mean(bnul >= bobs)
cat(sprintf("\nlibrary (i5) block nested in KO: within-line similarity of two DIFFERENT KOs, same block %.3f (n %d) vs different block %.3f (n %d); difference %.3f, permutation p = %s\n", mean(bs$rho[bs$same_block]), sum(bs$same_block), mean(bs$rho[!bs$same_block]), sum(!bs$same_block), bobs, fp(max(bp_, 1 / 5001))))
print(round(sapply(LINES, function(cl) c(same = mean(bs$rho[bs$line == cl & bs$same_block]), different = mean(bs$rho[bs$line == cl & !bs$same_block]))), 3))
kd <- read.csv(file.path(OUT, "spec74_knockdown_vs_agreement.csv")); KD <- tapply(kd$own_target_log2FC, list(kd$KO, kd$line), mean)[KOS, LINES]; vf <- read.csv(file.path(OUT, "spec74_variance_fractions.csv"))
cat("\nown-target knock-down, log2 NTC/KO (positive = knocked down); cells < 0.3 (< 19% reduction):", sum(KD < 0.3), "of 32\n"); print(round(KD, 2))
save_fig("fig30_where_reproducibility_is_lost", 15.5, 9, function() {
  layout(matrix(c(1, 1, 2, 3, 4, 5), nrow = 2), widths = c(1.2, 1, 1)); par(family = "sans", oma = c(0, 0, 3.4, 0.4))
  par(mar = c(4.6, 15, 3.6, 2)); n <- nrow(lad); cc <- c("#4a3aa7", "#1baf7a", "#eda100", "#eb6834", "#b3400f"); plot(NA, xlim = c(-0.3, 1), ylim = c(n + 0.5, 0.5), yaxt = "n", bty = "n", xlab = "Spearman of the KO-specific KO x Cm effects (160 values)", ylab = "", las = 1)
  abline(v = 0, col = INK, lwd = 0.7); for (i in seq_len(n)) { rect(0, i - 0.3, lad$value[i], i + 0.3, col = cc[i], border = NA); if (is.finite(lad$lo[i])) segments(lad$lo[i], i, lad$hi[i], i, lwd = 2, col = INK2)
    text(max(lad$value[i], lad$hi[i], 0, na.rm = TRUE) + 0.03, i, sprintf("%.2f", lad$value[i]), adj = 0, cex = 0.9, font = 2); text(-0.34, i, lad$step[i], adj = 1, xpd = NA, cex = 0.82) }
  mtext("A  where the KO-specific signal is lost", side = 3, adj = 0, line = 1.7, font = 2, cex = 0.82, at = -1.05); mtext("each data set reproduces itself; the signal is lost at the step to another culture, even of the same line (bar = mean, line = range)", side = 3, adj = 0, line = 0.5, cex = 0.55, col = INK2, at = -1.05)
  par(mar = c(4.6, 4.8, 3.6, 1)); boxplot(rho ~ same_block, bs, names = c("different library block", "same library block"), col = c("#d5d4cf", "#eda100"), border = INK2, outline = FALSE, las = 1, ylab = "similarity of two DIFFERENT KOs within a line", cex.axis = 0.75, xlab = "")
  points(jitter(as.numeric(bs$same_block) + 1, amount = 0.15), bs$rho, pch = 21, bg = adjustcolor(LINE_COL[bs$line], 0.8), col = INK2, cex = 0.8); legend("topleft", short, pt.bg = LINE_COL[LINES], pch = 21, col = INK2, bty = "n", cex = 0.62)
  mtext("B  library batch is nested in the KO", side = 3, adj = 0, line = 1.7, font = 2, cex = 0.82); mtext(sprintf("same - different = %.2f, permutation p = %s (labels permuted within line)", bobs, fp(max(bp_, 1 / 5001))), side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2)
  par(mar = c(4.6, 6, 3.6, 1)); image(1:4, 1:8, t(pmax(pmin(KD, 1.2), -1.2))[, 8:1], col = rev(DIV), zlim = c(-1.2, 1.2), axes = FALSE, xlab = "", ylab = ""); for (i in 1:8) for (j in 1:4) text(j, 9 - i, sprintf("%+.2f", KD[i, j]), cex = 0.75, font = ifelse(KD[i, j] >= 0.3, 2, 1))
  axis(1, 1:4, short, tick = FALSE, cex.axis = 0.8); text(0.4, 8:1, KOS, adj = 1, xpd = NA, font = 2, cex = 0.8, col = KO_COL[KOS]); mtext("C  target transcript after KO, log2 NTC / KO", side = 3, adj = 0, line = 1.7, font = 2, cex = 0.82)
  mtext(sprintf("blue = reduced; bold >= 0.3; %d of 32 KO x line show < 19%% reduction", sum(KD < 0.3)), side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2)
  par(mar = c(4.6, 4.8, 3.6, 1)); tm <- c("technical", "cell_cycle", "replicate", "KO", "residual"); V <- sapply(LINES, function(cl) vf$median_R2[match(tm, vf$term[vf$line == cl]) + min(which(vf$line == cl)) - 1])
  bp <- barplot(V, names.arg = short, col = c("#898781", "#e87ba4", "#d5d4cf", "#4a3aa7", "#efefef"), border = NA, las = 1, ylab = "median R2 per gene (sequential)", cex.names = 0.8, ylim = c(0, 1.25)); ch <- sapply(LINES, function(cl) vf$chance_R2[vf$line == cl & vf$term == "KO"])
  base <- colSums(V[1:3, ]); segments(bp - 0.5, base + ch, bp + 0.5, base + ch, col = "red", lwd = 2, lty = 2); legend("top", c("technical (% assigned, depth)", "cell cycle", "replicate", "KO group", "residual"), fill = c("#898781", "#e87ba4", "#d5d4cf", "#4a3aa7", "#efefef"), border = NA, bty = "n", cex = 0.58, ncol = 2)
  mtext("D  variance budget of a gene within a line", side = 3, adj = 0, line = 1.7, font = 2, cex = 0.82); mtext("red = KO term expected by chance with these group sizes (spec74)", side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2)
  k <- scr[grepl("^spec_", scr$what), ]; k$KO <- sub("spec_", "", k$what); par(mar = c(4.6, 4.8, 3.6, 1)); plot(k$n_cells, k$reliability, log = "x", pch = 21, bg = KO_COL[k$KO], col = INK2, cex = 2, ylim = c(-0.05, 1.05), bty = "n", las = 1, xlab = "cells per KO in the screen (log)", ylab = "screen reliability of the KO profile")
  text(k$n_cells, k$reliability, k$KO, pos = 3, cex = 0.68, font = 2, col = KO_COL[k$KO], xpd = NA); abline(h = 0.5, lty = 3, col = GRID); mtext("E  two screen KOs cannot be reproduced by anything", side = 3, adj = 0, line = 1.7, font = 2, cex = 0.82); mtext("split-half reliability (spec114): HMGA1 13 cells, SOX9 59 cells", side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2)
  mtext("Why the KO programmes reproduce weakly between CROP-seq and bulk RNA-seq", outer = TRUE, side = 3, adj = 0.01, line = 1.7, font = 2, cex = 1.02)
  mtext("KO-specific effect = KO minus the mean of the other 8 groups, on the 20 Cm module scores; neither normalisation nor covariates change this picture (spec115, 1,280 configurations)", outer = TRUE, side = 3, adj = 0.01, line = 0.5, cex = 0.62, col = INK2)
})
