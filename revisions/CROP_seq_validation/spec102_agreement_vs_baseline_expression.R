# Spec 102: does the agreement of a KO in a cell line depend on the baseline expression of the KO target in that line?
# per KO x line (32 points):
#   agreement  = mean Spearman of the KO profile of this line with the SAME KO in the other 3 lines, minus the mean with the OTHER KOs of those lines (primary configuration: S70, TMM + limma t of NTC - KO,
#                KO-centred, top-500 KO-signal genes; spec76 matrices), no covariate and cell cycle regressed
#   baseline   = mean TMM logCPM of the target gene in the NTC samples of the line;   knock-down = NTC - KO of the target gene (logCPM) in the line
# tests: Spearman over the 32 points; within KO (both variables centred per KO: is the line that expresses the TF more the one that agrees more?); within line (centred per line); lm agreement ~ KO + line + baseline
# H2052 only (8 points): own-KO Spearman CROP-seq vs bulk (spec100 matrix, genes variable in both, 9-group centring) vs baseline expression in bulk H2052 NTC
# needs: spec72 (S70), spec76 (SET=S70), spec100
# outputs: spec102_agreement_vs_baseline.csv, spec102_tests.csv, fig16_agreement_vs_baseline_expression.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
fab <- load_full_audited_bulk(); meta <- fab$meta; pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); ids <- pr$samples; SH <- function(l) sub("NCI-|-211H", "", l); fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
stopifnot(all(KOS %in% rownames(fab$cs))); y <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(fab$cs[, ids])), log = TRUE, prior.count = 2)[KOS, ]
D <- do.call(rbind, lapply(LINES, function(cl) do.call(rbind, lapply(KOS, function(k) { s <- ids[meta[ids, "Cell_Line"] == cl]; g <- meta[s, "Gene_Control"]; data.frame(KO = k, line = cl, baseline = mean(y[k, s[g == "NTC"]]), knockdown = mean(y[k, s[g == "NTC"]]) - mean(y[k, s[g == k]]), n_KO = sum(g == k)) }))))
agr <- function(cv) { z <- read.csv(file.path(OUT, sprintf("spec76_KOline_correlation_signal500_%s_S70.csv", cv)), check.names = FALSE); M <- as.matrix(z[, -1]); rn <- z[, 1]; ko <- sub(" .*", "", rn); ln <- sub(".* ", "", rn)
  t(sapply(seq_len(nrow(D)), function(i) { li <- SH(D$line[i]); r <- which(ko == D$KO[i] & ln == li); stopifnot(length(r) == 1); o <- ln != li; c(same = mean(M[r, o & ko == D$KO[i]]), other = mean(M[r, o & ko != D$KO[i]])) })) }
a <- agr("none"); D$same_none <- a[, "same"]; D$agree_none <- a[, "same"] - a[, "other"]; a <- agr("cc"); D$same_cc <- a[, "same"]; D$agree_cc <- a[, "same"] - a[, "other"]
C <- as.matrix(read.csv(file.path(OUT, "spec100_H2052_KOxKO_int2000.csv"), row.names = 1, check.names = FALSE)); D$crop_own <- ifelse(D$line == "NCI-H2052", diag(C)[D$KO], NA)
write.csv(D, file.path(OUT, "spec102_agreement_vs_baseline.csv"), row.names = FALSE)
cen <- function(v, g) v - ave(v, g); ct <- function(x, yv) { r <- suppressWarnings(cor.test(x, yv, method = "spearman", exact = FALSE)); c(rho = unname(r$estimate), p = r$p.value) }
TS <- do.call(rbind, lapply(c("baseline", "knockdown"), function(xv) do.call(rbind, lapply(c("agree_none", "agree_cc"), function(av) { x <- D[[xv]]; v <- D[[av]]; f <- summary(lm(v ~ D$KO + D$line + x))$coef
  rbind(data.frame(x = xv, agreement = av, test = "all 32 KO x line", t(ct(x, v))), data.frame(x = xv, agreement = av, test = "within KO (centred per KO)", t(ct(cen(x, D$KO), cen(v, D$KO)))),
        data.frame(x = xv, agreement = av, test = "within line (centred per line)", t(ct(cen(x, D$line), cen(v, D$line)))), data.frame(x = xv, agreement = av, test = "lm ~ KO + line + x (slope per log2 unit)", rho = f["x", 1], p = f["x", 4])) }))))
h <- D[D$line == "NCI-H2052", ]; TS <- rbind(TS, data.frame(x = c("baseline", "knockdown"), agreement = "CROP-seq vs bulk H2052 own-KO rho", test = "8 KOs", rbind(ct(h$baseline, h$crop_own), ct(h$knockdown, h$crop_own))))
write.csv(TS, file.path(OUT, "spec102_tests.csv"), row.names = FALSE); options(width = 200); print(TS, digits = 3, row.names = FALSE)
cat("\nbaseline logCPM (NTC) / knock-down (NTC - KO) / agreement (no covariate)\n"); for (v in c("baseline", "knockdown", "agree_none")) { cat("--", v, "\n"); print(round(t(sapply(KOS, function(k) setNames(D[[v]][D$KO == k][match(LINES, D$line[D$KO == k])], SH(LINES)))), 2)) }

save_fig("fig16_agreement_vs_baseline_expression", 14, 9.4, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4.2, 0.4), mar = c(4.6, 4.6, 3.6, 1), family = "sans"); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.8); mtext(b, side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2) }
  PCH <- setNames(c(21, 22, 24, 23), LINES); tt <- function(xv, av, tst) { r <- TS[TS$x == xv & TS$agreement == av & TS$test == tst, ]; sprintf("rho %.2f, p %s", r$rho, fp(r$p)) }
  sc <- function(x, v, xl, yl, lab = TRUE) { plot(x, v, pch = PCH[D$line], bg = KO_COL[D$KO], col = INK2, cex = 1.5, bty = "n", las = 1, xlab = xl, ylab = yl, cex.axis = 0.8); abline(h = 0, lty = 2, col = GRID); abline(lm(v ~ x), col = INK2, lty = 3) }
  sc(D$baseline, D$agree_none, "target expression in NTC (TMM logCPM)", "agreement with the other 3 lines (same KO - other KOs)"); legend("topleft", SH(LINES), pch = PCH, pt.bg = "grey70", bty = "n", cex = 0.65); legend("bottomright", KOS, pch = 21, pt.bg = KO_COL[KOS], bty = "n", cex = 0.55, ncol = 2)
  ttl("A  agreement vs baseline expression (no covariate)", sprintf("32 KO x line: %s; within line: %s", tt("baseline", "agree_none", "all 32 KO x line"), tt("baseline", "agree_none", "within line (centred per line)")))
  x <- cen(D$baseline, D$KO); v <- cen(D$agree_none, D$KO); sc(x, v, "target expression in NTC, relative to the KO mean over lines", "agreement, relative to the KO mean over lines"); abline(v = 0, lty = 2, col = GRID)
  ttl("B  within KO: does the line expressing the TF more agree more?", tt("baseline", "agree_none", "within KO (centred per KO)"))
  sc(D$baseline, D$agree_cc, "target expression in NTC (TMM logCPM)", "agreement with the other 3 lines, cell cycle regressed")
  ttl("C  agreement vs baseline expression (cell cycle regressed)", sprintf("32 KO x line: %s; within KO: %s", tt("baseline", "agree_cc", "all 32 KO x line"), tt("baseline", "agree_cc", "within KO (centred per KO)")))
  sc(D$knockdown, D$agree_none, "knock-down of the target, NTC - KO (log2)", "agreement with the other 3 lines (same KO - other KOs)"); abline(v = 0, lty = 2, col = GRID)
  ttl("D  agreement vs knock-down of the target transcript", sprintf("32 KO x line: %s; within KO: %s", tt("knockdown", "agree_none", "all 32 KO x line"), tt("knockdown", "agree_none", "within KO (centred per KO)")))
  plot(D$baseline, D$knockdown, pch = PCH[D$line], bg = KO_COL[D$KO], col = INK2, cex = 1.5, bty = "n", las = 1, xlab = "target expression in NTC (TMM logCPM)", ylab = "knock-down, NTC - KO (log2)", cex.axis = 0.8); abline(h = 0, lty = 2, col = GRID); r <- ct(D$baseline, D$knockdown)
  ttl("E  knock-down vs baseline expression", sprintf("rho %.2f, p %s", r[1], fp(r[2])))
  plot(h$baseline, h$crop_own, pch = 21, bg = KO_COL[h$KO], col = INK2, cex = 1.8, bty = "n", las = 1, xlab = "target expression in bulk H2052 NTC (TMM logCPM)", ylab = "own-KO Spearman, CROP-seq vs bulk H2052", cex.axis = 0.8, xlim = range(h$baseline) + c(-0.3, 0.5)); abline(h = 0, lty = 2, col = GRID); text(h$baseline, h$crop_own, h$KO, pos = 4, cex = 0.62, col = KO_COL[h$KO], font = 2)
  r <- TS[TS$agreement == "CROP-seq vs bulk H2052 own-KO rho" & TS$x == "baseline", ]; ttl("F  H2052: CROP-seq match vs baseline expression", sprintf("8 KOs: rho %.2f, p %s (genes variable in both, 9-group centring)", r$rho, fp(r$p)))
  mtext("Does the agreement of a KO in a cell line depend on how much that line expresses the target?", outer = TRUE, side = 3, adj = 0, line = 2.6, font = 2, cex = 1)
  mtext("agreement = mean Spearman with the same KO in the other 3 lines minus with the other KOs (S70, KO-centred limma t, top-500 KO-signal genes); colour = KO, symbol = line", outer = TRUE, side = 3, adj = 0, line = 1.4, cex = 0.65, col = INK2)
})
