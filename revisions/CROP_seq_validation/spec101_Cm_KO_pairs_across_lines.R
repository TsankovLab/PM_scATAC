# Spec 101: a Cm module in a given KO, CROP-seq and bulk across the 4 lines -- the fig14 (spec100) panels for any module / KO pairs.
# PAIRS_CM env: comma-separated Cm:KO (default Cm8:TWIST1,Cm9:TWIST1,Cm10:TCF3). Scores = AddModuleScore, top-50 module genes, z-scored within each data set (screen cells; bulk samples of a line).
# per data set: NTC - KO and other-KOs - KO (screen: Wilcoxon; bulk: Welch t, n/a with a single NTC or KO sample), rank of the KO among the 8 KOs for the largest drop; bulk lines combined: lm(z ~ line + KO)
# needs: spec72 (S70), spec75 (screen scores)
# outputs: spec101_Cm_KO_pairs.csv, spec101_Cm_KO_pairs_all_KOs.csv, fig15_Cm_KO_pairs_across_lines.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
PC <- strsplit(strsplit(Sys.getenv("PAIRS_CM", unset = "Cm8:TWIST1,Cm9:TWIST1,Cm10:TCF3"), ",")[[1]], ":"); GR <- c("NTC", KOS); fab <- load_full_audited_bulk(); meta <- fab$meta
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); ids <- pr$samples; sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); DS <- c("CROP-seq", LINES)
fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2)); SH <- function(l) sub("NCI-|-211H", "", l)
one <- function(cm, ko) { dat <- rbind(data.frame(ds = "CROP-seq", group = as.character(sc$group), v = sc$score[, cm]), data.frame(ds = meta[ids, "Cell_Line"], group = meta[ids, "Gene_Control"], v = pr$cm_scores$tmm[cm, ids])); dat <- dat[dat$group %in% GR, ]
  dat$z <- ave(dat$v, dat$ds, FUN = function(x) (x - mean(x)) / sd(x)); dat$cls <- ifelse(dat$group == "NTC", "NTC", ifelse(dat$group == ko, "KO", "other"))
  st <- do.call(rbind, lapply(DS, function(d) { a <- dat$z[dat$ds == d & dat$cls == "NTC"]; b <- dat$z[dat$ds == d & dat$cls == "KO"]; o <- dat$z[dat$ds == d & dat$cls == "other"]
    data.frame(Cm = cm, KO = ko, ds = d, n_NTC = length(a), n_KO = length(b), NTC_minus_KO = mean(a) - mean(b), p_vs_NTC = if (d == "CROP-seq") wilcox.test(a, b)$p.value else if (min(length(a), length(b)) >= 2) t.test(a, b)$p.value else NA,
               otherKO_minus_KO = mean(o) - mean(b), p_vs_otherKO = if (d == "CROP-seq") wilcox.test(o, b)$p.value else if (length(b) >= 2) t.test(o, b)$p.value else NA) }))
  bk <- dat[dat$ds != "CROP-seq", ]; f1 <- summary(lm(z ~ ds + I(cls == "KO"), bk[bk$cls != "other", ]))$coef; f2 <- summary(lm(z ~ ds + I(cls == "KO"), bk[bk$cls != "NTC", ]))$coef
  E <- do.call(rbind, lapply(DS, function(d) data.frame(Cm = cm, target = ko, ds = d, KO = KOS, NTC_minus_KO = sapply(KOS, function(k) mean(dat$z[dat$ds == d & dat$group == "NTC"]) - mean(dat$z[dat$ds == d & dat$group == k]))))); E$rank_most_down <- ave(-E$NTC_minus_KO, E$ds, FUN = rank)
  st$rank_most_down <- E$rank_most_down[E$KO == ko][match(st$ds, E$ds[E$KO == ko])]
  list(dat = dat, st = st, E = E, est = c(-f1[nrow(f1), 1], -f2[nrow(f2), 1]), p = c(f1[nrow(f1), 4], f2[nrow(f2), 4])) }
TAG <- ""
# PAIRS_CM=top:N -> rank all 160 Cm x KO pairs by agreement across the 4 bulk lines and plot the best N: pairs whose NTC - KO and other-KOs - KO both have the same sign in all 4 lines (KO-specific, not the shared any-KO response), ordered by the combined p (lm z ~ line + KO, NTC vs KO);
# the table (spec101_all_pairs_ranking.csv) also has the other-KOs - KO contrast and the CROP-seq effect
if (PC[[1]][1] == "top") { TAG <- "_top"; ALL <- do.call(rbind, lapply(paste0("Cm", 1:20), function(cm) do.call(rbind, lapply(KOS, function(ko) { r <- one(cm, ko); b <- r$st[r$st$ds != "CROP-seq", ]; s <- r$st[r$st$ds == "CROP-seq", ]
    data.frame(Cm = cm, KO = ko, bulk_NTC_minus_KO = r$est[1], bulk_p = r$p[1], n_lines_same_sign = max(sum(b$NTC_minus_KO > 0), sum(b$NTC_minus_KO < 0)), bulk_other_minus_KO = r$est[2], bulk_other_p = r$p[2],
               n_lines_same_sign_other = max(sum(b$otherKO_minus_KO > 0), sum(b$otherKO_minus_KO < 0)), screen_NTC_minus_KO = s$NTC_minus_KO, screen_p = s$p_vs_NTC, screen_same_direction = sign(s$NTC_minus_KO) == sign(r$est[1])) }))))
  ALL$bulk_fdr <- p.adjust(ALL$bulk_p, "BH"); ALL$consistent <- ALL$n_lines_same_sign == 4 & ALL$n_lines_same_sign_other == 4 & sign(ALL$bulk_NTC_minus_KO) == sign(ALL$bulk_other_minus_KO); ALL <- ALL[order(-ALL$consistent, ALL$bulk_p), ]; write.csv(ALL, file.path(OUT, "spec101_all_pairs_ranking.csv"), row.names = FALSE)
  options(width = 220); cat("pairs with the same sign in 4/4 lines for both NTC - KO and other KOs - KO:", sum(ALL$consistent), "of 160\n"); print(head(ALL, 20), digits = 2, row.names = FALSE); PC <- lapply(seq_len(as.integer(PC[[1]][2])), function(i) c(ALL$Cm[i], ALL$KO[i])) }
R <- lapply(PC, function(x) one(x[1], x[2])); write.csv(do.call(rbind, lapply(R, `[[`, "st")), file.path(OUT, sprintf("spec101_Cm_KO_pairs%s.csv", TAG)), row.names = FALSE); write.csv(do.call(rbind, lapply(R, `[[`, "E")), file.path(OUT, sprintf("spec101_Cm_KO_pairs_all_KOs%s.csv", TAG)), row.names = FALSE)
options(width = 200); for (r in R) { print(r$st, digits = 3, row.names = FALSE); cat(sprintf("bulk, 4 lines combined: NTC - KO = %+.2f SD, p = %.3g; other KOs - KO = %+.2f SD, p = %.3g\n\n", r$est[1], r$p[1], r$est[2], r$p[2])) }

save_fig(paste0("fig15_Cm_KO_pairs_across_lines", TAG), 13, 5.2 * length(R) + 0.6, function() {
  layout(matrix(seq_len(2 * length(R)), ncol = 2, byrow = TRUE), widths = c(1.25, 1)); par(oma = c(0, 0, 2.6, 0.4), family = "sans"); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.85); mtext(b, side = 3, adj = 0, line = 0.5, cex = 0.62, col = INK2) }
  for (r in R) { cm <- r$st$Cm[1]; ko <- r$st$KO[1]; dat <- r$dat; st <- r$st; E <- r$E; CC <- c(NTC = "grey45", KO = KO_COL[[ko]], other = "#d5d4cf"); par(mar = c(4.4, 4.6, 3.6, 1))
    plot(NA, xlim = c(0.4, 5.6), ylim = c(-3, 4.6), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = sprintf("%s score (z within data set)", cm), cex.axis = 0.8); abline(h = 0, lty = 2, col = GRID)
    for (i in seq_along(DS)) for (j in 1:3) { cl <- names(CC)[j]; v <- dat$z[dat$ds == DS[i] & dat$cls == cl]; at <- i + (j - 2) * 0.27
      if (length(v) >= 4) boxplot(v, at = at, add = TRUE, boxwex = 0.24, col = adjustcolor(CC[[cl]], 0.7), border = INK2, outline = FALSE, axes = FALSE)
      if (DS[i] != "CROP-seq") points(at + seq(-0.07, 0.07, length.out = length(v)), pmax(pmin(v, 3.2), -3), pch = 21, bg = CC[[cl]], col = INK2, cex = 0.9) else segments(at - 0.1, mean(v), at + 0.1, mean(v), lwd = 2, col = INK) }
    axis(1, at = 1:5, labels = c("CROP-seq\n(cells)", SH(LINES)), tick = FALSE, cex.axis = 0.78, padj = 0.5)
    text(1:5, 4.6, sprintf("NTC - %s %+.2f\np %s\nother - %s %+.2f\np %s", ko, st$NTC_minus_KO, ifelse(is.na(st$p_vs_NTC), "n/a (n = 1)", fp(st$p_vs_NTC)), ko, st$otherKO_minus_KO, ifelse(is.na(st$p_vs_otherKO), "n/a (n = 1)", fp(st$p_vs_otherKO))), cex = 0.56, adj = c(0.5, 1))
    legend("bottomright", c("NTC", paste(ko, "KO"), "other 7 KOs"), pt.bg = CC, pch = 22, pt.cex = 1.3, bty = "n", cex = 0.68)
    ttl(sprintf("%s in NTC, %s KO and the other KOs", cm, ko), sprintf("bulk, 4 lines combined: NTC - %s %+.2f SD, p = %s; other KOs - %s %+.2f SD, p = %s", ko, r$est[1], fp(r$p[1]), ko, r$est[2], fp(r$p[2])))
    plot(NA, xlim = c(0.5, 5.5), ylim = range(E$NTC_minus_KO) + c(-0.5, 0.4), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = sprintf("%s, NTC - KO (z units; > 0 = lower in the KO)", cm), cex.axis = 0.8); abline(h = 0, lty = 2, col = INK)
    for (i in seq_along(DS)) { e <- E[E$ds == DS[i], ]; points(rep(i, 8) + ifelse(e$KO == ko, 0, seq(-0.18, 0.18, length.out = 8)), e$NTC_minus_KO, pch = 21, bg = KO_COL[e$KO], col = INK2, cex = ifelse(e$KO == ko, 2, 1.1)); text(i, max(E$NTC_minus_KO) + 0.3, sprintf("%s #%g of 8", ko, e$rank_most_down[e$KO == ko]), cex = 0.62, font = 2, col = KO_COL[[ko]]) }
    axis(1, at = 1:5, labels = c("CROP-seq", SH(LINES)), tick = FALSE, cex.axis = 0.78); legend("bottomleft", KOS, pt.bg = KO_COL[KOS], pch = 21, pt.cex = 1.3, bty = "n", cex = 0.7, ncol = 4)
    ttl(sprintf("%s effect of every KO", cm), sprintf("large dot = %s; # = rank among the 8 KOs (1 = largest drop)", ko)) }
  mtext("Cm modules in their KO: CROP-seq and bulk RNA-seq across the 4 cell lines (AddModuleScore, top-50 module genes; S70)", outer = TRUE, side = 3, adj = 0, line = 1.2, font = 2, cex = 0.95)
})
