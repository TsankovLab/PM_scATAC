# Spec 15: do samples / guides / cells with STRONGER knockdown of their target gene reproduce the other dataset better?
#  (A) bulk, per replicate sample: knockdown = log2 CP10K of the target gene minus its NTC mean; the sample's cc-adjusted profile (replicate-level t, centred) is correlated with the
#      merged CROP-seq KO profiles; matched score = matched screen KO's r standardised across the 8 screen KOs (z), and its rank.
#  (B) CROP-seq per GUIDE (single-guide cells only, not merged per gene; guides with >= 30 cells): knockdown of the target (t vs NTC), guide-vs-guide reproducibility (same KO vs different KO),
#      and match of each guide's cc-regressed profile to the cc-adjusted bulk KOs (matched bulk KO's r standardised across the line's bulk KOs, averaged over lines).
#  (C) CROP-seq per CELL: within a KO, cells with the target gene undetected vs detected (low vs high expression); their cc-regressed profiles vs bulk, target gene excluded.
# All profiles: cc-regressed statistics (screen: Welch t after pooled cc regression; bulk: pooled-variance t after cc adjustment), shared response removed by centring.
# needs: spec1_prepare_stats.R
# outputs: spec15_bulk_samples.csv, spec15_guides.csv, spec15_guide_pairs.csv, spec15_cells.csv, fig19_knockdown_vs_match.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "spec_variants.R"))
set.seed(15); ctx <- load_spec_context(LINES); ss <- ctx$ss; bs <- ctx$bs; ALLC <- c(ORDER, EXTRA)
genes <- Reduce(intersect, c(list(ss$genes), lapply(bs[LINES], `[[`, "genes")))
srt <- load_screen(); md <- srt@meta.data; Yall <- GetAssayData(srt, layer = "data"); Y <- Yall[genes, ]
keepG1 <- md$Phase == "G1"; G1r <- screen_ko_stats(Y[, keepG1], md$merged_call[keepG1], md$cc[keepG1], TRUE)
zrow <- function(r, j) (r[j] - mean(r)) / sd(r)

## ================================================================ (A) bulk per replicate
REP <- list()
for (cl in LINES) { m <- ctx$meta[ctx$meta$Cell_Line == cl, ]; s <- rownames(m); y <- ctx$lg[genes, s, drop = FALSE]; cc <- ctx$ccb$cc[match(s, ctx$ccb$sample)]
  g <- factor(m$Gene_Control, levels = c("NTC", sort(setdiff(unique(m$Gene_Control), "NTC"), method = "radix"))); X <- model.matrix(~ cc + g); q <- qr(X)
  s2 <- pmax(colSums(qr.resid(q, t(y))^2) / (length(s) - q$rank), 1e-4); yadj <- y - outer(qr.coef(q, t(y))["cc", ], cc - mean(cc))
  ntc_s <- s[m$Gene_Control == "NTC"]; ko_s <- s[m$Gene_Control != "NTC"]; Tm <- (yadj[, ko_s, drop = FALSE] - rowMeans(yadj[, ntc_s, drop = FALSE])) / sqrt(s2 * (1 + 1 / length(ntc_s))); Tm[!is.finite(Tm)] <- 0
  REP[[cl]] <- list(t = Tm - rowMeans(Tm), ko = setNames(as.character(m[ko_s, "Gene_Control"]), ko_s), ntc = ntc_s) }
POW <- function(n) intersect(ORDER, names(n)[n >= 50]); n_g1 <- G1r$n
SV <- list(`all cells, cc regressed` = list(S = prep_dataset(ss$tcc[genes, ALLC], TRUE, FALSE)[, ORDER], n = ss$n), `G1 only, cc regressed` = list(S = prep_dataset(G1r$t[genes, ALLC], TRUE, FALSE)[, ORDER], n = n_g1))
bulk <- list()
for (cl in LINES) { ko <- REP[[cl]]$ko; ko <- ko[ko %in% ORDER]; for (smp in names(ko)) { tf <- ko[[smp]]; kd <- ctx$lg[tf, smp] - mean(ctx$lg[tf, REP[[cl]]$ntc])
  row <- data.frame(line = cl, TF = tf, sample = smp, knockdown_log2 = kd)
  for (sv in names(SV)) { r <- cor(REP[[cl]]$t[, smp], SV[[sv]]$S); names(r) <- ORDER; pw <- POW(SV[[sv]]$n)
    row[[paste0("z_", sv)]] <- zrow(r, tf); row[[paste0("r_", sv)]] <- r[[tf]]; row[[paste0("rank_", sv)]] <- sum(r > r[[tf]]) + 1
    row[[paste0("zpow_", sv)]] <- if (tf %in% pw) zrow(r[pw], tf) else NA }
  bulk[[length(bulk) + 1]] <- row } }
bulk <- do.call(rbind, bulk); write.csv(bulk, file.path(OUT, "spec15_bulk_samples.csv"), row.names = FALSE)
cat(sprintf("=== (A) bulk replicates: %d samples ===\n", nrow(bulk))); options(width = 200)
tfsum <- aggregate(cbind(knockdown_log2, `z_all cells, cc regressed`) ~ TF, bulk, mean); names(tfsum)[3] <- "mean_matched_z"; print(format(tfsum, digits = 3), row.names = FALSE)
cortab <- list()
for (sv in names(SV)) for (p in c("z", "zpow")) { z <- bulk[[paste0(p, "_", sv)]]; ok <- !is.na(z); ct_ <- cor.test(bulk$knockdown_log2[ok], z[ok], method = "spearman", exact = FALSE); strong <- bulk$knockdown_log2[ok] < -0.5
  cortab[[length(cortab) + 1]] <- data.frame(screen = sv, columns = ifelse(p == "z", "all 8 KOs", "KOs with >=50 cells"), n = sum(ok), rho_knockdown_vs_z = unname(ct_$estimate), p = ct_$p.value, n_strong = sum(strong), z_strong = mean(z[ok][strong]), z_weak = mean(z[ok][!strong]),
                                             wilcox_p = tryCatch(wilcox.test(z[ok][strong], z[ok][!strong])$p.value, error = function(e) NA)) }
cortab <- do.call(rbind, cortab)
# within (line, KO) replicate pairs: does the replicate with stronger knockdown match better?
pr <- do.call(rbind, lapply(split(bulk, paste(bulk$line, bulk$TF)), function(d) if (nrow(d) == 2) { d <- d[order(d$knockdown_log2), ]; data.frame(line = d$line[1], TF = d$TF[1], kd_stronger = d$knockdown_log2[1], kd_weaker = d$knockdown_log2[2], z_stronger = d$`z_all cells, cc regressed`[1], z_weaker = d$`z_all cells, cc regressed`[2]) }))
pr$stronger_better <- pr$z_stronger > pr$z_weaker; big <- abs(pr$kd_stronger - pr$kd_weaker) > 0.3
cat(sprintf("\nreplicate pairs: stronger-knockdown replicate matches better in %d of %d pairs (binomial p = %.2f); restricted to pairs whose knockdown differs by > 0.3 log2: %d of %d\n", sum(pr$stronger_better), nrow(pr), binom.test(sum(pr$stronger_better), nrow(pr))$p.value, sum(pr$stronger_better[big]), sum(big)))

## ================================================================ bulk reference profiles (cc-adjusted DESeq2 Wald, centred across the line's KOs)
BREF <- setNames(lapply(LINES, function(cl) { m <- bs[[cl]]$statcc[genes, ]; m[!is.finite(m)] <- 0; m - rowMeans(m) }), LINES)
base <- rowMeans(ss$tcc[genes, ALLC])                                         # shared screen response removed from guide / subgroup profiles
match_bulk <- function(prof, tf, gset = genes) { out <- sapply(LINES, function(cl) { B <- BREF[[cl]][gset, intersect(ORDER, colnames(BREF[[cl]])), drop = FALSE]; if (!tf %in% colnames(B) || ncol(B) < 6) return(NA_real_)
  r <- cor(prof[gset], B); zrow(r[1, ], match(tf, colnames(B))) }); out }

## ================================================================ (B) CROP-seq per guide
cc_all <- md$cc; is_single <- !grepl("\\|", md$crispr_calls); ntc_cells <- md$crispr_calls %in% c("NTC_1", "NTC_3")
gtab <- table(md$crispr_calls[is_single & !ntc_cells]); gids <- names(gtab)[gtab >= 30 & sub("_.*", "", names(gtab)) %in% ORDER]
sel <- (is_single & md$crispr_calls %in% gids) | ntc_cells; lab <- ifelse(ntc_cells, "NTC", md$crispr_calls)
GS <- screen_ko_stats(Y[, sel], lab[sel], cc_all[sel], TRUE)
guides <- do.call(rbind, lapply(gids, function(gd) { tf <- sub("_.*", "", gd); ix <- which(lab == gd)
  kd <- mean(Yall[tf, ix]) - mean(Yall[tf, which(lab == "NTC")]); prof <- GS$t[, gd] - base[rownames(GS$t)]; mz <- match_bulk(prof, tf, setdiff(genes, tf))
  data.frame(guide = gd, TF = tf, n_cells = length(ix), target_logexpr_diff = kd, target_t = if (tf %in% rownames(GS$t)) GS$t[tf, gd] else NA, detect_frac = mean(Yall[tf, ix] > 0), detect_frac_NTC = mean(Yall[tf, which(lab == "NTC")] > 0),
             matched_z_mean = mean(mz, na.rm = TRUE), setNames(as.list(mz), paste0("z_", LINES))) }))
write.csv(guides, file.path(OUT, "spec15_guides.csv"), row.names = FALSE); cat("\n=== (B) CROP-seq per guide ===\n"); print(format(guides[, c("guide", "n_cells", "target_logexpr_diff", "target_t", "detect_frac", "detect_frac_NTC", "matched_z_mean")], digits = 3), row.names = FALSE)
ct <- cor.test(guides$target_t, guides$matched_z_mean, method = "spearman", exact = FALSE); cat(sprintf("\nSpearman (guide target t vs matched z, n=%d): rho = %.2f (p = %.2f); negative t = stronger knockdown, so a NEGATIVE rho means stronger knockdown matches better\n", nrow(guides), ct$estimate, ct$p.value))
cat(sprintf("merged-per-gene screen KOs, same cc-regressed profile: mean matched z of the same KOs = %.2f | mean over these guides = %.2f\n", mean(sapply(unique(guides$TF), function(tf) mean(match_bulk(ss$tcc[genes, tf] - base, tf, setdiff(genes, tf)), na.rm = TRUE))), mean(guides$matched_z_mean)))
pg <- guides[guides$TF %in% names(which(table(guides$TF) >= 2)), ]; pairs_kd <- do.call(rbind, lapply(split(pg, pg$TF), function(d) { d <- d[order(d$target_t), ]; data.frame(TF = d$TF[1], strongest = d$guide[1], weakest = d$guide[nrow(d)], t_strong = d$target_t[1], t_weak = d$target_t[nrow(d)], z_strong = d$matched_z_mean[1], z_weak = d$matched_z_mean[nrow(d)]) }))
cat("\nwithin each KO: guide with the stronger knockdown (lowest target t) vs the weakest one\n"); print(format(pairs_kd, digits = 3), row.names = FALSE); cat(sprintf("stronger-knockdown guide matches better in %d of %d KOs\n", sum(pairs_kd$z_strong > pairs_kd$z_weak), nrow(pairs_kd)))
P <- GS$t[, gids] - base[rownames(GS$t)]; Cg <- cor(P); gp <- which(upper.tri(Cg), arr.ind = TRUE)
gpairs <- data.frame(g1 = gids[gp[, 1]], g2 = gids[gp[, 2]], r = Cg[gp]); gpairs$same_KO <- sub("_.*", "", gpairs$g1) == sub("_.*", "", gpairs$g2)
auc <- function(w, o) { r <- rank(c(w, o)); (sum(r[seq_along(w)]) - length(w) * (length(w) + 1) / 2) / (length(w) * length(o)) }
write.csv(gpairs, file.path(OUT, "spec15_guide_pairs.csv"), row.names = FALSE)
cat(sprintf("\nguide-vs-guide reproducibility (cc-regressed t, shared response removed): same-KO guide pairs n=%d median r = %.2f | different-KO pairs n=%d median r = %.2f | AUC = %.2f\n", sum(gpairs$same_KO), median(gpairs$r[gpairs$same_KO]), sum(!gpairs$same_KO), median(gpairs$r[!gpairs$same_KO]), auc(gpairs$r[gpairs$same_KO], gpairs$r[!gpairs$same_KO])))
print(format(gpairs[gpairs$same_KO, c("g1", "g2", "r")], digits = 3), row.names = FALSE)

## ================================================================ (C) CROP-seq per cell: target gene undetected vs detected
cells <- list(); ko_cells <- c("BPTF", "PITX1", "TCF3", "TEAD4", "TWIST1", "SOX9")
for (tf in ko_cells) { idx <- which(md$merged_call == tf); ntc_i <- which(md$merged_call == "NTC"); det <- Yall[tf, idx] > 0; nlo <- sum(!det); nhi <- sum(det); if (min(nlo, nhi) < 15) { cat(tf, ": too few cells in a group (undetected", nlo, ", detected", nhi, ")\n"); next }
  cells_i <- c(idx, ntc_i); labs <- c(ifelse(det, paste0(tf, "_hi"), paste0(tf, "_lo")), rep("NTC", length(ntc_i))); gs <- setdiff(genes, tf); St <- screen_ko_stats(Y[gs, cells_i], labs, cc_all[cells_i], TRUE)
  for (grp in c("lo", "hi")) { nm <- paste0(tf, "_", grp); prof <- St$t[, nm] - base[gs]; mz <- match_bulk(prof, tf, gs); cells[[length(cells) + 1]] <- data.frame(TF = tf, group = ifelse(grp == "lo", "target undetected", "target detected"), n_cells = sum(labs == nm), mean_target_expr = mean(Yall[tf, cells_i[labs == nm]]), matched_z_mean = mean(mz, na.rm = TRUE), setNames(as.list(mz), paste0("z_", LINES))) } }
cells <- do.call(rbind, cells); write.csv(cells, file.path(OUT, "spec15_cells.csv"), row.names = FALSE); cat("\n=== (C) CROP-seq per cell: target undetected (lower expression) vs detected ===\n"); print(format(cells[, c("TF", "group", "n_cells", "mean_target_expr", "matched_z_mean")], digits = 3), row.names = FALSE)
w <- reshape(cells[, c("TF", "group", "matched_z_mean")], idvar = "TF", timevar = "group", direction = "wide"); names(w) <- c("TF", "z_undetected", "z_detected"); cat(sprintf("\nundetected cells match better than detected in %d of %d KOs (mean z %.2f vs %.2f)\n", sum(w$z_undetected > w$z_detected), nrow(w), mean(w$z_undetected), mean(w$z_detected)))

## ================================================================ figure
save_fig("fig19_knockdown_vs_match", 15, 9.4, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4, 0.5), mar = c(4.4, 4.6, 3.6, 1), family = "sans")
  sv <- "all cells, cc regressed"; zc <- paste0("z_", sv); pc <- LINE_COL[bulk$line]
  plot(bulk$knockdown_log2, bulk[[zc]], pch = 21, bg = pc, col = "white", cex = 1.3, bty = "n", las = 1, xlab = "target-gene knockdown in the bulk sample (log2 CP10K vs NTC)", ylab = "matched score z (vs CROP-seq)", cex.lab = 0.8, cex.axis = 0.75); abline(h = 0, v = 0, col = GRID)
  cs <- cortab[cortab$screen == sv & cortab$columns == "all 8 KOs", ]; mtext("A  Bulk replicates: knockdown vs match to the screen", side = 3, adj = 0, line = 1.7, font = 2, cex = 0.8); mtext(sprintf("n=%d samples; Spearman rho = %+.2f (p = %.2f); rho < 0 = stronger knockdown matches better", cs$n, cs$rho_knockdown_vs_z, cs$p), side = 3, adj = 0, line = 0.6, cex = 0.55, col = INK2)
  legend("topleft", sub("NCI-", "", LINES), pt.bg = LINE_COL[LINES], pch = 21, col = "white", bty = "n", cex = 0.65)
  strong <- bulk$knockdown_log2 < -0.5; vv <- list(bulk[[zc]][strong], bulk[[zc]][!strong]); boxplot(vv, names = c(sprintf("strong knockdown\n(n=%d)", sum(strong)), sprintf("weak / none\n(n=%d)", sum(!strong))), col = c("#eb6834", GREY), border = INK2, medcol = INK, outline = FALSE, las = 1, ylab = "matched score z", cex.axis = 0.75, cex.lab = 0.8, boxwex = 0.6, whisklty = 1, staplelty = 1)
  for (i in 1:2) points(i + runif(length(vv[[i]]), -0.2, 0.2), vv[[i]], pch = 16, cex = 0.5, col = adjustcolor(INK, 0.5)); abline(h = 0, col = GRID)
  mtext("B  Strong (< -0.5 log2) vs weak knockdown", side = 3, adj = 0, line = 1.7, font = 2, cex = 0.8); mtext(sprintf("mean z %.2f vs %.2f; Wilcoxon p = %.2f", mean(vv[[1]]), mean(vv[[2]]), wilcox.test(vv[[1]], vv[[2]])$p.value), side = 3, adj = 0, line = 0.6, cex = 0.55, col = INK2)
  plot(NA, xlim = c(0.5, 2.5), ylim = range(c(pr$z_stronger, pr$z_weaker)) + c(-0.2, 0.2), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "matched score z", cex.axis = 0.75, cex.lab = 0.8); axis(1, at = 1:2, labels = c("stronger-knockdown\nreplicate", "weaker-knockdown\nreplicate"), tick = FALSE, cex.axis = 0.75, line = 0.6)
  segments(1, pr$z_stronger, 2, pr$z_weaker, col = ifelse(pr$stronger_better, "#eb6834", "#2a78d6"), lwd = 1.3); points(rep(1, nrow(pr)), pr$z_stronger, pch = 21, bg = LINE_COL[pr$line], col = "white", cex = 1.1); points(rep(2, nrow(pr)), pr$z_weaker, pch = 21, bg = LINE_COL[pr$line], col = "white", cex = 1.1)
  mtext("C  The two replicates of a KO", side = 3, adj = 0, line = 1.7, font = 2, cex = 0.8); mtext(sprintf("stronger replicate matches better in %d of %d pairs (orange)", sum(pr$stronger_better), nrow(pr)), side = 3, adj = 0, line = 0.6, cex = 0.55, col = INK2)
  tcol <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7", "#008300", "#b3400f"), ORDER)
  plot(guides$target_t, guides$matched_z_mean, pch = 21, bg = tcol[guides$TF], col = "white", cex = 2, bty = "n", las = 1, xlab = "guide: target-gene t vs NTC (negative = knockdown)", ylab = "matched score z (vs bulk)", cex.lab = 0.8, cex.axis = 0.75); abline(h = 0, v = 0, col = GRID); text(guides$target_t, guides$matched_z_mean, sub("_", "-", guides$guide), pos = 3, cex = 0.55)
  mtext("D  CROP-seq individual guides (not merged)", side = 3, adj = 0, line = 1.7, font = 2, cex = 0.8); mtext(sprintf("n=%d guides; Spearman rho = %+.2f (p = %.2f); rho < 0 = stronger knockdown matches better", nrow(guides), ct$estimate, ct$p.value), side = 3, adj = 0, line = 0.6, cex = 0.55, col = INK2)
  boxplot(list(gpairs$r[gpairs$same_KO], gpairs$r[!gpairs$same_KO]), names = c(sprintf("same KO\n(n=%d pairs)", sum(gpairs$same_KO)), sprintf("different KO\n(n=%d)", sum(!gpairs$same_KO))), col = c("#4a3aa7", GREY), border = INK2, medcol = INK, outline = FALSE, las = 1, ylab = "guide-vs-guide r", cex.axis = 0.75, cex.lab = 0.8, boxwex = 0.6, whisklty = 1, staplelty = 1)
  points(1 + runif(sum(gpairs$same_KO), -0.2, 0.2), gpairs$r[gpairs$same_KO], pch = 21, bg = "#4a3aa7", col = "white", cex = 1.1); points(2 + runif(sum(!gpairs$same_KO), -0.2, 0.2), gpairs$r[!gpairs$same_KO], pch = 16, cex = 0.4, col = adjustcolor(INK2, 0.4)); abline(h = 0, col = GRID)
  mtext("E  Guide-level reproducibility in the screen", side = 3, adj = 0, line = 1.7, font = 2, cex = 0.8); mtext(sprintf("AUC = %.2f (same-KO guide pair correlates higher than a different-KO pair)", auc(gpairs$r[gpairs$same_KO], gpairs$r[!gpairs$same_KO])), side = 3, adj = 0, line = 0.6, cex = 0.55, col = INK2)
  plot(NA, xlim = c(0.5, 2.5), ylim = range(cells$matched_z_mean) + c(-0.2, 0.2), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "matched score z (vs bulk)", cex.axis = 0.75, cex.lab = 0.8); axis(1, at = 1:2, labels = c("target undetected\n(lower expression)", "target detected"), tick = FALSE, cex.axis = 0.75, line = 0.6); abline(h = 0, col = GRID)
  for (tf in unique(cells$TF)) { a <- cells$matched_z_mean[cells$TF == tf & cells$group == "target undetected"]; b <- cells$matched_z_mean[cells$TF == tf & cells$group == "target detected"]; segments(1, a, 2, b, col = tcol[tf], lwd = 1.5); points(c(1, 2), c(a, b), pch = 21, bg = tcol[tf], col = "white", cex = 1.4); text(2.06, b, tf, adj = 0, cex = 0.6, xpd = NA) }
  mtext("F  CROP-seq cells within a KO", side = 3, adj = 0, line = 1.7, font = 2, cex = 0.8); mtext(sprintf("undetected cells match better in %d of %d KOs", sum(w$z_undetected > w$z_detected), nrow(w)), side = 3, adj = 0, line = 0.6, cex = 0.55, col = INK2)
  mtext("Do stronger knockdowns reproduce the other dataset better?", outer = TRUE, side = 3, adj = 0, line = 2.4, font = 2, cex = 1)
  mtext("A-C bulk vs merged CROP-seq; D-F CROP-seq guides / cells vs cc-adjusted bulk (avg of lines); all profiles cc adjusted and centred on the shared response; z = matched KO's r standardised across the 8 KOs (0 = no better than other KOs)", outer = TRUE, side = 3, adj = 0, line = 1.1, cex = 0.55, col = INK2)
})
