# Spec 83: adjust every cell line for the sarcomatoid score (scS = Cm17 module score: AXL, VIM, THY1, FBN1, SERPINE2, ...) and recompute the similarities. S70 set.
# Per line: TMM logCPM ~ group + covariate(s), limma moderated t of NTC - KO (agree_common.R ko_effects), covariates centred within the line:
#   none | Cm17 | cc | cc + Cm17      (Cm17 = the sample's AddModuleScore on TMM logCPM, top-50 genes, from spec72; cc = S + G2M index)
# NB: Cm17 differs between KO groups (R2 0.5-0.9 in bulk), so the adjustment also removes the part of a KO effect that lies along the sarcomatoid axis.
#   A across lines : agreement() of KO-centred profiles (Spearman, KO-label permutation, as spec76) on KO-signal genes (selected from the 2 other lines, recomputed per
#                    adjustment) | CROP-seq genes | Cm20_genes | Hallmark; per-TF p
#   B vs CROP-seq  : 8 KO x 20 Cm module effects per line vs the screen (cell module score, NTC - KO median; the screen gets the SAME adjustment: score regressed on the
#                    cell's Cm17 and/or cell-cycle score); uncentred and KO-centred; plus the 4-line mean
# needs: spec72_build_profiles.R (S70), spec75 (spec75_screen_scores_by_cell.rds), spec71 audit
# outputs: spec83_across_lines.csv, spec83_across_lines_per_KO.csv, spec83_vs_cropseq.csv, fig9E_Cm17_adjustment.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000")); perms <- make_perms(N_PERM)
fab <- load_full_audited_bulk(); meta <- fab$meta; cs <- fab$cs; q <- load_audit()
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); ids <- pr$samples; genes <- pr$genes; top50 <- cm_top50(); CM <- names(top50)
cm_sc <- pr$cm_scores$tmm; ss <- readRDS(file.path(OUT, "spec_screen_stats.rds"))
ADJ <- c(none = "no adjustment", Cm17 = "Cm17 (scS)", cc = "cell cycle", cc_Cm17 = "cell cycle + Cm17")
covX <- function(s, adj) { X <- switch(adj, none = NULL, Cm17 = cbind(Cm17 = cm_sc["Cm17", s]), cc = cbind(cc = q[s, "cc"]), cc_Cm17 = cbind(cc = q[s, "cc"], Cm17 = cm_sc["Cm17", s])); if (is.null(X)) NULL else scale(X, scale = FALSE) }

## ---------------------------------------------------------------- profiles per adjustment (genes and modules)
arr <- function(f) array(NA_real_, c(length(f), length(KOS), length(LINES)), dimnames = list(f, KOS, LINES))
PR <- lapply(setNames(names(ADJ), names(ADJ)), function(adj) { A <- arr(genes); Mo <- arr(CM); sig <- matrix(NA_real_, length(genes), length(LINES), dimnames = list(genes, LINES))
  for (cl in LINES) { s <- ids[meta[ids, "Cell_Line"] == cl]; g <- meta[s, "Gene_Control"]; X <- covX(s, adj); y <- normalise(cs[genes, s], "tmm")$y
    e <- ko_effects(y, g, X); A[, , cl] <- e$t; sig[, cl] <- matrixStats::rowMaxs(abs(e$t), na.rm = TRUE); Mo[, , cl] <- ko_effects(cm_sc[CM, s], g, X)$t }
  message("profiles: ", adj); list(A = A, M = Mo, sig = sig) })

## ---------------------------------------------------------------- A across lines
setmean <- function(A, sets) { sets <- lapply(sets, intersect, genes); sets <- sets[sapply(sets, length) >= 10]; out <- array(NA_real_, c(length(sets), dim(A)[2:3]), dimnames = c(list(names(sets)), dimnames(A)[2:3])); for (i in seq_along(sets)) out[i, , ] <- colMeans(A[sets[[i]], , , drop = FALSE]); out }
hall <- { suppressPackageStartupMessages(library(msigdbr)); h <- msigdbr(species = "Homo sapiens", category = "H"); split(h$gene_symbol, sub("HALLMARK_", "", h$gs_name)) }
scr <- intersect(unique(unlist(lapply(KOS, function(k) ss$genes[order(abs(ss$tcc[, k]), decreasing = TRUE)[1:100]]))), genes)
SPACES <- c(signal500 = "KO-signal genes", screen = "CROP-seq genes", Cm20_genes = "Cm modules", Hallmark = "Hallmark")
A_res <- list(); K_res <- list()
for (adj in names(ADJ)) { p <- PR[[adj]]
  for (sp in names(SPACES)) { CL <- switch(sp, signal500 = cl_from_array(p$A, "signal", p$sig), screen = cl_from_array(p$A, scr), Cm20_genes = cl_from_array(setmean(p$A, top50)), Hallmark = cl_from_array(setmean(p$A, hall)))
    a <- agreement(CL, perms); A_res[[length(A_res) + 1]] <- data.frame(adj = adj, space = sp, t(a$global)); K_res[[length(K_res) + 1]] <- cbind(adj = adj, space = sp, a$per_KO) } }
A_res <- do.call(rbind, A_res); K_res <- do.call(rbind, K_res)
write.csv(A_res, file.path(OUT, "spec83_across_lines.csv"), row.names = FALSE); write.csv(K_res, file.path(OUT, "spec83_across_lines_per_KO.csv"), row.names = FALSE)
options(width = 220); cat("\n=== A across lines: global permutation p (z) ===\n")
print(reshape(transform(A_res[, c("adj", "space")], v = sprintf("%.2g (z %.2f)", A_res$perm_p, A_res$perm_z)), idvar = "space", timevar = "adj", direction = "wide"), row.names = FALSE)
cat("\nper TF, KO-signal genes (perm p):\n"); x <- K_res[K_res$space == "signal500", ]; print(reshape(transform(x[, c("KO", "adj")], p = round(x$perm_p, 4)), idvar = "KO", timevar = "adj", direction = "wide"), row.names = FALSE)

## ---------------------------------------------------------------- B vs CROP-seq
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); S0 <- sc$score[, CM]; grp <- as.character(sc$group)
screen_eff <- function(adj) { Z <- switch(adj, none = NULL, Cm17 = cbind(S0[, "Cm17"]), cc = cbind(sc$cc), cc_Cm17 = cbind(sc$cc, S0[, "Cm17"]))
  S <- if (is.null(Z)) S0 else apply(S0, 2, function(y) resid(lm(y ~ Z)) + mean(y)); n <- apply(S[grp == "NTC", ], 2, median)
  sapply(KOS, function(k) SGN * (n - apply(S[grp == k, , drop = FALSE], 2, median))) }
B_res <- do.call(rbind, lapply(names(ADJ), function(adj) { E <- screen_eff(adj); Mb <- PR[[adj]]$M; keep <- CM[CM != "Cm17" | adj %in% c("none", "cc")]   # Cm17 itself is ~0 once adjusted
  do.call(rbind, lapply(c("uncentred", "KO-centred"), function(md) { f <- function(X) if (md == "uncentred") X else X - rowMeans(X, na.rm = TRUE)
    per <- do.call(rbind, lapply(LINES, function(cl) { x <- as.vector(f(E[keep, ])); y <- as.vector(f(Mb[keep, , cl])); ct <- cor.test(x, y, method = "spearman", exact = FALSE)
      data.frame(adj = adj, mode = md, line = cl, rho = unname(ct$estimate), p = ct$p.value, sign_agree = mean(sign(x) == sign(y))) }))
    x <- as.vector(f(E[keep, ])); y <- as.vector(f(apply(Mb[keep, , , drop = FALSE], c(1, 2), mean))); ct <- cor.test(x, y, method = "spearman", exact = FALSE)
    rbind(per, data.frame(adj = adj, mode = md, line = "mean of 4 lines", rho = unname(ct$estimate), p = ct$p.value, sign_agree = mean(sign(x) == sign(y)))) })) }))
write.csv(B_res, file.path(OUT, "spec83_vs_cropseq.csv"), row.names = FALSE)
cat("\n=== B bulk vs CROP-seq (KO x Cm module effects, Spearman rho) ===\n"); print(reshape(transform(B_res[, c("mode", "line", "adj")], v = sprintf("%.2f (p %.2g)", B_res$rho, B_res$p)), idvar = c("mode", "line"), timevar = "adj", direction = "wide"), row.names = FALSE)

## ---------------------------------------------------------------- figure
ACOL <- c(none = "grey60", Cm17 = "#e87ba4", cc = "#2a78d6", cc_Cm17 = "#4a3aa7")
save_fig("fig9E_Cm17_adjustment", 15, 5.6, function() {
  layout(matrix(1:3, nrow = 1), widths = c(1.1, 1.2, 1)); par(oma = c(0, 0, 3.4, 0.4), family = "sans"); ttl <- function(x) mtext(x, side = 3, adj = 0, line = 0.6, font = 2, cex = 0.8)
  par(mar = c(5.2, 4.2, 2.6, 0.6)); Z <- sapply(names(ADJ), function(a) A_res$perm_z[A_res$adj == a][match(names(SPACES), A_res$space[A_res$adj == a])]); P <- sapply(names(ADJ), function(a) A_res$perm_p[A_res$adj == a][match(names(SPACES), A_res$space[A_res$adj == a])])
  bp <- barplot(t(Z), beside = TRUE, col = ACOL, border = NA, las = 1, ylab = "global agreement across lines (permutation z)", names.arg = SPACES, cex.names = 0.75, cex.axis = 0.75, ylim = range(0, Z) * c(1.1, 1.2)); abline(h = qnorm(0.95), lty = 3, col = MUTED)
  text(bp, t(Z), ifelse(t(P) < 0.001, "***", ifelse(t(P) < 0.01, "**", ifelse(t(P) < 0.05, "*", ""))), pos = 3, cex = 0.8, xpd = NA); legend("topright", ADJ, fill = ACOL, border = NA, bty = "n", cex = 0.7); ttl("A  agreement across cell lines")
  par(mar = c(4.6, 4.2, 2.6, 0.6)); x <- K_res[K_res$space == "signal500", ]; o <- x$KO[x$adj == "none"][order(x$perm_p[x$adj == "none"])]
  Zk <- sapply(names(ADJ), function(a) x$perm_z[x$adj == a][match(o, x$KO[x$adj == a])]); Pk <- sapply(names(ADJ), function(a) x$perm_p[x$adj == a][match(o, x$KO[x$adj == a])])
  bp <- barplot(t(Zk), beside = TRUE, col = ACOL, border = NA, las = 1, ylab = "same-KO agreement (permutation z)", names.arg = rep("", length(o)), cex.axis = 0.75, ylim = range(Zk) * c(1.15, 1.25)); abline(h = c(0, qnorm(0.95)), lty = c(1, 3), col = c(INK, MUTED))
  text(bp, t(Zk), ifelse(t(Pk) < 0.01, "**", ifelse(t(Pk) < 0.05, "*", "")), pos = 3, cex = 0.75, xpd = NA); text(colMeans(bp), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), o, srt = 45, adj = 1, xpd = NA, cex = 0.78, font = 2, col = KO_COL[o]); ttl("B  per TF (KO-signal genes)")
  par(mar = c(5.2, 4.2, 2.6, 0.6)); d <- B_res[B_res$mode == "uncentred", ]; lv <- c(LINES, "mean of 4 lines"); R <- sapply(names(ADJ), function(a) d$rho[d$adj == a][match(lv, d$line[d$adj == a])]); Pp <- sapply(names(ADJ), function(a) d$p[d$adj == a][match(lv, d$line[d$adj == a])])
  bp <- barplot(t(R), beside = TRUE, col = ACOL, border = NA, las = 1, ylab = "Spearman rho, bulk vs CROP-seq (KO x Cm)", names.arg = sub("NCI-|-211H", "", sub("mean of 4 lines", "4-line mean", lv)), cex.names = 0.75, cex.axis = 0.75); abline(h = 0, col = INK)
  text(bp, pmax(t(R), 0), ifelse(t(Pp) < 0.05, "*", ""), pos = 3, cex = 0.9, xpd = NA); ttl("C  matching with CROP-seq (NTC - KO module effects)")
  dc <- B_res[B_res$mode == "KO-centred" & B_res$line == "mean of 4 lines", ]; mtext(sprintf("KO-centred, 4-line mean: %s", paste(sprintf("%s %.2f", ADJ[dc$adj], dc$rho), collapse = "; ")), side = 1, line = 3.6, cex = 0.5, col = INK2)
  mtext("Adjusting each cell line for the sarcomatoid score (scS = Cm17 module score) before computing KO profiles (S70)", outer = TRUE, side = 3, adj = 0, line = 1.8, font = 2, cex = 1)
  mtext("TMM logCPM ~ group + covariate(s), limma t of NTC - KO; across lines: KO-centred profiles, Spearman, KO-label permutation; vs CROP-seq: screen cells adjusted the same way; * p<0.05 ** p<0.01 *** p<0.001", outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.58, col = INK2)
})
