# Spec 126: the most reproducible TF KOs -- agreement with the CROP-seq screen, agreement across the bulk lines, and which Cm modules they move.
# All contrasts are KO vs the other KOs (no NTC): value = mean of the KOs minus the KO, so > 0 = lower in that KO (same direction as NTC - KO).
#   A  CROP-seq vs each bulk line, per TF: own-KO Spearman of the gene profiles (spec124: KO cells / KO samples, KO-centred, genes in the top 2000 variable of both), # = rank of the own KO among the 8 bulk KOs;
#      a dot marks the cells that are still ranked 1-2 after regressing the cell cycle out of both data sets
#   B  bulk across the 4 lines, per TF (S91, 11 KOs): Spearman of the same KO between two lines vs other KOs (limma t, KO-centred, top-500 KO-signal genes), KO-label permutation p; sorted
#   C  the two axes together: x = mean same-KO Spearman across lines, y = mean own-KO Spearman with the CROP-seq over the 4 lines
#   D  Cm modules of the TFs that do best on both axes: module effect in the CROP-seq and in each bulk line (AddModuleScore top-50 genes, z within data set; mean of the KOs minus the KO)
#   E  per TF and line: Spearman over the 20 modules between the CROP-seq and the bulk module effects
# needs: spec124 (per-KO table), spec72 (S70 lists), spec75 (screen module scores), agree_common.R::load_bulk_S91
# outputs: spec126_TF_summary.csv, spec126_Cm_effects.csv, spec126_Cm_match.csv, fig43_reproducible_TFs_and_Cm.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R")); suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3", width = 220)
N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000")); fp <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))); SH <- function(l) sub("NCI-|-211H", "", l); kcol <- function(k) ifelse(k %in% names(KO_COL), KO_COL[k], "grey30")
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); keep <- pr$genes; top50 <- cm_top50(); CM <- names(top50)

## ---------------------------------------------------------------- A  CROP-seq vs bulk per TF (spec124)
pk <- read.csv(file.path(OUT, "spec124_per_KO.csv"), stringsAsFactors = FALSE); g2 <- function(sm, bm, col) { x <- pk[pk$space == "int2000" & pk$screen == sm & pk$bulk == bm, ]; sapply(LINES, function(cl) x[[col]][x$line == cl][match(KOS, x$KO[x$line == cl])]) }
A_own <- g2("raw", "raw", "own"); A_rank <- g2("raw", "raw", "rank"); A_rank_cc <- g2("cc", "cc", "rank"); A_own_cc <- g2("cc", "cc", "own"); rownames(A_own) <- rownames(A_rank) <- rownames(A_rank_cc) <- rownames(A_own_cc) <- KOS

## ---------------------------------------------------------------- B  across bulk lines (S91, 11 KOs)
B <- load_bulk_S91(pr$samples); TT <- array(NA_real_, c(length(keep), length(KOS11), 4), dimnames = list(keep, KOS11, LINES))
for (cl in LINES) { s <- B$ids[B$meta[B$ids, "Cell_Line"] == cl]; TT[, , cl] <- ko_effects(normalise(B$cs[keep, s], "tmm")$y, B$meta[s, "Gene_Control"], NULL, kos = KOS11)$t }
sig <- apply(abs(TT), c(1, 3), max, na.rm = TRUE); ag <- agreement(cl_from_array(TT, "signal", sig), make_perms(N_PERM, K = length(KOS11))); kb <- ag$per_KO[order(ag$per_KO$perm_p), ]
cat("across lines (S91): global p", fp(ag$global[["perm_p"]]), "\n"); print(kb[, c("KO", "n_pairs", "mean_intra", "median_cross", "perm_p", "perm_fdr")], digits = 3, row.names = FALSE)

## ---------------------------------------------------------------- C  summary per TF
SUM <- data.frame(KO = KOS, crop_mean_own = rowMeans(A_own), crop_lines_rank12 = rowSums(A_rank <= 2), crop_lines_rank1 = rowSums(A_rank == 1), crop_mean_own_cc = rowMeans(A_own_cc), crop_lines_rank12_cc = rowSums(A_rank_cc <= 2),
                  lines_mean_same = ag$per_KO$mean_intra[match(KOS, ag$per_KO$KO)], lines_p = ag$per_KO$perm_p[match(KOS, ag$per_KO$KO)], lines_fdr = ag$per_KO$perm_fdr[match(KOS, ag$per_KO$KO)])
SUM <- SUM[order(-(rank(SUM$crop_mean_own) + rank(SUM$lines_mean_same))), ]; write.csv(SUM, file.path(OUT, "spec126_TF_summary.csv"), row.names = FALSE); cat("\nTF summary (sorted by the two axes together):\n"); print(SUM, digits = 2, row.names = FALSE)
TOP <- head(SUM$KO, 3)

## ---------------------------------------------------------------- D / E  Cm modules: KO vs the other KOs, screen and each line
B8 <- load_bulk_S91(pr$samples, extra_kos = FALSE); ids <- B8$ids; y <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(B8$cs[, ids])), log = TRUE, prior.count = 2); y <- y[rowSums(B8$cs[, ids]) > 0, ]; y <- y - min(y)
sb <- SetAssayData(CreateSeuratObject(counts = B8$cs[rownames(y), ids]), slot = "data", new.data = as.matrix(y)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1)); SCO <- t(as.matrix(sb@meta.data[ids, paste0("Cm_", 1:20)])); rownames(SCO) <- CM
eff <- function(z, grp) { gm <- sapply(KOS, function(k) rowMeans(z[, grp == k, drop = FALSE])); rowMeans(gm) - gm }                                  # modules x KO: mean of the KOs minus the KO
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); gs <- as.character(sc$group); ck <- gs %in% KOS; EF <- list(`CROP-seq` = eff(t(scale(sc$score[ck, CM])), gs[ck]))
for (cl in LINES) { s <- ids[B8$meta[ids, "Cell_Line"] == cl & B8$meta[ids, "Gene_Control"] %in% KOS]; EF[[cl]] <- eff(t(scale(t(SCO[, s]))), B8$meta[s, "Gene_Control"]) }
CE <- do.call(rbind, lapply(names(EF), function(n) data.frame(source = n, Cm = rep(CM, 8), KO = rep(KOS, each = 20), effect = as.vector(EF[[n]])))); write.csv(CE, file.path(OUT, "spec126_Cm_effects.csv"), row.names = FALSE)
MM <- sapply(LINES, function(cl) sapply(KOS, function(k) cor(EF[["CROP-seq"]][, k], EF[[cl]][, k], method = "spearman"))); MP <- sapply(LINES, function(cl) sapply(KOS, function(k) suppressWarnings(cor.test(EF[["CROP-seq"]][, k], EF[[cl]][, k], method = "spearman", exact = FALSE)$p.value)))
write.csv(data.frame(KO = KOS, MM, p = MP, check.names = FALSE), file.path(OUT, "spec126_Cm_match.csv"), row.names = FALSE); cat("\nmodule-level Spearman, CROP-seq vs bulk line (20 Cm), per TF:\n"); print(round(MM, 2)); cat("p:\n"); print(signif(MP, 2))
for (k in TOP) { S <- EF[["CROP-seq"]][, k]; Bm <- sapply(LINES, function(cl) EF[[cl]][, k]); agree <- rowSums(sign(Bm) == sign(S)); o <- order(-agree, -abs(S)); cat(sprintf("\n%s: modules with the screen direction in >= 3 of 4 bulk lines (screen effect; bulk H2052 MSTO H2452 H28):\n", k))
  for (i in o[agree[o] >= 3]) cat(sprintf("  %-5s %s in KO  screen %+.2f | bulk %s | %d/4 lines\n", CM[i], ifelse(S[i] > 0, "lower ", "higher"), S[i], paste(sprintf("%+.2f", Bm[i, ]), collapse = " "), agree[i])) }

## ---------------------------------------------------------------- figure
pal <- colorRampPalette(c("#2a78d6", "#f7f6f2", "#eb6834"))(101); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.8); mtext(b, side = 3, adj = 0, line = 0.5, cex = 0.56, col = INK2) }
save_fig("fig43_reproducible_TFs_and_Cm", 18, 11.4, function() {
  layout(rbind(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3), c(4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7)), heights = c(1, 1.15)); par(oma = c(0, 0, 4, 0.4), family = "sans")
  ## A
  o <- SUM$KO; M <- A_own[o, ]; R <- A_rank[o, ]; Rc <- A_rank_cc[o, ]; zl <- 0.6; par(mar = c(4, 6, 3.6, 4.4)); image(1:4, 1:8, t(pmax(pmin(M, zl), -zl)[8:1, ]), col = pal, zlim = c(-zl, zl), axes = FALSE, xlab = "", ylab = "")
  text(rep(1:4, each = 8), rep(8:1, 4), sprintf("%+.2f\n#%s", as.vector(M), format(as.vector(R))), cex = 0.68, font = ifelse(as.vector(R) <= 2, 2, 1)); w <- which(Rc <= 2, arr.ind = TRUE); points(w[, 2] + 0.36, 9 - w[, 1] + 0.3, pch = 16, cex = 0.9, col = "black")
  axis(1, at = 1:4, labels = SH(LINES), tick = FALSE, cex.axis = 0.9); text(0.4, 8:1, o, adj = 1, xpd = NA, cex = 0.9, font = 2, col = kcol(o)); yy <- seq(2, 7, length.out = 102); rect(4.75, yy[-102], 4.95, yy[-1], col = pal, border = NA, xpd = NA); text(5.02, c(2, 4.5, 7), sprintf("%.1f", c(-zl, 0, zl)), adj = 0, xpd = NA, cex = 0.62)
  ttl("A  match with the CROP-seq, per TF and bulk line", "own-KO Spearman of the gene profiles (KO vs other KOs); # = rank among the 8 bulk KOs; dot = still ranked 1-2 with the cell cycle regressed out")
  ## B
  ob <- kb$KO; n <- length(ob); xr <- range(c(unlist(ag$cross), ag$intra), na.rm = TRUE); par(mar = c(4, 5.4, 3.6, 6)); plot(NA, xlim = xr, ylim = c(n + 0.6, 0.4), yaxt = "n", bty = "n", las = 1, ylab = "", xlab = "Spearman between two bulk cell lines", cex.axis = 0.8); abline(v = 0, lty = 2, col = INK)
  for (i in seq_len(n)) { j <- match(ob[i], ag$per_KO$KO); v <- ag$intra[ob[i], ]; v <- v[!is.na(v)]; boxplot(ag$cross[[j]], at = i + 0.2, add = TRUE, horizontal = TRUE, boxwex = 0.32, col = "#d5d4cf", border = INK2, outline = FALSE, axes = FALSE)
    if (length(v) >= 3) boxplot(v, at = i - 0.2, add = TRUE, horizontal = TRUE, boxwex = 0.32, col = adjustcolor(kcol(ob[i]), 0.7), border = INK2, outline = FALSE, axes = FALSE); points(v, i - 0.2 + seq(-0.08, 0.08, length.out = length(v)), pch = 16, cex = 0.7, col = INK2)
    text(xr[2] + 0.03 * diff(xr), i, sprintf("p = %s%s", fp(kb$perm_p[i]), ifelse(kb$perm_fdr[i] < 0.05, " *", "")), adj = 0, xpd = NA, cex = 0.72, font = ifelse(kb$perm_p[i] < 0.05, 2, 1)) }
  text(xr[1] - 0.04 * diff(xr), seq_len(n), ob, adj = 1, xpd = NA, font = 2, cex = 0.85, col = kcol(ob)); ttl("B  reproducibility across the 4 bulk lines, per TF", sprintf("colour = same KO in two lines, grey = other KOs; KO-label permutation p, * FDR < 0.05; all KOs p = %s", fp(ag$global[["perm_p"]])))
  ## C
  par(mar = c(4.4, 4.8, 3.6, 1)); plot(SUM$lines_mean_same, SUM$crop_mean_own, pch = 21, bg = kcol(SUM$KO), col = INK2, cex = 2.4, bty = "n", las = 1, xlab = "across bulk lines: mean same-KO Spearman", ylab = "with the CROP-seq: mean own-KO Spearman over the 4 lines", cex.axis = 0.8, xlim = range(SUM$lines_mean_same) + c(-0.02, 0.04), ylim = range(SUM$crop_mean_own) + c(-0.03, 0.04))
  abline(h = 0, v = 0, lty = 2, col = GRID); text(SUM$lines_mean_same, SUM$crop_mean_own, SUM$KO, pos = 4, cex = 0.75, font = 2, col = kcol(SUM$KO), xpd = NA); ttl("C  the two kinds of reproducibility together", "top right = reproducible across bulk lines and matching the screen")
  ## D  module heatmaps of the top TFs
  for (k in TOP) { Mk <- cbind(`CROP-seq` = EF[["CROP-seq"]][, k], sapply(LINES, function(cl) EF[[cl]][, k])); colnames(Mk) <- c("CROP-seq", SH(LINES)); Z <- sweep(Mk, 2, apply(abs(Mk), 2, max), "/"); agree <- rowSums(sign(Mk[, -1]) == sign(Mk[, 1])); oo <- order(-Mk[, 1])
    par(mar = c(4, 4.2, 3.6, 3.2)); image(1:5, 1:20, t(Z[rev(oo), ]), col = pal, zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = ""); abline(v = 1.5, lwd = 2, col = "white"); axis(1, at = 1:5, labels = colnames(Mk), tick = FALSE, cex.axis = 0.72, las = 2); axis(2, at = 20:1, labels = CM[oo], tick = FALSE, las = 1, cex.axis = 0.66)
    text(5.75, 20:1, ifelse(agree[oo] >= 3, sprintf("%d/4", agree[oo]), ""), cex = 0.6, xpd = NA, font = 2); mtext(sprintf("D  %s: Cm modules", k), side = 3, adj = 0, line = 2, font = 2, cex = 0.8, col = kcol(k))
    mtext(sprintf("orange = lower in the KO; n/4 = bulk lines agreeing\nrho vs screen: %s", paste(sprintf("%s %.2f", SH(LINES), MM[k, ]), collapse = ", ")), side = 3, adj = 0, line = 0.2, cex = 0.5, col = INK2) }
  ## E
  o <- SUM$KO; M <- MM[o, ]; P <- MP[o, ]; par(mar = c(4, 6, 3.6, 4.4)); zl <- 0.8; image(1:4, 1:8, t(pmax(pmin(M, zl), -zl)[8:1, ]), col = pal, zlim = c(-zl, zl), axes = FALSE, xlab = "", ylab = ""); text(rep(1:4, each = 8), rep(8:1, 4), sprintf("%+.2f%s", as.vector(M), ifelse(as.vector(P) < 0.05, "*", "")), cex = 0.72, font = ifelse(as.vector(P) < 0.05, 2, 1))
  axis(1, at = 1:4, labels = SH(LINES), tick = FALSE, cex.axis = 0.9); text(0.4, 8:1, o, adj = 1, xpd = NA, cex = 0.9, font = 2, col = kcol(o)); ttl("E  module-level match with the CROP-seq", "Spearman over the 20 Cm modules; * p < 0.05")
  mtext("Most reproducible TF KOs: agreement with the CROP-seq, across the bulk cell lines, and the Cm modules they move", outer = TRUE, side = 3, adj = 0, line = 2.4, font = 2, cex = 1.05)
  mtext("all effects are KO vs the other KOs (no NTC): mean of the KOs minus the KO, so > 0 / orange = lower in that KO; corrected sample set (H2052_NTC_R3 as BPTF)", outer = TRUE, side = 3, adj = 0, line = 1.2, cex = 0.68, col = INK2)
})
