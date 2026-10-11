# Spec 124: KO-specific signal after removing the cell cycle from both data sets, without any NTC subtraction.
# CROP-seq: cells of the 8 shared KOs only (no NTC). Log-normalised expression regressed on S.Score + G2M.Score over those cells; KO profile = mean residual of the KO's cells.
# bulk    : per cell line, the KO samples only (no NTC; S73 labels: H2052_NTC_R3 = BPTF). TMM logCPM regressed on the cell-cycle score (S.Score + G2M.Score) over the KO samples of the line;
#           KO profile = mean residual of the KO's samples.
# both    : each gene centred on the mean of the 8 KO profiles of its data set (KO vs the other KOs; no control involved).
# 2 x 2   : screen as is / cell cycle regressed  x  bulk as is / cell cycle regressed, the CROP-seq compared with each of the 4 bulk lines.
# genes   : top 1000 / 2000 variable in both (screen: vst ranking of the KO cells; bulk: variance of the processed sample matrix of the line) and all shared genes.
# KO tiers per line (editing evidence, spec118 best-guide indel % in the KO samples of the line; guides detected in all, spec119): all 8 KOs | indel >= 5% | indel >= 10%
# metric  : Spearman KO x KO (screen rows, bulk columns) over the KOs of the tier; own minus other; exact permutation of the bulk KO labels; per KO own rho and rank among the 8 bulk KOs.
# needs: spec72 (S70 lists), spec118 (guide summary), agree_common.R::load_bulk_S91, spec_common.R
# outputs: spec124_summary.csv, spec124_per_KO.csv, spec124_indel_tiers.csv, fig41_cc_regressed_KO_profiles.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages(library(Seurat)); options(width = 230); fp <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))); SH <- function(l) sub("NCI-|-211H", "", l)
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); keep <- pr$genes; B <- load_bulk_S91(pr$samples, extra_kos = FALSE); m <- B$meta

## ---------------------------------------------------------------- editing tiers per KO x line
gsum <- read.csv(file.path(OUT, "spec118_guide_summary.csv"), stringsAsFactors = FALSE); gsum <- gsum[gsum$target %in% KOS & gsum$KO_pairs >= 10, ]
IND <- matrix(NA_real_, 8, 4, dimnames = list(KOS, LINES)); for (i in seq_len(nrow(gsum))) { v <- gsum$KO_indel_pct[i]; r <- gsum$target[i]; cl <- gsum$line[i]; if (is.na(IND[r, cl]) || v > IND[r, cl]) IND[r, cl] <- v }
write.csv(IND, file.path(OUT, "spec124_indel_tiers.csv")); cat("best-guide indel % per KO x line (NA = too few reads at the sites):\n"); print(round(IND, 1))
TIER <- list(all = function(cl) KOS, `indel>=5` = function(cl) KOS[!is.na(IND[, cl]) & IND[, cl] >= 5], `indel>=10` = function(cl) KOS[!is.na(IND[, cl]) & IND[, cl] >= 10])

## ---------------------------------------------------------------- screen profiles (KO cells only)
srt <- load_screen(); grp <- as.character(srt$merged_call); ck <- which(grp %in% KOS); Y <- GetAssayData(srt, layer = "data")[, ck]; g <- grp[ck]; md <- srt@meta.data[ck, ]
Y <- Y[Matrix::rowSums(Y > 0) >= 0.01 * ncol(Y), ]; Gs <- intersect(keep, rownames(Y)); Y <- Y[Gs, ]; cat("\nscreen: KO cells", ncol(Y), "| genes shared with bulk", length(Gs), "\n"); print(table(g)[KOS])
Xc <- scale(cbind(S = md$S.Score, G2M = md$G2M.Score), scale = FALSE); Bc <- as.matrix(Y %*% Xc) %*% solve(crossprod(Xc))                       # genes x 2 regression coefficients (intercept absorbed by centring)
gm <- sapply(KOS, function(k) Matrix::rowMeans(Y[, g == k, drop = FALSE])); gx <- sapply(KOS, function(k) colMeans(Xc[g == k, , drop = FALSE]))   # group means of expression and of the centred scores
SCR <- list(raw = gm - rowMeans(gm), cc = { r <- gm - Bc %*% gx; r - rowMeans(r) }); cat("cell-cycle score of the KO groups in the screen (S + G2M, centred):", paste(KOS, round(colSums(gx), 2), collapse = "  "), "\n")
hv <- intersect(VariableFeatures(FindVariableFeatures(srt[, ck], selection.method = "vst", nfeatures = nrow(srt), verbose = FALSE)), Gs)

## ---------------------------------------------------------------- bulk profiles per line (KO samples only)
BEST_REP <- Sys.getenv("BEST_REP", unset = "0") == "1"; TAG <- if (BEST_REP) "_bestrep" else ""      # BEST_REP=1: KO profile from the single replicate with the highest indel % (spec120 evidence table); the regression still uses all KO samples
ev <- read.csv(file.path(OUT, "spec120_sample_cut_evidence.csv"), stringsAsFactors = FALSE); ev_pct <- setNames(ev$best_pct, ev$sample)
BLK <- list(); BV <- list()
for (cl in LINES) { s <- B$ids[m[B$ids, "Cell_Line"] == cl & m[B$ids, "Gene_Control"] %in% KOS]; gl <- m[s, "Gene_Control"]; y <- normalise(B$cs[keep, s], "tmm")$y[Gs, ]; x <- B$cc[s] - mean(B$cc[s]); yr <- y - (y %*% x / sum(x^2)) %*% t(x)
  pick <- function(k) { w <- which(gl == k); if (BEST_REP && length(w) > 1) { e <- ev_pct[s[w]]; e[is.na(e)] <- -1; w <- w[which.max(e)] }; w }; prof <- function(M) { p <- sapply(KOS, function(k) rowMeans(M[, pick(k), drop = FALSE])); p - rowMeans(p) }; BLK[[cl]] <- list(raw = prof(y), cc = prof(yr)); BV[[cl]] <- list(raw = matrixStats::rowVars(y), cc = matrixStats::rowVars(yr))
  cat(sprintf("%s: %d KO samples; cell-cycle score explains %.0f%% of the variance of the top-2000 variable genes (median R2)\n", cl, length(s), 100 * median((cor(t(y[order(-BV[[cl]]$raw)[1:2000], ]), x))^2))) }

## ---------------------------------------------------------------- comparisons
perm_all <- function(n) { p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(seq_len(n))) }; PM <- lapply(1:8, perm_all)
RES <- list(); PK <- list(); MODES <- expand.grid(screen = c("raw", "cc"), bulk = c("raw", "cc"), stringsAsFactors = FALSE)
for (cl in LINES) for (i in seq_len(nrow(MODES))) { sm <- MODES$screen[i]; bm <- MODES$bulk[i]; vb <- Gs[order(-BV[[cl]][[bm]])]
  for (sp in c("int1000", "int2000", "all")) { f <- switch(sp, int1000 = intersect(hv[1:1000], vb[1:1000]), int2000 = intersect(hv[1:2000], vb[1:2000]), all = Gs); C <- cor(SCR[[sm]][f, ], BLK[[cl]][[bm]][f, ], method = "spearman")
    PK[[length(PK) + 1]] <- data.frame(line = cl, screen = sm, bulk = bm, space = sp, KO = KOS, own = diag(C), rank = sapply(1:8, function(k) rank(-C[k, ])[k]), indel = IND[, cl])
    for (tn in names(TIER)) { ks <- TIER[[tn]](cl); if (length(ks) < 3) next; Ct <- C[ks, ks]; dm <- mean(diag(Ct)); P <- PM[[length(ks)]]; nul <- apply(P, 1, function(p) mean(Ct[cbind(seq_along(ks), p)]))
      RES[[length(RES) + 1]] <- data.frame(line = cl, screen = sm, bulk = bm, space = sp, n_genes = length(f), tier = tn, n_KO = length(ks), KOs = paste(ks, collapse = " "), diag = dm, off = mean(Ct[row(Ct) != col(Ct)]), perm_p = mean(nul >= dm - 1e-12), top1 = sum(apply(Ct, 1, which.max) == seq_along(ks))) } } }
RES <- do.call(rbind, RES); RES$dmo <- RES$diag - RES$off; RES$mode <- paste0("screen ", RES$screen, " / bulk ", RES$bulk); PK <- do.call(rbind, PK); write.csv(RES, file.path(OUT, sprintf("spec124_summary%s.csv", TAG)), row.names = FALSE); write.csv(PK, file.path(OUT, sprintf("spec124_per_KO%s.csv", TAG)), row.names = FALSE)
for (tn in names(TIER)) { x <- RES[RES$space == "int2000" & RES$tier == tn, ]; cat(sprintf("\n=== own - other Spearman (perm p; KOs), genes in the top 2000 of both, tier: %s ===\n", tn)); x$v <- sprintf("%+.3f (p %s; %d)", x$dmo, fp(x$perm_p), x$n_KO); print(reshape(x[, c("line", "mode", "v")], idvar = "line", timevar = "mode", direction = "wide"), row.names = FALSE) }
x <- RES[RES$tier == "all" & RES$screen == RES$bulk, ]; cat("\n-- all 8 KOs, raw/raw vs cc/cc, by gene set\n"); x$v <- sprintf("%+.3f (p %s)", x$dmo, fp(x$perm_p)); print(reshape(x[, c("line", "space", "mode", "v")], idvar = c("line", "space"), timevar = "mode", direction = "wide"), row.names = FALSE)
for (md in list(c("raw", "raw"), c("cc", "cc"))) { x <- PK[PK$space == "int2000" & PK$screen == md[1] & PK$bulk == md[2], ]; cat(sprintf("\n-- per KO own rho (rank among 8), screen %s / bulk %s, top 2000 in both\n", md[1], md[2])); x$v <- sprintf("%+.2f (%s)", x$own, format(x$rank)); print(reshape(x[, c("KO", "line", "v")], idvar = "KO", timevar = "line", direction = "wide"), row.names = FALSE) }
x <- PK[PK$space == "int2000" & PK$screen == "cc" & PK$bulk == "cc" & !is.na(PK$indel), ]; r <- cor.test(x$own, x$indel, method = "spearman", exact = FALSE); cat(sprintf("\nown rho (cc/cc) vs indel %% over %d KO x line: Spearman %.2f, p %.2g | mean own rho: indel >= 10%% %.3f, < 10%% %.3f\n", nrow(x), r$estimate, r$p.value, mean(x$own[x$indel >= 10]), mean(x$own[x$indel < 10])))
x0 <- PK[PK$space == "int2000" & PK$screen == "raw" & PK$bulk == "raw", ]; x1 <- PK[PK$space == "int2000" & PK$screen == "cc" & PK$bulk == "cc", ]; cat("KO x line with the own KO ranked first: raw/raw", sum(x0$rank == 1), "of 32 | cc/cc", sum(x1$rank == 1), "of 32 | ranked 1-2: raw", sum(x0$rank <= 2), "| cc", sum(x1$rank <= 2), "(expected by chance: 4 and 8)\n")

## ---------------------------------------------------------------- figure
MC <- c(`screen raw / bulk raw` = "grey55", `screen cc / bulk raw` = "#8fb8ea", `screen raw / bulk cc` = "#f0a58f", `screen cc / bulk cc` = "#2a78d6")
save_fig(paste0("fig41_cc_regressed_KO_profiles", TAG), 16, 10.4, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4.2, 0.4), family = "sans"); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.8); mtext(b, side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2) }
  for (tn in c("all", "indel>=10")) { x <- RES[RES$space == "int2000" & RES$tier == tn, ]; M <- sapply(LINES, function(cl) x$dmo[x$line == cl][match(names(MC), x$mode[x$line == cl])]); P <- sapply(LINES, function(cl) x$perm_p[x$line == cl][match(names(MC), x$mode[x$line == cl])]); nk <- sapply(LINES, function(cl) x$n_KO[x$line == cl][1])
    par(mar = c(4.6, 4.6, 3.6, 1)); bp <- barplot(M, beside = TRUE, col = MC, border = NA, las = 1, names.arg = sprintf("%s\n(%s KOs)", SH(LINES), ifelse(is.na(nk), "<3", nk)), cex.names = 0.8, ylab = "own-KO minus other-KO Spearman, CROP-seq vs bulk", cex.axis = 0.8, ylim = range(c(0, M), na.rm = TRUE) + c(-0.03, 0.07)); abline(h = 0, col = INK)
    text(bp, ifelse(is.na(M), 0, pmax(M, 0)) + 0.004, ifelse(is.na(P), "", fp(P)), cex = 0.55, srt = 90, adj = c(0, 0.5), font = ifelse(!is.na(P) & P < 0.05, 2, 1)); if (tn == "all") legend("topright", names(MC), fill = MC, border = NA, bty = "n", cex = 0.62)
    ttl(if (tn == "all") "A  all 8 KOs" else "B  only KOs with >= 10% indels in the line", "genes in the top 2000 variable of both; numbers = exact permutation p; cc = cell cycle regressed") }
  x <- RES[RES$tier == "all" & RES$screen == "cc" & RES$bulk == "cc", ]; M <- sapply(LINES, function(cl) x$dmo[x$line == cl][match(c("int1000", "int2000", "all"), x$space[x$line == cl])]); P <- sapply(LINES, function(cl) x$perm_p[x$line == cl][match(c("int1000", "int2000", "all"), x$space[x$line == cl])])
  par(mar = c(4.6, 4.6, 3.6, 1)); gc <- c("#1baf7a", "#2a78d6", "#4a3aa7"); bp <- barplot(M, beside = TRUE, col = gc, border = NA, las = 1, names.arg = SH(LINES), cex.names = 0.85, ylab = "own-KO minus other-KO Spearman", cex.axis = 0.8, ylim = range(c(0, M)) + c(-0.03, 0.07)); abline(h = 0, col = INK)
  text(bp, pmax(M, 0) + 0.004, fp(P), cex = 0.55, srt = 90, adj = c(0, 0.5), font = ifelse(P < 0.05, 2, 1)); legend("topright", c("top 1000 in both", "top 2000 in both", "all shared genes"), fill = gc, border = NA, bty = "n", cex = 0.65); ttl("C  cell cycle regressed in both: gene sets", "all 8 KOs")
  pal <- colorRampPalette(c("#2a78d6", "#f7f6f2", "#eb6834"))(101)
  for (md in list(c("raw", "raw"), c("cc", "cc"))) { x <- PK[PK$space == "int2000" & PK$screen == md[1] & PK$bulk == md[2], ]; M <- sapply(LINES, function(cl) x$own[x$line == cl][match(KOS, x$KO[x$line == cl])]); R <- sapply(LINES, function(cl) x$rank[x$line == cl][match(KOS, x$KO[x$line == cl])]); rownames(M) <- KOS
    par(mar = c(4.6, 6, 3.6, 4.6)); zl <- 0.6; image(1:4, 1:8, t(pmax(pmin(M, zl), -zl)[8:1, ]), col = pal, zlim = c(-zl, zl), axes = FALSE, xlab = "", ylab = ""); text(rep(1:4, each = 8), rep(8:1, 4), sprintf("%.2f\n#%s", as.vector(M), format(as.vector(R))), cex = 0.62, font = ifelse(as.vector(R) == 1, 2, 1))
    lo <- which(is.na(IND) | IND < 10, arr.ind = TRUE); points(lo[, 2] + 0.38, 9 - lo[, 1] + 0.36, pch = 4, cex = 0.7, col = "grey30"); axis(1, at = 1:4, labels = SH(LINES), tick = FALSE, cex.axis = 0.85); text(0.4, 8:1, KOS, adj = 1, xpd = NA, cex = 0.82, font = 2, col = KO_COL[KOS])
    yy <- seq(2, 7, length.out = 102); rect(4.75, yy[-102], 4.95, yy[-1], col = pal, border = NA, xpd = NA); text(5.02, c(2, 4.5, 7), sprintf("%.1f", c(-zl, 0, zl)), adj = 0, xpd = NA, cex = 0.6)
    ttl(sprintf("%s  per KO: own-KO Spearman, %s", if (md[1] == "raw") "D" else "E", if (md[1] == "raw") "as is" else "cell cycle regressed in both"), "# = rank of the own KO among the 8 bulk KOs of the line; x = KO with < 10% indels (or not assessable) in that line") }
  x <- PK[PK$space == "int2000" & PK$screen == "cc" & PK$bulk == "cc" & !is.na(PK$indel), ]; par(mar = c(4.6, 4.6, 3.6, 1)); plot(x$indel, x$own, pch = c(21, 22, 24, 23)[match(x$line, LINES)], bg = KO_COL[x$KO], col = INK2, cex = 1.6, bty = "n", las = 1, xlab = "best-guide indel % in the KO samples of the line", ylab = "own-KO Spearman, cell cycle regressed in both", cex.axis = 0.8)
  abline(h = 0, lty = 2, col = GRID); abline(v = 10, lty = 3, col = MUTED); abline(lm(own ~ indel, x), lty = 3, col = INK2); r <- cor.test(x$own, x$indel, method = "spearman", exact = FALSE); legend("topleft", SH(LINES), pch = c(21, 22, 24, 23), pt.bg = "grey70", bty = "n", cex = 0.65); legend("bottomright", KOS, pch = 21, pt.bg = KO_COL[KOS], bty = "n", cex = 0.55, ncol = 2)
  ttl("F  does better editing give a better match?", sprintf("%d KO x line: Spearman %.2f, p %s", nrow(x), r$estimate, fp(r$p.value)))
  mtext("KO-specific signal, CROP-seq vs each bulk line, with the cell cycle regressed out and no NTC subtraction", outer = TRUE, side = 3, adj = 0, line = 2.6, font = 2, cex = 1)
  mtext("KO cells / KO samples only; profiles = KO mean minus the mean of the 8 KOs of the data set; screen: regression on S.Score + G2M.Score over the cells; bulk: on the cell-cycle score over the KO samples of the line", outer = TRUE, side = 3, adj = 0, line = 1.4, cex = 0.65, col = INK2)
})
