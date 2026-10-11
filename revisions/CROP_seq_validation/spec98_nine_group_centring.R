# Spec 98: correlations after subtracting, per gene, the mean of the 9 groups (NTC + 8 KOs) of each data set -- the centring of fig5E / fig5G step 3 -- instead of NTC subtraction.
#   bulk  : per line, group mean of TMM logCPM (S70) for NTC and the 8 KOs; profile of a KO = its group mean minus the mean of the 9 group means
#   screen: pseudobulk logCPM (counts summed per guide group) for NTC and the 8 KOs; same centring over its 9 groups
#   "c9" = centred over the 9 groups (requested) | "z9" = c9 divided by the gene's SD over the 9 groups | "c8" = centred over the 8 KOs only (= NTC - KO centred; the earlier analyses), sign flipped
#   to KO-minus-mean so all three are comparable (correlations are unaffected by flipping both sides)
# A  NCI-H2052, CROP-seq vs bulk: 8 x 8 KO Spearman, own vs other, exact permutation p of the bulk KO labels (8!), own-KO rank; gene sets: all shared | top 1000/2000/3000 H2052-variable | screen HVGs | 20 Cm modules (mean AddModuleScore per group)
# B  bulk across the 4 lines: agreement() (Spearman, KO-label permutation) on top-500 KO-signal genes (from the 2 lines not compared) | top 2000 variable | all genes | 20 Cm modules
# needs: spec72 (S70), spec88 (screen HVG list), spec_common.R, spec1
# outputs: spec98_H2052_crop_vs_bulk.csv, spec98_H2052_per_KO.csv, spec98_across_lines.csv, spec98_across_lines_per_KO.csv, fig11_nine_group_centring.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages(library(Seurat)); N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000")); permsL <- make_perms(N_PERM)
CL <- "NCI-H2052"; GR <- c("NTC", KOS); fab <- load_full_audited_bulk(); meta <- fab$meta; pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); ids <- pr$samples
srt <- load_screen(); grp <- as.character(srt$merged_call); X <- GetAssayData(srt, layer = "counts")
PB <- sapply(GR, function(k) Matrix::rowSums(X[, grp == k, drop = FALSE])); lcS <- edgeR::cpm(edgeR::DGEList(PB), log = TRUE, prior.count = 2); expS <- rownames(lcS)[rowSums(edgeR::cpm(edgeR::DGEList(PB)) >= 10) >= 1]
G <- intersect(pr$genes, expS); cat("genes:", length(G), "\n")
bulkG <- lapply(setNames(LINES, LINES), function(cl) { s <- ids[meta[ids, "Cell_Line"] == cl]; y <- normalise(fab$cs[pr$genes, s], "tmm")$y; g <- meta[s, "Gene_Control"]; sapply(GR, function(k) rowMeans(y[, g == k, drop = FALSE])) })   # genes x 9 groups
norm9 <- function(M, mode) { M <- M[, GR]; switch(mode, c9 = (M - rowMeans(M))[, KOS], z9 = { Z <- (M - rowMeans(M)) / matrixStats::rowSds(M); Z[!is.finite(Z)] <- 0; Z[, KOS] }, c8 = M[, KOS] - rowMeans(M[, KOS])) }
## 20 Cm modules (AddModuleScore, top-50 genes): mean score per group, appended as extra rows so that the same centring applies
CM <- paste0("Cm", 1:20); scS <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); lcS <- rbind(lcS, sapply(GR, function(k) colMeans(scS$score[as.character(scS$group) == k, CM, drop = FALSE])))
for (cl in LINES) { s <- ids[meta[ids, "Cell_Line"] == cl]; g <- meta[s, "Gene_Control"]; bulkG[[cl]] <- rbind(bulkG[[cl]], sapply(GR, function(k) rowMeans(pr$cm_scores$tmm[CM, s[g == k], drop = FALSE]))) }
MODES <- c(c9 = "minus mean of the 9 groups", z9 = "z over the 9 groups", c8 = "minus mean of the 8 KOs (earlier)"); MC <- c(c9 = "#eb6834", z9 = "#eda100", c8 = "grey55")
perms8 <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(1:8)) })

## ---------------------------------------------------------------- A  H2052 screen vs bulk
Lh <- read.csv(file.path(OUT, "spec88_screen_pseudobulk_logCPM.csv"), row.names = 1, check.names = FALSE); vH <- G[order(pr$gene_var[G, CL], decreasing = TRUE)]
SPA <- list(all = G, top1000 = vH[1:1000], top2000 = vH[1:2000], top3000 = vH[1:3000], screenHVG = intersect(rownames(Lh), G), Cm20 = CM); SPAL <- c(all = "all shared genes", top1000 = "top 1000 variable", top2000 = "top 2000 variable", top3000 = "top 3000 variable", screenHVG = "screen HVGs", Cm20 = "20 Cm modules")
rA <- list(); kA <- list()
for (sp in names(SPA)) for (md in names(MODES)) { f <- SPA[[sp]]; S <- norm9(lcS[f, ], md); B <- norm9(bulkG[[CL]][f, ], md); C <- cor(S, B, method = "spearman"); dm <- mean(diag(C)); nul <- apply(perms8, 1, function(p) mean(C[cbind(1:8, p)]))
  rA[[length(rA) + 1]] <- data.frame(space = sp, n_genes = length(f), mode = md, diag = dm, off = mean(C[row(C) != col(C)]), perm_p = mean(nul >= dm - 1e-12), top1 = sum(apply(C, 1, which.max) == 1:8))
  kA[[length(kA) + 1]] <- data.frame(space = sp, mode = md, KO = KOS, own = diag(C), median_other = sapply(1:8, function(i) median(C[i, -i])), rank = sapply(1:8, function(i) rank(-C[i, ])[i])) }
rA <- do.call(rbind, rA); rA$dmo <- rA$diag - rA$off; kA <- do.call(rbind, kA)
write.csv(rA, file.path(OUT, "spec98_H2052_crop_vs_bulk.csv"), row.names = FALSE); write.csv(kA, file.path(OUT, "spec98_H2052_per_KO.csv"), row.names = FALSE)
options(width = 220); cat("\n=== A  H2052 CROP-seq vs bulk ===\n"); print(rA[, c("space", "mode", "diag", "off", "dmo", "perm_p", "top1")], digits = 3, row.names = FALSE)
for (md in names(MODES)) { x <- kA[kA$mode == md, ]; cat("\n-- per TF own rho (rank),", MODES[[md]], "\n"); print(reshape(transform(x[, c("KO", "space")], o = sprintf("%.2f (%s)", x$own, format(x$rank))), idvar = "KO", timevar = "space", direction = "wide"), row.names = FALSE) }

## ---------------------------------------------------------------- B  across bulk lines
SPB <- c("signal500", "top2000", "all", "Cm20"); vAll <- pr$genes[order(rowMeans(pr$gene_var), decreasing = TRUE)]; rB <- list(); kB <- list(); AG <- list()
for (md in names(MODES)) { A <- array(NA_real_, c(length(pr$genes), 8, 4), dimnames = list(pr$genes, KOS, LINES)); for (cl in LINES) A[, , cl] <- norm9(bulkG[[cl]][pr$genes, ], md); Am <- array(NA_real_, c(20, 8, 4), dimnames = list(CM, KOS, LINES)); for (cl in LINES) Am[, , cl] <- norm9(bulkG[[cl]][CM, ], md)
  for (sp in SPB) { CLs <- switch(sp, Cm20 = cl_from_array(Am, NULL, centre = FALSE), signal500 = cl_from_array(A, "signal", pr$gene_signal, centre = FALSE), top2000 = cl_from_array(A, vAll[1:2000], centre = FALSE), all = cl_from_array(A, NULL, centre = FALSE))
    a <- agreement(CLs, permsL); AG[[paste(md, sp)]] <- a; rB[[length(rB) + 1]] <- data.frame(mode = md, space = sp, t(a$global)); kB[[length(kB) + 1]] <- cbind(mode = md, space = sp, a$per_KO) } }
rB <- do.call(rbind, rB); kB <- do.call(rbind, kB); write.csv(rB, file.path(OUT, "spec98_across_lines.csv"), row.names = FALSE); write.csv(kB, file.path(OUT, "spec98_across_lines_per_KO.csv"), row.names = FALSE)
cat("\n=== B  bulk across lines: global ===\n"); print(rB[, c("space", "mode", "diag_mean", "off_mean", "perm_p", "perm_z")], digits = 3, row.names = FALSE)
for (sp in SPB) { x <- kB[kB$space == sp, ]; cat("\n-- per TF perm p,", sp, "\n"); print(reshape(transform(x[, c("KO", "mode")], p = sprintf("%.4f (r %.2f)", x$perm_p, x$mean_intra)), idvar = "KO", timevar = "mode", direction = "wide"), row.names = FALSE) }

## ---------------------------------------------------------------- figure
fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
save_fig("fig11_nine_group_centring", 16, 10.4, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE), widths = c(1, 1.1, 1.1)); par(oma = c(0, 0, 4.4, 0.4), family = "sans"); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.8); mtext(b, side = 3, adj = 0, line = 0.4, cex = 0.58, col = INK2) }
  par(mar = c(6, 4.6, 3.4, 1)); D <- sapply(names(MODES), function(m) rA$dmo[rA$mode == m][match(names(SPA), rA$space[rA$mode == m])]); P <- sapply(names(MODES), function(m) rA$perm_p[rA$mode == m][match(names(SPA), rA$space[rA$mode == m])])
  bp <- barplot(t(D), beside = TRUE, col = MC, border = NA, las = 1, names.arg = rep("", length(SPA)), ylab = "own-KO minus other-KO Spearman", cex.axis = 0.8, ylim = range(c(0, D)) * c(1.3, 1.5)); abline(h = 0, col = INK)
  text(bp, pmax(t(D), 0) + 0.004, sprintf("p %s", fp(t(P))), cex = 0.55, srt = 90, adj = c(0, 0.5), font = ifelse(t(P) < 0.05, 2, 1)); text(colMeans(bp), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), SPAL, srt = 30, adj = 1, xpd = NA, cex = 0.72)
  legend("topright", MODES, fill = MC, border = NA, bty = "n", cex = 0.62); ttl("A  H2052, CROP-seq vs bulk: all 8 KOs", "exact permutation p of the bulk KO labels")
  for (md in c("c9", "c8")) { x <- kA[kA$mode == md & kA$space == "top2000", ]; o <- x$KO[order(-x$own)]; x <- x[match(o, x$KO), ]; par(mar = c(4.6, 4.6, 3.4, 1))
    bp <- barplot(x$own, col = adjustcolor(KO_COL[o], 0.8), border = NA, las = 1, ylim = c(min(x$own, x$median_other) - 0.08, max(x$own) + 0.2), ylab = "Spearman, CROP-seq KO vs bulk H2052 KO", names.arg = rep("", 8), cex.axis = 0.8); abline(h = 0, col = INK)
    segments(bp - 0.4, x$median_other, bp + 0.4, x$median_other, lty = 2, col = "grey25", lwd = 1.3); text(bp, pmax(x$own, 0) + 0.02, sprintf("#%s\np=%s", format(x$rank), fp(x$rank / 8)), cex = 0.62, adj = c(0.5, 0), font = ifelse(x$rank == 1, 2, 1))
    text(bp, par("usr")[3] - 0.02 * diff(par("usr")[3:4]), o, srt = 45, adj = 1, xpd = NA, cex = 0.8, font = 2, col = KO_COL[o]); g <- rA[rA$mode == md & rA$space == "top2000", ]
    ttl(sprintf("%s  H2052 per TF -- %s", if (md == "c9") "B" else "C", MODES[[md]]), sprintf("top 2000 variable genes; dash = median vs other bulk KOs; # = rank among 8; all KOs p = %s", fp(g$perm_p))) }
  o2 <- read.csv(file.path(OUT, "spec82_TF_agreement_ranked.csv"))$KO
  for (md in c("c9", "c8")) { a <- AG[[paste(md, "signal500")]]; k <- a$per_KO; yr <- range(c(unlist(a$cross), a$intra)); yr <- yr + c(-0.02, 0.17) * diff(yr); par(mar = c(4.6, 4.6, 3.4, 1))
    plot(NA, xlim = c(0.5, 8.5), ylim = yr, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "Spearman between two bulk cell lines", cex.axis = 0.8); abline(h = 0, lty = 2, col = INK)
    for (i in seq_along(o2)) { j <- match(o2[i], k$KO); boxplot(a$cross[[j]], at = i - 0.2, add = TRUE, boxwex = 0.32, col = "#d5d4cf", border = INK2, outline = FALSE, axes = FALSE)
      boxplot(a$intra[o2[i], ], at = i + 0.2, add = TRUE, boxwex = 0.32, col = adjustcolor(KO_COL[[o2[i]]], 0.65), border = INK2, outline = FALSE, axes = FALSE); points(i + 0.2 + seq(-0.1, 0.1, length.out = 6), a$intra[o2[i], ], pch = 16, cex = 0.75, col = INK2)
      text(i, yr[2], sprintf("p=%s\n%s", fp(k$perm_p[j]), ifelse(k$perm_fdr[j] < 0.01, "**", ifelse(k$perm_fdr[j] < 0.05, "*", ""))), cex = 0.66, adj = c(0.5, 1), font = ifelse(k$perm_p[j] < 0.05, 2, 1)) }
    text(1:8, par("usr")[3] - 0.03 * diff(par("usr")[3:4]), o2, srt = 45, adj = 1, xpd = NA, cex = 0.8, font = 2, col = KO_COL[o2])
    ttl(sprintf("%s  across the 4 bulk lines -- %s", if (md == "c9") "D" else "E", MODES[[md]]), sprintf("top-500 KO-signal genes; colour = same KO, grey = other KOs; global p = %s", fp(a$global[["perm_p"]]))) }
  par(mar = c(6, 4.6, 3.4, 1)); sps <- SPB; Z <- sapply(names(MODES), function(m) rB$perm_z[rB$mode == m][match(sps, rB$space[rB$mode == m])]); Pz <- sapply(names(MODES), function(m) rB$perm_p[rB$mode == m][match(sps, rB$space[rB$mode == m])])
  bp <- barplot(t(Z), beside = TRUE, col = MC, border = NA, las = 1, names.arg = c("top-500\nKO-signal", "top 2000\nvariable", "all genes", "20 Cm\nmodules"), ylab = "global agreement across lines (permutation z)", cex.axis = 0.8, cex.names = 0.75, ylim = range(c(0, Z)) * c(1.1, 1.3)); abline(h = qnorm(0.95), lty = 3, col = MUTED)
  text(bp, pmax(t(Z), 0) + 0.05, sprintf("p %s", fp(t(Pz))), cex = 0.58, adj = c(0.5, 0), font = ifelse(t(Pz) < 0.05, 2, 1)); ttl("F  across bulk lines: all KOs", "KO-label permutation within lines")
  mtext("Profiles centred on the mean of the 9 groups (NTC + 8 KOs) of each data set.   Top: NCI-H2052, CROP-seq vs bulk.   Bottom: bulk, across the 4 cell lines", outer = TRUE, side = 3, adj = 0, line = 2.8, font = 2, cex = 1)
  mtext("bulk: group mean TMM logCPM per line (S70); screen: pseudobulk logCPM per guide group; no NTC subtraction; grey = the earlier centring over the 8 KOs only, for reference", outer = TRUE, side = 3, adj = 0, line = 1.5, cex = 0.65, col = INK2)
})
