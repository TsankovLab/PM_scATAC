# Spec 95: CROP-seq vs bulk NCI-H2052 (same cell line).
# Part 1 -- shared feature space: genes in the top N of BOTH data sets, N = 1000 / 2000 / 3000, ranked by (a) variance across the KO groups (screen: pseudobulk logCPM over NTC + 8 KOs;
#           bulk: H2052 sample variance) or (b) KO signal (max |t| over KOs; screen Welch t vs NTC, bulk limma t). Views: NTC - KO centred over KOs | t vs other KOs (spec94).
#           Metric: own-KO minus other-KO Spearman, exact permutation p of bulk KO labels (8!), own-KO rank.
# Part 2 -- the screen axes that split the TF clusters: PCA of the screen KO profiles (pseudobulk NTC - KO, centred over the 8 KOs, screen HVGs); top 150 genes at each end of PC1 / PC2.
#           Bulk H2052 KO score on each axis = mean centred NTC - KO over the + genes minus over the - genes; agreement with the screen PC scores over the 8 KOs (Spearman, exact
#           permutation p); also the same top genes as a feature space for the KO x KO match.
# needs: spec88 (pseudobulk logCPM incl. all genes via load_screen), spec72 (S70), spec1 (screen stats)
# outputs: spec95_shared_space.csv, spec95_shared_space_per_KO.csv, spec95_screen_axes.csv, spec95_screen_axis_genes.csv, fig9Q_H2052_shared_space.{png,pdf}, fig9R_H2052_screen_axes.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages({ library(Seurat); library(ComplexHeatmap); library(circlize); library(grid) })
CL <- "NCI-H2052"; pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); ss <- readRDS(file.path(OUT, "spec_screen_stats.rds"))
srt <- load_screen(); grp <- as.character(srt$merged_call); X <- GetAssayData(srt, layer = "counts"); Y <- GetAssayData(srt, layer = "data")
PB <- sapply(c("NTC", KOS), function(k) Matrix::rowSums(X[, grp == k, drop = FALSE])); lc <- edgeR::cpm(edgeR::DGEList(PB), log = TRUE, prior.count = 2)
G <- Reduce(intersect, list(pr$genes, rownames(lc)[rowSums(edgeR::cpm(edgeR::DGEList(PB)) >= 10) >= 1], ss$genes)); cat("genes expressed in both:", length(G), "\n")
cen <- function(P) P - rowMeans(P)
welch_rest <- function(f) { Yf <- Y[f, ]; cells <- grp %in% KOS; sapply(KOS, function(k) { a <- which(grp == k); b <- which(cells & grp != k)
  ma <- Matrix::rowMeans(Yf[, a]); mb <- Matrix::rowMeans(Yf[, b]); va <- (Matrix::rowSums(Yf[, a]^2) - length(a) * ma^2) / (length(a) - 1); vb <- (Matrix::rowSums(Yf[, b]^2) - length(b) * mb^2) / (length(b) - 1)
  t <- (mb - ma) / sqrt(pmax(va, 1e-8) / length(a) + pmax(vb, 1e-8) / length(b)); t[!is.finite(t)] <- 0; t }) }
S_fc <- cen(SGN * (lc[G, "NTC"] - lc[G, KOS])); B_fc <- cen(pr$P$tmm$none$coef[G, , CL]); S_tr <- welch_rest(G); B_tr <- pr$P$tmm$none$t_rest[G, , CL]
perms <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(1:8)) })
score <- function(S, B) { C <- cor(S, B, method = "spearman"); dm <- mean(diag(C)); nul <- apply(perms, 1, function(p) mean(C[cbind(1:8, p)]))
  list(g = data.frame(diag = dm, off = mean(C[row(C) != col(C)]), perm_p = mean(nul >= dm - 1e-12), top1 = sum(apply(C, 1, which.max) == 1:8)),
       k = data.frame(KO = KOS, own = diag(C), median_other = sapply(1:8, function(i) median(C[i, -i])), rank = sapply(1:8, function(i) rank(-C[i, ])[i]))) }

## ---------------------------------------------------------------- Part 1: shared feature spaces
rk_var <- list(screen = G[order(matrixStats::rowVars(lc[G, ]), decreasing = TRUE)], bulk = G[order(pr$gene_var[G, CL], decreasing = TRUE)])
s_sig <- apply(abs(-ss$t0[G, KOS]), 1, max); s_sig[!is.finite(s_sig)] <- 0; rk_sig <- list(screen = G[order(s_sig, decreasing = TRUE)], bulk = G[order(apply(abs(pr$P$tmm$none$t[G, , CL]), 1, max), decreasing = TRUE)])
SP <- list(); for (n in c(1000, 2000, 3000)) { SP[[sprintf("var_top%d", n)]] <- intersect(rk_var$screen[1:n], rk_var$bulk[1:n]); SP[[sprintf("signal_top%d", n)]] <- intersect(rk_sig$screen[1:n], rk_sig$bulk[1:n]) }
SP$all <- G
r1 <- list(); k1 <- list()
for (sp in names(SP)) { f <- SP[[sp]]; for (v in c("fc", "t")) { a <- if (v == "fc") score(S_fc[f, ], B_fc[f, ]) else score(S_tr[f, ], B_tr[f, ])
  r1[[length(r1) + 1]] <- cbind(space = sp, n_genes = length(f), view = v, a$g); k1[[length(k1) + 1]] <- cbind(space = sp, view = v, a$k) } }
r1 <- do.call(rbind, r1); r1$dmo <- r1$diag - r1$off; k1 <- do.call(rbind, k1)
write.csv(r1, file.path(OUT, "spec95_shared_space.csv"), row.names = FALSE); write.csv(k1, file.path(OUT, "spec95_shared_space_per_KO.csv"), row.names = FALSE)
options(width = 220); cat("\n=== Part 1: shared feature spaces ===\n"); print(r1[, c("space", "n_genes", "view", "diag", "off", "dmo", "perm_p", "top1")], digits = 3, row.names = FALSE)
for (v in c("fc", "t")) { x <- k1[k1$view == v, ]; cat("\n-- own rho (rank),", v, "\n"); print(reshape(transform(x[, c("KO", "space")], o = sprintf("%.2f (%d)", x$own, x$rank)), idvar = "KO", timevar = "space", direction = "wide"), row.names = FALSE) }

## ---------------------------------------------------------------- Part 2: screen axes that split the TF clusters
# Pseudobulk PCs are dominated by the two smallest groups (TEAD2 20 cells, HMGA1 13), so the axes are defined (i) directly from the screen clusters of spec89
# (A = TCF3 / TEAD4 / TWIST1, B = PITX1 / BPTF): genes ranked by the screen Welch t of A cells vs B cells; and (ii) by PCA of the noise-weighted screen t profiles (t vs other KOs).
CA <- c("TCF3", "TEAD4", "TWIST1"); CB <- c("PITX1", "BPTF"); NTOP <- 150; Hs <- G
tAB <- { Yf <- Y[Hs, ]; a <- which(grp %in% CA); b <- which(grp %in% CB); ma <- Matrix::rowMeans(Yf[, a]); mb <- Matrix::rowMeans(Yf[, b])
  va <- (Matrix::rowSums(Yf[, a]^2) - length(a) * ma^2) / (length(a) - 1); vb <- (Matrix::rowSums(Yf[, b]^2) - length(b) * mb^2) / (length(b) - 1); t <- (ma - mb) / sqrt(pmax(va, 1e-8) / length(a) + pmax(vb, 1e-8) / length(b)); t[!is.finite(t)] <- 0; t }
pct <- prcomp(t(S_tr[Hs, ]), center = TRUE); vet <- 100 * pct$sdev^2 / sum(pct$sdev^2)
AXES <- list(clusterAB = tAB, tPC1 = pct$rotation[, 1], tPC2 = pct$rotation[, 2]); AXL <- c(clusterAB = sprintf("screen cluster axis: %s vs %s", paste(CA, collapse = "/"), paste(CB, collapse = "/")),
  tPC1 = sprintf("screen t-profile PC1 (%.0f%%)", vet[1]), tPC2 = sprintf("screen t-profile PC2 (%.0f%%)", vet[2]))
Sfc <- S_fc[Hs, ]; Bfc <- B_fc[Hs, ]; Btr <- B_tr[Hs, ]; ax <- list(); genes_out <- list()
for (an in names(AXES)) { ld <- AXES[[an]]; up <- names(sort(ld, decreasing = TRUE))[1:NTOP]; dn <- names(sort(ld))[1:NTOP]
  sc_s <- colMeans(Sfc[up, ]) - colMeans(Sfc[dn, ]); sc_b <- colMeans(Bfc[up, ]) - colMeans(Bfc[dn, ])
  obs <- cor(sc_s, sc_b, method = "spearman"); nul <- apply(perms, 1, function(p) cor(sc_s, sc_b[p], method = "spearman"))
  sepB <- mean(sc_b[CA]) - mean(sc_b[CB]); nB <- apply(perms, 1, function(p) { x <- sc_b[p]; mean(x[match(CA, KOS)]) - mean(x[match(CB, KOS)]) })
  m <- score(S_fc[c(up, dn), ], B_fc[c(up, dn), ]); mt <- score(S_tr[c(up, dn), ], B_tr[c(up, dn), ])
  ax[[an]] <- data.frame(axis = an, KO = KOS, screen_score = sc_s[KOS], bulk_score = sc_b[KOS], axis_rho = obs, axis_perm_p = mean(nul >= obs - 1e-12),
                         bulk_clusterA_minus_B = sepB, bulk_sep_p = mean(abs(nB) >= abs(sepB) - 1e-12), KOxKO_dmo_fc = m$g$diag - m$g$off, KOxKO_p_fc = m$g$perm_p, KOxKO_dmo_t = mt$g$diag - mt$g$off, KOxKO_p_t = mt$g$perm_p)
  genes_out[[an]] <- data.frame(axis = an, end = rep(c("+", "-"), each = NTOP), gene = c(up, dn), screen_axis_value = ld[c(up, dn)]) }
ax <- do.call(rbind, ax); write.csv(ax, file.path(OUT, "spec95_screen_axes.csv"), row.names = FALSE); write.csv(do.call(rbind, genes_out), file.path(OUT, "spec95_screen_axis_genes.csv"), row.names = FALSE)
cat("\n=== Part 2: screen axes, scored in bulk H2052 ===\n"); print(unique(ax[, c("axis", "axis_rho", "axis_perm_p", "bulk_clusterA_minus_B", "bulk_sep_p", "KOxKO_dmo_fc", "KOxKO_p_fc", "KOxKO_dmo_t", "KOxKO_p_t")]), digits = 3, row.names = FALSE)
print(reshape(ax[, c("axis", "KO", "screen_score", "bulk_score")], idvar = "KO", timevar = "axis", direction = "wide"), digits = 2, row.names = FALSE)
hall <- { suppressPackageStartupMessages(library(msigdbr)); h <- msigdbr(species = "Homo sapiens", category = "H"); split(h$gene_symbol, sub("HALLMARK_", "", h$gs_name)) }
for (an in names(AXES)) { ld <- AXES[[an]]; idx <- lapply(hall, function(x) which(names(ld) %in% x)); idx <- idx[sapply(idx, length) >= 10]; cr <- limma::cameraPR(ld, idx)
  cat(sprintf("\nHallmark on %s (Up = + end):\n", an)); print(head(cr[, c("NGenes", "Direction", "PValue", "FDR")], 6), digits = 3) }

## ---------------------------------------------------------------- figures
VC <- c(fc = "#2a78d6", t = "#1baf7a"); VL <- c(fc = "NTC - KO, centred", t = "t vs other KOs")
save_fig("fig9Q_H2052_shared_space", 15, 7.8, function() {
  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE), heights = c(1, 1.1)); par(oma = c(0, 0, 3.2, 0.4), family = "sans")
  par(mar = c(6.6, 4.4, 2.6, 0.8)); sps <- names(SP); D <- sapply(names(VL), function(v) r1$dmo[r1$view == v][match(sps, r1$space[r1$view == v])]); Pp <- sapply(names(VL), function(v) r1$perm_p[r1$view == v][match(sps, r1$space[r1$view == v])])
  ng <- r1$n_genes[r1$view == "fc"][match(sps, r1$space[r1$view == "fc"])]
  bp <- barplot(t(D), beside = TRUE, col = VC, border = NA, las = 1, names.arg = rep("", length(sps)), ylab = "own-KO minus other-KO Spearman", cex.axis = 0.75, ylim = range(c(0, D)) * c(1.15, 1.35)); abline(h = 0, col = INK)
  text(bp, pmax(t(D), 0) + 0.004, sprintf("p %.2f", t(Pp)), cex = 0.55, adj = c(0.5, 0), font = ifelse(t(Pp) < 0.05, 2, 1)); text(colMeans(bp), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), sprintf("%s\n(%d genes)", sps, ng), srt = 30, adj = 1, xpd = NA, cex = 0.65)
  legend("topright", VL, fill = VC, border = NA, bty = "n", cex = 0.7); mtext("A  shared feature spaces (top N in both screen and bulk H2052, by variance or by KO signal)", side = 3, adj = 0, line = 0.6, font = 2, cex = 0.8)
  for (v in names(VL)) { x <- k1[k1$view == v, ]; M <- sapply(sps, function(s) x$own[x$space == s][match(KOS, x$KO[x$space == s])]); R <- sapply(sps, function(s) x$rank[x$space == s][match(KOS, x$KO[x$space == s])])
    par(mar = c(6.6, 5.4, 2.6, 1)); image(seq_along(sps), 1:8, t(pmax(pmin(M, 0.6), -0.6))[, 8:1], col = colorRampPalette(c("#2a78d6", "#f0efec", "#eb6834"))(101), zlim = c(-0.6, 0.6), axes = FALSE, xlab = "", ylab = "")
    for (i in 1:8) for (jj in seq_along(sps)) text(jj, 9 - i, sprintf("%.2f\n#%d", M[i, jj], R[i, jj]), cex = 0.5, font = ifelse(R[i, jj] == 1, 2, 1))
    text(seq_along(sps), 0.3, sps, srt = 30, adj = 1, xpd = NA, cex = 0.62); text(0.4, 8:1, KOS, adj = 1, xpd = NA, cex = 0.72, font = 2, col = KO_COL[KOS]); box(col = MUTED)
    mtext(sprintf("%s  per TF: own-KO rho, # = rank among 8 -- %s", c(fc = "B", t = "C")[v], VL[[v]]), side = 3, adj = 0, line = 0.6, font = 2, cex = 0.75) }
  mtext("CROP-seq vs bulk NCI-H2052 on feature spaces shared by both data sets", outer = TRUE, side = 3, adj = 0, line = 1.7, font = 2, cex = 1)
  mtext("var_topN = genes among the top N most variable in the screen pseudobulk AND in H2052 bulk; signal_topN = same with KO signal (max |t| over KOs); p = exact permutation of the bulk KO labels", outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.6, col = INK2)
})
save_fig("fig9R_H2052_screen_axes", 15, 5.6, function() {
  layout(matrix(1:4, nrow = 1), widths = c(1, 1, 1, 0.9)); par(oma = c(0, 0, 3.2, 0.4), mar = c(4.4, 4.4, 2.8, 1), family = "sans")
  for (j in seq_along(AXES)) { an <- names(AXES)[j]; d <- ax[ax$axis == an, ]; plot(d$screen_score, d$bulk_score, pch = 21, bg = KO_COL[d$KO], col = ifelse(d$KO %in% CA, "black", ifelse(d$KO %in% CB, "grey40", "white")), lwd = 2, cex = 2.2, bty = "n", las = 1, cex.axis = 0.75,
      xlab = "CROP-seq KO score on the axis genes", ylab = "bulk H2052 KO score on the same genes"); abline(h = 0, v = 0, col = GRID)
    text(d$screen_score, d$bulk_score, d$KO, pos = 3, cex = 0.7, col = KO_COL[d$KO], font = 2, xpd = NA)
    mtext(sprintf("%s  %s", LETTERS[j], AXL[[an]]), side = 3, adj = 0, line = 1.3, font = 2, cex = 0.66)
    mtext(sprintf("KO scores screen vs bulk rho %.2f (p %.2f); bulk cluster A - B %+.2f (p %.2f)", d$axis_rho[1], d$axis_perm_p[1], d$bulk_clusterA_minus_B[1], d$bulk_sep_p[1]), side = 3, adj = 0, line = 0.3, cex = 0.52, col = INK2) }
  par(mar = c(6, 4.4, 2.8, 1)); M <- rbind(fc = sapply(names(AXES), function(a) ax$KOxKO_dmo_fc[ax$axis == a][1]), t = sapply(names(AXES), function(a) ax$KOxKO_dmo_t[ax$axis == a][1]))
  Pm <- rbind(sapply(names(AXES), function(a) ax$KOxKO_p_fc[ax$axis == a][1]), sapply(names(AXES), function(a) ax$KOxKO_p_t[ax$axis == a][1]))
  bp <- barplot(M, beside = TRUE, col = VC, border = NA, names.arg = c("cluster\naxis", "t PC1", "t PC2"), las = 1, ylab = "own-KO minus other-KO Spearman", cex.axis = 0.75, ylim = range(c(0, M)) * c(1.2, 1.4)); abline(h = 0, col = INK)
  text(bp, pmax(M, 0) + 0.005, sprintf("p %.2f", Pm), cex = 0.6, adj = c(0.5, 0)); legend("topright", VL, fill = VC, border = NA, bty = "n", cex = 0.6); mtext("D  KO x KO match on the axis genes", side = 3, adj = 0, line = 0.6, font = 2, cex = 0.7)
  mtext("Do the CROP-seq axes that split the TF clusters also separate the KOs in bulk NCI-H2052?", outer = TRUE, side = 3, adj = 0, line = 1.7, font = 2, cex = 1)
  mtext(sprintf("top %d genes at each end of each screen axis; KO score = mean centred NTC - KO over + genes minus - genes; black outline = cluster A (%s), grey = cluster B (%s); p = exact permutation of the 8 KO labels", NTOP, paste(CA, collapse = "/"), paste(CB, collapse = "/")), outer = TRUE, side = 3, adj = 0, line = 0.5, cex = 0.6, col = INK2)
})
