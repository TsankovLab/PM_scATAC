# Spec 99: spec98 (profiles centred on the mean of the 9 groups, NTC + 8 KOs, of each data set) on the INTERSECTION of the variable genes of the data sets compared,
# and a check for a bulk cell line that consistently underperforms (removed and rerun).
# A  NCI-H2052, CROP-seq vs bulk: genes = top-N screen HVGs (Seurat vst on the cells) that are also in the top-N variable genes of bulk H2052 (variance of TMM logCPM over its samples); N = 1000 / 2000 / 3000 / 5000
# B  bulk across lines: genes = intersection of the top-N variable genes of every line compared (N = 2000 / 3000 / 5000); top-500 KO-signal genes and top-2000 mean variance kept as reference
# C  per line: mean own-minus-other Spearman of the 3 line pairs that contain the line, per feature space -> the underperforming line = lowest average rank
# D  leave-one-line-out: agreement of the remaining 3 lines (gene intersections recomputed on those 3 lines; KO-signal genes from the third line)
# modes: c9 = minus mean of the 9 groups | z9 = z over the 9 groups | c8 = minus mean of the 8 KOs (earlier)
# needs: spec72 (S70), spec_common.R
# outputs: spec99_H2052_intersect.csv, spec99_H2052_intersect_per_KO.csv, spec99_across_lines_intersect.csv, spec99_per_line.csv, spec99_leave_one_line_out.csv, spec99_leave_one_line_out_per_KO.csv, fig12_intersect_variable_genes_drop_line.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages(library(Seurat)); N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000")); permsL <- make_perms(N_PERM)
CL <- "NCI-H2052"; GR <- c("NTC", KOS); fab <- load_full_audited_bulk(); meta <- fab$meta; pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); ids <- pr$samples
srt <- load_screen(); grp <- as.character(srt$merged_call); X <- GetAssayData(srt, layer = "counts")
PB <- sapply(GR, function(k) Matrix::rowSums(X[, grp == k, drop = FALSE])); lcS <- edgeR::cpm(edgeR::DGEList(PB), log = TRUE, prior.count = 2); expS <- rownames(lcS)[rowSums(edgeR::cpm(edgeR::DGEList(PB)) >= 10) >= 1]
G <- intersect(pr$genes, expS); hvS <- intersect(VariableFeatures(FindVariableFeatures(srt, selection.method = "vst", nfeatures = nrow(srt), verbose = FALSE)), G)      # screen genes ranked by vst variance
bulkG <- lapply(setNames(LINES, LINES), function(cl) { s <- ids[meta[ids, "Cell_Line"] == cl]; y <- normalise(fab$cs[pr$genes, s], "tmm")$y; g <- meta[s, "Gene_Control"]; sapply(GR, function(k) rowMeans(y[, g == k, drop = FALSE])) })
norm9 <- function(M, mode) { M <- M[, GR]; switch(mode, c9 = (M - rowMeans(M))[, KOS], z9 = { Z <- (M - rowMeans(M)) / matrixStats::rowSds(M); Z[!is.finite(Z)] <- 0; Z[, KOS] }, c8 = M[, KOS] - rowMeans(M[, KOS])) }
MODES <- c(c9 = "minus mean of the 9 groups", z9 = "z over the 9 groups", c8 = "minus mean of the 8 KOs (earlier)"); MC <- c(c9 = "#eb6834", z9 = "#eda100", c8 = "grey55")
perms8 <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(1:8)) })
topv <- function(cl, n, pool = pr$genes) head(pool[order(pr$gene_var[pool, cl], decreasing = TRUE)], n); options(width = 220)

## ---------------------------------------------------------------- A  H2052 screen vs bulk, intersection of variable genes
NS <- c(1000, 2000, 3000, 5000); SPA <- setNames(lapply(NS, function(n) intersect(hvS[seq_len(n)], topv(CL, n, G))), paste0("int", NS)); rA <- list(); kA <- list()
for (sp in names(SPA)) for (md in names(MODES)) { f <- SPA[[sp]]; S <- norm9(lcS[f, ], md); B <- norm9(bulkG[[CL]][f, ], md); C <- cor(S, B, method = "spearman"); dm <- mean(diag(C)); nul <- apply(perms8, 1, function(p) mean(C[cbind(1:8, p)]))
  rA[[length(rA) + 1]] <- data.frame(space = sp, n_genes = length(f), mode = md, diag = dm, off = mean(C[row(C) != col(C)]), perm_p = mean(nul >= dm - 1e-12), top1 = sum(apply(C, 1, which.max) == 1:8))
  kA[[length(kA) + 1]] <- data.frame(space = sp, mode = md, KO = KOS, own = diag(C), median_other = sapply(1:8, function(i) median(C[i, -i])), rank = sapply(1:8, function(i) rank(-C[i, ])[i])) }
rA <- do.call(rbind, rA); rA$dmo <- rA$diag - rA$off; kA <- do.call(rbind, kA); write.csv(rA, file.path(OUT, "spec99_H2052_intersect.csv"), row.names = FALSE); write.csv(kA, file.path(OUT, "spec99_H2052_intersect_per_KO.csv"), row.names = FALSE)
cat("\n=== A  H2052 CROP-seq vs bulk, intersection of top-N variable genes ===\n"); print(rA[, c("space", "n_genes", "mode", "diag", "off", "dmo", "perm_p", "top1")], digits = 3, row.names = FALSE)
for (md in names(MODES)) { x <- kA[kA$mode == md, ]; cat("\n-- per TF own rho (rank),", MODES[[md]], "\n"); print(reshape(transform(x[, c("KO", "space")], o = sprintf("%.2f (%s)", x$own, format(x$rank))), idvar = "KO", timevar = "space", direction = "wide"), row.names = FALSE) }

## ---------------------------------------------------------------- B-D  bulk across lines
AR <- lapply(setNames(names(MODES), names(MODES)), function(md) { A <- array(NA_real_, c(length(pr$genes), 8, 4), dimnames = list(pr$genes, KOS, LINES)); for (cl in LINES) A[, , cl] <- norm9(bulkG[[cl]], md); A })
SPB <- c(signal500 = "top-500 KO-signal", meanvar2000 = "top 2000 mean variance", int2000 = "intersect top 2000", int3000 = "intersect top 3000", int5000 = "intersect top 5000")
run <- function(md, sp, lines) { prs <- Filter(function(p) all(p %in% lines), PAIRS); A <- AR[[md]]
  f <- switch(sp, signal500 = "signal", meanvar2000 = pr$genes[order(rowMeans(pr$gene_var[, lines]), decreasing = TRUE)][1:2000], Reduce(intersect, lapply(lines, topv, n = as.integer(sub("int", "", sp)))))
  CLs <- if (identical(f, "signal")) cl_from_array(A, "signal", pr$gene_signal, centre = FALSE, pairs = prs, lines_avail = lines) else cl_from_array(A, f, centre = FALSE, pairs = prs)
  a <- agreement(CLs, permsL, pairs = prs); a$n_genes <- if (identical(f, "signal")) 500L else length(f); a$pairs <- prs; a$pair_dmo <- sapply(CLs, function(C) mean(diag(C)) - mean(C[row(C) != col(C)])); a }
rB <- list(); pl <- list(); AG <- list()
for (md in names(MODES)) for (sp in names(SPB)) { a <- run(md, sp, LINES); AG[[paste(md, sp)]] <- a; rB[[length(rB) + 1]] <- data.frame(mode = md, space = sp, n_genes = a$n_genes, t(a$global))
  pl[[length(pl) + 1]] <- data.frame(mode = md, space = sp, line = LINES, pair_dmo = sapply(LINES, function(l) mean(a$pair_dmo[sapply(a$pairs, function(p) l %in% p)])), without = sapply(LINES, function(l) mean(a$pair_dmo[!sapply(a$pairs, function(p) l %in% p)]))) }
rB <- do.call(rbind, rB); pl <- do.call(rbind, pl); pl$rank <- ave(pl$pair_dmo, pl$mode, pl$space, FUN = rank); write.csv(rB, file.path(OUT, "spec99_across_lines_intersect.csv"), row.names = FALSE); write.csv(pl, file.path(OUT, "spec99_per_line.csv"), row.names = FALSE)
cat("\n=== B  bulk across the 4 lines ===\n"); print(rB[, c("space", "n_genes", "mode", "diag_mean", "off_mean", "perm_p", "perm_z")], digits = 3, row.names = FALSE)
cat("\n=== C  per line: mean own-minus-other of the 3 pairs containing the line (rank 1 = worst) ===\n"); x <- pl[pl$mode == "c9", ]; print(reshape(transform(x[, c("line", "space")], v = sprintf("%.3f (%s)", x$pair_dmo, format(x$rank))), idvar = "line", timevar = "space", direction = "wide"), row.names = FALSE)
mr <- tapply(pl$rank, pl$line, mean)[LINES]; cat("mean rank over spaces and modes:\n"); print(round(mr, 2)); WORST <- names(which.min(mr)); cat("underperforming line:", WORST, "(worst in", sum(pl$rank[pl$line == WORST] == 1), "of", sum(pl$line == WORST), "space x mode)\n")
rL <- list(); kL <- list(); LO <- list()
for (dl in LINES) for (md in names(MODES)) for (sp in names(SPB)) { a <- run(md, sp, setdiff(LINES, dl)); LO[[paste(dl, md, sp)]] <- a; rL[[length(rL) + 1]] <- data.frame(dropped = dl, mode = md, space = sp, n_genes = a$n_genes, t(a$global)); kL[[length(kL) + 1]] <- cbind(dropped = dl, mode = md, space = sp, a$per_KO) }
rL <- do.call(rbind, rL); kL <- do.call(rbind, kL); write.csv(rL, file.path(OUT, "spec99_leave_one_line_out.csv"), row.names = FALSE); write.csv(kL, file.path(OUT, "spec99_leave_one_line_out_per_KO.csv"), row.names = FALSE)
cat("\n=== D  leave one line out (c9): global diag - off, perm p ===\n"); x <- rL[rL$mode == "c9", ]; print(reshape(transform(x[, c("dropped", "space")], v = sprintf("%.3f p=%.4f", x$diag_mean - x$off_mean, x$perm_p)), idvar = "dropped", timevar = "space", direction = "wide"), row.names = FALSE)
cat("\n-- all 4 lines (c9) for reference\n"); x <- rB[rB$mode == "c9", ]; print(data.frame(space = x$space, v = sprintf("%.3f p=%.4f", x$diag_mean - x$off_mean, x$perm_p)), row.names = FALSE)
cat("\n-- per TF perm p (mean same-KO rho) without", WORST, "vs all 4 lines, c9\n")
for (sp in names(SPB)) { a <- AG[[paste("c9", sp)]]$per_KO; b <- LO[[paste(WORST, "c9", sp)]]$per_KO; cat(sp, "\n"); print(data.frame(KO = a$KO, all4 = sprintf("%.4f (r %.2f)", a$perm_p, a$mean_intra), dropped = sprintf("%.4f (r %.2f)", b$perm_p, b$mean_intra), fdr_dropped = signif(b$perm_fdr, 2)), row.names = FALSE) }

## ---------------------------------------------------------------- figure
fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2)); SH <- function(l) sub("NCI-|-211H", "", l)
save_fig("fig12_intersect_variable_genes_drop_line", 16, 10.4, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE), widths = c(1, 1.1, 1.1)); par(oma = c(0, 0, 4.4, 0.4), family = "sans"); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.8); mtext(b, side = 3, adj = 0, line = 0.4, cex = 0.58, col = INK2) }
  ## A
  par(mar = c(6, 4.6, 3.4, 1)); gm <- function(col) sapply(names(MODES), function(m) rA[[col]][rA$mode == m][match(names(SPA), rA$space[rA$mode == m])]); D <- gm("dmo"); P <- gm("perm_p")
  bp <- barplot(t(D), beside = TRUE, col = MC, border = NA, las = 1, names.arg = rep("", length(SPA)), ylab = "own-KO minus other-KO Spearman", cex.axis = 0.8, ylim = range(c(0, D)) + c(-0.01, 0.035)); abline(h = 0, col = INK)
  text(bp, pmax(t(D), 0) + 0.003, sprintf("p %s", fp(t(P))), cex = 0.55, srt = 90, adj = c(0, 0.5), font = ifelse(t(P) < 0.05, 2, 1))
  text(colMeans(bp), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), sprintf("top %d both\n(%d genes)", NS, lengths(SPA)), srt = 30, adj = 1, xpd = NA, cex = 0.7)
  legend("topright", MODES, fill = MC, border = NA, bty = "n", cex = 0.62); ttl("A  H2052, CROP-seq vs bulk: all 8 KOs", "genes variable in both the screen and bulk H2052; exact permutation p")
  ## B
  x <- kA[kA$mode == "c9" & kA$space == "int2000", ]; o <- x$KO[order(-x$own)]; x <- x[match(o, x$KO), ]; par(mar = c(4.6, 4.6, 3.4, 1))
  bp <- barplot(x$own, col = adjustcolor(KO_COL[o], 0.8), border = NA, las = 1, ylim = c(min(x$own, x$median_other) - 0.08, max(x$own) + 0.2), ylab = "Spearman, CROP-seq KO vs bulk H2052 KO", names.arg = rep("", 8), cex.axis = 0.8); abline(h = 0, col = INK)
  segments(bp - 0.4, x$median_other, bp + 0.4, x$median_other, lty = 2, col = "grey25", lwd = 1.3); text(bp, pmax(x$own, 0) + 0.02, sprintf("#%s\np=%s", format(x$rank), fp(x$rank / 8)), cex = 0.62, adj = c(0.5, 0), font = ifelse(x$rank == 1, 2, 1))
  text(bp, par("usr")[3] - 0.02 * diff(par("usr")[3:4]), o, srt = 45, adj = 1, xpd = NA, cex = 0.8, font = 2, col = KO_COL[o]); g <- rA[rA$mode == "c9" & rA$space == "int2000", ]
  ttl("B  H2052 per TF (minus mean of the 9 groups)", sprintf("%d genes in the top 2000 of both; dash = median vs other bulk KOs; # = rank among 8; all KOs p = %s", g$n_genes, fp(g$perm_p)))
  ## C  per line
  x <- pl[pl$mode == "c9", ]; M <- t(sapply(LINES, function(l) x$pair_dmo[x$line == l][match(names(SPB), x$space[x$line == l])])); par(mar = c(6, 4.6, 3.4, 1))
  bp <- barplot(M, beside = TRUE, col = LINE_COL[LINES], border = NA, las = 1, names.arg = rep("", length(SPB)), ylab = "own-KO minus other-KO Spearman, pairs containing the line", cex.axis = 0.8, ylim = range(c(0, M)) + c(-0.005, 0.02)); abline(h = 0, col = INK)
  text(colMeans(bp), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), SPB, srt = 30, adj = 1, xpd = NA, cex = 0.7); legend("topright", SH(LINES), fill = LINE_COL[LINES], border = NA, bty = "n", cex = 0.7)
  ttl("C  which bulk line underperforms?", sprintf("mean over the 3 line pairs that contain the line (9-group centring); lowest average rank: %s", SH(WORST)))
  ## D  leave one out
  x <- rL[rL$mode == "c9", ]; Z <- t(sapply(c("none", LINES), function(l) if (l == "none") rB$perm_z[rB$mode == "c9"][match(names(SPB), rB$space[rB$mode == "c9"])] else x$perm_z[x$dropped == l][match(names(SPB), x$space[x$dropped == l])]))
  Pz <- t(sapply(c("none", LINES), function(l) if (l == "none") rB$perm_p[rB$mode == "c9"][match(names(SPB), rB$space[rB$mode == "c9"])] else x$perm_p[x$dropped == l][match(names(SPB), x$space[x$dropped == l])]))
  par(mar = c(6, 4.6, 3.4, 1)); cc <- c("grey20", LINE_COL[LINES]); bp <- barplot(Z, beside = TRUE, col = cc, border = NA, las = 1, names.arg = rep("", length(SPB)), ylab = "global agreement across lines (permutation z)", cex.axis = 0.8, ylim = range(c(0, Z)) + c(-0.1, 1.1)); abline(h = c(0, qnorm(0.95)), lty = c(1, 3), col = c(INK, MUTED))
  text(bp, pmax(Z, 0) + 0.06, fp(Pz), cex = 0.5, srt = 90, adj = c(0, 0.5), font = ifelse(Pz < 0.05, 2, 1)); text(colMeans(bp), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), SPB, srt = 30, adj = 1, xpd = NA, cex = 0.7)
  legend("topright", c("all 4 lines", paste("without", SH(LINES))), fill = cc, border = NA, bty = "n", cex = 0.62); ttl("D  leave one bulk line out: all KOs", "numbers = permutation p; gene sets recomputed on the lines kept")
  ## E, F  per TF: all 4 lines vs without the worst line
  o2 <- read.csv(file.path(OUT, "spec82_TF_agreement_ranked.csv"))$KO
  for (w in 1:2) { a <- if (w == 1) AG[["c9 signal500"]] else LO[[paste(WORST, "c9", "signal500")]]; k <- a$per_KO; yr <- range(c(unlist(a$cross), a$intra)); yr <- yr + c(-0.02, 0.17) * diff(yr); par(mar = c(4.6, 4.6, 3.4, 1))
    plot(NA, xlim = c(0.5, 8.5), ylim = yr, xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "Spearman between two bulk cell lines", cex.axis = 0.8); abline(h = 0, lty = 2, col = INK)
    for (i in seq_along(o2)) { j <- match(o2[i], k$KO); v <- a$intra[o2[i], ]; boxplot(a$cross[[j]], at = i - 0.2, add = TRUE, boxwex = 0.32, col = "#d5d4cf", border = INK2, outline = FALSE, axes = FALSE)
      boxplot(v, at = i + 0.2, add = TRUE, boxwex = 0.32, col = adjustcolor(KO_COL[[o2[i]]], 0.65), border = INK2, outline = FALSE, axes = FALSE); points(i + 0.2 + seq(-0.1, 0.1, length.out = length(v)), v, pch = 16, cex = 0.75, col = INK2)
      text(i, yr[2], sprintf("p=%s\n%s", fp(k$perm_p[j]), ifelse(k$perm_fdr[j] < 0.01, "**", ifelse(k$perm_fdr[j] < 0.05, "*", ""))), cex = 0.66, adj = c(0.5, 1), font = ifelse(k$perm_p[j] < 0.05, 2, 1)) }
    text(1:8, par("usr")[3] - 0.03 * diff(par("usr")[3:4]), o2, srt = 45, adj = 1, xpd = NA, cex = 0.8, font = 2, col = KO_COL[o2])
    ttl(if (w == 1) "E  per TF, all 4 bulk lines" else sprintf("F  per TF, without %s (3 lines)", SH(WORST)), sprintf("top-500 KO-signal genes, 9-group centring; colour = same KO, grey = other KOs; global p = %s; * / ** FDR < 0.05 / 0.01", fp(a$global[["perm_p"]]))) }
  mtext("9-group centring on genes variable in BOTH data sets.   Top: NCI-H2052 CROP-seq vs bulk, and the bulk line that underperforms.   Bottom: bulk across lines with one line left out", outer = TRUE, side = 3, adj = 0, line = 2.8, font = 2, cex = 0.95)
  mtext("bulk: group mean TMM logCPM per line (S70), variance over the samples of the line; screen: pseudobulk logCPM per guide group, HVGs from Seurat vst on the cells; Spearman; 20,000 KO-label permutations", outer = TRUE, side = 3, adj = 0, line = 1.5, cex = 0.65, col = INK2)
})

## ---------------------------------------------------------------- Cm modules: top-25 genes per module, restricted to the variable genes shared by the data sets compared
# module genes = first 25 genes of each Cm that are in the intersected top-N variable genes (N = 2000 / 3000 / 5000); modules left with < MIN_G genes are dropped. Scores = AddModuleScore (screen: per cell;
# bulk: per sample on TMM logCPM of all S70 samples), mean per group, then the same 9-group centring. H2052: genes variable in screen and bulk H2052; across lines: genes variable in every line compared.
# outputs: spec99_Cm25_genes.csv, spec99_Cm25_H2052.csv, spec99_Cm25_H2052_per_KO.csv, spec99_Cm25_across_lines.csv, spec99_Cm25_across_lines_per_KO.csv, fig12B_Cm25_variable_genes.{png,pdf}
options(Seurat.object.assay.version = "v3"); MIN_G <- 3; cm25 <- lapply({ m <- readRDS(CM_RDS); m[order(as.integer(sub("Cm", "", names(m))))] }, head, 25); CMN <- names(cm25); NM <- c(2000, 3000, 5000)
yB <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(fab$cs[, ids])), log = TRUE, prior.count = 2); yB <- yB[rowSums(fab$cs[, ids]) > 0, ]; yB <- yB - min(yB); sb <- CreateSeuratObject(counts = fab$cs[rownames(yB), ids]); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(yB))
ms <- function(obj, sets) { o <- suppressWarnings(AddModuleScore(obj, features = sets, name = "Cm_", seed = 1)); x <- as.matrix(o@meta.data[, paste0("Cm_", seq_along(sets)), drop = FALSE]); colnames(x) <- names(sets); x }   # cells / samples x modules
bulk_groups <- function(sets, cl) { x <- ms(sb, sets); s <- ids[meta[ids, "Cell_Line"] == cl]; g <- meta[s, "Gene_Control"]; sapply(GR, function(k) colMeans(x[s[g == k], , drop = FALSE])) }
restrict <- function(V) Filter(function(g) length(g) >= MIN_G, lapply(cm25, intersect, V))
gl <- list(); rM <- list(); kM <- list()
for (n in NM) { sets <- restrict(intersect(hvS[seq_len(n)], topv(CL, n, G))); gl[[length(gl) + 1]] <- data.frame(analysis = "H2052 screen vs bulk", N = n, Cm = CMN, n_genes = sapply(cm25, function(g) length(intersect(g, intersect(hvS[seq_len(n)], topv(CL, n, G))))))
  xs <- ms(srt, sets); S9 <- sapply(GR, function(k) colMeans(xs[grp == k, , drop = FALSE])); B9 <- bulk_groups(sets, CL)
  for (md in names(MODES)) { C <- cor(norm9(S9, md), norm9(B9, md), method = "spearman"); dm <- mean(diag(C)); nul <- apply(perms8, 1, function(p) mean(C[cbind(1:8, p)]))
    rM[[length(rM) + 1]] <- data.frame(N = n, n_modules = length(sets), median_genes = median(lengths(sets)), mode = md, diag = dm, off = mean(C[row(C) != col(C)]), perm_p = mean(nul >= dm - 1e-12), top1 = sum(apply(C, 1, which.max) == 1:8))
    kM[[length(kM) + 1]] <- data.frame(N = n, mode = md, KO = KOS, own = diag(C), median_other = sapply(1:8, function(i) median(C[i, -i])), rank = sapply(1:8, function(i) rank(-C[i, ])[i])) } }
rM <- do.call(rbind, rM); rM$dmo <- rM$diag - rM$off; kM <- do.call(rbind, kM); write.csv(rM, file.path(OUT, "spec99_Cm25_H2052.csv"), row.names = FALSE); write.csv(kM, file.path(OUT, "spec99_Cm25_H2052_per_KO.csv"), row.names = FALSE)
cat("\n=== Cm (top 25, variable genes only): H2052 CROP-seq vs bulk ===\n"); print(rM[, c("N", "n_modules", "median_genes", "mode", "diag", "off", "dmo", "perm_p", "top1")], digits = 3, row.names = FALSE)
x <- kM[kM$mode == "c9", ]; cat("\n-- per TF own rho (rank), minus mean of the 9 groups\n"); print(reshape(transform(x[, c("KO", "N")], o = sprintf("%.2f (%s)", x$own, format(x$rank))), idvar = "KO", timevar = "N", direction = "wide"), row.names = FALSE)
rN <- list(); kN <- list(); AGM <- list(); LSETS <- c(list(all4 = LINES), setNames(lapply(LINES, function(l) setdiff(LINES, l)), paste0("no_", SH(LINES))))
for (n in NM) for (ls in names(LSETS)) { lines <- LSETS[[ls]]; V <- Reduce(intersect, lapply(lines, topv, n = n)); sets <- restrict(V); if (ls == "all4") gl[[length(gl) + 1]] <- data.frame(analysis = "bulk 4 lines", N = n, Cm = CMN, n_genes = sapply(cm25, function(g) length(intersect(g, V))))
  if (length(sets) < 5) next; x <- ms(sb, sets); B9 <- lapply(setNames(lines, lines), function(cl) { s <- ids[meta[ids, "Cell_Line"] == cl]; g <- meta[s, "Gene_Control"]; sapply(GR, function(k) colMeans(x[s[g == k], , drop = FALSE])) }); prs <- Filter(function(p) all(p %in% lines), PAIRS)
  for (md in names(MODES)) { A <- array(NA_real_, c(length(sets), 8, 4), dimnames = list(names(sets), KOS, LINES)); for (cl in lines) A[, , cl] <- norm9(B9[[cl]], md); a <- agreement(cl_from_array(A, NULL, centre = FALSE, pairs = prs), permsL, pairs = prs); AGM[[paste(n, ls, md)]] <- a
    rN[[length(rN) + 1]] <- data.frame(N = n, lines = ls, n_modules = length(sets), median_genes = median(lengths(sets)), mode = md, t(a$global)); kN[[length(kN) + 1]] <- cbind(N = n, lines = ls, mode = md, a$per_KO) } }
rN <- do.call(rbind, rN); kN <- do.call(rbind, kN); gl <- do.call(rbind, gl); write.csv(rN, file.path(OUT, "spec99_Cm25_across_lines.csv"), row.names = FALSE); write.csv(kN, file.path(OUT, "spec99_Cm25_across_lines_per_KO.csv"), row.names = FALSE); write.csv(gl, file.path(OUT, "spec99_Cm25_genes.csv"), row.names = FALSE)
cat("\n=== Cm (top 25, variable genes only): bulk across lines ===\n"); print(rN[, c("N", "lines", "n_modules", "median_genes", "mode", "diag_mean", "off_mean", "perm_p", "perm_z")], digits = 3, row.names = FALSE)
x <- kN[kN$mode == "c9" & kN$lines %in% c("all4", paste0("no_", SH(WORST))), ]; cat("\n-- per TF perm p (mean same-KO rho), c9\n"); print(reshape(transform(x[, c("KO", "N", "lines")], id = paste(N, lines), v = sprintf("%.3f (r %.2f)", x$perm_p, x$mean_intra))[, c("KO", "id", "v")], idvar = "KO", timevar = "id", direction = "wide"), row.names = FALSE)
cat("\n-- genes kept per module (of 25)\n"); print(reshape(transform(gl, id = paste(analysis, N))[, c("Cm", "id", "n_genes")], idvar = "Cm", timevar = "id", direction = "wide"), row.names = FALSE)

save_fig("fig12B_Cm25_variable_genes", 15, 9.6, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4.4, 0.4), family = "sans"); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.8); mtext(b, side = 3, adj = 0, line = 0.4, cex = 0.58, col = INK2) }
  par(mar = c(5, 4.6, 3.4, 1)); gm <- function(col) sapply(names(MODES), function(m) rM[[col]][rM$mode == m][match(NM, rM$N[rM$mode == m])]); D <- gm("dmo"); P <- gm("perm_p"); u <- rM[rM$mode == "c9", ][match(NM, rM$N[rM$mode == "c9"]), ]
  bp <- barplot(t(D), beside = TRUE, col = MC, border = NA, las = 1, names.arg = sprintf("top %d both\n%d Cm, median %g genes", NM, u$n_modules, u$median_genes), cex.names = 0.62, ylab = "own-KO minus other-KO Spearman", cex.axis = 0.8, ylim = range(c(0, D)) + c(-0.03, 0.09)); abline(h = 0, col = INK)
  text(bp, pmax(t(D), 0) + 0.006, sprintf("p %s", fp(t(P))), cex = 0.55, srt = 90, adj = c(0, 0.5), font = ifelse(t(P) < 0.05, 2, 1)); legend("topright", MODES, fill = MC, border = NA, bty = "n", cex = 0.6); ttl("A  H2052, CROP-seq vs bulk: all 8 KOs, Cm modules", "exact permutation p of the bulk KO labels")
  for (n in c(3000, 5000)) { x <- kM[kM$mode == "c9" & kM$N == n, ]; o <- x$KO[order(-x$own)]; x <- x[match(o, x$KO), ]; par(mar = c(4.6, 4.6, 3.4, 1))
    bp <- barplot(x$own, col = adjustcolor(KO_COL[o], 0.8), border = NA, las = 1, ylim = c(min(x$own, x$median_other) - 0.1, max(x$own) + 0.25), ylab = "Spearman over Cm, CROP-seq KO vs bulk H2052 KO", names.arg = rep("", 8), cex.axis = 0.8); abline(h = 0, col = INK)
    segments(bp - 0.4, x$median_other, bp + 0.4, x$median_other, lty = 2, col = "grey25", lwd = 1.3); text(bp, pmax(x$own, 0) + 0.02, sprintf("#%s", format(x$rank)), cex = 0.66, adj = c(0.5, 0), font = ifelse(x$rank == 1, 2, 1))
    text(bp, par("usr")[3] - 0.02 * diff(par("usr")[3:4]), o, srt = 45, adj = 1, xpd = NA, cex = 0.8, font = 2, col = KO_COL[o]); g <- rM[rM$mode == "c9" & rM$N == n, ]
    ttl(sprintf("%s  H2052 per TF, module genes in top %d of both", if (n == 3000) "B" else "C", n), sprintf("%d modules; 9-group centring; dash = median vs other bulk KOs; # = rank among 8; all KOs p = %s", g$n_modules, fp(g$perm_p))) }
  x <- rN[rN$mode == "c9", ]; lsn <- names(LSETS); Z <- sapply(NM, function(n) x$perm_z[x$N == n][match(lsn, x$lines[x$N == n])]); Pz <- sapply(NM, function(n) x$perm_p[x$N == n][match(lsn, x$lines[x$N == n])]); cc <- c("grey20", LINE_COL[LINES]); par(mar = c(5, 4.6, 3.4, 1))
  u <- x[x$lines == "all4", ][match(NM, x$N[x$lines == "all4"]), ]; bp <- barplot(Z, beside = TRUE, col = cc, border = NA, las = 1, names.arg = sprintf("top %d all lines\n%s Cm (4 lines)", NM, ifelse(is.na(u$n_modules), "<5", u$n_modules)), cex.names = 0.62, ylab = "global agreement across lines (permutation z)", cex.axis = 0.8, ylim = range(c(0, Z), na.rm = TRUE) + c(-0.1, 1.1)); abline(h = c(0, qnorm(0.95)), lty = c(1, 3), col = c(INK, MUTED))
  text(bp, pmax(Z, 0) + 0.06, ifelse(is.na(Pz), "", fp(Pz)), cex = 0.5, srt = 90, adj = c(0, 0.5), font = ifelse(!is.na(Pz) & Pz < 0.05, 2, 1)); legend("topright", c("all 4 lines", paste("without", SH(LINES))), fill = cc, border = NA, bty = "n", cex = 0.6); ttl("D  bulk across lines: all KOs, Cm modules", "9-group centring; numbers = permutation p")
  for (ls in c("all4", paste0("no_", SH(WORST)))) { a <- AGM[[paste(5000, ls, "c9")]]; k <- a$per_KO; o <- k$KO[order(k$perm_p)]; k <- k[match(o, k$KO), ]; par(mar = c(4.6, 4.6, 3.4, 1))
    bp <- barplot(k$mean_intra, col = adjustcolor(KO_COL[o], 0.8), border = NA, las = 1, ylim = range(c(0, k$mean_intra, k$median_cross)) + c(-0.08, 0.15), ylab = "mean same-KO Spearman over Cm, between bulk lines", names.arg = rep("", 8), cex.axis = 0.8); abline(h = 0, col = INK)
    segments(bp - 0.4, k$median_cross, bp + 0.4, k$median_cross, lty = 2, col = "grey25", lwd = 1.3); text(bp, pmax(k$mean_intra, 0) + 0.015, sprintf("p=%s%s", fp(k$perm_p), ifelse(k$perm_fdr < 0.05, "*", "")), cex = 0.62, adj = c(0.5, 0), font = ifelse(k$perm_p < 0.05, 2, 1))
    text(bp, par("usr")[3] - 0.02 * diff(par("usr")[3:4]), o, srt = 45, adj = 1, xpd = NA, cex = 0.8, font = 2, col = KO_COL[o])
    ttl(sprintf("%s  per TF across bulk lines, %s", if (ls == "all4") "E" else "F", if (ls == "all4") "all 4 lines" else paste("without", SH(WORST))), sprintf("module genes in top 5000 of every line; dash = median vs other KOs; global p = %s; * FDR < 0.05", fp(a$global[["perm_p"]]))) }
  mtext("Cm modules scored on the top-25 module genes that are variable in the data sets compared (AddModuleScore), profiles centred on the mean of the 9 groups", outer = TRUE, side = 3, adj = 0, line = 2.8, font = 2, cex = 0.95)
  mtext(sprintf("modules with fewer than %d genes left are dropped; top: NCI-H2052 CROP-seq vs bulk; bottom: bulk across lines (S70), 20,000 KO-label permutations", MIN_G), outer = TRUE, side = 3, adj = 0, line = 1.5, cex = 0.65, col = INK2)
})
