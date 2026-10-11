# Spec 127: does correcting the bulk for technical batch structure improve the match with the CROP-seq?
# Everything is KO vs the other KOs (KO samples / KO cells only, no NTC), corrected labels (H2052_NTC_R3 = BPTF), genes in the top 2000 variable of both data sets, Spearman, exact permutation of the bulk KO labels.
# corrections of the bulk (applied per cell line before the KO profiles are built):
#   none      TMM logCPM as is
#   decont    carry-over removed at the count level: counts - RHO x (sample depth) x (pooled expression fraction of the other samples of its i5 group, all 106 samples of the run), RHO = 0.7% (spec125)
#   tech      regression on technical covariates over the KO samples of the line: log10 depth, % assigned, % mitochondrial, log10 liver-gene CPM (carry-over marker)
#   cc        regression on the cell-cycle score
#   cc+tech   both
#   i5block   each KO centred on the KOs of its own i5 index group within the line (removes an additive i5 shift; the KOs of a line are split over 2-3 i5 groups)
#   ruv1-3    replicate-based: gene-space directions of the differences between replicates of the same KO (top k singular vectors) removed from all samples of the line
#   decont+cc+ruv2   combined
# screen: as is, or cell cycle regressed (S.Score + G2M.Score over the KO cells), to pair with the bulk corrections that include cc
# comparisons: (1) the CROP-seq vs each bulk line; (2) the CROP-seq vs one pooled bulk profile per KO from a joint model over the 4 lines:
#   pooled_mean   mean over the lines of the per-line KO-centred profiles          pooled_lm      limma  ~ line + KO (KO effect common to the lines)
#   pooled_block  limma ~ line:i5 block + KO (batch = i5 group within line)         pooled_block_cc  the same + cell-cycle score
# i7: unique within an i5 group, so it cannot be a batch term within a line; reported: are samples that share an i7 sequence more alike than others (residuals after line x KO)?
# needs: spec72 (S70 lists), spec71 audit, sample sheet, agree_common.R::load_bulk_S91, spec_common.R
# outputs: spec127_per_line.csv, spec127_pooled.csv, spec127_per_KO.csv, fig44_batch_corrections_vs_cropseq.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages({ library(Seurat); library(limma) }); options(width = 230); RHO <- as.numeric(Sys.getenv("RHO", unset = "0.007")); SH <- function(l) sub("NCI-|-211H", "", l)
fp <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2)))
pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); keep <- pr$genes; B <- load_bulk_S91(pr$samples, extra_kos = FALSE); m <- B$meta; q <- load_audit()
ks <- B$ids[m[B$ids, "Gene_Control"] %in% KOS]; ss <- readLines(file.path(PROC, "sample_sheet/SampleSheet.csv")); ss <- ss[(grep("^\\[BCLConvert_Data\\]", ss) + 2):length(ss)]; ss <- do.call(rbind, strsplit(ss[nzchar(ss)], ",")); i7 <- setNames(ss[, 3], ss[, 2]); i5 <- setNames(ss[, 4], ss[, 2])
cat("KO samples:", length(ks), "\ni5 group x line:\n"); print(table(i5[ks], m[ks, "Cell_Line"])); cat("KOs per i5 group within each line:\n"); for (cl in LINES) { s <- ks[m[ks, "Cell_Line"] == cl]; for (g in unique(i5[s])) cat(sprintf("  %-9s %s: %s\n", SH(cl), g, paste(sort(unique(m[s[i5[s] == g], "Gene_Control"])), collapse = " "))) }

## ---------------------------------------------------------------- carry-over removal at the count level (all samples of the run define the pools)
ct <- read.delim(file.path(PROC, "counts/gene_counts_matrix_stranded_reverse.tsv"), row.names = 1, check.names = FALSE); sy <- read.delim(file.path(PROC, "counts/gene_id_to_symbol.tsv"), header = FALSE, col.names = c("id", "symbol")); sy <- sy[sy$id %in% rownames(ct) & !is.na(sy$symbol) & sy$symbol != "", ]
csA <- rowsum(as.matrix(ct[sy$id, ]), group = sy$symbol)[rownames(B$cs), ]; dep <- colSums(csA); grp_tot <- sapply(unique(i5[colnames(csA)]), function(g) rowSums(csA[, i5[colnames(csA)] == g, drop = FALSE])); grp_dep <- colSums(grp_tot)
dec <- sapply(ks, function(s) { g <- i5[[s]]; pool <- (grp_tot[, g] - csA[, s]) / (grp_dep[[g]] - dep[[s]]); pmax(csA[, s] - RHO * dep[[s]] * pool, 0) }); rownames(dec) <- rownames(csA)
LIVG <- intersect(c("ALB", "APOA1", "APOA2", "APOB", "FGB", "FGG", "TF", "TTR", "SERPINA1", "AFP"), rownames(csA)); liv <- log10(1e6 * colSums(csA[LIVG, ks]) / dep[ks] + 1); liv_dec <- log10(1e6 * colSums(dec[LIVG, ]) / colSums(dec) + 1)
cat(sprintf("\ncarry-over removal (RHO = %.2f%%): reads removed per sample median %.2f%%; liver-gene CPM median %.1f -> %.1f in the two liver-containing i5 groups\n", 100 * RHO, 100 * median(1 - colSums(dec) / dep[ks]), median((10^liv - 1)[i5[ks] %in% c("ATATTCAC", "GCGCCTGT")]), median((10^liv_dec - 1)[i5[ks] %in% c("ATATTCAC", "GCGCCTGT")])))

## ---------------------------------------------------------------- screen profiles (KO cells only), as is and cell cycle regressed
srt <- load_screen(); grp <- as.character(srt$merged_call); ck <- which(grp %in% KOS); Y <- GetAssayData(srt, layer = "data")[, ck]; g <- grp[ck]; md <- srt@meta.data[ck, ]; Y <- Y[Matrix::rowSums(Y > 0) >= 0.01 * ncol(Y), ]; Gs <- intersect(keep, rownames(Y)); Y <- Y[Gs, ]
Xc <- scale(cbind(md$S.Score, md$G2M.Score), scale = FALSE); Bc <- as.matrix(Y %*% Xc) %*% solve(crossprod(Xc)); gm <- sapply(KOS, function(k) Matrix::rowMeans(Y[, g == k, drop = FALSE])); gx <- sapply(KOS, function(k) colMeans(Xc[g == k, , drop = FALSE]))
SCR <- list(raw = gm - rowMeans(gm), cc = { r <- gm - Bc %*% gx; r - rowMeans(r) }); hv <- intersect(VariableFeatures(FindVariableFeatures(srt[, ck], selection.method = "vst", nfeatures = nrow(srt), verbose = FALSE)), Gs)

## ---------------------------------------------------------------- bulk corrections per line
reg <- function(y, X) { X <- scale(X, scale = FALSE); X <- X[, apply(X, 2, sd) > 0, drop = FALSE]; y - y %*% X %*% solve(crossprod(X)) %*% t(X) }                    # unprotected regression over the KO samples of the line
ruv <- function(y, gl, k) { R <- y; for (kk in unique(gl)) { w <- gl == kk; R[, w] <- y[, w, drop = FALSE] - rowMeans(y[, w, drop = FALSE]) }; R <- R[, gl %in% names(which(table(gl) > 1)), drop = FALSE]
  U <- svd(R, nu = k, nv = 0)$u; yc <- y - rowMeans(y); y - U %*% crossprod(U, yc) }                                                                              # remove the replicate-noise directions from all samples
CORR <- c("none", "decont", "tech", "cc", "cc+tech", "i5block", "ruv1", "ruv2", "ruv3", "decont+cc+ruv2"); USES_CC <- c(cc = TRUE, `cc+tech` = TRUE, `decont+cc+ruv2` = TRUE)
MAT <- list(); for (cl in LINES) { s <- ks[m[ks, "Cell_Line"] == cl]; gl <- m[s, "Gene_Control"]; y0 <- normalise(B$cs[keep, s], "tmm")$y[Gs, ]; yd <- normalise(dec[keep, s], "tmm")$y[Gs, ]
  tech <- cbind(depth = log10(q[s, "assigned_M"]), pa = q[s, "pct_assigned"], mito = q[s, "pct_mito"], liver = liv[s]); cc <- cbind(cc = B$cc[s])
  MAT[[cl]] <- list(samples = s, gl = gl, i5 = i5[s], y = list(none = y0, decont = yd, tech = reg(y0, tech), cc = reg(y0, cc), `cc+tech` = reg(y0, cbind(cc, tech)), i5block = y0, ruv1 = ruv(y0, gl, 1), ruv2 = ruv(y0, gl, 2), ruv3 = ruv(y0, gl, 3), `decont+cc+ruv2` = ruv(reg(yd, cc), gl, 2))) }
prof <- function(cl, cr) { M <- MAT[[cl]]; y <- M$y[[cr]]; p <- sapply(KOS, function(k) rowMeans(y[, M$gl == k, drop = FALSE]))
  if (cr == "i5block") { blk <- sapply(KOS, function(k) names(which.max(table(M$i5[M$gl == k])))); for (b in unique(blk)) { w <- blk == b; p[, w] <- if (sum(w) > 1) p[, w, drop = FALSE] - rowMeans(p[, w, drop = FALSE]) else NA }; p } else p - rowMeans(p) }
perms8 <- local({ pf <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(pf(v[-i]), function(r) c(v[i], r)))); do.call(rbind, pf(1:8)) })
cmp <- function(S, Bp, f) { ok <- which(colSums(!is.na(Bp)) > 0); C <- cor(S[f, ok], Bp[f, ok], method = "spearman"); n <- length(ok); P <- if (n == 8) perms8 else unique(t(apply(perms8, 1, function(p) p[p <= n])))[, seq_len(n), drop = FALSE]
  dm <- mean(diag(C)); nul <- apply(P, 1, function(p) mean(C[cbind(seq_len(n), p)])); list(C = C, kos = colnames(Bp)[ok], dmo = dm - mean(C[row(C) != col(C)]), p = mean(nul >= dm - 1e-12), own = diag(C), rank = sapply(seq_len(n), function(i) rank(-C[i, ])[i])) }
vsel <- function(v) { vb <- Gs[order(-v)]; intersect(hv[1:2000], vb[1:2000]) }

## ---------------------------------------------------------------- 1  per line
R1 <- list(); PK <- list()
for (cl in LINES) for (cr in CORR) { sm <- if (isTRUE(USES_CC[cr])) "cc" else "raw"; Bp <- prof(cl, cr); f <- vsel(matrixStats::rowVars(MAT[[cl]]$y[[cr]])); r <- cmp(SCR[[sm]], Bp, f)
  R1[[length(R1) + 1]] <- data.frame(line = cl, correction = cr, screen = sm, n_genes = length(f), n_KO = length(r$kos), dmo = r$dmo, p = r$p, n_rank1 = sum(r$rank == 1)); PK[[length(PK) + 1]] <- data.frame(target = cl, correction = cr, KO = r$kos, own = r$own, rank = r$rank) }
R1 <- do.call(rbind, R1); write.csv(R1, file.path(OUT, "spec127_per_line.csv"), row.names = FALSE); cat("\n=== 1  CROP-seq vs each bulk line: own - other Spearman (p) ===\n"); R1$v <- sprintf("%+.3f (%s)", R1$dmo, fp(R1$p)); print(reshape(R1[, c("correction", "line", "v")], idvar = "correction", timevar = "line", direction = "wide"), row.names = FALSE)
cat("mean over the 4 lines:\n"); print(round(sort(tapply(R1$dmo, R1$correction, mean), decreasing = TRUE), 3))

## ---------------------------------------------------------------- 2  pooled bulk profile per KO
line <- factor(m[ks, "Cell_Line"], levels = LINES); ko <- factor(m[ks, "Gene_Control"], levels = KOS); blk <- factor(paste(line, i5[ks])); yall <- do.call(cbind, lapply(LINES, function(cl) MAT[[cl]]$y$none))[, ks]; ydec <- do.call(cbind, lapply(LINES, function(cl) MAT[[cl]]$y$decont))[, ks]
kocoef <- function(y, design) { ne <- nonEstimable(design); if (!is.null(ne)) cat("  not estimable:", ne, "\n"); fit <- lmFit(y, design); cf <- fit$coefficients; K <- cbind(0, cf[, grep("^ko", colnames(cf)), drop = FALSE]); colnames(K) <- KOS; K - rowMeans(K, na.rm = TRUE) }   # KO effects vs the first KO, then centred over the 8 KOs
ccv <- B$cc[ks] - ave(B$cc[ks], line); cat("\npooled designs (rank / columns): line + KO", qr(model.matrix(~ line + ko))$rank, "/", ncol(model.matrix(~ line + ko)), "| block + KO", qr(model.matrix(~ blk + ko))$rank, "/", ncol(model.matrix(~ blk + ko)), "\n")
POOL <- list(pooled_mean = Reduce(`+`, lapply(LINES, function(cl) prof(cl, "none"))) / 4, pooled_lm = kocoef(yall, model.matrix(~ line + ko)), pooled_block = kocoef(yall, model.matrix(~ blk + ko)), pooled_block_cc = kocoef(yall, model.matrix(~ blk + ko + ccv)),
             pooled_block_cc_decont = kocoef(ydec, model.matrix(~ blk + ko + ccv)), pooled_3lines_block_cc = { w <- line != "MSTO-211H"; kocoef(yall[, w], model.matrix(~ droplevels(blk[w]) + ko[w] + ccv[w]) -> dd) }, pooled_mean_ruv2 = Reduce(`+`, lapply(LINES, function(cl) prof(cl, "ruv2"))) / 4)
vpool <- matrixStats::rowVars(yall - t(apply(yall, 1, function(v) ave(v, line)))); f <- vsel(vpool); R2 <- list()
for (pn in names(POOL)) for (sm in c("raw", "cc")) { r <- cmp(SCR[[sm]], POOL[[pn]], f); R2[[length(R2) + 1]] <- data.frame(pooled = pn, screen = sm, n_genes = length(f), dmo = r$dmo, p = r$p, n_rank1 = sum(r$rank == 1), n_rank12 = sum(r$rank <= 2)); PK[[length(PK) + 1]] <- data.frame(target = paste0(pn, " / screen ", sm), correction = pn, KO = r$kos, own = r$own, rank = r$rank) }
R2 <- do.call(rbind, R2); PK <- do.call(rbind, PK); write.csv(R2, file.path(OUT, "spec127_pooled.csv"), row.names = FALSE); write.csv(PK, file.path(OUT, "spec127_per_KO.csv"), row.names = FALSE)
cat("\n=== 2  CROP-seq vs the pooled bulk profile (4 lines together) ===\n"); print(R2, digits = 3, row.names = FALSE)
for (tg in c("pooled_mean / screen raw", "pooled_block_cc / screen cc", "pooled_3lines_block_cc / screen cc")) { x <- PK[PK$target == tg, ]; cat(sprintf("-- per KO, %s: %s\n", tg, paste(sprintf("%s %+.2f(#%s)", x$KO, x$own, format(x$rank)), collapse = "  "))) }

## ---------------------------------------------------------------- i7: are samples sharing an i7 sequence more alike?
res <- yall; for (cl in LINES) for (k in KOS) { w <- line == cl & ko == k; if (sum(w) > 1) res[, w] <- yall[, w] - rowMeans(yall[, w, drop = FALSE]) else res[, w] <- NA }; okc <- colSums(is.na(res)) == 0; tv <- order(-matrixStats::rowVars(res[, okc]))[1:3000]; Cr <- cor(res[tv, okc], method = "spearman"); s2 <- ks[okc]
ut <- which(upper.tri(Cr), arr.ind = TRUE); same_i7 <- i7[s2][ut[, 1]] == i7[s2][ut[, 2]]; same_i5 <- i5[s2][ut[, 1]] == i5[s2][ut[, 2]]; mate <- paste(line, ko)[okc][ut[, 1]] == paste(line, ko)[okc][ut[, 2]]; v <- Cr[ut]
cat(sprintf("\nreplicate-residual similarity (Spearman, pairs that are not replicate mates): same i7 %.3f (n = %d) vs different i7 %.3f, Wilcoxon p %s | same i5 group %.3f (n = %d) vs different %.3f, p %s\n", mean(v[same_i7 & !mate]), sum(same_i7 & !mate), mean(v[!same_i7 & !mate]), fp(wilcox.test(v[same_i7 & !mate], v[!same_i7 & !mate])$p.value),
            mean(v[same_i5 & !mate]), sum(same_i5 & !mate), mean(v[!same_i5 & !mate]), fp(wilcox.test(v[same_i5 & !mate], v[!same_i5 & !mate])$p.value)))

## ---------------------------------------------------------------- figure
pal <- colorRampPalette(c("#2a78d6", "#f7f6f2", "#eb6834"))(101)
save_fig("fig44_batch_corrections_vs_cropseq", 16, 9.6, function() {
  layout(matrix(1:3, nrow = 1), widths = c(1.15, 1, 1.25)); par(oma = c(0, 0, 4, 0.4), family = "sans"); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.8); mtext(b, side = 3, adj = 0, line = 0.5, cex = 0.56, col = INK2) }
  M <- sapply(LINES, function(cl) R1$dmo[R1$line == cl][match(CORR, R1$correction[R1$line == cl])]); P <- sapply(LINES, function(cl) R1$p[R1$line == cl][match(CORR, R1$correction[R1$line == cl])]); M <- cbind(M, mean = rowMeans(M)); n <- length(CORR); zl <- 0.15
  par(mar = c(4, 9, 3.6, 1)); image(1:5, 1:n, t(pmax(pmin(M, zl), -zl)[n:1, ]), col = pal, zlim = c(-zl, zl), axes = FALSE, xlab = "", ylab = ""); text(rep(1:4, each = n), rep(n:1, 4), sprintf("%+.3f\np %s", as.vector(M[, 1:4]), fp(as.vector(P))), cex = 0.66, font = ifelse(as.vector(P) < 0.05, 2, 1))
  text(5, n:1, sprintf("%+.3f", M[, 5]), cex = 0.72, font = 2); axis(1, at = 1:5, labels = c(SH(LINES), "mean"), tick = FALSE, cex.axis = 0.85); axis(2, at = n:1, labels = CORR, tick = FALSE, las = 1, cex.axis = 0.8); abline(v = 4.5, lwd = 2, col = "white")
  ttl("A  CROP-seq vs each bulk line, by correction of the bulk", "own-KO minus other-KO Spearman and exact permutation p (8 KOs; i5block: KOs with a block-mate)")
  x <- R2; x$lab <- paste0(x$pooled, "\nscreen ", ifelse(x$screen == "cc", "cc regressed", "as is")); par(mar = c(4, 12, 3.6, 3)); bp <- barplot(rev(x$dmo), horiz = TRUE, col = rev(ifelse(x$screen == "cc", "#2a78d6", "grey55")), border = NA, las = 1, names.arg = rev(x$lab), cex.names = 0.6, xlab = "own-KO minus other-KO Spearman", cex.axis = 0.8, xlim = range(c(0, x$dmo)) + c(-0.02, 0.06)); abline(v = 0, col = INK)
  text(pmax(rev(x$dmo), 0) + 0.004, bp, sprintf("p %s; %d KO #1", fp(rev(x$p)), rev(x$n_rank1)), adj = 0, cex = 0.62, font = ifelse(rev(x$p) < 0.05, 2, 1), xpd = NA); ttl("B  CROP-seq vs one pooled bulk profile per KO", "joint models over the 4 lines; block = i5 group within line")
  tg <- c(paste0("pooled_mean / screen raw"), "pooled_block_cc / screen cc", "pooled_3lines_block_cc / screen cc", "pooled_mean_ruv2 / screen raw"); Mo <- sapply(tg, function(t) PK$own[PK$target == t][match(KOS, PK$KO[PK$target == t])]); Ro <- sapply(tg, function(t) PK$rank[PK$target == t][match(KOS, PK$KO[PK$target == t])])
  Ml <- sapply(LINES, function(cl) PK$own[PK$target == cl & PK$correction == "none"][match(KOS, PK$KO[PK$target == cl & PK$correction == "none"])]); Rl <- sapply(LINES, function(cl) PK$rank[PK$target == cl & PK$correction == "none"][match(KOS, PK$KO[PK$target == cl & PK$correction == "none"])])
  MM <- cbind(Ml, Mo); RR <- cbind(Rl, Ro); nc <- ncol(MM); zl2 <- 0.6; par(mar = c(8.5, 6, 3.6, 1)); image(seq_len(nc), 1:8, t(pmax(pmin(MM, zl2), -zl2)[8:1, ]), col = pal, zlim = c(-zl2, zl2), axes = FALSE, xlab = "", ylab = ""); text(rep(seq_len(nc), each = 8), rep(8:1, nc), sprintf("%+.2f\n#%s", as.vector(MM), format(as.vector(RR))), cex = 0.6, font = ifelse(as.vector(RR) <= 2, 2, 1))
  axis(1, at = seq_len(nc), labels = c(paste(SH(LINES), "(none)"), "pooled mean", "pooled block+cc", "3 lines block+cc", "pooled mean ruv2"), tick = FALSE, las = 2, cex.axis = 0.72); text(0.4, 8:1, KOS, adj = 1, xpd = NA, cex = 0.85, font = 2, col = KO_COL[KOS]); abline(v = 4.5, lwd = 2, col = "white")
  ttl("C  per TF: own-KO Spearman with the CROP-seq", "# = rank of the own KO among the 8 bulk KOs; left: each line uncorrected; right: pooled profiles")
  mtext("Batch corrections of the bulk and the match with the CROP-seq", outer = TRUE, side = 3, adj = 0, line = 2.4, font = 2, cex = 1.05)
  mtext("KO vs the other KOs, no NTC; genes in the top 2000 variable of both; decont = carry-over within the i5 group removed; tech = depth, % assigned, % mito, liver-gene signal; ruv = replicate-noise directions removed", outer = TRUE, side = 3, adj = 0, line = 1.2, cex = 0.66, col = INK2)
})

## ---------------------------------------------------------------- grid: corrections x number of top genes x editing tier x replicate choice
# top-N variable genes in both data sets (N = 500 ... 5000, or all shared genes); KO tiers by the best-guide indel % of the KO in the line (spec124_indel_tiers.csv): all | >= 5% | >= 10%;
# replicates: all, or only the best-edited replicate of each KO (spec120 evidence table). Pooled profile: mean over the lines in which the KO passes the tier.
# outputs: spec127_grid.csv, fig44b_batch_corrections_grid.{png,pdf}
IND <- as.matrix(read.csv(file.path(OUT, "spec124_indel_tiers.csv"), row.names = 1, check.names = FALSE))[KOS, LINES]; ev <- read.csv(file.path(OUT, "spec120_sample_cut_evidence.csv"), stringsAsFactors = FALSE); evp <- setNames(ev$best_pct, ev$sample)
THR <- c(all = -Inf, `indel>=5` = 5, `indel>=10` = 10); pass <- function(cl, tn) if (tn == "all") KOS else KOS[!is.na(IND[, cl]) & IND[, cl] >= THR[[tn]]]; NS <- c(500, 1000, 2000, 3000, 5000, Inf); GC <- c("none", "decont", "tech", "cc", "cc+tech", "ruv2", "decont+cc+ruv2")
profg <- function(cl, cr, best) { M <- MAT[[cl]]; y <- M$y[[cr]]; pick <- function(k) { w <- which(M$gl == k); if (best && length(w) > 1) { e <- evp[M$samples[w]]; e[is.na(e)] <- -1; w <- w[which.max(e)] }; w }; p <- sapply(KOS, function(k) rowMeans(y[, pick(k), drop = FALSE])); p - rowMeans(p) }
PMl <- lapply(1:8, function(n) { pf <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(pf(v[-i]), function(r) c(v[i], r)))); do.call(rbind, pf(seq_len(n))) })
cmpk <- function(S, Bp, f, kk) { C <- cor(S[f, kk], Bp[f, kk], method = "spearman"); n <- length(kk); dm <- mean(diag(C)); nul <- apply(PMl[[n]], 1, function(p) mean(C[cbind(seq_len(n), p)])); c(dmo = dm - mean(C[row(C) != col(C)]), p = mean(nul >= dm - 1e-12), n_rank1 = sum(apply(C, 1, which.max) == seq_len(n))) }
gsel <- function(v, N) if (is.infinite(N)) Gs else { vb <- Gs[order(-v)]; intersect(hv[seq_len(N)], vb[seq_len(N)]) }; GR <- list()
for (cr in GC) { sm <- if (isTRUE(USES_CC[cr])) "cc" else "raw"; for (best in c(FALSE, TRUE)) { PL <- lapply(setNames(LINES, LINES), function(cl) profg(cl, cr, best)); VL <- lapply(setNames(LINES, LINES), function(cl) matrixStats::rowVars(MAT[[cl]]$y[[cr]]))
  for (tn in names(THR)) { # pooled profile: per KO, mean over the lines where it passes the tier
    PP <- sapply(KOS, function(k) { ls <- LINES[sapply(LINES, function(cl) k %in% pass(cl, tn))]; if (!length(ls)) rep(NA_real_, length(Gs)) else rowMeans(sapply(ls, function(cl) PL[[cl]][, k])) }); rownames(PP) <- Gs; kp <- KOS[colSums(!is.na(PP)) > 0]
    for (N in NS) { for (cl in LINES) { kk <- pass(cl, tn); if (length(kk) < 4) next; f <- gsel(VL[[cl]], N); r <- cmpk(SCR[[sm]], PL[[cl]], f, kk); GR[[length(GR) + 1]] <- data.frame(target = cl, correction = cr, replicates = ifelse(best, "best-edited", "all"), tier = tn, N = N, n_genes = length(f), n_KO = length(kk), t(r)) }
      f <- gsel(Reduce(`+`, VL) / 4, N); r <- cmpk(SCR[[sm]], PP, f, kp); GR[[length(GR) + 1]] <- data.frame(target = "pooled 4 lines", correction = cr, replicates = ifelse(best, "best-edited", "all"), tier = tn, N = N, n_genes = length(f), n_KO = length(kp), t(r)) } } } }
GR <- do.call(rbind, GR); write.csv(GR, file.path(OUT, "spec127_grid.csv"), row.names = FALSE); nb <- sum(GR$p < 0.05)
cat(sprintf("\n=== grid: %d combinations; p < 0.05 in %d (%.1f%%; 5%% expected by chance); p < 0.01 in %d ===\n", nrow(GR), nb, 100 * nb / nrow(GR), sum(GR$p < 0.01))); cat("p < 0.05 by target:\n"); print(table(GR$target, GR$p < 0.05))
cat("\nmean own - other by number of genes:\n"); print(round(tapply(GR$dmo, list(GR$target, ifelse(is.infinite(GR$N), "all", GR$N)), mean)[, c("500", "1000", "2000", "3000", "5000", "all")], 3))
cat("\nby tier x replicates:\n"); print(round(tapply(GR$dmo, list(GR$target, paste(GR$tier, GR$replicates)), mean), 3)); cat("\nby correction:\n"); print(round(tapply(GR$dmo, list(GR$target, GR$correction), mean)[, GC], 3))
cat("\nbest 15 combinations:\n"); print(head(GR[order(GR$p, -GR$dmo), c("target", "correction", "replicates", "tier", "N", "n_genes", "n_KO", "dmo", "p", "n_rank1")], 15), digits = 3, row.names = FALSE)
cat("\nbest combination per target:\n"); print(do.call(rbind, lapply(split(GR, GR$target), function(d) d[order(d$p, -d$dmo)[1], c("target", "correction", "replicates", "tier", "N", "n_KO", "dmo", "p")])), digits = 3, row.names = FALSE)
save_fig("fig44b_batch_corrections_grid", 16, 5.8, function() {
  layout(matrix(1:3, nrow = 1)); par(oma = c(0, 0, 3, 0.4), family = "sans", mar = c(4.6, 4.8, 3.6, 1)); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.8); mtext(b, side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2) }; TG <- c(LINES, "pooled 4 lines"); tc <- c(LINE_COL[LINES], `pooled 4 lines` = "black")
  for (cr in c("none", "decont+cc+ruv2")) { d <- GR[GR$correction == cr, ]; plot(NA, xlim = c(1, 6), ylim = range(GR$dmo[GR$correction %in% c("none", "decont+cc+ruv2")]), xaxt = "n", bty = "n", las = 1, xlab = "top-N variable genes in both data sets", ylab = "own-KO minus other-KO Spearman", cex.axis = 0.8); abline(h = 0, lty = 2, col = INK); axis(1, at = 1:6, labels = c(500, 1000, 2000, 3000, 5000, "all"), cex.axis = 0.8)
    for (tg in TG) for (sty in list(c("all", "all", 1, 16), c("indel>=10", "best-edited", 2, 17))) { e <- d[d$target == tg & d$tier == sty[1] & d$replicates == sty[2], ]; e <- e[match(NS, e$N), ]; lines(1:6, e$dmo, col = tc[[tg]], lty = as.integer(sty[3]), lwd = 1.6); points(1:6, e$dmo, col = tc[[tg]], pch = ifelse(!is.na(e$p) & e$p < 0.05, as.integer(sty[4]), NA), cex = 1.5) }
    if (cr == "none") { legend("bottomright", sub("NCI-|-211H", "", TG), col = tc, lwd = 2, bty = "n", cex = 0.65); legend("bottomleft", c("all KOs, all replicates", "indel >= 10%, best-edited replicate", "symbol = p < 0.05"), lty = c(1, 2, NA), pch = c(NA, NA, 16), bty = "n", cex = 0.62) }
    ttl(if (cr == "none") "A  no correction" else "B  carry-over + cell cycle + replicate noise removed", "by number of top genes, target and KO / replicate selection") }
  hist(GR$p, breaks = seq(0, 1, 0.05), col = "grey75", border = "white", las = 1, main = "", xlab = "permutation p", ylab = "combinations", cex.axis = 0.8); abline(h = nrow(GR) / 20, lty = 2, col = "#c8322b"); ttl("C  p-values of all combinations", sprintf("%d combinations; %d with p < 0.05; dashed = expected if there were no signal", nrow(GR), nb))
  mtext("CROP-seq vs bulk over corrections, number of top genes, editing tiers and replicate choice", outer = TRUE, side = 3, adj = 0, line = 1.5, font = 2, cex = 1); mtext("KO vs the other KOs, no NTC; exact permutation of the bulk KO labels", outer = TRUE, side = 3, adj = 0, line = 0.3, cex = 0.65, col = INK2) })
