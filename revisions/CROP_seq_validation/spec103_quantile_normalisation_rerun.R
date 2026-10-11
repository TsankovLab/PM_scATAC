# Spec 103: the main analyses rerun with quantile normalisation applied after the depth normalisation, side by side with the current normalisation.
#   tmm  = TMM logCPM (current)            tmmq = TMM logCPM, then quantile normalisation across all S70 bulk samples (limma::normalizeQuantiles; every sample gets the same distribution)
#   screen (H2052 comparison only): pseudobulk logCPM of the 9 guide groups, quantile-normalised across the 9 groups in the tmmq run
# 1  NCI-H2052, CROP-seq vs bulk: 9-group centring; genes variable in both (top 1000 / 2000), top 2000 H2052-variable, all shared; own - other Spearman, exact permutation p (8!), per TF own rho and rank
# 2  bulk across the 4 lines: limma t of NTC - KO, KO-centred, Spearman, KO-label permutation; top-500 KO-signal genes | top 2000 mean variance | all genes; no covariate and cell cycle regressed
# 3  Cm modules (AddModuleScore, top-50 genes, scored on each normalisation): module x KO effects across lines as in spec101 -- the four pairs from the text and the number / identity of consistent pairs
# needs: spec72 (S70: sample and gene lists), spec71 audit (cc), spec75 (screen module scores), spec_common.R
# outputs: spec103_H2052.csv, spec103_H2052_per_KO.csv, spec103_across_lines.csv, spec103_across_lines_per_KO.csv, spec103_Cm_pairs.csv, fig17_quantile_normalisation.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3", width = 220); N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000")); permsL <- make_perms(N_PERM)
CL <- "NCI-H2052"; GR <- c("NTC", KOS); fab <- load_full_audited_bulk(); meta <- fab$meta; cs <- fab$cs; q <- load_audit(); pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); ids <- pr$samples; keep <- pr$genes
fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2)); NC <- c(tmm = "grey55", tmmq = "#2a78d6"); NL <- c(tmm = "TMM logCPM (current)", tmmq = "TMM logCPM + quantile")
qn <- function(y) { z <- limma::normalizeQuantiles(y); dimnames(z) <- dimnames(y); z }
yT <- do.call(cbind, lapply(LINES, function(cl) normalise(cs[keep, ids[meta[ids, "Cell_Line"] == cl]], "tmm")$y))[, ids]; Y <- list(tmm = yT, tmmq = qn(yT))          # TMM per line as in spec72, quantiles across all 70 samples
cat("sample quantiles of TMM logCPM before quantile normalisation (range over the 70 samples): median", round(range(apply(yT, 2, median)), 2), " 90%", round(range(apply(yT, 2, quantile, 0.9)), 2), "\n")

## ---------------------------------------------------------------- 1  H2052 CROP-seq vs bulk
srt <- load_screen(); grp <- as.character(srt$merged_call); X <- GetAssayData(srt, layer = "counts"); PB <- sapply(GR, function(k) Matrix::rowSums(X[, grp == k, drop = FALSE]))
lc <- edgeR::cpm(edgeR::DGEList(PB), log = TRUE, prior.count = 2); expS <- rownames(lc)[rowSums(edgeR::cpm(edgeR::DGEList(PB)) >= 10) >= 1]; G <- intersect(keep, expS); LS <- list(tmm = lc[G, ], tmmq = qn(lc[expS, ])[G, ])
hvS <- intersect(VariableFeatures(FindVariableFeatures(srt, selection.method = "vst", nfeatures = nrow(srt), verbose = FALSE)), G); c9 <- function(M) (M[, GR] - rowMeans(M[, GR]))[, KOS]
perms8 <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(1:8)) })
sH <- ids[meta[ids, "Cell_Line"] == CL]; rA <- list(); kA <- list()
for (nm in names(Y)) { yH <- Y[[nm]][, sH]; bH <- sapply(GR, function(k) rowMeans(yH[, meta[sH, "Gene_Control"] == k, drop = FALSE])); vH <- G[order(matrixStats::rowVars(yH[G, ]), decreasing = TRUE)]
  SP <- list(int1000 = intersect(hvS[1:1000], vH[1:1000]), int2000 = intersect(hvS[1:2000], vH[1:2000]), top2000 = vH[1:2000], all = G)
  for (sp in names(SP)) { f <- SP[[sp]]; C <- cor(c9(LS[[nm]][f, ]), c9(bH[f, ]), method = "spearman"); dm <- mean(diag(C)); nul <- apply(perms8, 1, function(p) mean(C[cbind(1:8, p)]))
    rA[[length(rA) + 1]] <- data.frame(norm = nm, space = sp, n_genes = length(f), diag = dm, off = mean(C[row(C) != col(C)]), perm_p = mean(nul >= dm - 1e-12), top1 = sum(apply(C, 1, which.max) == 1:8))
    kA[[length(kA) + 1]] <- data.frame(norm = nm, space = sp, KO = KOS, own = diag(C), rank = sapply(1:8, function(i) rank(-C[i, ])[i])) } }
rA <- do.call(rbind, rA); rA$dmo <- rA$diag - rA$off; kA <- do.call(rbind, kA); write.csv(rA, file.path(OUT, "spec103_H2052.csv"), row.names = FALSE); write.csv(kA, file.path(OUT, "spec103_H2052_per_KO.csv"), row.names = FALSE)
cat("\n=== 1  H2052 CROP-seq vs bulk (9-group centring) ===\n"); print(rA[order(rA$space), c("space", "norm", "n_genes", "diag", "off", "dmo", "perm_p", "top1")], digits = 3, row.names = FALSE)
x <- kA[kA$space == "int2000", ]; cat("\n-- per TF own rho (rank), genes in the top 2000 of both\n"); print(reshape(transform(x[, c("KO", "norm")], o = sprintf("%.2f (%s)", x$own, format(x$rank))), idvar = "KO", timevar = "norm", direction = "wide"), row.names = FALSE)

## ---------------------------------------------------------------- 2  bulk across lines
rB <- list(); kB <- list(); AG <- list()
for (nm in names(Y)) { TT <- list(); for (cv in c("none", "cc")) { A <- array(NA_real_, c(length(keep), 8, 4), dimnames = list(keep, KOS, LINES)); for (cl in LINES) { s <- ids[meta[ids, "Cell_Line"] == cl]; A[, , cl] <- ko_effects(Y[[nm]][, s], meta[s, "Gene_Control"], covar_matrix(q[s, ], cv))$t }; TT[[cv]] <- A }
  sig <- apply(abs(TT$none), c(1, 3), max); gv <- sapply(LINES, function(cl) matrixStats::rowVars(Y[[nm]][, ids[meta[ids, "Cell_Line"] == cl]])); vAll <- keep[order(rowMeans(gv), decreasing = TRUE)]
  for (cv in c("none", "cc")) for (sp in c("signal500", "top2000", "all")) { CLs <- switch(sp, signal500 = cl_from_array(TT[[cv]], "signal", sig), top2000 = cl_from_array(TT[[cv]], vAll[1:2000]), all = cl_from_array(TT[[cv]], NULL)); a <- agreement(CLs, permsL); AG[[paste(nm, cv, sp)]] <- a
    rB[[length(rB) + 1]] <- data.frame(norm = nm, cov = cv, space = sp, t(a$global), n_KO_p05 = sum(a$per_KO$perm_p < 0.05), n_KO_fdr05 = sum(a$per_KO$perm_fdr < 0.05)); kB[[length(kB) + 1]] <- cbind(norm = nm, cov = cv, space = sp, a$per_KO) } }
rB <- do.call(rbind, rB); kB <- do.call(rbind, kB); write.csv(rB, file.path(OUT, "spec103_across_lines.csv"), row.names = FALSE); write.csv(kB, file.path(OUT, "spec103_across_lines_per_KO.csv"), row.names = FALSE)
cat("\n=== 2  bulk across the 4 lines (KO-centred limma t) ===\n"); print(rB[order(rB$cov, rB$space), c("cov", "space", "norm", "diag_mean", "off_mean", "perm_p", "perm_z", "n_KO_p05", "n_KO_fdr05")], digits = 3, row.names = FALSE)
for (cv in c("none", "cc")) { x <- kB[kB$space == "signal500" & kB$cov == cv, ]; cat("\n-- per TF perm p (mean same-KO rho), top-500 KO-signal genes,", cv, "\n"); print(reshape(transform(x[, c("KO", "norm")], v = sprintf("%.4f (r %.2f)%s", x$perm_p, x$mean_intra, ifelse(x$perm_fdr < 0.05, " *", ""))), idvar = "KO", timevar = "norm", direction = "wide"), row.names = FALSE) }

## ---------------------------------------------------------------- 3  Cm modules
top50 <- cm_top50(); CM <- names(top50); det <- rownames(cs)[rowSums(cs[, ids]) > 0]; yA <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(cs[det, ids])), log = TRUE, prior.count = 2)
score <- function(y) { y <- y - min(y); sb <- CreateSeuratObject(counts = cs[rownames(y), ids]); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(y)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1)); x <- t(as.matrix(sb@meta.data[ids, paste0("Cm_", seq_along(CM))])); rownames(x) <- CM; x }
SCO <- list(tmm = pr$cm_scores$tmm[CM, ids], tmmq = score(qn(yA))); sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds"))
pairstat <- function(S, cm, ko) { d <- data.frame(line = meta[ids, "Cell_Line"], group = meta[ids, "Gene_Control"], v = S[cm, ids]); d <- d[d$group %in% GR, ]; d$z <- ave(d$v, d$line, FUN = function(x) (x - mean(x)) / sd(x)); d$cls <- ifelse(d$group == "NTC", "NTC", ifelse(d$group == ko, "KO", "other"))
  pl <- t(sapply(LINES, function(l) { e <- d[d$line == l, ]; c(mean(e$z[e$cls == "NTC"]) - mean(e$z[e$cls == "KO"]), mean(e$z[e$cls == "other"]) - mean(e$z[e$cls == "KO"])) }))
  f1 <- summary(lm(z ~ line + I(cls == "KO"), d[d$cls != "other", ]))$coef; f2 <- summary(lm(z ~ line + I(cls == "KO"), d[d$cls != "NTC", ]))$coef; e1 <- -f1[nrow(f1), 1]; e2 <- -f2[nrow(f2), 1]
  data.frame(Cm = cm, KO = ko, NTC_minus_KO = e1, p = f1[nrow(f1), 4], n_same = max(sum(pl[, 1] > 0), sum(pl[, 1] < 0)), other_minus_KO = e2, p_other = f2[nrow(f2), 4], n_same_other = max(sum(pl[, 2] > 0), sum(pl[, 2] < 0)),
             per_line = paste(sprintf("%+.2f", pl[, 1]), collapse = " "), consistent = all(pl[, 1] > 0) == TRUE & all(pl[, 2] > 0) | all(pl[, 1] < 0) & all(pl[, 2] < 0)) }
CP <- do.call(rbind, lapply(names(SCO), function(nm) cbind(norm = nm, do.call(rbind, lapply(CM, function(cm) do.call(rbind, lapply(KOS, function(ko) pairstat(SCO[[nm]], cm, ko)))))))); CP$pair <- paste(CP$Cm, CP$KO); write.csv(CP, file.path(OUT, "spec103_Cm_pairs.csv"), row.names = FALSE)
TXT <- c("Cm20 TEAD4", "Cm8 TWIST1", "Cm9 TWIST1", "Cm10 TCF3"); cons <- lapply(split(CP, CP$norm), function(d) d$pair[d$consistent]); SHOW <- unique(c(TXT, unlist(cons)))
cat("\n=== 3  Cm x KO across lines (z within line; NTC - KO combined over the 4 lines; per line = H2052 MSTO H2452 H28) ===\n"); x <- CP[CP$pair %in% SHOW, ]; x <- x[order(match(x$pair, SHOW), x$norm), ]
print(x[, c("pair", "norm", "NTC_minus_KO", "p", "n_same", "other_minus_KO", "p_other", "n_same_other", "consistent", "per_line")], digits = 2, row.names = FALSE)
cat("\nconsistent pairs (4/4 lines, both contrasts): tmm", length(cons$tmm), "| tmmq", length(cons$tmmq), "| shared", length(intersect(cons$tmm, cons$tmmq)), "\n  tmm only:", setdiff(cons$tmm, cons$tmmq), "\n  tmmq only:", setdiff(cons$tmmq, cons$tmm), "\n")
cat("correlation of module scores between the two normalisations (per module, over the 70 samples): median", round(median(sapply(CM, function(m) cor(SCO$tmm[m, ], SCO$tmmq[m, ]))), 3), " min", round(min(sapply(CM, function(m) cor(SCO$tmm[m, ], SCO$tmmq[m, ]))), 3), "\n")

## ---------------------------------------------------------------- figure
save_fig("fig17_quantile_normalisation", 16, 10.4, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4.2, 0.4), family = "sans", mar = c(5.2, 4.6, 3.6, 1)); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.8); mtext(b, side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2) }
  sps <- c("int1000", "int2000", "top2000", "all"); g <- function(col) sapply(names(NC), function(n) rA[[col]][rA$norm == n][match(sps, rA$space[rA$norm == n])]); Dm <- g("dmo"); P <- g("perm_p")
  bp <- barplot(t(Dm), beside = TRUE, col = NC, border = NA, las = 1, names.arg = c("top 1000\nin both", "top 2000\nin both", "top 2000\nH2052", "all shared"), cex.names = 0.7, ylab = "own-KO minus other-KO Spearman", cex.axis = 0.8, ylim = range(c(0, Dm)) + c(-0.01, 0.03)); abline(h = 0, col = INK)
  text(bp, pmax(t(Dm), 0) + 0.003, sprintf("p %s", fp(t(P))), cex = 0.6, font = ifelse(t(P) < 0.05, 2, 1)); legend("topright", NL, fill = NC, border = NA, bty = "n", cex = 0.7); ttl("A  H2052, CROP-seq vs bulk: all 8 KOs", "9-group centring; exact permutation p")
  o <- KOS[order(-kA$own[kA$norm == "tmm" & kA$space == "int2000"])]; M <- t(sapply(names(NC), function(n) { x <- kA[kA$norm == n & kA$space == "int2000", ]; x$own[match(o, x$KO)] })); R <- t(sapply(names(NC), function(n) { x <- kA[kA$norm == n & kA$space == "int2000", ]; x$rank[match(o, x$KO)] }))
  bp <- barplot(M, beside = TRUE, col = NC, border = NA, las = 1, names.arg = rep("", 8), ylab = "own-KO Spearman, CROP-seq vs bulk H2052", cex.axis = 0.8, ylim = range(M) + c(-0.05, 0.12)); abline(h = 0, col = INK); text(bp, pmax(M, 0) + 0.015, sprintf("#%s", format(R)), cex = 0.58, font = ifelse(R == 1, 2, 1))
  text(colMeans(bp), par("usr")[3] - 0.02 * diff(par("usr")[3:4]), o, srt = 45, adj = 1, xpd = NA, cex = 0.8, font = 2, col = KO_COL[o]); ttl("B  H2052 per TF", "genes in the top 2000 variable of both; # = rank of the own KO among 8")
  lab <- c(); Z <- NULL; Pz <- NULL; for (cv in c("none", "cc")) for (sp in c("signal500", "top2000", "all")) { x <- rB[rB$cov == cv & rB$space == sp, ]; Z <- cbind(Z, x$perm_z[match(names(NC), x$norm)]); Pz <- cbind(Pz, x$perm_p[match(names(NC), x$norm)]); lab <- c(lab, sprintf("%s\n%s", c(signal500 = "KO-signal 500", top2000 = "top 2000 var", all = "all genes")[sp], c(none = "no cov", cc = "cc regressed")[cv])) }
  bp <- barplot(Z, beside = TRUE, col = NC, border = NA, las = 1, names.arg = rep("", 6), ylab = "global agreement across lines (permutation z)", cex.axis = 0.8, ylim = range(c(0, Z)) + c(-0.1, 0.6)); abline(h = c(0, qnorm(0.95)), lty = c(1, 3), col = c(INK, MUTED))
  text(bp, pmax(Z, 0) + 0.08, fp(Pz), cex = 0.56, font = ifelse(Pz < 0.05, 2, 1)); text(colMeans(bp), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), lab, srt = 30, adj = 1, xpd = NA, cex = 0.66); ttl("C  bulk across the 4 lines: all KOs", "KO-centred limma t; numbers = permutation p")
  o2 <- read.csv(file.path(OUT, "spec82_TF_agreement_ranked.csv"))$KO
  for (cv in c("none", "cc")) { M <- t(sapply(names(NC), function(n) { k <- AG[[paste(n, cv, "signal500")]]$per_KO; -log10(k$perm_p[match(o2, k$KO)]) })); Fd <- t(sapply(names(NC), function(n) { k <- AG[[paste(n, cv, "signal500")]]$per_KO; k$perm_fdr[match(o2, k$KO)] }))
    bp <- barplot(M, beside = TRUE, col = NC, border = NA, las = 1, names.arg = rep("", 8), ylab = "-log10 permutation p, same KO across lines", cex.axis = 0.8, ylim = c(0, max(M) + 0.5)); abline(h = -log10(0.05), lty = 3, col = MUTED); text(bp, M + 0.1, ifelse(Fd < 0.01, "**", ifelse(Fd < 0.05, "*", "")), cex = 0.9)
    text(colMeans(bp), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), o2, srt = 45, adj = 1, xpd = NA, cex = 0.8, font = 2, col = KO_COL[o2]); gp <- sapply(names(NC), function(n) AG[[paste(n, cv, "signal500")]]$global[["perm_p"]])
    ttl(sprintf("%s  per TF across lines (%s)", if (cv == "none") "D" else "E", c(none = "no covariate", cc = "cell cycle regressed")[cv]), sprintf("top-500 KO-signal genes; global p: current %s, quantile %s; dotted = 0.05; * / ** FDR < 0.05 / 0.01", fp(gp[1]), fp(gp[2]))) }
  M <- t(sapply(names(NC), function(n) CP$NTC_minus_KO[CP$norm == n][match(SHOW, CP$pair[CP$norm == n])])); Pp <- t(sapply(names(NC), function(n) CP$p[CP$norm == n][match(SHOW, CP$pair[CP$norm == n])])); Cn <- t(sapply(names(NC), function(n) CP$consistent[CP$norm == n][match(SHOW, CP$pair[CP$norm == n])]))
  par(mar = c(6, 4.6, 3.6, 1)); bp <- barplot(M, beside = TRUE, col = NC, border = NA, las = 1, names.arg = rep("", length(SHOW)), ylab = "module, NTC - KO over the 4 lines (SD; > 0 = lower in KO)", cex.axis = 0.8, ylim = range(M) + c(-0.5, 0.3)); abline(h = 0, col = INK); abline(v = mean(c(bp[2, length(TXT)], bp[1, length(TXT) + 1])), lty = 3, col = MUTED)
  text(bp, ifelse(M > 0, M + 0.08, M - 0.08), ifelse(Cn, "4/4", ""), cex = 0.5); text(colMeans(bp), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), SHOW, srt = 45, adj = 1, xpd = NA, cex = 0.66, font = 2, col = KO_COL[sub(".* ", "", SHOW)])
  ttl("F  Cm modules in their KO across bulk lines", "left of the dotted line: pairs from the text; right: pairs consistent in 4/4 lines for both contrasts (marked 4/4)")
  mtext("Quantile normalisation after depth normalisation vs the current normalisation: main analyses", outer = TRUE, side = 3, adj = 0, line = 2.6, font = 2, cex = 1)
  mtext("grey = TMM logCPM (current); blue = TMM logCPM then quantile normalisation across the 70 bulk samples (screen pseudobulks quantile-normalised across the 9 groups in A-B); S70", outer = TRUE, side = 3, adj = 0, line = 1.4, cex = 0.65, col = INK2)
})
