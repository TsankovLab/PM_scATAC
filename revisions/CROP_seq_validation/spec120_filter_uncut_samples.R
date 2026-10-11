# Spec 120: do the main analyses improve when the bulk samples without evidence of Cas9 cutting are removed?
# evidence per KO sample (spec118 counts at the 3 guide sites of its target, cut +/- 10 bp):
#   cut            = at least one guide with >= 2 read pairs carrying an indel and >= 2% of the pairs at that site (background 99th percentile 1.9%)
#   no evidence    = no such guide although >= 40 read pairs cover the sites
#   not assessable = fewer than 40 pairs at the sites (target barely expressed in the line)
# sample sets:  S70 (all) | CUT = S70 minus "no evidence" | CUT2 = CUT minus "not assessable" | CUT10 = only KO samples whose best guide has >= 10% indel pairs (NTC always kept)
#               NTCfix = S70 with the H2052 NTC that carries BPTF guides replaced by the 3 clean H2052-genotype NTCs | NTCfix_CUT = NTCfix minus "no evidence"
# analyses per set (as in spec103):
#   1  NCI-H2052, CROP-seq vs bulk: profiles centred on the mean of the groups present (NTC + KOs left in H2052; the screen is centred on the same groups), genes variable in both (top 1000 / 2000),
#      own - other Spearman, exact permutation p over the KOs present, per-TF own rho and rank
#   2  bulk across the 4 lines: limma t of NTC - KO, KO-centred, top-500 KO-signal genes (signal recomputed in the set), Spearman, KO-label permutation; no covariate and cell cycle regressed;
#      a KO with no sample left in a line is missing there (NA-aware agreement, as for RUN3)
#   3  Cm modules in their KO across lines (z within line, NTC - KO, lm z ~ line + KO): Cm20-TEAD4 and the pairs consistent in 4/4 lines in S70
# needs: spec118 (site counts), spec119 (guide reads), spec72 (S70), spec71 audit, spec75 (screen scores), spec_common.R
# outputs: spec120_sample_cut_evidence.csv, spec120_H2052.csv, spec120_H2052_per_KO.csv, spec120_across_lines.csv, spec120_across_lines_per_KO.csv, spec120_Cm_pairs.csv, fig36_filter_uncut_samples.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3", width = 220); N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000")); permsL <- make_perms(N_PERM)
CL <- "NCI-H2052"; fab <- load_full_audited_bulk(); meta <- fab$meta; cs <- fab$cs; q <- load_audit(); pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); S70 <- pr$samples; keep <- pr$genes
fp <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))); SH <- function(l) sub("NCI-|-211H", "", l)

## ---------------------------------------------------------------- evidence of cutting per sample
C <- read.csv(file.path(OUT, "spec118_site_counts.csv")); St <- read.csv(file.path(OUT, "spec118_gRNA_sites.csv")); St <- St[!is.na(St$chr), ]; C <- C[paste(C$guide, C$site) %in% paste(St$guide, sprintf("%s:%d", St$chr, St$cut)), ]
C$KO <- q$KO[match(C$sample, q$sample)]; C$line <- q$line[match(C$sample, q$sample)]; C$target <- sub("_g.*", "", C$guide); d <- C[!is.na(C$KO) & C$KO == C$target & C$KO %in% KOS, ]
E <- do.call(rbind, lapply(split(d, d$sample), function(x) data.frame(sample = x$sample[1], line = x$line[1], KO = x$KO[1], S70 = x$sample[1] %in% S70, pairs = sum(x$n_cover), indel_pairs = sum(x$n_indel), best_pct = max(ifelse(x$n_cover >= 5, 100 * x$n_indel / x$n_cover, 0)),
                                                              n_guides_edited = sum(x$n_indel >= 2 & x$n_indel / pmax(1, x$n_cover) >= 0.02), indel_per_guide = paste(sprintf("%d/%d", x$n_indel, x$n_cover), collapse = " "))))
E$call <- ifelse(E$n_guides_edited >= 1, "cut", ifelse(E$pairs >= 40, "no evidence", "not assessable")); E <- E[order(match(E$line, LINES), match(E$KO, KOS)), ]; write.csv(E, file.path(OUT, "spec120_sample_cut_evidence.csv"), row.names = FALSE)
Es <- E[E$S70, ]; cat("S70 KO samples:", nrow(Es), "\n"); print(table(Es$call)); cat("no evidence:", Es$sample[Es$call == "no evidence"], "\nnot assessable:", Es$sample[Es$call == "not assessable"], "\nbest guide < 10% (among cut):", Es$sample[Es$call == "cut" & Es$best_pct < 10], "\n")
SETS <- list(S70 = S70, CUT = setdiff(S70, Es$sample[Es$call == "no evidence"]), CUT2 = setdiff(S70, Es$sample[Es$call != "cut"]), CUT10 = setdiff(S70, Es$sample[!(Es$call == "cut" & Es$best_pct >= 10)]))
# spec119: the only H2052 NTC of S70 (H2052_NTC_R3) carries BPTF guides (and BPTF indels); the 3 H2052-genotype NTCs labelled H28 carry no guide -> NTCfix replaces it by those three
BAD_NTC <- "H2052_NTC_R3"; EXTRA <- c("H28_NTC_R1", "H28_NTC_R2", "H28_NTC_R3"); SETS$NTCfix <- c(setdiff(S70, BAD_NTC), EXTRA); SETS$NTCfix_CUT <- c(setdiff(SETS$CUT, BAD_NTC), EXTRA)
SL <- c(S70 = "all 70 samples", CUT = "minus no evidence", CUT2 = "minus no evidence / not assessable", CUT10 = "only best guide >= 10% indel", NTCfix = "H2052 NTC replaced", NTCfix_CUT = "H2052 NTC replaced, minus no evidence"); SC <- c(S70 = "grey55", CUT = "#f0a58f", CUT2 = "#eb6834", CUT10 = "#a8380b", NTCfix = "#8fb8ea", NTCfix_CUT = "#2a78d6")
for (sn in names(SETS)) { m <- meta[SETS[[sn]], ]; tb <- table(factor(m$Cell_Line, LINES), factor(m$Gene_Control, c("NTC", KOS))); cat(sprintf("\n%s (%s): %d samples\n", sn, SL[[sn]], length(SETS[[sn]]))); print(tb) }

## ---------------------------------------------------------------- 1  H2052 CROP-seq vs bulk
srt <- load_screen(); grp <- as.character(srt$merged_call); X <- GetAssayData(srt, layer = "counts"); GR <- c("NTC", KOS); PB <- sapply(GR, function(k) Matrix::rowSums(X[, grp == k, drop = FALSE]))
lc <- edgeR::cpm(edgeR::DGEList(PB), log = TRUE, prior.count = 2); expS <- rownames(lc)[rowSums(edgeR::cpm(edgeR::DGEList(PB)) >= 10) >= 1]; G <- intersect(keep, expS)
hvS <- intersect(VariableFeatures(FindVariableFeatures(srt, selection.method = "vst", nfeatures = nrow(srt), verbose = FALSE)), G)
perm_all <- function(n) { p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(seq_len(n))) }
rA <- list(); kA <- list()
for (sn in names(SETS)) { sH <- SETS[[sn]][meta[SETS[[sn]], "Cell_Line"] == CL]; gH <- meta[sH, "Gene_Control"]; yH <- normalise(cs[keep, sH], "tmm")$y; gp <- intersect(GR, unique(gH)); Kp <- setdiff(gp, "NTC")
  bH <- sapply(gp, function(k) rowMeans(yH[, gH == k, drop = FALSE])); vH <- G[order(matrixStats::rowVars(yH[G, ]), decreasing = TRUE)]; cg <- function(M) (M[, gp] - rowMeans(M[, gp]))[, Kp]; P8 <- perm_all(length(Kp))
  for (n in c(1000, 2000)) { f <- intersect(hvS[seq_len(n)], vH[seq_len(n)]); Cm <- cor(cg(lc[f, ]), cg(bH[f, ]), method = "spearman"); dm <- mean(diag(Cm)); nul <- apply(P8, 1, function(p) mean(Cm[cbind(seq_along(Kp), p)]))
    rA[[length(rA) + 1]] <- data.frame(set = sn, space = paste0("int", n), n_genes = length(f), n_KO = length(Kp), n_samples = length(sH), diag = dm, off = mean(Cm[row(Cm) != col(Cm)]), perm_p = mean(nul >= dm - 1e-12), top1 = sum(apply(Cm, 1, which.max) == seq_along(Kp)))
    kA[[length(kA) + 1]] <- data.frame(set = sn, space = paste0("int", n), KO = Kp, own = diag(Cm), rank = sapply(seq_along(Kp), function(i) rank(-Cm[i, ])[i]), n_bulk = as.integer(table(gH)[Kp])) } }
rA <- do.call(rbind, rA); rA$dmo <- rA$diag - rA$off; kA <- do.call(rbind, kA); write.csv(rA, file.path(OUT, "spec120_H2052.csv"), row.names = FALSE); write.csv(kA, file.path(OUT, "spec120_H2052_per_KO.csv"), row.names = FALSE)
cat("\n=== 1  H2052 CROP-seq vs bulk ===\n"); print(rA[order(rA$space), c("space", "set", "n_samples", "n_KO", "n_genes", "diag", "off", "dmo", "perm_p", "top1")], digits = 3, row.names = FALSE)
x <- kA[kA$space == "int2000", ]; cat("\n-- per TF own rho (rank), top 2000 in both\n"); print(reshape(transform(x[, c("KO", "set")], o = sprintf("%.2f (%s)", x$own, format(x$rank))), idvar = "KO", timevar = "set", direction = "wide"), row.names = FALSE)

## ---------------------------------------------------------------- 2  bulk across lines
rB <- list(); kB <- list(); AG <- list()
for (sn in names(SETS)) { ids <- SETS[[sn]]; TT <- list()
  for (cv in c("none", "cc")) { A <- array(NA_real_, c(length(keep), 8, 4), dimnames = list(keep, KOS, LINES)); for (cl in LINES) { s <- ids[meta[ids, "Cell_Line"] == cl]; A[, , cl] <- ko_effects(normalise(cs[keep, s], "tmm")$y, meta[s, "Gene_Control"], covar_matrix(q[s, ], cv))$t }; TT[[cv]] <- A }
  sig <- apply(abs(TT$none), c(1, 3), max, na.rm = TRUE)
  for (cv in c("none", "cc")) { a <- agreement(cl_from_array(TT[[cv]], "signal", sig), permsL); AG[[paste(sn, cv)]] <- a
    rB[[length(rB) + 1]] <- data.frame(set = sn, cov = cv, n_samples = length(ids), t(a$global), n_KO_p05 = sum(a$per_KO$perm_p < 0.05, na.rm = TRUE), n_KO_fdr05 = sum(a$per_KO$perm_fdr < 0.05, na.rm = TRUE)); kB[[length(kB) + 1]] <- cbind(set = sn, cov = cv, a$per_KO) } }
rB <- do.call(rbind, rB); kB <- do.call(rbind, kB); write.csv(rB, file.path(OUT, "spec120_across_lines.csv"), row.names = FALSE); write.csv(kB, file.path(OUT, "spec120_across_lines_per_KO.csv"), row.names = FALSE)
cat("\n=== 2  bulk across the 4 lines (KO-centred limma t, top-500 KO-signal genes) ===\n"); print(rB[order(rB$cov), c("cov", "set", "n_samples", "diag_mean", "off_mean", "perm_p", "perm_z", "n_KO_p05", "n_KO_fdr05")], digits = 3, row.names = FALSE)
for (cv in c("none", "cc")) { x <- kB[kB$cov == cv, ]; cat("\n-- per TF perm p (mean same-KO rho; line pairs),", cv, "\n"); print(reshape(transform(x[, c("KO", "set")], v = sprintf("%s (r %.2f; %d)%s", fp(x$perm_p), x$mean_intra, x$n_pairs, ifelse(!is.na(x$perm_fdr) & x$perm_fdr < 0.05, " *", ""))), idvar = "KO", timevar = "set", direction = "wide"), row.names = FALSE) }

## ---------------------------------------------------------------- 3  Cm pairs
SCO <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70N4%s.rds", DSUF)))$cm_scores$tmm;   # module scores of the 73-sample set (S70 + the 3 extra H2052 NTCs)
 pairstat <- function(ids, cm, ko) { dd <- data.frame(line = meta[ids, "Cell_Line"], group = meta[ids, "Gene_Control"], v = SCO[cm, ids]); dd <- dd[dd$group %in% GR, ]; dd$z <- ave(dd$v, dd$line, FUN = function(x) (x - mean(x)) / sd(x)); dd$cls <- ifelse(dd$group == "NTC", "NTC", ifelse(dd$group == ko, "KO", "other"))
  pl <- sapply(LINES, function(l) { e <- dd[dd$line == l, ]; mean(e$z[e$cls == "NTC"]) - mean(e$z[e$cls == "KO"]) }); f1 <- summary(lm(z ~ line + I(cls == "KO"), dd[dd$cls != "other", ]))$coef; f2 <- summary(lm(z ~ line + I(cls == "KO"), dd[dd$cls != "NTC", ]))$coef
  data.frame(Cm = cm, KO = ko, n_KO = sum(dd$cls == "KO"), NTC_minus_KO = -f1[nrow(f1), 1], p = f1[nrow(f1), 4], other_minus_KO = -f2[nrow(f2), 1], p_other = f2[nrow(f2), 4], lines_same_sign = max(sum(pl > 0, na.rm = TRUE), sum(pl < 0, na.rm = TRUE)), lines_with_KO = sum(!is.na(pl)), per_line = paste(ifelse(is.na(pl), "  n/a", sprintf("%+.2f", pl)), collapse = " ")) }
PAIRS_SHOW <- list(c("Cm20", "TEAD4"), c("Cm8", "TWIST1"), c("Cm9", "TWIST1"), c("Cm10", "TCF3"), c("Cm4", "PITX1"), c("Cm7", "PITX1"), c("Cm8", "PITX1"), c("Cm14", "PITX1"), c("Cm11", "TCF3"), c("Cm19", "TCF3"), c("Cm20", "TCF3"))
CP <- do.call(rbind, lapply(names(SETS), function(sn) cbind(set = sn, do.call(rbind, lapply(PAIRS_SHOW, function(p) pairstat(SETS[[sn]], p[1], p[2])))))); CP$pair <- paste(CP$Cm, CP$KO); write.csv(CP, file.path(OUT, "spec120_Cm_pairs.csv"), row.names = FALSE)
cat("\n=== 3  Cm x KO across lines (NTC - KO, SD; per line = H2052 MSTO H2452 H28) ===\n"); print(CP[order(match(CP$pair, unique(CP$pair))), c("pair", "set", "n_KO", "NTC_minus_KO", "p", "other_minus_KO", "p_other", "lines_same_sign", "lines_with_KO", "per_line")], digits = 2, row.names = FALSE)

## ---------------------------------------------------------------- figure
save_fig("fig36_filter_uncut_samples", 16, 10.2, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4.2, 0.4), family = "sans", mar = c(5, 4.6, 3.6, 1)); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.8); mtext(b, side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2) }; sn <- names(SETS)
  g <- function(col, sp) sapply(sn, function(s) rA[[col]][rA$set == s & rA$space == sp]); Dm <- rbind(g("dmo", "int1000"), g("dmo", "int2000")); P <- rbind(g("perm_p", "int1000"), g("perm_p", "int2000"))
  bp <- barplot(t(Dm), beside = TRUE, col = SC[names(SETS)], border = NA, las = 1, names.arg = c("top 1000 in both", "top 2000 in both"), cex.names = 0.8, ylab = "own-KO minus other-KO Spearman", cex.axis = 0.8, ylim = range(c(0, Dm)) + c(-0.01, 0.05)); abline(h = 0, col = INK)
  text(bp, pmax(t(Dm), 0) + 0.004, sprintf("p %s", fp(t(P))), cex = 0.6, font = ifelse(t(P) < 0.05, 2, 1)); legend("topleft", sprintf("%s (%d)", SL[names(SETS)], lengths(SETS)), fill = SC[names(SETS)], border = NA, bty = "n", cex = 0.6); ttl("A  H2052, CROP-seq vs bulk: all KOs present", "centred on the groups present; exact permutation p")
  x <- kA[kA$space == "int2000", ]; o <- KOS[order(-x$own[x$set == "S70"][match(KOS, x$KO[x$set == "S70"])])]; M <- t(sapply(sn, function(s) x$own[x$set == s][match(o, x$KO[x$set == s])])); R <- t(sapply(sn, function(s) x$rank[x$set == s][match(o, x$KO[x$set == s])]))
  bp <- barplot(M, beside = TRUE, col = SC[names(SETS)], border = NA, las = 1, names.arg = rep("", 8), ylab = "own-KO Spearman, CROP-seq vs bulk H2052", cex.axis = 0.8, ylim = range(M, na.rm = TRUE) + c(-0.06, 0.12)); abline(h = 0, col = INK); text(bp, ifelse(is.na(M), 0, pmax(M, 0)) + 0.015, ifelse(is.na(R), "x", sprintf("#%s", format(R))), cex = 0.5, font = ifelse(!is.na(R) & R == 1, 2, 1))
  text(colMeans(bp), par("usr")[3] - 0.02 * diff(par("usr")[3:4]), o, srt = 45, adj = 1, xpd = NA, cex = 0.8, font = 2, col = KO_COL[o]); ttl("B  H2052 per TF", "genes in the top 2000 of both; # = rank of the own KO; x = KO removed from H2052 in that set")
  Z <- sapply(c("none", "cc"), function(cv) rB$perm_z[rB$cov == cv][match(sn, rB$set[rB$cov == cv])]); Pz <- sapply(c("none", "cc"), function(cv) rB$perm_p[rB$cov == cv][match(sn, rB$set[rB$cov == cv])])
  bp <- barplot(Z, beside = TRUE, col = SC[names(SETS)], border = NA, las = 1, names.arg = c("no covariate", "cell cycle regressed"), cex.names = 0.8, ylab = "global agreement across lines (permutation z)", cex.axis = 0.8, ylim = range(c(0, Z)) + c(-0.1, 0.6)); abline(h = c(0, qnorm(0.95)), lty = c(1, 3), col = c(INK, MUTED))
  text(bp, pmax(Z, 0) + 0.08, fp(Pz), cex = 0.6, font = ifelse(Pz < 0.05, 2, 1)); ttl("C  bulk across the 4 lines: all KOs", "top-500 KO-signal genes, KO-centred limma t; numbers = permutation p")
  o2 <- read.csv(file.path(OUT, "spec82_TF_agreement_ranked.csv"))$KO
  for (cv in c("none", "cc")) { M <- t(sapply(sn, function(s) { k <- AG[[paste(s, cv)]]$per_KO; -log10(k$perm_p[match(o2, k$KO)]) })); Fd <- t(sapply(sn, function(s) { k <- AG[[paste(s, cv)]]$per_KO; k$perm_fdr[match(o2, k$KO)] }))
    bp <- barplot(M, beside = TRUE, col = SC[names(SETS)], border = NA, las = 1, names.arg = rep("", 8), ylab = "-log10 permutation p, same KO across lines", cex.axis = 0.8, ylim = c(0, max(M, na.rm = TRUE) + 0.5)); abline(h = -log10(0.05), lty = 3, col = MUTED); text(bp, M + 0.1, ifelse(!is.na(Fd) & Fd < 0.01, "**", ifelse(!is.na(Fd) & Fd < 0.05, "*", "")), cex = 0.9)
    text(colMeans(bp), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), o2, srt = 45, adj = 1, xpd = NA, cex = 0.8, font = 2, col = KO_COL[o2]); ttl(sprintf("%s  per TF across lines (%s)", if (cv == "none") "D" else "E", c(none = "no covariate", cc = "cell cycle regressed")[cv]), "dotted = p 0.05; * / ** FDR < 0.05 / 0.01") }
  prs <- unique(CP$pair); M <- t(sapply(sn, function(s) CP$NTC_minus_KO[CP$set == s][match(prs, CP$pair[CP$set == s])])); par(mar = c(6, 4.6, 3.6, 1))
  bp <- barplot(M, beside = TRUE, col = SC[names(SETS)], border = NA, las = 1, names.arg = rep("", length(prs)), ylab = "module, NTC - KO over the lines (SD; > 0 = lower in KO)", cex.axis = 0.8, ylim = range(M) + c(-0.3, 0.3)); abline(h = 0, col = INK); abline(v = mean(c(bp[length(sn), 4], bp[1, 5])), lty = 3, col = MUTED)
  text(colMeans(bp), par("usr")[3] - 0.03 * diff(par("usr")[3:4]), prs, srt = 45, adj = 1, xpd = NA, cex = 0.66, font = 2, col = KO_COL[sub(".* ", "", prs)]); ttl("F  Cm modules in their KO across bulk lines", "left of the dotted line: pairs from the text; right: pairs consistent in 4/4 lines")
  mtext("Removing bulk samples without evidence of Cas9 cutting at their guide sites: effect on the main analyses", outer = TRUE, side = 3, adj = 0, line = 2.6, font = 2, cex = 1)
  mtext("evidence = indels in the RNA-seq reads at the cut sites of the target's guides (spec118); grey = all S70 samples, darker = stricter filter; NTC samples always kept", outer = TRUE, side = 3, adj = 0, line = 1.4, cex = 0.65, col = INK2)
})
