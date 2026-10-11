# Spec 122: main analyses rerun on the corrected sample set S91 (agree_common.R::load_bulk_S91):
#   - the H2052 controls are H28_NTC_R1-R3 (H2052 genotype, no guide reads); H2052_NTC_R3, which carries BPTF guides, is relabelled BPTF KO
#   - the RUNX1, RUNX2 and SNAI2 KO samples are included (11 KOs; not in the screen, so they only enter the bulk analyses)
# sets compared: S70      = the previous set and labels (8 KOs)
#                S73fix   = new H2052 controls + relabelled BPTF, 8 KOs
#                S91      = S73fix + RUNX1 / RUNX2 / SNAI2 samples as labelled (11 KOs)
#                S91b     = S91 with H2052_RUNX1_R1/R3 relabelled PITX1 (they carry PITX1 guides and PITX1 indels, spec119 / spec121)
# analyses:
#   1  NCI-H2052, CROP-seq vs bulk (the 8 shared KOs): profiles centred on the mean of NTC + the 8 KOs in both data sets, genes variable in both (top 1000 / 2000), own - other Spearman, exact permutation (8!), per TF
#   2  bulk across the 4 lines: limma t of NTC - KO, centred over the KOs of the line, top-500 KO-signal genes (signal recomputed per set from the 2 lines not compared), Spearman, KO-label permutation
#      (NA-aware: a KO absent from a line is missing there); no covariate and cell cycle regressed (cell-cycle score recomputed over the samples of the set)
#   3  Cm modules (AddModuleScore, top-50 genes, scored over the samples of the set): NTC - KO per line and combined, for the pairs of the text and the pairs consistent in 4/4 lines
#   4  NCI-H2052 KO x Cm heatmap (median NTC - median KO, 8 KOs x 20 Cm) vs the CROP-seq heatmap (spec81): full / shared (module means) / KO-specific Spearman
# needs: spec72 (S70 sample and gene lists), spec71 audit, spec81 (screen Cm x KO), spec_common.R
# outputs: spec122_samples.csv, spec122_H2052.csv, spec122_H2052_per_KO.csv, spec122_across_lines.csv, spec122_across_lines_per_KO.csv, spec122_Cm_pairs.csv, spec122_Cm_heatmap_H2052.csv
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3", width = 230); N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "20000"))
CL <- "NCI-H2052"; q <- load_audit(); pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); S70 <- pr$samples; keep <- pr$genes; fab0 <- load_full_audited_bulk(); top50 <- cm_top50(); CM <- names(top50)
fp <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))); SH <- function(l) sub("NCI-|-211H", "", l); kcol <- function(k) ifelse(k %in% names(KO_COL), KO_COL[k], "grey30")
score <- function(cs, ids) { y <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(cs[, ids])), log = TRUE, prior.count = 2); y <- y[rowSums(cs[, ids]) > 0, ]; y <- y - min(y); sb <- CreateSeuratObject(counts = cs[rownames(y), ids]); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(y))
  sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1)); x <- t(as.matrix(sb@meta.data[ids, paste0("Cm_", seq_along(CM))])); rownames(x) <- CM; x }
mk <- function(B, kos) { B$kos <- kos; B$cm <- score(B$cs, B$ids); B }
B70 <- list(cs = fab0$cs[, S70], meta = data.frame(fab0$meta[S70, c("Cell_Line", "Gene_Control")], label_original = fab0$meta[S70, "Gene_Control"]), ids = S70, cc = setNames(q[S70, "cc"], S70))
SETS <- list(S70 = mk(B70, KOS), S73fix = mk(load_bulk_S91(S70, extra_kos = FALSE), KOS), S91 = mk(load_bulk_S91(S70), KOS11), S91b = mk(load_bulk_S91(S70, runx1_h2052_as_pitx1 = TRUE), KOS11))
SL <- c(S70 = "previous set (70)", S73fix = "new H2052 controls, BPTF relabelled (73)", S91 = "+ RUNX1, RUNX2, SNAI2 (91)", S91b = "S91, H2052 RUNX1 samples as PITX1"); SC <- c(S70 = "grey55", S73fix = "#8fb8ea", S91 = "#2a78d6", S91b = "#1baf7a")
m <- SETS$S91$meta; write.csv(data.frame(sample = rownames(m), m, cc = SETS$S91$cc), file.path(OUT, "spec122_samples.csv"), row.names = FALSE); cat("S91:", nrow(m), "samples\n"); print(table(factor(m$Cell_Line, LINES), factor(m$Gene_Control, c("NTC", KOS11))))
cat("relabelled:", rownames(m)[m$Gene_Control != m$label_original], "| H2052 NTC:", rownames(m)[m$Cell_Line == CL & m$Gene_Control == "NTC"], "\n")

## ---------------------------------------------------------------- 1  H2052 CROP-seq vs bulk (8 shared KOs)
srt <- load_screen(); grp <- as.character(srt$merged_call); X <- GetAssayData(srt, layer = "counts"); GR <- c("NTC", KOS); PB <- sapply(GR, function(k) Matrix::rowSums(X[, grp == k, drop = FALSE]))
lc <- edgeR::cpm(edgeR::DGEList(PB), log = TRUE, prior.count = 2); expS <- rownames(lc)[rowSums(edgeR::cpm(edgeR::DGEList(PB)) >= 10) >= 1]; G <- intersect(keep, expS)
hvS <- intersect(VariableFeatures(FindVariableFeatures(srt, selection.method = "vst", nfeatures = nrow(srt), verbose = FALSE)), G); c9 <- function(M) (M[, GR] - rowMeans(M[, GR]))[, KOS]
perms8 <- local({ p <- function(v) if (length(v) <= 1) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(p(v[-i]), function(r) c(v[i], r)))); do.call(rbind, p(1:8)) }); rA <- list(); kA <- list()
for (sn in names(SETS)) { B <- SETS[[sn]]; sH <- B$ids[B$meta[B$ids, "Cell_Line"] == CL & B$meta[B$ids, "Gene_Control"] %in% GR]; gH <- B$meta[sH, "Gene_Control"]; yH <- normalise(B$cs[keep, sH], "tmm")$y
  bH <- sapply(GR, function(k) rowMeans(yH[, gH == k, drop = FALSE])); vH <- G[order(matrixStats::rowVars(yH[G, ]), decreasing = TRUE)]
  for (n in c(1000, 2000)) { f <- intersect(hvS[seq_len(n)], vH[seq_len(n)]); Cm <- cor(c9(lc[f, ]), c9(bH[f, ]), method = "spearman"); dm <- mean(diag(Cm)); nul <- apply(perms8, 1, function(p) mean(Cm[cbind(1:8, p)]))
    rA[[length(rA) + 1]] <- data.frame(set = sn, space = paste0("int", n), n_genes = length(f), n_samples = length(sH), n_NTC = sum(gH == "NTC"), diag = dm, off = mean(Cm[row(Cm) != col(Cm)]), perm_p = mean(nul >= dm - 1e-12), top1 = sum(apply(Cm, 1, which.max) == 1:8))
    kA[[length(kA) + 1]] <- data.frame(set = sn, space = paste0("int", n), KO = KOS, own = diag(Cm), rank = sapply(1:8, function(i) rank(-Cm[i, ])[i])) } }
rA <- do.call(rbind, rA); rA$dmo <- rA$diag - rA$off; kA <- do.call(rbind, kA); write.csv(rA, file.path(OUT, "spec122_H2052.csv"), row.names = FALSE); write.csv(kA, file.path(OUT, "spec122_H2052_per_KO.csv"), row.names = FALSE)
cat("\n=== 1  H2052 CROP-seq vs bulk (8 shared KOs) ===\n"); print(rA[order(rA$space), c("space", "set", "n_samples", "n_NTC", "n_genes", "diag", "off", "dmo", "perm_p", "top1")], digits = 3, row.names = FALSE)
x <- kA[kA$space == "int2000", ]; cat("\n-- per TF own rho (rank), top 2000 in both\n"); print(reshape(transform(x[, c("KO", "set")], o = sprintf("%.2f (%s)", x$own, format(x$rank))), idvar = "KO", timevar = "set", direction = "wide"), row.names = FALSE)

## ---------------------------------------------------------------- 2  bulk across lines
rB <- list(); kB <- list(); AG <- list()
for (sn in names(SETS)) { B <- SETS[[sn]]; ks <- B$kos; perms <- make_perms(N_PERM, K = length(ks)); TT <- list()
  for (cv in c("none", "cc")) { A <- array(NA_real_, c(length(keep), length(ks), 4), dimnames = list(keep, ks, LINES))
    for (cl in LINES) { s <- B$ids[B$meta[B$ids, "Cell_Line"] == cl]; Xc <- if (cv == "cc") scale(cbind(cc = B$cc[s]), scale = FALSE) else NULL; A[, , cl] <- ko_effects(normalise(B$cs[keep, s], "tmm")$y, B$meta[s, "Gene_Control"], Xc, kos = ks)$t }; TT[[cv]] <- A }
  sig <- apply(abs(TT$none), c(1, 3), max, na.rm = TRUE)
  for (cv in c("none", "cc")) { a <- agreement(cl_from_array(TT[[cv]], "signal", sig), perms); AG[[paste(sn, cv)]] <- a
    rB[[length(rB) + 1]] <- data.frame(set = sn, cov = cv, n_samples = length(B$ids), n_KO = length(ks), t(a$global), n_KO_p05 = sum(a$per_KO$perm_p < 0.05, na.rm = TRUE), n_KO_fdr05 = sum(a$per_KO$perm_fdr < 0.05, na.rm = TRUE)); kB[[length(kB) + 1]] <- cbind(set = sn, cov = cv, a$per_KO) } }
rB <- do.call(rbind, rB); kB <- do.call(rbind, kB); write.csv(rB, file.path(OUT, "spec122_across_lines.csv"), row.names = FALSE); write.csv(kB, file.path(OUT, "spec122_across_lines_per_KO.csv"), row.names = FALSE)
cat("\n=== 2  bulk across the 4 lines (KO-centred limma t, top-500 KO-signal genes) ===\n"); print(rB[order(rB$cov), c("cov", "set", "n_samples", "n_KO", "diag_mean", "off_mean", "perm_p", "perm_z", "n_KO_p05", "n_KO_fdr05")], digits = 3, row.names = FALSE)
for (cv in c("none", "cc")) { x <- kB[kB$cov == cv, ]; cat("\n-- per TF perm p (mean same-KO rho; line pairs),", cv, "\n"); print(reshape(transform(x[, c("KO", "set")], v = sprintf("%s (r %.2f; %d)%s", fp(x$perm_p), x$mean_intra, x$n_pairs, ifelse(!is.na(x$perm_fdr) & x$perm_fdr < 0.05, " *", ""))), idvar = "KO", timevar = "set", direction = "wide"), row.names = FALSE) }

## ---------------------------------------------------------------- 3  Cm pairs
pairstat <- function(B, cm, ko) { ids <- B$ids; dd <- data.frame(line = B$meta[ids, "Cell_Line"], group = B$meta[ids, "Gene_Control"], v = B$cm[cm, ids]); dd$z <- ave(dd$v, dd$line, FUN = function(x) (x - mean(x)) / sd(x)); dd$cls <- ifelse(dd$group == "NTC", "NTC", ifelse(dd$group == ko, "KO", "other"))
  pl <- sapply(LINES, function(l) { e <- dd[dd$line == l, ]; mean(e$z[e$cls == "NTC"]) - mean(e$z[e$cls == "KO"]) }); f1 <- summary(lm(z ~ line + I(cls == "KO"), dd[dd$cls != "other", ]))$coef; f2 <- summary(lm(z ~ line + I(cls == "KO"), dd[dd$cls != "NTC", ]))$coef
  data.frame(Cm = cm, KO = ko, n_KO = sum(dd$cls == "KO"), NTC_minus_KO = -f1[nrow(f1), 1], p = f1[nrow(f1), 4], other_minus_KO = -f2[nrow(f2), 1], p_other = f2[nrow(f2), 4], lines_same_sign = max(sum(pl > 0, na.rm = TRUE), sum(pl < 0, na.rm = TRUE)), lines_with_KO = sum(!is.na(pl)), per_line = paste(ifelse(is.na(pl), "  n/a", sprintf("%+.2f", pl)), collapse = " ")) }
PAIRS_SHOW <- list(c("Cm20", "TEAD4"), c("Cm8", "TWIST1"), c("Cm9", "TWIST1"), c("Cm10", "TCF3"), c("Cm1", "SOX9"), c("Cm1", "TEAD2"), c("Cm4", "PITX1"), c("Cm7", "PITX1"), c("Cm8", "PITX1"), c("Cm14", "PITX1"), c("Cm11", "TCF3"), c("Cm19", "TCF3"), c("Cm20", "TCF3"))
CP <- do.call(rbind, lapply(names(SETS), function(sn) cbind(set = sn, do.call(rbind, lapply(PAIRS_SHOW, function(p) pairstat(SETS[[sn]], p[1], p[2])))))); CP$pair <- paste(CP$Cm, CP$KO); write.csv(CP, file.path(OUT, "spec122_Cm_pairs.csv"), row.names = FALSE)
cat("\n=== 3  Cm x KO across lines (NTC - KO, SD within line; per line = H2052 MSTO H2452 H28) ===\n"); print(CP[order(match(CP$pair, unique(CP$pair))), c("pair", "set", "n_KO", "NTC_minus_KO", "p", "other_minus_KO", "p_other", "lines_same_sign", "lines_with_KO", "per_line")], digits = 2, row.names = FALSE)
# all 20 x 8 pairs: how many are consistent in 4/4 lines for both contrasts, per set
ALLP <- do.call(rbind, lapply(names(SETS), function(sn) cbind(set = sn, do.call(rbind, lapply(CM, function(cm) do.call(rbind, lapply(KOS, function(ko) pairstat(SETS[[sn]], cm, ko))))))))
cat("\npairs (20 Cm x 8 KOs) with NTC - KO of the same sign in 4/4 lines and p < 0.05:\n"); for (sn in names(SETS)) { d <- ALLP[ALLP$set == sn & ALLP$lines_same_sign == 4 & ALLP$p < 0.05, ]; d <- d[order(d$p), ]; cat(sprintf("  %-7s %2d: %s\n", sn, nrow(d), paste(sprintf("%s-%s(%+.1f)", d$Cm, d$KO, d$NTC_minus_KO), collapse = " "))) }

## ---------------------------------------------------------------- 4  H2052 KO x Cm heatmap vs CROP-seq
x81 <- read.csv(file.path(OUT, "spec81_Cm_KO_NTCminusKO_screen_and_bulk.csv")); d <- x81[x81$source == "screen", ]; S <- tapply(d$NTC_minus_KO, list(d$KO, d$Cm), mean)[KOS, CM]; HM <- list()
for (sn in names(SETS)) { B <- SETS[[sn]]; s <- B$ids[B$meta[B$ids, "Cell_Line"] == CL]; g <- B$meta[s, "Gene_Control"]; Bm <- t(sapply(KOS, function(k) apply(B$cm[, s[g == "NTC"], drop = FALSE], 1, median) - apply(B$cm[, s[g == k], drop = FALSE], 1, median)))[, CM]
  ct <- function(a, b) { r <- cor.test(as.vector(a), as.vector(b), method = "spearman", exact = FALSE); c(r$estimate, r$p.value) }; f <- ct(S, Bm); sh <- ct(colMeans(S), colMeans(Bm)); spc <- ct(sweep(S, 2, colMeans(S)), sweep(Bm, 2, colMeans(Bm)))
  HM[[sn]] <- data.frame(set = sn, n_NTC = sum(g == "NTC"), full_rho = f[1], full_p = f[2], shared_rho = sh[1], shared_p = sh[2], specific_rho = spc[1], specific_p = spc[2], frac_cells_same_sign = mean(sign(S) == sign(Bm))) }
HM <- do.call(rbind, HM); write.csv(HM, file.path(OUT, "spec122_Cm_heatmap_H2052.csv"), row.names = FALSE); cat("\n=== 4  H2052 KO x Cm heatmap (median NTC - median KO) vs CROP-seq, 160 cells ===\n"); print(HM, digits = 2, row.names = FALSE)

