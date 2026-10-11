# Spec 71: audit of the 76 bulk RNA-seq samples (load_full_audited_bulk(): NTC + the 8 CROP-seq KOs, 4 lines) -- read depth, QC, PCA, objective outlier flags,
# KO-label check and the candidate sample sets used by spec72-76. Replaces the eye-picked VISUAL_OUTLIERS with rules:
#   low_depth    : reads assigned to genes < 5 M
#   qc_outlier   : robust z (median/MAD within the line) beyond 3 on the bad side for log10 assigned reads, % assigned, genes detected, % mito, top-100-gene share
#   expr_outlier : robust z < -3 for the sample's median Spearman to the other samples of its line, or > 3 for its distance from the line centre in PC1-5
#                  (TMM logCPM, filterByExpr genes, cell cycle regressed)
#   KO_label     : own target gene NOT among the 2 lowest of the 8 targets (per-line z of logCPM) -> KO label unsupported by expression (weak KO or swap)
# PCA per line in 3 versions (raw | cc regressed | cc + log depth regressed) and R2 of PC1-5 with KO, replicate number, depth, % assigned, cc.
# Sample sets: S70 = load_bulk() screen KOs (baseline) | S76 = all | SQC = 76 minus low_depth / qc_outlier / expr_outlier / mixed samples. Compared on a criterion that
# does not use cross-line agreement: within-line replicate structure (AUC of same-KO vs different-KO sample pairs, centred Spearman, as spec11).
# needs: counts/*_ReadsPerGene.out.tab, logs/*_Log.final.out
# outputs: spec71_sample_audit.csv, spec71_PC_covariate_R2.csv, spec71_KO_label_check.csv, spec71_sample_sets.csv, spec71_sample_set_replicate_AUC.csv,
#          fig7A_sample_audit.{png,pdf}, fig7B_PCA_by_line_3versions.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages({ library(edgeR); library(Seurat) })
KOS <- ORDER; rz <- function(x) (x - median(x, na.rm = TRUE)) / (1.4826 * mad(x, constant = 1, na.rm = TRUE) + 1e-12)

fab <- load_full_audited_bulk(); meta <- fab$meta; cs <- fab$cs; S <- rownames(meta); stopifnot(length(S) == 76)
s70 <- rownames(load_bulk()$meta); s70 <- s70[s70 %in% S]; stopifnot(length(s70) == 70)

## ---------------------------------------------------------------- QC metrics
rpg <- t(sapply(S, function(s) { x <- read.delim(file.path(PROC, "counts", paste0(s, "_ReadsPerGene.out.tab")), header = FALSE)
  c(unmapped = x[1, 4], multi = x[2, 4], nofeature = x[3, 4], ambiguous = x[4, 4], assigned = sum(x[-(1:4), 4])) }))
lf <- t(sapply(S, function(s) { f <- file.path(PROC, "logs", paste0(s, "_Log.final.out")); if (!file.exists(f)) return(c(input = NA, uniq_pct = NA))
  l <- readLines(f); g <- function(p) as.numeric(sub("%", "", sub(".*\\|\\s*", "", grep(p, l, value = TRUE)[1])))
  c(input = g("Number of input reads"), uniq_pct = g("Uniquely mapped reads %")) }))
qc <- data.frame(sample = S, Sample_ID = meta$Sample_ID, line = meta$Cell_Line, line_as_named = meta$Cell_Line_original, KO = meta$Gene_Control, rep = meta$rep, sheet_flag = meta$flag,
                 in_S70 = S %in% s70, input_pairs_M = lf[, "input"] / 1e6, uniq_pct = lf[, "uniq_pct"], assigned_M = rpg[, "assigned"] / 1e6,
                 pct_assigned = 100 * rpg[, "assigned"] / rowSums(rpg), pct_unmapped = 100 * rpg[, "unmapped"] / rowSums(rpg), pct_nofeature = 100 * rpg[, "nofeature"] / rowSums(rpg),
                 genes_detected = colSums(cs > 0), pct_mito = 100 * colSums(cs[grepl("^MT-", rownames(cs)), ]) / colSums(cs),
                 top100_share = apply(cs, 2, function(v) sum(sort(v, decreasing = TRUE)[1:100]) / sum(v)), stringsAsFactors = FALSE)

## ---------------------------------------------------------------- expression: TMM logCPM, cell-cycle index
obj <- CreateSeuratObject(counts = cs, min.cells = 0, min.features = 0); obj <- NormalizeData(obj, verbose = FALSE)
obj <- CellCycleScoring(obj, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes, set.ident = FALSE)
qc$cc <- as.numeric(obj$S.Score + obj$G2M.Score)
resid_on <- function(y, X) { X <- cbind(1, scale(X, scale = FALSE)); y - t(X[, -1, drop = FALSE] %*% solve(crossprod(X), t(y %*% X))[-1, , drop = FALSE]) }   # genes x samples, keeps gene means
line_expr <- function(cl, ids = S) { s <- ids[meta[ids, "Cell_Line"] == cl]; d <- DGEList(cs[, s]); d <- d[filterByExpr(d, group = meta[s, "Gene_Control"]), , keep.lib.sizes = FALSE]
  cpm(calcNormFactors(d), log = TRUE, prior.count = 2) }
pca_of <- function(y, ntop = 5000) { v <- matrixStats::rowVars(y); yt <- y[order(v, decreasing = TRUE)[seq_len(min(ntop, nrow(y)))], ]; yt <- yt[matrixStats::rowSds(yt) > 0, ]
  p <- prcomp(t(yt), center = TRUE, scale. = TRUE); list(x = p$x, ve = p$sdev^2 / sum(p$sdev^2)) }
PCA <- list(); r2rows <- list(); qc$med_spearman <- NA; qc$pc_dist <- NA
for (cl in LINES) { y <- line_expr(cl); s <- colnames(y); q <- qc[match(s, qc$sample), ]
  V <- list(raw = y, cc = resid_on(y, cbind(q$cc)), cc_depth = resid_on(y, cbind(q$cc, log10(q$assigned_M))))
  PCA[[cl]] <- lapply(V, pca_of)
  for (v in names(V)) for (pc in 1:5) { z <- PCA[[cl]][[v]]$x[, pc]; r2 <- function(f) summary(lm(z ~ f))$r.squared
    r2rows[[length(r2rows) + 1]] <- data.frame(line = cl, version = v, PC = pc, var_explained = PCA[[cl]][[v]]$ve[pc], KO = r2(factor(q$KO)), NTC_vs_KO = r2(q$KO == "NTC"), replicate = r2(factor(q$rep)),
                                               log_depth = r2(log10(q$assigned_M)), pct_assigned = r2(q$pct_assigned), cc = r2(q$cc)) }
  yc <- V$cc; top <- order(matrixStats::rowVars(yc), decreasing = TRUE)[1:3000]; C <- cor(yc[top, ] - rowMeans(yc[top, ]), method = "spearman"); diag(C) <- NA
  qc$med_spearman[match(s, qc$sample)] <- apply(C, 1, median, na.rm = TRUE)
  X <- PCA[[cl]]$cc$x[, 1:5]; X <- sweep(X, 2, apply(X, 2, median)); X <- sweep(X, 2, apply(X, 2, mad) + 1e-12, "/"); qc$pc_dist[match(s, qc$sample)] <- sqrt(rowSums(X^2)) }
R2 <- do.call(rbind, r2rows); write.csv(R2, file.path(OUT, "spec71_PC_covariate_R2.csv"), row.names = FALSE)

## ---------------------------------------------------------------- flags
byl <- function(x, bad_low = TRUE) { z <- ave(x, qc$line, FUN = rz); if (bad_low) z < -3 else z > 3 }
qc$low_depth <- qc$assigned_M < 5
qc$qc_outlier <- byl(log10(qc$assigned_M)) | byl(qc$pct_assigned) | byl(qc$genes_detected) | byl(qc$pct_mito, FALSE) | byl(qc$top100_share, FALSE)
qc$qc_reason <- apply(cbind(depth = byl(log10(qc$assigned_M)), pct_assigned = byl(qc$pct_assigned), genes = byl(qc$genes_detected), mito = byl(qc$pct_mito, FALSE), top100 = byl(qc$top100_share, FALSE)), 1,
                      function(r) paste(names(r)[r], collapse = "+"))
qc$expr_outlier <- byl(qc$med_spearman) | byl(qc$pc_dist, FALSE)
qc$mixed <- grepl("^MIXED", qc$sample)
qc$visual_outlier_old <- qc$Sample_ID %in% VISUAL_OUTLIERS | qc$sample %in% VISUAL_OUTLIERS
qc$in_SQC <- !(qc$low_depth | qc$qc_outlier | qc$expr_outlier | qc$mixed)

## ---------------------------------------------------------------- KO-label check: per-line z of each target gene, every sample scored against all 8 targets
lcpm_all <- cpm(calcNormFactors(DGEList(cs)), log = TRUE, prior.count = 2)[KOS, ]
Z <- lcpm_all; for (cl in LINES) { i <- meta$Cell_Line == cl; Z[, i] <- t(scale(t(lcpm_all[, i]))) }
kl <- do.call(rbind, lapply(S[meta$Gene_Control != "NTC"], function(s) { k <- meta[s, "Gene_Control"]; z <- Z[, s]
  data.frame(sample = s, line = meta[s, "Cell_Line"], KO = k, own_z = z[k], own_rank = rank(z)[k], lowest_target = names(which.min(z)), lowest_z = min(z),
             own_logCPM = lcpm_all[k, s], NTC_logCPM = mean(lcpm_all[k, meta$Cell_Line == meta[s, "Cell_Line"] & meta$Gene_Control == "NTC"])) }))
kl$log2FC_vs_NTC <- kl$own_logCPM - kl$NTC_logCPM; kl$label_supported <- kl$own_rank <= 2
write.csv(kl, file.path(OUT, "spec71_KO_label_check.csv"), row.names = FALSE)
qc$KO_label_supported <- kl$label_supported[match(qc$sample, kl$sample)]; qc$own_target_log2FC <- kl$log2FC_vs_NTC[match(qc$sample, kl$sample)]
# NTC samples: is any target unusually low (an NTC that is really a KO)?
ntc <- S[meta$Gene_Control == "NTC"]; qc$NTC_lowest_target <- NA; qc$NTC_lowest_z <- NA
qc$NTC_lowest_target[match(ntc, qc$sample)] <- apply(Z[, ntc, drop = FALSE], 2, function(z) names(which.min(z))); qc$NTC_lowest_z[match(ntc, qc$sample)] <- apply(Z[, ntc, drop = FALSE], 2, min)
write.csv(qc, file.path(OUT, "spec71_sample_audit.csv"), row.names = FALSE)

## ---------------------------------------------------------------- sample sets + replicate-structure criterion (no cross-line information used)
sets <- list(S70 = s70, S76 = S, SQC = qc$sample[qc$in_SQC], SQC_in70 = qc$sample[qc$in_SQC & qc$in_S70])
write.csv(do.call(rbind, lapply(names(sets), function(n) data.frame(set = n, sample = sets[[n]]))), file.path(OUT, "spec71_sample_sets.csv"), row.names = FALSE)
auc <- function(a, b) { if (!length(a) || !length(b)) return(NA); mean(outer(a, b, ">") + 0.5 * outer(a, b, "==")) }
rep_auc <- do.call(rbind, lapply(names(sets), function(n) do.call(rbind, lapply(LINES, function(cl) { ids <- sets[[n]]; if (sum(meta[ids, "Cell_Line"] == cl) < 4) return(NULL)
  y <- line_expr(cl, ids); q <- qc[match(colnames(y), qc$sample), ]; y <- resid_on(y, cbind(q$cc)); top <- order(matrixStats::rowVars(y), decreasing = TRUE)[1:3000]
  C <- cor(y[top, ] - rowMeans(y[top, ]), method = "spearman"); ut <- upper.tri(C); same <- outer(q$KO, q$KO, "==")[ut]; isn <- (outer(q$KO == "NTC", q$KO == "NTC", "|"))[ut]; r <- C[ut]
  data.frame(set = n, line = cl, n_samples = ncol(y), n_NTC = sum(q$KO == "NTC"), n_same_pairs = sum(same), median_same = median(r[same]), median_diff = median(r[!same]), AUC_all = auc(r[same], r[!same]),
             AUC_KO_only = auc(r[same & !isn], r[!same & !isn])) }))))
write.csv(rep_auc, file.path(OUT, "spec71_sample_set_replicate_AUC.csv"), row.names = FALSE)

options(width = 250); cat("\n=== flagged samples (any rule) ===\n")
print(qc[qc$low_depth | qc$qc_outlier | qc$expr_outlier | qc$mixed | !qc$in_S70, c("sample", "line", "KO", "in_S70", "assigned_M", "pct_assigned", "genes_detected", "med_spearman", "pc_dist", "low_depth", "qc_outlier", "qc_reason", "expr_outlier", "visual_outlier_old", "in_SQC")], digits = 3, row.names = FALSE)
cat("\n=== KO samples whose own target is not among the 2 lowest targets ===\n"); print(kl[!kl$label_supported, ], digits = 3, row.names = FALSE)
cat("\n=== NTC samples: lowest target z ===\n"); print(qc[qc$KO == "NTC", c("sample", "line", "NTC_lowest_target", "NTC_lowest_z")], digits = 3, row.names = FALSE)
cat("\n=== R2 of PC1-3 with covariates (cc-regressed PCA) ===\n"); print(R2[R2$version == "cc" & R2$PC <= 3, ], digits = 2, row.names = FALSE)
cat("\n=== replicate structure per sample set ===\n"); print(rep_auc, digits = 3, row.names = FALSE)
cat("\nset sizes:", paste(names(sets), sapply(sets, length), collapse = " | "), "\n")

## ---------------------------------------------------------------- figures
flag_of <- function(q) ifelse(q$mixed, "mixed", ifelse(q$low_depth | q$qc_outlier, "QC", ifelse(q$expr_outlier, "expression", ifelse(!q$in_S70, "not in S70", ""))))
draw_pca <- function(cl, v, lab_all = FALSE) { r <- PCA[[cl]][[v]]; X <- r$x[, 1:2]; q <- qc[match(rownames(X), qc$sample), ]; fl <- flag_of(q)
  rg <- apply(X, 2, range); pad <- 0.15 * apply(rg, 2, diff)
  plot(X, xlim = rg[, 1] + c(-1, 1) * pad[1], ylim = rg[, 2] + c(-1, 1) * pad[2], type = "n", bty = "n", las = 1, cex.axis = 0.7, cex.lab = 0.8,
       xlab = sprintf("PC1 (%.0f%%)", 100 * r$ve[1]), ylab = sprintf("PC2 (%.0f%%)", 100 * r$ve[2])); abline(h = 0, v = 0, col = GRID)
  points(X, pch = ifelse(fl == "", 21, 24), bg = KO_COL[q$KO], col = ifelse(fl == "", "white", INK), cex = ifelse(fl == "", 1.6, 1.8), lwd = ifelse(fl == "", 1, 1.4))
  i <- if (lab_all) rep(TRUE, nrow(X)) else fl != "" | q$KO == "NTC"
  if (any(i)) text(X[i, 1], X[i, 2], paste(q$KO, q$rep)[i], pos = 3, cex = 0.5, offset = 0.4, col = INK2)
  mtext(cl, side = 3, adj = 0, line = 0.4, font = 2, cex = 0.8) }

save_fig("fig7A_sample_audit", 15, 9.2, function() {
  layout(rbind(rep(1, 4), rep(2, 4), 3:6), heights = c(1, 0.8, 1.35)); par(oma = c(0, 0, 2.6, 0.5), family = "sans")
  o <- order(match(qc$line, LINES), match(qc$KO, c("NTC", KOS)), qc$rep); q <- qc[o, ]; fl <- flag_of(q); x <- seq_len(nrow(q)) + 1.5 * (match(q$line, LINES) - 1)
  par(mar = c(0.6, 4.4, 2, 0.5)); plot(NA, xlim = range(x) + c(-1, 1), ylim = c(0.3, 200), log = "y", xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "reads assigned to genes (M)", cex.axis = 0.75, xaxs = "i")
  abline(h = c(1, 10, 100), col = GRID); abline(h = 5, col = INK2, lty = 2); rect(x - 0.38, 0.3, x + 0.38, q$assigned_M, col = KO_COL[q$KO], border = NA)
  points(x[fl != ""], q$assigned_M[fl != ""] * 1.35, pch = 25, bg = INK, col = INK, cex = 0.7)
  for (cl in LINES) text(mean(x[q$line == cl]), 190, sprintf("%s (n=%d)", cl, sum(q$line == cl)), cex = 0.85, font = 2, col = INK, xpd = NA)
  legend("bottomright", c("NTC", KOS), fill = KO_COL[c("NTC", KOS)], border = NA, bty = "n", ncol = 9, cex = 0.65, inset = c(0, -0.02)); mtext("A  read depth (dashed = 5 M; triangle = flagged)", side = 3, adj = 0, line = 0.6, font = 2, cex = 0.75)
  par(mar = c(6.4, 4.4, 1.6, 0.5)); plot(NA, xlim = range(x) + c(-1, 1), ylim = c(20, 100), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "% reads assigned", cex.axis = 0.75, xaxs = "i"); abline(h = c(40, 60, 80), col = GRID)
  points(x, q$pct_assigned, pch = ifelse(fl == "", 21, 24), bg = KO_COL[q$KO], col = ifelse(fl == "", "white", INK), cex = 1.3)
  text(x, 14, paste0(q$sample, ifelse(fl == "", "", paste0("  [", fl, "]"))), srt = 90, adj = 1, xpd = NA, cex = 0.42, col = ifelse(fl == "", INK2, INK), font = ifelse(fl == "", 1, 2))
  mtext("B  % of reads assigned to genes", side = 3, adj = 0, line = 0.3, font = 2, cex = 0.75)
  par(mar = c(4, 4.2, 2.2, 1)); for (cl in LINES) draw_pca(cl, "cc")
  mtext(sprintf("Audit of the %d bulk RNA-seq samples (NTC + 8 CROP-seq KOs, 4 lines; corrected cell line)", nrow(qc)), outer = TRUE, side = 3, adj = 0, line = 1.2, font = 2, cex = 1)
  mtext("C (bottom) PCA per line: TMM logCPM, top 5000 variable genes, cell cycle regressed; colour = KO; triangle = flagged (QC / expression outlier / mixed / not in the 70-sample set); NTC and flagged samples labelled", outer = TRUE, side = 3, adj = 0, line = 0.1, cex = 0.6, col = INK2)
})
save_fig("fig7B_PCA_by_line_3versions", 15, 11.5, function() {
  layout(matrix(1:12, nrow = 3, byrow = TRUE)); par(oma = c(0, 2.2, 2.6, 0.5), mar = c(4, 4.2, 2.2, 1), family = "sans")
  for (v in c("raw", "cc", "cc_depth")) for (cl in LINES) { draw_pca(cl, v, lab_all = TRUE); if (cl == LINES[1]) mtext(c(raw = "raw", cc = "cell cycle regressed", cc_depth = "cell cycle + depth regressed")[v], side = 2, line = 3.6, font = 2, cex = 0.75) }
  mtext("PCA per cell line, 3 versions (TMM logCPM, top 5000 variable genes); colour = KO, triangle = flagged", outer = TRUE, side = 3, adj = 0, line = 1, font = 2, cex = 1)
})
