# Spec 123: why the NCI-H2052 KO x Cm heatmap stops resembling the CROP-seq one when the clean H2052 controls are used (spec122).
# old control = H2052_NTC_R3 (carries BPTF guides; same i5 index group as the H2052 KO samples); new controls = H28_NTC_R1-R3 (H2052 genotype, no guides; i5 index groups shared with H28 / MSTO / liver-line samples)
# shared response of a data set = control minus the mean of the KO samples (8 KOs), per Cm module (AddModuleScore, top-50 genes); screen: spec81 (mean over the 8 KOs of NTC - KO)
#   A / B  screen shared response vs bulk H2052 shared response with the old / new control(s) (20 modules)
#   C      cell-cycle score (S.Score + G2M.Score) of the H2052 samples: old control, new controls, KO samples
#   D      Spearman of the bulk shared response with the screen's, old vs new control, before and after regressing the cell-cycle score out of both data sets; and old minus new control vs the screen
#   E      PCA of the H2052 samples (top 5000 variable genes), controls marked, colour = i5 index group
#   F      liver-gene signal (ALB, APOA1, APOA2, APOB, FGB, FGG, TF, TTR, SERPINA1, AFP; counts per million) by i5 index group: carry-over from the liver-line samples that share the index
# also printed: the same shared-response comparison for the other 3 lines (their own controls)
# needs: spec72 (S70), spec75 (screen scores), spec81, sample sheet (i5 = index2), agree_common.R::load_bulk_S91
# outputs: spec123_shared_response.csv, spec123_H2052_samples.csv, fig40_H2052_NTC_old_vs_new.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R")); suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3", width = 220)
CL <- "NCI-H2052"; pr <- readRDS(file.path(OUT, sprintf("spec72_profiles_S70%s.rds", DSUF))); B <- load_bulk_S91(pr$samples, extra_kos = FALSE); m <- B$meta; ids <- B$ids; keep <- pr$genes; top50 <- cm_top50(); CM <- names(top50)
sp <- function(a, b) { r <- suppressWarnings(cor.test(a, b, method = "spearman", exact = FALSE)); c(rho = unname(r$estimate), p = r$p.value) }; fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 0), formatC(p, format = "g", digits = 2))
y <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(B$cs[, ids])), log = TRUE, prior.count = 2); yy <- y[rowSums(B$cs[, ids]) > 0, ]; yy <- yy - min(yy)
sb <- SetAssayData(CreateSeuratObject(counts = B$cs[rownames(yy), ids]), slot = "data", new.data = as.matrix(yy)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1)); SCO <- t(as.matrix(sb@meta.data[ids, paste0("Cm_", 1:20)])); rownames(SCO) <- CM
h <- ids[m$Cell_Line == CL]; old <- S91_NTC_AS_BPTF; new <- S91_H2052_NTC; ko <- setdiff(h, c(old, new)); kom <- rowMeans(SCO[, ko])
x81 <- read.csv(file.path(OUT, "spec81_Cm_KO_NTCminusKO_screen_and_bulk.csv")); d <- x81[x81$source == "screen", ]; scr <- colMeans(tapply(d$NTC_minus_KO, list(d$KO, d$Cm), mean)[KOS, CM])
sh_old <- SCO[, old] - kom; sh_new <- rowMeans(SCO[, new]) - kom
sc <- readRDS(file.path(OUT, "spec75_screen_scores_by_cell.rds")); g <- as.character(sc$group); Sr <- apply(sc$score[, CM], 2, function(v) resid(lm(v ~ sc$cc))); scr_cc <- colMeans(Sr[g == "NTC", ]) - rowMeans(sapply(KOS, function(k) colMeans(Sr[g == k, , drop = FALSE])))
Br <- t(apply(SCO[, h], 1, function(v) resid(lm(v ~ B$cc[h])))); colnames(Br) <- h; sh_old_cc <- Br[, old] - rowMeans(Br[, ko]); sh_new_cc <- rowMeans(Br[, new]) - rowMeans(Br[, ko])
RES <- rbind(data.frame(comparison = "old control, as is", t(sp(scr, sh_old))), data.frame(comparison = "new controls, as is", t(sp(scr, sh_new))), data.frame(comparison = "old control, cell cycle regressed", t(sp(scr_cc, sh_old_cc))),
             data.frame(comparison = "new controls, cell cycle regressed", t(sp(scr_cc, sh_new_cc))), data.frame(comparison = "old control minus new controls", t(sp(scr, SCO[, old] - rowMeans(SCO[, new])))))
for (cl in setdiff(LINES, CL)) { s <- ids[m$Cell_Line == cl]; n <- s[m[s, "Gene_Control"] == "NTC"]; RES <- rbind(RES, data.frame(comparison = sprintf("%s (its %d controls)", cl, length(n)), t(sp(scr, rowMeans(SCO[, n, drop = FALSE]) - rowMeans(SCO[, setdiff(s, n)])))) )}
write.csv(RES, file.path(OUT, "spec123_shared_response.csv"), row.names = FALSE); cat("=== shared response (control - mean of KO samples, 20 modules) vs the CROP-seq shared response ===\n"); print(RES, digits = 2, row.names = FALSE)
ccS <- apply(sc$score[, CM], 2, function(v) cor(v, sc$cc)); cat("screen shared response vs how strongly each module follows the cell cycle in the cells: rho", round(sp(scr, ccS)[1], 2), "| most cell-cycle module:", names(which.max(ccS)), sprintf("(r %.2f)", max(ccS)), "\n")
# sample sheet: i5 index
ss <- readLines(file.path(PROC, "sample_sheet/SampleSheet.csv")); ss <- ss[(grep("^\\[BCLConvert_Data\\]", ss) + 2):length(ss)]; ss <- do.call(rbind, strsplit(ss[nzchar(ss)], ",")); i5 <- setNames(ss[, 4], ss[, 2])
HS <- data.frame(sample = h, role = ifelse(h == old, "old control", ifelse(h %in% new, "new control", "KO")), KO_label = m[h, "label_original"], i5 = i5[h], cc = B$cc[h]); yk <- normalise(B$cs[keep, h], "tmm")$y
vg <- keep[order(matrixStats::rowVars(yk), decreasing = TRUE)][1:5000]; pc <- prcomp(t(yk[vg, ] - rowMeans(yk[vg, ]))); HS$PC1 <- pc$x[, 1]; HS$PC2 <- pc$x[, 2]; ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2)); write.csv(HS, file.path(OUT, "spec123_H2052_samples.csv"), row.names = FALSE)
cat("\ncell-cycle score: old control", round(B$cc[old], 2), "| new controls", round(B$cc[new], 2), "| KO samples: mean", round(mean(B$cc[ko]), 2), ", max", round(max(B$cc[ko]), 2), "(rank of the old control among the 20 H2052 samples:", rank(-B$cc[h])[old], ")\n")
cen <- yk[vg, ] - rowMeans(yk[vg, ko]); Cc <- cor(cen, method = "spearman"); cat("similarity among the 3 new controls (genes centred on the KO mean):", round(Cc[new, new][upper.tri(Cc[new, new])], 2), "| new control vs KO samples: mean", round(mean(Cc[new, ko]), 2), "| KO vs KO: median", round(median(Cc[ko, ko][upper.tri(Cc[ko, ko])]), 2), "\n")
cat("i5 index groups: H2052 KO samples", paste(names(table(i5[ko])), table(i5[ko]), collapse = ", "), "| old control", i5[[old]], "| new controls", i5[new], "\n")
# liver-gene carry-over by i5 group, all samples of the count matrix
LIV <- c("ALB", "APOA1", "APOA2", "APOB", "FGB", "FGG", "TF", "TTR", "SERPINA1", "AFP"); ct <- read.delim(file.path(PROC, "counts/gene_counts_matrix_stranded_reverse.tsv"), row.names = 1, check.names = FALSE); sy <- read.delim(file.path(PROC, "counts/gene_id_to_symbol.tsv"), header = FALSE)
liv <- 1e6 * colSums(ct[rownames(ct) %in% sy[sy[, 2] %in% LIV, 1], ]) / colSums(ct); LV <- data.frame(sample = names(liv), cpm = liv, i5 = i5[names(liv)], meso = grepl("^(H2052|MSTO|H2452|H28)_", names(liv)) & !grepl("^MIXED", names(liv)))
cat("\nliver-gene counts per million, mesothelioma samples by i5 group (median):\n"); print(round(tapply(LV$cpm[LV$meso], LV$i5[LV$meso], median), 1)); cat("liver-line / mixed samples per i5 group:\n"); print(table(LV$i5[!LV$meso & !grepl("NR4A2|SCRAMBLE", LV$sample)]))
cat("old control", round(liv[old], 1), "| new controls", round(liv[new], 1), "| H2052 KO samples median", round(median(liv[ko]), 1), "\n")

I5C <- setNames(c("#2a78d6", "#1baf7a", "#eb6834", "#eda100", "#4a3aa7"), sort(unique(LV$i5[LV$meso])))
save_fig("fig40_H2052_NTC_old_vs_new", 16, 10, function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 4, 0.4), family = "sans", mar = c(4.8, 4.8, 3.6, 1)); ttl <- function(a, b) { mtext(a, side = 3, adj = 0, line = 1.6, font = 2, cex = 0.8); mtext(b, side = 3, adj = 0, line = 0.5, cex = 0.58, col = INK2) }
  scat <- function(b, lab, r, col) { plot(scr, b, pch = 21, bg = col, col = INK2, cex = 1.6, bty = "n", las = 1, xlab = "CROP-seq: NTC minus mean of the 8 KOs, per module", ylab = sprintf("bulk H2052: %s minus mean of the KO samples", lab), cex.axis = 0.8); abline(h = 0, v = 0, lty = 2, col = GRID); abline(lm(b ~ scr), col = INK2, lty = 3); text(scr, b, CM, pos = 3, cex = 0.55, col = INK2, xpd = NA) }
  scat(sh_old, "old control", RES[1, ], "grey45"); ttl("A  old control (H2052_NTC_R3, carries BPTF guides)", sprintf("20 modules: Spearman %.2f, p = %s", RES$rho[1], fp(RES$p[1])))
  scat(sh_new, "new controls", RES[2, ], "#2a78d6"); ttl("B  new controls (H28_NTC_R1-R3, no guides)", sprintf("20 modules: Spearman %.2f, p = %s", RES$rho[2], fp(RES$p[2])))
  rc <- c(`old control` = "grey20", `new control` = "#2a78d6", KO = "#d5d4cf"); xx <- c(`old control` = 1, `new control` = 2, KO = 3)[HS$role] + ifelse(HS$role == "KO", seq(-0.25, 0.25, length.out = sum(HS$role == "KO"))[cumsum(HS$role == "KO")], ifelse(HS$role == "new control", c(-0.12, 0, 0.12)[cumsum(HS$role == "new control")], 0))
  plot(xx, HS$cc, pch = 21, bg = rc[HS$role], col = INK2, cex = 1.7, xlim = c(0.5, 3.5), xaxt = "n", bty = "n", las = 1, xlab = "", ylab = "cell-cycle score (S.Score + G2M.Score)", cex.axis = 0.8); axis(1, at = 1:3, labels = c("old control", "new controls", "KO samples"), tick = FALSE, cex.axis = 0.85); abline(h = mean(B$cc[ko]), lty = 3, col = MUTED)
  ttl("C  the old control is the most proliferative H2052 sample", sprintf("dotted = mean of the KO samples; the screen's shared response follows the cell cycle (module-level rho %.2f)", sp(scr, ccS)[1]))
  bp <- barplot(RES$rho[1:5], col = c("grey45", "#2a78d6", "grey45", "#2a78d6", "#eda100"), border = NA, las = 1, names.arg = rep("", 5), ylab = "Spearman with the CROP-seq shared response (20 modules)", cex.axis = 0.8, ylim = range(c(0, RES$rho[1:5])) + c(-0.15, 0.2)); abline(h = 0, col = INK)
  text(bp, pmax(RES$rho[1:5], 0) + 0.03, sprintf("%.2f\np %s", RES$rho[1:5], fp(RES$p[1:5])), cex = 0.62, adj = c(0.5, 0)); text(bp, par("usr")[3] - 0.03 * diff(par("usr")[3:4]), c("old control", "new controls", "old, cell cycle\nregressed", "new, cell cycle\nregressed", "old minus new\ncontrol"), srt = 30, adj = 1, xpd = NA, cex = 0.7)
  ttl("D  where the earlier similarity came from", "regressing the cell-cycle score out of both data sets; and the difference between the two controls themselves")
  plot(HS$PC1, HS$PC2, pch = ifelse(HS$role == "KO", 21, ifelse(HS$role == "old control", 23, 24)), bg = I5C[HS$i5], col = INK2, cex = ifelse(HS$role == "KO", 1.4, 2.1), bty = "n", las = 1, xlab = sprintf("PC1 (%d%%)", ve[1]), ylab = sprintf("PC2 (%d%%)", ve[2]), cex.axis = 0.8)
  w <- HS$role != "KO"; text(HS$PC1[w], HS$PC2[w], HS$sample[w], pos = 4, cex = 0.6, xpd = NA); legend("topright", c(names(I5C), "KO sample", "old control", "new control"), pch = c(rep(22, length(I5C)), 21, 23, 24), pt.bg = c(I5C, rep("white", 3)), bty = "n", cex = 0.58, title = "i5 index group")
  ttl("E  PCA of the H2052 samples", "top 5000 variable genes; the new controls sit together, away from the KO samples, and have other i5 index groups")
  LVm <- LV[LV$meso, ]; LVm$i5 <- factor(LVm$i5, levels = names(I5C)); boxplot(log10(cpm + 0.1) ~ i5, LVm, col = adjustcolor(I5C, 0.5), border = INK2, outline = FALSE, las = 2, xlab = "", ylab = "liver genes, log10 counts per million", cex.axis = 0.62, ylim = range(log10(LVm$cpm + 0.1)))
  points(jitter(as.integer(LVm$i5), 0.6), log10(LVm$cpm + 0.1), pch = 16, cex = 0.5, col = "grey40"); w <- match(c(old, new), LVm$sample); points(as.integer(LVm$i5[w]), log10(LVm$cpm[w] + 0.1), pch = c(23, 24, 24, 24), bg = c("grey20", rep("#2a78d6", 3)), cex = 1.7)
  nl <- table(factor(LV$i5[!LV$meso & !grepl("NR4A2|SCRAMBLE", LV$sample)], levels = names(I5C))); mtext(sprintf("%d liver", nl), side = 3, at = seq_along(nl), line = -1, cex = 0.55, col = INK2)
  ttl("F  trace liver transcripts follow the i5 index group", "mesothelioma samples; 'n liver' = liver-line or mixed samples sharing the index; diamond = old control, triangles = new controls")
  mtext("Why the NCI-H2052 bulk no longer resembles the CROP-seq shared response with the clean controls", outer = TRUE, side = 3, adj = 0, line = 2.4, font = 2, cex = 1)
  mtext("shared response = control minus the mean of the KO samples, per Cm module (AddModuleScore, top-50 genes)", outer = TRUE, side = 3, adj = 0, line = 1.2, cex = 0.65, col = INK2)
})
