# Spec 55: does adding the 3 QC-dropped H2052 NTC samples (H28_NTC_R1/2/3: sequenced as "H28", genotyped NCI-H2052, PCA outliers with PC2>20) change the
# CROP-seq screen vs bulk NCI-H2052 agreement (same cell line)? Spearman, NTC - KO, per KO: own KO (screen KO vs H2052 bulk same KO) vs null (screen KO vs the
# OTHER H2052 bulk KOs), for the 20 Cm modules and the top 1000/2000/3000 variable genes, with H2052 NTC = 1 sample (curated sheet) vs 4 samples (+H28_NTC_R1/2/3).
# Also a sanity check that the 3 added samples look like H2052 (Spearman of their log-expression vs H2052_NTC_R3 and vs each other line's NTC).
# needs: git_repo/files/scRNA_PM_modules.rds, spec_common.R (load_screen()), repro_common.R (load_bulk(), load_full_audited_bulk())
# outputs: spec55_NTC_concordance.csv, spec55_screen_vs_H2052bulk_NTC1_vs_NTC4.csv, fig3M_screen_vs_H2052bulk_NTC1_vs_NTC4.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3")
SGN <- -1; TOPNS <- c(1000, 2000, 3000)   # NTC - KO
b <- load_bulk(); meta0 <- b$meta; lg0 <- log2(b$cp10k + 1); a <- load_full_audited_bulk()
add <- rownames(a$meta)[a$meta$Cell_Line == "NCI-H2052" & a$meta$Gene_Control == "NTC" & !rownames(a$meta) %in% rownames(meta0)]
g <- intersect(rownames(lg0), rownames(a$lg))
meta1 <- rbind(meta0[, c("Sample_ID", "Cell_Line", "Gene_Control")], a$meta[add, c("Sample_ID", "Cell_Line", "Gene_Control")]); lg1 <- cbind(lg0[g, ], a$lg[g, add])[, rownames(meta1)]
cnt1 <- cbind(b$counts[g, ], a$cs[g, add])[, rownames(meta1)]; cat("added:", paste(add, collapse = ", "), "\n")

## ---- sanity: do the added samples look like H2052? (top 3000 variable genes across all samples; Spearman on log2 CP10K)
vg <- names(sort(apply(lg1, 1, var), decreasing = TRUE))[1:3000]; ntc <- rownames(meta1)[meta1$Gene_Control == "NTC"]
cc <- cor(lg1[vg, ntc], method = "spearman"); ntc_line <- setNames(meta1[ntc, "Cell_Line"], ntc)
conc <- do.call(rbind, lapply(add, function(s) data.frame(sample = s, vs_H2052_NTC_R3 = cc[s, "H2052_NTC_R3"], vs_MSTO_NTC_mean = mean(cc[s, ntc[ntc_line == "MSTO-211H"]]),
       vs_H2452_NTC_mean = mean(cc[s, ntc[ntc_line == "NCI-H2452"]]), vs_H28_NTC_mean = mean(cc[s, ntc[ntc_line == "NCI-H28"]]), vs_other_added_mean = mean(cc[s, setdiff(add, s)]))))
write.csv(conc, file.path(OUT, "spec55_NTC_concordance.csv"), row.names = FALSE); cat("\nSpearman of added H2052 NTC samples vs NTCs of each line (top3000 vg):\n"); print(format(conc, digits = 3), row.names = FALSE)

## ---- screen profiles (NTC - KO, per KO) in the same feature spaces
cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds"); top50 <- lapply(cm_genes, head, 50); CM <- names(top50)
srt <- suppressWarnings(AddModuleScore(load_screen(), features = top50, name = "Cm_", seed = 1)); sgrp <- as.character(srt$merged_call)
scm <- as.matrix(srt@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(scm) <- CM; Y <- GetAssayData(srt, layer = "data")
KOS <- intersect(ORDER, unique(sgrp)); KOS <- KOS[KOS %in% meta0$Gene_Control[meta0$Cell_Line == "NCI-H2052"]]
prof <- function(M, grp, ko_names) t(sapply(ko_names, function(k) SGN * (apply(M[grp == k, , drop = FALSE], 2, median) - apply(M[grp == "NTC", , drop = FALSE], 2, median))))
lv <- sapply(LINES, function(cl) matrixStats::rowVars(as.matrix(lg1[, rownames(meta1)[meta1$Cell_Line == cl], drop = FALSE]))); vord <- intersect(rownames(lg1)[order(rowMeans(lv), decreasing = TRUE)], rownames(Y))
spaces <- c("20 Cm modules", sprintf("top %d variable genes", TOPNS))

bulk_prof <- function(meta, lg, sp, cnt) {   # H2052 bulk NTC-KO profile (KO x feature) for one sample set
  ids <- rownames(meta)[meta$Cell_Line == "NCI-H2052"]; grp <- meta[ids, "Gene_Control"]
  if (sp == "20 Cm modules") { sb <- CreateSeuratObject(counts = as.matrix(cnt[, rownames(meta)])); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg[, rownames(meta)])); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
    M <- as.matrix(sb@meta.data[ids, paste0("Cm_", seq_along(CM))]); colnames(M) <- CM }
  else { n <- as.integer(sub("top (\\d+) .*", "\\1", sp)); M <- t(lg[vord[1:n], ids]) }
  prof(M, grp, KOS) }
res <- list()
for (sp in spaces) {
  S <- if (sp == "20 Cm modules") prof(scm, sgrp, KOS) else { n <- as.integer(sub("top (\\d+) .*", "\\1", sp)); prof(t(as.matrix(Y[vord[1:n], ])), sgrp, KOS) }
  for (v in c("NTC=1", "NTC=4")) { m <- if (v == "NTC=1") meta0 else meta1; l <- if (v == "NTC=1") lg0 else lg1; cn <- if (v == "NTC=1") b$counts else cnt1
    B <- bulk_prof(m, l, sp, cn); S2 <- S[, colnames(B), drop = FALSE]
    for (k in KOS) res[[length(res) + 1]] <- data.frame(space = sp, H2052_NTC = v, KO = k, own = suppressWarnings(cor(S2[k, ], B[k, ], method = "spearman")),
        null_median = median(sapply(setdiff(KOS, k), function(k2) suppressWarnings(cor(S2[k, ], B[k2, ], method = "spearman"))))) } }
res <- do.call(rbind, res); res$separation <- res$own - res$null_median; write.csv(res, file.path(OUT, "spec55_screen_vs_H2052bulk_NTC1_vs_NTC4.csv"), row.names = FALSE)
options(width = 200)
for (sp in spaces) { cat(sprintf("\n== %s: screen vs H2052 bulk (Spearman, NTC-KO) ==\n", sp)); w <- reshape(res[res$space == sp, c("KO", "H2052_NTC", "own", "null_median")], idvar = "KO", timevar = "H2052_NTC", direction = "wide"); print(format(w, digits = 2), row.names = FALSE)
  r <- res[res$space == sp, ]; for (v in c("NTC=1", "NTC=4")) cat(sprintf("  %s: median own %.3f | median null %.3f | separation %+.3f\n", v, median(r$own[r$H2052_NTC == v]), median(r$null_median[r$H2052_NTC == v]), median(r$own[r$H2052_NTC == v] - r$null_median[r$H2052_NTC == v]))) }

save_fig("fig3M_screen_vs_H2052bulk_NTC1_vs_NTC4", 13, 8.6, function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.2, 0), family = "sans")
  for (sp in spaces) { r <- res[res$space == sp, ]; kos <- r$KO[r$H2052_NTC == "NTC=4"][order(-r$own[r$H2052_NTC == "NTC=4"])]; par(mar = c(5, 4.4, 2.6, 0.8))
    plot(NA, xlim = c(0.5, length(kos) + 0.5), ylim = c(-1, 1), xaxt = "n", xlab = "", ylab = "Spearman rho (screen vs H2052 bulk, NTC-KO)")
    axis(1, at = seq_along(kos), labels = FALSE); text(seq_along(kos), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), kos, srt = 45, adj = 1, xpd = NA, cex = 0.85); abline(h = 0, lty = 2, col = INK)
    for (i in seq_along(kos)) { o1 <- r[r$KO == kos[i] & r$H2052_NTC == "NTC=1", ]; o4 <- r[r$KO == kos[i] & r$H2052_NTC == "NTC=4", ]
      segments(i - 0.12, o1$own, i + 0.12, o4$own, col = "grey60"); points(i - 0.12, o1$own, pch = 1, cex = 1.5, col = KO_COL[[kos[i]]], lwd = 2); points(i + 0.12, o4$own, pch = 16, cex = 1.5, col = KO_COL[[kos[i]]])
      points(i - 0.12, o1$null_median, pch = 4, cex = 1, col = "grey40"); points(i + 0.12, o4$null_median, pch = 3, cex = 1, col = "grey40") }
    if (sp == spaces[1]) legend("bottomleft", legend = c("own KO, H2052 NTC = 1", "own KO, H2052 NTC = 4", "null median, NTC = 1", "null median, NTC = 4"), pch = c(1, 16, 4, 3), col = c("black", "black", "grey40", "grey40"), bty = "n", cex = 0.65)
    mtext(sp, side = 3, adj = 0, line = 0.4, font = 2, cex = 0.9) }
  mtext("CROP-seq screen vs bulk NCI-H2052 (same line): H2052 NTC = 1 sample vs 4 samples (+H28_NTC_R1/2/3), NTC-KO, Spearman", outer = TRUE, side = 3, adj = 0, line = 1.4, font = 2, cex = 0.95)
})
