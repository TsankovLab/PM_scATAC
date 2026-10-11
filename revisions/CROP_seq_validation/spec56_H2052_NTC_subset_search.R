# Spec 56: which subset of the 4 candidate NCI-H2052 NTC samples (H2052_NTC_R3 + the QC-dropped H28_NTC_R1/2/3) gives the best agreement between the CROP-seq screen and
# bulk H2052? All 15 non-empty subsets are scored (Spearman, NTC - KO, per KO: own = screen KO vs H2052 bulk same KO; null = screen KO vs the OTHER H2052 bulk KOs;
# separation = own - null) in 4 feature spaces (20 Cm modules, top 1000/2000/3000 variable genes); subsets ranked by mean separation across the spaces.
# SELECTION-BIAS CHECK: leave-one-KO-out -- pick the best subset using 7 KOs, evaluate on the held-out KO, compare with the fixed R3-only baseline and the all-4 set.
# Cm scores are computed once on the full 4-NTC object (AddModuleScore), then the NTC baseline (median) is taken over the chosen subset.
# needs: git_repo/files/scRNA_PM_modules.rds, spec_common.R, repro_common.R
# outputs: spec56_NTC_subset_scores.csv, spec56_NTC_subset_LOKO.csv, fig3N_H2052_NTC_subset_search.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3")
SGN <- -1; TOPNS <- c(1000, 2000, 3000)
b <- load_bulk(); meta0 <- b$meta; lg0 <- log2(b$cp10k + 1); a <- load_full_audited_bulk()
add <- rownames(a$meta)[a$meta$Cell_Line == "NCI-H2052" & a$meta$Gene_Control == "NTC" & !rownames(a$meta) %in% rownames(meta0)]
g <- intersect(rownames(lg0), rownames(a$lg)); meta1 <- rbind(meta0[, c("Sample_ID", "Cell_Line", "Gene_Control")], a$meta[add, c("Sample_ID", "Cell_Line", "Gene_Control")])
lg1 <- cbind(lg0[g, ], a$lg[g, add])[, rownames(meta1)]; cnt1 <- cbind(b$counts[g, ], a$cs[g, add])[, rownames(meta1)]
NTCS <- rownames(meta1)[meta1$Cell_Line == "NCI-H2052" & meta1$Gene_Control == "NTC"]; cat("candidate H2052 NTC samples:", paste(NTCS, collapse = ", "), "\n")
h_ko <- rownames(meta1)[meta1$Cell_Line == "NCI-H2052" & meta1$Gene_Control != "NTC"]

cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds"); top50 <- lapply(cm_genes, head, 50); CM <- names(top50)
srt <- suppressWarnings(AddModuleScore(load_screen(), features = top50, name = "Cm_", seed = 1)); sgrp <- as.character(srt$merged_call)
scm <- as.matrix(srt@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(scm) <- CM; Y <- GetAssayData(srt, layer = "data")
sb <- CreateSeuratObject(counts = as.matrix(cnt1)); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg1)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
bcm <- as.matrix(sb@meta.data[rownames(meta1), paste0("Cm_", seq_along(CM))]); colnames(bcm) <- CM
KOS <- intersect(ORDER, unique(sgrp)); KOS <- KOS[KOS %in% meta1$Gene_Control[meta1$Cell_Line == "NCI-H2052"]]
lv <- sapply(LINES, function(cl) matrixStats::rowVars(as.matrix(lg1[, rownames(meta1)[meta1$Cell_Line == cl], drop = FALSE]))); vord <- intersect(rownames(lg1)[order(rowMeans(lv), decreasing = TRUE)], rownames(Y))
spaces <- c("20 Cm modules", sprintf("top %d variable genes", TOPNS))
feat <- function(sp) if (sp == "20 Cm modules") list(S = scm, B = bcm) else { n <- as.integer(sub("top (\\d+) .*", "\\1", sp)); list(S = t(as.matrix(Y[vord[1:n], ])), B = t(lg1[vord[1:n], ])) }
med <- function(M, ids) apply(M[ids, , drop = FALSE], 2, median)
subsets <- unlist(lapply(seq_along(NTCS), function(k) combn(NTCS, k, simplify = FALSE)), recursive = FALSE); lab <- sapply(subsets, function(s) paste(sub("_NTC_", "", s), collapse = "+"))

rows <- list()
for (sp in spaces) { F <- feat(sp); Sd <- t(sapply(KOS, function(k) SGN * (med(F$S, which(sgrp == k)) - med(F$S, which(sgrp == "NTC")))))
  Bko <- t(sapply(KOS, function(k) med(F$B, rownames(meta1)[meta1$Cell_Line == "NCI-H2052" & meta1$Gene_Control == k])))
  for (i in seq_along(subsets)) { Bn <- med(F$B, subsets[[i]]); Bd <- SGN * sweep(Bko, 2, Bn, "-")
    for (k in KOS) rows[[length(rows) + 1]] <- data.frame(space = sp, subset = lab[i], n_ntc = length(subsets[[i]]), KO = k, own = suppressWarnings(cor(Sd[k, ], Bd[k, ], method = "spearman")),
        null = median(sapply(setdiff(KOS, k), function(k2) suppressWarnings(cor(Sd[k, ], Bd[k2, ], method = "spearman"))))) } }
R <- do.call(rbind, rows); R$sep <- R$own - R$null; write.csv(R, file.path(OUT, "spec56_NTC_subset_scores.csv"), row.names = FALSE)
sc <- aggregate(cbind(own, null, sep) ~ subset + n_ntc + space, R, median)
top <- do.call(rbind, lapply(split(sc, sc$subset), function(d) data.frame(subset = d$subset[1], n_ntc = d$n_ntc[1], mean_sep = mean(d$sep), mean_own = mean(d$own), mean_null = mean(d$null),
        sep_Cm = d$sep[d$space == "20 Cm modules"], sep_top1000 = d$sep[grepl("1000", d$space)], sep_top2000 = d$sep[grepl("2000", d$space)], sep_top3000 = d$sep[grepl("3000", d$space)])))
top <- top[order(-top$mean_sep), ]; options(width = 220); cat("\n=== all 15 H2052 NTC subsets, ranked by mean separation (own - null; median over KOs; mean over the 4 spaces) ===\n"); print(format(top, digits = 2), row.names = FALSE)

## ---- selection-bias check: leave-one-KO-out
per_ko <- aggregate(sep ~ subset + KO, R, mean)   # mean over spaces, per subset x KO
loko <- do.call(rbind, lapply(KOS, function(k) { tr <- aggregate(sep ~ subset, per_ko[per_ko$KO != k, ], mean); pick <- tr$subset[which.max(tr$sep)]; te <- per_ko[per_ko$KO == k, ]
  data.frame(held_out = k, picked = pick, heldout_sep_picked = te$sep[te$subset == pick], heldout_sep_R3only = te$sep[te$subset == "H2052+NA"][1], heldout_sep_all4 = te$sep[te$subset == lab[length(lab)]],
             heldout_sep_mean_all_subsets = mean(te$sep)) }))
r3 <- lab[sapply(subsets, function(s) identical(s, "H2052_NTC_R3"))]; loko$heldout_sep_R3only <- sapply(KOS, function(k) per_ko$sep[per_ko$KO == k & per_ko$subset == r3])
all4 <- lab[sapply(subsets, function(s) length(s) == 4)]; loko$heldout_sep_all4 <- sapply(KOS, function(k) per_ko$sep[per_ko$KO == k & per_ko$subset == all4])
write.csv(loko, file.path(OUT, "spec56_NTC_subset_LOKO.csv"), row.names = FALSE); cat("\n=== leave-one-KO-out: subset picked on 7 KOs, evaluated on held-out KO (separation, mean over spaces) ===\n"); print(format(loko, digits = 2), row.names = FALSE)
cat(sprintf("\nmean held-out separation: picked %.3f | R3 only %.3f | all 4 %.3f | average over all subsets %.3f\n", mean(loko$heldout_sep_picked), mean(loko$heldout_sep_R3only), mean(loko$heldout_sep_all4), mean(loko$heldout_sep_mean_all_subsets)))
cat(sprintf("in-sample best subset: %s (mean separation %.3f); range across subsets %.3f to %.3f\n", top$subset[1], top$mean_sep[1], min(top$mean_sep), max(top$mean_sep)))

save_fig("fig3N_H2052_NTC_subset_search", 12.5, 6.4, function() {
  layout(matrix(1:2, nrow = 1), widths = c(1.15, 1)); par(oma = c(0, 0, 3.2, 0), family = "sans")
  tt <- top[order(top$mean_sep), ]; par(mar = c(4.2, 12, 2.4, 1))
  cols <- ifelse(tt$n_ntc == 1 & tt$subset == sub("_NTC_", "", "H2052_NTC_R3") | tt$subset == "H2052R3", "#eda100", ifelse(tt$n_ntc == 4, "#4a3aa7", "grey60"))
  bp <- barplot(tt$mean_sep, horiz = TRUE, names.arg = rep("", nrow(tt)), col = cols, border = NA, xlab = "mean separation (own - null) across 4 feature spaces", xlim = range(c(0, tt$mean_sep)) + c(-0.02, 0.02))
  text(par("usr")[1], bp, tt$subset, adj = 1, xpd = NA, cex = 0.72); abline(v = 0, col = INK)
  mtext("all 15 NTC subsets (orange = H2052_NTC_R3 only; purple = all 4)", side = 3, adj = 0, line = 0.3, font = 2, cex = 0.8)
  par(mar = c(4.2, 4.2, 2.4, 1)); v <- rbind(picked = loko$heldout_sep_picked, R3_only = loko$heldout_sep_R3only, all4 = loko$heldout_sep_all4); colnames(v) <- loko$held_out
  bb <- barplot(v, beside = TRUE, col = c("#1baf7a", "#eda100", "#4a3aa7"), las = 2, ylab = "held-out KO separation (own - null)", border = NA, cex.names = 0.75); abline(h = 0, col = INK)
  legend("topright", legend = c("subset picked on other 7 KOs", "R3 only", "all 4"), fill = c("#1baf7a", "#eda100", "#4a3aa7"), border = NA, bty = "n", cex = 0.65)
  mtext("leave-one-KO-out check", side = 3, adj = 0, line = 0.3, font = 2, cex = 0.8)
  mtext("Which H2052 NTC set best reproduces the CROP-seq profile? (Spearman, NTC-KO, screen vs bulk H2052)", outer = TRUE, side = 3, adj = 0, line = 1.5, font = 2, cex = 0.95)
})
