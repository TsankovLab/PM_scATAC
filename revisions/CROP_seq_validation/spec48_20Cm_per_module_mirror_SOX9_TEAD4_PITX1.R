# Spec 48: for SOX9, TEAD4 and PITX1, do the individual Cm modules (e.g. Cm20 for TEAD4) go the same direction in the CROP-seq screen and in
# bulk? Per-module KO-NTC diff (median of AddModuleScore, top-50 genes per Cm from scRNA_PM_modules.rds -- same method as spec45), shown
# for the screen and for each bulk line (NCI-H2052 shown too, as the same-line reference; its NTC has n=1) plus bulk pooled across lines.
# needs: git_repo/files/scRNA_PM_modules.rds, spec_common.R (load_screen()), repro_common.R (load_bulk())
# outputs: spec48_per_module_diff_screen_vs_bulk.csv, fig3G_per_module_mirror_SOX9_TEAD4_PITX1.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat))
KOS <- c("SOX9", "TEAD4", "PITX1")

cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds")
top50 <- lapply(cm_genes, head, 50); CM <- names(top50)
options(Seurat.object.assay.version = "v3")

srt <- load_screen(); srt <- AddModuleScore(srt, features = top50, name = "Cm_", seed = 1)
sc <- as.matrix(srt@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(sc) <- CM; grp <- as.character(srt$merged_call)
med <- function(M, idx) apply(M[idx, , drop = FALSE], 2, median)
screen_diff <- t(sapply(KOS, function(k) med(sc, grp == k) - med(sc, grp == "NTC")))

b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1)
sb <- CreateSeuratObject(counts = as.matrix(b$counts)); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg))
sb <- AddModuleScore(sb, features = top50, name = "Cm_", seed = 1)
bs <- as.matrix(sb@meta.data[, paste0("Cm_", seq_along(CM))]); colnames(bs) <- CM; bs <- bs[rownames(meta), ]
bulk_diff <- function(cells) t(sapply(KOS, function(k) {
  ko <- intersect(rownames(meta)[meta$Gene_Control == k], cells); ntc <- intersect(rownames(meta)[meta$Gene_Control == "NTC"], cells)
  med(bs, ko) - med(bs, ntc) }))
rows <- list(Screen = screen_diff)
for (cl in LINES) rows[[cl]] <- bulk_diff(rownames(meta)[meta$Cell_Line == cl])
rows[["Bulk pooled"]] <- bulk_diff(rownames(meta))

long <- do.call(rbind, lapply(names(rows), function(s) do.call(rbind, lapply(KOS, function(k)
  data.frame(source = s, KO = k, Cm = CM, diff_KO_minus_NTC = rows[[s]][k, ])))))
write.csv(long, file.path(OUT, "spec48_per_module_diff_screen_vs_bulk.csv"), row.names = FALSE)

options(width = 220)
for (k in KOS) {
  cat(sprintf("\n=== %s: KO - NTC per Cm module (rows = source) ===\n", k))
  M <- do.call(rbind, lapply(names(rows), function(s) rows[[s]][k, ])); rownames(M) <- names(rows)
  print(round(M, 3))
  cat(sprintf("sign agreement screen vs bulk pooled: %d/20 modules; Pearson r = %.2f\n", sum(sign(M["Screen", ]) == sign(M["Bulk pooled", ])), cor(M["Screen", ], M["Bulk pooled", ])))
  cat(sprintf("Cm20: screen %+.3f | bulk pooled %+.3f | per line: %s\n", M["Screen", "Cm20"], M["Bulk pooled", "Cm20"],
              paste(sprintf("%s %+.3f", LINES, M[LINES, "Cm20"]), collapse = ", ")))
}

save_fig("fig3G_per_module_mirror_SOX9_TEAD4_PITX1", 13, 9, function() {
  layout(matrix(1:3, ncol = 1)); par(oma = c(0, 0, 3, 0), family = "sans")
  src <- names(rows); vmax <- max(abs(unlist(rows)))
  for (k in KOS) {
    M <- do.call(rbind, lapply(src, function(s) rows[[s]][k, ])); rownames(M) <- src
    ord <- order(M["Screen", ], decreasing = TRUE); M <- M[, ord]
    par(mar = c(3.5, 8, 2.2, 1))
    image(seq_len(ncol(M)), seq_len(nrow(M)), t(M)[, nrow(M):1], col = DIV, zlim = c(-vmax, vmax), axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(nrow(M))) for (j in seq_len(ncol(M))) text(j, nrow(M) - i + 1, sprintf("%.2f", M[i, j]), cex = 0.55, col = INK)
    text(seq_len(ncol(M)), 0.35, colnames(M), srt = 90, adj = 1, xpd = NA, cex = 0.75, col = ifelse(colnames(M) == "Cm20", "#d62728", INK2), font = ifelse(colnames(M) == "Cm20", 2, 1))
    axis(2, at = seq_len(nrow(M)), labels = rev(rownames(M)), las = 2, tick = FALSE, cex.axis = 0.75, col.axis = INK2)
    abline(h = nrow(M) - 0.5, col = "white", lwd = 3); box(col = MUTED, lwd = 0.6)
    mtext(k, side = 3, adj = 0, line = 0.3, font = 2, cex = 0.95, col = KO_COL[k])
  }
  mtext("Per-module KO-NTC diff (median AddModuleScore): CROP-seq screen vs each bulk line, modules sorted by screen effect (Cm20 in red)", outer = TRUE, side = 3, adj = 0, line = 1.2, font = 2, cex = 1.0)
})
