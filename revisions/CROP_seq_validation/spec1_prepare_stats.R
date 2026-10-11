# Spec 1: per-KO statistics vs NTC used by the specificity analysis (spec2, spec3).
#   screen: Welch t, mean difference vs NTC and absolute group mean of log-normalised expression per KO (all genes detected in >=1% cells and present in bulk),
#           without and with the cell-cycle index (cc = S.Score + G2M.Score) regressed out pooled over all cells (as Seurat ScaleData does)
#   bulk (RUN1, each cell line): DESeq2 Wald stat / log2FC per KO vs NTC, design ~condition and ~cc+condition (cc = bulk cell-cycle index)
# needs: screen_export_figs.R, bulk_cellcycle_RUN1.R
# outputs: spec_screen_stats.rds, spec_bulk_stats.rds, spec_ontarget_knockdown.csv
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R"))
suppressPackageStartupMessages(library(DESeq2))

## ------------------------------------------------ screen
b <- load_bulk(); meta <- b$meta
srt <- load_screen(); Y <- GetAssayData(srt, layer = "data")
keep <- rownames(Y) %in% rownames(b$cp10k) & Matrix::rowSums(Y > 0) >= 0.01 * ncol(Y); Y <- Y[keep, ]
message("screen genes kept (in bulk, detected in >=1% cells): ", nrow(Y), " | cells: ", ncol(Y))
S0 <- screen_ko_stats(Y, srt$merged_call, srt$cc, FALSE); S1 <- screen_ko_stats(Y, srt$merged_call, srt$cc, TRUE)
saveRDS(list(genes = rownames(Y), diff0 = S0$diff, t0 = S0$t, abs0 = S0$abs, diffcc = S1$diff, tcc = S1$t, abscc = S1$abs, n = S0$n), file.path(OUT, "spec_screen_stats.rds"))
cat("cor(t without vs with cc regression) per KO:\n"); print(round(sapply(colnames(S0$t), function(k) cor(S0$t[, k], S1$t[, k])), 3))

## ------------------------------------------------ bulk DESeq2 per cell line
ccb <- read.csv(file.path(OUT, "bulk_cellcycle_per_sample_RUN1.csv"))
out <- list(); ontarget <- list()
for (cl in LINES) {
  m <- meta[meta$Cell_Line == cl, ]; s <- rownames(m); x <- round(b$counts[, s]); storage.mode(x) <- "integer"
  x <- x[rowSums(x >= 10) >= 2, ]
  cond <- factor(m$Gene_Control, levels = c("NTC", sort(setdiff(unique(m$Gene_Control), "NTC"), method = "radix")))
  cc <- as.numeric(scale(ccb$cc[match(s, ccb$sample)]))
  fit <- function(des, cd) DESeq(DESeqDataSetFromMatrix(x, cd, des), quiet = TRUE)
  d0 <- fit(~ condition, DataFrame(condition = cond)); d1 <- fit(~ cc + condition, DataFrame(condition = cond, cc = cc))
  kos <- levels(cond)[-1]
  get <- function(d, what) sapply(setNames(kos, kos), function(k) results(d, contrast = c("condition", k, "NTC"))[[what]])
  o <- list(genes = rownames(x), stat0 = get(d0, "stat"), lfc0 = get(d0, "log2FoldChange"), padj0 = get(d0, "padj"),
            statcc = get(d1, "stat"), lfccc = get(d1, "log2FoldChange"), baseMean = rowMeans(counts(d0, normalized = TRUE)))
  for (a in c("stat0", "lfc0", "padj0", "statcc", "lfccc")) rownames(o[[a]]) <- o$genes
  out[[cl]] <- o
  ontarget[[cl]] <- do.call(rbind, lapply(kos, function(k) if (k %in% o$genes)
    data.frame(cell_line = cl, KO = k, n_KO = sum(m$Gene_Control == k), n_NTC = sum(m$Gene_Control == "NTC"), lfc = o$lfc0[k, k], stat = o$stat0[k, k], padj = o$padj0[k, k], lfc_cc = o$lfccc[k, k]) else NULL))
  message(cl, ": ", length(o$genes), " genes, ", length(kos), " KO groups")
}
saveRDS(out, file.path(OUT, "spec_bulk_stats.rds"))
ot <- do.call(rbind, ontarget); rownames(ot) <- NULL; write.csv(ot, file.path(OUT, "spec_ontarget_knockdown.csv"), row.names = FALSE)
print(format(ot, digits = 3), row.names = FALSE)
