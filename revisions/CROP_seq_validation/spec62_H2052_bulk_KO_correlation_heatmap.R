# Spec 62: sample-level correlation heatmap of the bulk NCI-H2052 samples only (no screen involved) -- Spearman correlation of log2 CP10K expression across
# all pairs of H2052 bulk samples (KOs + NTC), in the pipeline's usual feature-space sweep (top 1000/2000/3000 variable genes, ranked by bulk rowVar averaged across
# the 4 lines -- spec22/29/30/47 convention) plus the 20 Cm modules (AddModuleScore, top-50 genes/module). Samples ordered/coloured by KO group. Env NTC4=1 adds the
# 3 QC-dropped H2052 NTC samples (H28_NTC_R1/2/3); default is the curated sheet's samples only (1 NTC).
# Both RAW (plain log2 CP10K / module score) and GENE-CENTRED (each feature's mean across these H2052 samples subtracted first, as fig13c/spec11 does -- removes
# the shared-baseline structure that otherwise pushes every pairwise correlation toward 1 and swamps the KO-specific signal) versions are produced.
# needs: git_repo/files/scRNA_PM_modules.rds, repro_common.R
# outputs: spec62_H2052_sample_correlation{_raw,_centred}.csv, fig3X_H2052_bulk_sample_correlation_heatmap_{raw,centred}_{Cm,top1000,top2000,top3000}{_NTC4}.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat)); options(Seurat.object.assay.version = "v3"); TOPNS <- c(1000, 2000, 3000)
NTC4 <- Sys.getenv("NTC4", "0") == "1"; TAG <- if (NTC4) "_NTC4" else ""
b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1); cnt <- b$counts
if (NTC4) { a <- load_full_audited_bulk(); add <- rownames(a$meta)[a$meta$Cell_Line == "NCI-H2052" & a$meta$Gene_Control == "NTC" & !rownames(a$meta) %in% rownames(meta)]
  g <- intersect(rownames(lg), rownames(a$lg)); meta <- rbind(meta[, c("Sample_ID", "Cell_Line", "Gene_Control")], a$meta[add, c("Sample_ID", "Cell_Line", "Gene_Control")])
  lg <- cbind(lg[g, ], a$lg[g, add])[, rownames(meta)]; cnt <- cbind(cnt[g, ], a$cs[g, add])[, rownames(meta)] }
h <- rownames(meta)[meta$Cell_Line == "NCI-H2052" & meta$Gene_Control %in% c("NTC", ORDER)]; hm <- meta[h, ]; hm <- hm[order(factor(hm$Gene_Control, levels = c("NTC", ORDER))), ]; h <- rownames(hm)
cat(sprintf("H2052 bulk samples: %d (%s)\n", length(h), paste(sprintf("%s=%d", names(table(hm$Gene_Control)), table(hm$Gene_Control)), collapse = ", ")))

cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds"); top50 <- lapply(cm_genes, head, 50); CM <- names(top50)
sb <- CreateSeuratObject(counts = as.matrix(cnt[, h])); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg[, h])); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
bcm <- as.matrix(sb@meta.data[h, paste0("Cm_", seq_along(CM))]); colnames(bcm) <- CM
lv <- sapply(LINES, function(cl) matrixStats::rowVars(as.matrix(lg[, rownames(meta)[meta$Cell_Line == cl], drop = FALSE]))); vord <- rownames(lg)[order(rowMeans(lv), decreasing = TRUE)]
spaces <- c("20 Cm modules", sprintf("top %d variable genes", TOPNS)); sp_tag <- c("Cm", paste0("top", TOPNS)); names(sp_tag) <- spaces
feat_raw <- function(sp) if (sp == "20 Cm modules") t(bcm) else { n <- as.integer(sub("top (\\d+) .*", "\\1", sp)); lg[vord[1:n], h] }
centre <- function(M) M - rowMeans(M)   # subtract each feature's mean ACROSS THESE H2052 SAMPLES ONLY (fig13c/spec11 convention) -- removes the shared baseline

grp <- hm[h, "Gene_Control"]; options(width = 200)
KIND <- c(raw = "RAW (log2 CP10K / module score)", centred = "GENE-CENTRED (feature mean across these H2052 samples subtracted)")
for (kind in names(KIND)) {
  cors <- lapply(spaces, function(sp) { M <- feat_raw(sp); if (kind == "centred") M <- centre(M); cor(M, method = "spearman") }); names(cors) <- spaces
  long <- do.call(rbind, lapply(spaces, function(sp) { C <- cors[[sp]]; C[upper.tri(C, diag = TRUE)] <- NA
    d <- as.data.frame(as.table(C)); d <- d[!is.na(d$Freq), ]; colnames(d) <- c("sample_a", "sample_b", "rho"); d$space <- sp; d }))
  write.csv(long, file.path(OUT, paste0("spec62_H2052_sample_correlation_", kind, TAG, ".csv")), row.names = FALSE)
  cat(sprintf("\n=== %s ===\n", KIND[kind]))
  cat("mean within-KO-group vs cross-group Spearman correlation, per feature space (NTC excluded from 'within' since n=1 possible):\n")
  for (sp in spaces) { C <- cors[[sp]]; within <- c(); cross <- c()
    for (i in seq_along(h)) for (j in seq_along(h)) if (i < j) { if (grp[i] == grp[j] && grp[i] != "NTC") within <- c(within, C[i, j]) else cross <- c(cross, C[i, j]) }
    cat(sprintf("  %-24s within-KO mean rho = %.3f (n=%d) | cross-KO mean rho = %.3f (n=%d)\n", sp, mean(within), length(within), mean(cross), length(cross))) }

  for (sp in spaces) save_fig(paste0("fig3X_H2052_bulk_sample_correlation_heatmap_", kind, "_", sp_tag[[sp]], TAG), 8.5, 7.5, function() {
    C <- cors[[sp]]; n <- nrow(C); par(mar = c(6, 6.5, 3.6, 1), family = "sans")
    zlim <- range(C[row(C) != col(C)])   # scale to this panel's actual off-diagonal range, not the full -1..1 scale
    Cp <- C; diag(Cp) <- NA   # diagonal (self-correlation = 1.00, trivial) deliberately blanked instead of left to fall outside zlim by accident
    image(seq_len(n), seq_len(n), t(Cp)[, n:1], col = colorRampPalette(c("#2a78d6", "#f0efec", "#eb6834"))(100), zlim = zlim, axes = FALSE, xlab = "", ylab = "")
    for (i in seq_len(n)) for (j in seq_len(n)) text(j, n - i + 1, if (i == j) "" else sprintf("%.2f", C[i, j]), cex = 0.55, col = INK)
    gb <- which(diff(as.integer(factor(grp, levels = c("NTC", ORDER)))) != 0); for (g in gb) { abline(v = g + 0.5, col = "white", lwd = 2.2); abline(h = n - g + 0.5, col = "white", lwd = 2.2) }
    axis(1, at = seq_len(n), labels = FALSE); text(seq_len(n), 0.15, h, srt = 60, adj = 1, xpd = NA, cex = 0.6, col = KO_COL[grp])
    axis(2, at = n:1, labels = FALSE); text(0.2, n:1, h, adj = 1, xpd = NA, cex = 0.6, col = KO_COL[grp]); box(col = MUTED, lwd = 0.6)
    mtext(sprintf("NCI-H2052 bulk, sample-sample Spearman correlation (%s): %s", kind, sp), side = 3, adj = 0, line = 1.9, font = 2, cex = 0.9)
    mtext(sprintf("%s; colour scaled to this panel's off-diagonal range [%.2f, %.2f]; diagonal (self-corr.=1.00) blanked; samples grouped by KO (colour = KO)%s",
                  if (kind == "raw") "log2 CP10K (Cm: AddModuleScore top-50 genes/module)" else "each feature mean-centred across these H2052 samples first (fig13c/spec11 convention)",
                  zlim[1], zlim[2], if (NTC4) "; H2052 NTC = 4 samples (+H28_NTC_R1/2/3)" else ""), side = 3, adj = 0, line = 0.6, cex = 0.58, col = INK2)
  })
}
