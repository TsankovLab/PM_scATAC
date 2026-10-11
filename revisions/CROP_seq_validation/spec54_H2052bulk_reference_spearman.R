# Spec 54: Spearman screen-free cross-line check: bulk NCI-H2052 NTC-KO profile (per KO) used as the REFERENCE and correlated with the same KO in each OTHER bulk line
# (own-KO, coloured, dots = lines) vs the reference KO against every OTHER KO in that line (null, grey). Feature spaces: 20 Cm modules (AddModuleScore, top-50 genes,
# same as spec45) and the top 1000/2000/3000 most variable genes (bulk rowVar averaged across lines, spec22/29/30/47 convention). Direction NTC - KO (Spearman is sign-flip
# invariant when both sides are flipped). NCI-H2052 has only 1 NTC sample, so the reference is noisy (single-NTC median).
# env: H2052_NTC4=1 (add the 3 QC-dropped H2052 NTC samples, H28_NTC_R1/2/3 -> 4 NTC), EXCL_LINES (comma list of OTHER lines to drop, default none), DIFF_DIR (NTC_KO default)
# needs: git_repo/files/scRNA_PM_modules.rds, repro_common.R (load_bulk())
# outputs: spec54_H2052ref_spearman_{own,cross}{ESUF}.csv, spec54_H2052ref_summary{ESUF}.csv, fig3L_H2052bulk_ref_spearman_Cm_and_vg{ESUF}.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat))
DIFF_DIR <- Sys.getenv("DIFF_DIR", "NTC_KO"); SGN <- if (DIFF_DIR == "NTC_KO") -1 else 1; DLAB <- if (DIFF_DIR == "NTC_KO") "NTC-KO" else "KO-NTC"
EXCL <- setdiff(strsplit(Sys.getenv("EXCL_LINES", ""), ",")[[1]], ""); ESUF <- if (length(EXCL)) paste0("_no", paste(sub("NCI-|-211H", "", EXCL), collapse = "_")) else ""; if (Sys.getenv("H2052_NTC4", "0") == "1") ESUF <- paste0(ESUF, "_H2052NTC4")
REF <- "NCI-H2052"; OTHER <- setdiff(LINES, c(REF, EXCL)); TOPNS <- c(1000, 2000, 3000)
b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1); bcounts <- b$counts
# env H2052_NTC4=1: also use the 3 H2052 NTC samples dropped from the curated sheet (H28_NTC_R1/2/3: sequenced under an H28 label, genotyped as NCI-H2052 and
# registered as such in load_full_audited_bulk()/H28_NTC_NATIVE) -> H2052 gets 4 NTC (with H2052_NTC_R3) instead of 1. Note these 3 are PC2>20 outliers in the H2052 PCA.
ADD4 <- Sys.getenv("H2052_NTC4", "0") == "1"
if (ADD4) { a <- load_full_audited_bulk(); add <- rownames(a$meta)[a$meta$Cell_Line == "NCI-H2052" & a$meta$Gene_Control == "NTC" & !rownames(a$meta) %in% rownames(meta)]
  cat("adding H2052 NTC samples:", paste(add, collapse = ", "), "\n"); g <- intersect(rownames(lg), rownames(a$cs))
  meta <- rbind(meta[, c("Sample_ID", "Cell_Line", "Gene_Control")], a$meta[add, c("Sample_ID", "Cell_Line", "Gene_Control")])
  bcounts <- cbind(bcounts[g, ], a$cs[g, add]); lg <- log2(cbind(b$cp10k[g, ], a$cp10k[g, add]) + 1)[, rownames(meta)] }
KOS <- intersect(ORDER, unique(meta$Gene_Control[meta$Cell_Line == REF]))
cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds"); top50 <- lapply(cm_genes, head, 50); CM <- names(top50)
options(Seurat.object.assay.version = "v3"); sb <- CreateSeuratObject(counts = as.matrix(bcounts[rownames(lg), colnames(lg)])); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
bs <- as.matrix(sb@meta.data[rownames(meta), paste0("Cm_", seq_along(CM))]); colnames(bs) <- CM
lv <- sapply(LINES, function(cl) matrixStats::rowVars(as.matrix(lg[, rownames(meta)[meta$Cell_Line == cl], drop = FALSE]))); vord <- rownames(lg)[order(rowMeans(lv), decreasing = TRUE)]

## per-line NTC-KO profile matrices (KO x feature), features = matrix rows (genes) or Cm modules
diff_line <- function(cl, FM) { mm <- meta[meta$Cell_Line == cl, ]; ntc <- rownames(mm)[mm$Gene_Control == "NTC"]; if (!length(ntc)) return(NULL)
  t(sapply(KOS, function(k) { ko <- rownames(mm)[mm$Gene_Control == k]; if (!length(ko)) return(rep(NA_real_, nrow(FM)))
    SGN * (apply(FM[, ko, drop = FALSE], 1, median) - apply(FM[, ntc, drop = FALSE], 1, median)) })) }
spaces <- c(list("20 Cm modules" = t(bs)), setNames(lapply(TOPNS, function(n) lg[vord[1:n], , drop = FALSE]), sprintf("top %d variable genes", TOPNS)))
own <- list(); cross <- list()
for (sp in names(spaces)) { FMAT <- spaces[[sp]]; ref <- diff_line(REF, FMAT); dl <- lapply(setNames(OTHER, OTHER), function(cl) diff_line(cl, FMAT))
  for (cl in OTHER) { M <- dl[[cl]]; for (k in KOS) { if (all(is.na(M[k, ]))) next
      own[[length(own) + 1]] <- data.frame(space = sp, line = cl, KO = k, rho = suppressWarnings(cor(ref[k, ], M[k, ], method = "spearman", use = "complete.obs")))
      for (k2 in setdiff(KOS, k)) if (!all(is.na(M[k2, ]))) cross[[length(cross) + 1]] <- data.frame(space = sp, line = cl, ref_KO = k, line_KO = k2, rho = suppressWarnings(cor(ref[k, ], M[k2, ], method = "spearman", use = "complete.obs"))) } } }
own <- do.call(rbind, own); cross <- do.call(rbind, cross)
write.csv(own, file.path(OUT, paste0("spec54_H2052ref_spearman_own", ESUF, ".csv")), row.names = FALSE); write.csv(cross, file.path(OUT, paste0("spec54_H2052ref_spearman_cross", ESUF, ".csv")), row.names = FALSE)
summ <- do.call(rbind, lapply(names(spaces), function(sp) do.call(rbind, lapply(KOS, function(k) { o <- own$rho[own$space == sp & own$KO == k]; n <- cross$rho[cross$space == sp & cross$ref_KO == k]
  data.frame(space = sp, KO = k, own_median = median(o), null_median = median(n), separation = median(o) - median(n), n_lines = length(o)) }))))
write.csv(summ, file.path(OUT, paste0("spec54_H2052ref_summary", ESUF, ".csv")), row.names = FALSE); options(width = 200)
cat(sprintf("\nreference = bulk %s (%s); other lines: %s\n", REF, DLAB, paste(OTHER, collapse = ", ")))
for (sp in names(spaces)) { s <- summ[summ$space == sp, ]; cat(sprintf("\n== %s: pooled own-KO median rho %.3f | null %.3f | separation %+.3f ==\n", sp, median(own$rho[own$space == sp]), median(cross$rho[cross$space == sp]), median(own$rho[own$space == sp]) - median(cross$rho[cross$space == sp])))
  print(format(s[order(-s$own_median), c("KO", "own_median", "null_median", "separation")], digits = 2), row.names = FALSE) }

save_fig(paste0("fig3L_H2052bulk_ref_spearman_Cm_and_vg", ESUF), 13, 9.4, function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.4, 0), family = "sans")
  for (sp in names(spaces)) { o <- own[own$space == sp, ]; cr <- cross[cross$space == sp, ]; ko_ord <- names(sort(sapply(split(o$rho, o$KO), median), decreasing = TRUE)); n_ko <- length(ko_ord)
    par(mar = c(5, 4.4, 2.6, 0.8)); plot(NA, xlim = c(0.5, n_ko + 0.5), ylim = c(-1, 1), xaxt = "n", xlab = "", ylab = sprintf("Spearman rho (H2052 bulk vs line, %s)", DLAB))
    axis(1, at = seq_len(n_ko), labels = FALSE); text(seq_len(n_ko), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), ko_ord, srt = 45, adj = 1, xpd = NA, cex = 0.85)
    for (i in seq_len(n_ko)) boxplot(cr$rho[cr$ref_KO == ko_ord[i]], at = i - 0.18, add = TRUE, boxwex = 0.3, col = "grey85", border = "grey40", outline = FALSE, axes = FALSE)
    for (i in seq_len(n_ko)) boxplot(o$rho[o$KO == ko_ord[i]], at = i + 0.18, add = TRUE, boxwex = 0.3, col = adjustcolor(KO_COL[ko_ord[i]], 0.5), border = INK2, outline = FALSE, axes = FALSE)
    set.seed(1); for (i in seq_len(n_ko)) { d <- o[o$KO == ko_ord[i], ]; points(jitter(rep(i + 0.18, nrow(d)), amount = 0.06), d$rho, pch = 16, cex = 1, col = LINE_COL[d$line]) }
    abline(h = 0, col = INK, lty = 2)
    if (sp == names(spaces)[1]) legend("bottomleft", legend = c(OTHER, "own KO", "other KOs (null)"), col = c(LINE_COL[OTHER], "grey30", "grey60"), pch = c(rep(16, length(OTHER)), 15, 15), bty = "n", cex = 0.6)
    mtext(sprintf("%s (own %.3f vs null %.3f)", sp, median(o$rho), median(cr$rho)), side = 3, adj = 0, line = 0.4, font = 2, cex = 0.85) }
  mtext(sprintf("Bulk NCI-H2052 as reference vs other bulk lines, Spearman (%s)%s", DLAB, if (length(EXCL)) paste0("; ", paste(EXCL, collapse = ", "), " removed") else ""), outer = TRUE, side = 3, adj = 0, line = 1.9, font = 2, cex = 1.0)
  mtext(paste0("coloured = H2052 KO profile vs the SAME KO in that line; grey = H2052 KO vs every OTHER KO in that line (null); H2052 NTC n=", if (ADD4) "4 (+H28_NTC_R1/2/3)" else "1"), outer = TRUE, side = 3, adj = 0, line = 0.6, cex = 0.6, col = INK2)
})
