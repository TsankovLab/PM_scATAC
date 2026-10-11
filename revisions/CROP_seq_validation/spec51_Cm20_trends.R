# Spec 51: Cm20 trends across all shared KOs -- Cm20 module score (AddModuleScore, top-50 genes, same method as spec45/48) per KO group, in the CROP-seq
# screen and in each bulk line, shown as (A) NTC - KO median difference per KO (dots = sources), (B) raw score distributions per group and source.
# Direction: NTC - KO by default (positive = Cm20 higher in NTC, i.e. lower after KO); DIFF_DIR=KO_NTC for the reverse (files get _KOminusNTC).
# needs: git_repo/files/scRNA_PM_modules.rds, spec_common.R (load_screen()), repro_common.R (load_bulk())
# outputs: spec51_Cm20_diff_per_KO{DSUF}.csv, spec51_Cm20_scores_per_sample{DSUF}.csv, fig3H_Cm20_trends_diff{DSUF}, fig3I_Cm20_scores_by_group{DSUF}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat))
DIFF_DIR <- Sys.getenv("DIFF_DIR", "NTC_KO"); SGN <- if (DIFF_DIR == "NTC_KO") -1 else 1; DLAB <- if (DIFF_DIR == "NTC_KO") "NTC-KO" else "KO-NTC"; DSUF <- if (DIFF_DIR == "NTC_KO") "" else "_KOminusNTC"
MOD <- "Cm20"
cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds")
top50 <- lapply(cm_genes, head, 50); CM <- names(top50); mi <- match(MOD, CM); options(Seurat.object.assay.version = "v3")

srt <- load_screen(); srt <- AddModuleScore(srt, features = top50, name = "Cm_", seed = 1)
scr <- data.frame(source = "CROP-seq (H2052)", group = as.character(srt$merged_call), score = srt@meta.data[[paste0("Cm_", mi)]], stringsAsFactors = FALSE)
b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1)
sb <- CreateSeuratObject(counts = as.matrix(b$counts)); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg)); sb <- AddModuleScore(sb, features = top50, name = "Cm_", seed = 1)
bsc <- sb@meta.data[rownames(meta), paste0("Cm_", mi)]
blk <- data.frame(source = meta$Cell_Line, group = meta$Gene_Control, score = bsc, sample = rownames(meta), stringsAsFactors = FALSE)
scr$sample <- NA; all <- rbind(scr, blk); KOS <- intersect(ORDER, unique(scr$group)); all <- all[all$group %in% c("NTC", KOS), ]
write.csv(all, file.path(OUT, paste0("spec51_Cm20_scores_per_sample", DSUF, ".csv")), row.names = FALSE)

SRC <- c("CROP-seq (H2052)", LINES); n_ntc <- sapply(SRC, function(s) length(unique(all$sample[all$source == s & all$group == "NTC"])))
med <- function(s, g) { x <- all$score[all$source == s & all$group == g]; if (length(x)) median(x) else NA }
D <- do.call(rbind, lapply(SRC, function(s) data.frame(source = s, KO = KOS, diff = SGN * sapply(KOS, function(k) med(s, k) - med(s, "NTC")),
      n_KO = sapply(KOS, function(k) sum(all$source == s & all$group == k)), n_NTC = sum(all$source == s & all$group == "NTC"))))
write.csv(D, file.path(OUT, paste0("spec51_Cm20_diff_per_KO", DSUF, ".csv")), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== %s %s (median AddModuleScore diff; positive = higher in %s) ===\n", MOD, DLAB, if (SGN < 0) "NTC" else "KO"))
W <- reshape(D[, c("source", "KO", "diff")], idvar = "KO", timevar = "source", direction = "wide"); colnames(W) <- sub("diff\\.", "", colnames(W)); print(format(W, digits = 3), row.names = FALSE)
cat("\nsign agreement with screen per bulk line (n KOs):\n"); for (l in LINES) cat(sprintf("  %-10s %d/%d\n", l, sum(sign(W[[l]]) == sign(W[["CROP-seq (H2052)"]]), na.rm = TRUE), sum(!is.na(W[[l]]))))
cat("\nKO order by screen diff (high to low): ", paste(W$KO[order(-W[["CROP-seq (H2052)"]])], collapse = " > "), "\n")
cat("Spearman screen vs each bulk line across KOs:\n"); for (l in LINES) cat(sprintf("  %-10s rho=%.2f\n", l, suppressWarnings(cor(W[["CROP-seq (H2052)"]], W[[l]], method = "spearman", use = "complete.obs"))))

## fig3H: per KO, screen diamond + bulk-line dots (H2052 bulk hollow -- 1 NTC, reference only), KOs sorted by screen diff
save_fig(paste0("fig3H_Cm20_trends_diff", DSUF), 9.5, 5.6, function() {
  ko_ord <- W$KO[order(-W[["CROP-seq (H2052)"]])]; par(mar = c(5, 4.8, 3.6, 1), family = "sans")
  yl <- range(D$diff, na.rm = TRUE); yl <- yl + c(-0.08, 0.08) * diff(yl)
  plot(NA, xlim = c(0.5, length(ko_ord) + 0.5), ylim = yl, xaxt = "n", xlab = "", ylab = sprintf("%s Cm20 score diff (median)", DLAB))
  axis(1, at = seq_along(ko_ord), labels = FALSE); text(seq_along(ko_ord), par("usr")[3] - 0.04 * diff(par("usr")[3:4]), ko_ord, srt = 45, adj = 1, xpd = NA)
  abline(h = 0, lty = 2, col = INK); abline(v = seq_along(ko_ord)[-1] - 0.5, col = "grey92")
  off <- c("MSTO-211H" = -0.24, "NCI-H2452" = -0.08, "NCI-H28" = 0.08, "NCI-H2052" = 0.24)
  for (i in seq_along(ko_ord)) { for (l in LINES) { v <- D$diff[D$source == l & D$KO == ko_ord[i]]
      if (l == "NCI-H2052") points(i + off[[l]], v, pch = 1, cex = 1.3, col = LINE_COL[[l]], lwd = 1.8) else points(i + off[[l]], v, pch = 16, cex = 1.3, col = LINE_COL[[l]]) }
    points(i, D$diff[D$source == "CROP-seq (H2052)" & D$KO == ko_ord[i]], pch = 23, bg = LINE_COL[["CROP-seq"]], col = INK, cex = 2, lwd = 1.2) }
  legend("topright", legend = c("CROP-seq screen", "MSTO-211H", "NCI-H2452", "NCI-H28", "NCI-H2052 bulk (1 NTC)"), pch = c(23, 16, 16, 16, 1), pt.bg = LINE_COL[["CROP-seq"]],
         col = c(INK, LINE_COL[c("MSTO-211H", "NCI-H2452", "NCI-H28", "NCI-H2052")]), bty = "n", cex = 0.65, ncol = 2)
  mtext(sprintf("%s trends across KOs: CROP-seq screen vs bulk lines (%s)", MOD, DLAB), side = 3, adj = 0, line = 2.2, font = 2, cex = 1.0)
  mtext(sprintf("median AddModuleScore diff per KO; positive = Cm20 higher in %s; KOs sorted by screen effect", if (SGN < 0) "NTC (lower after KO)" else "KO"), side = 3, adj = 0, line = 0.9, cex = 0.58, col = INK2)
})

## fig3I: raw Cm20 score per group, one panel per source (NTC first, then KOs sorted by screen diff)
save_fig(paste0("fig3I_Cm20_scores_by_group", DSUF), 13, 7.4, function() {
  grp_ord <- c("NTC", W$KO[order(-W[["CROP-seq (H2052)"]])]); layout(matrix(1:6, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3.2, 0), family = "sans")
  for (s in SRC) { d <- all[all$source == s, ]; par(mar = c(5, 4, 2.2, 0.8)); d$group <- factor(d$group, levels = grp_ord)
    boxplot(score ~ group, d, las = 2, xlab = "", ylab = "Cm20 score", outline = FALSE, col = c("grey80", adjustcolor(KO_COL[grp_ord[-1]], 0.55)), border = INK2, cex.axis = 0.7,
            ylim = if (s == "CROP-seq (H2052)") NULL else range(d$score))
    if (s != "CROP-seq (H2052)") { set.seed(1); stripchart(score ~ group, d, vertical = TRUE, add = TRUE, method = "jitter", jitter = 0.12, pch = 16, cex = 0.6, col = INK) }
    abline(h = median(d$score[d$group == "NTC"]), lty = 2, col = "#d62728")
    mtext(sprintf("%s (n NTC = %s)", s, if (s == "CROP-seq (H2052)") format(sum(d$group == "NTC"), big.mark = ",") else n_ntc[[s]]), side = 3, adj = 0, line = 0.3, font = 2, cex = 0.85) }
  mtext(sprintf("%s module score per group, screen and bulk (dashed red = NTC median); groups ordered NTC then KOs by screen %s effect", MOD, DLAB), outer = TRUE, side = 3, adj = 0, line = 1.2, font = 2, cex = 0.95)
})
