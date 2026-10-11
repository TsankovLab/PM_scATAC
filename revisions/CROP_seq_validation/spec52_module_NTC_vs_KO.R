# Spec 52: is a given Cm module higher in NTC than in a given KO, in the screen and in each bulk line? (default Cm1, SOX9; env MOD, KO)
# Cm module score = AddModuleScore on the top-50 genes (same as spec45/48/51). Direction NTC - KO (positive = higher in NTC).
# Tests: screen = Wilcoxon on cells (NTC vs KO); per bulk line = NTC-KO difference of sample means/medians with a Welch t (n=2-3, descriptive);
# bulk pooled over lines with a line-adjusted linear model score ~ line + group (H2052 excluded: only 1 NTC sample) and a sign test across lines.
# needs: git_repo/files/scRNA_PM_modules.rds, spec_common.R (load_screen()), repro_common.R (load_bulk())
# outputs: spec52_{MOD}_NTC_vs_{KO}.csv, fig3J_{MOD}_NTC_vs_{KO}.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat))
MOD <- Sys.getenv("MOD", "Cm1"); KO <- Sys.getenv("KO", "SOX9")
cm_genes <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/files/scRNA_PM_modules.rds")
top50 <- lapply(cm_genes, head, 50); CM <- names(top50); mi <- match(MOD, CM); stopifnot(!is.na(mi)); options(Seurat.object.assay.version = "v3")
srt <- suppressWarnings(AddModuleScore(load_screen(), features = top50, name = "Cm_", seed = 1))
scr <- data.frame(source = "CROP-seq (H2052)", group = as.character(srt$merged_call), score = srt@meta.data[[paste0("Cm_", mi)]], sample = NA_character_)
b <- load_bulk(); meta <- b$meta; lg <- log2(b$cp10k + 1)
sb <- CreateSeuratObject(counts = as.matrix(b$counts)); sb <- SetAssayData(sb, slot = "data", new.data = as.matrix(lg)); sb <- suppressWarnings(AddModuleScore(sb, features = top50, name = "Cm_", seed = 1))
blk <- data.frame(source = meta$Cell_Line, group = meta$Gene_Control, score = sb@meta.data[rownames(meta), paste0("Cm_", mi)], sample = rownames(meta))
all <- rbind(scr, blk); all <- all[all$group %in% c("NTC", KO), ]; SRC <- c("CROP-seq (H2052)", LINES)
res <- do.call(rbind, lapply(SRC, function(s) { d <- all[all$source == s, ]; n <- d$score[d$group == "NTC"]; k <- d$score[d$group == KO]
  p <- if (s == "CROP-seq (H2052)") wilcox.test(n, k)$p.value else if (length(n) > 1 && length(k) > 1) tryCatch(t.test(n, k)$p.value, error = function(e) NA) else NA
  data.frame(source = s, n_NTC = length(n), n_KO = length(k), mean_NTC = mean(n), mean_KO = mean(k), median_NTC = median(n), median_KO = median(k),
             diff_mean = mean(n) - mean(k), diff_median = median(n) - median(k), NTC_higher = mean(n) > mean(k), p = p, test = if (s == "CROP-seq (H2052)") "Wilcoxon (cells)" else "Welch t (samples)") }))
write.csv(res, file.path(OUT, sprintf("spec52_%s_NTC_vs_%s.csv", MOD, KO)), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== %s: NTC vs %s (NTC - KO; positive = higher in NTC) ===\n", MOD, KO)); print(format(res[, c("source", "n_NTC", "n_KO", "mean_NTC", "mean_KO", "diff_mean", "diff_median", "NTC_higher", "p", "test")], digits = 3), row.names = FALSE)
bl <- all[all$source %in% setdiff(LINES, "NCI-H2052"), ]; if (nrow(bl)) { bl$group <- relevel(factor(bl$group), ref = KO); m <- summary(lm(score ~ source + group, bl))$coefficients["groupNTC", ]
  cat(sprintf("\nbulk (MSTO/H2452/H28, line-adjusted lm): NTC - %s = %+.3f (t=%.2f, p=%.3g)\n", KO, m[1], m[3], m[4])) }
oth <- setdiff(LINES, "NCI-H2052"); nh <- sum(res$NTC_higher[res$source %in% oth]); cat(sprintf("lines with NTC > %s (excluding H2052 bulk): %d/%d; incl. H2052 bulk (1 NTC): %d/4\n", KO, nh, length(oth), sum(res$NTC_higher[res$source %in% LINES])))

save_fig(sprintf("fig3J_%s_NTC_vs_%s", MOD, KO), 12, 4.8, function() {
  layout(matrix(1:5, nrow = 1)); par(oma = c(0, 0, 3.2, 0), family = "sans")
  for (s in SRC) { d <- all[all$source == s, ]; d$group <- factor(d$group, levels = c("NTC", KO)); par(mar = c(4, 4.2, 2.6, 0.5))
    boxplot(score ~ group, d, outline = FALSE, col = c("grey80", adjustcolor(KO_COL[[KO]], 0.55)), border = INK2, xlab = "", ylab = sprintf("%s score", MOD), names = c("NTC", KO), ylim = range(d$score))
    set.seed(1); stripchart(score ~ group, d, vertical = TRUE, add = TRUE, method = "jitter", jitter = 0.12, pch = 16, cex = if (s == "CROP-seq (H2052)") 0.15 else 1, col = if (s == "CROP-seq (H2052)") adjustcolor(INK, 0.3) else INK)
    r <- res[res$source == s, ]; mtext(sprintf("%s%s", s, if (s == "NCI-H2052") " (1 NTC)" else ""), side = 3, adj = 0, line = 1.3, font = 2, cex = 0.8)
    mtext(sprintf("NTC-KO %+.3f%s", r$diff_mean, if (is.finite(r$p)) sprintf(", p=%.2g", r$p) else ""), side = 3, adj = 0, line = 0.2, cex = 0.65, col = if (r$NTC_higher) "#d62728" else INK2) }
  mtext(sprintf("%s score, NTC vs %s, screen and each bulk line (NTC-KO; red label = NTC higher)", MOD, KO), outer = TRUE, side = 3, adj = 0, line = 1.5, font = 2, cex = 1.0)
})
