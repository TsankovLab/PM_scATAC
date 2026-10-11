# Spec 38: pull a given cNMF module (MODULE_N env var, default 16) from a given k solution (CNMF_K env var, default 30) in the
# CROP-seq screen's own cNMF sweep run, rank its per-cell usage across all screen KO groups (cnmf.usages.k_<K>), then take its
# top-defining genes (cnmf.gene_spectra_score.k_<K>) and score them as a module (Seurat::AddModuleScore) in bulk RNA-seq, both
# pooled across lines and per line, to check whether the same KO group stands out in bulk too.
# needs: CRISPR_cropseq_analysis/.../cNMF_crispr_20_30_vf5000/cnmf/cnmf.usages.k_<K>.dt_0_3.consensus.txt,
#        cnmf.gene_spectra_score.k_<K>.dt_0_3.txt, spec_common.R (load_screen()), repro_common.R (load_full_audited_bulk())
# outputs: spec38_module<N>_top_genes_k<K>.csv, spec38_module<N>_usage_by_KO_screen_k<K>.csv, spec38_module<N>_bulk_scores_k<K>.csv,
#          fig2W_k<K>_module<N>_usage_screen.{png,pdf}, fig2W_k<K>_module<N>_AddModuleScore_bulk.{png,pdf},
#          fig2W_k<K>_module<N>_AddModuleScore_bulk_per_line.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "spec_common.R")); source(file.path(here, "profiles.R"))
suppressPackageStartupMessages(library(Seurat))
CNMF_D <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CRISPR_cropseq_analysis/_cellranger_raw_Filter_400_800_25/no_harmony/cNMF/cNMF_crispr_20_30_vf5000/cnmf"
MODULE <- Sys.getenv("MODULE_N", "16")
K <- Sys.getenv("CNMF_K", "30")
SUF <- paste0("_k", K); TAG <- paste0("module", MODULE)
cat(sprintf("=== cNMF module %s, k=%s ===\n", MODULE, K))

## ---------------------------------------------------------------- module's top-defining genes (gene_spectra_score, k=K)
gss <- read.delim(file.path(CNMF_D, sprintf("cnmf.gene_spectra_score.k_%s.dt_0_3.txt", K)), row.names = 1, check.names = FALSE)
mN <- sort(setNames(as.numeric(gss[MODULE, ]), colnames(gss)), decreasing = TRUE)
top_genes <- data.frame(gene = names(mN)[1:50], gene_spectra_score = mN[1:50])
write.csv(top_genes, file.path(OUT, paste0("spec38_", TAG, "_top_genes", SUF, ".csv")), row.names = FALSE)
cat(sprintf("\n=== cNMF module %s (k=%s) top 20 defining genes ===\n", MODULE, K)); print(head(top_genes, 20), row.names = FALSE)

## ---------------------------------------------------------------- usage of this module across screen KO groups
usages <- read.delim(file.path(CNMF_D, sprintf("cnmf.usages.k_%s.dt_0_3.consensus.txt", K)), row.names = 1, check.names = FALSE)
srt <- load_screen(); grp_scr <- setNames(as.character(srt$merged_call), colnames(srt))
common <- intersect(rownames(usages), names(grp_scr))
cat(sprintf("\ncNMF usages cells: %d | screen (filtered) cells: %d | overlap: %d\n", nrow(usages), length(grp_scr), length(common)))
uN <- usages[common, MODULE]; g <- grp_scr[common]
usage_by_ko <- sort(sapply(split(uN, g), mean), decreasing = TRUE)
write.csv(data.frame(group = names(usage_by_ko), mean_usage = usage_by_ko, n_cells = sapply(split(uN, g), length)[names(usage_by_ko)]),
          file.path(OUT, paste0("spec38_", TAG, "_usage_by_KO_screen", SUF, ".csv")), row.names = FALSE)
cat(sprintf("\n=== mean cNMF module %s (k=%s) usage per screen KO group (highest first) ===\n", MODULE, K)); print(round(usage_by_ko, 4))
top_ko_screen <- setdiff(names(usage_by_ko), "NTC")[1]
cat(sprintf("\ntop non-NTC group: %s (rank %d of %d)\n", top_ko_screen, match(top_ko_screen, names(usage_by_ko)), length(usage_by_ko)))

save_fig(paste0("fig2W", SUF, "_", TAG, "_usage_screen"), 8, 5, function() {
  ord <- names(usage_by_ko); cols <- ifelse(ord == top_ko_screen, KO_COL[top_ko_screen], ifelse(ord == "NTC", "grey60", adjustcolor(INK2, 0.5)))
  par(mar = c(6, 4.5, 3.4, 1), family = "sans")
  bp <- barplot(usage_by_ko[ord], names.arg = ord, col = cols, border = INK2, las = 2, cex.names = 0.85,
                ylab = sprintf("mean cNMF module %s usage", MODULE), ylim = c(0, max(usage_by_ko) * 1.15))
  mtext(sprintf("cNMF module %s (k=%s) usage per CROP-seq KO group", MODULE, K), side = 3, adj = 0, line = 1.8, font = 2, cex = 1.0)
  mtext(sprintf("bars = mean per-cell usage (cnmf.usages.k_%s consensus); top non-NTC group (%s) highlighted", K, top_ko_screen), side = 3, adj = 0, line = 0.6, cex = 0.6, col = INK2)
})

## ---------------------------------------------------------------- score this module's top genes as an AddModuleScore module in bulk RNA-seq
fab <- load_full_audited_bulk(); meta <- fab$meta; lg <- fab$lg
g_bulk <- top_genes$gene[top_genes$gene %in% rownames(lg)]
cat(sprintf("\n%d/%d module-%s top genes present in bulk\n", length(g_bulk), nrow(top_genes), MODULE))

options(Seurat.object.assay.version = "v3")   # AddModuleScore errors ("subscript out of bounds") on the default Assay5 class here; classic Assay works
srt_bulk <- CreateSeuratObject(counts = as.matrix(fab$cs)); srt_bulk <- SetAssayData(srt_bulk, slot = "data", new.data = as.matrix(lg))
srt_bulk <- AddModuleScore(srt_bulk, features = list(g_bulk), name = "mN", seed = 1)
score <- setNames(srt_bulk$mN1, colnames(srt_bulk))[rownames(meta)]

by_ko <- sort(sapply(split(score, meta$Gene_Control), mean), decreasing = TRUE)
write.csv(data.frame(sample = names(score), Sample_ID = meta$Sample_ID, Cell_Line = meta$Cell_Line, Gene_Control = meta$Gene_Control, module_score = score),
          file.path(OUT, paste0("spec38_", TAG, "_bulk_scores", SUF, ".csv")), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== mean module-%s (k=%s) AddModuleScore per bulk KO condition, all 4 lines pooled (highest first) ===\n", MODULE, K)); print(round(by_ko, 4))
top_ko_bulk <- setdiff(names(by_ko), "NTC")[1]
cat(sprintf("\ntop non-NTC group in bulk: %s (rank %d of %d) | screen's top non-NTC group (%s) ranks %d of %d in bulk\n",
            top_ko_bulk, match(top_ko_bulk, names(by_ko)), length(by_ko), top_ko_screen, match(top_ko_screen, names(by_ko)), length(by_ko)))

by_ko_line <- do.call(rbind, lapply(LINES, function(cl) {
  m <- meta[meta$Cell_Line == cl, ]; s <- score[rownames(m)]
  x <- sort(sapply(split(s, m$Gene_Control), mean), decreasing = TRUE)
  data.frame(line = cl, KO = names(x), mean_score = x, rank = seq_along(x))
}))
cat(sprintf("\n=== per line: %s (screen's top non-NTC group) rank ===\n", top_ko_screen))
print(by_ko_line[by_ko_line$KO == top_ko_screen, ], row.names = FALSE)

save_fig(paste0("fig2W", SUF, "_", TAG, "_AddModuleScore_bulk"), 8, 5, function() {
  ord <- names(by_ko); cols <- ifelse(ord == top_ko_screen, KO_COL[top_ko_screen], ifelse(ord == "NTC", "grey60", adjustcolor(INK2, 0.5)))
  par(mar = c(6, 4.5, 3.4, 1), family = "sans")
  ylim <- range(c(0, by_ko)) + c(-1, 1) * 0.15 * diff(range(c(0, by_ko)))
  bp <- barplot(by_ko[ord], names.arg = ord, col = cols, border = INK2, las = 2, cex.names = 0.85,
                ylab = "mean AddModuleScore", ylim = ylim)
  abline(h = 0, col = INK)
  mtext(sprintf("cNMF module %s (k=%s) top-50 genes, scored in bulk RNA-seq (AddModuleScore)", MODULE, K), side = 3, adj = 0, line = 1.8, font = 2, cex = 1.0)
  mtext(sprintf("bars = mean module score per bulk KO condition, 4 lines pooled; %s (screen's top group) highlighted", top_ko_screen), side = 3, adj = 0, line = 0.6, cex = 0.6, col = INK2)
})

## ---------------------------------------------------------------- same thing, one panel per cell line
save_fig(paste0("fig2W", SUF, "_", TAG, "_AddModuleScore_bulk_per_line"), 11, 8, function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE)); par(oma = c(0, 0, 3, 0), family = "sans")
  yr <- range(c(0, by_ko_line$mean_score))
  for (cl in LINES) {
    d <- by_ko_line[by_ko_line$line == cl, ]; d <- d[order(-d$mean_score), ]
    cols <- ifelse(d$KO == top_ko_screen, KO_COL[top_ko_screen], ifelse(d$KO == "NTC", "grey60", adjustcolor(INK2, 0.5)))
    par(mar = c(6, 4.5, 3, 1))
    barplot(d$mean_score, names.arg = d$KO, col = cols, border = INK2, las = 2, cex.names = 0.8,
            ylab = "mean AddModuleScore", ylim = yr + c(-1, 1) * 0.1 * diff(yr))
    abline(h = 0, col = INK)
    mtext(cl, side = 3, adj = 0, line = 0.6, font = 2, cex = 0.95, col = LINE_COL[[cl]])
  }
  mtext(sprintf("cNMF module %s (k=%s) top-50 genes, scored in bulk RNA-seq, per cell line", MODULE, K), outer = TRUE, side = 3, adj = 0, line = 1.6, font = 2, cex = 1.0)
  mtext(sprintf("bars = mean module score per bulk KO condition within that line; %s highlighted", top_ko_screen), outer = TRUE, side = 3, adj = 0, line = 0.4, cex = 0.6, col = INK2)
})
