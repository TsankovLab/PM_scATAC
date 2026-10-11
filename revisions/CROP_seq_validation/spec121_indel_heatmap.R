# Spec 121: heatmap of the fraction of read pairs with an indel at the gRNA cut sites, target gene x bulk sample.
# value = % of read pairs covering the Cas9 cut site (+/- 10 bp) that carry an insertion or deletion there (spec118), best guide of the gene (guides with >= MINP pairs in the sample; grey = no guide with enough pairs)
# rows = the 16 target genes with guides in "Oligos NR4A2+CROPseq.xlsx"; columns = bulk samples, split by genotyped cell line (name-derived for samples outside the spec71 audit), ordered by KO label
# black frame = the gene is the sample's KO label; top annotation = KO label and the gene whose guides were detected in the sample's reads (spec119, >= 300 guide + scaffold reads)
# needs: spec118 (site counts), spec119 (guide calls), spec71 audit
# outputs: spec121_indel_matrix.csv, fig37_indel_heatmap_gene_x_sample.{png,pdf}
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(here, "repro_common.R")); source(file.path(here, "profiles.R")); source(file.path(here, "agree_common.R"))
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(grid) }); MINP <- 10; MIN_GUIDE_READS <- 300
C <- read.csv(file.path(OUT, "spec118_site_counts.csv"), stringsAsFactors = FALSE); St <- read.csv(file.path(OUT, "spec118_gRNA_sites.csv")); St <- St[!is.na(St$chr), ]; C <- C[paste(C$guide, C$site) %in% paste(St$guide, sprintf("%s:%d", St$chr, St$cut)), ]
C$gene <- sub("_g.*", "", C$guide); C <- C[!C$gene %in% c("NTC", "NR4A2enh"), ]; q <- load_audit(); PFX <- c(H2052 = "NCI-H2052", MSTO = "MSTO-211H", H2452 = "NCI-H2452", H28 = "NCI-H28")
ids <- grep("^(H2052|MSTO|H2452|H28)_[A-Za-z0-9]+_R[0-9]+$", unique(C$sample), value = TRUE); SM <- data.frame(sample = ids, line = PFX[sub("_.*", "", ids)], KO = sub("^[^_]+_([^_]+)_.*", "\\1", ids), stringsAsFactors = FALSE)
i <- match(SM$sample, q$sample); SM$line[!is.na(i)] <- q$line[i[!is.na(i)]]; SM$KO[!is.na(i)] <- q$KO[i[!is.na(i)]]; gc <- read.csv(file.path(OUT, "spec119_guide_calls.csv"), stringsAsFactors = FALSE)
j <- match(SM$sample, gc$sample); SM$guides_detected <- ifelse(!is.na(j) & gc$guide_reads[j] >= MIN_GUIDE_READS & gc$top_share[j] >= 0.8, gc$top_target[j], "none")
GENES <- c(KOS, "RUNX1", "RUNX2", "SNAI2", "MEF2A", "MEF2D", "TEAD1", "TEAD3", "TFAP2A"); KOL <- c("NTC", KOS, "RUNX1", "RUNX2", "SNAI2"); SM <- SM[order(match(SM$line, LINES), match(SM$KO, KOL), SM$sample), ]
M <- sapply(SM$sample, function(s) sapply(GENES, function(g) { x <- C[C$sample == s & C$gene == g & C$n_cover >= MINP, ]; if (!nrow(x)) NA else max(100 * x$n_indel / x$n_cover) })); write.csv(M, file.path(OUT, "spec121_indel_matrix.csv"))
own <- outer(GENES, SM$KO, "=="); cat("own-target cells: median", round(median(M[own], na.rm = TRUE), 1), "% | all other cells: median", round(median(M[!own], na.rm = TRUE), 2), "%, 99th percentile", round(quantile(M[!own], 0.99, na.rm = TRUE), 1), "%\n")
off <- which(!own & !is.na(M) & M >= 5, arr.ind = TRUE); cat("cells >= 5% where the gene is not the sample label:\n"); print(data.frame(gene = GENES[off[, 1]], sample = SM$sample[off[, 2]], label = SM$KO[off[, 2]], guides_detected = SM$guides_detected[off[, 2]], indel_pct = round(M[off], 1)), row.names = FALSE)
kc <- c(KO_COL[intersect(names(KO_COL), c(KOL, GENES))], none = "white"); kc <- kc[!duplicated(names(kc))]; mx <- 50
save_fig("fig37_indel_heatmap_gene_x_sample", 19, 6.6, function() {
  ha <- HeatmapAnnotation(`KO label` = SM$KO, `guides detected` = SM$guides_detected, col = list(`KO label` = kc[unique(SM$KO)], `guides detected` = kc[unique(SM$guides_detected)]), annotation_name_gp = gpar(fontsize = 8), simple_anno_size = unit(3.2, "mm"),
                          annotation_legend_param = list(`KO label` = list(nrow = 2, title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)), `guides detected` = list(nrow = 2, title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7))), gp = gpar(col = "grey70", lwd = 0.3))
  ht <- Heatmap(pmin(M, mx), name = "% read pairs\nwith an indel", col = colorRamp2(c(0, 2, 15, mx), c("#f7f6f2", "#f7f6f2", "#f0a58f", "#a8380b")), na_col = "grey82", cluster_rows = FALSE, cluster_columns = FALSE, column_split = factor(sub("NCI-|-211H", "", SM$line), levels = sub("NCI-|-211H", "", LINES)),
                top_annotation = ha, row_names_side = "left", row_names_gp = gpar(fontsize = 9, fontface = "bold", col = ifelse(GENES %in% names(KO_COL), KO_COL[GENES], "grey20")), column_names_gp = gpar(fontsize = 6), column_title_gp = gpar(fontsize = 10, fontface = "bold"), rect_gp = gpar(col = "white", lwd = 0.5), border = TRUE,
                cell_fun = function(j, i, x, y, w, h, fill) { if (GENES[i] == SM$KO[j]) grid.rect(x, y, w, h, gp = gpar(col = "black", fill = NA, lwd = 1.1)); v <- M[i, j]; if (!is.na(v) && v >= 5) grid.text(sprintf("%.0f", v), x, y, gp = gpar(fontsize = 5.2, col = if (v > 25) "white" else "black")) },
                heatmap_legend_param = list(at = c(0, 10, 20, 30, 40, 50), labels = c("0", "10", "20", "30", "40", ">=50"), title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7), direction = "horizontal"))
  draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE, padding = unit(c(2, 4, 12, 4), "mm"))
  grid.text("Indels at the gRNA cut sites: target gene x bulk sample", x = unit(4, "mm"), y = unit(1, "npc") - unit(4, "mm"), just = "left", gp = gpar(fontsize = 12, fontface = "bold"))
  grid.text(sprintf("%% of read pairs covering the Cas9 cut site (+/- 10 bp) with an insertion or deletion, best guide of the gene (>= %d pairs; grey = too few reads); numbers shown from 5%%; black frame = the sample's KO label", MINP), x = unit(4, "mm"), y = unit(1, "npc") - unit(9, "mm"), just = "left", gp = gpar(fontsize = 7.5, col = "grey30"))
})

## ---------------------------------------------------------------- guide reads (spec119) above the indel heatmap, same samples and gene order
# guide reads = reads with the protospacer followed by the sgRNA scaffold among the unmapped read pairs, summed over the 3 guides of the gene, per million unmapped reads (log10 colour scale)
# output: spec121_guide_reads_matrix.csv, fig37b_guide_reads_and_indels.{png,pdf}
gr <- read.csv(file.path(OUT, "spec119_guide_reads.csv"), stringsAsFactors = FALSE); gr$gene <- sub("_g.*", "", gr$guide)
R <- sapply(SM$sample, function(s) { x <- gr[gr$sample == s, ]; sapply(GENES, function(g) 1e6 * sum(x$reads[x$gene == g]) / max(1, x$unmapped_reads[1])) }); write.csv(R, file.path(OUT, "spec121_guide_reads_matrix.csv"))
cat("guide reads per million unmapped: own-target median", round(median(R[own]), 1), "| other cells median", round(median(R[!own]), 2), ", 99th percentile", round(quantile(R[!own], 0.99), 1), "\n")
cat("Spearman between guide reads and indel % over the own-target cells:", round(cor(R[own], M[own], method = "spearman", use = "complete.obs"), 2), "\n")
save_fig("fig37b_guide_reads_and_indels", 19, 11.4, function() {
  lines_f <- factor(sub("NCI-|-211H", "", SM$line), levels = sub("NCI-|-211H", "", LINES)); rn <- gpar(fontsize = 9, fontface = "bold", col = ifelse(GENES %in% names(KO_COL), KO_COL[GENES], "grey20"))
  ha <- HeatmapAnnotation(`KO label` = SM$KO, col = list(`KO label` = kc[unique(SM$KO)]), annotation_name_gp = gpar(fontsize = 8), simple_anno_size = unit(3.2, "mm"), gp = gpar(col = "grey70", lwd = 0.3),
                          annotation_legend_param = list(`KO label` = list(nrow = 2, title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7))))
  frame <- function(Mx, thr, fmt) function(j, i, x, y, w, h, fill) { if (GENES[i] == SM$KO[j]) grid.rect(x, y, w, h, gp = gpar(col = "black", fill = NA, lwd = 1.1)); v <- Mx[i, j]; if (!is.na(v) && v >= thr) grid.text(fmt(v), x, y, gp = gpar(fontsize = 5, col = "black")) }
  h1 <- Heatmap(log10(R + 1), name = "guide reads per million\nunmapped reads", col = colorRamp2(c(0, 1, 2, 3.5), c("#f7f6f2", "#f7f6f2", "#a9c9ee", "#1f5fb0")), cluster_rows = FALSE, cluster_columns = FALSE, column_split = lines_f, top_annotation = ha, row_names_side = "left", row_names_gp = rn,
                show_column_names = FALSE, column_title_gp = gpar(fontsize = 10, fontface = "bold"), rect_gp = gpar(col = "white", lwd = 0.5), border = TRUE, row_title = "guide reads", row_title_gp = gpar(fontsize = 10, fontface = "bold"), cell_fun = frame(R, 30, function(v) sprintf("%.0f", v)),
                heatmap_legend_param = list(at = 0:3, labels = c("0", "10", "100", "1000"), title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7), direction = "horizontal"))
  h2 <- Heatmap(pmin(M, mx), name = "% read pairs\nwith an indel", col = colorRamp2(c(0, 2, 15, mx), c("#f7f6f2", "#f7f6f2", "#f0a58f", "#a8380b")), na_col = "grey82", cluster_rows = FALSE, cluster_columns = FALSE, column_split = lines_f, row_names_side = "left", row_names_gp = rn,
                column_names_gp = gpar(fontsize = 6), column_title = NULL, rect_gp = gpar(col = "white", lwd = 0.5), border = TRUE, row_title = "indels at the cut site", row_title_gp = gpar(fontsize = 10, fontface = "bold"), cell_fun = frame(M, 5, function(v) sprintf("%.0f", v)),
                heatmap_legend_param = list(at = c(0, 10, 20, 30, 40, 50), labels = c("0", "10", "20", "30", "40", ">=50"), title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7), direction = "horizontal"))
  draw(h1 %v% h2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE, padding = unit(c(2, 4, 12, 4), "mm"))
  grid.text("Guide reads detected and indels at the cut site: target gene x bulk sample", x = unit(4, "mm"), y = unit(1, "npc") - unit(4, "mm"), just = "left", gp = gpar(fontsize = 12, fontface = "bold"))
  grid.text("top: reads carrying a guide of the gene followed by the sgRNA scaffold, per million unmapped reads (numbers from 30); bottom: % of read pairs with an indel at the best guide of the gene (numbers from 5%; grey = too few reads); black frame = the sample's KO label",
            x = unit(4, "mm"), y = unit(1, "npc") - unit(9, "mm"), just = "left", gp = gpar(fontsize = 7.5, col = "grey30"))
})
