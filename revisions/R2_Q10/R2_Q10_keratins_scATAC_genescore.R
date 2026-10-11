###############################################################################
# R2_Q10 -- keratin GENE SCORES in epithelioid tumour cells compared with normal mesothelium (scATAC-seq)
#
# The scRNA analysis (R2_Q10_keratins_epithelioid_vs_mesothelium.R) repeated on the chromatin side.
# Data: main/scatac_ArchR, cells labelled Malignant and Mesothelium (celltype_lv1); keratin "expression" is the ArchR
# gene activity score (GeneScoreMatrix), i.e. accessibility around the gene, not RNA.
#
# Epithelioid call from the scS-score computed ON THE scATAC CELLS, with the recipe already used for scATAC in
# R2_Q14_sarcscore_scATAC.R / scatac_tumor_SOX9.R: ArchR addModuleScore over the gene scores for the 25 malignant
# cNMF programmes (top 30 genes, nBin 25, nBgd 100), z-scored across the programmes within each cell, cNMF20 = scS.
# Tumour call: epithelioid when the median scS of the tumour's malignant cells is below 0. No keratin is among the cNMF20 genes.
#
# Comparison as in the scRNA script: % of cells with a non-zero score, mean score, log2 fold change, Wilcoxon over
# cells (BH over the keratins), and per-patient pseudobulk.
# Two things specific to gene scores: (1) they scale with DNA copy number, so a gain or loss of 12q13 / 17q21 (where
# the keratin genes sit in two clusters) in a tumour moves all keratins of that cluster together; (2) neighbouring
# genes of a cluster share accessibility windows, so their scores are not independent.
# The main project is not modified (output directory redirected, no saveArchRProject).
#
# Output: R2_Q10_scATAC_scS_per_sample.csv, R2_Q10_scATAC_keratin_stats.csv, R2_Q10_scATAC_keratin_pseudobulk_per_sample.csv,
#         R2_Q10_keratin_RNA_vs_ATAC.csv, Plots/R2_Q10_keratins_scATAC_genescore.pdf
###############################################################################
suppressPackageStartupMessages({ library(ArchR); library(Matrix); library(ggplot2); library(patchwork) })
addArchRGenome("hg38"); addArchRThreads(4); options(width = 220); set.seed(1)
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"; OUT <- file.path(ROOT, "git_repo_claude", "R2_Q10")
setwd(OUT); dir.create("Plots", showWarnings = FALSE); dir.create("archr_scratch", showWarnings = FALSE)
SARC <- "cNMF20"; TOP_GENES <- 30; NBIN <- 25; NBGD <- 100; MIN_GENES <- 5; MIN_PCT <- 0.05; MIN_CELLS <- 10

proj <- loadArchRProject(file.path(ROOT, "main", "scatac_ArchR"), showLogo = FALSE); proj@projectMetadata$outputDirectory <- file.path(OUT, "archr_scratch")
cd0 <- getCellColData(proj); proj <- proj[rownames(cd0)[cd0$celltype_lv1 %in% c("Malignant", "Mesothelium")], ]
cat("cells:", nCells(proj), "\n"); print(table(proj$Sample, proj$celltype_lv1))

## ---- scS on the malignant cells ------------------------------------------------------------------
spectra <- readRDS(file.path(ROOT, "git_repo", "files", "malignant_cnmf_genelist_25_nfeat_5000.rds")); gnames <- getFeatures(proj, useMatrix = "GeneScoreMatrix")
spectra <- lapply(spectra, function(x) { x <- head(x, TOP_GENES); x[x %in% gnames] }); spectra <- spectra[lengths(spectra) >= MIN_GENES]; stopifnot(SARC %in% names(spectra))
cat(SARC, "genes:", paste(spectra[[SARC]], collapse = ", "), "\nkeratins among them:", paste(grep("^KRT", spectra[[SARC]], value = TRUE), collapse = ", "), "\n")
MS_RDS <- file.path("archr_scratch", "module_scores_25programmes.rds")                     # the module scores take ~15 min; cached
if (file.exists(MS_RDS)) { md <- readRDS(MS_RDS); stopifnot(identical(rownames(md), getCellNames(proj))) } else {
  proj@cellColData <- proj@cellColData[, !colnames(proj@cellColData) %in% names(spectra)]
  proj <- addModuleScore(proj, useMatrix = "GeneScoreMatrix", name = "", features = spectra, nBin = NBIN, nBgd = NBGD, seed = 1, threads = getArchRThreads(), logFile = createLogFile("addModuleScore"))
  colnames(proj@cellColData) <- gsub("^\\.", "", colnames(proj@cellColData)); md <- as.data.frame(getCellColData(proj)); saveRDS(md, MS_RDS) }
md$sample <- as.character(md$Sample)
md$scS <- t(scale(t(as.matrix(md[, names(spectra)]))))[, SARC]
mal <- md[md$celltype_lv1 == "Malignant", ]; ss <- data.frame(sample = names(table(mal$sample)), n_malignant = as.integer(table(mal$sample))); ss$median_scS <- tapply(mal$scS, mal$sample, median)[ss$sample]; ss <- ss[order(ss$median_scS), ]
## Tumour call: epithelioid when the median scS of the tumour's malignant cells is below SCS_CUT = 0. The score is a z across the 25
## programmes within each cell, so 0 is "the sarcomatoid programme is no more active than the cell's average programme".
## (The automatic break rule of the scRNA script is not used here: on these medians it splits after the two lowest tumours.)
SCS_CUT <- 0
rna <- read.csv("R2_Q10_scS_per_sample.csv"); ss$scRNA_median_scS <- rna$median_scS[match(ss$sample, rna$sample)]; ss$scRNA_call <- rna$call[match(ss$sample, rna$sample)]
ss$call <- ifelse(ss$median_scS < SCS_CUT, "epithelioid", "higher scS")
ss$n_mesothelium <- as.integer(table(factor(md$sample[md$celltype_lv1 == "Mesothelium"], levels = ss$sample)))
write.csv(ss, "R2_Q10_scATAC_scS_per_sample.csv", row.names = FALSE)
cat(sprintf("\nscS cut-off: %.2f (tumour median below it = epithelioid)\n", SCS_CUT)); print(ss, row.names = FALSE, digits = 3)
okc <- !is.na(ss$scRNA_call); cat(sprintf("same call as from the tumour's own scRNA: %d of %d tumours with both\n", sum(ss$call[okc] == ss$scRNA_call[okc]), sum(okc)))
ok <- !is.na(ss$scRNA_median_scS); cat(sprintf("scATAC vs scRNA tumour scS: Spearman %.2f over %d tumours\n", cor(ss$median_scS[ok], ss$scRNA_median_scS[ok], method = "spearman"), sum(ok)))
EPI <- ss$sample[ss$call == "epithelioid"]; GL <- c("Normal mesothelium", "Epithelioid tumour cells", "Higher-scS tumour cells")
md$group <- ifelse(md$celltype_lv1 == "Mesothelium", GL[1], ifelse(md$sample %in% EPI, GL[2], GL[3])); md$group_cell <- ifelse(md$celltype_lv1 == "Mesothelium", GL[1], ifelse(md$scS < SCS_CUT, GL[2], GL[3]))
cat("\ncells per group (tumour-level call):\n"); print(table(md$group)[GL]); cat("cells per group (cell-level call):\n"); print(table(md$group_cell)[GL])

## ---- keratin gene scores ---------------------------------------------------------------------------
se <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix"); X <- assay(se); rownames(X) <- rowData(se)$name; X <- X[, rownames(md)]
kr <- sort(grep("^KRT[0-9]", rownames(X), value = TRUE)); pct <- sapply(GL, function(g) Matrix::rowMeans(X[kr, md$group == g, drop = FALSE] > 0)); kr <- kr[apply(pct, 1, max) >= MIN_PCT]
cat("\nkeratin genes with a non-zero gene score in >=", 100 * MIN_PCT, "% of the cells of a group:", length(kr), "\n")
stat <- function(gcol) do.call(rbind, lapply(kr, function(g) { a <- X[g, md[[gcol]] == GL[2]]; b <- X[g, md[[gcol]] == GL[1]]
  data.frame(gene = g, pct_epithelioid = 100 * mean(a > 0), pct_mesothelium = 100 * mean(b > 0), mean_epithelioid = mean(a), mean_mesothelium = mean(b), log2FC = log2((mean(a) + 0.01) / (mean(b) + 0.01)), p_wilcox = suppressWarnings(wilcox.test(a, b)$p.value)) }))
ST <- stat("group"); ST$fdr <- p.adjust(ST$p_wilcox, "BH"); S2 <- stat("group_cell"); ST$log2FC_cell_level_call <- S2$log2FC; ST$fdr_cell_level_call <- p.adjust(S2$p_wilcox, "BH")
md$unit <- paste(md$sample, ifelse(md$celltype_lv1 == "Mesothelium", "Mesothelium", "Malignant"), sep = "|"); tu <- table(md$unit); units <- names(tu)[tu >= MIN_CELLS]
PBl <- log2(sapply(units, function(u) Matrix::rowMeans(X[kr, md$unit == u, drop = FALSE])) + 0.01)
ui <- data.frame(unit = units, sample = sub("\\|.*", "", units), compartment = sub(".*\\|", "", units), n_cells = as.integer(tu[units])); ui$call <- ifelse(ui$compartment == "Mesothelium", "Normal mesothelium", ss$call[match(ui$sample, ss$sample)]); ui$scS <- ss$median_scS[match(ui$sample, ss$sample)]
write.csv(cbind(ui, t(round(PBl, 3))), "R2_Q10_scATAC_keratin_pseudobulk_per_sample.csv", row.names = FALSE)
ue <- ui$unit[ui$call == "epithelioid"]; um <- ui$unit[ui$call == "Normal mesothelium"]
ST$pseudobulk_log2FC <- rowMeans(PBl[, ue, drop = FALSE]) - rowMeans(PBl[, um, drop = FALSE]); ST$n_epi_tumours_above_all_meso <- rowSums(PBl[, ue, drop = FALSE] > apply(PBl[, um, drop = FALSE], 1, max)); ST$n_epi_tumours_below_all_meso <- rowSums(PBl[, ue, drop = FALSE] < apply(PBl[, um, drop = FALSE], 1, min))
gr <- getGenes(proj); ST$chr <- as.character(seqnames(gr))[match(ST$gene, gr$symbol)]; ST <- ST[order(-ST$log2FC), ]; write.csv(ST, "R2_Q10_scATAC_keratin_stats.csv", row.names = FALSE)
cat(sprintf("\n=== keratin gene scores: epithelioid tumour cells (%d cells, %d tumours) vs normal mesothelium (%d cells, %d patients with >= %d cells) ===\n", sum(md$group == GL[2]), length(ue), sum(md$group == GL[1]), length(um), MIN_CELLS)); print(format(ST, digits = 3), row.names = FALSE)
cat(sprintf("\ncell level: %d of %d keratins higher in epithelioid cells (FDR < 0.05), %d lower\n", sum(ST$fdr < 0.05 & ST$log2FC > 0), nrow(ST), sum(ST$fdr < 0.05 & ST$log2FC < 0)))
cat("by chromosome (keratin clusters): higher / lower / n.s.\n"); print(table(ST$chr, ifelse(ST$fdr < 0.05, ifelse(ST$log2FC > 0, "higher", "lower"), "n.s.")))

## ---- agreement with the scRNA result ------------------------------------------------------------------
R <- read.csv("R2_Q10_keratin_stats.csv"); C <- merge(R[, c("gene", "log2FC", "fdr", "pseudobulk_log2FC")], ST[, c("gene", "log2FC", "fdr", "pseudobulk_log2FC")], by = "gene", suffixes = c("_RNA", "_ATAC"))
C$same_direction <- sign(C$log2FC_RNA) == sign(C$log2FC_ATAC); C$both_sig_same <- C$same_direction & C$fdr_RNA < 0.05 & C$fdr_ATAC < 0.05; write.csv(C, "R2_Q10_keratin_RNA_vs_ATAC.csv", row.names = FALSE)
cat(sprintf("\n=== scRNA vs scATAC, keratins tested in both: %d ; same direction %d ; Spearman of log2FC %.2f (p %.3g); significant in both with the same sign: %s ===\n", nrow(C), sum(C$same_direction), cor(C$log2FC_RNA, C$log2FC_ATAC, method = "spearman"),
            cor.test(C$log2FC_RNA, C$log2FC_ATAC, method = "spearman", exact = FALSE)$p.value, paste(sprintf("%s (%s)", C$gene[C$both_sig_same], ifelse(C$log2FC_RNA[C$both_sig_same] > 0, "up", "down")), collapse = ", "))); print(format(C[order(-C$log2FC_RNA), ], digits = 3), row.names = FALSE)

## ---- figure ------------------------------------------------------------------------------------------
th <- theme_classic(base_size = 9); mal$sample <- factor(mal$sample, levels = ss$sample); mal$call <- ss$call[match(mal$sample, ss$sample)]; gord <- ST$gene
pA <- ggplot(mal, aes(sample, scS, fill = call)) + geom_violin(scale = "width", linewidth = 0.25) + geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.25) + geom_hline(yintercept = SCS_CUT, linetype = 2, colour = "red") +
  scale_fill_manual(values = c(epithelioid = "plum4", `higher scS` = "grey45"), name = NULL) + th + theme(legend.position = "top", legend.key.size = unit(0.3, "cm")) + labs(x = NULL, y = "scS-score from gene scores (malignant cells)", title = "A  Tumours called by the scATAC scS-score", subtitle = sprintf("red = cut-off %.1f on the tumour median", SCS_CUT))
dd <- do.call(rbind, lapply(GL, function(g) data.frame(group = g, gene = kr, pct = 100 * Matrix::rowMeans(X[kr, md$group == g, drop = FALSE] > 0), mean = Matrix::rowMeans(X[kr, md$group == g, drop = FALSE])))); dd$gene <- factor(dd$gene, levels = rev(gord)); dd$group <- factor(dd$group, levels = GL)
dd$z <- ave(dd$mean, dd$gene, FUN = function(v) (v - mean(v)) / ifelse(sd(v) > 0, sd(v), 1))
pB <- ggplot(dd, aes(group, gene)) + geom_point(aes(size = pct, fill = z), shape = 21, colour = "grey30", stroke = 0.25) + scale_size(range = c(0.3, 4), name = "% cells\nscore > 0") + scale_fill_gradient2(low = "#2a78d6", mid = "white", high = "#b2182b", name = "mean score\n(z per gene)") +
  scale_x_discrete(labels = c("Normal\nmesothelium", "Epithelioid\ntumour", "Higher-scS\ntumour")) + th + theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 6), legend.key.size = unit(0.3, "cm")) + labs(x = NULL, y = NULL, title = "B  Keratin gene scores by group", subtitle = "genes ordered by log2FC epithelioid vs mesothelium")
ST$gene <- factor(ST$gene, levels = rev(gord)); ST$sig <- ifelse(ST$fdr < 0.05, ifelse(ST$log2FC > 0, "higher in epithelioid", "higher in mesothelium"), "n.s.")
pC <- ggplot(ST, aes(log2FC, gene, fill = sig)) + geom_col(width = 0.7) + geom_point(aes(x = pseudobulk_log2FC), shape = 23, fill = "white", size = 1.3) + geom_vline(xintercept = 0, linewidth = 0.3) + scale_fill_manual(values = c(`higher in epithelioid` = "plum4", `higher in mesothelium` = "olivedrab3", `n.s.` = "grey80"), name = "cell-level FDR < 0.05") +
  th + theme(legend.position = "top", legend.key.size = unit(0.3, "cm"), legend.title = element_text(size = 7), axis.text.y = element_text(size = 6)) + guides(fill = guide_legend(nrow = 2)) + labs(x = "log2 fold change of the gene score, epithelioid / mesothelium", y = NULL, title = "C  Epithelioid vs normal mesothelium", subtitle = "bar = all cells; diamond = mean of per-patient pseudobulks")
pE <- ggplot(C, aes(log2FC_RNA, log2FC_ATAC)) + geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") + geom_point(aes(fill = both_sig_same), shape = 21, size = 2.2) + ggrepel::geom_text_repel(aes(label = gene), size = 2.3, max.overlaps = Inf) +
  scale_fill_manual(values = c(`TRUE` = "plum4", `FALSE` = "grey85"), name = "significant in both,\nsame direction") + th + theme(legend.key.size = unit(0.3, "cm"), legend.title = element_text(size = 7)) + labs(x = "scRNA: log2FC expression", y = "scATAC: log2FC gene score", title = "D  scRNA vs scATAC", subtitle = sprintf("Spearman %.2f over %d keratins", cor(C$log2FC_RNA, C$log2FC_ATAC, method = "spearman"), nrow(C)))
ui <- ui[order(ui$call != "Normal mesothelium", ui$scS), ]; hd <- do.call(rbind, lapply(seq_len(nrow(ui)), function(i) data.frame(lab = sprintf("%s\n%s (n=%d)", ui$sample[i], ifelse(ui$compartment[i] == "Mesothelium", "meso", "tumour"), ui$n_cells[i]), call = ui$call[i], gene = kr, v = PBl[, ui$unit[i]])))
hd$lab <- factor(hd$lab, levels = unique(hd$lab)); hd$gene <- factor(hd$gene, levels = rev(gord)); hd$z <- ave(hd$v, hd$gene, FUN = function(v) (v - mean(v)) / ifelse(sd(v) > 0, sd(v), 1)); hd$call <- factor(ifelse(hd$call == "epithelioid", "Epithelioid tumours", ifelse(hd$call == "Normal mesothelium", "Normal mesothelium", "Higher-scS tumours")), levels = c("Normal mesothelium", "Epithelioid tumours", "Higher-scS tumours"))
pD <- ggplot(hd, aes(lab, gene, fill = z)) + geom_tile(colour = "white", linewidth = 0.3) + facet_grid(. ~ call, scales = "free_x", space = "free_x") + scale_fill_gradient2(low = "#2a78d6", mid = "white", high = "#b2182b", name = "pseudobulk\n(z per gene)") +
  th + theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), strip.background = element_blank(), strip.text = element_text(size = 7, face = "bold"), legend.key.size = unit(0.3, "cm")) + labs(x = NULL, y = NULL, title = "E  Per patient (pseudobulk gene score)", subtitle = sprintf("patients with >= %d cells of the compartment; tumours ordered by scS", MIN_CELLS))
pdf("Plots/R2_Q10_keratins_scATAC_genescore.pdf", width = 15, height = 11); print((pA | pB | pC | pE) / pD + plot_layout(heights = c(1, 1.15))); dev.off()
cat("\nDONE\n")
