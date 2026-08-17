#!/usr/bin/env Rscript
# R1_Q3_SOX9_module_P23_atac.R
#
# AddModuleScore of the SOX9 program (top-20 genes from P4 tumor scRNA cluster 32)
# on the P23 scATAC project, scored over the GeneScoreMatrix, then FeaturePlot
# on the P23 ATAC UMAP. Reference panels: SOX9_cluster (C9 high / C4 low) and a
# violin of the module score per cluster.

suppressMessages({ library(ArchR); library(ggplot2); library(patchwork) })
addArchRThreads(1)

PROJ <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
PDF  <- file.path(OUT, "Plots", "SOX9_program_module_P23_atac.pdf")

genes <- read.csv(file.path(OUT, "cluster32_top20_genes.csv"))$gene
message("SOX9 program (", length(genes), " genes): ", paste(genes, collapse = ", "))

proj <- loadArchRProject(PROJ, showLogo = FALSE)
message("cells: ", nCells(proj))

## which genes are present in the GeneScoreMatrix
gsm  <- getMatrixFromProject(proj, "GeneScoreMatrix")
avail <- intersect(genes, rowData(gsm)$name)
miss  <- setdiff(genes, avail)
message("genes in GeneScoreMatrix: ", length(avail), "/", length(genes),
        if (length(miss)) paste0("   (missing: ", paste(miss, collapse=","), ")") else "")

## ArchR module score over gene-activity
proj <- addModuleScore(proj, useMatrix = "GeneScoreMatrix",
                       name = "SOX9prog",
                       features = list(SOX9prog = avail), seed = 1)
sc_name <- grep("SOX9prog", colnames(getCellColData(proj)), value = TRUE)[1]
message("module column: ", sc_name)

## pick a UMAP embedding
embs <- names(proj@embeddings)
emb  <- grep("UMAP", embs, value = TRUE)
emb  <- if (length(emb)) emb[1] else embs[1]
message("embedding: ", emb)

## per-cluster summary
cd <- as.data.frame(getCellColData(proj, select = c("SOX9_cluster", sc_name)))
colnames(cd)[2] <- "module"
cat("\n#### SOX9 program module score by cluster ####\n")
print(tapply(cd$module, cd$SOX9_cluster, summary))
wt <- wilcox.test(module ~ SOX9_cluster, data = cd)
cat(sprintf("Wilcoxon C9(high) vs C4(low): p = %.2e\n", wt$p.value))

## plots
p_mod <- plotEmbedding(proj, colorBy = "cellColData", name = sc_name,
                       embedding = emb, plotAs = "points", size = 0.4) +
         ggtitle("SOX9 program module (gene-score) - P23 ATAC")
p_cl  <- plotEmbedding(proj, colorBy = "cellColData", name = "SOX9_cluster",
                       embedding = emb, plotAs = "points", size = 0.4) +
         ggtitle("SOX9 high (C9) vs low (C4)")
p_vln <- ggplot(cd, aes(SOX9_cluster, module, fill = SOX9_cluster)) +
         geom_violin(scale = "width", trim = TRUE) +
         geom_boxplot(width = 0.12, outlier.size = 0.2, fill = "white") +
         scale_fill_manual(values = c(SOX9_high_P23="#E41A1C", SOX9_low_P23="#377EB8")) +
         labs(title = "SOX9 program by cluster", x = NULL, y = "module score") +
         theme_classic(base_size = 11) + theme(legend.position = "none")

dir.create(dirname(PDF), showWarnings = FALSE, recursive = TRUE)
pdf(PDF, width = 7, height = 6)
print(p_mod); print(p_cl); print(p_vln)
dev.off()
message("\nWrote ", PDF)

# also dump per-cell scores
write.csv(data.frame(cell = rownames(cd), cd), file.path(OUT, "SOX9_program_module_P23_percell.csv"),
          row.names = FALSE)
