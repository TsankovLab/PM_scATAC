suppressMessages({ library(Seurat) })
set.seed(1234)
OUT <- Sys.getenv("CROP_OUT", unset = "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CROPseq_validation_bulkRNA/processed/de/perturbation_panel/figs_repro")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

srt <- readRDS("/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/CRISPR_cropseq_analysis/_cellranger_raw_Filter_400_800_25/no_harmony/srt_filtered.rds")
cat("cells loaded:", ncol(srt), "\n")
cat("unique_guide table:\n"); print(table(srt$unique_guide, useNA="ifany"))

# --- same filters as CRISPR_cropseq.R ---
srt <- srt[, srt$unique_guide %in% TRUE]
res2 <- as.data.frame(table(srt$crispr_calls)); res2 <- res2[order(-res2$Freq), ]
gene_left <- as.character(res2$Var1[res2$Freq > 9])
srt <- srt[, srt$crispr_calls %in% gene_left]

merged_guides <- list(
  MEF2D = c('MEF2D-1','MEF2D-3','MEF2D−3|NTC−3','MEF2D−1|MEF2D−3','MEF2D−1|NTC−3','MEF2D−1|MEF2D−3|NTC−3','MEF2D−1|NTC−1'),
  MEF2A = c('MEF2A−3','MEF2A−1','MEF2A−3|NTC−3','MEF2A−1|NTC−3','MEF2A−3|NTC−1'),
  PITX1 = c('PITX1−2','PITX1−3','PITX1−1','PITX1−2|PITX1−3','PITX1−1|PITX1−2','PITX1−1|PITX1−3','NTC-1|PITX1-2','NTC-3|PITX1-2','NTC-3|PITX1-3','NTC-1|PITX1-3'),
  BPTF = c('BPTF−3','BPTF−2'), HMGA1 = 'HMGA1−3', TEAD2 = 'TEAD2-1', NTC = c('NTC-1','NTC-3'),
  TCF3 = c('TCF3-1','TCF3-3'), SOX9 = 'SOX9-3', TWIST1 = c('TWIST1-1','TWIST1-3'), TEAD4 = c('TEAD4-2','TEAD4-3'))
merged_guides <- lapply(merged_guides, function(x) gsub("[−–-]", "_", x))
srt$merged_call <- srt$crispr_calls
for (i in names(merged_guides)) srt$merged_call[srt$merged_call %in% merged_guides[[i]]] <- i
srt <- srt[, !grepl('\\|', srt$merged_call)]
cat("merged_call counts:\n"); print(table(srt$merged_call))

# --- variable genes (nfeat = 5000, as in the original) ---
srt <- FindVariableFeatures(srt, nfeatures = 5000, verbose = FALSE)
vf <- VariableFeatures(srt)
cat("n variable features:", length(vf), "\n")
writeLines(vf, file.path(OUT, "screen_top5000_variable_genes.txt"))

avg <- AverageExpression(srt, group.by = "merged_call", features = vf, verbose = FALSE)[[1]]
avg <- log2(as.matrix(avg) + 1)
colnames(avg) <- gsub("-", "_", colnames(avg))
write.csv(avg, file.path(OUT, "screen_var5000_log2avg_by_KO.csv"))
write.csv(as.data.frame(table(srt$merged_call)), file.path(OUT, "screen_ncells_by_KO.csv"), row.names = FALSE)

# --- cell cycle index: cc = S.Score + G2M.Score ---
md <- srt@meta.data[, c("merged_call", "cc", "S.Score", "G2M.Score", "Phase")]
cat("cc == S+G2M ? max abs diff:", max(abs(md$cc - (md$S.Score + md$G2M.Score))), "\n")
write.csv(md, file.path(OUT, "screen_cellcycle_per_cell.csv"))

# which cc gene list reproduces stored scores?
tmp <- srt; tmp <- NormalizeData(tmp, verbose = FALSE)
for (nm in c("cc.genes.updated.2019", "cc.genes")) {
  g <- get(nm, envir = asNamespace("Seurat")) 
  tmp2 <- CellCycleScoring(tmp, s.features = g$s.genes, g2m.features = g$g2m.genes, set.ident = FALSE)
  cat(nm, ": cor S =", cor(tmp2$S.Score, srt$S.Score), " cor G2M =", cor(tmp2$G2M.Score, srt$G2M.Score), "\n")
}
