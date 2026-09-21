###############################################################################
# R2_Q14 -- per-tumour x compartment TF readouts from the MAIN scATAC ArchR project
# (all cell types), for testing whether BAP1-associated TF hits sit outside tumour cells.
#
# One project for every compartment, so motif activity and gene score are computed on the
# same peak set / gene model for malignant and non-malignant cells alike.
#   compartments (celltype_lv1): Malignant | TNK = T_cells + NK | Myeloid | pDC |
#     B_Plasma = B_cells + Plasma | Stroma = Fibroblasts + Endothelial + SmoothMuscle |
#     Mesothelium | Alveolar
# Group means per sample x compartment (getGroupSE, divideN):
#   MotifMatrix deviations and z; GeneScoreMatrix (scaleTo 1e4, then log2(x+1)) for TF genes.
# NOTHING is written back into the shared project: output directory redirected here, the
# project is never saved.
# Output: compartment_motif_dev.csv, compartment_motif_z.csv, compartment_genescore_TF.csv,
#         compartment_groups.csv (sample, compartment, n_cells)
###############################################################################
suppressPackageStartupMessages({ library(ArchR); library(SummarizedExperiment) })
addArchRThreads(4); addArchRGenome("hg38")
ROOT <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM"
OUT  <- file.path(ROOT, "git_repo_claude", "R2_Q14"); setwd(OUT)

proj <- loadArchRProject(file.path(ROOT, "main", "scatac_ArchR"), showLogo = FALSE)
proj@projectMetadata$outputDirectory <- file.path(OUT, "archr_scratch")
dir.create(file.path(OUT, "archr_scratch"), showWarnings = FALSE)
cd <- getCellColData(proj)
MAP <- c(Malignant = "Malignant", T_cells = "TNK", NK = "TNK", Myeloid = "Myeloid", pDCs = "pDC",
         B_cells = "B_Plasma", Plasma = "B_Plasma", Fibroblasts = "Stroma", Endothelial = "Stroma",
         SmoothMuscle = "Stroma", Mesothelium = "Mesothelium", Alveolar = "Alveolar")
comp <- unname(MAP[as.character(cd$celltype_lv1)])
cat("unmapped cell types:", paste(setdiff(unique(as.character(cd$celltype_lv1)), names(MAP)), collapse = ", "), "\n")
proj$R2Q14_group <- paste(as.character(cd$Sample), comp, sep = "__")
G <- as.data.frame(table(proj$R2Q14_group)); names(G) <- c("group", "n_cells")
G$sample <- sub("__.*", "", G$group); G$compartment <- sub(".*__", "", G$group)
write.csv(G[, c("group", "sample", "compartment", "n_cells")], "compartment_groups.csv", row.names = FALSE)
cat("cells per compartment:\n"); print(tapply(G$n_cells, G$compartment, sum))

## The group SE of a MotifMatrix holds deviations and z stacked in ONE assay; rowData
## `seqnames` says which block a row belongs to.  Split on it, not on assay names.
se <- getGroupSE(proj, useMatrix = "MotifMatrix", groupBy = "R2Q14_group", divideN = TRUE)
rd <- rowData(se); blk <- as.character(rd$seqnames)
cat("MotifMatrix group SE blocks:", paste(names(table(blk)), table(blk), collapse = " | "), "\n")
for (b in intersect(c("deviations", "z"), unique(blk))) {
  ix <- which(blk == b)
  M <- as.matrix(assay(se))[ix, , drop = FALSE]
  rownames(M) <- make.unique(sub("_[0-9]+$", "", rd$name[ix]))
  write.csv(M, sprintf("compartment_motif_%s.csv", if (b == "deviations") "dev" else "z"))
}
tf <- rownames(M)
if (!file.exists("compartment_genescore_TF.csv")) {
gse <- getGroupSE(proj, useMatrix = "GeneScoreMatrix", groupBy = "R2Q14_group", divideN = TRUE, scaleTo = 10000)
gn <- rowData(gse)$name
paper <- read.csv(file.path(OUT, "affinity_regression", "prepared", "paper_2B.csv"), stringsAsFactors = FALSE)
keep <- which(gn %in% union(sub("\\.[0-9]+$", "", tf), paper$TF) & !duplicated(gn))
GS <- log2(as.matrix(assay(gse))[keep, , drop = FALSE] + 1); rownames(GS) <- gn[keep]
write.csv(GS, "compartment_genescore_TF.csv")
}
cat(sprintf("motif matrix: %d motifs x %d groups\n", length(tf), ncol(se)))
cat("DONE\n")
