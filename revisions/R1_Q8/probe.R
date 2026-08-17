suppressMessages({
  library(ArchR)
  library(Seurat)
  library(SummarizedExperiment)
})
addArchRThreads(threads = 1)
addArchRGenome("hg38")

base <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment"

cat("=== ArchR project ===\n")
archp <- loadArchRProject(file.path(base, "scatac_ArchR"), showLogo = FALSE)
cc <- as.data.frame(archp@cellColData)
cat("nCells:", nrow(cc), "\n")
cat("Has Sample3:", "Sample3" %in% colnames(cc), "\n")
if ("Sample3" %in% colnames(cc)) print(sort(table(cc$Sample3)))
cat("Available matrices:\n"); print(getAvailableMatrices(archp))

cat("\n=== srt object ===\n")
srt <- readRDS(file.path(base, "scrna", "srt.rds"))
cat("class:", class(srt), "\n")
cat("nCells:", ncol(srt), "  nGenes:", nrow(srt), "\n")
cat("meta cols:\n"); print(grep("ample", colnames(srt@meta.data), value = TRUE))
if ("sampleID3" %in% colnames(srt@meta.data)) print(sort(table(srt$sampleID3)))
cat("Assays:\n"); print(Assays(srt))
cat("DONE\n")
