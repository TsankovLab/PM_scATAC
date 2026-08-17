suppressMessages({
  library(ArchR)
  library(SummarizedExperiment)
})
addArchRThreads(threads = 1)
addArchRGenome("hg38")

base <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment"
archp <- loadArchRProject(file.path(base, "scatac_ArchR"), showLogo = FALSE)
cc <- as.data.frame(archp@cellColData)

cat("nCells:", nrow(cc), "\n")
cat("colnames of cellColData:\n"); print(colnames(cc))

cat("\n=== Sample column ===\n")
if ("Sample" %in% colnames(cc)) print(sort(table(cc$Sample)))
cat("\n=== Sample3 column ===\n")
if ("Sample3" %in% colnames(cc)) print(sort(table(cc$Sample3)))

cat("\n=== Available matrices ===\n")
print(getAvailableMatrices(archp))

# Is BAP1 in GeneScore rownames?
gs <- getMatrixFromProject(archp, useMatrix = "GeneScoreMatrix")
gnames <- rowData(gs)$name
cat("\nBAP1 in GeneScoreMatrix:", "BAP1" %in% gnames, "\n")

# Motif names sample
mm <- getMatrixFromProject(archp, useMatrix = "MotifMatrix")
mnames <- gsub('_.*','', rowData(mm)$name)
cat("N motifs:", length(mnames), "\n")
cat("RUNX motifs present:\n"); print(grep("^RUNX", mnames, value = TRUE))
cat("MotifMatrix assay names:\n"); print(assayNames(mm))
cat("DONE\n")
