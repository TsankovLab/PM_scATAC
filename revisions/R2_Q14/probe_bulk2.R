D <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/bulkRNA_meso"
bl   <- readRDS(file.path(D, "bulk_RNA_studies.rds"))
meta <- readRDS(file.path(D, "bulk_RNA_studies_metadata.rds"))

mm  <- meta[["mesomics"]]
cat("mesomics meta dim:", dim(mm), "\n")
cat("has Sample col:", "Sample" %in% colnames(mm), "\n")
cat("IHC.BAP1 table:\n"); print(table(mm$IHC.BAP1, useNA = "always"))
cat("\nfirst Sample vs expr colnames:\n")
print(head(mm$Sample)); print(head(colnames(bl[["mesomics"]])))
cat("expr cols found in meta$Sample:",
    sum(colnames(bl[["mesomics"]]) %in% mm$Sample), "of", ncol(bl[["mesomics"]]), "\n")

# any BAP1 columns anywhere in tcga/bueno we missed (broaden search)
for (s in c("tcga","bueno")) {
  cat("\n", s, "cols with 'BAP' (any case):\n")
  print(grep("bap", colnames(meta[[s]]), ignore.case = TRUE, value = TRUE))
}
cat("\nDONE\n")
