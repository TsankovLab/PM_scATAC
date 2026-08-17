D <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/bulkRNA_meso"
bl   <- readRDS(file.path(D, "bulk_RNA_studies.rds"))
meta <- readRDS(file.path(D, "bulk_RNA_studies_metadata.rds"))

genes <- c("BAP1","RUNX1","RUNX2","RUNX3","CBFB")
for (s in names(bl)) {
  m <- bl[[s]]
  cat("\n==== ", s, " : ", nrow(m), "genes x", ncol(m), "samples ====\n")
  cat("range of values:", paste(round(range(m, na.rm=TRUE),2), collapse=" .. "), "\n")
  present <- genes[genes %in% rownames(m)]
  cat("target genes present:", paste(present, collapse=", "), "\n")
  cat("meta rows:", nrow(meta[[s]]), " ; ncol:", ncol(meta[[s]]), "\n")
  cat("colnames(meta) matching BAP/mut/alter/status/subtype:\n")
  print(grep("BAP|mut|Mut|alter|Alter|status|Status|subtype|driver|CDKN|NF2",
             colnames(meta[[s]]), value = TRUE))
  # sample name alignment
  cat("all expr cols in meta rownames? ",
      all(colnames(m) %in% rownames(meta[[s]])), "\n")
}
cat("\nDONE\n")
