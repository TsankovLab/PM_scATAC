# add_myeloid_compartment_col.R  (R1_Q4)
#
# Adds celltype_main = "Myeloid" to every cell in the myeloid ArchR project
# so that the reference-mapping projection can select all cells in one call
# (SUBSET_COL=celltype_main SUBSET_VALUE=Myeloid) instead of looping per
# celltype and copying Arrow files 7 times.
#
# Run via: bsub < submit_add_myeloid_compartment.sh

source(local({a <- grep("^--file=", commandArgs(FALSE), value = TRUE); for (d in c(if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))), getwd())) { while (!file.exists(file.path(d, "config.R")) && d != dirname(d)) d <- dirname(d); if (file.exists(file.path(d, "config.R"))) break }; Sys.getenv("MESO_CONFIG", file.path(d, "config.R"))}))   # data / results / reference roots: see config.R at the repo root
suppressPackageStartupMessages(library(ArchR))

addArchRThreads(1)
addArchRGenome("hg19")

proj_dir <- paste0(DATA_DIR, "/myeloid_cells/scatac_ArchR")

message("Loading: ", proj_dir)
archp <- loadArchRProject(proj_dir, showLogo = FALSE)
message("  Cells: ", nCells(archp))

archp <- addCellColData(
  ArchRProj = archp,
  data      = rep("Myeloid", nCells(archp)),
  name      = "celltype_main",
  cells     = getCellNames(archp),
  force     = TRUE
)

message("celltype_main added — unique values: ", paste(unique(archp$celltype_main), collapse = ", "))

saveArchRProject(archp, outputDirectory = proj_dir, overwrite = TRUE, load = FALSE)
message("Saved.")
