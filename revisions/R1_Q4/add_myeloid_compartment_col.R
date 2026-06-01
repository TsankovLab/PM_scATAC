# add_myeloid_compartment_col.R  (R1_Q4)
#
# Adds celltype_main = "Myeloid" to every cell in the myeloid ArchR project
# so that the reference-mapping projection can select all cells in one call
# (SUBSET_COL=celltype_main SUBSET_VALUE=Myeloid) instead of looping per
# celltype and copying Arrow files 7 times.
#
# Run via: bsub < submit_add_myeloid_compartment.sh

suppressPackageStartupMessages(library(ArchR))

addArchRThreads(1)
addArchRGenome("hg19")

proj_dir <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR"

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
