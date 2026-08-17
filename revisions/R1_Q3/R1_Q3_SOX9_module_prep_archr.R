#!/usr/bin/env Rscript
# R1_Q3_SOX9_module_prep_archr.R
#
# Prep for the purity-gated ChromBPNet iteration: subset the P23 tumor ArchR
# project to the module-gated SOX9+ / SOX9- cells (defined from the SOX9-program
# module score, top-20 P4 cluster-32 genes) and add a SOX9_module_gate column.
# Saves a new project the ChromBPNet pipeline groups on.
#
# SOX9pos = module >= 0.20 (335 cells) ; SOX9neg = module <= -0.10 (631 cells)
# (excluded middle band dropped)

suppressPackageStartupMessages(library(ArchR))
addArchRGenome("hg38"); addArchRThreads(4)

SRC  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_P23"
OUT  <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_module_P23"
GATE <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3/SOX9_module_gated_cells.csv"

message("Loading ArchR project: ", SRC)
archp <- loadArchRProject(SRC, showLogo = FALSE)
message(sprintf("  Loaded %d cells", nCells(archp)))

g <- read.csv(GATE, stringsAsFactors = FALSE)
g <- g[g$gate %in% c("SOX9pos", "SOX9neg"), c("cell", "gate")]
message(sprintf("  Gated cells: %d (pos+neg)", nrow(g)))

keep <- archp$cellNames %in% g$cell
message(sprintf("  Cells in project matching gate list: %d", sum(keep)))
archp_sub <- archp[keep]

archp_sub$SOX9_module_gate <- g$gate[match(archp_sub$cellNames, g$cell)]
message("SOX9_module_gate label counts:")
print(table(archp_sub$SOX9_module_gate))
stopifnot(!any(is.na(archp_sub$SOX9_module_gate)))

message("Saving subset ArchR project to: ", OUT)
saveArchRProject(ArchRProj = archp_sub, outputDirectory = OUT, overwrite = TRUE, load = FALSE)
message("Done.")
