#!/usr/bin/env Rscript
# R1_Q3_SOX9_module_footprint_archr.R
#
# Model-free Tn5 footprint cross-check of the ChromBPNet result: instead of the
# synthetic ChromBPNet marginal footprint, aggregate the OBSERVED Tn5 insertion
# (fragment) profile at SOX9 motif instances (the same motif positions used by
# chromVAR, Annotations/Motif-Positions-In-Peaks.rds) in the SOX9-high vs
# SOX9-low module-gated cells, using ArchR getFootprints / plotFootprints.
#
# Positive controls (CTCF, an AP-1 and an NFI motif) are footprinted alongside
# SOX9 to confirm the assay recovers a clear Tn5 footprint where one exists.
#
# ArchR writes coverages/plots to outputDirectory; we redirect that to scratch
# (never into the project tree) and copy the final PDFs to R1_Q3 afterwards.

suppressPackageStartupMessages(library(ArchR))
# attach the BSgenome explicitly: ArchR's k-mer Tn5-bias step does
# eval(parse(text="BSgenome.Hsapiens.UCSC.hg38")), which needs the object
# attached, not merely installed (addArchRGenome alone is insufficient).
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
addArchRGenome("hg38"); addArchRThreads(8)

PROJ    <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR_SOX9_module_P23"
SCRATCH <- "/sc/arion/scratch/giottb01/meso_SOX9_module_footprint"
OUTR1Q3 <- "/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q3"
GROUP   <- "SOX9_module_gate"
dir.create(SCRATCH, showWarnings = FALSE, recursive = TRUE)

message("Loading ArchR project: ", PROJ)
proj <- loadArchRProject(PROJ, showLogo = FALSE)

# redirect ALL ArchR output to scratch (arrow files keep their absolute paths)
proj@projectMetadata$outputDirectory <- SCRATCH
message("Output redirected to: ", getOutputDirectory(proj))
message("Gate counts:"); print(table(proj@cellColData[[GROUP]]))

# ---- motif positions used by chromVAR ----
positions <- getPositions(proj, name = "Motif")
allnm <- names(positions)
sox_idx <- grep("SOX9", allnm, ignore.case = TRUE, value = TRUE)
message("SOX9 motif entries found: ", paste(sox_idx, collapse = ", "))

pick1 <- function(pat) { h <- grep(pat, allnm, ignore.case = TRUE, value = TRUE); if (length(h)) h[1] else NULL }
controls <- unique(na.omit(c(pick1("^CTCF"), pick1("^FOS"), pick1("^JUN"), pick1("NFI"))))
message("Positive-control motifs: ", paste(controls, collapse = ", "))

sel <- unique(c(sox_idx, controls))
sel <- sel[sel %in% allnm]
stopifnot(length(sel) > 0)
saveRDS(sel, file.path(SCRATCH, "footprint_motifs_selected.rds"))

# ---- pseudobulk coverages for the two gate groups ----
message("addGroupCoverages on ", GROUP)
proj <- addGroupCoverages(proj, groupBy = GROUP, minCells = 25, maxCells = 1000,
                          minReplicates = 2, maxReplicates = 4, force = TRUE)

# ---- footprints (observed Tn5 insertions around motif sites) ----
message("getFootprints (flank=250) for: ", paste(sel, collapse = ", "))
seFoot <- getFootprints(ArchRProj = proj, positions = positions[sel],
                        groupBy = GROUP, flank = 250)
saveRDS(seFoot, file.path(SCRATCH, "seFoot_SOX9_module.rds"))

# ---- plots: Tn5-bias-subtracted and divided ----
for (nm in c("Subtract", "Divide")) {
  tryCatch(
    plotFootprints(seFoot = seFoot, ArchRProj = proj, normMethod = nm,
                   plotName = paste0("SOX9_module_footprints_", nm),
                   addDOC = FALSE, smoothWindow = 5),
    error = function(e) message("plotFootprints(", nm, ") failed: ", conditionMessage(e)))
}

# ---- dump raw SOX9 footprint matrices to CSV (for quantification in R1_Q3) ----
for (m in sox_idx) {
  tryCatch({
    if (!m %in% names(assays(seFoot))) stop("no assay")
    mat  <- assays(seFoot)[[m]]              # rows: positions (Bias + per-group); cols vary
    safe <- gsub("[^A-Za-z0-9]+", "_", m)
    write.csv(as.data.frame(mat),
              file.path(OUTR1Q3, paste0("SOX9_footprint_obsTn5_", safe, ".csv")),
              row.names = TRUE)
    message("wrote CSV for ", m)
  }, error = function(e) message("CSV dump for ", m, " skipped: ", conditionMessage(e)))
}

# ---- copy plot PDFs into R1_Q3 ----
pdfs <- list.files(file.path(SCRATCH, "Plots"),
                   pattern = "SOX9_module_footprints_.*\\.pdf$", full.names = TRUE)
for (p in pdfs) file.copy(p, file.path(OUTR1Q3, basename(p)), overwrite = TRUE)
message("Copied PDFs to R1_Q3: ", paste(basename(pdfs), collapse = ", "))
message("DONE")
