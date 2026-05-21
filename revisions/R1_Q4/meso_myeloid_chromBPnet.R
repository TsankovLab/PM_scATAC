# meso_myeloid_chromBPnet.R  (R1_Q4)
#
# Runs ChromBPNet for three myeloid populations:
#   1. momac         - monocyte-derived macrophages  (myeloid scatac_ArchR, archp$momac)
#   2. resident      - tissue-resident macrophages   (myeloid scatac_ArchR, archp$momac)
#   3. PBMC_CD14_Mono - blood CD14 monocytes         (meso_vs_pbmc ArchR)
#
# Uses Xmen pipeline utils (chromBPnet_call_peaks.R + chromBPnet_master.sh).
# Fragment extraction + peak calling runs here; training/contrib jobs are bsub'd.
# Output: /sc/arion/scratch/giottb01/scatac_meso/revisions/R1_Q4
#
# Run on login node:
#   conda activate /sc/arion/work/giottb01/conda/envs/meso_scatac
#   Rscript meso_myeloid_chromBPnet.R

set.seed(1234)
suppressPackageStartupMessages(library(ArchR))
addArchRThreads(threads = 8)
addArchRGenome("hg38")

repodir    <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/Xmen/Xmen'
grefdir    <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet'
biasdir    <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet/NKT_cells/bias_model'
chromBPdir <- '/sc/arion/scratch/giottb01/scatac_meso/revisions/R1_Q4'

myeloid_projdir <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR'
pbmc_projdir    <- '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/meso_vs_pbmc'

dir.create(chromBPdir, recursive = TRUE, showWarnings = FALSE)
system(paste0('chmod +x ', file.path(repodir, 'utils', 'chromBPnet_master.sh')))

submit_master <- function(celltype) {
  logpfx  <- file.path(chromBPdir, paste0('cBP_master_', celltype))
  command <- paste(
    "bsub",
    "-J",  paste0(celltype, '_CBPmaster'),
    "-P acc_Tsankov_Normal_Lung",
    "-q premium",
    "-n 8",
    "-W 144:00",
    "-R 'rusage[mem=64000]'",
    "-R 'span[hosts=1]'",
    "-o", paste0(logpfx, '.out'),
    "-e", paste0(logpfx, '.err'),
    file.path(repodir, 'utils', 'chromBPnet_master.sh'),
    chromBPdir, grefdir, repodir, celltype, biasdir
  )
  system(command)
  message(sprintf("  Submitted master job for: %s", celltype))
}

# ---- 1. momac / resident (myeloid scatac_ArchR) -------------------------------
message("=== Loading myeloid project ===")
archp <- loadArchRProject(myeloid_projdir, showLogo = FALSE)
archp <- archp[!is.na(archp$momac)]
message(sprintf("  Total cells: %d | momac: %d | resident: %d",
  nrow(archp@cellColData),
  sum(archp$momac == "momac"),
  sum(archp$momac == "resident")))

metaGroupName <- 'momac'
if (exists('fragments')) rm(fragments)
message("=== Extracting fragments and calling peaks: momac / resident ===")
source(file.path(repodir, 'utils', 'chromBPnet_call_peaks.R'))

message("=== Submitting ChromBPNet master jobs: momac / resident ===")
for (ct in c('momac', 'resident')) submit_master(ct)

# ---- 2. PBMC_CD14_Mono (meso_vs_pbmc) ----------------------------------------
message("=== Loading meso_vs_pbmc project ===")
archp <- loadArchRProject(pbmc_projdir, showLogo = FALSE)

ct_col   <- as.character(archp$celltype_unified)
tiss_col <- as.character(archp$tissue)
archp$tissue_ct <- ifelse(
  ct_col %in% c("CD14_Mono", "CD16_Mono", "DC"),
  paste0(tiss_col, "_", ct_col),
  ct_col
)
archp <- archp[archp$tissue_ct == "PBMC_CD14_Mono"]
message(sprintf("  PBMC_CD14_Mono cells: %d", nrow(archp@cellColData)))

metaGroupName <- 'tissue_ct'
if (exists('fragments')) rm(fragments)
message("=== Extracting fragments and calling peaks: PBMC_CD14_Mono ===")
source(file.path(repodir, 'utils', 'chromBPnet_call_peaks.R'))

message("=== Submitting ChromBPNet master job: PBMC_CD14_Mono ===")
submit_master('PBMC_CD14_Mono')

message("All jobs submitted. Done.")
