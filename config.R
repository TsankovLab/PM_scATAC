# config.R: root folders read and written by the analysis scripts. Every R script sources this file: it is found by searching upward from the
# working directory (run scripts from inside the repository), or set MESO_CONFIG=/path/to/config.R.
# Set the roots once per installation, in one of three ways (later ones win):
#   1. edit the placeholders below;  2. set the environment variables (MESO_DATA_DIR, ...);  3. put assignments in config.local.R next to this file
#      (git-ignored, sourced last).
#   DATA_DIR     project data: ArchR projects, Seurat objects, fragment files, bulk RNA-seq counts, in the subfolders the scripts name
#                (main/, tumor_compartment/, NKT_cells/, myeloid_cells/, stroma/, Endothelial/, B_cells/, CROPseq_validation_bulkRNA/, ...)
#   RESULTS_DIR  outputs of the revision analyses, one subfolder per analysis as in revisions/ (e.g. R1_Q4/output/)
#   REF_DIR      shared reference data: DBs/ (gene sets, motif databases), Public_data/ (published cohorts), chromBPnet/ (ChromBPNet reference
#                files and bias models), Xmen/ (external ArchR references, e.g. greenleaf_PBMC_BM = Granja et al. 2019, GSE139369), TCGA/, ...
#   GENOME_DIR   genome assemblies and annotation: human/hg38.fa, human/gencode.v46.annotation.gtf
#   SCRATCH_DIR  large intermediate files (alignments, model training, caches)
#   CONDA_ENVS   folder containing the conda environments (meso_scatac, chrombpnet, finemo, ...: see environment.yaml, software_used.csv)
#   CONDA_INIT   conda's etc/profile.d/conda.sh, sourced by job scripts that activate environments
#   LSF_ACCOUNT  LSF project / account for scripts that submit jobs with bsub (#BSUB headers carry the placeholder <LSF_ACCOUNT>)
REPO_DIR <- local({ f <- NULL; for (i in rev(seq_len(sys.nframe()))) if (!is.null(sys.frame(i)$ofile)) { f <- sys.frame(i)$ofile; break }
  if (is.null(f)) getwd() else dirname(normalizePath(f)) })
DATA_DIR    <- Sys.getenv("MESO_DATA_DIR",    "/path/to/data")
RESULTS_DIR <- Sys.getenv("MESO_RESULTS_DIR", "/path/to/results")
REF_DIR     <- Sys.getenv("MESO_REF_DIR",     "/path/to/references")
GENOME_DIR  <- Sys.getenv("MESO_GENOME_DIR",  "/path/to/genomes")
SCRATCH_DIR <- Sys.getenv("MESO_SCRATCH_DIR", "/path/to/scratch")
CONDA_ENVS  <- Sys.getenv("MESO_CONDA_ENVS",  "/path/to/conda/envs")
CONDA_INIT  <- Sys.getenv("MESO_CONDA_INIT",  "/path/to/conda/etc/profile.d/conda.sh")
LSF_ACCOUNT <- Sys.getenv("MESO_LSF_ACCOUNT", "<LSF_ACCOUNT>")
if (file.exists(file.path(REPO_DIR, "config.local.R"))) source(file.path(REPO_DIR, "config.local.R"))
local({ unset <- c(DATA_DIR = DATA_DIR, RESULTS_DIR = RESULTS_DIR, REF_DIR = REF_DIR, GENOME_DIR = GENOME_DIR, SCRATCH_DIR = SCRATCH_DIR)
  unset <- names(unset)[startsWith(unset, "/path/to/")]
  if (length(unset)) message("config.R: not set yet: ", paste(unset, collapse = ", "), " (see the header of config.R)") })
