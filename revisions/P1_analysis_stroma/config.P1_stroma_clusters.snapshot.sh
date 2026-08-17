#!/bin/bash
# config.sh — Fill in before running either pipeline.
# Sourced automatically by all pipeline scripts; do not run directly.
#
# Sections:
#   1. HPC settings
#   2. Conda environments
#   3. Genome reference & bias model
#   4. Pipeline I/O
#   5. Cell types & ArchR project

# ---- 1. HPC settings --------------------------------------------------------

# LSF project account (-P flag)
HPC_PROJECT="acc_Tsankov_Normal_Lung"

# Path to conda initialisation script on your cluster
CONDA_INIT="/hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh"

# GPU resource strings for LSF (-R flag)
GPU_TRAIN="a100"       # training jobs
GPU_CONTRIB="h100nvl"  # contribution score jobs
GPU_FINEMO="a100"      # FI-NeMo jobs

# ---- 2. Conda environments --------------------------------------------------

CONDA_ENV_CHROMBPNET="chrombpnet"   # chrombpnet + modisco + MACS2
CONDA_ENV_FINEMO="finemo"           # fi-nemo
CONDA_ENV_H5PY="h5py"              # HDF5 averaging (wiggletools / wigToBigWig)
CONDA_ENV_R="meso_scatac"           # R with ArchR + ggplot2

# ---- 3. Genome reference & bias model ---------------------------------------

# Directory with hg38 reference files. Expected structure:
#   genome_references/hg38.genome.fa
#   hg38.chrom.sizes
#   blacklist.bed.gz
#   folds/fold_{0..4}.json
GREFDIR="/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet"

# Pre-trained ATAC-seq bias model. Expected structure:
#   fold_{0..4}/models/bias.h5
BIASDIR="/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet/NKT_cells/bias_model"

# MEME motif database used by modisco report and TOMTOM.
# Auto-detect from chrombpnet env:
#   conda run -n chrombpnet python -c "import chrombpnet,os; print(os.path.join(os.path.dirname(chrombpnet.__file__),'data/motifs.meme.txt'))"
MEME_DB="/sc/arion/work/giottb01/conda/envs/chrombpnet/lib/python3.8/site-packages/chrombpnet/data/motifs.meme.txt"

# HDMA repo: path to HDMA/code/03-chrombpnet (contains chrombpnet_utils module)
# Clone from: https://github.com/GreenleafLab/HDMA
HDMA_UTILS="/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/HDMA/code/03-chrombpnet"

# ---- 4. Pipeline I/O --------------------------------------------------------

# Scratch directory where all ChromBPNet outputs are written (one sub-folder per cell type)
WORKDIR="/sc/arion/scratch/giottb01/scatac_meso/revisions/P1_stroma_chrombpnet"

# Path to this repository (set automatically — do not change)
REPODIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- 5. Cell types & ArchR project (ChromBPNet entry point) -----------------

# ArchR project directory
ARCHR_PROJECT="/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/P1_analysis_stroma/scatac_ArchR_P1_stroma"

# Column in cellColData that holds cell-type labels
CELLTYPE_COLUMN="Clusters"

# Cell types to run (leave empty to auto-detect all unique values in CELLTYPE_COLUMN)
CELLTYPES=("C1" "C2" "C3" "C4")
