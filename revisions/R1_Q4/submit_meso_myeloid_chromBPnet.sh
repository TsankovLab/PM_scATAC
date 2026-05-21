#!/bin/bash
#BSUB -J meso_myeloid_cBP_setup
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 8
#BSUB -W 12:00
#BSUB -R "rusage[mem=32000]"
#BSUB -R "span[hosts=1]"
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/meso_myeloid_cBP_setup_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/meso_myeloid_cBP_setup_%J.err

mkdir -p /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs

source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh
conda activate /sc/arion/work/giottb01/conda/envs/meso_scatac

export OPENBLAS_NUM_THREADS=1

Rscript /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/meso_myeloid_chromBPnet.R
