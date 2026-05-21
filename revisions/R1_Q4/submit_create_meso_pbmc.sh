#!/bin/bash
#BSUB -J meso_pbmc_archr
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 8
#BSUB -R "rusage[mem=16000] span[hosts=1]"
#BSUB -W 6:00
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/meso_pbmc_archr_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/meso_pbmc_archr_%J.err

mkdir -p /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs

export OPENBLAS_NUM_THREADS=1

/sc/arion/work/giottb01/conda/envs/meso_scatac/bin/Rscript \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/create_meso_pbmc_archr.R
