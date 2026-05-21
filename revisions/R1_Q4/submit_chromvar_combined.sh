#!/bin/bash
#BSUB -J chromvar_combined
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 36
#BSUB -R "rusage[mem=16000] span[hosts=1]"
#BSUB -W 12:00
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/chromvar_combined_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/chromvar_combined_%J.err

mkdir -p /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs

export OPENBLAS_NUM_THREADS=1

/sc/arion/work/giottb01/conda/envs/meso_scatac/bin/Rscript \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q4/run_chromvar_combined.R
