#!/bin/bash
#BSUB -J P1_stroma
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 4
#BSUB -R "rusage[mem=24000]"
#BSUB -W 2:00
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/P1_analysis_stroma/logs/P1_stroma_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/P1_analysis_stroma/logs/P1_stroma_%J.err

module purge
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

/sc/arion/work/giottb01/conda/envs/meso_scatac/bin/Rscript \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/P1_analysis_stroma/P1_stroma_recluster_WT1.R
