#!/bin/bash
#BSUB -J celltype_barplots
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 2
#BSUB -R "rusage[mem=32000]"
#BSUB -W 2:00
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/celltype_barplots_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/celltype_barplots_%J.err

mkdir -p /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs

module purge

/sc/arion/work/giottb01/conda/envs/meso_scatac/bin/Rscript \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q2/celltype_barplots.R
