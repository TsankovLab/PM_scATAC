#!/bin/bash
#BSUB -J myeloid_SE_hub_enrich
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 4
#BSUB -R "rusage[mem=128000]"
#BSUB -W 8:00
#BSUB -o /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/myeloid_SE_hub_enrich_%J.out
#BSUB -e /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs/myeloid_SE_hub_enrich_%J.err

mkdir -p /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/logs

module purge

/sc/arion/work/giottb01/conda/envs/meso_scatac/bin/Rscript \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R1_Q7/myeloid_SE_hub_enrichment.R
