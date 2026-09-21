#!/bin/bash
cd "$(dirname "$0")"
bsub -P acc_Tsankov_Normal_Lung -q premium -n 4 \
     -R "rusage[mem=16000] span[hosts=1]" -W 4:00 \
     -J R2Q14_sarc -o sarc_atac.out -e sarc_atac.err \
     "source /hpc/packages/minerva-rocky9/miniforge3/24.7.1-2/miniforge/etc/profile.d/conda.sh && conda activate meso_scatac && Rscript $PWD/R2_Q14_sarcscore_scATAC.R"
