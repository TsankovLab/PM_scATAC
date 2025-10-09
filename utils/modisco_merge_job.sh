#!/bin/bash
#BSUB -J modisco_merge
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 8
#BSUB -W 48:00
#BSUB -gpu num=1
#BSUB -R v100
#BSUB -R rusage[mem=8000]
#BSUB -R span[hosts=1]
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

ml anaconda3/2022.10

output_dir=${1}
model_head=${2}
cluster_key=${3}
modisco_dir=${4}
contribs_dir=${5}
batch=${6}

active activate chrombpnet 
#source /sc/arion/work/giottb01/conda/etc/profile.d/conda.sh
cd /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/HDMA
    
/sc/arion/work/giottb01/conda/envs/chrombpnet/bin/python -u ../git_repo/utils/03-merge_modisco.py \
                    --out-dir ${output_dir}/merged_motifs \
                    --model-head ${model_head} \
                    --cluster-key ${cluster_key} \
                    --modisco-dir ${modisco_dir} \
                    --contribs-dir ${contribs_dir} \
                    --batch ${batch}
