#!/bin/bash
#BSUB -J chrBP_model
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 8
#BSUB -W 48:00
#BSUB -gpu num=2
#BSUB -R h100nvl
#BSUB -R rusage[mem=32000]
#BSUB -R span[hosts=1]
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

ml anaconda3/2022.10
ml cuda/11.7.0
ml cudnn/8.9.5-11
ml proxies
ml java/11.0.2
ml tensorrt/8.5.3.1


source activate chrombpnet

chromBPdir=${1}
echo $chromBPdir
grefdir=${2}
echo $grefdir
celltype=${3}
echo $celltype
fold_number=${4}
echo $fold_number


#mkdir $chromBPdir
cd ${chromBPdir}/${celltype}

# Train chrombpnet bias model
# BIAS_MODEL_H5=${celltype}_bias_model/fold_${fold_number}/models/model_bias.h5
# if [ ! -f "${BIAS_MODEL_H5}" ]; then
#     rm -r ${celltype}_model/fold_$fold_number/
#     echo "run training bias model"

# Train bias model
chrombpnet bias pipeline \
        -ifrag ../fragments_${celltype}.tsv \
        -d "ATAC" \
        -g ${grefdir}/genome_references/hg38.genome.fa \
        -c ${grefdir}/hg38.chrom.sizes \
        -p ${celltype}_peakset_all_no_blacklist.bed \
        -n bias_model/output_negatives_f${fold_number}_negatives.bed \
        -fl ${grefdir}/folds/fold_${fold_number}.json \
        -b 0.5 \
        -o bias_model/fold_${fold_number}/ #\
        #-fp 
       
## Troubleshoot numpy (version installed should be 1.23.4)
# python
# import deepdish as dd
# import numpy as np
# d = {"test": np.array([0,0])}
# dd.io.save("test.h5", d, compression='blosc')


# !!!!!!ATTENTION!!!!!HOW TO USE THE NEW GPUS ADDED!!!!!ATTENTION!!!!!!!!!!!!!!!!!!!!!!!!!
# To submit jobs to NVLinked H100 GPU nodes, flag "-R h100nvl" is required.
# PLEASE ADD  #BSUB -R h100nvl to your LSF script or -R h100nvl to your LSF command line

# To submit jobs to L40S GPU nodes, flag "-R l40s" is required. 
# PLEASE ADD #BSUB -R l40s to your LSF script or -R l40s to your LSF command line. 