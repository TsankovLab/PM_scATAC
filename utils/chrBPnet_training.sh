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
biasdir=${5}
echo $biasdir

#mkdir $chromBPdir
mkdir $chromBPdir/$celltype
cd ${chromBPdir}/${celltype}


# Train chrombpnet model
MODEL_H5=no_bias_model/fold_${fold_number}/models/chrombpnet_nobias.h5

if [ ! -f "${MODEL_H5}" ]; then
    echo "chrombpnet_nobias.h5 file not found. Training chromBPnet with bias correction model..."
rm -r no_bias_model/fold_$fold_number/
chrombpnet pipeline \
    -ifrag ../fragments_${celltype}.tsv \
    -d "ATAC" \
    -g $grefdir/genome_references/hg38.genome.fa \
    -c $grefdir/hg38.chrom.sizes \
    -p ${celltype}_peakset_all_no_blacklist.bed \
    -n  no_bias_model/output_negatives_f${fold_number}_negatives.bed \
    -fl $grefdir/folds/fold_$fold_number.json \
    -b $biasdir/fold_${fold_number}/models/bias.h5 \
    -o no_bias_model/fold_$fold_number
else
    echo "chrombpnet_nobias.h5 file found!"
fi

# Compute contribution scores
REGIONS=${celltype}_peakset_all_no_blacklist.bed
GENOME=$grefdir/genome_references/hg38.genome.fa
CHROM_SIZES=$grefdir/hg38.chrom.sizes
count_scores_file=no_bias_model/fold_${fold_number}/contribution_scores.counts_scores.h5

if [ ! -f "${count_scores_file}" ]; then
    echo "Contribution scores file not found. Computing contribution scores..."
    chrombpnet contribs_bw -m $MODEL_H5 -r $REGIONS -g $GENOME -c $CHROM_SIZES -op no_bias_model/fold_${fold_number}/
else
    echo "Contribution scores file already exists: ${count_scores_file}"
fi

#chrombpnet contribs_bw -m $MODEL_H5 -r $REGIONS -g $GENOME -c $CHROM_SIZES -op $OUTPUT_PREFIX 




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