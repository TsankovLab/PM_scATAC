#!/bin/bash
#BSUB -J chrBP_model
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 8
#BSUB -W 48:00
#BSUB -gpu num=1
#BSUB -R v100
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


mkdir $chromBPdir
cd $chromBPdir

rm -r ${celltype}_model/fold_$fold_number/

chrombpnet pipeline \
    -ifrag fragments_$celltype.tsv \
    -d "ATAC" \
    -g $grefdir/genome_references/hg38.genome.fa \
    -c $grefdir/hg38.chrom.sizes \
    -p peakset_$celltype.bed \
    -n  output_negatives.bed \
    -fl $grefdir/folds/fold_$fold_number.json \
    -b bias_model/models/model_bias.h5 \
    -o ${celltype}_model/fold_$fold_number/

