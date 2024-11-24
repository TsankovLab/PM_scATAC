#!/bin/bash
#BSUB -J tfmodisco
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q premium
#BSUB -n 1
#BSUB -R rusage[mem=32000]
#BSUB -R span[hosts=1]
#BSUB -W 3:00
#BSUB -o tfmodisco.stdout
#BSUB -eo tfmodisco.stderr
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
celltype=${2}
echo $celltype
fold_number=${3}
echo $fold_number


# Get contribution score bigwigs
# Make sure to have version >=2.2.1  installed ny pulling from github like this: pip install git+https://github.com/jmschrei/tfmodisco-lite.git

mkdir ${chromBPdir}/${celltype}_model/fold_${fold_number}/modisco/
cd ${chromBPdir}/${celltype}_model/fold_${fold_number}/modisco/

modisco motifs -i ../${celltype}_contribution_scores.counts_scores.h5 -n 1000000 -o modisco_results.h5 #-w 2000

modisco report -i modisco_results.h5 -o report/ -s report/ -m /sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/HOCOMOCO_db/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme

