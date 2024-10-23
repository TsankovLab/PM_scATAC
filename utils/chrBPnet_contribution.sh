#!/bin/bash
#BSUB -J chrBP_model
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


source activate chrombpnet

# Get contribution score bigwigs
cd $chromBPdir

MODEL_H5=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet/${celltype}_model/fold_$fold_number/models/chrombpnet_nobias.h5
REGIONS=peakset_$celltype.bed
GENOME=$grefdir/genome_references/hg38.genome.fa
CHROM_SIZES=$grefdir/hg38.chrom.sizes
OUTPUT_PREFIX=${celltype}_model/fold_$fold_number/${celltype}_contribution_scores
chrombpnet contribs_bw -m $MODEL_H5 -r $REGIONS -g $GENOME -c $CHROM_SIZES -op $OUTPUT_PREFIX 
