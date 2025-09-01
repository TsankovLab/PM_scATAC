#!/bin/bash
#BSUB -J chrBP_model
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 8
#BSUB -W 48:00
#BSUB -gpu num=2
#BSUB -R v100
#BSUB -R rusage[mem=32000]
#BSUB -R span[hosts=1]
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

# ml anaconda3/2022.10
ml cuda/11.7.0
ml cudnn/8.9.5-11
ml proxies
ml java/11.0.2
ml tensorrt/8.5.3.1

#ml anaconda3/2020.11
source activate chrombpnet

chromBPdir=${1}
echo $chromBPdir
grefdir=${2}
echo $grefdir
repodir=${3}
echo $repodir
celltype=${4}
echo $celltype
MODEL_H5=${5}
echo $MODEL_H5
#mkdir $chromBPdir
#chromBPdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet
#celltype=low_P23
#MODEL_H5=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet/${celltype}/no_bias_model/fold_0/models/chrombpnet_nobias.h5
#/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet
grefdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet
motif_file=motif_footprints.txt

chromBPdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet
celltype=SOX9_low_P23
MODEL_H5=${chromBPdir}/${celltype}/no_bias_model/fold_0/models/chrombpnet_nobias.h5
grefdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet
OUTPUT_PREFIX=${chromBPdir}/${celltype}/footprints

cd $chromBPdir/$celltype

mkdir footprints
cd footprints
for fold_number in {0..4}; do
    # motif_file=motif_endothelial.TF.txt

    # mkdir -p $chromBPdir/$celltype   # only needed once, so you can keep it outside the loop if you like

    chrombpnet footprints -m $MODEL_H5 \
        -r ../${celltype}_peakset_all_no_blacklist.bed \
        -g $grefdir/genome_references/hg38.genome.fa \
        -fl $grefdir/folds/fold_${fold_number}.json \
        -op ${OUTPUT_PREFIX}/fold${fold_number} \
        -pwm_f ../../${motif_file} # [-bs BATCH_SIZE] [--ylim YLIM]
done



