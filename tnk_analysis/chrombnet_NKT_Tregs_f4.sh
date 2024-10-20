#!/bin/bash
#BSUB -J chrbp_model_Tregs_f4
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


cd /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet


chrombpnet pipeline \
    -ifrag fragments_Tregs.tsv \
    -d "ATAC" \
    -g ../../../../../genome_references/hg38.genome.fa \
    -c ../../../../../chromBPnet/hg38.chrom.sizes \
    -p peakset_Tregs.bed \
    -n  output_negatives.bed \
    -fl folds/fold_4.json \
    -b bias_model/models/tnk_bias.h5 \
    -o Tregs_model_f4/

