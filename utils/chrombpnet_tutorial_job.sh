#!/bin/bash
#BSUB -J bias_model_tutorial
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 1
#BSUB -W 24:00
#BSUB -gpu num=1
#BSUB -R v100
#BSUB -R rusage[mem=32000]
#BSUB -R span[hosts=1]
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

ml anaconda3/2022.10
ml cuda/11.7.0;
ml cudnn/8.9.5-11;
ml proxies;
ml java/11.0.2;


source activate chrombpnet


cd /sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet/tutorial/

rm -r bias_model

chrombpnet bias pipeline \
        -ibam data/downloads/merged.bam \
        -d "ATAC" \
        -g ../../genome_references/hg38.genome.fa \
        -c ../hg38.chrom.sizes \
        -p data/peaks_no_blacklist.bed \
        -n data/output_negatives.bed \
        -fl data/splits/fold_0.json \
        -b 0.5 \
        -o bias_model \
        -fp k562

