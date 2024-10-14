#!/bin/bash
#BSUB -J bias_model_tutorial
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 4
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


cd /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet


# Train bias model
chrombpnet bias pipeline \
        -ifrag fragments_NKT_cells.tsv \
        -d "ATAC" \
        -g ../../../../../genome_references/hg38.genome.fa \
        -c ../../../../../chromBPnet/hg38.chrom.sizes \
        -p peakset_all_no_blacklist.bed \
        -n output_negatives.bed \
        -fl splits/fold_0.json \
        -b 0.5 \
        -o bias_model/ \
        -fp tnk

