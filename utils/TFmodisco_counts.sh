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


# Get contribution score bigwigs
# Make sure to have version >=2.2.1  installed ny pulling from github like this: pip install git+https://github.com/jmschrei/tfmodisco-lite.git

mkdir ${chromBPdir}/${celltype}_model/modisco_counts/
cd ${chromBPdir}/${celltype}_model/modisco_counts/

if [ ! -f "modisco_results_counts.h5" ]; then
    echo "modisco_results_counts.h5 file not found. Running TFmodisco on contribution score counts ..."
modisco motifs -i ../averaged_contributions_counts.h5 -n 1000000 -o modisco_results_counts.h5 #-w 2000
modisco report -i modisco_results_counts.h5 -o report/ -s report/ -m /sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/HOCOMOCO_db/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
else
    echo "modisco_results_counts file found!"
fi

