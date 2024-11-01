#!/bin/bash
#BSUB -J finemo
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



chromBPdir=${1}
echo $chromBPdir
celltype=${2}
echo $celltype
fold_number=${3}
echo $fold_number


source activate chrombpnet

cd ${chromBPdir}/${celltype}_model/

# Take mean of contribution bw 
wiggletools mean fold_0/${celltype}_contribution_scores.counts_scores.bw fold_1/${celltype}_contribution_scores.counts_scores.bw fold_2/${celltype}_contribution_scores.counts_scores.bw fold_3/${celltype}_contribution_scores.counts_scores.bw fold_4/${celltype}_contribution_scores.counts_scores.bw > ${celltype}_cntr_score_counts_mean.bw


source activate finemo
finemo extract-regions-chrombpnet-h5 -c fold_0/${celltype}_contribution_scores.counts_scores.h5 fold_1/${celltype}_contribution_scores.counts_scores.h5 fold_2/${celltype}_contribution_scores.counts_scores.h5 fold_3/${celltype}_contribution_scores.counts_scores.h5 fold_4/${celltype}_contribution_scores.counts_scores.h5 -o motif_calls -w 2000

finemo call-hits -r motif_calls.npz -m fold_0/modisco/modisco_results.h5 -o finemo_out -p ../peakset_${celltype}.bed -J

finemo report -r motif_calls.npz -H finemo_out/hits.tsv -p ../peakset_${celltype}.bed -m fold_0/modisco/modisco_results.h5 -o finemo_out/report -W 2000