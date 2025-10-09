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
#output_dir=${3}
model_head=${3}
cd ${chromBPdir}/${celltype}/no_bias_model


# Take mean of contribution bw 
#wiggletools mean fold_0/${celltype}_contribution_scores.counts_scores.bw fold_1/${celltype}_contribution_scores.counts_scores.bw fold_2/${celltype}_contribution_scores.counts_scores.bw fold_3/${celltype}_contribution_scores.counts_scores.bw fold_4/${celltype}_contribution_scores.counts_scores.bw > ${celltype}_cntr_score_counts_mean.bw

source activate finemo

if [ ! -f "finemo_out_${model_head}/hits_counts.tsv" ]; then
    echo "hits_counts.tsv file not found. Running finemo on TFmodisco outputs from contribution score counts ..."
finemo extract-regions-chrombpnet-h5 -c ${celltype}_averaged_contribution_scores_${model_head}.h5 -o motif_calls_${model_head} -w 2000 #fold_$fold_number/${celltype}_contribution_scores.counts_scores.h5 fold_2/${celltype}_contribution_scores.counts_scores.h5 fold_3/${celltype}_contribution_scores.counts_scores.h5 fold_4/${celltype}_contribution_scores.counts_scores.h5 
finemo call-hits -r motif_calls_${model_head}.npz -m compiled_${model_head}/modisco_compiled.h5 -o finemo_out_${model_head} -p ../${celltype}_peakset_all_no_blacklist.bed -J
finemo report -r motif_calls_${model_head}.npz -H finemo_out_${model_head}/hits_counts.tsv -p ../${celltype}_peakset_all_no_blacklist.bed -m compiled_${model_head}/modisco_compiled.h5 -o finemo_out_${model_head}/report -W 2000
else
    echo "hits_counts.tsv file found!"
fi

# if [ ! -f "finemo_out/hits_profile.tsv" ]; then
#     echo "hits_profile.tsv file not found. Running finemo on TFmodisco outputs from contribution score profile ..."
# finemo extract-regions-chrombpnet-h5 -c ${celltype}_averaged_contribution_scores_profile.h5 -o motif_calls_profile -w 2000 #fold_$fold_number/${celltype}_contribution_scores.counts_scores.h5 fold_2/${celltype}_contribution_scores.counts_scores.h5 fold_3/${celltype}_contribution_scores.counts_scores.h5 fold_4/${celltype}_contribution_scores.counts_scores.h5 
# finemo call-hits -r motif_calls_profile.npz -m modisco_profile/modisco_results_profile.h5 -o finemo_out_profile -p ../${celltype}_peakset_all_no_blacklist.bed -J
# finemo report -r motif_calls_profile.npz -H finemo_out_profile/hits_profile.tsv -p ../${celltype}_peakset_all_no_blacklist.bed -m modisco_profile/modisco_results_profile.h5 -o finemo_out_profile/report -W 2000
# else
#     echo "hits_profile.tsv file found!"
# fi
