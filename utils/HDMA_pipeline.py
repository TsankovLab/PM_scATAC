conda activate chrombpnet

cd /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/HDMA

python code/03-chrombpnet/02-compendium/01-modisco_to_pfm.py \
  -c ../main/scatac_ArchR/chrombpnet_models_paths.tsv \
  -o ../main/scatac_ArchR/chromBPnet/pfms.txt \
  -p pos_patterns \
  -t 0.3 \
  -ml 3


python code/03-chrombpnet/02-compendium/01-modisco_to_pfm.py \
  -c ../main/scatac_ArchR/chrombpnet_models_paths.tsv \
  -o ../main/scatac_ArchR/chromBPnet/pfms_neg.txt \
  -p neg_patterns \
  -t 0.3 \
  -ml 3


# Positive motifs
gimme cluster \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/pfms.txt \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/cluster_pos \
  -t 0.8 \
  -N 8

# Negative motifs
gimme cluster \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/pfms_neg.txt \
  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/cluster_neg \
  -t 0.8 \
  -N 8

output_dir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/modisco_merged_counts
model_head=counts


### Merge modisco motifs per model head ###
#BATCH_SIZE=5
#bias_params="Heart_c0_thresh0.4"
#out_dir=${modisco_merged_dir}
#model_head="counts"

output_dir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/modisco_merged_${model_head}
mkdir $output_dir

model_head=counts

awk -F'\t' 'BEGIN{OFS="\t"} {n=int((NR-1)/3)+1; print $1, "pos_patterns", $2, n}' /sc/arion/scratch/giottb01/chromBPnet/cluster_pos/cluster_key.txt > /sc/arion/scratch/giottb01/chromBPnet/cluster_pos/cluster_key_modified_pos.txt
awk -F'\t' 'BEGIN{OFS="\t"} {n=int((NR-1)/3)+1; print $1, "neg_patterns", $2, n}' /sc/arion/scratch/giottb01/chromBPnet/cluster_neg/cluster_key.txt > /sc/arion/scratch/giottb01/chromBPnet/cluster_neg/cluster_key_modified_neg.txt
cat /sc/arion/scratch/giottb01/chromBPnet/cluster_pos/cluster_key_modified_pos.txt /sc/arion/scratch/giottb01/chromBPnet/cluster_neg/cluster_key_modified_neg.txt > ${output_dir}/cluster_key_combined.txt


cluster_key=${output_dir}/cluster_key_combined.txt


modisco_dir=/sc/arion/scratch/giottb01/chromBPnet
contribs_dir=/sc/arion/scratch/giottb01/chromBPnet

# set inputs
#cluster_key="${gimme_cluster_dir%/}/gimme_cluster_all_cluster_key.tsv"

# SUBMIT JOBS ------------------------------------------------------------------

# get batches
batches=$(cut -f 4 $cluster_key | sort | uniq)

chmod +x ../git_repo/utils/modisco_merge_job.sh

for batch in ${batches[@]}; do

bsub -J modisco_merge \
-P acc_Tsankov_Normal_Lung \
-q premium \
-n 8 \
-W 12:00 \
-R rusage[mem=80000] \
-R span[hosts=1] \
-o ${output_dir}/modisco_merge_${batch}.out \
-e ${output_dir}/modisco_merge_${batch}.err \
../git_repo/utils/modisco_merge_job.sh "$output_dir" "$model_head" "$cluster_key" "$modisco_dir" "$contribs_dir" "$batch"

done 

output_dir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/modisco_merged_${model_head}
cluster_key_path=${output_dir}/cluster_key_combined.txt
mkdir /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/compiled/
compiled_modisco_h5_path=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/compiled/modisco_compiled.h5
compiled_modisco_tsv_path=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/compiled/modisco_compiled.tsv
modisco_merged_dir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/modisco_merged_${model_head}

python ../git_repo/utils/04a-compile_modisco_obj.py \
  --output_dir ${output_dir} \
  --cluster_key_path ${cluster_key_path} \
  --compiled_modisco_h5_path ${compiled_modisco_h5_path} \
  --compiled_modisco_tsv_path ${compiled_modisco_tsv_path} \
  --modisco_merged_dir ${modisco_merged_dir}



# load conda env
conda activate chrombpnet
output_dir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/main/scatac_ArchR/chromBPnet/compiled_${model_head}
mkdir $output_dir

compiled_h5=$output_dir/modisco_compiled.h5
out_dir=$output_dir

echo $compiled_h5
echo $out_dir
motif_db=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/DBs/HOCOMOCO_db/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme

python -u code/03-chrombpnet/02-compendium/04b-get_tomtom_matches.py --modisco-h5 $compiled_h5 \
    --out-dir $out_dir \
    --meme-db $motif_db \
    --verbose True
                    
#echo "done."