cd /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scatac_PM/main/scatac_ArchR

python modisco_to_pfm.py \
  -c chrombpnet_models_paths.tsv \
  -o chromBPnet/pfms.txt \
  -p pos_patterns \
  -t 0.3 \
  -ml 3


python modisco_to_pfm.py \
  -c chrombpnet_models_paths.tsv \
  -o chromBPnet/pfms.txt \
  -p neg_patterns \
  -t 0.3 \
  -ml 3
