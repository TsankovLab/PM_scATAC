conda activate chrombpnet

### Create background regions file - Go in the output chromBPnet folder ####
cd /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet

# Remove regions overlapping black listed regions
bedtools slop -i /sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet/blacklist.bed.gz -g /sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet/hg38.chrom.sizes -b 1057 > temp.bed
bedtools intersect -v -a peakset_all.bed -b temp.bed  > peakset_all_no_blacklist.bed
wc -l peakset_all_no_blacklist.bed # # Make sure number of peaks is not more than 250K

# Generate training validation and test chromosome sets
head -n 23  ../../../../../chromBPnet/hg38.chrom.sizes >  hg38.chrom.subset.sizes

mkdir splits
chrombpnet prep splits -c hg38.chrom.subset.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op splits/fold_0

rm -r output_auxiliary
chrombpnet prep nonpeaks -g /sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet/genome_references/hg38.genome.fa -p peakset_all_no_blacklist.bed -c  /sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet/hg38.chrom.sizes -fl splits/fold_0.json -br /sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet/blacklist.bed.gz -o output


### Copy bias model from TNK compartment ####
cp -r /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet/bias_model .




