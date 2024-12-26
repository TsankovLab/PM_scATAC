ml anaconda3/2020.11
source activate chrombpnet

chromBPdir=${1}
echo $chromBPdir
grefdir=${2}
echo $grefdir
repodir=${3}
echo $repodir
celltype=${4}
echo $celltype

### Create background regions file - Go in the output chromBPnet folder ####
cd $chromBPdir

# Remove regions overlapping black listed regions
bedtools slop -i ${grefdir}/blacklist.bed.gz -g ${grefdir}/hg38.chrom.sizes -b 1057 > temp.bed
bedtools intersect -v -a peakset_all.bed -b temp.bed  > peakset_all_no_blacklist.bed
wc -l peakset_all_no_blacklist.bed # # Make sure number of peaks is not more than 250K

# Generate training validation and test chromosome sets
head -n 23  ${grefdir}/hg38.chrom.sizes >  hg38.chrom.subset.sizes

# mkdir splits
# chrombpnet prep splits -c hg38.chrom.subset.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op splits/fold_0

# Train chrombpnet bias model
BIAS_MODEL_H5=${celltype}_bias_model/models/model_bias.h5

if [ ! -f "${BIAS_MODEL_H5}" ]; then
    echo "chrombpnet_nobias.h5 file not found. Training chromBPnet bias model..."
rm -r output_auxiliary
chrombpnet prep nonpeaks \
	-g ${grefdir}/genome_references/hg38.genome.fa -p \
	peakset_all_no_blacklist.bed \
	-c ${grefdir}/hg38.chrom.sizes \
	-fl splits/fold_${fold_number}.json \
	-br ${grefdir}/blacklist.bed.gz \
	-o output_negatives_${fold_number}

rm -r ${celltype}_model/fold_$fold_number/

# Train bias model
chrombpnet bias pipeline \
        -ifrag fragments_${celltype}.tsv \
        -d "ATAC" \
        -g ${grefdir}/genome_references/hg38.genome.fa \
        -c ${grefdir}/hg38.chrom.sizes \
        -p peakset_all_no_blacklist.bed \
        -n output_negatives_${fold_number}.bed \
        -fl ${grefdir}/splits/fold_${foldr_number}.json \
        -b 0.5 \
        -o ${celltype}_bias_model/ \
        -fp $celltype

else
    echo "model_bias.h5 file found!"
fi
