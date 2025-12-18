#!/bin/bash
#BSUB -J chrBP_model
#BSUB -P acc_Tsankov_Normal_Lung
#BSUB -q gpu
#BSUB -n 8
#BSUB -W 48:00
#BSUB -gpu num=2
#BSUB -R v100
#BSUB -R rusage[mem=32000]
#BSUB -R span[hosts=1]
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

# ml anaconda3/2022.10
ml cuda/11.7.0
ml cudnn/8.9.5-11
ml proxies
ml java/11.0.2
ml tensorrt/8.5.3.1

#ml anaconda3/2020.11
source activate chrombpnet

chromBPdir=${1}
echo $chromBPdir
grefdir=${2}
echo $grefdir
repodir=${3}
echo $repodir
celltype=${4}
echo $celltype
MODEL_H5=${5}
echo $MODEL_H5
#mkdir $chromBPdir
#chromBPdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet
#celltype=low_P23
#MODEL_H5=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet/${celltype}/no_bias_model/fold_0/models/chrombpnet_nobias.h5
#/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet
grefdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet
motif_file=motif_footprints2.txt
motif_file=motif_footprints_IRF.txt

chromBPdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/myeloid_cells/scatac_ArchR/chromBPnet
#chromBPdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/NKT_cells/scatac_ArchR/chromBPnet
#chromBPdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/scatac_ArchR/chromBPnet
celltype=inflamed
celltype=NK_KLRC1
celltype=SOX9_high_P23
grefdir=/sc/arion/projects/Tsankov_Normal_Lung/Bruno/chromBPnet
OUTPUT_PREFIX=${chromBPdir}/${celltype}/footprints

cd $chromBPdir/$celltype

mkdir footprints
cd footprints
for fold_number in {0..4}; do
    # motif_file=motif_endothelial.TF.txt

    # mkdir -p $chromBPdir/$celltype   # only needed once, so you can keep it outside the loop if you like

    MODEL_H5=${chromBPdir}/${celltype}/no_bias_model/fold_${fold_number}/models/chrombpnet_nobias.h5
    chrombpnet footprints -m $MODEL_H5 \
        -r ../${celltype}_peakset_all_no_blacklist.bed \
        -g $grefdir/genome_references/hg38.genome.fa \
        -fl $grefdir/folds/fold_${fold_number}.json \
        -op ${OUTPUT_PREFIX}/fold${fold_number} \
        -pwm_f ../../${motif_file} # [-bs BATCH_SIZE] [--ylim YLIM]
done

export HDF5_PLUGIN_PATH=$CONDA_PREFIX/lib/hdf5/plugin

export FOOTPRINT_GROUPS=$(awk '{print $1}' ../../${motif_file} \
                         | sort -u \
                         | paste -sd "," -),control

echo "Using groups: $FOOTPRINT_GROUPS"

python3 - <<'EOF'
import h5py
import numpy as np
import glob
import os

# Files to average
files = sorted(glob.glob("fold*_footprints.h5"))
print("Found files:", files)

# Read groups from environment variable
groups = os.environ["FOOTPRINT_GROUPS"].split(",")
print("Using groups:", groups)

# Get datasets from the first file dynamically
with h5py.File(files[0], "r") as f:
    datasets = {grp: list(f[grp].keys()) for grp in groups}

# Output file
with h5py.File("average_footprints.h5", "w") as out:
    for grp in groups:
        print(f"Processing group: {grp}")
        grp_out = out.create_group(grp)
        for ds in datasets[grp]:
            arrays = []
            for fn in files:
                with h5py.File(fn, "r") as f:
                    arr = f[f"{grp}/{ds}"][:]
                    arrays.append(arr)
            mean_arr = np.mean(arrays, axis=0)
            grp_out.create_dataset(ds, data=mean_arr)

print("Saved average_footprints.h5")
EOF

python3 - <<'EOF'
import deepdish as dd
import matplotlib.pyplot as plt
import numpy as np
import os

# Load averaged footprints
footprints = dd.io.load("average_footprints.h5")

# Get groups from environment variable
groups = os.environ["FOOTPRINT_GROUPS"].split(",")
print("Using groups:", groups)

plt.figure(figsize=(6,4))

for grp in groups:
    fp = footprints[grp]['i0']  # assuming 'i0' contains the footprint
    center = fp.shape[0] // 2
    plt.plot(range(200), fp[center-100:center+100], label=grp)

plt.xlabel("200bp around motif")
plt.ylabel("Probability")
plt.xticks(ticks=[0,100,200], labels=[-100,0,100])
plt.legend()
plt.tight_layout()
plt.savefig("average_footprints_plot.pdf")  # <-- PDF instead of PNG
EOF

### Now compare composite NKFBI from inflamed and non-inflamed

#motif_file="motifs.tsv"

# Build FOOTPRINT_GROUPS from first column + control
export FOOTPRINT_GROUPS=$(awk '{print $1}' ../../${motif_file} | paste -sd "," -),control
echo "Using groups: $FOOTPRINT_GROUPS"


python3 - <<'EOF'
import deepdish as dd
import matplotlib.pyplot as plt
import numpy as np
import os

# Paths to average footprint files
file1 = "average_footprints.h5"
file2 = "../../inflamed/footprints/average_footprints.h5"

# Load footprints
footprints1 = dd.io.load(file1)
footprints2 = dd.io.load(file2)

# Get groups from environment variable
groups = os.environ["FOOTPRINT_GROUPS"].split(",")
print("Using groups:", groups)

plt.figure(figsize=(8,6))

# Define a color palette
colors = plt.cm.tab10.colors

for grp_idx, grp in enumerate(groups):
    fp_list = []
    sources = [footprints1, footprints2]
    
    for file_idx, fp_data in enumerate(sources):
        if grp in fp_data:
            fp_list.append(fp_data[grp]['i0'])  # assuming 'i0' contains the footprint
    
    if not fp_list:
        continue
    
    # Plot each instance separately
    center = fp_list[0].shape[0] // 2
    for inst_idx, fp in enumerate(fp_list):
        plt.plot(
            range(200),
            fp[center-100:center+100],
            label=f"{grp}_file{inst_idx+1}",
            color=colors[grp_idx % len(colors)],
            linestyle='-' if inst_idx==1 else '--',
            alpha=0.8
        )

plt.xlabel("200bp around motif", fontsize=12)
plt.ylabel("Probability", fontsize=12)
plt.xticks(ticks=[0,100,200], labels=[-100,0,100])
plt.legend(fontsize=8, ncol=2)
plt.tight_layout()
plt.savefig("average_footprints_combined_plot.pdf")
EOF
