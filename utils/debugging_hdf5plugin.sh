#!/bin/bash
#ml anaconda3/2020.11
conda env list
source activate h5py

# Debugging: Output environment variables
echo "Environment variables:"
env | grep -i conda

# Debugging: Verify Python path
python -c "import sys; print(sys.executable)"

# Debugging: Test importing hdf5plugin
python -c "import hdf5plugin; print('hdf5plugin imported successfully')"

# Run your script
python /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo/utils/average_CNT_scores.py
