import h5py
import numpy as np
import sys
import os

# Get variables from sys.argv
chromBPct_dir = sys.argv[1]  # First argument
#var2 = sys.argv[2]  # Second argument
#var3 = sys.argv[3]  # Third argument

os.chdir (chromBPct_dir)

# Define file paths
print ('average contribution count scores')
h5_files = ["fold_0/C2_test_contribution_scores.counts_scores.h5", 
"fold_1/C2_test_contribution_scores.counts_scores.h5",
"fold_2/C2_test_contribution_scores.counts_scores.h5",
"fold_3/C2_test_contribution_scores.counts_scores.h5",
"fold_4/C2_test_contribution_scores.counts_scores.h5"]  # List of H5 file paths

# Function to compute the average for a specific key
def average_h5_keys(files, key):
    total = None  # Initialize for accumulation
    count = 0     # Counter for the number of files
    
    for file in files:

        if not os.path.exists(file):
            print(f"Warning: File {file} does not exist.")
            continue
            
        with h5py.File(file, 'r') as h5_file:
            # Read the dataset
            data = h5_file[key][:]
            
            # Accumulate the sum
            if total is None:
                total = np.zeros_like(data)
            total += data
            count += 1
    
    if count > 0:
        return total / count
    else:
        raise ValueError("No valid files found to process.")

# Keys to average
keys = ["projected_shap/seq", "raw/seq", "shap/seq"]

# Create a new H5 file to store the averaged data
output_file = "averaged_count_scores.h5"
with h5py.File(output_file, 'w') as h5_output:
    for key in keys:
        # Compute the average for the current key
        averaged_data = average_h5_keys(h5_files, key)
        
        # Write the averaged data to the new file
        h5_output.create_dataset(key, data=averaged_data)

# Verify the output file structure
with h5py.File(output_file, 'r') as h5_output:
    print("Averaged H5 File Structure:")
    h5_output.visititems(lambda name, obj: print(f"{name} ({'Group' if isinstance(obj, h5py.Group) else 'Dataset'})"))


# Define file paths
print ('average contribution profile scores')
h5_files = ["fold_0/C2_test_contribution_scores.profile_scores.h5", 
"fold_1/C2_test_contribution_scores.profile_scores.h5",
"fold_2/C2_test_contribution_scores.profile_scores.h5",
"fold_3/C2_test_contribution_scores.profile_scores.h5",
"fold_4/C2_test_contribution_scores.profile_scores.h5"]  # List of H5 file paths

# Function to compute the average for a specific key
def average_h5_keys(files, key):
    total = None  # Initialize for accumulation
    count = 0     # Counter for the number of files
    
    for file in files:
        with h5py.File(file, 'r') as h5_file:
            # Read the dataset
            data = h5_file[key][:]
            
            # Accumulate the sum
            if total is None:
                total = np.zeros_like(data)
            total += data
            count += 1
    
    # Compute the average
    return total / count

# Keys to average
keys = ["projected_shap/seq", "raw/seq", "shap/seq"]

# Create a new H5 file to store the averaged data
output_file = "averaged_profile_scores.h5"
with h5py.File(output_file, 'w') as h5_output:
    for key in keys:
        # Compute the average for the current key
        averaged_data = average_h5_keys(h5_files, key)
        
        # Write the averaged data to the new file
        h5_output.create_dataset(key, data=averaged_data)

# Verify the output file structure
with h5py.File(output_file, 'r') as h5_output:
    print("Averaged H5 File Structure:")
    h5_output.visititems(lambda name, obj: print(f"{name} ({'Group' if isinstance(obj, h5py.Group) else 'Dataset'})"))


