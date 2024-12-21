import h5py
import numpy as np

# Define file paths
h5_files = ["file1.h5", "file2.h5", "file3.h5"]  # List of H5 file paths

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
output_file = "averaged_data.h5"
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
