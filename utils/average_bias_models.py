import h5py
import numpy as np
import sys
import os
import hdf5plugin

# Get variables from sys.argv
chromBPct_dir = sys.argv[1]  # First argument
celltype = sys.argv[2]
#var2 = sys.argv[2]  # Second argument
#var3 = sys.argv[3]  # Third argument

os.chdir (chromBPct_dir)

# Define file paths
# print ('average contribution count scores')


h5_counts_files = [
    f"bias_model/fold_0/models/bias.h5",
    f"bias_model/fold_1/modles/bias.h5",
    f"bias_model/fold_2/modles/bias.h5",
    f"bias_model/fold_3/modles/bias.h5",
    f"bias_model/fold_4/modles/bias.h5"
]  # List of H5 file paths

# List of input H5 files
#h5_files = ["file1.h5", "file2.h5", "file3.h5"]  # Replace with your file paths
output_file = "averaged_bias_model.h5"

# Open all files in parallel and compute the average
with h5py.File(h5_counts_files[0], 'r') as ref_file:  # Use the first file as a reference for structure
    with h5py.File(output_file, 'w') as out_file:
        for group_key in ref_file.keys():  # Iterate through groups
            ref_group = ref_file[group_key]
            out_group = out_file.create_group(group_key)  # Create group in the output file
            
            for subkey in ref_group.keys():  # Iterate through datasets in the group
                # Load all corresponding datasets from all files
                data_stack = np.array([h5py.File(f, 'r')[group_key][subkey][:] for f in h5_counts_files])
                
                # Compute the average across the stacked datasets
                averaged_data = np.mean(data_stack, axis=0)
                
                # Save the averaged data into the output file
                out_group.create_dataset(subkey, data=averaged_data)

print(f"Averaged data saved in: {output_file}")




# Function to compute averages for all keys across files
# def average_h5_keys_across_files(files):
#     averages = {}
#     subkey_sums = {}
#     subkey_counts = {}
    
#     # Iterate over all files
#     for file in files:
#         with h5py.File(file, 'r') as h5_file:
#             # Iterate over all groups in the file
#             for group_key in h5_file.keys():
#                 group = h5_file[group_key]
                
#                 # Iterate over all subkeys (datasets) within the group
#                 for subkey in group.keys():
#                     data = group[subkey][:]
                    
#                     # Accumulate the sum and count for each subkey
#                     if subkey not in subkey_sums:
#                         subkey_sums[subkey] = np.zeros_like(data, dtype=np.float64)
#                         subkey_counts[subkey] = 0
                    
#                     subkey_sums[subkey] += data
#                     subkey_counts[subkey] += 1
    
#     # Compute the averages for each subkey
#     for subkey in subkey_sums:
#         if subkey_counts[subkey] == 0:
#             raise ValueError(f"No data found for subkey '{subkey}'")
#         averages[subkey] = subkey_sums[subkey] / subkey_counts[subkey]
    
#     return averages

# # Define the list of H5 files


# # Create a new H5 file to store the averaged data
# output_file = "averaged_contributions_counts.h5"
# with h5py.File(output_file, 'w') as h5_output:
#     try:
#         # Compute averages for all keys across files
#         averages = average_h5_keys_across_files(h5_counts_files)
        
#         # Write the averaged data to the new file, preserving group and key structure
#         for subkey, avg_data in averages.items():
#             print(f"Processing group: {group_key}, subkey: {subkey}")  # Debug statement
#             # Determine which group the key belongs to by inspecting the first file
#             with h5py.File(h5_counts_files[0], 'r') as h5_file:
#                 for group_key in h5_file.keys():
#                     group = h5_file[group_key]
#                     if subkey in group:
#                         # If the subkey exists in the group, write it to the new file under the same structure
#                         if group_key not in h5_output:
#                             h5_output.create_group(group_key)  # Create the group if it doesn't exist
#                         h5_output.create_dataset(f"{group_key}/{subkey}", data=avg_data)
#                         break
        
#         print(f"Averaged data saved in: {output_file}")
#     except ValueError as e:
#         print("Error:", e)


# with h5py.File(output_file, 'r') as h5_output:
#     def print_structure(name, obj):
#         print(f"{name} ({'Group' if isinstance(obj, h5py.Group) else 'Dataset'})")
#     h5_output.visititems(print_structure)


# # Create a new H5 file to store the averaged data
# output_file = "averaged_contributions_profiles.h5"
# with h5py.File(output_file, 'w') as h5_output:
#     try:
#         # Compute averages for all keys across files
#         averages = average_h5_keys_across_files(h5_files_profiles_counts)
        
#         # Write the averaged data to the new file, preserving group and key structure
#         for subkey, avg_data in averages.items():
#             # Determine which group the key belongs to by inspecting the first file
#             with h5py.File(h5_files_profiles_counts[0], 'r') as h5_file:
#                 for group_key in h5_file.keys():
#                     group = h5_file[group_key]
#                     if subkey in group:
#                         # If the subkey exists in the group, write it to the new file under the same structure
#                         if group_key not in h5_output:
#                             h5_output.create_group(group_key)  # Create the group if it doesn't exist
#                         h5_output.create_dataset(f"{group_key}/{subkey}", data=avg_data)
#                         break
        
#         print(f"Averaged data saved in: {output_file}")
#     except ValueError as e:
#         print("Error:", e)


# # group_key = ["projected_shap","raw","shap"]

# # try:
# #     averages = average_h5_group(h5_files, group_key)
# #     for subkey, avg in averages.items():
# #         print(f"Average for subkey '{subkey}': {avg}")
# # except ValueError as e:
# #     print("Error:", e)

# # def list_h5_keys_recursive(h5_obj, prefix=''):
# #     for key in h5_obj.keys():
# #         full_key = f"{prefix}/{key}" if prefix else key
# #         print(full_key)
# #         if isinstance(h5_obj[key], h5py.Group):  # If the key is a group, recurse
# #             list_h5_keys_recursive(h5_obj[key], full_key)
# # with h5py.File(h5_files[4], 'r') as h5_file:
# #     list_h5_keys_recursive(h5_file)


# # # import h5py

# # # def inspect_h5_file(file_path, key):
# # #     try:
# # #         with h5py.File(file_path, 'r') as h5_file:
# # #             # Check if the key exists
# # #             if key in h5_file:
# # #                 item = h5_file[key]
# # #                 print(f"Key '{key}' found.")
# # #                 print(f"Type: {type(item)}")
# # #                 if isinstance(item, h5py.Dataset):
# # #                     print(f"Shape: {item.shape}")
# # #                     print(f"Dtype: {item.dtype}")
# # #                 elif isinstance(item, h5py.Group):
# # #                     print(f"'{key}' is a group. Available subkeys: {list(item.keys())}")
# # #                 else:
# # #                     print(f"Unknown object type: {type(item)}")
# # #             else:
# # #                 print(f"Key '{key}' not found in file.")
# # #     except Exception as e:
# # #         print(f"Error inspecting file: {e}")

# # # # Example usage
# # # inspect_h5_file('fold_4/C2_test_contribution_scores.counts_scores.h5', 'shap')

# # # Keys to average
# # keys = ["projected_shap", "raw", "shap"]
# # keys = ["shap"]

# # # Create a new H5 file to store the averaged data
# # output_file = "averaged_count_scores.h5"
# # with h5py.File(output_file, 'w') as h5_output:
# #     for key in keys:
# #         # Compute the average for the current key
# #         averaged_data = average_h5_keys(h5_files, key)
        
# #         # Write the averaged data to the new file
# #         h5_output.create_dataset(key, data=averaged_data)

# # # Verify the output file structure
# # with h5py.File(output_file, 'r') as h5_output:
# #     print("Averaged H5 File Structure:")
# #     h5_output.visititems(lambda name, obj: print(f"{name} ({'Group' if isinstance(obj, h5py.Group) else 'Dataset'})"))


# # # Define file paths
# # print ('average contribution profile scores')
# # h5_files = ["fold_0/C2_test_contribution_scores.profile_scores.h5", 
# # "fold_1/C2_test_contribution_scores.profile_scores.h5",
# # "fold_2/C2_test_contribution_scores.profile_scores.h5",
# # "fold_3/C2_test_contribution_scores.profile_scores.h5",
# # "fold_4/C2_test_contribution_scores.profile_scores.h5"]  # List of H5 file paths

# # # Function to compute the average for a specific key
# # def average_h5_keys(files, key):
# #     total = None  # Initialize for accumulation
# #     count = 0     # Counter for the number of files
    
# #     for file in files:
# #         with h5py.File(file, 'r') as h5_file:
# #             # Read the dataset
# #             data = h5_file[key][:]
            
# #             # Accumulate the sum
# #             if total is None:
# #                 total = np.zeros_like(data)
# #             total += data
# #             count += 1
    
# #     # Compute the average
# #     return total / count

# # # Keys to average
# # keys = ["projected_shap/seq", "raw/seq", "shap/seq"]

# # # Create a new H5 file to store the averaged data
# # output_file = "averaged_profile_scores.h5"
# # with h5py.File(output_file, 'w') as h5_output:
# #     for key in keys:
# #         # Compute the average for the current key
# #         averaged_data = average_h5_keys(h5_files, key)
        
# #         # Write the averaged data to the new file
# #         h5_output.create_dataset(key, data=averaged_data)

# # # Verify the output file structure
# # with h5py.File(output_file, 'r') as h5_output:
# #     print("Averaged H5 File Structure:")
# #     h5_output.visititems(lambda name, obj: print(f"{name} ({'Group' if isinstance(obj, h5py.Group) else 'Dataset'})"))


