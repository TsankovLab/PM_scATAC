

xargs -n 1 curl -O -L < /sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/tumor_compartment/ENCODE_normal/ENCODE_query_fragment_files.txt
for tarfile in *.tar.gz; do
    tar -xzf "$tarfile"
done

#!/bin/bash

#!/bin/bash

# Define the base directory containing the folders
BASE_DIR="encode_scatac_dcc_2/results"

# Define the output directory where compressed files will be stored
OUTPUT_DIR="compressed_fragments"
mkdir -p "$OUTPUT_DIR"

# Find and process each compressed fragment file
find "$BASE_DIR" -type f -name "fragments.tsv.gz" | while read -r fragment_file; do
    # Extract the sample name from the file path (corrected extraction)
    sample_name=$(basename "$(dirname "$(dirname "$fragment_file")")")

    # Temporary file for decompression
    temp_file="${fragment_file%.gz}"

    # Decompress the fragment file
    if gunzip -c "$fragment_file" > "$temp_file"; then
        echo "Decompressed: $fragment_file"
        
        # Recompress using bgzip
        if bgzip -c "$temp_file" > "$OUTPUT_DIR/${sample_name}_fragments.tsv.gz"; then
            echo "Recompressed: $temp_file -> $OUTPUT_DIR/${sample_name}_fragments.tsv.gz"
        else
            echo "Error recompressing: $temp_file" >&2
        fi
        
        # Clean up the temporary decompressed file
        rm -f "$temp_file"
    else
        echo "Error decompressing: $fragment_file" >&2
    fi
done

echo "All fragment files have been processed. Check $OUTPUT_DIR for results."
