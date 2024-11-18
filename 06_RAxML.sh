#!/bin/bash

# Directory containing .phy files
input_dir="./05_fasta2phylip/phylip" # Replace this with the actual directory

# Loop through all .phy files in the directory
for file in "$input_dir"/*.phy; do
    # Extract the base name of the file (e.g., AT1G01100 from AT1G01100.phy)
    base_name=$(basename "$file" .phy)
    
    # Run the RAxML command with the current file
    raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 500 -s "$file" -n "$base_name.tree"
done