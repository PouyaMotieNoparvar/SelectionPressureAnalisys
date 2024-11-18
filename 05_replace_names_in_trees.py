import os
import re

# Define the paths to the files and directories
txt_file_path = '05_fasta2phylip/sequence_names/sequence_name_log.txt'
newick_dir = '06_tree/02_plain_trees'
output_dir = '06_tree/03_truncated_names'  # New directory to save modified files

# Ensure the output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Step 1: Read the sequence_name_log.txt and create a mapping dictionary
name_mapping = {}

with open(txt_file_path, 'r') as f:
    for line in f:
        # Skip the header line
        if line.startswith("Source File"):
            continue
        # Split the line into columns
        columns = line.strip().split('\t')
        original_name = columns[1].strip()
        new_name = columns[2].strip()

        # Remove the '>' from the original name if present (for mapping)
        original_name = original_name.lstrip('>')  # Strip any leading '>' characters

        # Store the mapping in the dictionary
        name_mapping[original_name] = new_name

# Step 2: Iterate through the Newick files in the specified directory
newick_files = [f for f in os.listdir(newick_dir) if f.endswith('.newick')]

# Process each Newick file
for newick_file in newick_files:
    newick_file_path = os.path.join(newick_dir, newick_file)

    # Read the contents of the Newick file
    with open(newick_file_path, 'r') as f:
        tree_data = f.read()

    # Step 3: Replace the node names in the tree
    # We will use regular expressions to replace each original name with the new name
    for original_name, new_name in name_mapping.items():
        # Only replace if the original name exists in the tree (ignoring '>' sign)
        tree_data = re.sub(rf'\b{re.escape(original_name)}\b', new_name, tree_data)

    # Step 4: Save the modified tree to the new directory
    output_file_path = os.path.join(output_dir, newick_file)

    with open(output_file_path, 'w') as f:
        f.write(tree_data)

    print(f"Updated tree in file: {newick_file}")

print("All Newick files processed and saved to:", output_dir)
