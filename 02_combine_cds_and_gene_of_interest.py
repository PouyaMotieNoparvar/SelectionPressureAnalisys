import os

# Define folder paths
folder_genes_to_study = "01_genes_to_study"
folder_orthogroups_cds = "02_orthogroups_cds_unique"
output_folder = "02_orthogroups_cds_plus_genes_of_interest"

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

def read_fasta_as_single_line(filepath):
    """
    Reads a FASTA file and returns its content with each sequence as a single line.
    """
    with open(filepath, 'r') as file:
        lines = file.readlines()
    
    result = []
    current_sequence = []
    
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            # If there's an ongoing sequence, join it and add to the result
            if current_sequence:
                result.append("".join(current_sequence))
                current_sequence = []
            # Add the header line as is
            result.append(line)
        else:
            # Collect sequence lines
            current_sequence.append(line)
    
    # Add the last sequence if there is one
    if current_sequence:
        result.append("".join(current_sequence))
    
    return "\n".join(result)

# Process files
for filename in os.listdir(folder_genes_to_study):
    # Paths for each file
    genes_to_study_path = os.path.join(folder_genes_to_study, filename)
    orthogroups_cds_path = os.path.join(folder_orthogroups_cds, filename)
    output_file_path = os.path.join(output_folder, filename)
    
    # Check if the corresponding file exists in the orthogroups_cds folder
    if os.path.exists(orthogroups_cds_path):
        # Read both FASTA files with sequences in single-line format
        cds_content = read_fasta_as_single_line(orthogroups_cds_path)
        genes_content = read_fasta_as_single_line(genes_to_study_path)
        
        # Combine the contents
        combined_content = cds_content + "\n" + genes_content
        
        # Write the combined content to the output folder
        with open(output_file_path, 'w') as output_file:
            output_file.write(combined_content)
        
        print(f"Processed {filename} and saved to {output_file_path}")
    else:
        print(f"Warning: {filename} not found in {folder_orthogroups_cds}")
