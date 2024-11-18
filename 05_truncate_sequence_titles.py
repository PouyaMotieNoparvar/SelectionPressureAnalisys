import os
import random
import string

# Define paths
input_folder = "./05_pal2nal_results"  # Folder containing input FASTA files
output_folder = "./05_fasta2phylip/sequence_names"  # Folder for output files
log_file_path = os.path.join(output_folder, "sequence_name_log.txt")  # Path for log file

# Ensure the output directory exists
os.makedirs(output_folder, exist_ok=True)

def generate_unique_id(existing_ids, length=6):
    """Generate a unique ID of specified length."""
    while True:
        new_id = ''.join(random.choices(string.ascii_lowercase + string.digits, k=length))
        if new_id not in existing_ids:
            existing_ids.add(new_id)
            return new_id

def process_fasta_files(input_folder, output_folder, log_file_path):
    """Process FASTA files, truncate titles, and generate a log file."""
    log_entries = []
    
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".fasta"):
            input_file_path = os.path.join(input_folder, file_name)
            output_file_path = os.path.join(output_folder, file_name)
            
            existing_ids = set()
            with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
                for line in infile:
                    if line.startswith(">"):  # Title line
                        original_name = line.strip()[1:]  # Remove '>' for comparison
                        if original_name == "A_thaliana":
                            # Keep the title unchanged for A_thaliana
                            new_name = original_name
                        else:
                            new_name = generate_unique_id(existing_ids)
                        
                        log_entries.append(f"{file_name}\t>{original_name}\t>{new_name}")
                        outfile.write(f">{new_name}\n")  # Write title line
                    else:
                        outfile.write(line)  # Write sequence line
    
    # Write the log file
    with open(log_file_path, 'w') as log_file:
        log_file.write("Source File\tOriginal Name\tNew Name\n")
        log_file.write("\n".join(log_entries))

# Execute the script
if __name__ == "__main__":
    process_fasta_files(input_folder, output_folder, log_file_path)
    print(f"Processing complete. Log file saved at: {log_file_path}")
