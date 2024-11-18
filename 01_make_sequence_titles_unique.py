import os
from collections import defaultdict

# Directories
input_dir = "02_orthogroups_cds"
output_dir = "02_orthogroups_cds_unique"
log_file = os.path.join(output_dir, "sequence_rename_log.txt")

# Ensure the input directory exists
if not os.path.exists(input_dir):
    raise FileNotFoundError(f"The directory {input_dir} does not exist.")

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Prepare log content
log_entries = []

# Process each FASTA file in the directory
for filename in os.listdir(input_dir):
    if filename.endswith(".fasta") or filename.endswith(".fa"):
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)
        unique_names = defaultdict(int)
        output_lines = []

        with open(input_path, "r") as file:
            for line in file:
                if line.startswith(">"):
                    original_name = line.strip()
                    # Generate a unique name if necessary
                    unique_name = original_name
                    if unique_names[original_name] > 0:
                        unique_name = f"{original_name}_{unique_names[original_name]}"
                    unique_names[original_name] += 1

                    # Log the change
                    if unique_name != original_name:
                        log_entries.append(f"{filename}\t{original_name}\t{unique_name}")

                    output_lines.append(unique_name + "\n")
                else:
                    output_lines.append(line)

        # Write to the new FASTA file
        with open(output_path, "w") as file:
            file.writelines(output_lines)

# Write the log file
with open(log_file, "w") as log:
    log.write("Source File\tOriginal Name\tUnique Name\n")
    log.write("\n".join(log_entries))

print(f"Sequence names updated. Files saved in '{output_dir}' and log file created: {log_file}")
