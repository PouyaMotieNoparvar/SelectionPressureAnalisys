import os

def process_fasta(file_path, log_file):
    """Process a single FASTA file to make sequence titles unique and log the changes."""
    with open(file_path, 'r') as file:
        lines = file.readlines()

    titles = {}
    new_lines = []
    changes = []

    for line in lines:
        if line.startswith(">"):
            if line in titles:
                # Increment the counter for repeated titles
                titles[line] += 1
                unique_line = line.strip() + f"_{titles[line]}\n"
                changes.append(f"{line.strip()} -> {unique_line.strip()} in {os.path.basename(file_path)}")
            else:
                # Add the title to the dictionary with an initial count
                titles[line] = 1
                unique_line = line
            new_lines.append(unique_line)
        else:
            # Non-title lines remain unchanged
            new_lines.append(line)

    # Overwrite the file with the updated lines
    with open(file_path, 'w') as file:
        file.writelines(new_lines)

    # Log changes
    with open(log_file, 'a') as log:
        for change in changes:
            log.write(change + "\n")

def main():
    # Directory containing the FASTA files
    directory = "/mnt/e/Positive_Selection/New_Pipeline/04_align_mafft/cds_alignment"  # Update with your WSL path
    log_file = os.path.join(directory, "log.txt")

    # Clear the log file at the beginning
    with open(log_file, 'w') as log:
        log.write("Sequence Title Changes Log\n")
        log.write("=" * 30 + "\n")

    fasta_files = [
        "AT1G02780.fasta", "AT1G74230.fasta", "AT2G20820.fasta", "AT2G27720.fasta", "AT2G41840.fasta",
        "AT3G21420.fasta", "AT3G30720.fasta", "AT3G52590.fasta", "AT3G53730.fasta", "AT3G60245.fasta",
        "AT3G62250.fasta", "AT4G25740.fasta", "AT4G38680.fasta", "AT5G27850.fasta", "AT5G65360.fasta",
        "AT5G67600.fasta", "AT3G12145.fasta"
    ]

    for fasta_file in fasta_files:
        file_path = os.path.join(directory, fasta_file)
        if os.path.isfile(file_path):
            print(f"Processing {fasta_file}...")
            process_fasta(file_path, log_file)
        else:
            print(f"File {fasta_file} not found in the directory.")

if __name__ == "__main__":
    main()
