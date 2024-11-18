import os
import logging

def reformat_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sequence = ""
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    outfile.write(sequence + '\n')
                    sequence = ""
                outfile.write(line + '\n')
            else:
                sequence += line
        if sequence:
            outfile.write(sequence + '\n')

def reformat_fasta_files(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for filename in os.listdir(input_dir):
        if filename.endswith('.fasta'):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, filename.replace('.fasta', '.fasta'))
            
            if os.path.exists(output_file):
                base, ext = os.path.splitext(output_file)
                output_file = base + ext  # Example: Rename the file if it already exists
            
            try:
                reformat_fasta(input_file, output_file)
                logging.info(f"Processed {filename}")
            except Exception as e:
                logging.error(f"Error processing {filename}: {str(e)}")

# Configure logging
logging.basicConfig(level=logging.INFO)

# Define your input and output directories
input_dir = "./05_pal2nal_results"
output_dir = "./working_directory"

# Run the function to process all files
reformat_fasta_files(input_dir, output_dir)

