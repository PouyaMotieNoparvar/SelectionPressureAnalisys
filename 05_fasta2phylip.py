import os
from Bio import SeqIO

def fasta_to_phylip(fasta_file, phylip_file):
    """
    Convert a single FASTA file to PHYLIP format and save it.
    """
    # Read the FASTA file using SeqIO
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    # Open the output PHYLIP file for writing
    with open(phylip_file, "w") as out_file:
        # Write the header: number of sequences and length of the first sequence
        out_file.write(f"{len(records)} {len(records[0].seq)}\n")
        
        # Write each sequence in the PHYLIP format
        for record in records:
            # Each sequence should be written with a name and sequence data
            out_file.write(f"{record.id.ljust(10)} {str(record.seq)}\n")

def convert_fasta_to_phylip(input_dir, output_dir):
    """
    Convert all FASTA files in the input directory to PHYLIP format and save them to the output directory.
    """
    # Make sure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Iterate through all files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            # Full path to the current FASTA file
            fasta_file = os.path.join(input_dir, filename)
            
            # Output PHYLIP file path (same name but in the output directory)
            phylip_file = os.path.join(output_dir, filename.replace(".fasta", ".phy").replace(".fa", ".phy"))
            
            # Convert the current FASTA file to PHYLIP format
            fasta_to_phylip(fasta_file, phylip_file)
            print(f"Converted {filename} to PHYLIP format and saved as {phylip_file}")

# Define the input and output directories
input_dir = "./05_fasta2phylip/sequence_names"
output_dir = "./05_fasta2phylip/phylip"

# Call the function to perform the conversion
convert_fasta_to_phylip(input_dir, output_dir)
