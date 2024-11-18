#!/bin/bash

# Path to pal2nal.pl script (inside the current directory)
PAL2NAL_SCRIPT="./pal2nal.v14/pal2nal.pl"

# Define directories for protein and CDS alignments
PROTEIN_DIR="./04_align_mafft/protein_alignment"
CDS_DIR="./02_orthogroups_cds_plus_genes_of_interest"

# Output directory for results
OUTPUT_DIR="./05_pal2nal_results"
mkdir -p "$OUTPUT_DIR"

# Loop over all protein alignment files in the PROTEIN_DIR
for PROTEIN_FILE in "$PROTEIN_DIR"/*.fasta; do
  # Extract the base filename (e.g., 1234.fasta -> 1234)
  BASENAME=$(basename "$PROTEIN_FILE" .fasta)

  # Define corresponding CDS file
  CDS_FILE="$CDS_DIR/$BASENAME.fasta"

  # Define output file path
  OUTPUT_FILE="$OUTPUT_DIR/$BASENAME.fasta"

  # Check if the corresponding CDS file exists
  if [ -f "$CDS_FILE" ]; then
    # Run pal2nal.pl with the current pair of files
    perl "$PAL2NAL_SCRIPT" "$PROTEIN_FILE" "$CDS_FILE" -output fasta > "$OUTPUT_FILE"
    echo "Processed $BASENAME: Output saved to $OUTPUT_FILE"
  else
    echo "Warning: CDS file for $BASENAME not found, skipping..."
  fi
done
