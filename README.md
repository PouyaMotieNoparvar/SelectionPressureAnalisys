SelectionPressureAnalysis
This repository provides custom scripts (i.e., pseudopipeline) to analyze selection pressure across multiple genes simultaneously.

Pipeline Overview
The pipeline includes the following steps:

Retrieve CDS Sequences of Orthologous Genes
Extract CDS (coding sequence) data for a set of orthologous genes.

Combine the Fasta Sequences
Merge the individual fasta files of sequences into a single file.

Clean and Remove Duplicates
Clean the sequences and ensure only unique sequences are present in each fasta file.

Align Sequences Using MAFFT
Perform sequence alignment using the MAFFT algorithm.

Convert Fasta to PHYLIP Format
Convert the aligned sequences from Fasta to PHYLIP format (version 3.2).

Apply pal2nal
Use pal2nal to synchronize the coding sequences with the protein sequences.

Generate Phylogenetic Trees Using RAxML
Build phylogenetic trees based on the aligned sequences using RAxML.

Perform CODEML Analysis
Conduct CODEML analysis to estimate selection pressures on the gene sequences.

Requirements
To run these scripts, you will need to install the necessary tools and packages, including:

MAFFT
pal2nal
RAxML
CODEML (from PAML)
Please ensure all required dependencies are installed before running the pipeline.
