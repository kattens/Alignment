#!/bin/bash

# CONFIGURATION
FASTA_DIR="PATH TO FASTA FILES"       
OUTPUT_DIR="PATH TO THE RESULT CSV SAVE FOLDER"

#MODIFY THE HYPERPARAMETERS AS NEEDED
EVALUE="1e-5"
MAX_HITS=10
TAXID="5833"

# Load BLAST module
module load blast-plus/2.14.1

#PICK DATABASE TYPE, NR, PDB, SWISSPROT,...
NR_DB=$(blastdb_path -db nr -dbtype prot)

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Loop over all FASTA files
for fasta_file in "$FASTA_DIR"/*.fasta; do
    base_name=$(basename "$fasta_file" .fasta)
    output_file="$OUTPUT_DIR/${base_name}_nr_hits.csv" #CAN CHANGE THE NAME IF NEEDED

    echo "Running BLASTP for Protein Sequnece: $base_name"

    blastp \
        -query "$fasta_file" \
        -db "$NR_DB" \
        -taxids "$TAXID" \
        -evalue "$EVALUE" \
        -max_target_seqs "$MAX_HITS" \
        -outfmt "10 qseqid sseqid stitle evalue bitscore pident length" \
        -out "$output_file"

    echo "Saved results to: $output_file"
done

echo "All FASTA files processed."
