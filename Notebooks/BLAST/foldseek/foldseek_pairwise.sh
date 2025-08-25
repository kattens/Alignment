#!/bin/bash

# === [INPUT PARAMETERS] ===
PDB1_ID="1HHO"
PDB2_ID="4HHB"

# === [CUSTOM PATHS] ===
INPUT_DIR="/home/k_ensafitakaldani001_umb_edu/BLAST/foldseek/input
OUTPUT_DIR="/home/k_ensafitakaldani001_umb_edu/BLAST/foldseek/output"

# === [SETUP] ===
mkdir -p "$INPUT_DIR" "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit 1

# === [1. Download PDB Files] ===
echo "📥 Downloading PDB files to: $INPUT_DIR"
wget -q -O "${INPUT_DIR}/${PDB1_ID}.pdb" "https://files.rcsb.org/download/${PDB1_ID}.pdb"
wget -q -O "${INPUT_DIR}/${PDB2_ID}.pdb" "https://files.rcsb.org/download/${PDB2_ID}.pdb"

# === [2. Check Downloads] ===
if [[ ! -s "${INPUT_DIR}/${PDB1_ID}.pdb" || ! -s "${INPUT_DIR}/${PDB2_ID}.pdb" ]]; then
    echo "❌ Error: One or both PDB files failed to download or are empty."
    exit 1
fi

# === [3. Create Foldseek Databases] ===
echo "📦 Creating Foldseek databases..."
foldseek createdb "${INPUT_DIR}/${PDB1_ID}.pdb" queryDB
foldseek createdb "${INPUT_DIR}/${PDB2_ID}.pdb" targetDB

# === [4. Run Foldseek Search] ===
echo "🔍 Running Foldseek search..."
foldseek search queryDB targetDB results tmp --threads 4 --max-seqs 1

# === [5. Convert Alignments to TSV] ===
echo "📄 Formatting alignment output..."
foldseek convertalis queryDB targetDB results alignment_raw.txt

# === [6. Add Column Headers ===]
echo -e "query\ttarget\tfident\talnlen\tmismatch\tgapopen\tqstart\tqend\ttstart\ttend\tevalue\tbitscore" > alignment.txt
cat alignment_raw.txt >> alignment.txt

# === [7. Display Summary] ===
echo -e "\n✅ Foldseek alignment complete! Summary:"
column -t alignment.txt

echo -e "\n📁 Annotated result saved to: $OUTPUT_DIR/alignment.txt"
