import os
import csv
from chimera import runCommand

# Setup directories
base_dir = os.path.dirname(os.path.abspath(__file__))
malaria_dir = os.path.join(base_dir, "malariapdb_renamed")
target_dir = os.path.join(base_dir, "targetpdb_renamed")
log_dir = os.path.join(base_dir, "logs")
aligned_dir = os.path.join(base_dir, "aligned_pdbs")
os.makedirs(log_dir, exist_ok=True)
os.makedirs(aligned_dir, exist_ok=True)

# Path to the dataset
csv_path = os.path.join(base_dir, "Ultimate_Dataset.csv")

# Load and process CSV
with open(csv_path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for count, row in enumerate(reader, start=1):
        malaria_file = os.path.join(malaria_dir, row["malaria_renamed"])
        target_file = os.path.join(target_dir, row["target_renamed"])

        if not os.path.exists(malaria_file):
            print(f"Missing: {malaria_file}")
            continue
        if not os.path.exists(target_file):
            print(f"Missing: {target_file}")
            continue

        # Prepare names
        target_name = os.path.basename(target_file).replace('.pdb', '')
        malaria_name = os.path.basename(malaria_file).replace('.pdb', '')
        aligned_filename = f"aligned_{count}_{target_name}_vs_{malaria_name}.pdb"
        aligned_path = os.path.join(aligned_dir, aligned_filename)
        log_path = os.path.join(log_dir, f"output{count}.txt")

        try:
            # Run Chimera commands
            runCommand("close all")
            runCommand(f"open {target_file}")
            runCommand(f"open {malaria_file}")
            runCommand("matchmaker #1 #2")
            runCommand(f"write format pdb 0 {aligned_path}")
            runCommand(f"log save {log_path}")
        except Exception as e:
            print(f"Failed on row {count}: {e}")
