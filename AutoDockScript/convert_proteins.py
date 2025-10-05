"""
This script takes pdb files of proteins 
as input and converts them into pdbqt format
"""


import os
import subprocess


base_dir = os.path.dirname(os.path.abspath(__file__))


input_folder = os.path.join(base_dir,'malaria')

# Where to save the output PDBQT files
output_folder = os.path.join('malaria_pdbqt')
if not os.path.exists(output_folder):
    os.makedirs(output_folder, exist_ok=True)

log_file_path = os.path.join(base_dir, "failed_conversions_proteins.txt")
failed = []


prepare_script = os.path.join(base_dir,'prepare_receptor4.py')

for filename in os.listdir(input_folder):
    if filename.endswith(".pdb"):
        protein_path = os.path.join(input_folder, filename)
        output_name = os.path.splitext(filename)[0] + ".pdbqt"
        output_path = os.path.join(output_folder, output_name)

        cmd = [
            os.path.join(base_dir, "myenv", "Scripts", "python.exe"), prepare_script,
            "-r", protein_path,
            "-o", output_path,
            "-v"
        ]

        try:
            subprocess.run(cmd, cwd=base_dir, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Failed: {filename}")
            failed.append(filename)

if failed:
    with open(log_file_path, "w") as log:
        log.write("Failed protein conversions:\n")
        for name in failed:
            log.write(name + "\n")
    print(f"\nLogged {len(failed)} failures to {log_file_path}")
else:
    print("\n All ligands processed successfully.")