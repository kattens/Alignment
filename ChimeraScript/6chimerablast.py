import os
import csv
from chimerax.core.commands import run

def setup_dirs():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    malaria_dir = os.path.join(base_dir, "malariapdb_renamed")
    target_dir = os.path.join(base_dir, "targetpdb_renamed")
    aligned_dir = os.path.join(base_dir, "aligned_pdbs")
    log_dir = os.path.join(base_dir, "logs")

    if not os.path.exists(aligned_dir):
        os.makedirs(aligned_dir)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    return base_dir, malaria_dir, target_dir, aligned_dir, log_dir

def align_and_save(session,target_path, malaria_path, count, aligned_dir, log_dir):
    run(session,"close all")
    run(session,"log clear")

    run(session,f"open {target_path}")
    run(session,f"open {malaria_path}")

    run(session,"matchmaker #2 to #1")  # Align malaria to target

    aligned_filename = f"aligned_{count+1}_{os.path.basename(target_path).replace('.pdb','')}_vs_{os.path.basename(malaria_path).replace('.pdb','')}.pdb"
    aligned_path = os.path.join(aligned_dir, aligned_filename)

    run(session, f'save {aligned_path}')
    log_path = os.path.join(log_dir, f"output{count+1}.txt")
    run(session,f'log save {log_path}')

def main(session):
    base_dir, malaria_dir, target_dir, aligned_dir, log_dir = setup_dirs()
    csv_path = os.path.join(base_dir, "Ultimate_Dataset.csv")

    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for count, row in enumerate(reader, start=1):
            malaria_file = os.path.join(malaria_dir, row["malaria_renamed"])
            target_file = os.path.join(target_dir, row["target_renamed"])

            if not os.path.exists(malaria_file):
                print("Missing malaria file:", malaria_file)
                continue
            if not os.path.exists(target_file):
                print("Missing target file:", target_file)
                continue
            
            try:
                align_and_save(session, target_file, malaria_file, count, aligned_dir, log_dir)
            except Exception as e:
                print(f"Error during alignment #{count+1} ({os.path.basename(target_file)} vs {os.path.basename(malaria_file)}): {e}")

main(session)
