from chimerax.core.commands import run
import csv
import os

base_dir = os.path.dirname(os.path.abspath(__file__))

# folder to store log files
logs_path = os.path.join(base_dir,"logs")
if not os.path.exists(logs_path):
    logs_dir = os.makedirs(logs_path)

# folder to store aligned structures 
aligned_pdb_path = os.path.join(base_dir,"aligned_pdbs")
if not os.path.exists(aligned_pdb_path):
    aligned_pdb_dir = os.mkdir(aligned_pdb_path)

# structral alignment 
def run_alignment(session, target, malaria, count):
    run(session,"log clear")

    run(session, f'open {target}')
    run(session, f"open {malaria}")
    
    #performs structural alignment
    run(session,f"matchmaker #2 to #1")

    # saves the aligned pdb files
    run(session,f'save {aligned_pdb_path}/aligned_{target}_and_{malaria}.pdb')

    run(session, f"log save {logs_path}/output{count}.txt")

    run(session,"close all")

import re
import pandas as pd 
import os

# extracts the rmsd and alignment length
def extract_info_from_log(filepath):
    with open(filepath, 'r', encoding='utf-8') as file:
        html = file.read()

    match = re.findall(
        r'RMSD between (\d+) pruned atom pairs is ([\d.]+) angstroms; \(across all (\d+) pairs: ([\d.]+)\)',
        html
    )
    if match:
        pruned_pairs, rmsd, all_pairs, overall_rmsd = match[-1]
        return {
            "pruned_pairs": int(pruned_pairs),
            "rmsd": float(rmsd),
            "all_pairs": int(all_pairs),
            "overall_rmsd": float(overall_rmsd)
        }
    else:
        print("Could not find  info in the file.")
        return None

df = pd.read_csv('Ultimate_Dataset.csv')


base_dir = os.path.dirname(os.path.abspath(__file__))
folder_path = os.path.join(base_dir, "logs")

# adds rmsd and alignment length to the csv
for i in range(len(df)+1):
    file_path = os.path.join(folder_path, f"output{i+1}.txt")
    print(f"Looking for: {file_path}")
    if os.path.exists(file_path):
        result = extract_info_from_log(file_path)
        if result:
            df.at[i, '3d_rmsd'] = result["rmsd"]
            df.at[i, '3d_align_len'] = result["pruned_pairs"]
    else:
        print(f"File not found: {file_path}")


df.to_csv("Ultimate_with_rmsd.csv", index=False)

def main(session):
    input_csv = os.path.join(base_dir,"merged.csv")

    count = 0

    with open(input_csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            target = row['target_name'] 
            malaria = row['malaria_name'] 
            #target_chain = row['target_chain']
            #malaria_chain = row['malaria_chain']
            
            count += 1
            run_alignment(session, target, malaria,count)

main(session)