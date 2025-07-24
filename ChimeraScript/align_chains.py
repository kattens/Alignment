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