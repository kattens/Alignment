from chimerax.core.commands import run
import csv
import os

base_dir = os.path.dirname(os.path.abspath(__file__))


# method for structral alignment 
def run_alignment(session, target, malaria, count):
    run(session,"log clear")

    # Loads  PDBs at specified path
    run(session, f'open {target}')
    run(session, f"open {malaria}")
    

    run(session,f"matchmaker #2 to #1")

    aligned_pdb_path = os.path.join(base_dir,"aligned_pdbs")

    # saves the aligned pdb files
    run(session,f'save {aligned_pdb_path}/aligned_{target}_and_{malaria}.pdb')

    logs_path = os.path.join(base_dir,"logs")

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