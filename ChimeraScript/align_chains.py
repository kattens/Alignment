from chimerax.core.commands import run
import csv
import os
count=0

base_dir = os.path.dirname(os.path.abspath(__file__))
print(base_dir)


# method for structral alignment 
def run_alignment(session, target, malaria, target_chain, malaria_chain, count):
    run(session,"log clear")

    # Loads  PDBs at specified path
    run(session, f'open {base_dir}/pdbs/{target}')
    run(session, f"open {base_dir}/pdbs/{malaria}")
    
    # alignment by specific chain
    run(session, f"matchmaker #2/{malaria_chain} to #1/{target_chain}")
    
    run(session, f"log save {base_dir}/logs/output{count}.txt")

    run(session,"close all")


def main(session):
    input_csv = f"{base_dir}/malaria_mapping.csv"
    count = 0

    with open(input_csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            target = row['target'] + ".pdb"
            malaria = row['malaria'] + ".pdb"
            target_chain = row['target_chain']
            malaria_chain = row['malaria_chain']
            
            count += 1
            run_alignment(session, target, malaria, target_chain, malaria_chain, count)
