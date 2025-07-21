from chimerax.core.commands import run
import csv

count=0

# method for structral alignment 
def run_alignment(session,target,malaria,target_chain,malaria_chain):
    run(session,"log clear")

    # Loads  PDBs at specified path
    run(session, f'open C:/malaria_align_project/pdbs/{target}')
    run(session, f"open C:/malaria_align_project/pdbs/{malaria}")
    
    # alignment by specific chain
    run(session, f"matchmaker #2/{malaria_chain} to #1/{target_chain}")
    
    run(session, f"log save C:/malaria_align_project/logs/output{count}.txt")

    run(session,"close all")



with open("C:/malaria_align_project/malaria_mapping.csv", newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        # grabs the name of the pdb
        target = row['target']
        malaria = row['malaria']

        # grabs the chains that are being focused on
        target_chain = row['target_chain']
        malaria_chain = row['malaria_chain']

        count+=1

        run_alignment(session, f'{target}.pdb', f'{malaria}.pdb',target_chain,malaria_chain)

