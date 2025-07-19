import os
import pandas as pd
from Bio.PDB import PDBParser
from Preprocessing.Processing import ParsePDB
from Align.superimpose import iterative_alignment, get_aligned_coords, extract_sequence_and_coords

def run_preprocessing_pipeline(target_pdb, malaria_pdb, target_chain_id, malaria_chain_id, drug_resname, input_dir, output_dir):
    # Full paths
    target_pdb_path = os.path.join(input_dir, target_pdb)
    malaria_pdb_path = os.path.join(input_dir, malaria_pdb)

    drug_resname = row['drug_info'].split(':')[0].strip()
    
    # Output files-> prefixes are for making the output unique
    prefix = os.path.splitext(target_pdb)[0] + "_" + os.path.splitext(malaria_pdb)[0]
    out_chain_only = os.path.join(output_dir, f"{prefix}_target_chain_only.pdb")
    out_ligand = os.path.join(output_dir, f"{prefix}_target_chain_ligand.pdb")
    out_target = os.path.join(output_dir, f"{prefix}_target_renamed.pdb")
    out_malaria = os.path.join(output_dir, f"{prefix}_malaria_renamed.pdb")

    #debug ^_^
    print(out_chain_only,out_ligand,out_malaria,out_target)
    print(target_pdb_path,malaria_pdb_path)


    

    # Preprocessing
    print(f"\nProcessing: {target_pdb} vs {malaria_pdb}")
    parser = ParsePDB(target_pdb_path, target_chain_id, drug_resname)
    parser.keep_target(out_chain_only)
    parser.keep_drug(out_chain_only, out_ligand)
    parser.rename_target_ligand(out_ligand, target_chain_id, out_target)
    parser.keep_malaria(malaria_pdb_path, malaria_chain_id, out_malaria)

    # Alignment: rename chains as expected
    target_chain_id = "T"
    malaria_chain_id = "M"

    parser = PDBParser(QUIET=True)
    structure1 = parser.get_structure("struct1", out_target)
    structure2 = parser.get_structure("struct2", out_malaria)

    seq1, coords1, _ = extract_sequence_and_coords(structure1[0][target_chain_id])
    seq2, coords2, _ = extract_sequence_and_coords(structure2[0][malaria_chain_id])
    fixed_coords, moving_coords = get_aligned_coords(seq1, coords1, seq2, coords2)
    aligned_fixed, aligned_moving, rmsd = iterative_alignment(fixed_coords, moving_coords)
    
    


def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    input_dir = os.path.join(base_dir, "PDBs")
    output_dir = os.path.join(base_dir, "OutputPDBs")
    os.makedirs(output_dir, exist_ok=True)

    csv_path = os.path.join(base_dir, "pairs.csv")
    df = pd.read_csv(csv_path)

    for idx, row in df.iterrows():
       
        
        run_preprocessing_pipeline(
            target_pdb=row['target_pdb'].str[:4],
            malaria_pdb=row['malaria_pdb'].str[:4],
            drug_resname=row['ligand'],
            target_chain_id=row['target_chain'],
            malaria_chain_id=row['malaria_chain'],
            input_dir=input_dir,
            output_dir=output_dir
        )
        

if __name__ == "__main__":
    main()
