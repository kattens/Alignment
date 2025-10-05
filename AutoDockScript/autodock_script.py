import os
import subprocess
import pandas as pd

from Bio.PDB import PDBParser
import numpy as np


vina_executable = r"C:\Program Files (x86)\PyRx\vina.exe"  # AutoDock Vina executable path

base_dir = os.path.dirname(os.path.abspath(__file__))


output_folder = os.path.join(base_dir,"output")
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

log_folder = os.path.join(base_dir,"docking_logs")
if not os.path.exists(log_folder):
    os.mkdir(log_folder)

csv_folder = os.path.join(base_dir,"top10_affinities")
if not os.path.exists(csv_folder):
    os.mkdir(csv_folder)

lig_folder = os.path.join(base_dir,'ligands_pdbqt')
prot_folder = os.path.join(base_dir,'malaria_pdbqt')

csv_path = os.path.join('docking_dataset.csv')


def estimate_docking_box(pdb_file, buffer=5.0):
    """
    Estimates the center and size of the docking box from a PDBQT file.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", pdb_file)
    
    coords = np.array([atom.coord for atom in structure.get_atoms() if atom.element != "H"])
    
    # Calculate bounding box
    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    center = coords.mean(axis=0)
    size = max_coords - min_coords + buffer
    
    return tuple(center.tolist()), tuple(size.tolist())



def run_vina(size,center,log_folder,prefix,receptor,ligand):
    cmd = [
        vina_executable,
        "--receptor", receptor,
        "--ligand", ligand,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--out", os.path.join(output_folder,prefix+'.pdbqt'),
        "--log", os.path.join(log_folder,prefix)
    ]
    print("Running AutoDock Vina...")
    subprocess.run(cmd, check=True)



def parse_affinities(log_folder,prefix):
    affinities = []
    log_file = os.path.join(log_folder,prefix)
    with open(log_file, "r") as f:
        for line in f:
            if line.strip().startswith(tuple(str(i) for i in range(1,11))):
                
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        affinities.append(float(parts[1]))
                    except ValueError:
                        continue
    return affinities[:10]  # Top 10


def save_to_csv(affinities,prefix,csv_folder):
    df = pd.DataFrame({"Rank": range(1, len(affinities)+1),
                       "Binding_Affinity_kcal/mol": affinities})
    csv_file = os.path.join(csv_folder,prefix)
    df.to_csv(csv_file, index=False)
    print(f"Top 10 affinities saved to {csv_file}.csv")






def run_docking(lig_folder, prot_folder, log_folder,csv_folder):
    df = pd.read_csv(csv_path)
    prefix_col = "dock_helper"
    failed = []
    for prefix in df[prefix_col]:
        lig_path = os.path.join(lig_folder, f"{prefix}_ligand.pdbqt")
        prot_path = os.path.join(prot_folder, f"{prefix}_malaria.pdbqt")

        if os.path.exists(lig_path) and os.path.exists(prot_path):
            print(f"Found pair: {prefix}")

            try:
                center, size = estimate_docking_box(prot_path)
                size = (25.0,25.0,25.0) # we will keep each dimension at 25 angstroms for now

                print("Docking box center:", center, "size:", size)
                run_vina(size, center,log_folder,prefix,prot_path,lig_path)
                affinities = parse_affinities(log_folder,prefix)
                if affinities:
                    save_to_csv(affinities,prefix,csv_folder)
                else:
                    print("No binding affinities found. Check docking log.")
            except Exception as e:
                print(f"Error: {e}")


        else:
            print(f"Missing file(s) for: {prefix}")
            failed.append(prefix)

    failed_file_path = os.path.join(base_dir, "failed_dockings.txt")


    if failed:
        with open(failed_file_path, "w") as log:
            log.write("File not found:\n")
            for name in failed:
                log.write(name + "\n")
        print(f"\nLogged {len(failed)} failures to {failed_file_path}")
    else:
        print("\nAll dockings done successfully.")



if __name__ == "__main__":
    run_docking(lig_folder,prot_folder,log_folder,csv_folder)