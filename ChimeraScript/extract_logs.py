import re
import pandas as pd 
import os

# extracts the rmsd and pairs for the html
def extract_info_from_log(filepath):
    with open(filepath, 'r', encoding='utf-8') as file:
        html = file.read()

    # Regular expression to extract the line with RMSD
    match = re.findall(
        r'RMSD between (\d+) pruned atom pairs is ([\d.]+) angstroms; \(across all (\d+) pairs: ([\d.]+)\)',
        html
    )
    if match:
        # Get the last RMSD reported in the log
        pruned_pairs, rmsd, all_pairs, overall_rmsd = match
        return {
            "pruned_pairs": int(pruned_pairs),
            "rmsd": float(rmsd),
            "all_pairs": int(all_pairs),
            "overall_rmsd": float(overall_rmsd)
        }
    else:
        print("Could not find  info in the file.")
        return None

df = pd.read_csv('malaria_mapping.csv')

# create new columns
df['rmsd'] = None
df['pruned_pairs']= None
df['all_pairs'] = None

folder_path = "C:/malaria_align_project/logs"

for i in range(len(df)):
    file_path = os.path.join(folder_path, f"output{i+1}.txt")
    print(f"Looking for: {file_path}")
    if os.path.exists(file_path):
        result = extract_info_from_log(file_path)
        if result:
            df.at[i, 'rmsd'] = result["rmsd"]
            df.at[i, 'pruned_pairs'] = result["pruned_pairs"]
            df.at[i, 'all_pairs'] = result["all_pairs"]
    else:
        print(f"File not found: {file_path}")


df.to_csv("Malaria_Mapping_with_SMILES_&_sequences_columns_seperated.csv")