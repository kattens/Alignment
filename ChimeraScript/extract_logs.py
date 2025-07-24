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

df = pd.read_csv('merged(3).csv')

# create new columns
df['rmsd'] = None
df['alignment_length'] = None

base_dir = os.path.dirname(os.path.abspath(__file__))
folder_path = os.path.join(base_dir, "logs")
print(folder_path)

for i in range(len(df)):
    file_path = os.path.join(folder_path, f"output{i+1}.txt")
    print(f"Looking for: {file_path}")
    if os.path.exists(file_path):
        result = extract_info_from_log(file_path)
        if result:
            df.at[i, 'rmsd'] = result["rmsd"]
            df.at[i, 'alignment_length'] = result["pruned_pairs"]
    else:
        print(f"File not found: {file_path}")


df.to_csv("merged_with_rmsd.csv", index=False)
