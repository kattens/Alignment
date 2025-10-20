#DOWNLOAD CODE FOR ANY WEBSITE AND THE REQUEST SENDING
#CAN USE THIS TO DOWNLOAD ALPHAOFLD OR FASTA OR PDB OR CSV FILES AS LONG AS WE HAVE A LINK
import numpy as np
import pandas as pd
import os
import requests


#SET THE PATHS
save_dir = 'path to the save folder'
os.makedirs(save_dir, exist_ok=True)
log_file = os.path.join(save_dir, 'failed_downloads.log')

#KEEP TRACK OF THE ENTRIES
downloaded = []
failed = []

#CALL THE DATAFRAME
df = pd.read_csv("path to the csv file")


for i in df['column name']:
    fasta_name = f"{i}"
    fasta_path = os.path.join(save_dir, fasta_name)

    # Skip if already downloaded
    if os.path.exists(fasta_path):
        print(f"Already downloaded: {fasta_name}")
        downloaded.append(i)
        continue

    # Build UniProt query URL
    query_url = f'the url from the download link/{i}.type of the file like .pdb or .fasta or...'

    try:
        response = requests.get(query_url, timeout=10)
        if response.status_code == 200:
            with open(fasta_path, 'w') as f:
                f.write(response.text)
            print(f"Downloaded: {fasta_name}")
            downloaded.append(i)
        else:
            raise Exception("Invalid or empty response")

    except Exception as e:
        print(f"Failed: {i} ({e})")
        failed.append(i)
        with open(log_file, 'a') as log:
            log.write(f"{i}\t{e}\n")

print(f"\nFinished. {len(downloaded)} downloaded, {len(failed)} failed.")
