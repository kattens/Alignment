# This is just a simple query to get the pubchem_id and pdb_id as one query
import json
import requests
import os   

def pdb_query(pubchem_id, pdb_id):
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "text_chem",
                            "parameters": {
                                "attribute": "rcsb_chem_comp_related.resource_accession_code",
                                "operator": "exact_match",
                                "negation": False,
                                "value": pubchem_id  # string, not list
                            }
                        },
                        {
                            "type": "terminal",
                            "service": "text_chem",
                            "parameters": {
                                "attribute": "rcsb_chem_comp_related.resource_name",
                                "operator": "exact_match",
                                "value": "PubChem",
                                "negation": False
                            }
                        }
                    ],
                    "label": "nested-attribute"
                },
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": pdb_id  # string, not list
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 25
            },
            "results_content_type": [
                "experimental"
            ],
            "sort": [
                {
                    "sort_by": "score",
                    "direction": "desc"
                }
            ],
            "scoring_strategy": "combined"
        }
    }


def download_pdb(query, file_dir):
    # Step 1: Send query to RCSB -> this work the same as the regex url address modified by previous code
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, headers=headers, data=json.dumps(query))

    if response.status_code != 200:
        print(f"Query failed: {response.status_code}")
        print(response.text)
        return

    # Step 2: Extract PDB IDs
    result = response.json()
    pdb_ids = [entry["identifier"] for entry in result.get("result_set", [])]
    if not pdb_ids:
        print("No PDB entries found.")
        return

    print(f"Found PDB IDs: {pdb_ids}")

    # Step 3: Download PDB files
    os.makedirs(file_dir, exist_ok=True)
    for pdb_id in pdb_ids:
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        pdb_response = requests.get(pdb_url)
        if pdb_response.status_code == 200:
            with open(os.path.join(file_dir, f"{pdb_id}.pdb"), "w") as f:
                f.write(pdb_response.text)
            print(f"Downloaded: {pdb_id}")
        else:
            print(f"Failed to download: {pdb_id}")