"""
https://colab.research.google.com/github/kattens/PubChem-Data-Handler/blob/main/Blast_Run.ipynb#scrollTo=fSTLZjIcndWI

"""


import requests
from Bio.Blast import NCBIWWW, NCBIXML

def fetch_uniprot_sequence(accession_id):
    """
    Fetch a sequence from UniProt using the accession ID.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{accession_id}.fasta"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
        else:
            print(f"UniProt query failed for {accession_id} with status code: {response.status_code}")
            return None
    except Exception as e:
        print(f"Error querying UniProt for {accession_id}: {e}")
        return None

def fetch_and_blast_sequence(accession_id, taxonomy, blast_db="nr",
                             blast_type="blastp", expect=0.1, matrix_name="BLOSUM62",
                             alignments=50, hitlist_size=50, filter="F", gapcosts="11 1"):
    """
    Fetch a sequence from UniProt, run BLAST, print and return results as a dictionary.
    """
    sequence_data = fetch_uniprot_sequence(accession_id)
    if not sequence_data:
        print(f"No sequence data found for Accession ID: {accession_id}")
        return None

    entrez_query = f"txid{taxonomy}[ORGN]"
    result_handle = NCBIWWW.qblast(blast_type, blast_db, sequence_data,
                                   expect=expect, matrix_name=matrix_name,
                                   alignments=alignments, hitlist_size=hitlist_size,
                                   filter=filter, gapcosts=gapcosts,
                                   entrez_query=entrez_query)

    blast_record = NCBIXML.read(result_handle)
    blast_results = []
    print(f"Results for Accession ID: {accession_id}")
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            result_detail = {
                "hit_id": alignment.hit_id,
                "hit_def": alignment.hit_def,
                "e_value": hsp.expect,
                "score": hsp.score,
                "query_align": hsp.query[:100],  # First 100 characters for example
                "subject_align": hsp.sbjct[:100]
            }
            blast_results.append(result_detail)

            # Print each result
            print(f"  Hit ID: {alignment.hit_id}")
            print(f"  Hit Description: {alignment.hit_def}")
            print(f"    E-value: {hsp.expect}")
            print(f"    Score: {hsp.score}")
            print(f"    Query Alignment: {hsp.query[:100]}...")
            print(f"    Subject Alignment: {hsp.sbjct[:100]}...")
            print("-" * 80)

    if blast_results:
        return blast_results
    else:
        print(f"No BLAST results found for Accession ID: {accession_id}")
        return None
        if np.linalg.norm(fc[i] - mc[i]) < distance_cutoff:
                fixed_c.append(fc[i])
                moving_c.append(mc[i])
        else:
            print(f"Pruning residue {i} due to distance > {distance_cutoff} Å") 
            