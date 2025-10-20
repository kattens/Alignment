# pdb_utils.py

from __future__ import annotations
from typing import List, Dict, Tuple, Any, Optional, Iterable
from collections import defaultdict
import os
import re
import csv
import math
import time
import requests
import numpy as np
import pandas as pd

from Bio.PDB import (
    PDBParser, PDBIO, Select, PPBuilder, Selection, NeighborSearch, is_aa, Superimposer
)
from Bio.Blast import NCBIWWW, NCBIXML

# =============================================================================
#                               CORE PDB HELPERS
# =============================================================================

def get_structure(pdb_file: str):
    """Load a PDB structure."""
    parser = PDBParser(QUIET=True, PERMISSIVE=True)
    return parser.get_structure("struct", pdb_file)

def parse_all_entities(pdb_file: str) -> Dict[str, Dict[str, Any]]:
    """
    For each chain: protein AA sequence + categorized residues (amino_acids/ligands/waters/others).
    """
    structure = get_structure(pdb_file)
    ppb = PPBuilder()
    results: Dict[str, Dict[str, Any]] = {}

    for model in structure:
        for chain in model:
            chain_id = chain.id
            sequence = ''.join(str(pp.get_sequence()) for pp in ppb.build_peptides(chain))

            amino_acids, ligands, waters, others = [], [], [], []
            for residue in chain:
                hetfield, resseq, icode = residue.id
                resname = residue.get_resname()
                atoms = [{'atom_name': a.get_name(), 'element': a.element} for a in residue]
                info = {'residue_name': resname, 'resseq': resseq, 'insertion_code': icode, 'atom_sequence': atoms}

                if hetfield == ' ' and is_aa(residue, standard=True):
                    amino_acids.append(info)
                elif resname == 'HOH':
                    waters.append(info)
                elif hetfield != ' ':
                    ligands.append(info)
                else:
                    others.append(info)

            results[chain_id] = {
                'protein_sequence': sequence,
                'amino_acids': amino_acids,
                'ligands': ligands,
                'waters': waters,
                'others': others
            }
    return results

def get_atoms_from_chain(structure, chain_id: str, atom_names: Optional[Iterable[str]] = None,
                         aa_only: bool = True, exclude_h: bool = True):
    """Return atoms from a specific chain, optionally filtered to AA and atom names (e.g., ['CA'])."""
    atoms = []
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for residue in chain:
                if aa_only and not is_aa(residue):
                    continue
                for atom in residue:
                    if exclude_h and atom.element == 'H':
                        continue
                    if atom_names is None or atom.get_name() in atom_names:
                        atoms.append(atom)
    return atoms

def get_ligand_atoms(structure, ligand_resname: str, ligand_chain_id: Optional[str] = None,
                     exclude_h: bool = True):
    """Return atoms of residues matching ligand_resname (optionally only from a given chain)."""
    atoms = []
    for model in structure:
        for chain in model:
            if ligand_chain_id and chain.id != ligand_chain_id:
                continue
            for residue in chain:
                hetfield, _, _ = residue.id
                if hetfield != ' ' and residue.get_resname().strip() == ligand_resname:
                    for atom in residue:
                        if exclude_h and atom.element == 'H':
                            continue
                        atoms.append(atom)
    return atoms

# ---------------------------- #
# ---------- Filters --------- #
# ---------------------------- #

def filter_pdb_by_predicate(input_pdb: str, output_pdb: str, keep_line_fn):
    """
    Generic PDB filter. keep_line_fn(line) -> bool decides whether to keep the line.
    Writes only ATOM/HETATM that satisfy predicate; copies other records unchanged.
    """
    with open(input_pdb, "r") as fi, open(output_pdb, "w") as fo:
        for line in fi:
            if line.startswith(("ATOM", "HETATM")):
                if keep_line_fn(line):
                    fo.write(line)
            else:
                fo.write(line)

def keep_chain_and_ligand_predicate(protein_chain="M", ligand_chain="T", ligand_resname="D"):
    """
    Keep ATOMs from protein_chain and HETATMs for ligand_resname on ligand_chain.
    """
    def _pred(line: str) -> bool:
        rec = line[0:6].strip()
        chain = line[21].strip() if len(line) >= 22 else ""
        resn = line[17:20].strip() if len(line) >= 20 else ""
        if rec == "ATOM" and chain == protein_chain:
            return True
        if rec == "HETATM" and chain == ligand_chain and resn == ligand_resname:
            return True
        return False
    return _pred

def filter_chain_and_ligand(input_pdb: str, output_pdb: str,
                            protein_chain="M", ligand_chain="T", ligand_resname="D"):
    """Convenience wrapper: keep chain M and ligand D on chain T (defaults)."""
    pred = keep_chain_and_ligand_predicate(protein_chain, ligand_chain, ligand_resname)
    filter_pdb_by_predicate(input_pdb, output_pdb, pred)

def filter_chain_excluding(input_pdb: str, output_pdb: str, exclude_chain: str = "T"):
    """Remove ATOM lines from exclude_chain; keep everything else (incl. HETATM)."""
    def _pred(line: str) -> bool:
        rec = line[0:6].strip()
        chain = line[21].strip() if len(line) >= 22 else ""
        if rec == "ATOM" and chain == exclude_chain:
            return False
        return True
    filter_pdb_by_predicate(input_pdb, output_pdb, _pred)

def batch_filter(input_dir: str, output_dir: str, predicate_factory, *factory_args, **factory_kwargs):
    """Apply a predicate filter to all .pdb files in a folder."""
    os.makedirs(output_dir, exist_ok=True)
    for fname in os.listdir(input_dir):
        if not fname.lower().endswith(".pdb"):
            continue
        inp = os.path.join(input_dir, fname)
        out = os.path.join(output_dir, fname)
        pred = predicate_factory(*factory_args, **factory_kwargs)
        filter_pdb_by_predicate(inp, out, pred)

# ---------------------------- #
# ---------- RMSD ------------ #
# ---------------------------- #

def compute_rmsd_atoms(atoms1, atoms2) -> float:
    """RMSD with Bio.PDB.Superimposer; lists must be equal length."""
    if len(atoms1) != len(atoms2):
        raise ValueError(f"Atom lists have different lengths: {len(atoms1)} vs {len(atoms2)}")
    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    return sup.rms

def rmsd_between_ligand_and_chain(pdb_file: str, protein_chain_id: str = "M",
                                  ligand_resname: str = "D", ligand_chain_id: Optional[str] = None) -> Optional[float]:
    """
    Truncates to min length if atom counts differ, then computes RMSD.
    NOTE: heuristic (ligand vs protein atoms) — interpret cautiously.
    """
    structure = get_structure(pdb_file)
    lig_atoms = get_ligand_atoms(structure, ligand_resname, ligand_chain_id)
    prot_atoms = get_atoms_from_chain(structure, protein_chain_id, atom_names=None, aa_only=True, exclude_h=True)
    if not lig_atoms or not prot_atoms:
        return None
    n = min(len(lig_atoms), len(prot_atoms))
    return compute_rmsd_atoms(lig_atoms[:n], prot_atoms[:n])

# ---------------------------- #
# ---- Binding site search ----#
# ---------------------------- #

def find_ligand_binding_sites(pdb_file: str, distance_cutoff: float = 5.0) -> List[Dict[str, Any]]:
    """
    Find residues within distance_cutoff of each non-water hetero residue (ligand).
    Returns list of dicts per ligand with binding residue (id, name) sets.
    """
    structure = get_structure(pdb_file)
    model = next(structure.get_models())
    atoms = Selection.unfold_entities(model, 'A')
    ns = NeighborSearch(atoms)

    results: List[Dict[str, Any]] = []
    for chain in model:
        chain_id = chain.id
        for residue in chain:
            hetfield, resseq, icode = residue.get_id()
            if hetfield == ' ' or residue.get_resname() == 'HOH':
                continue
            ligand_atoms = list(residue.get_atoms())
            binding_residues = set()
            for atom in ligand_atoms:
                for neighbor in ns.search(atom.get_coord(), distance_cutoff, level='R'):
                    if neighbor.get_parent().id != chain_id:
                        continue
                    if neighbor.get_id()[0] == ' ':
                        binding_residues.add((neighbor.get_id()[1], neighbor.get_resname()))
            results.append({
                'chain_id': chain_id,
                'ligand_name': residue.get_resname(),
                'ligand_resseq': resseq,
                'binding_site_residues': sorted(binding_residues)
            })
    return results

# ---------------------------- #
# ------ Ligand–CA stats ------#
# ---------------------------- #

def ligand_ca_contacts(pdb_path: str, protein_chain_id: str = 'M', ligand_resname: str = 'D',
                       distance_threshold: float = 5.0) -> List[Dict[str, Any]]:
    """
    Compute contacts between all ligand atoms (any chain) and CA atoms in protein_chain_id.
    """
    structure = get_structure(pdb_path)
    lig_atoms = get_ligand_atoms(structure, ligand_resname, ligand_chain_id=None, exclude_h=False)
    ca_atoms = get_atoms_from_chain(structure, protein_chain_id, atom_names=["CA"], aa_only=True, exclude_h=False)

    contacts: List[Dict[str, Any]] = []
    for l in lig_atoms:
        lxyz = l.get_coord()
        lres = l.get_parent()
        for ca in ca_atoms:
            cxyz = ca.get_coord()
            dist = np.linalg.norm(lxyz - cxyz)
            if dist <= distance_threshold:
                cres = ca.get_parent()
                contacts.append({
                    'ligand_residue': lres.get_resname(),
                    'ligand_resseq': lres.id[1],
                    'ligand_atom': l.get_name(),
                    'protein_residue': cres.get_resname(),
                    'protein_resseq': cres.id[1],
                    'ca_atom': ca.get_name(),
                    'distance': float(round(dist, 3))
                })
    return contacts

def save_contacts_to_csv(contacts: List[Dict[str, Any]], output_csv: str):
    """Write contact dicts to CSV."""
    if not contacts:
        pd.DataFrame(columns=[
            "ligand_residue","ligand_resseq","ligand_atom","protein_residue",
            "protein_resseq","ca_atom","distance"
        ]).to_csv(output_csv, index=False)
        return
    pd.DataFrame(contacts).to_csv(output_csv, index=False)

def process_contacts_folder(input_dir: str, output_dir: str, protein_chain_id: str = 'M',
                            ligand_resname: str = 'D', distance_threshold: float = 5.0):
    """Run ligand_ca_contacts on all .pdb files in a folder and save CSVs (same basename)."""
    os.makedirs(output_dir, exist_ok=True)
    for fname in os.listdir(input_dir):
        if not fname.lower().endswith(".pdb"):
            continue
        pdb_path = os.path.join(input_dir, fname)
        out_csv = os.path.join(output_dir, os.path.splitext(fname)[0] + ".csv")
        contacts = ligand_ca_contacts(pdb_path, protein_chain_id, ligand_resname, distance_threshold)
        save_contacts_to_csv(contacts, out_csv)

# =============================================================================
#                         RCSB / UniProt / AlphaFold HELPERS
# =============================================================================

def fetch_json(url: str) -> Dict[str, Any]:
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.json()

def get_structure_summary_df(pdb_id: str) -> pd.DataFrame:
    """
    Return summary DataFrame for a PDB entry
    (molecules/chains/gene names/organisms/lengths/mutations/ligands).
    """
    pdb_id = pdb_id[:4].upper()
    entry = fetch_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
    polymer_ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", []) or []
    ligand_ids = entry.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_ids", []) or []

    molecule_list, chains_list, gene_names_list, organisms_list, lengths_list, mutations_list = [], [], [], [], [], []

    for ent_id in polymer_ids:
        poly = fetch_json(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{ent_id}")
        molecule = poly.get("entity", {}).get("pdbx_description", "N/A")
        chains = poly.get("rcsb_polymer_entity", {}).get("pdbx_strand_id", []) or []
        orgs = [o.get("ncbi_scientific_name", "N/A") for o in poly.get("rcsb_entity_source_organism", [])]
        gene_names = []
        if poly.get("rcsb_entity_source_organism"):
            gene_names = [g.get("value") for g in poly["rcsb_entity_source_organism"][0].get("rcsb_gene_name", [])]
        length = poly.get("entity_poly", {}).get("rcsb_sample_sequence_length", "N/A")
        mutations = poly.get("entity_poly", {}).get("rcsb_mutation_count", "N/A")

        molecule_list.append(molecule)
        chains_list.append(", ".join(chains))
        gene_names_list.append(", ".join(gene_names) if gene_names else "N/A")
        organisms_list.append(", ".join(orgs))
        lengths_list.append(length)
        mutations_list.append(mutations)

    ligand_names = []
    for ent_id in ligand_ids:
        lig = fetch_json(f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{ent_id}")
        nonpoly = lig.get("pdbx_entity_nonpoly", {})
        name = nonpoly.get("name", "N/A")
        comp_id = nonpoly.get("comp_id", "N/A")
        ligand_names.append(f"{comp_id}: {name}")

    return pd.DataFrame([{
        "PDB_ID": pdb_id,
        "Molecules": molecule_list,
        "Chains": chains_list,
        "Gene_Names": gene_names_list,
        "Organisms": organisms_list,
        "Sequence_Lengths": lengths_list,
        "Mutations": mutations_list,
        "Ligands": ligand_names
    }])

def annotate_df_with_rcsb(df: pd.DataFrame, pdb_col: str, out_gene_col: str, out_lig_col: str) -> pd.DataFrame:
    """
    For a dataframe with PDB IDs (possibly with chain suffix), add columns with Gene_Names and Ligands lists.
    """
    ids = df[pdb_col].astype(str).str[:4].str.upper().unique().tolist()
    info = []
    for pid in ids:
        try:
            info.append(get_structure_summary_df(pid))
        except Exception:
            info.append(pd.DataFrame([{
                "PDB_ID": pid,
                "Molecules": [], "Chains": [], "Gene_Names": ["N/A"],
                "Organisms": [], "Sequence_Lengths": [], "Mutations": [], "Ligands": ["N/A"]
            }]))
        # polite pacing for APIs
        time.sleep(0.1)

    combined = pd.concat(info, ignore_index=True)
    lookup = {r["PDB_ID"]: {"Gene_Names": r["Gene_Names"], "Ligands": r["Ligands"]} for _, r in combined.iterrows()}
    df[out_gene_col] = df[pdb_col].astype(str).str[:4].str.upper().map(lambda x: lookup.get(x, {}).get("Gene_Names", ["N/A"]))
    df[out_lig_col]  = df[pdb_col].astype(str).str[:4].str.upper().map(lambda x: lookup.get(x, {}).get("Ligands", ["N/A"]))
    return df

def query_rcsb_by_gene_name(gene_name: str, page_size: int = 100, timeout: int = 30) -> List[str]:
    """
    Text query RCSB for entries that match a gene name exactly (experimental entries).
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    payload = {
        "query": {
            "type": "terminal",
            "label": "text",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                "operator": "exact_match",
                "value": gene_name
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": page_size},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
            "scoring_strategy": "combined"
        }
    }
    r = requests.post(url, json=payload, timeout=timeout)
    if r.status_code != 200:
        return []
    return [e["identifier"] for e in r.json().get("result_set", [])]

def get_uniprot_id(raw_query: str, timeout: int = 30) -> Optional[str]:
    """Resolve a query (name/accession-like) to a single UniProt primary accession."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": raw_query, "format": "json", "fields": "accession", "size": 1}
    r = requests.get(url, params=params, timeout=timeout)
    if r.status_code != 200:
        return None
    results = r.json().get("results", [])
    return results[0]["primaryAccession"] if results else None

def get_alphafold_filename(uniprot_id: Optional[str]) -> Optional[str]:
    """Return expected AF2 v4 filename for a UniProt ID."""
    return f"AF-{uniprot_id}-F1-model_v4.pdb" if uniprot_id else None

def annotate_structure_info(df: pd.DataFrame, pdb_col_name: str, prefix: str) -> pd.DataFrame:
    """
    Convenience wrapper around annotate_df_with_rcsb using name + prefix convention.
    Adds: f'{prefix}_gene_names', f'{prefix}_ligands'
    """
    out_gene = f"{prefix}_gene_names"
    out_lig  = f"{prefix}_ligands"
    return annotate_df_with_rcsb(df.copy(), pdb_col=pdb_col_name, out_gene_col=out_gene, out_lig_col=out_lig)

# ---------------------------- #
# --------- Downloads -------- #
# ---------------------------- #

def download_pdb(pdb_id: str, save_dir: str) -> Optional[str]:
    """Download a PDB file from RCSB if missing. Returns path or None on failure."""
    os.makedirs(save_dir, exist_ok=True)
    pid = pdb_id.strip()[:4].upper()
    url = f'https://files.rcsb.org/download/{pid}.pdb'
    out = os.path.join(save_dir, f"{pid}.pdb")
    if os.path.exists(out):
        return out
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        with open(out, "w") as f:
            f.write(r.text)
        return out
    except Exception:
        return None

# ---------------------------- #
# --- File/CSV small utils --- #
# ---------------------------- #

name_re = re.compile(
    r"^(?P<target>[0-9A-Za-z]{4})_(?P<tchain>[A-Za-z0-9])_"
    r"(?P<malaria>[0-9A-Za-z]{4})_(?P<mchain>[A-Za-z0-9])_"
    r"(?P<ligand>[A-Za-z0-9]+)$"
)

def parse_filename(fname_no_ext: str) -> Optional[Tuple[str, str, str]]:
    """Parse 'target_targetChain_malaria_malariaChain_ligand' -> (target, malaria, ligand)."""
    m = name_re.match(fname_no_ext)
    if not m:
        return None
    return m.group("target"), m.group("malaria"), m.group("ligand")

def safe_unique_ligand_atom_count(df_csv: pd.DataFrame) -> int:
    """Try several reasonable columns to estimate unique ligand atom count in a contacts CSV."""
    if "ligand_coord" in df_csv.columns:
        return df_csv["ligand_coord"].nunique(dropna=True)
    if {"ligand_atom", "ligand_resseq"}.issubset(df_csv.columns):
        return (df_csv["ligand_atom"].astype(str) + ":" + df_csv["ligand_resseq"].astype(str)).nunique(dropna=True)
    if "ligand_atom" in df_csv.columns:
        return df_csv["ligand_atom"].astype(str).nunique(dropna=True)
    return len(df_csv.drop_duplicates())

def count_ligand_atoms_in_pdb(pdb_path: str, residue_name: str = "D"):
    """Count number of HETATM lines for a given residue_name (or return codes if not found/error)."""
    try:
        atom_count, found = 0, False
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("HETATM") and line[17:20].strip() == residue_name:
                    found = True
                    atom_count += 1
        return atom_count if found else "RESIDUE_NOT_FOUND"
    except Exception as e:
        return f"ERROR: {e}"

# =============================================================================
#                            UniProt + ONLINE BLAST
# =============================================================================

def fetch_uniprot_sequence(accession_id: str, timeout: int = 30) -> Optional[str]:
    """
    Fetch the FASTA for a UniProt accession ID. Returns FASTA text (str) or None.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{accession_id}.fasta"
    try:
        r = requests.get(url, timeout=timeout)
        if r.status_code == 200 and r.text.strip():
            return r.text
        return None
    except Exception:
        return None

def fetch_and_blast_sequence(
    accession_id: str,
    taxonomy: int,
    blast_db: str = "nr",
    blast_type: str = "blastp",
    expect: float = 0.1,
    matrix_name: str = "BLOSUM62",
    alignments: int = 50,
    hitlist_size: int = 50,
    filter: str = "F",
    gapcosts: str = "11 1",
    verbose: bool = False,
    sleep_seconds: float = 0.0
) -> Optional[List[Dict[str, Any]]]:
    """
    Fetch a sequence from UniProt and run NCBI online BLAST via Biopython.
    Returns a list of result dicts (id/def/e_value/score/query_align/subject_align) or None.

    NOTE: NCBI imposes usage limits—use sparingly and consider local BLAST for batches.
    """
    fasta_text = fetch_uniprot_sequence(accession_id)
    if not fasta_text:
        if verbose:
            print(f"[BLAST] No FASTA for UniProt: {accession_id}")
        return None

    entrez_query = f"txid{taxonomy}[ORGN]"
    if sleep_seconds:
        time.sleep(sleep_seconds)  # be gentle to NCBI

    handle = NCBIWWW.qblast(
        program=blast_type,
        database=blast_db,
        sequence=fasta_text,
        expect=expect,
        matrix_name=matrix_name,
        alignments=alignments,
        hitlist_size=hitlist_size,
        filter=filter,
        gapcosts=gapcosts,
        entrez_query=entrez_query,
    )
    record = NCBIXML.read(handle)

    out: List[Dict[str, Any]] = []
    for aln in record.alignments:
        for hsp in aln.hsps:
            out.append({
                "hit_id": aln.hit_id,
                "hit_def": aln.hit_def,
                "e_value": hsp.expect,
                "score": hsp.score,
                "query_align": getattr(hsp, "query", "")[:200],
                "subject_align": getattr(hsp, "sbjct", "")[:200],
            })
    if verbose:
        if out:
            print(f"[BLAST] {accession_id}: {len(out)} HSPs")
        else:
            print(f"[BLAST] {accession_id}: no hits")
    return out if out else None

# =============================================================================
#                       PIPELINED DATAFRAME OPERATIONS
# =============================================================================

def prepare_genename_table(df_gene_names: pd.DataFrame) -> pd.DataFrame:
    """
    Clean and augment a PlasmoDB gene table:
      - Extract a normalized 'Accession' (best-effort)
      - UniProt_ID (via search)
      - AlphaFold_File (expected AF2 file name)
      - Gene_PDB (comma-joined RCSB entries by gene name)
    """
    df = df_gene_names.copy()
    # Best-effort extraction: keep prefix up to a trailing digit block if present
    df["Accession"] = df["Accession"].astype(str).str.extract(r'^(.*?0?\d*)', expand=False)

    df["UniProt_ID"] = df["Accession"].apply(lambda q: get_uniprot_id(q) if pd.notna(q) and q.strip() else None)
    df["AlphaFold_File"] = df["UniProt_ID"].apply(get_alphafold_filename)

    def _rcsb(g: str) -> str:
        if pd.isna(g) or not str(g).strip():
            return "No Match"
        ids = query_rcsb_by_gene_name(str(g).strip())
        return ", ".join(ids) if ids else "No Match"

    df["Gene_PDB"] = df["Accession"].apply(_rcsb)
    return df

def annotate_targets_and_malaria(
    df_targets: pd.DataFrame,
    target_pdb_col: str = "target_chain_id",
    malaria_pdb_col: str = "malaria_match_id"
) -> pd.DataFrame:
    """
    Add RCSB annotations to a target/malaria mapping table.
    Produces new columns:
      - target_gene_names, target_ligands
      - malaria_gene_names, malaria_ligands
    """
    df = df_targets.copy()
    df = annotate_structure_info(df, pdb_col_name=target_pdb_col, prefix="target")
    df = annotate_structure_info(df, pdb_col_name=malaria_pdb_col, prefix="malaria")
    return df
