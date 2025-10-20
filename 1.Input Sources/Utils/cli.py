#!/usr/bin/env python3
"""
cli.py — Command-line interface for protein/ligand utilities powered by Biopython + RCSB/UniProt/NCBI.

Requires:
  - Python 3.9+
  - biopython, pandas, numpy, requests

Install:
  pip install biopython pandas numpy requests

Usage examples:

  # Parse entities and print JSON to stdout
  python cli.py parse-entities --pdb path/to/file.pdb

  # Filter to chain M and ligand D on chain T
  python cli.py filter --in input.pdb --out output.pdb --protein-chain M --ligand-chain T --ligand-resname D

  # Compute ligand–CA contacts and write CSV
  python cli.py contacts --pdb input.pdb --out-csv contacts.csv --protein-chain M --ligand-resname D --cutoff 5.0

  # Batch contacts for a folder of PDB files
  python cli.py contacts-batch --in-dir ./pdbs --out-dir ./contacts --protein-chain M --ligand-resname D --cutoff 5.0

  # RMSD (heuristic) between ligand D and protein chain M
  python cli.py rmsd --pdb input.pdb --protein-chain M --ligand-resname D

  # Ligand binding sites within cutoff
  python cli.py binding-sites --pdb input.pdb --cutoff 5.0 --out-json sites.json

  # Download a PDB by ID to a folder
  python cli.py download-pdb --pdb-id 1ABC --out-dir ./pdbs

  # Annotate a CSV with RCSB gene/ligand columns from a PDB ID column
  python cli.py annotate-rcsb --csv in.csv --pdb-col pdb --out out.csv --gene-col gene_names --lig-col ligands

  # Prepare/augment a PlasmoDB gene names CSV (UniProt/AF filename/Gene->PDB hits)
  python cli.py prepare-genenames --in gene_names.csv --out genenames_aug.csv

  # Annotate targets/malaria CSV with RCSB (adds *_gene_names, *_ligands)
  python cli.py annotate-targets --in targets.csv --out targets_aug.csv --target-col target_chain_id --malaria-col malaria_match_id

  # Online BLAST of a UniProt accession, restrict by taxonomy (e.g. 5833 for Plasmodium)
  python cli.py blast --accession Q8I7G5 --taxonomy 5833 --out-json blast.json --verbose
"""

import os
import sys
import json
import argparse
from typing import Any, Dict, List, Optional

# Ensure local import (pdb_utils.py in the same folder)
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

try:
    import pdb_utils as U
except ImportError as e:
    print("ERROR: Could not import pdb_utils.py. Make sure cli.py is in the same directory as pdb_utils.py.", file=sys.stderr)
    raise

import pandas as pd


def _print_or_save_json(data: Any, out_json: Optional[str]):
    if out_json:
        with open(out_json, "w") as f:
            json.dump(data, f, indent=2, default=str)
        print(f"Saved JSON -> {out_json}")
    else:
        print(json.dumps(data, indent=2, default=str))


def cmd_parse_entities(args: argparse.Namespace):
    res = U.parse_all_entities(args.pdb)
    _print_or_save_json(res, args.out_json)


def cmd_filter(args: argparse.Namespace):
    U.filter_chain_and_ligand(
        input_pdb=args.input,
        output_pdb=args.output,
        protein_chain=args.protein_chain,
        ligand_chain=args.ligand_chain,
        ligand_resname=args.ligand_resname
    )
    print(f"Filtered -> {args.output}")


def cmd_filter_exclude(args: argparse.Namespace):
    U.filter_chain_excluding(args.input, args.output, exclude_chain=args.exclude_chain)
    print(f"Filtered (exclude chain {args.exclude_chain}) -> {args.output}")


def cmd_batch_filter(args: argparse.Namespace):
    factory = U.keep_chain_and_ligand_predicate
    U.batch_filter(
        input_dir=args.in_dir,
        output_dir=args.out_dir,
        predicate_factory=factory,
        protein_chain=args.protein_chain,
        ligand_chain=args.ligand_chain,
        ligand_resname=args.ligand_resname
    )
    print(f"Batch filter complete -> {args.out_dir}")


def cmd_contacts(args: argparse.Namespace):
    contacts = U.ligand_ca_contacts(
        pdb_path=args.pdb,
        protein_chain_id=args.protein_chain,
        ligand_resname=args.ligand_resname,
        distance_threshold=args.cutoff
    )
    if args.out_csv:
        U.save_contacts_to_csv(contacts, args.out_csv)
        print(f"Saved contacts CSV -> {args.out_csv}")
    else:
        _print_or_save_json(contacts, None)


def cmd_contacts_batch(args: argparse.Namespace):
    U.process_contacts_folder(
        input_dir=args.in_dir,
        output_dir=args.out_dir,
        protein_chain_id=args.protein_chain,
        ligand_resname=args.ligand_resname,
        distance_threshold=args.cutoff
    )
    print(f"Batch contacts complete -> {args.out_dir}")


def cmd_rmsd(args: argparse.Namespace):
    rmsd = U.rmsd_between_ligand_and_chain(
        pdb_file=args.pdb,
        protein_chain_id=args.protein_chain,
        ligand_resname=args.ligand_resname,
        ligand_chain_id=args.ligand_chain if args.ligand_chain else None
    )
    if rmsd is None:
        print("No RMSD computed (missing ligand or protein atoms).")
        return
    print(f"RMSD: {rmsd:.3f} Å")


def cmd_binding_sites(args: argparse.Namespace):
    sites = U.find_ligand_binding_sites(args.pdb, distance_cutoff=args.cutoff)
    _print_or_save_json(sites, args.out_json)


def cmd_download_pdb(args: argparse.Namespace):
    path = U.download_pdb(args.pdb_id, args.out_dir)
    if path:
        print(f"Downloaded -> {path}")
    else:
        print("Download failed.")


def cmd_annotate_rcsb(args: argparse.Namespace):
    df = pd.read_csv(args.csv)
    df_out = U.annotate_df_with_rcsb(
        df=df,
        pdb_col=args.pdb_col,
        out_gene_col=args.gene_col,
        out_lig_col=args.lig_col
    )
    df_out.to_csv(args.out, index=False)
    print(f"Saved -> {args.out}")


def cmd_prepare_genenames(args: argparse.Namespace):
    df = pd.read_csv(args.input)
    df_out = U.prepare_genename_table(df)
    df_out.to_csv(args.out, index=False)
    print(f"Saved -> {args.out}")


def cmd_annotate_targets(args: argparse.Namespace):
    df = pd.read_csv(args.input)
    df_out = U.annotate_targets_and_malaria(
        df_targets=df,
        target_pdb_col=args.target_col,
        malaria_pdb_col=args.malaria_col
    )
    df_out.to_csv(args.out, index=False)
    print(f"Saved -> {args.out}")


def cmd_blast(args: argparse.Namespace):
    hits = U.fetch_and_blast_sequence(
        accession_id=args.accession,
        taxonomy=args.taxonomy,
        blast_db=args.db,
        blast_type=args.type,
        expect=args.expect,
        matrix_name=args.matrix,
        alignments=args.alignments,
        hitlist_size=args.hitlist_size,
        filter=args.filter,
        gapcosts=args.gapcosts,
        verbose=args.verbose,
        sleep_seconds=args.sleep
    )
    if args.out_json:
        _print_or_save_json(hits if hits is not None else [], args.out_json)
    else:
        _print_or_save_json(hits if hits is not None else [], None)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Protein/Ligand CLI over Biopython + RCSB/UniProt/NCBI helpers (see --help)."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # parse-entities
    sp = sub.add_parser("parse-entities", help="Parse chains/residues/sequence categories.")
    sp.add_argument("--pdb", required=True, help="Path to PDB file")
    sp.add_argument("--out-json", help="Optional JSON output path")
    sp.set_defaults(func=cmd_parse_entities)

    # filter
    sp = sub.add_parser("filter", help="Keep chain (protein) + ligand (on chain).")
    sp.add_argument("--in", dest="input", required=True, help="Input PDB")
    sp.add_argument("--out", dest="output", required=True, help="Output PDB")
    sp.add_argument("--protein-chain", default="M")
    sp.add_argument("--ligand-chain", default="T")
    sp.add_argument("--ligand-resname", default="D")
    sp.set_defaults(func=cmd_filter)

    # filter-exclude
    sp = sub.add_parser("filter-exclude", help="Remove ATOM lines from a chain (e.g., T).")
    sp.add_argument("--in", dest="input", required=True, help="Input PDB")
    sp.add_argument("--out", dest="output", required=True, help="Output PDB")
    sp.add_argument("--exclude-chain", default="T")
    sp.set_defaults(func=cmd_filter_exclude)

    # batch-filter
    sp = sub.add_parser("batch-filter", help="Batch filter PDBs in a folder (keep protein+ligand predicate).")
    sp.add_argument("--in-dir", required=True)
    sp.add_argument("--out-dir", required=True)
    sp.add_argument("--protein-chain", default="M")
    sp.add_argument("--ligand-chain", default="T")
    sp.add_argument("--ligand-resname", default="D")
    sp.set_defaults(func=cmd_batch_filter)

    # contacts
    sp = sub.add_parser("contacts", help="Ligand–CA contacts (JSON to stdout or CSV).")
    sp.add_argument("--pdb", required=True)
    sp.add_argument("--out-csv", help="Write contacts to CSV (otherwise prints JSON)")
    sp.add_argument("--protein-chain", default="M")
    sp.add_argument("--ligand-resname", default="D")
    sp.add_argument("--cutoff", type=float, default=5.0)
    sp.set_defaults(func=cmd_contacts)

    # contacts-batch
    sp = sub.add_parser("contacts-batch", help="Batch compute contacts for a folder.")
    sp.add_argument("--in-dir", required=True)
    sp.add_argument("--out-dir", required=True)
    sp.add_argument("--protein-chain", default="M")
    sp.add_argument("--ligand-resname", default="D")
    sp.add_argument("--cutoff", type=float, default=5.0)
    sp.set_defaults(func=cmd_contacts_batch)

    # rmsd
    sp = sub.add_parser("rmsd", help="Heuristic RMSD (ligand vs protein chain).")
    sp.add_argument("--pdb", required=True)
    sp.add_argument("--protein-chain", default="M")
    sp.add_argument("--ligand-resname", default="D")
    sp.add_argument("--ligand-chain", help="Optional ligand chain to restrict")
    sp.set_defaults(func=cmd_rmsd)

    # binding-sites
    sp = sub.add_parser("binding-sites", help="Ligand binding sites within cutoff.")
    sp.add_argument("--pdb", required=True)
    sp.add_argument("--cutoff", type=float, default=5.0)
    sp.add_argument("--out-json", help="Optional JSON output path")
    sp.set_defaults(func=cmd_binding_sites)

    # download-pdb
    sp = sub.add_parser("download-pdb", help="Download a PDB by 4-letter ID.")
    sp.add_argument("--pdb-id", required=True)
    sp.add_argument("--out-dir", required=True)
    sp.set_defaults(func=cmd_download_pdb)

    # annotate-rcsb
    sp = sub.add_parser("annotate-rcsb", help="Annotate CSV with RCSB gene/ligand lists from a PDB column.")
    sp.add_argument("--csv", required=True)
    sp.add_argument("--pdb-col", required=True, help="Column name containing PDB IDs (optionally with chain)")
    sp.add_argument("--out", required=True)
    sp.add_argument("--gene-col", default="gene_names", help="Output column for gene names list")
    sp.add_argument("--lig-col", default="ligands", help="Output column for ligands list")
    sp.set_defaults(func=cmd_annotate_rcsb)

    # prepare-genenames
    sp = sub.add_parser("prepare-genenames", help="Augment PlasmoDB gene table (UniProt/AF/RCSB by gene).")
    sp.add_argument("--in", dest="input", required=True)
    sp.add_argument("--out", required=True)
    sp.set_defaults(func=cmd_prepare_genenames)

    # annotate-targets
    sp = sub.add_parser("annotate-targets", help="Annotate targets/malaria CSV with RCSB (adds *_gene_names, *_ligands).")
    sp.add_argument("--in", dest="input", required=True)
    sp.add_argument("--out", required=True)
    sp.add_argument("--target-col", default="target_chain_id")
    sp.add_argument("--malaria-col", default="malaria_match_id")
    sp.set_defaults(func=cmd_annotate_targets)

    # blast
    sp = sub.add_parser("blast", help="Online BLAST of a UniProt accession (respect NCBI usage limits).")
    sp.add_argument("--accession", required=True)
    sp.add_argument("--taxonomy", type=int, required=True, help="NCBI Taxonomy ID (e.g., 5833 for Plasmodium)")
    sp.add_argument("--db", default="nr")
    sp.add_argument("--type", default="blastp")
    sp.add_argument("--expect", type=float, default=0.1)
    sp.add_argument("--matrix", default="BLOSUM62")
    sp.add_argument("--alignments", type=int, default=50)
    sp.add_argument("--hitlist-size", type=int, default=50)
    sp.add_argument("--filter", default="F")
    sp.add_argument("--gapcosts", default="11 1")
    sp.add_argument("--sleep", type=float, default=0.0, help="Sleep seconds before submitting to NCBI (be gentle)")
    sp.add_argument("--out-json", help="Write results to JSON")
    sp.add_argument("--verbose", action="store_true")
    sp.set_defaults(func=cmd_blast)

    return p


def main(argv: Optional[List[str]] = None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
