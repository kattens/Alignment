# Alignment

Pipeline for protein alignment and downstream analysis (BLAST, structural tools, docking helpers).
---

## 🔬 Drug Repurposing Pipeline Overview

This repository implements a **modular, stage-based pipeline** for drug repurposing research.
Each stage is isolated in its own folder (`packages/<stage>`), with a separate conda environment (`envs/<stage>.yml`), YAML config (`configs/<stage>.yaml`), and CLI entrypoint (`<stage>-stage`).

The pipeline can be run **end-to-end** or **stage-by-stage** (rerun only the part you need).
Intermediate data lives under `data/processed/`, results under `results/`. Raw data is kept in `data/raw/` and ignored by Git.

---

### 📑 Stages

| #      | Stage                              | Description                                                                                          | Inputs → Outputs                                                                                                                  |
| ------ | ---------------------------------- | ---------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------- |
| **1**  | **Ingest (PubChem preprocessing)** | Download and clean PubChem data (CSV/SMILES, QC).                                                    | `data/raw/pubchem/*.csv` → `data/processed/ingest/pubchem_clean.csv`                                                              |
| **2**  | **BLAST search**                   | Run BLAST queries against sequence DB.                                                               | `data/processed/sequences/query.fasta`, DB → `results/blast/out.tsv`                                                              |
| **3**  | **BLAST post-processing**          | Parse and filter BLAST results (top hits, merges).                                                   | `results/blast/out.tsv` → `data/processed/blast/blast_clean.csv`                                                                  |
| **4**  | **Gene–Ligand–PDB join**           | Map genes to ligands and PDBs; assemble metadata.                                                    | `data/processed/blast/blast_clean.csv` + lookups → `data/processed/join/gene_ligand_map.csv`                                      |
| **5**  | **Chain separation**               | Split target/malaria PDBs into individual chains.                                                    | `data/processed/pdb/*.pdb` → `data/processed/pdb_split/{chain}.pdb`                                                               |
| **6**  | **Chimera operations**             | Headless Chimera/ChimeraX for renaming + structural alignment.                                       | `data/processed/..._renamed/*.pdb` → `results/chimera/aligned_*.pdb`, logs                                                        |
| **7**  | **Alignment post-processing**      | Parse alignment logs; extract metrics (RMSD, TM-score, etc.).                                        | `results/chimera/*.pdb` + logs → `results/align/metrics.csv`                                                                      |
| **8**  | **Reporting (intermediate)**       | Generate summary plots and tables from processed alignment results.                                  | `data/processed/*` + `results/align/*` → `results/report/alignments/*`                                                            |
| **9**  | **AutoDock Vina**                  | Dock ligands into predicted binding sites; generate binding affinities.                              | `data/processed/docking/receptors/*.pdbqt`, ligands, grid boxes → `results/autodock/poses/*.pdbqt`, `results/autodock/scores.csv` |
| **10** | **AlphaFold 3 inference**          | Predict protein–ligand or protein–protein complexes with AlphaFold 3 (API/local).                    | sequences/structures + ligands → `results/af3/predictions/*`, `results/af3/metrics.csv`                                           |
| **11** | **Docking vs AF3 comparison**      | Compare Vina docking scores with AlphaFold 3 confidence metrics (e.g., pocket overlap, ipTM, pLDDT). | `results/autodock/scores.csv` + `results/af3/metrics.csv` → `results/compare/summary.csv`, plots                                  |

---

### ⚙️ How it works

* Each stage is an **installable mini-package** (`pyproject.toml`) with a Typer CLI.
* Each has its own **conda environment** (`envs/<stage>.yml`) for reproducibility.
* Parameters and file paths are controlled via **YAML configs** (`configs/<stage>.yaml`).
* A central runner script (`scripts/run_stage.sh`) lets you run one stage at a time:

  ```bash
  bash scripts/run_stage.sh preprocess configs/preprocess.yaml
  ```
* A `Makefile` provides shortcuts:

  ```bash
  make blast         # run stage 2
  make preprocess    # run stage 3
  make autodock      # run stage 9
  make af3           # run stage 10
  make compare       # run stage 11
  ```

---

### 📂 Folder structure

```
.
├── envs/            # per-stage conda envs
├── configs/         # per-stage YAML configs
├── data/            
│   ├── raw/         # input datasets (gitignored)
│   └── processed/   # cleaned outputs per stage (gitignored)
├── results/         # figures, logs, docking poses, AF3 preds (gitignored)
├── packages/        # each stage = installable Python package
│   ├── ingest_pubchem/
│   ├── blast/
│   ├── blast_preprocess/
│   ├── gene_ligand_join/
│   ├── chain_separation/
│   ├── chimera_ops/
│   ├── align/
│   ├── report/
│   ├── autodock/            # AutoDock Vina stage
│   ├── alphafold3_infer/    # AlphaFold 3 inference stage
│   └── compare_docking/     # Vina vs AF3 comparison
├── scripts/
│   └── run_stage.sh
└── tests/
 ```

## Typical run order:
make ingest_pubchem
make blast
make blast_preprocess
make gene_ligand_join
make chain_separation
make chimera_ops
make align
make report
make autodock
make af3
make compare
