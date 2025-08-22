# Alignment

Pipeline for protein alignment and downstream analysis (BLAST, structural tools, docking helpers).

## Quickstart
```bash
# clone
git clone git@github.com:kattens/Alignment.git
cd Alignment

# create env (conda or venv)
conda create -n align python=3.11 -y && conda activate align
pip install -r requirements.txt

# run an example
python run.py --help
