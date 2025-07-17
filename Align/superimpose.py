from Bio.PDB import PDBParser, Superimposer, Polypeptide
from Bio.Align import PairwiseAligner
import numpy as np
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.PDB.PDBIO import PDBIO


def extract_sequence_and_coords(chain):
    """Extract the sequence and CA coordinates from a chain."""
    seq = []
    coords = []
    res_ids = []
    residues=[]
    one_letter_seq=""
    for res in chain:
        if Polypeptide.is_aa(res, standard=True) and 'CA' in res:
            seq.append(res.resname)
            coords.append(res['CA'])
            res_ids.append(res.get_id()[1])  # sequence number
            residues.append(res)
    from Bio.Data.IUPACData import protein_letters_3to1
    try:
        one_letter_seq = ''.join([protein_letters_3to1.get(res.get_resname().capitalize(), 'X') for res in residues])
    except Exception as e:
        one_letter_seq = ""
    return one_letter_seq, coords, res_ids

def get_aligned_coords(seq1, coords1, seq2, coords2):
    blosum62 = substitution_matrices.load("BLOSUM62")
    """Align two sequences and return aligned coordinates."""
    alignments = pairwise2.align.localds(seq1, seq2, blosum62, -11, -1,one_alignment_only=True)
    aligned1, aligned2, _, _, _ = alignments[0]
    aligned_coords1 = []
    aligned_coords2 = []
    print(seq1)
    print(aligned2)
    i1 = i2 = 0
    for a1, a2 in zip(aligned1, aligned2):
        if a1 != '-' and a2 != '-':
            aligned_coords1.append(coords1[i1])
            aligned_coords2.append(coords2[i2])
        if a1 != '-':
            i1 += 1
        if a2 != '-':
            i2 += 1
    print(type(aligned_coords1))
    return aligned_coords1, aligned_coords2

def iterative_alignment(fixed_coords, moving_coords, max_iterations=10, distance_cutoff=4.0):
    """Iterative alignment with pruning based on 5Å distance cutoff."""
  #  assert fixed_coords.shape == moving_coords.shape
    fc = fixed_coords
    mc = moving_coords
  
    for iteration in range(max_iterations):
        indices = np.arange(0,len(fc))
        fixed_c = []
        moving_c = []
        si = Superimposer()
        si.set_atoms(fc, mc)
        si.apply(mc)
        rmsd = si.rms
        print(f"Final RMSD: {rmsd:.3f} Å with {len(indices)} residues after {iteration+1} iterations.")

        for i in indices:
           fixed_c.append(fc[i].get_coord())
           moving_c.append(mc[i].get_coord())
        distances = np.linalg.norm(np.array(fixed_c) - np.array(moving_c), axis=1)
        new_indices = np.where(distances <= distance_cutoff)[0]
        if len(new_indices) == len(indices):
            print(f"Converged after {iteration+1} iterations.")
            break
        if len(new_indices) < 3:
            print("Too few residues to align. Stopping.")
            break
        fcn = []
        mcn = []
        for i in new_indices:
           fcn.append(fc[i])
           mcn.append(mc[i])
        fc = fcn
        mc = mcn

    rmsd = si.rms
    print(f"Final RMSD: {rmsd:.3f} Å with {len(indices)} residues.")
    si.apply(moving_coords)
    return fixed_coords, moving_coords, rmsd


