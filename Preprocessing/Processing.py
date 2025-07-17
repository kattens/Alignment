from Bio.PDB import PDBParser, PDBIO, Select

class ParsePDB:
    def __init__(self, pdb_name, chain_id, drug_id):
        self.pdb_name = pdb_name
        self.chain_id = chain_id
        self.drug_id = drug_id
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure("struct", pdb_name)

    def accept_chain(self, chain):
        return chain.id == self.chain_id    

    def accept_residue(self, residue):
        hetfield, _, _ = residue.get_id()
        return (hetfield == ' ' or residue.get_resname() == self.drug_id)

    def save_pdb(self, out_file):
        io = PDBIO()
        io.set_structure(self.structure)
        io.save(out_file)

    def keep_target(self, out_file):
        class ChainSelector(Select):
            def accept_chain(inner_self, chain):
                return self.accept_chain(chain)
        io = PDBIO()
        io.set_structure(self.structure)
        io.save(out_file, ChainSelector())

    def keep_drug(self, input_file, out_file):
        struct = self.parser.get_structure("filtered", input_file)
        class ResidueSelector(Select):
            def accept_residue(inner_self, residue):
                return self.accept_residue(residue)
        io = PDBIO()
        io.set_structure(struct)
        io.save(out_file, ResidueSelector())

    def keep_malaria(self, input_file, chain_id, out_file):
        struct = self.parser.get_structure("malaria", input_file)
        for model in struct:
            chains_to_remove = []
            for chain in model:
                if chain.id == chain_id:
                    chain.id = 'M'
                else:
                    chains_to_remove.append(chain)
            for ch in chains_to_remove:
                model.detach_child(ch.id)
        io = PDBIO()
        io.set_structure(struct)
        io.save(out_file)

    def rename_target_ligand(self, input_file, chain_id, out_file):
        struct = self.parser.get_structure("target", input_file)
        for model in struct:
            for chain in model:
                if chain.id == chain_id:
                    chain.id = 'T'
                for residue in chain:
                    hetfield, _, _ = residue.get_id()
                    if hetfield != ' ' and residue.get_resname() == self.drug_id:
                        residue.resname = 'D'
        io = PDBIO()
        io.set_structure(struct)
        io.save(out_file)



