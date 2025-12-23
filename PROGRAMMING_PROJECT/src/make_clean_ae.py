from Bio.PDB import PDBParser, PDBIO, Select
import os

class ChainSelect(Select):
    def accept_chain(self, chain):
        return chain.id in ("A", "E")
    def accept_residue(self, residue):
        # keep only standard residues (remove HETATM waters/ions/ligands)
        return residue.id[0] == " "

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(base_dir, 'data')
input_pdb = os.path.join(data_dir, "6m0j_raw.pdb")
output_pdb = os.path.join(data_dir, "6m0j_clean.pdb")

parser = PDBParser(QUIET=True)
structure = parser.get_structure("6M0J", input_pdb)

io = PDBIO()
io.set_structure(structure)
io.save(output_pdb, ChainSelect())
print(f"Wrote {output_pdb} (chains A+E, no heteroatoms)")
