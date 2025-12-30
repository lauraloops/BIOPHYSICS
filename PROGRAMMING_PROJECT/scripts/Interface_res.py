### Identify interface residues between two chains in a PDB structure using Biopython
### 1. Load the PDB structure
### 2. Extract the two chains of interest
### 3. For each residue in chain A, check if any atom is within cutoff distance of any atom in chain B
### 4. Collect and print the interface residues for both chains


from Bio.PDB import *
import os

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
pdb_path = os.path.join(base_dir, 'data', "6m0j_clean_from_pymol.pdb")

pdb = PDBParser().get_structure("complex", pdb_path)
model = pdb[0]
chainA = model['A']   # RBD
chainB = model['E']   # ACE2

cutoff = 6.0  
### to avoid duplicates we use sets to store interface residues
interface_A = set()
interface_B = set()

for resA in chainA:
    for atomA in resA:
        for resB in chainB:
            for atomB in resB:
                if atomA - atomB < cutoff: 
                    interface_A.add(resA.get_id()[1])
                    interface_B.add(resB.get_id()[1])

print("Interface residues A:", sorted(interface_A))
print("Interface residues B:", sorted(interface_B))
