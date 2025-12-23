from biobb_structure_checking.structure_checking import StructureChecking
import os
import biobb_structure_checking

# Base directory
base_dir = os.path.dirname(biobb_structure_checking.__file__)
project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Input PDB
input_pdb = os.path.join(project_dir, 'data', "6m0j_clean.pdb")

# Arguments for StructureChecking
args_dict = {
    "input_structure_path": input_pdb,
    "non_interactive": True,
    "output_format": "pdb"
}

# Create StructureChecking object (automatically loads structure)
sc = StructureChecking(base_dir_path=base_dir, args=args_dict)

# Run all checks (missing atoms, altloc, clashes, etc.)
print("\nRunning checkall...")
sc.checkall()

# Fix everything that can be fixed
print("\nRunning fixall...")
sc.fixall()

# Add hydrogen atoms
print("\nAdding hydrogens...")
sc.add_hydrogen()

# Save the cleaned structure to a new file
output_pdb = "6m0j_prepared.pdb"
print(f"\nSaving repaired structure to {output_pdb}")
sc.save_structure(output_pdb)

print("\n✓ Finished cleaning PDB file.")
print("Next: Generate PDBQT using OpenBabel.")
