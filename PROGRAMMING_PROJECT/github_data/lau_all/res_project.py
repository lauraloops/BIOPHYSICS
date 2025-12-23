from Bio.PDB import PDBParser, NeighborSearch
import argparse

def get_interface_residues(pdb_file, chain1_id, chain2_id, cutoff=5.5):
    """
    Returns interface residues between chain1 and chain2 
    using atom–atom distance <= cutoff (Å).
    """

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)

    # Extract chains
    chain1 = structure[0][chain1_id]
    chain2 = structure[0][chain2_id]

    # Collect atoms
    atoms_chain1 = [atom for atom in chain1.get_atoms() if atom.element != 'H']
    atoms_chain2 = [atom for atom in chain2.get_atoms() if atom.element != 'H']

    # Build neighbor search
    ns = NeighborSearch(atoms_chain1 + atoms_chain2)

    interface1 = set()
    interface2 = set()

    # Search for contacts
    for atom1 in atoms_chain1:
        neighbors = ns.search(atom1.coord, cutoff)
        for atom2 in neighbors:
            # skip atoms of the same chain
            if atom2 in atoms_chain1:
                continue

            res1 = atom1.get_parent()
            res2 = atom2.get_parent()
            interface1.add((chain1_id, res1.get_resname(), res1.id[1]))
            interface2.add((chain2_id, res2.get_resname(), res2.id[1]))

    return sorted(interface1, key=lambda x: x[2]), sorted(interface2, key=lambda x: x[2])



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find interface residues between two chains")
    parser.add_argument("pdb", help="Input PDB file")
    parser.add_argument("chain1", help="First chain ID")
    parser.add_argument("chain2", help="Second chain ID")
    parser.add_argument("--cutoff", type=float, default=5.5, help="Distance cutoff (Å)")

    args = parser.parse_args()

    interface1, interface2 = get_interface_residues(args.pdb, args.chain1, args.chain2, args.cutoff)

    print(f"\nInterface residues in chain {args.chain1}:")
    for res in interface1:
        print(f"{res[0]} {res[1]} {res[2]}")

    print(f"\nInterface residues in chain {args.chain2}:")
    for res in interface2:
        print(f"{res[0]} {res[1]} {res[2]}")
