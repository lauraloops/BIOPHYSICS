#!/usr/bin/env python3
import argparse
from Bio.PDB import PDBParser, NeighborSearch, Selection

def get_interface_residues(pdb_file, chain1_id, chain2_id, cutoff):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)
    model = structure[0]

    try:
        chain1 = model[chain1_id]
        chain2 = model[chain2_id]
    except KeyError:
        raise ValueError(f"Chain IDs {chain1_id} or {chain2_id} not found in {pdb_file}")

    atoms_chain1 = list(chain1.get_atoms())
    atoms_chain2 = list(chain2.get_atoms())

    # NeighborSearch needs all atoms; we'll use it twice (both directions)
    all_atoms = atoms_chain1 + atoms_chain2
    ns = NeighborSearch(all_atoms)

    iface_residues_chain1 = set()
    iface_residues_chain2 = set()

    # For each atom in chain1, search neighbors in chain2
    for atom in atoms_chain1:
        neighbors = ns.search(atom.get_coord(), cutoff, level="A")
        for neigh in neighbors:
            parent_chain = neigh.get_parent().get_parent()  # residue -> chain
            if parent_chain.id == chain2_id:
                iface_residues_chain1.add(atom.get_parent())
                iface_residues_chain2.add(neigh.get_parent())

    return iface_residues_chain1, iface_residues_chain2

def residue_id(res):
    """Return a nice string for a Biopython residue object."""
    resname = res.get_resname()
    chain_id = res.get_parent().id
    seq_id = res.get_id()[1]  # (hetflag, resseq, icode) -> use resseq
    icode = res.get_id()[2].strip()
    if icode:
        return f"{resname} {chain_id}{seq_id}{icode}"
    else:
        return f"{resname} {chain_id}{seq_id}"

def main():
    parser = argparse.ArgumentParser(
        description="Find interface residues between two chains in a PDB file."
    )
    parser.add_argument("pdb", help="Input PDB file (clean structure)")
    parser.add_argument("chain1", help="First chain ID (e.g. A)")
    parser.add_argument("chain2", help="Second chain ID (e.g. E)")
    parser.add_argument(
        "-d", "--distance", type=float, default=5.0,
        help="Distance cutoff in angstroms (default: 5.0 Å)"
    )
    args = parser.parse_args()

    iface1, iface2 = get_interface_residues(
        args.pdb, args.chain1, args.chain2, args.distance
    )

    print(f"Interface residues on chain {args.chain1} (within {args.distance} Å of chain {args.chain2}):")
    for res in sorted(iface1, key=lambda r: r.get_id()[1]):
        print("  ", residue_id(res))

    print()
    print(f"Interface residues on chain {args.chain2} (within {args.distance} Å of chain {args.chain1}):")
    for res in sorted(iface2, key=lambda r: r.get_id()[1]):
        print("  ", residue_id(res))

if __name__ == "__main__":
    main()
