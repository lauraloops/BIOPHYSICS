#!/usr/bin/env python

"""
    Initial structure setup and annotation for interaction energy evaluation
    FIXED VERSION — compatible with OpenBabel-generated PDBQT files
"""

import argparse
import os

from Bio.PDB.PDBParser import PDBParser
from forcefield import VdwParamset

parser = argparse.ArgumentParser(
    prog='structure_setup',
    description='basic structure setup'
)


parser.add_argument(
    '--vdw',
    action='store',
    dest='vdwprm_file',
    default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'vdwprm'),
    help='Vdw parameter file'
)

parser.add_argument('pdb_file', help='Input PDB', type=open)
parser.add_argument('pdbqt_file', help='Input PDBQT', type=open)
parser_pdb = PDBParser(PERMISSIVE=1)

parser.add_argument('--asaA', dest='asa_chain_A', required=True, help='NACCESS .asa file for chain A')
parser.add_argument('--asaE', dest='asa_chain_E', required=True, help='NACCESS .asa file for chain E')

args = parser.parse_args()

print("PDB File:", args.pdb_file.name)
print("PDBQT File:", args.pdbqt_file.name)
print("Parameters File:", args.vdwprm_file)

# ---------------------------------------------------------
# Load VdW parameters
# ---------------------------------------------------------
ff_params = VdwParamset(args.vdwprm_file)

# ---------------------------------------------------------
# Load PDB structure
# ---------------------------------------------------------

print('Parsing PDB', args.pdb_file.name)
st = parser_pdb.get_structure('STR', args.pdb_file.name)

# ---------------------------------------------------------
# Parse PDBQT properly (store raw lines instead of broken fixed columns)
# ---------------------------------------------------------
print(f"Parsing PDBQT {args.pdbqt_file.name}")
params = [{}]  # so that atom serial numbers match index

for line in args.pdbqt_file:
    line = line.rstrip()

    # Only store ATOM/HETATM lines
    if line.startswith("ATOM") or line.startswith("HETATM"):
        params.append({'raw': line})

# ---------------------------------------------------------
# Fix atom serial numbers so they match PDBQT order
# ---------------------------------------------------------
i = 1
for at in st.get_atoms():
    at.serial_number = i
    i += 1

# ---------------------------------------------------------
# Assign charge, atom type, VDW parameters
# ---------------------------------------------------------
total_charge = 0.0

for at in st.get_atoms():

    # Extract raw PDBQT line for this atom
    raw = params[at.serial_number]['raw']
    parts = raw.split()

    # Charge = second-to-last column
    charge_str = parts[-2]
    try:
        charge = float(charge_str)
    except:
        charge = 0.0   # fallback for malformed values

    # Atom type = last column
    atom_type = parts[-1]
    if atom_type not in ff_params.at_types:
        atom_type = "C"  # safe fallback

    at.xtra['charge'] = charge
    at.xtra['atom_type'] = atom_type
    at.xtra['vdw'] = ff_params.at_types[atom_type]

    total_charge += charge

print(f"Total Charge: {total_charge:8.2f}")

# ---------------------------------------------------------
# NACCESS residue ASA parsing (from .asa files)
# ---------------------------------------------------------
def parse_naccess_asa(path):
    """
    Parse NACCESS .asa output and return dict:
    ASA[(chain_id, resseq, icode)] = asa_value

    Notes:
    - NACCESS formats can vary slightly. We robustly extract:
      chain (1 char), residue number (int), insertion code (optional), ASA (float)
    - We ignore non-residue summary lines.
    """
    asa = {}
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            # NACCESS residue lines typically start with "RES"
            # Example patterns seen in the wild:
            # "RES GLN A  42   12.34 ..."
            # "RES GLN A  42A  12.34 ..." (insertion code)
            if not line.startswith("RES"):
                continue

            parts = line.split()
            # Minimum expected fields: RES, RESNAME, CHAIN, RESNUM, ASA
            # parts[0]=RES, parts[1]=GLN, parts[2]=A, parts[3]=42, parts[4]=12.34
            if len(parts) < 5:
                continue

            chain_id = parts[2]
            resnum_token = parts[3]

            # resnum_token may contain insertion code appended, e.g. "42A"
            # separate numeric prefix from optional trailing insertion code
            num_str = ""
            icode = " "
            for ch in resnum_token:
                if ch.isdigit() or (ch == "-" and not num_str):
                    num_str += ch
                else:
                    icode = ch
                    break

            try:
                resseq = int(num_str)
            except ValueError:
                continue

            try:
                asa_val = float(parts[4])
            except ValueError:
                continue

            asa[(chain_id, resseq, icode)] = asa_val

    return asa


# Read ASA for each chain (NACCESS outputs)
print("Reading NACCESS ASA files...")
asa_A = parse_naccess_asa(args.asa_chain_A)
asa_E = parse_naccess_asa(args.asa_chain_E)

# Merge into one lookup dict
ASA_res = {}
ASA_res.update(asa_A)
ASA_res.update(asa_E)

print(f"Loaded ASA entries: chain A = {len(asa_A)}, chain E = {len(asa_E)}")

# Attach residue ASA to the structure (residue-level)
missing = 0
total_res = 0

for chain in st[0]:
    for res in chain:
        hetflag, resseq, icode = res.id
        # only standard residues (ignore hetero / waters; should already be removed)
        if hetflag != " ":
            continue

        total_res += 1
        key = (chain.id, resseq, icode if icode else " ")
        asa_val = ASA_res.get(key, None)

        if asa_val is None:
            # Allow fallback ignoring insertion code if needed
            key2 = (chain.id, resseq, " ")
            asa_val = ASA_res.get(key2, None)

        if asa_val is None:
            missing += 1
            res.xtra["ASA_NACCESS"] = 0.0
        else:
            res.xtra["ASA_NACCESS"] = asa_val

print(f"Residues in structure: {total_res} | Missing ASA entries: {missing}")

# Optional: attach residue ASA to atoms for compatibility with existing downstream code.
# This does NOT create real atom-level ASA; it simply propagates residue ASA to its atoms.
for chain in st[0]:
    for res in chain:
        hetflag, resseq, icode = res.id
        if hetflag != " ":
            continue
        asa_val = res.xtra.get("ASA_NACCESS", 0.0)
        for atom in res:
            atom.xtra["EXP_NACCESS"] = asa_val


# ---------------------------------------------------------
# Print example atom data
# ---------------------------------------------------------
print("ASA check example A353:", st[0]['A'][353].xtra.get("ASA_NACCESS"))
print(vars(st[0]['A'][42]['N']))
print(vars(st[0]['A'][42]['N'].xtra['vdw']))
