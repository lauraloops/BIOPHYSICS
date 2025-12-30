#!/usr/bin/env python3
import pandas as pd

DDG_CSV = "results/alanine_scanning/alanine_ddg.csv"
PDB_IN  = "data/6m0j_prepared.pdb"
PDB_OUT = "results/alanine_scanning/6m0j_ddg_bfactor.pdb"

df = pd.read_csv(DDG_CSV)
ddg_map = {(row["Chain"], int(row["Residue"])): float(row["ΔΔG (kcal/mol)"])
           for _, row in df.iterrows()}

def set_bfactor(line, b):
    # PDB B-factor is columns 61-66 (1-indexed); format width 6, 2 decimals
    b_str = f"{b:6.2f}"
    return line[:60] + b_str + line[66:]

with open(PDB_IN) as fin, open(PDB_OUT, "w") as fout:
    for line in fin:
        if line.startswith(("ATOM  ", "HETATM")):
            chain = line[21].strip()
            resi = int(line[22:26])
            b = ddg_map.get((chain, resi), 0.0)
            line = set_bfactor(line, b)
        fout.write(line)

print(f"Wrote: {PDB_OUT}")
