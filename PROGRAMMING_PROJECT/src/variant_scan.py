#!/usr/bin/env python3
import os
import csv
import subprocess
import tempfile
from pathlib import Path

# -----------------------
# Paths
# -----------------------
DATA = "data"
SRC = "src"
RESULTS = "results/variants"

WT_PDB = f"{DATA}/6m0j_prepared.pdb"
WT_PDBQT = f"{DATA}/6m0j_prepared.pdbqt"
WT_CSV = "results/WT/WT_interaction_energies.csv"
ENERGY_SCRIPT = f"{SRC}/int_energies_AE.py"

Path(f"{RESULTS}/mutants").mkdir(parents=True, exist_ok=True)
Path(f"{RESULTS}/energies").mkdir(parents=True, exist_ok=True)

# -----------------------
# Mutations (Variants of concern)
# Chain E is Spike RBD (correct mapping)
# -----------------------
MUTATIONS = [
    ("Alpha", "E", 501, "TYR"),  # N501Y
    ("Beta",  "E", 417, "ASN"),  # K417N
    ("Beta",  "E", 484, "LYS"),  # E484K
    ("Beta",  "E", 501, "TYR"),  # N501Y
    ("Delta", "E", 452, "ARG"),  # L452R
    ("Delta", "E", 478, "LYS"),  # T478K
]

AA3_TO_AA1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}

# -----------------------
# Helpers
# -----------------------
def extract_total_energy(csv_file: str) -> float:
    """
    Reads TOTAL(A–E) line from your CSV.
    Example: TOTAL(A–E),,elec,vdw,solv,total
    """
    with open(csv_file) as f:
        r = csv.reader(f)
        for row in r:
            if row and row[0].startswith("TOTAL"):
                return float(row[-1])
    raise RuntimeError(f"TOTAL energy not found in {csv_file}")

def pdb_to_pdbqt(pdb: str, pdbqt: str):
    subprocess.run(
        ["obabel", pdb, "-O", pdbqt, "--partialcharge", "gasteiger"],
        check=True
    )

def run_energy(pdb: str, pdbqt: str, out_csv: str):
    subprocess.run(
        ["python3", ENERGY_SCRIPT],
        check=True,
        env=dict(os.environ, PDB=pdb, PDBQT=pdbqt, OUTCSV=out_csv)
    )

def mutate_with_pymol(pdb_in: str, chain: str, resid: int, new_aa3: str, pdb_out: str):
    """
    Uses PyMOL mutagenesis wizard (headless) to rebuild sidechain.
    Requires 'pymol' executable in PATH.
    """
    selection = f"chain {chain} and resi {resid}"

    pml = f"""
load {pdb_in}, wt
remove hetatm
wizard mutagenesis
refresh_wizard

python
from pymol import cmd
w = cmd.get_wizard()
w.set_mode("{new_aa3}")
cmd.do("select mut_sel, {selection}")
w.do_select("mut_sel")
# Apply with default rotamer choice (PyMOL picks a reasonable one)
w.apply()
cmd.set_wizard()
python end

save {pdb_out}, wt
quit
"""
    with tempfile.NamedTemporaryFile("w", suffix=".pml", delete=False) as tf:
        tf.write(pml)
        pml_path = tf.name

    try:
        subprocess.run(["pymol", "-cq", pml_path], check=True)
    finally:
        try:
            os.remove(pml_path)
        except OSError:
            pass

# -----------------------
# MAIN
# -----------------------
def main():
    # WT reference (from your already computed WT CSV)
    wt_energy = extract_total_energy(WT_CSV)

    out_csv = f"{RESULTS}/variant_ddg.csv"
    rows = []

    for variant, chain, resid, newaa3 in MUTATIONS:
        tag = f"{variant}_{chain}{resid}_{newaa3}"
        mut_pdb = f"{RESULTS}/mutants/{tag}.pdb"
        mut_pdbqt = f"{RESULTS}/mutants/{tag}.pdbqt"
        e_csv = f"{RESULTS}/energies/{tag}.csv"

        print(f"\n=== Building {tag} ===")
        mutate_with_pymol(WT_PDB, chain, resid, newaa3, mut_pdb)
        pdb_to_pdbqt(mut_pdb, mut_pdbqt)

        print(f"=== Energy {tag} ===")
        run_energy(mut_pdb, mut_pdbqt, e_csv)

        mut_energy = extract_total_energy(e_csv)
        ddg = mut_energy - wt_energy

        rows.append([variant, chain, resid, newaa3, mut_energy, wt_energy, ddg])

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Variant", "Chain", "Residue", "Mut_AA3",
                    "ΔG_variant (kcal/mol)", "ΔG_WT (kcal/mol)", "ΔΔG_variant (kcal/mol)"])
        w.writerows(rows)

    print(f"\n✅ Variant scan done → {out_csv}")

if __name__ == "__main__":
    main()
