import os
import subprocess
import csv
from Bio.PDB import PDBParser, PDBIO

# -----------------------
# Paths
# -----------------------
DATA = "data"
SRC = "src"
RESULTS = "results/alanine_scanning"

WT_PDB = f"{DATA}/6m0j_prepared.pdb"
WT_PDBQT = f"{DATA}/6m0j_prepared.pdbqt"
WT_CSV = "results/WT/WT_interaction_energies.csv"

INTERFACE_A = "results/interface/interface_chain_A.txt"
INTERFACE_E = "results/interface/interface_chain_E.txt"

ENERGY_SCRIPT = f"{SRC}/int_energies_AE.py"

os.makedirs(f"{RESULTS}/mutants", exist_ok=True)
os.makedirs(f"{RESULTS}/energies", exist_ok=True)

BACKBONE_KEEP = {"N", "CA", "C", "O", "CB"}


# -----------------------
# Helpers
# -----------------------
def load_interface(path):
    with open(path) as f:
        return [int(x) for x in f.read().split() if x.isdigit()]


def mutate_to_alanine(pdb_in, chain_id, resid, pdb_out):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("mut", pdb_in)

    if chain_id not in structure[0]:
        raise KeyError(f"Chain {chain_id} not found in {pdb_in}")

    found = False
    for res in structure[0][chain_id]:
        if res.id[0] != " ":
            continue
        if res.id[1] == resid:
            found = True
            res.resname = "ALA"
            for atom in list(res):
                if atom.name not in BACKBONE_KEEP:
                    res.detach_child(atom.id)
            break

    if not found:
        raise ValueError(f"Residue {resid} not found in chain {chain_id}")

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_out)


def pdb_to_pdbqt(pdb, pdbqt):
    subprocess.run(
        ["obabel", pdb, "-O", pdbqt, "--partialcharge", "gasteiger"],
        check=True
    )


def extract_total_energy(csv_file):
    """
    Accepts either:
      TOTAL, ... , <value>
    or
      TOTAL(A–E),, <elec>, <vdw>, <solv>, <total>
    """
    with open(csv_file) as f:
        for row in csv.reader(f):
            if not row:
                continue
            key = row[0].strip()

            if key == "TOTAL":
                # old style
                return float(row[-1])

            if key.startswith("TOTAL(") or key.startswith("TOTAL"):
                # new style: TOTAL(A–E),,elec,vdw,solv,total
                try:
                    return float(row[-1])
                except Exception:
                    pass

    raise RuntimeError(f"Total energy not found in {csv_file}")


def run_energy(pdb, pdbqt, out_csv):
    env = dict(os.environ)
    env["PDB"] = pdb
    env["PDBQT"] = pdbqt
    env["OUTCSV"] = out_csv

    subprocess.run(
        ["python3", ENERGY_SCRIPT],
        check=True,
        env=env
    )

    if not os.path.exists(out_csv):
        raise FileNotFoundError(f"Energy script did not write OUTCSV: {out_csv}")


# -----------------------
# MAIN
# -----------------------
def main():
    iface_A = load_interface(INTERFACE_A)
    iface_E = load_interface(INTERFACE_E)

    # Safety filter (optional)
    iface_A = [r for r in iface_A if r >= 10]
    iface_E = [r for r in iface_E if r >= 10]

    WT_energy = extract_total_energy(WT_CSV)

    ddg_results = []

    # --- Chain A scan ---
    for resid in iface_A:
        tag = f"A_{resid}_ALA"
        pdb = f"{RESULTS}/mutants/{tag}.pdb"
        pdbqt = f"{RESULTS}/mutants/{tag}.pdbqt"
        csv_out = f"{RESULTS}/energies/{tag}.csv"

        try:
            mutate_to_alanine(WT_PDB, "A", resid, pdb)
            pdb_to_pdbqt(pdb, pdbqt)
            run_energy(pdb, pdbqt, csv_out)

            mut_energy = extract_total_energy(csv_out)
            ddg = mut_energy - WT_energy
            ddg_results.append(("A", resid, ddg))

        except Exception as e:
            print(f"Skipping {tag}: {e}")
            continue

    # --- Chain E scan ---
    for resid in iface_E:
        tag = f"E_{resid}_ALA"
        pdb = f"{RESULTS}/mutants/{tag}.pdb"
        pdbqt = f"{RESULTS}/mutants/{tag}.pdbqt"
        csv_out = f"{RESULTS}/energies/{tag}.csv"

        try:
            mutate_to_alanine(WT_PDB, "E", resid, pdb)
            pdb_to_pdbqt(pdb, pdbqt)
            run_energy(pdb, pdbqt, csv_out)

            mut_energy = extract_total_energy(csv_out)
            ddg = mut_energy - WT_energy
            ddg_results.append(("E", resid, ddg))

        except Exception as e:
            print(f"Skipping {tag}: {e}")
            continue

    # Save ΔΔG table
    out_ddg = f"{RESULTS}/alanine_ddg.csv"
    with open(out_ddg, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Chain", "Residue", "ΔΔG (kcal/mol)"])
        for row in ddg_results:
            w.writerow(row)

    print(f"Alanine scanning finished → {out_ddg}")


if __name__ == "__main__":
    main()
