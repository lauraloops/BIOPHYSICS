#!/usr/bin/env python3
import os
import csv
import math
from Bio.PDB import PDBParser
from forcefield import VdwParamset


# ---------------------------
# NACCESS RSA parsing (RES lines)
# Uses "Non-polar ABS" column (index 8)
# ---------------------------
def parse_naccess_rsa(path):
    asa = {}
    with open(path, "r") as f:
        for line in f:
            if not line.startswith("RES"):
                continue
            parts = line.split()
            if len(parts) < 10:
                continue

            chain_id = parts[2]
            resnum_token = parts[3]

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
                asa_val = float(parts[8])  # Non-polar ABS
            except ValueError:
                continue

            asa[(chain_id, resseq, icode)] = asa_val

    return asa


# ---------------------------
# Load and annotate structure: charges/types/vdw + ASA (bound & unbound)
# ---------------------------
def load_annotated_structure(pdb_file, pdbqt_file, vdw_file,
                            rsa_complex, rsa_chainA, rsa_chainE):
    ff_params = VdwParamset(vdw_file)

    parser = PDBParser(PERMISSIVE=1)
    st = parser.get_structure("STR", pdb_file)

    # Read PDBQT lines (ATOM/HETATM)
    params = [{}]
    with open(pdbqt_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("ATOM") or line.startswith("HETATM"):
                params.append({"raw": line})

    # Force serial numbers to match PDBQT order
    i = 1
    for at in st.get_atoms():
        at.serial_number = i
        i += 1

    # Assign charge/type/vdw
    total_charge = 0.0
    for at in st.get_atoms():
        raw = params[at.serial_number]["raw"]
        parts = raw.split()

        try:
            charge = float(parts[-2])
        except Exception:
            charge = 0.0

        atom_type = parts[-1]
        if atom_type not in ff_params.at_types:
            atom_type = "C"

        at.xtra["charge"] = charge
        at.xtra["atom_type"] = atom_type
        at.xtra["vdw"] = ff_params.at_types[atom_type]

        total_charge += charge

    print(f"Total charge (from PDBQT): {total_charge:.2f} e")

    # Load RSA: bound (complex) and unbound (chains)
    asa_complex = parse_naccess_rsa(rsa_complex)
    asa_A_unb   = parse_naccess_rsa(rsa_chainA)
    asa_E_unb   = parse_naccess_rsa(rsa_chainE)

    # Attach residue ASA
    missing_bound = 0
    missing_unb = 0
    total_res = 0

    model = st[0]
    for chain in model:
        for res in chain:
            hetflag, resseq, icode = res.id
            if hetflag != " ":
                continue

            total_res += 1
            key = (chain.id, resseq, icode if icode else " ")

            # bound ASA (complex RSA)
            asa_b = asa_complex.get(key)
            if asa_b is None:
                asa_b = asa_complex.get((chain.id, resseq, " "))
            if asa_b is None:
                missing_bound += 1
                asa_b = 0.0

            # unbound ASA (chain RSA)
            if chain.id == "A":
                asa_u = asa_A_unb.get(key)
                if asa_u is None:
                    asa_u = asa_A_unb.get((chain.id, resseq, " "))
            elif chain.id == "E":
                asa_u = asa_E_unb.get(key)
                if asa_u is None:
                    asa_u = asa_E_unb.get((chain.id, resseq, " "))
            else:
                asa_u = 0.0

            if asa_u is None:
                missing_unb += 1
                asa_u = 0.0

            res.xtra["ASA_BOUND"] = asa_b
            res.xtra["ASA_UNBOUND"] = asa_u

    print(f"ASA attached (Å^2). Residues={total_res}, missing bound={missing_bound}, missing unbound={missing_unb}")
    return st


# ---------------------------
# Energy terms
# ---------------------------
def coulomb_energy(atom_i, atom_j, k_elec=332.0636, eps_r=80.0):
    qi = atom_i.xtra["charge"]
    qj = atom_j.xtra["charge"]
    r = atom_i - atom_j
    if r == 0:
        return 0.0
    return k_elec * qi * qj / (eps_r * r)


def lennard_jones_energy(atom_i, atom_j):
    vdw_i = atom_i.xtra["vdw"]
    vdw_j = atom_j.xtra["vdw"]

    eps_ij = math.sqrt(vdw_i.eps * vdw_j.eps)
    sig_ij = 0.5 * (vdw_i.sig + vdw_j.sig)

    r = atom_i - atom_j
    if r == 0:
        return 0.0

    sr = sig_ij / r
    sr6 = sr**6
    sr12 = sr6**2
    return 4.0 * eps_ij * (sr12 - sr6)


def residue_pair_energy(resA, resE, cutoff=8.0):
    E_elec = 0.0
    E_vdw = 0.0
    for atA in resA:
        if not atA.element.strip():
            continue
        for atE in resE:
            if not atE.element.strip():
                continue
            r = atA - atE
            if r > cutoff:
                continue
            E_elec += coulomb_energy(atA, atE)
            E_vdw  += lennard_jones_energy(atA, atE)
    return E_elec, E_vdw


def solvation_energy_residue(res, asa_key):
    asa_res = res.xtra.get(asa_key, 0.0)
    atoms = [at for at in res if at.element.strip()]
    if not atoms or asa_res == 0.0:
        return 0.0

    sum_fsrf = 0.0
    for at in atoms:
        sum_fsrf += at.xtra["vdw"].fsrf

    asa_per_atom = asa_res / len(atoms)
    return sum_fsrf * asa_per_atom


def delta_solvation_residue(res):
    return solvation_energy_residue(res, "ASA_BOUND") - solvation_energy_residue(res, "ASA_UNBOUND")


# ---------------------------
# MAIN
# ---------------------------
def main():
    # Read inputs from environment (this is how your alanine scan calls it)
    pdb_file   = os.environ.get("PDB",   os.path.join("data", "6m0j_prepared.pdb"))
    pdbqt_file = os.environ.get("PDBQT", os.path.join("data", "6m0j_prepared.pdbqt"))
    out_csv    = os.environ.get("OUTCSV", "interaction_energies_RBD_ACE2.csv")

    vdw_file   = os.path.join("data", "vdwprm")
    rsa_complex = os.path.join("results", "naccess", "6m0j_prepared.rsa")
    rsa_A       = os.path.join("results", "naccess", "6m0j_chain_A.rsa")
    rsa_E       = os.path.join("results", "naccess", "6m0j_chain_E.rsa")

    cutoff_contact = 6.0  # Å (teacher)
    cutoff_energy  = 8.0  # Å

    st = load_annotated_structure(pdb_file, pdbqt_file, vdw_file, rsa_complex, rsa_A, rsa_E)
    model = st[0]
    chainA = model["A"]
    chainE = model["E"]

    # Sanity check
    try:
        print(f"Sanity check ASA_BOUND A353 (Å^2): {model['A'][353].xtra.get('ASA_BOUND', 0.0):.2f}")
    except Exception:
        print("Sanity check ASA_BOUND A353 (Å^2): 0.00 (residue missing?)")

    # Interface residues based on cutoff_contact (more correct breaks)
    interface_A = set()
    interface_E = set()

    for resA in chainA:
        for atA in resA:
            if not atA.element.strip():
                continue
            found = False
            for resE in chainE:
                for atE in resE:
                    if not atE.element.strip():
                        continue
                    if (atA - atE) < cutoff_contact:
                        interface_A.add(resA.id[1])
                        interface_E.add(resE.id[1])
                        found = True
                        break
                if found:
                    break

    print("Interface residues A:", sorted(interface_A))
    print("Interface residues E:", sorted(interface_E))

    res_energy_A = {}
    res_energy_E = {}

    # ---- Chain A per-residue
    for resA in chainA:
        if resA.id[1] not in interface_A:
            continue

        Ee_tot = 0.0
        Ev_tot = 0.0
        for resE in chainE:
            if resE.id[1] not in interface_E:
                continue
            Ee, Ev = residue_pair_energy(resA, resE, cutoff=cutoff_energy)
            Ee_tot += Ee
            Ev_tot += Ev

        dSolv = delta_solvation_residue(resA)
        Etot = Ee_tot + Ev_tot + dSolv
        res_energy_A[resA.id[1]] = (Ee_tot, Ev_tot, dSolv, Etot)

    # ---- Chain E per-residue
    for resE in chainE:
        if resE.id[1] not in interface_E:
            continue

        Ee_tot = 0.0
        Ev_tot = 0.0
        for resA in chainA:
            if resA.id[1] not in interface_A:
                continue
            Ee, Ev = residue_pair_energy(resA, resE, cutoff=cutoff_energy)
            Ee_tot += Ee
            Ev_tot += Ev

        dSolv = delta_solvation_residue(resE)
        Etot = Ee_tot + Ev_tot + dSolv
        res_energy_E[resE.id[1]] = (Ee_tot, Ev_tot, dSolv, Etot)

    # Print tables
    print("\nResidue interaction energies for chain A (RBD):")
    print("ResID   ΔG_elec   ΔG_vdw   ΔG_solv(Δ)   ΔG_total   [kcal/mol]")
    for resid in sorted(res_energy_A):
        Ee, Ev, Es, Et = res_energy_A[resid]
        print(f"{resid:4d}  {Ee:8.3f} {Ev:8.3f} {Es:11.3f} {Et:10.3f}")

    print("\nResidue interaction energies for chain E (ACE2):")
    print("ResID   ΔG_elec   ΔG_vdw   ΔG_solv(Δ)   ΔG_total   [kcal/mol]")
    for resid in sorted(res_energy_E):
        Ee, Ev, Es, Et = res_energy_E[resid]
        print(f"{resid:4d}  {Ee:8.3f} {Ev:8.3f} {Es:11.3f} {Et:10.3f}")

    # Totals (avoid double-counting elec/vdw)
    tot_elec = sum(v[0] for v in res_energy_A.values())
    tot_vdw  = sum(v[1] for v in res_energy_A.values())
    tot_dsolv = sum(v[2] for v in res_energy_A.values()) + sum(v[2] for v in res_energy_E.values())
    tot_dG = tot_elec + tot_vdw + tot_dsolv

    print("\nTOTAL interaction free energy (interface-based):")
    print(f"ΔG_elect(A–E) = {tot_elec: .3f} kcal/mol")
    print(f"ΔG_vdw(A–E)   = {tot_vdw: .3f} kcal/mol")
    print(f"ΔG_solv       = {tot_dsolv: .3f} kcal/mol   (G_complex - G_A - G_E, interface residues)")
    print(f"ΔG_total      = {tot_dG: .3f} kcal/mol")

    # Ensure OUTCSV directory exists
    out_dir = os.path.dirname(out_csv)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # Export CSV
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Chain", "Residue",
                    "ΔG_elec (kcal/mol)", "ΔG_vdw (kcal/mol)", "ΔG_solv (kcal/mol)", "ΔG_total (kcal/mol)"])

        for resid, (Ee, Ev, Es, Et) in sorted(res_energy_A.items()):
            w.writerow(["A", resid, Ee, Ev, Es, Et])

        for resid, (Ee, Ev, Es, Et) in sorted(res_energy_E.items()):
            w.writerow(["E", resid, Ee, Ev, Es, Et])

        w.writerow([])
        w.writerow(["TOTAL(A–E)", "", tot_elec, tot_vdw, tot_dsolv, tot_dG])

    print(f"\nCSV written: {out_csv}")


if __name__ == "__main__":
    main()
