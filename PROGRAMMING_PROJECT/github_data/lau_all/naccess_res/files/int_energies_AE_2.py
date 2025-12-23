from Bio.PDB import PDBParser
import os
from forcefield import VdwParamset
import math

# ---------------------------
# 1. Load structure and charges
# ---------------------------
def load_annotated_structure(pdb_file, pdbqt_file, vdw_file):
    ff_params = VdwParamset(vdw_file)
    parser = PDBParser(PERMISSIVE=1)
    st = parser.get_structure("STR", pdb_file)

    # Parse PDBQT charges
    params = [{}]
    with open(pdbqt_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("ATOM") or line.startswith("HETATM"):
                params.append({"raw": line})

    # Assign charges and vdW parameters
    i = 1
    total_charge = 0.0
    for at in st.get_atoms():
        at.serial_number = i
        i += 1
        raw = params[at.serial_number]["raw"]
        parts = raw.split()
        charge_str = parts[-2]
        try:
            charge = float(charge_str)
        except Exception:
            charge = 0.0
        atom_type = parts[-1]
        if atom_type not in ff_params.at_types:
            atom_type = "C"
        at.xtra["charge"] = charge
        at.xtra["atom_type"] = atom_type
        at.xtra["vdw"] = ff_params.at_types[atom_type]
        total_charge += charge

    print(f"Total charge (recomputed): {total_charge:8.2f}")
    return st

# ---------------------------
# 2. Energy functions
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

def residue_pair_energy(resA, resB, cutoff=8.0):
    E_elec = 0.0
    E_vdw = 0.0
    for atA in resA:
        if atA.element.strip() == "H":
            continue
        for atB in resB:
            if atB.element.strip() == "H":
                continue
            r = atA - atB
            if r > cutoff:
                continue
            E_elec += coulomb_energy(atA, atB)
            E_vdw += lennard_jones_energy(atA, atB)
    return E_elec, E_vdw

# ---------------------------
# 3. Solvation energies from NACCESS
# ---------------------------
def read_naccess_asa(naccess_file):
    """
    Reads ASA from NACCESS .asa file (per atom) using fixed-width columns.
    Returns dict: key = (chain, resnum, atom_name), value = ASA
    """
    asa_dict = {}
    with open(naccess_file) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                res_name  = line[17:20].strip()
                chain     = line[21].strip()
                res_num   = int(line[22:26].strip())
                asa       = float(line[54:60].strip())  # absolute ASA
                asa_dict[(chain, res_num, atom_name)] = asa
    return asa_dict


def compute_solvation_energy(st, asa_dict, coeff=0.005):
    """
    Compute solvation energy per residue:
    E_solv = sum_over_atoms(ASA_atom * coeff)
    """
    res_solv = {}
    for res in st.get_residues():
        chain = res.get_parent().id
        resnum = res.id[1]
        total_solv = 0.0
        for at in res:
            if at.element.strip() == "H":
                continue
            key = (chain, resnum, at.name)
            asa = asa_dict.get(key, 0.0)
            total_solv += asa * coeff
        res_solv[resnum] = total_solv
    return res_solv

# ---------------------------
# 4. Main: compute total ΔG including solvation
# ---------------------------
def main():
    pdb_file = "6m0j_prepared.pdb"
    pdbqt_file = "6m0j_prepared.pdbqt"
    vdw_file = "vdwprm"
    
    # NACCESS files
    naccess_chainA = "6m0j_chain_A.asa"
    naccess_chainE = "6m0j_chain_E.asa"
    naccess_complex = "6m0j_prepared.asa"
    
    st = load_annotated_structure(pdb_file, pdbqt_file, vdw_file)
    model = st[0]
    chainA = model["A"]
    chainE = model["E"]
    
    # Define interface residues
    cutoff_contact = 5.0
    interface_A = set()
    interface_E = set()
    for resA in chainA:
        for atA in resA:
            if atA.element == "H":
                continue
            for resE in chainE:
                for atE in resE:
                    if atE.element == "H":
                        continue
                    if atA - atE < cutoff_contact:
                        interface_A.add(resA.id[1])
                        interface_E.add(resE.id[1])
    
    print("Interface residues A:", sorted(interface_A))
    print("Interface residues E:", sorted(interface_E))
    
    # Per-residue electrostatic and vdW energies
    res_energy_A = {}
    res_energy_E = {}
    cutoff_energy = 8.0
    
    for resA in chainA:
        if resA.id[1] not in interface_A:
            continue
        E_elec_tot = 0.0
        E_vdw_tot = 0.0
        for resE in chainE:
            if resE.id[1] not in interface_E:
                continue
            E_elec, E_vdw = residue_pair_energy(resA, resE, cutoff=cutoff_energy)
            E_elec_tot += E_elec
            E_vdw_tot += E_vdw
        res_energy_A[resA.id[1]] = (E_elec_tot, E_vdw_tot, E_elec_tot + E_vdw_tot)
    
    for resE in chainE:
        if resE.id[1] not in interface_E:
            continue
        E_elec_tot = 0.0
        E_vdw_tot = 0.0
        for resA in chainA:
            if resA.id[1] not in interface_A:
                continue
            E_elec, E_vdw = residue_pair_energy(resA, resE, cutoff=cutoff_energy)
            E_elec_tot += E_elec
            E_vdw_tot += E_vdw
        res_energy_E[resE.id[1]] = (E_elec_tot, E_vdw_tot, E_elec_tot + E_vdw_tot)
    
    # Solvation energies
    asa_A = read_naccess_asa(naccess_chainA)
    asa_E = read_naccess_asa(naccess_chainE)
    asa_complex = read_naccess_asa(naccess_complex)
    
    solv_A = compute_solvation_energy(chainA, asa_A)
    solv_E = compute_solvation_energy(chainE, asa_E)
    solv_complex = compute_solvation_energy(model, asa_complex)
    
    # Compute ΔG interaction including solvation
    print("\nResidue interaction energies with solvation for chain A (RBD):")
    print("ResID   E_elec    E_vdw    E_solv   E_total [kcal/mol]")
    for resnum in sorted(interface_A):
        Ee, Ev, Et = res_energy_A[resnum]
        Es = solv_complex.get(resnum,0) - solv_A.get(resnum,0) - 0  # contribution from E neglected per-residue
        E_total = Et + Es
        print(f"{resnum:4d}  {Ee:8.3f} {Ev:8.3f} {Es:8.3f} {E_total:8.3f}")
    
    print("\nResidue interaction energies with solvation for chain E (ACE2):")
    print("ResID   E_elec    E_vdw    E_solv   E_total [kcal/mol]")
    for resnum in sorted(interface_E):
        Ee, Ev, Et = res_energy_E[resnum]
        Es = solv_complex.get(resnum,0) - 0 - solv_E.get(resnum,0)  # contribution from A neglected per-residue
        E_total = Et + Es
        print(f"{resnum:4d}  {Ee:8.3f} {Ev:8.3f} {Es:8.3f} {E_total:8.3f}")

if __name__ == "__main__":
    main()
