#!/usr/bin/env python3
"""
foldx_compare.py

Runs FoldX RepairPDB + AnalyseComplex for a Spike RBD–ACE2 complex and
parses:
- FoldX interaction energy between the two chains (A and E by default)
- FoldX interface residues (if present)
- A small hotspot list from FoldX interface residue energy contributions (optional heuristic)
- Optional overlap vs a user-provided hotspot list

Outputs CSVs into results/foldx/:
- foldx_summary.csv
- foldx_interface_residues.csv
- foldx_hotspots.csv
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ----------------------------
# Helpers
# ----------------------------

def die(msg: str, code: int = 1) -> None:
    raise RuntimeError(msg)

def log(msg: str) -> None:
    print(msg, flush=True)

def run_cmd(cmd: List[str], cwd: Path) -> None:
    p = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    # FoldX writes most info to stdout
    if p.returncode != 0:
        log(p.stdout)
        log(p.stderr)
        die(f"Command failed ({p.returncode}): {' '.join(cmd)}")
    # keep stdout visible (FoldX progress)
    if p.stdout.strip():
        print(p.stdout, end="")
    if p.stderr.strip():
        print(p.stderr, end="", file=sys.stderr)

def find_first(patterns: List[str], folder: Path) -> Optional[Path]:
    for pat in patterns:
        hits = sorted(folder.glob(pat))
        if hits:
            return hits[0]
    return None

def find_all(patterns: List[str], folder: Path) -> List[Path]:
    out: List[Path] = []
    for pat in patterns:
        out.extend(sorted(folder.glob(pat)))
    # unique preserving order
    seen = set()
    uniq = []
    for p in out:
        if p not in seen:
            uniq.append(p)
            seen.add(p)
    return uniq


# ----------------------------
# Parsers
# ----------------------------

_NUM_RE = re.compile(r"(-?\d+(?:\.\d+)?)")

@dataclass
class InteractionBlock:
    chain1: str
    chain2: str
    terms: Dict[str, float]
    total: float


def parse_key_value_terms(lines: List[str]) -> Dict[str, float]:
    """
    Parse lines of the form:
      Energy_VdW      =               -14.77
    Returns dict { "Energy_VdW": -14.77, ... }
    """
    terms: Dict[str, float] = {}
    for ln in lines:
        if "=" not in ln:
            continue
        left, right = ln.split("=", 1)
        key = left.strip()
        m = _NUM_RE.search(right)
        if not m:
            continue
        terms[key] = float(m.group(1))
    return terms


def parse_interaction_fxout(fxout: Path, chain_a: str = "A", chain_b: str = "E") -> InteractionBlock:
    """
    Robustly parse FoldX interaction energy from Interaction_*.fxout.

    FoldX has multiple output formats depending on version/build/options.
    This parser tries several strategies in order:

    (1) "interaction between A and E" block + the following "Total = ..."
    (2) Any line containing 'interaction' + a float number
    (3) Last occurrence of a "Total = <float>" in the whole file
    (4) Last float-like number in the file (very last resort)

    Returns InteractionBlock(chain1, chain2, terms, total).
    """
    lines = fxout.read_text(errors="ignore").splitlines()

    # --- Strategy 1: classic block header: "interaction between A and E"
    rx_header = re.compile(
        r"interaction\s+between\s+([A-Za-z0-9])\s+and\s+([A-Za-z0-9])",
        re.IGNORECASE,
    )

    matches: List[Tuple[int, str, str]] = []
    for i, ln in enumerate(lines):
        m = rx_header.search(ln)
        if m:
            matches.append((i, m.group(1), m.group(2)))

    target_pairs = {(chain_a.upper(), chain_b.upper()), (chain_b.upper(), chain_a.upper())}
    best: Optional[Tuple[int, str, str]] = None
    for i, c1, c2 in matches:
        if (c1.upper(), c2.upper()) in target_pairs:
            best = (i, c1, c2)
            break
    if best is None and matches:
        best = matches[0]

    if best is not None:
        start_idx, c1, c2 = best
        block_lines: List[str] = []
        total_val: Optional[float] = None

        for ln in lines[start_idx:]:
            block_lines.append(ln)
            if ln.strip().lower().startswith("total"):
                m = _NUM_RE.search(ln)
                if m:
                    total_val = float(m.group(1))
                    break

        if total_val is not None:
            terms = parse_key_value_terms(block_lines)
            return InteractionBlock(chain1=c1, chain2=c2, terms=terms, total=total_val)

    # --- Strategy 2: line mentions "interaction" and contains a float (many FoldX variants do this)
    rx_inter_line = re.compile(r"interaction", re.IGNORECASE)
    for ln in reversed(lines):
        if rx_inter_line.search(ln):
            m = _NUM_RE.search(ln)
            if m:
                # chains unknown in this format, but we keep your defaults
                return InteractionBlock(
                    chain1=chain_a, chain2=chain_b, terms={}, total=float(m.group(1))
                )

    # --- Strategy 3: last "Total = <float>" in the file (common in some Interaction_*.fxout)
    total_candidates: List[float] = []
    for ln in lines:
        if ln.strip().lower().startswith("total") and "=" in ln:
            m = _NUM_RE.search(ln)
            if m:
                total_candidates.append(float(m.group(1)))
    if total_candidates:
        return InteractionBlock(
            chain1=chain_a, chain2=chain_b, terms={}, total=total_candidates[-1]
        )

    # --- Strategy 4: absolute last resort: last float anywhere in file
    last_num: Optional[float] = None
    for ln in reversed(lines):
        m = _NUM_RE.search(ln)
        if m:
            last_num = float(m.group(1))
            break
    if last_num is not None:
        return InteractionBlock(chain1=chain_a, chain2=chain_b, terms={}, total=last_num)

    die(
        f"Could not extract interaction energy from {fxout.name}. "
        f"Tried: interaction header block, interaction lines, Total lines, and last numeric fallback."
    )



def parse_interface_residues_fxout(fxout: Path) -> List[str]:
    """
    Interface residues file sometimes includes residue identifiers/lines.
    We'll extract tokens like: 'A 353', 'E 505', or 'TYRE505' etc.
    Since FoldX formatting varies, we return raw unique residue-like tokens.
    """
    lines = fxout.read_text(errors="ignore").splitlines()
    residues: List[str] = []

    # Common patterns:
    # - chain + number: "A 353", "E 505"
    rx1 = re.compile(r"\b([A-Za-z0-9])\s*(\d{1,4})\b")
    # - 3-letter + chain + number: "TYRE505", "GLNE493"
    rx2 = re.compile(r"\b([A-Z]{3})([A-Za-z0-9])(\d{1,4})\b")

    for ln in lines:
        for m in rx2.finditer(ln.upper()):
            residues.append(f"{m.group(2)}{m.group(3)}")  # e.g., E505
        # rx1 can be too broad; only accept if line mentions interface/residue
        if "res" in ln.lower() or "interface" in ln.lower():
            for m in rx1.finditer(ln):
                residues.append(f"{m.group(1).upper()}{m.group(2)}")

    # unique preserve order
    seen = set()
    uniq = []
    for r in residues:
        if r not in seen:
            uniq.append(r)
            seen.add(r)
    return uniq


# ----------------------------
# Main pipeline
# ----------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True, help="Input PDB path (e.g., results/foldx/6m0j_prepared.pdb)")
    ap.add_argument("--foldx", default="results/foldx/foldx", help="FoldX executable path (default: results/foldx/foldx)")
    ap.add_argument("--workdir", default="results/foldx", help="Working directory where FoldX runs (default: results/foldx)")
    ap.add_argument("--chain-a", default="A", help="Chain ID for Spike RBD (default: A)")
    ap.add_argument("--chain-b", default="E", help="Chain ID for ACE2 (default: E)")
    ap.add_argument("--custom-hotspots", default="", help="Comma-separated residues like 'E505,E486,E493' for overlap statement")
    args = ap.parse_args()

    project_root = Path.cwd()
    workdir = (project_root / args.workdir).resolve()
    foldx_bin = (project_root / args.foldx).resolve()
    pdb_in = (project_root / args.pdb).resolve()

    if not workdir.exists():
        die(f"Workdir not found: {workdir}")
    if not foldx_bin.exists():
        die(f"FoldX executable not found: {foldx_bin}")
    if not os.access(foldx_bin, os.X_OK):
        die(f"FoldX is not executable: {foldx_bin} (run: chmod +x {foldx_bin})")
    if not pdb_in.exists():
        die(f"Input PDB not found: {pdb_in}")

    # Ensure PDB is in workdir (FoldX expects local file names often)
    pdb_local = workdir / pdb_in.name
    if pdb_in != pdb_local:
        shutil.copy2(pdb_in, pdb_local)
        log(f"[OK] Copied PDB into workdir: {pdb_local.name}")
    else:
        log(f"[OK] Using PDB in workdir: {pdb_local.name}")

    # 1) RepairPDB
    log("[RUN] FoldX RepairPDB ...")
    run_cmd([str(foldx_bin), "--command", "RepairPDB", "--pdb", pdb_local.name, "--output-dir", "."], cwd=workdir)
    log("[OK] RepairPDB finished")

    # Find repaired PDB (FoldX may produce Repair_<pdb>.pdb OR <pdb>_Repair.pdb)
    repaired_pdb = find_first(
        patterns=[
            f"Repair_{pdb_local.stem}.pdb",
            f"{pdb_local.stem}_Repair.pdb",
            "Repair_*.pdb",
            "*_Repair.pdb",
        ],
        folder=workdir,
    )
    if not repaired_pdb or not repaired_pdb.exists():
        die("RepairPDB did not produce a repaired PDB (Repair_*.pdb or *_Repair.pdb). Check FoldX output files.")

    log(f"[OK] Repaired PDB: {repaired_pdb.name}")

    # 2) AnalyseComplex
    log("[RUN] FoldX AnalyseComplex ...")
    run_cmd([str(foldx_bin), "--command", "AnalyseComplex", "--pdb", repaired_pdb.name, "--output-dir", "."], cwd=workdir)
    log("[OK] AnalyseComplex finished")

    # Find interaction fxout (naming varies)
    interaction_fx = find_first(
        patterns=[
            "Interaction_*_AC.fxout",
            "Interaction_*.fxout",
        ],
        folder=workdir,
    )
    if not interaction_fx:
        # sometimes files are generated but with different prefix; list candidates
        fxouts = sorted(workdir.glob("*.fxout"))
        die(f"AnalyseComplex output not found. Existing .fxout files: {[p.name for p in fxouts]}")

    # Parse interaction block robustly
    inter = parse_interaction_fxout(interaction_fx, chain_a=args.chain_a, chain_b=args.chain_b)

    # Optional: interface residues fxout
    interface_fx = find_first(
        patterns=[
            "Interface_Residues_*_AC.fxout",
            "Interface_Residues_*.fxout",
        ],
        folder=workdir,
    )
    interface_res = parse_interface_residues_fxout(interface_fx) if interface_fx else []

    # Build simple hotspot list heuristic (top residues can also be computed elsewhere;
    # here we just pass interface residue tokens if present)
    # You can later replace this with a better rule if needed.
    foldx_hotspots = interface_res[:10]  # keep small

    # Write summary CSV
    summary_csv = workdir / "foldx_summary.csv"
    with summary_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["interaction_chains", f"{inter.chain1}-{inter.chain2}"])
        w.writerow(["interaction_total_kcal_per_mol", f"{inter.total:.3f}"])
        w.writerow(["interaction_fxout", interaction_fx.name])
        w.writerow(["repaired_pdb", repaired_pdb.name])

    # Write interface residues CSV
    interface_csv = workdir / "foldx_interface_residues.csv"
    with interface_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["residue"])
        for r in interface_res:
            w.writerow([r])

    # Write hotspots CSV (heuristic)
    hotspots_csv = workdir / "foldx_hotspots.csv"
    with hotspots_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["residue"])
        for r in foldx_hotspots:
            w.writerow([r])

    log(f"[RESULT] FoldX interaction energy ({inter.chain1}-{inter.chain2}) Total = {inter.total:.2f} kcal/mol")
    log(f"[OK] Wrote: {summary_csv.relative_to(project_root)}")
    log(f"[OK] Wrote: {interface_csv.relative_to(project_root)}")
    log(f"[OK] Wrote: {hotspots_csv.relative_to(project_root)}")

    # Optional overlap statement
    if args.custom_hotspots.strip():
        custom = [x.strip().upper() for x in args.custom_hotspots.split(",") if x.strip()]
        foldx_set = set([x.upper() for x in foldx_hotspots])
        custom_set = set(custom)
        overlap = sorted(list(foldx_set.intersection(custom_set)))
        log("[OVERLAP] Custom hotspots vs FoldX (heuristic) hotspots")
        log(f"  Custom: {custom}")
        log(f"  FoldX : {foldx_hotspots}")
        log(f"  Overlap: {overlap if overlap else 'None'}")


if __name__ == "__main__":
    main()
