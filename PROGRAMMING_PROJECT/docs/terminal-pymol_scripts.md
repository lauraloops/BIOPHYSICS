# Terminal & PyMOL Scripts — Spike RBD–ACE2 (6M0J) Interface + Energy Workflow

This document consolidates **all terminal commands** and **PyMOL scripts** used in the project, in a single, readable format.

---

## Report structure (for reference)

- **< 20 pages**
- Introduction (≈2 pages)  
- Methods  
- Results  
- Discussion  
- Conclusions  
- Bibliography  
- Appendix  

---

# STEP 1 — Structure preparation (visual + cleaning + preparation)

## 1.1 Visual inspection in PyMOL (distance selection)

**Goal:** visually inspect interface contacts and choose a distance that includes direct contact residues.  
We tested **4 Å, 5 Å, 6 Å**, and decided that **6 Å** is the best compromise (direct contacts + adjacent residues).  
As requested, we then considered **6 Å + (1–2 Å) → 7 Å**.

### PyMOL script (interface distance tests)

```pymol
# STEP 1: Visual inspection in PyMOL

# Load structure (use your actual file name)
# load 6m0j_raw.pdb

# Clean for visualization
remove solvent
remove hetatm

# Color chains
color cyan, chain A
color pink, chain E

# ---- 4 Å ----
select test4, (chain A within 4 of chain E)
show sticks, test4
color red, test4

select test4b, (chain E within 4 of chain A)
show sticks, test4b
color hotpink, test4b

# ---- 5 Å ----
select test5, (chain A within 5 of chain E)
show sticks, test5
color red, test5

select test5b, (chain E within 5 of chain A)
show sticks, test5b
color hotpink, test5b

# ---- 6 Å ----
select test6, (chain A within 6 of chain E)
show sticks, test6
color red, test6

select test6b, (chain E within 6 of chain A)
show sticks, test6b
color blue, test6b

# Decision: optimal distance ~6 Å
# Add 1–2 Å => 7 Å
select test7, (chain A within 7 of chain E)
show sticks, test7
color red, test7

select test7b, (chain E within 7 of chain A)
show sticks, test7b
color blue, test7b

# Render
set ray_opaque_background, off
ray
png myimage.png, dpi=300
