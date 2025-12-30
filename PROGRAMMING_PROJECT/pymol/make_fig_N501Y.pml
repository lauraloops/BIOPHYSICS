# ==============================
# WT vs N501Y (Alpha) figure
# ==============================

reinitialize
set orthoscopic, on
set cartoon_fancy_helices, on
set cartoon_smooth_loops, on
set cartoon_sampling, 14
bg_color white
set ray_opaque_background, off
set antialias, 2
set depth_cue, 0
set stick_radius, 0.18
set sphere_scale, 0.25
set label_size, 20
set label_color, black
set ray_trace_mode, 1
set ray_shadows, 0

# ---- Paths (EDIT IF NEEDED) ----
load /home/laura/Desktop/BIOPHYSICS/PROGRAMMING_PROJECT/data/6m0j_prepared.pdb, wt
load /home/laura/Desktop/BIOPHYSICS/PROGRAMMING_PROJECT/results/variants/mutants/Alpha_E501_TYR.pdb, mut

# ---- Clean up ----
remove (wt and hetatm)
remove (mut and hetatm)

# ---- Define chains (your system: A=ACE2, E=RBD) ----
select wt_ACE2, wt and chain A
select wt_RBD,  wt and chain E
select mut_ACE2, mut and chain A
select mut_RBD,  mut and chain E

# ---- Global representations ----
hide everything, all
show cartoon, wt_ACE2 or wt_RBD
show cartoon, mut_ACE2 or mut_RBD

# Colors (keep consistent & legible)
color cyan, wt_ACE2
color lightpink, wt_RBD
color marine, mut_ACE2
color salmon, mut_RBD

# ---- Align mutant to WT (important for comparison) ----
align mut and chain A+E, wt and chain A+E

# ---- Interface region around residue 501 ----
# Residue 501 in RBD (chain E)
select wt_501,  wt and chain E and resi 501
select mut_501, mut and chain E and resi 501

# Nearby ACE2 residues around 501 (within 6 Å)
select wt_nearACE2_501,  wt and chain A and byres (wt and chain A within 6 of wt_501)
select mut_nearACE2_501, mut and chain A and byres (mut and chain A within 6 of mut_501)

# Nearby RBD residues around 501 (within 6 Å)
select wt_nearRBD_501,  wt and chain E and byres (wt and chain E within 6 of wt_501)
select mut_nearRBD_501, mut and chain E and byres (mut and chain E within 6 of mut_501)

# Show close neighborhood in sticks
show sticks, wt_nearACE2_501 or wt_nearRBD_501
show sticks, mut_nearACE2_501 or mut_nearRBD_501

# Highlight residue 501 strongly
color red, mut_501
color yelloworange, wt_501
show sticks, wt_501 or mut_501

# ---- “Crowding / steric” hint with spheres ----
# Show spheres for heavy atoms of the mutated residue + nearest ACE2 atoms
select mut501_heavy, mut_501 and not name H*
select mutACE2_contact, (mut and chain A) within 4.0 of mut501_heavy

show spheres, mut501_heavy
show spheres, mutACE2_contact
set sphere_transparency, 0.55, mut501_heavy
set sphere_transparency, 0.70, mutACE2_contact
color red, mut501_heavy
color gray70, mutACE2_contact

# ---- Labels (minimal, readable) ----
label wt_501 and name CA, "WT N501"
label mut_501 and name CA, "Mut Y501"

# Optional: label closest ACE2 residue(s) for context (pick one or two)
# (You can adjust these after you see which ones are near)
select mutACE2_close_res, byres (mutACE2_contact)
label mutACE2_close_res and name CA, "%s%s" % (resn, resi)

# ---- Views / Scenes ----
# Overview
orient wt_ACE2 or wt_RBD
zoom (wt_ACE2 or wt_RBD), 10
scene overview, store

# Close-up on 501 interface
zoom (wt_501 or mut_501 or mutACE2_contact), 10
scene closeup, store

# ---- Output images ----
# Make a nice close-up PNG for the report
scene closeup, recall
ray 2200, 1600
png /home/laura/Desktop/BIOPHYSICS/PROGRAMMING_PROJECT/results/figures/N501Y_closeup.png, dpi=300

# Optional overview image
scene overview, recall
ray 2200, 1600
png /home/laura/Desktop/BIOPHYSICS/PROGRAMMING_PROJECT/results/figures/N501Y_overview.png, dpi=300
