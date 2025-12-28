load results/alanine_scanning/6m0j_ddg_bfactor.pdb, complex

hide everything
show cartoon, complex
color gray80, complex

# Color by ΔΔG stored in b-factor
# Red = higher ΔΔG (destabilizing), Blue = lower ΔΔG (stabilizing)
spectrum b, blue_white_red, complex, minimum=-10, maximum=3

# Highlight interface chains
color marine, chain A
color orange, chain E

# Show hotspots as sticks
select hotspots, (chain E and resi 456+486+493+505)
show sticks, hotspots
set stick_radius, 0.2
color red, hotspots

# Nice view and render
bg_color white
set ray_opaque_background, off
orient chain A or chain E
ray 2000,1500
png results/alanine_scanning/ddg_pymol.png, dpi=300
