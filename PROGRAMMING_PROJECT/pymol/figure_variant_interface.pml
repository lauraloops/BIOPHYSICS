bg_color white
set ray_opaque_background, off
set cartoon_fancy_helices, on
set cartoon_smooth_loops, on
set cartoon_transparency, 0.0

# Load structure
load Alpha_E501_TYR.pdb, complex

# Define chains
select ACE2, chain A
select RBD, chain E

# Cartoon representation
show cartoon, ACE2
show cartoon, RBD
color marine, ACE2
color salmon, RBD

# Interface residues (from WT definition)
select interface_A, chain A and resi 19+20+23+24+25+27+28+29+30+31+34+35+37+38+41+42+45+79+82+83+325+326+330+352+353+354+355+357+386+393
select interface_E, chain E and resi 403+417+445+446+449+453+455+456+475+485+486+487+489+493+496+498+499+500+501+503+505

# Show interface sticks
show sticks, interface_A or interface_E
color black, interface_A or interface_E
set stick_radius, 0.2

# Highlight variant residue N501Y
select variant_501, chain E and resi 501
show sticks, variant_501
color red, variant_501
set stick_radius, 0.3

# Improve rendering
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 0

# Camera
orient interface_E
zoom interface_E, 12

# High-quality render
ray 2400,1800
png Figure_N501Y_interface.png, dpi=300
