reinitialize
fetch 6M0J, async=0

remove solvent
remove hetatm
remove not chain A+E

hide everything
show cartoon, chain A+E

color cyan, chain A
color pink, chain E

# --- Distance tests ---
# 4 Å
select test4,  (chain A within 4 of chain E)
show sticks, test4
color red, test4

select test4b, (chain E within 4 of chain A)
show sticks, test4b
color hotpink, test4b

# 5 Å
select test5,  (chain A within 5 of chain E)
show sticks, test5
color red, test5

select test5b, (chain E within 5 of chain A)
show sticks, test5b
color hotpink, test5b

# 6 Å  (your chosen “optimal”)
select test6,  (chain A within 6 of chain E)
show sticks, test6
color red, test6

select test6b, (chain E within 6 of chain A)
show sticks, test6b
color blue, test6b

# Optional: 7 Å “add 1–2 Å margin”
select test7,  (chain A within 7 of chain E)
show sticks, test7
color red, test7

select test7b, (chain E within 7 of chain A)
show sticks, test7b
color blue, test7b

# render
set ray_opaque_background, off
bg_color white
ray 2400,1800
png myimage.png, dpi=300
