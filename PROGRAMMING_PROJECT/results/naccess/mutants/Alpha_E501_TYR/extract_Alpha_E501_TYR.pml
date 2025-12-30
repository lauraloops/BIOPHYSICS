load /home/laura/Desktop/BIOPHYSICS/PROGRAMMING_PROJECT/results/variants/mutants/Alpha_E501_TYR.pdb, m

# sanity: show we loaded atoms (prints to stdout)
python
from pymol import cmd
print("ATOM_COUNT:", cmd.count_atoms("m"))
python end

save /home/laura/Desktop/BIOPHYSICS/PROGRAMMING_PROJECT/results/naccess/mutants/Alpha_E501_TYR/Alpha_E501_TYR_complex.pdb, m
save /home/laura/Desktop/BIOPHYSICS/PROGRAMMING_PROJECT/results/naccess/mutants/Alpha_E501_TYR/Alpha_E501_TYR_chain_A.pdb, m and chain A
save /home/laura/Desktop/BIOPHYSICS/PROGRAMMING_PROJECT/results/naccess/mutants/Alpha_E501_TYR/Alpha_E501_TYR_chain_E.pdb, m and chain E

quit
