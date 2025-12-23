# fix_pdbqt.py

import re

input_file = "6m0j_prepared.pdbqt"
output_file = "6m0j_fixed.pdbqt"

float_re = re.compile(r'^[-+]?[0-9]*\.?[0-9]+$')  # regex to detect floats

with open(input_file) as f, open(output_file, "w") as out:
    for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Extract charge (columns 70–76)
            raw_charge = line[69:76].strip()

            # If charge is not a valid float, replace with 0.0
            if not float_re.match(raw_charge):
                fixed_charge = " 0.000"
                line = line[:69] + fixed_charge + line[76:]

        out.write(line)

print("✓ Fixed PDBQT saved as:", output_file)
