# fix_pdbqt_strict.py

import re

input_file = "6m0j_prepared.pdbqt"
output_file = "6m0j_fixed2.pdbqt"

float_re = re.compile(r'^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$')

with open(input_file) as f, open(output_file, "w") as out:

    for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):

            # Extract full charge string
            raw_charge = line[69:76].strip()

            # If charge isn't a valid float → replace
            if not float_re.match(raw_charge):
                raw_charge = "0.000"

            # Rebuild line with fixed charge in the correct column
            charge_field = f"{float(raw_charge):7.3f}"
            line = line[:69] + charge_field + line[76:]

        out.write(line)

print("✓ Strict fixed PDBQT saved as:", output_file)
