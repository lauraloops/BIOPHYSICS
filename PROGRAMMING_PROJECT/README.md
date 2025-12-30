# SARS-CoV-2 RBD-ACE2 Interaction Energy Analysis

This project analyzes the interaction energies between the SARS-CoV-2 Spike RBD and human ACE2 receptor, including the effect of variants and alanine scanning.

## Project Structure

- **data/**: Input data files (PDB structures, parameters).
- **src/**: Core Python modules and helper scripts.
- **scripts/**: Executable scripts for running analyses (e.g., variant scanning, alanine scanning).
- **results/**: Output files (energies, structures, logs).
- **pymol/**: PyMOL scripts for visualization and mutagenesis.
- **foldx/**: FoldX plugin and configuration files.
- **BioPhysics/**: (Excluded from version control) Course materials and exercises.

## Requirements

- Python 3
- Biopython
- NumPy
- OpenBabel (command line tool `obabel`)
- PyMOL (command line tool `pymol`)

Install Python dependencies:
```bash
pip install -r requirements.txt
```

## Usage

Run scripts from the project root directory.

### Variant Scanning
To calculate interaction energies for specific variants:
```bash
python scripts/variant_scan.py
```

### Alanine Scanning
To perform alanine scanning on the interface:
```bash
python scripts/alanine_scan.py
```

## Notes
- Ensure `obabel` and `pymol` are in your system PATH.
- The `BioPhysics` folder contains supplementary material and is not part of the main analysis pipeline.
