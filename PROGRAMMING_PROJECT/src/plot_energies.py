import pandas as pd
import matplotlib.pyplot as plt
import os

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
results_dir = os.path.join(base_dir, 'results')
csv_path = os.path.join(results_dir, "interaction_energies_RBD_ACE2.csv")

# Load CSV
df = pd.read_csv(csv_path)

# Column name EXACTLY as in CSV
TOTAL_COL = "ΔG_total (kcal/mol)"

# Separate chains
df_A = df[df["Chain"] == "A"].copy()
df_E = df[df["Chain"] == "E"].copy()

# Sort by total interaction energy
df_A.sort_values(TOTAL_COL, inplace=True)
df_E.sort_values(TOTAL_COL, inplace=True)

# ---------------------------
# Plot RBD (Chain A)
# ---------------------------
plt.figure(figsize=(12,4))
plt.bar(
    df_A["Residue"].astype(str),
    df_A[TOTAL_COL]
)
plt.axhline(0, color="black", linewidth=0.8)
plt.title("RBD (Chain A) – Residue Interaction Free Energies")
plt.ylabel("ΔG_total (kcal/mol)")
plt.xlabel("Residue")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("RBD_interface_energies.png", dpi=300)
plt.close()

# ---------------------------
# Plot ACE2 (Chain E)
# ---------------------------
plt.figure(figsize=(12,4))
plt.bar(
    df_E["Residue"].astype(str),
    df_E[TOTAL_COL]
)
plt.axhline(0, color="black", linewidth=0.8)
plt.title("ACE2 (Chain E) – Residue Interaction Free Energies")
plt.ylabel("ΔG_total (kcal/mol)")
plt.xlabel("Residue")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("ACE2_interface_energies.png", dpi=300)
plt.close()

print("Plots saved: RBD_interface_energies.png, ACE2_interface_energies.png")

# ---------------------------
# Top stabilizers and destabilizers
# ---------------------------
print("\nTop 5 stabilizing residues – RBD (most negative ΔG):")
print(df_A.nsmallest(5, TOTAL_COL))

print("\nTop 5 destabilizing residues – RBD (most positive ΔG):")
print(df_A.nlargest(5, TOTAL_COL))

print("\nTop 5 stabilizing residues – ACE2:")
print(df_E.nsmallest(5, TOTAL_COL))

print("\nTop 5 destabilizing residues – ACE2:")
print(df_E.nlargest(5, TOTAL_COL))

# Optional: export per-chain CSVs
df_A.to_csv(os.path.join(results_dir, "RBD_interface_energies.csv"), index=False)
df_E.to_csv(os.path.join(results_dir, "ACE2_interface_energies.csv"), index=False)
