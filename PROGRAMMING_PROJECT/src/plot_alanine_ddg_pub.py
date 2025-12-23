#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

CSV_IN = "results/alanine_scanning/alanine_ddg.csv"
PNG_OUT = "results/alanine_scanning/alanine_ddg_pub.png"

df = pd.read_csv(CSV_IN)

# Ensure numeric
df["Residue"] = pd.to_numeric(df["Residue"])
df["ΔΔG (kcal/mol)"] = pd.to_numeric(df["ΔΔG (kcal/mol)"])

dfA = df[df["Chain"] == "A"].sort_values("Residue")
dfE = df[df["Chain"] == "E"].sort_values("Residue")

# Hotspots: destabilizing mutations
hot = df[df["ΔΔG (kcal/mol)"] > 1.0].sort_values("ΔΔG (kcal/mol)", ascending=False)

fig, axes = plt.subplots(
    2, 1, figsize=(11, 6), sharey=True, constrained_layout=True
)

# --- Chain A (RBD) ---
axes[0].bar(dfA["Residue"], dfA["ΔΔG (kcal/mol)"], width=1.0)
axes[0].axhline(0, lw=0.8, color="black")
axes[0].set_title("Alanine scanning ΔΔG at the RBD–ACE2 interface")
axes[0].set_ylabel("ΔΔG (kcal/mol)")
axes[0].text(0.01, 0.90, "RBD (Chain A)", transform=axes[0].transAxes)

# --- Chain E (ACE2) ---
axes[1].bar(dfE["Residue"], dfE["ΔΔG (kcal/mol)"], width=1.0)
axes[1].axhline(0, lw=0.8, color="black")
axes[1].set_xlabel("Residue number")
axes[1].set_ylabel("ΔΔG (kcal/mol)")
axes[1].text(0.01, 0.90, "ACE2 (Chain E)", transform=axes[1].transAxes)

# Annotate hotspots on their panel
for _, r in hot.iterrows():
    ax = axes[0] if r["Chain"] == "A" else axes[1]
    x = int(r["Residue"])
    y = float(r["ΔΔG (kcal/mol)"])
    ax.annotate(f"{r['Chain']}{x}",
                xy=(x, y),
                xytext=(x, y + 0.2),
                ha="center",
                fontsize=8,
                rotation=0)

# Nice y-limits (auto but padded)
ymin = df["ΔΔG (kcal/mol)"].min()
ymax = df["ΔΔG (kcal/mol)"].max()
pad = 0.5
axes[0].set_ylim(ymin - pad, ymax + pad)

plt.savefig(PNG_OUT, dpi=600)
print(f"Saved: {PNG_OUT}")
plt.show()
