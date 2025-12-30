#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

CSV_IN  = "results/alanine_scanning/alanine_ddg.csv"
PNG_OUT = "results/alanine_scanning/alanine_ddg_pub_v2.png"

df = pd.read_csv(CSV_IN)

# Ensure numeric
df["Residue"] = pd.to_numeric(df["Residue"])
df["ΔΔG (kcal/mol)"] = pd.to_numeric(df["ΔΔG (kcal/mol)"])

# IMPORTANT: Chain A = ACE2, Chain E = RBD
dfA = df[df["Chain"] == "A"].sort_values("Residue")  # ACE2
dfE = df[df["Chain"] == "E"].sort_values("Residue")  # RBD

# Use categorical x positions so residues are evenly spaced (no huge gaps)
dfA = dfA.reset_index(drop=True)
dfE = dfE.reset_index(drop=True)
dfA["xpos"] = range(len(dfA))
dfE["xpos"] = range(len(dfE))

# Hotspots (destabilizing)
hot = df[df["ΔΔG (kcal/mol)"] > 1.0].sort_values("ΔΔG (kcal/mol)", ascending=False)

fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharey=True, constrained_layout=True)

# ---- ACE2 (Chain A) ----
axes[0].bar(dfA["xpos"], dfA["ΔΔG (kcal/mol)"], width=0.85)
axes[0].axhline(0, lw=0.8, color="black")
axes[0].set_title("Alanine scanning ΔΔG at the RBD–ACE2 interface")
axes[0].set_ylabel("ΔΔG (kcal/mol)")
axes[0].text(0.01, 0.90, "ACE2 (Chain A)", transform=axes[0].transAxes)

axes[0].set_xticks(dfA["xpos"])
axes[0].set_xticklabels(dfA["Residue"].astype(int), rotation=90, fontsize=7)

# ---- RBD (Chain E) ----
axes[1].bar(dfE["xpos"], dfE["ΔΔG (kcal/mol)"], width=0.85)
axes[1].axhline(0, lw=0.8, color="black")
axes[1].set_xlabel("Interface residue (sequence number)")
axes[1].set_ylabel("ΔΔG (kcal/mol)")
axes[1].text(0.01, 0.90, "Spike RBD (Chain E)", transform=axes[1].transAxes)

axes[1].set_xticks(dfE["xpos"])
axes[1].set_xticklabels(dfE["Residue"].astype(int), rotation=90, fontsize=7)

# Annotate hotspots (top 6)
hot_top = hot.head(6)
for _, r in hot_top.iterrows():
    chain = r["Chain"]
    resid = int(r["Residue"])
    y = float(r["ΔΔG (kcal/mol)"])

    if chain == "A":
        row = dfA.index[dfA["Residue"] == resid]
        if len(row) == 0:
            continue
        x = int(dfA.loc[row[0], "xpos"])
        ax = axes[0]
        label = f"A{resid}"
    else:
        row = dfE.index[dfE["Residue"] == resid]
        if len(row) == 0:
            continue
        x = int(dfE.loc[row[0], "xpos"])
        ax = axes[1]
        label = f"E{resid}"

    ax.annotate(label, xy=(x, y), xytext=(x, y + 0.25),
                ha="center", fontsize=8)

# y-lims padded
ymin = df["ΔΔG (kcal/mol)"].min()
ymax = df["ΔΔG (kcal/mol)"].max()
pad = 0.5
axes[0].set_ylim(ymin - pad, ymax + pad)

plt.savefig(PNG_OUT, dpi=600)
print(f"Saved: {PNG_OUT}")
plt.show()
