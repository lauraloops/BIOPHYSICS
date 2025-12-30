#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt

CSV_IN  = "results/alanine_scanning/alanine_ddg.csv"
PNG_OUT = "results/alanine_scanning/alanine_ddg_overlay_v2.png"

df = pd.read_csv(CSV_IN)

df["Residue"] = pd.to_numeric(df["Residue"])
df["ΔΔG (kcal/mol)"] = pd.to_numeric(df["ΔΔG (kcal/mol)"])

# In your project:
# Chain A = ACE2, Chain E = Spike RBD
ace2 = df[df["Chain"] == "A"].copy().sort_values("Residue")
rbd  = df[df["Chain"] == "E"].copy().sort_values("Residue")

# --- Build a "stitched" x-axis ---
# Put RBD right after ACE2 with a small spacer
spacer = 15
ace2_min = ace2["Residue"].min()
ace2_max = ace2["Residue"].max()
rbd_min  = rbd["Residue"].min()

# ACE2 stays as-is
ace2["x"] = ace2["Residue"]

# RBD shifted so its first residue starts after ACE2 max + spacer
shift = (ace2_max + spacer) - rbd_min
rbd["x"] = rbd["Residue"] + shift

# --- Plot ---
plt.figure(figsize=(12, 4.2))
plt.axhline(0, lw=0.8, color="black")

# Use bars with narrow width (like your original)
plt.bar(ace2["x"], ace2["ΔΔG (kcal/mol)"], width=1.0, label="ACE2 (Chain A)")
plt.bar(rbd["x"],  rbd["ΔΔG (kcal/mol)"], width=1.0, label="Spike RBD (Chain E)")

plt.title("Alanine scanning ΔΔG for the RBD–ACE2 interface")
plt.ylabel("ΔΔG (kcal/mol)")
plt.xlabel("Residue number (stitched ACE2 then RBD)")

# --- Add tick labels: show real residue numbers in both regions ---
# Keep ticks sparse so it stays readable
ticks = []

# ACE2 ticks every ~25 residues
for t in range(int(ace2_min), int(ace2_max)+1, 25):
    ticks.append((t, str(t)))

# RBD ticks every ~10 residues (still sparse)
for t in range(int(rbd_min), int(rbd["Residue"].max())+1, 10):
    ticks.append((t + shift, str(t)))

xticks_pos = [p for p, lab in ticks]
xticks_lab = [lab for p, lab in ticks]
plt.xticks(xticks_pos, xticks_lab, rotation=0)

# Mark the split visually
plt.axvline(ace2_max + spacer/2, lw=0.8, color="gray")

# Hotspots label (ΔΔG > 1)
hot = df[df["ΔΔG (kcal/mol)"] > 1.0].copy()
for _, r in hot.iterrows():
    resid = int(r["Residue"])
    y = float(r["ΔΔG (kcal/mol)"])
    if r["Chain"] == "A":
        x = resid
        label = f"A{resid}"
    else:
        x = resid + shift
        label = f"E{resid}"
    plt.annotate(label, (x, y), xytext=(x, y+0.25), ha="center", fontsize=8)

plt.legend(frameon=True)

# y padding
ymin = df["ΔΔG (kcal/mol)"].min()
ymax = df["ΔΔG (kcal/mol)"].max()
plt.ylim(ymin - 0.5, ymax + 0.8)

plt.tight_layout()
plt.savefig(PNG_OUT, dpi=600)
print(f"Saved: {PNG_OUT}")
plt.show()
