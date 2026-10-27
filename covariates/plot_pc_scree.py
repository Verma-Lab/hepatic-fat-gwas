#!/usr/bin/env python3
"""
Scree plots of PMBB Release 2026-4.0 imputed genetic PC eigenvalues (combined,
EUR, AFR), to justify how many PCs are used as GWAS covariates.

Run with:
    python3 plot_pc_scree.py
"""
import glob
import sys
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# SETTINGS
# ============================================================================
PCA_DIR = "/static/PMBB/PMBB-Release-2026-4.0/Imputed/PCA"
PATTERNS = {
    "combined": f"{PCA_DIR}/combined/*no1KG.eigenvalues.tsv",
    "EUR": f"{PCA_DIR}/by_ancestry/*eigenvalues.EUR.tsv",
    "AFR": f"{PCA_DIR}/by_ancestry/*eigenvalues.AFR.tsv",
}
N_PCS = {"combined": 5, "EUR": 5, "AFR": 5}   # cutoffs used in the GWAS (decided 2026-10-12)
N_EXPECTED = 20
OUTPUT_DIR = Path("/project/verma_shared/projects/Liver_IDPs/covariates/pc_scree_plots")

USED_COLOR = "#2a6fb0"
UNUSED_COLOR = "#a8a8a8"


# ============================================================================
# READ EIGENVALUES
# ============================================================================
def find_file(group, pattern):
    matches = sorted(glob.glob(pattern))
    if len(matches) != 1:
        sys.exit(f"{group}: expected exactly 1 file for pattern\n  {pattern}\n"
                 f"found {len(matches)}: {matches}")
    return matches[0]


def read_eigenvalues(group, path):
    """Last non-empty line -> 20 floats, named PC1-PC20 (header ignored)."""
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]
    if not lines:
        sys.exit(f"{group}: {path} is empty")
    try:
        values = [float(v) for v in lines[-1].split()]
    except ValueError:
        sys.exit(f"{group}: last line of {path} is not all numbers:\n  {lines[-1][:200]}")
    if len(values) != N_EXPECTED:
        sys.exit(f"{group}: expected {N_EXPECTED} eigenvalues in {path}, got {len(values)}")
    return pd.Series(values, index=[f"PC{i}" for i in range(1, N_EXPECTED + 1)], name=group)


eigen = {}
for group, pattern in PATTERNS.items():
    path = find_file(group, pattern)
    eigen[group] = read_eigenvalues(group, path)
    print(f"{group}: {path}")
    print("  " + " ".join(f"{v:g}" for v in eigen[group]))


# ============================================================================
# TABLE: EIGENVALUE AND % OF THE SUM OF THE 20
# ============================================================================
table = pd.DataFrame({"PC": eigen["combined"].index})
for group, values in eigen.items():
    table[f"{group}_eigenvalue"] = values.values
    table[f"{group}_pct_of_20"] = (100 * values / values.sum()).round(2).values
print("\nEigenvalues and % of the sum of PC1-20 (cutoffs: "
      + ", ".join(f"{g} PC1-{n}" for g, n in N_PCS.items()) + ")")
print(table.to_string(index=False))


# ============================================================================
# SCREE PLOT: one panel per group, own y-axis
# ============================================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
x = list(range(1, N_EXPECTED + 1))
for ax, (group, values) in zip(axes, eigen.items()):
    n = N_PCS[group]
    y = values.values
    ax.plot(x[n - 1:], y[n - 1:], color=UNUSED_COLOR, linewidth=2, zorder=1)
    ax.plot(x[:n], y[:n], color=USED_COLOR, linewidth=2, zorder=2)
    ax.scatter(x[n:], y[n:], color=UNUSED_COLOR, s=40, zorder=3, edgecolors="white", linewidths=1)
    ax.scatter(x[:n], y[:n], color=USED_COLOR, s=40, zorder=4, edgecolors="white", linewidths=1)
    ax.axvline(n + 0.5, color="#555555", linestyle="--", linewidth=1)
    ax.text(n + 0.8, 0.95, f"PC1–{n} used", transform=ax.get_xaxis_transform(),
            ha="left", va="top", fontsize=10, color="#333333")
    ax.set_title(f"{'Combined' if group == 'combined' else group} ({n} PCs)")
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Eigenvalue")
    ax.set_xticks(x)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(axis="y", color="#e6e6e6", linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
fig.suptitle("PMBB Release 2026-4.0 imputed genetic PCs: scree plots")
fig.tight_layout()

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(OUTPUT_DIR / f"pc_scree.{ext}", dpi=300)
plt.close(fig)
table.to_csv(OUTPUT_DIR / "pc_eigenvalues.tsv", sep="\t", index=False)
print(f"\nSaved {OUTPUT_DIR / 'pc_scree.png'}, pc_scree.pdf and pc_eigenvalues.tsv")
