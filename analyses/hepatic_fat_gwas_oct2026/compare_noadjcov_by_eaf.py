#!/usr/bin/env python3
"""
Pairs each variant's test statistic between two SAIGE runs that differ only in is_noadjCov
(flag on vs off) and shows, by effect-allele frequency, how much the flag shrinks χ², next to
the shrinkage predicted by the issue #206 centering error. Also lists the strongest EAF > 0.5
variants in the fixed run with their p-value in the original run.

Run with:
    python3 compare_noadjcov_by_eaf.py <flag_on_run_dir> <flag_off_run_dir> [--out noadjcov_by_eaf]
"""
import argparse
import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2

BINS = np.round(np.arange(0.0, 1.0001, 0.05), 2)
COLS = ["variant_id", "effect_allele_frequency", "beta", "variance", "p_value"]
PNPLA3 = "chr22_43928847_C_G"  # rs738409, EAF 0.23: never flipped (negative control)

parser = argparse.ArgumentParser()
parser.add_argument("on_dir")
parser.add_argument("off_dir")
parser.add_argument("--out", default="noadjcov_by_eaf")
args = parser.parse_args()


def load(run_dir, anc):
    hits = [f for f in glob.glob(f"{run_dir}/**/PMBB_{anc}.hepatic_fat.gwas.saige.gz", recursive=True)
            if "/work/" not in f]
    if len(hits) != 1:
        raise SystemExit(f"expected one {anc} sumstats file under {run_dir}, found {hits}")
    print(f"  {hits[0]}")
    d = pd.read_csv(hits[0], sep="\t", usecols=COLS).set_index("variant_id")
    d["chi2"] = chi2.isf(d["p_value"], 1)
    return d


def predicted(p):
    # χ² on / χ² off when flipped genotypes (mean 2(1-p)) are centred at 2p: 2pq / (2pq + 4(1-2p)²)
    pq2 = 2 * p * (1 - p)
    return np.where(p > 0.5, pq2 / (pq2 + 4 * (1 - 2 * p) ** 2), 1.0)


tables = {}
for anc in ("EUR", "AFR"):
    print(f"\n########## {anc}")
    on, off = load(args.on_dir, anc), load(args.off_dir, anc)
    m = on.join(off, lsuffix="_on", rsuffix="_off", how="inner")
    eaf = m["effect_allele_frequency_off"]
    print(f"  paired variants: {len(m):,}  |  EAF range {eaf.min():.3f}–{eaf.max():.3f}  |  "
          f"EAF > 0.5: {(eaf > 0.5).sum():,} ({(eaf > 0.5).mean():.1%})")
    m["ratio"] = m["chi2_on"] / m["chi2_off"]
    # Same quantity from the variance column: not affected by p-value rounding near p = 1
    m["var_ratio"] = m["variance_off"] / m["variance_on"]
    m["bin"] = pd.cut(eaf, BINS, right=False)
    t = m.groupby("bin", observed=True).agg(n=("ratio", "size"), eaf_mid=("effect_allele_frequency_off", "median"),
                                            median_ratio=("ratio", "median"),
                                            median_var_ratio=("var_ratio", "median"),
                                            lambda_on=("chi2_on", lambda c: c.median() / 0.4549),
                                            lambda_off=("chi2_off", lambda c: c.median() / 0.4549))
    t["predicted_ratio"] = predicted(t["eaf_mid"].to_numpy())
    print("\n  χ² ratio (flag on / flag off) by EAF")
    print(t.round(3).to_string())
    tables[anc] = t

    top = m[eaf > 0.5].nsmallest(5, "p_value_off")
    print("\n  Strongest 5 EAF > 0.5 and 3 EAF > 0.8 variants in the fixed run (flipped inside SAIGE), "
          "then PNPLA3 (not flipped)")
    high = m[eaf > 0.8].drop(top.index, errors="ignore").nsmallest(3, "p_value_off")
    examples = pd.concat([top, high, m[m.index == PNPLA3]])
    for vid, x in examples.iterrows():
        p = x["effect_allele_frequency_off"]
        print(f"    {vid}  EAF {p:.4f}  χ² ratio observed {x['ratio']:.3f}  predicted {float(predicted(p)):.3f}")
        for run in ("on", "off"):
            print(f"      flag {run:3s}: var {x[f'variance_{run}']:.4g}  beta {x[f'beta_{run}']:.4f}  "
                  f"p {x[f'p_value_{run}']:.2e}")

out = pd.concat(tables, names=["ancestry"]).reset_index()
out["bin"] = out["bin"].astype(str)
out.to_csv(f"{args.out}.tsv", sep="\t", index=False)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
grid = np.linspace(0.01, 0.99, 200)
for ax, (anc, t) in zip(axes, tables.items()):
    ax.axhline(1, color="#9a9a96", lw=1, zorder=1)
    ax.axvline(0.5, color="#9a9a96", lw=1, ls=":", zorder=1)
    ax.plot(grid, predicted(grid), color="#3d3d3a", lw=2, label="predicted (centering error)", zorder=2)
    ax.scatter(t["eaf_mid"], t["median_ratio"], s=60, color="#2a78d6", edgecolor="white", linewidth=2,
               label="observed (median per EAF bin)", zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(0.003, 2)
    ax.set_xlim(0, 1)
    ax.set_title(anc, fontsize=11)
    ax.set_xlabel("Effect allele frequency (AF_Allele2)")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(color="#e6e6e3", lw=0.6)
axes[0].set_ylabel("χ² with flag on / χ² with flag off\n(same variant, log scale)")
axes[0].legend(frameon=False, fontsize=9, loc="lower left")
fig.suptitle("SL INT, SAIGE 1.5.1 PGEN: is_noadjCov=TRUE shrinks χ² only when the effect allele is major",
             fontsize=11)
fig.tight_layout()
fig.savefig(f"{args.out}.png", dpi=200)
print(f"\nsaved {args.out}.tsv and {args.out}.png")