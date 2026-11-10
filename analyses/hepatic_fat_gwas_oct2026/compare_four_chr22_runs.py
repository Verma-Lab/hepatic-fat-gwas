#!/usr/bin/env python3
"""
Compares chr22 SAIGE step 2 results from several runs for one ancestry: variant counts,
coverage and N; λ for EAF < 0.5 / ≥ 0.5 on each run's own variants and on the variants shared
by all runs (matched on chr:pos and allele pair, either orientation); PNPLA3 rs738409.
Each run is a pipeline sumstats file (.gz) or a folder of raw SAIGE chunk outputs.

Run with:
    python3 compare_four_chr22_runs.py EUR "PGEN 1.5.1=<path>" "pVCF 1.5.1=<path>" "BED 1.5.1=<path>" "BED 1.4.5=<path>"
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import chi2

COLS = {"CHR": "chromosome", "POS": "base_pair_location", "Allele1": "other_allele", "Allele2": "effect_allele",
        "AF_Allele2": "effect_allele_frequency", "BETA": "beta", "SE": "standard_error", "p.value": "p_value", "N": "N"}
PNPLA3_POS = 43928847  # rs738409, GRCh38

parser = argparse.ArgumentParser()
parser.add_argument("ancestry")
parser.add_argument("runs", nargs="+", help='"label=path" for each run')
args = parser.parse_args()


def load(path):
    if os.path.isdir(path):
        files = sorted(f for f in glob.glob(f"{path}/PMBB_*.txt*") if not f.endswith(".index"))
        d = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    else:
        d = pd.read_csv(path, sep="\t")
    d = d.rename(columns=COLS)[list(COLS.values())]
    d["chromosome"] = d["chromosome"].astype(str).str.replace("^chr", "", regex=True)
    d = d[d["chromosome"] == "22"].copy()
    alleles = np.sort(d[["other_allele", "effect_allele"]].to_numpy(dtype=str), axis=1)
    d["key"] = "22:" + d["base_pair_location"].astype(str) + ":" + alleles[:, 0] + ":" + alleles[:, 1]
    n_dup = d["key"].duplicated(keep=False).sum()
    d = d.drop_duplicates("key", keep=False).set_index("key")
    d["chi2"] = chi2.isf(d["p_value"], 1)
    return d, n_dup


def lam(c):
    return c.median() / 0.4549


runs = {}
for r in args.runs:
    label, path = r.split("=", 1)
    runs[label] = (path, *load(path))

print(f"########## {args.ancestry}\n\nResult files:")
for label, (path, _, _) in runs.items():
    print(f"  {label:12s} {path}")

print("\nA. Variant sets")
print(f"  {'run':12s} {'variants':>9s} {'dup keys dropped':>17s} {'min pos':>11s} {'max pos':>11s} {'N':>7s}")
for label, (_, d, n_dup) in runs.items():
    print(f"  {label:12s} {len(d):>9,} {n_dup:>17,} {d['base_pair_location'].min():>11,} "
          f"{d['base_pair_location'].max():>11,} {d['N'].median():>7.0f}")
shared = None
for _, d, _ in runs.values():
    shared = d.index if shared is None else shared.intersection(d.index)
print(f"  shared by all {len(runs)} runs: {len(shared):,}")

for title, subset in (("B1. λ by EAF, all variants in each run", None), ("B2. λ by EAF, shared variants", shared)):
    print(f"\n{title}  (EAF from each run's own effect allele)")
    print(f"  {'run':12s} {'λ EAF<0.5':>10s} {'n':>9s} {'λ EAF≥0.5':>10s} {'n':>9s} {'λ all':>7s}")
    for label, (_, d, _) in runs.items():
        s = d if subset is None else d.loc[subset]
        lo = s["effect_allele_frequency"] < 0.5
        print(f"  {label:12s} {lam(s.loc[lo, 'chi2']):>10.3f} {lo.sum():>9,} "
              f"{lam(s.loc[~lo, 'chi2']):>10.3f} {(~lo).sum():>9,} {lam(s['chi2']):>7.3f}")

print("\nC. PNPLA3 rs738409 (chr22:43928847, GRCh38)")
print(f"  {'run':12s} {'EA':>3s} {'EAF':>6s} {'beta':>8s} {'SE':>7s} {'p':>9s}")
for label, (_, d, _) in runs.items():
    for _, x in d[d["base_pair_location"] == PNPLA3_POS].iterrows():
        print(f"  {label:12s} {x['effect_allele']:>3s} {x['effect_allele_frequency']:>6.3f} {x['beta']:>8.4f} "
              f"{x['standard_error']:>7.4f} {x['p_value']:>9.2e}")