#!/usr/bin/env python3
"""
Prints λGC and mean χ² by effect-allele frequency for one or more SAIGE
summary stats files (optionally one chromosome). A well-calibrated GWAS gives
λ ≈ 1 in every EAF bin; the Oct 2026 PGEN runs collapse for EAF > 0.5.

Run with:
    python3 check_lambda_by_eaf.py <sumstats.gz> [<sumstats2.gz> ...] [--chr 22]
"""
import argparse

import pandas as pd
from scipy.stats import chi2

BINS = [(0, 0.5), (0.5, 1), (0.01, 0.05), (0.05, 0.2), (0.2, 0.5), (0.5, 0.8), (0.8, 0.95), (0.95, 0.99)]

parser = argparse.ArgumentParser()
parser.add_argument("files", nargs="+")
parser.add_argument("--chr", help="only this chromosome, e.g. 22")
args = parser.parse_args()

for path in args.files:
    d = pd.read_csv(path, sep="\t", usecols=["chromosome", "effect_allele_frequency", "p_value"],
                    dtype={"chromosome": str})
    if args.chr:
        d = d[d["chromosome"].str.replace("^chr", "", regex=True) == str(args.chr)]
    c = pd.Series(chi2.isf(d["p_value"], 1), index=d.index)
    eaf = d["effect_allele_frequency"]
    print(f"\n{path}" + (f"  [chr{args.chr}]" if args.chr else "") + f"  ({len(d):,} variants)")
    print(f"  {'all':14s} n={len(c):>10,}  lambda={c.median() / 0.4549:.3f}  mean chi2={c.mean():.3f}")
    for lo, hi in BINS:
        m = (eaf >= lo) & (eaf < hi)
        if m.any():
            print(f"  {f'EAF {lo}-{hi}':14s} n={m.sum():>10,}  lambda={c[m].median() / 0.4549:.3f}  mean chi2={c[m].mean():.3f}")
