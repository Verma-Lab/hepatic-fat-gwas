"""
Prints λGC and mean χ² by effect-allele frequency, split by imputation quality
(WEIGHTED_R2 from the release's sites-only VCFs), for one SAIGE summary stats file.
Tests whether the EAF > 0.5 deflation goes away among well-imputed variants.

Run with:
    python3 check_lambda_by_eaf_r2.py <sumstats.gz> "<sites_only_dir>/*.vcf.gz" [--chr 22] [--r2 0.7]
"""
import argparse
import glob
import gzip

import pandas as pd
from scipy.stats import chi2

BINS = [(0, 0.5), (0.5, 1), (0.01, 0.05), (0.05, 0.2), (0.2, 0.5), (0.5, 0.8), (0.8, 0.95), (0.95, 0.99)]

parser = argparse.ArgumentParser()
parser.add_argument("sumstats")
parser.add_argument("sites_glob", help='quoted glob of sites-only VCFs, e.g. "/static/.../chunked_pVCF_sites_only/*.vcf.gz"')
parser.add_argument("--chr", help="only this chromosome, e.g. 22")
parser.add_argument("--r2", type=float, default=0.7)
args = parser.parse_args()

d = pd.read_csv(args.sumstats, sep="\t", usecols=["chromosome", "variant_id", "effect_allele_frequency", "p_value"],
                dtype={"chromosome": str})
d["chromosome"] = d["chromosome"].str.replace("^chr", "", regex=True)
if args.chr:
    d = d[d["chromosome"] == str(args.chr)]
d = d.set_index("variant_id")

# WEIGHTED_R2 and "typed in any batch" for the variants in the sumstats only
wanted = set(d.index)
r2, typed = {}, {}
files = sorted(glob.glob(args.sites_glob))
if args.chr:
    files = [f for f in files if f"_chr{args.chr}_" in f]
print(f"Reading {len(files)} sites-only VCFs")
for f in files:
    with gzip.open(f, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            x = line.rstrip("\n").split("\t", 8)
            if x[2] not in wanted:
                continue
            info = x[7].split(";")
            r2[x[2]] = next((float(i[12:]) for i in info if i.startswith("WEIGHTED_R2=")), float("nan"))
            typed[x[2]] = any(i == "TYPED" or i.startswith("TYPED_V") for i in info)

d["r2"] = pd.Series(r2)
d["typed"] = pd.Series(typed)
print(f"{args.sumstats}" + (f"  [chr{args.chr}]" if args.chr else ""))
print(f"  matched to sites VCFs: {d['r2'].notna().sum():,} of {len(d):,}")

d["chi2"] = chi2.isf(d["p_value"], 1)
groups = [
    ("all", d["r2"].notna()),
    (f"R2 >= {args.r2}", d["r2"] >= args.r2),
    (f"R2 >= {args.r2}, not typed", (d["r2"] >= args.r2) & ~d["typed"].astype(bool)),
    (f"R2 < {args.r2}", d["r2"] < args.r2),
]
for name, g in groups:
    s = d[g]
    print(f"\n  [{name}]  n={len(s):,}  lambda={s['chi2'].median() / 0.4549:.3f}  mean chi2={s['chi2'].mean():.3f}")
    for lo, hi in BINS:
        m = (s["effect_allele_frequency"] >= lo) & (s["effect_allele_frequency"] < hi)
        if m.any():
            c = s.loc[m, "chi2"]
            print(f"    {f'EAF {lo}-{hi}':14s} n={m.sum():>10,}  lambda={c.median() / 0.4549:.3f}  mean chi2={c.mean():.3f}")
