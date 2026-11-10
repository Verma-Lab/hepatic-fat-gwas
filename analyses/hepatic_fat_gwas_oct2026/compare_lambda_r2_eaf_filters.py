#!/usr/bin/env python3
"""
Genome-wide λGC for EUR and AFR under three variant filters: none, imputation R² ≥ 0.7, and
R² ≥ 0.7 with 0.01 < EAF < 0.5. R² = WEIGHTED_R2 from the release's sites-only VCFs, matched
on variant_id as in check_lambda_by_eaf_r2.py. Also reports the EAF range of each results file.

Run with:
    python3 compare_lambda_r2_eaf_filters.py <run_dir> [--out lambda_filters.tsv] [--r2 0.7]
"""
import argparse
import glob
import gzip

import pandas as pd
from scipy.stats import chi2

SITES_GLOB = "/static/PMBB/PMBB-Release-2026-4.0/Imputed/chunked_pVCF_sites_only/*.vcf.gz"

parser = argparse.ArgumentParser()
parser.add_argument("run_dir", help="pipeline run folder holding PMBB_{EUR,AFR}.hepatic_fat.gwas.saige.gz")
parser.add_argument("--out", default="lambda_filters.tsv")
parser.add_argument("--r2", type=float, default=0.7)
parser.add_argument("--sites_glob", default=SITES_GLOB)
args = parser.parse_args()


def find_sumstats(anc):
    hits = [f for f in glob.glob(f"{args.run_dir}/**/PMBB_{anc}.hepatic_fat.gwas.saige.gz", recursive=True)
            if "/work/" not in f]
    if len(hits) != 1:
        raise SystemExit(f"expected one {anc} sumstats file under {args.run_dir}, found {hits}")
    return hits[0]


ss = {}
for anc in ("EUR", "AFR"):
    path = find_sumstats(anc)
    print(f"{anc}: {path}")
    ss[anc] = pd.read_csv(path, sep="\t", usecols=["variant_id", "effect_allele_frequency", "p_value"]
                          ).set_index("variant_id")

# WEIGHTED_R2 for the variants in either results file, one pass over the sites VCFs
wanted = set(ss["EUR"].index) | set(ss["AFR"].index)
r2 = {}
files = sorted(glob.glob(args.sites_glob))
print(f"Reading {len(files)} sites-only VCFs")
for f in files:
    with gzip.open(f, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            x = line.rstrip("\n").split("\t", 8)
            if x[2] not in wanted:
                continue
            r2[x[2]] = next((float(i[12:]) for i in x[7].split(";") if i.startswith("WEIGHTED_R2=")), float("nan"))
r2 = pd.Series(r2)

print("\nEAF range check")
for anc, d in ss.items():
    eaf = d["effect_allele_frequency"]
    d["r2"] = r2.reindex(d.index)
    d["chi2"] = chi2.isf(d["p_value"], 1)
    print(f"  {anc}: min {eaf.min():.5f}  max {eaf.max():.5f}  EAF < 0.01: {(eaf < 0.01).sum():,}  "
          f"EAF > 0.99: {(eaf > 0.99).sum():,}  matched to R²: {d['r2'].notna().sum():,} of {len(d):,}")

filters = {
    "no filter": lambda d: pd.Series(True, index=d.index),
    f"R2 >= {args.r2}": lambda d: d["r2"] >= args.r2,
    f"R2 >= {args.r2} & 0.01 < EAF < 0.5": lambda d: (d["r2"] >= args.r2) & (d["effect_allele_frequency"] > 0.01)
                                                    & (d["effect_allele_frequency"] < 0.5),
}
rows = []
for name, f in filters.items():
    row = {"filter": name}
    for anc, d in ss.items():
        c = d.loc[f(d), "chi2"]
        row[f"{anc}_lambda"] = round(c.median() / 0.4549, 3)
        row[f"{anc}_n"] = len(c)
    rows.append(row)
table = pd.DataFrame(rows)
print("\n" + table.to_string(index=False))
table.to_csv(args.out, sep="\t", index=False)
print(f"\nsaved {args.out}")