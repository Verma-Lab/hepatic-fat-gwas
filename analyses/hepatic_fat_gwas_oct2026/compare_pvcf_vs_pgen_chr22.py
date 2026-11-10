#!/usr/bin/env python3
"""
Combines the chr22 pVCF step 2 chunks (raw SAIGE output) for one ancestry and compares them
with the PGEN run: λ by effect-allele frequency for each run, on all variants and on the
variants shared by both (matched by position and allele pair), plus PNPLA3 rs738409.

Run with:
    python3 compare_pvcf_vs_pgen_chr22.py <pvcf_dir> <pgen_sumstats.gz>
"""
import argparse
import glob

import numpy as np
import pandas as pd
from scipy.stats import chi2

BINS = [(0, 0.5), (0.5, 1), (0.01, 0.05), (0.05, 0.2), (0.2, 0.5), (0.5, 0.8), (0.8, 0.95), (0.95, 0.99)]
RAW = {"CHR": "chromosome", "POS": "base_pair_location", "MarkerID": "variant_id", "Allele1": "other_allele",
       "Allele2": "effect_allele", "AF_Allele2": "effect_allele_frequency", "BETA": "beta",
       "SE": "standard_error", "p.value": "p_value", "N": "N"}
PNPLA3_POS = 43928847

parser = argparse.ArgumentParser()
parser.add_argument("pvcf_dir")
parser.add_argument("pgen")
args = parser.parse_args()


def add_key(d):
    d["chromosome"] = d["chromosome"].astype(str).str.replace("^chr", "", regex=True)
    alleles = np.sort(d[["other_allele", "effect_allele"]].to_numpy(dtype=str), axis=1)
    d["key"] = d["chromosome"] + ":" + d["base_pair_location"].astype(str) + ":" + alleles[:, 0] + ":" + alleles[:, 1]
    return d.drop_duplicates("key", keep=False).set_index("key")


def lambda_table(name, d):
    c = pd.Series(chi2.isf(d["p_value"], 1), index=d.index)
    eaf = d["effect_allele_frequency"]
    print(f"\n  [{name}]  n={len(d):,}  lambda={c.median() / 0.4549:.3f}  mean chi2={c.mean():.3f}")
    for lo, hi in BINS:
        m = (eaf >= lo) & (eaf < hi)
        if m.any():
            print(f"    {f'EAF {lo}-{hi}':14s} n={m.sum():>9,}  lambda={c[m].median() / 0.4549:.3f}  mean chi2={c[m].mean():.3f}")


files = sorted(f for f in glob.glob(f"{args.pvcf_dir}/PMBB_*.txt*") if not f.endswith(".index"))
vcf = pd.concat([pd.read_csv(f, sep=r"\s+") for f in files], ignore_index=True)
vcf = add_key(vcf.rename(columns=RAW)[list(RAW.values())])
pgen = pd.read_csv(args.pgen, sep="\t", usecols=list(RAW.values()), dtype={"chromosome": str})
pgen = add_key(pgen[pgen["chromosome"].str.replace("^chr", "", regex=True) == "22"])

print(f"{len(files)} pVCF chunk files | pVCF variants {len(vcf):,} | PGEN chr22 variants {len(pgen):,}")
print(f"N: pVCF {vcf['N'].median():.0f} | PGEN {pgen['N'].median():.0f}")
shared = vcf.index.intersection(pgen.index)
print(f"shared {len(shared):,} | pVCF only {len(vcf) - len(shared):,} | PGEN only {len(pgen) - len(shared):,}")
same_ea = (vcf.loc[shared, "effect_allele"] == pgen.loc[shared, "effect_allele"]).mean()
print(f"same effect allele in both runs: {same_ea:.1%}")

lambda_table("pVCF, all variants", vcf)
lambda_table("pVCF, shared variants", vcf.loc[shared])
lambda_table("PGEN, shared variants", pgen.loc[shared])

print("\nPNPLA3 rs738409 (chr22:43928847):")
for name, d in (("pVCF", vcf), ("PGEN", pgen)):
    r = d[d["base_pair_location"] == PNPLA3_POS]
    for _, x in r.iterrows():
        print(f"  {name}: {x['other_allele']}>{x['effect_allele']}  EAF {x['effect_allele_frequency']:.3f}  "
              f"beta {x['beta']:.4f}  SE {x['standard_error']:.4f}  p {x['p_value']:.2e}")