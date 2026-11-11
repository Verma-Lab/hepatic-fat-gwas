#!/usr/bin/env python3
"""
Writes EUR and AFR SAIGE summary stats restricted to imputation R² ≥ 0.7, with the R² added as a
last column (imputation_r2). R² = WEIGHTED_R2 from the release's sites-only VCFs, matched on
variant_id as in compare_lambda_r2_eaf_filters.py; genotyped (TYPED) variants with no R² are kept
with imputation_r2 = NA. Output goes next to each input file.

Run with:
    python3 filter_sumstats_r2.py <run_dir> [--r2 0.7]
"""
import argparse
import glob
import gzip

import pandas as pd

SITES_GLOB = "/static/PMBB/PMBB-Release-2026-4.0/Imputed/chunked_pVCF_sites_only/*.vcf.gz"

parser = argparse.ArgumentParser()
parser.add_argument("run_dir", help="pipeline run folder holding PMBB_{EUR,AFR}.hepatic_fat.gwas.saige.gz")
parser.add_argument("--r2", type=float, default=0.7)
parser.add_argument("--sites_glob", default=SITES_GLOB)
args = parser.parse_args()


def find_sumstats(anc):
    hits = [f for f in glob.glob(f"{args.run_dir}/**/PMBB_{anc}.hepatic_fat.gwas.saige.gz", recursive=True)
            if "/work/" not in f]
    if len(hits) != 1:
        raise SystemExit(f"expected one {anc} sumstats file under {args.run_dir}, found {hits}")
    return hits[0]


paths = {anc: find_sumstats(anc) for anc in ("EUR", "AFR")}
wanted = set()
for anc, path in paths.items():
    print(f"{anc}: {path}")
    wanted |= set(pd.read_csv(path, sep="\t", usecols=["variant_id"])["variant_id"])

# WEIGHTED_R2 and "typed in any batch" for the variants in either results file, one pass over the
# sites VCFs (TYPED flag parsed as in check_lambda_by_eaf_r2.py)
r2, typed = {}, set()
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
            info = x[7].split(";")
            r2[x[2]] = next((i[12:] for i in info if i.startswith("WEIGHTED_R2=")), "nan")
            if any(i == "TYPED" or i.startswith("TYPED_V") for i in info):
                typed.add(x[2])

# Copy each line unchanged (no reformatting of numbers) if its R² passes, or if it is genotyped
# with no R²
for anc, path in paths.items():
    out = path.replace(".saige.gz", f".saige.r2_{args.r2}.gz")
    n = dict.fromkeys(["in", "kept", "kept_typed_no_r2", "not_in_sites", "no_r2_not_typed",
                       "low_r2_typed", "low_r2_imputed"], 0)
    with gzip.open(path, "rt") as fin, gzip.open(out, "wt", compresslevel=5) as fout:
        header = fin.readline().rstrip("\n")
        id_col = header.split("\t").index("variant_id")
        fout.write(header + "\timputation_r2\n")
        for line in fin:
            n["in"] += 1
            vid = line.split("\t", id_col + 1)[id_col]
            val = r2.get(vid)
            if val is None:
                n["not_in_sites"] += 1
                continue
            if val == "nan":
                if vid in typed:
                    fout.write(line.rstrip("\n") + "\tNA\n")
                    n["kept"] += 1
                    n["kept_typed_no_r2"] += 1
                else:
                    n["no_r2_not_typed"] += 1
                continue
            if float(val) >= args.r2:
                fout.write(line.rstrip("\n") + f"\t{val}\n")
                n["kept"] += 1
            else:
                n["low_r2_typed" if vid in typed else "low_r2_imputed"] += 1
    print(f"\n  {anc}: kept {n['kept']:,} of {n['in']:,} variants ({n['kept'] / n['in']:.1%})")
    for label, key in [("kept, genotyped with no R² (imputation_r2 = NA)", "kept_typed_no_r2"),
                       ("dropped, not in sites VCFs", "not_in_sites"),
                       ("dropped, no R² and not genotyped", "no_r2_not_typed"),
                       (f"dropped, R² < {args.r2}, imputed only", "low_r2_imputed"),
                       (f"dropped, R² < {args.r2}, genotyped in some batch", "low_r2_typed")]:
        print(f"    {label + ':':<50} {n[key]:,}")
    print(f"    saved {out}")