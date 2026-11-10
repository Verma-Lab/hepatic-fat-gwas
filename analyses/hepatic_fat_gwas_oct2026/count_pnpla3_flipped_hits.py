#!/usr/bin/env python3
"""
Counts variants within ±500 kb of PNPLA3 rs738409 (chr22:43928847, GRCh38) with EAF > 0.5 and
p < 5e-8 in the flag-off run, and how many of those are still p < 5e-8 in the flag-on run.

Run with:
    python3 count_pnpla3_flipped_hits.py <flag_on_run_dir> <flag_off_run_dir>
"""
import argparse
import glob

import pandas as pd

POS, WINDOW, GWS = 43928847, 500_000, 5e-8
COLS = ["chromosome", "base_pair_location", "variant_id", "effect_allele_frequency", "p_value"]

parser = argparse.ArgumentParser()
parser.add_argument("on_dir")
parser.add_argument("off_dir")
args = parser.parse_args()


def load(run_dir, anc):
    hits = [f for f in glob.glob(f"{run_dir}/**/PMBB_{anc}.hepatic_fat.gwas.saige.gz", recursive=True)
            if "/work/" not in f]
    if len(hits) != 1:
        raise SystemExit(f"expected one {anc} sumstats file under {run_dir}, found {hits}")
    d = pd.read_csv(hits[0], sep="\t", usecols=COLS, dtype={"chromosome": str})
    d = d[(d["chromosome"].str.replace("^chr", "", regex=True) == "22")
          & (d["base_pair_location"].sub(POS).abs() <= WINDOW)]
    return d.set_index("variant_id")


for anc in ("EUR", "AFR"):
    m = load(args.on_dir, anc).join(load(args.off_dir, anc)[["p_value"]], lsuffix="_on", rsuffix="_off", how="inner")
    major = m["effect_allele_frequency"] > 0.5
    hit_off = m["p_value_off"] < GWS
    print(f"\n{anc}: {len(m):,} variants within ±500 kb of rs738409 ({major.sum():,} with EAF > 0.5)")
    print(f"  p < 5e-8, flag off:                 {hit_off.sum():,} total, {(hit_off & major).sum():,} with EAF > 0.5")
    print(f"  EAF > 0.5 hits still p < 5e-8 on:   {(hit_off & major & (m['p_value_on'] < GWS)).sum():,}")
    print(f"  EAF < 0.5 hits still p < 5e-8 on:   {(hit_off & ~major & (m['p_value_on'] < GWS)).sum():,}"
          f" of {(hit_off & ~major).sum():,}")