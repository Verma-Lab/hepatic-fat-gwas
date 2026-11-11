#!/usr/bin/env python3
"""
Formats SAIGE summary stats (GRCh38, no rsIDs) as FUMA input.

FUMA >= 2.1.0 accepts GRCh38 positions without rsIDs only if the GRCh38 option
is ticked on the submission page AND the file is tab-separated with exactly
this header, in this order (otherwise the job fails):
    chromosome  base_pair_location  effect_allele  other_allele  beta  p_value
FUMA looks up rsIDs (dbSNP v157) by chromosome/position, but its reference data
is GRCh37, so results come back in GRCh37 (PNPLA3 rs738409: chr22:43928847 in
GRCh38 -> chr22:44324727 in GRCh37). N is entered on the submission page.

Run with:
    python3 prepare_fuma_input.py                                  # SL_transformed, EUR + AFR
    python3 prepare_fuma_input.py --run SL_untransformed --ancestry EUR
    python3 prepare_fuma_input.py --suffix r2_0.7                  # R² ≥ 0.7 files from filter_sumstats_r2.py
"""
import argparse
import gzip
import io
import os
import sys
from pathlib import Path

import pandas as pd

# ============================================================================
# SETTINGS
# ============================================================================
BASE = Path("/project/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_oct2026")
IN_COLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "p_value"]
OUT_COLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "p_value"]  # FUMA layout (b)
AUTOSOMES = {str(i) for i in range(1, 23)}
CHUNKSIZE = 10_000_000
GWS = 5e-8
N_BY_ANCESTRY = {"EUR": 12867, "AFR": 4540}
PNPLA3 = ("22", "43928847")      # rs738409, GRCh38; expect effect G / other C, beta > 0


def sumstats_path(args, ancestry):
    name = f"PMBB_{ancestry}.hepatic_fat.gwas.saige{'.' + args.suffix if args.suffix else ''}.gz"
    if args.sumstats_dir:
        return Path(args.sumstats_dir) / name
    return BASE / args.run / f"PMBB_{ancestry}" / "Sumstats" / name


# ============================================================================
# FORMAT ONE FILE
# ============================================================================
def format_one(in_path, out_path):
    header = pd.read_csv(in_path, sep="\t", nrows=0).columns.tolist()
    missing = [c for c in IN_COLS if c not in header]
    if missing:
        sys.exit(f"{in_path}: missing column(s) {missing}\n  header: {header}")

    drops = {"not autosome 1-22": 0, "missing/non-numeric position": 0, "missing allele": 0,
             "missing/non-numeric beta": 0, "missing/non-numeric p": 0, "p not in (0, 1]": 0}
    n_in = n_out = n_gws = 0
    pnpla3 = None
    tmp_path = out_path.with_name(out_path.name + ".partial")

    # write to .partial, but store the final name inside the gzip (FUMA checks files weren't renamed)
    with open(tmp_path, "wb") as raw, \
            gzip.GzipFile(filename=out_path.name[:-3], mode="wb", fileobj=raw, compresslevel=6) as gz, \
            io.TextIOWrapper(gz, encoding="utf-8") as out:
        out.write("\t".join(OUT_COLS) + "\n")
        reader = pd.read_csv(in_path, sep="\t", usecols=IN_COLS, dtype=str, chunksize=CHUNKSIZE)
        for i, chunk in enumerate(reader, 1):
            n_in += len(chunk)
            chunk = chunk.apply(lambda s: s.str.strip())
            chunk["chromosome"] = chunk["chromosome"].str.replace(r"^chr", "", case=False, regex=True)

            # drop rows in this order; each row is counted under the first reason that applies
            keep = chunk["chromosome"].isin(AUTOSOMES)
            drops["not autosome 1-22"] += int((~keep).sum())
            pos = pd.to_numeric(chunk["base_pair_location"], errors="coerce")
            bad = keep & pos.isna()
            drops["missing/non-numeric position"] += int(bad.sum()); keep &= ~bad
            bad = keep & (chunk["effect_allele"].isna() | chunk["other_allele"].isna()
                          | (chunk["effect_allele"] == "") | (chunk["other_allele"] == ""))
            drops["missing allele"] += int(bad.sum()); keep &= ~bad
            beta = pd.to_numeric(chunk["beta"], errors="coerce")
            bad = keep & beta.isna()
            drops["missing/non-numeric beta"] += int(bad.sum()); keep &= ~bad
            p = pd.to_numeric(chunk["p_value"], errors="coerce")
            bad = keep & p.isna()
            drops["missing/non-numeric p"] += int(bad.sum()); keep &= ~bad
            bad = keep & ~((p > 0) & (p <= 1))
            drops["p not in (0, 1]"] += int(bad.sum()); keep &= ~bad

            kept = chunk.loc[keep, OUT_COLS].copy()
            kept["base_pair_location"] = pos[keep].astype("int64").astype(str)   # no "1234.0"
            n_out += len(kept)
            n_gws += int((p[keep] < GWS).sum())
            hit = kept[(kept["chromosome"] == PNPLA3[0]) & (kept["base_pair_location"] == PNPLA3[1])]
            if len(hit):
                pnpla3 = hit
            kept.to_csv(out, sep="\t", header=False, index=False)   # beta/p written as in the input
            print(f"  chunk {i}: {len(chunk):,} rows read, {len(kept):,} kept")

    os.replace(tmp_path, out_path)
    return n_in, n_out, drops, n_gws, pnpla3


# ============================================================================
# MAIN
# ============================================================================
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", default="SL_transformed", help="run folder (default: SL_transformed)")
    parser.add_argument("--ancestry", choices=["EUR", "AFR", "both"], default="both")
    parser.add_argument("--outdir", default=str(BASE / "fuma_input"))
    parser.add_argument("--overwrite", action="store_true", help="replace existing output files")
    parser.add_argument("--sumstats-dir", help="read PMBB_<anc>.hepatic_fat.gwas.saige.gz from this folder "
                                               "instead of <run>/PMBB_<anc>/Sumstats/ (e.g. a local copy)")
    parser.add_argument("--suffix", help="read PMBB_<anc>.hepatic_fat.gwas.saige.<suffix>.gz instead, e.g. r2_0.7 "
                                         "for the filter_sumstats_r2.py output; also added to the output name")
    args = parser.parse_args()
    ancestries = ["EUR", "AFR"] if args.ancestry == "both" else [args.ancestry]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # check every input/output before doing any work
    jobs = []
    for anc in ancestries:
        in_path = sumstats_path(args, anc)
        tag = f"{args.run}_{anc}_{args.suffix}" if args.suffix else f"{args.run}_{anc}"
        out_path = outdir / f"{tag}_fuma_GRCh38.txt.gz"
        if not in_path.exists():
            sys.exit(f"Input not found: {in_path}")
        if out_path.exists() and not args.overwrite:
            sys.exit(f"Output already exists: {out_path}\n  Use --overwrite to replace it.")
        jobs.append((anc, in_path, out_path))

    for anc, in_path, out_path in jobs:
        print(f"\n{'=' * 70}\n{args.run} {anc}\n{'=' * 70}")
        print(f"Input:  {in_path}\nOutput: {out_path}")
        n_in, n_out, drops, n_gws, pnpla3 = format_one(in_path, out_path)

        print(f"\nRows in: {n_in:,} | rows out: {n_out:,} | dropped: {n_in - n_out:,}")
        for reason, n in drops.items():
            print(f"  dropped, {reason}: {n:,}")
        print("Header written (tab-separated): " + "\t".join(OUT_COLS))
        if pnpla3 is None:
            print(f"PNPLA3 rs738409 ({PNPLA3[0]}:{PNPLA3[1]}): NOT FOUND in output")
        else:
            print(f"PNPLA3 rs738409 ({PNPLA3[0]}:{PNPLA3[1]}):")
            print(pnpla3.to_string(index=False))
            r = pnpla3.iloc[0]
            ok = r["effect_allele"] == "G" and r["other_allele"] == "C" and float(r["beta"]) > 0
            print("  -> G/C with positive beta: OK" if ok else "  -> CHECK: expected effect G, other C, beta > 0")
        print(f"Variants with p < {GWS:g}: {n_gws:,}")
        print(f"Output size: {out_path.stat().st_size / 1e6:,.1f} MB")

    print(f"\n{'=' * 70}\nFUMA submission settings\n{'=' * 70}")
    print("  - Tick the GRCh38 option.")
    print("  - Leave the column-name section empty (the header is already in FUMA's format).")
    for anc in ancestries:
        print(f"  - {anc}: N = {N_BY_ANCESTRY[anc]}; reference panel 1000G Phase 3 {anc}.")
    print("  - Keep the default lead-SNP settings.")
    print("  - Results come back in GRCh37 (e.g. PNPLA3 rs738409 at chr22:44324727).")


if __name__ == "__main__":
    main()
