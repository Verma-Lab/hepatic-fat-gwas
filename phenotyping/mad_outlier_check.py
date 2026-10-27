#!/usr/bin/env python3
"""
Sensitivity check: the current ±6 SD study-level outlier rule vs MAD-based rules
(median ± k × scaled MAD), per organ × phase. Study values are built as in
build_hepatic_fat_phenotype.py. Writes no phenotype files.

Run with:
    python3 mad_outlier_check.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# SETTINGS
# ============================================================================
CLEAN_CACHE_FILE = Path("/project/verma_shared/projects/Liver_IDPs/phenotyping/idp_liver_spleen_clean.pkl")
FAM_FILE = "/static/PMBB/PMBB-Release-2026-4.0/Imputed/chunked_bed_files/PMBB-Release-2026-4.0_genetic_imputed.463_chr8_chunk1_149284_2198252.fam"
OUTPUT_DIR = Path("/project/verma_shared/projects/Liver_IDPs/phenotyping/mad_check")

PHASES = ["native", "portal_venous"]
ORGANS = ["liver", "spleen"]
SD_CUTOFF = 6
MAD_KS = [3, 4, 5, 6, 8]
PLOT_MAD_KS = [4, 6, 8]
NATIVE_SPLEEN_MAX_HU = 100
NATIVE_LIVER_FAT_HU = 10          # native liver below this = likely real severe fat

EXPECTED_SD_REMOVED = {("liver", "native"): 8, ("liver", "portal_venous"): 6,
                       ("spleen", "native"): 117, ("spleen", "portal_venous"): 43}
EXPECTED_N_STUDIES = 141_894


def section(title):
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")


# ============================================================================
# STUDY VALUES (same as build_hepatic_fat_phenotype.py load_and_qc)
# ============================================================================
idp = pd.read_pickle(CLEAN_CACHE_FILE)
idp = idp[idp["phase"].isin(PHASES)].copy()
idp["person_id"] = idp["person_id"].astype(str)
studies = (
    idp.groupby(["person_id", "accession", "anatomical_structure", "phase"])["attenuation_Mean"]
    .median().rename("value").reset_index()
    .rename(columns={"anatomical_structure": "organ"})
)
status = "matches expected" if len(studies) == EXPECTED_N_STUDIES else f"EXPECTED {EXPECTED_N_STUDIES:,}"
print(f"Study values (median attenuation_Mean across series, clipping QC'd): {len(studies):,} ({status})")

geno_ids = set(pd.read_csv(FAM_FILE, sep=r"\s+", header=None, usecols=[1], dtype=str)[1])
studies["genotyped"] = studies["person_id"].isin(geno_ids)
print(f"Genotyped IDs in .fam: {len(geno_ids):,}; study values from genotyped people: {studies['genotyped'].sum():,}")


# ============================================================================
# BOUNDS AND REMOVAL MASKS PER RULE
# ============================================================================
def rule_bounds(x, rule, k):
    if rule == "SD":
        m, s = x.mean(), x.std()                                  # ddof = 1
    else:
        m, s = x.median(), median_abs_deviation(x, scale="normal")  # × 1.4826
    return m - k * s, m + k * s


RULES = [("SD", SD_CUTOFF)] + [("MAD", k) for k in MAD_KS]
removed = {r: pd.Series(False, index=studies.index) for r in RULES}   # removed by the rule itself
range_added = {r: pd.Series(False, index=studies.index) for r in RULES}
rows = []
bounds_for_plot = {}

for organ in ORGANS:
    for phase in PHASES:
        grp = studies[(studies["organ"] == organ) & (studies["phase"] == phase)]
        x = grp["value"]
        for rule, k in RULES:
            lo, hi = rule_bounds(x, rule, k)
            bounds_for_plot[(organ, phase, rule, k)] = (lo, hi)
            low, high = x < lo, x > hi
            out = low | high
            removed[(rule, k)].loc[grp.index] = out

            kept_native_spleen_high = (~out) & (x > NATIVE_SPLEEN_MAX_HU) if (organ, phase) == ("spleen", "native") \
                else pd.Series(False, index=grp.index)
            range_added[(rule, k)].loc[grp.index] = kept_native_spleen_high

            sd_out = removed[("SD", SD_CUTOFF)].loc[grp.index]
            row = {"organ": organ, "phase": phase, "rule": rule, "k": k,
                   "lower_HU": round(lo, 1), "upper_HU": round(hi, 1), "n_studies": len(x),
                   "n_low": int(low.sum()), "n_high": int(high.sum()), "n_total": int(out.sum()),
                   "both_with_sd6": int((out & sd_out).sum()),
                   "rule_only": int((out & ~sd_out).sum()),
                   "sd6_only": int((~out & sd_out).sum()),
                   "native_spleen_gt100_added": int(kept_native_spleen_high.sum()) if (organ, phase) == ("spleen", "native") else np.nan,
                   f"native_liver_lt{NATIVE_LIVER_FAT_HU}_removed": int((out & (x < NATIVE_LIVER_FAT_HU)).sum()) if (organ, phase) == ("liver", "native") else np.nan}
            rows.append(row)

summary = pd.DataFrame(rows)
for col in ["native_spleen_gt100_added", f"native_liver_lt{NATIVE_LIVER_FAT_HU}_removed"]:
    summary[col] = summary[col].astype("Int64")      # whole numbers, NA where not applicable

# People who lose all their studies in that organ × phase (genotyped only)
lose_all = []
for (organ, phase, rule, k) in summary[["organ", "phase", "rule", "k"]].itertuples(index=False):
    grp = studies[(studies["organ"] == organ) & (studies["phase"] == phase) & studies["genotyped"]]
    gone = removed[(rule, k)].loc[grp.index] | range_added[(rule, k)].loc[grp.index]
    kept_people = set(grp.loc[~gone, "person_id"])
    lose_all.append(grp.loc[gone, "person_id"].drop_duplicates().isin(kept_people).eq(False).sum())
summary["geno_people_lose_all"] = lose_all


# ============================================================================
# PRINT
# ============================================================================
section(f"Current rule check: ±{SD_CUTOFF} SD (should match 8 / 6 / 117 / 43)")
for _, r in summary[summary["rule"] == "SD"].iterrows():
    exp = EXPECTED_SD_REMOVED[(r["organ"], r["phase"])]
    tag = "matches expected" if r["n_total"] == exp else f"EXPECTED {exp}"
    print(f"  {r['organ']:6s} {r['phase']:14s} [{r['lower_HU']}, {r['upper_HU']}]  "
          f"low {r['n_low']:,} | high {r['n_high']:,} | total {r['n_total']:,}  ({tag})")

section("Per organ × phase: bounds (HU) and studies removed by each rule")
pd.set_option("display.width", 250)
show = summary.drop(columns=["n_studies"]).copy()
print(show.to_string(index=False, na_rep="-"))

section(f"Native liver studies below {NATIVE_LIVER_FAT_HU} HU (likely real severe fat) removed by each rule")
nl = summary[(summary["organ"] == "liver") & (summary["phase"] == "native")]
n_fat = ((studies["organ"] == "liver") & (studies["phase"] == "native") & (studies["value"] < NATIVE_LIVER_FAT_HU)).sum()
print(f"  native liver studies below {NATIVE_LIVER_FAT_HU} HU: {n_fat:,}")
for _, r in nl.iterrows():
    print(f"  {r['rule']} k={r['k']}: lower bound {r['lower_HU']} HU -> removes "
          f"{int(r[f'native_liver_lt{NATIVE_LIVER_FAT_HU}_removed']):,} of them")

section("Totals per rule (rule + native spleen > 100 HU check)")
tot_rows = []
for rule, k in RULES:
    any_gone = removed[(rule, k)] | range_added[(rule, k)]
    g = studies[studies["genotyped"]]
    gone_g = any_gone.loc[g.index]
    # genotyped people whose every study (any organ/phase) is removed
    kept_people = set(g.loc[~gone_g, "person_id"])
    lost_people = g.loc[gone_g, "person_id"].drop_duplicates()
    tot_rows.append({"rule": rule, "k": k,
                     "removed_by_rule": int(removed[(rule, k)].sum()),
                     "added_by_native_spleen_gt100": int(range_added[(rule, k)].sum()),
                     "total_removed": int(any_gone.sum()),
                     "pct_of_studies": round(100 * any_gone.sum() / len(studies), 3),
                     "geno_people_lose_all_studies": int((~lost_people.isin(kept_people)).sum())})
totals = pd.DataFrame(tot_rows)
print(f"(out of {len(studies):,} study values)")
print(totals.to_string(index=False))
print("\ngeno_people_lose_all (per organ × phase table): genotyped people left with no study in that organ × phase.")
print("geno_people_lose_all_studies (totals): genotyped people left with no study in any organ × phase.")


# ============================================================================
# SAVE
# ============================================================================
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
summary.to_csv(OUTPUT_DIR / "mad_vs_sd_summary.tsv", sep="\t", index=False, na_rep="NA")
totals.to_csv(OUTPUT_DIR / "mad_vs_sd_totals.tsv", sep="\t", index=False)

mad_colors = {4: "#e08a00", 6: "#7b3fb5", 8: "#1a9e77"}
for organ in ORGANS:
    for phase in PHASES:
        x = studies.loc[(studies["organ"] == organ) & (studies["phase"] == phase), "value"]
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.hist(x, bins=150, color="#8fb3d9", log=True)
        lo, hi = bounds_for_plot[(organ, phase, "SD", SD_CUTOFF)]
        ax.axvline(lo, color="#c0392b", linestyle="--", linewidth=1.5, label=f"±{SD_CUTOFF} SD (current)")
        ax.axvline(hi, color="#c0392b", linestyle="--", linewidth=1.5)
        for k in PLOT_MAD_KS:
            lo, hi = bounds_for_plot[(organ, phase, "MAD", k)]
            ax.axvline(lo, color=mad_colors[k], linestyle=":", linewidth=2, label=f"MAD k = {k}")
            ax.axvline(hi, color=mad_colors[k], linestyle=":", linewidth=2)
        if (organ, phase) == ("spleen", "native"):
            ax.axvline(NATIVE_SPLEEN_MAX_HU, color="#555555", linestyle="-.", linewidth=1,
                       label=f"native spleen > {NATIVE_SPLEEN_MAX_HU} HU check")
        ax.set_title(f"{organ}, {phase.replace('_', ' ')}: study values after clipping QC (n = {len(x):,})")
        ax.set_xlabel("attenuation_Mean, HU (median per study)")
        ax.set_ylabel("studies (log scale)")
        ax.legend(frameon=False)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        fig.tight_layout()
        fig.savefig(OUTPUT_DIR / f"mad_vs_sd_{organ}_{phase}.png", dpi=200)
        plt.close(fig)

print(f"\nSaved to {OUTPUT_DIR}: mad_vs_sd_summary.tsv, mad_vs_sd_totals.tsv, mad_vs_sd_<organ>_<phase>.png (4)")