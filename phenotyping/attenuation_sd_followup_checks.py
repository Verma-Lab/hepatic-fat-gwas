#!/usr/bin/env python3
"""
Follow-up checks on the ±6 SD cutoff, per anatomical structure and phase.

Not part of attenuation_sd_outlier_check.py. These are the extra checks to decide 
between ±4 and ±6 SD and to justify the native spleen HU range check (>100 HU).

For each organ + phase (study level = median attenuation_Mean across a
study's series, after clipping QC):
  1. Where ±4 and ±6 SD fall, and which studies ±4 SD would remove that
     ±6 SD keeps (low vs high side).
  2. Same-person check: are flagged studies one-offs (the person's other
     studies in that phase are normal) or consistent for that person?
  3. Median impact: how much one-off studies move a person's value
     (person value = median across their studies), in HU and percentile.
  4. Phase probability: are extreme studies likely mislabeled phases?
     (reads phase_probability from the raw IDP file; not in the clean cache)
  5. Dates: for flagged studies, how far away the person's nearest other
     study is, its value, and whether a contrast study is nearby.

Results from 10/7/26 (for reference):
  - Native liver: ±4 SD would remove 36 more, all low (-19 to +2 HU, ~37-49%
    fat) -> likely real severe fat -> keep ±6 SD.
  - PV liver: ±4 SD would remove 27 more (23 high, 4 low). High ones are people
    bright in every PV study (spleen bright too); phase_probability = 1.0.
    3 one-off low studies shift person values by only 0.1-3.3 HU -> keep ±6 SD.
  - Native spleen: ±6 SD cuts above ~155 HU; 64 studies at 100-153 HU are
    kept by ±6 SD. 31 of 48 with other native studies are normal (~30-60 HU)
    in those -> range check >100 HU.
  - PV spleen: ±6 SD cuts above ~325 HU; >300 check would add only 27 (not
    clear errors) -> no range check. ±4 SD would remove 188 more (~278 HU vs
    ~127 in the same people's other PV studies) -> consider ±4 SD if the PV
    phenotype uses spleen.

Run inside an interactive session (not the login node). Cells are marked
with "# %%" so they can be run one at a time.
"""

# %% Cell 0: imports and settings
import pandas as pd
import numpy as np
from pathlib import Path

IDP_FILE = "/static/PMBB/PMBB-Release-2026-4.0/Imaging/CT/idp/PMBB-Release-2026-4.0_imaging_ct_idp.txt"
CLEAN_CACHE = "idp_liver_spleen_clean.pkl"   # clipping-QC'd rows (from the clipping script)

OUT_DIR = Path("qc_clipping/attenuation_sd_followup")
OUT_DIR.mkdir(parents=True, exist_ok=True)

K_STRICT, K_LOOSE = 4, 6
COMBOS = [("liver", "native"), ("liver", "portal_venous"),
          ("spleen", "native"), ("spleen", "portal_venous")]
RANGE_CUTS = {("spleen", "native"): 100}       # HU range check (native spleen only)

pd.set_option("display.width", 200)


# %% Cell 1: load the clipping-QC'd cache
idp_clean = pd.read_pickle(CLEAN_CACHE)
idp_clean["study_date"] = pd.to_datetime(idp_clean["study_date"])
print(f"Loaded {len(idp_clean):,} rows")


# %% Cell 2: helper functions
def study_values(df, organ, phase):
    """One value per study: median attenuation_Mean across its series."""
    d = df[(df["anatomical_structure"] == organ) & (df["phase"] == phase)]
    return d.groupby(["person_id", "accession"])["attenuation_Mean"].median()


def sd_bounds(sm, k):
    m, s = sm.mean(), sm.std()
    return m - k * s, m + k * s


def cutoff_groups(sm, range_cut=None):
    """Studies removed by ±6 SD, extra studies ±4 SD would remove, and range-check studies."""
    lo4, hi4 = sd_bounds(sm, K_STRICT)
    lo6, hi6 = sd_bounds(sm, K_LOOSE)
    groups = {
        "removed_by_6sd": sm[(sm < lo6) | (sm > hi6)],
        "extra_low_4sd_only": sm[(sm < lo4) & (sm >= lo6)],
        "extra_high_4sd_only": sm[(sm > hi4) & (sm <= hi6)],
    }
    if range_cut is not None:
        groups[f"range_check_gt{range_cut}_inside_6sd"] = sm[(sm > range_cut) & (sm >= lo6) & (sm <= hi6)]
    return groups, (lo4, hi4, lo6, hi6)


def other_studies(sm, ids):
    """For each flagged study: the same person's OTHER studies in this phase (count + median)."""
    by_person = sm.reset_index()
    rows = []
    for pid, acc in ids:
        o = by_person[(by_person["person_id"] == pid) & (by_person["accession"] != acc)]["attenuation_Mean"]
        rows.append({"person_id": pid, "accession": acc, "value": sm.loc[(pid, acc)],
                     "n_other": len(o), "other_median": o.median() if len(o) else np.nan})
    return pd.DataFrame(rows)


def median_impact(sm, ids):
    """Person value (median across studies) with vs without the flagged studies, in HU and percentile."""
    ids = pd.MultiIndex.from_tuples(list(ids), names=["person_id", "accession"])
    people = ids.get_level_values("person_id").unique()
    with_ = sm.groupby(level="person_id").median()
    without = sm.drop(ids).groupby(level="person_id").median()
    res = pd.DataFrame({"with": with_.reindex(people), "without": without.reindex(people)})
    res["diff_HU"] = res["with"] - res["without"]
    res["pct_with"] = [100 * (with_ < v).mean() for v in res["with"]]
    res["pct_without"] = [100 * (with_ < v).mean() if pd.notna(v) else np.nan for v in res["without"]]
    return res


def phase_probability(organ, phase):
    """Study-level phase_probability from the raw IDP file (slow: scans the whole file)."""
    cols = ["person_id", "accession", "anatomical_structure", "phase", "phase_probability"]
    chunks = []
    for chunk in pd.read_csv(IDP_FILE, sep="\t", usecols=cols, chunksize=500_000):
        chunks.append(chunk[(chunk["anatomical_structure"] == organ) & (chunk["phase"] == phase)])
    pp = pd.concat(chunks, ignore_index=True)
    return pp.groupby(["person_id", "accession"])["phase_probability"].median()


def date_context(df, organ, phase, ids):
    """Nearest other study (days, value), other studies before/after, and nearby contrast studies."""
    d = df[(df["anatomical_structure"] == organ) & (df["phase"] == phase)]
    st = (d.groupby(["person_id", "accession"])
            .agg(value=("attenuation_Mean", "median"), date=("study_date", "min"))
            .reset_index())
    con = df[df["phase"] != "native"][["person_id", "accession", "study_date"]].drop_duplicates()
    rows = []
    for pid, acc in ids:
        f = st[(st["person_id"] == pid) & (st["accession"] == acc)].iloc[0]
        o = st[(st["person_id"] == pid) & (st["accession"] != acc)]
        c = con[con["person_id"] == pid]
        days_before = (f["date"] - c.loc[c["accession"] != acc, "study_date"]).dt.days
        rec = {"person_id": pid, "accession": acc, "value": f["value"], "n_other": len(o),
               "contrast_same_accession": (c["accession"] == acc).any(),
               "contrast_other_acc_0to3d_before": ((days_before >= 0) & (days_before <= 3)).any()}
        if len(o):
            gaps = (o["date"] - f["date"]).dt.days
            i = gaps.abs().idxmin()
            rec.update(nearest_days=gaps[i], nearest_value=o.loc[i, "value"],
                       before_median=o.loc[gaps < 0, "value"].median(),
                       after_median=o.loc[gaps > 0, "value"].median())
        rows.append(rec)
    return pd.DataFrame(rows)


def mask_ids(df):
    """Replace person_id / accession with P1, S1, ... for printing (same ID -> same code)."""
    out = df.copy()
    for col, prefix in [("person_id", "P"), ("accession", "S")]:
        if col in out.columns:
            out[col] = prefix + (out[col].astype("category").cat.codes + 1).astype(str)
    return out


# %% Cell 3: checks 1 + 2 for every organ and phase (±4 vs ±6 SD, same-person check)
study = {}
for organ, phase in COMBOS:
    sm = study_values(idp_clean, organ, phase)
    study[(organ, phase)] = sm
    groups, (lo4, hi4, lo6, hi6) = cutoff_groups(sm, RANGE_CUTS.get((organ, phase)))
    print(f"\n=== {organ} {phase}: n={len(sm):,}  mean {sm.mean():.1f}  SD {sm.std():.1f}")
    print(f"±4 SD: [{lo4:.1f}, {hi4:.1f}]   ±6 SD: [{lo6:.1f}, {hi6:.1f}]")
    for name, g in groups.items():
        if len(g) == 0:
            print(f"  {name}: 0")
            continue
        o = other_studies(sm, g.index)
        has = o["n_other"] > 0
        print(f"  {name}: {len(g)} studies, {g.min():.1f} to {g.max():.1f} HU (median {g.median():.1f}); "
              f"{has.sum()} have other studies, typical other-study median {o.loc[has, 'other_median'].median():.1f}")
        o.to_csv(OUT_DIR / f"{organ}_{phase}_{name}.tsv", sep="\t", index=False)   # real IDs, stays on cluster


# %% Cell 4: check 3, median impact of one-off studies (example: PV liver extra low, minus consistently low people)
# Edit `ids` to flagged studies to test.
sm = study[("liver", "portal_venous")]
groups, _ = cutoff_groups(sm)
o = other_studies(sm, groups["extra_low_4sd_only"].index)
one_off = o[o["other_median"] > 60]                 # person is normal in other studies -> one-off
ids = pd.MultiIndex.from_frame(one_off[["person_id", "accession"]])
res = median_impact(sm, ids)
res.index = [f"P{i + 1}" for i in range(len(res))]
print(res.round(1))


# %% Cell 5: check 4, phase probability of extreme studies (example: PV liver; slow)
sm = study[("liver", "portal_venous")]
pp = phase_probability("liver", "portal_venous")
lo4, hi4 = sd_bounds(sm, K_STRICT)
high_ids, low_ids = sm[sm > hi4].index, sm[sm < lo4].index
other_ids = sm.index.difference(high_ids.union(low_ids))
for label, idx in [("High (>4 SD)", high_ids), ("Low (<-4 SD)", low_ids), ("All others", other_ids)]:
    print(label, pp.reindex(idx).describe().round(3).to_dict())


# %% Cell 6: check 5, dates for the native spleen range-check studies (100 HU to +6 SD)
sm = study[("spleen", "native")]
cut = RANGE_CUTS[("spleen", "native")]
groups, _ = cutoff_groups(sm, cut)
flag_ids = groups[f"range_check_gt{cut}_inside_6sd"].index
if len(flag_ids) == 0:
    print("No studies in the range-check group")
else:
    dc = date_context(idp_clean, "spleen", "native", flag_ids)
    dc.to_csv(OUT_DIR / "spleen_native_range_check_dates.tsv", sep="\t", index=False)
    print(mask_ids(dc).sort_values("value").round(1).to_string(index=False))
    print("\nNo other native study:", (dc["n_other"] == 0).sum(), "of", len(dc))