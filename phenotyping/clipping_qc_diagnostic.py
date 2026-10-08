"""
Removes clipped liver/spleen attenuation values (CT brightness, in HU) from the
IDP data. A clipped value is one stuck at the image's storage limit rather than
real tissue. Saves the cleaned rows to idp_liver_spleen_clean.pkl.

Run cell by cell ("# %%"): Cells 1-2 build the raw cache once, Cell 10 applies
the filter.
"""

# %% Cell 0: imports and settings
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")  # remove this line in Jupyter if you want plots inline
import matplotlib.pyplot as plt

IDP_FILE = "/static/PMBB/PMBB-Release-2026-4.0/Imaging/CT/idp/PMBB-Release-2026-4.0_imaging_ct_idp.txt"
RAW_CACHE_FILE = Path("idp_liver_spleen_raw.pkl")
CLEAN_CACHE_FILE = Path("idp_liver_spleen_clean.pkl")

KEEP_COLS = [
    "person_id", "accession", "study_date", "anatomical_structure", "phase",
    "attenuation_Mean", "attenuation_Minimum", "attenuation_Maximum",
]
STRUCTURES = ["liver", "spleen"]

# Clipping cutoffs actually used by the filter -- update once you've decided.
MIN_CLIP_VALUES = [-32768]
MAX_CLIP_VALUES = [3071] #3071 was flagged by the data, sits at the CT scale ceiling (it makes the mean inflated about 5-6 HU higher)

# Known padding/saturation values, always reported in the diagnostic.
CANDIDATE_MIN_VALUES = [-32768]
CANDIDATE_MAX_VALUES = [32767] 

SPIKE_WINDOW = 15       # compare each value against this many neighbors on each side
SPIKE_MULTIPLIER = 15   # flag a value if its count > this x its LOCAL neighbor median
SPIKE_MIN_COUNT = 20    # ...and its count is at least this high

QC_OUTPUT_DIR = Path("qc_clipping")
QC_OUTPUT_DIR.mkdir(exist_ok=True)

RUN_SENSITIVITY_CHECK = False  # Cell 9b's 27-combination grid search -- turn on when you need it,
                                # off (default) to run the rest quickly without a batch job


# %% Cell 1 (run once): confirm column and structure names cheaply
header = pd.read_csv(IDP_FILE, sep="\t", nrows=0).columns.tolist()
print("Columns:", header)

missing = [c for c in KEEP_COLS if c not in header]
print("Missing from KEEP_COLS:", missing or "none")

names = pd.read_csv(IDP_FILE, sep="\t", usecols=["anatomical_structure", "phase"])
print(names["anatomical_structure"].value_counts().head(40))
print(names["phase"].value_counts())
# If liver/spleen are spelled differently (e.g. "Liver"), fix STRUCTURES in cell 0.


# %% Cell 2 (run once): read in chunks, keep liver/spleen only, save raw cache
# Skips the expensive read if the cache already exists -- safe to run as a
# plain top-to-bottom batch script without redoing this every time.
if not RAW_CACHE_FILE.exists():
    chunks = []
    for chunk in pd.read_csv(IDP_FILE, sep="\t", usecols=KEEP_COLS, chunksize=500_000):
        chunks.append(chunk[chunk["anatomical_structure"].isin(STRUCTURES)])
    idp = pd.concat(chunks, ignore_index=True)
    idp.to_pickle(RAW_CACHE_FILE)
    print(f"Saved {len(idp):,} liver/spleen rows to {RAW_CACHE_FILE}")
else:
    print(f"{RAW_CACHE_FILE} already exists -- skipping rebuild")


# %% Cell 3 (start here in later sessions): load the raw cache
idp = pd.read_pickle(RAW_CACHE_FILE)
print(f"Loaded {len(idp):,} rows")
print(idp.groupby(["anatomical_structure", "phase"]).size())


# %% Cell 4: helper functions
def detect_spikes(value_counts, window=SPIKE_WINDOW,
                   multiplier=SPIKE_MULTIPLIER, min_count=SPIKE_MIN_COUNT):
    """Flag values whose count is far higher than the LOCAL median count of
    nearby values (sorted numerically) -- isolates discrete artificial
    spikes (e.g. padding codes) from smooth natural variation in the bulk
    of the distribution, where moderately-high counts are normal."""
    if len(value_counts) < window * 2 + 1:
        return []
    vc_sorted = value_counts.sort_index()
    counts = vc_sorted.values
    values = vc_sorted.index.values
    n = len(counts)
    spikes = []
    for i in range(n):
        lo, hi = max(0, i - window), min(n, i + window + 1)
        neighbor_counts = np.concatenate([counts[lo:i], counts[i + 1:hi]])
        if len(neighbor_counts) == 0:
            continue
        local_med = np.median(neighbor_counts)
        threshold = max(local_med * multiplier, min_count)
        if counts[i] > threshold:
            spikes.append(values[i])
    return spikes


def summarize_mean_at_value(df, col, value, mean_col="attenuation_Mean"):
    """Compare attenuation_Mean for rows AT col==value vs all other rows."""
    def stats(s):
        if len(s) == 0:
            return dict(n=0, mean=np.nan, median=np.nan, sd=np.nan, iqr=np.nan)
        return dict(n=len(s), mean=s.mean(), median=s.median(), sd=s.std(),
                    iqr=s.quantile(0.75) - s.quantile(0.25))
    return stats(df.loc[df[col] == value, mean_col]), stats(df.loc[df[col] != value, mean_col])


def diagnose(idp, structure, phase, make_plots=True):
    """
    Clipping diagnostic for ONE structure/phase. Returns a list of summary
    rows, so you can run combinations one at a time and inspect each.
    """
    sub = idp[(idp["anatomical_structure"] == structure) & (idp["phase"] == phase)]
    print("=" * 70)
    print(f"{structure} / {phase}  (n={len(sub):,})")
    print("=" * 70)
    if len(sub) == 0:
        return []

    rows = []
    for edge, col, candidates in [
        ("min", "attenuation_Minimum", CANDIDATE_MIN_VALUES),
        ("max", "attenuation_Maximum", CANDIDATE_MAX_VALUES),
    ]:
        vc = sub[col].value_counts()
        print(f"\nTop 20 {col} values:")
        print(vc.head(20))

        spikes = detect_spikes(vc)
        check_values = sorted(set(spikes) | set(candidates))
        print(f"\nCandidate/spike values checked for {col}:")
        for value in check_values:
            count = int(vc.get(value, 0))
            pct = 100 * count / len(sub)
            tags = [t for t, f in [("SPIKE", value in spikes), ("candidate", value in candidates)] if f]
            print(f"  {value}: count={count:,} ({pct:.4f}%) [{', '.join(tags)}]")

            at, other = summarize_mean_at_value(sub, col, value)
            rows.append({
                "anatomical_structure": structure, "phase": phase, "edge": edge,
                "value": value, "count": count, "pct_of_rows": pct,
                "is_spike": value in spikes, "is_candidate": value in candidates,
                **{f"at_value_{k}": v for k, v in at.items()},
                **{f"other_{k}": v for k, v in other.items()},
            })
        #No collapsing done here these are raw rows
        if make_plots:
            fig, ax = plt.subplots(figsize=(8, 5))
            sub[col].plot.hist(bins=100, ax=ax, log=True)
            ax.set_title(f"{structure} / {phase} -- {col} (log y)")
            ax.set_xlabel(col)
            fig.tight_layout()
            fig.savefig(QC_OUTPUT_DIR / f"hist_{structure}_{phase}_{edge}.png", dpi=150)
            plt.close(fig)

            # Comparison: same columns, collapsed to one median value per
            # study (person_id, accession, checks whether this
            # shows up at the individual/study level too, or only in raw rows.
            study_level = sub.groupby(["person_id", "accession"])[col].median()
            fig, ax = plt.subplots(figsize=(8, 5))
            study_level.plot.hist(bins=100, ax=ax, log=True)
            ax.set_title(f"{structure} / {phase} -- {col} (log y, median per study)")
            ax.set_xlabel(col)
            fig.tight_layout()
            fig.savefig(QC_OUTPUT_DIR / f"hist_{structure}_{phase}_{edge}_study_median.png", dpi=150)
            plt.close(fig)

    return rows


# %% Cell 5: diagnostic, one combination at a time -- all phases, not just
# native/portal_venous, since clipping could affect any phase
summary_rows = []
summary_rows += diagnose(idp, "liver", "native")

# %% Cell 6
summary_rows += diagnose(idp, "spleen", "native")

# %% Cell 7
summary_rows += diagnose(idp, "liver", "portal_venous")

# %% Cell 8
summary_rows += diagnose(idp, "spleen", "portal_venous")

# %% Cell 8b
summary_rows += diagnose(idp, "liver", "arterial_early")
summary_rows += diagnose(idp, "spleen", "arterial_early")

# %% Cell 8c
summary_rows += diagnose(idp, "liver", "arterial_late")
summary_rows += diagnose(idp, "spleen", "arterial_late")

# %% Cell 9: save the diagnostic summary table
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(QC_OUTPUT_DIR / "clipping_qc_summary.tsv", sep="\t", index=False)
summary_df


# %% Cell 9b: sensitivity check -- does the spike result depend on the settings?
# Re-runs the spike detector across a grid of reasonable settings and counts how
# often each value gets flagged. Values flagged under nearly every setting are
# robust findings; values flagged only under loose settings are borderline.
# Gated by RUN_SENSITIVITY_CHECK (Cell 0) -- this is the expensive part (27x
# the spike detector), off by default so the rest of the script runs quickly.
if RUN_SENSITIVITY_CHECK:
    import itertools

    WINDOWS = [10, 15, 25]
    MULTIPLIERS = [10, 15, 25]
    MIN_COUNTS = [10, 20, 50]
    settings = list(itertools.product(WINDOWS, MULTIPLIERS, MIN_COUNTS))  # 27 combinations

    records = []
    for (structure, phase), sub in idp.groupby(["anatomical_structure", "phase"]):
        for edge, col in [("min", "attenuation_Minimum"), ("max", "attenuation_Maximum")]:
            vc = sub[col].value_counts()
            for window, mult, min_count in settings:
                for value in detect_spikes(vc, window=window, multiplier=mult, min_count=min_count):
                    records.append({
                        "anatomical_structure": structure, "phase": phase, "edge": edge,
                        "value": value, "count": int(vc[value]),
                        "window": window, "multiplier": mult, "min_count": min_count,
                    })

    sens = pd.DataFrame(records)

    if sens.empty:
        print("No spikes flagged under any setting.")
    else:
        stability = (
            sens.groupby(["anatomical_structure", "phase", "edge", "value", "count"])
            .size()
            .rename("n_settings_flagged")
            .reset_index()
        )
        stability["pct_settings_flagged"] = 100 * stability["n_settings_flagged"] / len(settings)
        stability = stability.sort_values("pct_settings_flagged", ascending=False)

        stability.to_csv(QC_OUTPUT_DIR / "spike_sensitivity.tsv", sep="\t", index=False)
        print(stability.head(40).to_string(index=False))

        robust = stability[stability["pct_settings_flagged"] >= 80]
        print(f"\n{len(robust)} value(s) flagged under >=80% of settings combinations (robust findings):")
        print(robust[["anatomical_structure", "phase", "edge", "value", "count"]].to_string(index=False))

        # Run the mean comparison on every robust value, so each one gets the confirmation step
        col_for_edge = {"min": "attenuation_Minimum", "max": "attenuation_Maximum"}
        check = []
        for _, r in robust.iterrows():
            sub = idp[(idp["anatomical_structure"] == r["anatomical_structure"]) &
                      (idp["phase"] == r["phase"])]
            at, other = summarize_mean_at_value(sub, col_for_edge[r["edge"]], r["value"])
            check.append({
                "anatomical_structure": r["anatomical_structure"], "phase": r["phase"],
                "edge": r["edge"], "value": r["value"],
                "at_value_n": at["n"], "at_value_median_mean": at["median"],
                "other_median_mean": other["median"],
            })
        robust_check = pd.DataFrame(check)
        robust_check.to_csv(QC_OUTPUT_DIR / "robust_spike_mean_check.tsv", sep="\t", index=False)
        print(robust_check.to_string(index=False))
else:
    print("RUN_SENSITIVITY_CHECK is False -- skipping Cell 9b.")


# %% Cell 10: apply the decided clipping filter, save a CLEANED cache for
# other scripts to consume (edit MIN/MAX_CLIP_VALUES in cell 0 first)
min_clip_mask = pd.Series(False, index=idp.index)
for v in MIN_CLIP_VALUES:
    min_clip_mask |= idp["attenuation_Minimum"] <= v

max_clip_mask = pd.Series(False, index=idp.index)
for v in MAX_CLIP_VALUES:
    max_clip_mask |= idp["attenuation_Maximum"] >= v

clip_mask = min_clip_mask | max_clip_mask
idp_clean = idp[~clip_mask]

n_dropped = clip_mask.sum()
print(f"Dropped {n_dropped:,} / {len(idp):,} rows ({100 * n_dropped / len(idp):.2f}%)")
print(idp.assign(min_edge=min_clip_mask, max_edge=max_clip_mask)
         .loc[clip_mask]
         .groupby(["anatomical_structure", "phase"])[["min_edge", "max_edge"]].sum())

idp_clean.to_pickle(CLEAN_CACHE_FILE)
print(f"Saved {len(idp_clean):,} QC-cleaned rows to {CLEAN_CACHE_FILE}")


# %% Cell 11: attenuation_Maximum AFTER clipping (confirms nothing above
# MAX_CLIP_VALUES remains) -- liver/spleen x native/portal_venous only
for structure in ["liver", "spleen"]:
    for phase in ["native", "portal_venous"]:
        sub_clean = idp_clean[
            (idp_clean["anatomical_structure"] == structure) & (idp_clean["phase"] == phase)
        ]
        if len(sub_clean) == 0:
            continue

        study_level_clean = sub_clean.groupby(["person_id", "accession"])["attenuation_Maximum"].median()

        fig, ax = plt.subplots(figsize=(8, 5))
        study_level_clean.plot.hist(bins=100, ax=ax, log=True)
        ax.set_title(f"{structure} / {phase} -- attenuation_Maximum, post-QC (log y, median per study)")
        ax.set_xlabel("attenuation_Maximum")
        fig.tight_layout()
        fig.savefig(QC_OUTPUT_DIR / f"hist_{structure}_{phase}_max_post_qc.png", dpi=150)
        plt.close(fig)
        print(f"{structure} / {phase}: max study-level attenuation_Maximum after QC = "
              f"{study_level_clean.max():.1f}")
