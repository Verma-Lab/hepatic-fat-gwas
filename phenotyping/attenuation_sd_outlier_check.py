#!/usr/bin/env python3
"""
Checks where a +/-N SD outlier filter on attenuation_Mean would actually
fall on the real distribution, before committing to a specific threshold

Every distribution here is collapsed to one median attenuation_Mean value
per study per phase (person_id, accession, phase) -- never raw rows.
Compares pre-QC (raw cache) vs post-QC (clipping-filtered cache), and
native/PV separately vs combined (pooled, illustration only). Also reports
skew after SD trimming and flags implausible studies (flag only, no exclusion).

Run with:
    python3 attenuation_sd_outlier_check.py
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

RAW_CACHE_FILE = Path("idp_liver_spleen_raw.pkl")
CLEAN_CACHE_FILE = Path("idp_liver_spleen_clean.pkl")
OUTPUT_DIR = Path("qc_clipping/attenuation_sd_outlier_check")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SD_THRESHOLDS = [3, 4, 5, 6, 7]

# Plausibility limits (rough, not published standards) -- flag only, nothing excluded.
NATIVE_SPLEEN_MAX_HU = 100   # native spleen study-level value above this is flagged
NATIVE_LIVER_MAX_HU = 100    # native liver study-level value above this is flagged
ANY_STUDY_MAX_HU = 300       # any organ/phase study-level value above this is flagged
AIR_MASK_MAX_HU = -900       # row attenuation_Maximum below this = mask entirely air
DENSE_MASK_MIN_HU = 0        # native spleen row attenuation_Minimum above this = whole mask dense

ZOOM_LOWER_Q, ZOOM_UPPER_Q = 0.001, 0.999   # zoomed plots: 0.1th-99.9th percentile
ZOOM_MARGIN_FRAC = 0.05                     # ...plus this fraction of the span on each side

PHASES = ["native", "portal_venous"]
STRUCTURES = ["liver", "spleen"]


def study_level_means(df, structure, phase=None):
    """One median attenuation_Mean per study per phase (person_id, accession, phase).
    phase=None pools native + portal_venous together (the "combined" view); a study
    with both phases contributes one value per phase, never a blended median."""
    sub = df[df["anatomical_structure"] == structure]
    sub = sub[sub["phase"] == phase] if phase is not None else sub[sub["phase"].isin(["native", "portal_venous"])]
    return (
        sub.groupby(["person_id", "accession", "phase"])["attenuation_Mean"]
        .median()
        .reset_index()["attenuation_Mean"]
    )


# --- Plot helpers: SD/mean lines, log-y and zoomed variants -------------------
def sd_bounds(study_means, n_sd):
    mean_val, sd_val = study_means.mean(), study_means.std()
    return mean_val - n_sd * sd_val, mean_val + n_sd * sd_val


def draw_mean_and_sd_lines(ax, mean_val, sd_val):
    ax.axvline(mean_val, color="black", linestyle="-", linewidth=1, label="mean")
    for n_sd, color in [(4, "orange"), (6, "red")]:
        ax.axvline(mean_val - n_sd * sd_val, color=color, linestyle="--", label=f"±{n_sd} SD")
        ax.axvline(mean_val + n_sd * sd_val, color=color, linestyle="--")


def zoom_xlim(study_means):
    """0.1th-99.9th percentile, widened if needed so the +/-6 SD lines stay visible."""
    lo, hi = study_means.quantile(ZOOM_LOWER_Q), study_means.quantile(ZOOM_UPPER_Q)
    lo6, hi6 = sd_bounds(study_means, 6)
    lo, hi = min(lo, lo6), max(hi, hi6)
    margin = ZOOM_MARGIN_FRAC * (hi - lo)
    return lo - margin, hi + margin


def save_histogram(study_means, title, path, log=False, zoom=False):
    mean_val, sd_val = study_means.mean(), study_means.std()
    fig, ax = plt.subplots(figsize=(8, 5))
    if zoom:
        xlim = zoom_xlim(study_means)
        n_outside = ((study_means < xlim[0]) | (study_means > xlim[1])).sum()
        study_means.plot.hist(bins=100, range=xlim, ax=ax, log=log)
        ax.set_xlim(xlim)
        title += f"\nzoomed: {n_outside:,} studies outside view"
    else:
        study_means.plot.hist(bins=100, ax=ax, log=log)
    draw_mean_and_sd_lines(ax, mean_val, sd_val)
    ax.set_title(title)
    ax.set_xlabel("attenuation_Mean (median per study)")
    ax.set_ylabel("studies (log scale)" if log else "studies")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def report_and_plot(study_means, label, filename_stub, combined=False):
    if len(study_means) == 0:
        print(f"{label}: no data")
        return

    mean_val, sd_val = study_means.mean(), study_means.std()
    skew_val = study_means.skew()
    if abs(skew_val) < 0.5:
        skew_tag = "approximately symmetric"
    elif abs(skew_val) < 1:
        skew_tag = "moderately skewed"
    else:
        skew_tag = "highly skewed"
    print("=" * 70)
    print(f"{label} (n={len(study_means):,}, mean={mean_val:.2f}, SD={sd_val:.2f}, "
          f"skew={skew_val:.3f} [{skew_tag}])")
    print("=" * 70)

    if combined:
        # Combined native+PV is bimodal (phase shift), so its SD is inflated by the
        # gap between phases -- no exclusion counts from this view.
        print("  NOTE: combined native+PV view is for illustration only. SD cutoffs should be")
        print("  applied per organ AND phase; no exclusion recommendations computed here.")
    else:
        for n_sd in SD_THRESHOLDS:
            lo, hi = mean_val - n_sd * sd_val, mean_val + n_sd * sd_val
            n_excluded = ((study_means < lo) | (study_means > hi)).sum()
            pct_excluded = 100 * n_excluded / len(study_means)
            print(f"  +/-{n_sd} SD -> range [{lo:.1f}, {hi:.1f}]: "
                  f"{n_excluded:,} studies excluded ({pct_excluded:.3f}%)")

    title = f"{label}\n(median per study, skew={skew_val:.3f} [{skew_tag}])"
    if combined:
        title += " -- SD lines illustrative only"
    save_histogram(study_means, title, OUTPUT_DIR / f"{filename_stub}.png")
    save_histogram(study_means, title, OUTPUT_DIR / f"{filename_stub}_logy.png", log=True)
    save_histogram(study_means, title, OUTPUT_DIR / f"{filename_stub}_zoom.png", zoom=True)
    print()


idp_raw = pd.read_pickle(RAW_CACHE_FILE)
idp_clean = pd.read_pickle(CLEAN_CACHE_FILE)


# --- Bug-fix check: study counts with old vs fixed groupby --------------------
# Old groupby was (person_id, accession); fixed adds phase. Per-phase n must be
# identical; combined n grows by the number of studies that have both phases.
def n_studies_old_groupby(df, structure, phase=None):
    sub = df[df["anatomical_structure"] == structure]
    sub = sub[sub["phase"] == phase] if phase is not None else sub[sub["phase"].isin(["native", "portal_venous"])]
    return sub.groupby(["person_id", "accession"]).ngroups


print("=" * 70)
print("Bug-fix check: n studies, old groupby (person_id, accession) vs fixed (+ phase)")
print("=" * 70)
for qc_label, idp_df in [("pre_qc", idp_raw), ("post_qc", idp_clean)]:
    for structure in STRUCTURES:
        for phase in PHASES:
            n_old = n_studies_old_groupby(idp_df, structure, phase)
            n_new = len(study_level_means(idp_df, structure, phase))
            status = "OK, unchanged" if n_old == n_new else "MISMATCH"
            print(f"  {qc_label:8s} {structure:6s} {phase:14s} before={n_old:,}  after={n_new:,}  [{status}]")
        n_old = n_studies_old_groupby(idp_df, structure, None)
        n_new = len(study_level_means(idp_df, structure, None))
        print(f"  {qc_label:8s} {structure:6s} {'combined':14s} before={n_old:,}  after={n_new:,}  "
              f"(+{n_new - n_old:,} = studies with both native and PV, previously blended)")
print()


for qc_label, idp_df in [("pre_qc", idp_raw), ("post_qc", idp_clean)]:
    for structure in ["liver", "spleen"]:
        for phase in ["native", "portal_venous"]:
            sm = study_level_means(idp_df, structure, phase)
            report_and_plot(
                sm, f"{structure} / {phase} ({qc_label})",
                f"sd_outlier_{structure}_{phase}_{qc_label}",
            )

        sm_combined = study_level_means(idp_df, structure, phase=None)
        report_and_plot(
            sm_combined, f"{structure} / combined native+PV ({qc_label})",
            f"sd_outlier_{structure}_combined_{qc_label}", combined=True,
        )

post_qc_study_means = {
    (structure, phase): study_level_means(idp_clean, structure, phase)
    for structure in STRUCTURES for phase in PHASES
}


# --- Slide figure: 2x2 post-QC grid, log y, per organ x phase -----------------
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
for i, structure in enumerate(STRUCTURES):
    for j, phase in enumerate(PHASES):
        ax = axes[i, j]
        sm = post_qc_study_means[(structure, phase)]
        if len(sm) == 0:
            ax.set_title(f"{structure.capitalize()}, {phase.replace('_', ' ')}: no data")
            continue
        sm.plot.hist(bins=100, ax=ax, log=True, color="steelblue")
        draw_mean_and_sd_lines(ax, sm.mean(), sm.std())
        ax.set_title(f"{structure.capitalize()}, {phase.replace('_', ' ')} (n={len(sm):,})")
        ax.set_xlabel("attenuation_Mean, HU (median per study)")
        ax.set_ylabel("studies (log scale)")
handles, labels = axes[0, 0].get_legend_handles_labels()
handles, labels = zip(*[(h, l) for h, l in zip(handles, labels) if l != "attenuation_Mean"])
fig.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False)
fig.suptitle("Study-level attenuation, post clipping QC, with ±4 / ±6 SD cutoffs")
fig.tight_layout(rect=[0, 0.05, 1, 0.96])
fig.savefig(OUTPUT_DIR / "sd_outlier_grid_post_qc.png", dpi=200)
plt.close(fig)


# --- Skewness after SD trimming (post-QC, per organ x phase) -----------------
# SD bounds are computed once from all studies, same as the exclusion counts above.
print("=" * 70)
print("Skew after SD trimming (post-QC, median per study)")
print("=" * 70)
skew_rows = []
for (structure, phase), sm in post_qc_study_means.items():
    if len(sm) == 0:
        continue
    for threshold in ["all", 6, 4]:
        if threshold == "all":
            kept = sm
        else:
            lo, hi = sd_bounds(sm, threshold)
            kept = sm[(sm >= lo) & (sm <= hi)]
        skew_rows.append({
            "organ": structure, "phase": phase,
            "threshold": "all" if threshold == "all" else f"+/-{threshold}SD",
            "n_kept": len(kept), "n_removed": len(sm) - len(kept),
            "mean": kept.mean(), "sd": kept.std(), "skew": kept.skew(),
        })
skew_df = pd.DataFrame(skew_rows)
print(skew_df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
skew_df.to_csv(OUTPUT_DIR / "skew_by_threshold.tsv", sep="\t", index=False)
print()


# --- Plausibility flags (post-QC, flag only -- nothing excluded) --------------
# One row per (study, phase, rule); a study can be flagged by more than one rule.
sub = idp_clean[idp_clean["phase"].isin(PHASES)].copy()
sub["row_air"] = sub["attenuation_Maximum"] < AIR_MASK_MAX_HU
sub["row_dense"] = sub["attenuation_Minimum"] > DENSE_MASK_MIN_HU
studies = (
    sub.groupby(["person_id", "accession", "anatomical_structure", "phase"])
    .agg(study_value=("attenuation_Mean", "median"),
         n_rows=("attenuation_Mean", "size"),
         n_rows_air=("row_air", "sum"),
         n_rows_dense=("row_dense", "sum"))
    .reset_index()
    .rename(columns={"anatomical_structure": "organ"})
)

outside_6sd = pd.Series(False, index=studies.index)
for (structure, phase), sm in post_qc_study_means.items():
    if len(sm) == 0:
        continue
    lo, hi = sd_bounds(sm, 6)
    in_group = (studies["organ"] == structure) & (studies["phase"] == phase)
    outside_6sd |= in_group & ((studies["study_value"] < lo) | (studies["study_value"] > hi))
studies["outside_6sd"] = outside_6sd

is_native = studies["phase"] == "native"
rules = {
    f"native_spleen_gt_{NATIVE_SPLEEN_MAX_HU}HU":
        is_native & (studies["organ"] == "spleen") & (studies["study_value"] > NATIVE_SPLEEN_MAX_HU),
    f"native_liver_gt_{NATIVE_LIVER_MAX_HU}HU":
        is_native & (studies["organ"] == "liver") & (studies["study_value"] > NATIVE_LIVER_MAX_HU),
    f"any_gt_{ANY_STUDY_MAX_HU}HU":
        studies["study_value"] > ANY_STUDY_MAX_HU,
    f"air_mask_row_max_lt_{AIR_MASK_MAX_HU}HU":
        studies["n_rows_air"] > 0,
    f"native_spleen_dense_mask_row_min_gt_{DENSE_MASK_MIN_HU}HU":
        is_native & (studies["organ"] == "spleen") & (studies["n_rows_dense"] > 0),
}
flags = pd.concat(
    [studies[mask].assign(rule=rule) for rule, mask in rules.items()], ignore_index=True
)
flag_cols = ["person_id", "accession", "organ", "phase", "study_value", "rule",
             "n_rows", "n_rows_air", "n_rows_dense", "outside_6sd"]
flags = flags[flag_cols]
flags.to_csv(OUTPUT_DIR / "plausibility_flags.tsv", sep="\t", index=False)

print("=" * 70)
print("Plausibility flags (post-QC, median per study; flag only, nothing excluded)")
print("=" * 70)
print("Air-mask / dense-mask rules flag a study if ANY of its rows meets the condition")
print("(see n_rows_air / n_rows_dense vs n_rows in the TSV).\n")
print("Per rule (studies flagged | also outside +/-6 SD | caught ONLY by this rule's limits):")
for rule in rules:
    r = flags[flags["rule"] == rule]
    print(f"  {rule:45s} {len(r):6,} | {r['outside_6sd'].sum():6,} | {(~r['outside_6sd']).sum():6,}")
print()
print("Per rule x organ/phase:")
if len(flags):
    print(flags.groupby(["rule", "organ", "phase"]).size().rename("n_flagged").to_string())
else:
    print("  none flagged")
print()

# Unique studies (any rule) vs +/-6 SD, per organ/phase
unique_flagged = flags.drop_duplicates(["person_id", "accession", "organ", "phase"])
print("Unique flagged studies (any rule) vs +/-6 SD, per organ/phase:")
for structure in STRUCTURES:
    for phase in PHASES:
        u = unique_flagged[(unique_flagged["organ"] == structure) & (unique_flagged["phase"] == phase)]
        n_sd6 = ((studies["organ"] == structure) & (studies["phase"] == phase) & studies["outside_6sd"]).sum()
        print(f"  {structure:6s} {phase:14s} flagged={len(u):,}  also outside 6SD={u['outside_6sd'].sum():,}  "
              f"rules only (6SD misses)={(~u['outside_6sd']).sum():,}  [total outside 6SD={n_sd6:,}]")
print(f"\nTotal unique flagged: {len(unique_flagged):,} "
      f"(also 6SD: {unique_flagged['outside_6sd'].sum():,}, "
      f"rules only: {(~unique_flagged['outside_6sd']).sum():,})")


