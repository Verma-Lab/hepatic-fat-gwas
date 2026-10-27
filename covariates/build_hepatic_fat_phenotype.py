#!/usr/bin/env python3
"""
Builds the SAIGE phenotype + covariate file for the hepatic fat GWAS from
PMBB Release 2026-4.0 CT IDPs: final QC, one value per person, native + PV
combined with a phase covariate. Replaces create_saige_covariate_file.py.

Writes one file per phenotype x transform x ancestry (EUR, AFR).

Run with (writes all 8 files):
    python3 build_hepatic_fat_phenotype.py
Rebuild one phenotype and/or transform (still writes both EUR and AFR):
    python3 build_hepatic_fat_phenotype.py --phenotype SL --transform untransformed
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# SETTINGS
# ============================================================================
CLEAN_CACHE_FILE = Path("/project/verma_shared/projects/Liver_IDPs/phenotyping/idp_liver_spleen_clean.pkl")
FAM_FILE = "/static/PMBB/PMBB-Release-2026-4.0/Imputed/chunked_bed_files/PMBB-Release-2026-4.0_genetic_imputed.463_chr8_chunk1_149284_2198252.fam"
COVARIATE_FILE = "/static/PMBB/PMBB-Release-2026-4.0/Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_covariates.txt"
CONDITION_FILE = "/static/PMBB/PMBB-Release-2026-4.0/Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_condition_occurrence.txt"
PCA_DIR = "/static/PMBB/PMBB-Release-2026-4.0/Imputed/PCA/by_ancestry"
ANCESTRIES = ["EUR", "AFR"]
EIGENVECTOR_FILES = {anc: f"{PCA_DIR}/PMBB-Release-2026-4.0_genetic_imputed.no1KG.eigenvectors.{anc}.tsv"
                     for anc in ANCESTRIES}
OUTPUT_DIR = Path("/project/verma_shared/projects/Liver_IDPs/covariates")

PHASE_PRIORITY = "native"            # person has both phases -> use this one
PHASES = ["native", "portal_venous"]
SD_CUTOFF = 6
NATIVE_SPLEEN_MAX_HU = 100
RANKNORM_K = 0.375
# Genetic PCs: PC1-5 from the ancestry-specific imputed PCAs, for EUR and AFR (decided 2026-10-12).
# Each GWAS runs in one ancestry; scree elbows at ~PC5 in both; PMBB README advises <= 5
# imputed PCs; Park et al. 2022 used ancestry-specific PCs. See CLAUDE.md "Genetic PC decision".
PCS = [f"PC{i}" for i in range(1, 6)]
# Covariate file columns (2026-4.0 uses lowercase); written out as Sample_age, Sex, batch.
# batch = genotyping batch (1/2/3), a categorical covariate (decided 2026-10-12).
COVARIATE_COLS = ["person_id", "batch", "sequenced_gender", "sample_age"]
PHENO_COL = "hepatic_fat"            # must match the SAIGE config
RAW_COL = "hepatic_fat_raw_HU"

EXCLUSION_CODES = [
    'B18.0', 'B18.1', 'B18.2', '070.32', '070.21', '070.22', '070.23',
    '070.31', '070.33', '070.54', '571.0', 'K70.0', '571.1', 'K70.1',
    '571.2', 'K70.3', 'K70.2', '571.3', 'K70.4', 'K70.40', 'K70.41',
    'K70.9', '303.0', '303.9', 'F10.229', 'F10.20'
]

# Expected QC counts (confirmed on cluster, Oct 2026)
EXPECTED_SD_REMOVED = {("liver", "native"): 8, ("liver", "portal_venous"): 6,
                       ("spleen", "native"): 117, ("spleen", "portal_venous"): 43}
EXPECTED_RANGE_REMOVED = 64
# Genotyped people with >=1 study (native or PV), before one phase per person
SANITY_TARGETS = {"SL": 19355, "liver": 19447}
SANITY_TOLERANCE = 0.05


# ============================================================================
# RANK-BASED INVERSE NORMAL TRANSFORMATION
# ============================================================================
def rank_norm(x, k=RANKNORM_K):
    """Same formula as RNOmni::RankNorm and the lab's RankNorm in generate_qc_idps.R (confirmed)."""
    r = rankdata(x, method="average")
    return norm.ppf((r - k) / (len(r) - 2 * k + 1))


def check_rank_norm():
    # Reference from R: qnorm((rank(x, ties.method="average") - 0.375) / (length(x) - 2*0.375 + 1))
    x = np.array([3.2, 1.5, 1.5, 7.0, -2.1, 4.4, 4.4, 4.4, 0.0, 10.0])
    r_values = np.array([-0.1225808439, -0.5104216427, -0.5104216427, 1.0004905456, -1.5466352714,
                         0.3754617702, 0.3754617702, 0.3754617702, -1.0004905456, 1.5466352714])
    max_diff = np.abs(rank_norm(x) - r_values).max()
    if max_diff > 1e-8:
        sys.exit(f"RankNorm check FAILED: max abs diff vs R = {max_diff:.2e}")
    print(f"RankNorm check vs R (10 values, with ties): max abs diff {max_diff:.1e} -- OK\n")


# ============================================================================
# HELPERS
# ============================================================================
def section(title):
    print(f"\n{'=' * 70}\n{title}\n{'=' * 70}")


def count_rows(step, df, person_level):
    """Studies and people per phase (person-level: studies = sum of n_studies)."""
    rows = []
    for phase in PHASES + ["all"]:
        d = df if phase == "all" else df[df["phase"] == phase]
        if person_level:
            rows.append({"step": step, "phase": phase, "n_studies": int(d["n_studies"].sum()),
                         "n_people": d["person_id"].nunique()})
        else:
            rows.append({"step": step, "phase": phase, "n_studies": len(d),
                         "n_people": d["person_id"].nunique()})
    return rows


# ============================================================================
# STEP 1: COLLAPSE SERIES -> STUDY, +/-6 SD, NATIVE SPLEEN > 100 HU
# ============================================================================
def load_and_qc():
    section("Step 1: load clipping-QC'd rows, collapse to study, +/-6 SD, native spleen check")
    idp = pd.read_pickle(CLEAN_CACHE_FILE)
    idp = idp[idp["phase"].isin(PHASES)].copy()
    idp["person_id"] = idp["person_id"].astype(str)
    print(f"Clipping-QC'd rows (native + PV): {len(idp):,}")

    studies = (
        idp.groupby(["person_id", "accession", "anatomical_structure", "phase"])["attenuation_Mean"]
        .median().rename("value").reset_index()
        .rename(columns={"anatomical_structure": "organ"})
    )
    print(f"Study values (median attenuation_Mean across series): {len(studies):,}\n")

    keep = pd.Series(True, index=studies.index)
    print(f"+/-{SD_CUTOFF} SD per organ/phase (all studies, ddof=1):")
    for (organ, phase), grp in studies.groupby(["organ", "phase"]):
        m, s = grp["value"].mean(), grp["value"].std()
        lo, hi = m - SD_CUTOFF * s, m + SD_CUTOFF * s
        out = (grp["value"] < lo) | (grp["value"] > hi)
        keep.loc[grp.index] = ~out
        expected = EXPECTED_SD_REMOVED.get((organ, phase))
        status = "matches expected" if out.sum() == expected else f"EXPECTED {expected}"
        print(f"  {organ:6s} {phase:14s} [{lo:.1f}, {hi:.1f}]  removed {out.sum():,} / {len(grp):,}  ({status})")

    range_out = (keep & (studies["organ"] == "spleen") & (studies["phase"] == "native")
                 & (studies["value"] > NATIVE_SPLEEN_MAX_HU))
    status = "matches expected" if range_out.sum() == EXPECTED_RANGE_REMOVED else f"EXPECTED {EXPECTED_RANGE_REMOVED}"
    print(f"Native spleen > {NATIVE_SPLEEN_MAX_HU} HU (after SD): removed {range_out.sum():,}  ({status})")

    qc = studies[keep & ~range_out].reset_index(drop=True)
    n_removed = len(studies) - len(qc)
    print(f"Total study values removed: {n_removed:,} of {len(studies):,} ({100 * n_removed / len(studies):.2f}%)")
    return studies, qc


# ============================================================================
# STEPS 2-3: STUDY LEVEL, THEN ONE VALUE PER PERSON PER PHASE
# ============================================================================
def study_level(qc, phenotype):
    """SL: accessions with both liver and spleen in that phase. liver: liver studies."""
    if phenotype == "SL":
        wide = qc.pivot_table(index=["person_id", "accession", "phase"], columns="organ", values="value").reset_index()
        return wide.dropna(subset=["liver", "spleen"])
    liver = qc[qc["organ"] == "liver"].rename(columns={"value": "liver"})
    return liver[["person_id", "accession", "phase", "liver"]]


def person_level(study, phenotype):
    """Median across a person's studies per phase, per organ; S-L from the person-level medians."""
    g = study.groupby(["person_id", "phase"])
    person = g["liver"].median().rename("liver_median").to_frame()
    person["n_studies"] = g.size()
    if phenotype == "SL":
        person["spleen_median"] = g["spleen"].median()
        person["phenotype_raw"] = person["spleen_median"] - person["liver_median"]
    else:
        person["phenotype_raw"] = -person["liver_median"]   # sign flip: higher = more fat
    return person.reset_index()


# ============================================================================
# STEP 5: ONE PHASE PER PERSON
# ============================================================================
def choose_phase(person):
    phases_by_person = person.groupby("person_id")["phase"].agg(set)
    category = phases_by_person.map(
        lambda s: "both" if len(s) == 2 else ("native_only" if "native" in s else "pv_only"))
    has_priority = set(person.loc[person["phase"] == PHASE_PRIORITY, "person_id"])
    chosen = person[(person["phase"] == PHASE_PRIORITY) | ~person["person_id"].isin(has_priority)].copy()
    chosen["phase_pv"] = (chosen["phase"] == "portal_venous").astype(int)
    chosen["phase_category"] = chosen["person_id"].map(category)
    return chosen


# ============================================================================
# COVARIATES AND EXCLUSIONS (loaded once)
# ============================================================================
def load_covariates():
    section("Covariates")
    available = pd.read_csv(COVARIATE_FILE, sep="\t", nrows=0).columns.tolist()
    missing = [c for c in COVARIATE_COLS if c not in available]
    if missing:
        sys.exit(f"Covariate column(s) {missing} not found in {COVARIATE_FILE}.\nColumns: {available}")
    covar = pd.read_csv(COVARIATE_FILE, sep="\t", usecols=COVARIATE_COLS,
                        dtype={"person_id": str, "batch": str})
    print("sequenced_gender values:", covar["sequenced_gender"].value_counts(dropna=False).to_dict())
    covar["Sex"] = (covar["sequenced_gender"].astype(str).str.strip().str.lower()
                    .map({"male": 1, "m": 1, "female": 2, "f": 2}))
    print("Sex coding (Male = 1, Female = 2):", covar["Sex"].value_counts(dropna=False).to_dict())
    covar["batch"] = covar["batch"].str.strip()
    print("batch values (categorical, written as-is):", covar["batch"].value_counts(dropna=False).to_dict())
    covar = covar.rename(columns={"sample_age": "Sample_age"})
    return covar[["person_id", "Sample_age", "Sex", "batch"]]


def load_ancestry_pcs(ancestry, geno_ids):
    """PC1-5 from one ancestry's eigenvector file, read by column position, with alignment checks."""
    path = EIGENVECTOR_FILES[ancestry]
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        first_row = f.readline().rstrip("\n").split("\t")
    print(f"{ancestry}: {path}")
    print(f"  header: {len(header)} names (last 3: {header[-3:]}) | first data row: {len(first_row)} fields")
    expected = ["person_id"] + [f"imputed_PC_{ancestry}{i}" for i in range(1, len(PCS) + 1)]
    if header[:len(expected)] != expected:
        sys.exit(f"{ancestry}: header doesn't start with {expected}\n  got {header[:len(expected)]}")

    # read by position: column 0 = ID, columns 1-5 = PC1-5 (doesn't depend on the header lining up)
    pcs = pd.read_csv(path, sep="\t", header=None, skiprows=1, usecols=range(len(PCS) + 1), dtype=str)
    pcs.columns = ["person_id"] + PCS
    pcs[PCS] = pcs[PCS].apply(pd.to_numeric, errors="coerce")

    problems = []
    n_bad_pc = pcs[PCS].isna().any(axis=1).sum()
    if n_bad_pc:
        problems.append(f"{n_bad_pc:,} rows have a non-numeric value in PC1-{len(PCS)}")
    if pcs[PCS].abs().max().max() >= 1:
        problems.append("PC values >= 1 in absolute value (eigenvector entries should be small)")
    frac_in_fam = pcs["person_id"].isin(geno_ids).mean()
    if frac_in_fam < 0.9:
        problems.append(f"only {frac_in_fam:.1%} of IDs are in the .fam file")
    n_dup = pcs["person_id"].duplicated().sum()
    if n_dup:
        problems.append(f"{n_dup:,} duplicated IDs")
    if problems:
        sys.exit(f"{ancestry}: ID and PC columns don't line up as expected:\n  - " + "\n  - ".join(problems)
                 + f"\n  first row: {first_row[:len(PCS) + 1]}")
    print(f"  {len(pcs):,} people; {frac_in_fam:.1%} of IDs in .fam; PC1-{len(PCS)} numeric -- OK")
    return pcs


def load_exclusion_ids():
    section("Liver-disease exclusion codes")
    path = Path(CONDITION_FILE)
    if not path.exists():
        candidates = sorted(p.name for p in path.parent.glob("*condition*"))
        sys.exit(f"Condition file not found: {path}\nCandidates in {path.parent}: {candidates}")
    ids = set()
    for chunk in pd.read_csv(path, sep="\t", usecols=["person_id", "condition_source_value"],
                             dtype=str, chunksize=2_000_000):
        ids.update(chunk.loc[chunk["condition_source_value"].isin(EXCLUSION_CODES), "person_id"])
    print(f"People with any exclusion code (whole release): {len(ids):,}")
    return ids


# ============================================================================
# BUILD FILES FOR ONE PHENOTYPE
# ============================================================================
def sanity_check(studies_clip_only, study, geno_ids, phenotype):
    target = SANITY_TARGETS[phenotype]
    n_final = study.loc[study["person_id"].isin(geno_ids), "person_id"].nunique()
    clip_study = study_level(studies_clip_only, phenotype)
    n_clip = clip_study.loc[clip_study["person_id"].isin(geno_ids), "person_id"].nunique()
    print(f"Sanity: genotyped people with >=1 study, native or PV (target ~{target:,}): "
          f"final QC {n_final:,} | clipping QC only {n_clip:,}")
    if all(abs(n - target) / target > SANITY_TOLERANCE for n in (n_final, n_clip)):
        print(f"  WARNING: both differ from the target by more than {SANITY_TOLERANCE:.0%}")


def plot_histograms(final, tag, transform, plot_dir):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    unit = "HU" if transform == "untransformed" else "inverse-normal units"
    for ax, (label, d) in zip(axes, [("native", final[final["phase_pv"] == 0]),
                                     ("portal_venous", final[final["phase_pv"] == 1]),
                                     ("pooled", final)]):
        if len(d):
            ax.hist(d[PHENO_COL], bins=80, color="steelblue")
        ax.set_title(f"{label} (n={len(d):,})")
        ax.set_xlabel(f"{PHENO_COL} ({unit})")
        ax.set_ylabel("people")
    fig.suptitle(f"{tag}: one value per person, final file")
    fig.tight_layout()
    fig.savefig(plot_dir / f"hist_{tag}.png", dpi=150)
    plt.close(fig)


def write_file(final, phenotype, transform, ancestry, plot_dir, n_missing_batch):
    """Transform within phase (this ancestry only), write the CSV, log the summary."""
    tag = f"{phenotype}_{transform}_{ancestry}"
    out = final.copy()
    if transform == "transformed":
        out["phenotype"] = out.groupby("phase")["phenotype_raw"].transform(rank_norm)
    else:
        out["phenotype"] = out["phenotype_raw"]
    out_cols = ["person_id", PHENO_COL, RAW_COL, "phase_pv", "phase_category", "n_studies",
                "Sample_age", "Sex", "batch"] + PCS
    out = out.rename(columns={"phenotype": PHENO_COL,
                              "phenotype_raw": RAW_COL})[out_cols]
    out_path = OUTPUT_DIR / f"hepatic_fat_pheno_{tag}.csv"
    out.to_csv(out_path, sep=",", index=False, na_rep="NA")

    print(f"\n--- {out_path.name} ---")
    print(f"Dropped for missing batch: {n_missing_batch:,}")
    print(f"Final n: {len(out):,}  |  native (phase_pv=0): {(out['phase_pv'] == 0).sum():,}  |  "
          f"PV (phase_pv=1): {(out['phase_pv'] == 1).sum():,}")
    print("phase_category:", out["phase_category"].value_counts().to_dict())
    print("batch:", out["batch"].value_counts().sort_index().to_dict())
    print(f"Females (Sex=2): {(out['Sex'] == 2).sum():,}  |  Males (Sex=1): {(out['Sex'] == 1).sum():,}")
    if transform == "transformed":
        stats = out.groupby("phase_pv")[PHENO_COL].agg(["mean", "std"]).round(3)
        print("Transformed mean/SD per phase (expect ~0 / ~1):\n" + stats.to_string())
    plot_histograms(out, tag, transform, plot_dir)
    return out


def build(phenotype, transforms, studies_clip_only, qc, geno_ids, covar, exclusion_ids, ancestry_pcs, plot_dir):
    section(f"Phenotype {phenotype}")
    counts = []

    study = study_level(qc, phenotype)
    counts += count_rows("2. studies after final QC" + (" (liver + spleen)" if phenotype == "SL" else ""),
                         study, person_level=False)
    sanity_check(studies_clip_only, study, geno_ids, phenotype)

    person = person_level(study, phenotype)
    counts += count_rows("3. person x phase medians", person, person_level=True)

    person = person[person["person_id"].isin(geno_ids)]
    if len(person) == 0:
        sys.exit("No person_id matched the .fam PMBB IDs -- the ID formats differ.\n"
                 f"  First 5 IDP person_id: {study['person_id'].drop_duplicates().head(5).tolist()}\n"
                 f"  First 5 .fam PMBB_ID:  {sorted(geno_ids)[:5]}")
    counts += count_rows("4. genotyped", person, person_level=True)

    chosen = choose_phase(person)
    counts += count_rows(f"5. one phase per person ({PHASE_PRIORITY} first)", chosen, person_level=True)

    merged = chosen.merge(covar, on="person_id", how="inner")
    counts += count_rows("7a. covariates merged", merged, person_level=True)
    n_before = len(merged)
    merged = merged[~merged["person_id"].isin(exclusion_ids)]
    print(f"Excluded for liver-disease codes: {n_before - len(merged):,}")
    counts += count_rows("7b. after exclusions", merged, person_level=True)

    # Ancestry = membership in that ancestry's eigenvector file
    in_anc = {anc: merged["person_id"].isin(set(pcs["person_id"])) for anc, pcs in ancestry_pcs.items()}
    n_both = (in_anc["EUR"] & in_anc["AFR"]).sum()
    n_neither = (~in_anc["EUR"] & ~in_anc["AFR"]).sum()
    print(f"\nAncestry of people with a phenotype (after exclusions): EUR {in_anc['EUR'].sum():,} | "
          f"AFR {in_anc['AFR'].sum():,} | neither {n_neither:,} | both {n_both:,}")
    if n_both:
        sys.exit(f"{n_both:,} people are in both the EUR and AFR eigenvector files -- expected 0")

    print("\nCounts at each step (person-level steps: n_studies = studies behind each person's value):")
    print(pd.DataFrame(counts).to_string(index=False))

    for ancestry, pcs in ancestry_pcs.items():
        sub = merged[in_anc[ancestry]].merge(pcs, on="person_id", how="inner")
        anc_counts = count_rows(f"8. {ancestry} + PCs merged", sub, person_level=True)
        n_missing_batch = sub["batch"].isna().sum()
        # drop before the transform, so untransformed and transformed keep the same people
        final = sub.dropna(subset=["phenotype_raw", "phase_pv", "phase_category", "n_studies",
                                   "Sample_age", "Sex", "batch"] + PCS)
        anc_counts += count_rows(f"9. {ancestry} after missing-value drop (final)", final, person_level=True)
        print(f"\n{ancestry}:\n" + pd.DataFrame(anc_counts).to_string(index=False))

        written = {t: write_file(final, phenotype, t, ancestry, plot_dir, n_missing_batch) for t in transforms}
        if len(written) == 2:
            u, t = written["untransformed"], written["transformed"]
            same = len(u) == len(t) and u["person_id"].tolist() == t["person_id"].tolist()
            print(f"\nUntransformed vs transformed ({phenotype}, {ancestry}): n {len(u):,} vs {len(t):,}, "
                  f"IDs {'identical' if same else 'DIFFER'}")
            if not same:
                sys.exit("Untransformed and transformed files differ in people -- this should not happen")


# ============================================================================
# MAIN
# ============================================================================
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--phenotype", choices=["SL", "liver"], help="only this phenotype (default: both)")
    parser.add_argument("--transform", choices=["untransformed", "transformed"], help="only this transform (default: both)")
    args = parser.parse_args()
    phenotypes = [args.phenotype] if args.phenotype else ["SL", "liver"]
    transforms = [args.transform] if args.transform else ["untransformed", "transformed"]

    check_rank_norm()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_dir = OUTPUT_DIR / "hepatic_fat_pheno_plots"
    plot_dir.mkdir(exist_ok=True)

    studies_clip_only, qc = load_and_qc()
    geno_ids = set(pd.read_csv(FAM_FILE, sep=r"\s+", header=None, usecols=[1], dtype=str)[1])
    print(f"\nGenotyped IDs in .fam: {len(geno_ids):,}")
    covar = load_covariates()
    section(f"Ancestry-specific PCs (PC1-{len(PCS)})")
    ancestry_pcs = {anc: load_ancestry_pcs(anc, geno_ids) for anc in ANCESTRIES}
    exclusion_ids = load_exclusion_ids()

    for phenotype in phenotypes:
        build(phenotype, transforms, studies_clip_only, qc, geno_ids, covar, exclusion_ids, ancestry_pcs, plot_dir)


if __name__ == "__main__":
    main()

