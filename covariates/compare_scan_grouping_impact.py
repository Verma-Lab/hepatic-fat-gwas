#!/usr/bin/env python3
"""
Quantifies the actual impact of the AccessionNumber_StudyUID overcounting bug
(confirmed 2026-09-29 via cross-reference against the consolidated IDP file:
~4x overcount in 2/2 people checked) on the real hepatic_fat phenotype used
by create_saige_covariate_file.py -- the file that fed the actual GWAS.

Runs the exact same liver/spleen -> scan_df -> patient_df pipeline as
create_saige_covariate_file.py, twice: once grouped by the original
AccessionNumber_StudyUID (buggy), once by true_accession (the numeric
segment before the first underscore -- the fix). Same population as the
real pipeline: ALL genotyped patients, both sexes, no menopause/female
restriction (this is about the core phenotype, not the menopause side-check).

Reports: how much the scan-row count drops, and -- the number that actually
matters -- how much the final patient-level hepatic_fat value changes
between the two groupings, and whether bigger overcounting for a patient
predicts a bigger change in their phenotype value.

Run with:
    python3 compare_scan_grouping_impact.py
"""
import pandas as pd
from pathlib import Path

DATA_DIR = Path("/project/vermalab_radar/data/processed/volume_attenuation")
MAPPING_FILE = Path("/project/vermalab_radar/PMBB_RAD_ID_to_PMBB_ID_map_20250407.csv")
GENOTYPE = Path("/static/PMBB/PMBB-Release-2024-3.0/Imputed/chunked_bed_files")


def map_rad_to_pmbb_id(df, mapping_df, rad_id_column='id'):
    df = df.copy()
    df["rad_id"] = df[rad_id_column].str.replace('PMBB', '', regex=False)
    return df.merge(mapping_df, left_on='rad_id', right_on='PMBB_RAD_ID', how='left')


def build_patient_hepatic_fat(liver_mapped, spleen_mapped, group_col, genotyped_ids):
    """Reproduces create_saige_covariate_file.py's pipeline exactly, parameterized
    by which column defines 'one scan' -- AccessionNumber_StudyUID (original) or
    true_accession (fixed)."""
    liver_scan = (
        liver_mapped.rename(columns={"mean_attenuation": "liver_hu"})
        .groupby(["PMBB_ID", group_col])["liver_hu"].median().reset_index()
    )
    spleen_scan = (
        spleen_mapped.rename(columns={"mean_attenuation": "spleen_hu"})
        .groupby(["PMBB_ID", group_col])["spleen_hu"].median().reset_index()
    )
    liver_scan = liver_scan[liver_scan['PMBB_ID'].isin(genotyped_ids)]
    spleen_scan = spleen_scan[spleen_scan['PMBB_ID'].isin(genotyped_ids)]

    scan_df = liver_scan.merge(spleen_scan, on=["PMBB_ID", group_col], how="inner").dropna()
    scan_df["hepatic_fat"] = scan_df["spleen_hu"] - scan_df["liver_hu"]

    patient_df = scan_df.groupby("PMBB_ID")["hepatic_fat"].median().reset_index()
    scans_per_person = scan_df.groupby("PMBB_ID").size().rename("n_scans")

    return patient_df, scans_per_person


print("=" * 70)
print("Loading data (same as create_saige_covariate_file.py)")
print("=" * 70)
mapping_df = pd.read_csv(MAPPING_FILE)
liver = pd.read_csv(DATA_DIR / "liver_cleaned.csv")
spleen = pd.read_csv(DATA_DIR / "spleen_cleaned.csv")
geno = pd.read_csv(
    GENOTYPE / "PMBB-Release-2024-3.0_genetic_imputed.chr8_chunk31_80698638_83631804.fam",
    header=None, sep=r'\s+', names=["IID", "PMBB_ID", "father", "mother", "sex", "pheno"]
)
genotyped_ids = set(geno["PMBB_ID"])

liver_mapped = map_rad_to_pmbb_id(liver, mapping_df)
spleen_mapped = map_rad_to_pmbb_id(spleen, mapping_df)

liver_mapped['true_accession'] = liver_mapped['AccessionNumber_StudyUID'].str.split('_').str[0]
spleen_mapped['true_accession'] = spleen_mapped['AccessionNumber_StudyUID'].str.split('_').str[0]

print("\n" + "=" * 70)
print("Building patient_df TWICE: original grouping vs. corrected grouping")
print("=" * 70)
patient_orig, scans_orig = build_patient_hepatic_fat(
    liver_mapped, spleen_mapped, 'AccessionNumber_StudyUID', genotyped_ids)
patient_fixed, scans_fixed = build_patient_hepatic_fat(
    liver_mapped, spleen_mapped, 'true_accession', genotyped_ids)

print(f"\nOriginal (AccessionNumber_StudyUID) grouping:")
print(f"  Total scan-rows: {scans_orig.sum():,} | Patients: {len(patient_orig):,}")
print(f"  Scans/person -- mean: {scans_orig.mean():.2f} | median: {scans_orig.median():.1f} | "
      f"max: {scans_orig.max():,}")

print(f"\nCorrected (true_accession) grouping:")
print(f"  Total scan-rows: {scans_fixed.sum():,} | Patients: {len(patient_fixed):,}")
print(f"  Scans/person -- mean: {scans_fixed.mean():.2f} | median: {scans_fixed.median():.1f} | "
      f"max: {scans_fixed.max():,}")

print(f"\nOverall scan-row inflation factor: {scans_orig.sum() / scans_fixed.sum():.2f}x")

print("\n" + "=" * 70)
print("Comparing patient-level hepatic_fat: original vs. corrected")
print("=" * 70)
comparison = (
    patient_orig.rename(columns={'hepatic_fat': 'hepatic_fat_original'})
    .merge(patient_fixed.rename(columns={'hepatic_fat': 'hepatic_fat_fixed'}), on='PMBB_ID', how='inner')
)
comparison['abs_diff'] = (comparison['hepatic_fat_original'] - comparison['hepatic_fat_fixed']).abs()
comparison = comparison.merge(scans_orig.rename('n_scans_original'), on='PMBB_ID')
comparison = comparison.merge(scans_fixed.rename('n_scans_fixed'), on='PMBB_ID')
comparison['inflation_ratio'] = comparison['n_scans_original'] / comparison['n_scans_fixed']

print(f"\nPatients compared: {len(comparison):,}")
print(f"Correlation (original vs. fixed hepatic_fat): {comparison['hepatic_fat_original'].corr(comparison['hepatic_fat_fixed']):.4f}")
print("\nAbsolute difference in hepatic_fat (original - fixed):")
print(comparison['abs_diff'].describe())

for thresh in [0.5, 1.0, 2.0]:
    n_large = (comparison['abs_diff'] > thresh).sum()
    print(f"  Patients with |difference| > {thresh}: {n_large:,} "
          f"({100 * n_large / len(comparison):.1f}%)")

print(f"\nCorrelation (per-patient inflation ratio vs. |hepatic_fat difference|): "
      f"{comparison['inflation_ratio'].corr(comparison['abs_diff']):.4f}")
print("(a positive correlation here means patients whose scan count was more inflated "
      "also saw a bigger change in their phenotype value -- i.e., the bug isn't just "
      "noise, it's systematically distorting the most-affected patients' values)")

print("\nDone.")
