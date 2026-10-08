"""
Checks whether the released, already-PMBB_ID-mapped IDP file
(PMBB-Release-2024-3.0_imaging_ct_idps.txt) produces a hepatic_fat phenotype
close enough to the existing production covariate file
(hepatic_fat_all_covariates.csv) that downstream GWAS/FUMA/coloc/LDSC would
NOT need to be rerun -- or whether it diverges enough that they would.

Builds patient-level hepatic_fat from the released IDP file (grouping by
person_id + accession, since accession is supposed to be one real study per
your advisor's explanation -- no RAD ID mapping needed, this file already
has person_id), then compares it directly against the PMBB_ID/hepatic_fat
columns already sitting in the production covariate file.

Run with:
    python3 compare_idp_vs_production_phenotype.py
"""
import pandas as pd
from pathlib import Path

IDP_FILE = Path("/static/PMBB/PMBB-Release-2024-3.0/Imaging/IDP/PMBB-Release-2024-3.0_imaging_ct_idps.txt")
PRODUCTION_COVARIATE_FILE = Path("/home/agaro/verma_shared/projects/Liver_IDPs/covariates/hepatic_fat_all_covariates.csv")


def build_patient_hepatic_fat_from_idp(idp):
    liver = idp[idp["anatomical_structure"] == "liver"]
    spleen = idp[idp["anatomical_structure"] == "spleen"]

    liver_scan = (
        liver.groupby(["person_id", "accession"])["attenuation_Mean"]
        .median()
        .rename("liver_hu")
        .reset_index()
    )
    spleen_scan = (
        spleen.groupby(["person_id", "accession"])["attenuation_Mean"]
        .median()
        .rename("spleen_hu")
        .reset_index()
    )

    scan_df = liver_scan.merge(spleen_scan, on=["person_id", "accession"], how="inner").dropna()
    scan_df["hepatic_fat"] = scan_df["spleen_hu"] - scan_df["liver_hu"]

    patient_df = (
        scan_df.groupby("person_id")["hepatic_fat"]
        .median()
        .reset_index()
        .rename(columns={"person_id": "PMBB_ID"})
    )
    scans_per_person = scan_df.groupby("person_id").size().rename("n_accessions_idp")

    return patient_df, scans_per_person


print("=" * 70)
print("Loading released IDP file and production covariate file")
print("=" * 70)
idp = pd.read_csv(IDP_FILE, sep="\t")
production = pd.read_csv(PRODUCTION_COVARIATE_FILE)

print(f"IDP file rows: {len(idp):,} | unique person_id: {idp['person_id'].nunique():,}")
print(f"Production covariate file patients: {len(production):,}")

print("\n" + "=" * 70)
print("Building patient-level hepatic_fat from the IDP file")
print("=" * 70)
patient_idp, scans_idp = build_patient_hepatic_fat_from_idp(idp)
print(f"Patients with hepatic_fat from IDP file: {len(patient_idp):,}")
print(f"Accessions/person -- mean: {scans_idp.mean():.2f} | median: {scans_idp.median():.1f} | "
      f"max: {scans_idp.max():,}")

print("\n" + "=" * 70)
print("Comparing IDP-derived hepatic_fat vs. existing production phenotype")
print("=" * 70)
comparison = (
    production[["PMBB_ID", "hepatic_fat"]]
    .rename(columns={"hepatic_fat": "hepatic_fat_production"})
    .merge(
        patient_idp.rename(columns={"hepatic_fat": "hepatic_fat_idp"}),
        on="PMBB_ID",
        how="inner",
    )
)
comparison["abs_diff"] = (comparison["hepatic_fat_production"] - comparison["hepatic_fat_idp"]).abs()

print(f"\nPatients compared: {len(comparison):,} (out of {len(production):,} in production file)")
print(f"Correlation (production vs. IDP-derived hepatic_fat): "
      f"{comparison['hepatic_fat_production'].corr(comparison['hepatic_fat_idp']):.4f}")

print("\nAbsolute difference in hepatic_fat (production - IDP-derived):")
print(comparison["abs_diff"].describe())

for thresh in [0.5, 1.0, 2.0]:
    n_large = (comparison["abs_diff"] > thresh).sum()
    print(f"  Patients with |difference| > {thresh}: {n_large:,} "
          f"({100 * n_large / len(comparison):.1f}%)")

print("\nDone.")
print("If correlation is very high (e.g. >0.98) and few patients cross the")
print("thresholds above, the released IDP file likely does NOT require")
print("rerunning downstream GWAS/FUMA/coloc/LDSC. If it diverges meaningfully")
print("(similar to or worse than the AccessionNumber_StudyUID vs. true_accession")
print("comparison), that's evidence the production phenotype needs to be")
print("regenerated from this file before those results can be trusted.")
