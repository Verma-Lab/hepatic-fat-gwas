#!/usr/bin/env python3
"""
Diagnostic: for each female's liver/spleen scan, was it taken before or
after menopause? Relevant because the existing SAIGE covariate
pipeline medians hepatic fat across all of a patient's scans, which would
blur a pre/post-menopause signal if scans span that transition.

Two classifications, since GYN_HX_MENOPAUS_AGE only covers a small subset:
  - Age-proxy: <45 vs >=55 
  - Precise: age_at_scan >= GYN_HX_MENOPAUS_AGE, for females who have that value

Grouping key: (PMBB_ID, AccessionNumber_StudyUID) = one scan.

Final step reports the same breakdowns restricted to genotyped
individuals only

Run with:
    python3 check_menopause_scan_timing.py
"""
import pandas as pd
from pathlib import Path

DATA_DIR = Path("/project/vermalab_radar/data/processed/volume_attenuation")
MAPPING_FILE = Path("/project/vermalab_radar/PMBB_RAD_ID_to_PMBB_ID_map_20250407.csv")
MENOPAUSE_FILE = Path("/static/PMBB/PMBB-Release-2024-3.0/Phenotypes/3.0/PMBB-Release-2024-3.1_phenotype_menopause.txt")
GENOTYPE_FAM = Path("/static/PMBB/PMBB-Release-2024-3.0/Imputed/chunked_bed_files/"
                     "PMBB-Release-2024-3.0_genetic_imputed.chr8_chunk31_80698638_83631804.fam")


def map_rad_to_pmbb_id(df, mapping_df, rad_id_column='id'):
    df = df.copy()
    df["rad_id"] = df[rad_id_column].str.replace('PMBB', '', regex=False)
    return df.merge(mapping_df, left_on='rad_id', right_on='PMBB_RAD_ID', how='left')


print("=" * 70)
print("STEP 1: Load and map to PMBB_ID")
print("=" * 70)
mapping_df = pd.read_csv(MAPPING_FILE)
liver = pd.read_csv(DATA_DIR / "liver_cleaned.csv")
spleen = pd.read_csv(DATA_DIR / "spleen_cleaned.csv")
menopause = pd.read_csv(MENOPAUSE_FILE, sep="\t")

liver_mapped = map_rad_to_pmbb_id(liver, mapping_df)
spleen_mapped = map_rad_to_pmbb_id(spleen, mapping_df)
print(f"  Liver rows: {len(liver_mapped):,} | matched PMBB_ID: {liver_mapped['PMBB_ID'].notna().sum():,}")
print(f"  Spleen rows: {len(spleen_mapped):,} | matched PMBB_ID: {spleen_mapped['PMBB_ID'].notna().sum():,}")

print("\n" + "=" * 70)
print("STEP 2: Collapse to one row per scan (PMBB_ID, AccessionNumber_StudyUID)")
print("=" * 70)
liver_scan = (
    liver_mapped[liver_mapped['sex'] == 'Female']
    .groupby(["PMBB_ID", "AccessionNumber_StudyUID"])
    .agg(age=("age", "first"), liver_hu=("mean_attenuation", "median"))
    .reset_index()
)
spleen_scan = (
    spleen_mapped[spleen_mapped['sex'] == 'Female']
    .groupby(["PMBB_ID", "AccessionNumber_StudyUID"])
    .agg(spleen_hu=("mean_attenuation", "median"))
    .reset_index()
)
print(f"  Unique female liver scans: {len(liver_scan):,} ({liver_scan['PMBB_ID'].nunique():,} people)")
print(f"  Unique female spleen scans: {len(spleen_scan):,} ({spleen_scan['PMBB_ID'].nunique():,} people)")

scan_df = liver_scan.merge(spleen_scan, on=["PMBB_ID", "AccessionNumber_StudyUID"], how="inner").dropna()
print(f"  Scans with both liver+spleen: {len(scan_df):,} ({scan_df['PMBB_ID'].nunique():,} people)")

print("\n" + "=" * 70)
print("STEP 3: Age-proxy classification (<45 vs >=55)")
print("=" * 70)
scan_df['age_group'] = pd.cut(scan_df['age'], bins=[0, 45, 55, 150],
                                labels=['<45', '45-55 (excluded)', '>=55'], right=False)
print("Scan-level counts:")
print(scan_df['age_group'].value_counts())

people_with_55plus_scan = scan_df.loc[scan_df['age'] >= 55, 'PMBB_ID'].nunique()
print(f"\nPeople with >=1 scan at age >=55: {people_with_55plus_scan:,} (of {scan_df['PMBB_ID'].nunique():,} total)")

print("\n" + "=" * 70)
print("STEP 4: Precise classification using GYN_HX_MENOPAUS_AGE")
print("=" * 70)
scan_df = scan_df.merge(menopause[['person_id', 'GYN_HX_MENOPAUS_AGE']],
                          left_on='PMBB_ID', right_on='person_id', how='left')
has_meno_age = scan_df['GYN_HX_MENOPAUS_AGE'].notna()
print(f"  Scans with a known GYN_HX_MENOPAUS_AGE for that person: {has_meno_age.sum():,} "
      f"({scan_df.loc[has_meno_age, 'PMBB_ID'].nunique():,} people)")

scan_df['post_menopause_precise'] = pd.NA
scan_df.loc[has_meno_age, 'post_menopause_precise'] = (
    scan_df.loc[has_meno_age, 'age'] >= scan_df.loc[has_meno_age, 'GYN_HX_MENOPAUS_AGE']
)
print("\nScan-level pre/post breakdown (precise, subset with known menopause age):")
print(scan_df.loc[has_meno_age, 'post_menopause_precise'].value_counts())

# Positive = years after menopause at that scan; negative = years before (until) menopause.
scan_df['years_since_menopause'] = pd.NA
scan_df.loc[has_meno_age, 'years_since_menopause'] = (
    scan_df.loc[has_meno_age, 'age'] - scan_df.loc[has_meno_age, 'GYN_HX_MENOPAUS_AGE']
)
print("\nyears_since_menopause distribution (subset with known menopause age, pre+post blended):")
print(scan_df.loc[has_meno_age, 'years_since_menopause'].astype(float).describe())

# Split into two sign-consistent tables instead of one blended, confusingly-signed one:
# post-menopause scans are already positive ("years post-menopause"); pre-menopause
# scans get flipped to positive too, relabeled "years until menopause" (a countdown,
# not a negative "years since" something that hasn't happened yet).
post_mask = has_meno_age & (scan_df['post_menopause_precise'] == True)
pre_mask = has_meno_age & (scan_df['post_menopause_precise'] == False)

print(f"\nyears_post_menopause distribution ({post_mask.sum():,} post-menopause scans):")
print(scan_df.loc[post_mask, 'years_since_menopause'].astype(float).describe())

years_until_menopause = -scan_df.loc[pre_mask, 'years_since_menopause'].astype(float)
print(f"\nyears_until_menopause distribution ({pre_mask.sum():,} pre-menopause scans):")
print(years_until_menopause.describe())

print("\n" + "=" * 70)
print("STEP 5: People whose scans span the menopause transition")
print("=" * 70)
print("(relevant to whether median-across-scans blurs a pre/post signal)")
spans = (
    scan_df.loc[has_meno_age]
    .groupby('PMBB_ID')['post_menopause_precise']
    .nunique()
)
n_spanning = (spans == 2).sum()
print(f"  People with >=1 pre-menopause AND >=1 post-menopause scan: {n_spanning:,} "
      f"(of {spans.shape[0]:,} people with known menopause age and >=1 dated scan)")

print("\n" + "=" * 70)
print("STEP 6: Same breakdowns, restricted to genotyped individuals")
print("=" * 70)
fam = pd.read_csv(GENOTYPE_FAM, header=None, sep=r'\s+',
                   names=["FID", "PMBB_ID", "father", "mother", "sex", "pheno"])
genotyped_ids = set(fam['PMBB_ID'].astype(str))

geno_df = scan_df[scan_df['PMBB_ID'].astype(str).isin(genotyped_ids)]
print(f"  Scans before genotype filter: {len(scan_df):,} ({scan_df['PMBB_ID'].nunique():,} people)")
print(f"  Scans after genotype filter: {len(geno_df):,} ({geno_df['PMBB_ID'].nunique():,} people)")

print("\nAge-proxy, genotyped only:")
print(geno_df['age_group'].value_counts())

people_with_55plus_scan_geno = geno_df.loc[geno_df['age'] >= 55, 'PMBB_ID'].nunique()
print(f"\nPeople with >=1 scan at age >=55 (genotyped): {people_with_55plus_scan_geno:,} "
      f"(of {geno_df['PMBB_ID'].nunique():,} total genotyped)")

people_with_under45_scan_geno = geno_df.loc[geno_df['age'] < 45, 'PMBB_ID'].nunique()
print(f"People with >=1 scan at age <45 (genotyped): {people_with_under45_scan_geno:,} "
      f"(of {geno_df['PMBB_ID'].nunique():,} total genotyped)")

has_meno_age_geno = geno_df['GYN_HX_MENOPAUS_AGE'].notna()
print(f"\nGenotyped scans with known GYN_HX_MENOPAUS_AGE: {has_meno_age_geno.sum():,} "
      f"({geno_df.loc[has_meno_age_geno, 'PMBB_ID'].nunique():,} people)")
print("Scan-level pre/post breakdown (precise, genotyped only):")
print(geno_df.loc[has_meno_age_geno, 'post_menopause_precise'].value_counts())

people_with_post_scan_geno = geno_df.loc[
    has_meno_age_geno & (geno_df['post_menopause_precise'] == True), 'PMBB_ID'].nunique()
people_with_pre_scan_geno = geno_df.loc[
    has_meno_age_geno & (geno_df['post_menopause_precise'] == False), 'PMBB_ID'].nunique()
n_known_meno_geno = geno_df.loc[has_meno_age_geno, 'PMBB_ID'].nunique()
print(f"People with >=1 post-menopause scan (precise, genotyped): {people_with_post_scan_geno:,} "
      f"(of {n_known_meno_geno:,} with known menopause age)")
print(f"People with >=1 pre-menopause scan (precise, genotyped): {people_with_pre_scan_geno:,} "
      f"(of {n_known_meno_geno:,} with known menopause age)")

print("\nyears_since_menopause distribution (genotyped only, pre+post blended):")
print(geno_df.loc[has_meno_age_geno, 'years_since_menopause'].astype(float).describe())

post_mask_geno = has_meno_age_geno & (geno_df['post_menopause_precise'] == True)
pre_mask_geno = has_meno_age_geno & (geno_df['post_menopause_precise'] == False)

print(f"\nyears_post_menopause distribution, genotyped ({post_mask_geno.sum():,} post-menopause scans):")
print(geno_df.loc[post_mask_geno, 'years_since_menopause'].astype(float).describe())

years_until_menopause_geno = -geno_df.loc[pre_mask_geno, 'years_since_menopause'].astype(float)
print(f"\nyears_until_menopause distribution, genotyped ({pre_mask_geno.sum():,} pre-menopause scans):")
print(years_until_menopause_geno.describe())

spans_geno = (
    geno_df.loc[has_meno_age_geno]
    .groupby('PMBB_ID')['post_menopause_precise']
    .nunique()
)
n_spanning_geno = (spans_geno == 2).sum()
print(f"\n  Genotyped people with >=1 pre AND >=1 post-menopause scan: {n_spanning_geno:,} "
      f"(of {spans_geno.shape[0]:,} genotyped people with known menopause age and >=1 dated scan)")

print("\n" + "=" * 70)
print("STEP 7: Sanity checks (genotyped only)")
print("=" * 70)

print("\n--- Check 1: Scans per person ---")
scans_per_person = geno_df.groupby('PMBB_ID').size()
print(f"Mean: {scans_per_person.mean():.2f} | Median: {scans_per_person.median():.1f} | "
      f"Max: {scans_per_person.max():,}")

top5_ids = scans_per_person.sort_values(ascending=False).head(5).index
print("\nTop 5 people by scan count -- ages and study IDs (showing up to 5 scans each):")
for pid in top5_ids:
    rows = (geno_df.loc[geno_df['PMBB_ID'] == pid, ['AccessionNumber_StudyUID', 'age']]
            .sort_values('age'))
    print(f"\n  PMBB_ID={pid} ({len(rows)} scans total, showing first 5):")
    print(rows.head(5).to_string(index=False))

# No real calendar date exists in this data (confirmed 2026-09-28 that PI_Time is
# post-injection contrast timing, not a date) -- use exact-age match within a person
# as the best available proxy for "same scan date", since age is computed to high
# decimal precision from the real date and two different real dates matching to
# that precision by chance is effectively impossible.
print("\nStudy IDs sharing the same age within a person (proxy for same scan date -- "
      "no real date column exists in this data):")
dup_groups = geno_df.groupby(['PMBB_ID', 'age'])['AccessionNumber_StudyUID'].apply(lambda x: sorted(x.unique()))
dup_groups = dup_groups[dup_groups.apply(len) > 1]
print(f"  (person, age) pairs with >1 distinct study ID: {len(dup_groups):,}")
if len(dup_groups) > 0:
    print(dup_groups.head(20))

print("\n--- Check 2: Menopause age distribution ---")
# Person-level, not scan-level, so people with many scans don't get overweighted.
meno_person = (geno_df.loc[has_meno_age_geno, ['PMBB_ID', 'GYN_HX_MENOPAUS_AGE']]
               .drop_duplicates('PMBB_ID'))
print(f"N people: {len(meno_person):,}")
print(meno_person['GYN_HX_MENOPAUS_AGE'].astype(float).describe())

n_meno_under40 = (meno_person['GYN_HX_MENOPAUS_AGE'].astype(float) < 40).sum()
n_meno_over60 = (meno_person['GYN_HX_MENOPAUS_AGE'].astype(float) > 60).sum()
print(f"\nWomen with menopause age <40: {n_meno_under40:,}")
print(f"Women with menopause age >60: {n_meno_over60:,}")

print("\nScans with years_since_menopause > 40 (scan age vs. menopause age):")
years_since_menopause_num = pd.to_numeric(geno_df['years_since_menopause'], errors='coerce')
extreme_mask = has_meno_age_geno & (years_since_menopause_num > 40)
extreme = geno_df.loc[extreme_mask, ['PMBB_ID', 'age', 'GYN_HX_MENOPAUS_AGE', 'years_since_menopause']]
print(f"  {len(extreme):,} scans (showing up to 5)")
print(extreme.head(5).to_string(index=False))

print("\n--- Check 3: Scans within +/-1 year of recorded menopause age ---")
near_transition_geno = has_meno_age_geno & (years_since_menopause_num.abs() <= 1)
print(f"Scans within +/-1 year of menopause: {near_transition_geno.sum():,}")

print("\nPre/post breakdown -- WITH near-transition scans (original, from Step 6):")
print(geno_df.loc[has_meno_age_geno, 'post_menopause_precise'].value_counts())

excl_near_geno = has_meno_age_geno & ~near_transition_geno
print("\nPre/post breakdown -- EXCLUDING near-transition scans:")
print(geno_df.loc[excl_near_geno, 'post_menopause_precise'].value_counts())

spans_excl_near = (
    geno_df.loc[excl_near_geno]
    .groupby('PMBB_ID')['post_menopause_precise']
    .nunique()
)
n_spanning_excl_near = (spans_excl_near == 2).sum()
print(f"\nSpanning-patient count -- WITH near-transition scans (original, from Step 6): "
      f"{n_spanning_geno:,} (of {spans_geno.shape[0]:,})")
print(f"Spanning-patient count -- EXCLUDING near-transition scans: "
      f"{n_spanning_excl_near:,} (of {spans_excl_near.shape[0]:,})")

print("\nDone.")
