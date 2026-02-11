import pandas as pd
from pathlib import Path

# ============================================================================
# PATHS
# ============================================================================
DATA_DIR = Path("/project/vermalab_radar/data/processed/volume_attenuation")
MAPPING_FILE = Path("/project/vermalab_radar/PMBB_RAD_ID_to_PMBB_ID_map_20250407.csv")
GENOTYPE = Path("/static/PMBB/PMBB-Release-2024-3.0/Imputed/chunked_bed_files")
PHENOTYPE = Path("/static/PMBB/PMBB-Release-2024-3.0/Phenotypes/3.0/PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
COVARIATE_FILE = Path("/static/PMBB/PMBB-Release-2024-3.0/Phenotypes/3.0")
COVARIATE_OUTPUT = Path("/home/agaro/verma_shared/projects/Liver_IDPs/covariates")
COVARIATE_OUTPUT.mkdir(exist_ok=True, parents=True)

# ============================================================================
# LOAD DATA
# ============================================================================

mapping_df = pd.read_csv(MAPPING_FILE)
liver = pd.read_csv(DATA_DIR / "liver_cleaned.csv")
spleen = pd.read_csv(DATA_DIR / "spleen_cleaned.csv")
geno = pd.read_csv(
    GENOTYPE / "PMBB-Release-2024-3.0_genetic_imputed.chr8_chunk31_80698638_83631804.fam", 
    header=None, 
    sep=r'\s+', 
    names=["IID", "PMBB_ID", "father", "mother", "sex", "pheno"]
)
covar = pd.read_csv(COVARIATE_FILE / "PMBB-Release-2024-3.0_covariates.txt", sep="\t")
covar.rename(columns={"person_id": "PMBB_ID"}, inplace=True)

# ============================================================================
# FUNCTION: MAP RAD_ID TO PMBB_ID
# ============================================================================
def map_rad_to_pmbb_id(df, mapping_df, rad_id_column='id', verbose=True):
    """Map RAD IDs to PMBB IDs for imaging data."""
    df = df.copy()
    df["rad_id"] = df[rad_id_column].str.replace('PMBB', '', regex=False)
    df_mapped = df.merge(
        mapping_df,
        left_on='rad_id',
        right_on='PMBB_RAD_ID',
        how='left'
    )
    if verbose:
        print(f"  Scans before mapping: {len(df)}")
        print(f"  Scans after mapping: {len(df_mapped)}")
        print(f"  Successful matches: {df_mapped['PMBB_ID'].notna().sum()}")
    return df_mapped

# ============================================================================
# MAP LIVER AND SPLEEN DATA
# ============================================================================
print("\nMapping liver data...")
liver_mapped = map_rad_to_pmbb_id(liver, mapping_df)

print("\nMapping spleen data...")
spleen_mapped = map_rad_to_pmbb_id(spleen, mapping_df)

# ============================================================================
# CREATE SCAN-LEVEL DATA (MEDIAN ACROSS SCANS)
# ============================================================================
print("\nCreating scan-level data...")
liver_scan = (
    liver_mapped
    .rename(columns={"mean_attenuation": "liver_hu"})
    .groupby(["PMBB_ID", "AccessionNumber_StudyUID"])["liver_hu"]
    .median()
    .reset_index()
)

spleen_scan = (
    spleen_mapped
    .rename(columns={"mean_attenuation": "spleen_hu"})
    .groupby(["PMBB_ID", "AccessionNumber_StudyUID"])["spleen_hu"]
    .median()
    .reset_index()
)

# ============================================================================
# FILTER TO GENOTYPED INDIVIDUALS
# ============================================================================
print("\nFiltering to genotyped individuals...")
print(f"  Before filter - Liver: {liver_scan['PMBB_ID'].nunique()} unique individuals")
print(f"  Before filter - Spleen: {spleen_scan['PMBB_ID'].nunique()} unique individuals")

liver_scan = liver_scan[liver_scan['PMBB_ID'].isin(geno["PMBB_ID"])]
spleen_scan = spleen_scan[spleen_scan['PMBB_ID'].isin(geno["PMBB_ID"])]

print(f"  After filter - Liver: {liver_scan['PMBB_ID'].nunique()} unique individuals")
print(f"  After filter - Spleen: {spleen_scan['PMBB_ID'].nunique()} unique individuals")

# ============================================================================
# CALCULATE HEPATIC FAT (SCAN-LEVEL)
# ============================================================================
print("\nCalculating hepatic fat...")
scan_df = (
    liver_scan
    .merge(
        spleen_scan,
        on=["PMBB_ID", "AccessionNumber_StudyUID"],
        how="inner"
    )
    .dropna()
)

scan_df["hepatic_fat"] = scan_df["spleen_hu"] - scan_df["liver_hu"]

# ============================================================================
# PATIENT-LEVEL HEPATIC FAT (MEDIAN ACROSS SCANS)
# ============================================================================
print("\nCreating patient-level hepatic fat (median across scans)...")
patient_df = (
    scan_df
    .groupby("PMBB_ID")["hepatic_fat"]
    .median()
    .reset_index()
)
print(f"  Total patients with hepatic fat: {len(patient_df)}")

# ============================================================================
# MERGE WITH COVARIATES
# ============================================================================
print("\nMerging with covariates (age, sex, PCs)...")
sex_patient_df = (
    patient_df
    .merge(
        covar,
        on=["PMBB_ID"],
        how="inner"
    )
)
print(f"  Patients after merging covariates: {len(sex_patient_df)}")

# ============================================================================
# CREATE COVARIATE FILE FOR SAIGE
# ============================================================================
print("\nCreating SAIGE covariate file...")

# Select columns needed
covariate_df = sex_patient_df[['PMBB_ID', 'Sample_age', 'Sequenced_gender', 
                                 'hepatic_fat', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']].copy()

# Code sex: Male = 1, Female = 2
covariate_df['Sex'] = covariate_df['Sequenced_gender'].map({'Male': 1, 'Female': 2})

# Verify coding
print("\nSex coding verification:")
print(covariate_df[['Sequenced_gender', 'Sex']].value_counts())

saige_covariates = covariate_df.drop(columns=['Sequenced_gender'], axis =1)

# Remove any rows with missing data
saige_covariates_clean = saige_covariates.dropna()

print(f"\nRows before removing missing: {len(saige_covariates_clean)}")
print(f"Rows after removing missing: {len(saige_covariates_clean)}")


# ============================================================================
# EXCLUDE INDIVIDUALS WITH SPECIFIC DIAGNOSIS CODES
# ============================================================================
phenotype = pd.read_csv(PHENOTYPE, sep="\t")

n_before_exclusion = len(saige_covariates_clean) 

# Define exclusion codes
exclusion_codes = [
    'B18.0', 'B18.1', 'B18.2', '070.32', '070.21', '070.22', '070.23', 
    '070.31', '070.33', '070.54', '571.0', 'K70.0', '571.1', 'K70.1', 
    '571.2', 'K70.3', 'K70.2', '571.3', 'K70.4', 'K70.40', 'K70.41', 
    'K70.9', '303.0', '303.9', 'F10.229', 'F10.20'
]

# Find individuals with any of these codes
individuals_to_exclude = phenotype[
    phenotype['condition_source_value'].isin(exclusion_codes)
]['person_id'].unique()

print(f"  Individuals with exclusion codes: {len(individuals_to_exclude)}")

# Exclude these individuals from the covariate file
saige_covariates_clean_exclusion = saige_covariates_clean[
    ~saige_covariates_clean['PMBB_ID'].isin(individuals_to_exclude)
]

# Calculate and report the difference
n_dropped = n_before_exclusion - len(saige_covariates_clean_exclusion)

print(f"  Participants before exclusion: {n_before_exclusion}")
print(f"  Participants dropped: {n_dropped}")
print(f"  Participants after exclusion: {len(saige_covariates_clean_exclusion)}")

# ============================================================================
# SAVE COVARIATE FILE
# ============================================================================
output_file = COVARIATE_OUTPUT / "hepatic_fat_all_covariates.csv"
saige_covariates_clean_exclusion.to_csv(output_file, sep=',', index=False, na_rep='NA')

print(f"\n{'='*70}")
print(f"✓ COVARIATE FILE CREATED SUCCESSFULLY")
print(f"{'='*70}")
print(f"Location: {output_file}")
print(f"Total participants: {len(saige_covariates_clean_exclusion)}")
print(f"  Females (sex=2): {(saige_covariates_clean_exclusion['Sex'] == 2).sum()}")
print(f"  Males (sex=1): {(saige_covariates_clean_exclusion['Sex'] == 1).sum()}")

# Show first few rows
print("\nFirst 5 rows of covariate file:")
print(saige_covariates_clean_exclusion.head())


