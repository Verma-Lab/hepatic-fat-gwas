import pandas as pd
import numpy as np

#Load liver, spleen, menopause, and mapping data
liver = pd.read_csv("/project/vermalab_radar/data/processed/volume_attenuation/liver_cleaned.csv")
spleen = pd.read_csv("/project/vermalab_radar/data/processed/volume_attenuation/spleen_cleaned.csv")
menopause = pd.read_csv("/static/PMBB/PMBB-Release-2024-3.0/Phenotypes/3.0/PMBB-Release-2024-3.1_phenotype_menopause.txt", sep = "\t")
id_mapping = pd.read_csv("/project/vermalab_radar/PMBB_RAD_ID_to_PMBB_ID_map_20250407.csv", sep = ",")
fam = pd.read_csv("/static/PMBB/PMBB-Release-2024-3.0/Imputed/chunked_bed_files/PMBB-Release-2024-3.0_genetic_imputed.chr18_chunk22_66391159_69245318.fam", sep = "\t")

#Format fam file 
fam.columns = ["FID", "IID", "father", "mother", "sex", "pheno"]
genotype_ids = set(fam['IID'].astype(str))


pd.set_option('display.max_columns', None)

print("="*60)
print("STEP 0: ID Mapping")
print("="*60)
print(f"Mapping file has {len(id_mapping)} entries")
print(f"Sample mapping:\n{id_mapping.head()}")

# Extract RAD_ID from the IDP 'id' column (remove 'PMBB' prefix)
liver['rad_id'] = liver['id'].str.replace('PMBB', '', regex=False)
spleen['rad_id'] = spleen['id'].str.replace('PMBB', '', regex=False)

print(f"\nSample liver RAD IDs after extraction:")
print(liver[['id', 'rad_id']].head())

# Map RAD_ID to PMBB_ID using the mapping file 
liver_mapped = liver.merge(
    id_mapping,
    left_on='rad_id',
    right_on='PMBB_RAD_ID',
    how='left'
)

spleen_mapped = spleen.merge(
    id_mapping,
    left_on='rad_id',
    right_on='PMBB_RAD_ID',
    how='left'
)

print(f"\nLiver scans before mapping: {len(liver)}") #472777
print(f"Liver scans after mapping: {len(liver_mapped)}") #472777
print(f"Liver scans with successful PMBB_ID match: {liver_mapped['PMBB_ID'].notna().sum()}") #472777

print(f"\nSpleen scans before mapping: {len(spleen)}") #472590
print(f"Spleen scans after mapping: {len(spleen_mapped)}") #472590
print(f"Spleen scans with successful PMBB_ID match: {spleen_mapped['PMBB_ID'].notna().sum()}") #472590

# 1. Check total participants with liver and spleen IDPs
print("="*60)
print("STEP 1: Overall IDP Availability")
print("="*60)
print(f"Unique participants with liver: {liver_mapped['id'].nunique()}") #59904
print(f"Unique participants with spleen: {spleen_mapped['id'].nunique()}") #59,906

# Find participants with BOTH liver and spleen measurements
liver_pmbb_ids = set(liver_mapped['PMBB_ID'].dropna())
spleen_pmbb_ids = set(spleen_mapped['PMBB_ID'].dropna())
both_ids = liver_pmbb_ids & spleen_pmbb_ids
print(f"\nParticipants with BOTH liver and spleen (mapped IDs): {len(both_ids)}") #59899

# 2. Filter for females
print("\n" + "="*60)
print("STEP 2: Female Participants")
print("="*60)

# Keep one row per participant 
liver_per_person = liver_mapped.groupby('PMBB_ID').first().reset_index()
spleen_per_person = spleen_mapped.groupby('PMBB_ID').first().reset_index()

merged = liver_per_person.merge(
    spleen_per_person[['PMBB_ID', 'mean_attenuation', 'volume_ml']], 
    on='PMBB_ID', 
    suffixes=('_liver', '_spleen'),
    how='inner'
)

# Filter for females 
female_mask = merged['sex'].isin(['Female'])
females = merged[female_mask]

print(f"\nFemale participants with both IDPs: {len(females)}") #32,039 femles with both IDP's
print(f"Percentage female: {len(females)/len(merged)*100:.1f}%") #53.5%

# 3. Age distribution in females
print("\n" + "="*60)
print("STEP 3: Age Distribution in Females")
print("="*60)
print(females['age'].describe())

# 4. Check overlap with menopause data
print("\n" + "="*60)
print("STEP 4: Overlap with Menopause Data")
print("="*60)

# Get female IDs with both IDPs
female_ids_with_idps = set(females['PMBB_ID'].dropna())
print(f"Female participants with both IDPs: {len(female_ids_with_idps)}")

# Get IDs with menopause data
menopause_ids = set(menopause['person_id'].unique())
print(f"Participants with menopause data: {len(menopause_ids)}") #2440

# Find overlap
overlap_ids = female_ids_with_idps & menopause_ids
print(f"\nFemales with BOTH liver/spleen IDPs AND menopause data: {len(overlap_ids)}") #1379
print(f"Percentage of female IDP participants with menopause data: {len(overlap_ids)/len(female_ids_with_idps)*100:.1f}%") #4.3% 

# Create a dataframe of the women who have both IDPs AND menopause data
women_with_idps_and_menopause = females[females['PMBB_ID'].isin(overlap_ids)].copy()

# How many have mean_attenuation values from liver and spleen?
print("\nLiver mean_attenuation availability:")
print(women_with_idps_and_menopause['mean_attenuation_liver'].notna().value_counts())

print("\nSpleen mean_attenuation availability:")
print(women_with_idps_and_menopause['mean_attenuation_spleen'].notna().value_counts())

women_with_idps_and_menopause.head()

# Merge menopause variables into the final dataframe
women_with_idps_and_menopause = women_with_idps_and_menopause.merge(
    menopause,
    left_on='PMBB_ID',
    right_on='person_id',
    how='left'
)

print("\nColumns after merging menopause data:")
print(women_with_idps_and_menopause.columns)

# After creating women_with_idps_and_menopause, filter for genotypes
women_with_idps_and_menopause_genotyped = women_with_idps_and_menopause[
    women_with_idps_and_menopause['PMBB_ID'].isin(genotype_ids)
]

print(f"Women with liver + spleen + menopause + genotypes: {len(women_with_idps_and_menopause_genotyped)}")

# Count women by age groups
under_51 = women_with_idps_and_menopause[women_with_idps_and_menopause['age'] < 51]
over_51 = women_with_idps_and_menopause[women_with_idps_and_menopause['age'] >= 51]

print(f"Number of women under 51: {len(under_51)}")
print(f"Number of women 51 or older: {len(over_51)}")

# How many people have IDP's and imputed data 
idp_ids = set(merged['PMBB_ID'])
idp_with_genotype = idp_ids & genotype_ids
len(idp_with_genotype)

# How many women have IDP's and imputed data 
idp_female_ids = set(females['PMBB_ID'])
female_idp_with_genotype = idp_female_ids & genotype_ids
len(female_idp_with_genotype)

# How many women with IDP's and menopause data have imputed data 
idp_meno_ids = set(women_with_idps_and_menopause['PMBB_ID'])
meno_idp_with_genotype = idp_meno_ids & genotype_ids
len(meno_idp_with_genotype)

#Filter to males
male_mask = merged['sex'].isin(['Male'])
males = merged[male_mask]

idp_male_ids = set(males['PMBB_ID'])
male_idps_with_genotype = idp_male_ids & genotype_ids
len(male_idps_with_genotype)