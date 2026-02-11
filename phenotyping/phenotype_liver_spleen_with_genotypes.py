import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import skew, kurtosis, probplot
from pathlib import Path

# Paths
COVARIATE_INPUT = Path("/home/agaro/verma_shared/projects/Liver_IDPs/covariates/hepatic_fat_all_covariates.csv")
ANCESTRY = Path("/static/PMBB/PMBB-Release-2024-3.0/Phenotypes/3.0/PMBB-Release-2024-3.0_covariates.txt")
FIG_DIR = Path("/home/agaro/verma_shared/projects/Liver_IDPs/figures")
MENOPAUSE = Path("/static/PMBB/PMBB-Release-2024-3.0/Phenotypes/3.0/PMBB-Release-2024-3.1_phenotype_menopause.txt")
# Load data 
df = pd.read_csv(COVARIATE_INPUT, sep = ",")
ancestry = pd.read_csv(ANCESTRY, sep = "\t")
menopause = pd.read_csv(MENOPAUSE, sep = "\t")

# Rename person_id to PMBB_ID in ancestry file and menopause file to match
ancestry.rename(columns={"person_id": "PMBB_ID"}, inplace=True)
menopause.rename(columns={"person_id": "PMBB_ID"}, inplace=True)


# Merge to get Class (ancestry) column - merge on PMBB_ID!
df = df.merge(ancestry[['PMBB_ID', 'Class']], on='PMBB_ID', how='left')

df = df.merge(menopause[['PMBB_ID', 'GYN_HX_MENOPAUS_AGE']], on='PMBB_ID', how='left')


print(f"\nMerge check:")
print(f"  Total rows: {len(df)}")
print(f"  Rows with ancestry data: {df['Class'].notna().sum()}")

# When creating parameters ask myself what might I want to change the next time I use it
# include a parameter without an equals sign (default) when its always different and has no default, in this case df 
# Does this value change between function calls?
# ├─ YES → Make it a parameter
# │   └─ Does it usually stay the same?
# │       ├─ YES → Give it a default value
# │       └─ NO → Make it required (no default)

# ============================================================================
# DESCRIPTIVE STATISTICS
# ============================================================================
print("\n" + "="*70)
print("COHORT CHARACTERISTICS")
print("="*70)

# Overall statistics
print(f"\nTotal participants: {len(df)}")

# Sex breakdown
n_female = (df['Sex'] == 2).sum()
n_male = (df['Sex'] == 1).sum()
print(f"\nSex distribution:")
print(f"  Females (Sex=2): {n_female} ({n_female/len(df)*100:.1f}%)")
print(f"  Males (Sex=1): {n_male} ({n_male/len(df)*100:.1f}%)")

# Age statistics
print(f"\nAge statistics:")
print(f"  Overall mean age: {df['Sample_age'].mean():.1f} years")
print(f"  Overall median age: {df['Sample_age'].median():.1f} years")
print(f"  Age range: {df['Sample_age'].min():.0f} - {df['Sample_age'].max():.0f} years")

# Split by sex for detailed stats
female_df = df[df['Sex'] == 2].copy()
male_df = df[df['Sex'] == 1].copy()

print(f"\nFemale age statistics:")
print(f"  Mean: {female_df['Sample_age'].mean():.1f} years")
print(f"  Median: {female_df['Sample_age'].median():.1f} years")
print(female_df['Sample_age'].describe())

print(f"\nMale age statistics:")
print(f"  Mean: {male_df['Sample_age'].mean():.1f} years")
print(f"  Median: {male_df['Sample_age'].median():.1f} years")
print(male_df['Sample_age'].describe())

#Steatosis prevalence
n_steatosis_overall = (df['hepatic_fat'] >= 10).sum()
print(f"  % with mild steatosis (≥10 ΔHU): {(df['hepatic_fat'] >= 10).sum() / len(df) * 100:.2f}%")

n_steatosis_male = (male_df['hepatic_fat'] >= 10).sum()
print(f"  % with mild steatosis (≥10 ΔHU): {(male_df['hepatic_fat'] >= 10).sum() / len(male_df) * 100:.2f}%")

n_steatosis_female = (female_df['hepatic_fat'] >= 10).sum()
print(f"  % with mild steatosis (≥10 ΔHU): {(female_df['hepatic_fat'] >= 10).sum() / len(female_df) * 100:.2f}%")


# Ancestry breakdown
print(f"\n" + "="*70)
print("ANCESTRY BREAKDOWN")
print("="*70)

print(f"\nOverall ancestry:")
for ancestry_group in ['EUR', 'AFR', 'SAS', 'AMR', 'EAS']:
    n = (df['Class'] == ancestry_group).sum()
    print(f"  {ancestry_group}: {n} ({n/len(df)*100:.1f}%)")

print(f"\nFemale ancestry:")
for ancestry_group in ['EUR', 'AFR', 'SAS', 'AMR', 'EAS']:
    n = (female_df['Class'] == ancestry_group).sum()
    print(f"  {ancestry_group}: {n} ({n/len(female_df)*100:.1f}%)")

print(f"\nMale ancestry:")
for ancestry_group in ['EUR', 'AFR', 'SAS', 'AMR', 'EAS']:
    n = (male_df['Class'] == ancestry_group).sum()
    print(f"  {ancestry_group}: {n} ({n/len(male_df)*100:.1f}%)")

# Menopausal breakdown 
print(f"\n" + "="*70)
print("MENOPAUSE CHARACTERISTICS (FEMALES ONLY)")
print("="*70)

# Filter to females with menopause data
female_menopause_df = female_df[female_df['GYN_HX_MENOPAUS_AGE'].notna()].copy()

print(f"\nMenopause data availability:")
print(f"  Total females in cohort: {len(female_df)}")
print(f"  Females with menopause data: {len(female_menopause_df)} ({len(female_menopause_df)/len(female_df)*100:.1f}%)")
print(f"  Females without menopause data: {female_df['GYN_HX_MENOPAUS_AGE'].isna().sum()}")

# Menopause age statistics
print(f"\nAge at menopause statistics:")
print(f"  Mean age at menopause: {female_menopause_df['GYN_HX_MENOPAUS_AGE'].mean():.1f} ± {female_menopause_df['GYN_HX_MENOPAUS_AGE'].std():.1f} years")
print(f"  Median age at menopause: {female_menopause_df['GYN_HX_MENOPAUS_AGE'].median():.1f} years")
print(f"  Range: {female_menopause_df['GYN_HX_MENOPAUS_AGE'].min():.0f} - {female_menopause_df['GYN_HX_MENOPAUS_AGE'].max():.0f} years")

# Current age distribution for women with menopause data
print(f"\nCurrent age distribution (women with menopause data):")
print(f"  Mean current age: {female_menopause_df['Sample_age'].mean():.1f} years")
print(f"  Median current age: {female_menopause_df['Sample_age'].median():.1f} years")

# Create age groups for women with menopause data
female_menopause_df['age_group'] = pd.cut(
    female_menopause_df['Sample_age'],
    bins=[0, 40, 50, 60, 70, 80, 100],
    labels=['<40', '40-49', '50-59', '60-69', '70-79', '80+'],
    right=False
)

print(f"\nAge group distribution (women with menopause data):")
age_group_counts = female_menopause_df['age_group'].value_counts().sort_index()
for age_group, count in age_group_counts.items():
    pct = count / len(female_menopause_df) * 100
    print(f"  {age_group}: {count} ({pct:.1f}%)")

# Hepatic fat by age group for menopausal women
print(f"\nHepatic fat by age group (women with menopause data):")
hepatic_fat_by_age = female_menopause_df.groupby('age_group', observed=True)['hepatic_fat'].agg([
    ('mean', 'mean'),
    ('median', 'median'),
    ('std', 'std'),
    ('count', 'count')
]).reset_index()

print(hepatic_fat_by_age.to_string(index=False))



def plot_histogram(
    x,
    bins,
    xlim,
    xlabel,
    ylabel,
    title,
    outpath,
    edgecolor="black"
):
    plt.figure(figsize=(8, 6))
    plt.hist(
        x,
        bins=bins,
        edgecolor=edgecolor,
        linewidth=0.5,
        alpha=0.85
    )
    plt.xlim(*xlim)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()



plot_histogram(
    df["hepatic_fat"],
    bins=80,
    xlim=(-80, 140),
    xlabel="Hepatic Fat (ΔHU)",
    ylabel="Number of Patients",
    title="Patient-level Distribution of Median Hepatic Fat in Genotyped Individuals in PMBB",
    outpath=FIG_DIR / "hepatic_fat_patient_level_median_histogram_with_genotype_data.jpg"
)

#The distribution of median hepatic fat (liver-spleen attenuation difference) among genotyped PMBB participants (N=~19,000) 
# shows a right-skewed distribution centered near 0 ΔHU. The majority of individuals (peak ~4,000 patients) have hepatic fat values 
# around 0 ΔHU, indicating normal liver fat content. However, there is a substantial right tail extending beyond +25 ΔHU, 
# representing individuals with hepatic steatosis (fatty liver). 

#Now I want to plot two different histograms, one for men and one for women


plot_histogram(
    female_df["hepatic_fat"],
    bins=80,
    xlim=(-80, 140),
    xlabel="Hepatic Fat (ΔHU)",
    ylabel="Number of Patients",
    title="Female Patient-level Distribution of Median Hepatic Fat in Genotyped Individuals in PMBB",
    outpath=FIG_DIR / "hepatic_fat_female_patient_level_median_histogram_with_genotype_data.jpg"
)


plot_histogram(
    male_df["hepatic_fat"],
    bins=80,
    xlim=(-80, 140),
    xlabel="Hepatic Fat (ΔHU)",
    ylabel="Number of Patients",
    title="Male Patient-level Distribution of Median Hepatic Fat in Genotyped Individuals in PMBB",
    outpath=FIG_DIR / "hepatic_fat_male_patient_level_median_histogram_with_genotype_data.jpg"
)


#Create overlay plot between male and female 
plt.figure(figsize=(10, 6))

plt.hist(
    female_df["hepatic_fat"],
    bins=80,
    range=(-80, 140),
    edgecolor="black",
    linewidth=0.5,
    alpha=0.6,
    label='Female',
    color='#E91E63'  # Pink
)

plt.hist(
    male_df["hepatic_fat"],
    bins=80,
    range=(-80, 140),
    edgecolor="black",
    linewidth=0.5,
    alpha=0.6,
    label='Male',
    color='#2196F3'  # Blue
)

plt.xlim(-80, 125)
plt.xlabel("Hepatic Fat (ΔHU)", fontsize=12)
plt.ylabel("Number of Patients", fontsize=12)
plt.title("Patient-level Distribution of Median Hepatic Fat in Genotyped Individuals in PMBB by Sex", fontsize=13)
plt.legend(loc='upper right', fontsize=11)
plt.tight_layout()
plt.savefig(FIG_DIR / "hepatic_fat_patient_level_by_sex_overlaid.jpg", dpi=300)
plt.close()

#Plot age against hepatic fat
# Create age bins for visualization
# Define age bins (adjust ranges as needed based on your data)
# Create age bins for both sexes
# Create age bins for both sexes
# female_patient_df_copy = female_patient_df.copy()
# female_patient_df_copy['age_bin'] = pd.cut(
#     female_patient_df_copy['Sample_age'], 
#     bins=[0, 30, 40, 50, 60, 70, 80, 100],
#     labels=['<30', '30-39', '40-49', '50-59', '60-69', '70-79', '80+'],
#     right=False
# )

# male_patient_df_copy = male_patient_df.copy()
# male_patient_df_copy['age_bin'] = pd.cut(
#     male_patient_df_copy['Sample_age'], 
#     bins=[0, 30, 40, 50, 60, 70, 80, 100],
#     labels=['<30', '30-39', '40-49', '50-59', '60-69', '70-79', '80+'],
#     right=False
# )

# Calculate MEDIAN hepatic fat by age bin for each sex
# female_age_summary = (
#     female_patient_df_copy
#     .groupby('age_bin', observed=True)['hepatic_fat']
#     .agg(['median', 'std', 'count'])  # Changed from 'mean' to 'median'
#     .reset_index()
# )

# male_age_summary = (
#     male_patient_df_copy
#     .groupby('age_bin', observed=True)['hepatic_fat']
#     .agg(['median', 'std', 'count'])  # Changed from 'mean' to 'median'
#     .reset_index()
# )

# # Print the summaries to verify
# print("\n=== HEPATIC FAT BY AGE AND SEX ===")
# print("\nFEMALE:")
# print(female_age_summary)
# print("\nMALE:")
# print(male_age_summary)

# Create the overlaid bar plot
# Create the overlaid bar plot
fig, ax = plt.subplots(figsize=(12, 8))

# Scatter plot for females
# ax.scatter(
#     female_patient_df['hepatic_fat'],
#     female_patient_df['Sample_age'],
#     alpha=0.4,
#     s=20,
#     color='#E91E63',
#     label='Female',
#     edgecolors='none'
# )

# Scatter plot for males
# ax.scatter(
#     male_patient_df['hepatic_fat'],
#     male_patient_df['Sample_age'],
#     alpha=0.4,
#     s=20,
#     color='#2196F3',
#     label='Male',
#     edgecolors='none'
# )

# ax.set_xlabel('Hepatic Fat (ΔHU)', fontsize=12)
# ax.set_ylabel('Age (years)', fontsize=12)
# ax.set_title('Hepatic Fat Distribution by Age and Sex in PMBB', fontsize=13)
# ax.legend(loc='upper right', fontsize=11)
# ax.axvline(x=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
# ax.grid(True, alpha=0.3)
# ax.set_xlim(-80, 140)
# ax.set_ylim(15, 95)

# plt.tight_layout()
# plt.savefig(FIG_DIR / "hepatic_fat_vs_age_scatter.jpg", dpi=300)
# plt.close()

# print("\n=== SCATTER PLOT SAVED ===")

################
##### AGE BIN BOXPLOT 

def plot_hepatic_fat_by_age_bin_and_sex(
    df,
    age_col="Sample_age",
    fat_col="hepatic_fat",
    sex_col="Sex",
    bins=(0, 40, 50, 60, 70, 80, 100),
    labels=("<40", "40–49", "50–59", "60–69", "70–79", "80+"),
    outpath=FIG_DIR / "hepatic_fat_by_age_bin_and_sex_boxplot.jpg",
    ylims=(-40, 120),
    show_n=True
):
    """
    Plot age-binned hepatic fat distributions stratified by sex using boxplots.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe containing age, hepatic fat, and sex.
    age_col : str
        Column name for age.
    fat_col : str
        Column name for hepatic fat (ΔHU).
    sex_col : str
        Column name for sex (1=Male, 2=Female).
    bins : tuple
        Age bin edges.
    labels : tuple
        Labels corresponding to age bins.
    outpath : Path or str
        Where to save the figure. If None, the plot is shown.
    ylims : tuple
        Y-axis limits.
    show_n : bool
        Whether to annotate sample size per age bin.
    """

    plot_df = df.copy()

    # Create age bins
    plot_df["age_bin"] = pd.cut(
        plot_df[age_col],
        bins=bins,
        labels=labels,
        right=False
    )

    # Map sex labels
    plot_df["Sex_label"] = plot_df[sex_col].map({1: "Male", 2: "Female"})

    # Drop missing values
    plot_df = plot_df.dropna(subset=["age_bin", fat_col, "Sex_label"])

    # Plot
    sns.set_style("whitegrid")
    g = sns.catplot(
        data=plot_df,
        x="age_bin",
        y=fat_col,
        col="Sex_label",
        kind="box",
        sharey=True,
        height=5,
        aspect=1.2,
        showfliers=False
    )

    g.set_axis_labels("Age group (years)", "Hepatic fat (ΔHU)")
    g.set_titles("{col_name}")
    g.set(ylim=ylims)

    # Steatosis threshold
    for ax in g.axes.flat:
        ax.axhline(10, linestyle="--", color="gray", alpha=0.6, linewidth=1)

    # Optional N annotations
    if show_n:
        for ax, sex in zip(g.axes.flat, ["Female", "Male"]):
            sub = plot_df[plot_df["Sex_label"] == sex]
            counts = sub.groupby("age_bin", observed=True).size()

            for i, age_bin in enumerate(labels):
                if age_bin in counts:
                    ax.text(
                        i,
                        ylims[0] + 5,
                        f"n={counts[age_bin]}",
                        ha="center",
                        va="bottom",
                        fontsize=9,
                        color="dimgray"
                    )

    plt.tight_layout()

    if outpath:
        plt.savefig(outpath, dpi=300)
        plt.close()
    else:
        plt.show()


