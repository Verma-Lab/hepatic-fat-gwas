import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import skew, kurtosis, probplot
from pathlib import Path

# Paths
DATA_DIR = Path("/project/vermalab_radar/data/processed/volume_attenuation")
FIG_DIR = Path("/home/agaro/verma_shared/projects/Liver_IDPs/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

def summarize_distributions(x, label): 
    print(f"{label} skewness: {skew(x): .3f}")
    print(f"{label} kurtosis: {kurtosis(x):.3f}")

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

#SERIES LEVEL EDA (LIVER)
liver_series = (
    pd.read_csv(DATA_DIR / "liver_cleaned.csv", usecols=["mean_attenuation"])
    .dropna()
    ["mean_attenuation"]
)

summarize_distributions(liver_series, "Liver series-level") #(skewness:4.301, kurtosis:115.628)

plot_histogram(
    liver_series,
    bins=80,
    xlim=(10, 180),
    xlabel="Mean Liver Attenuation (HU)",
    ylabel="Count",
    title="Series-level Liver Mean Attenuation",
    outpath=FIG_DIR / "liver_mean_attenuation_series_level_histogram.jpg"
)

#SERIES LEVEL EDA (SPLEEN)
spleen_series = (
    pd.read_csv(DATA_DIR / "spleen_cleaned.csv", usecols=["mean_attenuation"])
    .dropna()
    ["mean_attenuation"]
)

summarize_distributions(spleen_series, "Spleen series-level") #(skewness:2.114, kurtosis:30.729)

plot_histogram(
    spleen_series,
    bins=80,
    xlim=(0, 200),
    xlabel="Mean Spleen Attenuation (HU)",
    ylabel="Count",
    title="Series-level Spleen Mean Attenuation",
    outpath=FIG_DIR / "spleen_mean_attenuation_series_level_histogram.jpg"
)

#SCAN-LEVEL EDA (median across scans) 
#LIVER
liver_scan = (
    pd.read_csv(
        DATA_DIR / "liver_cleaned.csv",
        usecols=["PMBB_FINAL", "AccessionNumber_StudyUID", "mean_attenuation"]
    )
    .rename(columns={"mean_attenuation": "liver_hu"})
    .groupby(["PMBB_FINAL", "AccessionNumber_StudyUID"])["liver_hu"]
    .median()
    .reset_index()
)

summarize_distributions(liver_scan["liver_hu"], "Liver scan-level") #(skewness:3.930, kurtosis:104.652)

plot_histogram(
    liver_scan["liver_hu"],
    bins=80,
    xlim=(0, 200),
    xlabel="Scan-level Liver Mean Attenuation (HU)",
    ylabel="Number of CT Scans",
    title="Distribution of Scan-level Liver Attenuation (Median Across Series)",
    outpath=FIG_DIR / "liver_scan_level_median_across_series_histogram.jpg"
)

#SPLEEN 
spleen_scan = (
    pd.read_csv(
        DATA_DIR / "spleen_cleaned.csv",
        usecols=["PMBB_FINAL", "AccessionNumber_StudyUID", "mean_attenuation"]
    )
    .rename(columns={"mean_attenuation": "spleen_hu"})
    .groupby(["PMBB_FINAL", "AccessionNumber_StudyUID"])["spleen_hu"]
    .median()
    .reset_index()
)

summarize_distributions(spleen_scan["spleen_hu"], "Spleen scan-level") #(skewness:2.000, kurtosis:28.385)

plot_histogram(
    spleen_scan["spleen_hu"],
    bins=80,
    xlim=(0, 200),
    xlabel="Scan-level Spleen Mean Attenuation (HU)",
    ylabel="Number of CT Scans",
    title="Distribution of Scan-level Spleen Attenuation (Median Across Series)",
    outpath=FIG_DIR / "spleen_scan_level_median_across_series_histogram.jpg"
)

#HEPATIC FAT (SCAN-LEVEL)
scan_df = (
    liver_scan
    .merge(
        spleen_scan,
        on=["PMBB_FINAL", "AccessionNumber_StudyUID"],
        how="inner"
    )
    .dropna()
)

scan_df["hepatic_fat"] = scan_df["spleen_hu"] - scan_df["liver_hu"]

summarize_distributions(scan_df["hepatic_fat"], "Hepatic fat scan-level")



#HEPATIC FAT (PATIENT LEVEL)
patient_df = (
    scan_df
    .groupby("PMBB_FINAL")["hepatic_fat"]
    .median()
    .reset_index()
)

summarize_distributions(patient_df["hepatic_fat"], "Hepatic fat patient-level")

plot_histogram(
    patient_df["hepatic_fat"],
    bins=80,
    xlim=(-80, 140),
    xlabel="Hepatic Fat (ΔHU)",
    ylabel="Number of Patients",
    title="Patient-level Distribution of Median Hepatic Fat in PMBB",
    outpath=FIG_DIR / "hepatic_fat_patient_level_median_histogram.jpg"
)


# # Load data for liver 
# df = pd.read_csv("/project/vermalab_radar/data/processed/volume_attenuation/liver_cleaned.csv")

# # Keep only rows with valid mean attenuation
# df_liver = df[['mean_attenuation']].dropna()

# print(df_liver.shape)

# plt.figure()
# plt.hist(df_liver['mean_attenuation'], bins=80)
# plt.xlim(10, 180)
# plt.xlabel("Mean Liver Attenuation (HU)")
# plt.ylabel("Count")
# plt.title("Distribution of Liver Mean Attenuation")
# plt.show()


# plt.savefig("/home/agaro/verma_shared/projects/Liver_IDPs/figures/liver_mean_attenuation_histogram.jpg", dpi=300)
# plt.close()


# att_liver = df_liver['mean_attenuation']

# print("Skewness:", skew(att_liver)) #skewness is 4.3 which means the data is positively skewed 
# print("Kurtosis:", kurtosis(att_liver)) #kurtosis is 115, extremely leptokurtic distribution, meaning the data has significantly heavy tails and a very high probability of extreme outliers compared to a normal distribution


# #Load the data for the spleen 
# df_spleen = pd.read_csv("/project/vermalab_radar/data/processed/volume_attenuation/spleen_cleaned.csv")

# #Keep only the rows with valid mean attenuation
# spleen_att = df_spleen[["mean_attenuation"]].dropna()

# print(spleen_att.shape)

# plt.figure()
# plt.hist(spleen_att['mean_attenuation'], bins=80)
# plt.xlim(0, 200)
# plt.xlabel("Mean Spleen Attenuation (HU)")
# plt.ylabel("Count")
# plt.title("Distribution of Spleen Mean Attenuation")
# plt.show()


# plt.savefig("/home/agaro/verma_shared/projects/Liver_IDPs/figures/spleen_mean_attenuation_histogram.jpg", dpi=300)
# plt.close()

# print("Skewness:", skew(spleen_att)) #2.11 skewness
# print("Kurtosis:", kurtosis(spleen_att)) #30.728 kurtosis 

# #SCAN LEVEL EDA

# # Liver
# liver = pd.read_csv(
#     "/project/vermalab_radar/data/processed/volume_attenuation/liver_cleaned.csv",
#     usecols=["PMBB_FINAL", "AccessionNumber_StudyUID", "mean_attenuation"]
# ).rename(columns={"mean_attenuation": "liver_hu"})

# #Find all rows that belong to the same patient AND the same CT scan, and treat them as a group, and then find the median value of the mean attenuation
# liver_scan = (
#     liver
#     .groupby(["PMBB_FINAL", "AccessionNumber_StudyUID"])["liver_hu"]
#     .median()
#     .reset_index()
# )

# print("Liver scan-level skewness:", skew(liver_scan["liver_hu"])) #3.9 
# print("Liver scan-level kurtosis:", kurtosis(liver_scan["liver_hu"])) #104.65

# #Sanity check
# liver_scan.duplicated(
#     subset=["PMBB_FINAL", "AccessionNumber_StudyUID"]
# ).sum() #0

# #Find all rows that belong to the same patient AND the same CT scan, and treat them as a group, and then find the median value of the mean attenuation
# # Spleen
# spleen = pd.read_csv(
#     "/project/vermalab_radar/data/processed/volume_attenuation/spleen_cleaned.csv",
#     usecols=["PMBB_FINAL", "AccessionNumber_StudyUID", "mean_attenuation"]
# ).rename(columns={"mean_attenuation": "spleen_hu"})

# spleen_scan = (
#     spleen
#     .groupby(["PMBB_FINAL", "AccessionNumber_StudyUID"])["spleen_hu"]
#     .median()
#     .reset_index()
# )

# print("Spleen scan-level skewness:", skew(spleen_scan["spleen_hu"])) #1.99
# print("Spleen scan-level kurtosis:", kurtosis(spleen_scan["spleen_hu"])) #28.38

# #Sanity check 
# spleen_scan.duplicated(
#     subset=["PMBB_FINAL", "AccessionNumber_StudyUID"]
# ).sum() #0 

# plt.figure()
# plt.hist(spleen_scan['spleen_hu'], bins=80)
# plt.xlim(10, 200)
# plt.xlabel("Mean Spleen Attenuation (HU)")
# plt.ylabel("Number of CT Scans")
# plt.title("Scan-level Spleen Mean Attenuation (Median Across Series)")
# plt.show()


# plt.savefig("/home/agaro/verma_shared/projects/Liver_IDPs/figures/spleen_mean_attenuation_scan_level_median_across_series_histogram.jpg", dpi=300)
# plt.close()



# # Now we want to evaluate the hepatic fat distribution using the mean difference of both the spleen and the liver

# # In liver_scan and spleen_scan we collapsed multiple attenuation measurements from different reconstructed series into a single, scan-level HU value per patient per CT scan by taking the median across series.
# scan_df = liver_scan.merge(
#     spleen_scan,
#     on=["PMBB_FINAL", "AccessionNumber_StudyUID"],
#     how="inner"
# )


# print(scan_df.shape)

# scan_df = scan_df.dropna(subset=["liver_hu", "spleen_hu"])

# scan_df.duplicated(
#     subset=["PMBB_FINAL", "AccessionNumber_StudyUID"]
# ).sum()

# scan_df["hepatic_fat"] = scan_df["spleen_hu"] - scan_df["liver_hu"]

# print("Skewness:", skew(scan_df["hepatic_fat"])) #3.50
# print("Kurtosis:", kurtosis(scan_df["hepatic_fat"])) #68.59

# # plt.figure()
# # plt.hist(scan_df['hepatic_fat'], bins=180)
# # plt.xlim(-80, 140)
# # plt.xlabel("Hepatic Fat ( HU)")
# # plt.ylabel("Number of CT Scans")
# # plt.title("Scan-level Distribution of Median Hepatic Fat in PMBB")
# # plt.show()


# # plt.savefig("/home/agaro/verma_shared/projects/Liver_IDPs/figures/hepatic_fat_scan_level_from_series_median_histogram", dpi=300)
# # plt.close()

# # Patient level aggregation 
# patient_df = (
#     scan_df
#     .groupby("PMBB_FINAL")["hepatic_fat"]
#     .median()
#     .reset_index()
# )

# print("Patient-level skewness:", skew(patient_df["hepatic_fat"])) #3.166
# print("Patient-level kurtosis:", kurtosis(patient_df["hepatic_fat"])) #40.522


# plt.figure(figsize=(8, 6))
# plt.hist(
#     patient_df["hepatic_fat"],
#     bins=80,
#     edgecolor="black",   
#     linewidth=0.5,
#     alpha=0.85          
# )

# plt.xlim(-80, 140)
# plt.xlabel("Hepatic Fat (△HU)")
# plt.ylabel("Number of Patients")
# plt.title("Patient-level Distribution of Median Hepatic Fat in PMBB")

# plt.tight_layout()
# plt.show()

# plt.savefig("/home/agaro/verma_shared/projects/Liver_IDPs/figures/hepatic_fat_patient_level_from_series_median_histogram", dpi=300)