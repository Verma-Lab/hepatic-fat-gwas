#!/usr/bin/env python3
"""
For studies that have BOTH native and portal_venous phases, quantifies how
well they agree (ICC, Bland-Altman) -- informs the combined native/PV
phenotype decision, not a native-vs-PV "which one to use" choice.

Loads the QC-cleaned cache produced by clipping_qc_diagnostic.py -- run
that script first (through cell 10) to produce idp_liver_spleen_clean.pkl.
"""

# %% Cell 0: imports and settings
import pandas as pd
from pathlib import Path

CLEAN_CACHE_FILE = Path("idp_liver_spleen_clean.pkl")


# %% Cell 1: load the QC-cleaned cache
idp_clean = pd.read_pickle(CLEAN_CACHE_FILE)
print(f"Loaded {len(idp_clean):,} QC-cleaned rows")
print(idp_clean.groupby(["anatomical_structure", "phase"]).size())


# %% Cell 2: collapse series -> one value per study, structure, phase
target = idp_clean[idp_clean["phase"].isin(["native", "portal_venous"])]
collapsed = (
    target.groupby(["person_id", "accession", "anatomical_structure", "phase"])["attenuation_Mean"]
    .median()
    .reset_index()
)
pivoted = collapsed.pivot_table(
    index=["person_id", "accession"],
    columns=["anatomical_structure", "phase"],
    values="attenuation_Mean",
)
pivoted.columns = [f"{organ}_{phase}" for organ, phase in pivoted.columns]
pivoted = pivoted.reset_index()
pivoted.head()


# %% Cell 3: paired native vs portal venous correlations
required_cols = ["liver_native", "liver_portal_venous", "spleen_native", "spleen_portal_venous"]
paired = pivoted.dropna(subset=required_cols).copy()
print(f"Studies with both phases for both organs: {len(paired):,}")

paired["hepatic_fat_native"] = paired["spleen_native"] - paired["liver_native"]
paired["hepatic_fat_pv"] = paired["spleen_portal_venous"] - paired["liver_portal_venous"]

for label, a, b in [
    ("Liver", "liver_native", "liver_portal_venous"),
    ("Spleen", "spleen_native", "spleen_portal_venous"),
    ("Hepatic fat", "hepatic_fat_native", "hepatic_fat_pv"),
]:
    r = paired[a].corr(paired[b])
    print(f"{label:12s} native vs PV: r = {r:.4f}, r^2 = {r**2:.4f}")
    rho = paired[a].corr(paired[b], method="spearman")
    print(f"{label:12s} native vs PV: Spearman rho = {rho:.4f}")


# %% Cell 4: compare IQRs directly -- r^2 doesn't tell you whether the two
# phases' distributions actually overlap, just whether they move together
for label, a, b in [
    ("Liver", "liver_native", "liver_portal_venous"),
    ("Spleen", "spleen_native", "spleen_portal_venous"),
    ("Hepatic fat", "hepatic_fat_native", "hepatic_fat_pv"),
]:
    iqr_a = paired[a].quantile(0.75) - paired[a].quantile(0.25)
    iqr_b = paired[b].quantile(0.75) - paired[b].quantile(0.25)
    print(f"{label:12s} -- native IQR: {iqr_a:.2f} | PV IQR: {iqr_b:.2f} | "
          f"native median: {paired[a].median():.2f} | PV median: {paired[b].median():.2f}")


# %% Cell 5: ICC -- unlike r^2, this penalizes a systematic shift between
# native and PV, not just rank disagreement. Needs: pip install pingouin
import pingouin as pg

def compute_icc(df, study_id_col, native_col, pv_col, label):
    long_df = pd.concat([
        pd.DataFrame({"target": df[study_id_col], "rater": "native", "rating": df[native_col]}),
        pd.DataFrame({"target": df[study_id_col], "rater": "portal_venous", "rating": df[pv_col]}),
    ], ignore_index=True)
    icc_result = pg.intraclass_corr(data=long_df, targets="target", raters="rater", ratings="rating")
    print(f"\nICC -- {label}:")
    print(icc_result.to_string(index=False))
    return icc_result

paired["study_id"] = paired["person_id"].astype(str) + "_" + paired["accession"].astype(str)

for label, a, b in [
    ("Liver", "liver_native", "liver_portal_venous"),
    ("Spleen", "spleen_native", "spleen_portal_venous"),
    ("Hepatic fat", "hepatic_fat_native", "hepatic_fat_pv"),
]:
    compute_icc(paired, "study_id", a, b, label)
# Compare ICC2 (single, absolute agreement: penalizes a constant shift) with
# ICC3 (single, consistency: ignores a constant shift). ICC3 >> ICC2 means PV
# is shifted but preserves order; both low means PV reorders participants.


# %% Cell 6: Bland-Altman plots -- shows the mean difference (bias) and
# limits of agreement, and whether the difference is constant or changes
# with the magnitude of the measurement
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

BA_OUTPUT_DIR = Path("qc_clipping/native_vs_pv_icc_bland_altman")
BA_OUTPUT_DIR.mkdir(exist_ok=True)

def bland_altman_plot(native_vals, pv_vals, label, output_path):
    diff = native_vals - pv_vals
    mean_vals = (native_vals + pv_vals) / 2
    mean_diff = diff.mean()
    sd_diff = diff.std()
    loa_upper = mean_diff + 1.96 * sd_diff
    loa_lower = mean_diff - 1.96 * sd_diff

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(mean_vals, diff, alpha=0.4, s=15)
    ax.axhline(mean_diff, color="red", linestyle="-", label=f"Mean diff = {mean_diff:.2f}")
    ax.axhline(loa_upper, color="gray", linestyle="--", label=f"+1.96 SD = {loa_upper:.2f}")
    ax.axhline(loa_lower, color="gray", linestyle="--", label=f"-1.96 SD = {loa_lower:.2f}")
    ax.set_xlabel("Mean of native and portal venous")
    ax.set_ylabel("Native - Portal venous")
    ax.set_title(f"Bland-Altman: {label}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"{label:12s} -- mean diff: {mean_diff:.2f} | "
          f"limits of agreement: [{loa_lower:.2f}, {loa_upper:.2f}]")

for label, a, b, fname in [
    ("Liver", "liver_native", "liver_portal_venous", "bland_altman_liver.png"),
    ("Spleen", "spleen_native", "spleen_portal_venous", "bland_altman_spleen.png"),
    ("Hepatic fat", "hepatic_fat_native", "hepatic_fat_pv", "bland_altman_hepatic_fat.png"),
]:
    bland_altman_plot(paired[a], paired[b], label, BA_OUTPUT_DIR / fname)


# %% Cell 7: rerun agreement after final QC (study level, post-clipping):
# +/-6 SD per organ + phase on all studies, then native spleen > 100 HU check
SD_CUTOFF = 6
NATIVE_SPLEEN_MAX_HU = 100
keep = pd.Series(True, index=collapsed.index)
print(f"+/-{SD_CUTOFF} SD per organ/phase (all studies in that organ/phase):")
for (organ, phase), grp in collapsed.groupby(["anatomical_structure", "phase"]):
    m, s = grp["attenuation_Mean"].mean(), grp["attenuation_Mean"].std()
    lo, hi = m - SD_CUTOFF * s, m + SD_CUTOFF * s
    out = (grp["attenuation_Mean"] < lo) | (grp["attenuation_Mean"] > hi)
    keep.loc[grp.index] = ~out
    print(f"  {organ:6s} {phase:14s} range [{lo:.1f}, {hi:.1f}]  removed {out.sum():,} / {len(grp):,}")
range_out = (keep & (collapsed["anatomical_structure"] == "spleen") & (collapsed["phase"] == "native")
             & (collapsed["attenuation_Mean"] > NATIVE_SPLEEN_MAX_HU))
print(f"Native spleen > {NATIVE_SPLEEN_MAX_HU} HU (after SD): removed {range_out.sum():,} more")
collapsed_finalqc = collapsed[keep & ~range_out]
# rebuild pivot and paired set (liver + spleen, native + PV, same accession)
pivoted_finalqc = collapsed_finalqc.pivot_table(
    index=["person_id", "accession"], columns=["anatomical_structure", "phase"], values="attenuation_Mean")
pivoted_finalqc.columns = [f"{organ}_{phase}" for organ, phase in pivoted_finalqc.columns]
paired_finalqc = pivoted_finalqc.reset_index().dropna(subset=required_cols).copy()
paired_finalqc["hepatic_fat_native"] = paired_finalqc["spleen_native"] - paired_finalqc["liver_native"]
paired_finalqc["hepatic_fat_pv"] = paired_finalqc["spleen_portal_venous"] - paired_finalqc["liver_portal_venous"]
paired_finalqc["study_id"] = paired_finalqc["person_id"].astype(str) + "_" + paired_finalqc["accession"].astype(str)
print(f"\nPaired studies: before QC {len(paired):,} ({paired['person_id'].nunique():,} people) | "
      f"after QC {len(paired_finalqc):,} ({paired_finalqc['person_id'].nunique():,} people)")
# newer pingouin labels ICC2/ICC3 as ICC(A,1)/ICC(C,1); accept either
ICC_TYPE_ALIASES = {"ICC(A,1)": "ICC2", "ICC(C,1)": "ICC3"}
def icc_values(icc_result):
    return icc_result.set_index("Type")["ICC"].rename(index=ICC_TYPE_ALIASES)
summary_rows = []
for label, a, b, fname in [
    ("Liver", "liver_native", "liver_portal_venous", "bland_altman_liver_finalqc.png"),
    ("Spleen", "spleen_native", "spleen_portal_venous", "bland_altman_spleen_finalqc.png"),
    ("Hepatic fat", "hepatic_fat_native", "hepatic_fat_pv", "bland_altman_hepatic_fat_finalqc.png"),
]:
    icc = icc_values(compute_icc(paired_finalqc, "study_id", a, b, f"{label} (final QC)"))
    bland_altman_plot(paired_finalqc[a], paired_finalqc[b], f"{label} (final QC)", BA_OUTPUT_DIR / fname)
    diff = paired_finalqc[a] - paired_finalqc[b]
    summary_rows.append({
        "measure": label,
        "pearson_r": paired_finalqc[a].corr(paired_finalqc[b]),
        "spearman_rho": paired_finalqc[a].corr(paired_finalqc[b], method="spearman"),
        "ICC2_absolute": icc["ICC2"], "ICC3_consistency": icc["ICC3"],
        "BA_mean_diff": diff.mean(),
        "BA_loa_lower": diff.mean() - 1.96 * diff.std(), "BA_loa_upper": diff.mean() + 1.96 * diff.std(),
    })
print("\nAgreement after final QC (native - PV):")
print(pd.DataFrame(summary_rows).to_string(index=False, float_format=lambda x: f"{x:.3f}"))


# %% Cell 8: final QC (Cell 7) -> genotyped people only -> paired studies ->
# one median per person (as in the GWAS phenotype) -> agreement
FAM_FILE = "/static/PMBB/PMBB-Release-2026-4.0/Imputed/chunked_bed_files/PMBB-Release-2026-4.0_genetic_imputed.463_chr8_chunk1_149284_2198252.fam"
geno_ids = set(pd.read_csv(FAM_FILE, sep=r"\s+", header=None, usecols=[1], dtype=str)[1])
geno_mask = collapsed_finalqc["person_id"].astype(str).isin(geno_ids)
print(f"Genotyped IDs in .fam: {len(geno_ids):,} | people in final-QC data: "
      f"{collapsed_finalqc['person_id'].nunique():,} -> genotyped: "
      f"{collapsed_finalqc.loc[geno_mask, 'person_id'].nunique():,}")
collapsed_geno = collapsed_finalqc[geno_mask]
pivoted_geno = collapsed_geno.pivot_table(
    index=["person_id", "accession"], columns=["anatomical_structure", "phase"], values="attenuation_Mean")
pivoted_geno.columns = [f"{organ}_{phase}" for organ, phase in pivoted_geno.columns]
paired_geno = pivoted_geno.reset_index().dropna(subset=required_cols)
n_studies_geno = len(paired_geno)
# median across each person's paired studies, then hepatic fat from person-level values
person_geno = paired_geno.groupby("person_id")[required_cols].median().reset_index()
person_geno["person_id"] = person_geno["person_id"].astype(str)
person_geno["hepatic_fat_native"] = person_geno["spleen_native"] - person_geno["liver_native"]
person_geno["hepatic_fat_pv"] = person_geno["spleen_portal_venous"] - person_geno["liver_portal_venous"]
person_geno["neg_liver_pv"] = -person_geno["liver_portal_venous"]
print(f"Paired studies (genotyped): {n_studies_geno:,} -> people after per-person median: {len(person_geno):,}")
def agreement_row(df, id_col, a, b, label, ba_path=None, zscore=False):
    if zscore:
        df = df.assign(**{a: (df[a] - df[a].mean()) / df[a].std(), b: (df[b] - df[b].mean()) / df[b].std()})
    icc = icc_values(compute_icc(df, id_col, a, b, label))
    row = {"pearson_r": df[a].corr(df[b]), "spearman_rho": df[a].corr(df[b], method="spearman"),
           "ICC2_absolute": float("nan"), "ICC3_consistency": icc["ICC3"],
           "BA_mean_diff": float("nan"), "BA_loa_lower": float("nan"), "BA_loa_upper": float("nan")}
    if ba_path is not None:
        bland_altman_plot(df[a], df[b], label, ba_path)
        diff = df[a] - df[b]
        row.update({"ICC2_absolute": icc["ICC2"], "BA_mean_diff": diff.mean(),
                    "BA_loa_lower": diff.mean() - 1.96 * diff.std(), "BA_loa_upper": diff.mean() + 1.96 * diff.std()})
    return row
HF_VS_PV_LIVER = "Hepatic fat native vs -PV liver (z)"
cell8_rows = []
for label, a, b, fname in [
    ("Liver", "liver_native", "liver_portal_venous", "bland_altman_liver_finalqc_geno.png"),
    ("Spleen", "spleen_native", "spleen_portal_venous", "bland_altman_spleen_finalqc_geno.png"),
    ("Hepatic fat", "hepatic_fat_native", "hepatic_fat_pv", "bland_altman_hepatic_fat_finalqc_geno.png"),
]:
    row = agreement_row(person_geno, "person_id", a, b, f"{label} (final QC, genotyped, per person)",
                        ba_path=BA_OUTPUT_DIR / fname)
    cell8_rows.append({"measure": label, **row})
# Native spleen - liver vs -PV liver: different measures in different units, so
# absolute agreement (ICC2) and Bland-Altman differences are meaningless; z-score
# within each phase and report only correlation and ICC3 (consistency).
row = agreement_row(person_geno, "person_id", "hepatic_fat_native", "neg_liver_pv", HF_VS_PV_LIVER, zscore=True)
cell8_rows.append({"measure": HF_VS_PV_LIVER, **row})
cell8_df = pd.DataFrame(cell8_rows).assign(
    set="Cell 8: genotyped, per person", n_studies=n_studies_geno, n_people=len(person_geno))
# Cell 7 comparison (all studies, study level, not genotype-filtered); the -PV liver row wasn't in Cell 7, computed here
paired_finalqc["neg_liver_pv"] = -paired_finalqc["liver_portal_venous"]
row = agreement_row(paired_finalqc, "study_id", "hepatic_fat_native", "neg_liver_pv", f"{HF_VS_PV_LIVER}, Cell 7 data", zscore=True)
cell7_df = pd.DataFrame(summary_rows + [{"measure": HF_VS_PV_LIVER, **row}]).assign(
    set="Cell 7: all studies, per study", n_studies=len(paired_finalqc), n_people=paired_finalqc["person_id"].nunique())
order = ["Liver", "Spleen", "Hepatic fat", HF_VS_PV_LIVER]
summary_cols = ["measure", "set", "n_studies", "n_people", "pearson_r", "spearman_rho", "ICC2_absolute",
                "ICC3_consistency", "BA_mean_diff", "BA_loa_lower", "BA_loa_upper"]
summary_all = pd.concat([cell7_df, cell8_df], ignore_index=True)
summary_all["measure"] = pd.Categorical(summary_all["measure"], categories=order, ordered=True)
summary_all = summary_all.sort_values(["measure", "set"])[summary_cols]
summary_all["measure"] = summary_all["measure"].astype(str).replace("Hepatic fat", "Hepatic fat (spleen - liver)")
print("\nAgreement, Cell 7 vs Cell 8 (native vs PV; NaN = not applicable):")
print(summary_all.to_string(index=False, float_format=lambda x: f"{x:.3f}", na_rep="-"))


# %% Cell 9: per person, ALL people (no genotype filter) -- the middle step
# between Cell 7 (per study, all people) and Cell 8 (per person, genotyped)
import io
paired_all = pivoted_finalqc.reset_index().dropna(subset=required_cols)
person_all = paired_all.groupby("person_id")[required_cols].median().reset_index()
person_all["person_id"] = person_all["person_id"].astype(str)
person_all["hepatic_fat_native"] = person_all["spleen_native"] - person_all["liver_native"]
person_all["hepatic_fat_pv"] = person_all["spleen_portal_venous"] - person_all["liver_portal_venous"]
person_all["neg_liver_pv"] = -person_all["liver_portal_venous"]
print(f"Paired studies (all people): {len(paired_all):,} -> people after per-person median: {len(person_all):,}")
# agreement_row only returns ICC2/Bland-Altman when given a plot path; send the plot to an in-memory buffer (not saved)
cell9_rows = []
for label, a, b in [
    ("Liver", "liver_native", "liver_portal_venous"),
    ("Spleen", "spleen_native", "spleen_portal_venous"),
    ("Hepatic fat", "hepatic_fat_native", "hepatic_fat_pv"),
]:
    row = agreement_row(person_all, "person_id", a, b, f"{label} (final QC, all people, per person)",
                        ba_path=io.BytesIO())
    cell9_rows.append({"measure": label, **row})
row = agreement_row(person_all, "person_id", "hepatic_fat_native", "neg_liver_pv", HF_VS_PV_LIVER, zscore=True)
cell9_rows.append({"measure": HF_VS_PV_LIVER, **row})
cell9_df = pd.DataFrame(cell9_rows).assign(
    set="Cell 9: per person, all people", n_studies=len(paired_all), n_people=len(person_all))
set_order = ["Cell 7: all studies, per study", "Cell 9: per person, all people", "Cell 8: genotyped, per person"]
summary_3way = pd.concat([cell7_df, cell9_df, cell8_df], ignore_index=True)
summary_3way["measure"] = pd.Categorical(summary_3way["measure"], categories=order, ordered=True)
summary_3way["set"] = pd.Categorical(summary_3way["set"], categories=set_order, ordered=True)
summary_3way = summary_3way.sort_values(["measure", "set"])[summary_cols]
summary_3way["measure"] = summary_3way["measure"].astype(str).replace("Hepatic fat", "Hepatic fat (spleen - liver)")
summary_3way["set"] = summary_3way["set"].astype(str).replace({
    "Cell 7: all studies, per study": "Cell 7: per study, all people",
    "Cell 8: genotyped, per person": "Cell 8: per person, genotyped only"})
print("\nAgreement, Cell 7 -> Cell 9 -> Cell 8 (native vs PV; - = not applicable):")
print(summary_3way.to_string(index=False, float_format=lambda x: f"{x:.3f}", na_rep="-"))

