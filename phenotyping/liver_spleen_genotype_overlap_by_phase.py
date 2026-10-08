
#!/usr/bin/env python3
"""
Reports how many genotyped patients have usable liver and spleen IDP data,
broken out by contrast phase (native / any contrast / portal venous).
 
Run with:
    python3 liver_spleen_genotype_overlap_by_phase.py
"""
import pandas as pd
 
IDP_FILE = "/static/PMBB/PMBB-Release-2026-4.0/Imaging/CT/idp/PMBB-Release-2026-4.0_imaging_ct_idp.txt"
GENOTYPE = "/static/PMBB/PMBB-Release-2026-4.0/Imputed/chunked_bed_files/PMBB-Release-2026-4.0_genetic_imputed.463_chr8_chunk1_149284_2198252.fam"
 
CONTRAST_PHASES = {"portal_venous", "arterial_early", "arterial_late"}
VALID_PHASES = {"native"} | CONTRAST_PHASES  # excludes "missing"
 
idp = pd.read_csv(IDP_FILE, sep="\t")
idp = idp[idp["phase"].isin(VALID_PHASES)]
 
geno = pd.read_csv(
    GENOTYPE, header=None, sep=r'\s+', names=["IID", "PMBB_ID", "father", "mother", "sex", "pheno"]
)
genotyped_ids = set(geno["PMBB_ID"])
 
idp_geno = idp[idp["person_id"].isin(genotyped_ids)]
liver_idp = idp_geno[idp_geno["anatomical_structure"] == "liver"]
spleen_idp = idp_geno[idp_geno["anatomical_structure"] == "spleen"]
 
# --- Liver only ---
liver_phase_sets = liver_idp.groupby("person_id")["phase"].apply(set)
liver_native = liver_phase_sets.apply(lambda p: "native" in p)
liver_contrast = liver_phase_sets.apply(lambda p: bool(p & CONTRAST_PHASES))
print(f"Liver -- genotyped patients w/ native: {liver_native.sum():,} / {len(liver_phase_sets):,}")
print(f"Liver -- genotyped patients w/ contrast: {liver_contrast.sum():,} / {len(liver_phase_sets):,}")
 
# --- Spleen only ---
spleen_phase_sets = spleen_idp.groupby("person_id")["phase"].apply(set)
spleen_native = spleen_phase_sets.apply(lambda p: "native" in p)
spleen_contrast = spleen_phase_sets.apply(lambda p: bool(p & CONTRAST_PHASES))
print(f"Spleen -- genotyped patients w/ native: {spleen_native.sum():,} / {len(spleen_phase_sets):,}")
print(f"Spleen -- genotyped patients w/ contrast: {spleen_contrast.sum():,} / {len(spleen_phase_sets):,}")
 
# --- The number that matters: same accession, same phase, BOTH organs present ---
liver_pairs = liver_idp[["person_id", "accession", "phase"]].drop_duplicates()
spleen_pairs = spleen_idp[["person_id", "accession", "phase"]].drop_duplicates()
matched = liver_pairs.merge(spleen_pairs, on=["person_id", "accession", "phase"])
 
matched_native = matched[matched["phase"] == "native"]["person_id"].nunique()
matched_contrast = matched[matched["phase"].isin(CONTRAST_PHASES)]["person_id"].nunique()
matched_pv = matched[matched["phase"] == "portal_venous"]["person_id"].nunique()
matched_any = matched["person_id"].nunique()
 
print(f"\nPatients with >=1 study where BOTH liver+spleen are NATIVE: {matched_native:,}")
print(f"Patients with >=1 study where BOTH liver+spleen are ANY CONTRAST: {matched_contrast:,}")
print(f"Patients with >=1 study where BOTH liver+spleen are PORTAL VENOUS: {matched_pv:,}")
print(f"Total genotyped patients with any usable liver+spleen match: {matched_any:,}")
 
# --- Native + portal venous combined ----
nat_ids = set(matched.loc[matched["phase"] == "native", "person_id"])
pv_ids = set(matched.loc[matched["phase"] == "portal_venous", "person_id"])
print(f"\nNative OR PV (combined, unique people): {len(nat_ids | pv_ids):,}")
print(f"Both native AND PV:                      {len(nat_ids & pv_ids):,}")
print(f"PV only (gained by adding PV):           {len(pv_ids - nat_ids):,}")
 
# Liver alone (the Derstine paper favors liver alone on contrast)
liv_nat = set(liver_idp.loc[liver_idp["phase"] == "native", "person_id"])
liv_pv = set(liver_idp.loc[liver_idp["phase"] == "portal_venous", "person_id"])
print(f"\nLiver alone -- native OR PV combined:    {len(liv_nat | liv_pv):,}")
print(f"Liver alone -- PV only (gained):          {len(liv_pv - liv_nat):,}")
 
 
# --- Native + portal venous combined with QCed file ---
# Same counts as above, repeated on the clipping-QC'd cache from the clipping
# diagnostic script (cell 10). Note: clipping QC only (max >= 3071 or
# min <= -32768 removed); the +/-6 SD and plausibility QC have NOT been
# applied yet, so final counts may drop slightly.
print("\n" + "=" * 70)
print("AFTER CLIPPING QC (idp_liver_spleen_clean.pkl)")
print("=" * 70)
 
clean = pd.read_pickle("idp_liver_spleen_clean.pkl")
clean = clean[clean["phase"].isin(VALID_PHASES)]
clean_geno = clean[clean["person_id"].isin(genotyped_ids)]
clean_liver = clean_geno[clean_geno["anatomical_structure"] == "liver"]
clean_spleen = clean_geno[clean_geno["anatomical_structure"] == "spleen"]
 
# Same accession, same phase, BOTH organs present
clean_liver_pairs = clean_liver[["person_id", "accession", "phase"]].drop_duplicates()
clean_spleen_pairs = clean_spleen[["person_id", "accession", "phase"]].drop_duplicates()
clean_matched = clean_liver_pairs.merge(clean_spleen_pairs, on=["person_id", "accession", "phase"])
 
clean_nat_ids = set(clean_matched.loc[clean_matched["phase"] == "native", "person_id"])
clean_pv_ids = set(clean_matched.loc[clean_matched["phase"] == "portal_venous", "person_id"])
 
print(f"Patients with >=1 study where BOTH liver+spleen are NATIVE: {len(clean_nat_ids):,}")
print(f"Patients with >=1 study where BOTH liver+spleen are PORTAL VENOUS: {len(clean_pv_ids):,}")
print(f"\nNative OR PV (combined, unique people): {len(clean_nat_ids | clean_pv_ids):,}")
print(f"Both native AND PV:                      {len(clean_nat_ids & clean_pv_ids):,}")
print(f"PV only (gained by adding PV):           {len(clean_pv_ids - clean_nat_ids):,}")
 
# Liver alone
clean_liv_nat = set(clean_liver.loc[clean_liver["phase"] == "native", "person_id"])
clean_liv_pv = set(clean_liver.loc[clean_liver["phase"] == "portal_venous", "person_id"])
print(f"\nLiver alone -- native OR PV combined:    {len(clean_liv_nat | clean_liv_pv):,}")
print(f"Liver alone -- PV only (gained):          {len(clean_liv_pv - clean_liv_nat):,}")


# --- Venn diagrams: native vs PV, before and after clipping QC ---
# Uses only matplotlib (no extra packages). Saves two PNGs to the current folder.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
 
 
def draw_venn(native_set, pv_set, title, output_path):
    native_only = len(native_set - pv_set)
    both = len(native_set & pv_set)
    pv_only = len(pv_set - native_set)
    combined = len(native_set | pv_set)
 
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.add_patch(Circle((-0.5, 0), 1.0, color="#1D9E75", alpha=0.45, label="Native"))
    ax.add_patch(Circle((0.5, 0), 1.0, color="#7F77DD", alpha=0.45, label="Portal venous"))
 
    # Counts in each region
    ax.text(-1.0, 0, f"{native_only:,}\nNative only", ha="center", va="center", fontsize=11)
    ax.text(0, 0, f"{both:,}\nBoth", ha="center", va="center", fontsize=11)
    ax.text(1.0, 0, f"{pv_only:,}\nPV only", ha="center", va="center", fontsize=11)
 
    # Circle totals above, combined total below
    ax.text(-0.75, 1.2, f"Native: {len(native_set):,}", ha="center", fontsize=12, fontweight="bold")
    ax.text(0.75, 1.2, f"Portal venous: {len(pv_set):,}", ha="center", fontsize=12, fontweight="bold")
    ax.text(0, -1.3, f"Combined (native or PV): {combined:,}", ha="center", fontsize=12, fontweight="bold")
 
    ax.set_title(title, fontsize=12)
    ax.set_xlim(-1.8, 1.8)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Saved {output_path}")
 
 
print()
draw_venn(
    nat_ids, pv_ids,
    "Genotyped people with liver + spleen in the same study\nBefore clipping QC",
    "venn_native_pv_before_clipping.png",
)
draw_venn(
    clean_nat_ids, clean_pv_ids,
    "Genotyped people with liver + spleen in the same study\nAfter clipping QC",
    "venn_native_pv_after_clipping.png",
)
 