import argparse, glob, os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from manhattan_plot import ManhattanPlot

# Usage: python make_manhattan_plot_oct2026.py SL_untransformed PMBB_EUR
ap = argparse.ArgumentParser()
ap.add_argument("run", help="SL_untransformed, liver_untransformed, SL_transformed, liver_transformed")
ap.add_argument("cohort", help="PMBB_EUR or PMBB_AFR")
args = ap.parse_args()

pheno = "hepatic_fat"
cohort = args.cohort
base = "/project/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_oct2026"
run_dir = os.path.join(base, args.run)

sumstats_file = os.path.join(run_dir, cohort, "Sumstats", f"{cohort}.{pheno}.gwas.saige.gz")
annot_file = os.path.join(run_dir, "Annotations", "GWAS_biofilter_genes_rsids.csv")
output_dir = os.path.join(run_dir, "Plots_custom")

def newest(pattern):
    hits = glob.glob(pattern)
    if not hits:
        raise SystemExit(f"Not found: {pattern}")
    return max(hits, key=os.path.getmtime)

colnames_file = newest(os.path.join(run_dir, "work", "*", "*", "colnames.txt"))
pheno_table_file = newest(os.path.join(run_dir, "work", "*", "*", "pheno_summaries.csv"))
print("sumstats:", sumstats_file)
print("colnames:", colnames_file)
print("pheno table:", pheno_table_file)

# Load column mapping (same as January script)
columns_map = {}
with open(colnames_file) as fh:
    for line in fh:
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        columns_map[key.strip()] = value.strip().strip("'")
print("Column mappings:", columns_map)

# Phenotype info for the title
pheno_df = pd.read_csv(pheno_table_file)
pheno_df = pheno_df[(pheno_df["PHENO"] == pheno) & (pheno_df["COHORT"] == cohort)]
if pheno_df.empty:
    raise SystemExit(f"No row for {pheno} / {cohort} in {pheno_table_file}")
trait_type = "bin" if pheno_df["Cases"].count() != 0 else "quant"

os.makedirs(output_dir, exist_ok=True)
output_manhattan = os.path.join(output_dir, f"{cohort}.{pheno}.manhattan_vertical.png")
output_qq = os.path.join(output_dir, f"{cohort}.{pheno}.qq.png")

plot_title = f'GWAS {args.run.replace("_", " ")}, {cohort}: {pheno.replace("_", " ")}'
if trait_type == "bin":
    plot_title += f"\nCases = {pheno_df.iloc[0].Cases:,.0f}, Controls = {pheno_df.iloc[0].Controls:,.0f}"
else:
    plot_title += f"\nN = {pheno_df.iloc[0].N:,.0f}"
print("Title:", plot_title)

# Load and clean (same as January script)
mp = ManhattanPlot(sumstats_file, test_rows=None, title=plot_title)
mp.load_data()
mp.df["chromosome_noCHR"] = mp.df["chromosome"].astype(str).str.replace("chr", "").astype(int)
mp.df["variant_id"] = mp.df["variant_id"].str.replace("_", ":")
mp.clean_data(col_map={
    "chromosome_noCHR": "#CHROM",
    columns_map["POS"]: "POS",
    columns_map["MarkerID"]: "ID",
    columns_map["p.value"]: "P",
})
print("Data shape:", mp.df.shape)

# Annotations
if os.path.exists(annot_file):
    annot_df = pd.read_csv(annot_file)
    annot_df["ID"] = annot_df["Gene"]
    annot_df.iloc[:, 0] = annot_df.iloc[:, 0].str.replace("_", ":")
    mp.add_annotations(annot_df, extra_cols=["RSID"])
    print("Annotations added:", len(annot_df))
else:
    print("No annotation file found, skipping:", annot_file)

# Thin, threshold, plot (same as January script)
mp.get_thinned_data()
annot_thresh = 5E-8 if np.any(mp.thinned["P"].min() < 5E-8) else np.nanquantile(mp.thinned["P"], 10 / len(mp.thinned))
mp.update_plotting_parameters(
    vertical=True,
    sig=annot_thresh if not np.any(mp.thinned["P"] < 5E-8) else 5E-8,
    sug=annot_thresh,
    annot_thresh=annot_thresh,
    merge_genes=True,
)
mp.full_plot(
    save=output_manhattan, rep_boost=False,
    extra_cols={"beta": "BETA", "RSID": "RSID"}, number_cols=["BETA"],
    keep_chr_pos=False, with_table_bg=False, with_table_grid=False,
)
plt.clf()
print("Manhattan saved:", output_manhattan)

mp.qq_plot(save=output_qq)
print("QQ saved:", output_qq)
