# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 19:33:35 2026

@author: ammar
"""

import pandas as pd
import numpy as np
import os

# =============================================================================
# PATHS – adjust if your output folder differs
# =============================================================================
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results"
output_dir = os.path.join(base_dir, "outputs")
figure_dir = os.path.join(base_dir, "figures", "talib_manuscript")
os.makedirs(figure_dir, exist_ok=True)
table_output = os.path.join(figure_dir, "Table_S1_Q1_summary.csv")

# =============================================================================
# LOAD DATA
# =============================================================================
# Panel A summary (means and CI half-width)
summary_file = os.path.join(output_dir, "Q1_final_near_normal_summary_by_difference_type.csv")
df_summary = pd.read_csv(summary_file)

# GAM smooth terms (edf and p-value)
edf_file = os.path.join(output_dir, "Q1_final_smooth_gam_smooth_terms.csv")
edf_df = pd.read_csv(edf_file)

# =============================================================================
# FILTER ONLY T:ET SMOOTH TERMS AND PARSE ECOSYSTEM / DIFFERENCE
# =============================================================================
# Keep only terms that contain "s(Trans_ratio)"
edf_df = edf_df[edf_df["term"].str.contains(r"s\(Trans_ratio\)", regex=True)].copy()

def parse_smooth_term(term):
    # Example: "s(Trans_ratio):diff_ecosystemWUE_ET_minus_WUE_T__Upland"
    clean = term.replace("s(Trans_ratio):diff_ecosystem", "")
    diff_type, ecosystem = clean.split("__")
    diff_map = {
        "WUE_ET_minus_WUE_T": "WUE_ET - WUE_T",
        "WUE_ET_minus_WUE_E": "WUE_ET - WUE_E",
        "WUE_E_minus_WUE_T": "WUE_E - WUE_T",
    }
    return pd.Series({
        "difference_label": diff_map[diff_type],
        "ecosystem_label": ecosystem
    })

edf_parsed = edf_df["term"].apply(parse_smooth_term)
edf_df = pd.concat([edf_df, edf_parsed], axis=1)

# Keep relevant columns: edf and p‑value
p_col = [c for c in edf_df.columns if "p" in c.lower()][0]
edf_df = edf_df[["ecosystem_label", "difference_label", "edf", p_col]]
edf_df.rename(columns={p_col: "p_value"}, inplace=True)

# =============================================================================
# MERGE WITH SUMMARY DATA
# =============================================================================
df_summary = df_summary.rename(columns={"water_class": "ecosystem_label"})
table = df_summary.merge(edf_df, on=["ecosystem_label", "difference_label"], how="left")

# Build 95% CI interval as "lower to upper"
table["ci_lower"] = table["mean_difference"] - table["ci95_difference"]
table["ci_upper"] = table["mean_difference"] + table["ci95_difference"]
table["95% CI"] = table.apply(
    lambda r: f"{r['ci_lower']:.2f} to {r['ci_upper']:.2f}",
    axis=1
)

# Select and rename columns for the final table
table = table[[
    "ecosystem_label",
    "difference_label",
    "mean_difference",
    "95% CI",
    "edf",
    "p_value"
]]
table.rename(columns={
    "ecosystem_label": "Ecosystem",
    "difference_label": "WUE metric difference",
    "mean_difference": "Mean difference",
    "edf": "GAM edf",
    "p_value": "Smooth p-value"
}, inplace=True)

# Round numeric columns
table["Mean difference"] = table["Mean difference"].round(2)
table["GAM edf"] = table["GAM edf"].round(2)
table["Smooth p-value"] = table["Smooth p-value"].apply(
    lambda x: f"{x:.3f}" if x >= 0.001 else "<0.001"
)

# =============================================================================
# SORT IN A LOGICAL ORDER
# =============================================================================
ecosystem_order = ["Upland", "Freshwater", "Saline"]
diff_order = ["WUE_ET - WUE_T", "WUE_ET - WUE_E", "WUE_E - WUE_T"]

table["Ecosystem"] = pd.Categorical(table["Ecosystem"], categories=ecosystem_order, ordered=True)
table["WUE metric difference"] = pd.Categorical(table["WUE metric difference"], categories=diff_order, ordered=True)
table = table.sort_values(["Ecosystem", "WUE metric difference"]).reset_index(drop=True)

# =============================================================================
# SAVE CSV
# =============================================================================
table.to_csv(table_output, index=False)
print(f"✅ Supplementary table saved to: {table_output}")