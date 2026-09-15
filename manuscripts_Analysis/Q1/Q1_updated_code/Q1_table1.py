# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 19:33:35 2026

@author: ammar
"""

import pandas as pd
import numpy as np
import os

# =============================================================================
# PATHS – unchanged
# =============================================================================
base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results"
output_dir = os.path.join(base_dir, "outputs")
figure_dir = os.path.join(base_dir, "figures", "talib_manuscript")
os.makedirs(figure_dir, exist_ok=True)
table_output = os.path.join(figure_dir, "Table_S1_Q1_summary.csv")

# =============================================================================
# LOAD DATA
# =============================================================================
summary_file = os.path.join(output_dir, "Q1_final_near_normal_summary_by_difference_type.csv")
edf_file = os.path.join(output_dir, "Q1_final_smooth_gam_smooth_terms.csv")

df_summary = pd.read_csv(summary_file)
edf_df = pd.read_csv(edf_file)

# =============================================================================
# FILTER SUMMARY TO SINGLE DIFFERENCE
# =============================================================================
SINGLE_DIFF_LABEL = "WUE_ET - WUE_T"
SINGLE_DIFF_RAW = "WUE_ET_minus_WUE_T"

df_summary = df_summary[df_summary["difference_label"] == SINGLE_DIFF_LABEL].copy()
if df_summary.empty:
    raise ValueError(
        f"No rows in summary CSV with difference_label == '{SINGLE_DIFF_LABEL}'. "
        f"Found: {sorted(df_summary['difference_label'].unique().tolist()) if 'difference_label' in df_summary.columns else 'no column'}"
    )

# =============================================================================
# FILTER AND PARSE THE THREE ECOSYSTEM-SPECIFIC T:ET SMOOTH TERMS
# =============================================================================
# Keep only s(Trans_ratio) terms for the WUE_ET_minus_WUE_T difference.
# Excludes s(site_name) and s(month_f) random-effect smooths.
edf_df = edf_df[
    edf_df["term"].astype(str).str.contains(r"s\(Trans_ratio\)", regex=True, na=False) &
    edf_df["term"].astype(str).str.contains(SINGLE_DIFF_RAW, regex=False, na=False)
].copy()

if edf_df.empty:
    raise ValueError(
        f"No s(Trans_ratio) smooth terms found for '{SINGLE_DIFF_RAW}'. "
        "Check Q1_final_smooth_gam_smooth_terms.csv."
    )

def parse_smooth_term(term):
    """
    Parse an ecosystem-specific smooth term of the form:
        s(Trans_ratio):diff_ecosystemWUE_ET_minus_WUE_T__<Ecosystem>
    Returns a Series with the ecosystem label.
    """
    clean = term.replace("s(Trans_ratio):diff_ecosystem", "")
    parts = clean.split("__")
    if len(parts) != 2:
        raise ValueError(f"Unexpected smooth term format: {term}")
    diff_type, ecosystem = parts
    if diff_type != SINGLE_DIFF_RAW:
        raise ValueError(
            f"Unexpected difference type in smooth term '{term}'. "
            f"Expected only '{SINGLE_DIFF_RAW}'."
        )
    return pd.Series({"ecosystem_label": ecosystem})

edf_parsed = edf_df["term"].apply(parse_smooth_term)
edf_df = pd.concat([edf_df, edf_parsed], axis=1)

# Locate F and p-value columns robustly
f_col_candidates = [c for c in edf_df.columns if c.strip().lower() in ("f", "f statistic", "f_statistic")]
if not f_col_candidates:
    raise ValueError(f"Could not find F column in smooth terms. Columns: {edf_df.columns.tolist()}")
f_col = f_col_candidates[0]

p_col_candidates = [c for c in edf_df.columns if "p" in c.lower() and "value" in c.lower()]
if not p_col_candidates:
    p_col_candidates = [c for c in edf_df.columns if c.lower().startswith("p")]
if not p_col_candidates:
    raise ValueError(f"Could not find p-value column in smooth terms. Columns: {edf_df.columns.tolist()}")
p_col = p_col_candidates[0]

edf_df = edf_df[["ecosystem_label", "edf", f_col, p_col]].rename(
    columns={f_col: "F_stat", p_col: "p_value"}
)

# =============================================================================
# MERGE WITH SUMMARY DATA
# =============================================================================
df_summary = df_summary.rename(columns={"water_class": "ecosystem_label"})

table = df_summary.merge(edf_df, on="ecosystem_label", how="left")

# Build 95% CI interval as "lower to upper"
table["ci_lower"] = table["mean_difference"] - table["ci95_difference"]
table["ci_upper"] = table["mean_difference"] + table["ci95_difference"]
table["95% CI"] = table.apply(
    lambda r: f"{r['ci_lower']:.2f} to {r['ci_upper']:.2f}",
    axis=1
)

# =============================================================================
# SELECT, RENAME, SORT
# =============================================================================
table = table[[
    "ecosystem_label",
    "n_months",
    "n_sites",
    "mean_difference",
    "95% CI",
    "edf",
    "F_stat",
    "p_value",
]].rename(columns={
    "ecosystem_label": "Ecosystem",
    "n_months":        "Site-months",
    "n_sites":         "Sites",
    "mean_difference": "Mean WUE_ET \u2212 WUE_T",
    "edf":             "GAM edf",
    "F_stat":          "GAM F",
    "p_value":         "Smooth p-value",
})

# Round numeric columns
table["Mean WUE_ET \u2212 WUE_T"] = table["Mean WUE_ET \u2212 WUE_T"].round(2)
table["GAM edf"] = table["GAM edf"].round(2)
table["GAM F"]   = table["GAM F"].round(2)
table["Smooth p-value"] = table["Smooth p-value"].apply(
    lambda x: "<0.001" if pd.notna(x) and x < 0.001 else (f"{x:.3f}" if pd.notna(x) else "NA")
)

# Sort by ecosystem order
ecosystem_order = ["Upland", "Freshwater", "Saline"]
table["Ecosystem"] = pd.Categorical(table["Ecosystem"], categories=ecosystem_order, ordered=True)
table = table.sort_values("Ecosystem").reset_index(drop=True)

# =============================================================================
# FINAL VALIDATION CHECKS
# =============================================================================
assert len(table) == 3, f"Expected 3 rows, got {len(table)}"
assert set(table["Ecosystem"]) == {"Upland", "Freshwater", "Saline"}, \
    f"Unexpected ecosystems: {set(table['Ecosystem'])}"
# No WUE_E text anywhere in the final table
table_as_str = table.astype(str).apply(lambda col: col.str.cat(sep=" ")).str.cat(sep=" ")
assert "WUE_E" not in table_as_str or "WUE_ET" in table_as_str.replace("WUE_E", "WUE_ET"), \
    "Final table contains a WUE_E reference"
# Simpler & stricter check: no substring "WUE_E " or "WUE_E-" or "WUE_E," etc.
for forbidden in ["WUE_ET - WUE_E", "WUE_E - WUE_T", "WUE_ET_minus_WUE_E", "WUE_E_minus_WUE_T"]:
    assert forbidden not in table_as_str, f"Forbidden WUE_E text found: {forbidden}"
assert table["GAM edf"].notna().all(), "Some GAM edf values are missing"
assert table["GAM F"].notna().all(),   "Some GAM F values are missing"
assert table["Smooth p-value"].notna().all(), "Some smooth p-values are missing"

# =============================================================================
# PRINT AND SAVE
# =============================================================================
print("=" * 90)
print("FINAL SUPPLEMENTARY TABLE (WUE_ET - WUE_T only)")
print("=" * 90)
print(table.to_string(index=False))
print("=" * 90)

table.to_csv(table_output, index=False)
print(f"\n✅ Supplementary table saved to: {table_output}")