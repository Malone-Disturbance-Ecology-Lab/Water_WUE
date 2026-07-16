# -*- coding: utf-8 -*-
"""
Q4_spatial_upscaling_summary_table.py

Creates a summary table for Q4 spatial upscaling from input CSVs.
No hard-coded values; all numbers computed from data.
Outputs: CSV only, with zeros formatted as 0.00 (or 0.000 for resistance).
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================================
# PATHS
# ============================================================================
base_dir = Path(r"M:\Research\WUE_CUE\WUE_manuscript_version6\upscaling")

monthly_file = base_dir / "EDI_response_v2_upland_monthly_coast_summary.csv"
annual_file = base_dir / "EDI_response_v2_upland_annual_coast_summary.csv"
impacted_file = base_dir / "EDI_response_communication_annual_impacted_area_by_coast_2000_2025.csv"
resistance_file = base_dir / "spatial_performance_metric" / "spatial_resistance_analogue_two_panel_summary.csv"

output_dir = base_dir / "talib_publication"
output_dir.mkdir(parents=True, exist_ok=True)
output_csv = output_dir / "Q4_spatial_upscaling_summary_table.csv"

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def format_value(val, decimals=2, threshold=0.01, force_zero_two_decimals=True):
    """Format with 2 decimals unless |val| < 0.01 and !=0, then 4 decimals.
       Zero values always get '0.00' (or '0.000' if decimals=3).
    """
    if pd.isna(val):
        return np.nan
    if val == 0 and force_zero_two_decimals:
        return f"{0:.{decimals}f}"
    if abs(val) < threshold and val != 0:
        return f"{val:.4f}"
    else:
        return f"{val:.{decimals}f}"

# ============================================================================
# LOAD DATA
# ============================================================================
print("Loading input files...")
df_monthly = pd.read_csv(monthly_file)
df_annual = pd.read_csv(annual_file)
df_impacted = pd.read_csv(impacted_file)
df_resistance = pd.read_csv(resistance_file) if resistance_file.exists() else None
if df_resistance is None:
    print("Warning: resistance file not found, will compute from dry WUE change.")

# ============================================================================
# COAST NAME MAPPING (long to short)
# ============================================================================
coast_mapping = {
    "AK Coast": "Alaska",
    "Pacific Coast": "Pacific",
    "Gulf Coast": "Gulf",
    "Atlantic Coast": "Atlantic"
}

def standardize_coast(df, col="coast_region"):
    if col in df.columns:
        df[col] = df[col].replace(coast_mapping)
    return df

df_monthly = standardize_coast(df_monthly)
df_annual = standardize_coast(df_annual)
df_impacted = standardize_coast(df_impacted)
if df_resistance is not None:
    df_resistance = standardize_coast(df_resistance)

# ============================================================================
# COMPUTE SUMMARY STATISTICS
# ============================================================================
# 1. Monthly summary: mean SPEI and cumulative SPEI
monthly_agg = df_monthly.groupby("coast_region").agg(
    mean_SPEI=("mean_SPEI3", "mean"),
    cumul_SPEI=("mean_SPEI3", "sum")
).reset_index()

# 2. Annual summary: dry area, dry-month WUE_T change, dry pixel crossings
annual_agg = df_annual.groupby("coast_region").agg(
    dry_area=("dry_pct_pixels", "mean"),
    dry_WUE_change=("dry_mean_WUE_T_pct_change", "mean"),
    dry_cross_pos=("dry_pct_pixels_strong_increase", "mean"),
    dry_cross_neg=("dry_pct_pixels_strong_decrease", "mean")
).reset_index()

# 3. Impacted area summary
impacted_agg = df_impacted.groupby("coast_region").agg(
    neg_area=("annual_mean_negative_impacted_area_pct", "mean"),
    pos_area=("annual_mean_positive_impacted_area_pct", "mean"),
    any_area=("annual_mean_any_impacted_area_pct", "mean")
).reset_index()

# 4. Resistance analogue – pull from pixel panel
resist_agg = None
if df_resistance is not None:
    if {"panel", "coast_region", "mean"}.issubset(df_resistance.columns):
        df_res_pixel = df_resistance[df_resistance["panel"].astype(str).str.lower() == "pixel"].copy()
        if not df_res_pixel.empty:
            resist_agg = df_res_pixel[["coast_region", "mean"]].copy()
            resist_agg.rename(columns={"mean": "resistance"}, inplace=True)
    elif "mean_resistance_area_weighted" in df_resistance.columns:
        resist_agg = df_resistance[["coast_region", "mean_resistance_area_weighted"]].copy()
        resist_agg.rename(columns={"mean_resistance_area_weighted": "resistance"}, inplace=True)
    elif "mean_resistance_unweighted" in df_resistance.columns:
        resist_agg = df_resistance[["coast_region", "mean_resistance_unweighted"]].copy()
        resist_agg.rename(columns={"mean_resistance_unweighted": "resistance"}, inplace=True)

if resist_agg is None:
    print("Warning: resistance summary not usable; calculating resistance from dry WUE change.")

# ============================================================================
# MERGE ALL
# ============================================================================
result = monthly_agg.merge(annual_agg, on="coast_region", how="outer")
result = result.merge(impacted_agg, on="coast_region", how="outer")
if resist_agg is not None:
    result = result.merge(resist_agg, on="coast_region", how="outer")
else:
    result["resistance"] = 1 + result["dry_WUE_change"] / 100

# ============================================================================
# ADD INTERPRETATION
# ============================================================================
def assign_interpretation(coast):
    if coast == "Pacific":
        return "Strong positive WUE_T adjustment during SPEI-3 dry exposure"
    elif coast == "Atlantic":
        return "Weaker positive response with localized +5% dry-month crossings"
    elif coast == "Alaska":
        return "Modest SPEI-3 dry-month decline within ±5% threshold band"
    elif coast == "Gulf":
        return "Dry exposure occurred, but modeled WUE_T response remained below threshold"
    else:
        return ""

result["Q4 interpretation"] = result["coast_region"].apply(assign_interpretation)

# ============================================================================
# SELECT AND RENAME
# ============================================================================
final_columns = [
    "coast_region",
    "mean_SPEI",
    "cumul_SPEI",
    "dry_area",
    "neg_area",
    "pos_area",
    "any_area",
    "dry_WUE_change",
    "dry_cross_pos",
    "dry_cross_neg",
    "resistance",
    "Q4 interpretation"
]
result = result[final_columns]

result.rename(columns={
    "coast_region": "Coastline",
    "mean_SPEI": "Mean growing-season SPEI-3",
    "cumul_SPEI": "Cumulative growing-season SPEI-3",
    "dry_area": "Dry area (%)",
    "neg_area": "Negative ≤−5% response area (%)",   # <--- FIXED HERE
    "pos_area": "Positive ≥5% response area (%)",
    "any_area": "Any ±5% response area (%)",
    "dry_WUE_change": "Dry-month WUE_T change (%)",
    "dry_cross_pos": "Dry pixels crossing +5% (%)",
    "dry_cross_neg": "Dry pixels crossing −5% (%)",
    "resistance": "Resistance analogue"
}, inplace=True)

# Set coast order
coast_order = ["Alaska", "Pacific", "Gulf", "Atlantic"]
result["Coastline"] = pd.Categorical(result["Coastline"], categories=coast_order, ordered=True)
result = result.sort_values("Coastline").reset_index(drop=True)

# Convert Coastline to string (for safe fillna later)
result["Coastline"] = result["Coastline"].astype(str)

# ============================================================================
# FORMAT NUMBERS – with explicit zero handling
# ============================================================================
for col in result.columns:
    if col in ["Coastline", "Q4 interpretation"]:
        continue
    if col == "Resistance analogue":
        result[col] = result[col].apply(lambda x: format_value(x, decimals=3, threshold=0.01, force_zero_two_decimals=True) if pd.notna(x) else "")
    else:
        result[col] = result[col].apply(lambda x: format_value(x, decimals=2, threshold=0.01, force_zero_two_decimals=True) if pd.notna(x) else "")

# Fill any remaining NaN with empty string
for col in result.columns:
    if result[col].dtype == "object":
        result[col] = result[col].fillna("")
    else:
        result[col] = result[col].fillna("")

# ============================================================================
# SAVE CSV
# ============================================================================
result.to_csv(output_csv, index=False, encoding='utf-8-sig')
print(f"CSV saved: {output_csv}")

# ============================================================================
# DISPLAY TABLE FOR CHECKING
# ============================================================================
print("\n" + "="*80)
print("Q4 SPATIAL UPSCALING SUMMARY TABLE")
print("="*80)
print(result.to_string(index=False))