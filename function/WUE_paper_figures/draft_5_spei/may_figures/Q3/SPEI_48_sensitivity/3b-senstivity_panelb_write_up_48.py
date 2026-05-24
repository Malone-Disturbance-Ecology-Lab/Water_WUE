# -*- coding: utf-8 -*-
"""
Created on Tue May 12 13:40:27 2026

@author: ammar
"""
# =============================================================================
# DIAGNOSTIC: SPEI-48 DISTRIBUTION BY COASTAL REGION
# Purpose:
#   Test whether coastal regions differ in their SPEI-48 exposure distributions.
#   This addresses the concern that regional WUE_T sensitivity comparisons may be
#   affected by unequal hydroclimatic exposure.
#
# Main tests:
#   1. Kruskal-Wallis tests across coastal regions for SPEI-48 summaries:
#      min, q25, median, q75, max, range, IQR, mean, dry/wet fractions.
#   2. Pairwise Mann-Whitney tests only if Kruskal-Wallis is significant.
#   3. Saves site-level diagnostic table and summary table.
# =============================================================================

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# FILE PATHS
# =============================================================================

BASE_DIR = Path(r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results")

SLOPES_FILE = BASE_DIR / "Q3_analysis" / "Q3_WUET_SPEI48_slopes.csv"
MONTHLY_FILE = BASE_DIR / "monthly_data_after_outlier_removal.csv"
NN_MEDIANS_FILE = BASE_DIR / "wue_site_level_NN_medians_SPEI_1.csv"

OUTPUT_DIR = BASE_DIR / "Q3_analysis" / "diagnostic_SPEI48_distribution_by_coast"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 90)
print("DIAGNOSTIC: SPEI-48 DISTRIBUTION BY COASTAL REGION")
print("=" * 90)

# =============================================================================
# STEP 1: LOAD SPEI-48 SLOPES / METADATA
# =============================================================================

df_slopes = pd.read_csv(SLOPES_FILE)

required_cols = ["site_name", "slope_theilsen", "coast_region_analysis"]
missing = [c for c in required_cols if c not in df_slopes.columns]
if missing:
    raise ValueError(f"Missing required columns in slopes file: {missing}")

print("\n[STEP 1] Loaded SPEI-48 slopes file")
print(f"  Sites: {len(df_slopes)}")
print(f"  Columns: {df_slopes.columns.tolist()}")

# =============================================================================
# STEP 2: LOAD MONTHLY DATA AND APPLY SAME FILTERS AS SPEI-48 WORKFLOW
# =============================================================================

df_monthly = pd.read_csv(MONTHLY_FILE)
nn_data = pd.read_csv(NN_MEDIANS_FILE)

# Strict triple intersection sites
wue_sites = set(nn_data[nn_data["WUE_Metric"] == "WUE"]["site_name"].dropna().unique())
eva_sites = set(nn_data[nn_data["WUE_Metric"] == "WUE_eva"]["site_name"].dropna().unique())
tra_sites = set(nn_data[nn_data["WUE_Metric"] == "WUE_tra"]["site_name"].dropna().unique())
shared_sites = wue_sites.intersection(eva_sites).intersection(tra_sites)

df_monthly = df_monthly[df_monthly["site_name"].isin(shared_sites)].copy()

# Strict month filter: same valid WUE metrics
strict_mask = (
    df_monthly["WUE"].notna() &
    df_monthly["WUE_eva"].notna() &
    df_monthly["WUE_tra"].notna()
)
df_monthly = df_monthly[strict_mask].copy()

# Trans_ratio filter
if "Trans_ratio" in df_monthly.columns:
    before = len(df_monthly)
    df_monthly = df_monthly[df_monthly["Trans_ratio"] > 0].copy()
    print(f"\n[STEP 2] Trans_ratio > 0 filter removed {before - len(df_monthly)} rows")

# SPEI-48 and WUE_T valid rows
df_monthly = df_monthly.dropna(subset=["SPEI_48", "WUE_tra"]).copy()

print("\n[STEP 2] Monthly data after filters")
print(f"  Rows: {len(df_monthly)}")
print(f"  Sites: {df_monthly['site_name'].nunique()}")

# =============================================================================
# STEP 3: COMPUTE SITE-LEVEL SPEI-48 DISTRIBUTION METRICS
# =============================================================================

def fraction_between(x, low, high):
    return ((x >= low) & (x <= high)).mean()

site_rows = []

for site, g in df_monthly.groupby("site_name"):
    x = g["SPEI_48"].dropna()
    
    if len(x) >= 5:
        q25 = x.quantile(0.25)
        q75 = x.quantile(0.75)
        
        site_rows.append({
            "site_name": site,
            "n_months": len(x),
            "spei48_min": x.min(),
            "spei48_q25": q25,
            "spei48_median": x.median(),
            "spei48_q75": q75,
            "spei48_max": x.max(),
            "spei48_range": x.max() - x.min(),
            "spei48_iqr": q75 - q25,
            "spei48_mean": x.mean(),
            "frac_dry_neg": (x < 0).mean(),
            "frac_wet_pos": (x > 0).mean(),
            "frac_extreme_dry": (x <= -1).mean(),
            "frac_extreme_wet": (x >= 1).mean(),
            "frac_near_normal": fraction_between(x, -0.5, 0.5),
        })

df_spei = pd.DataFrame(site_rows)

print("\n[STEP 3] Site-level SPEI-48 distribution metrics computed")
print(f"  Sites: {len(df_spei)}")

# =============================================================================
# STEP 4: MERGE WITH COASTAL REGION METADATA
# =============================================================================

df = df_slopes.merge(df_spei, on="site_name", how="inner", suffixes=("_slopesfile", "_fresh"))

# If duplicate columns exist, use freshly computed monthly SPEI metrics
for col in ["spei48_min", "spei48_max", "spei48_range"]:
    fresh_col = f"{col}_fresh"
    old_col = f"{col}_slopesfile"
    
    if fresh_col in df.columns:
        df[col] = df[fresh_col]
    elif old_col in df.columns:
        df[col] = df[old_col]

if "spei48_median_fresh" in df.columns:
    df["spei48_median"] = df["spei48_median_fresh"]
elif "spei48_median_slopesfile" in df.columns:
    df["spei48_median"] = df["spei48_median_slopesfile"]

if "n_months_fresh" in df.columns:
    df["n_months"] = df["n_months_fresh"]
elif "n_months_slopesfile" in df.columns:
    df["n_months"] = df["n_months_slopesfile"]

df = df.dropna(subset=["coast_region_analysis"]).copy()

print("\n[STEP 4] Merged site-level SPEI metrics with coastal region")
print(f"  Sites merged: {len(df)}")
print("\nCoastal region counts:")
print(df["coast_region_analysis"].value_counts().to_string())

# =============================================================================
# STEP 5: REGIONAL SUMMARY TABLE
# =============================================================================

metrics = [
    "spei48_min",
    "spei48_q25",
    "spei48_median",
    "spei48_q75",
    "spei48_max",
    "spei48_range",
    "spei48_iqr",
    "spei48_mean",
    "frac_dry_neg",
    "frac_wet_pos",
    "frac_extreme_dry",
    "frac_extreme_wet",
    "frac_near_normal",
    "n_months"
]

region_summary_rows = []

print("\n" + "=" * 90)
print("REGIONAL SPEI-48 DISTRIBUTION SUMMARY")
print("=" * 90)

for region, g in df.groupby("coast_region_analysis"):
    row = {
        "coast_region_analysis": region,
        "N_sites": len(g),
    }
    
    print(f"\n{region} (N={len(g)} sites):")
    
    for metric in metrics:
        med = g[metric].median()
        q25 = g[metric].quantile(0.25)
        q75 = g[metric].quantile(0.75)
        row[f"{metric}_median"] = med
        row[f"{metric}_q25"] = q25
        row[f"{metric}_q75"] = q75
        
        if metric in ["spei48_median", "spei48_range", "spei48_iqr", "frac_dry_neg", "frac_wet_pos"]:
            print(f"  {metric}: median={med:.3f}, IQR=[{q25:.3f}, {q75:.3f}]")
    
    region_summary_rows.append(row)

region_summary_df = pd.DataFrame(region_summary_rows)
region_summary_csv = OUTPUT_DIR / "regional_spei48_distribution_summary.csv"
region_summary_df.to_csv(region_summary_csv, index=False)

print(f"\nSaved regional summary table: {region_summary_csv}")

# =============================================================================
# STEP 6: KRUSKAL-WALLIS TESTS ACROSS COASTAL REGIONS
# =============================================================================

print("\n" + "=" * 90)
print("KRUSKAL-WALLIS TESTS: DO SPEI-48 DISTRIBUTIONS DIFFER BY COASTAL REGION?")
print("=" * 90)

kw_results = []

regions = df["coast_region_analysis"].dropna().unique().tolist()

for metric in metrics:
    group_values = [
        df[df["coast_region_analysis"] == region][metric].dropna().values
        for region in regions
    ]
    
    # Keep only groups with at least 2 values
    valid_group_values = [vals for vals in group_values if len(vals) >= 2]
    
    if len(valid_group_values) >= 2:
        h_stat, p_val = kruskal(*valid_group_values)
        epsilon_sq = h_stat / (len(df) - 1)
        
        kw_results.append({
            "metric": metric,
            "H": h_stat,
            "p_value": p_val,
            "epsilon_sq": epsilon_sq,
            "significant_p05": p_val < 0.05
        })
        
        sig = " *" if p_val < 0.05 else ""
        print(f"\n{metric}:")
        print(f"  H = {h_stat:.3f}")
        print(f"  p = {p_val:.4f}{sig}")
        print(f"  epsilon² = {epsilon_sq:.3f}")

kw_results_df = pd.DataFrame(kw_results)
kw_csv = OUTPUT_DIR / "regional_spei48_distribution_kruskal_results.csv"
kw_results_df.to_csv(kw_csv, index=False)

print(f"\nSaved Kruskal-Wallis results: {kw_csv}")

# =============================================================================
# STEP 7: PAIRWISE TESTS FOR SIGNIFICANT METRICS ONLY
# =============================================================================

print("\n" + "=" * 90)
print("PAIRWISE TESTS FOR SIGNIFICANT REGIONAL SPEI-48 METRICS")
print("=" * 90)

pairwise_rows = []

significant_metrics = kw_results_df[kw_results_df["significant_p05"]]["metric"].tolist()

if len(significant_metrics) == 0:
    print("\nNo SPEI-48 distribution metric differed significantly among regions at p < 0.05.")
else:
    print(f"\nSignificant metrics: {significant_metrics}")
    
    for metric in significant_metrics:
        raw_pairwise = []
        
        for i, r1 in enumerate(regions):
            for r2 in regions[i+1:]:
                vals1 = df[df["coast_region_analysis"] == r1][metric].dropna().values
                vals2 = df[df["coast_region_analysis"] == r2][metric].dropna().values
                
                if len(vals1) >= 2 and len(vals2) >= 2:
                    u, p = mannwhitneyu(vals1, vals2, alternative="two-sided")
                    
                    raw_pairwise.append({
                        "metric": metric,
                        "region1": r1,
                        "region2": r2,
                        "n1": len(vals1),
                        "n2": len(vals2),
                        "median1": np.median(vals1),
                        "median2": np.median(vals2),
                        "U": u,
                        "p_raw": p
                    })
        
        if len(raw_pairwise) > 0:
            p_values = [row["p_raw"] for row in raw_pairwise]
            reject, p_adj = fdrcorrection(p_values, alpha=0.05)
            
            for row, adj, rej in zip(raw_pairwise, p_adj, reject):
                row["p_fdr"] = adj
                row["significant_fdr"] = rej
                pairwise_rows.append(row)
                
                if rej:
                    print(f"\n{metric}: {row['region1']} vs {row['region2']}")
                    print(f"  median = {row['median1']:.3f} vs {row['median2']:.3f}")
                    print(f"  p_FDR = {adj:.4f}")

pairwise_df = pd.DataFrame(pairwise_rows)
pairwise_csv = OUTPUT_DIR / "regional_spei48_distribution_pairwise_results.csv"
pairwise_df.to_csv(pairwise_csv, index=False)

print(f"\nSaved pairwise results: {pairwise_csv}")

# =============================================================================
# STEP 8: COMPARE REGIONAL SPEI DISTRIBUTION WITH REGIONAL SENSITIVITY
# =============================================================================
# This checks whether regions with wetter/drier SPEI exposure also have different
# median WUE_T sensitivity.

print("\n" + "=" * 90)
print("REGIONAL SPEI EXPOSURE VS REGIONAL WUE_T SENSITIVITY")
print("=" * 90)

regional_sensitivity = (
    df.groupby("coast_region_analysis")
    .agg(
        N_sites=("site_name", "count"),
        median_slope=("slope_theilsen", "median"),
        median_abs_slope=("slope_theilsen", lambda x: np.median(np.abs(x))),
        median_spei48=("spei48_median", "median"),
        median_spei48_range=("spei48_range", "median"),
        median_frac_dry=("frac_dry_neg", "median"),
        median_frac_wet=("frac_wet_pos", "median"),
    )
    .reset_index()
)

print(regional_sensitivity.to_string(index=False))

regional_sensitivity_csv = OUTPUT_DIR / "regional_spei48_exposure_vs_sensitivity_summary.csv"
regional_sensitivity.to_csv(regional_sensitivity_csv, index=False)

print(f"\nSaved regional exposure vs sensitivity table: {regional_sensitivity_csv}")

# =============================================================================
# STEP 9: SAVE FULL SITE-LEVEL TABLE
# =============================================================================

site_output_csv = OUTPUT_DIR / "site_level_spei48_distribution_by_region.csv"
df.to_csv(site_output_csv, index=False)

print("\nSaved site-level diagnostic table:")
print(site_output_csv)

# =============================================================================
# FINAL GUIDE
# =============================================================================

print("\n" + "=" * 90)
print("FINAL INTERPRETATION GUIDE")
print("=" * 90)

print("""
Use this diagnostic to decide wording for Panel B:

1. If regional SPEI-48 distribution summaries are NOT significantly different:
   → Regional sensitivity differences, if present, are not likely explained by unequal SPEI-48 exposure.

2. If some SPEI-48 distribution summaries ARE significantly different:
   → Be cautious. Say that regional sensitivity comparisons were evaluated alongside regional differences
     in SPEI-48 exposure, and avoid saying region alone explains sensitivity.

3. If Panel B regional sensitivity is NOT significant, as in your current SPEI-48 result:
   → You do not need a strong regional-exposure correction in the Results. A short statement is enough:
     "Regional SPEI-48 ranges overlapped broadly, and WUE_T sensitivity did not differ significantly
      among coastal regions."

4. If a reviewer asks for more:
   → Put this diagnostic in supplement/reviewer response.
""")

print("\n" + "=" * 90)
print("DIAGNOSTIC COMPLETE")
print("=" * 90)
