# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 09:36:05 2026

@author: ammar
"""

# =============================================================================
# Q1: T:ET ratio dataset by coast and water class under near-normal SPEI-1
#
# Saves only two clean data CSVs:
#   1. Q1_NN_TET_monthly_by_coast_waterclass.csv
#      - all near-normal monthly T:ET observations
#      - useful for density plots
#
#   2. Q1_NN_TET_site_level_by_coast_waterclass.csv
#      - one row per site
#      - useful for boxplots, jitter plots, and site-level tests
#
# Console only:
#   - summaries by coast
#   - summaries by water_class
#   - summaries by coast × water_class
#   - tests for T:ET differences by coast
#   - tests for T:ET differences by water_class
#
# Does NOT write anything to Q3.
# =============================================================================

import os
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

monthly_data_path = r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\WUE_CUE_monthly_merged_indices_clean.csv"

output_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q1\Q1_updated_results\outputs\T_ET_ratio_coast"
os.makedirs(output_dir, exist_ok=True)

monthly_output_file = os.path.join(
    output_dir,
    "Q1_NN_TET_monthly_by_coast_waterclass.csv"
)

site_output_file = os.path.join(
    output_dir,
    "Q1_NN_TET_site_level_by_coast_waterclass.csv"
)

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

ECOSYSTEM_CLASSES = ["Upland", "Freshwater", "Saline"]
COAST_REGION_LEVELS = ["Atlantic Coast", "Pacific Coast", "Gulf Coast", "AK Coast"]

def classify_coast_region(lat, long):
    """Same coast-region rule used in Q3/Q4."""
    if lat > 50:
        return "AK Coast"
    elif long > -100:
        return "Atlantic Coast"
    elif long < -120:
        return "Pacific Coast"
    else:
        return "Gulf Coast"

def bh_adjust(pvals):
    """Benjamini-Hochberg FDR adjustment."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adjusted = np.empty(n, dtype=float)

    running_min = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        value = ranked[i] * n / rank
        running_min = min(running_min, value)
        adjusted[i] = running_min

    out = np.empty(n, dtype=float)
    out[order] = np.minimum(adjusted, 1.0)
    return out

def print_summary(site_level, group_cols, title):
    """Print site-level T:ET summary."""
    summary = (
        site_level
        .groupby(group_cols, observed=True)
        .agg(
            n_sites=("site_name", "nunique"),
            n_months=("n_months", "sum"),
            mean_TET=("mean_TET", "mean"),
            median_TET=("mean_TET", "median"),
            sd_TET=("mean_TET", "std"),
            min_TET=("mean_TET", "min"),
            max_TET=("mean_TET", "max")
        )
        .reset_index()
    )

    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)
    print(summary.to_string(index=False))

def test_group_difference(site_level, group_col, value_col="mean_TET", context_label=""):
    """
    Kruskal-Wallis + pairwise Mann-Whitney tests using site-level mean T:ET.
    Tests are printed only, not saved.
    """
    test_df = site_level[[group_col, value_col]].dropna().copy()

    group_counts = test_df.groupby(group_col, observed=True)[value_col].size()
    valid_groups = [g for g, n in group_counts.items() if n >= 2]

    print("\n" + "=" * 80)
    print(f"TEST: {context_label}")
    print("=" * 80)
    print("Group counts:")
    print(group_counts.to_string())

    if len(valid_groups) < 2:
        print("\nNot enough groups with at least 2 sites.")
        return

    groups = [
        test_df.loc[test_df[group_col] == g, value_col].values
        for g in valid_groups
    ]

    h_stat, p_val = kruskal(*groups)

    print("\nKruskal-Wallis overall test:")
    print(f"  statistic = {h_stat:.4f}")
    print(f"  p-value   = {p_val:.6g}")

    pair_rows = []
    pvals = []

    for i, g1 in enumerate(valid_groups):
        for g2 in valid_groups[i + 1:]:
            x = test_df.loc[test_df[group_col] == g1, value_col].values
            y = test_df.loc[test_df[group_col] == g2, value_col].values

            if len(x) >= 2 and len(y) >= 2:
                stat, p_pair = mannwhitneyu(x, y, alternative="two-sided")
                pvals.append(p_pair)
                pair_rows.append({
                    "comparison": f"{g1} vs {g2}",
                    "n_1": len(x),
                    "n_2": len(y),
                    "mean_1": np.mean(x),
                    "mean_2": np.mean(y),
                    "difference_1_minus_2": np.mean(x) - np.mean(y),
                    "p_value": p_pair
                })

    if pair_rows:
        p_adj = bh_adjust(pvals)

        for row, adj in zip(pair_rows, p_adj):
            row["p_value_BH"] = adj
            row["significant_BH_0.05"] = adj < 0.05

        pair_df = pd.DataFrame(pair_rows)

        print("\nPairwise tests:")
        print(pair_df.to_string(index=False))

# -----------------------------------------------------------------------------
# Load monthly data
# -----------------------------------------------------------------------------

monthly = pd.read_csv(monthly_data_path)

for col in monthly.select_dtypes(include="object").columns:
    monthly[col] = monthly[col].where(
        monthly[col].isna(),
        monthly[col].astype(str).str.strip()
    )
    monthly[col] = monthly[col].replace({"": np.nan, "nan": np.nan, "NaN": np.nan})

# -----------------------------------------------------------------------------
# Q1-style filter
# -----------------------------------------------------------------------------
# Q1 requires all three WUE metrics plus Trans_ratio and SPEI_1.
# This keeps the T:ET dataset aligned with Q1.

required_cols = [
    "site_name", "Year", "month", "water_class", "lat", "long",
    "Trans_ratio", "WUE", "WUE_tra", "WUE_eva", "SPEI_1"
]

missing_cols = [c for c in required_cols if c not in monthly.columns]
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

q1_base = monthly.dropna(subset=required_cols).copy()

q1_base = q1_base[
    (q1_base["site_name"] != "") &
    (q1_base["water_class"].isin(ECOSYSTEM_CLASSES)) &
    np.isfinite(q1_base["lat"]) &
    np.isfinite(q1_base["long"]) &
    np.isfinite(q1_base["Trans_ratio"]) &
    (q1_base["Trans_ratio"] >= 0) &
    (q1_base["Trans_ratio"] <= 1) &
    np.isfinite(q1_base["WUE"]) &
    np.isfinite(q1_base["WUE_tra"]) &
    np.isfinite(q1_base["WUE_eva"]) &
    np.isfinite(q1_base["SPEI_1"])
].copy()

q1_base["coast_region"] = q1_base.apply(
    lambda row: classify_coast_region(row["lat"], row["long"]),
    axis=1
)

q1_base["coast_region"] = pd.Categorical(
    q1_base["coast_region"],
    categories=COAST_REGION_LEVELS,
    ordered=True
)

q1_base["water_class"] = pd.Categorical(
    q1_base["water_class"],
    categories=ECOSYSTEM_CLASSES,
    ordered=True
)

# -----------------------------------------------------------------------------
# Near-normal SPEI-1 filter
# -----------------------------------------------------------------------------

q1_nn = q1_base[
    (q1_base["SPEI_1"] >= -1) &
    (q1_base["SPEI_1"] <= 1)
].copy()

# -----------------------------------------------------------------------------
# Monthly-level dataset
# -----------------------------------------------------------------------------
# Use this for density plots.
# Each row is a site-month observation.

monthly_level = q1_nn[
    [
        "site_name", "Year", "month",
        "coast_region", "water_class",
        "Trans_ratio", "SPEI_1",
        "lat", "long"
    ]
].copy()

monthly_level = monthly_level.rename(columns={"Trans_ratio": "TET_ratio"})

# -----------------------------------------------------------------------------
# Site-level dataset
# -----------------------------------------------------------------------------
# Use this for main coast/water-class comparison.
# Each site contributes one mean T:ET value.

site_level = (
    q1_nn
    .groupby(["site_name", "coast_region", "water_class"], observed=True)
    .agg(
        n_months=("Trans_ratio", "count"),
        mean_TET=("Trans_ratio", "mean"),
        median_TET=("Trans_ratio", "median"),
        sd_TET=("Trans_ratio", "std"),
        min_TET=("Trans_ratio", "min"),
        max_TET=("Trans_ratio", "max"),
        mean_SPEI_1=("SPEI_1", "mean"),
        min_SPEI_1=("SPEI_1", "min"),
        max_SPEI_1=("SPEI_1", "max"),
        lat=("lat", "first"),
        long=("long", "first")
    )
    .reset_index()
)

# -----------------------------------------------------------------------------
# Save only data CSVs
# -----------------------------------------------------------------------------

monthly_level.to_csv(monthly_output_file, index=False)
site_level.to_csv(site_output_file, index=False)

# -----------------------------------------------------------------------------
# Console checks and summaries
# -----------------------------------------------------------------------------

print("=" * 80)
print("Q1 near-normal SPEI-1 T:ET datasets created")
print("=" * 80)

print(f"Q1-style base rows: {len(q1_base)}")
print(f"Q1 near-normal SPEI-1 rows: {len(q1_nn)}")
print(f"Sites in Q1 NN dataset: {q1_nn['site_name'].nunique()}")

print("\nSaved data CSVs:")
print(monthly_output_file)
print(site_output_file)

print_summary(
    site_level=site_level,
    group_cols=["coast_region"],
    title="Site-level T:ET summary by coast"
)

print_summary(
    site_level=site_level,
    group_cols=["water_class"],
    title="Site-level T:ET summary by water class"
)

print_summary(
    site_level=site_level,
    group_cols=["coast_region", "water_class"],
    title="Site-level T:ET summary by coast × water class"
)

# -----------------------------------------------------------------------------
# Tests printed only
# -----------------------------------------------------------------------------

test_group_difference(
    site_level=site_level,
    group_col="coast_region",
    value_col="mean_TET",
    context_label="Does site-level T:ET differ by coast?"
)

test_group_difference(
    site_level=site_level,
    group_col="water_class",
    value_col="mean_TET",
    context_label="Does site-level T:ET differ by water class?"
)

# Optional: water-class differences within each coast
for coast in COAST_REGION_LEVELS:
    sub = site_level[site_level["coast_region"] == coast].copy()
    if sub["water_class"].nunique() >= 2:
        test_group_difference(
            site_level=sub,
            group_col="water_class",
            value_col="mean_TET",
            context_label=f"Does site-level T:ET differ by water class within {coast}?"
        )

# Optional: coast differences within each water class
for wc in ECOSYSTEM_CLASSES:
    sub = site_level[site_level["water_class"] == wc].copy()
    if sub["coast_region"].nunique() >= 2:
        test_group_difference(
            site_level=sub,
            group_col="coast_region",
            value_col="mean_TET",
            context_label=f"Does site-level T:ET differ by coast within {wc}?"
        )

print("\nDone. Only two data CSVs were saved. No Q3 files were modified.")