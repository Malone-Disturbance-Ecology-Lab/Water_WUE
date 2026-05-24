# -*- coding: utf-8 -*-
"""
Created on Mon May 18 11:11:03 2026

@author: ammar
"""

# -*- coding: utf-8 -*-
"""
Q3 CLIMATE–BIOME DATASET PREPARATION
SPEI-48 WUE_T sensitivity analysis

Uses EXACT SAME filtering and slope logic as the main Q3 workflow:
    • strict triple intersection
    • Trans_ratio > 0
    • SPEI-48 vs WUE_T
    • minimum 5 monthly observations
    • Theil–Sen slopes
    • linear regression p-values

Goal:
    Prepare climate–biome grouped sensitivity dataset
    BEFORE any inferential testing.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
from sklearn.linear_model import TheilSenRegressor
import warnings
warnings.filterwarnings("ignore")

# =============================================================================
# FILE PATHS
# =============================================================================

BASE_DIR = Path(
    r"\\corellia.environment.yale.edu\MaloneLab\Research\WUE_CUE\data_products\results"
)

MONTHLY_FILE = BASE_DIR / "monthly_data_after_outlier_removal.csv"

NN_MEDIANS_FILE = BASE_DIR / "wue_site_level_NN_medians_SPEI_1.csv"

CLIMATE_BIOME_FILE = BASE_DIR / "ClimateBiome_WUE_WUET_FinalDataset.csv"

OUTPUT_DIR = BASE_DIR / "Q3_analysis" / "ClimateBiome_Q3"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# LOAD DATA
# =============================================================================

print("=" * 90)
print("Q3 CLIMATE–BIOME DATASET PREPARATION")
print("SPEI-48 WUE_T sensitivity analysis")
print("=" * 90)

df_monthly = pd.read_csv(MONTHLY_FILE)
nn_data = pd.read_csv(NN_MEDIANS_FILE)
clim_df = pd.read_csv(CLIMATE_BIOME_FILE)

print(f"\nLoaded monthly data: {len(df_monthly):,} rows")
print(f"Loaded NN medians file: {len(nn_data):,} rows")
print(f"Loaded climate–biome file: {len(clim_df):,} rows")

# =============================================================================
# STEP 1: STRICT TRIPLE INTERSECTION
# =============================================================================

print("\n" + "=" * 90)
print("STEP 1: STRICT TRIPLE INTERSECTION")
print("=" * 90)

wue_sites = set(
    nn_data[nn_data["WUE_Metric"] == "WUE"]["site_name"]
    .dropna()
    .unique()
)

eva_sites = set(
    nn_data[nn_data["WUE_Metric"] == "WUE_eva"]["site_name"]
    .dropna()
    .unique()
)

tra_sites = set(
    nn_data[nn_data["WUE_Metric"] == "WUE_tra"]["site_name"]
    .dropna()
    .unique()
)

shared_sites = (
    wue_sites
    .intersection(eva_sites)
    .intersection(tra_sites)
)

df_monthly = df_monthly[
    df_monthly["site_name"].isin(shared_sites)
].copy()

print(f"\nShared sites retained: {len(shared_sites)}")
print(f"Rows after triple intersection: {len(df_monthly):,}")

# =============================================================================
# STEP 2: STRICT MONTH FILTER
# =============================================================================

print("\n" + "=" * 90)
print("STEP 2: STRICT MONTH FILTER")
print("=" * 90)

strict_mask = (
    df_monthly["WUE"].notna() &
    df_monthly["WUE_eva"].notna() &
    df_monthly["WUE_tra"].notna()
)

before_rows = len(df_monthly)

df_monthly = df_monthly[strict_mask].copy()

after_rows = len(df_monthly)

print(f"\nRows removed: {before_rows - after_rows}")
print(f"Rows retained: {after_rows:,}")

# =============================================================================
# STEP 3: TRANS_RATIO FILTER
# =============================================================================

print("\n" + "=" * 90)
print("STEP 3: TRANS_RATIO FILTER")
print("=" * 90)

before_ratio = len(df_monthly)

df_monthly = df_monthly[
    df_monthly["Trans_ratio"] > 0
].copy()

after_ratio = len(df_monthly)

print(f"\nRows removed: {before_ratio - after_ratio}")
print(f"Rows retained: {after_ratio:,}")
print(f"Sites retained: {df_monthly['site_name'].nunique()}")

# =============================================================================
# STEP 4: COMPUTE SPEI-48 SENSITIVITY SLOPES
# =============================================================================

print("\n" + "=" * 90)
print("STEP 4: COMPUTE SITE-LEVEL SPEI-48 SENSITIVITY")
print("=" * 90)

site_results = []

for site, g in df_monthly.groupby("site_name"):

    g = g.dropna(subset=["SPEI_48", "WUE_tra"]).copy()

    n_months = len(g)

    if n_months >= 5:

        X = g["SPEI_48"].values.reshape(-1, 1)
        y = g["WUE_tra"].values

        # -------------------------------------------------------------
        # THEIL-SEN SLOPE
        # -------------------------------------------------------------

        try:
            ts = TheilSenRegressor(random_state=42)
            ts.fit(X, y)

            slope = ts.coef_[0]

        except:
            slope = np.nan

        # -------------------------------------------------------------
        # LINEAR REGRESSION P-VALUE
        # -------------------------------------------------------------

        try:
            _, _, _, p_value, _ = stats.linregress(
                g["SPEI_48"].values,
                y
            )

        except:
            p_value = np.nan

        # -------------------------------------------------------------
        # SIGN FLAGS
        # -------------------------------------------------------------

        if slope > 0:
            slope_sign = "Positive"

        elif slope < 0:
            slope_sign = "Negative"

        else:
            slope_sign = "Zero"

        # -------------------------------------------------------------
        # SAVE
        # -------------------------------------------------------------

        site_results.append({

            "site_name": site,

            "slope_theilsen": slope,

            "p_value": p_value,

            "is_significant": (
                (p_value < 0.05)
                if pd.notna(p_value)
                else False
            ),

            "slope_sign": slope_sign,

            "n_months": n_months,

            "spei48_min": g["SPEI_48"].min(),

            "spei48_max": g["SPEI_48"].max(),

            "spei48_median": g["SPEI_48"].median(),

            "spei48_range": (
                g["SPEI_48"].max() -
                g["SPEI_48"].min()
            )
        })

df_q3 = pd.DataFrame(site_results)

print(f"\nSites with valid slopes: {len(df_q3)}")

print(f"\nPositive slopes: {(df_q3['slope_theilsen'] > 0).sum()}")
print(f"Negative slopes: {(df_q3['slope_theilsen'] < 0).sum()}")

print(f"\nSignificant sites (p < 0.05): "
      f"{df_q3['is_significant'].sum()} / {len(df_q3)} "
      f"({df_q3['is_significant'].mean()*100:.1f}%)")

# =============================================================================
# STEP 5: MERGE CLIMATE–BIOME GROUPS
# =============================================================================

print("\n" + "=" * 90)
print("STEP 5: MERGE CLIMATE–BIOME GROUPS")
print("=" * 90)

group_cols = [
    "site_name",
    "Final_Group",
    "Broad_Biome",
    "climate_clean"
]

clim_unique = clim_df[group_cols].drop_duplicates(
    subset=["site_name"]
)

df_q3 = df_q3.merge(
    clim_unique,
    on="site_name",
    how="left"
)

# =============================================================================
# CHECK UNASSIGNED
# =============================================================================

missing_group = df_q3[
    df_q3["Final_Group"].isna()
]

print(f"\nSites without climate–biome group: {len(missing_group)}")

if len(missing_group) > 0:

    for s in missing_group["site_name"].tolist():
        print(f"  - {s}")

# =============================================================================
# REMOVE UNASSIGNED SITES
# =============================================================================

df_q3 = df_q3[
    df_q3["Final_Group"].notna()
].copy()

# =============================================================================
# STEP 6: SUPPORT SUMMARY
# =============================================================================

print("\n" + "=" * 90)
print("STEP 6: SUPPORT SUMMARY BY CLIMATE–BIOME GROUP")
print("=" * 90)

summary_rows = []

for grp, g in df_q3.groupby("Final_Group"):

    summary_rows.append({

        "Final_Group": grp,

        "N_sites": len(g),

        "Positive_slopes": (
            g["slope_theilsen"] > 0
        ).sum(),

        "Negative_slopes": (
            g["slope_theilsen"] < 0
        ).sum(),

        "Significant_sites": (
            g["is_significant"]
        ).sum(),

        "Median_slope": g["slope_theilsen"].median(),

        "Slope_IQR": (
            g["slope_theilsen"].quantile(0.75) -
            g["slope_theilsen"].quantile(0.25)
        ),

        "Slope_min": g["slope_theilsen"].min(),

        "Slope_max": g["slope_theilsen"].max(),

        "Median_n_months": g["n_months"].median()
    })

summary_df = pd.DataFrame(summary_rows)

summary_df = summary_df.sort_values(
    "Median_slope"
)

print(summary_df.to_string(index=False))

# =============================================================================
# STEP 7: SAVE OUTPUTS
# =============================================================================

print("\n" + "=" * 90)
print("STEP 7: SAVE OUTPUTS")
print("=" * 90)

site_output = OUTPUT_DIR / "Q3_ClimateBiome_SPEI48_sitelevel.csv"

summary_output = OUTPUT_DIR / "Q3_ClimateBiome_SPEI48_group_summary.csv"

df_q3.to_csv(site_output, index=False)

summary_df.to_csv(summary_output, index=False)

print(f"\nSaved site-level dataset:")
print(site_output)

print(f"\nSaved group summary:")
print(summary_output)

# =============================================================================
# FINAL SUMMARY
# =============================================================================

print("\n" + "=" * 90)
print("FINAL SUMMARY")
print("=" * 90)

print(f"\nFinal sites retained: {len(df_q3)}")

print(f"\nClimate–biome groups:")

for grp in sorted(df_q3["Final_Group"].unique()):

    n = (df_q3["Final_Group"] == grp).sum()

    print(f"  {grp}: {n}")

print(f"\nOverall slope distribution:")
print(f"  Positive: {(df_q3['slope_theilsen'] > 0).sum()}")
print(f"  Negative: {(df_q3['slope_theilsen'] < 0).sum()}")

print(f"\nOverall significant sites:")
print(f"  {(df_q3['is_significant']).sum()} / {len(df_q3)}")

print("\n✓ Q3 climate–biome dataset preparation complete")