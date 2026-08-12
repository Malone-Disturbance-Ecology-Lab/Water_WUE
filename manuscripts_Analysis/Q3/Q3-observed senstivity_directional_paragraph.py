# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 14:21:26 2026

@author: ammar
"""

#!/usr/bin/env python3
"""
Q3_DIAGNOSTIC_SIGNED_DRY_RESPONSE.py

READ ONLY.
Does not save or overwrite anything.

Purpose:
Check observed signed WUE_T responses during drought by coastline
and SPEI timescale.

Positive = WUE_T increased relative to near-normal baseline.
Negative = WUE_T decreased relative to near-normal baseline.
"""

import os
import numpy as np
import pandas as pd

# ============================================================
# PATH
# ============================================================

base_dir = r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3"

input_file = os.path.join(
    base_dir,
    "Q3_WUE_T_SPEI_sensitivity_outputs",
    "Q3_WUE_T_SPEI_site_response_from_near_normal.csv"
)

# ============================================================
# SETTINGS
# ============================================================

ALL_TIMESCALES = [
    "SPEI_1",
    "SPEI_3",
    "SPEI_6",
    "SPEI_12",
    "SPEI_24",
    "SPEI_36",
    "SPEI_48"
]

DISPLAY_TIMESCALES = ["SPEI_1", "SPEI_3", "SPEI_48"]

# Require at least 3 drought observations for a site × timescale
MIN_DRY_MONTHS = 3

# ============================================================
# LOAD
# ============================================================

print("=" * 90)
print("OBSERVED SIGNED WUE_T RESPONSE DURING SPEI DROUGHT")
print("=" * 90)

df = pd.read_csv(input_file)

required = [
    "site_name",
    "coast_region",
    "water_class",
    "SPEI_timescale",
    "SPEI_value",
    "WUE_T_pct_change"
]

missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(
        f"Missing required columns: {missing}\n"
        f"Available columns: {list(df.columns)}"
    )

print(f"\nInput rows: {len(df):,}")
print(f"Sites: {df['site_name'].nunique()}")

# ============================================================
# KEEP DRY OBSERVATIONS
# ============================================================

dry = df[
    (df["SPEI_timescale"].isin(ALL_TIMESCALES)) &
    (df["SPEI_value"] <= -1) &
    np.isfinite(df["WUE_T_pct_change"])
].copy()

print(f"Dry site-month-timescale rows: {len(dry):,}")

# ============================================================
# CALCULATE SIGNED RESPONSE FOR EACH SITE × TIMESCALE
# ============================================================

site_signed = (
    dry
    .groupby(
        ["site_name", "coast_region", "water_class", "SPEI_timescale"],
        as_index=False
    )
    .agg(
        n_dry_months=("WUE_T_pct_change", "size"),
        mean_signed_dry_pct_change=("WUE_T_pct_change", "mean")
    )
)

site_signed = site_signed[
    site_signed["n_dry_months"] >= MIN_DRY_MONTHS
].copy()

print(
    f"Site × timescale estimates with >= {MIN_DRY_MONTHS} "
    f"dry months: {len(site_signed):,}"
)

# ============================================================
# COAST × TIMESCALE SUMMARY
# ============================================================

summary = (
    site_signed
    .groupby(["coast_region", "SPEI_timescale"])
    ["mean_signed_dry_pct_change"]
    .agg(
        n_sites="count",
        mean="mean",
        median="median",
        q25=lambda x: x.quantile(0.25),
        q75=lambda x: x.quantile(0.75),
        min="min",
        max="max",
        pct_positive=lambda x: (x > 0).mean() * 100,
        pct_negative=lambda x: (x < 0).mean() * 100
    )
    .reset_index()
)

# ============================================================
# PRINT ALL SEVEN TIMESCALES
# ============================================================

print("\n" + "=" * 90)
print("ALL SEVEN SPEI TIMESCALES")
print("Positive = observed WUE_T increase during drought")
print("Negative = observed WUE_T decrease during drought")
print("=" * 90)

for ts in ALL_TIMESCALES:

    print(f"\n{ts}")
    print("-" * 90)

    sub = summary[summary["SPEI_timescale"] == ts].copy()

    if sub.empty:
        print("No qualifying sites.")
        continue

    print(
        sub[
            [
                "coast_region",
                "n_sites",
                "mean",
                "median",
                "q25",
                "q75",
                "pct_positive",
                "pct_negative"
            ]
        ]
        .round(1)
        .to_string(index=False)
    )

# ============================================================
# MANUSCRIPT TIMESCALES ONLY
# ============================================================

print("\n" + "=" * 90)
print("SELECTED TIMESCALES FOR FIGURE/TABLE: SPEI-1, SPEI-3, SPEI-48")
print("=" * 90)

selected = summary[
    summary["SPEI_timescale"].isin(DISPLAY_TIMESCALES)
].copy()

print(
    selected[
        [
            "coast_region",
            "SPEI_timescale",
            "n_sites",
            "mean",
            "median",
            "q25",
            "q75",
            "pct_positive",
            "pct_negative"
        ]
    ]
    .round(1)
    .sort_values(["SPEI_timescale", "coast_region"])
    .to_string(index=False)
)

# ============================================================
# SIMPLE SIGN CHECK
# ============================================================

print("\n" + "=" * 90)
print("DOMINANT OBSERVED DIRECTION")
print("=" * 90)

for _, row in selected.sort_values(
    ["SPEI_timescale", "coast_region"]
).iterrows():

    if row["median"] > 0:
        direction = "INCREASE"
    elif row["median"] < 0:
        direction = "DECREASE"
    else:
        direction = "NEUTRAL"

    print(
        f"{row['coast_region']:15s} "
        f"{row['SPEI_timescale']:8s}: "
        f"{direction:8s} | "
        f"median={row['median']:.1f}% | "
        f"mean={row['mean']:.1f}% | "
        f"{row['pct_positive']:.1f}% sites positive, "
        f"{row['pct_negative']:.1f}% sites negative "
        f"(n={int(row['n_sites'])})"
    )

print("\nDONE. No files were written.")