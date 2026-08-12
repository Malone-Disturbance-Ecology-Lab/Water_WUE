# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 13:50:29 2026

@author: ammar
"""

import pandas as pd

f = (
    r"M:\Research\WUE_CUE\WUE_manuscript_version6\Q3"
    r"\Q3_WUE_T_SPEI_sensitivity_outputs"
    r"\Q3_WUE_T_SPEI_model_data_long.csv"
)

d = pd.read_csv(f)

# Use SPEI-3 only so each site-month appears once
p = d[
    (d["coast_region"] == "Pacific Coast") &
    (d["SPEI_timescale"] == "SPEI_3")
].copy()

site = (
    p.groupby(["water_class", "site_name"], observed=True)
    .agg(
        n_months=("SPEI_value", "size"),
        min_SPEI3=("SPEI_value", "min"),
        max_SPEI3=("SPEI_value", "max"),
        n_dry_months=("SPEI_value", lambda x: (x <= -1).sum())
    )
    .reset_index()
)

print("=" * 80)
print("PACIFIC COAST OBSERVED SITE SUPPORT")
print("=" * 80)

summary = (
    site.groupby("water_class", observed=True)
    .agg(
        n_sites=("site_name", "nunique"),
        sites_with_drought=("n_dry_months", lambda x: (x > 0).sum()),
        total_dry_months=("n_dry_months", "sum"),
        minimum_observed_SPEI3=("min_SPEI3", "min")
    )
    .reset_index()
)

print(summary.to_string(index=False))

print("\nPACIFIC SALINE SITES:")
print(
    site[site["water_class"] == "Saline"]
    .sort_values("site_name")
    .to_string(index=False)
)

print("\nPACIFIC FRESHWATER SITES:")
print(
    site[site["water_class"] == "Freshwater"]
    .sort_values("site_name")
    .to_string(index=False)
)

print("\nDONE. Nothing was saved.")